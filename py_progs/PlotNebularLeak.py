#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Interactive Plotly visualization of the DecomposeCleanSky.py /
    sky_nebular_leak_eval.py nebular-leak workflow for one exposure: the
    shared-systematic residual pattern across FLUX/SKY_EAST/SKY_WEST in a
    chosen wavelength band, plus zoomed before/after panels at each
    requested nebular line showing the raw residual, the scaled template
    being subtracted, the corrected residual, and the fitted Gaussian.

Command line usage (if any):

    usage: PlotNebularLeak.py [-h] [-lines LIST] [-targets LIST]
                              [-template_ext EXT] [-template_band LO,HI]
                              [-leak_file PATH] [-title TITLE]
                              [-outfile PATH]
                              decomp_file

    where

    decomp_file     a DecomposeCleanSky.py output FITS file (WAVE,
                    <EXT>, <EXT>_BESTFIT, <EXT>_RESID extensions).

    -lines LIST     comma-separated NEBULAR_LINES names to show zoomed
                    panels for (default: siii_a,siii_b).

    -targets LIST   comma-separated extensions to show zoom panels for
                    (default: FLUX,SKY_WEST).

    -template_ext EXT
                    extension whose RESID is treated as the shared-
                    systematic template (default: SKY_EAST) -- must
                    match what -leak_file (if given) was produced with,
                    or the overlaid scaled-template/corrected curves
                    won't match that table's numbers. Recomputed locally
                    with sky_nebular_leak_eval.fit_template_scale, not
                    read from -leak_file, so this plot works even
                    without one.

    -template_band LO,HI
                    wavelength range (Angstrom) for the top overview
                    panel and the template scale fit (default: 9000,9600
                    -- matches sky_nebular_leak_eval.py's -template_band
                    default, and the flag name matches it too).

    -leak_file PATH sky_nebular_leak_eval.py output (<root>_lines.fits).
                    If given, the fitted Gaussian curves in the zoom
                    panels are drawn from its LAM0_FIT/SIGMA0_FIT/
                    LEAK_FLUX columns (matching exactly what that table
                    reports); otherwise this script fits them itself
                    with the same fit_leak_line routine.

    -title TITLE    plot title (default: the input filename).

    -outfile PATH   output HTML path (default:
                    Overview_Plot/<stem>.nebleak.html).

Description:

    Top panel: <EXT>_RESID for every extension present in decomp_file,
    overlaid over -band, sharing one y-axis -- this is the panel that
    shows the shared-systematic correlation directly (same shape,
    different amplitude across extensions within one exposure).

    One row of zoom panels per -lines entry, one column per -targets
    extension: raw residual (black), the scaled template
    k*template_ext_RESID being subtracted (grey dashed), the corrected
    residual (firebrick), and the fitted Gaussian on the corrected
    residual (green), each panel sharing a wavelength x-axis around that
    line's Doppler-shifted (vel=0 assumed -- see sky_nebular_leak_eval.py
    for the -v/-lmc/-smc convention this does not currently expose)
    rest wavelength.

Primary routines:

    plot_nebular_leak   build the full figure for one decomp_file.

Notes::

    Reuses sky_gaussfit.resolve_nebular_lines and
    sky_nebular_leak_eval.{fit_leak_line, fit_template_scale,
    line_exclusion_mask} rather than reimplementing this logic --
    matches sky_nebular_leak_eval.py's numbers by construction.

History::

    260909  ksl  Coding begun.

'''

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from sky_gaussfit import NEBULAR_LINES, resolve_nebular_lines
from sky_nebular_leak_eval import fit_leak_line, fit_template_scale, line_exclusion_mask, _gaussian

TRACE_WIDTH = 2.0
FIT_WIDTH = 2.5
SUBPLOT_TITLE_SIZE = 15
ROW_HEIGHT = 380


def _load_resid(hdul, ext):
    if f'{ext}_RESID' in hdul:
        return np.asarray(hdul[f'{ext}_RESID'].data, dtype=float)
    if ext in hdul and f'{ext}_BESTFIT' in hdul:
        return (np.asarray(hdul[ext].data, dtype=float)
               - np.asarray(hdul[f'{ext}_BESTFIT'].data, dtype=float))
    return None


def plot_nebular_leak(decomp_file, lines=('siii_a', 'siii_b'), targets=('FLUX', 'SKY_WEST'),
                      template_ext='SKY_EAST', band=(9000., 9600.), leak_file=None,
                      title=None, vel=0.0):
    '''
    Build the overview + per-line zoom figure for one DecomposeCleanSky.py
    output file.

    Parameters
    ----------
    decomp_file : str or Path
    lines : sequence of str
        NEBULAR_LINES names to show zoom panels for.
    targets : sequence of str
        Extensions to show zoom panels for (columns).
    template_ext : str
        Extension used as the shared-systematic template.
    band : (float, float)
    leak_file : str or Path, optional
        sky_nebular_leak_eval.py output table; if given, fitted-Gaussian
        curves are drawn from its columns instead of being refit here.
    title : str, optional
    vel : float
        Nebular systemic velocity (km/s), same convention as
        sky_nebular_leak_eval.py -- must match what decomp_file and
        leak_file (if given) were produced with.

    Returns
    -------
    plotly.graph_objects.Figure
    '''
    line_lookup = {name: (name, center, wmin, wmax) for name, center, wmin, wmax in NEBULAR_LINES}
    resolved, dropped = resolve_nebular_lines(vel)
    resolved_names = {name for name, *_ in resolved}
    lines = [ln for ln in lines if ln in line_lookup]

    leak_tab = Table.read(leak_file) if leak_file else None

    with fits.open(decomp_file) as hdul:
        wave = np.asarray(hdul['WAVE'].data, dtype=float)
        ext_names = [n for n in ('FLUX', 'SKY_EAST', 'SKY_WEST')
                    if f'{n}_RESID' in hdul or (n in hdul and f'{n}_BESTFIT' in hdul)]
        resids = {ext: _load_resid(hdul, ext) for ext in ext_names}

    band_sel = (wave >= band[0]) & (wave <= band[1])
    excl = line_exclusion_mask(wave, resolved, vel)
    template_resid = resids.get(template_ext)

    # Per-target template scale/correlation and corrected residual, computed
    # ONCE over the whole band (not per line -- the fit only depends on the
    # target extension, so refitting it per line would just repeat the same
    # answer) so the overview-AFTER panel and every per-line "after" panel
    # use an identical correction.
    target_k, target_r, target_corrected = {}, {}, {}
    for tgt in ext_names:
        target_corrected[tgt] = resids[tgt]
        target_k[tgt] = target_r[tgt] = np.nan
        if template_resid is not None and tgt != template_ext:
            k, r = fit_template_scale(template_resid[band_sel], resids[tgt][band_sel],
                                      exclude=excl[band_sel])
            target_k[tgt], target_r[tgt] = k, r
            if np.isfinite(k):
                corrected = resids[tgt].copy()
                corrected[band_sel] = resids[tgt][band_sel] - k * template_resid[band_sel]
                target_corrected[tgt] = corrected

    n_target_cols = max(len(targets), 1)
    n_cols = 2 * n_target_cols
    n_rows = 2 + len(lines)
    specs = [[{'colspan': n_cols}] + [None] * (n_cols - 1),
            [{'colspan': n_cols}] + [None] * (n_cols - 1)] + [[{}] * n_cols for _ in lines]
    subplot_titles = [f'Overview BEFORE, {band[0]:.0f}-{band[1]:.0f} A (raw residuals)',
                      f'Overview AFTER, {band[0]:.0f}-{band[1]:.0f} A '
                      f'(template-corrected where applicable)']
    for ln in lines:
        for tgt in targets:
            k_str = f' (k={target_k.get(tgt, np.nan):.2f})' if np.isfinite(target_k.get(tgt, np.nan)) else ''
            subplot_titles.append(f'{ln} ({line_lookup[ln][1]:.1f} A) -- {tgt} BEFORE')
            subplot_titles.append(f'{ln} ({line_lookup[ln][1]:.1f} A) -- {tgt} AFTER{k_str}')

    # 0.09 gives generous room between rows for a few lines; clamp so
    # vertical_spacing*(n_rows-1) never approaches/exceeds 1 (plotly's hard
    # limit) if a user asks for many -lines at once.
    vspace = min(0.09, 0.9 / max(n_rows - 1, 1))
    fig = make_subplots(rows=n_rows, cols=n_cols, specs=specs,
                        subplot_titles=subplot_titles, vertical_spacing=vspace)

    colors = {'FLUX': 'black', 'SKY_EAST': 'steelblue', 'SKY_WEST': 'firebrick'}
    for ext in ext_names:
        fig.add_trace(go.Scatter(x=wave[band_sel], y=resids[ext][band_sel], mode='lines',
                                 name=f'{ext} raw', line=dict(color=colors.get(ext, 'grey'), width=TRACE_WIDTH)),
                      row=1, col=1)
        fig.add_trace(go.Scatter(x=wave[band_sel], y=target_corrected[ext][band_sel], mode='lines',
                                 name=f'{ext} (template-corrected)' if ext != template_ext else f'{ext} (template)',
                                 line=dict(color=colors.get(ext, 'grey'), width=TRACE_WIDTH)),
                      row=2, col=1)
    fig.add_hline(y=0, line=dict(color='grey', width=1, dash='dot'), row=1, col=1)
    fig.add_hline(y=0, line=dict(color='grey', width=1, dash='dot'), row=2, col=1)

    for i, ln in enumerate(lines):
        name, center, wmin, wmax = line_lookup[ln]
        zz = 1.0 + vel / 3e5
        wave0 = zz * center
        wlo, whi = zz * wmin, zz * wmax
        sel = (wave >= wlo) & (wave <= whi)
        ww = wave[sel]
        row = i + 3

        for j, tgt in enumerate(targets):
            col_before, col_after = 2 * j + 1, 2 * j + 2
            target_resid = resids.get(tgt)
            if target_resid is None:
                continue
            k, corrected = target_k[tgt], target_corrected[tgt]
            has_correction = template_resid is not None and tgt != template_ext and np.isfinite(k)
            in_band = band[0] <= wave0 <= band[1]

            # Fitted Gaussian curve: from leak_file if given (matches its numbers
            # exactly), else refit here on whichever residual (corrected if a
            # template was applied and this line is in-band, else raw) sky_
            # nebular_leak_eval.py itself would have used.
            fit_row = None
            if leak_tab is not None:
                sel_tab = (leak_tab['EXT'] == tgt) & (leak_tab['LINE_NAME'] == ln)
                if sel_tab.sum():
                    fit_row = leak_tab[sel_tab][0]
                    amp, lam0, sigma0 = fit_row['LEAK_FLUX'], fit_row['LAM0_FIT'], fit_row['SIGMA0_FIT']
                    bkg = 0.0
            if fit_row is None:
                use = corrected if (has_correction and in_band) else target_resid
                fitres = fit_leak_line(wave, use, wave0, wlo, whi)
                if fitres is not None:
                    amp, lam0, sigma0, bkg = fitres['amp'], fitres['lam0'], fitres['sigma0'], 0.0
                else:
                    amp = None

            # BEFORE: raw residual + the scaled template about to be removed.
            fig.add_trace(go.Scatter(x=ww, y=target_resid[sel], mode='lines',
                                     name='raw', line=dict(color='black', width=TRACE_WIDTH),
                                     showlegend=(i == 0 and j == 0)),
                          row=row, col=col_before)
            if has_correction:
                fig.add_trace(go.Scatter(x=ww, y=k * template_resid[sel], mode='lines',
                                         name=f'k*{template_ext}',
                                         line=dict(color='grey', width=TRACE_WIDTH, dash='dash'),
                                         showlegend=(i == 0 and j == 0)),
                              row=row, col=col_before)
            fig.add_hline(y=0, line=dict(color='grey', width=1, dash='dot'), row=row, col=col_before)

            # AFTER: corrected residual (raw if no correction applies here) + fit.
            after_y = corrected[sel] if (has_correction and in_band) else target_resid[sel]
            fig.add_trace(go.Scatter(x=ww, y=after_y, mode='lines',
                                     name='corrected' if (has_correction and in_band) else 'raw (uncorrected)',
                                     line=dict(color='firebrick', width=TRACE_WIDTH),
                                     showlegend=(i == 0 and j == 0)),
                          row=row, col=col_after)
            if amp is not None and np.isfinite(amp):
                yfit = _gaussian(ww, amp, lam0, sigma0, bkg)
                fig.add_trace(go.Scatter(x=ww, y=yfit, mode='lines',
                                         name='Gaussian fit', line=dict(color='seagreen', width=FIT_WIDTH),
                                         showlegend=(i == 0 and j == 0)),
                              row=row, col=col_after)
            fig.add_hline(y=0, line=dict(color='grey', width=1, dash='dot'), row=row, col=col_after)

            # Match BEFORE/AFTER y-ranges so the reduction in amplitude/scatter
            # is visually obvious rather than hidden by independent autoscaling.
            combined = np.concatenate([target_resid[sel], after_y])
            finite = combined[np.isfinite(combined)]
            if finite.size:
                pad = 0.1 * (finite.max() - finite.min() or abs(finite.max()) or 1.0)
                yr = [finite.min() - pad, finite.max() + pad]
                fig.update_yaxes(range=yr, row=row, col=col_before)
                fig.update_yaxes(range=yr, row=row, col=col_after)

    fig.update_layout(
        title=dict(text=title or str(decomp_file), y=0.995, yanchor='top'),
        height=ROW_HEIGHT * n_rows + 60, width=420 * n_cols,
        margin=dict(t=130),
        # Legend sits below the title (y=1.02 in a taller top margin), not
        # sharing the same line -- the bigger legend font otherwise
        # overlapped the title text at the default margin.
        legend=dict(orientation='h', yanchor='bottom', y=1.02,
                   font=dict(size=SUBPLOT_TITLE_SIZE)))
    # Subplot titles are plain annotations plotly renders small by default
    # (noticeably smaller than the axis titles/ticks) -- bump them to match.
    fig.update_annotations(font_size=SUBPLOT_TITLE_SIZE)
    resid_label = 'Residual flux (erg s⁻¹ cm⁻² Å⁻¹)'
    # Plotly's default SI-prefix tick format ('50f' for 5e-14, etc.) is not
    # how flux values this small are conventionally read -- force plain
    # scientific notation instead, on every y-axis in the figure.
    fig.update_yaxes(exponentformat='e')
    for r in (1, 2):
        fig.update_yaxes(title_text=resid_label, row=r, col=1)
        fig.update_xaxes(title_text='Wavelength (Å)', row=r, col=1)
    for i in range(len(lines)):
        row = i + 3
        for c in range(n_cols):
            fig.update_yaxes(title_text=resid_label, row=row, col=c + 1)
            fig.update_xaxes(title_text='Wavelength (Å)', row=row, col=c + 1)
    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Visualize the DecomposeCleanSky.py / sky_nebular_leak_eval.py '
                    'shared-systematic and template-correction workflow for one exposure.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('decomp_file', help='DecomposeCleanSky.py output FITS file')
    parser.add_argument('-lines', default='siii_a,siii_b',
                        help='Comma-separated NEBULAR_LINES names to show zoom panels for')
    parser.add_argument('-targets', default='FLUX,SKY_WEST',
                        help='Comma-separated extensions to show zoom panels for')
    parser.add_argument('-template_ext', default='SKY_EAST',
                        help='Extension used as the shared-systematic template')
    parser.add_argument('-template_band', default='9000,9600', help='LO,HI (Angstrom)')
    parser.add_argument('-leak_file', default=None,
                        help='sky_nebular_leak_eval.py output table (<root>_lines.fits)')
    parser.add_argument('-v', dest='vel', type=float, default=0.0,
                        help='Nebular systemic velocity (km/s)')
    parser.add_argument('-title', default=None)
    parser.add_argument('-outfile', default=None,
                        help='Output HTML path (default: Overview_Plot/<stem>.nebleak.html)')
    args = parser.parse_args()

    lines = [s.strip() for s in args.lines.split(',') if s.strip()]
    targets = [s.strip() for s in args.targets.split(',') if s.strip()]
    band = tuple(float(x) for x in args.template_band.split(','))

    fig = plot_nebular_leak(args.decomp_file, lines=lines, targets=targets,
                            template_ext=args.template_ext, band=band,
                            leak_file=args.leak_file, title=args.title, vel=args.vel)

    outpath = Path(args.outfile) if args.outfile else \
        Path('Overview_Plot') / (Path(args.decomp_file).stem + '.nebleak.html')
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(outpath))
    print(f"Wrote {outpath}")


if __name__ == '__main__':
    main()
