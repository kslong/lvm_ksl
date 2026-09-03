#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Master (stacked) residual spectrum -- median FLUX_OBS - FLUX_PRED at
    every wavelength pixel across a whole BatchPredictSky.py corpus --
    split into 3 moon-brightness terciles by percentile of each
    exposure's own observed blue-continuum flux (a direct, model-
    independent brightness proxy; no moon geometry/ephemeris needed).
    Method-agnostic in the same sense as sky_residual_eval.py: works on
    any (wave, flux_obs, flux_pred) batch, however flux_pred was made.

Command line usage (if any):

    usage: PlotSkyResiduals.py [-h] [-mask_file PATH] [-outfile PATH]
                               [-html PATH]
                               batch_fits

    where

    batch_fits      a BatchPredictSky.py output FITS file (WAVE,
                    FLUX_OBS, FLUX_PRED, META[row,expnum]).

    -mask_file PATH
                    palace_make_mask.py mask FITS (default:
                    data/sky_mask.fits, same default as PlotSpec.py/
                    sky_residual_eval.py).

    -outfile PATH   output FITS path for the 3 master residual spectra +
                    percentile bands (default: <stem>_master_resid.fits).

    -html PATH      output interactive Plotly HTML comparing the 3
                    tercile median residuals (default:
                    Overview_Plot/<stem>_master_resid.html).

    -py_progs_dir PATH
                    path to the lvm_ksl repo's py_progs/ directory, which
                    supplies PlotSpec.get_sky_mask/GetSkyCont (default:
                    ~/SDSS/lvm_ksl/py_progs).

Description:

    1. Per exposure, the blue-continuum brightness proxy is
       median(FLUX_OBS) over clean (non-sky-line) pixels in
       GetSkyCont.ARM_EVAL_RANGES' B band (3650-5775 A) -- the observed
       spectrum itself, not a model-derived quantity, so this split is
       available regardless of which candidate produced FLUX_PRED.
    2. Exposures are cut into 3 terciles (dark/medium/bright) at the
       33rd/67th percentile of that proxy.
    3. Within each tercile, the median (and 16/84th percentile band) of
       FLUX_OBS - FLUX_PRED is computed at every wavelength pixel across
       all its exposures -- reveals persistent systematic residual
       structure a per-band or per-named-line aggregate could dilute or
       miss, and whether it depends on how bright the sky is.

History::

    260902  ksl  Coding begun.
    260903  ksl  Renamed from MasterResidualByMoon.py -- name now matches
        the PlotSpec*.py/PlotPredictSky.py convention (verb-first,
        "what it does" rather than "how it does it").
    260903  ksl  Switched every option from double-dash (--mask-file) to
        single-dash (-mask_file), matching py_progs/'s convention -- see
        BatchPredictSkyESO.py's History for the fuller note.

'''

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
import plotly.graph_objects as go

DEFAULT_PY_PROGS_DIR = Path('~/SDSS/lvm_ksl/py_progs').expanduser()

_pre = argparse.ArgumentParser(add_help=False)
_pre.add_argument('-py_progs_dir', default=str(DEFAULT_PY_PROGS_DIR))
_pre_args, _ = _pre.parse_known_args()

sys.path.insert(0, _pre_args.py_progs_dir)
from PlotSpec import get_sky_mask  # noqa: E402
from GetSkyCont import _interp_mask_to_wave, ARM_EVAL_RANGES  # noqa: E402

TERCILE_LABELS = ('dark', 'medium', 'bright')


def compute_master_residuals(wave, flux_obs, flux_pred, line_pred, mask_file=None):
    '''
    Depends only on the three saved spectra (FLUX_OBS, FLUX_PRED,
    LINE_PRED) plus the existing sky mask used to pick which pixels are
    genuinely continuum -- no independent continuum fit of any kind.

    Per row:
      cont_pred = flux_pred - line_pred     (already-available predicted continuum)
      A, alpha  = log-log fit of flux_obs/cont_pred vs (wave/w_ref) at
                  clean (continuum) pixels, weighted by cont_pred**2 --
                  same color-term convention as sky_residual_eval.py/
                  GetSkyCont.py's CONT_ALPHA, solved as a genuine 2-term
                  fit (amplitude + color), not their small-deviation
                  linearisation
      scaled_continuum = cont_pred * A * (wave/w_ref)**alpha
      line_error = flux_obs - scaled_continuum - line_pred

    i.e. rescale the model's own predicted continuum by one smooth power
    law across the whole spectrum (an amplitude and a color/slope term,
    not per-arm, not a single flat number), subtract that from the
    observed spectrum to get an observed-line-like signal, then subtract
    the predicted lines -- what's left is the error attributable to the
    line prediction specifically. The left/raw view stays
    flux_obs - flux_pred, unchanged.

    Returns
    -------
    dict with keys 'proxy' (n_rows,), 'tercile_edges' (2,),
    'tercile_label' (n_rows,), 'clean' (n_wave,), 'log_amp'/'color_alpha'
    (n_rows,), 'w_ref' (float), and per-tercile 'median'/'p16'/'p84' (raw
    residual) + 'line_median'/'line_p16'/'line_p84' (line error)
    (n_wave,) arrays.
    '''
    sky_mask = get_sky_mask(mask_file)
    if sky_mask is None:
        raise RuntimeError(f'No sky mask found ({mask_file or "data/sky_mask.fits"}); '
                           'required for the blue-continuum brightness proxy.')
    mask_wave, mask_bool = sky_mask
    clean = _interp_mask_to_wave(mask_wave, mask_bool, wave)

    b_lo, b_hi = next((lo, hi) for name, lo, hi in ARM_EVAL_RANGES if name == 'B')
    b_clean = clean & (wave >= b_lo) & (wave <= b_hi)
    if not np.any(b_clean):
        raise RuntimeError('No clean pixels in the B band -- cannot form the brightness proxy.')

    proxy = np.nanmedian(flux_obs[:, b_clean], axis=1)

    edges = np.nanpercentile(proxy, [100.0 / 3.0, 200.0 / 3.0])
    tercile_idx = np.digitize(proxy, edges)  # 0=dark, 1=medium, 2=bright
    tercile_label = np.asarray(TERCILE_LABELS)[tercile_idx]

    resid = flux_obs - flux_pred

    # Rescale the predicted continuum by a smooth power law in wavelength
    # (amplitude A x (wave/w_ref)^alpha), not a single constant -- same
    # "CONT_ALPHA" color-term convention already used by
    # sky_residual_eval.py/GetSkyCont.py, just solved as a genuine 2-term
    # log-log fit (amplitude + color) rather than their small-deviation
    # linearisation, since our amplitude offsets aren't always small.
    # One continuous curve across the whole spectrum, not per-arm.
    cont_pred = flux_pred - line_pred
    w_ref = float(np.median(wave[clean]))
    log_w = np.log(wave / w_ref)
    n_rows = flux_obs.shape[0]
    log_amp = np.zeros(n_rows)
    color_alpha = np.zeros(n_rows)
    for i in range(n_rows):
        valid = clean & np.isfinite(flux_obs[i]) & (flux_obs[i] > 0) & (cont_pred[i] > 0)
        if valid.sum() < 3:
            continue
        weight = cont_pred[i, valid] ** 2
        yvar = np.log(flux_obs[i, valid] / cont_pred[i, valid])
        alpha, ln_a = np.polyfit(log_w[valid], yvar, 1, w=np.sqrt(weight))
        log_amp[i] = ln_a
        color_alpha[i] = alpha
    scale_curve = np.exp(log_amp[:, None] + color_alpha[:, None] * log_w[None, :])
    line_error = flux_obs - scale_curve * cont_pred - line_pred

    out = dict(proxy=proxy, tercile_edges=edges, tercile_label=tercile_label,
              clean=clean, log_amp=log_amp, color_alpha=color_alpha, w_ref=w_ref)
    for k, label in enumerate(TERCILE_LABELS):
        sel = tercile_idx == k
        out[f'median_{label}'] = np.nanmedian(resid[sel], axis=0)
        out[f'p16_{label}'] = np.nanpercentile(resid[sel], 16, axis=0)
        out[f'p84_{label}'] = np.nanpercentile(resid[sel], 84, axis=0)
        out[f'line_median_{label}'] = np.nanmedian(line_error[sel], axis=0)
        out[f'line_p16_{label}'] = np.nanpercentile(line_error[sel], 16, axis=0)
        out[f'line_p84_{label}'] = np.nanpercentile(line_error[sel], 84, axis=0)
        out[f'n_{label}'] = int(sel.sum())
    return out


RANGE_PERCENTILE = 99.0
"""Y-range ceiling is set from this percentile of |band edge|, not the raw
max -- a single outlier pixel (typically one strong OH line, or the worst
Z-band continuum point) otherwise sets the scale and flattens everything
else to invisibility. Same principle as this project's log-axis convention
(99th pct, not .max())."""


def _shared_range(result):
    '''99th-percentile |band edge| over the bright tercile (widest), raw
    and continuum-subtracted combined -- just the starting range;
    matches='y' below keeps all 6 panels synced to whatever range the
    viewer zooms/pans to afterward.'''
    vals = np.concatenate([
        result['p84_bright'], result['p16_bright'],
        result['line_p84_bright'], result['line_p16_bright'],
    ])
    vals = vals[np.isfinite(vals)]
    span = np.nanpercentile(np.abs(vals), RANGE_PERCENTILE) if vals.size else 1.0
    return [-1.15 * span, 1.15 * span]


def plot_master_residuals(wave, result, title=None):
    '''
    6 stacked panels, full width each (easier to read fine, ~100 A scale
    structure than a 3x2 grid): rows 1-3 = FLUX_OBS - FLUX_PRED (total
    error) for dark/medium/bright sky; rows 4-6 = FLUX_OBS -
    scale*(FLUX_PRED-LINE_PRED) - LINE_PRED (line-prediction error only)
    for the same three, in the same order -- see compute_master_residuals
    for the per-row continuum rescaling. All 6 panels' axes are linked
    (matches='x'/'y'): zooming or panning any one panel updates all six
    together. Per panel: explicit median, p16, and p84 lines (not just a
    shaded fill) -- no per-exposure spaghetti.
    '''
    from plotly.subplots import make_subplots

    row_titles = {'dark': 'Dark Sky', 'medium': 'Medium Sky', 'bright': 'Bright Sky'}
    panel_order = (
        [(label, '', 'Total error') for label in TERCILE_LABELS]
        + [(label, 'line_', 'Line error') for label in TERCILE_LABELS]
    )
    # Rows 1 and 5 are deliberately blank spacers (no traces) so the "Total
    # Error"/"Line Error" group headers have real allocated room -- centred
    # in their own spacer row's domain, exactly like the panel row titles
    # are positioned relative to their own domain, rather than guessed as
    # an offset from the title block (whose rendered height changes with
    # its text length/wrapping and isn't something the figure layout knows
    # about -- confirmed by inspecting fig.layout's actual y-domains before
    # this fix, same issue "Line Error" had before its own spacer row).
    panel_rows = [2, 3, 4, 6, 7, 8]
    subplot_titles = ['', row_titles['dark'], row_titles['medium'], row_titles['bright'], '',
                      row_titles['dark'], row_titles['medium'], row_titles['bright']]
    fig = make_subplots(
        rows=8, cols=1, vertical_spacing=0.022,
        row_heights=[0.6, 1, 1, 1, 0.6, 1, 1, 1],
        subplot_titles=subplot_titles,
    )

    yrange = _shared_range(result)
    clean = result['clean']

    for row, (label, prefix, _) in zip(panel_rows, panel_order):
        is_total = (prefix == '')
        for stat, color, width in (('median', 'black', 1.5),
                                   ('p16', 'grey', 1), ('p84', 'grey', 1)):
            if is_total:
                # Sky-line highlight drawn FIRST (underneath) and WIDER, so
                # it shows as a halo around the data line rather than being
                # hidden by it once the data line is drawn on top -- unlike
                # PlotSpec.py's -mask (which replaces the line color at
                # masked pixels), here the actual value must stay visible
                # everywhere, with line-affected pixels just flagged by the
                # grey band behind them. Total-error panels only.
                masked = np.where(~clean, result[f'{prefix}{stat}_{label}'], np.nan)
                fig.add_trace(go.Scatter(
                    x=wave, y=masked, mode='lines',
                    line=dict(color='lightgrey', width=width + 3, dash='solid'),
                    name=f'{stat} (line-affected)', legendgroup=label, showlegend=False,
                ), row=row, col=1)
            fig.add_trace(go.Scatter(
                x=wave, y=result[f'{prefix}{stat}_{label}'], mode='lines',
                line=dict(color=color, width=width, dash='solid'),
                name=stat, legendgroup=label, showlegend=False,
            ), row=row, col=1)
        fig.add_hline(y=0.0, line=dict(color='lightgrey', width=1, dash='dot'), row=row, col=1)

    # matches='x'/'y' link every panel's zoom/pan together live (not just a
    # one-time shared starting range) -- boxed axes on every panel
    # (plotly_white only draws bottom+left by default; mirror=True adds
    # the matching top/right line).
    fig.update_xaxes(matches='x', tickformat='.0f', showline=True, linecolor='black',
                     mirror=True, ticks='outside')
    fig.update_yaxes(matches='y', range=yrange, exponentformat='e', showline=True,
                     linecolor='black', mirror=True, ticks='outside')
    # Rows 1 and 5 are the blank spacers -- hide their axes entirely (no
    # box, no ticks).
    for spacer_row in (1, 5):
        fig.update_xaxes(showline=False, showticklabels=False, row=spacer_row, col=1)
        fig.update_yaxes(showline=False, showticklabels=False, matches=None, row=spacer_row, col=1)
    fig.update_xaxes(title_text='Wavelength (Å)', row=8, col=1)
    fig.update_yaxes(title_text='FLUX_OBS − FLUX_PRED', row=3, col=1)
    fig.update_yaxes(title_text='FLUX_OBS − FLUX_PRED', row=7, col=1)

    # "Total Error"/"Line Error" centred in their own spacer row's domain --
    # real allocated space, so neither can land on top of a plot panel or
    # the title block, regardless of how tall the title text renders.
    for header_text, spacer_row in (('Total Error', 1), ('Line Error', 5)):
        axis_name = 'yaxis' if spacer_row == 1 else f'yaxis{spacer_row}'
        lo_dom, hi_dom = fig.layout[axis_name].domain
        fig.add_annotation(
            text=f'<b>{header_text}</b>', xref='paper', yref='paper',
            x=0.5, y=(lo_dom + hi_dom) / 2, showarrow=False,
            font=dict(size=24), xanchor='center', yanchor='middle',
        )

    # Title block, top to bottom: (1) what's being differenced, (2) how to
    # read a panel (black/grey, what the 16/84th percentiles mean, Total
    # vs. Line Error), (3) provenance -- source file, spectrum count, and
    # the brightness boundary for each of the three panels -- last, since
    # it's about where this particular run's data came from, not about how
    # to read the plot.
    opening_line = 'The difference between a set of Observed Spectra and (Model) Sky Spectrum.'
    explain_line = (
        'Each panel: black = median of (observed - predicted) sky flux at that wavelength, across all '
        'exposures in that brightness/error group. Grey = the 16th and 84th percentile of that same '
        'distribution across exposures -- the band between them is where the middle 68% of exposures fall '
        '(a robust, non-parametric analogue of a 1-sigma range). Total Error = FLUX_OBS - FLUX_PRED; '
        'Line Error = the same, with a color-rescaled predicted continuum and the predicted lines removed first.'
    )
    lo, hi = result['tercile_edges']
    n_total = sum(result[f'n_{label}'] for label in TERCILE_LABELS)
    source_name = title or 'unknown file'
    data_line = f'Data from: {source_name} with {n_total} spectra. Continuum limits (erg/s/cm2/A):'
    bullet_lines = [
        f'• Dark Sky: continuum less than {lo:.2e}',
        f'• Medium Sky: continuum between {lo:.2e} and {hi:.2e}',
        f'• Bright Sky: continuum greater than {hi:.2e}',
    ]
    # No <sup> here -- it silently shrinks its contents regardless of the
    # font size set below, which is why the explanation/data lines stayed
    # small even after the font size was already increased once.
    full_title = '<br>'.join([opening_line, explain_line, data_line] + bullet_lines)

    fig.update_layout(template='plotly_white',
                      title=dict(text=full_title, font=dict(size=16)),
                      showlegend=False, height=2350, width=1400,
                      margin=dict(t=420))
    return fig


def main():
    p = argparse.ArgumentParser(
        parents=[_pre],
        description='Master residual spectrum split into 3 moon-brightness '
                    '(blue-continuum) terciles.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('batch_fits', help='BatchPredictSky.py output FITS file')
    p.add_argument('-mask_file', default=None, dest='mask_file',
                   help='palace_make_mask.py mask FITS')
    p.add_argument('-outfile', default=None, help='output FITS path')
    p.add_argument('-html', default=None, help='output interactive HTML path')
    args = p.parse_args()

    with fits.open(args.batch_fits) as hdul:
        wave = np.asarray(hdul['WAVE'].data, dtype=np.float64)
        flux_obs = np.asarray(hdul['FLUX_OBS'].data, dtype=np.float64)
        flux_pred = np.asarray(hdul['FLUX_PRED'].data, dtype=np.float64)
        line_pred = np.asarray(hdul['LINE_PRED'].data, dtype=np.float64)

    result = compute_master_residuals(wave, flux_obs, flux_pred, line_pred, mask_file=args.mask_file)
    lo, hi = result['tercile_edges']
    print(f'Blue-continuum brightness proxy tercile edges: {lo:.3g}, {hi:.3g} erg/s/cm^2/A')
    for label in TERCILE_LABELS:
        print(f'  {label}: n={result[f"n_{label}"]}')

    outfile = args.outfile
    if outfile is None:
        outfile = str(Path(args.batch_fits).with_suffix('').as_posix()) + '_master_resid.fits'
    cols = [fits.Column(name='WAVE', format='D', array=wave)]
    for label in TERCILE_LABELS:
        for prefix in ('', 'line_'):
            for stat in ('median', 'p16', 'p84'):
                key = f'{prefix}{stat}_{label}'
                cols.append(fits.Column(name=key.upper(), format='D', array=result[key]))
    fits.BinTableHDU.from_columns(cols).writeto(outfile, overwrite=True)
    print(f'Wrote {outfile}')

    fig = plot_master_residuals(wave, result, title=Path(args.batch_fits).name)
    html = args.html
    if html is None:
        outdir = Path('Overview_Plot')
        outdir.mkdir(exist_ok=True)
        html = str(outdir / (Path(args.batch_fits).stem + '_master_resid.html'))
    fig.write_html(html, include_plotlyjs=True)
    print(f'Wrote {html}')


if __name__ == '__main__':
    main()
