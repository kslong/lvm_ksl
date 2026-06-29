#!/usr/bin/env python
# coding: utf-8
'''
                   Space Telescope Science Institute

Synopsis:

    Evaluate continuum fits produced by GetSkyCont.py by plotting flux,
    continuum, components, and residuals over a chosen wavelength window.
    Reads a FITS file produced by GetSkyCont.py and writes an interactive
    Plotly HTML file.  Also computes per-spectrum fit-quality statistics
    and writes them back to the DRP_ALL table in the input FITS file.

Command line usage (if any):

    usage: GetSkyCont_eval.py skycont_file.fits [wmin wmax] [-num N] [-out outroot]

    Arguments:

    skycont_file.fits  FITS file from GetSkyCont.py (process_skyfile or process_many).
    wmin               minimum wavelength in Angstroms (default 3600).
    wmax               maximum wavelength in Angstroms (default 9800).

    Options:

    -num N             overlay N randomly-selected individual spectra on the band (default 20; 0 = band only).
    -out outroot       output filename root; default is <stem>_<wmin>_<wmax>.

Description:

    Writes a single HTML file containing three interactive Plotly figures and
    updates the DRP_ALL table in the input FITS file with per-spectrum
    fit-quality statistics.

    Figure 1 — four-panel spectral overview:

    Panel 1 (Flux + Continuum, log): observed FLUX with the fitted CONT median
    overlaid in red.  Blue vertical bands mark the clean (fit) regions;
    grey bands mark the masked (line-affected, excluded) regions.

    Panel 2 (Total Continuum, log): the total CONT band, showing the overall
    smoothness and level of the B-spline model.

    Panel 3 (Components, log): MOON (red; B-spline x solar) and DIFFUSE
    (orange; plain B-splines) plotted together with the total CONT median
    (purple) for reference.  The relative amplitude of the components shows
    how much of the continuum is Moon/zodiacal vs. diffuse airglow.

    Panel 4 (Residual, linear): RESID = FLUX - CONT; should be near zero in
    clean regions and show sky line emission in the masked regions.

    Both log panels share a floor at 1e-18 to suppress near-zero noise.

    Figure 2 — per-arm residual histograms:

    Three panels (Blue 3600-5900 A, Red 5900-7600 A, NIR 7600-9800 A) each
    showing the distribution of all residual flux values (all spectra, all
    clean pixels) in that wavelength range.  The x-axis is residual flux; the
    y-axis is count N.  A Gaussian with center = median and sigma = NMAD is
    overlaid in black.  An annotation box in each panel reports N, the median,
    NMAD, skewness, and the 10th/90th percentiles.  The histogram range is
    clipped to median +/- 5*NMAD to suppress extreme outliers.

    Figure 3 — per-spectrum fit quality:

    Row 1: three scatter plots (Blue, Red, NIR) of per-spectrum median
    residual (x) vs NMAD (y).  Each point is one spectrum; hover shows the
    spectrum row index and its median/NMAD.  A dashed vertical line at x=0
    marks the ideal median; a dashed horizontal line marks the ensemble NMAD
    for reference.

    Row 2: ranked NMAD plot showing all three arms on the same panel, with
    spectra sorted worst-to-best by their overall (all clean pixels) NMAD.
    The x-axis is rank (0 = worst); hover shows the original spectrum row
    index.  A horizontal dashed line marks the ensemble NMAD per arm.

    Per-spectrum statistics written to DRP_ALL:

    For each arm (blue, red, nir) and for all arms combined (all), four
    float32 columns are added (or updated) in the DRP_ALL table:
        resid_med_<arm>   — median residual in clean pixels
        resid_nmad_<arm>  — NMAD of residuals in clean pixels
        resid_rms_<arm>   — RMS of residuals in clean pixels
        resid_skew_<arm>  — skewness of residuals in clean pixels

Primary routines:

    plot_eval   build and write the Plotly HTML evaluation plot.

Notes:

    Requires plotly (pip install plotly).
    Wavelength window defaults to the full LVM range (3600-9800 A) when
    wmin and wmax are not given.
    The input FITS file is updated in place; a backup is not created.

History:

    260628  ksl  Written.
    260628  ksl  Redesigned to four panels: Flux, Total Cont, Components, Residual.
    260629  ksl  Added per-arm residual histogram figure (Figure 2).
    260629  ksl  Added per-spectrum fit-quality figure (Figure 3) and DRP_ALL update.
'''

import sys
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.table import Table
import plotly.graph_objects as go
from plotly.subplots import make_subplots

ARM_DEFS = [
    ('Blue', 3600.0, 5900.0, 'rgba(31,119,180,0.6)'),
    ('Red',  5900.0, 7600.0, 'rgba(214,39,40,0.6)'),
    ('NIR',  7600.0, 9800.0, 'rgba(44,160,44,0.6)'),
]

_ARM_SOLID = {
    'Blue': 'rgb(31,119,180)',
    'Red':  'rgb(214,39,40)',
    'NIR':  'rgb(44,160,44)',
}


# ──────────────────────────────────────────────────────────────
# Trace helpers
# ──────────────────────────────────────────────────────────────

def _band_traces(wave, arr, name, color_band, color_med, show_band_legend, legend_ref='legend'):
    '''Median line + shaded 10/90-pct band traces.'''
    p10, p90 = np.percentile(arr, [10, 90], axis=0)
    med      = np.median(arr, axis=0)
    lower = go.Scatter(x=wave, y=p10, mode='lines', line=dict(width=0),
                       showlegend=False, hoverinfo='skip')
    upper = go.Scatter(x=wave, y=p90, mode='lines', line=dict(width=0),
                       fill='tonexty', fillcolor=color_band,
                       name=f'{name} 10–90 pct', showlegend=show_band_legend,
                       legend=legend_ref, hoverinfo='skip')
    median = go.Scatter(x=wave, y=med, mode='lines',
                        line=dict(color=color_med, width=1.5),
                        name=f'Median {name}', showlegend=True, legend=legend_ref)
    return [lower, upper, median]


def _sample_traces(wave, arr, idx):
    '''Semi-transparent individual spectrum traces.'''
    n = len(idx)
    if n == 0:
        return []
    alpha = min(0.7, 3.0 / max(1.0, n ** 0.5))
    return [
        go.Scatter(x=wave, y=arr[i], mode='lines',
                   line=dict(color=f'rgba(120,120,120,{alpha:.2f})', width=0.6),
                   showlegend=False, hoverinfo='skip')
        for i in idx
    ]


def _log_range(arr):
    '''Return (log_lo, log_hi) for a log axis from array arr.

    Uses the 99th percentile of positive values as the ceiling (immune to
    extreme outlier pixels) and the median of positive values minus 2 decades
    as the floor, capped at -18.
    '''
    pos = arr[arr > 0]
    if not len(pos):
        return -18.0, -14.0
    med    = np.median(pos)
    p99    = np.percentile(pos, 99)
    log_lo = max(np.log10(med) - 2.0, -18.0)
    log_hi = np.ceil(np.log10(p99)) + 0.5
    return log_lo, log_hi


def _contiguous_regions(wave, flag):
    '''Return (wlo, whi) tuples for contiguous wavelength regions where flag is True.'''
    if not np.any(flag):
        return []
    change = np.diff(flag.astype(int))
    starts = np.where(change ==  1)[0]      # last clean pixel before masked region
    ends   = np.where(change == -1)[0] + 1  # first clean pixel after masked region
    if flag[0]:
        starts = np.r_[0, starts]
    if flag[-1]:
        ends = np.r_[ends, len(wave)]
    return [(float(wave[s]), float(wave[min(e, len(wave) - 1)]))
            for s, e in zip(starts, ends)]


def _arm_stats(wave, clean, resid, wlo, whi):
    '''Residual statistics for clean pixels within [wlo, whi] Angstroms.

    Parameters
    ----------
    wave  : 1-D array, wavelength of display pixels
    clean : 1-D bool array, True = clean (used in fit)
    resid : 2-D array (n_spec, n_wave), residual spectra
    wlo, whi : float, wavelength range

    Returns
    -------
    dict with keys n, med, nmad, skew, p10, p90, vals  — or None if no data.
    '''
    sel = (wave >= wlo) & (wave <= whi) & clean
    if not sel.any():
        return None
    vals = resid[:, sel].flatten()
    vals = vals[np.isfinite(vals)]
    if len(vals) < 10:
        return None
    med  = float(np.median(vals))
    nmad = float(1.4826 * np.median(np.abs(vals - med)))
    mn   = float(np.mean(vals))
    sd   = float(np.std(vals, ddof=0))
    skew = float(np.mean(((vals - mn) / sd) ** 3)) if sd > 0 else 0.0
    p10, p90 = (float(v) for v in np.percentile(vals, [10, 90]))
    return dict(n=len(vals), med=med, nmad=nmad, skew=skew,
                p10=p10, p90=p90, vals=vals)


def _per_spec_stats(wave, clean, resid, wlo, whi):
    '''Per-spectrum residual statistics within [wlo, whi] Angstroms.

    Uses full-resolution arrays so downsampling does not bias the statistics.

    Parameters
    ----------
    wave  : 1-D array (n_pix,), full-resolution wavelength
    clean : 1-D bool array (n_pix,), True = clean pixel
    resid : 2-D array (n_spec, n_pix), full-resolution residuals
    wlo, whi : float, wavelength bounds

    Returns
    -------
    dict with 1-D float arrays of length n_spec: med, nmad, rms, skew
    or None if no clean pixels fall within the range.
    '''
    sel = (wave >= wlo) & (wave <= whi) & clean
    if not sel.any():
        return None
    r = resid[:, sel].astype(float)
    r[~np.isfinite(r)] = np.nan
    med  = np.nanmedian(r, axis=1)
    nmad = 1.4826 * np.nanmedian(np.abs(r - med[:, np.newaxis]), axis=1)
    rms  = np.sqrt(np.nanmean(r ** 2, axis=1))
    mn   = np.nanmean(r, axis=1)
    sd   = np.nanstd(r, axis=1, ddof=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        skew = np.where(sd > 0,
                        np.nanmean(((r - mn[:, np.newaxis]) / sd[:, np.newaxis]) ** 3,
                                   axis=1),
                        np.zeros(len(sd)))
    return dict(med=med, nmad=nmad, rms=rms, skew=skew)


def _fmt(v):
    '''Format a flux value compactly in scientific notation.'''
    if not np.isfinite(v) or abs(v) < 1e-30:
        return '0'
    return f'{v:.2e}'


# ──────────────────────────────────────────────────────────────
# Main plot function
# ──────────────────────────────────────────────────────────────

def plot_eval(filename, wmin=3600.0, wmax=9800.0, n_sample=20, outroot=''):
    '''Build spectral evaluation plots and per-spectrum quality diagnostics
    from a GetSkyCont output FITS file.

    Writes a three-figure HTML file and updates the DRP_ALL table in the
    input FITS file with per-spectrum fit-quality statistics (med, nmad, rms,
    skew) for each spectrograph arm and for all arms combined.

    Parameters
    ----------
    filename : str or Path
        skycont_*.fits file produced by GetSkyCont.py.
    wmin : float
        Minimum wavelength in Angstroms.
    wmax : float
        Maximum wavelength in Angstroms.
    n_sample : int
        Number of randomly-selected individual spectra to overlay.  0 = band only.
    outroot : str
        Output filename root.  Default: <stem>_<wmin>_<wmax>.
    '''
    # ── Read FITS ─────────────────────────────────────────────
    # Build pixel index: wavelength window + downsample to ~2000 display points.
    # This lets each 2D array be loaded and sliced in one step, avoiding large
    # float64 intermediates (5 arrays × N_spec × 12401 before reduction).
    hdul      = fits.open(filename)
    wave_all  = hdul['WAVE'].data.astype(float)
    mask_all  = hdul['MASK'].data.astype(bool)

    wmask   = (wave_all >= wmin) & (wave_all <= wmax)
    ds      = max(1, int(wmask.sum()) // 2000)
    pix_idx = np.where(wmask)[0][::ds]

    wave  = wave_all[pix_idx]
    clean = mask_all[pix_idx]

    flux    = hdul['FLUX'].data[:,    pix_idx].astype(float)
    cont    = hdul['CONT'].data[:,    pix_idx].astype(float)
    moon    = hdul['MOON'].data[:,    pix_idx].astype(float)
    diffuse = hdul['DIFFUSE'].data[:, pix_idx].astype(float)
    resid   = hdul['RESID'].data[:,   pix_idx].astype(float)
    # Full-resolution residual for per-spectrum statistics (not downsampled)
    resid_all = hdul['RESID'].data.astype(float)
    drp_table = Table(hdul['DRP_ALL'].data)
    hdul.close()

    if flux.ndim == 1:
        flux      = flux[np.newaxis, :]
        cont      = cont[np.newaxis, :]
        moon      = moon[np.newaxis, :]
        diffuse   = diffuse[np.newaxis, :]
        resid     = resid[np.newaxis, :]
        resid_all = resid_all[np.newaxis, :]

    n_spec = len(flux)

    rng        = np.random.default_rng(42)
    sample_idx = (rng.choice(n_spec, size=min(n_sample, n_spec), replace=False)
                  if n_sample > 0 else [])

    # ── Per-spectrum statistics (full resolution) ──────────────
    per_spec = {name: _per_spec_stats(wave_all, mask_all, resid_all, wlo, whi)
                for name, wlo, whi, _ in ARM_DEFS}
    per_spec['All'] = _per_spec_stats(wave_all, mask_all, resid_all, 3600.0, 9800.0)
    del wave_all, mask_all, resid_all

    # ── Update DRP_ALL in the input FITS file ──────────────────
    _sfx = {'Blue': 'blue', 'Red': 'red', 'NIR': 'nir', 'All': 'all'}
    for arm_name, st in per_spec.items():
        if st is None:
            continue
        sfx = _sfx[arm_name]
        for stat_name, vals in [('med',  st['med']),  ('nmad', st['nmad']),
                                 ('rms',  st['rms']),  ('skew', st['skew'])]:
            col = f'resid_{stat_name}_{sfx}'
            if col in drp_table.colnames:
                drp_table.remove_column(col)
            drp_table[col] = vals.astype(np.float32)

    with fits.open(filename, mode='update') as hdul_upd:
        new_drp = fits.table_to_hdu(drp_table)
        new_drp.name = 'DRP_ALL'
        for i, h in enumerate(hdul_upd):
            if h.name == 'DRP_ALL':
                hdul_upd[i] = new_drp
                break
        hdul_upd.flush()
    print(f'Updated DRP_ALL in {filename}')

    # ── Clip to 1e-18 for log-scale display ───────────────────
    # Keeps traces contiguous with no NaN artefacts.
    _floor    = 1e-18
    flux_l    = np.maximum(flux,    _floor)
    cont_l    = np.maximum(cont,    _floor)
    moon_l    = np.maximum(moon,    _floor)
    diffuse_l = np.maximum(diffuse, _floor)
    cont_med  = np.median(cont_l, axis=0)

    # Per-arm ensemble residual statistics (for Figures 2 and 3 references).
    arm_stats = {name: _arm_stats(wave, clean, resid, wlo, whi)
                 for name, wlo, whi, _ in ARM_DEFS}

    stem = Path(filename).stem

    # ══════════════════════════════════════════════════════════
    # Figure 1 — spectral panels (4-row, 1-col, shared x-axis)
    # ══════════════════════════════════════════════════════════
    fig = make_subplots(
        rows=4, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.05,
        subplot_titles=[
            f'Flux + Total Continuum  ({n_spec} spectra)  — log',
            'Total Continuum  — log',
            'Components: Moon (red) + Diffuse (orange)  — log',
            'Residual  (Flux − Continuum)  — linear',
        ],
    )

    # ── Panel 1: FLUX + total CONT median (log) ───────────────
    for tr in _sample_traces(wave, flux_l, sample_idx):
        fig.add_trace(tr, row=1, col=1)
    for tr in _band_traces(wave, flux_l, 'Flux', 'rgba(31,119,180,0.20)', 'rgb(31,119,180)',
                           show_band_legend=True, legend_ref='legend'):
        fig.add_trace(tr, row=1, col=1)
    fig.add_trace(
        go.Scatter(x=wave, y=cont_med, mode='lines',
                   line=dict(color='rgb(214,39,40)', width=1.5),
                   name='Median Cont', showlegend=True, legend='legend'),
        row=1, col=1,
    )

    # ── Panel 2: Total CONT band (log) ────────────────────────
    for tr in _band_traces(wave, cont_l, 'Cont', 'rgba(214,39,40,0.20)', 'rgb(214,39,40)',
                           show_band_legend=True, legend_ref='legend2'):
        fig.add_trace(tr, row=2, col=1)

    # ── Panel 3: MOON + DIFFUSE + total median (log) ──────────
    for tr in _band_traces(wave, moon_l, 'Moon', 'rgba(214,39,40,0.20)', 'rgb(214,39,40)',
                           show_band_legend=True, legend_ref='legend3'):
        fig.add_trace(tr, row=3, col=1)
    for tr in _band_traces(wave, diffuse_l, 'Diffuse', 'rgba(255,127,14,0.20)', 'rgb(255,127,14)',
                           show_band_legend=True, legend_ref='legend3'):
        fig.add_trace(tr, row=3, col=1)
    fig.add_trace(
        go.Scatter(x=wave, y=cont_med, mode='lines',
                   line=dict(color='rgb(148,103,189)', width=1.5),
                   name='Median Total', showlegend=True, legend='legend3'),
        row=3, col=1,
    )

    # ── Panel 4: RESID (linear) ────────────────────────────────
    for tr in _sample_traces(wave, resid, sample_idx):
        fig.add_trace(tr, row=4, col=1)
    for tr in _band_traces(wave, resid, 'Resid', 'rgba(44,160,44,0.20)', 'rgb(44,160,44)',
                           show_band_legend=True, legend_ref='legend4'):
        fig.add_trace(tr, row=4, col=1)

    # ── Fit-region indicator on panel 1 ───────────────────────
    _f_lo, _ = _log_range(flux)
    fit_y_val = 10 ** (_f_lo + 0.25)
    fit_y     = np.where(clean, fit_y_val, np.nan)
    fig.add_trace(
        go.Scatter(x=wave, y=fit_y, mode='lines',
                   line=dict(color='rgb(0,150,0)', width=3),
                   name='Fit pixels', showlegend=True, legend='legend'),
        row=1, col=1,
    )

    # ── Grey shading for masked regions ───────────────────────
    masked_regions = _contiguous_regions(wave, ~clean)
    _yref = {1: 'y domain', 2: 'y2 domain', 3: 'y3 domain', 4: 'y4 domain'}
    shapes = [
        dict(type='rect', xref='x', yref=_yref[row],
             x0=w0, x1=w1, y0=0, y1=1,
             fillcolor='rgba(180,180,180,0.30)', layer='below',
             line=dict(width=0))
        for w0, w1 in masked_regions
        for row in [1, 2, 3, 4]
    ]

    _leg = dict(xanchor='right', yanchor='top', x=0.99,
                bgcolor='rgba(255,255,255,0.7)')
    fig.update_layout(
        title=f'{stem}  {int(wmin)}–{int(wmax)} Å',
        height=1100,
        template='simple_white',
        shapes=shapes,
        legend  = dict(**_leg, y=0.99),
        legend2 = dict(**_leg, y=0.73),
        legend3 = dict(**_leg, y=0.49),
        legend4 = dict(**_leg, y=0.24),
    )

    _ax = dict(showline=True, linewidth=1, linecolor='black',
               mirror=True, ticks='outside', ticklen=4, showticklabels=True)
    fig.update_xaxes(**_ax)
    fig.update_yaxes(**_ax, exponentformat='e', showexponent='all')

    fig.update_xaxes(range=[wmin, wmax])
    for row in [1, 2, 3]:
        fig.update_xaxes(showticklabels=False, row=row, col=1)
    fig.update_xaxes(title_text='Wavelength (Å)', row=4, col=1)

    _comp_range = list(_log_range(cont))
    fig.update_yaxes(title_text='Flux',       type='log',                    row=1, col=1)
    fig.update_yaxes(title_text='Continuum',  type='log', range=_comp_range, row=2, col=1)
    fig.update_yaxes(title_text='Components', type='log', range=_comp_range, row=3, col=1)
    fig.update_yaxes(title_text='Residual',                                   row=4, col=1)

    # ══════════════════════════════════════════════════════════
    # Figure 2 — per-arm residual histograms (separate figure)
    # ══════════════════════════════════════════════════════════
    fig2 = make_subplots(
        rows=1, cols=3,
        horizontal_spacing=0.08,
        subplot_titles=[
            f'Residual distribution  Blue (3600–5900 Å)',
            f'Residual distribution  Red (5900–7600 Å)',
            f'Residual distribution  NIR (7600–9800 Å)',
        ],
    )

    for col_i, (name, wlo, whi, color) in enumerate(ARM_DEFS, start=1):
        st = arm_stats[name]
        if st is None:
            continue
        lo = st['med'] - 5.0 * st['nmad']
        hi = st['med'] + 5.0 * st['nmad']
        if lo >= hi:
            lo, hi = st['p10'] - 1e-20, st['p90'] + 1e-20
        counts, edges = np.histogram(st['vals'], bins=60, range=(lo, hi))
        centers = 0.5 * (edges[:-1] + edges[1:])
        bw      = edges[1] - edges[0]

        fig2.add_trace(
            go.Bar(x=centers, y=counts, marker_color=color,
                   marker_line_width=0, showlegend=False),
            row=1, col=col_i,
        )

        xg = np.linspace(lo, hi, 300)
        yg = (st['n'] * bw / (np.sqrt(2 * np.pi) * st['nmad'])
              * np.exp(-0.5 * ((xg - st['med']) / st['nmad']) ** 2))
        fig2.add_trace(
            go.Scatter(x=xg, y=yg, mode='lines',
                       line=dict(color='black', width=1.5),
                       showlegend=False),
            row=1, col=col_i,
        )

        ax_sfx = '' if col_i == 1 else str(col_i)
        fig2.add_annotation(
            x=0.03, y=0.97,
            xref=f'x{ax_sfx} domain', yref=f'y{ax_sfx} domain',
            xanchor='left', yanchor='top',
            showarrow=False, align='left',
            font=dict(size=10, family='monospace'),
            bgcolor='rgba(255,255,255,0.75)',
            bordercolor='grey', borderwidth=1,
            text=(f'N = {st["n"]}<br>'
                  f'Median = {_fmt(st["med"])}<br>'
                  f'NMAD   = {_fmt(st["nmad"])}<br>'
                  f'Skew   = {st["skew"]:.2f}<br>'
                  f'p10    = {_fmt(st["p10"])}<br>'
                  f'p90    = {_fmt(st["p90"])}'),
        )

    fig2.update_layout(
        height=450,
        template='simple_white',
        bargap=0.02,
        showlegend=False,
    )
    fig2.update_xaxes(**_ax, title_text='Residual flux',
                      exponentformat='e', showexponent='all')
    fig2.update_yaxes(**_ax, title_text='N',
                      exponentformat='none', showexponent='none')

    # ══════════════════════════════════════════════════════════
    # Figure 3 — per-spectrum fit quality
    # ══════════════════════════════════════════════════════════

    fig3 = make_subplots(
        rows=2, cols=3,
        specs=[[{}, {}, {}], [{'colspan': 3}, None, None]],
        subplot_titles=[
            'Blue (3600–5900 Å)  Median vs NMAD',
            'Red (5900–7600 Å)  Median vs NMAD',
            'NIR (7600–9800 Å)  Median vs NMAD',
            f'Per-spectrum NMAD by spectrum number  ({n_spec} spectra)',
        ],
        vertical_spacing=0.18,
        horizontal_spacing=0.08,
    )

    # ── Row 1: scatter of median residual vs NMAD ─────────────
    for col_i, (name, wlo, whi, color) in enumerate(ARM_DEFS, start=1):
        st  = per_spec[name]
        ens = arm_stats[name]
        if st is None:
            continue

        # Percentile-based axis limits so a few extreme spectra do not
        # compress the scale for the bulk of the distribution.
        p2_med, p98_med = (float(v) for v in np.nanpercentile(st['med'], [2, 98]))
        p98_nmad        = float(np.nanpercentile(st['nmad'], 98))
        x_pad = max(abs(p98_med - p2_med) * 0.12, 1e-22)
        y_pad = max(p98_nmad * 0.12, 1e-22)
        x_lo, x_hi = p2_med - x_pad, p98_med + x_pad
        y_hi        = p98_nmad + y_pad

        in_range = ((st['med']  >= x_lo) & (st['med']  <= x_hi) &
                    (st['nmad'] >= 0)    & (st['nmad'] <= y_hi))
        n_off = int(np.sum(~in_range))

        # Replace out-of-range points with NaN so they are absent from the
        # trace data; autorange then cannot snap back to include them.
        med_plot  = np.where(in_range, st['med'],  np.nan)
        nmad_plot = np.where(in_range, st['nmad'], np.nan)

        hover = [f'spec {i}<br>med  = {_fmt(st["med"][i])}<br>'
                 f'nmad = {_fmt(st["nmad"][i])}'
                 for i in range(n_spec)]

        fig3.add_trace(
            go.Scatter(
                x=med_plot, y=nmad_plot,
                mode='markers',
                marker=dict(color=_ARM_SOLID[name], size=7, opacity=0.7),
                text=hover, hovertemplate='%{text}<extra></extra>',
                showlegend=False,
            ),
            row=1, col=col_i,
        )

        # x=0 reference line (ideal median)
        fig3.add_trace(
            go.Scatter(x=[0, 0], y=[0, y_hi],
                       mode='lines',
                       line=dict(color='grey', width=1, dash='dash'),
                       showlegend=False, hoverinfo='skip'),
            row=1, col=col_i,
        )

        # Horizontal reference at ensemble NMAD
        if ens is not None:
            fig3.add_trace(
                go.Scatter(x=[x_lo, x_hi],
                           y=[ens['nmad'], ens['nmad']],
                           mode='lines',
                           line=dict(color='grey', width=1, dash='dot'),
                           showlegend=False, hoverinfo='skip'),
                row=1, col=col_i,
            )

        # Annotation: ensemble NMAD and count of off-scale spectra
        ax_sfx   = '' if col_i == 1 else str(col_i)
        off_str  = f'<br>{n_off} off-scale' if n_off > 0 else ''
        ens_str  = f'ens. NMAD={_fmt(ens["nmad"])}<br>' if ens is not None else ''
        fig3.add_annotation(
            x=0.97, y=0.97,
            xref=f'x{ax_sfx} domain', yref=f'y{ax_sfx} domain',
            xanchor='right', yanchor='top',
            showarrow=False, align='right',
            font=dict(size=10, family='monospace'),
            bgcolor='rgba(255,255,255,0.75)',
            bordercolor='grey', borderwidth=1,
            text=f'{ens_str}N={n_spec}{off_str}',
        )

        # Apply percentile-based axis limits
        fig3.update_xaxes(range=[x_lo, x_hi], row=1, col=col_i)
        fig3.update_yaxes(range=[0, y_hi],    row=1, col=col_i)

    # ── Row 2: NMAD by spectrum number (log scale) ────────────
    # Log scale shows both the outlier tail and the bulk simultaneously.
    # Original spectrum order lets the user spot clusters of bad observations.
    spec_x = np.arange(n_spec)
    for name, _, _, _ in ARM_DEFS:
        st  = per_spec[name]
        ens = arm_stats[name]
        if st is None:
            continue

        hover_spec = [f'spec {i}<br>nmad = {_fmt(st["nmad"][i])}'
                      for i in range(n_spec)]

        fig3.add_trace(
            go.Scatter(
                x=spec_x, y=st['nmad'],
                mode='markers+lines',
                line=dict(color=_ARM_SOLID[name], width=1),
                marker=dict(color=_ARM_SOLID[name], size=5),
                name=name,
                text=hover_spec,
                hovertemplate='%{text}<extra></extra>',
            ),
            row=2, col=1,
        )

        # Horizontal reference at ensemble NMAD — extend past axis ends
        if ens is not None:
            fig3.add_trace(
                go.Scatter(
                    x=[-1, n_spec],
                    y=[ens['nmad'], ens['nmad']],
                    mode='lines',
                    line=dict(color=_ARM_SOLID[name], width=1, dash='dot'),
                    showlegend=False, hoverinfo='skip',
                ),
                row=2, col=1,
            )

    fig3.update_layout(
        height=750,
        template='simple_white',
        showlegend=True,
        legend=dict(x=0.99, xanchor='right', y=0.35, yanchor='top',
                    bgcolor='rgba(255,255,255,0.7)'),
    )
    fig3.update_xaxes(**_ax)
    fig3.update_yaxes(**_ax, exponentformat='e', showexponent='all')
    # Row 1 x-axes carry flux values — need scientific notation explicitly.
    for _ci in [1, 2, 3]:
        fig3.update_xaxes(title_text='Median residual flux',
                          exponentformat='e', showexponent='all',
                          row=1, col=_ci)
    fig3.update_xaxes(title_text='Spectrum number',
                      range=[-0.5, n_spec - 0.5], row=2, col=1)
    fig3.update_yaxes(title_text='NMAD', row=1, col=1)
    fig3.update_yaxes(title_text='NMAD  (log)', type='log', row=2, col=1)

    # ══════════════════════════════════════════════════════════
    # Write all three figures into a single self-contained HTML file
    # ══════════════════════════════════════════════════════════
    if outroot == '':
        outroot = f'{stem}_{int(wmin)}_{int(wmax)}'
    outfile = f'{outroot}.html'

    html_spec = fig.to_html(full_html=False, include_plotlyjs=True)
    html_hist = fig2.to_html(full_html=False, include_plotlyjs=False)
    html_qual = fig3.to_html(full_html=False, include_plotlyjs=False)

    with open(outfile, 'w') as fh:
        fh.write('<!DOCTYPE html>\n<html>\n<body>\n')
        fh.write(html_spec)
        fh.write('\n')
        fh.write(html_hist)
        fh.write('\n')
        fh.write(html_qual)
        fh.write('\n</body>\n</html>\n')

    print(f'Wrote {outfile}')


# ──────────────────────────────────────────────────────────────
# Command-line interface
# ──────────────────────────────────────────────────────────────

def steer(argv):
    filename        = ''
    positional_nums = []
    n_sample        = 20
    outroot         = ''

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-num':
            i += 1
            n_sample = int(argv[i])
        elif argv[i] == '-out':
            i += 1
            outroot = argv[i]
        elif argv[i][0] == '-':
            print('Error: cannot parse command line:', argv)
            return
        elif filename == '':
            filename = argv[i]
        else:
            try:
                positional_nums.append(float(argv[i]))
            except ValueError:
                print('Error: cannot parse command line:', argv)
                return
        i += 1

    if not filename:
        print(__doc__)
        return

    wmin = positional_nums[0] if len(positional_nums) > 0 else 3600.0
    wmax = positional_nums[1] if len(positional_nums) > 1 else 9800.0

    plot_eval(filename, wmin=wmin, wmax=wmax, n_sample=n_sample, outroot=outroot)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
