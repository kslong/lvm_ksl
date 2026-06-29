#!/usr/bin/env python
# coding: utf-8
'''
                    Space Telescope Science Institute

Synopsis:

    Evaluate continuum fits produced by GetSkyCont.py by plotting flux,
    continuum, components, and residuals over a chosen wavelength window.
    Reads a FITS file produced by GetSkyCont.py and writes an interactive
    Plotly HTML file.

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

    Creates a four-panel interactive HTML plot from a skycont_*.fits file:

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

Primary routines:

    plot_eval   build and write the Plotly HTML evaluation plot.

Notes:

    Requires plotly (pip install plotly).
    Wavelength window defaults to the full LVM range (3600-9800 A) when
    wmin and wmax are not given.

History:

    260628  ksl  Written.
    260628  ksl  Redesigned to four panels: Flux, Total Cont, Components, Residual.
'''

import sys
import numpy as np
from pathlib import Path
from astropy.io import fits
import plotly.graph_objects as go
from plotly.subplots import make_subplots


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
    starts = np.where(change ==  1)[0] + 1
    ends   = np.where(change == -1)[0] + 1
    if flag[0]:
        starts = np.r_[0, starts]
    if flag[-1]:
        ends = np.r_[ends, len(wave)]
    return [(float(wave[s]), float(wave[min(e, len(wave) - 1)]))
            for s, e in zip(starts, ends)]



# ──────────────────────────────────────────────────────────────
# Main plot function
# ──────────────────────────────────────────────────────────────

def plot_eval(filename, wmin=3600.0, wmax=9800.0, n_sample=20, outroot=''):
    '''Build a four-panel Plotly evaluation plot from a GetSkyCont output FITS file.

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
    # Build pixel index: wavelength window + downsample to ~2000 display points.
    # This lets each 2D array be loaded and sliced in one step, avoiding large
    # float64 intermediates (5 arrays × N_spec × 12401 before reduction).
    hdul     = fits.open(filename)
    wave_all = hdul['WAVE'].data.astype(float)
    mask_all = hdul['MASK'].data.astype(bool)

    wmask   = (wave_all >= wmin) & (wave_all <= wmax)
    ds      = max(1, int(wmask.sum()) // 2000)
    pix_idx = np.where(wmask)[0][::ds]

    wave  = wave_all[pix_idx]
    clean = mask_all[pix_idx]
    del wave_all, mask_all

    flux    = hdul['FLUX'].data[:,    pix_idx].astype(float)
    cont    = hdul['CONT'].data[:,    pix_idx].astype(float)
    moon    = hdul['MOON'].data[:,    pix_idx].astype(float)
    diffuse = hdul['DIFFUSE'].data[:, pix_idx].astype(float)
    resid   = hdul['RESID'].data[:,   pix_idx].astype(float)
    hdul.close()

    if flux.ndim == 1:
        flux    = flux[np.newaxis, :]
        cont    = cont[np.newaxis, :]
        moon    = moon[np.newaxis, :]
        diffuse = diffuse[np.newaxis, :]
        resid   = resid[np.newaxis, :]

    n_spec = len(flux)

    rng        = np.random.default_rng(42)
    sample_idx = (rng.choice(n_spec, size=min(n_sample, n_spec), replace=False)
                  if n_sample > 0 else [])

    # Clip to 1e-18 for log-scale display (keeps traces contiguous, no NaN artefacts).
    _floor    = 1e-18
    flux_l    = np.maximum(flux,    _floor)
    cont_l    = np.maximum(cont,    _floor)
    moon_l    = np.maximum(moon,    _floor)
    diffuse_l = np.maximum(diffuse, _floor)
    cont_med  = np.median(cont_l, axis=0)

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
    # Sit the line just below the panel 1 floor (log_lo from _log_range).
    _f_lo, _ = _log_range(flux)
    fit_y_val = 10 ** (_f_lo + 0.25)
    fit_y     = np.where(clean, fit_y_val, np.nan)
    fig.add_trace(
        go.Scatter(x=wave, y=fit_y, mode='lines',
                   line=dict(color='rgb(0,150,0)', width=3),
                   name='Fit pixels', showlegend=True, legend='legend'),
        row=1, col=1,
    )

    # ── Grey shading for masked regions (all panels) ───────────
    # Build all shapes at once and pass to update_layout in a single call.
    # add_vrect() goes through Plotly's schema validation on every call, making
    # it O(n_regions * n_panels) slow; batch dicts bypass that overhead.
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
    # ── Layout ─────────────────────────────────────────────────
    stem = Path(filename).stem
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

    # Shared range for panels 2 and 3: derived from CONT (= MOON + DIFFUSE),
    # so its 99th pct always exceeds either component.
    # 99th pct avoids the handful of extreme outlier pixels that inflate .max().
    _comp_range = list(_log_range(cont))

    fig.update_yaxes(title_text='Flux',       type='log',                              row=1, col=1)
    fig.update_yaxes(title_text='Continuum',  type='log', range=_comp_range,           row=2, col=1)
    fig.update_yaxes(title_text='Components', type='log', range=_comp_range,           row=3, col=1)
    fig.update_yaxes(title_text='Residual',                                             row=4, col=1)

    if outroot == '':
        outroot = f'{stem}_{int(wmin)}_{int(wmax)}'
    outfile = f'{outroot}.html'
    fig.write_html(outfile)
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
