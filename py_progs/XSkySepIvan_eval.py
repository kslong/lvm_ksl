#!/usr/bin/env python
# coding: utf-8
'''
                    Space Telescope Science Institute

Synopsis:

    Evaluate PALACE sky decomposition quality by plotting flux, residuals,
    and continuum over a chosen wavelength window.  Reads a FITS file
    produced by XSkySepIvan.py and writes an interactive Plotly HTML file.

Command line usage (if any):

    usage: XSkySepIvan_eval.py palace_file.fits [wmin wmax] [-num N] [-out outroot]

    Arguments:

    palace_file.fits  FITS file from XSkySepIvan.py (process_many or process_skyfile).
    wmin              minimum wavelength in Angstroms (default 3600).
    wmax              maximum wavelength in Angstroms (default 9800).

    Options:

    -num N            overlay N randomly-selected individual spectra on the band (default 20; 0 = band only).
    -out outroot      output filename root; default is <stem>_<wmin>_<wmax>.

Description:

    Creates a three-panel interactive HTML plot from a palace_*.fits file:

    Panel 1 (Flux): median spectrum with a shaded 10th/90th percentile band
    and N individual spectra overlaid in light grey.

    Panel 2 (Residual): same layout for the RESID extension.

    Panel 3 (Continuum): median CONT spectrum on a log-y scale with the
    10th/90th percentile band shaded.

    The HTML file can be opened in any browser for interactive exploration
    (zoom, pan, hover values).

Primary routines:

    plot_eval   build and write the Plotly HTML evaluation plot.

Notes:

    Requires plotly (pip install plotly).
    Wavelength window defaults to the full LVM range (3600-9800 A) when
    wmin and wmax are not given.

History:

    260628  ksl  Written.
'''

import sys
import numpy as np
from pathlib import Path
from astropy.io import fits
import plotly.graph_objects as go
from plotly.subplots import make_subplots


_USAGE = '''Usage: XSkySepIvan_eval.py palace_file.fits [wmin wmax] [-num N] [-out outroot]

Arguments:
  palace_file.fits  FITS file from XSkySepIvan.py
  wmin              minimum wavelength in Angstroms (default 3600)
  wmax              maximum wavelength in Angstroms (default 9800)

Options:
  -num N        number of individual spectra to overlay (default 20; 0 = band only)
  -out outroot  output filename root (default: <stem>_<wmin>_<wmax>)
'''


def _band_traces(wave, arr, name, color_band, color_med, show_band_legend, legend_ref='legend'):
    '''Return Plotly traces for a median line with a shaded 10/90-pct band.'''
    p10 = np.nanpercentile(arr, 10, axis=0)
    p90 = np.nanpercentile(arr, 90, axis=0)
    med = np.nanmedian(arr, axis=0)
    lower = go.Scatter(
        x=wave, y=p10,
        mode='lines', line=dict(width=0),
        showlegend=False, hoverinfo='skip',
    )
    upper = go.Scatter(
        x=wave, y=p90,
        mode='lines', line=dict(width=0),
        fill='tonexty', fillcolor=color_band,
        name='10–90 pct', showlegend=show_band_legend,
        legend=legend_ref, hoverinfo='skip',
    )
    median = go.Scatter(
        x=wave, y=med,
        mode='lines', line=dict(color=color_med, width=1.5),
        name=f'Median {name}', showlegend=True,
        legend=legend_ref,
    )
    return [lower, upper, median]


def _sample_traces(wave, arr, idx):
    '''Return semi-transparent traces for a random sample of spectra.'''
    n = len(idx)
    if n == 0:
        return []
    alpha = min(0.7, 3.0 / max(1.0, n ** 0.5))
    traces = []
    for i in idx:
        traces.append(go.Scatter(
            x=wave, y=arr[i],
            mode='lines',
            line=dict(color=f'rgba(120,120,120,{alpha:.2f})', width=0.6),
            showlegend=False, hoverinfo='skip',
        ))
    return traces


def plot_eval(filename, wmin=3600.0, wmax=9800.0, n_sample=20, outroot=''):
    '''Build a three-panel Plotly evaluation plot from a PALACE output FITS file.

    Parameters
    ----------
    filename : str or Path
        PALACE output FITS file produced by XSkySepIvan.py.
    wmin : float
        Minimum wavelength in Angstroms.
    wmax : float
        Maximum wavelength in Angstroms.
    n_sample : int
        Number of randomly-selected individual spectra to overlay on the
        median+band summary.  0 = band only.
    outroot : str
        Output filename root.  Default: <stem>_<wmin>_<wmax>.
    '''
    hdul = fits.open(filename)
    wave = hdul['WAVE'].data
    flux = hdul['FLUX'].data
    resid = hdul['RESID'].data
    cont = hdul['CONT'].data
    hdul.close()

    # Normalise to 2-D (N_spec, N_pix)
    if flux.ndim == 1:
        flux  = flux[np.newaxis, :]
        resid = resid[np.newaxis, :]
        cont  = cont[np.newaxis, :]

    mask  = (wave >= wmin) & (wave <= wmax)
    wave  = wave[mask]
    flux  = flux[:,  mask]
    resid = resid[:, mask]
    cont  = cont[:,  mask]

    n_spec = len(flux)

    rng        = np.random.default_rng(42)
    sample_idx = rng.choice(n_spec, size=min(n_sample, n_spec), replace=False) if n_sample > 0 else []

    fig = make_subplots(
        rows=3, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.06,
        subplot_titles=[
            f'Flux  ({n_spec} spectra)',
            'Residual',
            'Continuum  (log scale)',
        ],
    )

    panels = [
        (flux,  'Flux',  'rgba(31,119,180,0.2)',  'rgb(31,119,180)',  1, 'legend'),
        (resid, 'Resid', 'rgba(214,39,40,0.2)',   'rgb(214,39,40)',   2, 'legend2'),
        (cont,  'Cont',  'rgba(44,160,44,0.2)',   'rgb(44,160,44)',   3, 'legend3'),
    ]
    for arr, name, c_band, c_med, row, leg_ref in panels:
        for tr in _sample_traces(wave, arr, sample_idx):
            fig.add_trace(tr, row=row, col=1)
        for tr in _band_traces(wave, arr, name, c_band, c_med,
                               show_band_legend=(row == 1), legend_ref=leg_ref):
            fig.add_trace(tr, row=row, col=1)

    stem = Path(filename).stem

    _leg = dict(xanchor='right', yanchor='top', x=0.99, bgcolor='rgba(255,255,255,0.7)')
    fig.update_layout(
        title=f'{stem}  {int(wmin)}–{int(wmax)} Å',
        height=900,
        template='simple_white',
        legend  = dict(**_leg, y=0.99),
        legend2 = dict(**_leg, y=0.64),
        legend3 = dict(**_leg, y=0.29),
    )

    # Apply box borders and tick formatting to all axes in one sweep
    _ax = dict(showline=True, linewidth=1, linecolor='black',
               mirror=True, ticks='outside', ticklen=4, showticklabels=True)
    fig.update_xaxes(**_ax)
    fig.update_yaxes(**_ax, exponentformat='e', showexponent='all')

    # Per-panel x-axis settings
    fig.update_xaxes(range=[wmin, wmax])
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_xaxes(showticklabels=False, row=2, col=1)
    fig.update_xaxes(title_text='Wavelength (Å)', row=3, col=1)

    # Per-panel y-axis titles and log scale
    fig.update_yaxes(title_text='Flux',       row=1, col=1)
    fig.update_yaxes(title_text='Residual',   row=2, col=1)
    fig.update_yaxes(title_text='Continuum',  type='log', row=3, col=1)

    if outroot == '':
        outroot = f'{stem}_{int(wmin)}_{int(wmax)}'
    outfile = f'{outroot}.html'
    fig.write_html(outfile)
    print(f'Wrote {outfile}')


def steer(argv):
    filename        = ''
    positional_nums = []
    n_sample        = 20
    outroot         = ''

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(_USAGE)
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
        print(_USAGE)
        return

    wmin = positional_nums[0] if len(positional_nums) > 0 else 3600.0
    wmax = positional_nums[1] if len(positional_nums) > 1 else 9800.0

    plot_eval(filename, wmin=wmin, wmax=wmax, n_sample=n_sample, outroot=outroot)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
