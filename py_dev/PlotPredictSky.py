#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Interactive Plotly comparison of PredictSky.py's predicted vs observed
    sky spectrum for one exposure, with the ensemble's 1-sigma prediction
    band and (optionally) sky-line-contaminated pixels recolored, and a
    residual sub-panel.

Command line usage (if any):

    usage: PlotPredictSky.py [-h] [-mask] [-mask_file PATH]
                             [-linear] [-title TITLE] [-outfile PATH]
                             [-py_progs_dir PATH]
                             fits_file

    where

    fits_file       is a FITS file written by PredictSky.py (WAVE,
                    FLUX_OBS, FLUX_PRED, FLUX_PRED_LO, FLUX_PRED_HI, COEF).

    -mask           recolor sky-line-contaminated pixels of the observed
                    spectrum light grey, using py_progs/PlotSpec.py's
                    get_sky_mask() (data/sky_mask.fits by default) --
                    off by default, matching PlotSpec.py/PlotSpec3.py's
                    own -mask convention.

    -mask_file PATH override the default sky mask FITS file (implies
                    -mask).

    -linear         use a linear y-axis instead of the default log axis
                    (sky flux spans several orders of magnitude between
                    faint continuum and bright OH lines, so log is the
                    more useful default for a full-spectrum overview).

    -title TITLE    plot title (default: the input filename).

    -outfile PATH   output HTML path (default:
                    Overview_Plot/<stem>.predictsky.html).

    -py_progs_dir PATH
                    path to the lvm_ksl repo's py_progs/ directory, which
                    supplies PlotSpec.get_sky_mask and GetSkyCont.
                    _interp_mask_to_wave (default: ~/SDSS/lvm_ksl/py_progs).

Description:

    Top panel: FLUX_OBS (unbroken black line, sky-line-masked pixels
    redrawn in light grey on top if -mask is given -- the same
    "recolor over an unbroken spectrum" convention as PlotSpec.py's -mask),
    FLUX_PRED (firebrick), and the ensemble's FLUX_PRED_LO..FLUX_PRED_HI
    band as a shaded region.  Bottom panel: FLUX_OBS - FLUX_PRED residual,
    with a zero reference line, sharing the wavelength axis.

Notes::

    Log-axis handling follows this project's established rules: the axis
    ceiling is set from a percentile (not the raw max, which is dominated
    by whichever OH line happens to be brightest) and flux is clipped to a
    small positive floor (not left as NaN/negative) so the traces stay
    continuous rather than showing gaps at zero-crossings.

History::

    260902  ksl  Coding begun.
    260904  ksl  Promoted from niv/ to py_dev/. Switched every option
        from double-dash (--mask-file) to single-dash (-mask_file),
        matching py_progs/'s convention. PY_PROGS_DIR is now CLI-
        overridable (-py_progs_dir) instead of hardcoded, matching
        EvalFluxResiduals.py's -py_progs_dir precedent.

'''

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
import plotly.graph_objects as go
from plotly.subplots import make_subplots

DEFAULT_PY_PROGS_DIR = Path('~/SDSS/lvm_ksl/py_progs').expanduser()

_pre = argparse.ArgumentParser(add_help=False)
_pre.add_argument('-py_progs_dir', default=str(DEFAULT_PY_PROGS_DIR))
_pre_args, _ = _pre.parse_known_args()

sys.path.insert(0, _pre_args.py_progs_dir)
from PlotSpec import get_sky_mask  # noqa: E402
from GetSkyCont import _interp_mask_to_wave  # noqa: E402

LOG_FLOOR = 1e-16
CEIL_PERCENTILE = 99.0


def plot_predict_sky(fits_file, mask=False, mask_file=None, log_y=True, title=None):
    '''
    Build the interactive figure for one PredictSky.py output file.

    Returns
    -------
    plotly.graph_objects.Figure
    '''
    with fits.open(fits_file) as hdul:
        wave = np.asarray(hdul['WAVE'].data, dtype=float)
        obs = np.asarray(hdul['FLUX_OBS'].data, dtype=float)
        pred = np.asarray(hdul['FLUX_PRED'].data, dtype=float)
        pred_lo = np.asarray(hdul['FLUX_PRED_LO'].data, dtype=float)
        pred_hi = np.asarray(hdul['FLUX_PRED_HI'].data, dtype=float)

    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True, row_heights=[0.72, 0.28],
        vertical_spacing=0.06,
    )

    def _clip(arr):
        return np.clip(arr, LOG_FLOOR, None) if log_y else arr

    fig.add_trace(go.Scatter(
        x=np.concatenate([wave, wave[::-1]]),
        y=np.concatenate([_clip(pred_hi), _clip(pred_lo)[::-1]]),
        fill='toself', fillcolor='rgba(214,39,40,0.15)',
        line=dict(color='rgba(255,255,255,0)'),
        name='Predicted 1σ', hoverinfo='skip',
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=wave, y=_clip(obs), mode='lines',
        line=dict(color='black', width=1), name='Observed',
    ), row=1, col=1)

    if mask:
        sky_mask = get_sky_mask(mask_file)
        if sky_mask is not None:
            mask_wave, mask_bool = sky_mask
            clean = _interp_mask_to_wave(mask_wave, mask_bool, wave)
            masked_obs = np.where(~clean, obs, np.nan)
            fig.add_trace(go.Scatter(
                x=wave, y=_clip(masked_obs), mode='lines',
                line=dict(color='lightgrey', width=1.5),
                name='Observed (sky-line masked)',
            ), row=1, col=1)
        else:
            print(f'-mask given but no sky mask found '
                  f'({mask_file or "data/sky_mask.fits"}); skipping overlay.')

    fig.add_trace(go.Scatter(
        x=wave, y=_clip(pred), mode='lines',
        line=dict(color='firebrick', width=1), name='Predicted',
    ), row=1, col=1)

    resid = obs - pred
    fig.add_trace(go.Scatter(
        x=wave, y=resid, mode='lines',
        line=dict(color='steelblue', width=1), name='Observed − Predicted',
        showlegend=False,
    ), row=2, col=1)
    fig.add_hline(y=0.0, line=dict(color='grey', width=1, dash='dot'), row=2, col=1)

    # Plotly's default tick formatting abbreviates large/small numbers with
    # SI-style suffixes (e.g. "6k" for 6000, "50f" for 50e-15) -- disabled
    # on both axes so wavelength reads as plain "6000" and the residual
    # axis reads in explicit scientific notation instead of femto-prefix.
    fig.update_xaxes(tickformat='.0f')
    fig.update_xaxes(title_text='Wavelength (Å)', row=2, col=1)
    fig.update_yaxes(title_text='Flux (erg s⁻¹ cm⁻² Å⁻¹)',
                     row=1, col=1)
    fig.update_yaxes(title_text='Residual (erg s⁻¹ cm⁻² Å⁻¹)',
                     exponentformat='e', row=2, col=1)

    if log_y:
        finite_pos = np.concatenate([obs[obs > 0], pred[pred > 0]])
        ceil = float(np.nanpercentile(finite_pos, CEIL_PERCENTILE)) if finite_pos.size else 1.0
        fig.update_yaxes(type='log', exponentformat='e',
                         range=[np.log10(LOG_FLOOR), np.log10(3.0 * ceil)],
                         row=1, col=1)

    fig.update_layout(
        template='plotly_white',
        title=title or Path(fits_file).name,
        legend=dict(orientation='h', yanchor='bottom', y=1.02),
        height=700,
    )
    return fig


def main():
    p = argparse.ArgumentParser(
        parents=[_pre],
        description="Interactive Plotly comparison of PredictSky.py's "
                    "predicted vs observed sky spectrum.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('fits_file', help='PredictSky.py output FITS file')
    p.add_argument('-mask', action='store_true',
                   help='recolor sky-line-contaminated pixels light grey')
    p.add_argument('-mask_file', default=None, dest='mask_file',
                   help='override the default sky mask FITS file (implies -mask)')
    p.add_argument('-linear', action='store_true',
                   help='use a linear y-axis instead of the default log axis')
    p.add_argument('-title', default=None, help='plot title')
    p.add_argument('-outfile', default=None, help='output HTML path')
    args = p.parse_args()

    mask = args.mask or (args.mask_file is not None)
    fig = plot_predict_sky(
        args.fits_file, mask=mask, mask_file=args.mask_file,
        log_y=not args.linear, title=args.title,
    )

    outfile = args.outfile
    if outfile is None:
        outdir = Path('Overview_Plot')
        outdir.mkdir(exist_ok=True)
        outfile = str(outdir / f'{Path(args.fits_file).stem}.predictsky.html')
    fig.write_html(outfile, include_plotlyjs=True)
    print('Wrote', outfile)


if __name__ == '__main__':
    main()
