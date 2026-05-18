#!/usr/bin/env python
# coding: utf-8

'''
                   Space Telescope Science Institute

Synopsis:

Display sky_gaussfit.py fit results for all 18 sky (airglow) lines using
one page per fitted quantity.  Three pages are produced:

    <root>_wave.png  — centroid wavelength residuals (wave - median)
    <root>_flux.png  — flux residuals (flux - median)
    <root>_fwhm.png  — FWHM residuals (fwhm - median)

Each page is a 6-row × 3-column grid of scatter plots (fiber RA vs Dec
colored by the residual), one panel per sky line.  Subplot titles give
the line name, median value, and MAD-based standard deviation.


Command line usage::

    plot_sky_gaussfit.py [-out root] [-s size] filename [filename ...]

    where

    -out root    Root name for output PNG files when a single input file is
                 given.  Ignored when multiple files are supplied (each file
                 uses its own stem as the root).
    -s size      Scatter marker size in points^2 (default 30).
    filename     One or more ASCII tables written by sky_gaussfit.py.
                 Each file is plotted independently; no stacking is done.


Description:

Each input file produces its own three output PNGs.  Color limits for each
quantity are set to ±(median × VRANGE_FRAC), giving a fixed, comparable
scale across all sky lines on a page.  Edit VRANGE_FRAC and FWHM_VRANGE_KMS
near the top of the file to adjust.

Primary routines:

    plot_quantity  - create and save one full-page figure for one quantity
    plot_all       - call plot_quantity for wave, flux, and fwhm


History:

260516 ksl Coding begun

'''

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import ascii


FIG_DIR = 'Figs_gaussfit_sky'

# The 18 airglow lines in wavelength order (must match sky_gaussfit output).
SKY_LINES = [
    'sky5577', 'sky6300', 'sky6363', 'sky6533', 'sky6553', 'sky6577',
    'sky6912', 'sky6923', 'sky6939', 'sky7358', 'sky7392', 'sky7914',
    'sky8344', 'sky8399', 'sky8827', 'sky8988', 'sky9552', 'sky9719',
]

# Central wavelengths (A) for each sky line — used to convert FWHM to km/s.
SKY_LINE_WAVE = {
    'sky5577': 5577.34,  'sky6300': 6300.31,  'sky6363': 6363.78,
    'sky6533': 6533.04,  'sky6553': 6553.0,   'sky6577': 6577.2,
    'sky6912': 6912.62,  'sky6923': 6923.22,  'sky6939': 6939.52,
    'sky7358': 7358.68,  'sky7392': 7392.21,  'sky7914': 7913.72,
    'sky8344': 8344.61,  'sky8399': 8399.18,  'sky8827': 8827.11,
    'sky8988': 8988.38,  'sky9552': 9552.55,  'sky9719': 9719.84,
}

# Grid layout for 18 panels.
NROWS, NCOLS = 6, 3

# Per-quantity color range: ±(|median| × VRANGE_FRAC).
# wave: Angstrom shift equivalent to the given km/s.
# flux: fractional variation of median.
# fwhm: handled separately in velocity space (see FWHM_VRANGE_KMS).
VRANGE_FRAC = {'wave': 5.0 / 3e5, 'flux': 0.02}

# Absolute ±half-width for FWHM color limits (km/s).
# FWHM_obs^2 = FWHM_inst^2 + sigma^2  (all in km/s); this sets the display range.
FWHM_VRANGE_KMS = 30.0

# Title format for median and MAD-std values.
MEDIAN_FMT  = {'wave': '%.3f', 'fwhm': '%.1f', 'flux': '%.3e'}

# Axis labels for each quantity.
QUANTITY_LABEL = {
    'wave': r'$\lambda$ residual (A)',
    'flux': 'Flux residual',
    'fwhm': 'FWHM residual (km/s)',
}

# Human-readable color-limit description for the suptitle of each page.
VRANGE_LABEL = {
    'wave': u'±5 km/s',
    'fwhm': u'±%.0f km/s' % FWHM_VRANGE_KMS,
    'flux': u'±2%',
}


def plot_panel(ax, xtable, var, quantity, marker_size=30, line_wave=None):
    '''
    Draw one scatter panel for xtable[var] on ax.
    Subtracts the median, applies fixed color limits, and returns
    (median_val, mad_std) or (None, None).

    For quantity='fwhm', pass line_wave (A) to convert FWHM from Angstroms
    to km/s via FWHM_kms = (FWHM_AA / line_wave) * c before plotting.
    Returned median_val and mad_std are then in km/s.
    '''
    if var not in xtable.colnames:
        ax.text(0.5, 0.5, 'not in table', transform=ax.transAxes,
                ha='center', va='center', fontsize=7)
        ax.set_visible(True)
        return None, None

    col = np.array(xtable[var], dtype=float)
    finite = col[np.isfinite(col)]
    if len(finite) == 0:
        ax.text(0.5, 0.5, 'no data', transform=ax.transAxes,
                ha='center', va='center', fontsize=7)
        return None, None

    # Convert FWHM (A) → km/s so all sky lines share the same velocity scale.
    if quantity == 'fwhm' and line_wave is not None:
        col = col / line_wave * 3e5
        finite = finite / line_wave * 3e5

    median_val = np.median(finite)
    col = col - median_val
    finite = finite - median_val
    mad_std = 1.4826 * np.median(np.abs(finite))

    if quantity == 'fwhm' and line_wave is not None:
        vmin, vmax = -FWHM_VRANGE_KMS, FWHM_VRANGE_KMS
    else:
        vrange = VRANGE_FRAC.get(quantity)
        if vrange is not None:
            delta = abs(median_val) * vrange
            vmin, vmax = -delta, delta
        else:
            vmin = np.percentile(finite, 5)
            vmax = np.percentile(finite, 95)

    sc = ax.scatter(xtable['ra'], xtable['dec'], c=col, cmap='viridis',
                    vmin=vmin, vmax=vmax, s=marker_size, linewidths=0, alpha=1)
    plt.colorbar(sc, ax=ax)
    ax.set_xlabel('RA', fontsize=7)
    ax.set_ylabel('Dec', fontsize=7)
    ax.tick_params(labelsize=7)
    return median_val, mad_std


def plot_quantity(xtable, quantity, outroot, title='', marker_size=30):
    '''
    Create one 6×3 grid page showing all SKY_LINES for the given quantity
    (wave, flux, or fwhm) and save it to FIG_DIR/<outroot>_<quantity>.png.
    '''
    fig, axes = plt.subplots(NROWS, NCOLS, figsize=(16, 4 * NROWS),
                              squeeze=False)

    for idx, line in enumerate(SKY_LINES):
        row, col_idx = divmod(idx, NCOLS)
        ax = axes[row][col_idx]
        var = '%s_%s' % (quantity, line)
        line_wave = SKY_LINE_WAVE.get(line) if quantity == 'fwhm' else None
        median_val, mad_std = plot_panel(ax, xtable, var, quantity, marker_size,
                                         line_wave=line_wave)
        if median_val is not None:
            fmt = MEDIAN_FMT.get(quantity, '%.4g')
            ax.set_title('%s  [med=%s  σ=%s]' % (line, fmt % median_val,
                                                   fmt % mad_std), fontsize=8)
        else:
            ax.set_title(line, fontsize=8)

    # Hide any unused panels (shouldn't happen with 18 lines and 6×3=18 cells)
    for idx in range(len(SKY_LINES), NROWS * NCOLS):
        row, col_idx = divmod(idx, NCOLS)
        axes[row][col_idx].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    if title:
        plt.suptitle('%s — %s  [color: %s]' % (
            title, QUANTITY_LABEL[quantity], VRANGE_LABEL.get(quantity, '')),
            fontsize=16, y=0.99)

    os.makedirs(FIG_DIR, exist_ok=True)
    figname = os.path.join(FIG_DIR, '%s_%s.png' % (outroot, quantity))
    plt.savefig(figname, dpi=100, bbox_inches='tight')
    print('Saved %s' % figname)
    plt.close(fig)


def plot_all(xtable, outroot='sky_gaussfit2', title='', marker_size=30):
    '''
    Produce wave, flux, and fwhm pages for all 18 sky lines.
    '''
    for quantity in ('wave', 'flux', 'fwhm'):
        plot_quantity(xtable, quantity, outroot, title, marker_size)


def steer(argv):
    outroot = ''
    filenames = []
    marker_size = 30

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-out':
            i += 1
            outroot = argv[i]
        elif argv[i] == '-s':
            i += 1
            marker_size = int(argv[i])
        elif argv[i][0] == '-':
            print('Unknown option: %s' % argv[i])
            return
        else:
            filenames.append(argv[i])
        i += 1

    if not filenames:
        print(__doc__)
        return

    for fname in filenames:
        try:
            xtable = ascii.read(fname)
            print('Read %d rows from %s' % (len(xtable), fname))
        except Exception as e:
            print('Could not read %s: %s' % (fname, e))
            continue

        stem = os.path.splitext(os.path.basename(fname))[0]
        root = outroot if (outroot and len(filenames) == 1) else stem
        plot_all(xtable, root, stem, marker_size)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
