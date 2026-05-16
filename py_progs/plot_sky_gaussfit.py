#!/usr/bin/env python
# coding: utf-8

'''
                   Space Telescope Science Institute

Synopsis:

Display sky_gaussfit.py fit results as sky maps of wavelength, flux, and
FWHM for each fitted line.  Up to 6 lines are shown per page (one row of
three panels per line: centroid wavelength, integrated flux, FWHM), with
colors scaled to the 5th-95th percentile of each quantity.  Pages are
saved as PNG files to the directory Figs_gaussfit_sky/.


Command line usage::

    plot_sky_gaussfit.py [-out root] [-lines line1,line2,...] [-np npage]
                         filename [filename ...]

    where

    -out root    Root name for output PNG files.  Default: stem of the
                 first input filename.
    -lines list  Comma-separated list of line names to plot (e.g.
                 ha,sii_a,sky6300).  Default: all lines inferred from
                 flux_ columns in the table, in column order.
    -np npage    Number of lines per page (default 6, max 6).
    -s size      Scatter marker size in points^2 (default 30).
    filename     One or more ASCII tables written by sky_gaussfit.py.
                 If multiple files are given they are stacked into one
                 combined table before plotting.


Description:

For each group of up to 6 lines the script creates one figure with rows
of three panels (wavelength, flux, FWHM).  Each panel is a scatter plot
of fiber RA vs Dec with the quantity shown as color.  Colorbars reflect
the 5th-95th percentile range so outliers do not dominate the scale.

Output filenames follow the pattern::

    Figs_gaussfit_sky/<root>_p01.png
    Figs_gaussfit_sky/<root>_p02.png
    ...

Primary routines:

    plot_one   - scatter plot of one quantity on one Axes object
    plot_line  - fill one row (wave, flux, fwhm) for a single line
    plot_page  - create and save one multi-row figure
    plot_all   - iterate over all lines in groups, calling plot_page


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
LINES_PER_PAGE = 6


def plot_one(ax, xtable, var, marker_size=30):
    '''
    Scatter plot of xtable[var] vs RA/Dec on ax.
    Color range is clipped to the 5th-95th percentile of finite values.
    '''
    if var not in xtable.colnames:
        ax.text(0.5, 0.5, '%s\nnot in table' % var,
                transform=ax.transAxes, ha='center', va='center', fontsize=8)
        return
    col = np.array(xtable[var], dtype=float)
    finite = col[np.isfinite(col)]
    if len(finite) == 0:
        ax.text(0.5, 0.5, 'no finite data', transform=ax.transAxes,
                ha='center', va='center', fontsize=8)
        return
    vmin = np.percentile(finite, 5)
    vmax = np.percentile(finite, 95)
    sc = ax.scatter(xtable['ra'], xtable['dec'], c=col, cmap='viridis',
                    vmin=vmin, vmax=vmax, s=marker_size, linewidths=0, alpha=1)
    plt.colorbar(sc, ax=ax)
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')


def plot_line(axes_row, xtable, line, marker_size=30):
    '''
    Fill one row of three Axes with wave, flux, and fwhm maps for line.
    '''
    for ax, suffix, title in zip(axes_row,
                                  ['wave', 'flux', 'fwhm'],
                                  ['Wavelength (A)', 'Flux', 'FWHM (A)']):
        var = '%s_%s' % (suffix, line)
        plot_one(ax, xtable, var, marker_size)
        ax.set_title('%s  %s' % (line, title), fontsize=9)


def plot_page(xtable, lines, page_num, outroot, marker_size=30):
    '''
    Create one figure with one row per line (up to LINES_PER_PAGE rows of
    3 panels each) and save it to FIG_DIR.
    '''
    nrows = len(lines)
    fig, axes = plt.subplots(nrows, 3, figsize=(12, 4 * nrows),
                              squeeze=False)
    for i, line in enumerate(lines):
        plot_line(axes[i], xtable, line, marker_size)

    plt.tight_layout()

    os.makedirs(FIG_DIR, exist_ok=True)
    figname = os.path.join(FIG_DIR, '%s_p%02d.png' % (outroot, page_num))
    plt.savefig(figname, dpi=100)
    print('Saved %s' % figname)
    plt.close(fig)


def plot_all(xtable, lines=None, outroot='sky_gaussfit', nper=LINES_PER_PAGE,
             marker_size=30):
    '''
    Plot all lines in groups of nper, saving one PNG per group.

    If lines is None the line list is inferred from flux_ columns in xtable,
    preserving the column order (which matches the sky_gaussfit fit order:
    nebular lines first, then airglow).
    '''
    if lines is None:
        lines = [c[5:] for c in xtable.colnames if c.startswith('flux_')]

    nper = min(nper, LINES_PER_PAGE)
    pages = [lines[i:i + nper] for i in range(0, len(lines), nper)]
    for page_num, page_lines in enumerate(pages, start=1):
        plot_page(xtable, page_lines, page_num, outroot, marker_size)


def steer(argv):
    outroot = ''
    lines = None
    filenames = []
    nper = LINES_PER_PAGE
    marker_size = 30

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-out':
            i += 1
            outroot = argv[i]
        elif argv[i] == '-lines':
            i += 1
            lines = argv[i].split(',')
        elif argv[i] == '-np':
            i += 1
            nper = int(argv[i])
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

    tables = []
    for fname in filenames:
        try:
            t = ascii.read(fname)
            tables.append(t)
            print('Read %d rows from %s' % (len(t), fname))
        except Exception as e:
            print('Could not read %s: %s' % (fname, e))

    if not tables:
        print('No tables could be read.')
        return

    xtable = vstack(tables) if len(tables) > 1 else tables[0]

    if outroot == '':
        base = os.path.basename(filenames[0])
        outroot = os.path.splitext(base)[0]

    plot_all(xtable, lines=lines, outroot=outroot, nper=nper,
             marker_size=marker_size)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
