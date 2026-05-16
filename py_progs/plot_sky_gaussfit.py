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
    -neb         Plot only the nebular group (LINE_GROUPS entry 'neb').
    -lines list  Comma-separated list of line names, saved as _p01.png,
                 _p02.png, … in groups of -np.  Overrides LINE_GROUPS.
    -np npage    Lines per page when using -lines (default 6, max 6).
    -s size      Scatter marker size in points^2 (default 30).
    filename     One or more ASCII tables written by sky_gaussfit.py.
                 If multiple files are given they are stacked before plotting.


Description:

The default behaviour (no -neb, no -lines) iterates over every entry in
the LINE_GROUPS table defined near the top of this file and saves one PNG
per entry named Figs_gaussfit_sky/<root>_<suffix>.png.  Edit LINE_GROUPS
to control which lines appear on which page.

With -neb only the 'neb' entry in LINE_GROUPS is plotted.

With -lines a custom comma-separated list is plotted in pages of -np lines
named _p01.png, _p02.png, etc.

Output filenames (default mode)::

    Figs_gaussfit_sky/<root>_neb.png
    Figs_gaussfit_sky/<root>_sky1.png
    Figs_gaussfit_sky/<root>_sky2.png
    Figs_gaussfit_sky/<root>_sky3.png

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

# Named line groups: each entry is (suffix, [line, ...]).
# The default run plots every group and saves <root>_<suffix>.png.
# Edit this table to control which lines appear on which page.
LINE_GROUPS = [
    ('neb',  ['oi_a', 'ha', 'nii_b', 'sii_a', 'sii_b']),
    ('sky1', ['sky5577', 'sky6300', 'sky6363', 'sky6533', 'sky6553', 'sky6577']),
    ('sky2', ['sky6912', 'sky6923', 'sky6939', 'sky7358', 'sky7392', 'sky7914']),
    ('sky3', ['sky8344', 'sky8399', 'sky8827', 'sky8988', 'sky9552', 'sky9719']),
]


# Per-column color range: fraction of the median used as the ±half-width.
# wave/fwhm: 30 km/s → delta_lambda = median * 30/3e5
# flux:      5 % of median
VRANGE_FRAC = {'wave': 10.0 / 3e5, 'fwhm': 100.0 / 3e5, 'flux': 0.02}

# Format string for the median value shown in each subplot title.
MEDIAN_FMT  = {'wave': '%.3f A', 'fwhm': '%.3f A', 'flux': '%.3e'}


def plot_one(ax, xtable, var, marker_size=30, subtract_median=False,
             vrange_frac=None):
    '''
    Scatter plot of xtable[var] vs RA/Dec on ax.
    If subtract_median is True the column median is removed before plotting
    and returned.  If vrange_frac is also given the color limits are set to
    ±(|median| * vrange_frac); otherwise the 5th-95th percentile is used.
    Returns the median, or None if the column is absent or all-NaN.
    '''
    if var not in xtable.colnames:
        ax.text(0.5, 0.5, '%s\nnot in table' % var,
                transform=ax.transAxes, ha='center', va='center', fontsize=8)
        return None
    col = np.array(xtable[var], dtype=float)
    finite = col[np.isfinite(col)]
    if len(finite) == 0:
        ax.text(0.5, 0.5, 'no finite data', transform=ax.transAxes,
                ha='center', va='center', fontsize=8)
        return None
    median_val = None
    if subtract_median:
        median_val = np.median(finite)
        col = col - median_val
        finite = finite - median_val
    if vrange_frac is not None and median_val is not None:
        delta = abs(median_val) * vrange_frac
        vmin, vmax = -delta, delta
    else:
        vmin = np.percentile(finite, 5)
        vmax = np.percentile(finite, 95)
    sc = ax.scatter(xtable['ra'], xtable['dec'], c=col, cmap='viridis',
                    vmin=vmin, vmax=vmax, s=marker_size, linewidths=0, alpha=1)
    plt.colorbar(sc, ax=ax)
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    return median_val


def plot_line(axes_row, xtable, line, marker_size=30):
    '''
    Fill one row of three Axes with wave, flux, and fwhm maps for line.
    Each panel shows quantity - median; color range is ±(|median|*VRANGE_FRAC).
    The median value and its units appear in each subplot title.
    '''
    for ax, suffix, label in zip(axes_row,
                                  ['wave', 'flux', 'fwhm'],
                                  ['Wavelength (A)', 'Flux', 'FWHM (A)']):
        var = '%s_%s' % (suffix, line)
        median_val = plot_one(ax, xtable, var, marker_size,
                              subtract_median=True,
                              vrange_frac=VRANGE_FRAC.get(suffix))
        if median_val is not None:
            median_str = MEDIAN_FMT.get(suffix, '%.4g') % median_val
            ax.set_title('%s  %s  [median=%s]' % (line, label, median_str),
                         fontsize=9)
        else:
            ax.set_title('%s  %s' % (line, label), fontsize=9)


def plot_page(xtable, lines, outroot, suffix, marker_size=30, title=''):
    '''
    Create one figure with one row per line (up to LINES_PER_PAGE rows of
    3 panels each) and save it to FIG_DIR/<outroot>_<suffix>.png.
    title is shown as a suptitle above all panels.
    '''
    nrows = len(lines)
    fig, axes = plt.subplots(nrows, 3, figsize=(16, 4 * nrows),
                              squeeze=False)
    for i, line in enumerate(lines):
        plot_line(axes[i], xtable, line, marker_size)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    if title:
        plt.suptitle(title, fontsize=16, y=0.99)

    os.makedirs(FIG_DIR, exist_ok=True)
    figname = os.path.join(FIG_DIR, '%s_%s.png' % (outroot, suffix))
    plt.savefig(figname, dpi=100, bbox_inches='tight')
    print('Saved %s' % figname)
    plt.close(fig)


def plot_all(xtable, outroot='sky_gaussfit', marker_size=30, title=''):
    '''
    Plot every group defined in LINE_GROUPS, saving one PNG per group.
    '''
    for suffix, lines in LINE_GROUPS:
        plot_page(xtable, lines, outroot, suffix, marker_size, title)


def plot_custom(xtable, lines, outroot='sky_gaussfit', nper=LINES_PER_PAGE,
                marker_size=30, title=''):
    '''
    Plot a custom line list in pages of nper, saved as _p01.png, _p02.png, etc.
    '''
    nper = min(nper, LINES_PER_PAGE)
    pages = [lines[i:i + nper] for i in range(0, len(lines), nper)]
    for page_num, page_lines in enumerate(pages, start=1):
        plot_page(xtable, page_lines, outroot, 'p%02d' % page_num,
                  marker_size, title)


def steer(argv):
    outroot = ''
    lines = None
    filenames = []
    nper = LINES_PER_PAGE
    marker_size = 30
    neb_mode = False

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-out':
            i += 1
            outroot = argv[i]
        elif argv[i] == '-neb':
            neb_mode = True
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

    # Suptitle: basename of each input file, extension stripped
    title = '  '.join(
        os.path.splitext(os.path.basename(f))[0] for f in filenames
    )

    if neb_mode:
        suffix, neb_lines = next(g for g in LINE_GROUPS if g[0] == 'neb')
        plot_page(xtable, neb_lines, outroot, suffix, marker_size, title)
    elif lines is not None:
        plot_custom(xtable, lines, outroot, nper, marker_size, title)
    else:
        plot_all(xtable, outroot, marker_size, title)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
