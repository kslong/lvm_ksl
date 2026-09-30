#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Measure how flat the sky is across the LVM science IFU, and how
    reproducible any departures are from exposure to exposure, using the
    per-fiber airglow-line fluxes that sky_gaussfit.py fits on lvmCFrame
    (un-sky-subtracted) files.  Since the sky is essentially uniform over
    the IFU, fiber-to-fiber and spectrograph-to-spectrograph differences
    in sky-line flux measure the instrument (throughput / fiber-flat
    errors), independent of any sky-subtraction method.

Command line usage (if any):

    usage: SkyLineFlatness.py [-h] [-lines L] [-clip lo hi] [-out ROOT]
                              fitfile [fitfile ...]

    Arguments::

        fitfile     two or more sky_gaussfit.py output tables, one per
                    exposure, from lvmCFrame inputs (they must contain the
                    spectrographid and exposure columns sky_gaussfit.py
                    writes since 260930)

    Options::

        -h          print this help and exit
        -lines L    comma-separated sky lines to use (default: the nine
                    brightest, sky5577,sky6300,sky6363,sky7914,sky8344,
                    sky8399,sky8827,sky9552,sky9719)
        -clip lo hi normalized fluxes outside [lo, hi] are treated as
                    failed fits (default 0.8 1.2)
        -out ROOT   output root (default SkyLineFlatness)

    Example::

        sky_gaussfit.py -np 8 -out skyfit_CF_14405 lvmCFrame-00014405.fits
        ... (one per exposure)
        SkyLineFlatness.py skyfit_CF_*.txt

Description:

    For each line, every fiber's fitted flux is divided by that
    exposure's median over the IFU (so changes in sky brightness between
    exposures drop out); values outside -clip are dropped as failed fits.
    Exposures are ordered by exposure number (i.e. time).  Then:

    1. Spectrograph offsets: for each exposure and line, the median
       normalized flux of each spectrograph's fibers.  1.0 means that
       spectrograph sees the line exactly as bright as the IFU median.

    2. Fiber-pattern reproducibility: for each exposure and line, the
       Spearman rank correlation r between its normalized per-fiber
       pattern and the per-fiber median of all OTHER exposures
       (leave-one-out), and the MAD-based rms of (pattern - reference).
       Spearman is used so that outlier fibers do not dominate r,
       consistent with the MAD-based scatters in 3.

    3. For each line, the MAD-based fiber-to-fiber rms of a single
       exposure, the rms of the all-exposure median pattern, and their
       squared ratio -- the fraction of a single exposure's fiber-to-fiber
       variance that is the same in every exposure.

    Summaries are printed.  Figure: spectrograph ratios for the b/r-arm
    lines and for the z-arm lines (same y-scale), and r per exposure.

Outputs::

    <root>_spec.txt   line, exposure, mjd, sp1, sp2, sp3 (item 1)
    <root>_expo.txt   exposure, mjd, tile_id, line, r_vs_others,
                      rms_diff_pct (item 2)
    <root>.png        the figure

Primary routines:

    read_fits            read and normalize the fit tables
    spectrograph_ratios  item 1
    pattern_agreement    items 2 and 3
    plot_flatness        the figure
    do_all               run everything and write the outputs
    steer                command-line driver

Notes:

    Arm boundaries for grouping the lines are taken from the line
    wavelengths in sky_gaussfit.SKY_LINES (b/r below ARM_Z_START A,
    z above).  Runs in the ksl environment.

History::

    260930 ksl Moved into py_progs from a Vela project script
        (skyfit_consistency.py) used on 21 exposures on 260929: now reads
        exposure/spectrograph from the fit tables themselves, correlation
        switched from Pearson to Spearman, options for lines, clip range
        and output root.

'''

import sys
import re
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from sky_gaussfit import SKY_LINES


def _usage_from_doc(doc):
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc

_USAGE = _usage_from_doc(__doc__)

DEFAULT_LINES = ['sky5577', 'sky6300', 'sky6363', 'sky7914', 'sky8344', 'sky8399',
                 'sky8827', 'sky9552', 'sky9719']
LINE_WAVE = {name: center for name, center, wmin, wmax in SKY_LINES}
ARM_Z_START = 7550.0   # A: lines redder than this are z-arm
SPECS = (1, 2, 3)


def rstd(a, axis=None):
    '''MAD-based robust standard deviation, ignoring NaNs'''
    med = np.nanmedian(a, axis=axis, keepdims=True)
    return 1.4826 * np.nanmedian(np.abs(a - med), axis=axis)


def read_fits(files, lines, clip=(0.8, 1.2)):
    '''
    Read sky_gaussfit.py tables and normalize each line by its per-exposure
    median.

    Returns:
        dict with exps (exposure numbers, time order), info (exposure ->
        (mjd, tile_id)), fibers (common fiberids), spec (their
        spectrographid), cube (line -> array (nexp, nfiber) of normalized
        flux, NaN where missing or clipped); or None on failure
    '''
    tabs = {}
    for f in files:
        try:
            t = ascii.read(f)
        except Exception as e:
            print('Error: could not read %s (%s)' % (f, e))
            return None
        missing = [c for c in ['fiberid', 'spectrographid', 'exposure'] + ['flux_' + l for l in lines]
                   if c not in t.colnames]
        if missing:
            print('Error: %s lacks column(s) %s; rerun sky_gaussfit.py (260930 or later) on the CFrame'
                  % (f, ', '.join(missing)))
            return None
        e = int(t['exposure'][0])
        if e in tabs:
            print('Error: exposure %d appears in more than one input file' % e)
            return None
        tabs[e] = t
    if len(tabs) < 2:
        print('Error: need at least two exposures')
        return None

    exps = sorted(tabs)
    info = {e: (int(tabs[e]['mjd'][0]) if 'mjd' in tabs[e].colnames else -1,
                int(tabs[e]['tile_id'][0]) if 'tile_id' in tabs[e].colnames else -1) for e in exps}
    fibers = np.array(sorted(set.intersection(*[set(tabs[e]['fiberid']) for e in exps])))
    spec_of = dict(zip(tabs[exps[0]]['fiberid'], tabs[exps[0]]['spectrographid']))
    spec = np.array([spec_of[i] for i in fibers])

    cube = {}
    for l in lines:
        arr = np.full((len(exps), len(fibers)), np.nan)
        for j, e in enumerate(exps):
            t = tabs[e]
            v = np.array(t['flux_' + l], float)
            v[~np.isfinite(v) | (v <= 0)] = np.nan
            v = v / np.nanmedian(v)
            v[(v < clip[0]) | (v > clip[1])] = np.nan   # failed/wild fits
            lookup = dict(zip(t['fiberid'], v))
            arr[j] = [lookup.get(i, np.nan) for i in fibers]
        cube[l] = arr
    return dict(exps=exps, info=info, fibers=fibers, spec=spec, cube=cube)


def spectrograph_ratios(d, lines):
    '''Per exposure and line, median normalized flux in each spectrograph (item 1).'''
    rows = []
    for l in lines:
        for j, e in enumerate(d['exps']):
            row = d['cube'][l][j]
            rows.append([l, e, d['info'][e][0]] + [np.nanmedian(row[d['spec'] == k]) for k in SPECS])
    tab = Table(rows=rows, names=['line', 'exposure', 'mjd'] + ['sp%d' % k for k in SPECS])
    for k in SPECS:
        tab['sp%d' % k].format = '.4f'
    return tab


def pattern_agreement(d, lines):
    '''
    Leave-one-out Spearman r and rms difference per exposure and line
    (item 2), and per-line reproducible fraction (item 3).

    Returns:
        expo_tab, summary (line -> (rms single exposure, rms median pattern,
        reproducible fraction))
    '''
    rows = []
    summary = {}
    for l in lines:
        c = d['cube'][l]
        ref_all = np.nanmedian(c, axis=0)
        rms_one = np.nanmedian(rstd(c, axis=1))
        summary[l] = (rms_one, rstd(ref_all), rstd(ref_all) ** 2 / np.nanmedian(rstd(c, axis=1) ** 2))
        for j, e in enumerate(d['exps']):
            ref = np.nanmedian(np.delete(c, j, axis=0), axis=0)
            ok = np.isfinite(c[j]) & np.isfinite(ref)
            r = spearmanr(c[j][ok], ref[ok]).correlation
            rows.append([e, d['info'][e][0], d['info'][e][1], l, r, 100 * rstd(c[j][ok] - ref[ok])])
    tab = Table(rows=rows, names=['exposure', 'mjd', 'tile_id', 'line', 'r_vs_others', 'rms_diff_pct'])
    tab['r_vs_others'].format = '.3f'
    tab['rms_diff_pct'].format = '.2f'
    return tab, summary


def plot_flatness(d, lines, spec_tab, expo_tab, outname):
    '''Spectrograph ratios for b/r and z lines, and r per exposure.'''
    exps = d['exps']
    x = np.arange(len(exps))
    labels = [str(e) for e in exps]
    line_col = dict(zip(lines, plt.cm.tab10(np.arange(len(lines)) % 10)))
    sp_marker = {1: 'o', 2: 's', 3: '^'}
    groups = [('b/r-arm sky lines', [l for l in lines if LINE_WAVE.get(l, 0) < ARM_Z_START]),
              ('z-arm sky lines', [l for l in lines if LINE_WAVE.get(l, 0) >= ARM_Z_START])]
    groups = [g for g in groups if g[1]]

    fig, axes = plt.subplots(len(groups) + 1, 1, figsize=(11, 4.3 * (len(groups) + 1)))
    for ax, (arm, group) in zip(axes[:-1], groups):
        for k, m in sp_marker.items():
            for l in group:
                s = spec_tab[spec_tab['line'] == l]
                ax.plot(x, s['sp%d' % k], m, color=line_col[l], ms=5)
        ax.axhline(1, color='k', lw=0.5)
        ax.set_xticks(x, labels, rotation=60, fontsize=8)
        ax.set_ylabel('median flux in spectrograph /\nmedian flux in whole IFU')
        ax.set_title('%s: each spectrograph\'s median sky-line flux relative to the whole IFU, per exposure\n'
                     '(1.0 = same as IFU; each point = one line in one exposure)' % arm, fontsize=10)
        handles = [Line2D([], [], ls='', marker='o', color=line_col[l], label=l) for l in group]
        handles += [Line2D([], [], ls='', marker=m, color='k', mfc='none', label='sp%d' % k)
                    for k, m in sp_marker.items()]
        ax.legend(handles=handles, fontsize=7, ncol=3, loc='upper right')
    ymin = min(a.get_ylim()[0] for a in axes[:-1])
    ymax = max(a.get_ylim()[1] for a in axes[:-1])
    for a in axes[:-1]:
        a.set_ylim(ymin, ymax)   # same y-scale so the arms can be compared directly

    ax = axes[-1]
    for l in lines:
        s = expo_tab[expo_tab['line'] == l]
        ax.plot(x, s['r_vs_others'], 'o-', color=line_col[l], ms=4, lw=0.8, label=l)
    ax.set_xticks(x, labels, rotation=60, fontsize=8)
    ax.set_ylabel('Spearman r')
    ax.set_title('Fiber-by-fiber agreement: Spearman r between each exposure\'s per-fiber sky-line flux pattern\n'
                 'and the median pattern of all OTHER exposures (1 = identical pattern, 0 = unrelated)', fontsize=10)
    ax.legend(fontsize=7, ncol=3)
    fig.tight_layout()
    fig.savefig(outname, dpi=100)
    plt.close(fig)


def do_all(files, lines=None, clip=(0.8, 1.2), outroot='SkyLineFlatness'):
    '''Run the full analysis on a set of sky_gaussfit.py tables and write the outputs.'''
    lines = lines or DEFAULT_LINES
    d = read_fits(files, lines, clip)
    if d is None:
        return None

    spec_tab = spectrograph_ratios(d, lines)
    expo_tab, summary = pattern_agreement(d, lines)

    print('%d exposures, %d common fibers' % (len(d['exps']), len(d['fibers'])))
    print('\nPer-spectrograph ratio (median over fibers / IFU median): mean +- std across exposures')
    for l in lines:
        s = spec_tab[spec_tab['line'] == l]
        print('  %-8s  ' % l + '  '.join('sp%d=%.4f+-%.4f' % (k, np.nanmean(s['sp%d' % k]), np.nanstd(s['sp%d' % k]))
                                         for k in SPECS))
    print('\nFiber-to-fiber scatter per line:')
    for l in lines:
        rms_one, rms_ref, frac = summary[l]
        print('  %-8s  rms per exposure (median) %.2f%%   rms of median pattern %.2f%%   reproducible fraction %.2f'
              % (l, 100 * rms_one, 100 * rms_ref, frac))
    print('\nPer exposure: median over lines of Spearman r (pattern vs other exposures), rms diff %')
    for e in d['exps']:
        s = expo_tab[expo_tab['exposure'] == e]
        k = np.argmin(s['r_vs_others'])
        print('  %6d MJD %d tile %d:  r=%.3f  rms diff=%.2f%%   (min r %.3f at %s)'
              % (e, d['info'][e][0], d['info'][e][1], np.median(s['r_vs_others']),
                 np.median(s['rms_diff_pct']), s['r_vs_others'][k], s['line'][k]))

    spec_tab.write(outroot + '_spec.txt', format='ascii.fixed_width_two_line', overwrite=True)
    expo_tab.write(outroot + '_expo.txt', format='ascii.fixed_width_two_line', overwrite=True)
    plot_flatness(d, lines, spec_tab, expo_tab, outroot + '.png')
    print('\nWrote %s_spec.txt, %s_expo.txt, %s.png' % (outroot, outroot, outroot))
    return spec_tab, expo_tab


def steer(argv):
    files = []
    lines = None
    clip = (0.8, 1.2)
    outroot = 'SkyLineFlatness'
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '-h':
            print(_USAGE)
            return
        elif arg == '-lines':
            i += 1
            lines = argv[i].split(',')
        elif arg == '-clip':
            clip = (float(argv[i + 1]), float(argv[i + 2]))
            i += 2
        elif arg == '-out':
            i += 1
            outroot = argv[i]
        elif arg.startswith('-'):
            print('Error: unknown option "%s"' % arg)
            print(_USAGE)
            return
        else:
            files.append(arg)
        i += 1
    if not files:
        print(_USAGE)
        return
    if lines:
        bad = [l for l in lines if l not in LINE_WAVE]
        if bad:
            print('Error: unknown sky line(s) %s; available: %s' % (bad, ', '.join(LINE_WAVE)))
            return
    do_all(files, lines=lines, clip=clip, outroot=outroot)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
