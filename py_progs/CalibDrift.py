#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Show how far the instrument drifts from the long-term master
    calibrations over the nights they serve, from the corrections the
    DRP records in the headers of reduced science frames, and write an
    html report with plots and a table.

Command line usage (if any):

    usage: CalibDrift.py [-h] [-out outroot] dir [dir ...]
           CalibDrift.py [-h] [-out outroot] [-base dir] [-every N] -drpall file

    Arguments::

        dir          directories searched (recursively) for lvmSFrame
                     files; a version name such as 1.3.2 is taken to be
                     that version's reductions under $SAS_BASE_DIR

    Options::

        -h           print this help and exit
        -out outroot name of the html file (outroot.html) and of the
                     subdirectory with its plots and table (outroot/);
                     default CalibDrift
        -drpall FILE take the exposures from a drpall file (fits or h5)
                     instead of searching directories, using its
                     location column (relative to $SAS_BASE_DIR); its
                     calib_mjd, exptime, focus and K-mirror positions
                     and altitudes are joined to the output table
        -base DIR    resolve the drpall locations under DIR instead of
                     $SAS_BASE_DIR
        -every N     read only every Nth exposure of the drpall file,
                     for a quick trial

Description::

    For each science exposure the DRP applies two corrections on top of
    the master calibrations and records them in the lvmSFrame header:

        <CAM> FIBER SHIFT   shift in pixels applied to the master traces
                            at extraction (the median of the per-block
                            shifts of _fix_fiber_thermal_shifts)
        <CAM> WAVE SKYOFF   mean sky line offset in Angstrom removed from
                            the master wavelength solution
                            (shift_wave_skylines)

    These measure how far the instrument has moved from the masters on
    that night.  The epoch of the masters used is read from the COMMENT
    cards that list the calibration files.

    The report shows both corrections against MJD, with the master
    epochs marked, and against the number of nights since the epoch,
    with statistics per camera and per epoch.  A short summary is
    printed to the terminal.

Primary routines:

    doit
    steer

Notes::

    Only primary headers are read, so the ancillary files of the
    reductions are not needed.  Missing or unreadable files are counted,
    listed in outroot/unread_files.txt and skipped.

    The DRP limits the trace shifts to +-3 pixels; if fewer than half
    the fibers of a spectrograph give a sky line offset, no offset is
    applied, which shows up as SKYOFF = 0.

    Only numpy, astropy, matplotlib and xhtml.py are needed (not
    lvmdrp), so it can be run where the reductions are stored, e.g. at
    Utah:

        CalibDrift.py -out CalibDrift_1.3.2 -drpall drpall-1.3.2.fits

History::

    261101 ksl Coding begun, to follow the trace shifts and sky line
               offsets the DRP applies across the nights each set of
               master calibrations serves

'''

import os
import re
import sys
import time
import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xhtml


CAMERAS = ['b1', 'b2', 'b3', 'r1', 'r2', 'r3', 'z1', 'z2', 'z3']
SPEC_COLORS = {'1': 'C0', '2': 'C1', '3': 'C2'}
TRACE_RANGE = 3.0


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before the History: line, so -h stays short
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


_USAGE = _usage_from_doc(__doc__)


def find_sframes(dirs):
    '''
    The lvmSFrame files in a list of directories, searched recursively;
    a bare version name is taken to be that version's reductions
    '''
    files = []
    for d in dirs:
        d = os.path.expanduser(d)
        if not os.path.isdir(d):
            d = os.path.join(os.environ.get('SAS_BASE_DIR', ''), 'sdsswork/lvm/spectro/redux', d)
        if not os.path.isdir(d):
            print('CalibDrift: %s is not a directory' % d)
            continue
        for root, _, names in os.walk(d, followlinks=True):
            files += [os.path.join(root, n) for n in names if n.startswith('lvmSFrame-') and n.endswith('.fits')]
    return sorted(set(files))


DRPALL_COLUMNS = ['calib_mjd', 'tileid', 'exptime', 'sci_alt', 'sci_kmpos', 'sci_focpos', 'skye_focpos', 'skyw_focpos',
                  'sun_alt', 'moon_alt']


def read_drpall(filename, base='', every=1):
    '''
    The exposures of a drpall file (fits or h5): a table with the full
    path of each lvmSFrame and the columns in DRPALL_COLUMNS
    '''
    if filename.endswith('.h5'):
        import h5py
        with h5py.File(filename, 'r') as h:
            xtab = Table(h['summary/table'][:])
    else:
        xtab = Table.read(filename, hdu=1)
    for col in xtab.colnames:
        if xtab[col].dtype.kind == 'S':
            xtab[col] = np.char.decode(xtab[col]).astype(str)
    if every > 1:
        xtab = xtab[::every]
    base = base or os.environ.get('SAS_BASE_DIR', '')
    xtab['path'] = [os.path.join(base, str(loc).strip().lstrip('/')) for loc in xtab['location']]
    keep = ['expnum', 'path'] + [c for c in DRPALL_COLUMNS if c in xtab.colnames]
    return xtab[keep]


def calib_epoch(header, kind):
    '''
    The MJD of the calibration directory of the master of the given kind
    (e.g. mwave, mcentroids) listed in the COMMENT cards, or -1
    '''
    text = ''.join(str(c) for c in header.get('COMMENT', []))
    m = re.search(r'calib/(\d{5})/lvm-%s-' % kind, text)
    return int(m.group(1)) if m else -1


def read_sframe(filename):
    '''
    One row of the table: the exposure, the master epochs used and the
    corrections applied for each camera
    '''
    h = fits.getheader(filename)
    row = dict(file=os.path.basename(filename), drpver=h.get('DRPVER', ''), mjd=int(h.get('MJD', -1)),
               expnum=int(h.get('EXPOSURE', -1)), tileid=int(h.get('TILE_ID', -1) or -1),
               wave_epoch=calib_epoch(h, 'mwave'), trace_epoch=calib_epoch(h, 'mcentroids'))
    for cam in CAMERAS:
        row['shift_%s' % cam] = float(h.get('%s FIBER SHIFT' % cam.upper(), np.nan))
        row['skyoff_%s' % cam] = float(h.get('%s WAVE SKYOFF' % cam.upper(), np.nan))
    return row


def plot_quantity(xtab, prefix, ylabel, title, filename, xcol='mjd', epochs=None, limit=None):
    '''
    One panel per channel, the three spectrographs in different colors
    '''
    fig, axs = plt.subplots(3, 1, figsize=(14, 11), sharex=True)
    fig.suptitle(title)
    for ax, ch in zip(axs, 'brz'):
        for sp in '123':
            cam = ch + sp
            y = np.asarray(xtab['%s_%s' % (prefix, cam)], float)
            ax.plot(xtab[xcol], y, 'o', ms=4 if len(xtab) < 500 else 1.5, color=SPEC_COLORS[sp], label=cam)
        if epochs is not None:
            for e in epochs:
                ax.axvline(e, color='0.6', lw=0.8, ls='--')
        if limit is not None:
            for s in (-limit, limit):
                ax.axhline(s, color='r', lw=0.8)
        ax.axhline(0, color='k', lw=0.5)
        ax.set_ylabel('%s %s' % (ch, ylabel))
        ax.legend(loc='upper left', fontsize=8, ncol=3)
    axs[-1].set_xlabel('MJD (dashed: master epochs used)' if xcol == 'mjd' else 'nights since the epoch of the masters')
    fig.savefig(filename, dpi=80)
    plt.close(fig)
    return filename


def camera_stats(xtab, prefix):
    rows = []
    for cam in CAMERAS:
        y = np.asarray(xtab['%s_%s' % (prefix, cam)], float)
        y = y[np.isfinite(y)]
        if y.size == 0:
            rows.append([cam, 0] + ['--'] * 5)
            continue
        rows.append([cam, y.size, '%+.3f' % np.median(y), '%.3f' % np.std(y), '%+.3f' % y.min(), '%+.3f' % y.max(),
                     '%d' % np.sum(y == 0)])
    return rows


def doit(dirs, outroot='CalibDrift', drpall=None, base='', every=1):
    '''
    Collect the corrections recorded in the lvmSFrame headers found in
    dirs, or listed in a drpall file, and write outroot.html with its
    plots and tables in outroot/

    History::

        261101 ksl Coding begun
    '''
    dtab = None
    if drpall:
        dtab = read_drpall(drpall, base=base, every=every)
        files = list(dtab['path'])
        source = '%s (%d exposures%s)' % (os.path.basename(drpall), len(files), ', every %dth' % every if every > 1 else '')
    else:
        files = find_sframes(dirs)
        source = ', '.join(dirs)
    if not files:
        print('CalibDrift: no lvmSFrame files found in %s' % (drpall or dirs))
        return
    print('Reading %d lvmSFrame headers' % len(files))
    rows, missing, failed = [], [], []
    t0 = time.time()
    for i, f in enumerate(files):
        if i and i % 1000 == 0:
            print('  %d of %d read (%.0f s), %d missing, %d unreadable' % (i, len(files), time.time() - t0,
                                                                         len(missing), len(failed)))
        if not os.path.isfile(f):
            missing.append(f)
            continue
        try:
            rows.append(read_sframe(f))
        except Exception as e:
            failed.append('%s: %s' % (f, e))
    if not rows:
        print('CalibDrift: none of the %d files could be read (%d missing, %d unreadable)' % (len(files), len(missing), len(failed)))
        return
    xtab = Table(rows=rows)
    if dtab is not None:
        dtab = dtab[[c for c in dtab.colnames if c != 'path']]
        dtab.rename_column('tileid', 'drpall_tileid') if 'tileid' in dtab.colnames else None
        from astropy.table import join
        xtab = join(xtab, dtab, keys='expnum', join_type='left')
    xtab.sort(['mjd', 'expnum'])
    xtab['nights_since_epoch'] = xtab['mjd'] - xtab['wave_epoch']
    os.makedirs(outroot, exist_ok=True)
    tabfile = os.path.join(outroot, 'calib_drift.txt')
    xtab.write(tabfile, format='ascii.fixed_width_two_line', overwrite=True)
    problems = []
    if missing or failed:
        probfile = os.path.join(outroot, 'unread_files.txt')
        with open(probfile, 'w') as g:
            g.write('\n'.join(['missing %s' % f for f in missing] + ['unreadable %s' % f for f in failed]) + '\n')
        problems.append('%d of %d files were missing and %d unreadable; see %s' % (
            len(missing), len(files), len(failed), xhtml.link(os.path.basename(probfile), probfile)))
    if 'calib_mjd' in xtab.colnames:
        has = ~np.asarray(getattr(xtab['calib_mjd'], 'mask', np.zeros(len(xtab), bool)))
        differ = has & (np.asarray(xtab['calib_mjd'].filled(-1) if hasattr(xtab['calib_mjd'], 'filled')
                                   else xtab['calib_mjd']) != np.asarray(xtab['wave_epoch']))
        if differ.any():
            problems.append('%d exposures have a drpall calib_mjd different from the epoch of the mwave file in '
                            'their header' % differ.sum())

    epochs = sorted(set(int(e) for e in xtab['wave_epoch'] if e > 0))
    versions = sorted(set(str(v) for v in xtab['drpver']))
    figs = [plot_quantity(xtab, 'shift', 'trace shift (pix)', 'Trace shift applied at extraction (FIBER SHIFT)',
                          os.path.join(outroot, 'shift_mjd.png'), epochs=epochs, limit=TRACE_RANGE),
            plot_quantity(xtab, 'skyoff', 'sky line offset (A)', 'Sky line offset removed from the wavelengths (WAVE SKYOFF)',
                          os.path.join(outroot, 'skyoff_mjd.png'), epochs=epochs),
            plot_quantity(xtab, 'shift', 'trace shift (pix)', 'Trace shift against nights since the epoch of the masters',
                          os.path.join(outroot, 'shift_age.png'), xcol='nights_since_epoch'),
            plot_quantity(xtab, 'skyoff', 'sky line offset (A)', 'Sky line offset against nights since the epoch of the masters',
                          os.path.join(outroot, 'skyoff_age.png'), xcol='nights_since_epoch')]
    scols = ['camera', 'n', 'median', 'std', 'min', 'max', 'n exactly 0']
    shift_rows = camera_stats(xtab, 'shift')
    skyoff_rows = camera_stats(xtab, 'skyoff')

    # the exposures per master epoch
    erows = [['epoch of masters', 'n exposures', 'MJD range', 'nights since epoch', 'max |trace shift| (pix)',
              'max |sky offset| (A)']]
    for e in epochs:
        s = xtab['wave_epoch'] == e
        sh = np.abs(np.array([xtab['shift_%s' % c][s] for c in CAMERAS], float))
        so = np.abs(np.array([xtab['skyoff_%s' % c][s] for c in CAMERAS], float))
        erows.append([e, int(s.sum()), '%d-%d' % (xtab['mjd'][s].min(), xtab['mjd'][s].max()),
                      '%d-%d' % (xtab['nights_since_epoch'][s].min(), xtab['nights_since_epoch'][s].max()),
                      '%.3f' % np.nanmax(sh), '%.3f' % np.nanmax(so)])

    string = xhtml.begin('Drift from the master calibrations')
    string += xhtml.add_list(['%d lvmSFrame headers from %s' % (len(xtab), source),
                              'DRP versions: %s' % ', '.join(versions),
                              'MJD %d-%d; %d master epochs' % (xtab['mjd'].min(), xtab['mjd'].max(), len(epochs))]
                             + problems)
    string += xhtml.paragraph('''FIBER SHIFT is the shift applied to the master traces at extraction (the DRP allows
+-%.0f pixels, red lines); WAVE SKYOFF is the mean sky line offset removed from the master wavelength solution
of each spectrograph. Both are corrections the DRP makes because the instrument has moved since the
masters were made. A SKYOFF of exactly 0 means fewer than half the fibers gave a sky line offset and
none was applied.''' % TRACE_RANGE)
    string += xhtml.hline()
    string += xhtml.h2('Trace shift applied at extraction (pixels)')
    string += xhtml.table([scols] + shift_rows)
    string += xhtml.h2('Sky line offset removed from the wavelengths (Angstrom)')
    string += xhtml.table([scols] + skyoff_rows)
    string += xhtml.h2('By master epoch')
    string += xhtml.table(erows)
    string += xhtml.hline()
    for fig in figs:
        string += xhtml.image(fig, width=1100, height=864)
    string += xhtml.paragraph('Table: ' + xhtml.link(os.path.basename(tabfile), tabfile))
    string += xhtml.end()
    with open(outroot + '.html', 'w') as g:
        g.write(string)

    print('\nCalibDrift: %d exposures, MJD %d-%d, %d master epochs, versions %s' % (
        len(xtab), xtab['mjd'].min(), xtab['mjd'].max(), len(epochs), ', '.join(versions)))
    for p in problems:
        print('  ' + re.sub('<[^>]*>', '', p))
    for label, srows, unit in (('trace shift', shift_rows, 'pix'), ('sky offset', skyoff_rows, 'A')):
        print('  %s (%s): ' % (label, unit) + '; '.join('%s %s..%s' % (r[0], r[4], r[5]) for r in srows if r[1]))
    print('Report: %s.html (plots and table in %s/)' % (outroot, outroot))


def steer(argv):
    '''
    Just a steering routine
    '''
    outroot = 'CalibDrift'
    drpall = None
    base = ''
    every = 1
    dirs = []
    i = 1
    while i < len(argv):
        if argv[i][0:2] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-out':
            i += 1
            outroot = argv[i]
        elif argv[i] == '-drpall':
            i += 1
            drpall = argv[i]
        elif argv[i] == '-base':
            i += 1
            base = argv[i]
        elif argv[i] == '-every':
            i += 1
            every = int(argv[i])
        elif argv[i][0] == '-':
            print('Error: unknown option %s' % argv[i])
            print(_USAGE)
            return
        else:
            dirs.append(argv[i])
        i += 1
    if not dirs and not drpall:
        print(_USAGE)
        return
    doit(dirs, outroot=outroot, drpall=drpall, base=base, every=every)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
