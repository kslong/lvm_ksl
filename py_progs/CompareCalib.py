#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Compare two sets of LVM long-term master calibrations -- two
    processings of the same epoch with different versions of the DRP,
    or two calibration epochs -- and write an html report, with plots
    and tables, that one can scroll through.

Command line usage (if any):

    usage: CompareCalib.py [-h] [-out outroot] [-only bias,trace,wave,flat]
                           [-labels labA,labB] [-same] [-note text] A B

    Arguments::

        A, B         the two sets of masters; B is the reference, and
                     ratios and differences are A/B and A-B.  Each is
                     - a directory holding the masters, e.g.
                       ~/Projects/lvm_data/calib/60995; an epoch
                       directory of a version's reductions is searched
                       for its calib/ subdirectory
                     - an MJD, a directory in $LVM_MASTER_DIR, e.g. 60995
                     - version/MJD, the calib directory of that
                       version's reductions under $SAS_BASE_DIR, e.g.
                       ksl-wave1/60995

    Options::

        -h           print this help and exit
        -out outroot name of the html file (outroot.html) and of the
                     subdirectory with its plots and tables (outroot/);
                     default CompareCalib_<labA>_vs_<labB>
        -only LIST   compare only some products, a comma separated list
                     of bias, trace, wave and flat
        -labels A,B  short names for the two sets in plots and tables;
                     default derived from A and B
        -same        treat A and B as the same epoch whatever their MJD
                     keywords say (a master carries the MJD of its first
                     input frame, so e.g. two halves of one arc sequence
                     can carry different MJDs)
        -note text   a paragraph shown at the top of the report, e.g. to
                     say how A and B were made

Description::

    For each type of master present in both sets:

    - Provenance: DRPVER, COMMIT, MJD, NFRAMES and the input frames
      (FRAMES extension, where present) side by side, so differences
      in code or inputs show first.
    - Bias (lvm-mbias): pixel differences per camera, with the robust
      sigma and the number of pixels beyond 5 sigma per quadrant, since
      the amplifiers differ in noise.
    - Traces (lvm-mcentroids, -msigmas, -mmodel): centroid and sigma
      differences per camera, outlier fibers, the ratio of the fiber
      models and the quadrant gains of both sets, which set the scale
      of the models.
    - Wavelengths (lvm-mwave, -mlsf): wavelength differences in Angstrom
      and km/s per camera, outlier fibers, the ratio of the LSF FWHM and
      the fit RMS values recorded in the headers.
    - Twilight fiber flats (lvm-mfiberflat_twilight): the per-fiber ratio,
      split into spectrograph factors, a gradient across the science IFU
      and a fiber to fiber residual with lvmdrp.qa.fiberflats
      (bin_fiberflat, ratio_stats, decompose_ratio, _header_corrections),
      and the SKYCORR state, factors and gradient in both headers.

    Outlier fibers are found after subtracting the local trend of the
    difference, the running median over 15 neighbouring fibers, so that
    smooth changes along the slit and the spectrum (e.g. thermal shifts
    between epochs) are not flagged while a single bad fiber is.
    Wavelength outliers use 0.05 A inside the detector and 0.25 A within
    200 columns of its ends, where the polynomials extrapolate past the
    last arc lines.  For each outlier the wiggle, the scatter along the
    fiber of its offset from its neighbours, is given for both sets;
    for wavelengths the set with the larger wiggle is named as the
    suspect.  Trace differences are mostly constant offsets, which the
    wiggle cannot see, so no suspect is named for them.

    For traces and wavelengths the report also shows what the science
    reduction does not remove: the trace difference after one shift per
    block of 36 fibers (as _fix_fiber_thermal_shifts applies at
    extraction), and the wavelength difference after a slit-smooth
    per-fiber offset measured at the reference sky lines (as
    shift_wave_skylines applies).  That part propagates into the
    science spectra.

    A short summary is printed to the terminal.

Primary routines:

    doit
    steer

Notes::

    The html file is written to the current working directory and its
    plots and tables to the subdirectory outroot/.

    Twilight flats are compared on a common wavelength grid (50 A bins),
    so flats with different wavelength solutions can be compared; the
    other products are compared pixel by pixel.

    Twilight flats that have had the sky-line factor correction
    (drp calibrations fiberflat-corrections, FIBERFLAT SKYCORR = True)
    differ from uncorrected ones by the recorded spectrograph factors.
    In DRP 1.3.2 the IFU gradient recorded with them (GCOEFF1-4) was not
    applied to the flat (see FixFlatGradient.py), so it does not appear
    in the ratio.

    Requires lvmdrp (for the fiber flat routines and REF_SKYLINES).

History::

    261101 ksl Coding begun, to compare a rebuild of the 60995 masters
               with the DRP 1.3.0dev ones and with neighbouring epochs:
               provenance, bias, traces, wavelengths/LSF and twilight
               flats, outliers relative to the local trend, and what
               the science reduction does not remove
    261104 ksl Report explains its sources (image data; the header
               values used are marked), the ideas used, every table
               column and figure, and what the science reduction does
               not remove means; raw trace map on one fixed scale plus
               a map minus the camera median; colour bars per panel
               where the scale differs; ASCII-only html
    261104 ksl -same and -note options, for comparing the two halves
               of a split arc sequence

'''

import os
import re
import sys
import numpy as np
import bottleneck as bn
from astropy.io import fits
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import ndimage, interpolate
import xhtml
from lvmdrp.core.rss import RSS
from lvmdrp.core.constants import REF_SKYLINES
from lvmdrp.qa.fiberflats import (wavelength_edges, bin_fiberflat, ratio_stats, decompose_ratio,
                                  _header_corrections)


CAMERAS = ['b1', 'b2', 'b3', 'r1', 'r2', 'r3', 'z1', 'z2', 'z3']
CHANNELS = ['b', 'r', 'z']
PRODUCTS = ['bias', 'trace', 'wave', 'flat']
C_KMS = 299792.458
CALIB_TILE = '0011XX/11111'
# fibers whose largest difference exceeds these are listed as outliers
TRACE_LIMIT = 0.05   # pixels
WAVE_LIMIT = 0.05    # Angstrom, inside EDGE columns of the ends of the detector
WAVE_EDGE_LIMIT = 0.25   # Angstrom, within EDGE columns of the ends, where the polynomials are unconstrained
EDGE = 200
TRACE_RAW_SCALE = 1.0    # pixels, colour scale of the raw trace difference map
TRACE_MED_SCALE = 0.1    # pixels, colour scale of the trace difference minus the camera median
NWORST = 5


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before the History: line, so -h stays short
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


_USAGE = _usage_from_doc(__doc__)


def resolve(spec):
    '''
    Turn a command line specification of a set of masters into a
    directory and a short label

    Description:

        spec is a directory, an MJD (a directory in $LVM_MASTER_DIR),
        or version/MJD (the calib directory of that version's
        reductions under $SAS_BASE_DIR)

    History::

        261101 ksl Coding begun
    '''
    spec = os.path.expanduser(spec)
    if os.path.isdir(spec) and not has_masters(spec) and os.path.isdir(os.path.join(spec, 'calib')):
        # an epoch directory of a version's reductions; the masters are in its calib/
        spec = os.path.join(spec, 'calib')
    if os.path.isdir(spec):
        parts = os.path.normpath(os.path.abspath(spec)).split(os.sep)
        label = parts[-2] + '/' + parts[-1] if parts[-1] == 'calib' and len(parts) > 2 else parts[-1]
        if parts[-1] == 'calib' and 'redux' in parts:
            # .../redux/<version>/0011XX/11111/<mjd>/calib
            label = parts[parts.index('redux') + 1] + '/' + parts[-2]
        elif 'calib' in parts[:-1]:
            label = 'calib/' + parts[-1]
        return os.path.abspath(spec), label
    if spec.isdigit():
        master = os.environ.get('LVM_MASTER_DIR', '')
        return os.path.join(master, spec), 'calib/' + spec
    if '/' in spec:
        version, mjd = spec.split('/')[:2]
        base = os.environ.get('SAS_BASE_DIR', '')
        return os.path.join(base, 'sdsswork/lvm/spectro/redux', version, CALIB_TILE, mjd, 'calib'), spec
    return spec, spec


def has_masters(folder):
    return any(f.startswith('lvm-m') and f.endswith('.fits') for f in os.listdir(folder))


def local_trend(d, window=15, step=16):
    '''
    The local trend of a fibers x columns difference image: at every
    step-th column the running median over window neighbouring fibers
    (reflected at the ends of the slit), interpolated along the columns

    Notes:

        A single bad fiber does not affect the median of its window, so
        it stands out after the trend is subtracted, while smooth changes
        along the slit and the spectrum, such as thermal shifts between
        epochs, are removed
    '''
    nf, nc = d.shape
    cols = np.unique(np.r_[np.arange(0, nc, step), nc - 1])
    half = window // 2
    sub = np.pad(d[:, cols], ((half, half), (0, 0)), mode='reflect')
    win = np.lib.stride_tricks.sliding_window_view(sub, window, axis=0)
    med = bn.nanmedian(win, axis=-1)
    trend = np.full(d.shape, np.nan)
    x = np.arange(nc)
    for i in range(nf):
        ok = np.isfinite(med[i])
        if ok.sum() >= 2:
            trend[i] = np.interp(x, cols[ok], med[i][ok])
    return trend


def drp_trace_residual(d, block=36, columns=(450, 3050)):
    '''
    What is left of a fibers x columns centroid difference after the
    correction the science reduction applies at extraction: one constant
    shift per fiber block, measured over the columns where
    _fix_fiber_thermal_shifts cross-correlates (500-3000, 100 wide)
    '''
    res = np.full(d.shape, np.nan)
    shifts = []
    for i0 in range(0, d.shape[0], block):
        s = bn.nanmedian(d[i0:i0 + block, columns[0]:columns[1]])
        shifts.append(s)
        res[i0:i0 + block] = d[i0:i0 + block] - s
    return res, np.array(shifts)


def drp_wave_residual(d, wave, channel):
    '''
    What is left of a fibers x columns wavelength difference after the
    correction the science reduction applies with shift_wave_skylines: a
    constant offset per fiber, measured at the reference sky lines,
    median filtered over 8 fibers and smoothed with a spline along the
    slit, then subtracted from the zeroth coefficient

    Returns the residual and the per-fiber offset that was removed
    '''
    lines = np.asarray(REF_SKYLINES[channel])
    nf = d.shape[0]
    off = np.full(nf, np.nan)
    for f in range(nf):
        ok = np.isfinite(d[f]) & np.isfinite(wave[f])
        if ok.sum() < 2:
            continue
        w = wave[f][ok]
        inside = (lines > w.min()) & (lines < w.max())
        if inside.any():
            off[f] = np.nanmedian(np.interp(lines[inside], w, d[f][ok]))
    smooth = np.zeros(nf)
    good = np.isfinite(off)
    if good.sum() > 0.5 * nf:
        x = np.arange(nf)[good]
        med = ndimage.median_filter(off[good], 8)
        knots = np.linspace(x[len(x) // 10], x[-len(x) // 10], 10)
        try:
            tck = interpolate.splrep(x, med, task=-1, t=knots)
            smooth = interpolate.splev(np.arange(nf), tck)
        except Exception:
            smooth = np.interp(np.arange(nf), x, med)
    return d - smooth[:, None], smooth


def mpath(folder, kind, name):
    return os.path.join(folder, 'lvm-%s-%s.fits' % (kind, name))


def both_exist(a, b, kind, names):
    return [n for n in names if os.path.isfile(mpath(a, kind, n)) and os.path.isfile(mpath(b, kind, n))]


def camera_slitmap(path, camera):
    '''
    The rows of the SLITMAP of one camera, in the order of the fibers of
    a camera product, or None if there is no SLITMAP
    '''
    with fits.open(path) as h:
        if 'SLITMAP' not in h:
            return None
        sm = Table(h['SLITMAP'].data)
    return sm[np.asarray(sm['spectrographid']) == int(camera[1])]


def fiber_label(sm, i):
    if sm is None:
        return '%d' % i
    return '%d(%s)' % (i, str(sm['orig_ifulabel'][i]).strip())


def neighbour_wiggle(data, fiber):
    '''
    Scatter along the fiber of its offset from the mean of its two
    neighbours, after removing the median offset; fibers at the ends
    of the slit use their one neighbour
    '''
    n = data.shape[0]
    lo, hi = max(fiber - 1, 0), min(fiber + 1, n - 1)
    if lo == fiber:
        lo = hi
    if hi == fiber:
        hi = lo
    d = data[fiber] - 0.5 * (data[lo] + data[hi])
    d = d[np.isfinite(d)]
    if d.size == 0:
        return np.nan
    return float(np.std(d - np.median(d)))


def robust_sigma(x):
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def fmt(x, form='%.4f'):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return '--'
    try:
        return form % x
    except TypeError:
        return str(x)


def write_table(rows, colnames, outdir, name):
    '''
    Write a list of rows as an ascii table in outdir and return the
    file name relative to the html file
    '''
    xtab = Table(rows=[[str(v) for v in r] for r in rows], names=colnames) if rows else Table(names=colnames)
    filename = os.path.join(outdir, name)
    xtab.write(filename, format='ascii.fixed_width_two_line', overwrite=True)
    return filename


def camera_grid(title):
    fig, axs = plt.subplots(3, 3, figsize=(15, 10), sharex=True, sharey=True)
    fig.suptitle(title)
    return fig, axs


def save_grid(fig, axs, im, label, filename, ylabel='fiber'):
    if im is not None:
        fig.colorbar(im, ax=axs, label=label, shrink=0.8)
    for x in axs[-1]:
        x.set_xlabel('column')
    for x in axs[:, 0]:
        x.set_ylabel(ylabel)
    fig.savefig(filename, dpi=80)
    plt.close(fig)
    return filename


def compare_provenance(a, b, laba, labb):
    '''
    Header keywords and input frames of the first file of each kind of
    master in A and B

    Returns the rows of an html table and the list of files whose
    input frames differ
    '''
    first = [('mbias', 'b1'), ('mcentroids', 'b1'), ('mmodel', 'b1'), ('mwave', 'b1'), ('mlsf', 'b1'),
             ('mfiberflat_twilight', 'b')]
    keys = ['MJD', 'DRPVER', 'COMMIT', 'NFRAMES', 'EXPOSURE']
    rows = [['file', 'set'] + keys + ['FRAMES (expnums)']]
    notes = []
    for kind, name in first:
        frames = {}
        for lab, folder in ((laba, a), (labb, b)):
            p = mpath(folder, kind, name)
            if not os.path.isfile(p):
                rows.append(['lvm-%s-%s' % (kind, name), lab] + ['missing'] + [''] * len(keys))
                continue
            with fits.open(p) as h:
                hdr = h[0].header
                fr = None
                if 'FRAMES' in h:
                    fr = sorted(set(int(x) for x in h['FRAMES'].data['expnum']))
            frames[lab] = fr
            frtxt = 'not recorded' if fr is None else '%d: %d-%d' % (len(fr), fr[0], fr[-1]) if fr else 'empty'
            rows.append(['lvm-%s-%s' % (kind, name), lab] + [str(hdr.get(k, '--')) for k in keys] + [frtxt])
        if len(frames) == 2 and None not in frames.values() and frames[laba] != frames[labb]:
            only_a = sorted(set(frames[laba]) - set(frames[labb]))
            only_b = sorted(set(frames[labb]) - set(frames[laba]))
            notes.append('lvm-%s-%s: input frames differ; only in %s: %s; only in %s: %s' %
                         (kind, name, laba, only_a, labb, only_b))
    return rows, notes


def compare_bias(a, b, laba, labb, outdir):
    '''
    Pixel differences of the master biases
    '''
    cams = both_exist(a, b, 'mbias', CAMERAS)
    if not cams:
        return None
    rows = []
    fig, axs = camera_grid('Bias: %s - %s (e-)' % (laba, labb))
    im = None
    for cam, ax in zip(CAMERAS, axs.flat):
        ax.set_title(cam)
        if cam not in cams:
            ax.text(0.5, 0.5, 'missing', transform=ax.transAxes, ha='center')
            continue
        da, ha = fits.getdata(mpath(a, 'mbias', cam)), fits.getheader(mpath(a, 'mbias', cam))
        db, hb = fits.getdata(mpath(b, 'mbias', cam)), fits.getheader(mpath(b, 'mbias', cam))
        d = (da - db).astype(float)
        ok = np.isfinite(d)
        dmax = np.abs(d[ok]).max() if ok.any() else np.nan
        if np.array_equal(da, db, equal_nan=True):
            equal = 'exact'
        elif dmax < 1e-4:
            equal = 'rounding'
        else:
            equal = 'no'
        sig = robust_sigma(d[ok])
        # the amplifiers differ in noise, so sigma and the 5 sigma count are per quadrant
        ny, nx = d.shape
        qsig, nhigh = [], 0
        for rs in (slice(0, ny // 2), slice(ny // 2, ny)):
            for cs in (slice(0, nx // 2), slice(nx // 2, nx)):
                q = d[rs, cs][np.isfinite(d[rs, cs])]
                s = robust_sigma(q)
                qsig.append(s)
                if q.size:
                    nhigh += int(np.sum(np.abs(q - np.median(q)) > max(5 * s, 1e-4)))
        rows.append([cam, '%s/%s' % (ha.get('NFRAMES', '--'), hb.get('NFRAMES', '--')),
                     fmt(np.nanmedian(da), '%.3f'), fmt(np.nanmedian(db), '%.3f'), fmt(np.median(d[ok]), '%+.4f'),
                     fmt(sig, '%.4f'), ' '.join(fmt(s, '%.3f') for s in qsig), str(nhigh),
                     fmt(100.0 * nhigh / max(ok.sum(), 1), '%.2f'), fmt(dmax, '%.2e'), equal])
        lim = max(5 * sig, 1e-3) if np.isfinite(sig) else 1
        im = ax.imshow(d, aspect='auto', origin='lower', cmap='RdBu_r', vmin=-lim, vmax=lim, interpolation='nearest')
    colnames = ['camera', 'NFRAMES A/B', 'median A', 'median B', 'median diff', 'robust sigma diff',
                'robust sigma by quadrant (lower left, lower right, upper left, upper right)', 'n pixels >5 sigma',
                '% pixels >5 sigma', 'max |diff|', 'equal (exact, rounding <1e-4 e-, no)']
    fig_name = save_grid(fig, axs, im, 'e- (scale +-5 robust sigma per camera)', os.path.join(outdir, 'bias_diff.png'),
                         ylabel='row')
    tab_name = write_table(rows, colnames, outdir, 'bias.txt')
    nexact = sum(r[-1] == 'exact' for r in rows)
    nround = sum(r[-1] == 'rounding' for r in rows)
    summary = 'bias: %d cameras; %d exactly equal, %d equal to rounding (<1e-4 e-); robust sigma of the ' \
        'difference %s to %s e-; %d pixels beyond 5 sigma' % (
            len(rows), nexact, nround, fmt(min(float(r[5]) for r in rows), '%.4f'),
            fmt(max(float(r[5]) for r in rows), '%.4f'), sum(int(r[7]) for r in rows))
    return dict(rows=[colnames] + rows, figs=[fig_name], tables=[tab_name], summary=summary)


def gains(header, cam):
    return [header.get('%s AMP%d GAIN' % (cam.upper(), i), np.nan) for i in range(1, 5)]


def compare_traces(a, b, laba, labb, outdir):
    '''
    Differences of the fiber centroids, sigmas and models
    '''
    cams = both_exist(a, b, 'mcentroids', CAMERAS)
    if not cams:
        return None
    rows, outliers, drows = [], [], []
    figr, axr = camera_grid('Trace centroids: %s - %s (pixels)' % (laba, labb))
    figm, axm = camera_grid('Trace centroids: %s - %s minus the camera median (pixels)' % (laba, labb))
    figc, axc = camera_grid('Trace centroids: %s - %s minus the local trend (pixels)' % (laba, labb))
    figs, axs = camera_grid('Trace sigmas: %s - %s (pixels)' % (laba, labb))
    figf, axf = plt.subplots(3, 3, figsize=(15, 10), sharex=True)
    figf.suptitle('Largest |centroid difference - local trend| along each fiber: %s - %s' % (laba, labb))
    imc = ims = imr = imm = None
    for cam, a0, am, a1, a2, a3 in zip(CAMERAS, axr.flat, axm.flat, axc.flat, axs.flat, axf.flat):
        for x in (a0, am, a1, a2, a3):
            x.set_title(cam)
        if cam not in cams:
            continue
        ca, cb = fits.getdata(mpath(a, 'mcentroids', cam)), fits.getdata(mpath(b, 'mcentroids', cam))
        ha, hb = fits.getheader(mpath(a, 'mcentroids', cam)), fits.getheader(mpath(b, 'mcentroids', cam))
        sm = camera_slitmap(mpath(a, 'mcentroids', cam), cam)
        draw = (ca - cb).astype(float)
        dc = draw - local_trend(draw)
        dres, bshift = drp_trace_residual(draw)
        fmax = np.nanmax(np.abs(dc), axis=1)
        nout = int(np.sum(fmax > TRACE_LIMIT))
        worst = np.argsort(-np.nan_to_num(fmax))[:NWORST]
        ds = None
        if os.path.isfile(mpath(a, 'msigmas', cam)) and os.path.isfile(mpath(b, 'msigmas', cam)):
            sa, sb = fits.getdata(mpath(a, 'msigmas', cam)), fits.getdata(mpath(b, 'msigmas', cam))
            ds = (sa - sb).astype(float)
        mrel = mp99 = np.nan
        if os.path.isfile(mpath(a, 'mmodel', cam)) and os.path.isfile(mpath(b, 'mmodel', cam)):
            ma, mb = fits.getdata(mpath(a, 'mmodel', cam)), fits.getdata(mpath(b, 'mmodel', cam))
            s = mb > 1000
            r = np.abs(ma[s] / mb[s] - 1)
            mrel, mp99 = np.median(r), np.percentile(r, 99)
        ga, gb = gains(ha, cam), gains(hb, cam)
        drows.append([cam, fmt(np.nanmin(bshift), '%+.3f'), fmt(np.nanmax(bshift), '%+.3f'), fmt(np.nanstd(dres)),
                      fmt(np.nanpercentile(np.abs(dres), 99)), fmt(np.nanmax(np.abs(dres)), '%.3f'),
                      str(int(np.sum(np.nanmax(np.abs(dres), axis=1) > TRACE_LIMIT)))])
        rows.append([cam, fmt(np.nanmedian(draw), '%+.4f'), fmt(np.nanstd(draw)), fmt(np.nanstd(dc)),
                     fmt(np.nanpercentile(np.abs(dc), 99)), fmt(np.nanmax(np.abs(dc)), '%.3f'),
                     str(nout), ' '.join('%s:%.2f' % (fiber_label(sm, f), fmax[f]) for f in worst if fmax[f] > TRACE_LIMIT) or '--',
                     fmt(np.nanstd(ds) if ds is not None else np.nan), fmt(mrel, '%.2e'), fmt(mp99, '%.2e'),
                     ' '.join('%.2f/%.2f' % (x, y) for x, y in zip(ga, gb))])
        for f in np.where(fmax > TRACE_LIMIT)[0]:
            cols = np.where(np.abs(dc[f]) > TRACE_LIMIT)[0]
            wa, wb = neighbour_wiggle(ca, f), neighbour_wiggle(cb, f)
            outliers.append([cam, f, fiber_label(sm, f).split('(')[-1].rstrip(')') if sm is not None else '--',
                             str(sm['fibstatus'][f]) if sm is not None else '--', '%.3f' % fmax[f],
                             '%d (%d-%d)' % (len(cols), cols.min(), cols.max()),
                             '%.3f' % np.nanmax(np.abs(ds[f])) if ds is not None else '--', '%.3f' % wa, '%.3f' % wb])
        # one fixed scale for all cameras, so the colours mean the same everywhere
        imr = a0.imshow(draw, aspect='auto', origin='lower', cmap='RdBu_r', vmin=-TRACE_RAW_SCALE,
                        vmax=TRACE_RAW_SCALE, interpolation='nearest')
        dmed = np.nanmedian(draw)
        imm = am.imshow(draw - dmed, aspect='auto', origin='lower', cmap='RdBu_r', vmin=-TRACE_MED_SCALE,
                        vmax=TRACE_MED_SCALE, interpolation='nearest')
        am.set_title('%s: median %+.3f' % (cam, dmed))
        imc = a1.imshow(dc, aspect='auto', origin='lower', cmap='RdBu_r', vmin=-0.02, vmax=0.02, interpolation='nearest')
        if ds is not None:
            ims = a2.imshow(ds, aspect='auto', origin='lower', cmap='RdBu_r', vmin=-0.02, vmax=0.02, interpolation='nearest')
        a3.semilogy(np.arange(len(fmax)), np.maximum(fmax, 1e-5), '.', ms=2)
        a3.axhline(TRACE_LIMIT, color='r', lw=0.8)
    for x in axf[-1]:
        x.set_xlabel('fiber')
    for x in axf[:, 0]:
        x.set_ylabel('max |dcen - trend| (pix)')
    f0 = save_grid(figr, axr, imr, 'pixels', os.path.join(outdir, 'trace_centroid_raw.png'))
    fm = save_grid(figm, axm, imm, 'pixels', os.path.join(outdir, 'trace_centroid_minus_median.png'))
    f1 = save_grid(figc, axc, imc, 'pixels', os.path.join(outdir, 'trace_centroid_diff.png'))
    f2 = save_grid(figs, axs, ims, 'pixels', os.path.join(outdir, 'trace_sigma_diff.png'))
    figf.savefig(os.path.join(outdir, 'trace_fiber_max.png'), dpi=80)
    plt.close(figf)
    colnames = ['camera', 'median dcen', 'rms dcen', 'rms dcen-trend', 'p99 |dcen-trend|', 'max |dcen-trend|',
                'n fibers >%.2f' % TRACE_LIMIT, 'worst fibers (max |dcen-trend|)',
                'rms dsig', 'model median |A/B-1|', 'model p99 |A/B-1|', 'gains A/B (amp1-4)']
    ocols = ['camera', 'fiber', 'ifulabel', 'fibstatus', 'max |dcen-trend|', 'columns >limit', 'max |dsig|',
             'wiggle A', 'wiggle B']
    outliers.sort(key=lambda r: -float(r[4]))
    dcols = ['camera', 'block shift min', 'block shift max', 'rms left', 'p99 |left|', 'max |left|',
             'n fibers >%.2f left' % TRACE_LIMIT]
    t3 = write_table(drows, dcols, outdir, 'trace_drp_residual.txt')
    t1 = write_table(rows, colnames, outdir, 'trace.txt')
    t2 = write_table(outliers, ocols, outdir, 'trace_outliers.txt')
    drp_summary = 'traces after the per-block shift the science reduction applies: rms left %s to %s pix; ' \
        '%d fibers still off by more than %.2f pix' % (fmt(min(float(r[3]) for r in drows)), fmt(max(float(r[3]) for r in drows)),
                                                     sum(int(r[6]) for r in drows), TRACE_LIMIT)
    summary = 'traces: %d cameras; median offsets %s to %s pix; %d fibers depart from the local trend by more ' \
        'than %.2f pix; largest %s pix' % (len(rows), fmt(min(float(r[1]) for r in rows), '%+.4f'),
                                           fmt(max(float(r[1]) for r in rows), '%+.4f'), len(outliers), TRACE_LIMIT,
                                           fmt(max(float(r[5]) for r in rows), '%.3f'))
    return dict(rows=[colnames] + rows, outliers=[ocols] + outliers, drows=[dcols] + drows, drp_summary=drp_summary,
                figs=[f0, fm, f1, f2, os.path.join(outdir, 'trace_fiber_max.png')], tables=[t1, t2, t3], summary=summary)


def rms_keywords(header, cam):
    '''
    The fit RMS keywords of a wave master: the current ones, or the
    PIPE DISP RMS ones of older versions
    '''
    new = '%s WAVE DISP RMS MEDIAN' % cam
    if new in header:
        return ('%.3f/%.3f' % (header[new], header['%s WAVE DISP RMS MAX' % cam]),
                '%.3f/%.3f' % (header['%s WAVE FWHM RMS MEDIAN' % cam], header['%s WAVE FWHM RMS MAX' % cam]))
    if 'PIPE DISP RMS MEDIAN' in header:
        return ('--', 'PIPE %.3f/%.3f' % (header['PIPE DISP RMS MEDIAN'], header['PIPE DISP RMS MAX']))
    return ('--', '--')


def compare_wave(a, b, laba, labb, outdir):
    '''
    Differences of the wavelength solutions and the LSF
    '''
    cams = both_exist(a, b, 'mwave', CAMERAS)
    if not cams:
        return None
    rows, outliers, drows = [], [], []
    figd, axd = camera_grid('Wavelength: %s - %s after the sky-line offset the science reduction applies (Angstrom)'
                            % (laba, labb))
    imd = None
    figr, axr = camera_grid('Wavelength: %s - %s (Angstrom, scale per camera)' % (laba, labb))
    figw, axw = camera_grid('Wavelength: %s - %s minus the local trend (Angstrom)' % (laba, labb))
    figl, axl = camera_grid('LSF FWHM: %s / %s - 1' % (laba, labb))
    figf, axf = plt.subplots(3, 3, figsize=(15, 10), sharex=True)
    figf.suptitle('Largest |wavelength difference - local trend| along each fiber: %s - %s' % (laba, labb))
    imw = iml = None
    for cam, a0, a1, a2, a3 in zip(CAMERAS, axr.flat, axw.flat, axl.flat, axf.flat):
        for x in (a0, a1, a2, a3):
            x.set_title(cam)
        if cam not in cams:
            continue
        wa, wb = fits.getdata(mpath(a, 'mwave', cam)).astype(float), fits.getdata(mpath(b, 'mwave', cam)).astype(float)
        ha, hb = fits.getheader(mpath(a, 'mwave', cam)), fits.getheader(mpath(b, 'mwave', cam))
        sm = camera_slitmap(mpath(a, 'mcentroids', cam), cam) if os.path.isfile(mpath(a, 'mcentroids', cam)) else None
        bad = np.zeros(wa.shape, bool)
        for folder in (a, b):
            with fits.open(mpath(folder, 'mwave', cam)) as h:
                if 'BADPIX' in h:
                    bad |= h['BADPIX'].data.astype(bool)
        d = np.where(bad, np.nan, wa - wb)
        med = np.nanmedian(d)
        v = d / wb * C_KMS
        dd = d - local_trend(d)
        dres, soff = drp_wave_residual(d, wb, cam[0])
        blue = bn.nanmedian(dres[:, :EDGE], axis=1)
        red = bn.nanmedian(dres[:, -EDGE:], axis=1)
        drows.append([cam, fmt(np.nanmin(soff), '%+.4f'), fmt(np.nanmax(soff), '%+.4f'),
                      fmt(np.nanstd(dres[:, EDGE:-EDGE]), '%.4f'), fmt(np.nanpercentile(np.abs(dres[:, EDGE:-EDGE]), 99), '%.4f'),
                      fmt(np.nanmedian(blue), '%+.4f'), fmt(np.nanmedian(red), '%+.4f'),
                      fmt(np.nanmedian(blue) / np.nanmedian(wb[:, EDGE // 2]) * C_KMS, '%+.2f'),
                      fmt(np.nanmedian(red) / np.nanmedian(wb[:, -EDGE // 2]) * C_KMS, '%+.2f'),
                      str(int(np.sum(np.nanmax(np.abs(dres[:, EDGE:-EDGE]), axis=1) > WAVE_LIMIT)))])
        dpanel = axd.flat[CAMERAS.index(cam)]
        dpanel.set_title(cam)
        imd = dpanel.imshow(dres, aspect='auto', origin='lower', cmap='RdBu_r', vmin=-0.1, vmax=0.1, interpolation='nearest')
        fmax = np.nanmax(np.abs(dd), axis=1)
        imax = np.nanmax(np.abs(dd[:, EDGE:-EDGE]), axis=1)
        emax = np.fmax(np.nanmax(np.abs(dd[:, :EDGE]), axis=1), np.nanmax(np.abs(dd[:, -EDGE:]), axis=1))
        frms = np.sqrt(np.nanmean(dd ** 2, axis=1))
        flag = (imax > WAVE_LIMIT) | (emax > WAVE_EDGE_LIMIT)
        nout = int(np.sum(flag))
        worst = [f for f in np.argsort(-np.nan_to_num(fmax)) if flag[f]][:NWORST]
        lrel = lp99 = np.nan
        lr = None
        if os.path.isfile(mpath(a, 'mlsf', cam)) and os.path.isfile(mpath(b, 'mlsf', cam)):
            la, lb = fits.getdata(mpath(a, 'mlsf', cam)), fits.getdata(mpath(b, 'mlsf', cam))
            lr = np.where(bad, np.nan, la / lb - 1)
            lrel, lp99 = np.nanmedian(lr), np.nanpercentile(np.abs(lr), 99)
        ra, rb = rms_keywords(ha, cam), rms_keywords(hb, cam)
        rows.append([cam, fmt(med, '%+.4f'), fmt(np.nanmedian(v), '%+.2f'), fmt(np.nanstd(d), '%.4f'),
                     fmt(np.median(frms), '%.4f'),
                     fmt(np.nanpercentile(np.abs(dd[:, EDGE:-EDGE]), 99), '%.4f'), fmt(np.nanmax(imax), '%.3f'),
                     fmt(np.nanmax(emax), '%.3f'), str(nout),
                     ' '.join('%s:%.3f' % (fiber_label(sm, f), fmax[f]) for f in worst) or '--',
                     fmt(lrel, '%+.4f'), fmt(lp99, '%.4f'), ra[0], rb[0], ra[1], rb[1]])
        for f in np.where(flag)[0]:
            col = int(np.nanargmax(np.abs(dd[f])))
            wga, wgb = neighbour_wiggle(wa, f), neighbour_wiggle(wb, f)
            outliers.append([cam, f, fiber_label(sm, f).split('(')[-1].rstrip(')') if sm is not None else '--',
                             '%.3f' % imax[f], '%.3f' % emax[f], col, '%.3f' % frms[f], '%.3f' % wga, '%.3f' % wgb,
                             laba if wga > wgb else labb])
        lim = max(0.05, np.nanpercentile(np.abs(d - med), 99))
        im0 = a0.imshow(d - med, aspect='auto', origin='lower', cmap='RdBu_r', vmin=-lim, vmax=lim, interpolation='nearest')
        a0.set_title('%s: median %+.3f, scale +-%.2f' % (cam, med, lim))
        figr.colorbar(im0, ax=a0, label='Angstrom')   # each camera has its own scale, so its own colour bar
        imw = a1.imshow(dd, aspect='auto', origin='lower', cmap='RdBu_r', vmin=-0.05, vmax=0.05, interpolation='nearest')
        if lr is not None:
            iml = a2.imshow(lr, aspect='auto', origin='lower', cmap='RdBu_r', vmin=-0.05, vmax=0.05, interpolation='nearest')
        a3.semilogy(np.arange(len(fmax)), np.maximum(fmax, 1e-5), '.', ms=2)
        a3.axhline(WAVE_LIMIT, color='r', lw=0.8)
    for x in axf[-1]:
        x.set_xlabel('fiber')
    for x in axf[:, 0]:
        x.set_ylabel('max |dwave - trend| (A)')
    f0 = save_grid(figr, axr, None, 'Angstrom', os.path.join(outdir, 'wave_raw.png'))
    f1 = save_grid(figw, axw, imw, 'Angstrom', os.path.join(outdir, 'wave_diff.png'))
    f2 = save_grid(figl, axl, iml, 'relative', os.path.join(outdir, 'lsf_ratio.png'))
    figf.savefig(os.path.join(outdir, 'wave_fiber_max.png'), dpi=80)
    plt.close(figf)
    colnames = ['camera', 'median dwave (A)', 'median dv (km/s)', 'rms dwave (A)', 'median per-fiber rms about trend',
                'p99 |dwave-trend| inside', 'max |dwave-trend| inside', 'max |dwave-trend| edges',
                'n fibers >%.2f A inside or >%.2f A at edges' % (WAVE_LIMIT, WAVE_EDGE_LIMIT),
                'worst flagged fibers (max |dwave-trend|)', 'LSF median A/B-1',
                'LSF p99 |A/B-1|', 'wave RMS A (pix med/max)', 'wave RMS B', 'FWHM RMS A', 'FWHM RMS B']
    ocols = ['camera', 'fiber', 'ifulabel', 'max |dwave-trend| inside (A)', 'max at edges (A)', 'column of max',
             'rms about trend (A)', 'wiggle A (A)', 'wiggle B (A)',
             'suspect']
    outliers.sort(key=lambda r: (-float(r[3]), -float(r[4])))
    fd = save_grid(figd, axd, imd, 'Angstrom', os.path.join(outdir, 'wave_drp_residual.png'))
    dcols = ['camera', 'sky offset removed min (A)', 'max (A)', 'rms left inside (A)', 'p99 |left| inside (A)',
             'median left, blue end (A)', 'median left, red end (A)', 'blue end (km/s)', 'red end (km/s)',
             'n fibers >%.2f A left inside' % WAVE_LIMIT]
    t3 = write_table(drows, dcols, outdir, 'wave_drp_residual.txt')
    t1 = write_table(rows, colnames, outdir, 'wave.txt')
    t2 = write_table(outliers, ocols, outdir, 'wave_outliers.txt')
    ends = [abs(float(r[c])) for r in drows for c in (5, 6) if r[c] != '--']
    drp_summary = 'wavelengths after the sky-line offset the science reduction applies: rms left inside %s to %s A; ' \
        'largest median left at a detector end %s A; %d fibers still off by more than %.2f A inside' % (
            fmt(min(float(r[3]) for r in drows)), fmt(max(float(r[3]) for r in drows)), fmt(max(ends) if ends else np.nan),
            sum(int(r[9]) for r in drows), WAVE_LIMIT)
    summary = 'wavelengths: %d cameras; median offsets %s to %s A; %d fibers depart from the local trend by more ' \
        'than %.2f A inside or %.2f A within %d columns of the edges; largest inside %s A, at edges %s A' % (
        len(rows), fmt(min(float(r[1]) for r in rows), '%+.4f'), fmt(max(float(r[1]) for r in rows), '%+.4f'),
        len(outliers), WAVE_LIMIT, WAVE_EDGE_LIMIT, EDGE, fmt(max(float(r[6]) for r in rows), '%.3f'),
        fmt(max(float(r[7]) for r in rows), '%.3f'))
    return dict(rows=[colnames] + rows, outliers=[ocols] + outliers, drows=[dcols] + drows, drp_summary=drp_summary,
                dfigs=[fd], figs=[f0, f1, f2, os.path.join(outdir, 'wave_fiber_max.png')], tables=[t1, t2, t3],
                summary=summary)


def ifu_panel(ax, x, y, values, title, lim, center=1.0):
    sc = ax.scatter(x, y, c=values, s=9, cmap='RdBu_r', vmin=center - lim, vmax=center + lim, marker='h')
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    return sc


def compare_flats(a, b, laba, labb, outdir, kind='mfiberflat_twilight', bin_width=50.0):
    '''
    Ratio of the twilight fiber flats, split into spectrograph factors,
    a gradient across the science IFU and a fiber to fiber residual
    with the routines of lvmdrp.qa.fiberflats
    '''
    chans = both_exist(a, b, kind, CHANNELS)
    if not chans:
        return None
    rows, hrows = [], []
    fig, axs = plt.subplots(len(chans), 3, figsize=(15, 5 * len(chans)), squeeze=False)
    fig.suptitle('Twilight fiber flats on the science IFU: %s / %s per fiber' % (laba, labb))
    fig2, axs2 = plt.subplots(len(chans), 1, figsize=(15, 4 * len(chans)), squeeze=False)
    fig2.suptitle('Fiber flat ratio %s / %s per fiber, along the slit' % (laba, labb))
    for ich, ch in enumerate(chans):
        fa, fb = RSS.from_file(mpath(a, kind, ch)), RSS.from_file(mpath(b, kind, ch))
        edges = wavelength_edges(fb, bin_width=bin_width)
        ba, bb = bin_fiberflat(fa, edges), bin_fiberflat(fb, edges)
        sm = fb._slitmap
        tel = np.char.strip(np.asarray(sm['telescope']).astype(str))
        stat_a = np.asarray(fa._slitmap['fibstatus']) if 'fibstatus' in fa._slitmap.colnames else np.zeros(len(sm))
        stat_b = np.asarray(sm['fibstatus']) if 'fibstatus' in sm.colnames else np.zeros(len(sm))
        good = np.isin(tel, ('Sci', 'SkyE', 'SkyW')) & (stat_a == 0) & (stat_b == 0)
        with np.errstate(invalid='ignore', divide='ignore'):
            ratio = ba / bb
        fiber_ratio = bn.nanmedian(ratio, axis=1)
        fiber_ratio[~good] = np.nan
        stats = ratio_stats(ratio[good])
        parts = decompose_ratio(sm, fiber_ratio, good, deg=1)
        scale = parts['level'] * parts['model']
        with np.errstate(invalid='ignore', divide='ignore'):
            corr = fiber_ratio / scale
            corr_stats = ratio_stats((ratio / scale[:, None])[good])
        sci = parts['gradients'].get('Sci', {})
        fstd = np.nanstd(fiber_ratio[good] / parts['level'])
        rows.append([ch, stats['nvalues'], fmt(stats['median']), fmt(stats['sigma']), fmt(stats['frac_off'], '%.3f'),
                     fmt(fstd), fmt(parts['level']),
                     ' '.join(fmt(parts['offsets'].get(i, np.nan), '%+.4f') for i in (1, 2, 3)),
                     fmt(sci.get('gx', np.nan), '%+.4f'), fmt(sci.get('gy', np.nan), '%+.4f'), fmt(sci.get('amp', np.nan)),
                     fmt(np.nanstd(corr[good])), fmt(corr_stats['sigma']), fmt(corr_stats['frac_off'], '%.3f')])
        for lab, f in ((laba, fa), (labb, fb)):
            hc = _header_corrections(getattr(f, '_header', None), ch, f._slitmap)
            hdr = getattr(f, '_header', None)
            skycorr = hdr.get('%s FIBERFLAT SKYCORR' % ch, '--') if hdr is not None else '--'
            hrows.append([ch, lab, str(skycorr), ' '.join(fmt(hc['hdr_factor%d' % i], '%.4f') for i in (1, 2, 3)),
                          ' '.join(fmt(hc['hdr_offset_sp%d' % i], '%+.4f') for i in (1, 2, 3)),
                          fmt(hc['hdr_grad_x'], '%+.4f'), fmt(hc['hdr_grad_y'], '%+.4f'),
                          ' '.join(fmt(hc['hdr_gcoeff%d' % k], '%.4g') for k in (1, 2, 3, 4))])
        s = (tel == 'Sci') & good
        x, y = np.asarray(sm['xpmm'], float)[s], np.asarray(sm['ypmm'], float)[s]
        r = fiber_ratio[s] / parts['level']
        lim = max(3 * np.nanstd(r), 0.01)
        sc = ifu_panel(axs[ich, 0], x, y, r, '%s: ratio / median' % ch, lim)
        ifu_panel(axs[ich, 1], x, y, parts['model'][s], '%s: spectrograph factors x gradient' % ch, lim)
        ifu_panel(axs[ich, 2], x, y, corr[s], '%s: residual (fiber to fiber)' % ch, lim)
        fig.colorbar(sc, ax=axs[ich, :].tolist(), shrink=0.8)
        specid = np.asarray(sm['spectrographid']).astype(int)
        for sp, color in zip((1, 2, 3), ('C0', 'C1', 'C2')):
            m = good & (specid == sp)
            axs2[ich, 0].plot(np.where(m)[0], fiber_ratio[m] / parts['level'], '.', ms=2, color=color, label='sp%d' % sp)
        axs2[ich, 0].set_ylabel('%s ratio / median' % ch)
        axs2[ich, 0].set_ylim(1 - 4 * fstd, 1 + 4 * fstd)
        axs2[ich, 0].legend(loc='upper right', markerscale=4)
    axs2[-1, 0].set_xlabel('fiber (row of the flat)')
    f1 = os.path.join(outdir, 'flat_ifu.png')
    fig.savefig(f1, dpi=80)
    plt.close(fig)
    f2 = os.path.join(outdir, 'flat_slit.png')
    fig2.savefig(f2, dpi=80)
    plt.close(fig2)
    colnames = ['channel', 'n values', 'median ratio', 'robust sigma', 'frac off >1%', 'per-fiber std', 'level',
                'spectrograph offsets sp1 sp2 sp3', 'IFU gradient gx', 'gy', 'amplitude', 'per-fiber std after model',
                'robust sigma after model', 'frac off >1% after model']
    hcols = ['channel', 'set', 'SKYCORR', 'header factors sp1 sp2 sp3', 'relative offsets', 'header gradient gx', 'gy',
             'header GCOEFF1-4']
    t1 = write_table(rows, colnames, outdir, 'flat.txt')
    t2 = write_table(hrows, hcols, outdir, 'flat_header_corrections.txt')
    summary = 'twilight flats: ' + '; '.join('%s per-fiber std %s -> %s after factors+gradient' % (r[0], r[5], r[11])
                                            for r in rows)
    return dict(rows=[colnames] + rows, hrows=[hcols] + hrows, figs=[f1, f2], tables=[t1, t2], summary=summary)


# ---------------------------------------------------------------------------------------------------------
# Explanations written into the html report.  The report declares ISO-8859-1 (xhtml.begin), so everything
# here must be plain ASCII: special characters are written as html entities.
# ---------------------------------------------------------------------------------------------------------

INTRO = '''<h2>What this report is</h2>
<p>This page compares two sets of LVM long-term master calibrations, A and B, made by
<code>CompareCalib.py</code> (in <code>lvm_ksl/py_progs</code>). B is the reference: every difference is
A&minus;B and every ratio is A/B. The two sets can be</p>
<ul>
<li><b>the same epoch processed twice</b> (e.g. our rebuild against the team's masters). The inputs are
(nearly) the same raw frames, so any difference comes from the code, its settings or the choice of input
frames, and ideally everything is close to zero;</li>
<li><b>two different epochs</b>. Then the differences also contain real changes of the instrument between the
two calibration sequences (the spectra and traces move by fractions of a pixel between epochs). What matters
is not the raw difference but what remains after the corrections the DRP makes to every science exposure, and
single fibers that behave differently from their neighbours.</li>
</ul>

<h2>Where the numbers come from</h2>
<p>Almost everything is measured from the <b>image data</b> of the masters, compared pixel by pixel (bias,
traces, wavelengths, LSF) or on a common wavelength grid (twilight flats). The fiber labels and status come
from the <code>SLITMAP</code> table of each master, and pixels marked in the <code>BADPIX</code> masks are left
out. <b>Header keywords are used only</b> in these places, each marked as such below:</p>
<ul>
<li>the <b>Provenance</b> table (MJD, DRPVER, COMMIT, NFRAMES, EXPOSURE and the <code>FRAMES</code>
extension);</li>
<li>the number of input frames (NFRAMES) in the bias table;</li>
<li>the amplifier <b>gains</b> in the trace table;</li>
<li>the fit <b>RMS</b> values in the wavelength table;</li>
<li>the <b>corrections recorded in the flat headers</b> (SKYCORR, spectrograph factors, gradient).</li>
</ul>
<p>Header values are what the DRP says about itself and are not always right: the DRP 1.3.0dev wave masters,
for example, record the LSF RMS where the wavelength RMS belongs.</p>

<h2>Ideas used throughout</h2>
<ul>
<li><b>Local trend.</b> For traces and wavelengths the difference A&minus;B is an array of 648 fibers &times;
4086 columns per camera. At every 16th column the running median over 15 neighbouring fibers (along the slit)
is taken and interpolated along the columns: that is the local trend. Smooth changes along the slit and along
the spectrum, such as a thermal shift between epochs, are part of the trend; a single fiber that differs from
its neighbours is not, because one fiber cannot move the median of 15. So <i>difference minus trend</i>
isolates single bad fibers.</li>
<li><b>Outlier fibers</b> are fibers whose difference minus trend exceeds a limit somewhere along their length:
%(trace)s pixels for traces; %(wave)s &Aring; for wavelengths inside the detector, and %(edge_lim)s &Aring;
within %(edge)d columns of either end, where the polynomials extrapolate past the last arc lines and two
independent fits always differ most.</li>
<li><b>Wiggle.</b> For a fiber, its value minus the mean of its two neighbouring fibers, along the fiber
(one neighbour for the first and last fiber); the median of that is removed, so a constant offset does not
count, and the scatter that remains is the wiggle. It is computed separately for A and for B, so it says which
of the two has a fiber that bends differently from its neighbours.</li>
<li><b>What the science reduction does not remove.</b> Every science exposure is corrected for some changes
since the masters were made: the traces by one shift per block of 36 fibers, the wavelengths by one offset per
fiber, smooth along the slit, measured from sky lines. Differences of that form between A and B do not matter
for science data. Each of those corrections is applied here to A&minus;B, and what is left is what would reach
the science spectra.</li>
</ul>
''' % dict(trace='%.2f' % TRACE_LIMIT, wave='%.2f' % WAVE_LIMIT, edge_lim='%.2f' % WAVE_EDGE_LIMIT, edge=EDGE)

PROVENANCE_EXPLAIN = '''<p><b>From the headers.</b> For the first file of each kind of master (camera b1, or
channel b for the twilight flat; the other cameras are not checked) the table gives the header keywords MJD
(the epoch), DRPVER (the version <i>label</i> the masters were made under, i.e. the environment variable
<code>LVMDRP_VERSION</code>, not the installed package version), COMMIT (the DRP git commit), NFRAMES (the
number of input exposures per camera) and EXPOSURE (one exposure number, which is not used consistently by the
different products). The last column summarizes the <code>FRAMES</code> extension, the list of input exposures,
as count: first&ndash;last; masters from older DRP versions have no such extension (&quot;not recorded&quot;),
so their inputs cannot be checked. When both sets list their frames and the lists differ, the exposures found
in only one set are listed below the table.</p>
'''

SECTION_EXPLAIN = {}

SECTION_EXPLAIN['bias'] = '''<p><b>What is compared.</b> <code>lvm-mbias-&lt;camera&gt;.fits</code> is the
combined bias image of each camera (normally about 95 frames from the 11 nights starting at the epoch). The
two images are subtracted pixel by pixel. The image is split into its four quadrants (the two halves in each
direction), taken to be the four amplifiers, which differ in noise, so the noise of the difference and the
count of discrepant pixels are worked out per quadrant.</p>
<p><b>Columns of the table.</b></p>
<ul>
<li><b>NFRAMES A/B</b>: number of input frames, from the headers.</li>
<li><b>median A, median B</b>: median of each master image (e<sup>&minus;</sup>).</li>
<li><b>median diff</b>: median of A&minus;B over the image: a change of the overall bias level.</li>
<li><b>robust sigma diff</b>: 1.4826 &times; the median absolute deviation of A&minus;B over the whole image:
how different the two images are pixel to pixel.</li>
<li><b>robust sigma by quadrant</b>: the same for each quadrant (lower left, lower right, upper left, upper
right).</li>
<li><b>n pixels &gt;5 sigma</b>, <b>%% pixels</b>: pixels whose difference departs from its quadrant's median
by more than 5 of that quadrant's robust sigmas: structure (columns, glow, edges) that is in one master and
not the other.</li>
<li><b>max |diff|</b>: the largest absolute difference.</li>
<li><b>equal</b>: <i>exact</i> if the two images are identical, <i>rounding</i> if they differ by less than
10<sup>&minus;4</sup>&nbsp;e<sup>&minus;</sup> everywhere (single-precision rounding), <i>no</i> otherwise.</li>
</ul>
<p><b>How to read it.</b> Two processings of the same bias frames with the same code should be equal to
rounding. Between epochs the difference has real structure (in particular a glow-like pattern in r1, a bad
column and the edges in b2), so a few percent of pixels beyond 5 sigma is normal.</p>
'''

SECTION_EXPLAIN['trace'] = '''<p><b>What is compared.</b> For each camera the trace masters describe where
each of the 648 fibers falls on the detector: <code>lvm-mcentroids</code> gives the centre of each fiber's
profile (a detector row) at every column, <code>lvm-msigmas</code> its width (Gaussian sigma), and
<code>lvm-mmodel</code> a model image of the whole flat built from the fitted profiles. All are fitted by the
DRP to the dome flats of the epoch's calibration sequence. The centroid and sigma arrays (fibers &times;
columns, pixels) of A and B are subtracted element by element; <i>dcen</i> and <i>dsig</i> below are those
differences. The model images are compared as a ratio, A/B, over the pixels where B's model exceeds 1000
counts, i.e. on the fibers.</p>
<p><b>Columns of the first table.</b></p>
<ul>
<li><b>median dcen</b>: median of A&minus;B over all fibers and columns: the overall shift of the traces
(for two epochs mostly a real shift of the instrument).</li>
<li><b>rms dcen</b>: the scatter of A&minus;B, including any smooth pattern.</li>
<li><b>rms dcen-trend</b>, <b>p99 |dcen-trend|</b>, <b>max |dcen-trend|</b>: the scatter, 99th percentile and
largest value of the difference after the local trend is removed: fiber-to-fiber differences only.</li>
<li><b>n fibers &gt;%(trace)s</b> and <b>worst fibers</b>: the fibers whose |dcen&minus;trend| exceeds
%(trace)s pixels somewhere, and the worst five as fiber(IFU label):largest value. Fibers are numbered by their
row in the camera, 0&ndash;647.</li>
<li><b>rms dsig</b>: the scatter of the difference in profile width.</li>
<li><b>model median |A/B-1|</b>, <b>model p99 |A/B-1|</b>: how different the model images are, as a fraction.
The model is in counts, so it scales with the gain used to convert the flats.</li>
<li><b>gains A/B (amp1-4)</b>: <i>from the headers</i>, the gain of each of the four amplifiers in A and in
B. A change of gain between the two sets shows up directly in the model ratio, but does not move the
centroids.</li>
</ul>
'''

SECTION_EXPLAIN['wave'] = '''<p><b>What is compared.</b> <code>lvm-mwave-&lt;camera&gt;.fits</code> gives the
wavelength (&Aring;) of every pixel of every fiber (648 fibers &times; 4086 columns), evaluated from the
polynomial the DRP fits to the arc lines of each fiber; <code>lvm-mlsf-&lt;camera&gt;.fits</code> gives the line
spread function width (FWHM) at the same places. Both come from the arc exposures of the epoch's calibration
sequence. The arrays of A and B are subtracted element by element (<i>dwave</i>, &Aring;; also as a velocity,
<i>dv</i> = dwave / &lambda; &times; c), with the pixels masked in either set's <code>BADPIX</code> left out; the
LSF is compared as the ratio A/B&minus;1. &quot;Inside&quot; means columns %(edge)d&ndash;%(edge_hi)d; &quot;edges&quot;
means the first and last %(edge)d columns.</p>
<p><b>Columns of the first table.</b></p>
<ul>
<li><b>median dwave</b>, <b>median dv</b>: the overall shift of the wavelength scale, in &Aring; and km/s.</li>
<li><b>rms dwave</b>: the scatter of A&minus;B, including smooth patterns (for two epochs mostly a real shift and
stretch of the spectra).</li>
<li><b>median per-fiber rms about trend</b>: for each fiber the rms of dwave&minus;trend along the fiber; the
median over fibers. The typical fiber-to-fiber disagreement.</li>
<li><b>p99 / max |dwave-trend| inside</b>, <b>max |dwave-trend| edges</b>: the 99th percentile and the
largest departure from the local trend, inside the detector and at its ends.</li>
<li><b>n fibers &gt;%(wave)s A inside or &gt;%(edge_lim)s A at edges</b> and <b>worst flagged fibers</b>: the outlier fibers (more than %(wave)s &Aring; from
the trend inside, or %(edge_lim)s &Aring; at the edges), and the worst five as fiber(IFU label):largest
value.</li>
<li><b>LSF median A/B-1</b>, <b>LSF p99 |A/B-1|</b>: the overall change of the LSF width, and its spread, as a
fraction.</li>
<li><b>wave RMS A, wave RMS B, FWHM RMS A, FWHM RMS B</b>: <i>from the headers</i>, the RMS of the DRP's fits of
the wavelength and of the LSF width to the arc lines, in pixels, as median/maximum over fibers (keywords
<code>&lt;CAM&gt; WAVE DISP RMS</code> and <code>&lt;CAM&gt; WAVE FWHM RMS</code>). Older versions wrote only
<code>PIPE DISP RMS</code>, shown as &quot;PIPE&quot; under FWHM, because it actually held the LSF RMS: the
wavelength RMS was overwritten, so for those masters it is unknown.</li>
</ul>
'''

SECTION_EXPLAIN['flat'] = '''<p><b>What is compared.</b> <code>lvm-mfiberflat_twilight-&lt;channel&gt;.fits</code>
is one file per channel (b, r, z), all three spectrographs stacked (1944 fibers): the relative throughput of
every fiber at every pixel, about 1, made from the dusk twilight exposures. The two flats can have different
wavelength solutions, so they are not compared pixel by pixel: each fiber of each flat is resampled onto a
common 0.5&nbsp;&Aring; grid and medianed in 50&nbsp;&Aring; bins (lvmdrp's <code>bin_fiberflat</code>), and the
ratio A/B is taken bin by bin. Its median over the bins gives one ratio per fiber. Only fibers of the science
and the two sky IFUs (not the standard-star fibers) with good status (<code>fibstatus</code> 0) in both flats
are used.</p>
<p>The per-fiber ratios, divided by their median (the <i>level</i>), are then fitted with lvmdrp's
<code>decompose_ratio</code>, the DRP's own model for the fiber flat corrections
(<code>fit_sky_gradient</code>): one factor per spectrograph times a linear gradient across each IFU. The fit
splits the ratio into large-scale parts, which the DRP's <code>fiberflat-corrections</code> step adjusts
anyway, and a fiber-to-fiber residual, which is the real change of the flat field.</p>
<p><b>Columns of the first table.</b></p>
<ul>
<li><b>n values</b>: number of fiber &times; bin ratios used.</li>
<li><b>median ratio</b>, <b>robust sigma</b>, <b>frac off &gt;1%%</b>: median and robust scatter of all those
ratios, and the fraction that differ from 1 by more than 1%%.</li>
<li><b>per-fiber std</b>: the scatter over fibers of the per-fiber ratio divided by the level.</li>
<li><b>level</b>: the median per-fiber ratio, i.e. the overall normalization of A relative to B.</li>
<li><b>spectrograph offsets sp1 sp2 sp3</b>: the fitted factor of each spectrograph, relative to the mean (a
fraction; 0.01 = 1%%).</li>
<li><b>IFU gradient gx, gy, amplitude</b>: the fitted linear gradient across the science IFU, its two components
and its size at the edge of the IFU relative to the mean.</li>
<li><b>per-fiber std after model</b>, <b>robust sigma after model</b>, <b>frac off &gt;1%% after model</b>: the
same statistics after dividing by the fitted spectrograph factors and gradient: the fiber-to-fiber change.</li>
</ul>
<p><b>Second table, from the headers.</b> What each flat says was applied to it by
<code>fiberflat-corrections</code>: <b>SKYCORR</b> (whether the step was run), the <b>spectrograph factors</b>
(and the same relative to their mean), and the recorded <b>IFU gradient</b> (gx, gy and the coefficients
GCOEFF1&ndash;4). If only one set has SKYCORR = True, the spectrograph offsets in the ratio should match its
factors, with the sign reversed if it is B. In DRP 1.3.2 the gradient was recorded but never applied to the
flat (lvmdrp issue #268), so it does not appear in the ratio.</p>
'''
for _k in ('trace', 'wave', 'flat'):
    SECTION_EXPLAIN[_k] = SECTION_EXPLAIN[_k] % dict(trace='%.2f' % TRACE_LIMIT, wave='%.2f' % WAVE_LIMIT,
                                                     edge_lim='%.2f' % WAVE_EDGE_LIMIT, edge=EDGE,
                                                     edge_hi=4085 - EDGE)
SECTION_EXPLAIN['bias'] = SECTION_EXPLAIN['bias'].replace('%%', '%')

DRP_MEANING = {
    'trace': '''<div style="border-left: 4px solid #888; padding-left: 10px;">
<p><b>What this section answers.</b> Suppose the instrument on some night was in exactly the state that A
describes, but its science frames were extracted with B's traces (or the other way round). The DRP would shift
B's traces by one constant per 36-fiber block, measured from the science frame, and extract. The numbers and
the map below are how far the extraction positions would then still be off, in pixels: the error in trace
position that would reach the science spectra.</p>
<p><b>How to read it.</b></p>
<ul>
<li><b>It is a disagreement between two masters, not an error against the truth.</b> It assumes one of the two
describes the night correctly. In practice the instrument on a given night is in some third state, close to but
not identical with either.</li>
<li><b>It is a best case for the correction.</b> In a real exposure the block shifts are measured by
cross-correlating the traces with the science frame at 11 columns, with noise, and blocks whose shift looks
pathological are replaced by the median. Here each block shift is the median of the noiseless difference of
the two masters. Real residuals will be somewhat larger.</li>
<li><b>It covers only what the correction cannot do:</b> single fibers that moved, shifts that differ within a
block, and changes along the fibers (tilt or curvature). A shift of a whole block, however large, is
removed.</li>
</ul></div>
''',
    'wave': '''<div style="border-left: 4px solid #888; padding-left: 10px;">
<p><b>What this section answers.</b> Suppose the instrument on some night was in exactly the state that A
describes, but its science frames were reduced with B's wavelength solution (or the other way round). The DRP
would measure each fiber's offset from the sky lines in the science frame and shift that fiber's wavelengths by
a constant. The numbers and the map below are the wavelength error that would then still be left in the
science spectra, in &Aring; (and km/s at the detector ends).</p>
<p><b>How to read it.</b></p>
<ul>
<li><b>It is a disagreement between two masters, not an error against the truth.</b> It assumes one of the two
describes the night correctly. In practice the instrument on a given night is in some third state, close to but
not identical with either.</li>
<li><b>It is a best case for the correction.</b> In a real exposure the per-fiber offsets are measured from sky
lines with noise; fibers with S/N below 10 or failed fits are skipped, and in b a single line (5577&nbsp;&Aring;)
is used. Here the offsets are read exactly off the noiseless difference of the two masters at the sky-line
wavelengths, then smoothed as the DRP does. Real residuals will be larger, by the sky-line measurement
noise.</li>
<li><b>It covers only what the correction cannot do.</b> The DRP shifts each fiber by a constant, so what
remains is a stretch or bend along the spectrum, or a single fiber whose fit went wrong. In b, where the offset
is pinned at 5577&nbsp;&Aring;, anything that grows away from that line stays.</li>
<li><b>The ends of the detector are reported separately.</b> The rms leaves out the first and last %d columns,
where the polynomials extrapolate; the blue-end and red-end columns of the table give what is left there. A
stretch of the wavelength scale shows up as opposite signs at the two ends.</li>
</ul></div>
''' % EDGE}

DRP_EXPLAIN = {
    'trace': '''<p><b>The correction.</b> When a science frame is extracted, the DRP cross-correlates the
master traces with the frame at 11 columns between 500 and 3000 and shifts each block of 36 fibers by one
constant, the median of its measured shifts (<code>_fix_fiber_thermal_shifts</code>). Here the same is done to
A&minus;B: for each of the 18 blocks the median of the difference over columns 450&ndash;3050 is subtracted.
A shift of a whole block, however large, is removed; what is left (single fibers, shifts that differ within a
block, changes along the fibers) would change the extraction of science frames.</p>
<p><b>Columns.</b> <b>block shift min/max</b>: the range of the 18 block shifts removed (pixels).
<b>rms left</b>, <b>p99 |left|</b>, <b>max |left|</b>: the scatter, 99th percentile and largest value of what is
left, over all fibers and columns. <b>n fibers &gt;%(trace)s left</b>: fibers still off by more than %(trace)s
pixels somewhere.</p>
''' % dict(trace='%.2f' % TRACE_LIMIT),
    'wave': '''<p><b>The correction.</b> In every science exposure the DRP measures each fiber's offset from
the reference sky lines (b: 5577&nbsp;&Aring;; r: 6364, 7359, 7392&nbsp;&Aring;; z: 8399, 8988, 9553,
9720&nbsp;&Aring;), median filters the offsets over 8 fibers, fits a spline with 10 knots along the slit, and
subtracts the result from the wavelengths as a constant per fiber (<code>shift_wave_skylines</code>). Here the
same is done to A&minus;B: the difference at those lines gives the per-fiber offset, smoothed in the same way
(if fewer than half the fibers have an offset nothing is removed, as in the DRP). A constant, slit-smooth shift
is removed; what is left (a stretch or bend along the spectrum, a single bad fiber) stays in the science
spectra. In b, where the offset comes from one line, any change away from 5577&nbsp;&Aring; is left.</p>
<p><b>Columns.</b> <b>sky offset removed min/max</b>: the range of the per-fiber offset removed (&Aring;).
<b>rms left inside</b>, <b>p99 |left| inside</b>: the scatter and 99th percentile of what is left, over columns
%(edge)d&ndash;%(edge_hi)d. <b>median left, blue/red end</b>: for each fiber the median of what is left over the
first (blue) or last (red) %(edge)d columns, and the median of that over fibers, in &Aring; and km/s: a
stretch of the wavelength scale shows up as opposite signs at the two ends. <b>n fibers &gt;%(wave)s A left
inside</b>: fibers still off by more than %(wave)s &Aring; inside. The map shows what is left, on a scale of
&plusmn;0.1&nbsp;&Aring;.</p>
''' % dict(wave='%.2f' % WAVE_LIMIT, edge=EDGE, edge_hi=4085 - EDGE)}

OUTLIER_EXPLAIN = {
    'trace': '''<p>Every fiber whose |dcen&minus;trend| exceeds %(trace)s pixels somewhere along its length,
worst first. <b>fiber</b>: row in the camera (0&ndash;647); <b>ifulabel</b>, <b>fibstatus</b>: from the
<code>SLITMAP</code> of A (fibstatus 0 is a good fiber). <b>max |dcen-trend|</b>: the largest departure
(pixels). <b>columns &gt;limit</b>: how many columns exceed the limit, and the first&ndash;last of them; a
number near 4086 means the whole fiber is offset. <b>max |dsig|</b>: the largest width difference.
<b>wiggle A</b>, <b>wiggle B</b>: the wiggle of the fiber's centroid in each set. Trace differences are usually
near-constant offsets along the whole fiber, which the wiggle cannot see, so no suspect is named: deciding
which set is right needs a comparison with the flat itself.</p>
''' % dict(trace='%.2f' % TRACE_LIMIT),
    'wave': '''<p>Every fiber that departs from the local trend by more than %(wave)s &Aring; inside, or
%(edge_lim)s &Aring; at the edges, worst first. <b>fiber</b>: row in the camera (0&ndash;647);
<b>ifulabel</b>: from the <code>SLITMAP</code> of A's trace master. <b>max |dwave-trend| inside</b>, <b>max at
edges</b>: the largest departure in each region (&Aring;). <b>column of max</b>: where along the fiber the
largest departure is. <b>rms about trend</b>: the rms of the departure along the fiber. <b>wiggle A</b>,
<b>wiggle B</b>: the wiggle of the fiber's wavelength solution in each set. <b>suspect</b>: the set with the
larger wiggle, i.e. the one whose solution bends differently from its neighbouring fibers, and so probably
the one whose fit went wrong. A wiggle of a few hundredths of an &Aring;ngstr&ouml;m is normal; a few tenths
is a failed fit.</p>
''' % dict(wave='%.2f' % WAVE_LIMIT, edge_lim='%.2f' % WAVE_EDGE_LIMIT)}

# captions for the figures, by file name
FIG_CAPTION = {
    'bias_diff.png': '''Bias A&minus;B, one panel per camera (rows against columns), each on a scale of
&plusmn;5 times its own robust sigma, so structure shows whatever the noise level.''',
    'trace_centroid_raw.png': '''Trace centroid difference A&minus;B (pixels), i.e. the change of each fiber's
position across the dispersion, fibers against columns: D(f,x) = c<sub>A</sub>(f,x) &minus; c<sub>B</sub>(f,x),
where c is the centroid (detector row) of fiber f at column x. All cameras share one fixed scale of
&plusmn;1 pixel and one colour bar: white is no change, red means the fiber is at a higher row in A, blue
a lower row. For two epochs this is mostly a real shift of the whole camera, and a camera that has shifted as a
whole is a uniform colour.''',
    'trace_centroid_minus_median.png': '''The same difference minus the camera's median shift,
D(f,x) &minus; median(D), the median (over all fibers and columns of the camera) being given in each panel
title. All cameras share one fixed scale of &plusmn;0.1 pixels. With the overall shift removed, this shows
how the fibers vary around it: whole 36-fiber blocks that moved more or less than the rest appear as bands,
and noisy or mis-traced fibers as thin lines. The DRP removes block shifts at extraction (see &quot;What the
science reduction does not remove&quot;), but not differences within a block.''',
    'trace_centroid_diff.png': '''The trace centroid difference A&minus;B minus the local trend (pixels), fibers
against columns, on a fixed scale of &plusmn;0.02 pixels for all cameras, so one colour bar serves all panels.
<br><b>What &quot;minus the local trend&quot; means.</b> The trend is worked out separately for each camera
from the difference image itself (648 fibers &times; 4086 columns). At every 16th column, each fiber's trend
value is the <i>median</i> of the difference over the 15 fibers centred on it (itself and 7 on each side, with
the fibers reflected at the ends of the slit); these values are interpolated linearly along the columns and
subtracted. So each fiber is compared with its 15 nearest neighbours on the slit, and anything the
neighbours share disappears: a whole block shifting, a smooth change along the slit, a slow change along the
columns such as a tilt. What is left is how much a single fiber differs from its neighbours. One bad fiber
cannot move the median of 15, so it stays and shows up as a thin horizontal line. At a block boundary, the
last fiber of a block still has 8 of its 15 window fibers in its own block, so the median follows its own block
and clean steps between blocks leave no false residual. The limitation: a group of more than about 7 adjacent
fibers that all changed together looks like part of the trend.''',
    'trace_sigma_diff.png': '''Profile width difference A&minus;B (pixels), on a fixed scale of
&plusmn;0.02 pixels.''',
    'trace_fiber_max.png': '''For each fiber (row in the camera), the largest |dcen&minus;trend| along its
length, on a logarithmic scale. The red line is the %.2f pixel limit for outliers.''' % TRACE_LIMIT,
    'wave_drp_residual.png': '''What is left of the wavelength difference after the sky-line correction
(&Aring;), fibers against columns, on a scale of &plusmn;0.1&nbsp;&Aring;.''',
    'wave_raw.png': '''Wavelength difference A&minus;B minus the camera's median (&Aring;), fibers against
columns; the median is in the panel title. As for the first trace map, the colours are centred on zero but each
camera has its own range, shown by its own colour bar and in its title, so the same colour can mean different
amounts in different panels. For two epochs this is mostly a real shift and stretch of the spectra.''',
    'wave_diff.png': '''The wavelength difference A&minus;B minus the local trend (&Aring;), fibers against
columns, on a fixed scale of &plusmn;0.05&nbsp;&Aring; for all cameras. The strongest differences normally sit at
the two ends, where the polynomials extrapolate.
<br><b>What &quot;minus the local trend&quot; means.</b> The trend is worked out separately for each camera
from the difference image itself (648 fibers &times; 4086 columns). At every 16th column, each fiber's trend
value is the <i>median</i> of the difference over the 15 fibers centred on it (itself and 7 on each side, with
the fibers reflected at the ends of the slit); these values are interpolated linearly along the columns and
subtracted. So each fiber is compared with its 15 nearest neighbours on the slit, and anything the
neighbours share disappears: a whole block shifting, a smooth change along the slit, a slow change along the
columns such as a tilt. What is left is how much a single fiber differs from its neighbours. One bad fiber
cannot move the median of 15, so it stays and shows up as a thin horizontal line. At a block boundary, the
last fiber of a block still has 8 of its 15 window fibers in its own block, so the median follows its own block
and clean steps between blocks leave no false residual. The limitation: a group of more than about 7 adjacent
fibers that all changed together looks like part of the trend.''',
    'lsf_ratio.png': '''LSF width ratio A/B&minus;1, fibers against columns, on a scale of &plusmn;5%.''',
    'wave_fiber_max.png': '''For each fiber (row in the camera), the largest |dwave&minus;trend| anywhere along
its length, on a logarithmic scale. The red line is the %.2f&nbsp;&Aring; limit used inside the detector; the
limit at the edges is %.2f&nbsp;&Aring;, so points between the two may come from the edges only.''' % (
        WAVE_LIMIT, WAVE_EDGE_LIMIT),
    'flat_ifu.png': '''The science IFU, one hexagon per fiber at its position, one row per channel. Left: the
per-fiber ratio A/B divided by the level. Middle: the fitted model, spectrograph factors times gradient. Right:
the residual, ratio divided by model, i.e. the fiber-to-fiber change. All three panels of a row share one
colour scale (&plusmn;3 times the scatter of the ratio, at least 1%), centred on 1.''',
    'flat_slit.png': '''The per-fiber ratio A/B divided by the level against the fiber's row in the flat (the
order along the slit), coloured by spectrograph: a step between spectrographs is a spectrograph factor; steps
within one spectrograph can come from the amplifiers (a gain change).''',
}


def doit(spec_a, spec_b, outroot='', only=None, labels=None, same=False, note=''):
    '''
    Compare two sets of masters and write outroot.html with its plots
    and tables in outroot/

    Description:

        spec_a and spec_b are resolved with resolve(); B is the
        reference.  only is a list of the products to compare, by
        default all of bias, trace, wave and flat.  same forces the
        same-epoch wording; note is a paragraph shown at the top

    History::

        261101 ksl Coding begun
    '''
    a, laba = resolve(spec_a)
    b, labb = resolve(spec_b)
    if labels:
        laba, labb = labels
    for folder in (a, b):
        if not os.path.isdir(folder):
            print('CompareCalib: %s is not a directory' % folder)
            return
        if not has_masters(folder):
            print('CompareCalib: no master calibrations (lvm-m*.fits) in %s' % folder)
            return
    if outroot == '':
        outroot = 'CompareCalib_%s_vs_%s' % (laba.replace('/', '_'), labb.replace('/', '_'))
    outdir = outroot
    os.makedirs(outdir, exist_ok=True)
    only = only or PRODUCTS

    mjds = []
    for folder in (a, b):
        mjd = '--'
        for f in sorted(os.listdir(folder)):
            if f.startswith('lvm-m') and f.endswith('.fits'):
                mjd = fits.getheader(os.path.join(folder, f)).get('MJD', '--')
                break
        mjds.append(mjd)
    if same and mjds[0] != mjds[1]:
        mode = 'the same epoch (the MJD keywords say %s and %s, the MJDs of the first input frames)' % tuple(mjds)
    elif mjds[0] == mjds[1]:
        mode = 'the same epoch (%s)' % mjds[0]
    else:
        mode = 'two different epochs (%s and %s)' % tuple(mjds)

    string = xhtml.begin('Calibration comparison: %s vs %s' % (laba, labb))
    string += xhtml.add_list(['A = %s: %s' % (laba, a), 'B = %s: %s (reference)' % (labb, b),
                              'Comparing %s; ratios and differences are A/B and A-B' % mode])
    if note:
        string += '<p style="border-left:4px solid #888;padding-left:8px"><b>Note:</b> %s</p>\n' % note
    string += INTRO
    string += xhtml.hline()
    string += xhtml.h2('Provenance')
    string += PROVENANCE_EXPLAIN
    prov_rows, prov_notes = compare_provenance(a, b, laba, labb)
    string += xhtml.table(prov_rows)
    if prov_notes:
        string += xhtml.add_list(prov_notes)

    def add_figure(fig):
        shape = plt.imread(fig).shape
        caption = FIG_CAPTION.get(os.path.basename(fig))
        out = '<p><i>%s</i></p>\n' % caption if caption else ''
        return out + xhtml.image(fig, width=1200, height=1200 * shape[0] // shape[1])

    summaries = []
    sections = [('bias', 'Bias', compare_bias), ('trace', 'Traces', compare_traces),
                ('wave', 'Wavelengths and LSF', compare_wave), ('flat', 'Twilight fiber flats', compare_flats)]
    for key, title, func in sections:
        if key not in only:
            continue
        print('Comparing %s' % title.lower())
        string += xhtml.hline()
        string += xhtml.h2(title)
        try:
            result = func(a, b, laba, labb, outdir)
        except Exception as e:
            string += xhtml.paragraph('Comparison failed: %s: %s' % (type(e).__name__, e))
            summaries.append('%s: comparison failed: %s: %s' % (key, type(e).__name__, e))
            continue
        if result is None:
            string += xhtml.paragraph('Not present in both sets')
            summaries.append('%s: not present in both sets' % key)
            continue
        summaries.append(result['summary'])
        string += xhtml.paragraph('<b>Summary:</b> ' + result['summary'])
        string += SECTION_EXPLAIN.get(key, '')
        string += xhtml.table(result['rows'])
        if 'hrows' in result:
            string += xhtml.h3('Corrections recorded in the flat headers')
            string += xhtml.table(result['hrows'])
        if 'drows' in result:
            string += xhtml.h3('What the science reduction does not remove')
            string += DRP_MEANING.get(key, '')
            string += xhtml.paragraph('<b>Summary:</b> ' + result['drp_summary'])
            string += DRP_EXPLAIN.get(key, '')
            string += xhtml.table(result['drows'])
            for fig in result.get('dfigs', []):
                string += add_figure(fig)
            summaries.append('  ' + result['drp_summary'])
        if result['figs']:
            string += xhtml.h3('Maps and plots')
        for fig in result['figs']:
            string += add_figure(fig)
        if 'outliers' in result:
            string += xhtml.h3('Outlier fibers (%d), largest first' % (len(result['outliers']) - 1))
            string += OUTLIER_EXPLAIN.get(key, '')
            string += xhtml.table(result['outliers'][:41])
            if len(result['outliers']) > 41:
                string += xhtml.paragraph('First 40 shown; all are in the table file below')
        string += xhtml.paragraph('Tables: ' + ', '.join(xhtml.link(os.path.basename(t), t) for t in result['tables']))

    string += xhtml.hline()
    string += xhtml.end()
    with open(outroot + '.html', 'w') as g:
        g.write(string)

    print('\nCompareCalib: %s vs %s, %s' % (laba, labb, mode))
    for note in prov_notes:
        print('  provenance: %s' % note)
    for s in summaries:
        print('  ' + s)
    print('Report: %s.html (plots and tables in %s/)' % (outroot, outdir))


def steer(argv):
    '''
    Just a steering routine
    '''
    outroot = ''
    only = None
    labels = None
    same = False
    note = ''
    specs = []
    i = 1
    while i < len(argv):
        if argv[i][0:2] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-out':
            i += 1
            outroot = argv[i]
        elif argv[i] == '-only':
            i += 1
            only = argv[i].split(',')
            unknown = [x for x in only if x not in PRODUCTS]
            if unknown:
                print('CompareCalib: unknown products %s; choose from %s' % (unknown, PRODUCTS))
                return
        elif argv[i] == '-labels':
            i += 1
            labels = argv[i].split(',')
            if len(labels) != 2:
                print('CompareCalib: -labels needs two names separated by a comma')
                return
        elif argv[i] == '-same':
            same = True
        elif argv[i] == '-note':
            i += 1
            note = argv[i]
        elif argv[i][0] == '-':
            print('Error: unknown option %s' % argv[i])
            print(_USAGE)
            return
        else:
            specs.append(argv[i])
        i += 1
    if len(specs) != 2:
        print(_USAGE)
        return
    doit(specs[0], specs[1], outroot=outroot, only=only, labels=labels, same=same, note=note)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
