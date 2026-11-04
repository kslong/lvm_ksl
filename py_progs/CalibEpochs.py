#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Take a census of the long-term master calibrations of all the
    calibration epochs, reduce each epoch and camera to a few numbers,
    and find the epochs that are outliers or where the calibrations
    changed in a step, with an html report of interactive plots and a
    ranked table.

Command line usage (if any):

    usage: CalibEpochs.py [-h] [-out outroot] [-calib DIR] [-epochs FILE] [-only LIST] [-thresh T]
           CalibEpochs.py [-h] [-out outroot] [-thresh T] -table calib_epochs.txt

    Options::

        -h           print this help and exit
        -out outroot name of the html file (outroot.html) and of the
                     subdirectory with its plots and tables (outroot/);
                     default CalibEpochs
        -calib DIR   directory holding one subdirectory of masters per
                     epoch (default $LVM_MASTER_DIR)
        -epochs FILE lvmcore calibration-epochs.yaml (default
                     $LVMCORE_DIR/calibrations/calibration-epochs.yaml)
        -only LIST   measure only some products, a comma separated list
                     of bias, trace, wave and flat
        -thresh T    flag a change larger than T times the typical
                     change between consecutive epochs (default 4)
        -table FILE  make the report from a calib_epochs.txt table
                     written by an earlier run, without reading masters

Description::

    For each epoch and camera the masters are reduced to these numbers:

    Traces (lvm-mcentroids, -msigmas)
        trace_pos     median fiber position (row) at the central columns
        trace_tilt    median change of fiber position from column 500 to
                      column 3500
        trace_width   median fiber sigma at the central columns
        trace_nmasked number of fibers masked along their whole length
        trace_left    rms of the centroid difference from the previous
                      epoch after one shift per block of 36 fibers (the
                      correction the DRP applies at extraction), i.e. the
                      part of the change the science reduction keeps
    Wavelengths (lvm-mwave, -mlsf)
        wave_c        median wavelength at the central column
        wave_disp     median dispersion (A/pixel) between columns 200
                      and 3886
        lsf_c         median LSF FWHM at the central column
        wave_nwiggle  number of fibers whose solution wiggles relative to
                      its neighbouring fibers (more than 5 times the
                      camera's median wiggle, and more than 0.03 A)
        wave_maxwiggle the largest such wiggle, flagged (with the fiber)
                      if above 0.1 A whatever the other epochs
        wave_left     rms of the wavelength difference from the previous
                      epoch after the slit-smooth per-fiber sky line
                      offset the DRP applies (shift_wave_skylines)
    Bias (lvm-mbias)
        bias_level    median level
        bias_qrange   range of the four quadrant medians
        bias_sigma    median over the quadrants of the robust sigma
        bias_nhot     pixels more than 10 sigma above their quadrant
    Twilight fiber flats (lvm-mfiberflat_twilight), per channel
        flat_sp1..3   spectrograph offsets of the ratio to the previous
                      epoch's flat (lvmdrp.qa.fiberflats.decompose_ratio)
        flat_grad     gradient across the science IFU of that ratio
        flat_left     per-fiber scatter of that ratio after the
                      spectrograph factors and gradient
        flat_nodd     fibers of good status whose flat is below 0.2 or
                      above 2, flagged (with the fibers) whatever the
                      other epochs

    Positions shift between epochs for all the spectrographs together,
    so trace_pos, trace_tilt, wave_c, wave_disp and lsf_c are also given
    relative to the mean of the other two spectrographs of the same
    channel (the _rel quantities), which brings out a problem in one
    spectrograph.

    Each series is then searched for outliers and steps, with the
    typical change between consecutive epochs (robust sigma of the
    differences) as the scale:

        outlier  an epoch that differs from both neighbours by more than
                 thresh times the scale, in opposite directions, while
                 the two neighbours agree with each other
        step     a change from the previous epoch larger than thresh
                 times the scale that is not part of an outlier

    For the quantities that compare an epoch with the previous one
    (trace_left, wave_left, flat_*), values more than thresh robust
    sigmas above the median are flagged as large changes.  The changes
    between epochs are heavy tailed, so every flag must also exceed a
    minimum size for its quantity (MINIMUM), near the 97th percentile
    of the changes between consecutive epochs; bias_nhot is not
    flagged.  A step in a position is usually a real change of the
    instrument; if it happened during an epoch, the CalibDrift.py
    report shows large corrections in that epoch.

    The report shows each quantity against epoch, one panel per camera,
    with the flagged epochs marked, and a table of all flags ranked by
    their size.  CompareCalib.py on a flagged epoch and its neighbours
    gives the details.

Primary routines:

    doit
    steer

Notes::

    Requires lvmdrp (fiber flat routines, REF_SKYLINES) through
    CompareCalib.py.  Reading the masters of 43 epochs takes several
    minutes, mostly the bias frames; the table is written to
    outroot/calib_epochs.txt, so the report can be remade with -table.

History::

    261101 ksl Coding begun, to find outlier epochs and steps in the
               long-term master calibrations of all the epochs
    261104 ksl Report explains where the numbers come from (image
               data, no headers) and how each quantity is measured;
               axis labels on the plots; ASCII-only html

'''

import os
import re
import sys
import time
import numpy as np
import bottleneck as bn
from astropy.io import fits
from astropy.table import Table
import xhtml
from CompareCalib import drp_trace_residual, drp_wave_residual, robust_sigma
from FixFlatGradient import default_epochs_file, load_epochs


def _usage_from_doc(doc):
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


_USAGE = _usage_from_doc(__doc__)

CAMERAS = ['b1', 'b2', 'b3', 'r1', 'r2', 'r3', 'z1', 'z2', 'z3']
CHANNELS = ['b', 'r', 'z']
PRODUCTS = ['bias', 'trace', 'wave', 'flat']
CENTER = slice(1993, 2094)
EDGE = 200

# name: (unit, group, description); the _rel versions are added for the RELATIVE ones
QUANTITIES = {
    'trace_pos': ('pix', 'trace', 'median fiber position at the central columns'),
    'trace_tilt': ('pix', 'trace', 'median change of fiber position from column 500 to 3500'),
    'trace_width': ('pix', 'trace', 'median fiber sigma at the central columns'),
    'trace_nmasked': ('fibers', 'trace', 'fibers masked along their whole length'),
    'trace_left': ('pix', 'trace', 'rms change from the previous epoch left after one shift per fiber block'),
    'wave_c': ('A', 'wave', 'median wavelength at the central column'),
    'wave_disp': ('A/pix', 'wave', 'median dispersion between columns 200 and 3886'),
    'lsf_c': ('A', 'wave', 'median LSF FWHM at the central column'),
    'wave_nwiggle': ('fibers', 'wave', 'fibers whose solution wiggles relative to their neighbours'),
    'wave_maxwiggle': ('A', 'wave', 'largest wiggle of a fiber solution relative to its neighbours'),
    'wave_left': ('A', 'wave', 'rms change from the previous epoch left after the sky line offset'),
    'bias_level': ('e-', 'bias', 'median bias level'),
    'bias_qrange': ('e-', 'bias', 'range of the four quadrant medians'),
    'bias_sigma': ('e-', 'bias', 'median over the quadrants of the robust sigma'),
    'bias_nhot': ('pixels', 'bias', 'pixels more than 10 sigma above their quadrant'),
}
RELATIVE = ['trace_pos', 'trace_tilt', 'wave_c', 'wave_disp', 'lsf_c']
FLAT_QUANTITIES = {
    'flat_sp1': ('', 'flat', 'spectrograph 1 offset of the ratio to the previous epoch'),
    'flat_sp2': ('', 'flat', 'spectrograph 2 offset of the ratio to the previous epoch'),
    'flat_sp3': ('', 'flat', 'spectrograph 3 offset of the ratio to the previous epoch'),
    'flat_grad': ('', 'flat', 'gradient across the science IFU of the ratio to the previous epoch'),
    'flat_left': ('', 'flat', 'per-fiber scatter of the ratio after spectrograph factors and gradient'),
    'flat_nodd': ('fibers', 'flat', 'fibers of good status whose flat is below 0.2 or above 2'),
}
# quantities that compare an epoch with the previous one: flagged by size, not by the outlier/step test
CHANGE_QUANTITIES = ['trace_left', 'wave_left', 'flat_sp1', 'flat_sp2', 'flat_sp3', 'flat_grad', 'flat_left']
# smallest scale for counts, so a series that is nearly constant does not flag a change of one
FLOOR = {'trace_nmasked': 1.0, 'wave_nwiggle': 1.0, 'bias_nhot': 10.0}
# smallest change that is flagged, whatever the statistics: the changes between epochs are heavy tailed,
# so a purely statistical threshold flags hundreds of small ones; these are near the 97th percentile of the
# changes between consecutive epochs (DRP 1.3.0dev masters, 43 epochs) and of a size that matters for the data
MINIMUM = {'trace_pos_rel': 0.5, 'trace_tilt_rel': 0.05, 'trace_width': 0.03, 'trace_nmasked': 2,
           'trace_left': 0.03, 'wave_c_rel': 0.25, 'wave_disp_rel': 2e-5, 'lsf_c_rel': 0.03, 'wave_nwiggle': 2,
           'wave_maxwiggle': 0.05, 'wave_left': 0.02, 'bias_level': 0.15, 'bias_qrange': 0.3, 'bias_sigma': 0.25,
           'flat_sp1': 0.01, 'flat_sp2': 0.01, 'flat_sp3': 0.01, 'flat_grad': 0.01, 'flat_left': 0.005}
# measured and plotted, but too unstable between epochs to be flagged
NOT_FLAGGED = ['bias_nhot', 'wave_maxwiggle', 'flat_nodd']


def mpath(folder, kind, name):
    return os.path.join(folder, 'lvm-%s-%s.fits' % (kind, name))


def wiggles(wave, step=8):
    '''
    Per fiber, the scatter along the fiber of its wavelength offset from
    the mean of its two neighbours, after removing the median offset
    '''
    w = wave[:, ::step].astype(float)
    nb = np.empty_like(w)
    nb[1:-1] = 0.5 * (w[:-2] + w[2:])
    nb[0], nb[-1] = w[1], w[-2]
    d = w - nb
    d = d - bn.nanmedian(d, axis=1)[:, None]
    return bn.nanstd(d, axis=1)


def trace_numbers(folder, cam, prev_cent):
    cent = fits.getdata(mpath(folder, 'mcentroids', cam)).astype(float)
    bad = fits.getdata(mpath(folder, 'mcentroids', cam), 'BADPIX').astype(bool)
    sig = fits.getdata(mpath(folder, 'msigmas', cam)).astype(float)
    c = np.where(bad, np.nan, cent)
    out = dict(trace_pos=bn.nanmedian(c[:, CENTER]),
               trace_tilt=bn.nanmedian(bn.nanmedian(c[:, 3450:3551], axis=1) - bn.nanmedian(c[:, 450:551], axis=1)),
               trace_width=bn.nanmedian(np.where(bad, np.nan, sig)[:, CENTER]),
               trace_nmasked=int(np.sum(bad.all(axis=1))), trace_left=np.nan)
    if prev_cent is not None:
        res, _ = drp_trace_residual(c - prev_cent)
        out['trace_left'] = bn.nanstd(res)
    return out, c


def wave_numbers(folder, cam, prev_wave):
    wave = fits.getdata(mpath(folder, 'mwave', cam)).astype(float)
    lsf = fits.getdata(mpath(folder, 'mlsf', cam)).astype(float)
    with fits.open(mpath(folder, 'mwave', cam)) as h:
        bad = h['BADPIX'].data.astype(bool) if 'BADPIX' in h else np.zeros(wave.shape, bool)
    w = np.where(bad, np.nan, wave)
    wig = wiggles(w)
    med = bn.nanmedian(wig)
    out = dict(wave_c=bn.nanmedian(w[:, CENTER]),
               wave_disp=bn.nanmedian((w[:, -EDGE] - w[:, EDGE]) / (wave.shape[1] - 2 * EDGE)),
               lsf_c=bn.nanmedian(np.where(bad, np.nan, lsf)[:, CENTER]),
               wave_nwiggle=int(np.sum((wig > 5 * med) & (wig > 0.03))), wave_maxwiggle=bn.nanmax(wig),
               wave_wigglefiber=int(np.nanargmax(wig)) if np.isfinite(wig).any() else -1, wave_left=np.nan)
    if prev_wave is not None:
        res, _ = drp_wave_residual(w - prev_wave, prev_wave, cam[0])
        out['wave_left'] = bn.nanstd(res[:, EDGE:-EDGE])
    return out, w


def bias_numbers(folder, cam):
    d = fits.getdata(mpath(folder, 'mbias', cam)).astype(np.float32)
    ny, nx = d.shape
    levels, sigmas, nhot = [], [], 0
    for rs in (slice(0, ny // 2), slice(ny // 2, ny)):
        for cs in (slice(0, nx // 2), slice(nx // 2, nx)):
            q = d[rs, cs]
            q = q[np.isfinite(q)]
            m = np.median(q)
            s = robust_sigma(q[::7])
            levels.append(m)
            sigmas.append(s)
            nhot += int(np.sum(q > m + 10 * s))
    return dict(bias_level=float(np.nanmedian(d)), bias_qrange=float(np.ptp(levels)),
                bias_sigma=float(np.median(sigmas)), bias_nhot=nhot)


def flat_numbers(folder, ch, prev):
    from lvmdrp.core.rss import RSS
    from lvmdrp.qa.fiberflats import wavelength_edges, bin_fiberflat, decompose_ratio
    flat = RSS.from_file(mpath(folder, 'mfiberflat_twilight', ch))
    out = dict(flat_sp1=np.nan, flat_sp2=np.nan, flat_sp3=np.nan, flat_grad=np.nan, flat_left=np.nan)
    # fibers of good status whose flat is far from 1: dividing by it would wreck their spectra
    data = np.where(flat._mask.astype(bool), np.nan, flat._data) if flat._mask is not None else flat._data
    level = bn.nanmedian(data, axis=1)
    status = np.asarray(flat._slitmap['fibstatus']) if 'fibstatus' in flat._slitmap.colnames else np.zeros(len(level))
    odd = np.where((status == 0) & ((level < 0.2) | (level > 2)))[0]
    out['flat_nodd'] = len(odd)
    out['flat_oddfibers'] = ' '.join('%d(%s)' % (int(flat._slitmap['fiberid'][i]), str(flat._slitmap['orig_ifulabel'][i]).strip())
                                     for i in odd) or '--'
    if prev is None:
        return out, flat
    edges = wavelength_edges(prev, bin_width=50.0)
    ratio = bin_fiberflat(flat, edges) / bin_fiberflat(prev, edges)
    sm = prev._slitmap
    tel = np.char.strip(np.asarray(sm['telescope']).astype(str))
    good = np.isin(tel, ('Sci', 'SkyE', 'SkyW'))
    for f in (flat, prev):
        if 'fibstatus' in f._slitmap.colnames:
            good &= np.asarray(f._slitmap['fibstatus']) == 0
    fiber_ratio = bn.nanmedian(ratio, axis=1)
    fiber_ratio[~good] = np.nan
    parts = decompose_ratio(sm, fiber_ratio, good, deg=1)
    sci = parts['gradients'].get('Sci', {})
    out.update({'flat_sp%d' % i: parts['offsets'].get(i, np.nan) for i in (1, 2, 3)})
    out['flat_grad'] = sci.get('amp', np.nan)
    out['flat_left'] = np.nanstd(fiber_ratio[good] / (parts['level'] * parts['model'][good]))
    return out, flat


def measure(epochs, calib, only):
    '''
    The census: one row per epoch and camera (and per epoch and channel
    for the flats), measured in epoch order so each epoch can be
    compared with the previous one
    '''
    rows = []
    prev_c, prev_w, prev_f = {}, {}, {}
    t0 = time.time()
    for k, e in enumerate(epochs):
        folder = os.path.join(calib, str(e))
        print('  epoch %d (%d of %d, %.0f s)' % (e, k + 1, len(epochs), time.time() - t0))
        for cam in CAMERAS:
            row = dict(epoch=int(e), camera=cam)
            try:
                if 'trace' in only:
                    out, prev_c[cam] = trace_numbers(folder, cam, prev_c.get(cam))
                    row.update(out)
                if 'wave' in only:
                    out, prev_w[cam] = wave_numbers(folder, cam, prev_w.get(cam))
                    row.update(out)
                if 'bias' in only:
                    row.update(bias_numbers(folder, cam))
            except Exception as err:
                print('CalibEpochs: epoch %d %s: %s: %s' % (e, cam, type(err).__name__, err))
            rows.append(row)
        if 'flat' in only:
            for ch in CHANNELS:
                row = dict(epoch=int(e), camera=ch)
                try:
                    out, prev_f[ch] = flat_numbers(folder, ch, prev_f.get(ch))
                    row.update(out)
                except Exception as err:
                    print('CalibEpochs: epoch %d %s flat: %s: %s' % (e, ch, type(err).__name__, err))
                rows.append(row)
    names = ['epoch', 'camera'] + list(QUANTITIES) + ['wave_wigglefiber'] + list(FLAT_QUANTITIES)
    xtab = Table(names=names + ['flat_oddfibers'], dtype=[int, 'U2'] + [float] * (len(names) - 2) + ['U200'])
    for r in rows:
        xtab.add_row([r.get(n, np.nan) for n in names] + [r.get('flat_oddfibers', '--')])
    return xtab


def add_relative(xtab):
    '''
    For the RELATIVE quantities, the value minus the mean of the other two
    spectrographs of the same channel at the same epoch
    '''
    for q in RELATIVE:
        rel = np.full(len(xtab), np.nan)
        for i, r in enumerate(xtab):
            if r['camera'] not in CAMERAS:
                continue
            others = [o for o in CAMERAS if o[0] == r['camera'][0] and o != r['camera']]
            s = (xtab['epoch'] == r['epoch']) & np.isin(xtab['camera'], others)
            rel[i] = r[q] - np.nanmean(xtab[q][s])
        xtab['%s_rel' % q] = rel
    return xtab


def series(xtab, q, cam):
    s = xtab['camera'] == cam
    return np.asarray(xtab['epoch'][s]), np.asarray(xtab[q][s], float)


def classify(epochs, x, q, thresh):
    '''
    Outliers and steps in one series; returns a list of
    (epoch, kind, value, previous, next, score)
    '''
    flags = []
    if q in NOT_FLAGGED:
        return flags
    ok = np.isfinite(x)
    e, v = epochs[ok], x[ok]
    if len(v) < 4:
        return flags
    minimum = MINIMUM.get(q, 0.0)
    if q in CHANGE_QUANTITIES:
        scale = max(robust_sigma(v), FLOOR.get(q, 0.0), 1e-12)
        med = np.median(v)
        for i in range(len(v)):
            dev = abs(v[i] - med) if q.startswith('flat_sp') or q == 'flat_grad' else v[i] - med
            score = dev / scale
            if score > thresh and dev > minimum:
                flags.append((int(e[i]), 'large change', v[i], v[i - 1] if i else np.nan,
                              v[i + 1] if i + 1 < len(v) else np.nan, score))
        return flags
    diffs = np.diff(v)
    scale = max(robust_sigma(diffs), FLOOR.get(q, 0.0), 1e-12)
    outlier = np.zeros(len(v), bool)
    for i in range(1, len(v) - 1):
        dp, dn = v[i] - v[i - 1], v[i + 1] - v[i]
        small = min(abs(dp), abs(dn))
        if dp * dn < 0 and small > max(thresh * scale, minimum) and abs(v[i + 1] - v[i - 1]) < 0.5 * small:
            outlier[i] = True
            flags.append((int(e[i]), 'outlier', v[i], v[i - 1], v[i + 1], small / scale))
    for i in range(1, len(v)):
        dp = v[i] - v[i - 1]
        if abs(dp) > max(thresh * scale, minimum) and not outlier[i] and not outlier[i - 1]:
            flags.append((int(e[i]), 'step', v[i], v[i - 1], v[i + 1] if i + 1 < len(v) else np.nan, abs(dp) / scale))
    return flags


# absolute limits: an epoch is flagged whatever its neighbours if a value exceeds these
ABSOLUTE = {'wave_maxwiggle': 0.1, 'flat_nodd': 0.5}


def absolute_flags(xtab):
    '''
    Epochs whose masters are bad in themselves: a fiber whose wavelength
    solution wiggles by more than ABSOLUTE['wave_maxwiggle'] relative to
    its neighbours, or good-status fibers whose twilight flat is below
    0.2 or above 2
    '''
    rows = []
    for r in xtab:
        q = 'wave_maxwiggle'
        if r['camera'] in CAMERAS and np.isfinite(r[q]) and r[q] > ABSOLUTE[q]:
            rows.append((int(r['epoch']), r['camera'], q, 'bad fiber %d' % int(r['wave_wigglefiber']), r[q], np.nan,
                         np.nan, r[q] / ABSOLUTE[q]))
        q = 'flat_nodd'
        if r['camera'] in CHANNELS and np.isfinite(r[q]) and r[q] > ABSOLUTE[q]:
            rows.append((int(r['epoch']), r['camera'], q, 'bad flat %s' % r['flat_oddfibers'], r[q], np.nan,
                         np.nan, 10.0 * r[q]))
    return rows


def find_flags(xtab, thresh):
    rows = absolute_flags(xtab)
    camera_q = [q for q in QUANTITIES if q not in RELATIVE] + ['%s_rel' % q for q in RELATIVE]
    for q in camera_q:
        for cam in CAMERAS:
            e, x = series(xtab, q, cam)
            rows += [(f[0], cam, q) + f[1:] for f in classify(e, x, q, thresh)]
    for q in FLAT_QUANTITIES:
        for ch in CHANNELS:
            e, x = series(xtab, q, ch)
            rows += [(f[0], ch, q) + f[1:] for f in classify(e, x, q, thresh)]
    rows.sort(key=lambda r: -r[-1])
    return rows


def _unit(q):
    base = q[:-4] if q.endswith('_rel') else q
    return (QUANTITIES.get(base) or FLAT_QUANTITIES.get(base))[0]


def _describe(q):
    base = q[:-4] if q.endswith('_rel') else q
    d = (QUANTITIES.get(base) or FLAT_QUANTITIES.get(base))[2]
    return d + (', relative to the other two spectrographs' if q.endswith('_rel') else '')


def plot_quantity(xtab, q, flags, cams):
    '''
    One panel per camera (or channel) of a quantity against epoch, minus
    the camera's median over the epochs, with the flagged epochs marked
    '''
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    ncol = 3 if len(cams) == 9 else 1
    nrow = 3
    fig = make_subplots(rows=nrow, cols=ncol, shared_xaxes='all', shared_yaxes='all', vertical_spacing=0.06,
                        horizontal_spacing=0.03, subplot_titles=cams)
    colors = {'outlier': 'red', 'step': 'darkorange', 'large change': 'purple'}
    for i, cam in enumerate(cams):
        row, col = i // ncol + 1, i % ncol + 1
        e, x = series(xtab, q, cam)
        med = np.nanmedian(x) if np.isfinite(x).any() else 0.0
        fig.add_trace(go.Scatter(x=e, y=x - med, mode='lines+markers', line=dict(color='#1f77b4', width=1),
                                 marker=dict(size=5), customdata=x, name=cam, showlegend=False,
                                 hovertemplate='epoch %{x}<br>value %{customdata:.5g}<br>minus median %{y:.4g}'
                                               '<extra></extra>'), row=row, col=col)
        for kind, color in colors.items():
            fe = [f[0] for f in flags if f[1] == cam and f[2] == q and f[3] == kind]
            if fe:
                k = np.isin(e, fe)
                fig.add_trace(go.Scatter(x=e[k], y=x[k] - med, mode='markers',
                                         marker=dict(size=11, color=color, symbol='circle-open', line=dict(width=2)),
                                         name=kind, legendgroup=kind, showlegend=False,
                                         hovertemplate=kind + ' at epoch %{x}<extra></extra>'), row=row, col=col)
        fig.add_hline(y=0, line=dict(color='#888888', width=1), row=row, col=col)
    for kind, color in colors.items():
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', name=kind,
                                 marker=dict(size=11, color=color, symbol='circle-open', line=dict(width=2))))
    fig.update_xaxes(tickformat='d')
    # axis labels: y on the left-hand panels (the y axes are shared), x on the bottom row
    unit = _unit(q) or 'fraction'
    for r in range(1, nrow + 1):
        fig.update_yaxes(title_text='%s - median (%s)' % (q, unit), title_font=dict(size=11), row=r, col=1)
    for c in range(1, ncol + 1):
        fig.update_xaxes(title_text='epoch (MJD)', row=nrow, col=c)
    fig.update_layout(title='%s (%s): %s' % (q, _unit(q), _describe(q)), height=750 if ncol == 3 else 650,
                      margin=dict(l=60, r=20, t=80, b=40), legend=dict(orientation='h', y=1.08))
    return fig.to_html(full_html=False, include_plotlyjs=False)


def _plotly_js(outdir):
    jsfile = os.path.join(outdir, 'plotly.min.js')
    if not os.path.isfile(jsfile):
        from plotly.offline import get_plotlyjs
        with open(jsfile, 'w') as g:
            g.write(get_plotlyjs())
    return '<script src="%s"></script>\n' % jsfile


PROVENANCE = """<h2>Where these numbers come from</h2>
<p>Every number on this page was measured by <code>CalibEpochs.py</code> (in <code>lvm_ksl/py_progs</code>)
from the <b>image data</b> of the team's long-term master calibrations: one directory of masters per
calibration epoch in <code>$LVM_MASTER_DIR</code>, built with DRP 1.3.0dev and used for the DRP 1.3.2 science
reductions. <b>No header keywords are used.</b> Besides the data arrays, the program reads only the
<code>BADPIX</code> mask extensions, the <code>WAVE_TRACE</code> and <code>SLITMAP</code> tables of the twilight
flats (their wavelengths, and each fiber's status and label), and the list of epochs from lvmcore's
<code>calibration-epochs.yaml</code>. So the values do not depend on what the DRP wrote into the headers, some
of which are known to be wrong (e.g. the 1.3.0dev wave masters record the LSF RMS where the wavelength RMS
belongs).</p>
<p>The masters themselves were made by the team, not by us. How each quantity is measured, and what it is
sensitive to, is explained at the head of its section of plots below.</p>
"""

# detailed explanation at the head of each section of plots
SECTION_EXPLAIN = {}

SECTION_EXPLAIN['Traces'] = """
<p><b>What the masters contain.</b> For each camera the trace masters describe where each of the 648 fibers
falls on the detector. <code>lvm-mcentroids-&lt;camera&gt;.fits</code> gives the centre of each fiber's
profile, as a detector row, at every column (an array of 648 fibers &times; 4086 columns);
<code>lvm-msigmas-&lt;camera&gt;.fits</code> gives the width (Gaussian sigma, pixels) of the profile at the same
places. Both are fitted by the DRP to the dome flats of the epoch's calibration sequence. Fibers or pixels
marked in the centroid master's <code>BADPIX</code> mask are left out of everything below.</p>
<p><b>What is measured, per camera and epoch:</b></p>
<ul>
<li><b>trace_pos</b> (pixels): the median centroid over all fibers in the 101 central columns (1993&ndash;2093),
i.e. where the slit as a whole falls on the detector. Its absolute value is not meaningful; changes between
epochs show how far the image of the slit moved along the detector columns (perpendicular to the
dispersion).</li>
<li><b>trace_tilt</b> (pixels): for each fiber, the median centroid in columns 3450&ndash;3550 minus that in
columns 450&ndash;550; the median over fibers. It is the rise of the traces over about 3000 columns; a change
means the slit image rotated or bent relative to the detector.</li>
<li><b>trace_width</b> (pixels): the median sigma over all fibers in the central columns. A change means a
change of focus, or of the fitted profile shape.</li>
<li><b>trace_nmasked</b> (fibers): fibers masked along their whole length (dead or blocked fibers, or fibers
the trace fit gave up on). A jump means fibers were lost, or the fit failed.</li>
<li><b>trace_pos_rel, trace_tilt_rel</b> (pixels): the same, minus the mean of the other two spectrographs of
the same channel at the same epoch. When the whole instrument shifts, all three spectrographs move together
and these stay near zero; they change when <i>one</i> spectrograph moves on its own, which points to a problem
in, or work on, a single spectrograph. These are the versions plotted and flagged.</li>
<li><b>trace_left</b> (pixels): how much of the change since the previous epoch the science reduction would
<i>not</i> remove. The centroid array minus the previous epoch's is taken fiber by fiber and column by column;
the fibers are split into the 18 blocks of 36, and from each block the median of its difference over columns
450&ndash;3050 is subtracted. That is the one constant shift per block that the DRP applies when it extracts a
science frame (<code>_fix_fiber_thermal_shifts</code>, which cross-correlates at 11 columns between 500 and
3000). The rms of what is left, over all fibers and columns, is trace_left. A whole-block shift, however large,
gives trace_left &asymp; 0, because the DRP corrects it; what remains comes from single fibers that moved,
shifts that differ within a block, or changes along the fibers (tilt or curvature), none of which the
extraction corrects, so it goes straight into the science extraction. Typical values are 0.006&ndash;0.017
pixels; r1 is persistently higher (median 0.046 pixels), because its traces are noisier from fiber to fiber
(seen also in the <code>CompareCalib.py</code> comparisons of 60995 with 60985 and 61065).</li>
</ul>
<p><b>What is not done.</b> The traces are compared with each other, not with the flats they came from, so a
fiber mis-traced in the same way in two epochs would not show up. Checking the fits against the flats
themselves is done separately (<code>CompareCalib.py</code>, and our rebuild of 60995).</p>
"""

SECTION_EXPLAIN['Wavelengths and LSF'] = """
<p><b>What the masters contain.</b> For each camera, <code>lvm-mwave-&lt;camera&gt;.fits</code> gives the
wavelength (&Aring;) of every pixel of every fiber (648 fibers &times; 4086 columns), evaluated from the
polynomial the DRP fits to the arc lines measured in each fiber; <code>lvm-mlsf-&lt;camera&gt;.fits</code> gives
the line spread function width (FWHM) at the same places. Both come from the arc lamp exposures of the epoch's
calibration sequence (in practice the twelve 10&nbsp;s arcs). Pixels marked in the wave master's
<code>BADPIX</code> mask are left out.</p>
<p><b>What is measured, per camera and epoch:</b></p>
<ul>
<li><b>wave_c</b> (&Aring;): the median wavelength over all fibers in the 101 central columns. Changes between
epochs show how far the spectra moved along the dispersion.</li>
<li><b>wave_disp</b> (&Aring;/pixel): for each fiber, the wavelength at column 3886 minus that at column 200,
divided by the 3686 columns between them; the median over fibers. A change means the spectrum was stretched
or compressed on the detector.</li>
<li><b>lsf_c</b>: the median LSF FWHM over all fibers in the central columns. A change means a change of focus
or resolution.</li>
<li><b>wave_c_rel, wave_disp_rel, lsf_c_rel</b>: the same, minus the mean of the other two spectrographs of the
same channel at the same epoch, as for the traces; these are the versions plotted and flagged. A step here
means one spectrograph changed on its own.</li>
<li><b>wave_maxwiggle</b> (&Aring;) and <b>wave_nwiggle</b> (fibers): a check on single fibers whose solution has
gone wrong. Every 8th column, each fiber's wavelength minus the mean of its two neighbouring fibers is taken
(a single neighbour for the first and last fiber); the median of that difference along the fiber is removed,
so a constant offset does not count, and the scatter along the fiber that remains is the fiber's
<i>wiggle</i>. Smooth changes across the slit cancel, so a fiber whose polynomial bends differently from its
neighbours stands out. wave_maxwiggle is the largest wiggle in the camera (the fiber's row is given in the
table of bad masters), and any value above 0.1&nbsp;&Aring; marks the epoch as bad in itself, whatever the other
epochs. wave_nwiggle counts the fibers whose wiggle is more than 5 times the camera's median wiggle and more
than 0.03&nbsp;&Aring;. Such fibers keep their error in the science spectra, since the sky line correction only
removes a constant offset, smoothed along the slit.</li>
<li><b>wave_left</b> (&Aring;): how much of the change since the previous epoch the science reduction would
<i>not</i> remove. The wavelength array minus the previous epoch's is taken; for each fiber its value at the
DRP's reference sky lines (b: 5577&nbsp;&Aring;; r: 6364, 7359, 7392&nbsp;&Aring;; z: 8399, 8988, 9553,
9720&nbsp;&Aring;) is averaged (median) into one offset, the offsets are median filtered over 8 fibers and
smoothed along the slit with a 10-knot spline, and the result is subtracted from every column, as
<code>shift_wave_skylines</code> does to a science frame (if fewer than half the fibers have an offset, nothing
is subtracted, as in the DRP). The rms of what is left, over all fibers and columns 200&ndash;3886 (the
200 columns at each end, where the polynomials are least constrained, are left out), is wave_left. What
remains is a change of dispersion, single fibers, or structure along the spectrum; in b, where the offset is
measured at a single line, any change away from 5577&nbsp;&Aring; remains. Typical values are about
0.02&nbsp;&Aring; in b and 0.003&ndash;0.007&nbsp;&Aring; in r and z.</li>
</ul>
<p><b>What is not done.</b> The solutions are compared with each other and with neighbouring fibers, not with
the measured arc lines, so their absolute accuracy is not checked here; nor is their zero point against the
sky lines (the <code>CalibDrift.py</code> report measures that in the science frames).</p>
"""

SECTION_EXPLAIN['Bias'] = """
<p><b>What the masters contain.</b> <code>lvm-mbias-&lt;camera&gt;.fits</code> is the combined bias image of
each camera, made from the biases of (normally) the 11 nights starting at the epoch, about 95 frames. For the
numbers below the image is split into its four quadrants (the two halves in each direction), which are taken
to be the four amplifiers.</p>
<p><b>What is measured, per camera and epoch:</b></p>
<ul>
<li><b>bias_level</b> (e<sup>&minus;</sup>): the median of the whole image.</li>
<li><b>bias_qrange</b> (e<sup>&minus;</sup>): the range (largest minus smallest) of the four quadrant
medians, i.e. how unequal the amplifiers' residual levels are.</li>
<li><b>bias_sigma</b> (e<sup>&minus;</sup>): for each quadrant, the robust sigma of its pixels (from every 7th
pixel); the median over the quadrants. Since the master averages many frames, this is mostly fixed structure
in the bias (patterns, columns, glow), not the read noise of a single frame. r1 is much higher than the other
cameras (median 1.8 against 0.07&ndash;0.15&nbsp;e<sup>&minus;</sup>), from the glow-like pattern in its bias.</li>
<li><b>bias_nhot</b> (pixels): pixels more than 10 sigma above their quadrant's median. Shown but not flagged,
because it varies too much from epoch to epoch.</li>
</ul>
"""

SECTION_EXPLAIN['Twilight fiber flats (ratio to the previous epoch)'] = """
<p><b>What the masters contain.</b> <code>lvm-mfiberflat_twilight-&lt;channel&gt;.fits</code> is one file per
channel (b, r, z) with all three spectrographs stacked (1944 fibers &times; 4086 pixels): the relative
throughput of every fiber at every pixel, about 1, made from the 12 dusk twilight exposures of the epoch's
sequence. The masters in <code>$LVM_MASTER_DIR</code> have also had the DRP's <code>fiberflat-corrections</code>
step applied (spectrograph factors from sky lines in science frames). The wavelength of each pixel is taken
from the file's <code>WAVE_TRACE</code> table, and each fiber's status, telescope and label from its
<code>SLITMAP</code>; pixels marked in <code>BADPIX</code>, and non-positive values, are left out.</p>
<p><b>What is measured, per channel and epoch:</b></p>
<ul>
<li><b>flat_nodd</b> (fibers): fibers of good status (<code>fibstatus</code> 0) whose median flat over
wavelength is below 0.2 or above 2. Dividing a spectrum by such a value multiplies or divides it by a large
factor, so any such fiber marks the epoch as bad in itself, whatever the other epochs (the fibers are listed,
as fiberid(label), in the table of bad masters).</li>
</ul>
<p>The other quantities compare the flat with the previous epoch's. Both flats are resampled onto a common
0.5&nbsp;&Aring; grid and medianed in 50&nbsp;&Aring; bins (lvmdrp's <code>bin_fiberflat</code>), so fibers with
different wavelength solutions are compared at the same wavelengths; the ratio this epoch / previous epoch is
taken bin by bin, and its median over the bins gives one ratio per fiber. Only fibers of the science and the
two sky IFUs (not the standard-star fibers) with good status in both epochs are used. The per-fiber ratios,
normalized by their median, are then fitted with lvmdrp's <code>decompose_ratio</code>, which is the DRP's own
model for the fiber flat corrections (<code>fit_sky_gradient</code>): one factor per spectrograph times a
linear gradient across each IFU.</p>
<ul>
<li><b>flat_sp1, flat_sp2, flat_sp3</b>: the fitted offset of each spectrograph (a fraction; 0.01 = 1%). A
change in the overall throughput of one spectrograph, or in the spectrograph factors applied by
<code>fiberflat-corrections</code>, shows up here. Typically below 0.05%, as most consecutive epochs carry the
same factors.</li>
<li><b>flat_grad</b>: the amplitude of the fitted linear gradient across the science IFU (a fraction): a change
in how evenly the twilight sky illuminated the IFU, or in the gradient correction.</li>
<li><b>flat_left</b>: the scatter over fibers of the ratio after dividing by the fitted spectrograph factors
and gradient (a fraction): the genuine change of the fiber-to-fiber flat field. Typically 0.3&ndash;0.4%.</li>
</ul>
<p><b>What is not done.</b> Only the median over wavelength of each fiber's ratio is used, so a change in the
shape of a fiber's flat with wavelength is not measured. How flat the sky comes out across fibers in reduced
science frames, which is the real test of a fiber flat, is not checked here.</p>
"""

EXPLAIN = """<h2>How epochs are flagged</h2>
<p>Each epoch's long-term masters are reduced to a few numbers per camera (per channel for the
twilight flats), listed below, and each series of numbers over the epochs is searched for epochs that do
not fit. The scale for a series is the typical change between consecutive epochs (the robust sigma of the
differences), and with a threshold T (here %s):</p>
<ul>
<li><b>outlier</b> (red): an epoch that differs from both neighbours by more than T times the scale, in
opposite directions, while the two neighbours agree with each other: something went wrong with that epoch.</li>
<li><b>step</b> (orange): a change from the previous epoch larger than T times the scale that is not part of
an outlier: the calibration, or the instrument, changed and stayed changed.</li>
<li><b>large change</b> (purple): for the quantities that already compare an epoch with the previous one
(trace_left, wave_left, flat_*), a value more than T robust sigmas from the median.</li>
<li><b>bad fiber</b> / <b>bad flat</b>: an epoch whose masters are bad in themselves, whatever its neighbours:
a fiber whose wavelength solution wiggles by more than %s A relative to its neighbouring fibers (the fiber is
named), or fibers of good status whose twilight flat is below 0.2 or above 2 (named), which would multiply or
divide their spectra by a large factor.</li>
</ul>
<p>Changes between epochs are heavy tailed, so a change is flagged only if it is also larger than a minimum
size for its quantity (the last column of the table below), close to the 97th percentile of the changes
between consecutive epochs and of a size that matters for the data. bias_nhot is shown but not flagged; it
varies too much between epochs.</p>
<p>A step in the positions (trace_pos_rel, wave_c_rel) is usually a real change of the instrument, e.g. after
work on a spectrograph. If the change happened during an epoch, the nights after it were reduced with masters
made before it, and need large corrections: compare with the trace shifts and sky line offsets in the
CalibDrift.py report for the epoch before the step.</p>
<p>Positions move between epochs for all the spectrographs together (the instrument shifts as a whole), so
for the positions the plots and flags use the value relative to the mean of the other two spectrographs of
the same channel (the _rel quantities); a problem in one spectrograph shows up there. In the plots each
quantity is shown minus the camera's median over the epochs, so that all the panels share their axes; hover
for the value itself. A flagged epoch can be examined in detail with CompareCalib.py against its
neighbours.</p>"""


def doit(outroot='CalibEpochs', calib='', epochs_file='', only=None, table=None, thresh=4.0):
    '''
    Take the census (or read it from table), find the outliers and steps,
    and write outroot.html with its plots and tables in outroot/

    History::

        261101 ksl Coding begun
    '''
    os.makedirs(outroot, exist_ok=True)
    if table:
        xtab = Table.read(table, format='ascii.fixed_width_two_line')
        source = 'the table %s, measured from the masters by an earlier run' % table
    else:
        calib = calib or os.getenv('LVM_MASTER_DIR', '')
        epochs_file = epochs_file or default_epochs_file()
        epochs = load_epochs(epochs_file) if epochs_file else None
        if epochs is None:
            print('CalibEpochs: no epochs file: set $LVMCORE_DIR or use -epochs')
            return
        epochs = [e for e in epochs if os.path.isdir(os.path.join(calib, str(e)))]
        if not epochs:
            print('CalibEpochs: no epoch directories found in %s' % calib)
            return
        print('CalibEpochs: measuring the masters of %d epochs in %s' % (len(epochs), calib))
        xtab = measure(epochs, calib, only or PRODUCTS)
        source = 'the masters in %s (%d epochs)' % (calib, len(epochs))
    xtab = add_relative(xtab)
    tabfile = os.path.join(outroot, 'calib_epochs.txt')
    xtab.write(tabfile, format='ascii.fixed_width_two_line', overwrite=True)

    flags = find_flags(xtab, thresh)
    ftab_rows = [['epoch', 'camera', 'quantity', 'kind', 'value', 'previous epoch', 'next epoch', 'score', 'unit']]
    for f in flags:
        ftab_rows.append([f[0], f[1], f[2], f[3], '%.5g' % f[4], '%.5g' % f[5] if np.isfinite(f[5]) else '--',
                          '%.5g' % f[6] if np.isfinite(f[6]) else '--', '%.1f' % f[7], _unit(f[2])])
    flagfile = os.path.join(outroot, 'calib_epochs_flags.txt')
    Table(rows=[[str(v) for v in r] for r in ftab_rows[1:]], names=ftab_rows[0]).write(
        flagfile, format='ascii.fixed_width_two_line', overwrite=True) if flags else None

    # how often each epoch is flagged, as a summary
    epochs_all = sorted(set(int(e) for e in xtab['epoch']))
    erows = [['epoch', 'flags', 'outliers', 'steps', 'large changes', 'bad fibers or flats', 'cameras', 'quantities']]
    for e in epochs_all:
        fe = [f for f in flags if f[0] == e]
        if fe:
            erows.append([e, len(fe), sum(f[3] == 'outlier' for f in fe), sum(f[3] == 'step' for f in fe),
                          sum(f[3] == 'large change' for f in fe), sum(f[3].startswith('bad') for f in fe),
                          ' '.join(sorted(set(f[1] for f in fe))), ' '.join(sorted(set(f[2] for f in fe)))])
    erows[1:] = sorted(erows[1:], key=lambda r: -r[1])
    brows = [['epoch', 'camera', 'problem', 'value']] + [[f[0], f[1], f[3], '%.3g' % f[4]] for f in
                                                          sorted(flags, key=lambda f: (f[0], f[1])) if f[3].startswith('bad')]

    qrows = [['quantity', 'unit', 'what it is', 'smallest change flagged']]
    for q, (unit, _, desc) in list(QUANTITIES.items()) + list(FLAT_QUANTITIES.items()):
        key = '%s_rel' % q if q in RELATIVE else q
        limit = ('absolute: > %g' % ABSOLUTE[key]) if key in ABSOLUTE else \
            ('not flagged' if q in NOT_FLAGGED else '%g' % MINIMUM.get(key, 0))
        qrows.append([key, unit, desc + (', relative to the other two spectrographs' if q in RELATIVE else ''), limit])

    js = _plotly_js(outroot)
    string = xhtml.begin('Census of the calibration epochs')
    string += xhtml.add_list(['numbers from %s' % source, '%d epochs, %d to %d' % (
        len(epochs_all), epochs_all[0], epochs_all[-1]), '%d flags on %d epochs' % (len(flags), len(erows) - 1)])
    string += PROVENANCE
    string += EXPLAIN % (thresh, ABSOLUTE['wave_maxwiggle'])
    string += xhtml.table(qrows)
    string += xhtml.hline()
    string += xhtml.h2('Epochs with the most flags')
    string += xhtml.table(erows)
    string += xhtml.h2('Masters that are bad in themselves (%d)' % (len(brows) - 1))
    string += xhtml.paragraph('Fibers are given as row in the camera (wavelengths) or fiberid(label) (flats).')
    string += xhtml.table(brows)
    string += xhtml.h2('All flags, largest first')
    string += xhtml.paragraph('The score is the change in units of the scale of the series.')
    string += xhtml.table(ftab_rows[:61])
    if len(ftab_rows) > 61:
        string += xhtml.paragraph('First 60 shown; all are in %s' % xhtml.link(os.path.basename(flagfile), flagfile))
    string += xhtml.hline()
    string += js
    groups = [('Traces', [q for q in QUANTITIES if QUANTITIES[q][1] == 'trace']),
              ('Wavelengths and LSF', [q for q in QUANTITIES if QUANTITIES[q][1] == 'wave']),
              ('Bias', [q for q in QUANTITIES if QUANTITIES[q][1] == 'bias'])]
    for title, qs in groups:
        qs = [('%s_rel' % q if q in RELATIVE else q) for q in qs]
        qs = [q for q in qs if q in xtab.colnames and np.isfinite(np.asarray(xtab[q], float)).any()]
        if not qs:
            continue
        string += xhtml.h2(title)
        string += SECTION_EXPLAIN.get(title, '')
        for q in qs:
            string += plot_quantity(xtab, q, flags, CAMERAS)
    fq = [q for q in FLAT_QUANTITIES if np.isfinite(np.asarray(xtab[q], float)).any()]
    if fq:
        title = 'Twilight fiber flats (ratio to the previous epoch)'
        string += xhtml.h2(title)
        string += SECTION_EXPLAIN.get(title, '')
        for q in fq:
            string += plot_quantity(xtab, q, flags, CHANNELS)
    string += xhtml.hline()
    string += xhtml.paragraph('Tables: ' + xhtml.link(os.path.basename(tabfile), tabfile) +
                              (', ' + xhtml.link(os.path.basename(flagfile), flagfile) if flags else ''))
    string += xhtml.end()
    with open(outroot + '.html', 'w') as g:
        g.write(string)

    print('\nCalibEpochs: %d epochs, %d flags on %d epochs' % (len(epochs_all), len(flags), len(erows) - 1))
    for r in erows[1:11]:
        print('  epoch %s: %d flags (%d outliers, %d steps, %d large changes, %d bad) in %s' % tuple(r[:7]))
    for r in brows[1:]:
        print('  bad: epoch %s %s %s' % (r[0], r[1], r[2]))
    print('Report: %s.html (plots and tables in %s/)' % (outroot, outroot))


def steer(argv):
    outroot, calib, epochs_file, only, table, thresh = 'CalibEpochs', '', '', None, None, 4.0
    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-out':
            i += 1
            outroot = argv[i]
        elif argv[i] == '-calib':
            i += 1
            calib = argv[i]
        elif argv[i] == '-epochs':
            i += 1
            epochs_file = argv[i]
        elif argv[i] == '-only':
            i += 1
            only = argv[i].split(',')
            unknown = [x for x in only if x not in PRODUCTS]
            if unknown:
                print('CalibEpochs: unknown products %s; choose from %s' % (unknown, PRODUCTS))
                return
        elif argv[i] == '-table':
            i += 1
            table = argv[i]
        elif argv[i] == '-thresh':
            i += 1
            thresh = float(argv[i])
        else:
            print('Error: unknown option %s' % argv[i])
            print(_USAGE)
            return
        i += 1
    doit(outroot=outroot, calib=calib, epochs_file=epochs_file, only=only, table=table, thresh=thresh)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
