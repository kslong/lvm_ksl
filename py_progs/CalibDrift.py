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

    usage: CalibDrift.py [-h] [-out outroot] [-minexp S] dir [dir ...]
           CalibDrift.py [-h] [-out outroot] [-minexp S] [-base dir] [-every N] -drpall file
           CalibDrift.py [-h] [-out outroot] [-minexp S] [-drpall file] -table calib_drift.txt

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
        -table FILE  make the report from a calib_drift.txt table written
                     by an earlier run, without reading any headers (e.g.
                     a table made at Utah); if the table predates the use
                     of SMJD as the night, give -drpall too, and the night
                     is taken from the drpall mjd column
        -minexp S    shortest exposure time in seconds used for the
                     statistics (default 300); shorter exposures are
                     shown in grey in the plots

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
    cards that list the calibration files.  The night is the SDSS MJD
    (SMJD), not the MJD keyword: the MJD changes at 0h UT, during the
    night at LCO, so early evening exposures would otherwise be counted
    with the previous night.

    The report explains both corrections and shows, as interactive
    plotly figures with one panel per camera: each correction against
    the night, one point per exposure, with the median of each epoch;
    and the nightly medians against the nights since the epoch of the
    masters.  Tables give per camera statistics, where the variation
    comes from (within a night, between nights of an epoch, between
    epochs), the epochs in which one spectrograph departs from the
    other two, and the nights with failed or large corrections.  A short
    summary is printed to the terminal.

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

    plotly.min.js is written once into outroot/, so the report works
    offline.  Only numpy, astropy, plotly and xhtml.py are needed (not
    lvmdrp), so it can be run where the reductions are stored, e.g. at
    Utah:

        CalibDrift.py -out CalibDrift_1.3.2 -drpall drpall-1.3.2.fits

History::

    261101 ksl Coding begun, to follow the trace shifts and sky line
               offsets the DRP applies across the nights each set of
               master calibrations serves
    261101 ksl Report rebuilt as interactive plotly figures with an
               explanation of each plot and table: one panel per camera
               with shared axes, each exposure at its observation time
               (MJD + 0.3), the calibration epochs marked; tables of
               where the variation comes from, of epochs in which one
               spectrograph departs from the others and of nights with
               failed corrections.  The night is now the SDSS MJD
               (SMJD), not the MJD keyword, which had put the early
               evening exposures of 299 of 876 nights with the previous
               night.  Added -table and -minexp

'''

import os
import re
import sys
import time
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import xhtml


CAMERAS = ['b1', 'b2', 'b3', 'r1', 'r2', 'r3', 'z1', 'z2', 'z3']
SPEC_COLORS = {'1': '#1f77b4', '2': '#ff7f0e', '3': '#2ca02c'}
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
    # the night is the SDSS MJD (SMJD); MJD is the plain MJD, which changes at 0h UT, in the
    # middle of the night at LCO, so early evening exposures carry the previous day's MJD
    night = h.get('SMJD', h.get('MJD', -1))
    try:
        mjd_obs = Time(h['OBSTIME']).mjd
    except Exception:
        mjd_obs = np.nan
    row = dict(file=os.path.basename(filename), drpver=h.get('DRPVER', ''), mjd=int(night),
               header_mjd=int(h.get('MJD', -1)), mjd_obs=mjd_obs,
               expnum=int(h.get('EXPOSURE', -1)), tileid=int(h.get('TILE_ID', -1) or -1),
               wave_epoch=calib_epoch(h, 'mwave'), trace_epoch=calib_epoch(h, 'mcentroids'))
    for cam in CAMERAS:
        row['shift_%s' % cam] = float(h.get('%s FIBER SHIFT' % cam.upper(), np.nan))
        row['skyoff_%s' % cam] = float(h.get('%s WAVE SKYOFF' % cam.upper(), np.nan))
    return row


EXPLAIN_INTRO = """<p>For every science exposure the DRP starts from the long-term master calibrations of
an epoch (traces, wavelength solution, LSF, fiber flat) and then makes two small corrections for how far
the instrument has moved since those masters were made. Both are recorded, per camera, in the header of the
lvmSFrame, and this page collects them for every exposure read.</p>
<ul>
<li><b>Trace shift</b> (<code>&lt;CAM&gt; FIBER SHIFT</code>, pixels) is a shift <b>across the
fibers</b>, i.e. perpendicular to the dispersion, along the detector rows. At extraction the master
fiber traces are cross-correlated with the science image at several columns, a shift is found for each
block of 36 fibers, and the traces are moved by it before the spectra are extracted; the header gives the
median over the blocks. A positive value means the fibers sit at larger row numbers than in the masters.</li>
<li><b>Sky line offset</b> (<code>&lt;CAM&gt; WAVE SKYOFF</code>, Angstrom) is a shift <b>along the
dispersion</b>. After extraction the positions of bright sky lines (b 5577 A; r and z three or four lines)
are measured in each fiber, smoothed along the slit, and subtracted from the wavelength solution as a
constant per fiber; the header gives the mean over the fibers of each spectrograph. It is how far the
master wavelength solution, which comes from arc lamps, is from the sky on that exposure. If fewer than half
the fibers of a spectrograph give a measurement, no offset is applied and the value is exactly 0.</li>
</ul>
<p>The two are orthogonal on the detector and corrected independently. Small values mean the masters
describe the instrument well on that night. The night is the SDSS MJD (SMJD), which is the same for all the
exposures of a night; the MJD keyword changes at 0h UT, during the night at LCO.</p>"""

EXPLAIN_TIME = """<p>One panel per camera: rows are the b, r and z channels and columns spectrographs 1, 2
and 3, so the three spectrographs no longer hide each other. Each dot is one exposure (hover for the
exposure number, time, epoch and exposure time). Each exposure is placed at its observation time, MJD + 0.3,
whose integer part is the SDSS night (SMJD) and whose fraction runs through the night (about 0.2 at dusk to
0.75 at dawn), so zooming in on a night shows its exposures in time order. All nine panels share their axes, so zooming or panning
one moves them all. The line is the median of each master epoch, drawn at the
median night of the epoch. The dotted vertical lines, with a triangle at the top of each panel, mark the
nights of the master calibrations (the epochs); hover over a triangle for its MJD. Each exposure uses the
masters of the latest epoch before it. Exposures shorter than
%d s, which often have too little sky for a good measurement, are shown in grey. The y range is limited to
+-%s; %s exposures lie outside it (see the table of large corrections). Zoom by dragging, double click to
reset, and click a legend entry to hide it.</p>"""

EXPLAIN_AGE = """<p>The same corrections against the number of nights since the epoch of the masters
used, one point per night (the median of the night's exposures of at least %d s), coloured by epoch. If the
masters went stale, the points would move away from the epoch's starting value as the nights go by.</p>"""

EXPLAIN_STATS = """<p>Per camera, for exposures of at least %d s with the failures removed (|trace shift|
&gt; %s pix, |sky offset| &gt; %s A, or a sky offset of exactly 0): the median, and half the 16-84 percentile
range as a robust scatter. The last columns count the failures and large values among all exposures.</p>"""

EXPLAIN_LAYERS = """<p>Where the variation comes from, per camera: the robust scatter of the exposures of a
night about the night's median; the scatter of the nightly medians within an epoch about the epoch's median
(the median over epochs); and the scatter of the epoch medians. The last column is the median change from the
first to the last exposure of a night. If the between-epoch column is the largest, the corrections are set
mostly by which masters are used, not by drift during an epoch.</p>"""

EXPLAIN_DEPART = """<p>Epochs in which one spectrograph's median correction differs from the mean of the
other two by more than %s pix (trace shift) or %s A (sky offset). The three spectrographs share the
enclosure, so such a difference points to that spectrograph's masters for the epoch rather than to the
instrument moving.</p>"""

EXPLAIN_LARGE = """<p>The nights with the most exposures with a failed or large correction (|trace shift|
&gt; 1 pix, |sky offset| &gt; 0.5 A, or a sky offset of exactly 0 in any camera).</p>"""


def _good(xtab, prefix, cam, min_exptime, shift_lim, sky_lim):
    v = np.asarray(xtab['%s_%s' % (prefix, cam)], float)
    lim = shift_lim if prefix == 'shift' else sky_lim
    g = np.isfinite(v) & (np.abs(v) < lim) & (v != 0)
    if 'exptime' in xtab.colnames:
        g &= np.asarray(xtab['exptime'], float) >= min_exptime
    return v, g


def _half68(x):
    return 0.5 * (np.percentile(x, 84) - np.percentile(x, 16)) if len(x) > 2 else np.nan


def _robust(x):
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) > 2 else np.nan


def _plotly_js(outdir):
    """
    Write plotly.min.js into outdir once, so the report works offline
    """
    jsfile = os.path.join(outdir, 'plotly.min.js')
    if not os.path.isfile(jsfile):
        from plotly.offline import get_plotlyjs
        with open(jsfile, 'w') as g:
            g.write(get_plotlyjs())
    return '<script src="%s"></script>\n' % jsfile


def plot_time(xtab, prefix, title, ylabel, yrange, epochs, min_exptime, shift_lim, sky_lim):
    """
    A 3x3 grid of panels (channels x spectrographs) of a correction against
    the night, one point per exposure, with the epoch medians
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=3, cols=3, shared_xaxes='all', shared_yaxes='all', vertical_spacing=0.05,
                        horizontal_spacing=0.03, subplot_titles=CAMERAS)
    exptime = np.asarray(xtab['exptime'], float) if 'exptime' in xtab.colnames else np.full(len(xtab), np.nan)
    long_ = ~(exptime < min_exptime)
    # compact types keep the html small: the hover data are repeated in every panel
    custom = np.column_stack([xtab['expnum'], xtab['wave_epoch'], np.nan_to_num(exptime, nan=-1),
                              xtab['nights_since_epoch']]).astype(np.int32)
    # each exposure at its time: MJD(obstime) + 0.3, whose integer part is the SDSS night, so a night's
    # exposures spread out in time order instead of stacking at the integer night
    if 'mjd_obs' in xtab.colnames:
        night = np.where(np.isfinite(np.asarray(xtab['mjd_obs'], float)), np.asarray(xtab['mjd_obs'], float) + 0.3,
                         np.asarray(xtab['mjd'], float))
    else:
        night = np.asarray(xtab['mjd'], float)
    hover = ('exposure %{customdata[0]:.0f}<br>SMJD %{x:.3f}<br>epoch %{customdata[1]:.0f} (+%{customdata[3]:.0f} nights)'
             '<br>exptime %{customdata[2]:.0f} s<br>%{y:.3f}<extra></extra>')
    for i, cam in enumerate(CAMERAS):
        row, col = i // 3 + 1, i % 3 + 1
        v, g = _good(xtab, prefix, cam, min_exptime, shift_lim, sky_lim)
        color = SPEC_COLORS[cam[1]]
        for sel, c, name in ((long_, color, '>= %d s' % min_exptime), (~long_, '#aaaaaa', '< %d s' % min_exptime)):
            if sel.any():
                fig.add_trace(go.Scattergl(x=night[sel], y=v[sel].astype(np.float32), mode='markers',
                                           marker=dict(size=3, color=c, opacity=0.5), customdata=custom[sel],
                                           hovertemplate=hover, name='%s %s' % (cam, name), legendgroup=cam,
                                           showlegend=False), row=row, col=col)
        em, ex = [], []
        for e in epochs:
            k = g & (np.asarray(xtab['wave_epoch']) == e)
            if k.sum() >= 3:
                em.append(np.median(v[k]))
                ex.append(np.median(np.asarray(xtab['mjd'])[k]))
        fig.add_trace(go.Scatter(x=ex, y=np.round(em, 4), mode='lines+markers', line=dict(color='black', width=1.5),
                                 marker=dict(size=5), name='epoch median', legendgroup='median',
                                 showlegend=(i == 0), hovertemplate='epoch median %{y:.3f}<extra></extra>'),
                      row=row, col=col)
        # the dates of the master calibrations: a dotted line at each epoch, and a marker at the top of the
        # panel that names the epoch on hover (text labels would overlap: epochs are ~3 weeks apart)
        for e in epochs:
            fig.add_vline(x=e, line=dict(color='rgba(80,80,80,0.55)', width=1, dash='dot'), layer='above',
                          row=row, col=col)
        fig.add_trace(go.Scatter(x=epochs, y=[0.97 * yrange[1]] * len(epochs), mode='markers',
                                 marker=dict(symbol='triangle-down', size=7, color='#555555'),
                                 name='calibration epoch', legendgroup='epochs', showlegend=(i == 0),
                                 hovertemplate='master calibrations of epoch %{x}<extra></extra>'),
                      row=row, col=col)
        fig.add_hline(y=0, line=dict(color='#888888', width=1), row=row, col=col)
    fig.update_yaxes(range=yrange)
    for r, ch in enumerate('brz'):
        fig.update_yaxes(title_text='%s %s' % (ch, ylabel), row=r + 1, col=1)
    # plain numbers on the night axis, not plotly's default SI abbreviation (60.8k); integers when zoomed
    # out, fractions of the night when zoomed in
    fig.update_xaxes(tickformat='.8~g')
    for c in range(3):
        fig.update_xaxes(title_text='SMJD (MJD + 0.3; the fraction is the time in the night)', row=3, col=c + 1)
    fig.update_layout(title=title, height=950, margin=dict(l=70, r=20, t=80, b=50),
                      legend=dict(orientation='h', y=1.06))
    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_age(xtab, prefix, title, ylabel, yrange, epochs, min_exptime, shift_lim, sky_lim):
    """
    A 3x3 grid of panels of the nightly median correction against the
    nights since the epoch of the masters, coloured by epoch
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plotly.express as px
    palette = px.colors.qualitative.Dark24
    fig = make_subplots(rows=3, cols=3, shared_xaxes='all', shared_yaxes='all', vertical_spacing=0.05,
                        horizontal_spacing=0.03, subplot_titles=CAMERAS)
    night = np.asarray(xtab['mjd'])
    epoch = np.asarray(xtab['wave_epoch'])
    for i, cam in enumerate(CAMERAS):
        row, col = i // 3 + 1, i % 3 + 1
        v, g = _good(xtab, prefix, cam, min_exptime, shift_lim, sky_lim)
        for j, e in enumerate(epochs):
            k = g & (epoch == e)
            if not k.any():
                continue
            nights = np.unique(night[k])
            med = [np.median(v[k & (night == n)]) for n in nights]
            fig.add_trace(go.Scatter(x=nights - e, y=np.round(med, 4), mode='lines+markers',
                                     line=dict(color=palette[j % len(palette)], width=1), marker=dict(size=4),
                                     name='epoch %d' % e, legendgroup='e%d' % e, showlegend=(i == 0),
                                     customdata=nights,
                                     hovertemplate='epoch %d<br>night %%{customdata} (+%%{x})<br>%%{y:.3f}<extra></extra>' % e),
                          row=row, col=col)
        fig.add_hline(y=0, line=dict(color='#888888', width=1), row=row, col=col)
    fig.update_yaxes(range=yrange)
    for r, ch in enumerate('brz'):
        fig.update_yaxes(title_text='%s %s' % (ch, ylabel), row=r + 1, col=1)
    for c in range(3):
        fig.update_xaxes(title_text='nights since the epoch', row=3, col=c + 1)
    fig.update_layout(title=title, height=950, margin=dict(l=70, r=20, t=80, b=50))
    return fig.to_html(full_html=False, include_plotlyjs=False)


def camera_stats(xtab, min_exptime, shift_lim, sky_lim):
    rows = [['camera', 'trace shift median (pix)', 'robust scatter', 'n |shift| > 1', 'n |shift| > 3',
             'sky offset median (A)', 'robust scatter', 'n |offset| > 0.5', 'n offset = 0']]
    for cam in CAMERAS:
        vs, gs = _good(xtab, 'shift', cam, min_exptime, shift_lim, sky_lim)
        vo, go_ = _good(xtab, 'skyoff', cam, min_exptime, shift_lim, sky_lim)
        rows.append([cam, '%+.3f' % np.median(vs[gs]), '%.3f' % _half68(vs[gs]), int(np.sum(np.abs(vs) > 1)),
                     int(np.sum(np.abs(vs) > 3)), '%+.3f' % np.median(vo[go_]), '%.3f' % _half68(vo[go_]),
                     int(np.sum(np.abs(vo) > 0.5)), int(np.sum(vo == 0))])
    return rows


def variation_layers(xtab, prefix, min_exptime, shift_lim, sky_lim):
    rows = [['camera', 'within a night', 'nights within an epoch', 'between epochs', 'first to last of a night']]
    night = np.asarray(xtab['mjd'])
    epoch = np.asarray(xtab['wave_epoch'])
    expnum = np.asarray(xtab['expnum'])
    for cam in CAMERAS:
        v, g = _good(xtab, prefix, cam, min_exptime, shift_lim, sky_lim)
        within, firstlast, per_epoch = [], [], {}
        for n in np.unique(night[g]):
            k = g & (night == n)
            vv = v[k][np.argsort(expnum[k])]
            if len(vv) >= 4:
                within.append(vv - np.median(vv))
                firstlast.append(np.median(vv[-2:]) - np.median(vv[:2]))
            per_epoch.setdefault(epoch[k][0], []).append(np.median(vv))
        nights_in_epoch = np.median([_robust(np.array(x)) for x in per_epoch.values() if len(x) >= 5])
        between = np.std([np.median(x) for x in per_epoch.values()])
        rows.append([cam, '%.3f' % _robust(np.concatenate(within)), '%.3f' % nights_in_epoch, '%.3f' % between,
                     '%+.3f' % np.median(firstlast)])
    return rows


def departures(xtab, epochs, min_exptime, shift_lim, sky_lim, shift_thr=0.2, sky_thr=0.1):
    rows = [['epoch', 'n exposures', 'nights', 'camera', 'correction', 'this spectrograph', 'other two', 'difference']]
    epoch = np.asarray(xtab['wave_epoch'])
    for e in epochs:
        s = epoch == e
        for prefix, thr, unit in (('shift', shift_thr, 'pix'), ('skyoff', sky_thr, 'A')):
            for ch in 'brz':
                med = {}
                for sp in '123':
                    v, g = _good(xtab, prefix, ch + sp, min_exptime, shift_lim, sky_lim)
                    k = g & s
                    med[sp] = np.median(v[k]) if k.sum() > 5 else np.nan
                for sp in '123':
                    others = np.nanmean([med[o] for o in '123' if o != sp])
                    d = med[sp] - others
                    if np.isfinite(d) and abs(d) > thr:
                        rows.append([e, int(s.sum()), len(set(np.asarray(xtab['mjd'])[s])), ch + sp,
                                     'trace shift (pix)' if prefix == 'shift' else 'sky offset (A)',
                                     '%+.3f' % med[sp], '%+.3f' % others, '%+.3f' % d])
    return rows


def large_nights(xtab, nmax=20):
    bad = np.zeros(len(xtab), bool)
    for cam in CAMERAS:
        s = np.asarray(xtab['shift_%s' % cam], float)
        o = np.asarray(xtab['skyoff_%s' % cam], float)
        bad |= (np.abs(s) > 1) | (np.abs(o) > 0.5) | (o == 0)
    night = np.asarray(xtab['mjd'])
    rows = [['night', 'exposures with a failed or large correction', 'exposures that night', 'epoch',
             'median exptime of those (s)']]
    nights, counts = np.unique(night[bad], return_counts=True)
    for k in np.argsort(-counts)[:nmax]:
        n = nights[k]
        sel = bad & (night == n)
        et = np.median(np.asarray(xtab['exptime'], float)[sel]) if 'exptime' in xtab.colnames else np.nan
        rows.append([int(n), int(counts[k]), int(np.sum(night == n)), int(np.asarray(xtab['wave_epoch'])[sel][0]),
                     '%.0f' % et])
    return rows, int(bad.sum()), len(nights)


def read_headers(files):
    """
    Read the lvmSFrame headers; returns the table rows and the lists of
    missing and unreadable files
    """
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
    return rows, missing, failed


def doit(dirs, outroot='CalibDrift', drpall=None, base='', every=1, table=None, min_exptime=300):
    """
    Collect the corrections recorded in the lvmSFrame headers found in
    dirs or listed in a drpall file, or read them from a table written
    earlier, and write outroot.html with its plots and tables in outroot/

    History::

        261101 ksl Coding begun
    """
    shift_lim, sky_lim = 1.0, 0.5
    problems = []
    dtab = None
    if drpall:
        dtab = read_drpall(drpall, base=base, every=every)
    if table:
        xtab = Table.read(table, format='ascii.fixed_width_two_line')
        source = table
        if 'header_mjd' not in xtab.colnames:
            # written before the night was taken from SMJD: the mjd column is the MJD keyword
            if dtab is None:
                problems.append('the table has no header_mjd column, so its mjd column is the MJD keyword, not '
                                'the night; give -drpall to take the night from the drpall file')
            else:
                full = Table.read(drpall, hdu=1) if not drpall.endswith('.h5') else None
                if full is not None and 'mjd' in full.colnames:
                    night = dict(zip(np.asarray(full['expnum']), np.asarray(full['mjd'])))
                    xtab.rename_column('mjd', 'header_mjd')
                    xtab['mjd'] = [night.get(e, m) for e, m in zip(xtab['expnum'], xtab['header_mjd'])]
                    problems.append('nights taken from the mjd column of %s (the table predates SMJD)' %
                                    os.path.basename(drpall))
        if 'mjd_obs' not in xtab.colnames and dtab is not None:
            full = Table.read(drpall, hdu=1) if not drpall.endswith('.h5') else None
            if full is not None and 'obstime' in full.colnames:
                obs = dict(zip(np.asarray(full['expnum']),
                               Time([str(o).strip() for o in full['obstime']], format='isot').mjd))
                xtab['mjd_obs'] = [obs.get(e, np.nan) for e in xtab['expnum']]
                problems.append('observation times taken from the obstime column of %s' % os.path.basename(drpall))
        xtab['nights_since_epoch'] = xtab['mjd'] - xtab['wave_epoch']
    else:
        if dtab is not None:
            files = list(dtab['path'])
            source = '%s (%d exposures%s)' % (os.path.basename(drpall), len(files),
                                              ', every %dth' % every if every > 1 else '')
        else:
            files = find_sframes(dirs)
            source = ', '.join(dirs)
        if not files:
            print('CalibDrift: no lvmSFrame files found in %s' % (drpall or dirs))
            return
        rows, missing, failed = read_headers(files)
        if not rows:
            print('CalibDrift: none of the %d files could be read (%d missing, %d unreadable)' % (
                len(files), len(missing), len(failed)))
            return
        xtab = Table(rows=rows)
        if dtab is not None:
            dtab = dtab[[c for c in dtab.colnames if c != 'path']]
            if 'tileid' in dtab.colnames:
                dtab.rename_column('tileid', 'drpall_tileid')
            from astropy.table import join
            xtab = join(xtab, dtab, keys='expnum', join_type='left')
        xtab['nights_since_epoch'] = xtab['mjd'] - xtab['wave_epoch']
        os.makedirs(outroot, exist_ok=True)
        if missing or failed:
            probfile = os.path.join(outroot, 'unread_files.txt')
            with open(probfile, 'w') as g:
                g.write('\n'.join(['missing %s' % f for f in missing] + ['unreadable %s' % f for f in failed]) + '\n')
            problems.append('%d of %d files were missing and %d unreadable; see %s' % (
                len(missing), len(files), len(failed), xhtml.link(os.path.basename(probfile), probfile)))
    xtab.sort(['mjd', 'expnum'])
    os.makedirs(outroot, exist_ok=True)
    tabfile = os.path.join(outroot, 'calib_drift.txt')
    xtab.write(tabfile, format='ascii.fixed_width_two_line', overwrite=True)
    if 'calib_mjd' in xtab.colnames:
        cm = np.ma.filled(np.ma.asarray(xtab['calib_mjd']), -1)
        differ = (cm != -1) & (cm != np.asarray(xtab['wave_epoch']))
        if differ.any():
            problems.append('%d exposures have a drpall calib_mjd different from the epoch of the mwave file in '
                            'their header' % differ.sum())

    epochs = sorted(set(int(e) for e in xtab['wave_epoch'] if e > 0))
    versions = sorted(set(str(v) for v in xtab['drpver']))
    n_out_shift = int(sum(np.sum(np.abs(np.asarray(xtab['shift_%s' % c], float)) > 1.2) for c in CAMERAS))
    n_out_sky = int(sum(np.sum(np.abs(np.asarray(xtab['skyoff_%s' % c], float)) > 0.6) for c in CAMERAS))
    stats = camera_stats(xtab, min_exptime, shift_lim, sky_lim)
    layers_shift = variation_layers(xtab, 'shift', min_exptime, shift_lim, sky_lim)
    layers_sky = variation_layers(xtab, 'skyoff', min_exptime, shift_lim, sky_lim)
    depart = departures(xtab, epochs, min_exptime, shift_lim, sky_lim)
    large, nbad, nbadnights = large_nights(xtab)
    erows = [['epoch of masters', 'n exposures', 'nights', 'night range', 'nights since epoch']]
    for e in epochs:
        s = np.asarray(xtab['wave_epoch']) == e
        n = np.asarray(xtab['mjd'])[s]
        erows.append([e, int(s.sum()), len(set(n)), '%d-%d' % (n.min(), n.max()), '%d-%d' % (n.min() - e, n.max() - e)])

    js = _plotly_js(outroot)
    string = xhtml.begin('Drift from the master calibrations')
    string += xhtml.add_list(['%d exposures from %s' % (len(xtab), source),
                              'DRP versions: %s' % ', '.join(versions),
                              'nights %d-%d (%d nights), %d master epochs' % (
                                  xtab['mjd'].min(), xtab['mjd'].max(), len(set(xtab['mjd'])), len(epochs))]
                             + problems)
    string += EXPLAIN_INTRO
    string += xhtml.hline()
    string += xhtml.h2('Corrections per camera')
    string += EXPLAIN_STATS % (min_exptime, shift_lim, sky_lim)
    string += xhtml.table(stats)
    string += xhtml.hline()
    string += xhtml.h2('Trace shift across the fibers, against night')
    string += EXPLAIN_TIME % (min_exptime, '1.2 pix', n_out_shift)
    string += js
    string += plot_time(xtab, 'shift', 'Trace shift applied at extraction (FIBER SHIFT)', 'shift (pix)', [-1.2, 1.2],
                        epochs, min_exptime, shift_lim, sky_lim)
    string += xhtml.h2('Sky line offset along the dispersion, against night')
    string += EXPLAIN_TIME % (min_exptime, '0.6 A', n_out_sky)
    string += plot_time(xtab, 'skyoff', 'Sky line offset removed from the wavelengths (WAVE SKYOFF)', 'offset (A)',
                        [-0.6, 0.6], epochs, min_exptime, shift_lim, sky_lim)
    string += xhtml.hline()
    string += xhtml.h2('Against nights since the epoch of the masters')
    string += EXPLAIN_AGE % min_exptime
    string += plot_age(xtab, 'shift', 'Nightly median trace shift against nights since the epoch', 'shift (pix)',
                       [-1.2, 1.2], epochs, min_exptime, shift_lim, sky_lim)
    string += plot_age(xtab, 'skyoff', 'Nightly median sky line offset against nights since the epoch', 'offset (A)',
                       [-0.6, 0.6], epochs, min_exptime, shift_lim, sky_lim)
    string += xhtml.hline()
    string += xhtml.h2('Where the variation comes from')
    string += EXPLAIN_LAYERS
    string += xhtml.h3('Trace shift (pix)')
    string += xhtml.table(layers_shift)
    string += xhtml.h3('Sky line offset (A)')
    string += xhtml.table(layers_sky)
    string += xhtml.hline()
    string += xhtml.h2('Epochs where one spectrograph departs from the others (%d)' % (len(depart) - 1))
    string += EXPLAIN_DEPART % ('0.2', '0.1')
    string += xhtml.table(depart)
    string += xhtml.hline()
    string += xhtml.h2('Nights with failed or large corrections')
    string += xhtml.paragraph('%d exposures on %d nights.' % (nbad, nbadnights))
    string += EXPLAIN_LARGE
    string += xhtml.table(large)
    string += xhtml.hline()
    string += xhtml.h2('Master epochs')
    string += xhtml.table(erows)
    string += xhtml.paragraph('Table of all exposures: ' + xhtml.link(os.path.basename(tabfile), tabfile))
    string += xhtml.end()
    with open(outroot + '.html', 'w') as g:
        g.write(string)

    print('\nCalibDrift: %d exposures, nights %d-%d, %d master epochs, versions %s' % (
        len(xtab), xtab['mjd'].min(), xtab['mjd'].max(), len(epochs), ', '.join(versions)))
    for p_ in problems:
        print('  ' + re.sub('<[^>]*>', '', p_))
    print('  trace shift median (pix): ' + '; '.join('%s %s' % (r[0], r[1]) for r in stats[1:]))
    print('  sky offset median (A):    ' + '; '.join('%s %s' % (r[0], r[5]) for r in stats[1:]))
    print('  %d epochs x cameras where one spectrograph departs from the others; %d exposures with failed or '
          'large corrections on %d nights' % (len(depart) - 1, nbad, nbadnights))
    print('Report: %s.html (plots and table in %s/)' % (outroot, outroot))


def steer(argv):
    '''
    Just a steering routine
    '''
    outroot = 'CalibDrift'
    drpall = None
    base = ''
    every = 1
    table = None
    min_exptime = 300
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
        elif argv[i] == '-table':
            i += 1
            table = argv[i]
        elif argv[i] == '-minexp':
            i += 1
            min_exptime = float(argv[i])
        elif argv[i][0] == '-':
            print('Error: unknown option %s' % argv[i])
            print(_USAGE)
            return
        else:
            dirs.append(argv[i])
        i += 1
    if not dirs and not drpall and not table:
        print(_USAGE)
        return
    doit(dirs, outroot=outroot, drpall=drpall, base=base, every=every, table=table, min_exptime=min_exptime)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
