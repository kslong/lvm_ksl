#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Analyse the per-exposure fiber-flat planes measured by
    SkyLinePlane.py: check them against the gradient the DRP recorded in
    the master flats, model what is left with the telescope focus, the
    K-mirror angle and the time since the ffactor window, and write an
    html report with plots and tables.

Command line usage (if any):

    usage: SkyLinePlaneAnalysis.py [-h] [-out outroot] [-minexp S] [-minglat D] [-all] table.txt

    Arguments::

        table.txt    a table written by SkyLinePlane.py

    Options::

        -h           print this help and exit
        -out outroot name of the html file (outroot.html) and of the
                     subdirectory with its plots and tables (outroot/);
                     default SkyLinePlaneAnalysis
        -minexp S    shortest exposure time used (default 300 s)
        -minglat D   smallest |galactic latitude| used (default 15 deg),
                     to keep [OI] 6364 from the target out of the r arm
        -all         keep exposures with DRPQUAL != 0 (default: drop them)

Description::

    For each exposure and arm SkyLinePlane.py gives the plane fitted to
    the sky line (grad, dir) and the plane predicted from the GCOEFF of
    the master flat (pgrad, pdir).  In a DRP 1.3.2 CFrame the predicted
    plane is the gradient the DRP measured but did not apply, so

        residual = measured plane - predicted plane

    (as x,y vectors, % per 10 mm) is what is left after a perfect
    application of the recorded correction, i.e. after FixFlatGradient.py
    (for a .gradfix table the predicted plane is not subtracted).  The
    report:

    1.  checks the method: inside the ffactor science windows the mean
        residual should be about zero
    2.  shows the measured and predicted gradient, and the residual
        x and y components, against MJD, with the ffactor epochs and
        calibration epochs marked
    3.  fits each residual component with nested linear models:
            c0 + c_f*(focus - 35.5) + c_s*sin(K) + c_c*cos(K) + c_t*dt/100
        where focus is TESCIFO, K the K-mirror angle TESCIKM and dt the
        days outside the ffactor science window (0 inside), plus a
        version with one intercept per ffactor epoch; and gives the
        coefficients, their errors and the rms before and after
    4.  separates focus from season: the same focus fit to the values
        minus their nightly means, which uses only the focus changes
        within a night
    5.  compares dark (Moon below the horizon or illumination < 0.3)
        with bright time
    6.  splits the scatter left after the full model into within a
        night, between nights of a calibration epoch, and between
        calibration epochs, and tabulates the median residual of each
        calibration epoch, to look for steps
    7.  shows the spectrograph imbalance left (measured minus
        predicted sp1 and sp3) against focus and MJD
    8.  at bright Moon only, compares the continuum plane with the line
        plane (a check that the effect is a throughput change)

    Exposures are dropped if their DRPQUAL is not 0, they are shorter
    than -minexp, |b| < -minglat, fewer than 500 fibers were fitted, or
    the residual is more than 6 MAD from the median of its arm (listed
    in outroot/outliers.txt).

Primary routines:

    doit
    steer

Notes::

    Only numpy, astropy, matplotlib and xhtml.py are needed.  The
    ffactor science windows are read from lvmcore
    fiberflat-factor-epochs.yaml if $LVMCORE_DIR is set, else the copy
    in FF_WINDOWS below is used.

    The report text is ASCII only (xhtml.begin declares ISO-8859-1).

History::

    261105 ksl Coding begun, to analyse the Utah run of SkyLinePlane.py
               over the 1.3.2 science exposures
    261106 ksl Text columns holding numbers (None for empty header
               keywords) are converted to float, and placeholders
               such as -999 set to NaN

'''

import os
import re
import sys
import numpy as np
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # noqa: E402
import xhtml                      # noqa: E402

ARMS = 'brz'
FOCUS0 = 35.5
# ffactor epoch: (first, last) MJD of its science window (lvmcore, 2026-11)
FF_WINDOWS = {60177: (60255, 60284), 60339: (60339, 60368), 60465: (60465, 60494),
              60573: (60573, 60602), 60744: (60750, 60779), 60985: (60985, 61014)}
MODELS = [('focus', ['f']), ('K-mirror', ['s', 'c']), ('focus + K-mirror', ['f', 's', 'c']),
          ('focus + K-mirror + dt', ['f', 's', 'c', 't']),
          ('focus + K-mirror + dt + epoch intercepts', ['f', 's', 'c', 't', 'e'])]


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before the History: line, so -h stays short
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


_USAGE = _usage_from_doc(__doc__)


def ff_windows():
    '''
    The ffactor epochs and their science windows
    '''
    core = os.environ.get('LVMCORE_DIR', '')
    fname = os.path.join(core, 'calibrations', 'fiberflat-factor-epochs.yaml')
    if core and os.path.isfile(fname):
        try:
            import yaml
            with open(fname) as f:
                ep = yaml.safe_load(f)['epochs']
            out = {}
            for k, v in ep.items():
                lo, hi = (int(s) for s in str(v['flavors']['science']).split('-'))
                out[int(k)] = (lo, hi)
            return out
        except Exception as err:
            print('SkyLinePlaneAnalysis: could not read %s (%s), using the built-in windows' % (fname, err))
    return FF_WINDOWS


def _robust_sigma(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) else np.nan


def to_numeric(xtab, keep=('file', 'drpver', 'fluxcal', 'gradfix', 'why')):
    '''
    Convert text columns that hold numbers (with None or empty entries,
    as written by earlier SkyLinePlane.py versions) to float, NaN for
    the entries that are not numbers
    '''
    for c in xtab.colnames:
        if c in keep or xtab[c].dtype.kind not in 'US':
            continue
        v = np.full(len(xtab), np.nan)
        nbad = 0
        for i, x in enumerate(np.asarray(xtab[c]).astype(str)):
            try:
                v[i] = float(x)
            except ValueError:
                nbad += 1
        if nbad < len(xtab):
            xtab[c] = v
            if nbad:
                print('SkyLinePlaneAnalysis: column %s: %d entries are not numbers, set to NaN' % (c, nbad))
    return xtab


def prepare(xtab, minexp, minglat, keep_all):
    '''
    Add the night, ffactor epoch, dt, focus and the residual components;
    returns the selected table and the dropped exposures with reasons
    '''
    xtab = xtab.copy()
    night = xtab['smjd'] if 'smjd' in xtab.colnames else xtab['mjd']
    xtab['night'] = np.where(np.asarray(night) > 0, night, xtab['mjd'])
    win = ff_windows()
    ffs = np.array(sorted(win))
    i = np.clip(np.searchsorted(ffs, xtab['mjd'], side='right') - 1, 0, len(ffs) - 1)
    xtab['ff'] = ffs[i]
    lo = np.array([win[k][0] for k in xtab['ff']])
    hi = np.array([win[k][1] for k in xtab['ff']])
    m = np.asarray(xtab['mjd'])
    xtab['dt'] = np.where(m < lo, m - lo, np.where(m > hi, m - hi, 0))
    foc = np.asarray(xtab['focus'], float) if 'focus' in xtab.colnames else np.full(len(xtab), np.nan)
    if 'drpall_sci_focpos' in xtab.colnames:
        alt = np.asarray(xtab['drpall_sci_focpos'], float)
        foc = np.where(np.isfinite(foc), foc, alt)
    xtab['foc'] = foc
    gradfix = 'gradfix' in xtab.colnames and np.all(np.asarray(xtab['gradfix']).astype(str) == 'True')
    for arm in ARMS:
        g, d = np.asarray(xtab[arm + '_grad'], float), np.radians(np.asarray(xtab[arm + '_dir'], float))
        pg, pd = np.asarray(xtab[arm + '_pgrad'], float), np.radians(np.asarray(xtab[arm + '_pdir'], float))
        xtab[arm + '_mx'], xtab[arm + '_my'] = g * np.cos(d), g * np.sin(d)
        xtab[arm + '_px'], xtab[arm + '_py'] = pg * np.cos(pd), pg * np.sin(pd)
        if gradfix:
            xtab[arm + '_dx'], xtab[arm + '_dy'] = xtab[arm + '_mx'], xtab[arm + '_my']
            xtab[arm + '_dsp1'], xtab[arm + '_dsp3'] = xtab[arm + '_sp1'], xtab[arm + '_sp3']
        else:
            xtab[arm + '_dx'] = xtab[arm + '_mx'] - xtab[arm + '_px']
            xtab[arm + '_dy'] = xtab[arm + '_my'] - xtab[arm + '_py']
            xtab[arm + '_dsp1'] = np.asarray(xtab[arm + '_sp1'], float) - np.asarray(xtab[arm + '_psp1'], float)
            xtab[arm + '_dsp3'] = np.asarray(xtab[arm + '_sp3'], float) - np.asarray(xtab[arm + '_psp3'], float)
    why = [''] * len(xtab)

    def flag(mask, reason):
        for k in np.where(mask)[0]:
            if not why[k]:
                why[k] = reason

    if not keep_all and 'drpqual' in xtab.colnames:
        flag(np.asarray(xtab['drpqual']) != 0, 'drpqual')
    if 'exptime' in xtab.colnames:
        flag(np.asarray(xtab['exptime'], float) < minexp, 'exptime')
    flag(~(np.abs(np.asarray(xtab['glat'], float)) >= minglat), 'glat')
    for arm in ARMS:
        flag(~(np.asarray(xtab[arm + '_n'], float) >= 500), 'nfib_' + arm)
        flag(~np.isfinite(np.asarray(xtab[arm + '_dx'], float)), 'nan_' + arm)
    ok = np.array([w == '' for w in why])
    for arm in ARMS:
        for c in ('_dx', '_dy'):
            v = np.asarray(xtab[arm + c], float)
            med, sig = np.median(v[ok]), _robust_sigma(v[ok])
            flag(ok & (np.abs(v - med) > 6 * sig), 'outlier_' + arm + c[1:])
    xtab['why'] = why
    sel = xtab[np.array([w == '' for w in why])]
    return sel, xtab[np.array([w != '' for w in why])], gradfix


def design(xtab, terms):
    '''
    Design matrix for the nested models (always with an intercept, or
    one intercept per ffactor epoch if 'e' is in terms)
    '''
    km = np.radians(np.asarray(xtab['kmir'], float))
    cols, names = [], []
    if 'e' in terms:
        for ff in sorted(set(xtab['ff'])):
            cols.append((np.asarray(xtab['ff']) == ff).astype(float))
            names.append('c0[%d]' % ff)
    else:
        cols.append(np.ones(len(xtab)))
        names.append('c0')
    if 'f' in terms:
        cols.append(np.asarray(xtab['foc'], float) - FOCUS0)
        names.append('c_f')
    if 's' in terms:
        cols.append(np.sin(km))
        names.append('c_s')
    if 'c' in terms:
        cols.append(np.cos(km))
        names.append('c_c')
    if 't' in terms:
        cols.append(np.asarray(xtab['dt'], float) / 100.)
        names.append('c_t')
    return np.vstack(cols).T, names


def linfit(A, v):
    '''
    Least squares with errors; returns coefficients, errors, residuals
    '''
    good = np.all(np.isfinite(A), axis=1) & np.isfinite(v)
    c, *_ = np.linalg.lstsq(A[good], v[good], rcond=None)
    res = np.full(len(v), np.nan)
    res[good] = v[good] - A[good] @ c
    dof = max(1, good.sum() - A.shape[1])
    s2 = np.nansum(res ** 2) / dof
    try:
        err = np.sqrt(np.diag(np.linalg.inv(A[good].T @ A[good])) * s2)
    except np.linalg.LinAlgError:
        err = np.full(len(c), np.nan)
    return c, err, res


def nightly_demean(xtab, v):
    '''
    v minus its nightly mean, NaN for nights with one exposure
    '''
    v = np.asarray(v, float)
    out = np.full(len(v), np.nan)
    nights = np.asarray(xtab['night'])
    for n in np.unique(nights):
        m = nights == n
        if m.sum() > 1:
            out[m] = v[m] - np.nanmean(v[m])
    return out


def variance_split(xtab, res):
    '''
    rms of the residual within a night, of nightly means about their
    calibration epoch's mean, and of calibration-epoch means
    '''
    res = np.asarray(res, float)
    nights, epochs = np.asarray(xtab['night']), np.asarray(xtab['epoch'])
    within, nmeans, nep = [], [], []
    for n in np.unique(nights):
        m = (nights == n) & np.isfinite(res)
        if m.sum() > 1:
            within += list(res[m] - res[m].mean())
        if m.sum():
            nmeans.append(res[m].mean())
            nep.append(epochs[m][0])
    nmeans, nep = np.array(nmeans), np.array(nep)
    between_n, emeans = [], []
    for e in np.unique(nep):
        m = nep == e
        if m.sum() > 1:
            between_n += list(nmeans[m] - nmeans[m].mean())
        emeans.append(nmeans[m].mean())
    q = lambda a: np.sqrt(np.mean(np.square(a))) if len(a) else np.nan  # noqa: E731
    return q(within), len(within), q(between_n), len(between_n), np.std(emeans), len(emeans)


def _mark_epochs(ax, xtab):
    for ff in sorted(set(xtab['ff'])):
        ax.axvline(ff, color='k', ls='--', lw=0.8)
    for e in sorted(set(np.asarray(xtab['epoch']))):
        if e > 0:
            ax.axvline(e, color='0.75', ls=':', lw=0.6)


def plot_time(xtab, outdir):
    '''
    Measured and predicted gradient, and the residual components, vs MJD
    '''
    fig, ax = plt.subplots(3, 3, figsize=(16, 11), sharex=True, layout='constrained')
    for j, arm in enumerate(ARMS):
        m = np.asarray(xtab['mjd'])
        ax[0, j].plot(m, xtab[arm + '_grad'], '.', ms=2, color='C0', label='measured')
        ax[0, j].plot(m, xtab[arm + '_pgrad'], '.', ms=2, color='C3', label='predicted (GCOEFF)')
        ax[0, j].set_title('%s arm' % arm)
        for i, (c, lab) in enumerate((('_dx', 'residual x'), ('_dy', 'residual y'))):
            sc = ax[i + 1, j].scatter(m, xtab[arm + c], c=xtab['foc'], s=3, cmap='viridis')
            ax[i + 1, j].axhline(0, color='k', lw=0.5)
            ax[i + 1, j].set_ylabel('%s (%% / 10 mm)' % lab)
        for a in ax[:, j]:
            _mark_epochs(a, xtab)
            a.grid(alpha=0.3)
        ax[2, j].set_xlabel('MJD')
    ax[0, 0].set_ylabel('|gradient| (% / 10 mm)')
    ax[0, 0].legend(fontsize=8, markerscale=4)
    fig.colorbar(sc, ax=ax[1:, :].ravel().tolist(), label='focus TESCIFO')
    fig.suptitle('dashed: ffactor epochs; dotted: calibration epochs')
    fig.savefig(os.path.join(outdir, 'time.png'), dpi=90)
    plt.close(fig)


def plot_drivers(xtab, outdir):
    '''
    Residual y against focus, residual x against K-mirror angle, coloured by dt
    '''
    fig, ax = plt.subplots(2, 3, figsize=(16, 9), layout='constrained')
    for j, arm in enumerate(ARMS):
        sc = ax[0, j].scatter(xtab['foc'], xtab[arm + '_dy'], c=xtab['dt'], s=3, cmap='plasma')
        ax[0, j].set_xlabel('focus TESCIFO')
        ax[0, j].set_ylabel('residual y (% / 10 mm)')
        ax[0, j].set_title('%s arm' % arm)
        ax[1, j].scatter(xtab['kmir'], xtab[arm + '_dx'], c=xtab['dt'], s=3, cmap='plasma')
        ax[1, j].set_xlabel('K-mirror angle TESCIKM (deg)')
        ax[1, j].set_ylabel('residual x (% / 10 mm)')
        for a in ax[:, j]:
            a.axhline(0, color='k', lw=0.5)
            a.grid(alpha=0.3)
    fig.colorbar(sc, ax=ax.ravel().tolist(), label='days outside the ffactor window')
    fig.savefig(os.path.join(outdir, 'drivers.png'), dpi=90)
    plt.close(fig)


def plot_after(xtab, res_full, outdir):
    '''
    Nightly medians of the residual after the full model, vs MJD
    '''
    fig, ax = plt.subplots(2, 3, figsize=(16, 8), sharex=True, layout='constrained')
    nights = np.asarray(xtab['night'])
    un = np.unique(nights)
    for j, arm in enumerate(ARMS):
        for i, c in enumerate(('dx', 'dy')):
            r = res_full[(arm, c)]
            med = [np.nanmedian(r[nights == n]) for n in un]
            ax[i, j].plot(un, med, '.', ms=3)
            ax[i, j].axhline(0, color='k', lw=0.5)
            _mark_epochs(ax[i, j], xtab)
            ax[i, j].set_ylabel('nightly median, residual %s\nafter full model (%% / 10 mm)' % c[1])
            ax[i, j].grid(alpha=0.3)
        ax[0, j].set_title('%s arm' % arm)
        ax[1, j].set_xlabel('MJD')
    fig.savefig(os.path.join(outdir, 'after_model.png'), dpi=90)
    plt.close(fig)


def plot_balance(xtab, outdir):
    '''
    Spectrograph imbalance left (sp1, sp3) vs focus and MJD
    '''
    fig, ax = plt.subplots(2, 3, figsize=(16, 8), layout='constrained')
    for j, arm in enumerate(ARMS):
        for s, col in (('1', 'C0'), ('3', 'C2')):
            ax[0, j].plot(xtab['foc'], xtab['%s_dsp%s' % (arm, s)], '.', ms=2, color=col, label='sp' + s)
            ax[1, j].plot(xtab['mjd'], xtab['%s_dsp%s' % (arm, s)], '.', ms=2, color=col, label='sp' + s)
        ax[0, j].set_xlabel('focus TESCIFO')
        ax[1, j].set_xlabel('MJD')
        _mark_epochs(ax[1, j], xtab)
        ax[0, j].set_title('%s arm' % arm)
        for a in ax[:, j]:
            a.axhline(0, color='k', lw=0.5)
            a.set_ylabel('imbalance left (%)')
            a.grid(alpha=0.3)
    ax[0, 0].legend(fontsize=8, markerscale=4)
    fig.savefig(os.path.join(outdir, 'balance.png'), dpi=90)
    plt.close(fig)


def _f(x, fmt='%.2f'):
    return fmt % x if np.isfinite(x) else '-'


def doit(table, outroot='SkyLinePlaneAnalysis', minexp=300., minglat=15., keep_all=False):
    '''
    Description:
        Read a SkyLinePlane.py table, analyse it and write outroot.html
        with its plots and tables in outroot/

    Notes:
        A short summary is printed to the terminal

    History:
        261105 ksl Coding begun
    '''
    xtab = Table.read(table, format='ascii.fixed_width_two_line')
    xtab = to_numeric(xtab)
    for c in ('focus', 'skye_focus', 'skyw_focus', 'kmir', 'pa', 'alt', 'labtemp', 'ccdtemp', 'moon_fli',
              'moon_alt', 'moon_sep', 'sun_alt', 'drpall_sci_focpos', 'drpall_sci_kmpos'):
        if c in xtab.colnames and xtab[c].dtype.kind == 'f':
            bad = np.asarray(xtab[c]) <= -900
            if bad.any():
                xtab[c][bad] = np.nan
                print('SkyLinePlaneAnalysis: column %s: %d placeholder values (<= -900) set to NaN' % (c, bad.sum()))
    sel, dropped, gradfix = prepare(xtab, minexp, minglat, keep_all)
    outdir = outroot
    os.makedirs(outdir, exist_ok=True)
    if len(dropped):
        dropped['file', 'expnum', 'mjd', 'why'].write(os.path.join(outdir, 'outliers.txt'),
                                                     format='ascii.fixed_width_two_line', overwrite=True)
    else:
        with open(os.path.join(outdir, 'outliers.txt'), 'w') as f:
            f.write('# no exposures dropped\n')
    nd = {}
    for w in dropped['why']:
        nd[w.split('_')[0]] = nd.get(w.split('_')[0], 0) + 1
    print('SkyLinePlaneAnalysis: %d exposures, %d used, dropped %s' % (len(xtab), len(sel), nd))
    if len(sel) < 10:
        print('SkyLinePlaneAnalysis: too few exposures to analyse')
        return

    page = xhtml.begin('Per-exposure fiber-flat planes from sky lines')
    page += xhtml.paragraph(
        'Table %s: %d exposures, %d used (MJD %d-%d, %d nights, %d calibration epochs). '
        'Dropped: %s; the list is in <a href="%s/outliers.txt">outliers.txt</a>. %s'
        % (os.path.basename(table), len(xtab), len(sel), np.min(sel['mjd']), np.max(sel['mjd']),
           len(set(sel['night'])), len(set(sel['epoch'])),
           ', '.join('%s %d' % kv for kv in sorted(nd.items())) or 'none', outdir,
           'The CFrames are FixFlatGradient.py copies, so the residual is the measured plane itself.'
           if gradfix else
           'Residual = measured plane minus the plane predicted from the GCOEFF of the master flat, '
           'i.e. what FixFlatGradient.py (or the DRP fix) would leave.'))
    page += xhtml.paragraph(
        'Planes are vectors in the IFU (x = xpmm, y = ypmm) in %% per 10 mm. sp1/sp3 = mean of the fitted '
        'model over the Sci fibers of spectrograph 1/3, minus 1. focus = TESCIFO, K = K-mirror angle '
        'TESCIKM, dt = days outside the ffactor science window (0 inside). Focus is referred to %.1f.'
        % FOCUS0)

    # 1. validation
    rows = [['arm', 'N in window', 'mean residual x', 'mean residual y', 'rms x', 'rms y',
             'median |measured| - |predicted|', 'median formal error']]
    inw = sel[sel['dt'] == 0]
    for arm in ARMS:
        rows.append([arm, len(inw), _f(np.mean(inw[arm + '_dx'])), _f(np.mean(inw[arm + '_dy'])),
                     _f(np.std(inw[arm + '_dx'])), _f(np.std(inw[arm + '_dy'])),
                     _f(np.median(inw[arm + '_grad'] - inw[arm + '_pgrad'])), _f(np.median(inw[arm + '_egrad']))])
    page += xhtml.hline() + xhtml.h2('1. Does the measurement reproduce the DRP inside the ffactor windows?')
    page += xhtml.paragraph('Inside the science windows from which the DRP measured GCOEFF, the mean residual '
                            'should be close to zero. The rms is the exposure-to-exposure scatter, which the '
                            'DRP averages over.')
    page += xhtml.table(rows)

    # 2. time series
    plot_time(sel, outdir)
    page += xhtml.hline() + xhtml.h2('2. Gradient and residual against time')
    page += xhtml.paragraph('Top: measured and predicted |gradient|. Middle and bottom: residual x and y, '
                            'coloured by focus.')
    page += xhtml.image('%s/time.png' % outdir, width=1400, height=960)

    # 3. models
    page += xhtml.hline() + xhtml.h2('3. What drives the residual: nested linear models')
    plot_drivers(sel, outdir)
    page += xhtml.image('%s/drivers.png' % outdir, width=1400, height=790)
    res_full, summary = {}, []
    mrows = [['arm', 'comp', 'model', 'coefficients (+- error)', 'rms before', 'rms after']]
    for arm in ARMS:
        for c in ('dx', 'dy'):
            v = np.asarray(sel['%s_%s' % (arm, c)], float)
            for name, terms in MODELS:
                A, names = design(sel, terms)
                coef, err, res = linfit(A, v)
                txt = ', '.join('%s %.2f+-%.2f' % (n, a, e) for n, a, e in zip(names, coef, err)
                                if not n.startswith('c0[') and n != 'c0')
                if 'e' in terms:
                    txt += '; epoch intercepts ' + ' '.join('%.1f' % a for n, a in zip(names, coef) if n.startswith('c0['))
                mrows.append([arm, c[1], name, txt, _f(np.nanstd(v)), _f(np.nanstd(res))])
                if name == 'focus + K-mirror + dt':
                    res_full[(arm, c)] = res
                    summary.append((arm, c, dict(zip(names, coef)), np.nanstd(v), np.nanstd(res)))
    page += xhtml.paragraph('Each residual component is fitted with c0 + c_f*(focus - %.1f) + c_s*sin(K) + '
                            'c_c*cos(K) + c_t*dt/100, and nested subsets. c_f is in %% per 10 mm per unit of '
                            'focus, c_t in %% per 10 mm per 100 days.' % FOCUS0)
    page += xhtml.table(mrows)

    # 4. focus vs season
    frows = [['arm', 'comp', 'N (nights with >1 exposure)', 'c_f from within-night changes', 'c_f from all',
              'rms of within-night focus changes']]
    for arm in ARMS:
        for c in ('dx', 'dy'):
            v = nightly_demean(sel, sel['%s_%s' % (arm, c)])
            f = nightly_demean(sel, sel['foc'])
            g = np.isfinite(v) & np.isfinite(f)
            if g.sum() > 10 and np.std(f[g]) > 0:
                cf, ef, _ = linfit(f[g][:, None], v[g])
                cw = '%.2f+-%.2f' % (cf[0], ef[0])
            else:
                cw = '-'
            A, names = design(sel, ['f'])
            ca, ea, _ = linfit(A, np.asarray(sel['%s_%s' % (arm, c)], float))
            frows.append([arm, c[1], g.sum(), cw, '%.2f+-%.2f' % (ca[1], ea[1]), _f(np.std(f[g]) if g.any() else np.nan)])
    page += xhtml.hline() + xhtml.h2('4. Focus or season?')
    page += xhtml.paragraph('Focus follows the temperature, so it also follows the season. Within a night the '
                            'season does not change, so a focus slope fitted to the values minus their nightly '
                            'means measures the focus effect alone. If it agrees with the slope from all the '
                            'data, the effect is the focus.')
    page += xhtml.table(frows)

    # 5. dark vs bright
    moon_alt = np.asarray(sel['moon_alt'], float)
    fli = np.asarray(sel['moon_fli'], float)
    dark = (moon_alt < 0) | (fli < 0.3)
    drows = [['arm', 'comp', 'time', 'N', 'mean', 'c_f (focus only)', 'c_s (K-mirror only)']]
    for arm in ARMS:
        for c in ('dx', 'dy'):
            for lab, m in (('dark', dark), ('bright', ~dark)):
                s = sel[m]
                if len(s) < 10:
                    continue
                v = np.asarray(s['%s_%s' % (arm, c)], float)
                cf, ef, _ = linfit(design(s, ['f'])[0], v)
                cs, es, _ = linfit(design(s, ['s', 'c'])[0], v)
                drows.append([arm, c[1], lab, len(s), _f(np.mean(v)), '%.2f+-%.2f' % (cf[1], ef[1]),
                              '%.2f+-%.2f' % (cs[1], es[1])])
    page += xhtml.hline() + xhtml.h2('5. Dark and bright time')
    page += xhtml.paragraph('Dark = Moon below the horizon or illumination below 0.3. The sky line is a '
                            'difference from the local continuum, so moonlight should not matter; this checks it.')
    page += xhtml.table(drows)

    # 6. what is left
    vrows = [['arm', 'comp', 'rms after full model', 'within a night (N)', 'nightly means within a calibration epoch (N)',
              'calibration-epoch means (N)']]
    for arm in ARMS:
        for c in ('dx', 'dy'):
            w, nw, b, nb, e, ne = variance_split(sel, res_full[(arm, c)])
            vrows.append([arm, c[1], _f(np.nanstd(res_full[(arm, c)])), '%s (%d)' % (_f(w), nw),
                          '%s (%d)' % (_f(b), nb), '%s (%d)' % (_f(e), ne)])
    erows = [['calibration epoch', 'ffactor epoch', 'N'] + ['%s %s' % (a, c) for a in ARMS for c in ('x', 'y')]]
    for e in sorted(set(sel['epoch'])):
        m = np.asarray(sel['epoch']) == e
        erows.append([e, sel['ff'][m][0], m.sum()] +
                     [_f(np.nanmedian(res_full[(a, c)][m])) for a in ARMS for c in ('dx', 'dy')])
    plot_after(sel, res_full, outdir)
    page += xhtml.hline() + xhtml.h2('6. What is left after focus + K-mirror + dt')
    page += xhtml.paragraph('Where the remaining scatter lives. Within a night it is mostly measurement noise; '
                            'between nights or epochs it is an effect the model misses. Steps at calibration '
                            'epochs would show in the table of epoch medians and in the plot of nightly medians.')
    page += xhtml.table(vrows)
    page += xhtml.image('%s/after_model.png' % outdir, width=1400, height=700)
    page += xhtml.h3('Median residual after the model, per calibration epoch (% per 10 mm)')
    page += xhtml.table(erows)

    # 7. spectrograph balance
    plot_balance(sel, outdir)
    brows = [['arm', 'median sp1 left (%)', 'median sp3 left (%)', 'rms sp1', 'rms sp3', 'r(sp1, focus)', 'r(sp3, focus)']]
    for arm in ARMS:
        s1, s3 = np.asarray(sel[arm + '_dsp1'], float), np.asarray(sel[arm + '_dsp3'], float)
        brows.append([arm, _f(np.median(s1)), _f(np.median(s3)), _f(np.std(s1)), _f(np.std(s3)),
                      _f(np.corrcoef(s1, sel['foc'])[0, 1]), _f(np.corrcoef(s3, sel['foc'])[0, 1])])
    page += xhtml.hline() + xhtml.h2('7. Spectrograph imbalance left')
    page += xhtml.paragraph('The sky-line asymmetry between spectrographs a user sees after the recorded '
                            'correction is applied.')
    page += xhtml.table(brows)
    page += xhtml.image('%s/balance.png' % outdir, width=1400, height=700)

    # 8. continuum check at bright Moon
    bm = (moon_alt > 20) & (fli > 0.8)
    crows = [['arm', 'N bright Moon', 'r(x) continuum vs line', 'r(y) continuum vs line', 'median |cont| - |line|']]
    for arm in ARMS:
        s = sel[bm]
        if len(s) < 10:
            continue
        cg, cd = np.asarray(s[arm + '_cgrad'], float), np.radians(np.asarray(s[arm + '_cdir'], float))
        g = np.isfinite(cg)
        crows.append([arm, g.sum(), _f(np.corrcoef(cg[g] * np.cos(cd[g]), s[arm + '_mx'][g])[0, 1]),
                      _f(np.corrcoef(cg[g] * np.sin(cd[g]), s[arm + '_my'][g])[0, 1]),
                      _f(np.median(cg[g] - s[arm + '_grad'][g]))])
    page += xhtml.hline() + xhtml.h2('8. Continuum check (bright Moon only)')
    page += xhtml.paragraph('Moon above 20 deg and illumination above 0.8, where scattered moonlight '
                            'dominates the continuum. If the continuum plane follows the line plane, the '
                            'effect is a throughput change, not a property of the airglow. Moonlight has a weak '
                            'gradient of its own, so exact agreement is not expected.')
    page += xhtml.table(crows)
    page += xhtml.end()
    with open(outroot + '.html', 'w', encoding='ascii', errors='xmlcharrefreplace') as f:
        f.write(page)
    sel.write(os.path.join(outdir, 'selected.fits'), overwrite=True)

    print('SkyLinePlaneAnalysis: full model (focus + K-mirror + dt), rms before -> after, % per 10 mm:')
    for arm, c, coef, before, after in summary:
        print('  %s %s  c_f %6.2f  c_s %6.2f  c_c %6.2f  c_t %6.2f   %.2f -> %.2f'
              % (arm, c[1], coef['c_f'], coef['c_s'], coef['c_c'], coef['c_t'], before, after))
    print('SkyLinePlaneAnalysis: wrote %s.html' % outroot)


def steer(argv):
    '''
    Parse the command line and run doit
    '''
    table, outroot, minexp, minglat, keep_all = '', 'SkyLinePlaneAnalysis', 300., 15., False
    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-out':
            i += 1
            outroot = argv[i]
        elif argv[i] == '-minexp':
            i += 1
            minexp = float(argv[i])
        elif argv[i] == '-minglat':
            i += 1
            minglat = float(argv[i])
        elif argv[i] == '-all':
            keep_all = True
        elif argv[i].startswith('-'):
            print('SkyLinePlaneAnalysis: unknown option %s' % argv[i])
            print(_USAGE)
            return
        else:
            table = argv[i]
        i += 1
    if not table:
        print(_USAGE)
        return
    doit(table, outroot, minexp, minglat, keep_all)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
