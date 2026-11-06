#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Fit the large-scale flat-field model (spectrograph factors and a
    science-IFU plane per arm, as functions of the telescope focus, the
    K-mirror angle and the date) to a SkyLinePlane.py table, report how
    much each stage of the model improves the flat field, and write the
    coefficients for ApplyFlatModel.py.

Command line usage (if any):

    usage: FitFlatModel.py [-h] [-out coeffs.ecsv] [-calib DIR] [-minexp S] [-drpall FILE] table.txt

    Arguments::

        table.txt    a table written by SkyLinePlane.py from DRP 1.3.2
                     lvmCFrames (not .gradfix copies)

    Options::

        -h           print this help and exit
        -out FILE    the coefficients (default flat_model_coeffs.ecsv); the
                     report goes to the same name with .txt
        -calib DIR   directory holding <epoch>/lvm-mfiberflat_twilight-<arm>.fits,
                     for the factors 1.3.2 applied and the fiber positions
                     (default $LVM_MASTER_DIR, else
                     $SAS_BASE_DIR/sdsswork/lvm/sandbox/calib)
        -minexp S    shortest exposure time used (default 300 s)
        -drpall FILE a drpall file with obstime, sci_ra, sci_dec and sci_pa,
                     to test an airglow term (diagnostic only, not part of
                     the model)

Description::

    The target is the large-scale throughput of the science IFU relative
    to the stage-1 (twilight, large scale removed) master flat, per arm
    and exposure:

        T = F_s * (1 + (gx*xpmm + gy*ypmm)/1000)

    In a 1.3.2 CFrame the flat used was the stage-1 flat x the applied
    FACTOR1-3 of the master-flat header (no gradient, sdss/lvmdrp#268),
    so F_s = applied factor x measured residual factor (normalized to a
    mean of 1) and gx, gy = the measured plane (% per 10 mm).  Five
    components per arm are modelled: F1-1, F2-1, F3-1 (%), gx, gy.

    Stages, each tested by 5-fold cross-validation over whole nights (the
    model never sees the night it predicts):

        0   stage-1 flat only (no large-scale correction)
        A0  1.3.2 as released: the applied factors
        A   1.3.2 as intended: the applied factors + the GCOEFF plane
        B   one constant per ffactor epoch
        C   B + focus (TESCIFO - 35.5, its square, max(0, 35.5 - TESCIFO))
            + K-mirror (sin, cos of K and 2K)
        D   C + Legendre time terms (degree 1-3) within each ffactor epoch;
            this is the model written out
        D'  D with one constant per calibration epoch
        E   D + a nightly offset from the other exposures of the night
        F   a direct fit in the exposure itself (measurement noise only)

    The error of each stage is the rms over the Sci fibers of the
    large-scale pattern error (%), using the fiber positions and
    spectrographs, with the overall normalization removed; "true" removes
    the measurement noise (the covariance of the five components from
    consecutive exposures with almost the same focus and K-mirror angle)
    in quadrature, and "median" is the median per-exposure error.

    Exposures are used if DRPQUAL = 0, at least -minexp long, with focus
    and K-mirror angle in the header, and with planes (relative to the
    GCOEFF prediction) within 6 MAD of the median in every arm (the
    selection of SkyLinePlaneAnalysis.py); per arm, those with fewer than
    85% of the usual number of fitted fibers (part of a spectrograph
    missing) are left out.  These outliers, about 1%, are what the
    per-exposure guard of ApplyFlatModel.py is for.  All latitudes are used: the b and z sky lines show no
    target contamination near the Galactic plane.

    The coefficients of model D fitted to all the data are written as an
    ecsv table (arm, component, term, value), whose metadata hold what is
    needed to evaluate the model: focus0, the MJD range of each ffactor
    epoch and a description of the terms.

Primary routines:

    doit
    fit_arm

Notes::

    Needs numpy, astropy, pyyaml and SkyLinePlaneAnalysis.py (in this
    directory), not lvmdrp.

    The fiber positions are read from the SLITMAP of a master flat in
    -calib.

    Typical results (2026-11-06, 21,800 exposures): fiber rms of the
    large-scale error, b / r / z, 1.3.2 as released 3.3 / 3.5 / 2.5%,
    1.3.2 as intended 1.0 / 1.1 / 0.9%, model D 0.47 / 0.59 / 0.45%
    (median 0.24 / 0.34 / 0.24%), direct fit 0.11 / 0.26 / 0.13%.

History::

    261106 ksl Coding begun, from the prototype FlatCorr_analysis/flat_model_fit.py

'''

import os
import re
import sys
import numpy as np
from numpy.polynomial import legendre
from astropy.io import fits
from astropy.table import Table, join

from SkyLinePlaneAnalysis import to_numeric, prepare

COMP = ['F1', 'F2', 'F3', 'gx', 'gy']
FOCUS0 = 35.5


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before the History: line, so -h stays short
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


_USAGE = _usage_from_doc(__doc__)


def fiber_basis(calib, epoch):
    '''
    The 5x5 matrix C such that the rms (%) over the Sci fibers of the
    pattern of a component error vector e is sqrt(e C e)
    '''
    p = os.path.join(calib, str(epoch), 'lvm-mfiberflat_twilight-b.fits')
    sm = fits.getdata(p, 'SLITMAP')
    m = np.char.strip(np.asarray(sm['telescope']).astype(str)) == 'Sci'
    x = np.asarray(sm['xpmm'], float)[m]
    y = np.asarray(sm['ypmm'], float)[m]
    s = np.asarray(sm['spectrographid'], int)[m]
    B = np.vstack([(s == 1).astype(float), (s == 2).astype(float), (s == 3).astype(float), x / 10, y / 10]).T
    B = B - B.mean(axis=0)
    return B.T @ B / len(B)


def applied_factors(calib, epochs, arm):
    '''
    {epoch: FACTOR1-3 of the master flat, normalized to a mean of 1}
    '''
    out = {}
    for e in epochs:
        p = os.path.join(calib, str(e), 'lvm-mfiberflat_twilight-%s.fits' % arm)
        try:
            h = fits.getheader(p)
            f = np.array([h['%s FIBERFLAT FACTOR%d' % (arm.upper(), k)] for k in (1, 2, 3)], float)
            out[e] = f / f.mean()
        except (OSError, KeyError):
            out[e] = np.full(3, np.nan)
    return out


def targets(t, arm, fac):
    '''
    The 5 components per exposure (F in %, plane in % per 10 mm), and the
    predictions of 1.3.2 as released (A0) and as intended (A)
    '''
    ep = np.asarray(t['epoch'])
    app = np.array([fac.get(e, np.full(3, np.nan)) for e in ep])
    meas = np.vstack([np.asarray(t['%s_f%d' % (arm, k)], float) for k in (1, 2, 3)]).T
    F = app * meas
    F = F / F.mean(axis=1)[:, None]
    V = np.column_stack([100 * (F - 1), np.asarray(t[arm + '_mx'], float), np.asarray(t[arm + '_my'], float)])
    A0 = np.column_stack([100 * (app - 1), np.zeros(len(t)), np.zeros(len(t))])
    A = np.column_stack([100 * (app - 1), np.asarray(t[arm + '_px'], float), np.asarray(t[arm + '_py'], float)])
    return V, A0, A


def design(t, level, extra=None, groups='ff'):
    '''
    Design matrix and term names of the nested models (level 1 = B,
    3 = C, 4 = D), with constants per ffactor epoch (groups='ff') or per
    calibration epoch (groups='epoch')
    '''
    ff = np.asarray(t['ff'])
    ffs = sorted(set(ff))
    cols, names = [], []
    grp = np.asarray(t[groups])
    for e in sorted(set(grp)):
        cols.append((grp == e).astype(float))
        names.append('c0_%d' % e)
    if level >= 3:
        f = np.asarray(t['foc'], float) - FOCUS0
        k = np.radians(np.asarray(t['kmir'], float))
        cols += [f, f * f, np.maximum(0, -f), np.sin(k), np.cos(k), np.sin(2 * k), np.cos(2 * k)]
        names += ['focus', 'focus2', 'focus_hinge', 'sinK', 'cosK', 'sin2K', 'cos2K']
    if level >= 4:
        m = np.asarray(t['mjd'], float)
        for e in ffs:
            sel = ff == e
            lo, hi = m[sel].min(), m[sel].max()
            u = np.where(sel, 2 * (m - lo) / max(1., hi - lo) - 1, 0.)
            for deg in (1, 2, 3):
                cc = np.zeros(deg + 1)
                cc[deg] = 1
                cols.append(np.where(sel, legendre.legval(u, cc), 0.))
                names.append('t%d_%d' % (deg, e))
    if extra is not None:
        for n, c in extra:
            cols.append(c)
            names.append(n)
    return np.vstack(cols).T, names


def cv(t, V, A, folds=5, seed=1):
    '''
    Cross-validated predictions, the folds made of whole nights
    '''
    nights = np.asarray(t['night'])
    un = np.unique(nights)
    fold_of = dict(zip(un, np.random.default_rng(seed).integers(0, folds, len(un))))
    fold = np.array([fold_of[n] for n in nights])
    P = np.full(V.shape, np.nan)
    for k in range(folds):
        tr, te = fold != k, fold == k
        for j in range(V.shape[1]):
            g = tr & np.isfinite(V[:, j])
            c, *_ = np.linalg.lstsq(A[g], V[g, j], rcond=None)
            P[te, j] = A[te] @ c
    return P


def nightly(t, V, P):
    '''
    P plus the mean residual of the other exposures of the same night
    '''
    nights = np.asarray(t['night'])
    out = np.full(V.shape, np.nan)
    R = V - P
    for n in np.unique(nights):
        idx = np.where(nights == n)[0]
        if len(idx) < 2:
            continue
        s = np.nansum(R[idx], axis=0)
        out[idx] = P[idx] + (s[None, :] - R[idx]) / (len(idx) - 1)
    return out


def pattern_rms(E, C, median=False):
    '''
    rms over exposures (or the median) of the per-exposure fiber-pattern
    rms (%) of the error vectors E
    '''
    g = np.all(np.isfinite(E), axis=1)
    per = np.einsum('ij,jk,ik->i', E[g], C, E[g])
    return np.sqrt(np.median(per) if median else np.mean(per))


def noise_cov(t, V):
    '''
    Covariance of the measurement noise of the 5 components, from
    consecutive exposures of a night with almost the same focus and
    K-mirror angle; and the number of pairs
    '''
    o = np.argsort(np.asarray(t['mjd'] * 1e6 + t['expnum']))
    n, e = np.asarray(t['night'])[o], np.asarray(t['expnum'])[o]
    f, k = np.asarray(t['foc'], float)[o], np.asarray(t['kmir'], float)[o]
    ok = (n[1:] == n[:-1]) & (np.abs(f[1:] - f[:-1]) < 0.03) & (np.abs(k[1:] - k[:-1]) < 5) & (e[1:] - e[:-1] < 5)
    D = (V[o][1:] - V[o][:-1])[ok]
    D = D[np.all(np.isfinite(D), axis=1)]
    mad = 1.4826 * np.median(np.abs(D - np.median(D, axis=0)), axis=0)
    D = D[np.all(np.abs(D - np.median(D, axis=0)) < 5 * mad, axis=1)]
    return np.cov(D.T) / 2, ok.sum()


def airglow_terms(t, drpall):
    '''
    The exposures matched to a drpall file, and the zenith-direction x
    cot(alt) terms (an airglow gradient) for them
    '''
    from astropy.time import Time
    from astropy.coordinates import EarthLocation, Longitude
    import astropy.units as u
    d = Table.read(drpall, hdu=1)['expnum', 'obstime', 'sci_ra', 'sci_dec', 'sci_pa']
    for c in d.colnames:
        if d[c].dtype.kind == 'S':
            d[c] = np.char.decode(d[c]).astype(str)
    t = join(t, d, keys='expnum')
    t = t[np.asarray(t['sci_ra'], float) > -900]
    lco = EarthLocation(lat=-29.0146 * u.deg, lon=-70.6926 * u.deg, height=2380 * u.m)
    tm = Time(list(t['obstime']), format='isot', scale='utc')
    ha = (tm.sidereal_time('apparent', longitude=lco.lon) -
          Longitude(np.asarray(t['sci_ra'], float) * u.deg)).wrap_at(180 * u.deg).rad
    dec, phi = np.radians(np.asarray(t['sci_dec'], float)), lco.lat.rad
    q = np.arctan2(np.sin(ha), np.tan(phi) * np.cos(dec) - np.sin(dec) * np.cos(ha))
    alt = np.arcsin(np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.cos(ha))
    pa = np.radians(np.asarray(t['sci_pa'], float))
    cot = 1 / np.tan(alt)
    extra = [('Zc_s-', cot * np.sin(q - pa)), ('Zc_c-', cot * np.cos(q - pa)),
             ('Zc_s+', cot * np.sin(q + pa)), ('Zc_c+', cot * np.cos(q + pa))]
    return t, extra


def fit_arm(t, arm, calib, C, drpall, say):
    '''
    Report the stages for one arm; returns the coefficient rows of model D
    '''
    epochs = sorted(set(np.asarray(t['epoch'])))
    fac = applied_factors(calib, epochs, arm)
    V, A0p, Ap = targets(t, arm, fac)
    nfib = np.asarray(t[arm + '_n'], float)
    broken = nfib < 0.85 * np.median(nfib)
    good = np.all(np.isfinite(V), axis=1) & ~broken
    tt, V, A0p, Ap = t[good], V[good], A0p[good], Ap[good]
    Sig, npair = noise_cov(tt, V)
    noise = np.sqrt(np.trace(C @ Sig))
    preds = {'0  stage-1 flat only': np.zeros_like(V), 'A0 1.3.2 as released (factors)': A0p,
             'A  1.3.2 as intended (factors + GCOEFF)': Ap}
    for lab, lev in (('B  constant per ffactor epoch', 1), ('C  B + focus + K-mirror', 3),
                     ('D  C + smooth time (the model)', 4)):
        preds[lab] = cv(tt, V, design(tt, lev)[0])
    preds["D' D, constants per calibration epoch"] = cv(tt, V, design(tt, 4, groups='epoch')[0])
    preds['E  D + nightly offset'] = nightly(tt, V, preds['D  C + smooth time (the model)'])
    say()
    say('== %s arm: N = %d (%d exposures with < 85%% of the usual fibers left out); noise from %d pairs: %.2f%% '
        'fiber rms (components: %s)' % (arm, len(tt), broken.sum(), npair, noise,
                                        ', '.join('%s %.2f' % (c, np.sqrt(Sig[i, i])) for i, c in enumerate(COMP))))
    say('   %-42s %8s %8s %8s | %12s %12s' % ('stage', 'error', 'true', 'median', 'true |b|<15', 'true |b|>=15'))
    gl = np.abs(np.asarray(tt['glat'], float))

    def tr(v):
        return np.sqrt(max(0, v * v - noise * noise))

    for lab, P in preds.items():
        E = P - V
        say('   %-42s %8.2f %8.2f %8.2f | %12.2f %12.2f'
            % (lab, pattern_rms(E, C), tr(pattern_rms(E, C)), pattern_rms(E, C, median=True),
               tr(pattern_rms(E[gl < 15], C)), tr(pattern_rms(E[gl >= 15], C))))
    say('   %-42s %8s %8.2f %8s |' % ('F  direct fit (noise)', '', noise, ''))
    E = preds['D  C + smooth time (the model)'] - V
    say('   model D, per component: rms of data -> CV residual: ' +
        ', '.join('%s %.2f->%.2f' % (c, np.nanstd(V[:, i]), np.nanstd(E[:, i])) for i, c in enumerate(COMP)))
    if drpall:
        ta, extra = airglow_terms(tt, drpall)
        Va, _, _ = targets(ta, arm, fac)
        g2 = np.all(np.isfinite(Va), axis=1)
        ta, Va = ta[g2], Va[g2]
        extra = [(n_, c[g2]) for n_, c in extra]
        Az, _ = design(ta, 4, extra)
        amp = []
        for j in range(len(COMP)):
            co, *_ = np.linalg.lstsq(Az, Va[:, j], rcond=None)
            amp.append(np.sqrt(np.sum(co[-4:] ** 2)))
        say('   airglow term (matched to %s, N = %d): CV error D %.2f -> D + airglow %.2f; amplitude per unit '
            'cot(alt): %s' % (os.path.basename(drpall), len(ta), pattern_rms(cv(ta, Va, design(ta, 4)[0]) - Va, C),
                              pattern_rms(cv(ta, Va, Az) - Va, C), ', '.join('%s %.2f' % (c, a) for c, a in zip(COMP, amp))))
    A, names = design(tt, 4)
    rows = []
    for j, c in enumerate(COMP):
        co, *_ = np.linalg.lstsq(A, V[:, j], rcond=None)
        rows += [(arm, c, n_, float(v)) for n_, v in zip(names, co)]
    return rows, tt


def doit(table, out='flat_model_coeffs.ecsv', calib='', minexp=300., drpall=''):
    '''
    Description:
        Fit the model to a SkyLinePlane.py table, print and write the
        report, and write the coefficients

    History:
        261106 ksl Coding begun
    '''
    if not calib:
        calib = os.environ.get('LVM_MASTER_DIR', '') or \
            os.path.join(os.environ.get('SAS_BASE_DIR', ''), 'sdsswork/lvm/sandbox/calib')
    lines = []

    def say(s=''):
        print(s)
        lines.append(s)

    t = Table.read(table, format='ascii.fixed_width_two_line')
    t = to_numeric(t)
    for c in ('focus', 'kmir', 'drpall_sci_focpos', 'drpall_sci_kmpos'):
        if c in t.colnames and t[c].dtype.kind == 'f':
            t[c][np.asarray(t[c]) <= -900] = np.nan
    t, _, _ = prepare(t, minexp, 0., False, clip=True)
    t = t[np.isfinite(np.asarray(t['foc'], float)) & np.isfinite(np.asarray(t['kmir'], float))]
    if len(t) < 100:
        print('FitFlatModel: only %d usable exposures, stopping' % len(t))
        return
    C = fiber_basis(calib, sorted(set(np.asarray(t['epoch'])))[0])
    say('Large-scale flat model fitted to %s: %d exposures, %d nights, all latitudes'
        % (os.path.basename(table), len(t), len(set(t['night']))))
    say('Error = rms over the Sci fibers of the large-scale pattern error (%), relative to the stage-1 flat;')
    say('"true" removes the measurement noise (from consecutive exposures) in quadrature.')
    rows, ranges = [], {}
    for arm in 'brz':
        r, tt = fit_arm(t, arm, calib, C, drpall, say)
        rows += r
        if arm == 'b':
            m, ff = np.asarray(tt['mjd'], float), np.asarray(tt['ff'])
            ranges = {int(e): [float(m[ff == e].min()), float(m[ff == e].max())] for e in sorted(set(ff))}
    ct = Table(rows=rows, names=('arm', 'component', 'term', 'value'))
    ct.meta['description'] = ('Large-scale flat model: per arm, F1-F3 in % (spectrograph factors relative to the '
                              'stage-1 twilight flat, mean 1) and the Sci-IFU plane gx, gy in % per 10 mm '
                              '(G = 1 + (gx*xpmm + gy*ypmm)/1000); fitted to ' + os.path.basename(table))
    ct.meta['focus0'] = FOCUS0
    ct.meta['ff_ranges'] = ranges
    ct.meta['terms'] = ('c0_<ff>: 1 in ffactor epoch <ff>; focus = TESCIFO - focus0; focus2 = focus^2; '
                        'focus_hinge = max(0, -focus); sinK, cosK, sin2K, cos2K of TESCIKM; t<d>_<ff> = Legendre '
                        'P_d(u) in epoch <ff>, u = 2 (MJD - lo)/(hi - lo) - 1 with [lo, hi] = ff_ranges[<ff>]')
    ct.meta['created'] = 'FitFlatModel.py'
    ct.write(out, overwrite=True)
    say()
    say('FitFlatModel: model D coefficients (%d) written to %s' % (len(ct), out))
    with open(os.path.splitext(out)[0] + '.txt', 'w') as f:
        f.write('\n'.join(lines) + '\n')


def steer(argv):
    '''
    Parse the command line and run doit
    '''
    table, out, calib, minexp, drpall = '', 'flat_model_coeffs.ecsv', '', 300., ''
    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-out':
            i += 1
            out = argv[i]
        elif argv[i] == '-calib':
            i += 1
            calib = argv[i]
        elif argv[i] == '-minexp':
            i += 1
            minexp = float(argv[i])
        elif argv[i] == '-drpall':
            i += 1
            drpall = argv[i]
        elif argv[i].startswith('-'):
            print('FitFlatModel: unknown option %s' % argv[i])
            print(_USAGE)
            return
        else:
            table = argv[i]
        i += 1
    if not table:
        print(_USAGE)
        return
    doit(table, out, calib, minexp, drpall)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
