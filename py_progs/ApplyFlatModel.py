#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Apply the large-scale flat-field model (spectrograph factors and a
    science-IFU plane predicted from the telescope focus, the K-mirror
    angle and the date) to DRP 1.3.2 lvmCFrames, optionally followed by a
    direct fit of what is left, and write corrected copies.

Command line usage (if any):

    usage: ApplyFlatModel.py [-h] -model coeffs.ecsv [-outdir DIR] [-calib DIR]
                             [-direct] [-guard X] [-list FILE] cframe ...

    Arguments::

        cframe       DRP 1.3.2 lvmCFrame files (not FixFlatGradient.py
                     .gradfix copies: the model replaces that correction)

    Options::

        -h           print this help and exit
        -model FILE  the model coefficients (flat_model_coeffs.ecsv, from
                     FlatCorr_analysis/flat_model_fit.py); required
        -outdir DIR  where the corrected copies go (default .), written as
                     <name>.flatmodel.fits
        -calib DIR   directory holding <epoch>/lvm-mfiberflat_twilight-<arm>.fits,
                     for the factors 1.3.2 applied (default
                     $LVM_MASTER_DIR, else $SAS_BASE_DIR/sdsswork/lvm/sandbox/calib)
        -direct      after the model, fit the plane and factors left in the
                     sky lines of the corrected exposure and apply them if
                     the guard passes
        -guard X     largest allowed difference, in % per 10 mm, between
                     the residual planes from 5577 (b) and 8399 (z) for the
                     direct step to be applied (default 1.5)
        -list FILE   read CFrame names from FILE (one per line), in
                     addition to those on the command line

Description::

    In 1.3.2 each CFrame was flat-fielded with the stage-1 (twilight,
    large scale removed) master flat times the spectrograph factors
    FACTOR1-3 of its header; the IFU gradient GCOEFF was not applied
    (sdss/lvmdrp#268).  The model gives, for each arm, the large-scale
    throughput relative to the stage-1 flat:

        T = F_s * G,   G = 1 + (gx*xpmm + gy*ypmm)/1000

    with F_s the factor of spectrograph s (normalized to a mean of 1) and
    gx, gy the plane in % per 10 mm, each a linear function of:

        c0_<ff>      a constant per ffactor epoch <ff>
        focus        TESCIFO - focus0, its square, max(0, -focus)
        K            sin and cos of TESCIKM and of 2*TESCIKM
        time         Legendre polynomials of degree 1-3 of the MJD within
                     the ffactor epoch, clipped at the ends of the range
                     the model was fitted over

    For each CFrame:

    1.  the model is evaluated from the header (TESCIFO, TESCIKM, MJD);
        no data from the exposure itself are used
    2.  the factors 1.3.2 applied are read from the master flat of the
        exposure's calibration epoch (the mwave master named in the
        header, else lvmcore calibration-epochs.yaml) and normalized
    3.  per arm the Sci fibers are divided by (F_s/A_s) x G, with G
        normalized to a mean of 1 over the Sci fibers, and IVAR is
        multiplied by its square; the arms are joined with the linear
        ramps across the overlaps of FixFlatGradient.py
    4.  with -direct, the plane and factors left in the corrected
        exposure are fitted to 5577 (b) and 8399 (z) as in
        SkyLinePlane.py; if the two planes agree to within -guard, b and z
        get their own residual and r the mean of the two (the changing
        part is achromatic, and 6364 is affected by [OI] from targets);
        otherwise the model alone is kept
    5.  the model values, the applied factors removed, and the direct
        step (whether applied, its guard value and its planes) are
        recorded in the primary header

Primary routines:

    doit
    apply_model

Notes::

    Only the Sci fibers' FLUX and IVAR change.  The factors 1.3.2 applied
    to the sky and standard-star fibers, the SKY_EAST/SKY_WEST spectra and
    the flux calibration are left as they were: the model is normalized
    to a mean of 1, so the mean level of the Sci fibers, and hence the
    sky subtraction, is kept to the precision of the change of the
    factors (a few 0.1%).

    Exposures without TESCIFO or TESCIKM are skipped; the model cannot be
    evaluated for them.

    Needs numpy, scipy, astropy, pyyaml, FixFlatGradient.py and
    SkyLinePlane.py (in this directory), not lvmdrp.

    The model is a prototype (2026-11-06), fitted to the Utah run of
    SkyLinePlane.py; see Obsidian DRP-Flat-Field-Gradient-2611.md.

History::

    261106 ksl Coding begun

'''

import os
import re
import sys
import numpy as np
from astropy.io import fits
from astropy.table import Table
from numpy.polynomial import legendre

from FixFlatGradient import arm_weights, load_epochs, default_epochs_file
from SkyLinePlane import (plane_fit, line_flux, calib_epoch, LINES, HALF, SIDES, CONT, CONT_HALF, STAR_Q,
                          _num, _cols)

COMP = ['F1', 'F2', 'F3', 'gx', 'gy']


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before the History: line, so -h stays short
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


_USAGE = _usage_from_doc(__doc__)


def read_model(filename):
    '''
    The coefficients as {arm: {component: {term: value}}}, the focus
    reference and the MJD range of each ffactor epoch
    '''
    t = Table.read(filename)
    coef = {}
    for r in t:
        coef.setdefault(str(r['arm']), {}).setdefault(str(r['component']), {})[str(r['term'])] = float(r['value'])
    ranges = {int(k): tuple(v) for k, v in t.meta['ff_ranges'].items()}
    return coef, float(t.meta.get('focus0', 35.5)), ranges


def model_terms(mjd, focus, kmir, focus0, ranges):
    '''
    The values of the terms for one exposure, as {term: value}
    '''
    ffs = sorted(ranges)
    ff = max([e for e in ffs if e <= mjd], default=ffs[0])
    f = focus - focus0
    k = np.radians(kmir)
    terms = {'c0_%d' % ff: 1., 'focus': f, 'focus2': f * f, 'focus_hinge': max(0., -f),
             'sinK': np.sin(k), 'cosK': np.cos(k), 'sin2K': np.sin(2 * k), 'cos2K': np.cos(2 * k)}
    lo, hi = ranges[ff]
    u = float(np.clip(2 * (mjd - lo) / max(1., hi - lo) - 1, -1, 1))
    for deg in (1, 2, 3):
        cc = np.zeros(deg + 1)
        cc[deg] = 1
        terms['t%d_%d' % (deg, ff)] = float(legendre.legval(u, cc))
    return terms, ff


def evaluate(coef, terms):
    '''
    {arm: {component: value}} for one exposure
    '''
    return {arm: {c: sum(v * terms.get(term, 0.) for term, v in coef[arm][c].items()) for c in COMP}
            for arm in coef}


def applied_factors(calib, epoch, arm):
    '''
    FACTOR1-3 of the master flat (normalized to a mean of 1), or None
    '''
    p = os.path.join(calib, str(epoch), 'lvm-mfiberflat_twilight-%s.fits' % arm)
    try:
        h = fits.getheader(p)
        f = np.array([h['%s FIBERFLAT FACTOR%d' % (arm.upper(), k)] for k in (1, 2, 3)], float)
    except (OSError, KeyError):
        return None
    return f / f.mean()


def arm_correction(F, gx, gy, x, y, spec):
    '''
    Per-fiber multiplicative throughput (mean 1 over the given fibers) for
    factors F (fractions, 3) and a plane gx, gy (% per 10 mm)
    '''
    G = 1 + (gx * x + gy * y) / 1000.
    k = F[spec - 1] * G
    return k / np.mean(k)


def residual_fit(flux, wave, x, y, spec):
    '''
    The plane (gx, gy) and factors left in the sky lines of 5577 and 8399
    of corrected Sci-fiber spectra, {arm: dict or None}
    '''
    out = {}
    for arm in 'bz':
        z = line_flux(flux, wave, (LINES[arm] - HALF, LINES[arm] + HALF), SIDES[arm])
        k0, k1 = _cols(wave, CONT[arm] - CONT_HALF, CONT[arm] + CONT_HALF)
        cont = np.nanmedian(flux[:, k0:k1].astype(float), axis=1)
        ok = np.isfinite(cont)
        qlo, qhi = np.percentile(cont[ok], STAR_Q)
        star = ~ok | (cont < qlo) | (cont > qhi)
        r = plane_fit(x, y, np.where(star, np.nan, z), spec)
        if r is not None:
            d = np.radians(r['dir'])
            r['gx'], r['gy'] = r['grad'] * np.cos(d), r['grad'] * np.sin(d)
        out[arm] = r
    return out


def apply_model(filename, coef, focus0, ranges, outdir='.', calib='', epochs=None, direct=False, guard=1.5,
                modelname=''):
    '''
    Write a model-corrected copy of one CFrame; returns the output name or None
    '''
    try:
        x = fits.open(filename)
    except Exception as err:
        print('ApplyFlatModel: could not open %s (%s)' % (filename, err))
        return None
    hdr = x['PRIMARY'].header
    base = os.path.basename(filename)
    if hdr.get('FGRADFIX', False) or hdr.get('FMODEL', False):
        print('ApplyFlatModel: %s already corrected (FGRADFIX or FMODEL), skipped' % base)
        return None
    if not str(hdr.get('DRPVER', '')).startswith('1.3.2'):
        print('ApplyFlatModel: warning: %s has DRPVER %s; the model was made for 1.3.2' % (base, hdr.get('DRPVER')))
    mjd, focus, kmir = _num(hdr, 'MJD'), _num(hdr, 'TESCIFO'), _num(hdr, 'TESCIKM')
    if not (np.isfinite(mjd) and np.isfinite(focus) and np.isfinite(kmir) and focus > -900 and kmir > -900):
        print('ApplyFlatModel: %s lacks MJD, TESCIFO or TESCIKM, skipped' % base)
        return None
    epoch = calib_epoch(hdr, epochs)
    if epoch < 0:
        print('ApplyFlatModel: no calibration epoch for %s, skipped' % base)
        return None
    app = {arm: applied_factors(calib, epoch, arm) for arm in 'brz'}
    if any(a is None for a in app.values()):
        print('ApplyFlatModel: applied factors of epoch %d not found in %s, skipped' % (epoch, calib))
        return None

    terms, ff = model_terms(mjd, focus, kmir, focus0, ranges)
    val = evaluate(coef, terms)

    st = x['SLITMAP'].data
    sci = np.where(np.char.strip(np.asarray(st['telescope']).astype(str)) == 'Sci')[0]
    xp, yp = np.asarray(st['xpmm'], float)[sci], np.asarray(st['ypmm'], float)[sci]
    spec = np.asarray(st['spectrographid'], int)[sci]
    wave = np.asarray(x['WAVE'].data, float)
    w = arm_weights(wave)

    kmod = {}
    for arm in 'brz':
        F = 1 + np.array([val[arm]['F1'], val[arm]['F2'], val[arm]['F3']]) / 100.
        F = F / F.mean()
        kmod[arm] = arm_correction(F / app[arm], val[arm]['gx'], val[arm]['gy'], xp, yp, spec)
    corr = sum(kmod[arm][:, None] * w[arm][None, :] for arm in 'brz')

    flux, ivar = x['FLUX'].data, x['IVAR'].data
    fs = (flux[sci] / corr).astype(np.float32)

    used_direct, gval, res = False, np.nan, None
    if direct:
        good = np.isfinite(fs)
        fs_clean = np.where(good, fs, np.nan)
        res = residual_fit(fs_clean, wave, xp, yp, spec)
        if res['b'] is not None and res['z'] is not None:
            gval = float(np.hypot(res['b']['gx'] - res['z']['gx'], res['b']['gy'] - res['z']['gy']))
            nmin = 0.85 * np.sum(np.isfinite(xp))
            if gval <= guard and min(res['b']['n'], res['z']['n']) >= 0.7 * nmin:
                used_direct = True
                kres = {}
                for arm in 'bz':
                    r = res[arm]
                    kres[arm] = arm_correction(np.array([r['f1'], r['f2'], r['f3']]), r['gx'], r['gy'], xp, yp, spec)
                rm = {k: 0.5 * (res['b'][k] + res['z'][k]) for k in ('f1', 'f2', 'f3', 'gx', 'gy')}
                kres['r'] = arm_correction(np.array([rm['f1'], rm['f2'], rm['f3']]), rm['gx'], rm['gy'], xp, yp, spec)
                cres = sum(kres[arm][:, None] * w[arm][None, :] for arm in 'brz')
                corr = corr * cres
                fs = (flux[sci] / corr).astype(np.float32)

    flux[sci] = fs.astype(flux.dtype)
    ivar[sci] = (ivar[sci] * corr ** 2).astype(ivar.dtype)

    hdr['FMODEL'] = (True, 'large-scale flat model applied to Sci fibers')
    hdr['FMODFILE'] = (os.path.basename(modelname)[:60], 'flat model coefficients')
    hdr['FMODFFE'] = (ff, 'ffactor epoch of the model terms')
    hdr['FMODCAL'] = (epoch, 'calibration epoch of the 1.3.2 factors removed')
    for arm in 'brz':
        A = arm.upper()
        for c in COMP:
            hdr['FM%s%s' % (A, c.upper())] = (round(val[arm][c], 4), '%s model %s (%s)' % (arm, c, '%' if c[0] == 'F' else '%/10mm'))
        for k in range(3):
            hdr['FM%sA%d' % (A, k + 1)] = (round(float(app[arm][k]), 5), '%s 1.3.2 factor %d removed (normalized)' % (arm, k + 1))
    hdr['FMDIRECT'] = (used_direct, 'direct residual fit applied')
    if direct:
        hdr['FMDGUARD'] = (round(gval, 3) if np.isfinite(gval) else -1., '|5577 - 8399| residual plane (%/10mm)')
        for arm in 'bz':
            r = res[arm] if res else None
            hdr['FMD%sGX' % arm.upper()] = (round(r['gx'], 3) if r else -999., '%s residual plane x (%%/10mm)' % arm)
            hdr['FMD%sGY' % arm.upper()] = (round(r['gy'], 3) if r else -999., '%s residual plane y (%%/10mm)' % arm)
    hdr.add_history('ApplyFlatModel: Sci FLUX/IVAR corrected with the large-scale flat model (ffactor epoch %d)%s'
                    % (ff, '; direct residual fit applied' if used_direct else ''))

    root = base[:-5] if base.endswith('.fits') else base
    outname = os.path.join(outdir, root + '.flatmodel.fits')
    x.writeto(outname, overwrite=True)
    x.close()
    print('ApplyFlatModel: %s -> %s  ff %d cal %d  model planes b %.1f,%.1f r %.1f,%.1f z %.1f,%.1f  direct %s%s'
          % (base, outname, ff, epoch, val['b']['gx'], val['b']['gy'], val['r']['gx'], val['r']['gy'],
             val['z']['gx'], val['z']['gy'], used_direct, (' (guard %.2f)' % gval) if direct else ''))
    return outname



def doit(files, model, outdir='.', calib='', direct=False, guard=1.5):
    '''
    Description:
        Apply the model to each CFrame in files

    History:
        261106 ksl Coding begun
    '''
    coef, focus0, ranges = read_model(model)
    if not calib:
        calib = os.environ.get('LVM_MASTER_DIR', '') or \
            os.path.join(os.environ.get('SAS_BASE_DIR', ''), 'sdsswork/lvm/sandbox/calib')
    try:
        epochs = load_epochs(default_epochs_file())
    except Exception:
        epochs = None
    os.makedirs(outdir, exist_ok=True)
    n = 0
    for fn in files:
        if apply_model(fn, coef, focus0, ranges, outdir, calib, epochs, direct, guard, model):
            n += 1
    print('ApplyFlatModel: corrected %d of %d files' % (n, len(files)))


def steer(argv):
    '''
    Parse the command line and run doit
    '''
    files, model, outdir, calib, direct, guard, listfile = [], '', '.', '', False, 1.5, ''
    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-model':
            i += 1
            model = argv[i]
        elif argv[i] == '-outdir':
            i += 1
            outdir = argv[i]
        elif argv[i] == '-calib':
            i += 1
            calib = argv[i]
        elif argv[i] == '-direct':
            direct = True
        elif argv[i] == '-guard':
            i += 1
            guard = float(argv[i])
        elif argv[i] == '-list':
            i += 1
            listfile = argv[i]
        elif argv[i].startswith('-'):
            print('ApplyFlatModel: unknown option %s' % argv[i])
            print(_USAGE)
            return
        else:
            files.append(argv[i])
        i += 1
    if listfile:
        with open(listfile) as f:
            files += [ln.split()[0] for ln in f if ln.strip() and not ln.startswith('#')]
    if not files or not model:
        print(_USAGE)
        return
    doit(files, model, outdir, calib, direct, guard)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
