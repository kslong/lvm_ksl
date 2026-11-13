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

    usage: ApplyFlatModel.py [-h] [-model coeffs.ecsv] [-outdir DIR | -outroot DIR]
                             [-calib DIR] [-direct] [-guard X] [-rguard X]
                             [-list FILE] [-drpall FILE] [-base DIR] [-mjd LO-HI]
                             [-np N] [-part I/N] [-redo] [-summary FILE] [cframe ...]
           ApplyFlatModel.py -skysub [-drpall FILE | cframe ...] [-outroot DIR] [-mjd LO-HI]
                             [-np N] [-part I/N] [-redo] [-summary FILE]
           ApplyFlatModel.py -report summary.txt [summary.txt ...]

    Arguments::

        cframe       DRP 1.3.2 lvmCFrame files (not FixFlatGradient.py
                     .gradfix copies: the model replaces that correction)

    Options::

        -h           print this help and exit
        -model FILE  the model coefficients (an ecsv from FitFlatModel.py);
                     a name not found as given is looked for in the lvm_ksl
                     data directory (../data from this script).  Default:
                     the 1.3.2 model, data/flat_model_coeffs_1.3.2.ecsv
        -outdir DIR  where the corrected copies go (default .), written as
                     <name>.flatmodel.fits
        -outroot DIR write the corrected copies under DIR with the same
                     <tilegrp>/<tileid>/<mjd>/lvmCFrame-<expnum>.fits path
                     and name as the originals have under their redux
                     version directory; default with -drpall:
                     $SAS_BASE_DIR/sdsswork/lvm/spectro/redux/knox.  An
                     original is never overwritten
        -drpall FILE take the exposures from a drpall file (fits or h5):
                     the lvmCFrame next to each lvmSFrame in its location
                     column (relative to $SAS_BASE_DIR, or to -base)
        -base DIR    resolve the drpall locations under DIR instead of
                     $SAS_BASE_DIR
        -mjd LO-HI   with -drpall, only the exposures with LO <= MJD <= HI
        -np N        correct N files at a time (default 1)
        -part I/N    correct only part I of N (I = 1..N) of the files:
                     every Nth file of the sorted list, starting at the
                     Ith, so N jobs (e.g. one per node) share a run.
                     The summary default becomes
                     ApplyFlatModel_summary_I_of_N.txt
        -redo        with -outroot, redo exposures whose output exists
                     (default: skip them, so a run can be restarted)
        -summary FILE one row per exposure: status, output, model planes,
                     direct step and the quality planes (default
                     ApplyFlatModel_summary.txt)
        -report      only print the quality report of the summary files
                     given as arguments (e.g. those of all the -part jobs)
        -calib DIR   directory holding <epoch>/lvm-mfiberflat_twilight-<arm>.fits,
                     for the factors 1.3.2 applied (default
                     $LVM_MASTER_DIR, else $SAS_BASE_DIR/sdsswork/lvm/sandbox/calib)
        -direct      after the model, fit the plane and factors left in the
                     sky lines of the corrected exposure and apply them if
                     the guard passes
        -guard X     largest allowed difference, in % per 10 mm, between
                     the residual planes from 5577 (b) and 8399 (z) for the
                     direct step to be applied (default 1.5)
        -rguard X    largest allowed difference, in % per 10 mm, between
                     the r OH residual plane and the mean of b and z for
                     the r OH lines to be used in r (default 1.5)
        -list FILE   read CFrame names from FILE (one per line), in
                     addition to those on the command line
        -skysub      instead of correcting, redo the sky subtraction of
                     corrected CFrames already written: the same files
                     and options as the correction run (-drpall, -outroot,
                     -mjd, -part, -np, -redo) select the corrected CFrames
                     (<outroot>/.../lvmCFrame-<expnum>.fits, or the
                     CFrames given when there is no -outroot), and an
                     lvmSFrame-<expnum>.fits is written next to each,
                     with lvmdrp's quick_sky_subtraction.  Needs lvmdrp
                     (any version from 1.3.2).  The summary default is
                     ApplyFlatModel_skysub_summary[_I_of_N].txt.  About
                     40 s and 3.7 GB of memory per file, so choose -np
                     by the memory of the node (e.g. 128 GB -> -np 30)

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
        exposure are fitted to 5577 (b), 8399 (z) and the sum of the r OH
        lines 6864, 6923, 7316 and 7341 (r), as in SkyLinePlane.py.  If
        the b and z planes agree to within -guard, b and z get their own
        residual.  r gets the average of its r OH residual and the mean of
        the b and z residuals, if the r OH plane agrees with that mean to
        within -rguard; otherwise r gets the mean of b and z alone.  6364
        is not used: in the v2 held-out tests it was no better than the
        model, while r OH and the b+z mean each left about 0.3-0.4 %/10 mm
        and their average less.  If the b and z guard fails, the model
        alone is kept in all arms
    5.  the model values, the applied factors removed, and the direct
        step (whether applied, its guard value and its planes) are
        recorded in the primary header, with the original file and the
        lvm_ksl commit
    6.  quality: the planes of 5577 and of the summed held-out r OH and
        z OH lines are measured before and after the correction, as
        SkyLinePlane.py measures them (Sci fibers with fibstatus 0,
        stars rejected by their continuum), and recorded as their
        amplitudes in %/10 mm (FMQB0/FMQB1 for 5577 before/after,
        FMQR0/FMQR1 for r OH, FMQZ0/FMQZ1 for z OH) and in the summary.
        The held-out lines are never fitted by the direct step, so r OH
        and z OH after the correction are the check; 5577 is fitted by
        it.  A quality report (medians before and after, fraction below
        1 %/10 mm, direct step rejected, worst exposures) ends each run,
        and -report gives it for any set of summary files

    Sky subtraction (-skysub): each corrected CFrame is given to lvmdrp's
    quick_sky_subtraction, as in the DRP.  lvmdrp writes the sky table to
    the ancillary directory of the reduction version it runs as, so each
    worker runs as the version of the tree the CFrames are in (for
    redux/knox, LVMDRP_VERSION=knox and LVM_SPECTRO_REDUX=redux) and
    creates the ancillary and qa directories, which 1.3.2 does not.  The
    drpall is not updated.

    For a large run, e.g. of all of 1.3.2 at Utah:

        ApplyFlatModel.py -direct -drpall drpall-1.3.2.fits -np 8

    reads the 1.3.2 CFrames where the drpall says they are and writes
    the corrected ones to redux/knox/<tilegrp>/<tileid>/<mjd>/.  To
    spread it over several nodes, run the same command in N jobs with
    -part 1/N ... -part N/N (e.g. a Slurm job array with
    -part ${SLURM_ARRAY_TASK_ID}/N); each writes its own summary.

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
    261107 ksl r direct step from the r OH lines averaged with the b+z
               mean, instead of the b+z mean alone; -rguard
    261113 ksl -skysub: the sky subtraction of the corrected CFrames
    261112 ksl -drpall, -outroot (same tree layout as the originals,
               default redux/knox), -mjd, -np, -part, -redo, -summary;
               original file and lvm_ksl commit in the header; -model
               defaults to the 1.3.2 model in the lvm_ksl data directory;
               quality planes before/after in the header and summary,
               a quality report at the end of a run, and -report

'''

import multiprocessing
import os
import re
import subprocess
import sys
import numpy as np
from astropy.io import fits
from astropy.table import Table
from numpy.polynomial import legendre

from FixFlatGradient import arm_weights, load_epochs, default_epochs_file
from SkyLinePlane import (plane_fit, line_flux, calib_epoch, LINES, HALF, SIDES, CONT, CONT_HALF, STAR_Q,
                          HOLDOUT, OH_GROUPS, _num, _cols, read_drpall)

COMP = ['F1', 'F2', 'F3', 'gx', 'gy']
DEFAULT_OUTROOT = 'sdsswork/lvm/spectro/redux/knox'   # under $SAS_BASE_DIR
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
DEFAULT_MODEL = 'flat_model_coeffs_1.3.2.ecsv'


def find_model(name=''):
    '''
    The model file: name as given if it exists, else name in the lvm_ksl
    data directory; with no name, the 1.3.2 model there.  '' if not found
    '''
    name = name or DEFAULT_MODEL
    for p in (name, os.path.join(DATA_DIR, os.path.basename(name))):
        if os.path.isfile(p):
            return os.path.abspath(p)
    return ''


def code_version():
    '''
    The lvm_ksl commit of this script (with -dirty if modified), or unknown
    '''
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        h = subprocess.run(['git', '-C', here, 'rev-parse', '--short=10', 'HEAD'],
                           capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(['git', '-C', here, 'diff', '--quiet', 'HEAD', '--', os.path.basename(__file__)]).returncode
        return h + ('-dirty' if dirty else '')
    except Exception:
        return 'unknown'


def tree_path(filename, outroot):
    '''
    The path of filename under outroot with the same last four parts,
    <tilegrp>/<tileid>/<mjd>/<name>, as the original under its redux
    version directory
    '''
    parts = os.path.normpath(os.path.abspath(filename)).split(os.sep)
    return os.path.join(outroot, *parts[-4:])


QUALITY = ('5577', 'r_oh', 'z_oh')
QKEY = {'5577': 'B', 'r_oh': 'R', 'z_oh': 'Z'}


def quality(flux, wave, x, y, spec):
    '''
    The planes (gx, gy) in %/10 mm of 5577 and of the summed held-out r OH
    and z OH lines in the Sci fibers given (fibers x wavelength), measured
    as SkyLinePlane.measure does; (nan, nan) where the fit fails
    '''
    stars = {}
    for arm in 'brz':
        k0, k1 = _cols(wave, CONT[arm] - CONT_HALF, CONT[arm] + CONT_HALF)
        with np.errstate(all='ignore'):
            cont = np.nanmedian(flux[:, k0:k1].astype(float), axis=1)
        ok = np.isfinite(cont)
        qlo, qhi = np.percentile(cont[ok], STAR_Q)
        stars[arm] = ~ok | (cont < qlo) | (cont > qhi)
    with np.errstate(all='ignore'):
        z = {'5577': np.where(stars['b'], np.nan,
                              line_flux(flux, wave, (LINES['b'] - HALF, LINES['b'] + HALF), SIDES['b']))}
        for name, members in OH_GROUPS.items():
            z[name] = np.sum([np.where(stars[HOLDOUT[m][0]], np.nan, line_flux(flux, wave, HOLDOUT[m][1], HOLDOUT[m][2]))
                              for m in members], axis=0)
    out = {}
    for name in QUALITY:
        r = plane_fit(x, y, z[name], spec)
        if r is None:
            out[name] = (np.nan, np.nan)
        else:
            d = np.radians(r['dir'])
            out[name] = (r['grad'] * np.cos(d), r['grad'] * np.sin(d))
    return out


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
    The plane (gx, gy) and factors left in corrected Sci-fiber spectra:
    in 5577 (b), 8399 (z) and the sum of the r OH lines (r), as
    {arm: dict or None}
    '''
    out = {}
    for arm in 'bzr':
        if arm == 'r':
            z = np.sum([line_flux(flux, wave, HOLDOUT[m][1], HOLDOUT[m][2]) for m in OH_GROUPS['r_oh']], axis=0)
        else:
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
                modelname='', rguard=1.5, outname='', version=''):
    '''
    Write a model-corrected copy of one CFrame, to outname if given, else
    to outdir as <name>.flatmodel.fits.  Returns a summary row (dict) with
    'status' 'ok' or the reason it was skipped
    '''
    info = {'file': filename, 'out': '', 'status': ''}
    try:
        x = fits.open(filename)
    except Exception as err:
        print('ApplyFlatModel: could not open %s (%s)' % (filename, err))
        info['status'] = 'unreadable'
        return info
    hdr = x['PRIMARY'].header
    base = os.path.basename(filename)
    if hdr.get('FGRADFIX', False) or hdr.get('FMODEL', False):
        print('ApplyFlatModel: %s already corrected (FGRADFIX or FMODEL), skipped' % base)
        info['status'] = 'already_corrected'
        return info
    if not str(hdr.get('DRPVER', '')).startswith('1.3.2'):
        print('ApplyFlatModel: warning: %s has DRPVER %s; the model was made for 1.3.2' % (base, hdr.get('DRPVER')))
    mjd, focus, kmir = _num(hdr, 'MJD'), _num(hdr, 'TESCIFO'), _num(hdr, 'TESCIKM')
    if not (np.isfinite(mjd) and np.isfinite(focus) and np.isfinite(kmir) and focus > -900 and kmir > -900):
        print('ApplyFlatModel: %s lacks MJD, TESCIFO or TESCIKM, skipped' % base)
        info['status'] = 'no_focus_or_kmirror'
        return info
    epoch = calib_epoch(hdr, epochs)
    if epoch < 0:
        print('ApplyFlatModel: no calibration epoch for %s, skipped' % base)
        info['status'] = 'no_calib_epoch'
        return info
    app = {arm: applied_factors(calib, epoch, arm) for arm in 'brz'}
    if any(a is None for a in app.values()):
        print('ApplyFlatModel: applied factors of epoch %d not found in %s, skipped' % (epoch, calib))
        info['status'] = 'no_factors'
        return info

    terms, ff = model_terms(mjd, focus, kmir, focus0, ranges)
    val = evaluate(coef, terms)

    st = x['SLITMAP'].data
    sci = np.where(np.char.strip(np.asarray(st['telescope']).astype(str)) == 'Sci')[0]
    xp, yp = np.asarray(st['xpmm'], float)[sci], np.asarray(st['ypmm'], float)[sci]
    spec = np.asarray(st['spectrographid'], int)[sci]
    wave = np.asarray(x['WAVE'].data, float)
    w = arm_weights(wave)
    # quality fibers, as in SkyLinePlane.py: Sci with fibstatus 0, as positions within sci
    qsel = np.where(np.asarray(st['fibstatus'])[sci] == 0)[0]
    q_before = quality(np.asarray(x['FLUX'].data[sci][qsel], np.float32), wave, xp[qsel], yp[qsel], spec[qsel])

    kmod = {}
    for arm in 'brz':
        F = 1 + np.array([val[arm]['F1'], val[arm]['F2'], val[arm]['F3']]) / 100.
        F = F / F.mean()
        kmod[arm] = arm_correction(F / app[arm], val[arm]['gx'], val[arm]['gy'], xp, yp, spec)
    corr = sum(kmod[arm][:, None] * w[arm][None, :] for arm in 'brz')

    flux, ivar = x['FLUX'].data, x['IVAR'].data
    fs = (flux[sci] / corr).astype(np.float32)

    used_direct, gval, res, rval, rmode = False, np.nan, None, np.nan, 'none'
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
                keys = ('f1', 'f2', 'f3', 'gx', 'gy')
                rm = {k: 0.5 * (res['b'][k] + res['z'][k]) for k in keys}
                rmode = 'bz'
                ro = res['r']
                if ro is not None:
                    rval = float(np.hypot(ro['gx'] - rm['gx'], ro['gy'] - rm['gy']))
                    if rval <= rguard and ro['n'] >= 0.7 * nmin:
                        rm = {k: 0.5 * (ro[k] + rm[k]) for k in keys}
                        rmode = 'oh+bz'
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
        for arm in 'bzr':
            r = res[arm] if res else None
            what = 'r OH' if arm == 'r' else arm
            hdr['FMD%sGX' % arm.upper()] = (round(r['gx'], 3) if r else -999., '%s residual plane x (%%/10mm)' % what)
            hdr['FMD%sGY' % arm.upper()] = (round(r['gy'], 3) if r else -999., '%s residual plane y (%%/10mm)' % what)
        hdr['FMDRGRD'] = (round(rval, 3) if np.isfinite(rval) else -1., '|r OH - mean(b,z)| residual plane (%/10mm)')
        hdr['FMDRMODE'] = (rmode, 'r direct: oh+bz, bz (b+z mean only) or none')
    q_after = quality(np.asarray(fs[qsel], np.float32), wave, xp[qsel], yp[qsel], spec[qsel])
    for name in QUALITY:
        for tag, q in (('0', q_before), ('1', q_after)):
            a = float(np.hypot(*q[name]))
            hdr['FMQ%s%s' % (QKEY[name], tag)] = (round(a, 3) if np.isfinite(a) else -1.,
                                                  '%s plane %s correction (%%/10mm)' % (name, 'before' if tag == '0' else 'after'))
    hdr['FMODCODE'] = ((version or code_version())[:40], 'lvm_ksl commit of ApplyFlatModel.py')
    hdr['FMODORIG'] = (os.path.abspath(filename), 'original lvmCFrame')
    hdr.add_history('ApplyFlatModel: Sci FLUX/IVAR corrected with the large-scale flat model (ffactor epoch %d)%s'
                    % (ff, '; direct residual fit applied' if used_direct else ''))

    if not outname:
        root = base[:-5] if base.endswith('.fits') else base
        outname = os.path.join(outdir, root + '.flatmodel.fits')
    if os.path.exists(outname) and os.path.samefile(outname, filename):
        print('ApplyFlatModel: %s would overwrite the original, skipped' % base)
        x.close()
        info['status'] = 'would_overwrite'
        return info
    os.makedirs(os.path.dirname(os.path.abspath(outname)), exist_ok=True)
    x.writeto(outname, overwrite=True)
    x.close()
    print('ApplyFlatModel: %s -> %s  ff %d cal %d  model planes b %.1f,%.1f r %.1f,%.1f z %.1f,%.1f  direct %s%s'
          % (base, outname, ff, epoch, val['b']['gx'], val['b']['gy'], val['r']['gx'], val['r']['gy'],
             val['z']['gx'], val['z']['gy'], used_direct,
             (' (guard %.2f, r %s %.2f)' % (gval, rmode, rval)) if direct else ''))
    info.update({'out': outname, 'status': 'ok', 'expnum': int(_num(hdr, 'EXPOSURE')) if np.isfinite(_num(hdr, 'EXPOSURE')) else -1,
                 'mjd': int(mjd), 'ff': ff, 'calib': epoch, 'direct': used_direct,
                 'guard': round(gval, 3) if np.isfinite(gval) else -1., 'rmode': rmode})
    for arm in 'brz':
        info['%s_gx' % arm] = round(val[arm]['gx'], 3)
        info['%s_gy' % arm] = round(val[arm]['gy'], 3)
    for name in QUALITY:
        k = QKEY[name].lower()
        for tag, q in (('0', q_before), ('1', q_after)):
            a = float(np.hypot(*q[name]))
            info['q%s%s' % (k, tag)] = round(a, 3) if np.isfinite(a) else -999.
        info['q%s1x' % k] = round(float(q_after[name][0]), 3) if np.isfinite(q_after[name][0]) else -999.
        info['q%s1y' % k] = round(float(q_after[name][1]), 3) if np.isfinite(q_after[name][1]) else -999.
    return info



_WORKER = {}


def _init_worker(model, calib, direct, guard, rguard, version):
    '''
    Read the model and the epochs once per worker process
    '''
    coef, focus0, ranges = read_model(model)
    try:
        epochs = load_epochs(default_epochs_file())
    except Exception:
        epochs = None
    _WORKER.update(coef=coef, focus0=focus0, ranges=ranges, calib=calib, epochs=epochs, direct=direct,
                   guard=guard, rguard=rguard, model=model, version=version)


def _work(job):
    '''
    Correct one file; job is (input, output name or '', outdir)
    '''
    fn, outname, outdir = job
    w = _WORKER
    try:
        return apply_model(fn, w['coef'], w['focus0'], w['ranges'], outdir, w['calib'], w['epochs'], w['direct'],
                           w['guard'], w['model'], w['rguard'], outname=outname, version=w['version'])
    except Exception as err:
        print('ApplyFlatModel: %s failed (%s)' % (fn, err))
        return {'file': fn, 'out': '', 'status': 'error: %s' % str(err)[:80]}


def write_summary(rows, summary):
    '''
    One row per exposure, fixed-width
    '''
    cols = ['status', 'expnum', 'mjd', 'ff', 'calib', 'direct', 'guard', 'rmode',
            'b_gx', 'b_gy', 'r_gx', 'r_gy', 'z_gx', 'z_gy',
            'qb0', 'qb1', 'qr0', 'qr1', 'qz0', 'qz1', 'qb1x', 'qb1y', 'qr1x', 'qr1y', 'qz1x', 'qz1y',
            'file', 'out']
    data = {c: [r.get(c, '' if c in ('status', 'rmode', 'file', 'out') else -999) for r in rows] for c in cols}
    data['direct'] = [str(v) for v in data['direct']]
    Table(data).write(summary, format='ascii.fixed_width_two_line', overwrite=True)


def quality_report(rows, nworst=10):
    '''
    Print the quality of the corrected exposures in rows (summary dicts or
    table rows): median |plane| of 5577, r OH and z OH before and after,
    the fraction below 1 %/10 mm, the direct steps rejected, and the
    exposures with the largest held-out plane after the correction
    '''
    ok = [r for r in rows if str(r['status']) == 'ok']
    status = {}
    for r in rows:
        status[str(r['status'])] = status.get(str(r['status']), 0) + 1
    print('\nApplyFlatModel quality report: %d exposures; status %s' % (len(rows), status))
    if not ok:
        return

    def col(k):
        v = np.array([float(r[k]) for r in ok])
        v[v <= -900] = np.nan
        return v

    print('  median |plane| %%/10mm (fraction < 1) [84th percentile], %d corrected exposures' % len(ok))
    for tag, lab in (('0', 'before'), ('1', 'after ')):
        txt = []
        for name, k in (('5577', 'b'), ('r OH', 'r'), ('z OH', 'z')):
            v = col('q%s%s' % (k, tag))
            txt.append('%s %5.2f (%.2f) [%5.2f]' % (name, np.nanmedian(v), np.nanmean(v < 1), np.nanpercentile(v, 84)))
        print('  %s  %s' % (lab, '   '.join(txt)))
    print('  (5577 is fitted by the direct step; r OH and z OH are held out)')
    direct = np.array([str(r['direct']) == 'True' for r in ok])
    print('  direct step applied in %d of %d (%.1f%% rejected by the guard)' % (direct.sum(), len(ok), 100 * (1 - direct.mean())))
    worst = np.fmax(col('qr1'), col('qz1'))
    nbad = int(np.sum(worst > 1))
    print('  %d exposures with r OH or z OH above 1 %%/10mm after the correction%s'
          % (nbad, '; the largest:' if nbad else ''))
    for i in np.argsort(-np.nan_to_num(worst, nan=-1))[:min(nworst, nbad)]:
        r = ok[i]
        print('    %6s  mjd %s  r OH %.2f  z OH %.2f  5577 %.2f  direct %s  guard %s'
              % (r['expnum'], r['mjd'], col('qr1')[i], col('qz1')[i], col('qb1')[i], r['direct'], r['guard']))


def report(summaries):
    '''
    The quality report of one or more summary files (e.g. of all -part jobs)
    '''
    rows = []
    for f in summaries:
        t = Table.read(f, format='ascii.fixed_width_two_line')
        rows += [dict(zip(t.colnames, r)) for r in t]
        print('ApplyFlatModel: %s, %d rows' % (f, len(t)))
    seen, uniq = set(), []
    for r in rows:   # a file may be listed in a later summary as done_before
        key = str(r['file'])
        if str(r['status']) == 'ok' or key not in seen:
            uniq.append(r)
        seen.add(key)
    quality_report(uniq)


def doit(files, model, outdir='.', calib='', direct=False, guard=1.5, rguard=1.5, outroot='', nproc=1,
         redo=False, summary='ApplyFlatModel_summary.txt'):
    '''
    Description:
        Apply the model to each CFrame in files, writing to outdir as
        <name>.flatmodel.fits, or with outroot to the same tree layout
        and names as the originals

    Notes:
        With outroot, files whose output exists are skipped unless redo,
        so a stopped run can be restarted.  The summary is written at
        the end, and every 500 files during a long run.

    History:
        261106 ksl Coding begun
        261112 ksl outroot, nproc, redo, summary
    '''
    if not calib:
        calib = os.environ.get('LVM_MASTER_DIR', '') or \
            os.path.join(os.environ.get('SAS_BASE_DIR', ''), 'sdsswork/lvm/sandbox/calib')
    version = code_version()
    print('ApplyFlatModel: model %s, calib %s, lvm_ksl %s, direct %s' % (model, calib, version, direct))

    rows, jobs = [], []
    for fn in files:
        if outroot:
            outname = tree_path(fn, outroot)
            if os.path.abspath(outname) == os.path.abspath(fn):
                print('ApplyFlatModel: %s is under -outroot; it would be overwritten, skipped' % fn)
                rows.append({'file': fn, 'out': outname, 'status': 'would_overwrite'})
                continue
            if os.path.exists(outname) and not redo:
                rows.append({'file': fn, 'out': outname, 'status': 'done_before'})
                continue
        else:
            outname = ''
        if not os.path.isfile(fn):
            rows.append({'file': fn, 'out': '', 'status': 'missing'})
            continue
        jobs.append((fn, outname, outdir))
    if not outroot:
        os.makedirs(outdir, exist_ok=True)
    print('ApplyFlatModel: %d files to correct, %d skipped beforehand' % (len(jobs), len(rows)))

    args = (model, calib, direct, guard, rguard, version)
    if nproc > 1:
        with multiprocessing.Pool(nproc, initializer=_init_worker, initargs=args) as pool:
            for k, r in enumerate(pool.imap_unordered(_work, jobs, chunksize=4)):
                rows.append(r)
                if (k + 1) % 500 == 0:
                    write_summary(rows, summary)
    else:
        _init_worker(*args)
        for k, job in enumerate(jobs):
            rows.append(_work(job))
            if (k + 1) % 500 == 0:
                write_summary(rows, summary)
    write_summary(rows, summary)
    n = sum(r['status'] == 'ok' for r in rows)
    print('ApplyFlatModel: corrected %d of %d files; summary in %s' % (n, len(files), summary))
    quality_report(rows)


def _init_skysub(version, redux):
    '''
    Run lvmdrp in this worker as reduction version `version` under `redux`,
    so its sky tables go to <redux>/<version>/.../ancillary
    '''
    os.environ['LVMDRP_VERSION'] = version
    os.environ['LVM_SPECTRO_REDUX'] = redux
    from lvmdrp.functions.skyMethod import quick_sky_subtraction
    from lvmdrp import path, __version__
    # importing lvmdrp (sdss_access) resets LVM_SPECTRO_REDUX from its own tree; set it again
    os.environ['LVM_SPECTRO_REDUX'] = redux
    _WORKER.update(skysub=quick_sky_subtraction, path=path, drpver=__version__, tree=os.path.join(redux, version))


def _skywork(job):
    '''
    Sky-subtract one corrected CFrame; job is (cframe, sframe)
    '''
    import time
    cframe, sframe = job
    w = _WORKER
    row = {'cframe': cframe, 'sframe': sframe, 'status': 'ok', 'expnum': -999, 'mjd': -999, 'seconds': -999.}
    t0 = time.time()
    try:
        hdr = fits.getheader(cframe)
        row['expnum'], row['mjd'] = int(hdr['EXPOSURE']), int(hdr['SMJD'])
        anc = w['path'].full('lvm_anc', mjd=row['mjd'], tileid=hdr['TILE_ID'], drpver=w['drpver'],
                             kind='sky', camera='brz', imagetype='table', expnum=row['expnum'])
        if not os.path.abspath(anc).startswith(w['tree'] + os.sep):
            raise RuntimeError('lvmdrp would write the sky table outside %s: %s' % (w['tree'], anc))
        os.makedirs(os.path.dirname(anc), exist_ok=True)
        os.makedirs(os.path.join(os.path.dirname(sframe), 'qa'), exist_ok=True)
        w['skysub'](in_cframe=cframe, out_sframe=sframe)
    except Exception as err:
        print('ApplyFlatModel: sky subtraction of %s failed (%s)' % (cframe, err))
        row['status'] = 'error: %s' % str(err)[:80]
    row['seconds'] = round(time.time() - t0, 1)
    return row


def write_skysub_summary(rows, summary):
    '''
    One row per CFrame of a -skysub run, fixed-width
    '''
    cols = ['status', 'expnum', 'mjd', 'seconds', 'cframe', 'sframe']
    data = {c: [r.get(c, '' if c in ('status', 'cframe', 'sframe') else -999) for r in rows] for c in cols}
    Table(data).write(summary, format='ascii.fixed_width_two_line', overwrite=True)


def doit_skysub(files, outroot='', nproc=1, redo=False, summary='ApplyFlatModel_skysub_summary.txt'):
    '''
    Description:
        Redo the sky subtraction of the corrected CFrames: for each file,
        the corrected CFrame (under outroot, as the correction run wrote
        it, or the file itself without outroot) is sky-subtracted into
        lvmSFrame-<expnum>.fits in the same directory

    Notes:
        SFrames that exist are skipped unless redo.  The CFrames must all
        be in one tree <redux>/<version>/<tilegrp>/<tileid>/<mjd>/, which
        sets the version lvmdrp runs as.

    History:
        261113 ksl Coding begun
    '''
    rows, jobs = [], []
    for fn in files:
        cframe = tree_path(fn, outroot) if outroot else fn
        sframe = os.path.join(os.path.dirname(cframe), os.path.basename(cframe).replace('lvmCFrame', 'lvmSFrame'))
        if sframe == cframe:
            rows.append({'cframe': cframe, 'sframe': '', 'status': 'not_a_cframe'})
        elif not os.path.isfile(cframe):
            rows.append({'cframe': cframe, 'sframe': sframe, 'status': 'missing'})
        elif os.path.exists(sframe) and not redo:
            rows.append({'cframe': cframe, 'sframe': sframe, 'status': 'done_before'})
        else:
            jobs.append((cframe, sframe))
    print('ApplyFlatModel: %d CFrames to sky-subtract, %d skipped beforehand' % (len(jobs), len(rows)))
    if jobs:
        # <redux>/<version>/<tilegrp>/<tileid>/<mjd>/lvmCFrame-<expnum>.fits
        trees = {os.path.dirname(os.path.abspath(c)).rsplit(os.sep, 3)[0] for c, _ in jobs}
        if len(trees) > 1:
            print('ApplyFlatModel: -skysub needs the CFrames in one reduction tree, not %s' % sorted(trees)[:3])
            return
        redux, version = os.path.split(trees.pop())
        print('ApplyFlatModel: sky subtraction as lvmdrp version %s under %s' % (version, redux))
        args = (version, redux)
        if nproc > 1:
            with multiprocessing.Pool(nproc, initializer=_init_skysub, initargs=args) as pool:
                for k, r in enumerate(pool.imap_unordered(_skywork, jobs, chunksize=1)):
                    rows.append(r)
                    if (k + 1) % 500 == 0:
                        write_skysub_summary(rows, summary)
        else:
            _init_skysub(*args)
            for k, job in enumerate(jobs):
                rows.append(_skywork(job))
                if (k + 1) % 500 == 0:
                    write_skysub_summary(rows, summary)
    write_skysub_summary(rows, summary)
    status = {}
    for r in rows:
        status[r['status']] = status.get(r['status'], 0) + 1
    t = [r['seconds'] for r in rows if r['status'] == 'ok']
    print('ApplyFlatModel: sky subtraction of %d files: %s%s; summary in %s'
          % (len(rows), status, ', median %.0f s per file' % np.median(t) if t else '', summary))


def steer(argv):
    '''
    Parse the command line and run doit
    '''
    files, model, outdir, calib, direct, guard, rguard, listfile = [], '', '.', '', False, 1.5, 1.5, ''
    outroot, drpall, base, mjd_range, nproc, redo, summary = '', '', '', None, 1, False, ''
    part, do_report, skysub = None, False, False
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
        elif argv[i] == '-outroot':
            i += 1
            outroot = argv[i]
        elif argv[i] == '-calib':
            i += 1
            calib = argv[i]
        elif argv[i] == '-direct':
            direct = True
        elif argv[i] == '-guard':
            i += 1
            guard = float(argv[i])
        elif argv[i] == '-rguard':
            i += 1
            rguard = float(argv[i])
        elif argv[i] == '-list':
            i += 1
            listfile = argv[i]
        elif argv[i] == '-drpall':
            i += 1
            drpall = argv[i]
        elif argv[i] == '-base':
            i += 1
            base = argv[i]
        elif argv[i] == '-mjd':
            i += 1
            lo, hi = argv[i].split('-')
            mjd_range = (int(lo), int(hi))
        elif argv[i] == '-np':
            i += 1
            nproc = int(argv[i])
        elif argv[i] == '-part':
            i += 1
            try:
                ip, npart = (int(v) for v in argv[i].split('/'))
                assert 1 <= ip <= npart
            except (ValueError, AssertionError):
                print('ApplyFlatModel: -part needs I/N with 1 <= I <= N, not %s' % argv[i])
                return
            part = (ip, npart)
        elif argv[i] == '-redo':
            redo = True
        elif argv[i] == '-report':
            do_report = True
        elif argv[i] == '-skysub':
            skysub = True
        elif argv[i] == '-summary':
            i += 1
            summary = argv[i]
        elif argv[i].startswith('-'):
            print('ApplyFlatModel: unknown option %s' % argv[i])
            print(_USAGE)
            return
        else:
            files.append(argv[i])
        i += 1
    if do_report:
        if not files:
            print('ApplyFlatModel: -report needs one or more summary files')
            return
        report(files)
        return
    if listfile:
        with open(listfile) as f:
            files += [ln.split()[0] for ln in f if ln.strip() and not ln.startswith('#')]
    if drpall:
        xtab = read_drpall(drpall, base=base)
        if mjd_range is not None:
            # the MJD is the directory of each file: .../<tilegrp>/<tileid>/<mjd>/lvmCFrame-<expnum>.fits
            m = np.array([int(os.path.basename(os.path.dirname(str(p)))) for p in xtab['path']])
            xtab = xtab[(m >= mjd_range[0]) & (m <= mjd_range[1])]
        files += [str(p) for p in xtab['path']]
        if not outroot:
            outroot = os.path.join(os.environ.get('SAS_BASE_DIR', ''), DEFAULT_OUTROOT)
    if not files:
        print(_USAGE)
        return
    if skysub:
        if part is not None:
            files = sorted(set(files))[part[0] - 1::part[1]]
            print('ApplyFlatModel: part %d of %d: %d files' % (part[0], part[1], len(files)))
        if not summary:
            summary = 'ApplyFlatModel_skysub_summary%s.txt' % ('_%d_of_%d' % part if part else '')
        doit_skysub(files, outroot, nproc, redo, summary)
        return
    found = find_model(model)
    if not found:
        print('ApplyFlatModel: model %s not found (as given or in %s)' % (model or DEFAULT_MODEL, DATA_DIR))
        return
    model = found
    if part is not None:
        files = sorted(set(files))[part[0] - 1::part[1]]
        print('ApplyFlatModel: part %d of %d: %d files' % (part[0], part[1], len(files)))
    if not summary:
        summary = 'ApplyFlatModel_summary%s.txt' % ('_%d_of_%d' % part if part else '')
    doit(files, model, outdir, calib, direct, guard, rguard, outroot, nproc, redo, summary)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
