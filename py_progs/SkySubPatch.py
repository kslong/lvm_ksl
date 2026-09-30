#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Sky-subtract an LVM lvmCFrame using a sky taken from a PATCH of the
    science field itself -- the compact group of fibers with the least
    nebular emission -- transformed for every fiber to that fiber's own
    wavelength offset, line-spread function and throughput, and write an
    lvmSFrame-layout file with each fiber's own sky in its SKY extension.

Command line usage (if any):

    usage: SkySubPatch.py [-h] [-mode groups|center|slit] [-lsf broaden|none]
                          [-nbkg N] [-bkg FILE] [-ngroup N] [-mask FILE]
                          [-np N] [-out ROOT] cframe [cframe ...]

    Arguments::

        cframe      one or more lvmCFrame FITS files, each processed
                    independently

    Options::

        -h          print this help and exit
        -mode M     calibration unit (default groups):
                      groups  non-overlapping groups of ~ngroup neighbouring
                              fibers on the sky; each fiber gets an
                              inverse-distance average of its 3 nearest
                              groups' calibrations
                      center  every fiber calibrated on itself + its
                              ngroup-1 nearest neighbours (slowest)
                      slit    blocks of ngroup consecutive fiberids within a
                              spectrograph
        -lsf L      broaden (default): at each wavelength broaden whichever
                    of fiber / sky is sharper (fiber FLUX, IVAR and LSF are
                    updated where the fiber is broadened); none: no LSF
                    matching
        -nbkg N     fibers in the background patch (default 20)
        -bkg FILE   use these fiberids (one per line, or a table with a
                    fiberid column) as the background instead of the
                    automatic choice
        -ngroup N   fibers per calibration group (default 7: a fiber and
                    its surrounding hexagon)
        -mask FILE  sky-line mask (palace_make_mask.py output; default
                    lvm_ksl data/sky_mask.fits)
        -np N       parallel worker processes (default 8)
        -out ROOT   output root (default lvmSFrame-<exposure>.patch, in the
                    current directory); with several inputs ROOT_<exp>

Description::

    1. Background patch: every fiber's nebular lines (Halpha, [NII]6583,
       [SII]6716/6731, [OIII]5007; plus Hbeta when the Moon is below the
       horizon, SKY MOON_ALT < 0; never [OII]) and bright sky lines are
       fitted with Gaussians.  After screening (masked pixels, continuum
       above median + 1 sigma = stars, anomalous sky-line throughput), the
       patch is the compact group of nbkg fibers in ONE spectrograph with
       the lowest emission score that is faint in every score line; odd
       fibers (residual scatter, continuum shape, stellar absorption) are
       replaced by the next-nearest candidate when the Moon is up, and only
       reported at dark time.  Background = mean of the patch spectra.

    2. Calibration (per unit, see -mode), from the unit's combined
       spectrum against the background, on sky-line pixels only (nebular
       windows, including [OI]6300/6364, and the b/r and r/z boundary
       regions excluded):
         - wavelength shift: constant in b, quadratic in wavelength through
           per-segment grid-search shifts in r and z;
         - relative LSF: two-sided kernel matching (lsf_kernel.py, adapted
           from Ivan Katkov's lsf_surface_iterative): a kernel in each
           direction, and at each wavelength the sharper spectrum is
           broadened;
         - throughput: one factor per arm, fitted after LSF matching.

    3. Each science fiber's sky = throughput x (background continuum +
       background sky lines shifted and, where the background is the
       sharper, broadened).  Where the FIBER is the sharper it is broadened
       instead: its FLUX and IVAR are convolved and its LSF extension
       updated.  Fibers that are not good science fibers (SkyE, SkyW,
       standards, fibstatus != 0) get the uncorrected background as sky.

    4. Output: FLUX = (broadened) CFrame FLUX - SKY, IVAR (incl. the sky
       variance), MASK, WAVE, LSF (updated), SKY, SKY_IVAR, FLUXCAL_*,
       SLITMAP, plus BACKGROUND (the 1-D background spectrum), BKGFIBERS
       (the patch) and CALIB (per fiber: unit, shift, throughput,
       fraction broadened per arm).  PRIMARY gains SSP* keywords.

    FLUX is the emission in EXCESS of the background patch: any emission
    the patch itself contains (in extended nebulae, a "floor") is
    subtracted from every fiber.

Primary routines::

    select_background  choose the patch and build the background spectrum
    make_units         define the calibration units for -mode
    calibrate_unit     fit shift / LSF kernels / throughput for one unit
    fiber_sky          build one fiber's sky (and broadened fiber) from a
                       calibration
    do_one             process one CFrame and write the SFrame
    steer              command-line driver

Notes:

    Runs in the ksl environment.  Numerical libraries are limited to one
    thread per process (parallelism is over calibration units).  The
    output name contains "SFrame", but QualSFrame.py's maps go through
    GetTelData, which swaps an SFrame name for its CFrame, so they show
    unsubtracted data.  Needs data/sky_mask.fits and
    data/lvm_sky_lines_all.dat from this repository and lsf_kernel.py
    alongside it in py_progs/.

History::

    260930 ksl Coding begun, from the Vela background-sky tests
        (select_background.py, lsf_group_test.py, fiber_test.py):
        background patch = 20 least-emission fibers in one spectrograph
        (lvm_gaussfit fits); calibration units by k-means groups (default),
        per-fiber neighbourhoods or slit blocks; per unit a wavelength
        shift, two-sided LSF kernels (lsf_kernel.py) and per-arm
        throughput fitted on sky lines only, [OI]6300/6364 excluded as
        they can contain source emission.  Moved to py_progs/ on branch
        sky_patch, with repository-relative paths.

'''

import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys
import re
import time
import warnings
import multiprocessing
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table
from scipy.optimize import curve_fit
from scipy.ndimage import median_filter, uniform_filter1d
from scipy.spatial import cKDTree

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
import lsf_kernel as LK
from GetSkyCont import load_mask, _interp_mask_to_wave

warnings.simplefilter('ignore')


def _usage_from_doc(doc):
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc

_USAGE = _usage_from_doc(__doc__)

DEFAULT_MASK = os.path.join(_REPO, 'data', 'sky_mask.fits')
LINELIST = os.path.join(_REPO, 'data', 'lvm_sky_lines_all.dat')
# score lines: centre and fit window (A), the windows sky_gaussfit.py uses (a fixed +-5 A window
# is too narrow to pin down the background under [SII]; it changed the patch chosen)
NEB_FIT = {'ha': (6562.80, 6555., 6570.), 'nii_b': (6583.46, 6578., 6593.), 'sii_a': (6716.44, 6706., 6726.),
           'sii_b': (6730.815, 6721., 6741.), 'oiii_b': (5006.843, 5001., 5012.), 'hb': (4861.325, 4855., 4870.)}
LINES_MOON = ['ha', 'nii_b', 'sii_a', 'sii_b', 'oiii_b']
LINES_DARK = LINES_MOON + ['hb']
SKY_FIT = {'sky5577': (5577.34, 5572., 5582.), 'sky6300': (6300.30, 6295., 6305.), 'sky8344': (8344.61, 8339., 8349.),
           'sky8399': (8399.18, 8394., 8404.), 'sky8827': (8827.11, 8822., 8832.)}
# nebular lines excluded from the calibration fits; [OI]6300/6364 are also bright sky lines, but the
# source (shocks) can add [OI] to some fibers and bias the r-arm throughput, so they are excluded too
NEB_MASK = [3726.03, 3728.82, 4861.33, 4958.91, 5006.84, 5875.6, 6300.30, 6363.78, 6548.05, 6562.80, 6583.45,
            6716.44, 6730.82, 9068.6, 9530.6]
NEB_HALF = 5.0
EXCLUDE = [(5750.0, 5810.0), (7450.0, 7650.0)]   # b/r and r/z boundaries: problematic flux calibration
CONT_CUT = 1.0
THRU_CUT = 2.0
MAX_MASKFRAC = 0.05
MAX_RADIUS_NN = 3.5
PASS_PCT = 20
MAX_ITER = 10
ODD_RMS_CONT, ODD_RMS_SKY, ODD_CONT_FRAC, ODD_STAR, ODD_NSIG = 1.5, 2.0, 0.05, 0.93, 4.0
STAR_BANDS = [(3925.0, 3940.0), (3960.0, 3975.0), (5165.0, 5185.0)]
CONT_BANDS = {'b': (4000, 5500), 'r': (6000, 7400), 'z': (7700, 9500)}
SHIFT_GRID = np.arange(-0.30, 0.3001, 0.005)
SEG_WIDTH = {'R': 200.0, 'Z': 150.0}
SHIFT_POLY_DEG = 2
DIR_SMOOTH = 51
DLAM = 0.5


def rstd(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return 1.4826 * np.median(np.abs(a - np.median(a))) if a.size else np.nan


def continuum(w, spec, clean):
    idx = np.where(clean & np.isfinite(spec))[0]
    sm = median_filter(spec[idx], size=51, mode='nearest')
    return np.interp(w, w[idx], sm)


def arm_of(w):
    out = np.empty(w.size, dtype=object)
    for ch, lo, hi in LK.ARMS:
        out[LK.arm_mask(w, lo, hi)] = ch
    return out


def shift(w, spec, d):
    '''value at w taken from w - d (positive d moves features redward); d may vary with w'''
    return np.interp(w - d, w, spec)


def _gauss(x, a, c, s, b):
    return a * np.exp(-0.5 * ((x - c) / s) ** 2) + b


def fit_line_fluxes(w, flux, ivar, mask, rows, lines):
    '''
    Gaussian flux of each line (name: (centre, wmin, wmax)) in each of rows,
    using lvm_gaussfit.fit_gaussian_to_spectrum exactly as sky_gaussfit.py
    does (same windows, unweighted), so the background choice reproduces the
    validated sky_gaussfit-based selection.  NaN where the fit fails.
    '''
    from lvm_gaussfit import fit_gaussian_to_spectrum
    out = {k: np.full(len(rows), np.nan) for k in lines}
    for name, (wl, wmin, wmax) in lines.items():
        win = (w >= wmin - 2) & (w <= wmax + 2)
        for k, r in enumerate(rows):
            y = flux[r, win]
            if not np.all(np.isfinite(y)):
                continue
            try:
                res, _ = fit_gaussian_to_spectrum(Table([w[win], y], names=['WAVE', 'FLUX']), line=name,
                                                  init_wavelength=wl, init_fwhm=1., wavelength_min=wmin,
                                                  wavelength_max=wmax)
                out[name][k] = float(res['flux_' + name][0])
            except Exception:
                pass
    return out


# --------------------------------------------------------------------------- background
def select_background(d, nbkg, bkg_fiberids=None):
    '''
    Choose the background patch.  Returns dict with rows (CFrame row indices
    of the patch), spec (spectrograph), lines (score lines), moon_alt,
    passes, n_replaced, odd_kept, and the background spectrum (bkg, bkg_err,
    bkg_lsf).
    '''
    w, flux, ivar, mask, slit = d['wave'], d['flux'], d['ivar'], d['mask'], d['slit']
    sci = d['sci']
    moon_alt = float(d['hdr'].get('SKY MOON_ALT', np.nan))
    lines = LINES_DARK if moon_alt < 0 else LINES_MOON
    info = dict(moon_alt=moon_alt, lines=lines, passes=None, n_replaced=0, odd_kept=[])
    if bkg_fiberids is not None:
        pos = {f: i for i, f in enumerate(slit['fiberid'])}
        rows = np.array([pos[f] for f in bkg_fiberids if f in pos])
        info.update(rows=rows, spec=int(np.bincount(np.asarray(slit['spectrographid'])[rows]).argmax()),
                    source='user')
    else:
        t0 = time.time()
        neb = fit_line_fluxes(w, flux, ivar, mask, sci, {k: NEB_FIT[k] for k in lines})
        sky = fit_line_fluxes(w, flux, ivar, mask, sci, SKY_FIT)
        clean = d['clean']
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            fm = np.where(mask[sci] == 0, flux[sci], np.nan)
            cont = np.nanmedian(fm[:, clean], axis=1)
        maskfrac = np.mean(mask[sci] != 0, axis=1)
        thru = np.nanmean([sky[k] / np.nanmedian(sky[k]) for k in SKY_FIT], axis=0)
        score = np.nanmean([(neb[k] - np.nanmedian(neb[k])) / rstd(neb[k]) for k in lines], axis=0)
        cand = ((maskfrac <= MAX_MASKFRAC) & np.isfinite(cont) & (cont <= np.nanmedian(cont) + CONT_CUT * rstd(cont))
                & np.isfinite(thru) & (np.abs(thru - np.nanmedian(thru)) <= THRU_CUT * rstd(thru)) & np.isfinite(score))
        xy = np.column_stack([np.asarray(slit['ra'])[sci] * np.cos(np.radians(np.nanmean(slit['dec'][sci]))),
                              np.asarray(slit['dec'])[sci]])
        nn = np.median(cKDTree(xy).query(xy, k=2)[0][:, 1])
        spec = np.asarray(slit['spectrographid'])[sci]
        pct = {k: np.nanpercentile(neb[k], PASS_PCT) for k in lines}
        best = None
        for k in (1, 2, 3):
            idx = np.where(cand & (spec == k))[0]
            if len(idx) < nbkg:
                continue
            tree = cKDTree(xy[idx])
            for s in range(len(idx)):
                dist, nbi = tree.query(xy[idx[s]], k=nbkg)
                if dist.max() > MAX_RADIUS_NN * nn:
                    continue
                mem = idx[nbi]
                passes = all(np.nanmedian(neb[l][mem]) < pct[l] for l in lines)
                trial = (not passes, np.median(score[mem]), idx[s], k)
                if best is None or trial[:2] < best[:2]:
                    best = trial
        if best is None:
            print('Error: no compact background patch of %d candidate fibers found' % nbkg)
            return None
        failed, _, seed, kspec = best
        replace = moon_alt >= 0
        poolidx = np.where(cand & (spec == kspec))[0]
        order = poolidx[np.argsort(np.hypot(*(xy[poolidx] - xy[seed]).T))]
        excluded = set()
        for _ in range(MAX_ITER):
            members = np.array([i for i in order if i not in excluded][:nbkg])
            odd = _oddity(d, sci[members])
            if not odd.any():
                break
            if not replace:
                info['odd_kept'] = [int(v) for v in np.asarray(slit['fiberid'])[sci[members[odd]]]]
                break
            excluded.update(members[odd])
        info.update(rows=sci[members], spec=int(kspec), passes=not failed, n_replaced=len(excluded),
                    source='auto', radius=float(np.hypot(*(xy[members] - xy[seed]).T).max() / nn),
                    ratio={k: float(np.nanmedian(neb[k][members]) / np.nanmedian(neb[k])) for k in lines},
                    t_select=time.time() - t0)
    rows = info['rows']
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        f = np.where(mask[rows] == 0, flux[rows], np.nan)
        info['bkg'] = np.nanmean(f, axis=0)
        info['bkg_err'] = np.nanstd(f, axis=0) / np.sqrt(np.maximum(np.sum(np.isfinite(f), axis=0), 1))
        info['bkg_lsf'] = np.nanmean(d['lsf'][rows], axis=0)
    return info


def _oddity(d, rows):
    w, clean = d['wave'], d['clean']
    f = np.where(d['mask'][rows] == 0, d['flux'][rows], np.nan).astype(float)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        med = np.nanmedian(f, axis=0)
        sig = np.where(d['ivar'][rows] > 0, 1 / np.sqrt(d['ivar'][rows]), np.nan)
        res = (f - med) / sig
        inr = (w > 3700) & (w < 9600)
        rms_cont = np.sqrt(np.nanmean(res[:, clean & inr] ** 2, axis=1))
        rms_sky = np.sqrt(np.nanmean(res[:, ~clean & inr] ** 2, axis=1))

        def band_ratio(pix):
            ref = np.nanmedian(med[pix])
            r = np.nanmedian(f[:, pix], axis=1) / ref
            e = 1.2533 * np.sqrt(np.nanmean(sig[:, pix] ** 2, axis=1) / np.sum(np.isfinite(f[:, pix]), axis=1)) / abs(ref)
            return r, e
        odd = (rms_cont > ODD_RMS_CONT) | (rms_sky > ODD_RMS_SKY)
        for lo, hi in CONT_BANDS.values():
            r, e = band_ratio((w > lo) & (w < hi) & clean)
            odd |= (np.abs(r - 1) > ODD_CONT_FRAC) & (np.abs(r - 1) > ODD_NSIG * e)
        sr = [band_ratio((w > lo) & (w < hi)) for lo, hi in STAR_BANDS]
        star = np.nanmean([r for r, e in sr], axis=0)
        star_err = np.sqrt(np.nansum([e ** 2 for r, e in sr], axis=0)) / len(sr)
        odd |= (star < ODD_STAR) & ((1 - star) > ODD_NSIG * star_err)
    return odd


# --------------------------------------------------------------------------- calibration units
def make_units(d, mode, ngroup):
    '''
    Returns (units, assign): units = list of arrays of CFrame rows combined for
    each calibration; assign = {row: [(unit index, weight), ...]} for every
    good science fiber.
    '''
    slit = d['slit']
    sci = d['sci']
    xy = np.column_stack([np.asarray(slit['ra'])[sci] * np.cos(np.radians(np.nanmean(slit['dec'][sci]))),
                          np.asarray(slit['dec'])[sci]])
    nn = np.median(cKDTree(xy).query(xy, k=2)[0][:, 1])
    units, assign = [], {}
    if mode == 'center':
        tree = cKDTree(xy)
        for k, r in enumerate(sci):
            _, nb = tree.query(xy[k], k=ngroup)
            units.append(sci[np.atleast_1d(nb)])
            assign[r] = [(k, 1.0)]
    elif mode == 'slit':
        fid = np.asarray(slit['fiberid'])[sci]
        spec = np.asarray(slit['spectrographid'])[sci]
        for k in (1, 2, 3):
            idx = np.where(spec == k)[0]
            idx = idx[np.argsort(fid[idx])]
            for s in range(0, len(idx), ngroup):
                blk = idx[s:s + ngroup]
                if len(blk) < max(3, ngroup // 2) and units:
                    units[-1] = np.concatenate([units[-1], sci[blk]])
                    for i in blk:
                        assign[sci[i]] = [(len(units) - 1, 1.0)]
                    continue
                units.append(sci[blk])
                for i in blk:
                    assign[sci[i]] = [(len(units) - 1, 1.0)]
    else:   # groups: compact clusters of ~ngroup neighbouring fibers on the sky (k-means on position)
        from scipy.cluster.vq import kmeans2
        ncl = max(1, int(round(len(sci) / ngroup)))
        _, lab = kmeans2(xy / nn, ncl, minit='++', seed=1, iter=50)
        members = [np.where(lab == i)[0] for i in range(ncl) if np.any(lab == i)]
        # merge very small clusters into the nearest larger one
        centers = np.array([xy[g].mean(axis=0) for g in members])
        big = [i for i, g in enumerate(members) if len(g) >= max(3, ngroup // 2)]
        final = {i: list(members[i]) for i in big}
        for i, g in enumerate(members):
            if i in final:
                continue
            j = big[int(np.argmin(np.hypot(*(centers[big] - centers[i]).T)))]
            final[j] += list(g)
        keys = sorted(final)
        units = [sci[np.array(final[i])] for i in keys]
        cxy = np.array([xy[final[i]].mean(axis=0) for i in keys])
        ctree = cKDTree(cxy)
        for k, r in enumerate(sci):
            dist, nb = ctree.query(xy[k], k=min(3, len(keys)))
            wgt = 1.0 / (np.atleast_1d(dist) + nn)
            assign[r] = list(zip(np.atleast_1d(nb).tolist(), (wgt / wgt.sum()).tolist()))
    return units, assign


_SH = {}


def _init(shared):
    _SH.update(shared)


def calibrate_unit(task):
    '''fit shift, two-sided LSF kernels and per-arm throughput for one unit's combined spectrum'''
    uid, flux, ivar, mask = task
    w, clean, excl, near, arms = _SH['wave'], _SH['clean'], _SH['excl'], _SH['near'], _SH['arms']
    bkg_l, bkg_err = _SH['bkg_l'], _SH['bkg_err']
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fm = np.where(mask == 0, flux, np.nan)
        f = np.nanmean(fm, axis=0)
        nm = np.sum(np.isfinite(fm) & (ivar > 0), axis=0)
        var = np.nansum(np.where(ivar > 0, 1 / ivar, np.nan), axis=0) / np.maximum(nm, 1) ** 2
    good = np.isfinite(f) & (nm > 0)
    f_l = np.where(good, f - continuum(w, f, clean), 0.0)
    iv = np.where(good & ~excl, 1.0 / (var + bkg_err ** 2), 0.0)
    iv[~np.isfinite(iv)] = 0.0
    wts = np.sqrt(iv)
    skypix = (iv > 0) & ~clean & near

    def best(pix):
        chi = []
        for dd in SHIFT_GRID:
            tpl = shift(w, bkg_l, dd)
            sc = np.sum(wts[pix] ** 2 * f_l[pix] * tpl[pix]) / max(np.sum(wts[pix] ** 2 * tpl[pix] ** 2), 1e-300)
            chi.append(np.sum((wts[pix] * (f_l[pix] - sc * tpl[pix])) ** 2))
        return float(SHIFT_GRID[int(np.argmin(chi))])

    dfun = np.full(w.size, best(skypix & (arms == 'B')) if np.any(skypix & (arms == 'B')) else 0.0)
    for a in ('R', 'Z'):
        inarm = arms == a
        ok = inarm & (iv > 0)
        if not ok.any():
            continue
        lo, hi = w[ok].min(), w[ok].max()
        cen, dd, wt = [], [], []
        for s0 in np.arange(lo, hi, SEG_WIDTH[a]):
            pix = skypix & inarm & (w >= s0) & (w < s0 + SEG_WIDTH[a])
            if pix.sum() < 30:
                continue
            cen.append(s0 + 0.5 * SEG_WIDTH[a])
            dd.append(best(pix))
            wt.append(np.sqrt(np.sum(wts[pix] ** 2 * bkg_l[pix] ** 2)))
        if cen:
            deg = min(SHIFT_POLY_DEG, len(cen) - 1)
            coef = np.polyfit(np.array(cen) - 8000.0, dd, deg, w=wt)
            dfun[inarm] = np.clip(np.polyval(coef, w[inarm] - 8000.0), SHIFT_GRID[0], SHIFT_GRID[-1])
    bsh = shift(w, bkg_l, dfun)

    def arm_scales(target, source):
        s = np.ones(3)
        for k, a in enumerate(('B', 'R', 'Z')):
            pix = skypix & (arms == a)
            den = np.sum(wts[pix] ** 2 * source[pix] ** 2)
            if den > 0:
                s[k] = np.sum(wts[pix] ** 2 * target[pix] * source[pix]) / den
        return s

    armidx = np.select([arms == 'B', arms == 'R'], [0, 1], 2)
    out = dict(uid=uid, dfun=dfun.astype(np.float32))
    if _SH['lsf_mode'] == 'broaden':
        s0 = arm_scales(f_l, bsh)
        m = LK.two_sided_match(w, f_l, bsh, iv, fit_source=s0[armidx] * bsh, smooth_pix=DIR_SMOOTH)
        thr = arm_scales(m['target_matched'], m['source_matched'])
        out.update(k1=m['k1']['surface'].astype(np.float32), k2=m['k2']['surface'].astype(np.float32),
                   sig1=m['sigma1'].astype(np.float32), sig2=m['sigma2'].astype(np.float32))
    else:
        thr = arm_scales(f_l, bsh)
    out['thr'] = thr
    return out


def fiber_sky(d, bkginfo, cal, row):
    '''
    One fiber's sky (and possibly broadened FLUX/variance/LSF) from its
    (interpolated) calibration.  Returns flux_out, var_out, lsf_out, sky, skyvar, frac_broadened (per arm)
    '''
    w, arms = d['wave'], d['arms']
    armidx = np.select([arms == 'B', arms == 'R'], [0, 1], 2)
    bkg_c, bkg_l, bkg_err = bkginfo['bkg_c'], bkginfo['bkg_l'], bkginfo['bkg_err']
    thr = cal['thr'][armidx]
    bsh = shift(w, bkg_l, cal['dfun'])
    f = d['flux'][row].astype(float)
    var = np.where(d['ivar'][row] > 0, 1 / d['ivar'][row], np.nan)
    lsf = d['lsf'][row].astype(float).copy()
    fiber_broadened = np.zeros(w.size, bool)
    if 'k1' in cal:
        src_sharper = uniform_filter1d(cal['sig1'] - cal['sig2'], size=DIR_SMOOTH) >= 0
        lines = np.where(src_sharper, LK.apply_kernel(w, bsh, cal['k1']), bsh)
        fiber_broadened = ~src_sharper
        fb = LK.apply_kernel(w, np.nan_to_num(f), cal['k2'])
        vb = LK.apply_kernel(w, np.nan_to_num(var), cal['k2'] ** 2)
        f = np.where(fiber_broadened, fb, f)
        var = np.where(fiber_broadened, vb, var)
        lsf = np.where(fiber_broadened, np.sqrt(lsf ** 2 + (2.3548 * cal['sig2'] * DLAM) ** 2), lsf)
    else:
        lines = bsh
    sky = thr * (bkg_c + lines)
    skyvar = (thr * bkg_err) ** 2
    frac = [float(np.mean(fiber_broadened[arms == a])) for a in ('B', 'R', 'Z')]
    return f, var, lsf, sky, skyvar, frac


def _combine_cals(cals, weights):
    out = {'thr': sum(wt * c['thr'] for c, wt in zip(cals, weights)),
           'dfun': sum(wt * c['dfun'].astype(float) for c, wt in zip(cals, weights))}
    if 'k1' in cals[0]:
        for key in ('k1', 'k2', 'sig1', 'sig2'):
            out[key] = sum(wt * c[key].astype(float) for c, wt in zip(cals, weights))
    return out


# --------------------------------------------------------------------------- driver
def load(filename, mask_wave, mask_bool):
    x = fits.open(filename)
    slit = Table(x['SLITMAP'].data)
    w = x['WAVE'].data.astype(float)
    d = dict(x=x, hdr=x['PRIMARY'].header, slit=slit, wave=w, flux=x['FLUX'].data, ivar=x['IVAR'].data,
             mask=x['MASK'].data, lsf=x['LSF'].data.astype(float),
             sci=np.where((slit['telescope'] == 'Sci') & (slit['fibstatus'] == 0))[0],
             clean=_interp_mask_to_wave(mask_wave, mask_bool, w), arms=arm_of(w))
    excl = np.min(np.abs(w[:, None] - np.array(NEB_MASK)[None, :]), axis=1) <= NEB_HALF
    for lo, hi in EXCLUDE:
        excl |= (w > lo) & (w < hi)
    d['excl'] = excl
    return d


def label_near(w, lsf_fwhm):
    '''pixels within 3 sigma of a listed sky line'''
    t = ascii.read(LINELIST)
    t = t[(t['wave'] > w[0]) & (t['wave'] < w[-1])]
    near = np.zeros(w.size, bool)
    for wl in t['wave']:
        j = np.searchsorted(w, wl)
        sig = lsf_fwhm[min(j, w.size - 1)] / 2.355
        sel = slice(max(j - 12, 0), min(j + 12, w.size))
        near[sel] |= np.abs(w[sel] - wl) <= 3 * sig
    return near


def do_one(filename, mode='groups', lsf_mode='broaden', nbkg=20, bkg_fiberids=None, ngroup=7,
           mask_wave=None, mask_bool=None, nproc=8, outroot=''):
    t_start = time.time()
    try:
        d = load(filename, mask_wave, mask_bool)
    except Exception as e:
        print('Error: could not read %s (%s)' % (filename, e))
        return None
    expnum = int(d['hdr'].get('EXPOSURE', -1))
    bk = select_background(d, nbkg, bkg_fiberids)
    if bk is None:
        return None
    bk['bkg_c'] = continuum(d['wave'], bk['bkg'], d['clean'])
    bk['bkg_l'] = np.where(np.isfinite(bk['bkg']), bk['bkg'] - bk['bkg_c'], 0.0)
    print('%s: background patch %d fibers in spectrograph %d (%s); Moon alt %.1f -> score lines %s'
          % (os.path.basename(filename), len(bk['rows']), bk['spec'], bk['source'], bk['moon_alt'], ','.join(bk['lines'])))
    if bk['source'] == 'auto':
        print('    patch/field line ratio: ' + '  '.join('%s %.2f' % kv for kv in bk['ratio'].items()) +
              '; faint in every line: %s; replaced %d, odd kept %d' % (bk['passes'], bk['n_replaced'], len(bk['odd_kept'])))

    units, assign = make_units(d, mode, ngroup)
    near = label_near(d['wave'], np.nanmedian(d['lsf'][d['sci']], axis=0))
    shared = dict(wave=d['wave'], clean=d['clean'], excl=d['excl'], near=near, arms=d['arms'], bkg_l=bk['bkg_l'],
                  bkg_err=bk['bkg_err'], lsf_mode=lsf_mode)
    tasks = [(i, d['flux'][u], d['ivar'][u], d['mask'][u]) for i, u in enumerate(units)]
    t0 = time.time()
    if nproc > 1:
        with multiprocessing.get_context('spawn').Pool(nproc, initializer=_init, initargs=(shared,)) as pool:
            cals = pool.map(calibrate_unit, tasks, chunksize=1)
    else:
        _init(shared)
        cals = [calibrate_unit(t) for t in tasks]
    cals = sorted(cals, key=lambda c: c['uid'])
    t_cal = time.time() - t0
    print('    %d calibration units (mode %s, lsf %s) in %.1f s' % (len(units), mode, lsf_mode, t_cal))

    # build every fiber's sky
    nrow, npix = d['flux'].shape
    flux_out = d['flux'].astype(float).copy()
    var_out = np.where(d['ivar'] > 0, 1.0 / np.where(d['ivar'] > 0, d['ivar'], 1), np.nan)
    lsf_out = d['lsf'].copy()
    sky2d = np.tile(bk['bkg'], (nrow, 1))
    skyvar2d = np.tile(bk['bkg_err'] ** 2, (nrow, 1))
    calrows = []
    unit_of = {}
    for i, u in enumerate(units):
        for r in u:
            unit_of.setdefault(int(r), i)
    for r in range(nrow):
        if r in assign:
            pairs = assign[r]
            cal = _combine_cals([cals[i] for i, _ in pairs], [wt for _, wt in pairs])
            f, v, l, s, sv, frac = fiber_sky(d, bk, cal, r)
            flux_out[r], var_out[r], lsf_out[r], sky2d[r], skyvar2d[r] = f, v, l, s, sv
            k5577 = np.argmin(np.abs(d['wave'] - 5577.3))
            k8400 = np.argmin(np.abs(d['wave'] - 8400.0))
            calrows.append([int(d['slit']['fiberid'][r]), unit_of.get(r, -1), 1, float(cal['dfun'][k5577]),
                            float(cal['dfun'][k8400])] + [float(v) for v in cal['thr']] + frac)
        else:
            calrows.append([int(d['slit']['fiberid'][r]), -1, 0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out_flux = flux_out - sky2d
        tot_var = var_out + skyvar2d
        out_ivar = np.where(np.isfinite(tot_var) & (tot_var > 0), 1.0 / tot_var, 0.0)
        sky_ivar = np.where(skyvar2d > 0, 1.0 / skyvar2d, 0.0)
    out_mask = d['mask'].copy()
    out_mask[~np.isfinite(sky2d) | ~np.isfinite(out_flux)] = 1

    root = outroot or 'lvmSFrame-%08d.patch' % expnum
    write_output(d, root + '.fits', out_flux, out_ivar, out_mask, lsf_out, sky2d, sky_ivar, bk, calrows,
                 dict(mode=mode, lsf_mode=lsf_mode, ngroup=ngroup, nunits=len(units), source=filename))
    d['x'].close()
    print('    wrote %s  (total %.1f s)' % (root + '.fits', time.time() - t_start))
    return root + '.fits'


def write_output(d, outfits, out_flux, out_ivar, out_mask, lsf_out, sky2d, sky_ivar, bk, calrows, opts):
    x = d['x']
    hdr = x['PRIMARY'].header.copy()
    cards = [('SSPMODE', opts['mode'], 'SkySubPatch: calibration unit mode'),
             ('SSPLSF', opts['lsf_mode'], 'SkySubPatch: LSF matching'),
             ('SSPNGRP', opts['ngroup'], 'SkySubPatch: fibers per calibration unit'),
             ('SSPNUNIT', opts['nunits'], 'SkySubPatch: number of calibration units'),
             ('SSPNBKG', len(bk['rows']), 'SkySubPatch: background patch fibers'),
             ('SSPBSPEC', bk['spec'], 'SkySubPatch: background patch spectrograph'),
             ('SSPBSRC', bk['source'], 'SkySubPatch: patch auto-selected or user'),
             ('SSPLINES', ','.join(bk['lines']), 'SkySubPatch: score lines'),
             ('SSPSRC', os.path.basename(opts['source']), 'SkySubPatch: input CFrame')]
    for k, v, c in cards:
        hdr[k] = (v, c)

    def img(data, name, dtype, like='FLUX'):
        return fits.ImageHDU(data=np.asarray(data).astype(dtype), header=x[like].header, name=name)

    hdus = [fits.PrimaryHDU(header=hdr), img(out_flux, 'FLUX', np.float32), img(out_ivar, 'IVAR', np.float32),
            img(out_mask, 'MASK', np.uint8, 'MASK'), x['WAVE'].copy(), img(lsf_out, 'LSF', np.float32, 'LSF'),
            img(sky2d, 'SKY', np.float32), img(sky_ivar, 'SKY_IVAR', np.float32)]
    for ext in ('FLUXCAL_STD', 'FLUXCAL_SCI', 'FLUXCAL_MOD'):
        if ext in x:
            hdus.append(x[ext].copy())
    hdus.append(x['SLITMAP'].copy())
    bt = Table([d['wave'], bk['bkg'], bk['bkg_err'], bk['bkg_lsf']], names=['WAVE', 'BKG', 'BKG_ERR', 'BKG_LSF'])
    hdus.append(fits.BinTableHDU(bt, name='BACKGROUND'))
    bf = Table()
    bf['fiberid'] = np.asarray(d['slit']['fiberid'])[bk['rows']]
    bf['odd_kept'] = np.isin(bf['fiberid'], bk.get('odd_kept', []))
    hdus.append(fits.BinTableHDU(bf, name='BKGFIBERS'))
    ct = Table(rows=calrows, names=['fiberid', 'unit', 'corrected', 'shift_5577', 'shift_8400', 'thr_b', 'thr_r',
                                    'thr_z', 'fbroad_b', 'fbroad_r', 'fbroad_z'])
    hdus.append(fits.BinTableHDU(ct, name='CALIB'))
    fits.HDUList(hdus).writeto(outfits, overwrite=True)


def steer(argv):
    files = []
    mode, lsf_mode, nbkg, ngroup, nproc = 'groups', 'broaden', 20, 7, 8
    bkg_file, mask_file, outroot = '', DEFAULT_MASK, ''
    i = 1
    while i < len(argv):
        a = argv[i]
        if a == '-h':
            print(_USAGE)
            return
        elif a == '-mode':
            i += 1
            mode = argv[i]
        elif a == '-lsf':
            i += 1
            lsf_mode = argv[i]
        elif a == '-nbkg':
            i += 1
            nbkg = int(argv[i])
        elif a == '-bkg':
            i += 1
            bkg_file = argv[i]
        elif a == '-ngroup':
            i += 1
            ngroup = int(argv[i])
        elif a == '-mask':
            i += 1
            mask_file = argv[i]
        elif a == '-np':
            i += 1
            nproc = int(argv[i])
        elif a == '-out':
            i += 1
            outroot = argv[i]
        elif a.startswith('-'):
            print('Error: unknown option "%s"' % a)
            print(_USAGE)
            return
        else:
            files.append(a)
        i += 1
    if not files:
        print(_USAGE)
        return
    if mode not in ('groups', 'center', 'slit'):
        print('Error: -mode must be groups, center or slit')
        return
    if lsf_mode not in ('broaden', 'none'):
        print('Error: -lsf must be broaden or none')
        return
    if not os.path.exists(mask_file):
        print('Error: mask file not found: %s' % mask_file)
        return
    bkg_fiberids = None
    if bkg_file:
        try:
            t = ascii.read(bkg_file)
            col = 'fiberid' if 'fiberid' in t.colnames else t.colnames[0]
            bkg_fiberids = [int(v) for v in t[col]]
        except Exception as e:
            print('Error: could not read -bkg file %s (%s)' % (bkg_file, e))
            return
    mask_wave, mask_bool = load_mask(mask_file)
    for f in files:
        root = outroot
        if outroot and len(files) > 1:
            m = re.search(r'(\d{5,8})', os.path.basename(f))
            root = '%s_%s' % (outroot, int(m.group(1)) if m else os.path.basename(f).replace('.fits', ''))
        do_one(f, mode, lsf_mode, nbkg, bkg_fiberids, ngroup, mask_wave, mask_bool, nproc, root)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
