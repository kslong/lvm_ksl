#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Measure, in each science exposure, the large-scale fiber-flat error
    left in the lvmCFrames: three spectrograph factors plus a plane
    across the science IFU, fitted to the flux of one airglow line per
    arm, and write one table row per exposure.

Command line usage (if any):

    usage: SkyLinePlane.py [-h] [-out table.txt] [-calib DIR] [-notracers] [-list FILE] cframe ...
           SkyLinePlane.py [-h] [-out table.txt] [-calib DIR] [-notracers] [-base DIR] [-every N] -drpall file

    Arguments::

        cframe       lvmCFrame files (a FixFlatGradient.py .gradfix copy
                     works as well)

    Options::

        -h           print this help and exit
        -out FILE    output table (default SkyLinePlane.txt); during a
                     long run it is rewritten every 500 exposures
        -calib DIR   directory holding <epoch>/lvm-mfiberflat_twilight-<arm>.fits,
                     used for the predicted gradient (default
                     $LVM_MASTER_DIR, else $SAS_BASE_DIR/sdsswork/lvm/sandbox/calib)
        -list FILE   read CFrame names from FILE (one per line), in
                     addition to those on the command line
        -drpall FILE take the exposures from a drpall file (fits or h5):
                     the lvmCFrame next to each lvmSFrame in its location
                     column (relative to $SAS_BASE_DIR); its calib_mjd,
                     focus, K-mirror and altitude columns are joined to
                     the output table
        -base DIR    resolve the drpall locations under DIR instead of
                     $SAS_BASE_DIR
        -every N     read only every Nth exposure of the drpall file,
                     for a quick trial
        -notracers   skip the nebular tracers (which need the SKY_EAST
                     and SKY_WEST extensions, about 3x the reading)

Description::

    The DRP's stage-2 fiber-flat correction
    (run_calseq.create_fiberflats_corrections) fits spectrograph factors
    and an IFU gradient to sky lines in the science exposures of one
    window per "ffactor epoch" and records them in the master twilight
    flat (FACTOR1-3, GCOEFF1-4).  In DRP 1.3.2 the factors were applied
    to the flats but the gradient was not (sdss/lvmdrp#268).  So in a
    1.3.2 CFrame the plane fitted here is the missing gradient plus
    whatever changed since the window was measured; in a .gradfix
    CFrame it is the change alone.

    For each CFrame and arm (b 5577, r 6364, z 8399 A):

    1.  only the columns needed are read (memmap): +-4 A around the
        line, two clean continuum side bands next to it, and +-25 A
        around the DRP's star-rejection continuum (4600/6760/8170 A)
    2.  line flux of each Sci fiber (fibstatus 0) = the sum over the
        window of FLUX minus the median of the side bands
    3.  stars are rejected as the DRP does: fibers whose continuum is
        outside the 5-97 percentiles of the Sci fibers
    4.  z_i/median(z) = f_s(i) * (1 + a*y_i + b*x_i) is fitted, with
        x,y = xpmm,ypmm and s the spectrograph, clipping 4 x MAD
        outliers; the factors are normalized to a mean of 1
    5.  written: grad = 1000*hypot(a,b) (% per 10 mm), dir = atan2(a,b)
        (deg; 0 = +x, 90 = +y), their formal errors, the residual rms
        per fiber (%), the factors, and sp1-3 = the mean of the model
        over each spectrograph's Sci fibers minus 1 (%)
    6.  predicted (p prefix): the same fit applied to the GCOEFF
        gradient of the master flat of the exposure's calibration epoch,
        i.e. what a 1.3.2 CFrame shows if nothing changed
    7.  the same fit to the continuum (c prefix), only a check, and only
        meaningful at bright Moon (otherwise stars and nebulae dominate)

    8.  held-out sky lines, not used in the fit, to test it: the same
        plane fit, per line, to OH lines in r (6864, 6923, 7316, 7341)
        and z (8344, 8430, 8465, 8827, 8886, 8920), to their sums per arm
        (r_oh, z_oh); b has no other usable airglow line (Hg 5461 has a
        per-fiber rms of ~70%), so its held-out test is z_oh, through the
        achromatic changing part; the r OH lines sit in the telluric O2 B
        and H2O bands, so their side bands are poorer than those in z;
        written as <line>_gx, _gy (the plane as x, y components, % per
        10 mm), _egrad, _rms, _n
    9.  nebular tracers (unless -notracers): line flux of the science
        fibers minus that of the sky (SKYEW*SKY_EAST + SKYWW*SKY_WEST, as
        the DRP weights them) for Hb, [OIII] 5007, [OI] 6300, Ha,
        [NII] 6583, [SII] 6717+6731 and, as a direct measure of target
        emission in the fitted lines, 5577, 6364 and 8399 themselves
        (t5577, t6364, t8399); written as <tracer>_sci, _sky (medians
        over the Sci fibers), _exc (median of science - sky), _srms
        (robust rms over the fibers of science - sky) and _gx, _gy (a
        plain linear fit of science - sky across the IFU, flux per 10
        mm), in the units of the line fluxes
    10. spectrograph factors from the sky fibers: for SkyE and SkyW
        separately, the line flux (5577, 6364, 8399 and the r_oh and
        z_oh sums) of each good sky fiber, fibers whose continuum is an
        outlier (4 x MAD, e.g. a star) rejected; the median per
        spectrograph divided by the mean of the three medians.  The sky
        is uniform across each small sky IFU, so these are pure
        spectrograph factors with no IFU gradient, to compare with the
        factors of the science-fiber fit (which can trade off against
        the IFU plane).  Like those, they are relative to the flat the
        CFrame was reduced with (in 1.3.2 the stage-2 factors were
        applied to all fibers).  Written as <ifu>_<line>_f1, _f2, _f3,
        _ef (the mean error of a factor) and _n (fibers used), with ifu
        skye or skyw and line b, r, z, r_oh or z_oh

    Per exposure it also records MJD, tile, DRPQUAL, FGRADFIX, the
    calibration epoch, the field PA (POSCIPA), the K-mirror angle
    (TESCIKM), the focus of the three telescopes (TESCIFO etc.), the
    altitude, lab and CCD temperatures, Moon illumination, altitude and
    separation, Sun altitude and galactic latitude.

Primary routines:

    doit
    measure
    plane_fit

Notes::

    Needs numpy, scipy, astropy, pyyaml and FixFlatGradient.py (in this
    directory), not lvmdrp, so it can run where the reductions are
    stored, e.g. at Utah:

        SkyLinePlane.py -out SkyLinePlane_1.3.2.txt -drpall drpall-1.3.2.fits

    Each row of a CFrame spans only a few disk pages, so reading any set
    of windows reads nearly all of FLUX (89 MB); the tracers also read
    SKY_EAST and SKY_WEST.  About 0.2 s per CFrame on a local disk
    without the tracers, roughly 3x that with them.  The CFrames must not
    be compressed, or each file has to be decompressed in full.

    The calibration epoch is taken from the mwave master named in the
    COMMENT cards of the header, else from lvmcore
    calibration-epochs.yaml (the rule of FixFlatGradient.py).

    [OI] 6364 can come from the target as well as the sky, so the r arm
    has a few outliers in fields with nebular emission.

    The 2026-11-05 pilot on 66 exposures found the residual plane is
    achromatic, constant within a night, its y slope following the
    science telescope focus TESCIFO (about -4 % per 10 mm per mm) and its
    x slope sin(K-mirror angle).  Obsidian DRP-Flat-Field-Gradient-2611.md.

History::

    261105 ksl Coding begun (pilot as FlatCorr_analysis/skyline_plane.py),
               moved to lvm_ksl so it can run at Utah; -drpall added
    261106 ksl Header keywords with no value are written as NaN, not None
    261106 ksl Held-out sky lines (OH in r and z) and nebular
               tracers (science minus sky at Hb, [OIII], [OI] 6300, Ha,
               [NII], [SII] and the fitted lines), to test the direct
               fit on all observations; -notracers; the Sci rows of
               FLUX are now read once per file
    261108 ksl Spectrograph factors from the SkyE and SkyW sky fibers

'''

import os
import re
import sys
import numpy as np
from astropy.io import fits
from astropy.table import Table, join
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import least_squares

from FixFlatGradient import load_epochs, find_epoch, default_epochs_file

LINES = {'b': 5577.35, 'r': 6363.78, 'z': 8399.18}
HALF = 4.0
SIDES = {'b': ((5566.0, 5573.0), (5582.0, 5586.0)),
         'r': ((6350.0, 6355.0), (6366.5, 6378.0)),
         'z': ((8393.0, 8396.5), (8402.5, 8412.5))}
CONT = {'b': 4600.0, 'r': 6760.0, 'z': 8170.0}
CONT_HALF = 25.0
STAR_Q = (5, 97)
CLIP = 4.0
# held-out sky lines: name -> (arm whose star mask is used, line window, side bands)
HOLDOUT = {'oh6864': ('r', (6861.0, 6867.0), ((6852.0, 6858.0), (6874.0, 6878.0))),
           'oh6923': ('r', (6920.0, 6926.5), ((6916.0, 6918.5), (6929.0, 6935.0))),
           'oh7316': ('r', (7313.0, 7319.5), ((7307.0, 7310.0), (7321.0, 7325.0))),
           'oh7341': ('r', (7337.5, 7344.0), ((7332.0, 7335.5), (7347.0, 7354.0))),
           'oh8344': ('z', (8341.5, 8347.5), ((8335.0, 8339.0), (8357.0, 8362.0))),
           'oh8430': ('z', (8427.0, 8433.5), ((8420.5, 8424.0), (8436.0, 8442.0))),
           'oh8465': ('z', (8462.0, 8469.0), ((8456.0, 8459.5), (8471.0, 8478.0))),
           'oh8827': ('z', (8824.0, 8830.5), ((8815.0, 8821.0), (8840.0, 8846.0))),
           'oh8886': ('z', (8882.5, 8889.0), ((8874.0, 8880.0), (8892.0, 8899.0))),
           'oh8920': ('z', (8916.0, 8923.0), ((8907.0, 8913.0), (8926.0, 8933.0)))}
OH_GROUPS = {'r_oh': ['oh6864', 'oh6923', 'oh7316', 'oh7341'],
             'z_oh': ['oh8344', 'oh8430', 'oh8465', 'oh8827', 'oh8886', 'oh8920']}
# nebular tracers, measured as science minus sky: name -> (arm whose star mask is used, window, side bands)
TRACERS = {'hb': ('b', (4857.5, 4865.0), ((4845.0, 4852.0), (4870.0, 4877.0))),
           'oiii': ('b', (5003.0, 5011.0), ((4990.0, 4997.0), (5015.0, 5022.0))),
           't5577': ('b', (5577.35 - 4.0, 5577.35 + 4.0), ((5566.0, 5573.0), (5582.0, 5586.0))),
           'oi6300': ('r', (6296.5, 6304.5), ((6288.0, 6294.0), (6307.0, 6313.0))),
           't6364': ('r', (6363.78 - 4.0, 6363.78 + 4.0), ((6350.0, 6355.0), (6366.5, 6378.0))),
           'ha': ('r', (6559.5, 6566.0), ((6528.0, 6538.0), (6608.0, 6618.0))),
           'nii': ('r', (6580.0, 6587.0), ((6528.0, 6538.0), (6608.0, 6618.0))),
           'sii': ('r', (6712.5, 6735.0), ((6684.0, 6692.0), (6752.0, 6760.0))),
           't8399': ('z', (8399.18 - 4.0, 8399.18 + 4.0), ((8393.0, 8396.5), (8402.5, 8412.5)))}
SKY_IFUS = {'skye': 'SkyE', 'skyw': 'SkyW'}
SKY_MIN = 5   # fewest good fibers per spectrograph for a sky-fiber factor
DRPALL_COLUMNS = ['calib_mjd', 'tileid', 'exptime', 'sci_alt', 'sci_kmpos', 'sci_focpos', 'skye_focpos',
                  'skyw_focpos', 'sun_alt', 'moon_alt']


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before the History: line, so -h stays short
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


_USAGE = _usage_from_doc(__doc__)


def plane_fit(x, y, z, spec):
    '''
    Fit z/median(z) = f_s * (1 + a*y + b*x) with iterative clipping;
    a dict of the results, or None if fewer than 100 good fibers
    '''
    good = np.isfinite(z) & (z > 0)
    if good.sum() < 100:
        return None
    zn = z / np.median(z[good])
    idx = np.searchsorted(np.array([1, 2, 3]), spec)

    def resid(p, m):
        return zn[m] - p[:3][idx[m]] * (1 + p[3] * y[m] + p[4] * x[m])

    p = np.array([1., 1., 1., 0., 0.])
    m = good.copy()
    for _ in range(6):
        p = least_squares(resid, p, args=(m,)).x
        res = zn - p[:3][idx] * (1 + p[3] * y + p[4] * x)
        mad = 1.4826 * np.median(np.abs(res[m] - np.median(res[m])))
        mnew = good & (np.abs(res) < CLIP * mad)
        if np.array_equal(mnew, m):
            break
        m = mnew
    r = least_squares(resid, p, args=(m,))
    p = r.x
    sig2 = np.sum(r.fun ** 2) / max(1, m.sum() - 5)
    try:
        cov = np.linalg.inv(r.jac.T @ r.jac) * sig2
        ea, eb, cab = np.sqrt(cov[3, 3]), np.sqrt(cov[4, 4]), cov[3, 4]
    except np.linalg.LinAlgError:
        ea = eb = cab = np.nan
    fac = p[:3] / np.mean(p[:3])
    a, b = p[3], p[4]
    h = np.hypot(a, b)
    if h > 0:
        egrad = 1000 * np.sqrt(a * a * ea ** 2 + b * b * eb ** 2 + 2 * a * b * cab) / h
        edir = np.degrees(np.sqrt(b * b * ea ** 2 + a * a * eb ** 2 - 2 * a * b * cab) / h ** 2)
    else:
        egrad = edir = np.nan
    model = fac[idx] * (1 + a * y + b * x)
    spm = [100 * (np.mean(model[spec == s]) / np.mean(model) - 1) for s in (1, 2, 3)]
    return dict(n=int(m.sum()), f1=fac[0], f2=fac[1], f3=fac[2], grad=1000 * h, egrad=egrad,
                dir=np.degrees(np.arctan2(a, b)), edir=edir, rms=100 * np.sqrt(sig2),
                sp1=spm[0], sp2=spm[1], sp3=spm[2])


def _cols(wave, lo, hi):
    i = np.where((wave >= lo) & (wave <= hi))[0]
    return i[0], i[-1] + 1


def line_flux(arr, wave, window, sides):
    '''
    Line flux of each row of arr (fibers x wavelength): the sum over the
    window minus the median of the side bands, NaN if the window has a
    non-finite pixel
    '''
    i0, i1 = _cols(wave, *window)
    line = arr[:, i0:i1].astype(float)
    side = np.hstack([arr[:, slice(*_cols(wave, lo, hi))] for lo, hi in sides]).astype(float)
    cont = np.nanmedian(side, axis=1)
    z = np.nansum(line - cont[:, None], axis=1) * np.median(np.diff(wave))
    z[~np.all(np.isfinite(line), axis=1)] = np.nan
    return z


def _rob(a):
    a = a[np.isfinite(a)]
    return 1.4826 * np.median(np.abs(a - np.median(a))) if len(a) > 10 else np.nan


def _plane_cols(row, name, r):
    '''
    Write the held-out plane fit r (from plane_fit) as x, y components
    '''
    if r:
        d = np.radians(r['dir'])
        row.update({name + '_gx': r['grad'] * np.cos(d), name + '_gy': r['grad'] * np.sin(d),
                    name + '_egrad': r['egrad'], name + '_rms': r['rms'], name + '_n': r['n']})
    else:
        row.update({name + '_gx': np.nan, name + '_gy': np.nan, name + '_egrad': np.nan,
                    name + '_rms': np.nan, name + '_n': 0})


def sky_factors(z, cont, spec):
    '''
    Spectrograph factors from the fibers of one sky IFU: z the line
    flux, cont their continuum (for rejecting outliers such as stars),
    spec the spectrograph of each fiber.  The median of z per
    spectrograph divided by the mean of the three; a dict with f1-f3,
    ef (the mean error of a factor) and n, NaN if a spectrograph has
    fewer than SKY_MIN usable fibers
    '''
    out = {'f1': np.nan, 'f2': np.nan, 'f3': np.nan, 'ef': np.nan, 'n': 0}
    ok = np.isfinite(z) & np.isfinite(cont) & (z > 0)
    if ok.sum() < 3 * SKY_MIN:
        return out
    c = cont[ok]
    mad = 1.4826 * np.median(np.abs(c - np.median(c)))
    if mad > 0:
        ok[ok] = np.abs(c - np.median(c)) < CLIP * mad
    med, err = [], []
    for s in (1, 2, 3):
        zs = z[ok & (spec == s)]
        if len(zs) < SKY_MIN:
            return out
        med.append(np.median(zs))
        err.append(1.2533 * 1.4826 * np.median(np.abs(zs - np.median(zs))) / np.sqrt(len(zs)))
    mean = np.mean(med)
    out.update(f1=med[0] / mean, f2=med[1] / mean, f3=med[2] / mean, ef=float(np.mean(err)) / mean, n=int(ok.sum()))
    return out


def read_gcoeff(calib, epoch, arm):
    '''
    GCOEFF1-4 of the master twilight flat of an epoch and arm, or None
    '''
    p = os.path.join(calib, str(epoch), 'lvm-mfiberflat_twilight-%s.fits' % arm)
    if not os.path.isfile(p):
        return None
    h = fits.getheader(p)
    c = [h.get('%s FIBERFLAT GCOEFF%d' % (arm.upper(), k)) for k in (1, 2, 3, 4)]
    return None if any(v is None for v in c) else c


def calib_epoch(header, epochs):
    '''
    The calibration epoch: the directory of the mwave master in the
    COMMENT cards, else the rule of FixFlatGradient.py; -1 if neither
    '''
    text = ''.join(str(c) for c in header.get('COMMENT', []))
    m = re.search(r'calib/(\d{5})/lvm-mwave-', text)
    if m:
        return int(m.group(1))
    mjd = int(header.get('MJD', header.get('SMJD', 0)))
    e = find_epoch(mjd, epochs) if epochs is not None else None
    return -1 if e is None else e


def _num(header, key):
    '''
    A header value as a float, NaN if missing, empty or not a number
    '''
    try:
        return float(header[key])
    except (KeyError, TypeError, ValueError):
        return np.nan


def measure(fname, calib, epochs, gcache, tracers=True):
    '''
    One table row (a dict) for one CFrame
    '''
    row = {'file': os.path.basename(fname)}
    with fits.open(fname, memmap=True) as f:
        h = f[0].header
        wave = f['WAVE'].data.astype(float)
        sm = f['SLITMAP'].data
        sci = np.where((sm['telescope'] == 'Sci') & (sm['fibstatus'] == 0))[0]
        x = np.asarray(sm['xpmm'], float)[sci]
        y = np.asarray(sm['ypmm'], float)[sci]
        spec = np.asarray(sm['spectrographid'], int)[sci]
        flux = np.asarray(f['FLUX'].data[sci], dtype=np.float32)
        row.update(expnum=int(h.get('EXPOSURE', -1)), mjd=int(h.get('MJD', -1)), smjd=int(h.get('SMJD', -1)),
                   tile=int(h.get('TILE_ID', -1)), drpqual=int(h.get('DRPQUAL', -1)),
                   drpver=str(h.get('DRPVER', '')), gradfix=bool(h.get('FGRADFIX', False)),
                   fluxcal=str(h.get('FLUXCAL', '')), exptime=_num(h, 'EXPTIME'),
                   pa=_num(h, 'POSCIPA'), kmir=_num(h, 'TESCIKM'),
                   focus=_num(h, 'TESCIFO'), skye_focus=_num(h, 'TESKYEFO'),
                   skyw_focus=_num(h, 'TESKYWFO'), alt=_num(h, 'SCIALT'),
                   labtemp=_num(h, 'LABTEMP'), ccdtemp=_num(h, 'CCDTEMP1'),
                   moon_fli=_num(h, 'SKY MOON_FLI'), moon_alt=_num(h, 'SKY MOON_ALT'),
                   moon_sep=_num(h, 'SKY SCI_MOON_SEP'), sun_alt=_num(h, 'SKY SUN_ALT'))
        try:
            row['glat'] = SkyCoord(h['SCIRA'] * u.deg, h['SCIDEC'] * u.deg).galactic.b.deg
        except (KeyError, ValueError, TypeError):
            row['glat'] = np.nan
        epoch = calib_epoch(h, epochs)
        row['epoch'] = epoch
        stars = {}
        for arm in 'brz':
            z = line_flux(flux, wave, (LINES[arm] - HALF, LINES[arm] + HALF), SIDES[arm])
            k0, k1 = _cols(wave, CONT[arm] - CONT_HALF, CONT[arm] + CONT_HALF)
            cont = np.nanmedian(flux[:, k0:k1].astype(float), axis=1)
            ok = np.isfinite(cont)
            qlo, qhi = np.percentile(cont[ok], STAR_Q)
            star = ~ok | (cont < qlo) | (cont > qhi)
            stars[arm] = star
            zl = np.where(star, np.nan, z)
            row['%s_cont' % arm] = np.median(cont[ok])
            row['%s_line' % arm] = np.nanmedian(zl)
            for tag, val in (('', zl), ('c', np.where(star, np.nan, cont))):
                r = plane_fit(x, y, val, spec)
                for k in ('n', 'f1', 'f2', 'f3', 'grad', 'egrad', 'dir', 'edir', 'rms', 'sp1', 'sp2', 'sp3'):
                    if tag == 'c' and k in ('f1', 'f2', 'f3', 'egrad', 'edir'):
                        continue
                    row['%s_%s%s' % (arm, tag, k)] = r[k] if r else np.nan
            key = (epoch, arm)
            if key not in gcache:
                gc = read_gcoeff(calib, epoch, arm) if epoch > 0 else None
                if gc is None:
                    gcache[key] = None
                else:
                    G = gc[0] + gc[1] * y + gc[2] * x + gc[3] * x * y
                    gcache[key] = plane_fit(x, y, G, spec)
            r = gcache[key]
            for k in ('grad', 'dir', 'sp1', 'sp3'):
                row['%s_p%s' % (arm, k)] = r[k] if r else np.nan
        # held-out sky lines, singly and summed per arm
        hz = {}
        for name, (arm, window, sides) in HOLDOUT.items():
            hz[name] = np.where(stars[arm], np.nan, line_flux(flux, wave, window, sides))
            _plane_cols(row, name, plane_fit(x, y, hz[name], spec))
        for name, members in OH_GROUPS.items():
            _plane_cols(row, name, plane_fit(x, y, np.sum([hz[m] for m in members], axis=0), spec))
        # spectrograph factors from the sky fibers of each sky IFU
        for ifu, tel in SKY_IFUS.items():
            sk = np.where((sm['telescope'] == tel) & (sm['fibstatus'] == 0))[0]
            sflux = np.asarray(f['FLUX'].data[sk], dtype=np.float32) if len(sk) else None
            sspec = np.asarray(sm['spectrographid'], int)[sk]
            for name in ('b', 'r', 'z', 'r_oh', 'z_oh'):
                if sflux is None:
                    r = sky_factors(np.array([]), np.array([]), np.array([], int))
                else:
                    arm = name[0]
                    k0, k1 = _cols(wave, CONT[arm] - CONT_HALF, CONT[arm] + CONT_HALF)
                    scont = np.nanmedian(sflux[:, k0:k1].astype(float), axis=1)
                    if name in LINES:
                        zz = line_flux(sflux, wave, (LINES[name] - HALF, LINES[name] + HALF), SIDES[name])
                    else:
                        zz = np.sum([line_flux(sflux, wave, HOLDOUT[m][1], HOLDOUT[m][2]) for m in OH_GROUPS[name]], axis=0)
                    r = sky_factors(zz, scont, sspec)
                for k in ('f1', 'f2', 'f3', 'ef', 'n'):
                    row['%s_%s_%s' % (ifu, name, k)] = r[k]
        # nebular tracers: science minus sky
        if tracers:
            we, ww = _num(h, 'SKYEW'), _num(h, 'SKYWW')
            if not (np.isfinite(we) and np.isfinite(ww)) or we + ww <= 0:
                we = ww = 0.5
            sky = we * np.asarray(f['SKY_EAST'].data[sci], dtype=np.float32) + \
                ww * np.asarray(f['SKY_WEST'].data[sci], dtype=np.float32)
            for name, (arm, window, sides) in TRACERS.items():
                ls = line_flux(flux, wave, window, sides)
                lk = line_flux(sky, wave, window, sides)
                neb = np.where(stars[arm], np.nan, ls - lk)
                g = np.isfinite(neb)
                row[name + '_sci'] = np.nanmedian(np.where(stars[arm], np.nan, ls))
                row[name + '_sky'] = np.nanmedian(np.where(stars[arm], np.nan, lk))
                row[name + '_exc'] = np.nanmedian(neb) if g.any() else np.nan
                row[name + '_srms'] = _rob(neb)
                if g.sum() > 100:
                    A = np.vstack([np.ones(g.sum()), x[g], y[g]]).T
                    c, *_ = np.linalg.lstsq(A, neb[g], rcond=None)
                    row[name + '_gx'], row[name + '_gy'] = 10 * c[1], 10 * c[2]
                else:
                    row[name + '_gx'] = row[name + '_gy'] = np.nan
    return row


def read_drpall(filename, base='', every=1):
    '''
    The exposures of a drpall file (fits or h5): expnum, the path of the
    lvmCFrame next to each lvmSFrame, and the columns in DRPALL_COLUMNS
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
    xtab['path'] = [os.path.join(base, str(loc).strip().lstrip('/')).replace('lvmSFrame-', 'lvmCFrame-')
                    for loc in xtab['location']]
    keep = ['expnum', 'path'] + [c for c in DRPALL_COLUMNS if c in xtab.colnames]
    return xtab[keep]


def _write(rows, out, drpall_tab):
    xtab = Table(rows=rows)
    if drpall_tab is not None and len(xtab):
        extra = drpall_tab.copy()
        extra.remove_column('path')
        for c in extra.colnames:
            if c != 'expnum':
                extra.rename_column(c, 'drpall_' + c)
        xtab = join(xtab, extra, keys='expnum', join_type='left')
        xtab.sort('expnum')
    for n in xtab.colnames:
        if xtab[n].dtype.kind == 'f':
            xtab[n].format = '.4g'
    xtab.write(out, format='ascii.fixed_width_two_line', overwrite=True)
    return xtab


def doit(files, out='SkyLinePlane.txt', calib='', drpall_tab=None, tracers=True):
    '''
    Description:
        Measure each CFrame in files and write the table to out

    Notes:
        Files that are missing or fail are reported and skipped; the
        table is rewritten every 500 exposures so a long run can be
        inspected (or survive being stopped)

    History:
        261105 ksl Coding begun
    '''
    if not calib:
        calib = os.environ.get('LVM_MASTER_DIR', '') or \
            os.path.join(os.environ.get('SAS_BASE_DIR', ''), 'sdsswork/lvm/sandbox/calib')
    try:
        epochs = load_epochs(default_epochs_file())
    except Exception:
        epochs = None
    rows, failed, gcache = [], [], {}
    for k, fn in enumerate(files):
        if not os.path.isfile(fn):
            failed.append((fn, 'missing'))
            continue
        try:
            rows.append(measure(fn, calib, epochs, gcache, tracers))
        except Exception as err:
            failed.append((fn, str(err)))
            print('SkyLinePlane: failed on %s: %s' % (fn, err))
        if (k + 1) % 500 == 0:
            _write(rows, out, drpall_tab)
            print('SkyLinePlane: %d / %d (%d failed)' % (k + 1, len(files), len(failed)))
            sys.stdout.flush()
    if not rows:
        print('SkyLinePlane: nothing measured')
        return None
    xtab = _write(rows, out, drpall_tab)
    if failed:
        with open(os.path.splitext(out)[0] + '_failed.txt', 'w') as f:
            for fn, why in failed:
                f.write('%s %s\n' % (fn, why))
    print('SkyLinePlane: wrote %d rows to %s; %d files missing or failed' % (len(xtab), out, len(failed)))
    return xtab


def steer(argv):
    '''
    Parse the command line and run doit
    '''
    files, out, calib, listfile, drpall, base, every, tracers = [], 'SkyLinePlane.txt', '', '', '', '', 1, True
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
        elif argv[i] == '-list':
            i += 1
            listfile = argv[i]
        elif argv[i] == '-drpall':
            i += 1
            drpall = argv[i]
        elif argv[i] == '-base':
            i += 1
            base = argv[i]
        elif argv[i] == '-notracers':
            tracers = False
        elif argv[i] == '-every':
            i += 1
            every = int(argv[i])
        elif argv[i].startswith('-'):
            print('SkyLinePlane: unknown option %s' % argv[i])
            print(_USAGE)
            return
        else:
            files.append(argv[i])
        i += 1
    if listfile:
        with open(listfile) as f:
            files += [ln.split()[0] for ln in f if ln.strip() and not ln.startswith('#')]
    drpall_tab = None
    if drpall:
        drpall_tab = read_drpall(drpall, base, every)
        files += list(drpall_tab['path'])
    if not files:
        print(_USAGE)
        return
    doit(files, out, calib, drpall_tab, tracers)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
