#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Perform sky subtraction on an XCframe summary file using
    sky_decomp.lsf_surface_iterative.SkyDecompLSFSurfaceIterative (the
    same physically-motivated, per-row-LSF-refined decomposition
    DecomposeCleanSky.py uses) for the continuum/line separation,
    instead of the B-spline fit in SkySubDev1.py or the PALACE model in
    SkySubDev2.py.  In all other respects the algorithm is the same as
    SkySubDev1.py: for each row the continuum and line residuals of the
    science and sky spectra are separated, a scale factor for the sky
    lines is found by bisection, and the scaled sky is subtracted.

    Two methods are supported::

        nearest           continuum and lines both from the nearest
                          sky telescope
        farlines_nearcont lines from the far sky, continuum from the
                          near sky (default)

    Produces a FITS file with WAVE, FLUX (sky-subtracted), SKY, and
    DRP_ALL extensions -- the same layout as SkySubOrig/Drp/Dev1/Dev2.py,
    so SkySub_eval.py reads and compares this file exactly like those,
    with no changes needed there.

Command line usage (if any):

    usage: SkySubDev3.py [-method METHOD] [-delta N] [-v VEL] [-lmc] [-smc]
                         [-out ROOT] [-lvmsky_skysub PATH] filename

    Arguments::

        filename    XCframe FITS file to process

    Options::

        -method METHOD   sky subtraction method: nearest |
                         farlines_nearcont  (default: farlines_nearcont)
        -delta N         process every N-th row; useful for quick tests
                         (default: 1 = all rows)
        -v VEL           nebular systemic velocity (km/s) used to
                         Doppler-shift the exclusion windows masked out
                         of every spectrum's continuum/line fit (default:
                         0 -- see Description).
        -lmc / -smc      shortcuts for the LMC (~262 km/s) / SMC (~146
                         km/s) systemic velocity, same convention as
                         sky_gaussfit.py.  -v overrides these if both
                         are given.
        -out ROOT        output filename root; default is
                         <stem>_dev3_<method>
        -lvmsky_skysub PATH
                         path to the lvmsky repo's skysub/ directory,
                         which supplies the sky_decomp package this
                         script imports (default: ~/SDSS/lvmsky/skysub).

    Unlike SkySubDev1.py/SkySubDev2.py, no -mask argument exists here --
    SkyDecomp needs no external palace_mask file; it fits its own OH/
    atomic-line families directly against the data.

Description:

    The SkyDecompLSFSurfaceIterative instance is built once (from the
    full wavelength grid) and reused for every row -- see
    DecomposeCleanSky.build_decomp, imported directly rather than
    duplicated.  For each row:

    1. Science and sky spectra are read for the row.
    2. Near/far sky is determined from sci_ra/dec, skye_ra/dec,
       skyw_ra/dec in DRP_ALL (identical to SkySubDev1.py).
    3. Each of the three spectra (science, near sky, far sky) is fit
       with SkyDecomp, with sky_gaussfit.resolve_nebular_lines' windows
       (Doppler-shifted by -v/-lmc/-smc) excluded from the fit via
       ivar=0 -- the same mask-and-wrap approach DecomposeCleanSky.py
       uses. This matters for the two sky-telescope spectra specifically:
       DecomposeCleanSky.py/sky_nebular_leak_eval.py's validation work
       found real nebular-line leak in SKY_WEST on at least one tested
       exposure, so masking it out of the *continuum* fit (rather than
       assuming, as SkySubDev1/Dev2/Drp implicitly do, that the sky
       telescopes are nebula-free) keeps that leak from biasing the
       fitted continuum.  cont = bestfit_lsf - (oh+atom+orc+o2 line
       components); lines = spectrum - cont (same PredictSky.py
       LINE_KEYS convention, and the same "lines = observed - continuum"
       definition SkySubDev1.py uses -- kept identical so the two are
       comparable apples-to-apples).
    4. A global line scale factor r is found by ksl_bisection (from
       SkySubOrig.py), exactly as in SkySubDev1.py.
    5. For farlines_nearcont: sky = cont_near + r * lines_far
       For nearest:           sky = cont_near + r * lines_near
    6. sky-subtracted science = flux_sci - sky

    QA flag bits stored in DRP_ALL['QA_FLAGS']::

        0x01  NANDATA    NaN/inf found in input flux or sky data
        0x02  ZEROSKY    sky line vector is all-zero; scale unreliable
        0x04  NOTSOLVED  at least one of the three per-row SkyDecomp
                         fits reported a status other than Solved/
                         AlmostSolved (see SkyDecompResult.fit_status)
        0x08  FAILED     row raised an exception; spectrum filled with NaN

Notes:

    Output filename is ``<ROOT>.fits``.  If -out is omitted the name is
    derived as ``<stem>_dev3_<method>.fits``.

    Imports sky_decomp from an external lvmsky checkout (-lvmsky_skysub,
    default ~/SDSS/lvmsky/skysub) via sys.path, same as
    DecomposeCleanSky.py; only works when that checkout's tree contains
    the package.

History::

    260909  ksl  Coding begun.

'''

import sys
import os

# ensure py_progs siblings are importable when running directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

from SkySubOrig import fit_func, ksl_bisection, obstime_to_mjd
from sky_gaussfit import resolve_nebular_lines

# ──────────────────────────────────────────────────────────────
# DecomposeCleanSky.py (py_dev/) supplies build_decomp -- importing it
# triggers its own module-level -lvmsky_skysub argv parsing and sys.path
# setup for the lvmsky sky_decomp package, so nothing further is needed
# here beyond making py_dev/ importable.  DEFAULT_LVMSKY_SKYSUB below
# must match DecomposeCleanSky.DEFAULT_LVMSKY_SKYSUB -- both default to
# the same path, and a -lvmsky_skysub given on this script's own command
# line is parsed identically by both (same flag spelling).
# ──────────────────────────────────────────────────────────────

DEFAULT_LVMSKY_SKYSUB = os.path.expanduser('~/SDSS/lvmsky/skysub')

_PY_DEV_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py_dev')
sys.path.insert(0, _PY_DEV_DIR)

from DecomposeCleanSky import build_decomp, FACTOR  # noqa: E402

# ──────────────────────────────────────────────────────────────
# QA flag bits
# ──────────────────────────────────────────────────────────────
QA_NANDATA   = 1   # NaN or inf values found in input flux or sky data
QA_ZEROSKY   = 2   # sky line vector is all-zero; scale factor is unreliable
QA_NOTSOLVED = 4   # at least one of the three per-row SkyDecomp fits was not Solved/AlmostSolved
QA_FAILED    = 8   # row failed entirely; FLUX and SKY are NaN

_QA_FLAG_NAMES = {
    QA_NANDATA:   'NANDATA',
    QA_ZEROSKY:   'ZEROSKY',
    QA_NOTSOLVED: 'NOTSOLVED',
    QA_FAILED:    'FAILED',
}

_SOLVED_STATUSES = {'Solved', 'AlmostSolved'}
_LINE_KEYS = ('oh', 'atom', 'orc', 'o2')  # matches PredictSky.py's own convention

LMC_VEL = 262.
SMC_VEL = 146.

_USAGE = '''Usage:
  SkySubDev3.py [-method METHOD] [-delta N] [-v VEL] [-lmc] [-smc]
                [-out ROOT] [-lvmsky_skysub PATH] filename

Arguments:
  filename         XCframe FITS file to process

Options:
  -method METHOD   nearest | farlines_nearcont
                   (default: farlines_nearcont)
  -delta N         step size through rows for quick tests (default: 1)
  -v VEL           nebular systemic velocity (km/s, default: 0)
  -lmc / -smc      LMC (~262 km/s) / SMC (~146 km/s) shortcuts for -v
  -out ROOT        output filename root (default: <stem>_dev3_<method>)
  -lvmsky_skysub PATH
                   path to lvmsky's skysub/ dir (default: ~/SDSS/lvmsky/skysub)
'''


def resolve_velocity(vel=None, lmc=False, smc=False):
    '''Resolve the nebular systemic velocity from -v/-lmc/-smc (see
    DecomposeCleanSky.py/sky_nebular_leak_eval.py's identical helper).'''
    if vel is not None:
        return float(vel)
    if lmc:
        return LMC_VEL
    if smc:
        return SMC_VEL
    return 0.0


def nebular_mask(wave, vel=0.0):
    '''Boolean mask, True where a pixel is clean of every genuinely-
    nebular line's Doppler-shifted window (see DecomposeCleanSky.py's
    identical function for the sky-line-coincidence rationale).'''
    wave = np.asarray(wave, dtype=float)
    zz = 1.0 + vel / 3e5
    resolved, _dropped = resolve_nebular_lines(vel)
    clean = np.ones(wave.shape, dtype=bool)
    for _name, _center, wmin, wmax in resolved:
        clean &= ~((wave >= zz * wmin) & (wave <= zz * wmax))
    return clean


# ──────────────────────────────────────────────────────────────
# Per-row sky subtraction
# ──────────────────────────────────────────────────────────────

def _decompose(decomp, spec, ivar):
    '''
    Fit spec with decomp (nebular windows already reflected in ivar),
    return (cont, lines, status).  cont/lines are in the spectrum's
    native (non-FACTOR-scaled) units; lines = spec - cont.
    '''
    fit = decomp.fit(spec * FACTOR, ivar, verbose=False)
    line_total = np.zeros_like(spec, dtype=np.float64)
    for key in _LINE_KEYS:
        arr = fit.components.get(key)
        if arr is not None:
            line_total = line_total + np.asarray(arr, dtype=np.float64)
    bestfit = np.asarray(fit.bestfit_lsf, dtype=np.float64) / FACTOR
    cont = bestfit - line_total / FACTOR
    lines = spec - cont
    return cont, lines, fit.fit_status


def one_dev3(xfits, drp_all, row, wave, decomp, clean_mask,
            method='farlines_nearcont'):
    '''
    Sky-subtract a single row from an open XCframe FITS object.

    Uses the shared SkyDecompLSFSurfaceIterative instance (decomp) and
    nebular-line exclusion mask (clean_mask) to separate continuum and
    lines in sci, near-sky, and far-sky spectra, then scales and
    subtracts the sky.

    Returns (result_table, qa_flags, line_scale).  result_table has
    columns WAVE, SCI_FLUX, SKY.  line_scale is the ksl_bisection factor
    r applied to the sky lines (SKY = cont_near + r*use_lines; the
    continuum itself is not scaled).

    Returns (None, QA_FAILED, nan) on error.
    '''
    qa_flags = 0

    flux      = np.array(xfits['FLUX'].data[row],     dtype=float)
    skye_flux = np.array(xfits['SKY_EAST'].data[row], dtype=float)
    skyw_flux = np.array(xfits['SKY_WEST'].data[row], dtype=float)

    if not (np.all(np.isfinite(flux)) and
            np.all(np.isfinite(skye_flux)) and
            np.all(np.isfinite(skyw_flux))):
        qa_flags |= QA_NANDATA
        flux      = np.nan_to_num(flux,      nan=0.0, posinf=0.0, neginf=0.0)
        skye_flux = np.nan_to_num(skye_flux, nan=0.0, posinf=0.0, neginf=0.0)
        skyw_flux = np.nan_to_num(skyw_flux, nan=0.0, posinf=0.0, neginf=0.0)

    # determine near/far sky from angular separation
    sci_coord  = SkyCoord(ra=drp_all['sci_ra'][row]  * u.degree,
                          dec=drp_all['sci_dec'][row] * u.degree)
    skye_coord = SkyCoord(ra=drp_all['skye_ra'][row]  * u.degree,
                          dec=drp_all['skye_dec'][row] * u.degree)
    skyw_coord = SkyCoord(ra=drp_all['skyw_ra'][row]  * u.degree,
                          dec=drp_all['skyw_dec'][row] * u.degree)
    if sci_coord.separation(skye_coord) < sci_coord.separation(skyw_coord):
        sky_near, sky_far = skye_flux, skyw_flux
    else:
        sky_near, sky_far = skyw_flux, skye_flux

    ivar = np.where(clean_mask, 1.0, 0.0)

    cont_sci,  lines_sci,  st_sci  = _decompose(decomp, flux,     ivar)
    cont_near, lines_near, st_near = _decompose(decomp, sky_near, ivar)
    cont_far,  lines_far,  st_far  = _decompose(decomp, sky_far,  ivar)

    if not all(st in _SOLVED_STATUSES for st in (st_sci, st_near, st_far)):
        qa_flags |= QA_NOTSOLVED

    if method == 'farlines_nearcont':
        use_lines = lines_far
    elif method == 'nearest':
        use_lines = lines_near
    else:
        print('Error: unknown method "%s"' % method)
        return None, QA_FAILED, np.nan

    if np.dot(use_lines, use_lines) == 0:
        qa_flags |= QA_ZEROSKY

    r = ksl_bisection(fit_func, 0.5, 1.5, tol=0.001, maxiter=8,
                      args=(lines_sci, use_lines))

    sky          = cont_near + r * use_lines
    sci_flux_sub = flux - sky

    result = Table([wave, sci_flux_sub, sky], names=['WAVE', 'SCI_FLUX', 'SKY'])
    return result, qa_flags, r


# ──────────────────────────────────────────────────────────────
# Batch processing
# ──────────────────────────────────────────────────────────────

def do_all(filename, method='farlines_nearcont', idelta=1, vel=0.0,
          outroot='', lvmsky_skysub=DEFAULT_LVMSKY_SKYSUB):
    '''
    Process every row of an XCframe file and write sky-subtracted output.

    Parameters
    ----------
    filename : str
        Path to the input XCframe FITS file.
    method : str
        Sky subtraction method (nearest or farlines_nearcont).
    idelta : int
        Row step size (1 = all rows).
    vel : float
        Nebular systemic velocity (km/s) for the exclusion mask.
    outroot : str
        Output filename root; defaults to <stem>_dev3_<method>.
    lvmsky_skysub : str
        Path to lvmsky's skysub/ directory.
    '''
    x       = fits.open(filename)
    drp_all = Table(x['DRP_ALL'].data)
    wave    = np.array(x['WAVE'].data, dtype=float)

    clean_mask = nebular_mask(wave, vel=vel)
    print('Nebular-excluded pixels: %d / %d (vel=%.1f km/s)'
          % (int((~clean_mask).sum()), len(wave), vel))

    decomp = build_decomp(wave, lvmsky_skysub)

    nan_spectrum = np.full(len(wave), np.nan)

    final_flux      = []
    final_sky       = []
    select          = []
    qa_flags_list   = []
    line_scale_list = []

    i = 0
    while i < len(drp_all):
        try:
            ftab, row_flags, line_scale = one_dev3(
                x, drp_all, i, wave, decomp, clean_mask, method=method)
        except Exception as e:
            print('Row %d: exception (%s)' % (i, e))
            ftab = None
            row_flags = 0
            line_scale = np.nan

        if ftab is None:
            row_flags |= QA_FAILED
            final_flux.append(nan_spectrum.copy())
            final_sky.append(nan_spectrum.copy())
        else:
            final_flux.append(np.array(ftab['SCI_FLUX']))
            final_sky.append(np.array(ftab['SKY']))

        select.append(i)
        qa_flags_list.append(row_flags)
        line_scale_list.append(line_scale)
        i += idelta
        if i % 100 == 0:
            print('Completed %6d of %d in steps of %d' % (i, len(drp_all), idelta))

    n_failed = sum(1 for f in qa_flags_list if f & QA_FAILED)
    n_warned = sum(1 for f in qa_flags_list if f != 0 and not (f & QA_FAILED))
    print('\nProcessed %d rows: %d failed (NaN fill), %d with warnings'
          % (len(select), n_failed, n_warned))

    flagged = [(select[j], qa_flags_list[j])
               for j in range(len(select)) if qa_flags_list[j] != 0]
    if flagged:
        print('Rows with QA flags:')
        for orig_row, flags in flagged:
            active = [name for bit, name in _QA_FLAG_NAMES.items() if flags & bit]
            print('  Row %6d  flags=0x%02x  (%s)' % (orig_row, flags, ', '.join(active)))

    hdu1 = fits.PrimaryHDU(data=None)
    hdu1.header['Title']  = 'SkySubDev3'
    hdu1.header['METHOD'] = method
    hdu1.header['NEBVEL'] = vel
    hdu2 = fits.ImageHDU(data=wave.astype(np.float32),         name='WAVE')
    hdu3 = fits.ImageHDU(data=np.array(final_flux),            name='FLUX')
    hdu4 = fits.ImageHDU(data=np.array(final_sky),             name='SKY')

    xtab = drp_all[select].copy()
    xtab['QA_FLAGS']   = np.array(qa_flags_list, dtype=np.int32)
    xtab['LINE_SCALE'] = np.array(line_scale_list, dtype=float)
    if 'obstime' in xtab.colnames and 'mjd' in xtab.colnames:
        xtab['mjd'] = obstime_to_mjd(xtab['obstime'])
    hdu5 = fits.BinTableHDU(xtab, name='DRP_ALL')

    from astropy.wcs import WCS
    dwave = float(wave[1] - wave[0]) if len(wave) > 1 else 0.5
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [1, 1]
    wcs.wcs.crval = [float(wave[0]), 0]
    wcs.wcs.cdelt = [dwave, 1]
    wcs.wcs.ctype = ['WAVE', 'LINE']
    hdu3.header.update(wcs.to_header())
    hdu4.header.update(wcs.to_header())

    hdul = fits.HDUList([hdu1, hdu2, hdu3, hdu4, hdu5])

    if outroot == '':
        stem    = os.path.splitext(os.path.basename(filename))[0]
        outroot = '%s_dev3_%s' % (stem, method)
    outfile = '%s.fits' % outroot
    hdul.writeto(outfile, overwrite=True)
    print('Wrote results to %s' % outfile)
    x.close()


# ──────────────────────────────────────────────────────────────
# Command-line entry point
# ──────────────────────────────────────────────────────────────

if __name__ == '__main__':
    argv = sys.argv[1:]
    if not argv or '-h' in argv or '--help' in argv:
        print(_USAGE)
        sys.exit(0)

    method        = 'farlines_nearcont'
    idelta        = 1
    vel_arg       = None
    lmc           = False
    smc           = False
    outroot       = ''
    lvmsky_skysub = DEFAULT_LVMSKY_SKYSUB
    filename      = None

    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == '-method':
            i += 1
            method = argv[i]
        elif arg == '-delta':
            i += 1
            idelta = int(argv[i])
        elif arg == '-v':
            i += 1
            vel_arg = float(argv[i])
        elif arg == '-lmc':
            lmc = True
        elif arg == '-smc':
            smc = True
        elif arg == '-out':
            i += 1
            outroot = argv[i]
        elif arg == '-lvmsky_skysub':
            i += 1
            lvmsky_skysub = argv[i]
        elif arg.startswith('-'):
            print('Error: unknown option "%s"' % arg)
            print(_USAGE)
            sys.exit(1)
        else:
            if filename is not None:
                print('Error: unexpected argument "%s" (filename already set to "%s")'
                      % (arg, filename))
                print(_USAGE)
                sys.exit(1)
            filename = arg
        i += 1

    if filename is None:
        print('Error: no filename supplied')
        print(_USAGE)
        sys.exit(1)

    if not os.path.exists(filename):
        print('Error: file not found: %s' % filename)
        sys.exit(1)

    valid_methods = {'nearest', 'farlines_nearcont'}
    if method not in valid_methods:
        print('Error: -method must be one of: %s' % ', '.join(sorted(valid_methods)))
        sys.exit(1)

    vel = resolve_velocity(vel_arg, lmc, smc)

    do_all(filename=filename, method=method, idelta=idelta, vel=vel,
          outroot=outroot, lvmsky_skysub=lvmsky_skysub)
