#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Perform sky subtraction on an XCframe file using the ESO Sky Model
    to separate the sky continuum into its physical MOON, ZODI, and
    DIFFUSE components, then scale sky lines against science lines
    exactly as in SkySubOrig.py.  Three methods are supported:

        nearest          continuum and lines both from the nearest
                         sky telescope
        farthest         continuum and lines both from the farthest
                         sky telescope
        farlines_nearcont  scale lines from the far sky, continuum from
                           the near sky (default)

    Produces a FITS file with WAVE, FLUX (sky-subtracted), SKY, MOON,
    ZODI, DIFFUSE, and DRP_ALL extensions.  The DRP_ALL table carries
    the fitted ESO-model coefficients, the sky-line scale factor, a
    QA_FLAGS column, and (if sky_mask.fits is found) per-arm
    continuum-fit-quality stats against the raw pre-subtraction
    science/sky spectra.

Command line usage (if any):

    usage: SkySepESO.py [-method METHOD] [-delta N] [-out ROOT] filename

    Arguments::

        filename    XCframe FITS file to process

    Options::

        -method METHOD   sky subtraction method: nearest | farthest |
                         farlines_nearcont  (default: farlines_nearcont)
        -delta N         process every N-th row; useful for quick tests
                         (default: 1 = all rows)
        -out ROOT        output filename root; default is <stem>_eso_<method>

Description:

    For each row (fiber):

    1. RA/Dec and obstime for the science fiber and both sky telescopes
       are read from DRP_ALL, and the nearer/farther sky telescope is
       identified from their angular separation from the science fiber
       (as in SkySubOrig.py).

    2. The ESO sky model is fetched for the science fiber's coordinates
       and time (EsoSkyObs.run_sky_obs, engine='auto': local ESO SM-01
       model first, falling back to the SkyCalc web service if that
       fails), giving MOON/ZODI/DIFFUSE continuum templates on the
       model's own wavelength grid.  These are interpolated onto the instrument
       wavelength grid and fit to the observed flux as a non-negative
       3-component linear combination::

           CONT = a*MOON + b*ZODI + c*DIFFUSE   (a, b, c >= 0)

       using an iteratively downweighted, sigma-clipped least-squares
       fit so that sky/airglow line pixels do not pull the continuum
       fit upward.  LINES = FLUX - CONT.

    3. The same model-fetch-and-fit is repeated for whichever sky
       fiber(s) the chosen method requires.

    4. A line scale factor r is found by ksl_bisection (imported from
       SkySubOrig.py) minimising the same objective used by the other
       SkySub* scripts::

           sum |sci_lines * (sci_lines - r*sky_lines)| / ||sky_lines||^2

    5. SKY = sky_CONT + r*sky_LINES; sky-subtracted FLUX = FLUX - SKY.

    MOON/ZODI/DIFFUSE in the output file are the *scaled* components
    (a*MOON, b*ZODI, c*DIFFUSE) of whichever sky fit produced sky_CONT
    (the near-sky fit for nearest/farlines_nearcont, the far-sky fit
    for farthest) -- MOON+ZODI+DIFFUSE reconstructs the continuum part
    of SKY exactly; the remaining line term r*sky_LINES can be
    recovered as SKY - (MOON+ZODI+DIFFUSE).

    DRP_ALL gains these per-row columns::

        LINE_SCALE                     the bisection factor r
        SCI_MOON, SCI_ZODI, SCI_DIFFUSE  science-side fit coefficients
                                        (used only to derive sci_lines
                                        for the bisection; not part of
                                        the subtracted SKY)
        SKY_MOON, SKY_ZODI, SKY_DIFFUSE  coefficients of the sky fit
                                        that produced the MOON/ZODI/
                                        DIFFUSE extensions
        ERROR_MSG                      empty for a good row; otherwise
                                        a description (truncated to 200
                                        characters) of why the row
                                        failed -- see Notes.

    QA flag bits stored in DRP_ALL['QA_FLAGS']::

        0x01  NANDATA    NaN/inf found in input flux or sky data
        0x02  ZEROSKY    sky line vector is all-zero; scale unreliable
        0x04  MODELFAIL  the ESO sky model (local SM-01 and the SkyCalc
                         web-service fallback) could not be fetched for
                         this row
        0x08  FAILED     row raised an exception; spectrum filled with
                         NaN (always set alongside MODELFAIL, since a
                         model-fetch failure also fails the row)

Notes:

    Every row with a non-zero QA_FLAGS value is listed at the end of
    the run with the specific reason it failed (e.g. which engine
    EsoSkyObs.run_sky_obs tried and why each failed); this report is
    printed to the terminal and also written to ``<ROOT>_errors.txt`` next to the
    output FITS file.  The same reason string (truncated to 200
    characters) is stored per row in DRP_ALL['ERROR_MSG'], so the
    output FITS file is self-documenting even without the text file.

    Each row requires two or three ESO sky model fetches (one for the
    science fiber, one or two for the sky fiber(s) needed by -method),
    each of which writes a small per-call FITS file to the current
    working directory; that file is read and deleted immediately, so
    nothing accumulates on disk.  This makes SkySepESO.py far slower
    per row than the other SkySub* scripts -- use -delta for quick
    tests, e.g. -delta 50.

    EsoSkyObs.run_sky_obs's local engine requires the ESO SM-01 sky
    model binary, resolved dynamically from the ESO_SKY_MODEL
    environment variable, and only works on machines where it is set
    up; SkySepESO.py falls back automatically to the SkyCalc web
    service (EsoSkyObs.run_sky_obs's remote engine) if the local call
    fails, so it can still run (more slowly, and needing network
    access) on machines without the local model installed.

    Output filename is ``<ROOT>.fits``.  If -out is omitted the name
    is derived as ``<stem>_eso_<method>.fits``.

    Ported from the lvm_sky2506 prototype (SkySepMod.py's
    separate_continuum/fit_continuum plus SkySubModDev250624.ipynb's
    one_drp).  The MOON/ZODI/DIFFUSE output extensions and the DRP_ALL
    coefficient/scale-factor columns are new; the prototype only wrote
    ASCII tables for one row at a time and did not record fit
    coefficients anywhere.

History:

    260705 ksl  Ported from lvm_sky2506/SkySepMod.py +
                SkySubModDev250624.ipynb.
    260706 ksl  Added per-row failure reason tracking: DRP_ALL['ERROR_MSG'],
                a printed end-of-run summary table, and a companion
                <ROOT>_errors.txt file, so failures (e.g. ESO model
                rejecting a geometry, or a broken SkyCalc fallback) are
                explained, not just flagged.
    260706 ksl  DRP_ALL['mjd'] now recomputed precisely from 'obstime' via
                SkySubOrig.obstime_to_mjd(), instead of the truncated
                integer carried through from the input file.
    260706 ksl  Added per-arm continuum-fit-quality columns (SCI_*/SKY_* for
                med/nmad/rms/skew x b/r/z), evaluated against the raw
                pre-subtraction science and sky spectra using
                GetSkyCont.arm_continuum_stats(); requires sky_mask.fits.
    260711 ksl  Migrated _get_sky_model from SkyModelObs.do_one/SkyCalcObs.py
                (both retired) to EsoSkyObs.run_sky_obs(engine='auto'), which
                already implements the same local-then-remote fallback.
                Picks up a real fix as a side effect: the local engine now
                resolves the historical solar flux (GetSolar.get_flux) instead
                of SkyModelObs.py's old hardcoded msolflux=101.

'''

import sys
import os
import io
import contextlib
from pathlib import Path

# ensure py_progs siblings are importable when running directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
from scipy.optimize import minimize
from scipy.stats import sigmaclip

import EsoSkyObs
from SkySubOrig import fit_func, ksl_bisection, obstime_to_mjd

try:
    from GetSkyCont import (load_mask, _interp_mask_to_wave,
                            arm_continuum_stats, flatten_arm_stats)
    _HAVE_MASK = True
except ImportError:
    _HAVE_MASK = False

# ──────────────────────────────────────────────────────────────
# QA flag bits
# ──────────────────────────────────────────────────────────────
QA_NANDATA   = 1   # NaN or inf values found in input flux/sky data
QA_ZEROSKY   = 2   # sky line vector is all-zero; scale factor is unreliable
QA_MODELFAIL = 4   # ESO sky model (local + SkyCalc fallback) both failed
QA_FAILED    = 8   # row failed entirely; FLUX and SKY are NaN

_QA_FLAG_NAMES = {
    QA_NANDATA:   'NANDATA',
    QA_ZEROSKY:   'ZEROSKY',
    QA_MODELFAIL: 'MODELFAIL',
    QA_FAILED:    'FAILED',
}

_USAGE = '''Usage:
  SkySepESO.py [-method METHOD] [-delta N] [-out ROOT] filename

Arguments:
  filename         XCframe FITS file to process

Options:
  -method METHOD   nearest | farthest | farlines_nearcont
                   (default: farlines_nearcont)
  -delta N         step size through rows for quick tests (default: 1)
  -out ROOT        output filename root (default: <stem>_eso_<method>)
'''


# ──────────────────────────────────────────────────────────────
# ESO sky model fetch
# ──────────────────────────────────────────────────────────────

def _run_captured(func, *args, **kwargs):
    '''
    Call func, capturing anything it prints to stdout, echoing it back
    to the terminal unchanged (so live progress monitoring is
    unaffected), and returning (result, captured_text) so the printed
    diagnostics can also be folded into an error message if needed.
    '''
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        result = func(*args, **kwargs)
    captured = buf.getvalue()
    if captured:
        print(captured, end='' if captured.endswith('\n') else '\n')
    return result, captured.strip()


def _get_sky_model(ra, dec, obstime):
    '''
    Fetch an ESO sky model spectrum for the given coordinates/time.

    Uses EsoSkyObs.run_sky_obs(engine='auto'), which tries the local ESO
    SM-01 model first and falls back to the SkyCalc web service if that
    fails.  The per-call FITS file it writes is deleted immediately after
    being read, so nothing accumulates on disk across a full run.

    Returns (model_tab, err_msg).  model_tab is an astropy Table with
    WAVE, MOON, ZODI, DIFFUSE columns, and err_msg is ''.  If both the
    local model and the web fallback fail, model_tab is None and err_msg
    is whatever run_sky_obs printed while trying each engine.
    '''
    outroot, captured = _run_captured(EsoSkyObs.run_sky_obs,
                                      ra=ra, dec=dec, xtime=obstime, engine='auto')

    if outroot == '':
        err_msg = captured or 'EsoSkyObs.run_sky_obs failed, no message'
        return None, err_msg

    modelfile = '%s.fits' % outroot
    try:
        model_hdul = fits.open(modelfile)
        model_tab = Table(model_hdul[1].data)
        model_hdul.close()
    except Exception as e:
        return None, 'could not read sky model file %s (%s)' % (modelfile, e)
    finally:
        if os.path.exists(modelfile):
            os.remove(modelfile)
    return model_tab, ''


def _fit_eso_continuum(wave, flux, model_tab, max_iter=5, clip_sigma=3.0,
                       downweight_factor=5.0, initial_coeffs=(1.0, 1.0, 1.0)):
    '''
    Fit a non-negative 3-component (MOON, ZODI, DIFFUSE) linear
    combination of ESO sky model templates to an observed spectrum.

    Iteratively downweights positive residuals (sky/airglow lines) and
    sigma-clips outliers so the fit follows the continuum rather than
    being pulled upward by emission lines.  Ported from lvm_sky2506's
    SkySepMod.fit_continuum.

    The fit is done on flux rescaled by 1e17 (order-unity values) for
    numerical conditioning of the optimizer; this does not change the
    best-fit coefficients, only their numerical stability.

    Returns (continuum, coefficients, components):
      continuum   : array, same length as wave/flux
      coefficients: (a, b, c)
      components  : (a*MOON, b*ZODI, c*DIFFUSE), each same length as wave
    '''
    factor = 1e17
    flux_s  = flux * factor
    moon    = np.interp(wave, model_tab['WAVE'], model_tab['MOON'])    * factor
    zodi    = np.interp(wave, model_tab['WAVE'], model_tab['ZODI'])    * factor
    diffuse = np.interp(wave, model_tab['WAVE'], model_tab['DIFFUSE']) * factor

    weights = np.ones_like(flux_s)
    popt = np.array(initial_coeffs, dtype=float)

    def model_func(params):
        a, b, c = params
        return a * moon + b * zodi + c * diffuse

    def cost_function(params):
        residuals = flux_s - model_func(params)
        return np.sum(weights * residuals ** 2)

    bounds = [(0, None)] * 3
    for _ in range(max_iter):
        old_popt = popt.copy()
        result = minimize(cost_function, popt, method='L-BFGS-B', bounds=bounds)
        if not result.success:
            break
        popt = result.x

        residuals = flux_s - model_func(popt)
        clipped = sigmaclip(residuals, low=clip_sigma, high=clip_sigma).clipped

        new_weights = weights.copy()
        positive = residuals > 0
        new_weights[positive] = weights[positive] / downweight_factor
        outlier = ~np.isin(residuals, clipped)
        new_weights[outlier] = weights[outlier] / downweight_factor

        param_change  = (np.sum(np.abs(popt - old_popt))
                         / (np.sum(np.abs(old_popt)) + 1e-10))
        weight_change = (np.sum(np.abs(new_weights - weights))
                         / (np.sum(np.abs(weights)) + 1e-10))
        weights = new_weights
        if param_change < 0.01 and weight_change < 0.01:
            break

    a, b, c = popt
    components = (a * moon / factor, b * zodi / factor, c * diffuse / factor)
    continuum  = components[0] + components[1] + components[2]
    return continuum, (float(a), float(b), float(c)), components


# ──────────────────────────────────────────────────────────────
# Per-row sky subtraction
# ──────────────────────────────────────────────────────────────

def one_drp(xfits, row, method='farlines_nearcont', clean_mask=None):
    '''
    Sky-subtract a single spectrum (row) using the ESO sky model.

    Returns (scitab, qa_flags, coeffs, err_msg).

    scitab has columns WAVE, FLUX, SKY, SCI_FLUX, MOON, ZODI, DIFFUSE
    -- MOON/ZODI/DIFFUSE are the scaled components of the sky
    continuum that sum to the continuum part of SKY.

    coeffs is a dict with sci_moon/sci_zodi/sci_diffuse (science-side
    fit coefficients), sky_moon/sky_zodi/sky_diffuse (coefficients of
    the sky fit that produced MOON/ZODI/DIFFUSE), line_scale (the
    bisection factor r), and -- if clean_mask is given -- per-arm
    continuum-fit-quality stats (GetSkyCont.arm_continuum_stats/
    flatten_arm_stats) evaluated on the RAW pre-subtraction science
    spectrum (sci_lines) and the sky spectrum whose CONT went into the
    final SKY (near for nearest/farlines_nearcont, far for farthest;
    keys 'sci_<stat>_<arm>'/'sky_<stat>_<arm>').  This tests the ESO
    model continuum fit itself, not the final sky-subtracted result.

    err_msg is '' on success, else a description of what failed.

    On error, returns (None, QA_FAILED, {}, err_msg) or
    (None, QA_MODELFAIL, {}, err_msg).
    '''
    qa_flags = 0
    try:
        wave     = np.array(xfits['WAVE'].data, dtype=float)
        flux     = np.array(xfits['FLUX'].data[row],     dtype=float)
        skye     = np.array(xfits['SKY_EAST'].data[row], dtype=float)
        skyw     = np.array(xfits['SKY_WEST'].data[row], dtype=float)
        drp_all  = Table(xfits['DRP_ALL'].data)
        obstime  = drp_all['obstime'][row]
        sci_ra   = drp_all['sci_ra'][row]
        sci_dec  = drp_all['sci_dec'][row]
        skye_ra  = drp_all['skye_ra'][row]
        skye_dec = drp_all['skye_dec'][row]
        skyw_ra  = drp_all['skyw_ra'][row]
        skyw_dec = drp_all['skyw_dec'][row]
    except Exception as e:
        msg = 'could not read row data: %s' % e
        print('Row %d: %s' % (row, msg))
        return None, QA_FAILED, {}, msg

    if not (np.all(np.isfinite(flux)) and np.all(np.isfinite(skye))
            and np.all(np.isfinite(skyw))):
        qa_flags |= QA_NANDATA
        # nan_to_num's default replaces +-inf with +-1.8e308 (float64 max),
        # not 0 -- that "poison" value overflows through the ESO-model
        # continuum fit's matrix operations, corrupting the row.  Zero all
        # non-finite pixels explicitly instead.
        flux = np.nan_to_num(flux, nan=0.0, posinf=0.0, neginf=0.0)
        skye = np.nan_to_num(skye, nan=0.0, posinf=0.0, neginf=0.0)
        skyw = np.nan_to_num(skyw, nan=0.0, posinf=0.0, neginf=0.0)

    sci_coord  = SkyCoord(ra=sci_ra  * u.degree, dec=sci_dec  * u.degree)
    skye_coord = SkyCoord(ra=skye_ra * u.degree, dec=skye_dec * u.degree)
    skyw_coord = SkyCoord(ra=skyw_ra * u.degree, dec=skyw_dec * u.degree)
    de = sci_coord.separation(skye_coord)
    dw = sci_coord.separation(skyw_coord)
    if de < dw:
        sky_near, near_ra, near_dec = skye, skye_ra, skye_dec
        sky_far,  far_ra,  far_dec  = skyw, skyw_ra, skyw_dec
    else:
        sky_near, near_ra, near_dec = skyw, skyw_ra, skyw_dec
        sky_far,  far_ra,  far_dec  = skye, skye_ra, skye_dec

    # science continuum -- used only to derive sci_lines for the
    # bisection scale-factor fit; not part of the subtracted SKY
    sci_model, sci_err = _get_sky_model(sci_ra, sci_dec, obstime)
    if sci_model is None:
        return None, QA_MODELFAIL, {}, 'science fiber: %s' % sci_err
    sci_cont, sci_coeffs, _sci_components = _fit_eso_continuum(wave, flux, sci_model)
    sci_lines = flux - sci_cont

    def _sky_fit(ra, dec, spec):
        model, err = _get_sky_model(ra, dec, obstime)
        if model is None:
            return None, err
        return _fit_eso_continuum(wave, spec, model), ''

    if method == 'nearest':
        fit, err = _sky_fit(near_ra, near_dec, sky_near)
        if fit is None:
            return None, QA_MODELFAIL, {}, 'near sky fiber: %s' % err
        sky_cont, sky_coeffs, sky_components = fit
        sky_lines = sky_near - sky_cont
        sky_cont_resid = sky_lines

    elif method == 'farthest':
        fit, err = _sky_fit(far_ra, far_dec, sky_far)
        if fit is None:
            return None, QA_MODELFAIL, {}, 'far sky fiber: %s' % err
        sky_cont, sky_coeffs, sky_components = fit
        sky_lines = sky_far - sky_cont
        sky_cont_resid = sky_lines

    elif method == 'farlines_nearcont':
        far_fit,  far_err  = _sky_fit(far_ra,  far_dec,  sky_far)
        near_fit, near_err = _sky_fit(near_ra, near_dec, sky_near)
        if far_fit is None or near_fit is None:
            msgs = []
            if far_fit is None:
                msgs.append('far sky fiber: %s' % far_err)
            if near_fit is None:
                msgs.append('near sky fiber: %s' % near_err)
            return None, QA_MODELFAIL, {}, '; '.join(msgs)
        far_cont, _far_coeffs, _far_components = far_fit
        near_cont, sky_coeffs, sky_components  = near_fit
        sky_cont  = near_cont
        sky_lines = sky_far - far_cont
        # CONT used in the final SKY comes from the NEAR fit; evaluate
        # continuum-fit quality against the near fiber's own raw flux,
        # not sky_lines above (which is the far fiber's residual, used
        # only for the line-scale target).
        sky_cont_resid = sky_near - near_cont

    else:
        msg = 'unknown method "%s"' % method
        print('Error: %s' % msg)
        return None, QA_FAILED, {}, msg

    if np.dot(sky_lines, sky_lines) == 0:
        qa_flags |= QA_ZEROSKY

    r = ksl_bisection(fit_func, 0.5, 1.5, tol=0.001, maxiter=8,
                      args=(sci_lines, sky_lines))

    sky      = sky_cont + r * sky_lines
    sci_flux = flux - sky

    scitab = Table(
        [wave, flux, sky, sci_flux,
         sky_components[0], sky_components[1], sky_components[2]],
        names=['WAVE', 'FLUX', 'SKY', 'SCI_FLUX', 'MOON', 'ZODI', 'DIFFUSE'])

    coeffs = dict(
        sci_moon=sci_coeffs[0], sci_zodi=sci_coeffs[1], sci_diffuse=sci_coeffs[2],
        sky_moon=sky_coeffs[0], sky_zodi=sky_coeffs[1], sky_diffuse=sky_coeffs[2],
        line_scale=r,
    )
    if clean_mask is not None:
        coeffs.update(flatten_arm_stats(
            'sci', arm_continuum_stats(wave, sci_lines, clean_mask)))
        coeffs.update(flatten_arm_stats(
            'sky', arm_continuum_stats(wave, sky_cont_resid, clean_mask)))

    return scitab, qa_flags, coeffs, ''


# ──────────────────────────────────────────────────────────────
# Batch processing
# ──────────────────────────────────────────────────────────────

def do_all(filename, method='farlines_nearcont', idelta=1, outroot=''):
    '''
    Process every row of an XCframe file and write sky-subtracted output.

    Parameters
    ----------
    filename : str
        Path to the input XCframe FITS file.
    method : str
        Sky subtraction method (nearest, farthest, farlines_nearcont).
    idelta : int
        Row step size (1 = all rows).
    outroot : str
        Output filename root; defaults to <stem>_eso_<method>.
    '''
    x = fits.open(filename)
    drp_all    = Table(x['DRP_ALL'].data)
    final_wave = np.array(x['WAVE'].data)
    nan_spectrum = np.full(len(final_wave), np.nan)

    # Load sky-line mask once (same convention as SkySub_eval.py) purely for
    # the continuum-fit-quality evaluation below; the ESO-model fit itself
    # is unaffected if the mask is unavailable.
    clean_mask = None
    if _HAVE_MASK:
        _data_dir = Path(__file__).parent.parent / 'data'
        for _candidate in [Path('sky_mask.fits'), _data_dir / 'sky_mask.fits']:
            if _candidate.exists():
                try:
                    _mask_wave, _mask_arr = load_mask(str(_candidate))
                    clean_mask = _interp_mask_to_wave(_mask_wave, _mask_arr, final_wave)
                    print('Loaded sky mask for continuum-quality evaluation: %s' % _candidate)
                except Exception as _e:
                    print('Warning: could not load mask %s (%s)' % (_candidate, _e))
                break
    if clean_mask is None:
        print('Warning: sky_mask.fits not found; skipping continuum-quality columns')

    final_flux    = []
    final_sky     = []
    final_moon    = []
    final_zodi    = []
    final_diffuse = []
    select        = []
    qa_flags_list = []
    coeffs_list   = []
    err_msgs      = []

    i = 0
    while i < len(drp_all):
        err_msg = ''
        try:
            ftab, row_flags, coeffs, err_msg = one_drp(
                xfits=x, row=i, method=method, clean_mask=clean_mask)
        except Exception as e:
            print('Row %d: exception (%s)' % (i, e))
            ftab = None
            row_flags = 0
            coeffs = {}
            err_msg = 'exception: %s' % e

        if ftab is None:
            row_flags |= QA_FAILED
            final_flux.append(nan_spectrum.copy())
            final_sky.append(nan_spectrum.copy())
            final_moon.append(nan_spectrum.copy())
            final_zodi.append(nan_spectrum.copy())
            final_diffuse.append(nan_spectrum.copy())
        else:
            final_flux.append(np.array(ftab['SCI_FLUX']))
            final_sky.append(np.array(ftab['SKY']))
            final_moon.append(np.array(ftab['MOON']))
            final_zodi.append(np.array(ftab['ZODI']))
            final_diffuse.append(np.array(ftab['DIFFUSE']))

        select.append(i)
        qa_flags_list.append(row_flags)
        coeffs_list.append(coeffs)
        err_msgs.append(err_msg)
        i += idelta
        if i % 100 == 0:
            print('Completed %6d of %d in steps of %d' % (i, len(drp_all), idelta))

    n_failed    = sum(1 for f in qa_flags_list if f & QA_FAILED)
    n_modelfail = sum(1 for f in qa_flags_list if f & QA_MODELFAIL)
    n_warned    = sum(1 for f in qa_flags_list
                      if f != 0 and not (f & (QA_FAILED | QA_MODELFAIL)))
    print('\nProcessed %d rows: %d failed (NaN fill), %d model-fetch '
          'failures, %d with warnings'
          % (len(select), n_failed, n_modelfail, n_warned))

    if outroot == '':
        stem = os.path.splitext(os.path.basename(filename))[0]
        outroot = '%s_eso_%s' % (stem, method)

    flagged = [(select[j], qa_flags_list[j], err_msgs[j])
               for j in range(len(select)) if qa_flags_list[j] != 0]
    if flagged:
        report_lines = ['Rows with QA flags:',
                        '  %-8s %-22s %s' % ('Row', 'Flags', 'Reason')]
        for orig_row, flags, msg in flagged:
            active = [name for bit, name in _QA_FLAG_NAMES.items() if flags & bit]
            report_lines.append('  %-8d %-22s %s'
                                % (orig_row, ','.join(active), msg or '(no message)'))
        report = '\n'.join(report_lines)
        print(report)

        errfile = '%s_errors.txt' % outroot
        with open(errfile, 'w') as fh:
            fh.write(report + '\n')
        print('Wrote error summary to %s' % errfile)

    out_wave = np.array(x['WAVE'].data)

    hdu1 = fits.PrimaryHDU(data=None)
    hdu1.header['Title']  = 'SkySepESO'
    hdu1.header['METHOD'] = method
    hdu2 = fits.ImageHDU(data=out_wave,             name='WAVE')
    hdu3 = fits.ImageHDU(data=np.array(final_flux), name='FLUX')
    hdu4 = fits.ImageHDU(data=np.array(final_sky),  name='SKY')
    hdu5 = fits.ImageHDU(data=np.array(final_moon), name='MOON')
    hdu6 = fits.ImageHDU(data=np.array(final_zodi), name='ZODI')
    hdu7 = fits.ImageHDU(data=np.array(final_diffuse), name='DIFFUSE')

    xtab = Table(x['DRP_ALL'].data)
    xtab = xtab[select]
    xtab['QA_FLAGS']    = np.array(qa_flags_list, dtype=np.int32)
    xtab['LINE_SCALE']  = np.array([c.get('line_scale',  np.nan) for c in coeffs_list])
    xtab['SCI_MOON']    = np.array([c.get('sci_moon',    np.nan) for c in coeffs_list])
    xtab['SCI_ZODI']    = np.array([c.get('sci_zodi',    np.nan) for c in coeffs_list])
    xtab['SCI_DIFFUSE'] = np.array([c.get('sci_diffuse', np.nan) for c in coeffs_list])
    xtab['SKY_MOON']    = np.array([c.get('sky_moon',    np.nan) for c in coeffs_list])
    xtab['SKY_ZODI']    = np.array([c.get('sky_zodi',    np.nan) for c in coeffs_list])
    xtab['SKY_DIFFUSE'] = np.array([c.get('sky_diffuse', np.nan) for c in coeffs_list])
    # Continuum-fit-quality columns (raw pre-subtraction sci/sky spectra,
    # not the final sky-subtracted result -- see one_drp() docstring).
    _known_keys = {'line_scale', 'sci_moon', 'sci_zodi', 'sci_diffuse',
                  'sky_moon', 'sky_zodi', 'sky_diffuse'}
    _cont_keys = sorted({k for d in coeffs_list for k in d} - _known_keys)
    for _key in _cont_keys:
        xtab[_key.upper()] = np.array(
            [d.get(_key, np.nan) for d in coeffs_list], dtype=np.float32)
    xtab['ERROR_MSG']   = np.array([m[:200] for m in err_msgs])
    if 'obstime' in xtab.colnames and 'mjd' in xtab.colnames:
        xtab['mjd'] = obstime_to_mjd(xtab['obstime'])
    hdu8 = fits.BinTableHDU(xtab, name='DRP_ALL')

    dwave = out_wave[1] - out_wave[0] if len(out_wave) > 1 else 0.5
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [1, 1]
    wcs.wcs.crval = [float(out_wave[0]), 0]
    wcs.wcs.cdelt = [float(dwave), 1]
    wcs.wcs.ctype = ['WAVE', 'LINE']
    for hdu in (hdu3, hdu4, hdu5, hdu6, hdu7):
        hdu.header.update(wcs.to_header())

    hdul = fits.HDUList([hdu1, hdu2, hdu3, hdu4, hdu5, hdu6, hdu7, hdu8])

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

    method  = 'farlines_nearcont'
    idelta  = 1
    outroot = ''
    filename = None

    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == '-method':
            i += 1
            method = argv[i]
        elif arg == '-delta':
            i += 1
            idelta = int(argv[i])
        elif arg == '-out':
            i += 1
            outroot = argv[i]
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

    valid_methods = {'nearest', 'farthest', 'farlines_nearcont'}
    if method not in valid_methods:
        print('Error: -method must be one of: %s' % ', '.join(sorted(valid_methods)))
        sys.exit(1)

    do_all(filename=filename, method=method, idelta=idelta, outroot=outroot)
