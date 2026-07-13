#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Evaluate how well the ESO Sky Model reproduces REAL observed sky
    spectra -- a model-quality tool, not a sky-subtraction tool (compare
    SkySepESO.py, which subtracts the model from a science spectrum).
    For each observed sky spectrum, fits the ESO Sky Model's MOON, ZODI,
    DIFFUSE, and LINES templates -- all four, unlike SkySepESO.py, which
    only ever fits the three continuum templates -- against the real
    flux as a single 4-parameter non-negative linear combination, so the
    full model (continuum AND lines) can be judged against the real sky
    at once.

    Produces a FITS file that is the input file, extended: WAVE and FLUX
    are passed through unchanged; MOON, ZODI, DIFFUSE, LINES (the four
    scaled fitted model components) and RESID (FLUX minus their sum --
    the actual discrepancy between the real spectrum and the model's
    best 4-parameter fit) are added; DRP_ALL gains the fitted scaling
    factors, a QA_FLAGS column, and ERROR_MSG.

Command line usage (if any):

    usage: EsoSkyFit.py [-engine local|remote|auto] [-delta N] [-out ROOT] filename

    Arguments::

        filename    Sky_<name>.fits produced by GetSky_from_CFrame_sum.py's
                    extraction mode (WAVE, FLUX, DRP_ALL with ra/dec/obstime)

    Options::

        -engine E        'local' (default) | 'remote' | 'auto' -- which
                         EsoSkyObs.run_sky_obs engine to fetch the ESO
                         model from; see Notes for why 'local' is the
                         default here (unlike SkySepESO.py's 'auto')
        -delta N         process every N-th row; useful for quick tests
                         (default: 1 = all rows)
        -out ROOT        output filename root; default is <stem>_esofit

Description:

    For each row (sky observation):

    1. ra/dec/obstime are read from DRP_ALL (as merged by
       GetSky_from_CFrame_sum.py -- one physical sky pointing/time per
       row, not a science/sky pair).

    2. The ESO sky model is fetched for that position/time
       (SkySepESO._get_sky_model, reused unchanged; engine='local' by
       default here), giving MOON/ZODI/DIFFUSE/LINES templates on the
       model's own wavelength grid -- LINES here is the model's own
       predicted airglow emission-line spectrum, a genuine column
       SkySepESO._get_sky_model already returns but SkySepESO.py itself
       never uses (see Notes for why).

    3. _fit_full_model (this module's own function, NOT reused from
       SkySepESO.py -- see Notes) fits all four templates to the
       observed flux as a single non-negative 4-component linear
       combination::

           MODEL = a*MOON + b*ZODI + c*DIFFUSE + d*LINES   (a,b,c,d >= 0)

       using the same iteratively downweighted, sigma-clipped
       least-squares scheme SkySepESO._fit_eso_continuum uses (see that
       function and _fit_full_model's own docstring for why downweighting
       outliers during the fit doesn't suppress the final reported
       residual).  RESID = FLUX - MODEL is the real discrepancy between
       the observed spectrum and the model's best 4-parameter fit --
       continuum and lines together.

    Unlike SkySepESO.py, there is no science fiber, no near/far sky
    telescope selection, no line-scale bisection, and no subtraction --
    this is purely fitting the model's full predicted spectrum against a
    real sky spectrum, one row at a time.

    DRP_ALL gains these per-row columns::

        MOON_COEFF, ZODI_COEFF, DIFFUSE_COEFF, LINES_COEFF
                                                 the fitted a, b, c, d
        ERROR_MSG                               empty for a good row;
                                                 otherwise a description
                                                 (truncated to 200
                                                 characters) of why the
                                                 row failed

    QA flag bits stored in DRP_ALL['QA_FLAGS'] (a subset of
    SkySepESO.py's scheme -- no ZEROSKY, since there is no line-scale
    bisection here)::

        0x01  NANDATA    NaN/inf found in the input flux
        0x02  MODELFAIL  the ESO sky model could not be fetched for this
                         row with the requested -engine
        0x04  FAILED     row raised an exception; spectrum filled with
                         NaN (always set alongside MODELFAIL, since a
                         model-fetch failure also fails the row)

Notes:

    -engine defaults to 'local' here, unlike SkySepESO.py's 'auto'
    default: since this tool's whole purpose is judging ESO model
    quality, silently letting some rows fall back to the remote SkyCalc
    web service (a different vintage/resolution engine) while others use
    the local SM-01 model would contaminate that judgement.  A row whose
    local model call fails is flagged MODELFAIL/NaN-filled rather than
    silently substituting the remote engine.  Pass -engine auto (or
    remote) explicitly if that tradeoff is what you want instead.

    Moon-below-horizon geometries are a legitimate ESO Sky Model input,
    not an error: EsoSkyObs.py's local-engine config template documents
    altmoon's valid range as [-90,90] degrees, so a below-horizon Moon
    should simply produce a near-zero MOON template, which
    _fit_full_model handles without issue (the fitted a coefficient
    becomes numerically degenerate in that case, since a*0 contributes
    nothing to the residual regardless of a's value, but this has no
    effect on the fit quality or on zodi/diffuse/lines).

    Each row requires exactly one ESO sky model fetch (compare
    SkySepESO.py's two or three per row), which writes a small per-call
    FITS file to the current working directory; that file is read and
    deleted immediately, so nothing accumulates on disk.

    Output filename is ``<ROOT>.fits``.  If -out is omitted the name is
    derived as ``<stem>_esofit.fits``.

    Reuses SkySepESO._get_sky_model directly rather than duplicating it
    -- importing an underscore-prefixed helper across scripts has
    precedent elsewhere in this codebase (SkySepPalace.py imports
    _load_ref_lsf/_interp_lsf_to_wave from SkySubDev2.py; SkySubDev2.py
    imports _get_decomposer from XSkySepIvan.py).  The fit itself
    (_fit_full_model) is NOT reused from SkySepESO._fit_eso_continuum,
    deliberately: SkySepESO.py's real sky subtraction only ever fits the
    three continuum templates and intentionally does not trust the
    model's own absolute line predictions (instead scaling real
    sky-fiber line residuals against real science-fiber lines via
    bisection) -- a real, considered design choice for subtraction, not
    an oversight.  EsoSkyFit.py's purpose is different: judging the full
    model (including its own line predictions) against real data, so it
    needed its own 4-parameter version rather than changing the function
    SkySepESO.py depends on for production sky subtraction.

History:

    260711 ksl Coding begun.  First cut fit only the three continuum
        templates (MOON/ZODI/DIFFUSE) via SkySepESO._fit_eso_continuum,
        reused unchanged, with a LINES column defined as FLUX minus the
        fitted continuum -- realized this was wrong: since LINES was
        just whatever the real data minus fitted continuum, MOON+ZODI+
        DIFFUSE+LINES was a tautology (always equals FLUX exactly, telling
        you nothing about model quality), and LINES wasn't actually a
        model prediction at all, just an empirical leftover conflating
        real airglow lines with model error.  Corrected: fit all four
        ESO model templates (MOON/ZODI/DIFFUSE, and LINES -- the model's
        own predicted airglow line spectrum, a real column
        SkySepESO._get_sky_model already returns) jointly via a new
        _fit_full_model (4-parameter version of _fit_eso_continuum,
        added here rather than modifying SkySepESO.py -- see Notes), and
        added a genuine RESID = FLUX - (MOON+ZODI+DIFFUSE+LINES) column,
        which is now the real discrepancy between the model and the data.

'''

import sys
import os

# ensure py_progs siblings are importable when running directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import minimize
from scipy.stats import sigmaclip

from SkySepESO import _get_sky_model

# ──────────────────────────────────────────────────────────────
# QA flag bits
# ──────────────────────────────────────────────────────────────
QA_NANDATA   = 1   # NaN or inf values found in input flux
QA_MODELFAIL = 2   # ESO sky model could not be fetched for this row
QA_FAILED    = 4   # row failed entirely; spectrum filled with NaN

_QA_FLAG_NAMES = {
    QA_NANDATA:   'NANDATA',
    QA_MODELFAIL: 'MODELFAIL',
    QA_FAILED:    'FAILED',
}

_USAGE = '''Usage:
  EsoSkyFit.py [-engine local|remote|auto] [-delta N] [-out ROOT] filename

Arguments:
  filename         Sky_<name>.fits from GetSky_from_CFrame_sum.py

Options:
  -engine E        local (default) | remote | auto
  -delta N         step size through rows for quick tests (default: 1)
  -out ROOT        output filename root (default: <stem>_esofit)
'''


# ──────────────────────────────────────────────────────────────
# Full 4-component model fit (MOON, ZODI, DIFFUSE, LINES)
# ──────────────────────────────────────────────────────────────

def _fit_full_model(wave, flux, model_tab, max_iter=5, clip_sigma=3.0,
                    downweight_factor=5.0, initial_coeffs=(1.0, 1.0, 1.0, 1.0)):
    '''
    Fit a non-negative 4-component (MOON, ZODI, DIFFUSE, LINES) linear
    combination of ESO sky model templates to an observed spectrum.

    Unlike SkySepESO._fit_eso_continuum (which fits only the three
    continuum templates -- SkySepESO.py's real subtraction deliberately
    does not trust the model's own absolute line predictions, instead
    scaling real sky-fiber line residuals against real science-fiber
    lines), this fits the model's own predicted LINES column too, as a
    4th free amplitude -- since EsoSkyFit.py's purpose is judging how
    well the model as a whole (continuum AND lines) reproduces a real
    observed sky spectrum, not subtraction.

    Same iteratively downweighted, sigma-clipped least-squares scheme as
    _fit_eso_continuum (see that function's docstring): this keeps a few
    outlier pixels (cosmic rays, bad pixels, wavelength-mismatched line
    cores) from distorting the fitted amplitudes, without artificially
    shrinking the final reported residual, which is recomputed from the
    real data at every pixel regardless of what weight it had during the
    fit -- exactly the discrepancy this tool is meant to expose.

    Returns (model_total, coefficients, components):
      model_total : array, same length as wave/flux -- a*MOON+b*ZODI+c*DIFFUSE+d*LINES
      coefficients: (a, b, c, d)
      components  : (a*MOON, b*ZODI, c*DIFFUSE, d*LINES), each same length as wave
    '''
    factor = 1e17
    flux_s  = flux * factor
    moon    = np.interp(wave, model_tab['WAVE'], model_tab['MOON'])    * factor
    zodi    = np.interp(wave, model_tab['WAVE'], model_tab['ZODI'])    * factor
    diffuse = np.interp(wave, model_tab['WAVE'], model_tab['DIFFUSE']) * factor
    lines   = np.interp(wave, model_tab['WAVE'], model_tab['LINES'])   * factor

    weights = np.ones_like(flux_s)
    popt = np.array(initial_coeffs, dtype=float)

    def model_func(params):
        a, b, c, d = params
        return a * moon + b * zodi + c * diffuse + d * lines

    def cost_function(params):
        residuals = flux_s - model_func(params)
        return np.sum(weights * residuals ** 2)

    bounds = [(0, None)] * 4
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

    a, b, c, d = popt
    components = (a * moon / factor, b * zodi / factor,
                 c * diffuse / factor, d * lines / factor)
    model_total = components[0] + components[1] + components[2] + components[3]
    return model_total, (float(a), float(b), float(c), float(d)), components


# ──────────────────────────────────────────────────────────────
# Per-row model fit
# ──────────────────────────────────────────────────────────────

def one_row(xfits, row, engine='local'):
    '''
    Fit the ESO sky model (continuum + lines) to a single observed sky
    spectrum (row).

    Returns (result, qa_flags, coeffs, err_msg).

    result is a dict with keys MOON, ZODI, DIFFUSE, LINES (the four
    scaled fitted model components) and RESID (FLUX minus their sum --
    the discrepancy between the real spectrum and the model's best
    4-parameter fit), each an array the same length as WAVE.

    coeffs is a dict with moon/zodi/diffuse/lines (the fitted a, b, c, d
    scaling factors).

    err_msg is '' on success, else a description of what failed.

    On error, returns (None, QA_FAILED, {}, err_msg) or
    (None, QA_MODELFAIL, {}, err_msg).
    '''
    qa_flags = 0
    try:
        wave    = np.array(xfits['WAVE'].data, dtype=float)
        flux    = np.array(xfits['FLUX'].data[row], dtype=float)
        drp_all = Table(xfits['DRP_ALL'].data)
        ra      = drp_all['ra'][row]
        dec     = drp_all['dec'][row]
        obstime = drp_all['obstime'][row]
    except Exception as e:
        msg = 'could not read row data: %s' % e
        print('Row %d: %s' % (row, msg))
        return None, QA_FAILED, {}, msg

    if not np.all(np.isfinite(flux)):
        qa_flags |= QA_NANDATA
        # see SkySepESO.one_drp for why this is zeroed explicitly rather
        # than left to nan_to_num's default (+-inf -> +-1.8e308, which
        # would overflow through the fit's matrix operations)
        flux = np.nan_to_num(flux, nan=0.0, posinf=0.0, neginf=0.0)

    model_tab, err = _get_sky_model(ra, dec, obstime, engine=engine)
    if model_tab is None:
        return None, QA_MODELFAIL, {}, err

    model_total, (a, b, c, d), (moon_c, zodi_c, diffuse_c, lines_c) = _fit_full_model(
        wave, flux, model_tab)
    resid = flux - model_total

    result = dict(MOON=moon_c, ZODI=zodi_c, DIFFUSE=diffuse_c, LINES=lines_c, RESID=resid)
    coeffs = dict(moon=a, zodi=b, diffuse=c, lines=d)

    return result, qa_flags, coeffs, ''


# ──────────────────────────────────────────────────────────────
# Batch processing
# ──────────────────────────────────────────────────────────────

def do_all(filename, engine='local', idelta=1, outroot=''):
    '''
    Fit the ESO sky model's continuum to every row of a Sky_<name>.fits
    file (from GetSky_from_CFrame_sum.py) and write an extended copy.

    Parameters
    ----------
    filename : str
        Path to the input Sky_<name>.fits file.
    engine : str
        'local' (default), 'remote', or 'auto' -- passed to
        SkySepESO._get_sky_model / EsoSkyObs.run_sky_obs.
    idelta : int
        Row step size (1 = all rows).
    outroot : str
        Output filename root; defaults to <stem>_esofit.
    '''
    x = fits.open(filename)
    drp_all      = Table(x['DRP_ALL'].data)
    final_wave   = np.array(x['WAVE'].data)
    nan_spectrum = np.full(len(final_wave), np.nan)

    final_moon    = []
    final_zodi    = []
    final_diffuse = []
    final_lines   = []
    final_resid   = []
    select        = []
    qa_flags_list = []
    coeffs_list   = []
    err_msgs      = []

    i = 0
    while i < len(drp_all):
        try:
            result, row_flags, coeffs, err_msg = one_row(
                xfits=x, row=i, engine=engine)
        except Exception as e:
            print('Row %d: exception (%s)' % (i, e))
            result = None
            row_flags = 0
            coeffs = {}
            err_msg = 'exception: %s' % e

        if result is None:
            row_flags |= QA_FAILED
            final_moon.append(nan_spectrum.copy())
            final_zodi.append(nan_spectrum.copy())
            final_diffuse.append(nan_spectrum.copy())
            final_lines.append(nan_spectrum.copy())
            final_resid.append(nan_spectrum.copy())
        else:
            final_moon.append(result['MOON'])
            final_zodi.append(result['ZODI'])
            final_diffuse.append(result['DIFFUSE'])
            final_lines.append(result['LINES'])
            final_resid.append(result['RESID'])

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
        outroot = '%s_esofit' % stem

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

    # Output = the input file, extended: WAVE/FLUX/primary header keys
    # pass through unchanged; MOON/ZODI/DIFFUSE/LINES/RESID and new
    # DRP_ALL columns are added.
    header = x['PRIMARY'].header.copy()
    header['Title']  = 'EsoSkyFit'
    header['ENGINE'] = engine

    hdu1 = fits.PrimaryHDU(data=None, header=header)
    hdu2 = fits.ImageHDU(data=final_wave,                     name='WAVE')
    hdu3 = fits.ImageHDU(data=np.array(x['FLUX'].data)[select], name='FLUX')
    hdu4 = fits.ImageHDU(data=np.array(final_moon),    name='MOON')
    hdu5 = fits.ImageHDU(data=np.array(final_zodi),    name='ZODI')
    hdu6 = fits.ImageHDU(data=np.array(final_diffuse), name='DIFFUSE')
    hdu7 = fits.ImageHDU(data=np.array(final_lines),   name='LINES')
    hdu9 = fits.ImageHDU(data=np.array(final_resid),   name='RESID')

    xtab = drp_all[select]
    xtab['QA_FLAGS']      = np.array(qa_flags_list, dtype=np.int32)
    xtab['MOON_COEFF']    = np.array([c.get('moon',    np.nan) for c in coeffs_list])
    xtab['ZODI_COEFF']    = np.array([c.get('zodi',    np.nan) for c in coeffs_list])
    xtab['DIFFUSE_COEFF'] = np.array([c.get('diffuse', np.nan) for c in coeffs_list])
    xtab['LINES_COEFF']   = np.array([c.get('lines',   np.nan) for c in coeffs_list])
    xtab['ERROR_MSG']     = np.array([m[:200] for m in err_msgs])
    hdu8 = fits.BinTableHDU(xtab, name='DRP_ALL')

    hdul = fits.HDUList([hdu1, hdu2, hdu3, hdu4, hdu5, hdu6, hdu7, hdu9, hdu8])

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

    engine   = 'local'
    idelta   = 1
    outroot  = ''
    filename = None

    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == '-engine':
            i += 1
            engine = argv[i]
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

    valid_engines = {'local', 'remote', 'auto'}
    if engine not in valid_engines:
        print('Error: -engine must be one of: %s' % ', '.join(sorted(valid_engines)))
        sys.exit(1)

    do_all(filename=filename, engine=engine, idelta=idelta, outroot=outroot)
