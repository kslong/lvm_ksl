#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Perform alternative sky subtraction on an XCframe file by separating
    continuum and sky lines, then fitting a scale factor for the lines
    independently.  Three methods are supported::

        nearest          subtract the nearest sky telescope spectrum
        farthest         subtract the farthest sky telescope spectrum
        farlines_nearcont  scale lines from the far sky, continuum from the
                           near sky (default)

    Produces a FITS file with WAVE, FLUX (sky-subtracted), SKY, and DRP_ALL
    extensions.  The DRP_ALL table carries a QA_FLAGS column that records
    per-row quality issues, and a LINE_SCALE column recording the
    bisection factor r used to scale the sky lines (see Description).

Command line usage (if any):

    usage: SkySubOrig.py [-method METHOD] [-delta N] [-out ROOT] filename

    Arguments::

        filename    XCframe FITS file to process

    Options::

        -method METHOD   sky subtraction method: nearest | farthest |
                         farlines_nearcont  (default: farlines_nearcont)
        -delta N         process every N-th row; useful for quick tests
                         (default: 1 = all rows)
        -out ROOT        output filename root; default is <stem>_orig_<method>

Description:

    Each science spectrum is decomposed into a continuum (fitted with a
    degree-4 polynomial and sigma clipping) and line residuals.  The
    nearest and farthest sky telescopes are identified from the RA/Dec
    separation stored in DRP_ALL.  A scale factor for the sky lines is
    found by minimising::

        sum | sci_lines * (sci_lines - r * sky_lines) | / ||sky_lines||^2

    using a custom 4-point bisection search.  Only the sky lines are
    scaled by r; the continuum (CONT) is used exactly as fitted, with
    no additional scaling.  SKY = CONT + r*LINES, so LINE_SCALE (r) is
    the one factor needed to recover the sky-line contribution.

    Since the continuum itself is never rescaled, its fit quality
    against the RAW (pre-subtraction) science and sky spectra is
    evaluated directly and recorded per row: if sky_mask.fits is found
    (same convention as SkySub_eval.py), GetSkyCont.arm_continuum_stats
    computes, per spectrograph arm (B/R/Z), the median/NMAD/RMS/skew of
    FLUX-CONT in clean (sky-line-free) pixels -- for the science
    spectrum's own fit (SCI_MED_B etc.) and for whichever sky fiber's
    fit produced the CONT used in SKY (SKY_MED_B etc.).  This is
    unrelated to the final sky-subtracted SCI_FLUX -- see
    SkySub_eval.py's Figures 5/6 for that (post-subtraction) residual,
    which measures leftover source continuum rather than fit quality.

    QA flag bits stored in DRP_ALL['QA_FLAGS']::

        0x01  NANDATA   NaN/inf found in input flux or sky data
        0x02  ZEROSKY   sky line vector is all-zero; scale unreliable
        0x04  POORFIT   continuum polyfit poorly conditioned
        0x08  FAILED    row raised an exception; spectrum filled with NaN

Notes:

    Output filename is ``<ROOT>.fits``.  If -out is omitted the name
    is derived as ``<stem>_orig_<method>.fits`` where ``<stem>`` is the
    input filename without extension (e.g.
    ``XCframe_1.2.1_7325_48860_1_50_orig_farlines_nearcont.fits``).

    The RankWarning emitted by numpy.polyfit is caught portably across
    NumPy 1.x and 2.x.

    Continuum-quality columns (SCI_*/SKY_* per arm) require sky_mask.fits;
    searched for in the current directory, then in the lvm_ksl data/
    directory.  If not found, do_all() prints a warning and those columns
    are omitted (everything else is unaffected).

History::

    250604 ksl Notebook SkySubDev250604.ipynb — initial development
    260630 ksl Converted to executable script; fixed np.RankWarning for
               NumPy 2.x; removed debug prints from ksl_bisection
    260630 ksl Renamed SkySubDev → SkySubOrig; options may precede filename;
               default output name is <stem>_<method>.fits
    260630 ksl Default wavelength range changed to full spectrum (wmin/wmax
               now None by default in one_drp and do_all)
    260706 ksl Added obstime_to_mjd(); DRP_ALL['mjd'] is now recomputed
               precisely from 'obstime' instead of carrying through the
               truncated-integer mjd from the input file.
    260706 ksl one_drp() now also returns the bisection line-scale factor;
               do_all() records it as DRP_ALL['LINE_SCALE'].
    260706 ksl Added per-arm continuum-fit-quality columns (SCI_*/SKY_* for
               med/nmad/rms/skew x b/r/z), evaluated against the raw
               pre-subtraction science and sky spectra using
               GetSkyCont.arm_continuum_stats(); requires sky_mask.fits.
    260708 ksl Default output name is now <stem>_orig_<method>.fits (was
               <stem>_<method>.fits) -- brings it in line with
               SkySubDrp/Dev1/Dev2/SepESO, which all include their own
               method tag (_drp_/_dev1_/_dev2_/_eso_) in the default name.

'''

import sys
import os
import warnings as _warnings
from pathlib import Path

# ensure py_progs siblings are importable when running directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clip
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

try:
    from GetSkyCont import (load_mask, _interp_mask_to_wave,
                            arm_continuum_stats, flatten_arm_stats)
    _HAVE_MASK = True
except ImportError:
    _HAVE_MASK = False

# ──────────────────────────────────────────────────────────────
# NumPy 2.0 compatibility: RankWarning moved to numpy.exceptions
# ──────────────────────────────────────────────────────────────
try:
    _RankWarning = np.exceptions.RankWarning
except AttributeError:
    _RankWarning = np.RankWarning  # NumPy < 2.0

# ──────────────────────────────────────────────────────────────
# QA flag bits
# ──────────────────────────────────────────────────────────────
QA_NANDATA = 1   # NaN or inf values found in input flux/sky data
QA_ZEROSKY = 2   # sky line vector is all-zero; scale factor is unreliable
QA_POORFIT = 4   # polyfit poorly conditioned (too few unclipped points)
QA_FAILED  = 8   # row failed entirely; FLUX and SKY are NaN for this row

_QA_FLAG_NAMES = {
    QA_NANDATA: 'NANDATA',
    QA_ZEROSKY: 'ZEROSKY',
    QA_POORFIT: 'POORFIT',
    QA_FAILED:  'FAILED',
}

_USAGE = '''Usage:
  SkySubOrig.py [-method METHOD] [-delta N] [-out ROOT] filename

Arguments:
  filename         XCframe FITS file to process

Options:
  -method METHOD   nearest | farthest | farlines_nearcont
                   (default: farlines_nearcont)
  -delta N         step size through rows for quick tests (default: 1)
  -out ROOT        output filename root (default: <stem>_orig_<method>)
'''


# ──────────────────────────────────────────────────────────────
# Core functions
# ──────────────────────────────────────────────────────────────

def obstime_to_mjd(obstime):
    '''
    Precise MJD from OBSTIME (ISO-format timestamp string, or an array
    of them).

    DRP_ALL tables carry a truncated-integer MJD column alongside
    OBSTIME; that integer loses sub-day precision (all exposures on
    the same night share one value).  This recomputes MJD from OBSTIME
    directly via astropy so it can replace that column exactly.
    '''
    return Time(np.asarray(obstime, dtype=str), format='isot', scale='utc').mjd


def fit_func(r, sci, sky):
    '''
    Objective function minimised to find the sky scale factor.

    Computes  sum|sci * (sci - r*sky)| / ||sky||^2 , which weights
    bright sky-line regions more heavily.
    '''
    delta = sci - r * sky
    xresult = np.sum(np.abs(sci * delta))
    xnorm = np.dot(sky, sky)
    if xnorm == 0:
        return np.inf
    return xresult / xnorm


def ksl_bisection(func, a, b, tol=1e-2, args=(), maxiter=30):
    '''
    4-point bisection search for the minimum of func on [a, b].

    Returns the x-coordinate of the estimated minimum.
    '''
    for _ in range(maxiter):
        interval = [a, (3*a + b) / 4, (a + b) / 2, (a + 3*b) / 4, b]
        values = [func(x, *args) for x in interval]
        idx = values.index(min(values))
        if idx <= 1:
            a, b = interval[0], interval[2]
        elif idx == 2:
            a, b = interval[1], interval[3]
        else:
            a, b = interval[2], interval[4]
        if abs(a - b) < tol:
            break
    return (a + b) / 2


def polynomial_fit_with_outliers(spectrum_table, degree=3,
                                 sigma_lower=3, sigma_upper=3,
                                 grow=0, max_iter=10):
    '''
    Iterative sigma-clipping polynomial fit to a WAVE/FLUX table.

    Returns (output_table, qa_flags).  output_table gains CONT and MASK
    columns.  qa_flags ORs QA_NANDATA and/or QA_POORFIT as needed.
    '''
    output_table = spectrum_table.copy()
    x = spectrum_table['WAVE']
    x_centered = x - np.mean(x)
    qa_flags = 0

    flux = np.array(spectrum_table['FLUX'], dtype=float)
    finite = np.isfinite(flux)
    if not np.all(finite):
        qa_flags |= QA_NANDATA

    if np.sum(finite) <= degree:
        output_table['CONT'] = np.zeros(len(spectrum_table))
        output_table['MASK'] = ~finite
        return output_table, qa_flags

    coefficients = np.polyfit(x_centered[finite], flux[finite], degree)
    fitted_flux = np.polyval(coefficients, x_centered)

    mask = ~finite
    for _ in range(max_iter):
        residuals = flux - fitted_flux
        finite_res = np.isfinite(residuals)
        residuals_safe = np.where(finite_res, residuals, 0.0)

        with _warnings.catch_warnings():
            _warnings.simplefilter('ignore')
            clipped = sigma_clip(residuals_safe,
                                 sigma_lower=sigma_lower,
                                 sigma_upper=sigma_upper,
                                 grow=grow)
        mask = clipped.mask | ~finite_res

        n_good = int(np.sum(~mask))
        if n_good <= degree:
            qa_flags |= QA_POORFIT
            break

        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter('always')
            coefficients = np.polyfit(x_centered[~mask], flux[~mask], degree)
        if any(issubclass(w.category, _RankWarning) for w in caught):
            qa_flags |= QA_POORFIT
        fitted_flux = np.polyval(coefficients, x_centered)

    output_table['CONT'] = fitted_flux
    output_table['MASK'] = mask
    return output_table, qa_flags


def one_drp(xfits, row=300, wmin=None, wmax=None,
            method='farlines_nearcont', do_plot=False, clean_mask=None):
    '''
    Sky-subtract a single spectrum (row) from an open FITS object.

    Returns (scitab, qa_flags, line_scale, cont_stats).  scitab has
    columns WAVE, FLUX, CONT, MASK, LINES, SKY, SCI_FLUX.  line_scale
    is the ksl_bisection factor r applied to the sky lines (SKY = CONT
    + r*LINES; CONT itself is not scaled).

    cont_stats is a flat dict of per-arm continuum-fit-quality stats
    (from GetSkyCont.arm_continuum_stats/flatten_arm_stats), evaluated
    on the RAW (pre-subtraction) science spectrum and the raw sky
    spectrum whose CONT went into the final SKY -- i.e. this tests the
    continuum fit itself, not the final sky-subtracted result.  Keys
    are 'sci_<stat>_<arm>' and 'sky_<stat>_<arm>' for stat in
    med/nmad/rms/skew and arm in b/r/z.  Empty dict if clean_mask is
    None or GetSkyCont is unavailable.

    clean_mask : 1-D bool array (n_wave,) or None
        True = sky-line-free pixel (from sky_mask.fits), used only for
        cont_stats; the fit itself is unaffected.

    Returns (None, QA_FAILED, nan, {}) on error.
    '''
    qa_flags = 0
    try:
        wave  = xfits['WAVE'].data
        flux  = xfits['FLUX'].data[row]
        skye  = xfits['SKY_EAST'].data[row]
        skyw  = xfits['SKY_WEST'].data[row]
        drp_all  = Table(xfits['DRP_ALL'].data)
        sci_ra   = drp_all['sci_ra'][row]
        sci_dec  = drp_all['sci_dec'][row]
        skye_ra  = drp_all['skye_ra'][row]
        skye_dec = drp_all['skye_dec'][row]
        skyw_ra  = drp_all['skyw_ra'][row]
        skyw_dec = drp_all['skyw_dec'][row]
    except Exception as e:
        print('Row %d: could not read data (%s)' % (row, e))
        return None, QA_FAILED, np.nan, {}

    sci_coord  = SkyCoord(ra=sci_ra  * u.degree, dec=sci_dec  * u.degree)
    skye_coord = SkyCoord(ra=skye_ra * u.degree, dec=skye_dec * u.degree)
    skyw_coord = SkyCoord(ra=skyw_ra * u.degree, dec=skyw_dec * u.degree)
    de = sci_coord.separation(skye_coord)
    dw = sci_coord.separation(skyw_coord)
    if de < dw:
        sky_near, sky_far = skye, skyw
    else:
        sky_near, sky_far = skyw, skye

    # Select wavelength range (trim only if limits are given)
    if wmin is not None or wmax is not None:
        lo = wmin if wmin is not None else -np.inf
        hi = wmax if wmax is not None else  np.inf
        wcut     = (wave > lo) & (wave < hi)
        wave     = wave[wcut]
        flux     = flux[wcut]
        sky_near = sky_near[wcut]
        sky_far  = sky_far[wcut]
        if clean_mask is not None:
            clean_mask = clean_mask[wcut]

    from astropy.table import Table as _Table
    scitab = _Table([wave, flux], names=['WAVE', 'FLUX'])
    scitab, f = polynomial_fit_with_outliers(scitab, degree=4,
                                             sigma_lower=3, sigma_upper=1, grow=5)
    qa_flags |= f

    if method == 'nearest':
        skytab = _Table([wave, sky_near], names=['WAVE', 'FLUX'])
        skytab, f = polynomial_fit_with_outliers(skytab, degree=4,
                                                 sigma_lower=3, sigma_upper=1, grow=5)
        qa_flags |= f
        sky_cont_resid = np.asarray(skytab['FLUX'] - skytab['CONT'])

    elif method == 'farthest':
        skytab = _Table([wave, sky_far], names=['WAVE', 'FLUX'])
        skytab, f = polynomial_fit_with_outliers(skytab, degree=4,
                                                 sigma_lower=3, sigma_upper=1, grow=5)
        qa_flags |= f
        sky_cont_resid = np.asarray(skytab['FLUX'] - skytab['CONT'])

    elif method == 'farlines_nearcont':
        skytab_far = _Table([wave, sky_far], names=['WAVE', 'FLUX'])
        xfar, f = polynomial_fit_with_outliers(skytab_far, degree=4,
                                               sigma_lower=3, sigma_upper=1, grow=5)
        qa_flags |= f
        skytab_near = _Table([wave, sky_near], names=['WAVE', 'FLUX'])
        xnear, f = polynomial_fit_with_outliers(skytab_near, degree=4,
                                                sigma_lower=3, sigma_upper=1, grow=5)
        qa_flags |= f
        combined_flux = xfar['FLUX'] - xfar['CONT'] + xnear['CONT']
        skytab = _Table([wave, combined_flux, xnear['CONT']],
                        names=['WAVE', 'FLUX', 'CONT'])
        # CONT used in the final SKY comes from the NEAR fit; evaluate
        # continuum-fit quality against the near fiber's own raw flux,
        # not skytab['LINES'] below (which ends up being the far fiber's
        # residual, used only for the line-scale target).
        sky_cont_resid = np.asarray(xnear['FLUX'] - xnear['CONT'])

    else:
        print('Error: unknown method "%s"' % method)
        return None, QA_FAILED, np.nan, {}

    scitab['LINES'] = scitab['FLUX'] - scitab['CONT']
    if 'LINES' not in skytab.colnames:
        skytab['LINES'] = skytab['FLUX'] - skytab['CONT']

    if np.dot(skytab['LINES'], skytab['LINES']) == 0:
        qa_flags |= QA_ZEROSKY

    minimum = ksl_bisection(fit_func, 0.5, 1.5, tol=0.001, maxiter=8,
                            args=(scitab['LINES'], skytab['LINES']))

    scitab['SKY']      = skytab['CONT'] + minimum * skytab['LINES']
    scitab['SCI_FLUX'] = scitab['FLUX'] - scitab['SKY']

    cont_stats = {}
    if clean_mask is not None and _HAVE_MASK:
        sci_cont_resid = np.asarray(scitab['LINES'])
        cont_stats.update(flatten_arm_stats(
            'sci', arm_continuum_stats(wave, sci_cont_resid, clean_mask)))
        cont_stats.update(flatten_arm_stats(
            'sky', arm_continuum_stats(wave, sky_cont_resid, clean_mask)))

    return scitab, qa_flags, minimum, cont_stats


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
        Output filename root; defaults to the method name.
    '''
    x = fits.open(filename)
    drp_all    = Table(x['DRP_ALL'].data)
    final_wave = np.array(x['WAVE'].data)
    nan_spectrum = np.full(len(final_wave), np.nan)

    # Load sky-line mask once (same convention as SkySub_eval.py) purely for
    # the continuum-fit-quality evaluation below; the fit itself is unaffected
    # if the mask is unavailable.
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

    final_flux     = []
    final_sky      = []
    select         = []
    qa_flags_list  = []
    line_scale_list = []
    cont_stats_list = []

    i = 0
    while i < len(drp_all):
        try:
            ftab, row_flags, line_scale, cont_stats = one_drp(
                xfits=x, row=i, method=method, do_plot=False, clean_mask=clean_mask)
        except Exception as e:
            print('Row %d: exception (%s)' % (i, e))
            ftab = None
            row_flags = 0
            line_scale = np.nan
            cont_stats = {}

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
        cont_stats_list.append(cont_stats)
        i += idelta
        if i % 100 == 0:
            print('Completed %6d of %d in steps of %d' % (i, len(drp_all), idelta))

    n_failed  = sum(1 for f in qa_flags_list if f & QA_FAILED)
    n_warned  = sum(1 for f in qa_flags_list if f != 0 and not (f & QA_FAILED))
    print('\nProcessed %d rows: %d failed (NaN fill), %d with warnings'
          % (len(select), n_failed, n_warned))

    flagged = [(select[j], qa_flags_list[j])
               for j in range(len(select)) if qa_flags_list[j] != 0]
    if flagged:
        print('Rows with QA flags:')
        for orig_row, flags in flagged:
            active = [name for bit, name in _QA_FLAG_NAMES.items() if flags & bit]
            print('  Row %6d  flags=0x%02x  (%s)' % (orig_row, flags, ', '.join(active)))

    out_wave = np.array(x['WAVE'].data)

    hdu1 = fits.PrimaryHDU(data=None)
    hdu1.header['Title'] = 'SkySubOrig'
    hdu1.header['METHOD'] = method
    hdu2 = fits.ImageHDU(data=out_wave, name='WAVE')
    hdu3 = fits.ImageHDU(data=np.array(final_flux), name='FLUX')
    hdu4 = fits.ImageHDU(data=np.array(final_sky),  name='SKY')

    xtab = Table(x['DRP_ALL'].data)
    xtab = xtab[select]
    xtab['QA_FLAGS']   = np.array(qa_flags_list, dtype=np.int32)
    xtab['LINE_SCALE'] = np.array(line_scale_list, dtype=float)
    # Continuum-fit-quality columns (raw pre-subtraction sci/sky spectra,
    # not the final sky-subtracted result -- see one_drp() docstring).
    _cont_keys = sorted({k for d in cont_stats_list for k in d})
    for _key in _cont_keys:
        xtab[_key.upper()] = np.array(
            [d.get(_key, np.nan) for d in cont_stats_list], dtype=np.float32)
    if 'obstime' in xtab.colnames and 'mjd' in xtab.colnames:
        xtab['mjd'] = obstime_to_mjd(xtab['obstime'])
    hdu5 = fits.BinTableHDU(xtab, name='DRP_ALL')

    dwave = out_wave[1] - out_wave[0] if len(out_wave) > 1 else 0.5
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [1, 1]
    wcs.wcs.crval = [float(out_wave[0]), 0]
    wcs.wcs.cdelt = [float(dwave), 1]
    wcs.wcs.ctype = ['WAVE', 'LINE']
    hdu3.header.update(wcs.to_header())
    hdu4.header.update(wcs.to_header())

    hdul = fits.HDUList([hdu1, hdu2, hdu3, hdu4, hdu5])

    if outroot == '':
        stem = os.path.splitext(os.path.basename(filename))[0]
        outroot = '%s_orig_%s' % (stem, method)
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
