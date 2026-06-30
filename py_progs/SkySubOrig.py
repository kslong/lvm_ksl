#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Perform alternative sky subtraction on an XCframe file by separating
    continuum and sky lines, then fitting a scale factor for the lines
    independently.  Three methods are supported:

        nearest          subtract the nearest sky telescope spectrum
        farthest         subtract the farthest sky telescope spectrum
        farlines_nearcont  scale lines from the far sky, continuum from the
                           near sky (default)

    Produces a FITS file with WAVE, FLUX (sky-subtracted), SKY, and DRP_ALL
    extensions.  The DRP_ALL table carries a QA_FLAGS column that records
    per-row quality issues.

Command line usage (if any):

    usage: SkySubOrig.py [-method METHOD] [-delta N] [-out ROOT] filename

    Arguments:

    filename    XCframe FITS file to process

    Options:

    -method METHOD   sky subtraction method: nearest | farthest |
                     farlines_nearcont  (default: farlines_nearcont)
    -delta N         process every N-th row; useful for quick tests
                     (default: 1 = all rows)
    -out ROOT        output filename root; default is <stem>_<method>

Description:

    Each science spectrum is decomposed into a continuum (fitted with a
    degree-4 polynomial and sigma clipping) and line residuals.  The
    nearest and farthest sky telescopes are identified from the RA/Dec
    separation stored in DRP_ALL.  A scale factor for the sky lines is
    found by minimising::

        sum | sci_lines * (sci_lines - r * sky_lines) | / ||sky_lines||^2

    using a custom 4-point bisection search.

    QA flag bits stored in DRP_ALL['QA_FLAGS']:

        0x01  NANDATA   NaN/inf found in input flux or sky data
        0x02  ZEROSKY   sky line vector is all-zero; scale unreliable
        0x04  POORFIT   continuum polyfit poorly conditioned
        0x08  FAILED    row raised an exception; spectrum filled with NaN

Notes:

    Output filename is ``<ROOT>.fits``.  If -out is omitted the name
    is derived as ``<stem>_<method>.fits`` where ``<stem>`` is the
    input filename without extension (e.g.
    ``XCframe_1.2.1_7325_48860_1_50_farlines_nearcont.fits``).

    The RankWarning emitted by numpy.polyfit is caught portably across
    NumPy 1.x and 2.x.

History:

    250604 ksl Notebook SkySubDev250604.ipynb — initial development
    260630 ksl Converted to executable script; fixed np.RankWarning for
               NumPy 2.x; removed debug prints from ksl_bisection
    260630 ksl Renamed SkySubDev → SkySubOrig; options may precede filename;
               default output name is <stem>_<method>.fits
    260630 ksl Default wavelength range changed to full spectrum (wmin/wmax
               now None by default in one_drp and do_all)

'''

import sys
import os
import warnings as _warnings

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clip
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

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
  -out ROOT        output filename root (default: <stem>_<method>)
'''


# ──────────────────────────────────────────────────────────────
# Core functions
# ──────────────────────────────────────────────────────────────

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
            method='farlines_nearcont', do_plot=False):
    '''
    Sky-subtract a single spectrum (row) from an open FITS object.

    Returns (scitab, qa_flags).  scitab has columns WAVE, FLUX, CONT,
    MASK, LINES, SKY, SCI_FLUX.  Returns (None, QA_FAILED) on error.
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
        return None, QA_FAILED

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

    elif method == 'farthest':
        skytab = _Table([wave, sky_far], names=['WAVE', 'FLUX'])
        skytab, f = polynomial_fit_with_outliers(skytab, degree=4,
                                                 sigma_lower=3, sigma_upper=1, grow=5)
        qa_flags |= f

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

    else:
        print('Error: unknown method "%s"' % method)
        return None, QA_FAILED

    scitab['LINES'] = scitab['FLUX'] - scitab['CONT']
    if 'LINES' not in skytab.colnames:
        skytab['LINES'] = skytab['FLUX'] - skytab['CONT']

    if np.dot(skytab['LINES'], skytab['LINES']) == 0:
        qa_flags |= QA_ZEROSKY

    minimum = ksl_bisection(fit_func, 0.5, 1.5, tol=0.001, maxiter=8,
                            args=(scitab['LINES'], skytab['LINES']))

    scitab['SKY']      = skytab['CONT'] + minimum * skytab['LINES']
    scitab['SCI_FLUX'] = scitab['FLUX'] - scitab['SKY']
    return scitab, qa_flags


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

    final_flux    = []
    final_sky     = []
    select        = []
    qa_flags_list = []

    i = 0
    while i < len(drp_all):
        try:
            ftab, row_flags = one_drp(xfits=x, row=i,
                                      method=method, do_plot=False)
        except Exception as e:
            print('Row %d: exception (%s)' % (i, e))
            ftab = None
            row_flags = 0

        if ftab is None:
            row_flags |= QA_FAILED
            final_flux.append(nan_spectrum.copy())
            final_sky.append(nan_spectrum.copy())
        else:
            final_flux.append(np.array(ftab['SCI_FLUX']))
            final_sky.append(np.array(ftab['SKY']))

        select.append(i)
        qa_flags_list.append(row_flags)
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
    xtab['QA_FLAGS'] = np.array(qa_flags_list, dtype=np.int32)
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
        outroot = '%s_%s' % (stem, method)
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
