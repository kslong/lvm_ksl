#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Perform sky subtraction on an XCframe summary file using the
    B-spline continuum separation method from GetSkyCont.py instead
    of the polynomial fit used by SkySubOrig.py.  In all other
    respects the algorithm is the same: for each row the continuum
    and line residuals of the science and sky spectra are separated,
    a scale factor for the sky lines is found by bisection, and the
    scaled sky is subtracted.

    Two methods are supported:

        nearest           continuum and lines both from the nearest
                          sky telescope
        farlines_nearcont lines from the far sky, continuum from the
                          near sky (default)

    Produces a FITS file with WAVE, FLUX (sky-subtracted), SKY, and
    DRP_ALL extensions.  The DRP_ALL table carries a QA_FLAGS column
    that records per-row quality issues.

Command line usage (if any):

    usage: SkySubDev1.py [-method METHOD] [-delta N] [-kstep N]
                         [-out ROOT] -mask mask.fits filename

    Arguments:

    filename    XCframe FITS file to process

    Options:

    -mask file       palace_mask FITS file from palace_make_mask.py
                     (MASK extension: 1=clean).  If omitted the script
                     searches for sky_mask.fits in the current directory
                     then in the lvm_ksl data/ directory.
    -method METHOD   sky subtraction method: nearest |
                     farlines_nearcont  (default: farlines_nearcont)
    -delta N         process every N-th row; useful for quick tests
                     (default: 1 = all rows)
    -kstep N         B-spline knot spacing in Angstroms (default 100)
    -out ROOT        output filename root; default is <stem>_dev1_<method>

Description:

    The design matrix (two-component B-spline: DIFFUSE + MOON) is
    built once from the full wavelength grid and the palace mask,
    then reused for every row.  For each row:

    1. Science and sky spectra are read for the row.
    2. Near/far sky is determined from sci_ra/dec, skye_ra/dec,
       skyw_ra/dec in DRP_ALL.
    3. fit_continuum (from GetSkyCont.py) decomposes each spectrum
       into continuum and line residuals using non-negative least
       squares on the clean (mask=1) pixels.
    4. A global line scale factor r is found by ksl_bisection
       (from SkySubOrig.py) minimising::

           sum |sci_lines * (sci_lines - r * sky_lines)| / ||sky_lines||^2

    5. For farlines_nearcont: sky = cont_near + r * lines_far
       For nearest:           sky = cont_near + r * lines_near
    6. sky-subtracted science = flux_sci - sky

    QA flag bits stored in DRP_ALL['QA_FLAGS']:

        0x01  NANDATA   NaN/inf found in input flux or sky data
        0x02  ZEROSKY   sky line vector is all-zero; scale unreliable
        0x08  FAILED    row raised an exception; spectrum filled with NaN

Notes:

    Output filename is ``<ROOT>.fits``.  If -out is omitted the name
    is derived as ``<stem>_dev1_<method>.fits`` where ``<stem>`` is
    the input filename without extension.

    If -mask is omitted the script looks for ``sky_mask.fits`` first
    in the current working directory, then in the lvm_ksl ``data/``
    directory (sibling of ``py_progs/``).

History:

    260630 ksl Coding begun; imports continuum separation from
               GetSkyCont.py and bisection from SkySubOrig.py

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

from GetSkyCont import (build_design_matrix, fit_continuum,
                        load_mask, load_solar, _interp_mask_to_wave)
from SkySubOrig import fit_func, ksl_bisection

# ──────────────────────────────────────────────────────────────
# QA flag bits
# ──────────────────────────────────────────────────────────────
QA_NANDATA = 1   # NaN or inf values found in input flux or sky data
QA_ZEROSKY = 2   # sky line vector is all-zero; scale factor is unreliable
QA_FAILED  = 8   # row failed entirely; FLUX and SKY are NaN

_QA_FLAG_NAMES = {
    QA_NANDATA: 'NANDATA',
    QA_ZEROSKY: 'ZEROSKY',
    QA_FAILED:  'FAILED',
}

_USAGE = '''Usage:
  SkySubDev1.py [-method METHOD] [-delta N] [-kstep N] [-out ROOT]
                -mask mask.fits filename

Arguments:
  filename         XCframe FITS file to process

Options:
  -mask file       palace_mask FITS file (default: sky_mask.fits
                   searched in cwd then data/)
  -method METHOD   nearest | farlines_nearcont
                   (default: farlines_nearcont)
  -delta N         step size through rows for quick tests (default: 1)
  -kstep N         B-spline knot spacing in Angstroms (default: 100)
  -out ROOT        output filename root (default: <stem>_dev1_<method>)
'''


# ──────────────────────────────────────────────────────────────
# Per-row sky subtraction
# ──────────────────────────────────────────────────────────────

def one_drp(xfits, drp_all, row, wave, clean, A, n_b,
            method='farlines_nearcont'):
    '''
    Sky-subtract a single row from an open XCframe FITS object.

    Uses the pre-built B-spline design matrix (A, n_b) and clean-pixel
    mask to separate continuum and lines in sci, near-sky, and far-sky
    spectra, then scales and subtracts the sky.

    Returns (result_table, qa_flags).  result_table has columns WAVE,
    SCI_FLUX, SKY.  Returns (None, QA_FAILED) on error.
    '''
    qa_flags = 0

    flux      = np.array(xfits['FLUX'].data[row],     dtype=float)
    skye_flux = np.array(xfits['SKY_EAST'].data[row], dtype=float)
    skyw_flux = np.array(xfits['SKY_WEST'].data[row], dtype=float)

    if not (np.all(np.isfinite(flux)) and
            np.all(np.isfinite(skye_flux)) and
            np.all(np.isfinite(skyw_flux))):
        qa_flags |= QA_NANDATA
        flux      = np.nan_to_num(flux)
        skye_flux = np.nan_to_num(skye_flux)
        skyw_flux = np.nan_to_num(skyw_flux)

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

    # decompose science and sky spectra into continuum + lines
    cont_sci,  _, _, _, _ = fit_continuum(flux,     clean, A, n_b)
    cont_near, _, _, _, _ = fit_continuum(sky_near, clean, A, n_b)
    cont_far,  _, _, _, _ = fit_continuum(sky_far,  clean, A, n_b)

    lines_sci  = flux     - cont_sci
    lines_near = sky_near - cont_near
    lines_far  = sky_far  - cont_far

    if method == 'farlines_nearcont':
        use_lines = lines_far
    elif method == 'nearest':
        use_lines = lines_near
    else:
        print('Error: unknown method "%s"' % method)
        return None, QA_FAILED

    if np.dot(use_lines, use_lines) == 0:
        qa_flags |= QA_ZEROSKY

    r = ksl_bisection(fit_func, 0.5, 1.5, tol=0.001, maxiter=8,
                      args=(lines_sci, use_lines))

    sky         = cont_near + r * use_lines
    sci_flux_sub = flux - sky

    result = Table([wave, sci_flux_sub, sky], names=['WAVE', 'SCI_FLUX', 'SKY'])
    return result, qa_flags


# ──────────────────────────────────────────────────────────────
# Batch processing
# ──────────────────────────────────────────────────────────────

def do_all(filename, mask_file, method='farlines_nearcont',
           idelta=1, knot_step=100.0, outroot=''):
    '''
    Process every row of an XCframe file and write sky-subtracted output.

    Parameters
    ----------
    filename : str
        Path to the input XCframe FITS file.
    mask_file : str
        Path to the palace_mask FITS file (MASK extension: 1=clean).
    method : str
        Sky subtraction method (nearest or farlines_nearcont).
    idelta : int
        Row step size (1 = all rows).
    knot_step : float
        B-spline knot spacing in Angstroms.
    outroot : str
        Output filename root; defaults to <stem>_dev1_<method>.
    '''
    x       = fits.open(filename)
    drp_all = Table(x['DRP_ALL'].data)
    wave    = np.array(x['WAVE'].data, dtype=float)

    # load mask and build design matrix once
    mask_wave, mask_arr = load_mask(mask_file)
    clean = _interp_mask_to_wave(mask_wave, mask_arr, wave)
    print('Clean pixels: %d / %d' % (int(clean.sum()), len(wave)))

    solar  = load_solar(wave)
    A, n_b = build_design_matrix(wave, knot_step=knot_step, solar=solar)
    print('Design matrix: %d columns (%d diffuse + %d moon)'
          % (A.shape[1], n_b, A.shape[1] - n_b))

    nan_spectrum = np.full(len(wave), np.nan)

    final_flux    = []
    final_sky     = []
    select        = []
    qa_flags_list = []

    i = 0
    while i < len(drp_all):
        try:
            ftab, row_flags = one_drp(x, drp_all, i, wave, clean, A, n_b,
                                      method=method)
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
    hdu1.header['Title']  = 'SkySubDev1'
    hdu1.header['METHOD'] = method
    hdu1.header['MASK']   = os.path.basename(mask_file)
    hdu1.header['KSTEP']  = knot_step
    hdu2 = fits.ImageHDU(data=wave.astype(np.float32),         name='WAVE')
    hdu3 = fits.ImageHDU(data=np.array(final_flux),            name='FLUX')
    hdu4 = fits.ImageHDU(data=np.array(final_sky),             name='SKY')

    xtab = drp_all[select].copy()
    xtab['QA_FLAGS'] = np.array(qa_flags_list, dtype=np.int32)
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
        outroot = '%s_dev1_%s' % (stem, method)
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

    method    = 'farlines_nearcont'
    idelta    = 1
    knot_step = 100.0
    outroot   = ''
    mask_file = ''
    filename  = None

    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == '-method':
            i += 1
            method = argv[i]
        elif arg == '-delta':
            i += 1
            idelta = int(argv[i])
        elif arg == '-kstep':
            i += 1
            knot_step = float(argv[i])
        elif arg == '-out':
            i += 1
            outroot = argv[i]
        elif arg == '-mask':
            i += 1
            mask_file = argv[i]
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

    if not mask_file:
        _data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 '..', 'data')
        for _candidate in [os.path.join(os.getcwd(), 'sky_mask.fits'),
                            os.path.join(_data_dir,   'sky_mask.fits')]:
            if os.path.exists(_candidate):
                mask_file = _candidate
                print('Using default mask: %s' % mask_file)
                break
        if not mask_file:
            print('Error: no mask file supplied and sky_mask.fits not found '
                  'in the current directory or data/')
            print(_USAGE)
            sys.exit(1)
    elif not os.path.exists(mask_file):
        print('Error: mask file not found: %s' % mask_file)
        sys.exit(1)

    valid_methods = {'nearest', 'farlines_nearcont'}
    if method not in valid_methods:
        print('Error: -method must be one of: %s' % ', '.join(sorted(valid_methods)))
        sys.exit(1)

    do_all(filename=filename, mask_file=mask_file, method=method,
           idelta=idelta, knot_step=knot_step, outroot=outroot)
