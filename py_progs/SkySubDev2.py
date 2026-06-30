#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Perform sky subtraction on an XCframe summary file using the PALACE
    (Ivan Katkov / XSkySepIvan) spectral decomposition to separate the
    sky into emission lines and continuum.  For each row the near- and
    far-sky spectra are decomposed with the PALACE model; the resulting
    line and continuum components are then combined and subtracted from
    the science spectrum without any additional scaling.

    Two methods are supported:

        nearest           continuum and lines both from the nearest
                          sky telescope
        farlines_nearcont lines from the far sky, continuum from the
                          near sky (default)

    Produces a FITS file with WAVE, FLUX (sky-subtracted), SKY, and
    DRP_ALL extensions.  The DRP_ALL table carries a QA_FLAGS column
    that records per-row quality issues.

Command line usage (if any):

    usage: SkySubDev2.py [-method METHOD] [-delta N] [-lsf FWHM]
                         [-out ROOT] filename

    Arguments:

    filename    XCframe FITS file to process

    Options:

    -method METHOD   sky subtraction method: nearest |
                     farlines_nearcont  (default: farlines_nearcont)
    -delta N         process every N-th row; useful for quick tests
                     (default: 1 = all rows)
    -lsf FWHM        LSF FWHM in Angstroms (default: 1.3)
    -out ROOT        output filename root; default is <stem>_dev2_<method>

Description:

    The PALACE decomposer (SkyDecomp) is built once from the full
    wavelength grid and LSF, then reused for every row.  For each row:

    1. Science and sky spectra are read for the row.
    2. Near/far sky is determined from sci_ra/dec, skye_ra/dec,
       skyw_ra/dec in DRP_ALL.
    3. The near-sky and far-sky spectra are each decomposed by the
       PALACE model into::

           LINES = oh + atom + orc + o2
           CONT  = moon + diffuse

    4. The sky model is assembled without scaling:

       farlines_nearcont:  sky = CONT_near + LINES_far
       nearest:            sky = CONT_near + LINES_near

    5. sky-subtracted science = flux_sci - sky

    QA flag bits stored in DRP_ALL['QA_FLAGS']:

        0x01  NANDATA   NaN/inf found in input flux or sky data
        0x08  FAILED    row raised an exception; spectrum filled with NaN

Notes:

    Output filename is ``<ROOT>.fits``.  If -out is omitted the name
    is derived as ``<stem>_dev2_<method>.fits`` where ``<stem>`` is the
    input filename without extension.

    Requires the PALACE library (lvmsky/skysub/sky_decomp) and the
    palace data files; the paths are taken from XSkySepIvan.py.

History:

    260630 ksl Coding begun; imports PALACE decomposer from XSkySepIvan.py

'''

import sys
import os

# ensure py_progs siblings are importable when running directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u

from XSkySepIvan import _get_decomposer, estimate_ivar, DEFAULT_BASE_DIR

# ──────────────────────────────────────────────────────────────
# QA flag bits
# ──────────────────────────────────────────────────────────────
QA_NANDATA = 1   # NaN or inf values found in input flux or sky data
QA_FAILED  = 8   # row failed entirely; FLUX and SKY are NaN

_QA_FLAG_NAMES = {
    QA_NANDATA: 'NANDATA',
    QA_FAILED:  'FAILED',
}

_USAGE = '''Usage:
  SkySubDev2.py [-method METHOD] [-delta N] [-lsf FWHM] [-out ROOT] filename

Arguments:
  filename         XCframe FITS file to process

Options:
  -method METHOD   nearest | farlines_nearcont
                   (default: farlines_nearcont)
  -delta N         step size through rows for quick tests (default: 1)
  -lsf FWHM        LSF FWHM in Angstroms (default: 1.3)
  -out ROOT        output filename root (default: <stem>_dev2_<method>)
'''


# ──────────────────────────────────────────────────────────────
# PALACE decomposition of a raw flux array
# ──────────────────────────────────────────────────────────────

def _decompose(flux, wave, decomposer):
    '''
    Run the PALACE decomposer on a single flux array and return
    (lines, cont) where lines = oh+atom+orc+o2 and cont = moon+diffuse.
    '''
    flux = np.asarray(flux, float)
    ivar = estimate_ivar(flux, wave)

    finite = np.isfinite(flux)
    flux_scale = float(np.sqrt(np.nanmean(flux[finite] ** 2))) if finite.any() else 1.0
    flux_scale = max(flux_scale, 1e-30)

    result = decomposer.fit(flux / flux_scale, ivar * flux_scale ** 2)
    c      = result.components
    lines  = (c['oh'] + c['atom'] + c['orc'] + c['o2']) * flux_scale
    cont   = (c['moon'] + c['diffuse']) * flux_scale
    return lines, cont


# ──────────────────────────────────────────────────────────────
# Per-row sky subtraction
# ──────────────────────────────────────────────────────────────

def one_drp(xfits, drp_all, row, wave, decomposer,
            method='farlines_nearcont'):
    '''
    Sky-subtract a single row using PALACE decomposition with no scaling.

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
    sci_coord  = SkyCoord(ra=drp_all['sci_ra'][row]   * u.degree,
                          dec=drp_all['sci_dec'][row]  * u.degree)
    skye_coord = SkyCoord(ra=drp_all['skye_ra'][row]  * u.degree,
                          dec=drp_all['skye_dec'][row] * u.degree)
    skyw_coord = SkyCoord(ra=drp_all['skyw_ra'][row]  * u.degree,
                          dec=drp_all['skyw_dec'][row] * u.degree)
    if sci_coord.separation(skye_coord) < sci_coord.separation(skyw_coord):
        sky_near, sky_far = skye_flux, skyw_flux
    else:
        sky_near, sky_far = skyw_flux, skye_flux

    # PALACE decomposition — no scaling applied
    lines_near, cont_near = _decompose(sky_near, wave, decomposer)

    if method == 'farlines_nearcont':
        lines_far, _ = _decompose(sky_far, wave, decomposer)
        sky = cont_near + lines_far
    elif method == 'nearest':
        sky = cont_near + lines_near
    else:
        print('Error: unknown method "%s"' % method)
        return None, QA_FAILED

    result = Table([wave, flux - sky, sky], names=['WAVE', 'SCI_FLUX', 'SKY'])
    return result, qa_flags


# ──────────────────────────────────────────────────────────────
# Batch processing
# ──────────────────────────────────────────────────────────────

def do_all(filename, method='farlines_nearcont', idelta=1,
           fwhm_lsf=1.3, outroot=''):
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
    fwhm_lsf : float
        LSF FWHM in Angstroms.
    outroot : str
        Output filename root; defaults to <stem>_dev2_<method>.
    '''
    x       = fits.open(filename)
    drp_all = Table(x['DRP_ALL'].data)
    wave    = np.array(x['WAVE'].data, dtype=float)

    lsf_sigma  = fwhm_lsf / 2.355
    decomposer = _get_decomposer(wave, lsf_sigma, base_dir=DEFAULT_BASE_DIR)

    nan_spectrum = np.full(len(wave), np.nan)

    final_flux    = []
    final_sky     = []
    select        = []
    qa_flags_list = []

    i = 0
    while i < len(drp_all):
        try:
            ftab, row_flags = one_drp(x, drp_all, i, wave, decomposer,
                                      method=method)
        except Exception as e:
            print('Row %d: exception (%s)' % (i, e))
            ftab      = None
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
    hdu1.header['Title']   = 'SkySubDev2'
    hdu1.header['METHOD']  = method
    hdu1.header['FWHMLSF'] = fwhm_lsf
    hdu2 = fits.ImageHDU(data=wave.astype(np.float32),       name='WAVE')
    hdu3 = fits.ImageHDU(data=np.array(final_flux),          name='FLUX')
    hdu4 = fits.ImageHDU(data=np.array(final_sky),           name='SKY')

    xtab = drp_all[select].copy()
    xtab['QA_FLAGS'] = np.array(qa_flags_list, dtype=np.int32)
    hdu5 = fits.BinTableHDU(xtab, name='DRP_ALL')

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
        outroot = '%s_dev2_%s' % (stem, method)
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

    method   = 'farlines_nearcont'
    idelta   = 1
    fwhm_lsf = 1.3
    outroot  = ''
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
        elif arg == '-lsf':
            i += 1
            fwhm_lsf = float(argv[i])
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

    valid_methods = {'nearest', 'farlines_nearcont'}
    if method not in valid_methods:
        print('Error: -method must be one of: %s' % ', '.join(sorted(valid_methods)))
        sys.exit(1)

    do_all(filename=filename, method=method, idelta=idelta,
           fwhm_lsf=fwhm_lsf, outroot=outroot)
