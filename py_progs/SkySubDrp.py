#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Perform sky subtraction on an XCframe summary file using the
    lvmdrp create_skysub_spectrum routine.  For each row the script
    builds the SCI/SKYE/SKYW BinTable HDU that create_skysub_spectrum
    expects, using the per-row FLUX, SKY_EAST, and SKY_WEST spectra
    together with the RA/Dec information from the DRP_ALL table, then
    calls the DRP function to determine the sky model.

    Two methods are supported:

        nearest           continuum and lines both from the nearest
                          sky telescope
        farlines_nearcont lines from the far sky, continuum from the
                          near sky (default)

    Produces a FITS file with WAVE, FLUX (sky-subtracted), SKY, and
    DRP_ALL extensions.  The DRP_ALL table carries a QA_FLAGS column
    that records per-row quality issues.

Command line usage (if any):

    usage: SkySubDrp.py [-method METHOD] [-delta N] [-out ROOT] filename

    Arguments:

    filename    XCframe FITS file to process

    Options:

    -method METHOD   sky subtraction method: nearest |
                     farlines_nearcont  (default: farlines_nearcont)
    -delta N         process every N-th row; useful for quick tests
                     (default: 1 = all rows)
    -out ROOT        output filename root; default is <stem>_<method>

Description:

    For each row the following HDU is constructed and passed to
    create_skysub_spectrum::

        PRIMARY  (empty)
        SCI      BinTable  WAVE, FLUX, ERROR   RA/DEC from sci_ra/sci_dec
        SKYE     BinTable  WAVE, FLUX, ERROR   RA/DEC from skye_ra/skye_dec
        SKYW     BinTable  WAVE, FLUX, ERROR   RA/DEC from skyw_ra/skyw_dec

    Errors are estimated as sqrt(|flux|) since the XCframe format
    does not carry IVAR.  They affect only the propagated sky error
    stored internally; the sky model itself does not depend on them.

    QA flag bits stored in DRP_ALL['QA_FLAGS']:

        0x01  NANDATA   NaN/inf found in input flux or sky data
        0x08  FAILED    row raised an exception; spectrum filled with NaN

Notes:

    Output filename is ``<ROOT>.fits``.  If -out is omitted the name
    is derived as ``<stem>_drp_<method>.fits`` where ``<stem>`` is the
    input filename without extension.

    Requires the lvmdrp26 conda environment.

History:

    260630 ksl Coding begun, modelled on SkySubOrig.py; fakes the
               sky_hdu per-row from XCframe FLUX/SKY_EAST/SKY_WEST
               and DRP_ALL RA/Dec columns

'''

import sys
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from lvmdrp.functions.skyMethod import create_skysub_spectrum

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
  SkySubDrp.py [-method METHOD] [-delta N] [-out ROOT] filename

Arguments:
  filename         XCframe FITS file to process

Options:
  -method METHOD   nearest | farlines_nearcont
                   (default: farlines_nearcont)
  -delta N         step size through rows for quick tests (default: 1)
  -out ROOT        output filename root (default: <stem>_drp_<method>)
'''


# ──────────────────────────────────────────────────────────────
# Core functions
# ──────────────────────────────────────────────────────────────

def _build_sky_hdu(wave, sci_flux, skye_flux, skyw_flux,
                   sci_ra, sci_dec, skye_ra, skye_dec, skyw_ra, skyw_dec):
    '''
    Build the SCI/SKYE/SKYW BinTable HDUList expected by
    create_skysub_spectrum from per-row XCframe data.
    '''
    def _make_ext(name, flux, ra, dec):
        err = np.sqrt(np.abs(flux))
        tab = Table([wave, flux, err], names=['WAVE', 'FLUX', 'ERROR'])
        hdu = fits.BinTableHDU(tab, name=name)
        hdu.header['RA']  = float(ra)
        hdu.header['DEC'] = float(dec)
        return hdu

    return fits.HDUList([
        fits.PrimaryHDU(),
        _make_ext('SCI',  sci_flux,  sci_ra,  sci_dec),
        _make_ext('SKYE', skye_flux, skye_ra, skye_dec),
        _make_ext('SKYW', skyw_flux, skyw_ra, skyw_dec),
    ])


def one_drp(xfits, drp_all, row, wave, wmin=None, wmax=None,
            method='farlines_nearcont'):
    '''
    Sky-subtract a single row from an open XCframe FITS object.

    Returns (result_table, qa_flags).  result_table has columns
    WAVE, SCI_FLUX, SKY.  Returns (None, QA_FAILED) on error.
    '''
    qa_flags = 0

    flux      = np.array(xfits['FLUX'].data[row],      dtype=float)
    skye_flux = np.array(xfits['SKY_EAST'].data[row],  dtype=float)
    skyw_flux = np.array(xfits['SKY_WEST'].data[row],  dtype=float)

    if not (np.all(np.isfinite(flux)) and
            np.all(np.isfinite(skye_flux)) and
            np.all(np.isfinite(skyw_flux))):
        qa_flags |= QA_NANDATA
        # replace NaN so the DRP function can still run
        flux      = np.nan_to_num(flux)
        skye_flux = np.nan_to_num(skye_flux)
        skyw_flux = np.nan_to_num(skyw_flux)

    sci_ra   = drp_all['sci_ra'][row]
    sci_dec  = drp_all['sci_dec'][row]
    skye_ra  = drp_all['skye_ra'][row]
    skye_dec = drp_all['skye_dec'][row]
    skyw_ra  = drp_all['skyw_ra'][row]
    skyw_dec = drp_all['skyw_dec'][row]

    if wmin is not None or wmax is not None:
        lo = wmin if wmin is not None else -np.inf
        hi = wmax if wmax is not None else  np.inf
        wcut     = (wave > lo) & (wave < hi)
        wave_cut = wave[wcut]
        flux_cut = flux[wcut]
        skye_cut = skye_flux[wcut]
        skyw_cut = skyw_flux[wcut]
    else:
        wave_cut = wave
        flux_cut = flux
        skye_cut = skye_flux
        skyw_cut = skyw_flux

    sky_hdu = _build_sky_hdu(wave_cut, flux_cut, skye_cut, skyw_cut,
                             sci_ra, sci_dec, skye_ra, skye_dec, skyw_ra, skyw_dec)

    sky, _sky_err = create_skysub_spectrum(sky_hdu, tel='sci', method=method)

    result = Table([wave_cut, flux_cut - sky, sky],
                   names=['WAVE', 'SCI_FLUX', 'SKY'])
    return result, qa_flags


def do_all(filename, method='farlines_nearcont', idelta=1, outroot=''):
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
    outroot : str
        Output filename root; defaults to <stem>_<method>.
    '''
    x        = fits.open(filename)
    drp_all  = Table(x['DRP_ALL'].data)
    wave     = np.array(x['WAVE'].data)
    out_wave = wave

    nan_spectrum = np.full(len(out_wave), np.nan)

    final_flux    = []
    final_sky     = []
    select        = []
    qa_flags_list = []

    i = 0
    while i < len(drp_all):
        try:
            ftab, row_flags = one_drp(x, drp_all, i, wave, method=method)
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

    out_wave = wave

    hdu1 = fits.PrimaryHDU(data=None)
    hdu1.header['Title']  = 'SkySubDrp'
    hdu1.header['METHOD'] = method
    hdu2 = fits.ImageHDU(data=out_wave,               name='WAVE')
    hdu3 = fits.ImageHDU(data=np.array(final_flux),   name='FLUX')
    hdu4 = fits.ImageHDU(data=np.array(final_sky),    name='SKY')

    xtab = drp_all[select].copy()
    xtab['QA_FLAGS'] = np.array(qa_flags_list, dtype=np.int32)
    hdu5 = fits.BinTableHDU(xtab, name='DRP_ALL')

    dwave = float(out_wave[1] - out_wave[0]) if len(out_wave) > 1 else 0.5
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [1, 1]
    wcs.wcs.crval = [float(out_wave[0]), 0]
    wcs.wcs.cdelt = [dwave, 1]
    wcs.wcs.ctype = ['WAVE', 'LINE']
    hdu3.header.update(wcs.to_header())
    hdu4.header.update(wcs.to_header())

    hdul = fits.HDUList([hdu1, hdu2, hdu3, hdu4, hdu5])

    if outroot == '':
        stem    = os.path.splitext(os.path.basename(filename))[0]
        outroot = '%s_drp_%s' % (stem, method)
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

    do_all(filename=filename, method=method, idelta=idelta, outroot=outroot)
