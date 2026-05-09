#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Fit Gaussians to a fixed set of airglow lines in each row of a
SummarizeCframe output file, producing one output row per exposure.


Command line usage:

    gauss_offset.py [-ext FLUX] [-out root] filename [filename ...]

    where

    -ext extname  Extension to fit: FLUX (default), SKY_EAST, or SKY_WEST
    -out root     Root name for the output file.  If omitted the output is
                  named <input_stem>_gauss.<ext>.fits, e.g.
                  XCframe_1.2.0_10000_20000_10_50_gauss.FLUX.fits
    filename      One or more SummarizeCframe FITS files


Description:

Each row of the chosen flux extension in a SummarizeCframe file contains the
spectrum for one exposure (median science-fiber flux for FLUX, sky-telescope
spectra for SKY_EAST and SKY_WEST).  This script fits single Gaussians to the
following airglow lines in every row:

    sky5577   5577.34 A
    sky6300   6300.3  A
    sky6363   6363.8  A
    sky6533   6533.04 A
    sky6553   6553.0  A
    sky6577   6577.2  A
    sky6912   6912.6  A
    sky6923   6923.2  A
    sky6939   6939.5  A

For each line the fit returns the flux, centroid wavelength, FWHM, background
level, and RMSE.

Primary routines:

    do_one  - fit all airglow lines for a single spectrum table
    do_all  - iterate over all rows in a SummarizeCframe FITS file


Output:

A FITS binary table with one row per exposure.  Columns cover the fit
parameters (flux, wave, fwhm, back, rmse) for each line together with expnum
and mjd drawn from the drp_all extension.

The automatic output filename is <input_stem>_gauss.<ext>.fits, where <ext>
is the extension that was fit (e.g. FLUX, SKY_EAST, SKY_WEST).  This ensures
the output file is never confused with the input and that FLUX, SKY_EAST, and
SKY_WEST runs on the same input file do not overwrite each other.

If -out root is supplied the output is written to root.fits (adding .fits if
not already present).


Notes:

Because the input spectra are from CFrame (before sky subtraction), the
airglow lines are bright and well-suited to measuring centroid wavelengths
and fluxes as a function of exposure, which can reveal wavelength-calibration
drifts or changes in sky brightness over a night or survey.

History:

240604 ksl Coding begun (as sky_gaussfit.py)
260504 ksl Adapted for SummarizeCframe input; restricted to airglow lines only
260504 ksl FITS output; selectable extension (-ext)

'''

import sys
import os
from astropy.table import Table, vstack, hstack
from astropy.io import fits
from astropy.time import Time

from lvm_ksl.lvm_gaussfit import fit_gaussian_to_spectrum, check_for_nan


# Wavelengths taken from
# https://www.eso.org/observing/dfo/quality/UVES/pipeline/sky_spectrum.html
AIRGLOW_LINES = [
    ('sky5577', 5577.34668,  5572., 5582.),
    ('sky6300', 6300.308594, 6295., 6305.),
    ('sky6363', 6363.782715, 6358., 6368.),
    ('sky7358', 7358.680176, 7353., 7363.),
    ('sky7392', 7392.209961, 7387., 7397.),
    ('sky7914', 7913.717773, 7908., 7918.),
    ('sky8344', 8344.613281, 8339., 8349.),
    ('sky8399', 8399.175781, 8394., 8404.),
    ('sky8827', 8827.112305, 8822., 8832.),
    ('sky8988', 8988.383789, 8983., 8993.),
    ('sky9552', 9552.546875, 9547., 9557.),
    ('sky9719', 9719.838867, 9714., 9724.),
]


def do_one(spectrum_table):
    '''
    Fit all airglow lines for a single spectrum supplied as an astropy Table
    with WAVE and FLUX columns.  Returns an hstacked single-row Table of fit
    results, or an empty list if nothing could be fit.
    '''
    records = []

    for xname, xwcen, xwmin, xwmax in AIRGLOW_LINES:
        try:
            results, xspec = fit_gaussian_to_spectrum(
                spectrum_table, line=xname,
                init_wavelength=xwcen, init_fwhm=1.,
                wavelength_min=xwmin, wavelength_max=xwmax)
            records.append(results)
        except Exception as e:
            print('Fitting %s: An exception occurred: %s' % (xname, e))

    if len(records) == 0:
        print('Nothing fit for this spectrum')
        return []

    return hstack(records)


def do_all(filename, ext='FLUX', outname=''):
    '''
    Fit airglow lines for every exposure row in a SummarizeCframe FITS file.

    Parameters
    ----------
    filename : str
        Path to a SummarizeCframe output FITS file.
    ext : str
        Name of the flux extension to fit (FLUX, SKY_EAST, or SKY_WEST).
    outname : str
        Output filename.  If empty, auto-generated from the input stem.

    Returns
    -------
    results : astropy Table or None
    '''
    try:
        x = fits.open(filename)
    except Exception:
        print('Error: Could not open %s' % filename)
        return

    wave = x['WAVE'].data
    try:
        flux = x[ext].data
    except KeyError:
        print('Error: extension %s not found in %s' % (ext, filename))
        return

    drp = Table(x['drp_all'].data)
    n_exp = flux.shape[0]

    records = []
    for i in range(n_exp):
        one_spec = Table([wave, flux[i]], names=['WAVE', 'FLUX'])
        if check_for_nan(flux[i]):
            print('Too many NaNs for row %d (expnum %d), skipping'
                  % (i, drp['expnum'][i]))
            continue
        rtab = do_one(one_spec)
        if len(rtab) > 0:
            rtab['expnum'] = int(drp['expnum'][i])
            rtab['mjd']    = Time(drp['obstime'][i], format='isot', scale='utc').mjd
            records.append(rtab)
        else:
            print('Nothing fit for row %d (expnum %d)' % (i, drp['expnum'][i]))

    if len(records) == 0:
        print('No spectra successfully fit in %s' % filename)
        return

    results = vstack(records)

    if outname == '':
        stem = os.path.basename(filename).replace('.fits', '')
        outname = '%s_gauss.%s.fits' % (stem, ext)
    else:
        if not outname.endswith('.fits'):
            outname = outname + '.fits'

    results.write(outname, format='fits', overwrite=True)
    print('Wrote results to %s' % outname)
    return results


def steer(argv):
    outname = ''
    ext = 'FLUX'
    fitsfiles = []

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-ext':
            i += 1
            ext = argv[i]
        elif argv[i][0:4] == '-out':
            i += 1
            outname = argv[i]
        elif argv[i][0] == '-':
            print('Unknown option: %s' % argv[i])
            return
        elif argv[i].endswith('.fits'):
            fitsfiles.append(argv[i])
        else:
            print('Unrecognised argument: %s' % argv[i])
            return
        i += 1

    for one_file in fitsfiles:
        do_all(one_file, ext, outname)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
