#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Create a summary of CFrame (flux-calibrated, before sky subtraction) data
where each row contains the percentile spectrum for each of the three LVM
spectrographs from a single exposure.  Useful for evaluating
spectrograph-to-spectrograph variations before and after sky subtraction.

Command line usage::

    SummarizeSpec.py [-h] [-sf] [-ver drp_ver] [-percent 50] [-emin 900]
                     [-out name] exp_start exp_stop delta

Description:

This script processes multiple CFrame files (the default) and computes the
percentile spectrum for science fibers belonging to each of the three LVM
spectrographs (spectrographid = 1, 2, 3) as recorded in the SLITMAP
extension.  Only fibers with telescope == 'Sci' and fibstatus == 0 are
included.  The -sf switch causes SFrame files to be read instead.

File paths are taken from the location column of the drpall file, which
records SFrame paths.  For CFrame files the path is adjusted automatically
by replacing 'SFrame' with 'CFrame' in the location string.

An exposure is included in the output only if science fibers are found
for all three spectrographs.  Exposures where one or more spectrographs
are missing are skipped and recorded in a separate ASCII table.

Two output files are produced.  The main FITS file contains one row per
valid exposure in each flux extension.  The skipped-exposure table lists
every exposure that was rejected, with columns indicating which
spectrographs were present and which were absent.

Options: -h prints this documentation; -sf reads SFrame instead of CFrame
files; -ver drp_ver sets the DRP version (default 1.2.0); -percent N sets
the percentile to compute (default 50 = median); -emin N sets the minimum
exposure time in seconds to include (default 900); -out name sets the root
name of the output files.

Positional arguments: exp_start is the first exposure number to consider;
exp_stop is the last exposure number; delta processes every Nth exposure.

Primary routines:

    doit

Notes:

The three flux extensions are named FLUX1, FLUX2, and FLUX3, corresponding
to spectrographid values 1, 2, and 3 respectively.  Each is a 2-D array
with one row per valid exposure and one column per wavelength pixel.

The skipped-exposure table is always written, even if no exposures were
skipped, so there is always a record that the completeness check was
performed.

History:

260312 ksl Coding begun, based on SummarizeRings.py

'''

import sys
import os
import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u


def augment_drp_all(xtab):
    '''
    Add survey classification, near/far sky telescope, and redshift columns
    to a drpall table.

    Parameters:
        xtab (astropy.table.Table): drpall table with sci_ra, sci_dec,
            skye_ra, skye_dec, skyw_ra, skyw_dec columns.

    Returns:
        astropy.table.Table: Input table with Survey, Near, Far, and
        Redshift columns added.
    '''
    drp_all = xtab
    sci_ra = drp_all['sci_ra']
    sci_dec = drp_all['sci_dec']
    skye_ra = drp_all['skye_ra']
    skye_dec = drp_all['skye_dec']
    skyw_ra = drp_all['skyw_ra']
    skyw_dec = drp_all['skyw_dec']

    sci_coord = SkyCoord(ra=sci_ra * u.degree, dec=sci_dec * u.degree, frame='icrs')
    skye_coord = SkyCoord(ra=skye_ra * u.degree, dec=skye_dec * u.degree, frame='icrs')
    skyw_coord = SkyCoord(ra=skyw_ra * u.degree, dec=skyw_dec * u.degree, frame='icrs')
    de = sci_coord.separation(skye_coord).degree
    dw = sci_coord.separation(skyw_coord).degree
    galactic_lat = sci_coord.galactic.b.degree

    lmc_pos = SkyCoord(ra=80.89416666666668 * u.degree, dec=-69.75618 * u.degree, frame='icrs')
    lmc_sep = lmc_pos.separation(sci_coord).degree
    smc_pos = SkyCoord(ra=13.1875 * u.degree, dec=-72.8286 * u.degree, frame='icrs')
    smc_sep = smc_pos.separation(sci_coord).degree

    xlocal = np.select([np.fabs(galactic_lat) > 10], ['HighLat'], default='Plane')
    xlocal = np.select([lmc_sep < 8], ['LMC'], default=xlocal)
    xlocal = np.select([smc_sep < 5], ['SMC'], default=xlocal)
    drp_all['Survey'] = xlocal

    near = np.select([de < dw], ['SKY_EAST'], default='SKY_WEST')
    far = np.select([de >= dw], ['SKY_EAST'], default='SKY_WEST')
    drp_all['Near'] = near
    drp_all['Far'] = far
    drp_all['Redshift'] = np.select(
        [drp_all['Survey'] == 'LMC', drp_all['Survey'] == 'SMC'],
        [262., 146.], default=0.0)
    return drp_all


def read_drpall(drp_ver='1.2.0'):
    '''
    Read the drpall FITS file and return an augmented metadata table.

    Looks for the file locally first, then falls back to the Utah path.

    Parameters:
        drp_ver (str): DRP version string (default '1.2.0').

    Returns:
        astropy.table.Table: drpall table with survey classification
        columns added, or an empty list if the file cannot be found
        or read.
    '''
    DRPFILE = 'drpall-%s.fits' % drp_ver
    if os.path.isfile(DRPFILE):
        xfile = DRPFILE
    else:
        BASEDIR = '/uufs/chpc.utah.edu/common/home/sdss51/sdsswork/lvm/spectro/redux/%s/' % drp_ver
        xfile = '%s/%s' % (BASEDIR, DRPFILE)
        if not os.path.isfile(xfile):
            print('Error: Could not locate: %s' % xfile)
            return []

    print('Opening: %s' % xfile)
    try:
        drpall = fits.open(xfile)
    except Exception:
        print('Error: Located but could not read: %s' % xfile)
        return []

    drp_tab = Table(drpall[1].data)
    drp_tab = augment_drp_all(drp_tab)
    return drp_tab


def select(ztab, exp_start=4000, exp_stop=8000, delta=5, exp_min=900.):
    '''
    Select every Nth exposure between exp_start and exp_stop.

    Parameters:
        ztab (astropy.table.Table): drpall table.
        exp_start (int): First exposure number to include.
        exp_stop (int): Last exposure number to include.
        delta (int): Step size; selects every Nth row.
        exp_min (float): Minimum exposure time in seconds.

    Returns:
        astropy.table.Table: Filtered table.
    '''
    xtab = ztab[ztab['expnum'] >= exp_start]
    xtab = xtab[xtab['expnum'] <= exp_stop]
    if exp_min > 0:
        xtab = xtab[xtab['exptime'] >= exp_min]
    if delta > 1:
        xtab = xtab[::delta]
    return xtab


XTOP = '/uufs/chpc.utah.edu/common/home/sdss51/'
XRAINBOW = '/Users/long/Projects/lvm_data/sas'
XMUSKIE = '/home/long/Projects/lvm_data/sas'


def find_top():
    '''
    Determine the top-level data directory for the current machine.

    Returns:
        str: Path to the top-level data directory, or an empty string
        if the location cannot be determined.
    '''
    if os.path.isdir(XTOP):
        loc = 'Utah'
        topdir = XTOP
    elif os.path.isdir(XRAINBOW):
        loc = 'Rainbow'
        topdir = XRAINBOW
    elif os.path.isdir(XMUSKIE):
        loc = 'Muskie'
        topdir = XMUSKIE
    else:
        print('Error: I do not know where I am: %s' % os.getcwd())
        return ''

    print('We are on: %s' % loc)
    return topdir


def get_spec(xtab, spectrograph_id):
    '''
    Select good science fibers for a given spectrograph from a SLITMAP table.

    Filters to rows where telescope == 'Sci', fibstatus == 0, and
    spectrographid == spectrograph_id.

    Parameters:
        xtab (astropy.table.Table): SLITMAP table from an SFrame file.
        spectrograph_id (int): Spectrograph identifier (1, 2, or 3).

    Returns:
        astropy.table.Table: Filtered table of matching fibers.
    '''
    good = xtab[xtab['telescope'] == 'Sci']
    good = good[good['fibstatus'] == 0]
    good = good[good['spectrographid'] == spectrograph_id]
    return good


def get_all_spec_specs(filename, percentile=50):
    '''
    Compute percentile flux spectra for each spectrograph from one SFrame file.

    Opens the SFrame FITS file and checks whether science fibers are
    present for each of the three spectrographs (spectrographid 1, 2, 3).
    If any spectrograph has no valid science fibers the function returns
    None for all arrays together with a status dictionary that records
    which spectrographs were found.

    Parameters:
        filename (str): Path to the SFrame FITS file.
        percentile (int): Percentile to compute across fibers (default 50).

    Returns:
        tuple: (wav, flux1, flux2, flux3, status) where wav is the
        wavelength array, flux1/flux2/flux3 are 1-D percentile spectra
        for spectrographs 1/2/3, and status is a dict with keys 'SP1',
        'SP2', 'SP3' set to True if fibers were found and False otherwise.
        If the file cannot be opened, returns (None, None, None, None, None).
        If any spectrograph is missing, returns (None, None, None, None, status).
    '''
    try:
        x = fits.open(filename)
    except Exception:
        print('get_all_spec_specs: Could not open %s' % filename)
        return None, None, None, None, None

    xtab = Table(x['SLITMAP'].data)
    wav = x['WAVE'].data

    status = {'SP1': False, 'SP2': False, 'SP3': False}
    fiber_sets = {}
    for sp_id, key in zip([1, 2, 3], ['SP1', 'SP2', 'SP3']):
        fibers = get_spec(xtab, sp_id)
        if len(fibers) > 0:
            status[key] = True
            fiber_sets[sp_id] = fibers
        else:
            status[key] = False

    if not all(status.values()):
        x.close()
        return None, None, None, None, status

    flux_out = []
    for sp_id in [1, 2, 3]:
        fibers = fiber_sets[sp_id]
        sp_flux = x['FLUX'].data[fibers['fiberid'] - 1]
        sp_mask = x['MASK'].data[fibers['fiberid'] - 1]
        sp_flux = np.ma.masked_array(sp_flux, sp_mask)

        if percentile == 50:
            flux_percentile = np.ma.median(sp_flux, axis=0)
        else:
            sp_flux = np.ma.filled(sp_flux, np.nan)
            flux_percentile = np.nanpercentile(sp_flux, percentile, axis=0)

        flux_out.append(flux_percentile)

    x.close()
    return wav, flux_out[0], flux_out[1], flux_out[2], status


def make_spec_specs(xtab, data_dir, outfile='', percentile=50, file_type='CFrame'):
    '''
    Process multiple CFrame or SFrame files and write a FITS summary by spectrograph.

    Iterates over the exposures in xtab, calling get_all_spec_specs for
    each file that exists on disk.  Exposures for which not all three
    spectrographs have valid science fibers are skipped and recorded in a
    separate ASCII table written alongside the main output file.

    The drpall location column records SFrame paths.  When file_type is
    'CFrame' (the default), 'SFrame' is replaced with 'CFrame' in the path
    before checking for the file.

    The main FITS file contains extensions WAVE, FLUX1, FLUX2, FLUX3, and
    drp_all.  FLUX1/FLUX2/FLUX3 are 2-D arrays (n_exposures x n_wavelengths)
    containing the percentile spectrum for spectrographs 1, 2, and 3.

    The skipped-exposure table is always written and contains one row per
    rejected exposure with columns expnum, mjd, tileid, SP1, SP2, SP3
    showing which spectrographs were present (True) or absent (False).

    Parameters:
        xtab (astropy.table.Table): drpall table of selected exposures.
        data_dir (str): Top-level data directory prepended to the location
            column of xtab to form full file paths.
        outfile (str): Output filename root (default 'test_spec.fits').
        percentile (int): Percentile to compute across fibers (default 50).
        file_type (str): 'CFrame' (default) or 'SFrame'.

    Returns:
        None
    '''
    # Identify files that exist on disk
    select_idx = []
    xfiles = []
    for i in range(len(xtab)):
        xfile = '%s/%s' % (data_dir, xtab['location'][i])
        if file_type == 'CFrame':
            xfile = xfile.replace('SFrame', 'CFrame')
        if os.path.isfile(xfile):
            select_idx.append(i)
            xfiles.append(xfile)
    print('There are %d files to process' % len(select_idx))

    if len(select_idx) == 0:
        print('Warning: No files found to process. Check exposure range and data directory.')
        return

    xtab = xtab[select_idx]

    xflux1 = []
    xflux2 = []
    xflux3 = []
    wav = None
    good_idx = []

    skipped_expnum = []
    skipped_mjd = []
    skipped_tileid = []
    skipped_sp1 = []
    skipped_sp2 = []
    skipped_sp3 = []

    for i in range(len(xfiles)):
        xwav, flux1, flux2, flux3, status = get_all_spec_specs(xfiles[i], percentile)

        if xwav is None and status is None:
            # File could not be opened — skip silently
            continue

        if flux1 is None:
            # One or more spectrographs missing — record in skipped table
            row = xtab[i]
            skipped_expnum.append(int(row['expnum']))
            skipped_mjd.append(int(row['mjd']))
            skipped_tileid.append(int(row['tileid']))
            skipped_sp1.append(status['SP1'])
            skipped_sp2.append(status['SP2'])
            skipped_sp3.append(status['SP3'])
            print('  Skipped expnum %d: SP1=%s SP2=%s SP3=%s' % (
                row['expnum'], status['SP1'], status['SP2'], status['SP3']))
            continue

        wav = xwav
        xflux1.append(flux1)
        xflux2.append(flux2)
        xflux3.append(flux3)
        good_idx.append(i)

        if i % 10 == 0:
            print('Finished %d of %d' % (i, len(xfiles)))

    # Always write the skipped table
    if outfile == '' or outfile == 'test_spec.fits':
        skipped_file = 'test_spec.skipped.txt'
    else:
        skipped_file = outfile.replace('.fits', '') + '.skipped.txt'

    skipped_tab = Table([skipped_expnum, skipped_mjd, skipped_tileid,
                         skipped_sp1, skipped_sp2, skipped_sp3],
                        names=['expnum', 'mjd', 'tileid', 'SP1', 'SP2', 'SP3'])
    if len(skipped_tab) > 0:
        skipped_tab.write(skipped_file, format='ascii.fixed_width_two_line', overwrite=True)
        print('Skipped exposure table written to %s (%d exposures)' % (skipped_file, len(skipped_tab)))
    else:
        print('No exposures were skipped.')

    if wav is None or len(xflux1) == 0:
        print('Warning: No valid data extracted from any files.')
        return

    wav = np.array(wav)
    xflux1 = np.array(xflux1)
    xflux2 = np.array(xflux2)
    xflux3 = np.array(xflux3)

    good_tab = xtab[good_idx]

    print('SP1 shape: %s  SP2 shape: %s  SP3 shape: %s' % (
        xflux1.shape, xflux2.shape, xflux3.shape))

    # Build WCS for spectral axis
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [1, 1]
    wcs.wcs.crval = [wav[0], 0]
    wcs.wcs.cdelt = [0.5, 1]
    wcs.wcs.ctype = ['WAVE', 'LINE']

    hdu1 = fits.PrimaryHDU(data=None)
    hdu1.header['Title'] = '%s_Spec_Summary' % file_type
    hdu1.header['PERCENT'] = (percentile, 'Percentile used for flux')
    hdu1.header['SP1'] = ('spectrographid=1', 'Spectrograph 1 science fibers')
    hdu1.header['SP2'] = ('spectrographid=2', 'Spectrograph 2 science fibers')
    hdu1.header['SP3'] = ('spectrographid=3', 'Spectrograph 3 science fibers')

    hdu2 = fits.ImageHDU(data=wav, name='WAVE')
    hdu3 = fits.ImageHDU(data=xflux1, name='FLUX1')
    hdu4 = fits.ImageHDU(data=xflux2, name='FLUX2')
    hdu5 = fits.ImageHDU(data=xflux3, name='FLUX3')
    hdu6 = fits.BinTableHDU(good_tab, name='drp_all')

    hdu3.header.update(wcs.to_header())
    hdu4.header.update(wcs.to_header())
    hdu5.header.update(wcs.to_header())

    hdul = fits.HDUList([hdu1, hdu2, hdu3, hdu4, hdu5, hdu6])

    if outfile == '':
        outfile = 'test_spec.fits'
    else:
        if not outfile.endswith('.fits'):
            outfile = outfile + '.fits'

    hdul.writeto(outfile, overwrite=True)
    print('Wrote results to %s (%d exposures)' % (outfile, len(good_tab)))


def doit(exp_start=4000, exp_stop=8000, delta=5, exp_min=900., out_name='',
         drp_ver='1.2.0', percentile=50, file_type='CFrame'):
    '''
    Top-level driver: read drpall, select exposures, compute spectrograph spectra.

    Parameters:
        exp_start (int): First exposure number to process.
        exp_stop (int): Last exposure number to process.
        delta (int): Process every Nth exposure.
        exp_min (float): Minimum exposure time in seconds.
        out_name (str): Output file root name.
        drp_ver (str): DRP version string.
        percentile (int): Percentile to compute across fibers.
        file_type (str): 'CFrame' (default) or 'SFrame'.

    Returns:
        None
    '''
    xtop = find_top()
    xtab = read_drpall(drp_ver)
    if len(xtab) == 0:
        return
    ztab = select(xtab, exp_start, exp_stop, delta, exp_min)

    if out_name == '':
        out_name = 'XSpec_%s_%s_%d_%d_%d_%d.fits' % (file_type, drp_ver, exp_start, exp_stop, delta, percentile)

    make_spec_specs(ztab, xtop, outfile=out_name, percentile=percentile, file_type=file_type)


def steer(argv):
    '''
    Parse command line arguments and run the spectrograph summary.

    Parameters:
        argv (list): Command line argument list (typically sys.argv).

    Returns:
        None
    '''
    exp_start = -1
    exp_stop = -1
    delta = -1
    exp_min = 900
    percent = 50
    out_name = ''
    ver = '1.2.0'
    file_type = 'CFrame'

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-sf':
            file_type = 'SFrame'
        elif argv[i] == '-emin':
            i += 1
            exp_min = int(argv[i])
        elif argv[i] == '-out':
            i += 1
            out_name = argv[i]
        elif argv[i] == '-ver':
            i += 1
            ver = argv[i]
        elif argv[i][:5] == '-perc':
            i += 1
            percent = eval(argv[i])
        elif argv[i][0] == '-':
            print('Unknown option: %s' % argv[i])
            return
        elif exp_start < 0:
            exp_start = int(argv[i])
        elif exp_stop < 0:
            exp_stop = int(argv[i])
        elif delta < 0:
            delta = int(argv[i])
        i += 1

    if exp_start < 0 or exp_stop < 0:
        print('Usage: SummarizeSpec.py [-sf] [-ver drp_ver] [-percent 50] [-emin 900] [-out name] exp_start exp_stop delta')
        return

    if delta < 0:
        delta = 1

    doit(exp_start, exp_stop, delta, exp_min, out_name, drp_ver=ver, percentile=percent,
         file_type=file_type)


# Next lines permit one to run the routine from the command line
if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
