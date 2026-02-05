#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Create a summary of SFrame (sky-subtracted) data where each row contains the
percentile spectrum for 3 ring sets (inner, middle, outer) from a single exposure.
Useful for evaluating radial variations in sky subtraction quality.

Command line usage (if any):

    usage: SummarizeRings.py [-h] [-ver drp_ver] [-percent 50] [-emin 900] [-out whatever]
                             [-inner 1 9] [-middle 10 19] [-outer 20 25]
                             exp_start exp_stop delta

Description:

    This script processes multiple SFrame files and computes the percentile
    spectrum for three ring sets (inner, middle, outer) defined by ringnum
    in the slitmap. The output contains one row per exposure, with the
    sky-subtracted flux for each ring set and the median sky spectrum.

    Options: -h prints this documentation; -ver drp_ver sets the DRP version to
    use (default 1.2.0); -percent N sets the percentile to use (default 50 =
    median); -emin sets minimum exposure time to include (default 900); -out
    whatever sets the name or root name of output fits file; -inner min max sets
    ring range for inner set (default 1 9); -middle min max sets ring range for
    middle set (default 10 19); -outer min max sets ring range for outer set
    (default 20 25).

    Positional arguments: exp_start is the starting exposure number to consider;
    exp_stop is the exposure number to stop on; delta skips every Nth exposure.

Primary routines:

    doit

Notes:

History:

250128 ksl Coding begun, based on SummarizeSFrame.py

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.table import join, Table
import shutil
from datetime import datetime
from astropy.wcs import WCS


from astropy.coordinates import SkyCoord,  Galactocentric
import astropy.units as u

def augment_drp_all(xtab):

    drp_all=xtab
    sci_ra=drp_all['sci_ra']
    sci_dec=drp_all['sci_dec']
    skye_ra=drp_all['skye_ra']
    skye_dec=drp_all['skye_dec']
    skyw_ra=drp_all['skyw_ra']
    skyw_dec=drp_all['skyw_dec']

    sci_coord = SkyCoord(ra=sci_ra*u.degree, dec=sci_dec*u.degree, frame='icrs')
    skye_coord = SkyCoord(ra=skye_ra*u.degree, dec=skye_dec*u.degree, frame='icrs')
    skyw_coord = SkyCoord(ra=skyw_ra*u.degree, dec=skyw_dec*u.degree, frame='icrs')
    de = sci_coord.separation(skye_coord).degree
    dw = sci_coord.separation(skyw_coord).degree
    galactic_lat=sci_coord.galactic.b.degree
    lmc=262.
    smc=146.
    lmc_pos=SkyCoord(ra=80.89416666666668*u.degree,dec=-69.75618*u.degree,frame='icrs')
    lmc_sep=lmc_pos.separation(sci_coord)
    lmc_sep=lmc_sep.degree
    smc_pos=SkyCoord(ra=13.1875*u.degree,dec=-72.8286*u.degree,frame='icrs')
    smc_sep=smc_pos.separation(sci_coord)
    smc_sep=smc_sep.degree


    xlocal=np.select([np.fabs(galactic_lat)>10],['HighLat'],default='Plane')

    xlocal=np.select([lmc_sep<8],['LMC'],default=xlocal)
    xlocal=np.select([smc_sep<5],['SMC'],default=xlocal)
    drp_all['Survey']=xlocal
    near=np.select([de<dw],['SKY_EAST'],default='SKY_WEST')
    far=np.select([de>=dw],['SKY_EAST'],default='SKY_WEST')
    drp_all['Near']=near
    drp_all['Far']=far
    lmc=262.
    smc=146.
    drp_all['Redshift']=np.select([drp_all['Survey']=='LMC',drp_all['Survey']=='SMC'],[lmc,smc],default=0.0)
    return drp_all


def read_drpall(drp_ver='1.0.3'):
    DRPFILE='drpall-%s.fits' % (drp_ver)
    # First try to locate the DRP file locally, otherwise
    if os.path.isfile(DRPFILE):
        xfile=DRPFILE
    else:
        BASEDIR='/uufs/chpc.utah.edu/common/home/sdss51/sdsswork/lvm/spectro/redux/%s/' % (drp_ver)
        xfile='%s/%s' % (BASEDIR,DRPFILE)
        if os.path.isfile(xfile)==False:
            print('Error: Could not locate : ', xfile)
            return []

    print('Opening : ', xfile)
    try:
        drpall=fits.open(xfile)
    except:
        print('Error: Located but could not read  : ', xfile)
        return []

    drp_tab=Table(drpall[1].data)
    drp_tab=augment_drp_all(drp_tab)
    return  drp_tab


def select(ztab,exp_start=4000,exp_stop=8000,delta=5,exp_min=900.):
    '''
    select every nth exposure between two start and stop times
    '''
    xtab=ztab[ztab['expnum']>=exp_start]
    xtab=xtab[xtab['expnum']<=exp_stop]
    if exp_min>0:
        xtab=xtab[xtab['exptime']>=exp_min]
    if delta>1:
        xtab=xtab[::delta]
    return xtab


XTOP='/uufs/chpc.utah.edu/common/home/sdss51/'
XRAINBOW='/Users/long/Projects/lvm_data/sas'
XMUSKIE='/home/long/Projects/lvm_data/sas'

def find_top():
    if os.path.isdir(XTOP):
        loc='Utah'
        topdir=XTOP
    elif os.path.isdir(XRAINBOW):
        loc='Rainbow'
        topdir=XRAINBOW
    elif os.path.isdir(XMUSKIE):
        loc='Muskie'
        topdir=XMUSKIE
    else:
        print('Error: I do not know where I am:', os.getcwd())
        return ''

    print('We are on : ', loc)
    return topdir



def get_ring(xtab, ring_min=1, ring_max=3):
    '''
    Select good science fibers within a ring range from the slitmap table
    '''
    good = xtab[xtab['telescope'] == 'Sci']
    good = good[good['fibstatus'] == 0]
    good = good[(good['ringnum'] >= ring_min) & (good['ringnum'] <= ring_max)]
    return good


def get_ring_spec(filename, ring_min=1, ring_max=25, percentile=50):
    '''
    Get the percentile spectrum for fibers in a specified ring range

    Returns wavelength array and percentile flux spectrum
    '''
    try:
        x = fits.open(filename)
    except:
        print('get_ring_spec: Could not open %s' % filename)
        return None, None

    xtab = Table(x['SLITMAP'].data)
    ring_fibers = get_ring(xtab, ring_min=ring_min, ring_max=ring_max)

    if len(ring_fibers) == 0:
        print('get_ring_spec: No fibers found for ring %d-%d in %s' % (ring_min, ring_max, filename))
        return None, None

    wav = x['WAVE'].data
    ring_flux = x['FLUX'].data[ring_fibers['fiberid'] - 1]
    ring_mask = x['MASK'].data[ring_fibers['fiberid'] - 1]
    ring_flux = np.ma.masked_array(ring_flux, ring_mask)

    if percentile == 50:
        flux_percentile = np.ma.median(ring_flux, axis=0)
    else:
        ring_flux = np.ma.filled(ring_flux, np.nan)
        flux_percentile = np.nanpercentile(ring_flux, percentile, axis=0)

    x.close()
    return wav, flux_percentile


def get_all_ring_specs(filename, ring_sets, percentile=50):
    '''
    Get percentile spectra for all ring sets plus median sky

    Parameters:
        filename: path to SFrame file
        ring_sets: list of tuples [(min1, max1), (min2, max2), (min3, max3)]
        percentile: percentile to compute (default 50 = median)

    Returns:
        wav: wavelength array
        flux_inner: percentile flux for inner ring set
        flux_middle: percentile flux for middle ring set
        flux_outer: percentile flux for outer ring set
        sky_med: median sky spectrum across all science fibers
    '''
    try:
        x = fits.open(filename)
    except:
        print('get_all_ring_specs: Could not open %s' % filename)
        return None, None, None, None, None

    xtab = Table(x['SLITMAP'].data)
    wav = x['WAVE'].data

    # Get flux for each ring set
    flux_rings = []
    for ring_min, ring_max in ring_sets:
        ring_fibers = get_ring(xtab, ring_min=ring_min, ring_max=ring_max)

        if len(ring_fibers) == 0:
            print('get_all_ring_specs: No fibers found for ring %d-%d in %s' % (ring_min, ring_max, filename))
            flux_rings.append(np.full(len(wav), np.nan))
            continue

        ring_flux = x['FLUX'].data[ring_fibers['fiberid'] - 1]
        ring_mask = x['MASK'].data[ring_fibers['fiberid'] - 1]
        ring_flux = np.ma.masked_array(ring_flux, ring_mask)

        if percentile == 50:
            flux_percentile = np.ma.median(ring_flux, axis=0)
        else:
            ring_flux = np.ma.filled(ring_flux, np.nan)
            flux_percentile = np.nanpercentile(ring_flux, percentile, axis=0)

        flux_rings.append(flux_percentile)

    # Get median sky from all science fibers
    all_sci = get_ring(xtab, ring_min=1, ring_max=25)
    sci_sky = x['SKY'].data[all_sci['fiberid'] - 1]
    sci_mask = x['MASK'].data[all_sci['fiberid'] - 1]
    sci_sky = np.ma.masked_array(sci_sky, sci_mask)
    sky_med = np.ma.median(sci_sky, axis=0)

    x.close()
    return wav, flux_rings[0], flux_rings[1], flux_rings[2], sky_med


def make_ring_specs(xtab, data_dir, outfile='', percentile=50,
                    ring_sets=[(1, 9), (10, 19), (20, 25)]):
    '''
    Process multiple SFrame files and create output FITS with ring percentile spectra
    '''
    i = 0
    select_idx = []
    xfiles = []
    while i < len(xtab):
        xfile = '%s/%s' % (data_dir, xtab['location'][i])
        if os.path.isfile(xfile):
            select_idx.append(i)
            xfiles.append(xfile)
        i += 1
    print('There are %d files to process' % (len(select_idx)))

    if len(select_idx) == 0:
        print('Warning: No files found to process. Check exposure range and data directory.')
        return

    xtab = xtab[select_idx]

    i = 0
    xflux_inner = []
    xflux_middle = []
    xflux_outer = []
    xsky = []
    wav = None

    while i < len(xfiles):
        xwav, flux_inner, flux_middle, flux_outer, sky_med = get_all_ring_specs(
            xfiles[i], ring_sets, percentile)

        if xwav is None:
            i += 1
            continue

        wav = xwav
        xflux_inner.append(flux_inner)
        xflux_middle.append(flux_middle)
        xflux_outer.append(flux_outer)
        xsky.append(sky_med)

        if i % 10 == 0:
            print('Finished %d of %d' % (i, len(xfiles)))

        i += 1

    if wav is None or len(xflux_inner) == 0:
        print('Warning: No valid data extracted from any files.')
        return

    wav = np.array(wav)
    xflux_inner = np.array(xflux_inner)
    xflux_middle = np.array(xflux_middle)
    xflux_outer = np.array(xflux_outer)
    xsky = np.array(xsky)

    print('Inner shape: %s, Middle shape: %s, Outer shape: %s' %
          (xflux_inner.shape, xflux_middle.shape, xflux_outer.shape))

    # Create FITS file
    hdu1 = fits.PrimaryHDU(data=None)
    hdu1.header['Title'] = 'SFrame_Ring_Summary'
    hdu1.header['PERCENT'] = (percentile, 'Percentile used for flux')
    hdu1.header['INNER'] = ('%d-%d' % ring_sets[0], 'Inner ring range')
    hdu1.header['MIDDLE'] = ('%d-%d' % ring_sets[1], 'Middle ring range')
    hdu1.header['OUTER'] = ('%d-%d' % ring_sets[2], 'Outer ring range')

    hdu2 = fits.ImageHDU(data=wav, name='WAVE')
    hdu3 = fits.ImageHDU(data=xflux_inner, name='FLUX_INNER')
    hdu4 = fits.ImageHDU(data=xflux_middle, name='FLUX_MIDDLE')
    hdu5 = fits.ImageHDU(data=xflux_outer, name='FLUX_OUTER')
    hdu6 = fits.ImageHDU(data=xsky, name='SKY')
    hdu7 = fits.BinTableHDU(xtab, name='drp_all')

    wmin = wav[0]
    dwave = 0.5

    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [1, 1]
    wcs.wcs.crval = [wmin, 0]
    wcs.wcs.cdelt = [dwave, 1]
    wcs.wcs.ctype = ['WAVE', 'LINE']

    hdu3.header.update(wcs.to_header())
    hdu4.header.update(wcs.to_header())
    hdu5.header.update(wcs.to_header())
    hdu6.header.update(wcs.to_header())

    hdul = fits.HDUList([hdu1, hdu2, hdu3, hdu4, hdu5, hdu6, hdu7])

    if outfile == '':
        outfile = 'test_rings.fits'
    else:
        if outfile.count('.fits') == 0:
            outfile = outfile + '.fits'
    hdul.writeto(outfile, overwrite=True)
    print('Wrote results to %s' % outfile)
    return


def doit(exp_start=4000, exp_stop=8000, delta=5, exp_min=900., out_name='',
         drp_ver='1.1.0', percentile=50, ring_sets=[(1, 9), (10, 19), (20, 25)]):
    xtop = find_top()
    xtab = read_drpall(drp_ver)
    ztab = select(xtab, exp_start, exp_stop, delta)

    if out_name == '':
        out_name = 'XRings_%s_%d_%d_%d_%d.fits' % (drp_ver, exp_start, exp_stop, delta, percentile)
    make_ring_specs(xtab=ztab, data_dir=xtop, outfile=out_name, percentile=percentile,
                    ring_sets=ring_sets)


def steer(argv):
    '''
    SummarizeRings.py  exp_start expstop delta
    '''
    exp_start = -1
    exp_stop = -1
    delta = -1

    exp_min = 900
    percent = 50
    out_name = ''

    ver = '1.2.0'

    # Default ring sets
    inner = (1, 9)
    middle = (10, 19)
    outer = (20, 25)

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-emin':
            i += 1
            exp_min = int(argv[i])
        elif argv[i] == '-out':
            i += 1
            out_name = (argv[i])
        elif argv[i] == '-ver':
            i += 1
            ver = argv[i]
        elif argv[i][:5] == '-perc':
            i += 1
            percent = eval(argv[i])
        elif argv[i] == '-inner':
            i += 1
            ring_min = int(argv[i])
            i += 1
            ring_max = int(argv[i])
            inner = (ring_min, ring_max)
        elif argv[i] == '-middle':
            i += 1
            ring_min = int(argv[i])
            i += 1
            ring_max = int(argv[i])
            middle = (ring_min, ring_max)
        elif argv[i] == '-outer':
            i += 1
            ring_min = int(argv[i])
            i += 1
            ring_max = int(argv[i])
            outer = (ring_min, ring_max)
        elif argv[i][0] == '-':
            print('Unknown option : ', argv)
        elif exp_start < 0:
            exp_start = int(argv[i])
        elif exp_stop < 0:
            exp_stop = int(argv[i])
        elif delta < 0:
            delta = int(argv[i])
        i += 1

    if delta < 0:
        delta = 1

    ring_sets = [inner, middle, outer]

    doit(exp_start, exp_stop, delta, exp_min, out_name, drp_ver=ver,
         percentile=percent, ring_sets=ring_sets)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
