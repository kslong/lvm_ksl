#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Create a summary of SFrame (sky-subtracted) data where each row contains the
median spectrum (across science fibers) from a single exposure. Useful for
evaluating sky subtraction quality and comparing residuals across many exposures.

Command line usage (if any):

    usage: SummarizeSFrame.py [-h] [-ver drp_ver] [-percent 50] [-emin 900] [-out whatever]
                              [-by pixel|fiber] [-navg 10] [-sigma 3.0] [-maxiters 5]
                              [-mask FILE] exp_start exp_stop delta

Description:

    This script processes multiple SFrame files and computes a summary spectrum
    across all science fibers for each exposure. The output contains one row per
    exposure, with the sky-subtracted flux, the sky spectrum, and the inverse
    variance.

    This is the SFrame equivalent of SummarizeCframe.py. While SummarizeCframe
    works on CFrame files (before sky subtraction) and includes SKY_EAST/SKY_WEST,
    this script works on SFrame files (after sky subtraction) and includes the
    combined SKY and IVAR extensions.

    Options: -h prints this documentation; -ver drp_ver sets the DRP version to
    use (default 1.2.1); -percent N sets the percentile to use (default 50);
    -emin sets minimum exposure time to include (default 900); -out whatever
    sets the name or root name of output fits file.

    -by pixel|fiber selects how the summary spectrum is formed (default pixel):
    'pixel' computes the percentile independently at each wavelength pixel across
    fibers, as before -- the result is a per-pixel statistical composite, not any
    single fiber's real spectrum. 'fiber' instead ranks whole science fibers by
    sky-line-masked continuum flux, then averages the -navg fibers nearest the
    -percent rank via a sigma-clipped mean, applying that same fiber window to
    FLUX, SKY, IVAR, and LSF so the output reflects a real, consistent set of
    fibers. In fiber mode: -navg N sets the number of nearest-rank fibers
    combined (default 10); -sigma S and -maxiters K set the sigma-clipping
    threshold and iteration limit for the robust mean (defaults 3.0, 5); -mask
    FILE gives the sky-line mask (default: sky_mask.fits searched in cwd then
    data/). Auto-generated output filenames get a '_fiber' suffix in fiber mode;
    an explicit -out name is used as given in either mode.

    Positional arguments: exp_start is the starting exposure number to consider;
    exp_stop is the exposure number to stop on; delta skips every Nth exposure.

Primary routines:

    doit

Notes:
                                       
History:

240726 ksl Coding begun

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
from GetSkyCont import load_mask, _interp_mask_to_wave
from SummarizeCframe import _rank_window, _robust_mean


from astropy.coordinates import SkyCoord,  Galactocentric
import astropy.units as u


import re


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before a line consisting of "History:" (or
    "History::"/"Version History" -- whitespace/colon-insensitive), so
    -h stays short even as that section grows -- without hand-
    duplicating the Synopsis/Options text in a second string.  Anchored
    to a whole line (not a bare substring search) so it can't misfire on
    "History:" appearing mid-sentence, and returns doc unchanged if no
    such line is present.
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


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


def read_drpall(drp_ver='1.2.1'):
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
        print('Error: Locared but could not read  : ', xfile)
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
    # print(len(xtab))
    # print(np.unique(xtab['exptime'],return_counts=True))
    if exp_min>0:
        xtab=xtab[xtab['exptime']>=exp_min]
    # print(len(xtab))
    # print(np.unique(xtab['exptime'],return_counts=True))
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
    elif ox.path.isdir(XMUSKIE):
        loc='Muskie'
        topdir=XMUSKIE
    else:
        print('Error: I donot know where I am:', os.pwd())
        return ''

    print('We am on : ', loc)
    return topdir

        

def scifib(xtab,select='all',telescope=''):
    '''
    Select good fibers from a telescope, of a spefic
    type or all from a telescope from the slitmap table
    of a calbrated file
    '''
    # print(np.unique(xtab['fibstatus']))
    # print(np.unique(xtab['targettype']))
    ztab=xtab[xtab['fibstatus']==0]
    if select=='all':
        ztab=ztab[ztab['targettype']!='standard']
    else:
        ztab=ztab[ztab['targettype']==select]

    if telescope!='' and telescope!='all':
        ztab=ztab[ztab['telescope']==telescope]


    # print('Found %d fibers' % len(ztab))
    return ztab


def get_med_spec(filename= '/Users/long/Projects/lvm_data/sas/sdsswork/lvm/spectro/redux/1.0.3/0011XX/11111/60192/lvmSFrame-00004336.fits',percentile=50):
    try:
        x=fits.open(filename)
    except:
        print('gets_spec: Could not open %s' % filename)
        return

    xtab=Table(x['SLITMAP'].data)

    science_fibers=scifib(xtab,select='science',telescope='Sci')
    # skye_fibers=scifib(xtab,select='SKY',telescope='SkyE')
    # skyw_fibers=scifib(xtab,select='SKY',telescope='SkyW')

    wav=x['WAVE'].data
    sci_flux=x['FLUX'].data[science_fibers['fiberid']-1]
    sci_sky=x['SKY'].data[science_fibers['fiberid']-1]
    sci_var=x['IVAR'].data[science_fibers['fiberid']-1]
    sci_lsf=x['LSF'].data[science_fibers['fiberid']-1]

    sci_mask=x['MASK'].data[science_fibers['fiberid']-1]
    sci_flux=np.ma.masked_array(sci_flux,sci_mask)
    sci_sky=np.ma.masked_array(sci_sky,sci_mask)
    sci_var=np.ma.masked_array(sci_var,sci_mask)
    sci_lsf=np.ma.masked_array(sci_lsf,sci_mask)

    # foo=np.ma.median(sci_sky,axis=0)
    # print(foo.shape)
    # print(foo)
    if percentile==50:
        # print('using median')
        sci_flux_med=np.ma.median(sci_flux,axis=0)
        sci_sky_med=np.ma.median(sci_sky,axis=0)
        sci_lsf_med=np.ma.median(sci_lsf,axis=0)
    else:
        # print('Using percentile %d' % percentile)
        sci_flux = np.ma.filled(sci_flux, np.nan)
        sci_sky  = np.ma.filled(sci_sky, np.nan)
        sci_lsf  = np.ma.filled(sci_lsf, np.nan)
        sci_flux_med=np.nanpercentile(sci_flux,percentile,axis=0)
        sci_sky_med=np.nanpercentile(sci_sky,percentile,axis=0)
        sci_lsf_med=np.nanpercentile(sci_lsf,percentile,axis=0)
        # print(sci_flux_med.shape)
        # print(sci_flux_med)

    # The inverse variance has to be multipled by the number of science fibers/(1.253)**2 to get the 
    # inverse variance of the mediane
    sci_var_med=np.ma.median(sci_var,axis=0)*len(science_fibers['fiberid'])*0.63694


    # skye_flux= x['FLUX'].data[skye_fibers['fiberid']-1]
    # skye_sky=x['SKY'].data[skye_fibers['fiberid']-1]
    # skye_mask=x['MASK'].data[skye_fibers['fiberid']-1]
    # skye_flux=np.ma.masked_array(skye_flux,skye_mask)
    # skye_sky=np.ma.masked_array(skye_sky,skye_mask)
    # skye_flux_med=np.ma.median(skye_flux,axis=0)
    # skye_sky_med=np.ma.median(skye_sky,axis=0)

    # skyw_flux= x['FLUX'].data[skyw_fibers['fiberid']-1]
    # skyw_sky=x['SKY'].data[skyw_fibers['fiberid']-1]
    # skyw_mask=x['MASK'].data[skyw_fibers['fiberid']-1]
    # skyw_flux=np.ma.masked_array(skyw_flux,skyw_mask)
    # skyw_sky=np.ma.masked_array(skyw_sky,skyw_mask)
    # skyw_flux_med=np.ma.median(skyw_flux,axis=0)
    # skyw_sky_med=np.ma.median(skyw_sky,axis=0)
    return wav, sci_flux_med, sci_sky_med, sci_var_med, sci_lsf_med


def get_fiber_spec(filename, percent=50, navg=10, sigma=3.0, maxiters=5,
                   mask_wave=None, mask_bool=None, stat='median'):
    '''
    Fiber-based (not per-pixel) analogue of get_med_spec(): rank science
    fibers by sky-line-masked continuum flux, then combine the -navg
    fibers nearest the target percentile rank via a sigma-clipped mean.

    Same ranking metric and combination method as SummarizeCframe.py's
    get_fiber_spec(), but pulls FLUX/SKY/IVAR/LSF (SFrame extensions)
    instead of FLUX/SKY_EAST/SKY_WEST/LSF.  IVAR is combined with the
    same sigma-clipped-mean window as the other arrays; the
    "* n_fibers * 0.63694" inverse-variance-of-the-median correction
    used in get_med_spec()'s percentile==50 branch does not apply here
    and is not attempted.

    Returns (wav, sci_flux, sci_sky, sci_var, sci_lsf, meta), or None if
    the file could not be processed.
    '''
    try:
        x = fits.open(filename)
    except Exception:
        print('get_fiber_spec: Could not open %s' % filename)
        return None

    try:
        xtab = Table(x['SLITMAP'].data)
        sci = scifib(xtab, select='science', telescope='Sci')
        if len(sci) < 10:
            print('get_fiber_spec: only %d science fibers found in %s, skipping.'
                 % (len(sci), filename))
            return None

        wav      = x['WAVE'].data.astype(np.float64)
        sci_flux = x['FLUX'].data[sci['fiberid'] - 1].astype(np.float64)
        sci_sky  = x['SKY'].data[sci['fiberid'] - 1].astype(np.float64)
        sci_var  = x['IVAR'].data[sci['fiberid'] - 1].astype(np.float64)
        sci_lsf  = x['LSF'].data[sci['fiberid'] - 1].astype(np.float64)
        bad = x['MASK'].data[sci['fiberid'] - 1] != 0
        sci_flux[bad] = np.nan
        sci_sky[bad]  = np.nan
        sci_var[bad]  = np.nan
        sci_lsf[bad]  = np.nan

        clean = _interp_mask_to_wave(mask_wave, mask_bool, wav)
        if stat == 'mean':
            cont = np.nanmean(sci_flux[:, clean], axis=1)
        else:
            cont = np.nanmedian(sci_flux[:, clean], axis=1)

        good = np.isfinite(cont)
        if good.sum() < 10:
            print('get_fiber_spec: too few fibers with valid continuum flux in %s, skipping.'
                 % filename)
            return None
        sci_tab  = sci[good]
        cont     = cont[good]
        sci_flux = sci_flux[good]
        sci_sky  = sci_sky[good]
        sci_var  = sci_var[good]
        sci_lsf  = sci_lsf[good]

        order    = np.argsort(cont)
        n        = len(order)
        i_target = int(round(percent / 100.0 * (n - 1)))
        win      = _rank_window(order, i_target, navg)

        sci_flux_out = _robust_mean(sci_flux[win], sigma=sigma, maxiters=maxiters)
        sci_sky_out  = _robust_mean(sci_sky[win],  sigma=sigma, maxiters=maxiters)
        sci_var_out  = _robust_mean(sci_var[win],  sigma=sigma, maxiters=maxiters)
        sci_lsf_out  = _robust_mean(sci_lsf[win],  sigma=sigma, maxiters=maxiters)

        def _mode_int(arr):
            arr = np.asarray(arr, int)
            return int(np.bincount(arr).argmax())

        meta = dict(
            n_sci_fibers         = n,
            n_avg                = len(win),
            fiberid_list         = ','.join(str(v) for v in sci_tab['fiberid'][win]),
            ra_fiber             = float(np.mean(sci_tab['ra'][win])),
            dec_fiber            = float(np.mean(sci_tab['dec'][win])),
            spectrographid_fiber = _mode_int(sci_tab['spectrographid'][win]),
            contflux_fiber       = float(np.mean(cont[win])),
        )

        return wav, sci_flux_out, sci_sky_out, sci_var_out, sci_lsf_out, meta
    finally:
        x.close()


def make_med_spec(xtab,data_dir,outfile='',percentile=50,exp_start=None,
                  exp_stop=None,delta=None,exp_min=None,drp_ver=None,
                  by='pixel',navg=10,sigma=3.0,maxiters=5,mask_file=''):
    i=0
    select_idx=[]
    xfiles=[]
    while i < len(xtab):
        xfile='%s/%s' % (data_dir,xtab['location'][i])
        if os.path.isfile(xfile):
            select_idx.append(i)
            xfiles.append(xfile)
        i+=1
    print('There are %d files to process' % (len(select_idx)))

    xtab=xtab[select_idx]
    # print(xfiles)

    mask_wave = mask_bool = None
    if by=='fiber':
        mask_wave, mask_bool = load_mask(mask_file)

    i=0
    xsci_flux=[]
    xsci_sky=[]
    xsci_var=[]
    xsci_lsf=[]
    meta_list=[]
    good_rows=[]
    while i<len(xfiles):
        if by=='fiber':
            result=get_fiber_spec(xfiles[i],percent=percentile,navg=navg,
                                  sigma=sigma,maxiters=maxiters,
                                  mask_wave=mask_wave,mask_bool=mask_bool)
            if result is None:
                i+=1
                continue
            wav,sci_flux,sci_sky,sci_var,sci_lsf,meta=result
            meta_list.append(meta)
        else:
            wav,sci_flux,sci_sky,sci_var,sci_lsf=get_med_spec(xfiles[i],percentile)
        xsci_flux.append(sci_flux)
        xsci_sky.append(sci_sky)
        xsci_var.append(sci_var)
        xsci_lsf.append(sci_lsf)
        good_rows.append(i)
        if i%10==0:
            print('Finished %d of %d' % (i,len(xfiles)))

        i+=1

    if len(xsci_flux)==0:
        print('Warning: No valid data extracted from any files.')
        return

    if by=='fiber':
        xtab=xtab[good_rows]
        for col in meta_list[0]:
            xtab[col]=[m[col] for m in meta_list]

    # dtype=np.float32 here too: get_fiber_spec() upconverts wav to float64
    # internally for mask-interpolation precision, but that shouldn't leak
    # into the output -- pixel mode's wav is float32 already.
    wav=np.array(wav,dtype=np.float32)
    # dtype=np.float32 guards against np.nanpercentile/sigma_clipped_stats
    # silently upcasting to float64 (numpy 1.26 does this for any
    # percentile computation, even though the FITS extensions are float32);
    # IVAR needs it separately since the "* 0.63694" factor in get_med_spec()
    # is a Python float and upcasts regardless of -percent.
    xsci_flux=np.array(xsci_flux,dtype=np.float32)
    xsci_sky=np.array(xsci_sky,dtype=np.float32)
    xsci_var=np.array(xsci_var,dtype=np.float32)
    xsci_lsf=np.array(xsci_lsf,dtype=np.float32)
    print(xsci_flux.shape,xsci_sky.shape)
    hdu1 = fits.PrimaryHDU(data=None)
    hdu1.header['Title'] = 'SFrame_Suummary'
    hdu1.header['ROUTINE'] = ('SummarizeSframe', 'Script that produced this file')
    hdu1.header['PERCENT'] = (percentile, 'Percentile used for flux')
    hdu1.header['DRPVER'] = (drp_ver, 'DRP version used for drpall lookup')
    hdu1.header['EMIN'] = (exp_min, 'Minimum exposure time (s)')
    hdu1.header['EXPSTART'] = (exp_start, 'First exposure number selected')
    hdu1.header['EXPSTOP'] = (exp_stop, 'Last exposure number selected')
    hdu1.header['DELTA'] = (delta, 'Exposure-number stride')
    hdu1.header['SELECT'] = (by, 'pixel or fiber selection mode')
    if by=='fiber':
        hdu1.header['NAVG'] = (navg, 'Nearest-rank fibers combined per exposure')
        hdu1.header['SIGCLIP'] = (sigma, 'Sigma-clipping threshold for robust mean')
        hdu1.header['MAXITER'] = (maxiters, 'Sigma-clipping iteration limit')
        hdu1.header['MASKFILE'] = (os.path.basename(mask_file), 'Sky-line mask used for fiber ranking')
    hdu2= fits.ImageHDU(data=wav,name='WAVE')
    hdu3=fits.ImageHDU(data=xsci_flux,name='FLUX')
    hdu4=fits.ImageHDU(data=xsci_sky,name='SKY')
    hdu5=fits.ImageHDU(data=xsci_var,name='IVAR')
    hdu5a=fits.ImageHDU(data=xsci_lsf,name='LSF')
    hdu6 = fits.BinTableHDU(xtab, name='drp_all')

    wmin=wav[0]
    dwave=0.5

    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [1, 1]  # Reference pixel (1-based index)
    wcs.wcs.crval = [wmin, 0]  # Coordinates at the reference pixel: minimum wavelength and line number 0
    wcs.wcs.cdelt = [dwave, 1]  # Pixel scale: 0.5 Angstroms per pixel in wavelength, 1 per pixel in line number
    wcs.wcs.ctype = ['WAVE', 'LINE']  # Coordinate types: wavelength and line number

    hdu3.header.update(wcs.to_header())
    hdu4.header.update(wcs.to_header())
    hdu5.header.update(wcs.to_header())
    hdu5a.header.update(wcs.to_header())

    hdul = fits.HDUList([hdu1, hdu2, hdu3,hdu4,hdu5,hdu5a,hdu6])

    if outfile=='':
        outfile='test.fits'
    else:
        if outfile.count('.fits')==0:
            outfile=outfile+'.fits'
    hdul.writeto(outfile,overwrite=True)
    print('Wrote results to %s' % outfile)
    return



def doit(exp_start=4000,exp_stop=8000,delta=5,exp_min=900.,out_name='',drp_ver='1.2.1',
        percentile=50,by='pixel',navg=10,sigma=3.0,maxiters=5,mask_file=''):
    xtop=find_top()
    xtab=read_drpall(drp_ver)
    ztab=select(xtab,exp_start,exp_stop,delta)

    if by=='fiber' and not mask_file:
        _data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
        for _candidate in [os.path.join(os.getcwd(), 'sky_mask.fits'),
                           os.path.join(_data_dir,   'sky_mask.fits')]:
            if os.path.exists(_candidate):
                mask_file = _candidate
                print('Using default mask: %s' % mask_file)
                break
        if not mask_file:
            print('Error: -by fiber requires a mask file and sky_mask.fits was not '
                 'found in the current directory or data/')
            return

    if out_name=='':
        out_name='XSFrame_%s_%d_%d_%d_%d.fits' % (drp_ver,exp_start,exp_stop,delta,percentile)
        if by=='fiber':
            out_name=out_name.replace('.fits','_fiber.fits')
    make_med_spec(xtab=ztab,data_dir=xtop,outfile=out_name,percentile=percentile,
                 exp_start=exp_start,exp_stop=exp_stop,delta=delta,
                 exp_min=exp_min,drp_ver=drp_ver,by=by,navg=navg,sigma=sigma,
                 maxiters=maxiters,mask_file=mask_file)

def steer(argv):
    '''
    SummarizeSFrame.py  exp_start expstop delta
    '''
    exp_start=-1
    exp_stop=-1
    delta=-1

    exp_min=900
    percent=50
    out_name=''

    ver='1.2.1'
    by='pixel'
    navg=10
    sigma=3.0
    maxiters=5
    mask_file=''

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(_usage_from_doc(__doc__))
            return
        elif argv[i]=='-emin':
            i+=1
            exp_min=int(argv[i])
        elif argv[i]=='-out':
            i+=1
            out_name=(argv[i])
        elif argv[i]=='-ver':
            i+=1
            ver=argv[i]
        elif argv[i][:5]=='-perc':
            i+=1
            percent=eval(argv[i])
        elif argv[i]=='-by':
            i+=1
            by=argv[i]
        elif argv[i]=='-navg':
            i+=1
            navg=int(argv[i])
        elif argv[i]=='-sigma':
            i+=1
            sigma=float(argv[i])
        elif argv[i]=='-maxiters':
            i+=1
            maxiters=int(argv[i])
        elif argv[i]=='-mask':
            i+=1
            mask_file=argv[i]
        elif argv[i][0]=='-':
            print('Unknown option : ',argv)
        elif exp_start<0:
            exp_start=int(argv[i])
        elif exp_stop<0:
            exp_stop=int(argv[i])
        elif delta<0:
            delta=int(argv[i])
        i+=1

    if delta<0:
        delta=1

    if by not in ('pixel','fiber'):
        print('Error: -by must be pixel or fiber')
        return

    doit(exp_start,exp_stop,delta,exp_min,out_name,drp_ver=ver,percentile=percent,
        by=by,navg=navg,sigma=sigma,maxiters=maxiters,mask_file=mask_file)




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
