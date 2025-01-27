#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Look at how well skysubtraction is working on a large
number of files


Command line usage (if any):

    usage: SummarizeSFrame.py  [-h] [-ver drp_ver][-percent 50] [-emin 900] [-out whatever] exp_start expstop delta

Description:  

    where:

    -h prints this documentatikon
    -ver drp_ver where drp_ver is a version run of the drp, e.g the current default 1.1.1
    -percent is the percentile to use in the calculation, if missing use median or 50
    -emin sets the mimimum exposure time to include; the default is 900
    -out is name or rootname of output fits file
    exp_start is the number is the exposure to start with
    exp_stop is the number of the exposure to stop with
    delta is the number of exposures to skip in creating the output image

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

    sci_mask=x['MASK'].data[science_fibers['fiberid']-1]
    sci_flux=np.ma.masked_array(sci_flux,sci_mask)
    sci_sky=np.ma.masked_array(sci_sky,sci_mask)
    sci_var=np.ma.masked_array(sci_var,sci_mask)

    # foo=np.ma.median(sci_sky,axis=0)
    # print(foo.shape)
    # print(foo)
    if percentile==50:
        # print('using median')
        sci_flux_med=np.ma.median(sci_flux,axis=0)
        sci_sky_med=np.ma.median(sci_sky,axis=0)
    else:
        # print('Using percentile %d' % percentile)
        sci_flux = np.ma.filled(sci_flux, np.nan)
        sci_sky  = np.ma.filled(sci_sky, np.nan)
        sci_flux_med=np.nanpercentile(sci_flux,percentile,axis=0)
        sci_sky_med=np.nanpercentile(sci_sky,percentile,axis=0)
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
    return wav, sci_flux_med, sci_sky_med, sci_var_med


def make_med_spec(xtab,data_dir,outfile='',percentile=50):
    i=0
    select=[]
    xfiles=[]
    while i < len(xtab):
        xfile='%s/%s' % (data_dir,xtab['location'][i])
        if os.path.isfile(xfile):
            select.append(i)
            xfiles.append(xfile)
        i+=1
    print('There are %d files to process' % (len(select)))

    xtab=xtab[select]
    # print(xfiles)
    i=0
    xsci_flux=[]
    xsci_sky=[]
    xsci_var=[]
    while i<len(xfiles):
        wav,sci_flux,sci_sky,sci_var=get_med_spec(xfiles[i],percentile)
        xsci_flux.append(sci_flux)
        xsci_sky.append(sci_sky)
        xsci_var.append(sci_var)
        if i%10==0:
            print('Finished %d of %d' % (i,len(xfiles)))
                                          
        i+=1
    wav=np.array(wav)
    xsci_flux=np.array(xsci_flux)
    xsci_sky=np.array(xsci_sky)
    print(xsci_flux.shape,xsci_sky.shape)
    hdu1 = fits.PrimaryHDU(data=None)
    hdu1.header['Title'] = 'SFrame_Suummary'
    hdu2= fits.ImageHDU(data=wav,name='WAVE')
    hdu3=fits.ImageHDU(data=xsci_flux,name='FLUX')
    hdu4=fits.ImageHDU(data=xsci_sky,name='SKY')
    hdu5=fits.ImageHDU(data=xsci_var,name='IVAR')
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

    hdul = fits.HDUList([hdu1, hdu2, hdu3,hdu4,hdu5,hdu6])

    if outfile=='':
        outfile='test.fits'
    else:
        if outfile.count('.fits')==0:
            outfile=outfile+'.fits'
    hdul.writeto(outfile,overwrite=True)
    print('Wrote results to %s' % outfile)
    return



def doit(exp_start=4000,exp_stop=8000,delta=5,exp_min=900.,out_name='',drp_ver='1.1.0',percentile=50):
    xtop=find_top()
    xtab=read_drpall(drp_ver)
    ztab=select(xtab,exp_start,exp_stop,delta)
    
    if out_name=='':
        out_name='XSFrame_%s_%d_%d_%d_%d.fits' % (drp_ver,exp_start,exp_stop,delta,percentile)
    make_med_spec(xtab=ztab,data_dir=xtop,outfile=out_name,percentile=percentile)

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

    ver='1.1.0'

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-emin':
            i+=1
            exp_min=int(argv[i])
        elif argv[i]=='-out':
            i+=1
            out_name=(argv[i])
        elif argv[i]=='-ver']:
            i+=1
            ver=argv[i]
        elif argv[i][:5]=='-perc':
            i+=1
            percent=eval(argv[i])
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
                

    doit(exp_start,exp_stop,delta,exp_min,out_name,drp_ver=ver,percentile=percent)




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
