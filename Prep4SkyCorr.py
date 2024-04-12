#!/usr/bin/env python
# coding: utf-8


'''
                    Space Telescope Science Institute

Synopsis:  

Read one or more calibrated (but un-sky subtraced) lvmCframe and generate
files that can be used with SkyCorr.  The routine constructs 
mean spectra for the science, SkyE, and SkyW data that are in 
the frame or frames


Command line usage (if any):

    usage: Prep4SkyCorr.py [-h] [-dir whatever]  -all filename1 filename2 ..

Description:  

    where filenames is a list of files to be processed, and the optional
    parameters

    -h  will print this help, and quit
    -dir whatever will cause the directory whatever to be serched for the files that are indcatdd
    -all will cause all lvmCframe files to be processed
    -

Primary routines:

    Prep4SkyMean
    steer

Notes:
                                       
History:

240328 ksl Coding begun

'''



import os
from astropy.stats import biweight_location, biweight_scale, biweight_midvariance,mad_std
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from glob import glob




def Prep4SkyCorrSingle(filename='data/lvmCFrame-00006661.fits',fiber_id=10):
    '''
    Create a fits file that can be used as an input for skycorr
    from a single fiber
    '''
    x=fits.open(filename)
    # determine which telescope we are discussing
    slitmap=Table(x['SLITMAP'].data)
    one_slit=slitmap[slitmap['fiberid']==fiber_id]
    # print(one_slit)
    telescope=one_slit['telescope']
    print(telescope)
    
    if telescope=='Sci':
        ra=x['PRIMARY'].header['TESCIRA']
        dec=x['PRIMARY'].header['TESCIDE']
        airmass=x['PRIMARY'].header['TESCIAM']
        xtel='Sci'
    elif telescope=='SkyW':
        ra=x['PRIMARY'].header['TESKYWRA']
        dec=x['PRIMARY'].header['TESKYWDE']
        airmass=x['PRIMARY'].header['TESKYWAM']
        xtel='SkyW'
    elif telescope=='SkyE':
        ra=x['PRIMARY'].header['TESKYERA']
        dec=x['PRIMARY'].header['TESKYEDE']
        airmass=x['PRIMARY'].header['TESKYEAM']  
        xtel='SkyE'
    else:
        print('Error: No telescope identified')
        return
    
    # print(ra,dec,airmass)
    alt=90-np.arccos(1./airmass)*57.29578
    # print(alt)
    
    xtime=x['PRIMARY'].header['OBSTIME']
    print(xtime)
    word=xtime.split('T')
    word=word[-1].split(':')
    utc=eval(word[0])+eval(word[1])/60.+eval(word[2])/3600
    # print(utc)
    
    wave=x['WAVE'].data
    xflux=flux=x['FlUX'].data[fiber_id-1,:]
    xsky=x['SKY'].data[fiber_id-1,:]
    flux+=xsky
    error=x['ERROR'].data[fiber_id-1,:]
    # print(error.shape)
    
    xtab=Table([wave,flux,error,xflux,xsky], names=['WAVE','FLUX','ERROR','XFLUX','XSKY'])
    table_hdu = fits.BinTableHDU(xtab, name='') 
    
    new_primary_hdu = fits.PrimaryHDU(header=x[0].header)
    new = fits.HDUList([new_primary_hdu, table_hdu])
    
    new.info()
    
    new['Primary'].header['UTC']=utc
    new['Primary'].header['RA']=ra
    new['Primary'].header['DEC']=dec
    new['Primary'].header['ALT']=alt
    new['Primary'].header['TELE']=xtel
    
    
    words=filename.split('-')
    goo=words[-1].replace('.fits','')
    # print(goo)
    outfile='w%s_%04d.fits' % (goo,fiber_id)
    print('Writing: %s' % outfile)
    new.writeto(outfile,overwrite=True)
    
   


def scifib(xtab,select='all',telescope=''):
    # print(np.unique(xtab['fibstatus']))
    # print(np.unique(xtab['targettype']))
    ztab=xtab[xtab['fibstatus']==0]
    if select=='all':
        ztab=ztab[ztab['targettype']!='standard']
    else:
        ztab=ztab[ztab['targettype']==select]

    if telescope!='' and telescope!='all':
        ztab=ztab[ztab['telescope']==telescope]

    print('Found %d fibers' % len(ztab))
    return ztab


def get_spec(filename='NoSky/lvmCFrame-00004852.fits',select='science',telescope='',extension='FLUX'):
    '''
    '''
    try:
        xfits=fits.open(filename)
    except:
        print('Error: eval_qual: Could not open %s' % filename)
        return
    xtab=Table(xfits['SLITMAP'].data)
    xgood=scifib(xtab,select=select,telescope=telescope)
    wav=xfits['WAVE'].data
    xflux=xfits[extension].data[xgood['fiberid']-1]
    flux=biweight_location(xflux,axis=0,ignore_nan=True)
    std=mad_std(xflux,axis=0,ignore_nan=True)


    results=Table([wav,flux,std],names=['WAVE','FLUX','STD'])
    results['WAVE'].format='.2f'
    results['FLUX'].format='.4e'
    results['STD'].format='.4e'
    expo=xfits[0].header['EXPOSURE']
    outname='Spec_%05d' % expo
    if telescope!='':
        outname='%s_%s' % (outname,telescope)
    if extension!='FLUX':
        outname='%s_%s' % (outname,extension)

    results.write('%s.txt' % outname,format='ascii.fixed_width_two_line',overwrite=True)

    xfits.close()
    return outname








def Prep4SkyCorrMean(filename='data/lvmCFrame-00006661.nosky_sub.fits'):
    '''
    Create mean spectra for the science, skyE and skyW fibers
    in a calibrated data frame and write results to 
    fits files that can be used by SkyCorr
    '''
    
    try:
        x=fits.open(filename)
    except:
        print('Error: Could not open %s' % filename)
        return
    #
    xroot=x['PRIMARY'].header['EXPOSURE']
    
    # determine which telescope we are discussing
    slitmap=Table(x['SLITMAP'].data)     
    
    good_fibers=slitmap[slitmap['fibstatus']==0]
    
    sci=good_fibers[good_fibers['telescope']=='Sci']
    sky_e=good_fibers[good_fibers['telescope']=='SkyE']
    sky_w=good_fibers[good_fibers['telescope']=='SkyW']
    
    print(len(sci),len(sky_e),len(sky_w))
    
    xsci=x['FLUX'].data[sci['fiberid']-1]
    esci=x['Error'].data[sci['fiberid']-1]

    xsky=x['SKY'].data[sci['fiberid']-1]
    esky=x['SKY_Error'].data[sci['fiberid']-1]
    
    xsky_e=x['FLUX'].data[sky_e['fiberid']-1]
    esky_e=x['Error'].data[sky_e['fiberid']-1]
    
    xsky_w=x['FLUX'].data[sky_w['fiberid']-1]
    esky_w=x['Error'].data[sky_w['fiberid']-1]    


    # xflux=xfits[extension].data[xgood['fiberid']-1]
    # flux=biweight_location(xflux,axis=0,ignore_nan=True)
    # std=mad_std(xflux,axis=0,ignore_nan=True)


    sci_flux=biweight_location(xsci,axis=0,ignore_nan=True)
    sci_err=mad_std(xsci,axis=0,ignore_nan=True)

    sky=biweight_location(xsky,axis=0,ignore_nan=True)
    sky_e=biweight_location(xsky,axis=0,ignore_nan=True)
    
    sky_e_flux=biweight_location(xsky_e,axis=0,ignore_nan=True)
    sky_e_err=mad_std(xsky_e,axis=0,ignore_nan=True)
    
    
    sky_w_flux=biweight_location(xsky_w,axis=0,ignore_nan=True)
    sky_w_err=mad_std(xsky_w,axis=0,ignore_nan=True)

    print('a',np.sum(sky_e_flux-sky_w_flux))
    
    # sci_med=biweight_location(xsci,axis=0,ignore_nan=True))
    # esci_med=np.median(esci,axis=0)
    
    # sky_e_med=np.median(xsky_e,axis=0)
    # esky_e_med=np.median(esky_e,axis=0)
    
    # sky_w_med=np.median(xsky_w,axis=0)
    # esky_w_med=np.median(esky_w,axis=0)
    
    wave=x['WAVE'].data
    
    # plt.plot(wave,sci_med,label='Sci')
    # plt.plot(wave,sky_e_med,label='SkyE')
    # plt.plot(wave,sky_w_med,label='SkyW')
    
    # plt.plot(wave,sci_flux,label='Sci')
    # plt.plot(wave,sky_e_flux,label='SkyE')
    # plt.plot(wave,sky_w_flux,label='SkyW')
    
    # Now create the 3 sets of data
    xtime=x['PRIMARY'].header['OBSTIME']
    print(xtime)
    word=xtime.split('T')
    word=word[-1].split(':')
    utc=int(word[0])+int(word[1])/60.+eval(word[2])/3600
    
    # First do the science data
    
    ra=x['PRIMARY'].header['TESCIRA']
    dec=x['PRIMARY'].header['TESCIDE']
    airmass=x['PRIMARY'].header['TESCIAM']
    xtel='Sci'
    alt=90-np.arccos(1./airmass)*57.29578
    
    # xtab=Table([wave,sci_med,esci_med], names=['WAVE','FLUX','ERROR'])
    xtab=Table([wave,sci_flux,sci_err], names=['WAVE','FLUX','ERROR'])

    table_hdu = fits.BinTableHDU(xtab, name='') 
    new_primary_hdu = fits.PrimaryHDU(header=x[0].header)
    new = fits.HDUList([new_primary_hdu, table_hdu])
    
    new['Primary'].header['UTC']=utc
    new['Primary'].header['RA']=ra
    new['Primary'].header['DEC']=dec
    new['Primary'].header['ALT']=alt
    new['Primary'].header['TELE']=xtel
    
    Sci_file=outfile='Sci_%04d.fits' % xroot
    print('Writing: %s' % outfile)
    new.writeto(outfile,overwrite=True)      

    #now write a median super sky; header to this file is the same
    # so the rest should be easy

    xtab['FLUX'] = sky
    xtab['ERROR'] = sky_e

    table_hdu = fits.BinTableHDU(xtab, name='') 
    new_primary_hdu = fits.PrimaryHDU(header=x[0].header)
    new = fits.HDUList([new_primary_hdu, table_hdu])
    
    new['Primary'].header['UTC']=utc
    new['Primary'].header['RA']=ra
    new['Primary'].header['DEC']=dec
    new['Primary'].header['ALT']=alt
    new['Primary'].header['TELE']=xtel
    

    Sky_file=outfile='SkyM_%04d.fits' % xroot
    print('Writing: %s' % outfile)
    new.writeto(outfile,overwrite=True)      

    
    # Now do Sky E
    
    ra=x['PRIMARY'].header['TESKYERA']
    dec=x['PRIMARY'].header['TESKYEDE']
    airmass=x['PRIMARY'].header['TESKYEAM']  
    xtel='SkyE'
    alt=90-np.arccos(1./airmass)*57.29578
    
    # xtab=Table([wave,sky_e_med,esky_e_med], names=['WAVE','FLUX','ERROR'])
    xtab_sky_e=Table([wave,sky_e_flux,sky_e_err], names=['WAVE','FLUX','ERROR'])

    table_hdu = fits.BinTableHDU(xtab_sky_e, name='') 
    new_primary_hdu = fits.PrimaryHDU(header=x[0].header)
    new = fits.HDUList([new_primary_hdu, table_hdu])
    
    new['Primary'].header['UTC']=utc
    new['Primary'].header['RA']=ra
    new['Primary'].header['DEC']=dec
    new['Primary'].header['ALT']=alt
    new['Primary'].header['TELE']=xtel
    
    SkyE_file=outfile='SkyE_%04d.fits' % xroot
    print('Writing: %s' % outfile)
    new.writeto(outfile,overwrite=True)      

    # Finally do Sky W
    
    ra=x['PRIMARY'].header['TESKYWRA']
    dec=x['PRIMARY'].header['TESKYWDE']
    airmass=x['PRIMARY'].header['TESKYWAM']  
    xtel='SkyE'
    alt=90-np.arccos(1./airmass)*57.29578
    
    # xtab=Table([wave,sky_w_med,esky_w_med], names=['WAVE','FLUX','ERROR'])
    xtab_sky_w=Table([wave,sky_w_flux,sky_w_err], names=['WAVE','FLUX','ERROR'])

    print('b',np.sum(xtab_sky_e['FLUX']-xtab_sky_w['FLUX']))

    table_hdu = fits.BinTableHDU(xtab_sky_w, name='') 
    new_primary_hdu = fits.PrimaryHDU(header=x[0].header)
    new = fits.HDUList([new_primary_hdu, table_hdu])
    
    new['Primary'].header['UTC']=utc
    new['Primary'].header['RA']=ra
    new['Primary'].header['DEC']=dec
    new['Primary'].header['ALT']=alt
    new['Primary'].header['TELE']=xtel
    
    SkyW_file=outfile='SkyW_%04d.fits' % xroot
    print('Writing: %s' % outfile)
    new.writeto(outfile,overwrite=True)      

    Sci_file=Sci_file.replace('.fits','')
    SkyE_file=SkyE_file.replace('.fits','')
    SkyW_file=SkyW_file.replace('.fits','')
    return Sci_file,SkyE_file,SkyW_file
    

def steer(argv):
    '''
    This is mainly a steering routine
    '''
    xdir='.'
    files=[]
    do_all=False
    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__) 
            return
        elif argv[i]=='-dir':
            i+=1
            xdir=argv[i]
        elif argv[i]=='-all':
            do_all=True
        elif argv[i][0]=='-':
            print('Error: Problem with command line: ',argv)
            return
        else:
            files.append(argv[i])
        i+=1

    if do_all:
        files=glob('%s/*CFrame*.fits' % xdir )
    else:
        j=0
        while j<len(files):
            files[j]='%s/%s' % (xdir,files[j])
            j+=1

    if len(files)==0:
        print('Error: No files to process')
        return

    for one in files:
        Prep4SkyCorrMean(filename=one)
    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
