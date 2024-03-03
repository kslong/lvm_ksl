#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Create an image or cube from a LVM exposure


Command line usage (if any):

    usage: kslmap.py [-no_back]  [-image_type filter] filename

    where image_type indicates a predefined filter to plot.  The
    currenly allowed bands are 

    * ha
    * sii

    -no_back means not to subtract background from the image

Description:  

Primary routines:

    doit

Notes:

    This version creates images with the standard orientation for ds9
                                       
History:

231216 ksl Coding begun
w40303 ksl Added redshift corrections if RA and DEC are near the LMC or
    SMC

'''


import os
import sys
import numpy as np
import astropy.io.fits as fits
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from scipy.interpolate import griddata
from lvm_ksl import fib2radec




RSS_DIR='./data'
RADIAN=57.29578

def distance(r1,d1,r2,d2):
    '''
    distance(r1,d1,r2,d2)
    Return the angular offset between two ra,dec positions
    All variables are expected to be in degrees.
    Output is in degrees

    Note - This routine could easily be made more general
    '''
#    print 'distance',r1,d1,r2,d2
    r1=r1/RADIAN
    d1=d1/RADIAN
    r2=r2/RADIAN
    d2=d2/RADIAN
    xlambda=np.sin(d1)*np.sin(d2)+np.cos(d1)*np.cos(d2)*np.cos(r1-r2)
#    print 'xlambda ',xlambda
    if xlambda>=1.0:
        xlambda=0.0
    else:
        xlambda=np.arccos(xlambda)

    xlambda=xlambda*RADIAN
#    print 'angle ',xlambda
    return xlambda

def get_redshift(r,d):
    '''
    Get a redshift correction if the object is close to the 
    SMC or LMC
    '''

    lmc=[80.8942,-69.7561,6.,262.]
    smc=[13.1583,-72.8003,3.,146]

    if distance(r,d,lmc[0],lmc[1])<lmc[2]:
        print('Field is in the LMC')
        return lmc[3]/3e5
    elif distance(r,d,smc[0],smc[1])<smc[2]:
        print('Field is in the SMC')
        return smc[3]/3e5
    
    print('No redshift correct applied')
    return 0.0



def doit(filename,wrange=[6560,6566], out_label='',
         crange=None,do_astrom=True,cube=False,do_mask=True,interpolate=False):
    
    
    # Read RSS QL file
    # xname='%s/%s' % (RSS_DIR,filename)
    xname=fib2radec.locate_file(filename)
    print('\nReading file',xname)
    try:
        rss=fits.open(xname)
    except:
        print('Error": could not read row stack spectra at %s' % xname)
        return
    mjd=rss['PRIMARY'].header['MJD']
    expnum=rss['PRIMARY'].header['EXPOSURE']


    posang=0
    RAobs=0
    DECobs=0

    try:
        best=rss['PRIMARY'].header['Best']
        print('Found Best astrometry')
        posang=rss['PRIMARY'].header['BESTPA']
        RAobs=rss['PRIMARY'].header['BESTRA']
        DECobs=rss['PRIMARY'].header['BESTDEC']
        do_astrom=False

    except:
        print('The is no pre-existing astrometry')

    if do_astrom:
        LVMAGCAM_DIR = os.environ.get('LVMAGCAM_DIR')
        procscifile = f"{LVMAGCAM_DIR}/{mjd}/coadds/lvm.sci.coadd_s{expnum:0>8}.fits"
        if os.path.exists(procscifile):
            print('Found astrometry file')
            
            # Read WCS info from coadd
            agfile    = os.path.basename(procscifile)
            mfheader  = fits.getheader(procscifile, 1)
            outw      = WCS(mfheader)
            CDmatrix  = outw.pixel_scale_matrix
            posangrad = -1*np.arctan(CDmatrix[1,0]/CDmatrix[0,0])
            posang    = posangrad*180/np.pi
            IFUcen    = outw.pixel_to_world(2500,1000)
            RAobs     = IFUcen.ra.value
            DECobs    = IFUcen.dec.value
            print(f'Using center and PA of {RAobs:.6f}, {DECobs:.6f}, {posang:.6f} from {agfile}')
            # For recent data, the header keywords should be correct
            altRAobs  = mfheader['RAMEAS']
            altDECobs = mfheader['DECMEAS']
            altposang = mfheader['PAMEAS'] - 180.
            print(f'Header values are {altRAobs:.6f} and {altDECobs:.6f} with PA {altposang:.6f}')
        else:
            print(f'{procscifile} does not exist, skipping astrometry')
            do_astrom = False
            RAobs  = mfheader['RAMEAS']
            DECobs = mfheader['DECMEAS']
            posang = mfheader['POSCIPA']
    
    # Read fibermap and get x,y coordinates of fibers
    slittab = rss['SLITMAP'].data
    targettype=slittab['targettype']
    spectrograph=slittab['spectrographid']
    telescope=slittab['telescope']
    
    
    selsci=(telescope=='Sci')
    x=slittab['xpmm'][selsci]
    y=slittab['ypmm'][selsci]
    fibid=slittab['fiberid'][selsci]
    sciflux = rss['FLUX'].data[selsci]
    scimask = rss['MASK'].data[selsci]
    print('Selected',sciflux.shape[0],'science fibers')
    
    
    # select wavelength range of interest
    z=get_redshift(RAobs,DECobs)
    wrange[0]*=(1+z)
    wrange[1]*=(1+z)
    if crange!=None:
        crange[0]*=(1+z)
        wrange[1]*=(1+z)


    nfibers,nchans = rss['FLUX'].data.shape
    wcsrss=WCS(rss['FLUX'].header)
    wave=np.array(wcsrss.pixel_to_world(np.arange(nchans),0)[0]*1e10)
    selwave=(wave>=wrange[0])*(wave<=wrange[1])
    print('Requested wavelength range is',wrange[0],'to',wrange[1])
    print(f'Actual wavelength range is {wave[selwave][0]} to {wave[selwave][-1]}')
    if crange:
        print('Requested ontinuum wavelength range is',crange[0],'to',crange[1])
        cselwave=(wave>=crange[0])*(wave<=crange[1])
        print('Continuum wavelength range is',crange[0],'to',crange[1])


    # Fill the flux array
    if do_mask:
        sciflux[scimask==1] = np.nan
    if not cube:
        flux = np.nanmean(sciflux*selwave, axis=1)
        print(flux)
        print('Averages',np.nanmean(flux),np.nanmedian(flux))
        # Optional continuum subtraction
        if crange:
            cflux = np.nanmean(sciflux*cselwave,axis=1)
            print('Cont Averages',np.nanmean(cflux),np.nanmedian(cflux))
            flux  = flux - cflux
            print('Final Averages',np.nanmean(flux),np.nanmedian(flux))
    else:  # In my version, never reduce the number of channels
        flux = sciflux[:,selwave]
        crval3 = wave[selwave.argmax()]
        cdelt3 = rss['FLUX'].header['CDELT1'] * 1e10  # Angstroms

    xmean=np.nanmean(flux)
    xmed=np.nanmedian(flux)
    xstd=np.nanstd(flux)
    print('Final stats %.2e %.2e %.2e ' % (xmean,xmed,xstd))

    # Calculate image geometry
    platescale=112.36748321030637 # Focal plane platescale in "/mm
    print("Using Fiducial Platescale = ", platescale)
    pscale=0.01 # IFU image pixel scale in mm/pix
    rspaxel=35.3/platescale/2 # spaxel radius in mm assuming 35.3" diameter chromium mask
    npix=int(1800) # size of IFU image
    print(f'The output image is {npix} pixels square')

    # Create IFU Image
    xima, yima = np.meshgrid(np.arange(npix)-npix/2,np.arange(npix)-npix/2)
    xima=xima*pscale # x coordinate in mm of each pixel in image
    yima=yima*pscale # y coordinate in mm of each pixel in image
    if not cube:
        print('Making image')
        ima=np.full((npix,npix),np.nan)
        if not interpolate:
            for i in range(len(x)):
                sel=(xima-x[i])**2+(yima-y[i])**2<=rspaxel**2
                ima[sel]=flux[i]
        else:
            ima = griddata((x,y), flux, (xima,yima), method='linear')
    else:
        print('Making cube')
        nwave = flux.shape[1]
        ima=np.full((nwave,npix,npix),np.nan)
        if not interpolate:
            for i in range(len(x)):
                sel=(xima-x[i])**2+(yima-y[i])**2<=rspaxel**2
                for k in range(nwave):
                    ima[k, sel]=flux[i, k]
        else:
            for k in range(nwave):
                ima[k,:] = griddata((x,y), flux[:,k], (xima,yima), method='linear')
                
    # Create WCS for IFU image
    w = WCS(naxis=2)
    w.wcs.crpix = [int(npix/2)+1, int(npix/2)+1]
    skypscale=pscale*platescale/3600 # IFU image pixel scale in deg/pix
    # The extra 180 degrees is to get the image in the corect orientation.
    posangrad=(posang+180.)*np.pi/180
    w.wcs.cd=np.array([[skypscale*np.cos(posangrad),    -1*skypscale*np.sin(posangrad)],
                       [-1*skypscale*np.sin(posangrad), -1*skypscale*np.cos(posangrad)]])
    w.wcs.crval = [RAobs,DECobs]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    header = w.to_header()
    
    
    # Flip axes so that one produces an image with the standard orientation

    ima= np.flipud(np.fliplr(ima))

    # Save IFU image as fits
    if cube == False:
        hdu = fits.PrimaryHDU(ima.astype(np.float32), header=header)
    else:
        header['WCSAXES']=3
        header['CRVAL3']=crval3
        header['CDELT3']=cdelt3
        header['CUNIT3']='Angstrom'
        header['CRPIX3']=1.
        header['CTYPE3']='WAVE'
        hdu = fits.PrimaryHDU(ima.astype(np.float32), header=header)
            
    hdul=fits.HDUList([hdu])

    if out_label=='':
        if cube==False:
            outfile='x%05d_%04d_%04d.fits' % (expnum,wrange[0],wrange[1])
        else:
            outfile='c%05d_%04d_%04d.fits' % (expnum,wrange[0],wrange[1])
    else:
        if cube==False:
            outfile='x%05d_%s.fits' % (expnum,out_label)
        else:
            outfile='c%05d_%0s.fits' % (expnum,out_label)

    if crange==None:
        outfile=outfile.replace('.fits','_no_back.fits')

    hdul.writeto(outfile, overwrite=True)

    hdul.info()

    print('Suggested LoadFrame commdand')
    print('LoadFreme %s 1 %.2e %.2e' % (outfile,xmed-xstd,xmed+3*xstd))

def steer(argv):
    '''
    Get the inputs

    '''
    image_type='ha'
    xfits=''
    sub_back=True

    i=1
    while i<len(argv):
        print(argv[i])
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-band':
            i+=1
            image_type=argv[i]
        elif argv[i]=='-no_back':
            sub_back=False
        elif argv[i][0]=='-':
            print('Error: Unknown switch %s ' % argv[i])
            return
        elif xfits=='':
            xfits=argv[i]
        print(argv[i])
        i+=1


    xha=['ha',[6560,6566],[6590,6630]]
    xs2=['sii',[6710,6735],[6740,6760]]

    if xfits=='':
        print('Error: not enough arguments: ', argv)
        return

    xbands=[xha,xs2]

    good=False
    for one_band in xbands:
        if image_type==one_band[0]:
            print(one_band[1],one_band[2],image_type)
            if sub_back:
                doit(xfits,wrange=one_band[1],crange=one_band[2],out_label=image_type)
            else:
                doit(xfits,wrange=one_band[1],crange=None,out_label=image_type)
            good=True
            break
        if good==False:
            print('Error: The filter %s did not match a known band' % image_type)
            print(__doc__)

        return

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
