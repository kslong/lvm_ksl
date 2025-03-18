#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Create an image or cube from a LVM exposure


Command line usage (if any):

    usage: kslmap.py [-no_back]  [-band filter] filename

    where image_type indicates a predefined filter to plot.  The
    currenly allowed bands are 

    * ha
    * sii
    * x

    -no_back means not to subtract background from the image

    for band x, there are aditional parameters to be entered, e.g.

    -band x 6700 6800  will produce an image of the average flux 
        between 6700 and 6800 A. For this no backgournd is
        obtained

    The routine produces and output fits file, which starts with x
    normally, but starts with z if a CFrame file is provided.

Description:  

Primary routines:

    doit

Notes:

    This version creates images with the standard orientation for ds9
                                       
History:

231216 ksl Coding begun
240303 ksl Added redshift corrections if RA and DEC are near the LMC or
    SMC
240630 ksl Modified to use ra and dec's recorded in the SLITMAP extension

'''


import os
import sys
import numpy as np
import astropy.io.fits as fits
import astropy.units as u
from astropy.wcs.utils import fit_wcs_from_points
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from scipy.interpolate import griddata
from lvm_ksl import fib2radec
from astropy.coordinates import SkyCoord




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



def calculate_percentiles(arr, percentiles):
    # Calculate specified percentiles, handling NaN values
    result = np.nanpercentile(arr, percentiles)

    return result

platescale=112.36748321030637 # Focal plane platescale in "/mm
pscale=0.01
skypscale=pscale*platescale/3600 # IFU image pixel scale in deg/pix

def make_wcs(slitab):
    '''
    Create a WCS from the ra's and dec's in 
    the slitmap extension
    '''

    tab=slitab[slitab['telescope']=='Sci']

    # This choice for xcoords and ycoords is specific to the image  created
    # and flips to get a wcs with RA in creasing to the left.
    xcoords=900-tab['xpmm']/pscale
    ycoords=900-tab['ypmm']/pscale
    stars=np.array([xcoords,ycoords])
    # print('stars\n',stars)
    ra=np.array(tab['ra'])
    dec=np.array(tab['dec'])
    radec=np.transpose([ra,dec])
    # print('radec\n',radec)
    known_coords=SkyCoord(radec,frame='icrs',unit='deg')
    xwcs=fit_wcs_from_points(xy = stars, world_coords = known_coords, projection='TAN')

    xstars=np.transpose(stars)
    # print(xstars.shape)
    ra_dec = xwcs.all_pix2world(xstars, 1)
    # print(ra_dec)
    return xwcs

def doit(filename,out_label='',wrange=[6560,6566],
         crange=None,do_astrom=True,do_mask=True,interpolate=False):
    
    
    print('\nReading file',filename)
    try:
        rss=fits.open(filename)
    except:
        print('Error": could not read row stack spectra at %s' % filename)
        return
    mjd=rss['PRIMARY'].header['MJD']
    expnum=rss['PRIMARY'].header['EXPOSURE']

    xstart='x'
    if filename.count('lvmCFrame'):
        xstart='z'


    try:
        RAobs  = rss['PRIMARY'].header['POSCIRA']
    except:
        RAobs  = rss['PRIMARY'].header['TESCIRA']

    try:
        DECobs = rss['PRIMARY'].header['POSCIDE']
    except:
        DECobs = rss['PRIMARY'].header['TESCIDE']

    try:
        posang = rss['PRIMARY'].header['POSCIPA']
    except:
        print('Error: POSCIPA is missing, assuming 0')
        posang=0
    
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
    # print('Selected',sciflux.shape[0],'science fibers')
    
    
    # select wavelength range of interest
    z=get_redshift(RAobs,DECobs)
    wrange[0]=float(wrange[0])
    wrange[1]=float(wrange[1])
    wrange[0]*=(1+z)
    wrange[1]*=(1+z)
    if crange!=None:
        crange[0]=float(crange[0])
        crange[1]=float(crange[1])
        crange[0]*=(1+z)
        crange[1]*=(1+z)


    nfibers,nchans = rss['FLUX'].data.shape
    wcsrss=WCS(rss['FLUX'].header)
    wave=np.array(wcsrss.pixel_to_world(np.arange(nchans),0)[0]*1e10)
    selwave=(wave>=wrange[0])*(wave<=wrange[1])
    # print('Requested wavelength range is',wrange[0],'to',wrange[1])
    # print(f'Actual wavelength range is {wave[selwave][0]} to {wave[selwave][-1]}')
    if crange:
        # print('Requested continuum wavelength range is',crange[0],'to',crange[1])
        cselwave=(wave>=crange[0])*(wave<=crange[1])
        # print('Continuum wavelength range is',crange[0],'to',crange[1])


    # Fill the flux array
    if do_mask:
        sciflux[scimask==1] = np.nan
    flux = np.nanmean(sciflux*selwave, axis=1)
    # print(flux)
    print('Averages',np.nanmean(flux),np.nanmedian(flux))
    # Optional continuum subtraction
    if crange:
        cflux = np.nanmean(sciflux*cselwave,axis=1)
        print('Cont Averages',np.nanmean(cflux),np.nanmedian(cflux))
        flux  = flux - cflux
        print('Final Averages',np.nanmean(flux),np.nanmedian(flux))

    xmean=np.nanmean(flux)
    xmed=np.nanmedian(flux)
    xstd=np.nanstd(flux)

    percentiles=calculate_percentiles(flux,[5,95])

    print('Final stats %.2e %.2e %.2e ' % (xmean,xmed,xstd))
    # print('percentiles ',percentiles)


    # Calculate image geometry
    platescale=112.36748321030637 # Focal plane platescale in "/mm
    # print("Using Fiducial Platescale = ", platescale)
    pscale=0.01 # IFU image pixel scale in mm/pix
    rspaxel=35.3/platescale/2 # spaxel radius in mm assuming 35.3" diameter chromium mask
    npix=int(1800) # size of IFU image
    # print(f'The output image is {npix} pixels square')

    xima, yima = np.meshgrid(np.arange(npix)-npix/2,np.arange(npix)-npix/2)
    xima=xima*pscale # x coordinate in mm of each pixel in image
    yima=yima*pscale # y coordinate in mm of each pixel in image

    # print('Making image')
    ima=np.full((npix,npix),np.nan)
    if not interpolate:
        for i in range(len(x)):
            sel=(xima-x[i])**2+(yima-y[i])**2<=rspaxel**2
            ima[sel]=flux[i]
    else:
        ima = griddata((x,y), flux, (xima,yima), method='linear')
                
    # Make the header from the ra's and dec's in the slitmap if possible
    try:
        ww=make_wcs(slittab)
        header = ww.to_header()
    except:
        print('Error: could not create wcs from RA and Dec in SLITTAB, so using header keywords')
        # Create WCS for IFU image from keywords in the header
        w = WCS(naxis=2)
        w.wcs.crpix = [int(npix/2)+1, int(npix/2)+1]
        skypscale=pscale*platescale/3600 # IFU image pixel scale in deg/pix
        # The extra 180 degrees is to get the image in the corect orientation.
        posangrad=(posang+180.)*np.pi/180
        w.wcs.cd=np.array([[skypscale*np.cos(posangrad),    -1*skypscale*np.sin(posangrad)],
                       [-1*skypscale*np.sin(posangrad), -1*skypscale*np.cos(posangrad)]])
        w.wcs.crval = [RAobs,DECobs]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        header=w.to_header

    
    # Flip axes so that one produces an image with the standard orientation

    ima= np.flipud(np.fliplr(ima))

    # Save IFU image as fits
    hdu = fits.PrimaryHDU(ima.astype(np.float32), header=header)
    hdul=fits.HDUList([hdu])

    if os.path.isdir('./qdata/')==False:
        os.mkdir('./qdata/')

    if out_label=='':
            outfile='./qdata/%s%05d_%04d_%04d.fits' % (xstart,expnum,wrange[0],wrange[1])
    else:
            outfile='./qdata/%s%05d_%0s.fits' % (xstart,expnum,out_label)


    hdul.writeto(outfile, overwrite=True)

    # hdul.info()

    print('Suggested LoadFrame command')
    print('LoadFrame %s 1 %.2e %.2e' % (outfile,percentiles[0],percentiles[1]))
    return outfile

def steer(argv):
    '''
    Get the inputs

    '''
    image_type='ha'
    xfits=''
    sub_back=True

    xha=['ha',[6560,6566],[6590,6630]]
    xs2=['sii',[6710,6735],[6740,6760]]
    xxx=['x',[8000,9000],[6000,7000]]


    i=1
    while i<len(argv):
        print(argv[i])
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-band':
            i+=1
            image_type=argv[i]
            if image_type=='x':
                i+=1
                xxx[1][0]=eval(argv[i])
                i+=1
                xxx[1][1]=eval(argv[i])
                sub_back=False
        elif argv[i]=='-no_back':
            sub_back=False
        elif argv[i][0]=='-':
            print('Error: Unknown switch %s ' % argv[i])
            return
        elif xfits=='':
            xfits=argv[i]
        # print(argv[i])
        i+=1

    if xfits=='':
        print('Error: not enough arguments: ', argv)
        return

    xbands=[xha,xs2,xxx]

    good=False
    for one_band in xbands:
        # print(image_type,one_band[0])
        if image_type==one_band[0]:
            # print(one_band[1],one_band[2],image_type)
            if sub_back:
                doit(xfits,out_label=image_type,wrange=one_band[1],crange=one_band[2])
            else:
                doit(xfits,out_label=image_type,wrange=one_band[1],crange=None)
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
