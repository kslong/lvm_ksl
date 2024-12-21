#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Generate a fits file condaining an image of an RSS file
over a limited wavelength range

It should handle CFrame or SFrame data, but was
written to handle frames writtne with rss_combine


Command line usage (if any):

    usage: rss2image.py [-no_back]  [-band filter] filename

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



Description:  

Primary routines:

    do_one

Notes:
                                       
History:

241130 ksl Coding begun.  The routine was developed from the
    routine quickmap.py which only handles SFrame or CFrame
    files.  

'''

# # Create a new version of QuickMaP that is more like what I want to do with BigMaps
import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt
import os
import astropy.io.fits as fits
import astropy.units as u
from astropy.wcs.utils import fit_wcs_from_points
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from scipy.interpolate import griddata
from astropy.coordinates import SkyCoord
from astropy.table import Table,join,vstack
import numpy as np

RADIAN=57.29578

def distance(r1,d1,r2,d2):
    '''
    distance(r1,d1,r2,d2)
    Return the angular offset between two ra,dec positions
    All variables are expected to be in degrees.
    Output is in degrees

    Note - This routine could easily be made more general
    '''
    r1=r1/RADIAN
    d1=d1/RADIAN
    r2=r2/RADIAN
    d2=d2/RADIAN
    xlambda=np.sin(d1)*np.sin(d2)+np.cos(d1)*np.cos(d2)*np.cos(r1-r2)
    if xlambda>=1.0:
        xlambda=0.0
    else:
        xlambda=np.arccos(xlambda)

    xlambda=xlambda*RADIAN
    return xlambda


def get_redshift_correction(r,d):
    '''
    Get a redshift correction if the object is close to the
    SMC or LMC

    If there is no correction this returns 1., so this should
    be multiplied by a wavelength
    '''

    lmc=[80.8942,-69.7561,6.,262.]
    smc=[13.1583,-72.8003,3.,146]

    if distance(r,d,lmc[0],lmc[1])<lmc[2]:
        print('Field is in the LMC')
        return 1.+lmc[3]/3e5
    elif distance(r,d,smc[0],smc[1])<smc[2]:
        print('Field is in the SMC')
        return 1.+smc[3]/3e5

    print('No redshift correct applied')
    return 1.0

def calculate_percentiles(arr, percentiles):
    # Calculate specified percentiles, handling NaN values
    result = np.nanpercentile(arr, percentiles)

    return result

def clean_slitmap(slit_tab):
    '''
    Clean slit_tab and add a column for the row number
    '''
    nstart=len(slit_tab)
    slit_tab['Row']=np.arange(nstart)
    colnames=slit_tab.colnames
    for one in colnames:
        if one=='telescope':
            slit_tab=slit_tab[slit_tab['telescope']=='Sci']
        if one=='fibstatus':
            slit_tab=slit_tab[slit_tab['fibstatus']==0]

    nstop=len(slit_tab)
    print('clean_slitmap: %d --> %d so %d lost' % (nstart,nstop,nstart-nstop))
    
    return slit_tab



def create_wcs_manually(ra_center,dec_center,scale,size):
    '''
    Manually create a header, where scale is the pixel
    sale in arrsec, and size is the size in degrees
    '''



    # print('manual',ra_center,dec_center,scale,size)
    pixel_scale = scale / 3600.
    center=size/2*3600
    wcs = WCS(naxis=2)  # Define a 2D WCS object
    wcs.wcs.crpix = [center, center]  # Reference pixel
    wcs.wcs.crval = [ra_center, dec_center]  # Reference world coordinates (RA, Dec)
    wcs.wcs.cdelt = [-pixel_scale, pixel_scale]  # Pixel scale: negative for RA (increasing to the east), positive for Dec
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']  # Tangential projection for RA and Dec
    # wcs.naxis = [size, size] 
    # print(wcs)
    return wcs


def create_wcs_from_points(xpos,ypos,ra_deg,dec_deg):
    '''
    Create an wcs, from a series of points, and RA's a Dec's
    where scale represents to conversion 
    '''
    pix_coords=(xpos,ypos)
    sky_coords = SkyCoord(ra=ra_deg, dec=dec_deg, unit="deg", frame="icrs")
    wcs = fit_wcs_from_points(pix_coords, sky_coords, proj_point="center")
    # Fit the WCS using fit_wcs_from_points
    try:
        wcs = fit_wcs_from_points(pix_coords, sky_coords, proj_point="center")
        print("WCS successfully generated:")
        print(wcs)
        success=True
    except Exception as e:
        print(f"Failed to generate WCS: {e}")
        success=False

    if success==True:
        return wcs
    else:
        return Null






def get_flux(rss,wmin=6000,wmax=6300,xtype='sum',ext='FLUX'):
    '''
    get the mean, madian, or sum of the 'flux'in 
    a wavelenght interval  a 'standard rss file.
    This returns values for all of the spectra, and so can
    include rows of poor fibers and from all telescopes

    This was originally set always to get the FLUS
    but can obtain any extension that exits

    A table containg the row number, and the value
    '''
    wave=rss['WAVE'].data
    wtab=Table([wave,np.arange(len(wave))],names=['WAVE','i'])
    wtab=wtab[wtab['WAVE']>wmin]
    wtab=wtab[wtab['WAVE']<wmax]
    imin=wtab['i'][0]
    imax=wtab['i'][-1]
    try:
        zarray=rss[ext].data
    except:
        print('Could not obatin data from extension %s' % (ext))
        raise IOError
    farray=zarray[:,imin:imax]
    if xtype=='sum':
        flux=np.nansum(farray,axis=1)
    elif xtype=='med':
        flux=np.nanmedian(farray,axis=1)
    else:
        flux=np.nanmean(farray,axis=1)

    ztab=Table([np.arange(len(flux)),flux],names=['Row','FLUX'])
    nwave=imax-imin
    return nwave,ztab



def get_wcs(header):
    try:
        wcs = WCS(header)
        if wcs.has_celestial:
            print("The header contains a valid celestial WCS.")
            print(wcs)
            return wcs
        else:
            print("The header has WCS, but it does not define celestial coordinates.")
    except Exception as e:
        print(f"The header does not contain a valid WCS. Error: {e}")
    return Null




def do_one(filename='data/lvmSFrame-00007373.fits',wrange=[6550,6570],crange=None,outroot='',xscale=1,xsize=0.5,pos_type='slit',ext='FLUX'):
    '''
    This is the primary driving routine
    '''


    x=fits.open(filename)

    phead=x['PRIMARY'].header

    xtab=Table(x['SLITMAP'].data)
    xtab['Row']=np.arange(len(xtab))

    qtab=clean_slitmap(xtab)

    # print('After clean ',len(qtab))

    # Find out if the X and Y columns exist

    slit_names=qtab.colnames

    have_wcs=False
    try:
        xpos=qtab['X']
        ypos=qtab['Y']
        wcs=create_wcs_from_points(xpos,ypos,qtab['ra'],qtab['dec'])
        print('Boundaries ',np.min(xpos), np.max(xpos))
        print('Boundaries ',np.min(ypos), np.max(ypos))
        have_wcs=True
    except:
        print('Could not generate WCS using for slittab')

    if have_wcs==False:
        try:
            xra=phead['POSCIRA']
            xdec=phead['POSCIDE']
            wcs=create_wcs_manually(xra,xdec,xscale,xsize)
            print('Generated WCS from primary header')
            have_wcs=True
        except:
            print('Could not generate wcs from primary header :')

    if have_wcs==False:
        print('There is no point in contuing, so returning')


    xra=wcs.wcs.crval[0]
    xdec=wcs.wcs.crval[1]

    redshift_correction=get_redshift_correction(xra,xdec)

    # print(wrange)

    wmin=wrange[0]
    wmax=wrange[1]

    wmin*=redshift_correction
    wmax*=redshift_correction

    print('# wmin and wmax for image: %.1f %.1f' % (wmin,wmax))

    try:
        if ext=='FLUX':
            nwave,xflux=get_flux(x,wmin,wmax,xtype='sum',ext=ext)
        else:
            nwave,xflux=get_flux(x,wmin,wmax,xtype='med',ext=ext)
    except:
        print('Failed to retrieve data')
        return


    if crange!=None:
        cmin=crange[0]
        cmax=crange[1]
        print('# wmin and wmax for background: %.1f %.1f' % (cmin,cmax))
        nback,xback=get_flux(x,cmin,cmax,xtype='med')
        xflux['FLUX']-=nwave*xback['FLUX']

    xmean=np.nanmean(xflux['FLUX'])
    xmed=np.nanmedian(xflux['FLUX'])
    xstd=np.nanstd(xflux['FLUX'])

    percentiles=calculate_percentiles(xflux['FLUX'],[5,95])

    qtab=join(qtab,xflux,keys='Row',join_type='left')

    # So at this point I have the final WCS, and I have the fluxes, I need to calculate the xy postions of the fibers that I still have in terms of the WCS that I have generated
    print('This is the WCS we are using')
    print(wcs)
    x_pixels, y_pixels = wcs.world_to_pixel_values(qtab['ra'], qtab['dec'])
    qtab['X']=x_pixels
    qtab['Y']=y_pixels

    
    print('Reproject X:  %8.1f  %8.1f' % (np.min(qtab['X']),np.max(qtab['X'])))
    print('Reproject Y:  %8.1f  %8.1f' % (np.min(qtab['Y']),np.max(qtab['Y'])) )

    # Now we need to generate the final image

    npix=int(xsize*3600/xscale)
    xmax=np.max(qtab['X'])
    ymax=np.max(qtab['Y'])
    npix=int(np.max([xmax,ymax]))
    print('Making image with %d rows/columns' %  (npix))
    
    # xima, yima = np.meshgrid(np.arange(npix)-npix/2,np.arange(npix)-npix/2)
    xima, yima = np.meshgrid(np.arange(npix),np.arange(npix))

    ima=np.full((npix,npix),np.nan)
    rspaxel=35./2./xscale
    # print('npix, rspaxel ',npix,rspaxel)
    for one_row in qtab:
        sel=(xima-one_row['X'])**2+(yima-one_row['Y'])**2<=rspaxel**2
        ima[sel]=one_row['FLUX']

    # Now we havwe the image, and can create the ouput fits file

    primary_hdu = fits.PrimaryHDU(header=phead)
    hdu1 = fits.ImageHDU(data=ima, header=wcs.to_header(), name="IMAGE1")
    hdulist = fits.HDUList([primary_hdu])
    hdulist.append(hdu1)

    if outroot=='':
        outroot='test'
    output_file=outroot.replace('.fits','')+'.fits'
    if ext!='FLUX':
        output_file=output_file.replace('.fits','.%s.fits' % ext)
    hdulist.writeto(output_file, overwrite=True)
    print('\nFinal stats: mean %.2e median %.2e STD %.2e ' % (xmean,xmed,xstd))

    print('\nSuggested LoadFrame command')
    print('LoadFrame %s 1 %.2e %.2e' % (output_file,percentiles[0],percentiles[1]))
    return

    
def steer(argv):
    '''
    Get the inputs

    '''
    image_type='ha'
    xfits=''
    sub_back=True
    ext='FLUX'

    xha=['ha',[6560,6566],[6590,6630]]
    xs2=['sii',[6710,6735],[6740,6760]]
    xxx=['x',[8000,9000],[6000,7000]]

    files=[]


    i=1
    while i<len(argv):
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
        elif argv[i]=='-ext':
            i+=1
            ext=argv[i]
        elif argv[i]=='-no_back':
            sub_back=False
        elif argv[i][0]=='-':
            print('Error: Unknown switch %s ' % argv[i])
            return
        elif argv[i].count('fits'):
            files.append(argv[i])
        else:
            print('Error: Could not interpret command line :',argv)
        i+=1

    if len(files)==0:
        print('Error: not enough arguments: ', argv)
        return

    xbands=[xha,xs2,xxx]
    # xbands contains 3 possible intervals, and nomially I chose one of thes

    found_band=False
    for one_band in xbands:
        if image_type==one_band[0]:
            qband=one_band
            good=found_band=True
            break
    if found_band==False:
        print('It is not obvious a defined band exists for image_type %s' % (image_type))
        return

            

    for xfits in files:
        words=xfits.split('/')
        froot=words[-1].replace('.fits','')
        froot=froot.replace('lvmSFrame-000','')
        xoutroot='image_%s_%s' % (froot,qband[0])


        if sub_back | (ext !='FLUX'):
            do_one(xfits,wrange=qband[1],crange=qband[2],outroot=xoutroot,xscale=1,xsize=0.5,pos_type='slit',ext=ext)
        else:
            xoutroot='%s_no_back' % xoutroot
            do_one(xfits,wrange=qband[1],crange=None,outroot=xoutroot,xscale=1,xsize=0.5,pos_type='slit',ext=ext)

    return









# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)   
    else:
        print (__doc__)
