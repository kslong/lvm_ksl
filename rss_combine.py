#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Combine multiple LVM row-stacked (SFrame) spectra into a 
sincle row-staked spectrum, by creating a grid 
of fibers that cover the entire area, and then
assigning portions of the spectra of various fibers 
to the final one.  Alternatively, combine the
exposures for a single tile into one rss image with
"fibers" centered on the mean postion of the exposures.


Command line usage (if any):

    usage: rss_combine.py [-orig] [-sum] [-med] [-outroot xxxx] filenames

    where 
        -outroot xxxx   sets the name for the output file
        -orig           instead of creating set of output fibers on
                        a regular grid, here
                        the 'output' fiber positions are taken to
                        be the average of the positions of the input
                        fibers; the fractional contributions are still
                        calculated.
        -sum            instead of calculating fractional contributions
                        to the various 'output' fibers and dividing 
                        here the entire input fiber spectra  is added to one
                        of the 'output' fibers (assuming there is a
                        fiber which is near enough.  As in -orig, the
                        output fiber locations are taken from the input
                        fiber locations.  
        -med            By default the routine takes all of the
                        images, apportions the fluxes from each
                        of the images onto a new grid, and ultimately
                        averages the fluxes created.  With this
                        switch, the median is created.
                   

Description:  

    The routine sums multiple RSS images to create a final output image,
    with various options

    With no switchhes the routine reads the input SFrames, and constructs
    'output' fiber postions on a regular grid on the sky. Given the 
    positions of these fibers it apportions flux from the original images
    to the outuput fibers, based on the area of overlap.  

    With the -orig switch, the output grid is constructed from the input 
    grid, based on the average postion (ra, and dec) of each of the input
    fibers.  The original spectra are still apportion based ont he overlaps.
    This may be a better approach for a single tile than the regular grid
    approach

    With the -sum switch the positions of the rss fibers in the outuput image
    are taken as above, but additionally instead of splitting the spectra
    between output fibers, all of the flux from an imput fiber is assinged
    to one ouput fiber.  This may be the best approach if the input SFrame
    files were not dithered.

    A single fits file is produced.  

    * The primary header is simply taken from one of the imnput images
    * FLUX  - the flux is stored here, note that the number of rows will
        now depend on how many artificial fibers were created
    * IVAR  - the final inverse variance as calculated by converting the 
        inverse variance of the individual exposures to variance and
        the combining the varianance depending on what fraction of
        the flux from an input fiber was written into an output fiber.
        We did not use inverse variance weighting in
        creating the fluxes so we cannot simply sum the inverse variances.
    * MASK - Set to 0 if this fiber and waveleenght has anny valid EXPOSURE, 
        1 otherwise
    * WAVE - the standard set of wavelengtths
    * SLITMAP - a version of the slitmap file, that gives the ra's and dec's of
        the fibers, and enough other information so that the DAP should run.
        The SLITMAP also contains the X, Y positions of the fibers on
        the WCS defined below.
    * WCS_INFO - an extension that contains the WCS that contains the
        WCS that was created to redifine the fiber positions.  It is possile
        to use this to create an "image" on the sky, but one can equally ge
    * EXPOSURE - the effective number of exposures that went into
        each fiber. 

Primary routines:

    do_combine

Notes:

    The final image is not weighted by the errors.  All of the imput
    images are treated identically.

    There is not attempt to scale the images to one another globally,
    though one might argue that would be a good idea.

    It is apparent the when one uses the mean, and that exposure is
    somehow different, the effects of bad pixels are more apparent
    than when one uses the median.


                                       
History:

241201 ksl Coding begun

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from astropy.table import Table, join,vstack
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.wcs.utils import fit_wcs_from_points
from astropy.wcs import WCS
from scipy.interpolate import griddata
from astropy.coordinates import SkyCoord
import os
import psutil

def xcheck(xfiles):
    '''
    Check the input files for obvious problemss
    '''
    drp=[]
    commit=[]
    fluxcal=[]
    ra=[]
    dec=[]
    tile=[]
    obj=[]
    mjd=[]
    exptime=[]
    unit=[]
    helio=[]
    for one in xfiles:
        try:
            f=fits.open(one)
            xhead=f['PRIMARY'].header
            # print('Got ;',one)
        except:
            print('Could not open ',one)
            pass

        try:
            drp.append(xhead['drpver'])
        except:
            drp.append('Unknown')
        try:
            commit.append(xhead['COMMIT'])
        except:
            commit.append('Unknown')

        try:
            fluxcal.append(xhead['FLUXCAL'])
        except:
            fluxcal.append('Unknown')

        try:
            mjd.append(xhead['MJD'])
        except:
            mjd.append(-99)

        try:
            ra.append(xhead['POSCIRA'])
        except:
            ra.append(-37.0)

        try:
            dec.append(xhead['POSCIDE'])
        except:
            dec.append(-99.0)


        try:
            tile.append(xhead['TILE_ID'])
        except:
            tile.append('Unknown')

        try:
            obj.append(xhead['OBJECT'])
        except:
            obj.append('Unknown')

        try:
            unit.append(xhead['BUNIT'])
        except:
            unit.append('Unknown')


        try:
            exptime.append(xhead['EXPTIME'])
        except:
            exptime.append('Unknown')

        try:
            helio.append(xhead['WAVE HELIORV_SCI'])
        except:
            helio.append(-99.0)


    # print(len(xfiles),len(drp),len(fluxcal),len(ra),len(dec),len(tile),len(obj))

    xtab=Table([xfiles,mjd,drp,commit,fluxcal,unit,helio,ra,dec,tile,exptime,obj],names=['Filename','MJD','DRP','Commit','FluxCal','BUNIT','RA','HeleoV','Dec','Tile_ID','EXPTIME','Source_name'])
    xtab.sort((['Filename']))
    print(xtab)
    return xtab

def check_array_info(arr):
    array_type = "MaskedArray" if isinstance(arr, np.ma.MaskedArray) else "ndarray"
    print( {"shape": arr.shape, "dtype": arr.dtype, "type": array_type})

def check_memory(percentage=80.):
    '''
    Calculate the percentage of memory bening used and report a problem
    if it is greater than a given percentage
    '''
    mem=psutil.virtual_memory()
    frac=mem.used/mem.total*100
    if frac>percentage:
        print(f'Warning: Total Memory: {mem.total/1e9:.2f} GB')
        print(f'Warning:  Used Memory: {mem.used/1e9:.2f} GB')
        print(f'Warning %.2f of memor is being used' % (frac*100.))
    return frac


def get_size(xfiles,rad=0.25):
    '''
    From some SFRame images calculate the center
    and size of an image that would encompass them all
    '''
    ra=[]
    dec=[]
    for one_file in xfiles:
        x=fits.open(one_file)
        ra.append(x['PRIMARY'].header['TESCIRA'])
        dec.append(x['PRIMARY'].header['TESCIDE'])
    xtab=Table([xfiles,ra,dec],names=['Filename','RA','Dec'])
    ra_max=np.max(xtab['RA'])
    ra_min=np.min(xtab['RA'])
    dec_max=np.max(xtab['Dec'])
    dec_min=np.min(xtab['Dec'])

    ra_cen=0.5*(ra_max+ra_min)
    dec_cen=0.5*(dec_max+dec_min)

    delta_dec=dec_max-dec_min
    delta_ra=(ra_max-ra_min)*np.cos(dec_cen/57.29578)
    size=np.max([delta_ra,delta_dec])+2*rad
    # print(xtab)
    #old ra_cen=np.average(xtab['RA'])
    #old dec_cen=np.average(xtab['Dec'])
    # delta_dec=np.max(np.fabs(xtab['Dec']-dec_cen))
    # delta_ra=np.max(np.fabs(xtab['RA']-ra_cen))*np.cos(dec_cen/57.29578)
    # print(delta_dec,delta_ra)
    # size=2*rad+np.max([delta_ra,delta_dec])
    print('RA: %10.6f Dec: %10.6f size: %.1f' % (ra_cen,dec_cen,size))
    return ra_cen,dec_cen,size
    



def create_wcs(ra_deg, dec_deg,pos_ang, size_deg):
    """
    Create a WCS object centered on a given RA and DEC.
    
    Parameters:
        ra_deg (float): Right Ascension in degrees (center of the image).
        dec_deg (float): Declination in degrees (center of the image).
    
    Returns:
        WCS: Configured WCS object.
    """
    theta = np.deg2rad(pos_ang)  # Convert angle to radians
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
   
    # Define the WCS
    wcs = WCS(naxis=2)
    
    # Reference pixel is the center of the image
    image_size_pixels = 3600 * size_deg  # 1 degree = 3600 arcseconds
    wcs.wcs.crpix = [image_size_pixels / 2, image_size_pixels / 2]  # Center pixel

    
    
    # Set pixel scale: 1 arcsec/pixel -> 1/3600 degrees/pixel
    # wcs.wcs.cdelt = [-1/3600, 1/3600]  # Negative for RA to account for axis direction
    pixel_scale=1./3600.
    
    # Set reference coordinates (RA and DEC) at the center
    wcs.wcs.crval = [ra_deg, dec_deg]
    wcs.wcs.cd = np.array([
        [-pixel_scale * cos_theta, -pixel_scale * sin_theta],
        [pixel_scale * sin_theta,  pixel_scale * cos_theta]
        ])
    
    # Define the projection type
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]  # Tangential (gnomonic) projection
    
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
        # print(wcs)
        success=True
    except Exception as e:
        print(f"Failed to generate WCS: {e}")
        success=False

    if success==True:
        return wcs
    else:
        return Null



def generate_grid(wcs, spacing_arcsec):
    """
    Generate a regular grid of pixel positions centered on the WCS reference pixel,
    and compute the corresponding RA and Dec for each position.

    Parameters:
        wcs (astropy.wcs.WCS): The WCS object defining the coordinate system.
        spacing_arcsec (float): Spacing between grid points in arcseconds.

    Returns:
        tuple: Two 2D numpy arrays containing RA and Dec values for the grid.
    """
    # Get the WCS reference pixel and pixel scale
    crpix = wcs.wcs.crpix  # Reference pixel coordinates (1-based)

    cd_matrix = wcs.pixel_scale_matrix  # The CD matrix
    pixel_spacing_x = np.sqrt(cd_matrix[0, 0]**2 + cd_matrix[0, 1]**2)  # Scale along x-axis
    pixel_spacing_y = np.sqrt(cd_matrix[1, 0]**2 + cd_matrix[1, 1]**2)  # Scale along y-axis'

    # print(pixel_spacing_x,pixel_spacing_y) # This the pixel scale in arcsec
    
    # cdelt = wcs.wcs.cdelt  # Pixel scale in degrees
    spacing_deg = spacing_arcsec / 3600.0  # Convert spacing to degrees
    print('deg',spacing_deg)

    # Calculate pixel spacing in terms of pixels
    # pixel_spacing_x = spacing_deg / abs(cdelt[0])
    # pixel_spacing_y = spacing_deg / abs(cdelt[1])

    # Define the array size based on CRPIX
    array_size_x = int(2 * crpix[0])  # Assuming array size is 2 * CRPIX
    array_size_y = int(2 * crpix[1])  # Assuming array size is 2 * CRPIX
    # print(crpix[0],crpix[1])
    # print(array_size_x,array_size_y)

    # Define grid extents
    x_center, y_center = crpix  # Reference pixel is the center (1-based)
    x = np.arange(0, array_size_x, spacing_deg / pixel_spacing_x) - (array_size_x / 2) + x_center
    y = np.arange(0, array_size_y, spacing_deg /  pixel_spacing_y) - (array_size_y / 2) + y_center
    # print(len(x),len(y))
    # print(np.min(x),np.max(x))
    xx, yy = np.meshgrid(x, y)

    # Convert the pixel coordinates to RA and Dec
    ra_dec = wcs.wcs_pix2world(xx, yy, 1)  # 1 for 1-based indexing
    ra, dec = ra_dec
    
    x_flat = xx.flatten()
    y_flat = yy.flatten()
    ra_flat = ra.flatten()
    dec_flat = dec.flatten()

    fiber_id=np.arange(len(ra_flat))+1
    xtab=Table([fiber_id,x_flat,y_flat,ra_flat,dec_flat],names=['fiberid','X','Y','ra','dec'])
    xtab['ra'].format='.6f'
    xtab['dec'].format='.6f'
    
    

    return xtab



def circle_overlap_area(r, x):
    """
    Calculate the area of overlap between two circles of radius r,
    positioned x units apart. Supports scalar or array inputs for x.
    
    Parameters:
        r (float): Radius of the circles.
        x (float or np.ndarray): Distance(s) between the centers of the circles.
    
    Returns:
        float or np.ndarray: Area of the overlapping region(s).
    """
    x = np.asarray(x)  # Ensure x is a NumPy array for consistency
    
    # Initialize the overlap area array
    overlap_area = np.zeros_like(x, dtype=float)
    
    # Case 1: No overlap
    no_overlap = x >= 2 * r
    overlap_area[no_overlap] = 0
    
    # Case 2: Complete overlap
    complete_overlap = x == 0
    overlap_area[complete_overlap] = np.pi * r**2
    
    # Case 3: Partial overlap
    partial_overlap = (x > 0) & (x < 2 * r)
    x_partial = x[partial_overlap]
    
    # Calculate partial overlap areas
    alpha = 2 * np.arccos(x_partial / (2 * r))
    sector_area = 0.5 * r**2 * alpha
    triangle_area = 0.5 * r**2 * np.sin(alpha)
    overlap_area[partial_overlap] = 2 * (sector_area - triangle_area)
    
    return np.array(overlap_area)


def copy_fits_file(file_path):
    """
    Read a FITS file and create an internal copy of all its extensions.

    Parameters:
        file_path (str): Path to the FITS file to be copied.

    Returns:
        fits.HDUList: A copy of the FITS file as an HDUList object.
    """
    # Open the original FITS file
    with fits.open(file_path) as hdulist:
        # Create and return a copy of the HDUList
        hdulist_copy = fits.HDUList([hdu.copy() for hdu in hdulist])
    return hdulist_copy

def add_extension_with_same_shape(hdulist, extension_index, new_extension_name):
    """
    Add a new extension to an HDUList with the same shape as an existing extension.
    
    Parameters:
        hdulist (fits.HDUList): The HDUList to modify.
        extension_index (int): Index of the extension to copy the shape from.
        new_extension_name (str): Name of the new extension.

    Returns:
        fits.HDUList: The updated HDUList with the new extension added.
    """
    # Get the data shape from the specified extension
    original_data = hdulist[extension_index].data
    if original_data is None:
        raise ValueError(f"Extension at index {extension_index} has no data.")

    # Create a new data array with the same shape, initialized to zeros
    new_data = np.zeros_like(original_data)

    # Create a new ImageHDU with the data and header
    new_hdu = fits.ImageHDU(data=new_data, name=new_extension_name)

    # Append the new extension to the HDUList
    hdulist.append(new_hdu)

    return hdulist

def replace_header(header,new_wcs):

    for key in list(header.keys()):
        if key.startswith("WCS") or key in ["CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", 
                                            "CDELT1", "CDELT2", "CTYPE1", "CTYPE2", 
                                            "CD1_1", "CD1_2", "CD2_1", "CD2_2",
                                            "PC1_1", "PC1_2", "PC2_1", "PC2_2"]:
            header.remove(key, ignore_missing=True, remove_all=True)

    # Update the header with the new WCS
    header.update(new_wcs.to_header())
    return header



def get_wcs(header):
    try:
        wcs = WCS(header)
        if wcs.has_celestial:
            print("The header contains a valid celestial WCS.")
            print(wcs)
        else:
            print("The header has WCS, but it does not define celestial coordinates.")
    except Exception as e:
        print(f"The header does not contain a valid WCS. Error: {e}")



def prep_tables_square(wcs,filenames):
    '''
    This reads all of the files and finds the xy positions in a new WCS.  
    This version uses a wcs that has already been created to get x,y positions
    in that space

    This returns a list of tables one for each file
    '''
    ra=[]
    dec=[]
    pos_ang=[]
    slit=[]
    object=[]
    ok=True
    for one_file in filenames:
        x=fits.open(one_file)
        xhead=x['PRIMARY'].header
        # print('test',xhead['POSCIRA'])
        # print('test',xhead['POSCIDE'])
        # print('test',xhead['POSCIPA'])
        try:
            ra.append(xhead['POSCIRA'])
        except:
            print('no ra for %s' %filename)
            ra.append(0.0)
            ok=False
        try:
            dec.append(xhead['POSCIDE'])
        except:
            print('no dec for %s' % filename)
            dec.append(0.0)    
            ok=False
        try:
            pos_ang.append(xhead['POSCIPA'])
        except:
            print('Position angle missing, assuming 0.0 for %s' % one_file)
            pos_ang.append(0.0)
        try: 
            object.append(xhead['OBJECT'])
        except:
            object.append('Unknown')

        slit.append(Table(x['SLITMAP'].data))
        x.close()
    if ok==False:
        print('Something is badly wrong with headers, fix before proceeding')
        raise KeyError
    qtab=Table([filenames,object,ra,dec,pos_ang],names=['Filenames','OBJECT','RA','Dec','PA'])
    qtab.sort(['Filenames'])
    print(wcs)
    for one_tab in slit:
        x_pixels, y_pixels = wcs.world_to_pixel_values(one_tab['ra'], one_tab['dec'])
        one_tab['X']=x_pixels
        one_tab['Y']=y_pixels
        j=100
    for one_tab in slit:
        one_tab=one_tab[one_tab['telescope']=='Sci']
        one_tab=one_tab[one_tab['fibstatus']==0]
    return slit

def get_nearest(main_table,xtab,dd=35):
    '''
    Assign the closest pixel to each fiber in xtab to
    one in the final rss file.  But do not
    add if there is no overlap at all

    This is written to be analogous to frac_cale2
    '''

    dd2=dd*dd
    xtables=[]
    ff=np.array([1.])  # We are going to assign the flull vlaue of the flux to 1 fiber

    for one_row in main_table:
        xtab['d2']=(xtab['X']-one_row['X'])**2+(xtab['Y']-one_row['Y'])**2
        qtab=xtab[xtab['d2']<dd2]
        if len(qtab)>0:
            qtab.sort(['d2'])
            # So at this point we have the nearest fitst
            qtab=Table(qtab[0]) # We only want one fiber
            qtab['fib_master']=one_row['fiberid']
            qtab['frac']=ff
            xtables.append(qtab)
    ztab=vstack(xtables) 
    ztab.sort(['fib_master'])
    return ztab




def frac_calc2(main_table,xtab,dd=35.):
    '''
    xtab is the one that needs updating

    Note that it is very important when intializing new columns to make sure they
    are in the right format; otherwise one will run into the issue of format conversions
    
    '''
    dd2=dd*dd
    xtables=[]

    for one_row in main_table:
        xtab['d2']=(xtab['X']-one_row['X'])**2+(xtab['Y']-one_row['Y'])**2
        qtab=xtab[xtab['d2']<dd2]
        #print(qtab['line_no'])
        # print (ff.shape)
        if len(qtab)>0:
            qtab.sort(['d2'])
            ff=circle_overlap_area(dd,np.sqrt(qtab['d2']))
            ff=np.array(ff/np.sum(ff))
            # print(qtab['d2'],ff)
            qtab['fib_master']=one_row['fiberid']
            qtab['frac']=ff
            
            xtables.append(qtab)
    ztab=vstack(xtables) 
    ztab.sort(['fib_master'])
    return ztab

def gen_average_slit_tab(wcs,filenames):
    '''
    Create a slit_tab that contions the averagy ra's
    and decs of each fiber which can be used to 
    define a set of output fibers 
    '''

    ra=[]
    dec=[]
    base_tab=[]
    for one_file in filenames:
        x=fits.open(one_file)
        ztab=Table(x['SLITMAP'].data)
        ztab=ztab[ztab['telescope']=='Sci']
        if len(base_tab)==0:
            base_tab=ztab.copy()
            base_tab['orig_fiberid']=base_tab['fiberid']
        ra.append(np.array(ztab['ra']))
        dec.append(np.array(ztab['dec']))
    ra=np.array(ra)
    xra=np.nanmedian(ra,axis=0)
    xdec=np.nanmedian(dec,axis=0)
    base_tab['ra']=xra
    base_tab['dec']=xdec

    ra_med=np.nanmedian(ra,axis=1)
    dec_med=np.nanmedian(dec,axis=1)
    ddd=np.nanmedian(dec_med)
    ctheta=np.cos(np.deg2rad(ddd))
    # print(ctheta)

    delta_ra=3600*(np.max(ra_med)-np.min(ra_med))
    delta_dec=3600*(np.max(dec_med)-np.min(dec_med))



    print('Standard dev of offsets RA %.1f Dec %.1f' % (3600*np.std(ra_med)*ctheta,3600*np.std(dec_med)))
    print('Max-min of offsets      RA %.1f Dec %.1f' % (delta_ra*ctheta,delta_dec))

    ra_deg=np.array(base_tab['ra'])
    dec_deg=np.array(base_tab['dec'])

    x,y = wcs.world_to_pixel(SkyCoord(ra=ra_deg, dec=dec_deg, unit="deg"))

    base_tab['X']=x
    base_tab['Y']=y

    center_x, center_y = wcs.wcs.crpix


    base_tab=base_tab[base_tab['X']>0]
    base_tab=base_tab[base_tab['Y']>0]



    base_tab=base_tab[base_tab['X']< 2 * center_x]
    base_tab=base_tab[base_tab['Y']< 2 * center_y]

    # Having eliminated some fibers, we need to renumber the master fibers
    base_tab['fiberid']=np.arange(len(base_tab))+1

    
    # print(x)

    # Now we need to use the base_tab and calculate
    return base_tab['fiberid','orig_fiberid','X','Y','ra','dec']


def check_for_nan(flux,max_frac=0.5):
    '''
    Check an array for nans. Some number of which are allowed
    '''
    unique_elements, counts = np.unique(np.isnan(flux), return_counts=True)
    if True in unique_elements:
        nan_count = counts[unique_elements == True][0]
    else:
        return False

    if nan_count>max_frac*len(flux):
        return  True
    else:
        return False

def remap_one(filename,q,final_slitmap,shape):
    '''
    where filename is one file to read
    1 is the associated row of the frac table
    shape is the shape of the output 
    array
    '''

    # q.info()
    # final_slitmap.info()

    xtab=final_slitmap.copy()
    xtab['frac']=0.0
    qtab=q.copy()

    if os.path.isdir('xtmp')==False:
        os.makedirs('xtmp')
    
    VERY_BIG = 1e50
    VERY_SMALL = 1e-50

    xzero_array = np.zeros(shape, dtype=np.float32)

    print('\nStart    %s' % (filename))
    x=fits.open(filename)
    print('Flux  %10.3e %10.3e %10.e %10.3e'  % (np.nanmedian(x['FLUX'].data),np.nanstd(x['FLUX'].data),np.nanmin(x['FLUX'].data),np.nanmax(x['FLUX'].data)))
    zmax=np.nanmax(x['IVAR'].data)
    print('Ivar  %10.3e %10.3e %10.e %10.3e'  % (np.nanmedian(x['IVAR'].data),zmax*np.nanstd(x['IVAR'].data/zmax),np.nanmin(x['IVAR'].data),np.nanmax(x['IVAR'].data)))


    xflux=xzero_array.copy()
    xvar=xzero_array.copy()
    xexp=xzero_array.copy()
    # q contains the fraction allocation 
    for one_row in q:
        xmask=x['MASK'].data[one_row['fiberid']-1]
        one_flux=np.select([xmask==0],[x['FLUX'].data[one_row['fiberid']-1]],default=0)
        one_invar=np.select([xmask==0],[x['IVAR'].data[one_row['fiberid']-1]],default=VERY_SMALL)
        one_exp=np.select([xmask==0],[1.],default=0.0)
        non_finite_mask=~np.isfinite(one_flux)
        non_finite_sum=np.sum(non_finite_mask)
        if non_finite_sum>0:
            print('Found non_finite %d values' % non_finite_sum)
        one_vflux=x['IVAR'].data[one_row['fiberid']-1]
        # so one_flux, one_exp, and one_var represent one row in the flux, var, and expososure arrays
        xflux[one_row['fib_master']-1]+=one_flux*one_row['frac']
        xvar[one_row['fib_master']-1]+=one_row['frac']*one_invar  
        xexp[one_row['fib_master']-1]+=one_row['frac']*one_exp
        xtab['frac'][one_row['fib_master']-1]+=one_row['frac']
            
    # Set lines in ivar with no spectra to nan
    n=0
    while n<len(xtab):
        if xtab['frac'][n]==0:
            xflux[n]=np.nan
            xvar[n]=np.nan
        n+=1

    x['FLUX'].data=xflux
    x['IVAR'].data=xvar
    n=0
    x=add_extension_with_same_shape(x, 'FLUX','EXPOSURE')
    x['EXPOSURE'].data=xexp

    # We need to process the final_slit map table

    # print('BEFORE')
    # xtab.info()
    xtab['fibstatus']=np.select([xtab['EXPOSURE']==0,xtab['frac']==0],[99,98],0)
    # xtab.info()
    # print('AFTER')

    foo=fits.BinTableHDU(xtab)
    x['SLITMAP'].data=foo.data 

    xfilename=filename.split('/')[-1]
    x.writeto('xtmp/%s' % xfilename,overwrite=True)
    x.close()

    return 

def process_remapped_images(file_list, extension='FLUX', xproc='med',memory_limit=1_000_000_000):
    """
    Computes the median,average of an exension of  a list of FITS images while ignoring NaN and Inf values.
    and limiting the total amount of memory used

    Parameters:
    - file_list (list of str): List of FITS filenames to process.
    - extension name to process
    - type of processing med=median, sum=sum, ave or anything else average
    - memory_limit (int): Maximum memory usage in bytes (default: 1GB).
    """

    if not file_list:
        raise ValueError("No FITS files provided.")

    # Open first image to determine shape & dtype
    with fits.open(file_list[0], memmap=True) as hdul:
        shape = hdul[extension].data.shape  # (rows, cols)
        dtype = hdul[extension].data.dtype  # Data type

    # Estimate memory usage per row (in bytes)
    row_size = np.prod(shape[1:]) * np.dtype(dtype).itemsize  # Size of one row
    batch_size = max(1, memory_limit // (len(file_list) * row_size))  # Rows per batch

    print(f"Processing extension {extension}in batches of {batch_size} rows...")

    # Placeholder for median-filtered output
    xximage = np.zeros(shape, dtype=np.float32)

    # Process in batches
    for i in range(0, shape[0], batch_size):
        batch_stack = []

        # Read a batch of rows from each image
        for filename in file_list:

            with fits.open(filename, memmap=True) as hdul:
                batch_stack.append(hdul[extension].data[i:i+batch_size, :])  # Read batch

        # Convert to NumPy array and mask invalid values (NaN, Inf)
        batch_stack = np.array(batch_stack,dtype=np.float64)
        # print("  Shape:", batch_stack.shape)
        # print("  Dtype:", batch_stack.dtype)
        # print("  Contains NaN:", np.isnan(batch_stack).any())
        # print("  Contains Inf:", np.isinf(batch_stack).any())

        batch_stack[~np.isfinite(batch_stack)] = np.nan

        # Note that median for an even number of values aveages the two values closest to the medaina
        if xproc=='med':
            # print('Gotcha med')
            xximage[i:i+batch_size, :] = np.nanmedian(batch_stack, axis=0)
        elif xproc=='sum':
            # print('Gotcha sum')
            xximage[i:i+batch_size, :] = np.nansum(batch_stack, axis=0)
        else:
            # print('Gotcha ave')
            xximage[i:i+batch_size, :] = np.nanmean(batch_stack, axis=0)

        boogle=xximage[i:i+batch_size, :]
        # print('GGG %s  %e %e %e' % (extension,np.nanmedian(boogle),np.nanmin(boogle),np.nanmax(boogle)))



    # print('FINAL %s %e %e %e' %  (extension,np.nanmedian(xximage),np.nanmin(xximage),np.nanmax(xximage)))

    # ximage is a normal array
    return xximage


def do_combine(filenames,outroot='',fib_type='xy',c_type='ave'):
    ''' 
    The main routine that carries out the entire process.

    The basic steps are 
    * to create a WCS for tho output fibers, 
    * to calculate the positions of the output fibers on this WCS
    * to calculate the positions of all of the input fibers on this WCS
    * to calculate how much of each input fiber should be apportion to
    the output fibers, by one of several methods
    * to create intermediate images that contain the apportion fluxes
    from the input rss file to the output file
    * to sum or median filter the flux results
    * to calculated an new IVAR, and Mask extensions
    * to writhe everthing to and output fits file

    '''
    print('\nCombining images using %s for fluxes' %  (c_type))
    # First get the new WCS
     
    ra_center,dec_center,diam=get_size(xfiles=filenames,rad=0.25)
    print('ra,dec,size : %.2f %.2f %.2f' % (ra_center,dec_center,diam))

    # Now make a wcs
    wcs = create_wcs(ra_center, dec_center,0., diam)

    # Now create the new apetures; we need to do this here to we can set image sizes

    if fib_type=='xy':
        new_slitmap_table=generate_grid(wcs,35.)
    elif fib_type=='orig' or fib_type=='sum':
        new_slitmap_table=gen_average_slit_tab(wcs,filenames)
    else:
        print('Error: Option for creating output fibers is unknwon : ',fib_type)
        return

    # new_slitmap represents the beginning of the final SLIPMAP table

    # Now create a list of tables that contains the positions of all of the fibeers in the wc that was created

    slit=prep_tables_square(wcs, filenames)

    # Now calculate the fractions

    print('\nApportioning fractional contributions from individual to final virtual fiber positions')

    zslit=[]
    for one_tab in slit:
        if fib_type=='sum':
            one_tab=get_nearest(new_slitmap_table,one_tab)
        else:
            one_tab=frac_calc2(new_slitmap_table,one_tab)
        zslit.append(one_tab)

    # zslit is a list of tables, now with added columns indicating how to approtion the data from indvidual files into the final image

    print('The number of files being combined:',len(slit))


    foo=fits.open(filenames[0])
    final = fits.HDUList()
    final.append(foo['PRIMARY'])
    final.append(foo['FLUX'])
    final.append(foo['IVAR'])
    final.append(foo['MASK'])
    final.append(foo['WAVE'])
    final.append(foo['SLITMAP'])
    print('Adding WCS')
    # print(wcs)

    xheader = wcs.to_header()
    data = np.random.rand(6, 6)
    image_hdu = fits.ImageHDU(data=data, header=xheader, name="WCS_INFO")
    final.append(image_hdu)
    
    shape=[len(new_slitmap_table),12401]
    xzero_array = np.zeros(shape, dtype=np.float32)

    final['FLUX'].data=xzero_array.copy()
    final['IVAR'].data=xzero_array.copy()

    final=add_extension_with_same_shape(final, 'FLUX','EXPOSURE')
    
    zero_array=np.zeros(shape,dtype=np.int8)
    final['MASK'].data=zero_array
    
    

    eflux=np.zeros_like(final['FLUX'].data)
    new_slitmap_table['EXPOSURE']=0.0

    # now we want to crate various arrays that indicate how each input file would be split into the new filters system
    VERY_BIG = 1e50
    VERY_SMALL = 1e-50

    print('Begin accumulating data')

    i=0
    while i < len(filenames):
        print('starting to rmap  %s' % (filenames[i]),flush=True)
        q=zslit[i]
        for one_row in q:
            new_slitmap_table['EXPOSURE'][one_row['fib_master']-1]+=one_row['frac']
        remap_one(filenames[i],q,new_slitmap_table,shape)
        print('Finished remapping %s\n' % (filenames[i]),flush=True)
        i+=1


    # Now read them all back
    value=check_memory()
    mem=psutil.virtual_memory()
    print(f'Total Memory: {mem.total/1e9:.2f} GB')
    print(f' Used Memory: {mem.used/1e9:.2f} GB')
    print('This is %.2f percent of the available memory' % value)

    print('\nNow read the data back, and construct the final rss spectra')

    xfiles=[]
    for one_file in filenames:
        xfile=one_file.split('/')[-1]
        xfile='xtmp/%s' % xfile
        xfiles.append(xfile)

    xfl=process_remapped_images(file_list=xfiles, extension='FLUX', xproc=c_type,memory_limit=1_000_000_000)
    xexp=process_remapped_images(file_list=xfiles, extension='EXPOSURE', xproc='sum',memory_limit=1_000_000_000)
    xvar=process_remapped_images(file_list=xfiles, extension='IVAR', xproc='sum',memory_limit=1_000_000_000)
    final['FLUX'].data=xfl
    final['EXPOSURE'].data=xexp
    final['IVAR'].data=xvar


    print('Finished creating output arrays.\n')
    value=check_memory()
    mem=psutil.virtual_memory()
    print(f'Total Memory: {mem.total/1e9:.2f} GB')
    print(f' Used Memory: {mem.used/1e9:.2f} GB')
    print('This is %.2f pecent of the available memory' % value)


    new_slitmap_table['fibstatus']=0
    new_slitmap_table['targettype']='science'

    nbad=0
    for row in new_slitmap_table:
        xflux=final['FLUX'].data[row['fiberid']-1]
        if check_for_nan(xflux):
            row['fibstatus']=95
            nbad+=1
    ngood=len(new_slitmap_table)-nbad
    print('Of %d virtual fibers, %d were filled, and the %d were unfilled  ' % (len(new_slitmap_table),ngood,nbad))

    
    final['MASK'].data=np.select([final['EXPOSURE'].data==0],[1],default=0)

    hdu = fits.BinTableHDU(new_slitmap_table.as_array(), name='SLITMAP')
    final['SLITMAP']=hdu


    # Before writing everythin out, we need to add the size to the hdr

    nxmax=np.max(new_slitmap_table['X'])
    nymax=np.max(new_slitmap_table['Y'])
    final['WCS_INFO'].header['XNAXIS1']=int(nxmax)
    final['WCS_INFO'].header['XNAXIS2']=int(nymax)
    
    if outroot=='':
        outroot='test_square'

    print ('\n Final Stats for output image : %s.fits' % outroot)
    print('Combined  FLUX  %10.3e %10.3e %10.e %10.3e'  % (np.nanmedian(final['FLUX'].data),np.nanstd(final['FLUX'].data),np.nanmin(final['FLUX'].data),np.nanmax(final['FLUX'].data)))
    zmax=np.nanmax(final['IVAR'].data)
    print('Combined  IVAR  %10.3e %10.3e %10.e %10.3e\n'  % (np.nanmedian(final['IVAR'].data),zmax*np.nanstd(final['IVAR'].data/zmax),np.nanmin(final['IVAR'].data),np.nanmax(final['IVAR'].data)))


    new_slitmap_table.write('%s.tab' % outroot,format='ascii.fixed_width_two_line',overwrite=True)
    final.writeto(outroot+'.fits',overwrite=True)
    print('Wrote new RSS file: %s.fits' % (outroot))
    return
    

def steer(argv):
    xfiles=[]
    outroot='xtest'
    fib_type='xy'
    c_type='ave'
    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-outroot':
            i+=1
            outroot=argv[i]
        elif argv[i][:5]=='-orig' and fib_type!='sum':
            fib_type='orig'
        elif argv[i]=='-sum':
            fib_type=='sum'
        elif argv[i][0:4]=='-med':
            c_type='med'
        elif argv[i][0]=='-':
            print('Could not parse command line:',argv)
            return
        else:
            xfiles.append(argv[i])
        i+=1

    files=[]
    for one in xfiles:
        if one.count('*'):
            qfiles=glob(one)
            for one_file in qfiles:
                files.appen
        else:
            files.append(one)

    if len(files)==0:
        print('There are no files to combine: ',argv)
        return
    else:
        print('Creating new rss file from  %d files' % len(files))
        for one_file in files:
            print(one_file)
        print('\n')

    # Do some checks on the files
    xtab=xcheck(files)

    ftab=xtab[xtab['FluxCal']=='False']
    if len(ftab)>0:
        print('!!! There %d  are files that have not been flux calibrated' % (len(ftab)))
        print(ftab)
        print('!!! Fix this before proeeeding')
        return

    xver,counts=np.unique(xtab['DRP'],return_counts=True)
    if len(xver)>1:
        print('Warning - Multiple DRP versions are being processed')
        k=0
        while k<len(xver):
            print('Ver %20s  Number %3d' % (xver[k],counts[k]))
            k+=1
        print('Continuing')
    else:
        print('All files seem to be from the save DRP version and all are fluxe calibrated')

    do_combine(files,outroot,fib_type,c_type)   
    





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
