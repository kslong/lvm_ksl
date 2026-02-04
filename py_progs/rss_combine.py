#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Combine multiple LVM row-stacked (SFrame) spectra into a
single row-stacked spectrum, by creating a grid
of fibers that cover the entire area, and then
assigning portions of the spectra of various fibers
to the final one.  Alternatively, combine the
exposures for a single tile into one RSS image with
"fibers" centered on the mean position of the exposures.


Command line usage::

    rss_combine.py [-orig] [-sum] [-med] [-keep] [-outroot xxxx] filenames

The ``-outroot xxxx`` option sets the name for the output file.

The ``-orig`` option uses the average of the input fiber positions as
output fiber positions instead of creating a regular grid. The fractional
contributions are still calculated.

The ``-sum`` option assigns each input fiber's entire spectrum to the
nearest output fiber (no fractional splitting). Output fiber locations
are taken from the average input fiber locations as with -orig.

The ``-med`` option uses the median instead of the mean when combining
remapped images from multiple input files. This is more robust to outliers
and bad pixels. Note that -med only affects the final combination, not
how flux is apportioned (that is controlled by -orig and -sum).

The ``-keep`` option keeps the temporary xtmp/ directory containing
intermediate remapped files. By default, this directory is deleted
after processing.


Description:

The routine combines multiple RSS images to create a final output image,
with various options.

With no switches the routine reads the input SFrames, and constructs
output fiber positions on a regular grid on the sky. Given the positions
of these fibers it apportions flux from the original images to the output
fibers, based on the area of overlap.

With the -orig switch, the output grid is constructed from the input grid,
based on the average position (RA and Dec) of each of the input fibers.
The original spectra are still apportioned based on the overlaps. This may
be a better approach for a single tile than the regular grid approach.

With the -sum switch the positions of the RSS fibers in the output image
are taken as above, but additionally instead of splitting the spectra
between output fibers, all of the flux from an input fiber is assigned to
one output fiber. This may be the best approach if the input SFrame files
were not dithered.

A single FITS file is produced with the following extensions: PRIMARY (the
primary header from one of the input images), FLUX (the flux array, with
rows depending on how many artificial fibers were created), IVAR (the final
inverse variance), MASK (set to 0 if this fiber and wavelength has any valid
exposure, 1 otherwise), WAVE (the standard set of wavelengths), SLITMAP (a
version of the slitmap file with RAs, Decs, and X/Y positions of fibers),
WCS_INFO (the WCS created to define fiber positions), and EXPOSURE (the
effective number of exposures per fiber).

Primary routines:

The main routine is do_combine.

Notes:

The final image is not weighted by the errors. All of the input images are
treated identically.

There is no attempt to scale the images to one another globally, though one
might argue that would be a good idea.

It is apparent that when one uses the mean, and one exposure is somehow
different, the effects of bad pixels are more apparent than when one uses
the median.



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
import shutil
import psutil

def xcheck(xfiles):
    '''
    Check the input files for obvious problems.

    Validates that each file can be opened, has flux calibration applied,
    and extracts metadata including DRP version, MJD, coordinates, tile ID,
    exposure time, and heliocentric velocity.

    Parameters:
        xfiles (list): List of FITS filenames to check.

    Returns:
        astropy.table.Table: Table with metadata for each file including a Good column indicating whether the file passed quality checks.
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
    good=[]
    for one in xfiles:
        try:
            f=fits.open(one)
            xhead=f['PRIMARY'].header
            # print('Got ;',one)
        except:
            print('Could not open ',one)
            pass

        ok='Yes'
        try:
            drp.append(xhead['drpver'])
        except:
            drp.append('Unknown')
        try:
            commit.append(xhead['COMMIT'])
        except:
            commit.append('Unknown')

        try:
            xfluxcal=xhead['FLUXCAL']
            if xfluxcal=='FALSE' or xfluxcal == 'NONE':
                ok='No'
        except:
            xfluxcal='Unknown'
            ok='No'
        fluxcal.append(xfluxcal)


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
        good.append(ok)



    xtab=Table([xfiles,mjd,drp,commit,fluxcal,unit,helio,ra,dec,tile,exptime,obj,good],names=['Filename','MJD','DRP','Commit','FluxCal','BUNIT','RA','HeleoV','Dec','Tile_ID','EXPTIME','Source_name','Good'])
    xtab.sort((['Filename']))




    ftab=xtab[xtab['Good']=='No']
    if len(ftab)>0:
        print('!!! There %d of %d files that should not be used'% (len(ftab),len(xtab)))
        print(ftab)

    xver,counts=np.unique(xtab['DRP'],return_counts=True)
    if len(xver)>1:
        print('Warning - Multiple DRP versions are being processed')
        k=0
        while k<len(xver):
            print('Ver %20s  Number %3d' % (xver[k],counts[k]))
            k+=1
        print('Continuing')
    else:
        print('All files seem to be from the same DRP version and all are flux calibrated')


    # print(xtab)
    return xtab

def check_array_info(arr):
    '''
    Print diagnostic information about a numpy array.

    Parameters:
        arr (np.ndarray or np.ma.MaskedArray): The array to inspect.

    Returns:
        None: Prints shape, dtype, and array type to stdout.
    '''
    array_type = "MaskedArray" if isinstance(arr, np.ma.MaskedArray) else "ndarray"
    print( {"shape": arr.shape, "dtype": arr.dtype, "type": array_type})

def check_memory(percentage=80.):
    '''
    Calculate the percentage of memory being used and report a warning
    if it is greater than a given percentage.

    Parameters:
        percentage (float): Threshold percentage above which to print warnings. Default is 80%.

    Returns:
        float: The current memory usage as a percentage.
    '''
    mem=psutil.virtual_memory()
    frac=mem.used/mem.total*100
    if frac>percentage:
        print(f'Warning: Total Memory: {mem.total/1e9:.2f} GB')
        print(f'Warning:  Used Memory: {mem.used/1e9:.2f} GB')
        print(f'Warning: {frac:.2f}% of memory is being used')
    return frac


def get_size(xfiles,rad=0.25):
    '''
    From SFrame images calculate the center and size of an image
    that would encompass them all.

    Reads the telescope pointing coordinates (TESCIRA, TESCIDE) from
    each file's primary header to determine the bounding box.

    Parameters:
        xfiles (list): List of SFrame FITS filenames.
        rad (float): Additional radius in degrees to add as padding around the computed extent. Default is 0.25 degrees.

    Returns:
        tuple: (ra_center, dec_center, size) where ra_center and dec_center are the center coordinates in degrees, and size is the diameter of the bounding region in degrees.
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
    



def create_wcs(ra_deg, dec_deg, pos_ang, size_deg):
    """
    Create a WCS object centered on a given RA and Dec.

    Creates a tangent plane projection (TAN) with 1 arcsec/pixel scale,
    optionally rotated by the specified position angle.

    Parameters:
        ra_deg (float): Right Ascension in degrees (center of the image).
        dec_deg (float): Declination in degrees (center of the image).
        pos_ang (float): Position angle in degrees for rotation of the WCS.
        size_deg (float): Size of the image in degrees. The image will have 3600 times size_deg pixels on each side.

    Returns:
        astropy.wcs.WCS: Configured WCS object with TAN projection.
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

def create_wcs_from_points(xpos, ypos, ra_deg, dec_deg):
    '''
    Create a WCS from a series of pixel positions and their corresponding
    sky coordinates.

    Uses astropy's fit_wcs_from_points to derive a WCS transformation
    that maps the given pixel coordinates to the given RA/Dec values.

    Parameters:
        xpos (array-like): X pixel coordinates.
        ypos (array-like): Y pixel coordinates.
        ra_deg (array-like): Right Ascension values in degrees.
        dec_deg (array-like): Declination values in degrees.

    Returns:
        astropy.wcs.WCS or None: The fitted WCS object, or None if fitting fails.
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
        return None



def generate_grid(wcs, spacing_arcsec):
    """
    Generate a regular grid of pixel positions centered on the WCS reference pixel,
    and compute the corresponding RA and Dec for each position.

    Creates virtual fiber positions on a regular grid that can be used as
    output fiber locations when combining multiple RSS exposures.

    Parameters:
        wcs (astropy.wcs.WCS): The WCS object defining the coordinate system.
        spacing_arcsec (float): Spacing between grid points in arcseconds.

    Returns:
        astropy.table.Table: Table with columns fiberid, X, Y, ra, dec containing the grid positions in both pixel and sky coordinates.
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

def replace_header(header, new_wcs):
    '''
    Replace the WCS information in a FITS header with a new WCS.

    Removes all existing WCS-related keywords from the header and
    replaces them with the keywords from the new WCS object.

    Parameters:
        header (astropy.io.fits.Header): The FITS header to modify.
        new_wcs (astropy.wcs.WCS): The new WCS to insert into the header.

    Returns:
        astropy.io.fits.Header: The modified header with updated WCS.
    '''
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
    '''
    Extract and validate a WCS from a FITS header.

    Attempts to create a WCS object from the header and checks whether
    it contains valid celestial coordinates. Prints diagnostic messages
    about the WCS status.

    Parameters:
        header (astropy.io.fits.Header): The FITS header to parse.

    Returns:
        None: Prints WCS information or error messages to stdout.
    '''
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

def get_nearest(main_table, xtab, dd=35):
    '''
    Assign each input fiber to the nearest output fiber.

    For each output fiber in main_table, finds the nearest input fiber
    in xtab and assigns the full flux value (fraction=1.0) to it. This
    is used with the -sum option where no fractional splitting is done.

    This is written to be analogous to frac_calc2.

    Parameters:
        main_table (astropy.table.Table): Table of output fiber positions with fiberid, X, Y columns.
        xtab (astropy.table.Table): Table of input fiber positions with X, Y columns.
        dd (float): Maximum distance in pixels for a fiber to be considered overlapping. Default is 35 (matching fiber diameter).

    Returns:
        astropy.table.Table or list: Table with input fibers assigned to output fibers with fib_master and frac columns added, or empty list if no matches found.
    '''

    if len(main_table)==0:
        print('Error: get_nearest called with no data')
        return []

    dd2=dd*dd
    xtables=[]
    ff=np.array([1.])  # We are going to assign the full value of the flux to 1 fiber



    for one_row in main_table:
        xtab['d2']=(xtab['X']-one_row['X'])**2+(xtab['Y']-one_row['Y'])**2
        qtab=xtab[xtab['d2']<dd2]
        if len(qtab)>0:
            qtab.sort(['d2'])
            # So at this point we have the nearest first
            qtab=Table(qtab[0]) # We only want one fiber
            qtab['fib_master']=one_row['fiberid']
            qtab['frac']=ff
            xtables.append(qtab)

    try:
        ztab=vstack(xtables) 
        ztab.sort(['fib_master'])
    except ValueError:
        print('Error: get_nearest: nothing to stack')
        return []
    return ztab




def frac_calc2(main_table, xtab, dd=35.):
    '''
    Calculate fractional flux contributions from input fibers to output fibers.

    For each output fiber in main_table, finds all nearby input fibers in xtab
    and calculates what fraction of each input fiber's flux should be assigned
    to the output fiber, based on the area of overlap between circular apertures.

    Note that it is very important when initializing new columns to make sure they
    are in the right format; otherwise one will run into the issue of format conversions.

    Parameters:
        main_table (astropy.table.Table): Table of output fiber positions with fiberid, X, Y columns.
        xtab (astropy.table.Table): Table of input fiber positions with X, Y columns.
        dd (float): Fiber diameter in pixels used to calculate overlap area. Default is 35.

    Returns:
        astropy.table.Table or list: Table with input fibers and their fractional contributions to output fibers (fib_master and frac columns added), or empty list if no matches found.
    '''
    dd2=dd*dd
    xtables=[]

    if len(main_table)==0:
        print('Error: frac_calc2 received main_table with no rows')

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
    try:
        ztab=vstack(xtables) 
        ztab.sort(['fib_master'])
    except ValueError:
        print('Error: frac_calc2: stacking failed with nothing to stack')
        return[]
    return ztab

def gen_average_slit_tab(wcs, filenames):
    '''
    Create a slitmap table containing the average positions of fibers.

    Reads the SLITMAP from each input file and computes the median RA and Dec
    for each fiber across all exposures. This creates output fiber positions
    that represent the average pointing, which is useful for combining
    exposures of the same tile.

    Parameters:
        wcs (astropy.wcs.WCS): The WCS object for converting sky to pixel coordinates.
        filenames (list): List of SFrame FITS filenames to process.

    Returns:
        astropy.table.Table: Table with columns fiberid, orig_fiberid, X, Y, ra, dec for the averaged fiber positions.
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


def check_for_nan(flux, max_frac=0.5):
    '''
    Check an array for NaN values.

    Determines whether the fraction of NaN values in the array exceeds
    a specified threshold.

    Parameters:
        flux (np.ndarray): The array to check.
        max_frac (float): Maximum allowed fraction of NaN values. Default is 0.5 (50%).

    Returns:
        bool: True if NaN fraction exceeds max_frac, False otherwise.
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

def remap_one(filename, q, final_slitmap, shape):
    '''
    Remap flux from one input file to the output fiber grid.

    Reads an input SFrame file and redistributes its flux to the output
    fiber positions according to the fractional contributions in table q.
    Writes an intermediate file to the 'xtmp/' directory.

    Parameters:
        filename (str): Path to the input SFrame FITS file.
        q (astropy.table.Table): Table containing the fractional mapping with fiberid, fib_master, and frac columns.
        final_slitmap (astropy.table.Table): The output fiber position table.
        shape (tuple): Shape of the output arrays (n_fibers, n_wavelengths).

    Returns:
        None: Writes remapped data to xtmp directory.
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
        # so one_flux, one_exp, and one_var represent one row in the flux, var, and exposure arrays
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

def process_remapped_images(file_list, extension='FLUX', xproc='med', memory_limit=1_000_000_000):
    """
    Compute the median, sum, or average of a FITS extension across multiple images.

    Processes the images in batches to limit memory usage, ignoring NaN and Inf
    values in the computation.

    Parameters
    ----------
    file_list : list of str
        List of FITS filenames to process.
    extension : str
        Name of the FITS extension to process. Default is 'FLUX'.
    xproc : str
        Type of processing to perform. Options are 'med' for median,
        'sum' for sum, or 'ave' (or anything else) for mean.
    memory_limit : int
        Maximum memory usage in bytes. Default is 1GB.

    Returns
    -------
    np.ndarray
        The combined image array with the specified statistic computed
        along the file axis.
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

    print(f"Processing extension {extension} in batches of {batch_size} rows...")

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

        # Note that median for an even number of values averages the two values closest to the median
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


def do_combine(filenames, outroot='', fib_type='xy', c_type='ave', keep_tmp=False):
    '''
    Main routine that carries out the entire RSS combining process.

    The basic steps are: (1) create a WCS for the output fibers,
    (2) calculate the positions of output and input fibers on this WCS,
    (3) calculate how much of each input fiber should be apportioned to
    the output fibers, (4) create intermediate images with apportioned fluxes,
    (5) average or median filter the flux results, (6) calculate new IVAR
    and MASK extensions, and (7) write everything to an output FITS file.

    Parameters
    ----------
    filenames : list
        List of input SFrame FITS filenames to combine.
    outroot : str
        Root name for the output file. Default creates 'test_square.fits'.
    fib_type : str
        Method for creating output fiber positions and apportioning flux.
        This controls whether frac_calc2 or get_nearest is used.
        Options are 'xy' (regular grid, uses frac_calc2), 'orig' (average
        input positions, uses frac_calc2), or 'sum' (average input positions,
        uses get_nearest with no fractional splitting).
    c_type : str
        Method for combining the remapped images from multiple input files.
        This is applied after the flux apportionment step. Options are
        'ave' for mean (default) or 'med' for median.
    keep_tmp : bool
        If True, keep the temporary xtmp/ directory containing intermediate
        remapped files. Default is False (deletes temp files).

    Returns
    -------
    None
        Writes output to '<outroot>.fits' and '<outroot>.tab'.
    '''
    # Do Quality checks on the files
    xtab=xcheck(filenames)
    print(xtab)

    xgood=xtab[xtab['Good']=='Yes']
    if len(xgood)<len(xtab):
        print('do_combine: %d files of %d were eliminated due to quality checks' % (len(xtab)-len(xgood),len(xtab)))
        xbad=xtab[xtab['Good']=='No']
        for one in xbad:
            print('Eliminated %s' % one['Filename'])
        filenames=xgood['Filename']


    print('\nCombining images using %s for fluxes' %  (c_type))
    # First get the new WCS
     
    ra_center,dec_center,diam=get_size(xfiles=filenames,rad=0.25)
    print('ra,dec,size : %.2f %.2f %.2f' % (ra_center,dec_center,diam))

    # Now make a wcs
    wcs = create_wcs(ra_center, dec_center,0., diam)

    # Now create the new apertures; we need to do this here so we can set image sizes

    if fib_type=='xy':
        new_slitmap_table=generate_grid(wcs,35.)
    elif fib_type=='orig' or fib_type=='sum':
        new_slitmap_table=gen_average_slit_tab(wcs,filenames)
    else:
        print('Error: Option for creating output fibers is unknown: ',fib_type)
        return

    # new_slitmap represents the beginning of the final SLITMAP table

    # Now create a list of tables that contains the positions of all of the fibers in the WCS that was created

    slit=prep_tables_square(wcs, filenames)

    # Now calculate the fractions

    print('\nApportioning fractional contributions from individual to final virtual fiber positions')

    zslit=[]
    for one_tab in slit:
        if fib_type=='sum':
            one_tab=get_nearest(new_slitmap_table,one_tab)
        else:
            one_tab=frac_calc2(new_slitmap_table,one_tab)
        if len(one_tab)>0:
            zslit.append(one_tab)

    # zslit is a list of tables, now with added columns indicating how to apportion the data from individual files into the final image

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
    
    zero_array=np.zeros(shape,dtype=np.int64)
    final['MASK'].data=zero_array
    
    

    eflux=np.zeros_like(final['FLUX'].data)
    new_slitmap_table['EXPOSURE']=0.0

    # now we want to create various arrays that indicate how each input file would be split into the new fibers system
    VERY_BIG = 1e50
    VERY_SMALL = 1e-50

    print('Begin accumulating data')

    i=0
    while i < len(filenames):
        print('Starting to remap %s' % (filenames[i]),flush=True)
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
    print('This is %.2f percent of the available memory' % value)


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
    # Explicitly set BITPIX to a positive integer value
    final['MASK'].header['BITPIX'] = 64 
    # print(final['EXPOSURE'].data[5500])
    # print(final['MASK'].data[5500])

    hdu = fits.BinTableHDU(new_slitmap_table.as_array(), name='SLITMAP')
    final['SLITMAP']=hdu


    # Before writing everything out, we need to add the size to the header

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

    # Clean up temporary directory unless keep_tmp is True
    if keep_tmp==False and os.path.isdir('xtmp'):
        shutil.rmtree('./xtmp')
        print('Removed temporary directory xtmp/')
    elif keep_tmp:
        print('Keeping temporary directory xtmp/')

    return
    

def steer(argv):
    '''
    Parse command line arguments and run the RSS combining process.

    Parameters:
        argv (list): Command line arguments (sys.argv).

    Returns:
        None: Calls do_combine with parsed arguments.
    '''
    xfiles=[]
    outroot='xtest'
    fib_type='xy'
    c_type='ave'
    keep_tmp=False
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
            fib_type='sum'
        elif argv[i][0:4]=='-med':
            c_type='med'
        elif argv[i]=='-keep':
            keep_tmp=True
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
                files.append(one_file)
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

    do_combine(files,outroot,fib_type,c_type,keep_tmp)   
    





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
