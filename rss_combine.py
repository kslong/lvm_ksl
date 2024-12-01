#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Combine multiple LVM row-stacked spectra into a 
sincle row-staked spectrum, by creating a grid 
of fibegs that cover the entire area, and then
assigning portions of the spectra of various fibers 
to the final one.


Command line usage (if any):

    usage: rss_combine.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
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
    # print(xtab)
    ra_cen=np.average(xtab['RA'])
    dec_cen=np.average(xtab['Dec'])
    delta_dec=np.max(np.fabs(xtab['Dec']-dec_cen))
    delta_ra=np.max(np.fabs(xtab['RA']-ra_cen))*np.cos(dec_cen/57.29578)
    # print(delta_dec,delta_ra)
    size=2*rad+np.max([delta_ra,delta_dec])
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
        [pixel_scale * cos_theta, -pixel_scale * sin_theta],
        [pixel_scale * sin_theta,  pixel_scale * cos_theta]
        ])
    
    # Define the projection type
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]  # Tangential (gnomonic) projection
    
    return wcs


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

    print(pixel_spacing_x,pixel_spacing_y) # This the pixel scale in arcsec
    
    # cdelt = wcs.wcs.cdelt  # Pixel scale in degrees
    spacing_deg = spacing_arcsec / 3600.0  # Convert spacing to degrees
    print('deg',spacing_deg)

    # Calculate pixel spacing in terms of pixels
    # pixel_spacing_x = spacing_deg / abs(cdelt[0])
    # pixel_spacing_y = spacing_deg / abs(cdelt[1])

    # Define the array size based on CRPIX
    array_size_x = int(2 * crpix[0])  # Assuming array size is 2 * CRPIX
    array_size_y = int(2 * crpix[1])  # Assuming array size is 2 * CRPIX
    print(crpix[0],crpix[1])
    print(array_size_x,array_size_y)

    # Define grid extents
    x_center, y_center = crpix  # Reference pixel is the center (1-based)
    x = np.arange(0, array_size_x, spacing_deg / pixel_spacing_x) - (array_size_x / 2) + x_center
    y = np.arange(0, array_size_y, spacing_deg /  pixel_spacing_y) - (array_size_y / 2) + y_center
    print(len(x),len(y))
    print(np.min(x),np.max(x))
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
    This version uses a wcs that has already been creaed to get x,y positions
    in that apace

    This returns a list of tables
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
        print('XXX',one_tab['ra'][j],one_tab['dec'][j],one_tab['X'][j],one_tab['Y'][j])
        # print('stat ',np.nanmin(one_tab['X']),np.nanmax(one_tab['X']),np.nanmin(one_tab['Y']),np.nanmax(one_tab['Y']))
    for one_tab in slit:
        one_tab=one_tab[one_tab['telescope']=='Sci']
        one_tab=one_tab[one_tab['fibstatus']==0]
    return slit


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



def do_square(filenames,outroot=''):
    ''' 
    This is supposed to do the whole process
    '''
    # First get the new WCS
     
    ra_center,dec_center,diam=get_size(xfiles=filenames,rad=0.25)

    # Now make a wcs
    wcs = create_wcs(ra_center, dec_center,0., diam)

    # Now create the new apetures; we need to do this here to we can set image sizes

    new_slitmap_table=generate_grid(wcs,35.)

    new_slitmap_table.info()

    # Now create a table that contains the positions of all of the fibeers in the wc that was created

    slit=prep_tables_square(wcs, filenames)

    # Now calculate the fractions

    slit[0].info()


    zslit=[]
    for one_tab in slit:
        one_tab=frac_calc2(new_slitmap_table,one_tab)
        zslit.append(one_tab)
    print('test',len(slit))

    zslit[0].info()

    # At this point contributions calculated for each final 'fiber'

    print('The number of files being combined:',len(slit))

    # Now create a new fits file, we can modify

    foo=fits.open(filenames[0])
    final = fits.HDUList()
    final.append(foo['PRIMARY'])
    final.append(foo['FLUX'])
    final.append(foo['IVAR'])
    final.append(foo['MASK'])
    final.append(foo['WAVE'])
    final.append(foo['SLITMAP'])
    print('Adding WCS')
    print(wcs)
    # replace_header(foo['PRIMARY'].header,wcs)
    # final['PRIMARY'].header['MOSAIC']='TRUE'
    # 

    xheader = wcs.to_header()
    data = np.random.rand(6, 6)
    image_hdu = fits.ImageHDU(data=data, header=xheader, name="WCS_INFO")
    final.append(image_hdu)
    
                 
                 
    

    # final=copy_fits_file(filenames[0])
    shape=[len(new_slitmap_table),12401]
    zero_array = np.zeros(shape, dtype=np.float32)

    final['FLUX'].data=zero_array
    final['IVAR'].data=zero_array

    final=add_extension_with_same_shape(final, 'FLUX','EXPOSURE')
    
    zero_array=np.zeros(shape,dtype=np.int8)
    final['MASK'].data=zero_array
    
    
    final.info()
 
    # Now fill the arrays

    eflux=np.zeros_like(final['FLUX'].data)
    new_slitmap_table['EXPOSURE']=0.0



    # final.info()
    # zslit[0].info()
    i=0
    while i<len(filenames):
        x=fits.open(filenames[i])
        q=zslit[i]       
        for one_row in q:
            new_slitmap_table['EXPOSURE'][one_row['fib_master']-1]+=one_row['frac']
            flux=x['FLUX'].data[one_row['fiberid']-1]
            # sky=x['SKY'].data[one_row['fiberid']-1]
            vflux=x['IVAR'].data[one_row['fiberid']-1]
            # vsflux=x['SKY_IVAR'].data[one_row['fiberid']-1]
            if i==0:
                final['FLUX'].data[one_row['fib_master']-1]+=one_row['frac']*flux
                # final['SKY'].data[one_row['fib_master']-1]+=one_row['frac']*sky
                final['EXPOSURE'].data[one_row['fib_master']-1]+=one_row['frac']
                eflux[one_row['fib_master']-1]+=one_row['frac']/vflux
                # esflux[one_row['fib_master']-1]+=one_row['frac']/vsflux
                if one_row['fib_master']==103:
                    print('got for %d %f %e %e -> %e' % (i,one_row['frac'],np.nanmedian(final['FLUX'].data[one_row['fib_master']-1]),
                                                np.nanmedian(vflux), np.nanmedian(eflux[one_row['fib_master']-1])
                                               ))
            else:
                # print('ok',one_row['fiberid']-1,one_row['fibstatus'])
                if np.all(np.isnan(flux))== False:
                    valid_mask=~np.isnan(flux)
                    final['FLUX'].data[one_row['fib_master']-1][valid_mask]+=one_row['frac']*flux[valid_mask]
                    # final['SKY'].data[one_row['fib_master']-1][valid_mask]+=one_row['frac']*sky[valid_mask]
                    final['EXPOSURE'].data[one_row['fib_master']-1][valid_mask]+=one_row['frac']
                    eflux[one_row['fib_master']-1][valid_mask]+=one_row['frac']*1./vflux[valid_mask]
                    # esflux[one_row['fib_master']-1][valid_mask]+=one_row['frac']/vsflux[valid_mask]
                    if one_row['fib_master']==103:
                        print('got for %d %f %e %e -> %e' % (i,one_row['frac'],np.nanmedian(final['FLUX'].data[one_row['fib_master']-1]),
                                                np.nanmedian(vflux), np.nanmedian(eflux[one_row['fib_master']-1])
                                               ))
        i+=1
   
    plt.hist(new_slitmap_table['EXPOSURE'].data,1000,range=(0,len(filenames)),cumulative=-1,density=True)

    
    # plt.ylim(0,10)  
    # print('At end  %e' % (np.nanmedian(final['FLUX'].data[102])))
    # print('At exp  %e' % (np.nanmean(final['EXPOSURE'].data[102])))
    final['EXPOSURE'].data[final['EXPOSURE'].data==0]= 1.
    final['FLUX'].data/=final['EXPOSURE'].data
    # final['SKY'].data/=final['EXPOSURE'].data
    final['IVAR'].data=final['EXPOSURE'].data*final['EXPOSURE'].data/(eflux)
    # final['IVAR'].data=1./(eflux)
    # final['SKY_IVAR'].data=final['EXPOSURE'].data*final['EXPOSURE'].data/(esflux)

    hdu = fits.BinTableHDU(new_slitmap_table.as_array(), name='SLITMAP')
    final['SLITMAP']=hdu

    # Before writing everythin out, we need to add the size to the hdr

    nxmax=np.max(new_slitmap_table['X'])
    nymax=np.max(new_slitmap_table['Y'])
    final['WCS_INFO'].header['XNAXIS1']=int(nxmax)
    final['WCS_INFO'].header['XNAXIS2']=int(nymax)
    

    # final['SLITMAP'].data=new_slitmap_table

    print('At end  %e %e -> %e' % (np.nanmedian(final['FLUX'].data[102]),np.nanmedian(eflux[102]),np.nanmedian(final['IVAR'].data[102])))

    if outroot=='':
        outroot='test_square'

    final.writeto(outroot+'.fits',overwrite=True)
    

def steer(argv):
    xfiles=[]
    outroot='xtest'
    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-outroot':
            i+=1
            outroot=argv[i]
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
    else:
        print('Combining %d files' % len(files))
        for one_file in files:
            print(one_file)

    do_square(files,outroot)   
    





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
