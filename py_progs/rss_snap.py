#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Create one or more subimages of specific regions observed with 
lvm using a file derived from DRP_All and a set of postions 
containing the exposures asscoated with each object.


Command line usage::

    rss_snap.py [-h] [-keep] [-redo] [-all] xfile source_name

Arguments: xfile is a version of an expanded master file containing
the source names and associated exposures (one row per exposure to be
combined). source_name is one or more of the source_names in xfile.

Options:

-h
    Prints this documentation.

-keep
    Retains the temporary files from individual exposures in the
    directory xtmp. Without this switch the temporary files are deleted.

-redo
    Recreates the combined fits files, even if they exist in the Snap
    directory.

-all
    Causes all of the sources in xfile to be done.

Description:

Primary routines:

    do_one

Notes:
                                       
History:

250704 ksl Coding begun

'''

import os
import sys
from lvm_ksl import rss_combine_pos
from lvm_ksl import lvm_gaussfit
from lvm_ksl.get_vel import get_vel


from astropy.io import ascii, fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np


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



def get_files(xsum='lmc.out',source_name='whatever'):
    xtop=find_top()
    xtab=ascii.read(xsum)
    xtab=xtab[xtab['Source_name']==source_name]

    filenames=[]
    if len(xtab)==0:
        print('No obs found for %s' % source_name)
        return -99.,-99.,[]
    else:
        ra=xtab['RA'][0]
        dec=xtab['Dec'][0]
    for one in xtab:
        xfile='%s/%s' % (xtop,one['location'])
        if os.path.isfile(xfile):
            filenames.append(xfile)
        else:
            print('Could not find ',xfile)
        
    return ra, dec, filenames






def plot_one(xtab, var='flux_ha', ymin=0, ymax=0, label='', ax=None, marker_size=240):
    '''
    Make one panel of a larger plot
    '''

    ra = xtab['ra']
    dec = xtab['dec']
    if isinstance(var, str):
        value = xtab[var]
    else:
        value=var
        
    qval = value.copy()
    
    if ymin < ymax:
        value = np.select([value < ymin, value > ymax], [ymin, ymax], default=value)
    
    # Calculate RA spacing with cos(dec) correction
    mean_dec_rad = np.radians(np.mean(dec))
    cos_dec_factor = np.cos(mean_dec_rad)
    
    # Use provided axis or get current axis
    if ax is None:
        ax = plt.gca()
    
    # Apply the correction factor to RA values for plotting
    ra_scaled = ra * cos_dec_factor
    
    sc = ax.scatter(ra_scaled, dec, s=marker_size, c=value, cmap='viridis', alpha=0.8)
    
    # Invert RA axis so that RA increases to the left
    ax.invert_xaxis()
    
    # Label the axes
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    
    # Set aspect ratio to be equal
    ax.set_aspect('equal')
    
    # Calculate auto-range if needed
    if ymax == 0.0:
        ymin = np.percentile(value, 20)
        ymax = np.percentile(value, 80)
        sc.set_clim(ymin, ymax)
    
    # Add colorbar for flux values
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = plt.colorbar(sc, cax=cax)
    
    if label == '':
        label = 'Flux'
    cbar.set_label(label)
    
    # Properly set tick positions and labels to account for cos(dec) scaling
    # First, get the current x-ticks
    original_tick_locs = ax.get_xticks()
    
    # Prepare both positions and labels
    tick_positions = original_tick_locs.copy()
    tick_labels = [f"{tick/cos_dec_factor:.1f}" for tick in original_tick_locs]
    
    # Set both the tick positions and labels
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    
    return sc


# Add new routine
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_one_interpolated(xtab, var='flux_ha', ymin=0, ymax=0, label='', ax=None, 
                         grid_resolution=100, interpolation_method='linear', 
                         fill_value=np.nan, smooth_factor=1.0):
    '''
    Make one panel of a larger plot using interpolation to create an image
    
    Parameters:
    -----------
    xtab : astropy.table.Table
        Table containing 'ra', 'dec', and variable columns
    var : str or array-like
        Variable name in table or array of values to plot
    ymin, ymax : float
        Color scale limits (auto-calculated if ymax=0)
    label : str
        Label for colorbar
    ax : matplotlib.axes.Axes
        Axis to plot on (uses current axis if None)
    grid_resolution : int
        Number of grid points along each axis for interpolation
    interpolation_method : str
        Interpolation method: 'linear', 'nearest', or 'cubic'
    fill_value : float
        Value to use for points outside convex hull (np.nan for transparent)
    smooth_factor : float
        Factor to smooth the grid resolution (1.0 = no smoothing)
    '''
    
    ra = xtab['ra']
    dec = xtab['dec']
    
    if isinstance(var, str):
        value = xtab[var]
    else:
        value = var
    
    # Remove any NaN values
    mask = ~(np.isnan(ra) | np.isnan(dec) | np.isnan(value))
    ra = ra[mask]
    dec = dec[mask]
    value = value[mask]
    
    # Apply value clipping if specified
    qval = value.copy()
    if ymin < ymax:
        value = np.select([value < ymin, value > ymax], [ymin, ymax], default=value)
    
    # Calculate RA spacing with cos(dec) correction
    mean_dec_rad = np.radians(np.mean(dec))
    cos_dec_factor = np.cos(mean_dec_rad)
    
    # Use provided axis or get current axis
    if ax is None:
        ax = plt.gca()
    
    # Apply the correction factor to RA values for interpolation
    ra_scaled = ra * cos_dec_factor
    
    # Create regular grid for interpolation
    ra_min, ra_max = ra_scaled.min(), ra_scaled.max()
    dec_min, dec_max = dec.min(), dec.max()
    
    # Add small padding to grid
    ra_padding = (ra_max - ra_min) * 0.05
    dec_padding = (dec_max - dec_min) * 0.05
    
    # Create grid with specified resolution
    grid_ra = np.linspace(ra_min - ra_padding, ra_max + ra_padding, 
                         int(grid_resolution * smooth_factor))
    grid_dec = np.linspace(dec_min - dec_padding, dec_max + dec_padding, 
                          int(grid_resolution * smooth_factor))
    
    # Create meshgrid
    grid_ra_mesh, grid_dec_mesh = np.meshgrid(grid_ra, grid_dec)
    
    # Interpolate values onto regular grid
    points = np.column_stack((ra_scaled, dec))
    grid_values = griddata(points, value, (grid_ra_mesh, grid_dec_mesh), 
                          method=interpolation_method, fill_value=fill_value)
    
    # Create the image plot
    extent = [ra_max + ra_padding, ra_min - ra_padding, dec_min - dec_padding, dec_max + dec_padding]
    grid_values=np.fliplr(grid_values)
    im = ax.imshow(grid_values, extent=extent, origin='lower', 
                   cmap='viridis', alpha=0.8, aspect='equal')
    
    # Set labels
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    
    # Calculate auto-range if needed
    if ymax == 0.0:
        ymin = np.nanpercentile(grid_values, 20)
        ymax = np.nanpercentile(grid_values, 80)
    
    # Set color limits
    im.set_clim(ymin, ymax)
    
    # Add colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = plt.colorbar(im, cax=cax)
    
    if label == '':
        label = 'Flux'
    cbar.set_label(label)
    
    # Set up x-axis ticks to show true RA values
    current_ticks = ax.get_xticks()
    # Filter ticks to only those within the plot range
    valid_ticks = current_ticks[(current_ticks >= extent[1]) & (current_ticks <= extent[0])]
    tick_labels = [f"{tick/cos_dec_factor:.1f}" for tick in valid_ticks]
    ax.set_xticks(valid_ticks)
    ax.set_xticklabels(tick_labels)
    
    return im

# end of new routine



def fig1(xtab,title=''):
    ymin=np.percentile(xtab['flux_ha'],5)
    ymax=np.percentile(xtab['flux_ha'],95)
    xtab['flux_sii']=xtab['flux_sii_a']+xtab['flux_sii_b']
    xtab['s2:ha']=xtab['flux_sii']/xtab['flux_ha']
    plt.figure(1,(12,8))
    plt.clf()
    # plt.figure(1,calculate_figure_size(width=12, rows=2, cols=2, xtab=xtab))
    plt.subplot(2,2,1)
    plot_one_interpolated(xtab,'flux_ha',ymin,ymax,r'H$\alpha$ Flux')
    plt.subplot(2,2,2)
    ymin/=2.
    ymax/=2.
    ymin=np.percentile(xtab['flux_sii'],5)
    ymax=np.percentile(xtab['flux_sii'],95)
    plot_one_interpolated(xtab,'flux_sii',ymin,ymax,r'[SII] Flux')
    plt.subplot(2,2,3)
    ymin=0.0
    ymax=0.8
    ymin=np.percentile(xtab['s2:ha'],5)
    ymax=np.percentile(xtab['s2:ha'],95)
    plot_one_interpolated(xtab,'s2:ha',ymin,ymax,r'[SII]:h$\alpha$')
    plt.subplot(2,2,4)
    plot_one_interpolated(xtab,'fwhm_ha',1.0,2.5,label=r'H$\alpha$  FWHM')
    if title=='':
        title='test'
    plt.suptitle(title, fontsize=16)
    plt.tight_layout()
    os.makedirs('./Snap_fig',exist_ok=True)
    plt.savefig('./Snap_fig/%s.png' % (title))
    return

lmc=262
smc=146
galaxy=0


def one_snapshot(xsum,source_name,size_arcmin=10.,keep_tmp=False,redo=True):
    
    print('OK sports fans: ', redo)
    root_fit='Snap/%s'  % source_name
    print('OK',root_fit)
    # This has to be here to get the ra and dec
    ra,dec,filenames=get_files(xsum=xsum,source_name=source_name)
    xprocess=True
    if redo==False:
        if os.path.isfile('%s.fits' % root_fit)==True:
            print('We have created  %s.fits previously' % root_fit)
            xprocess=False
        else:
            print('Could not find %s.fits'  % root_fit)


    os.makedirs('Snap',exist_ok=True)

    if xprocess==True:
        # ra,dec,filenames=get_files(xsum=xsum,source_name=source_name)
        # This only return files that are available locally.

        if len(filenames)==0:
            print('No files for one_snapshot to use for source %s at %.3f %.3f' % (source_name,ra,dec))
            return

        rss_combine_pos.do_fixed(filenames,ra, dec, pa=0, size=size_arcmin/60.,c_type='ave',outroot=root_fit,keep_tmp=keep_tmp)

    os.makedirs('./Snap_gauss',exist_ok=True)
    root_spec='./Snap_gauss/%s' % source_name
    vel=get_vel(ra,dec)

    results=lvm_gaussfit.do_all('%s.fits' % root_fit, vel=vel,outname=root_spec,xplot=False) 

    results['flux_sii']=results['flux_sii_a']+results['flux_sii_b']
    results['s2:ha']=results['flux_sii']/results['flux_ha']

    fig1(results,title=source_name)
    return





def steer(argv):
    '''
    Usage: rss_snap.py [-h] [-keep] [-redo] [-all] file  source_name
    '''

    sources=[]
    xfile=''
    keep_tmp=False
    size=10.
    xall=False
    redo=False

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-keep':
            keep_tmp=True
        elif argv[i]=='-redo':
            redo=True
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-size':
            i+=1
            size=eval(argv[i])
        elif argv[i][0]=='-':
            print('Error: Could not parse command line',argv)
            return
        elif xfile=='':
            xfile=argv[i]
        else:
            sources.append(argv[i])
        i+=1

    if xall==True:
        ztab=ascii.read(xfile)
        sources=np.unique(ztab['Source_name'])
        print('Processing all %d source names in %s' % (len(sources),xfile))

    i=0
    while i<len(sources):

        source_name=sources[i]
        print('!! Beginning: %s: %d of %d sources to process' % (source_name,i+1,len(sources)))
        print('What ', redo)

        one_snapshot(xfile,source_name,size_arcmin=size,keep_tmp=keep_tmp,redo=redo)

        print('!! Finished : %s: %d of %d sources to process' % (source_name,i+1,len(sources)))

        i+=1
        






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
