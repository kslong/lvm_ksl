#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Create one or more subimages of specific regions observed with 
lvm using a file derived from DRP_All and a set of postions 
containing the exposures asscoated with each object.


Command line usage (if any):

    Usage: rss_snap.py [-h] [-keep] xfile  source_name ...

Description:  

    where:
        -h prints this documentation
        -keep retains the temporary files from individual exposures in the 
            directory xtmp. Without this switch the temporary files are deleted
        xfile is a file containing the sources names and associated exposures
        source_name is one or more of the source_names in xfile

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


def fig1(xtab,marker_size=50,plot_root=''):
    ymin=np.percentile(xtab['flux_ha'],5)
    ymax=np.percentile(xtab['flux_ha'],75)
    plt.figure(1,(12,8))
    # plt.figure(1,calculate_figure_size(width=12, rows=2, cols=2, xtab=xtab))
    plt.subplot(2,2,1)
    plot_one(xtab,'flux_ha',ymin,ymax,r'H$\alpha$ Flux',marker_size=marker_size)
    plt.subplot(2,2,2)
    plot_one(xtab,'flux_sii',ymin,ymax,r'[SII] Flux',marker_size=marker_size)
    plt.subplot(2,2,3)
    plot_one(xtab,'s2:ha',0.2,0.6,r'[SII]:h$\alpha$',marker_size=marker_size)
    plt.subplot(2,2,4)
    plot_one(xtab,'fwhm_ha',1.4,2.5,label=r'H$\alpha$  FWHM',marker_size=marker_size)
    plt.tight_layout()
    if plot_root=='':
        plot_root='test'
    plt.savefig('%s.png' % (plot_root))

def one_snapshot(xsum,source_name,vel,keep_tmp=False):

        ra,dec,filenames=get_files(xsum=xsum,source_name=source_name)
        if len(filenames)==0:
            print('No files for one_snapshot to use for source %s' % (source_name))
            return


        os.makedirs('Snap',exist_ok=True)

        root_fit='Snap/%s'  % source_name
        rss_combine_pos.do_fixed(filenames,ra, dec, pa=0, size=20./60.,c_type='ave',outroot=root_fit,keep_tmp=keep_tmp)

        os.makedirs('./Snap_gauss',exist_ok=True)
        root_spec='./Snap_gauss/%s' % source_name

        results=lvm_gaussfit.do_all('%s.fits' % root_fit, vel=lmc,outname=root_spec,xplot=False) 

        results['flux_sii']=results['flux_sii_a']+results['flux_sii_b']
        results['s2:ha']=results['flux_sii']/results['flux_ha']

        os.makedirs('./Snap_fig',exist_ok=True)
        root_plot='./Snap_fig/%s' % source_name
        fig1(results,plot_root=root_plot)
        return



lmc=262
smc=146
galaxy=0


def steer(argv):
    '''
    Usage: rss_snap.py [-h] [-keep] xfile  source_name
    '''

    sources=[]
    xfile=''
    vel=lmc
    keep_tmp=False

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-keep':
            keep_tmp=True
        elif argv[i][0]=='-':
            print('Error: Could not parse command line',argv)
            return
        elif xfile=='':
            xfile=argv[i]
        else:
            sources.append(argv[i])
        i+=1

    i=0
    while i<len(sources):
        source_name=sources[i]
        one_snapshot(xfile,source_name,vel,keep_tmp)

        i+=1
        






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
