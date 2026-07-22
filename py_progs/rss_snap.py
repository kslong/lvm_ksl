#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Create one or more subimages of specific regions observed with 
lvm using a file derived from DRP_All and a set of postions 
containing the exposures asscoated with each object.


Command line usage::

    rss_snap.py [-h] [-keep] [-redo] [-all] [-med] [-size arcmin] [-lmc] [-smc] [-vel vel] xfile source_name

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

-med
    Use median for combining remapped images. By default, mean is used.

-lmc
    Applies the LMC radial velocity (~262 km/s) when fitting nebular lines.

-smc
    Applies the SMC radial velocity (~146 km/s) when fitting nebular lines.

-vel vel
    Applies an arbitrary radial velocity (km/s) when fitting nebular lines.
    The default is 0.

Description:

Primary routines:

    do_one

Notes:
                                       
History:

250704 ksl Coding begun

'''

import os
import re
import sys
from lvm_ksl import rss_combine_pos
from lvm_ksl import lvm_gaussfit
from lvm_ksl import radec_plot


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
    elif os.path.isdir(XMUSKIE):
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






def _usage_from_doc(doc):
    '''
    __doc__ truncated just before a line consisting of "History:" (or
    "History::"/"Version History" -- whitespace/colon-insensitive), so
    -h stays short even as that section grows -- without hand-
    duplicating the Synopsis/Options text in a second string.  Anchored
    to a whole line (not a bare substring search) so it can't misfire on
    "History:" appearing mid-sentence, and returns doc unchanged if no
    such line is present.
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


def plot_one(xtab, var='flux_ha', ymin=0, ymax=0, label='', ax=None, marker_size=None):
    '''
    Thin wrapper around radec_plot.plot_scatter(), kept for backward
    compatibility with call sites in this module (fig1) and any
    notebooks already using rss_snap.plot_one() -- see
    radec_plot.plot_scatter() for the general-purpose, RA/Dec-column-
    agnostic implementation this now delegates to.
    '''
    return radec_plot.plot_scatter(xtab, var, ymin=ymin, ymax=ymax, label=label,
                                    ax=ax, marker_size=marker_size)


def plot_one_interpolated(xtab, var='flux_ha', ymin=0, ymax=0, label='', ax=None,
                           grid_resolution=100, interpolation_method='linear',
                           mask_radius=None):
    '''
    Thin wrapper around radec_plot.plot_interpolated(), kept for
    backward compatibility with call sites in this module (fig1) and
    any notebooks already using rss_snap.plot_one_interpolated() -- see
    radec_plot.plot_interpolated() for the general-purpose, RA/Dec-
    column-agnostic implementation this now delegates to.

    Note: the old fill_value/smooth_factor parameters are gone --
    fill_value was always np.nan in practice (masking is now handled
    by mask_radius instead), and smooth_factor was just a multiplier
    on grid_resolution a caller can apply directly.
    '''
    return radec_plot.plot_interpolated(xtab, var, ymin=ymin, ymax=ymax, label=label,
                                         ax=ax, grid_resolution=grid_resolution,
                                         interpolation_method=interpolation_method,
                                         mask_radius=mask_radius)



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


def one_snapshot(xsum,source_name,size_arcmin=10.,keep_tmp=False,redo=True,c_type='ave',vel=0.):

    print('OK sports fans: ', redo)
    root_fit='Snap/%s.%s'  % (source_name,c_type)
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

        rss_combine_pos.do_fixed(filenames,ra, dec, pa=0, size=size_arcmin/60.,c_type=c_type,outroot='Snap/%s' % source_name,keep_tmp=keep_tmp)

    os.makedirs('./Snap_gauss',exist_ok=True)
    root_spec='./Snap_gauss/%s' % source_name

    results=lvm_gaussfit.do_all('%s.fits' % root_fit, vel=vel,outname=root_spec,xplot=False)

    results['flux_sii']=results['flux_sii_a']+results['flux_sii_b']
    results['s2:ha']=results['flux_sii']/results['flux_ha']

    fig1(results,title=source_name)
    return





def steer(argv):
    '''
    Usage: rss_snap.py [-h] [-keep] [-redo] [-all] [-med] [-size arcmin] [-lmc] [-smc] [-vel vel] xfile source_name
    '''

    sources=[]
    xfile=''
    keep_tmp=False
    size=10.
    xall=False
    redo=False
    c_type='ave'
    vel=0.

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(_usage_from_doc(__doc__))
            return
        elif argv[i]=='-keep':
            keep_tmp=True
        elif argv[i]=='-redo':
            redo=True
        elif argv[i]=='-all':
            xall=True
        elif argv[i][0:4]=='-med':
            c_type='med'
        elif argv[i]=='-size':
            i+=1
            size=eval(argv[i])
        elif argv[i]=='-lmc':
            vel=lmc
        elif argv[i]=='-smc':
            vel=smc
        elif argv[i]=='-vel':
            i+=1
            vel=eval(argv[i])
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

        one_snapshot(xfile,source_name,size_arcmin=size,keep_tmp=keep_tmp,redo=redo,c_type=c_type,vel=vel)

        print('!! Finished : %s: %d of %d sources to process' % (source_name,i+1,len(sources)))

        i+=1
        






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
