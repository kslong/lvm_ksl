#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Rerun the quicksky portion of the pipeline on an
CFrame file


Command line usage (if any):

    usage: RunSky.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

240419 ksl Coding begun

'''

import sys
from astropy.io import ascii,fits
import matplotlib.pyplot as plt
from astropy.io import fits,ascii
from astropy.table import Table
import numpy as np
from glob import glob
from lvmdrp.functions.skyMethod import quick_sky_subtraction
import os
import pathlib
from lvm_ksl import sky_plot



def find_skytab(exp):
    '''
    Locate a skytable-brz file in the lvm directory
    structure
    '''
    xfile = 'lvm-skytable-brz-%08d.fits' % exp
    sas_base_dir = os.environ.get('SAS_BASE_DIR')
    if sas_base_dir is None:
        print("SAS_BASE_DIR environment variable is not set.")
        return None
    files=list(pathlib.Path(sas_base_dir).rglob('*'+xfile))

    return files[0]
    

def get_sky_tab(filename):
    '''
    Retrive a sky table and split it into 
    tables
    '''
    x=fits.open(filename)
    xsci=Table(x['SCI'].data)
    xsky_e=Table(x['SKYE'].data)
    xsky_w=Table(x['SKYW'].data)
    xsky_e_super=Table(x['SKYE_SUPER'].data)
    xsky_w_super=Table(x['SKYW_SUPER'].data)
    return xsci,xsky_e,xsky_w,xsky_e_super,xsky_w_super


def plot_one(xtab):
    '''
    Plot a single skytab table
    '''


    plt.figure(1,(6,6))
    plt.clf()
    plt.subplot(3,1,1)
    plt.semilogy(xtab['WAVE'],xtab['FLUX'])
    plt.tight_layout()
    try:
        plt.semilogy(xtab['WAVE'],xtab['CONT'])
    except:
        print('No Cont')
        return
    ymin,ymax=plt.ylim()
    plt.xlim(3500,9800)
    plt.subplot(3,1,2)
    plt.semilogy(xtab['WAVE'],xtab['LINES'],label='Lines')
    plt.legend()
    plt.ylim(ymin,ymax)
    plt.xlim(3500,9800)
    plt.tight_layout()
    plt.subplot(3,1,3)
    plt.plot(xtab['WAVE'],xtab['FLUX']-xtab['CONT'])
    plt.plot([3600,9800],[0,0],'r')
    plt.ylim(-.1e-13,0.1e-13)
    plt.xlim(3500,9800)
    plt.tight_layout()
                       


def steer(argv):
    '''
    This is mainly a steering routine
    '''

    i=1
    exposure=-1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i][0]=='-':
            print('Error: Poorly formated command line: ',argv)
        elif exposure==-1.:
            exposure=eval(argv[i])
        else:
            print('Error: Poorly formated command line: ',argv)
        i+=1


    if exposure==-1:
        return
    in_cframe='data/lvmCFrame-%08d.fits' % (exposure)
    out_sframe='foo.fits'
    sky_method='farlines_nearcont'

    print('Working on exposure %d' % exposure)
    filename = find_skytab(exposure)
    x=fits.open(filename)
    x.info()


    quick_sky_subtraction(in_cframe, out_sframe,
        skip_subtraction=False, skymethod=sky_method)


    sci,sky_e,sky_w,sky_e_super,sky_w_super=get_sky_tab(filename)

    plot_one(sci)
    plot_one(sky_e)
    plot_one(sky_w)
    plot_one(sky_e_super)
    plot_one(sky_w_super)

    sky_plot.eval_qual('foo.fits')
    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
