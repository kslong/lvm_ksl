#!/usr/bin/env python
# coding: utf-8


'''
                    Space Telescope Science Institute

Synopsis:  

Run SkyCorr on one set of data


Command line usage (if any):

    usage: RunSkyCorr.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

240328 ksl Coding begun

'''
# # Run and then evaluate the results of SkyCorr on any two files
# 
# This is a simple little notebook that contains a routine to run SkyCorr on any two files, assuming you have edited the
# SkyCorr parameter file lvm_base_par file to match your directory structure.
# 
# **These are the lines in lvm_base_par that need to be edited**
# 
# ````
# INST_DIR=/Users/long/Projects/sdss/skycorr/
# INPUT_OBJECT_SPECTRUM=/Users/long/Projects/sdss/skycorr_test/data/XXOBJECT.fits
# INPUT_SKY_SPECTRUM=/Users/long/Projects/sdss/skycorr_test/data/XXSKY.fits
# OUTPUT_DIR=/Users/long/Projects/sdss/skycorr_demo/output
# ````
# 
# All of the PATHS need to be full paths.  (SkyCorr can use relative paths, but they are relative to where SkyCorr is installed, not the directory in which SkyCorr is being run).  (Note The words XXOBJECT and XXSKY 
# should not be changed as they are used in the function update par below.)  
# 
# In this case, skycorr is beeing run from the directory skycorr_demo, the data are located in skycorr_test/data/
# 
# **Additionally, to make the script work, one must alter one line in run_sky below to point at the executable for sky_corr**
# 
# Note: Aside from providing a simple way to run sky_corr for the first time, this scipt is intended to be easily adaptable to run sky_corr multiple times with various choices for the science and sky data.  The script assumes that the individula fits files have 



import os

from astropy.io import fits,ascii
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import subprocess



def update_par(sci='sctest_57454_11_52',sky='sctest_57458_10_53',outroot=''):
    '''
    Read a parameter file with everything but the science and sky images
    defined, and write out a file that can be run by SkyCorr.
    
    The inputs are the rootnames of the data and sky files that you want to
    run sky corr on.
    '''
    x=open('lvm_base.par')
    lines=x.readlines()
    x.close()
    
    # make sure root names have been entered
    
    sci=sci.replace('.fits','')
    sky=sky.replace('.fits','')
    
    if outroot=='':
        outroot='%s_%s' % (sci,sky)
    
    outname=outroot+'.par'
    x=open(outname,'w')
    for line in lines:
        line=line.replace('XXOBJECT',sci)
        line=line.replace('XXSKY',sky)
        line=line.replace('XXOUTPUT',outroot)
        x.write(line)
    x.close()

    return outname
    
    
def run_sky(sci='sctest_57454_11_52',sky='sctest_57458_10_53',outroot='',print_output=False):
    '''
    Update the parameter file and run skycorr, sci and sky are the rootnames of the files
    that one wants to run skycorr on, and then actually run sky_coor
    '''
    parfile=update_par(sci,sky,outroot)
    print(parfile)
    # The next line must be edited to point to the binary that runs sky corr
    result=subprocess.run(['skycorr', parfile],capture_output=True,text=True)
    if print_output:
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)
    outroot=parfile.replace('.par','')
    return outroot
    

def summary_plots(outroot):

    filename='./output/%s_sci.fits' % outroot
    try:
        x=fits.open(filename)
    except:
        print('Error: Could not open %s' % filename)
        return
        
    sci=Table(x[1].data)

    try:
        xfile=filename.replace('sci.fits','sky.fits')
        x=fits.open(xfile)
    except:
        print('Error: Could not open %s' % xfile)
        return
    
    sky=Table(x[1].data)
    x.close()
    
    try:
        xfile=filename.replace('sci.fits','fit.fits')
        x=fits.open(xfile )
    except:
        print('Error: Could not open %s' % xfile)
        return    
    fit=Table(x[1].data)
    x.close()

    sci.info()
    
    plt.figure(1,(6,8))
    plt.clf()
    plt.subplot(311)
    plt.plot(sci['WAVE'],sci['FLUX'],label='Sci')
    plt.legend()
    plt.tight_layout()
    plt.subplot(312)
    plt.plot(sky['WAVE'],sky['FLUX'],label='Sky')
    plt.legend()
    plt.tight_layout()
    plt.subplot(313)
    plt.plot(sci['WAVE'],sci['FLUX']-sky['FLUX'],label='Sci-Sky')
    plt.legend()
    plt.tight_layout()
    
    
    plt.figure(2,(6,3))
    plt.clf()
    plt.plot(fit['lambda'],fit['scflux'],label='Skycorr-Subtracted')
    plt.legend()
    plt.tight_layout()
    
    plt.figure(3,(6,3))
    plt.clf()
    plt.plot(fit['lambda'],fit['cflux'],label='Sci-cont')
    plt.plot(fit['lambda'],fit['mcflux'],label='Sky-cont')
    plt.legend()
    plt.tight_layout()
    
    plt.figure(4,(6,8))
    plt.clf()
    plt.subplot(311)
    plt.plot(fit['lambda'],fit['lflux'],label='Sci-lines')
    plt.legend()
    plt.tight_layout()
    plt.subplot(312)
    plt.plot(fit['lambda'],fit['mlflux'],label='Sky-lines')
    plt.legend()
    plt.tight_layout()
    plt.subplot(313)
    plt.plot(fit['lambda'],fit['lflux']-fit['mlflux'],label='Diff-lines')
    plt.legend()
    plt.tight_layout()
    plt.savefig(outroot+'.png')
    
    
def do_one(xsci='Sci_4852',xsky='SkyE_4852'):
    outroot=run_sky(sci=xsci,sky=xsky,outroot='',print_output=False)
    summary_plots(outroot=outroot)

def steer(argv):
    '''
    The simple directs the flow
    '''

    science=''
    sky=''

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i][0]=='-':
            print('Error: Cannot parse command line: ',argv)
            return
        elif science=='':
            science=argv[i]
        elif sky=='':
            sky=argv[i]
        else:
            print('Error: Too many arguments in command line: ',argv)
            return
        i+=1

    if sky=='':
        print('Error: not enoubh argmeents: ',argv)

    do_one(science.replace('.fits',''),sky.replace('.fits',''))
    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)         
    else:
        print (__doc__)
