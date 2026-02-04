#!/usr/bin/env python
# coding: utf-8


'''
                    Space Telescope Science Institute

Synopsis:  

Run and then evaluate the results of SkyCorr on any two files


Command line usage (if any):

    usage: RunSkyCorr.py sci.fits sky.fits

    where sci.fits contains the science data, and sky.fits contains the sky

Description:  

    This routine will run SkyCorr on properly prepared input files, based
    on parameters given in a file called lvm_base.par.  If this file
    does not exist, one will be written, but then the program will
    exist, because certain lines in this file need to be edited as
    indicated below.  Most of the parameters in lvm_base.par can be 
    modified if desired.  

Primary routines:

    doit

Notes:

SkyCoor uses absolute paths in the 
These are the lines in lvm_base_par that (may) need to be edited::

    INST_DIR=/Users/long/Projects/sdss/skycorr/
    INPUT_OBJECT_SPECTRUM=/Users/long/Projects/sdss/skycorr_test/data/XXOBJECT.fits
    INPUT_SKY_SPECTRUM=/Users/long/Projects/sdss/skycorr_test/data/XXSKY.fits
    OUTPUT_DIR=/Users/long/Projects/sdss/skycorr_demo/output

All of the PATHS need to be full paths.  (SkyCorr can use relative paths, but they are relative to where SkyCorr is installed, not the directory in which SkyCorr is being run).  (Note The words XXOBJECT and XXSKY
should not be changed as they are used in the function update par below.)

In this case, skycorr is beeing run from the directory skycorr_demo, the data are located in skycorr_test/data/

For the routine to work, the skycorr (the executable) must be in one's path, or alternatively (but not recommened), one
one must alter one line in run_sky below to point at the executable for sky_corr.

                                       
History:

240328 ksl Coding begun

'''


import os

from astropy.io import fits,ascii
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import subprocess


xbase='''
# ----------------------------------------------------------------------------
# -------------------- INPUT PARAMETER FILE FOR SKYCORR ----------------------
# ----------------------------------------------------------------------------

# ---------------------------DIRECTORIES + FILES------------------------------

# Absolute path of skycorr installation directory
INST_DIR=/Users/long/SDSS/skycorr-1.1.2/   

# Absolute or relative (with respect to INST_DIR) path and filename of input
# object spectrum
INPUT_OBJECT_SPECTRUM=QQQ/XXOBJECT.fits

# Absolute or relative (with respect to INST_DIR) path and filename of input
# sky spectrum
INPUT_SKY_SPECTRUM=QQQ/XXSKY.fits

# Absolute or relative (with respect to INST_DIR) path and filename of output
# directory (will be created if not present; default: <INST_DIR>/output/)
OUTPUT_DIR=QQQ/output/

# Main name of diagnostic output files, extensions will be added
OUTPUT_NAME=XXOUTPUT

#------------------------------INPUT STRUCTURE--------------------------------

# Names of file columns (table) or extensions (image)
# A list of 4 labels has to be provided:
# 1: wavelength [image: NONE if dedicated extension does not exist]
# 2: flux [image: NONE if in zeroth, unnamed extension]
# 3: flux error [NONE if not present]
# 4: mask (integer: 1 = selected, 0 = rejected;
#          float:   0. = selected, otherwise rejected) [NONE if not present]
COL_NAMES=WAVE FLUX ERROR NONE

# Error relative to mean if no error column is provided (default: 0.01)
DEFAULT_ERROR=0.01

# Multiplicative factor to convert wavelength to micron
# e.g.: wavelength unit = A -> WLG_TO_MICRON = 1e-4
WLG_TO_MICRON=1e-4

# Wavelengths in vacuum (= vac) or air (= air)
VAC_AIR=air


# ----------------------------------------------------------------------------
# ------------------------- EXPERT MODE PARAMETERS ---------------------------
# ----------------------------------------------------------------------------

# ------------------------------FITS KEYWORDS---------------------------------

# FITS keyword of sky spectrum for Modified Julian Day (MJD) or date in years
# (default: MJD-OBS; optional parameter for value: DATE_VAL)
DATE_KEY=MJD

# FITS keyword of sky spectrum for UTC time in s
# (default: TM-START; optional parameter for value: TIME_VAL)
TIME_KEY=UTC

# FITS keyword of sky spectrum for telescope altitude angle in deg
# (default: ESO TEL ALT; optional parameter for value: TELALT_VAL)
TELALT_KEY=ALT

# ---------------------------REQUIRED INPUT DATA------------------------------

# Airglow line list
# Required directory: <INST_DIR>/sysdata/
LINETABNAME=airglow_groups.dat

# File for airglow scaling parameters
# Required directory: <INST_DIR>/sysdata/
VARDATNAME=airglow_var.dat

# FTP address (supplemented by "ftp://") for folder with monthly averages of
# solar radio flux at 10.7 cm
SOLDATURL=ftp.geolab.nrcan.gc.ca/data/solar_flux/monthly_averages

# File with monthly averages of solar radio flux at 10.7 cm
# Required directory: SOLDATURL or <INST_DIR>/sysdata/
SOLDATNAME=solflux_monthly_average.txt

# Solar radio flux at 10.7 cm:
# Positive value in sfu (= 0.01 MJy) or -1 [default] for corresponding monthly
# average from http://www.spaceweather.gc.ca. Download only if local file in
# <INST_DIR>/sysdata/ does not contain required month.
SOLFLUX=-1

# ---------------------------LINE IDENTIFICATION------------------------------

# Initial estimate of line FWHM [pixel]
FWHM=5.0

# Variable line width (linear increase with wavelength)? -- 1 = yes; 0 = no
VARFWHM=0

# Relative FWHM convergence criterion (default: 1e-2)
LTOL=1e-2

# Minimum distance to neighbouring lines for classification as isolated line:
# <MIN_LINE_DIST> * <FWHM> [pixel]
MIN_LINE_DIST=2.5

# Minimum line peak flux for consideration of lines from airglow line list:
# <FLUXLIM> * <median flux of identified lines>
# Automatic search -> FLUXLIM = -1 (default)
FLUXLIM=-1

# ---------------------------FITTING OF SKY LINES-----------------------------

# Relative chi^2 MPFIT convergence criterion (default: 1e-3)
FTOL=1e-3

# Relative parameter MPFIT convergence criterion (default: 1e-3)
XTOL=1e-3

# Relative chi^2 convergence criterion for iterative improvement of
# wavelength grid (default: 1e-3)
WTOL=1e-3

# Maximum degree of Chebyshev polynomial for wavelength grid correction:
# -1 = no correction
#  0 = linear term (coef. = 1) is also considered but not fitted
#  7 = default
CHEBY_MAX=7

# Minimum degree of Chebyshev polynomial for wavelength grid correction.
# CHEBY_MIN <= CHEBY_MAX:
# - Iterative increase of polynomial degree at least until CHEBY_MIN
#   (default: 3).
# - Procedure stops if chi^2 gets worse or CHEBY_MAX is reached.
# - Results of degree with best chi^2 are taken.
# CHEBY_MIN > CHEBY_MAX:
# - Iterative increase of polynomial degree until CHEBY_MAX is reached.
# - Results of degree CHEBY_MAX are taken.
CHEBY_MIN=3

# Initial constant term for wavelength grid correction (shift relative to half
# wavelength range)
CHEBY_CONST=0.

# Type of rebinning:
# 0 = simple rebinning (summation of pixel fractions)
# 1 = convolution with asymmetric, damped sinc kernel [default]
REBINTYPE=1

# Minimum relative weight of the strongest line group of a pixel for
# including a pixel in the line fitting procedure (default: 0.67)
WEIGHTLIM=0.67

# Sigma limit for excluding outliers (e.g. object emission lines) from
# estimate of group flux correction factors (default: 15.)
SIGLIM=15.

# Lower relative uncertainty limit for the consideration of a line group for
# the fitting procedure. The value is compared to the sigma-to-mean ratio of
# the group-specific flux correction factors of the initial estimate
# (default: 0. -> include all fittable line groups).
FITLIM=0.

# ---------------------------------PLOTTING-----------------------------------

# Diagnostic gnuplot plots:
# Options for output on screen:
# W - wxt terminal
# X - x11 terminal
# N - no screen output [default]
# NOTE: An illustration of the sky subtraction quality is plotted into a PS
#       file in the OUTPUT_DIR folder in any case.
PLOT_TYPE=N
'''

def write_par():
    '''
    Write out a base parameter file
    '''
    cwd=os.getcwd()
    print('lvm_base.par does not seem to exist, so writing one, but this needs to be edited to have correct paths')
    x=open('lvm_base.par','w')
    zbase=xbase.replace('QQQ',cwd)
    x.write(zbase)
    x.close()
    print(__doc__)
    exit()


def update_par(sci='sctest_57454_11_52',sky='sctest_57458_10_53',outroot=''):
    '''
    Read a parameter file with everything but the science and sky images
    defined, and write out a file that can be run by SkyCorr.
    
    The inputs are the rootnames of the data and sky files that you want to
    run sky corr on.
    '''
    try:
        x=open('lvm_base.par')
    except:
        write_par()

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
    if parfile==None:
        return None
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
    
    plt.figure(1,(12,6))
    plt.clf()
    plt.subplot(2,2,1)
    plt.semilogy(sci['WAVE'],sci['FLUX'],label='Sci')
    plt.xlim(3600,9800)
    plt.ylim(2e-15,2e-12)
    plt.text(4000,3e-15,outroot)
    plt.legend()
    plt.tight_layout()

    plt.subplot(2,2,2)
    plt.semilogy(sky['WAVE'],sky['FLUX'],label='Sky')
    plt.xlim(3600,9800)
    plt.ylim(2e-15,2e-12)
    plt.legend()
    plt.tight_layout()

    plt.subplot(2,2,3)
    plt.plot(sci['WAVE'],sci['FLUX']-sky['FLUX'],label='Simple Sci-Sky')
    plt.plot([3600,9800],[6e-15,6e-15],'k:')
    plt.plot([3600,9800],[-6e-15,-6e-15],'k:')
    plt.xlim(3600,9800)
    plt.ylim(-1e-14,1e-14)
    plt.legend()
    plt.tight_layout()

    plt.subplot(2,2,4)
    plt.plot(10000*fit['lambda'],fit['scflux'],label='Skycorr-Subtracted')
    plt.plot([3600,9800],[6e-15,6e-15],'k:')
    plt.plot([3600,9800],[-6e-15,-6e-15],'k:')
    plt.xlim(3600,9800)
    plt.ylim(-1e-14,1e-14)
    plt.plot([3600,9800],[6e-15,6e-15],'k:')
    plt.plot([3600,9800],[-6e-15,-6e-15],'k:')
    plt.legend()
    plt.tight_layout()

    plt.savefig('%s_skycorr.png' % outroot)

    plt.figure(2,(12,6))
    plt.clf()
    plt.subplot(2,2,1)
    plt.semilogy(sci['WAVE'],sci['FLUX'],label='Sci')
    plt.semilogy(10000*fit['lambda'],fit['cflux'],label='Sci-cont')
    plt.xlim(3600,9800)
    plt.ylim(2e-15,2e-12)
    plt.text(4000,3e-15,outroot)
    plt.legend()
    plt.tight_layout()



    plt.subplot(2,2,2)
    plt.semilogy(sky['WAVE'],sky['FLUX'],label='Sky')
    plt.semilogy(10000*fit['lambda'],fit['mcflux'],label='Sky-cont')
    plt.xlim(3600,9800)
    plt.ylim(2e-15,2e-12)
    plt.legend()
    plt.tight_layout()

    plt.subplot(2,2,3)

    plt.semilogy(10000*fit['lambda'],fit['cflux'],label='Sci-cont')
    plt.semilogy(10000*fit['lambda'],fit['mcflux'],label='Sky-cont')
    plt.semilogy(10000*fit['lambda'],np.abs(fit['cflux']-fit['mcflux']),label='|Sci-Sky-cont|')
    plt.plot([3600,9800],[6e-15,6e-15],'k:')
    plt.ylim(1e-16,1e-12)
    plt.legend()
    plt.tight_layout()

    plt.subplot(2,2,4)
    plt.plot(10000*fit['lambda'],fit['lflux']-fit['mlflux'],label='Diff-lines')
    plt.ylim(-1e-14,1e-14)
    plt.plot([3600,9800],[6e-15,6e-15],'k:')
    plt.plot([3600,9800],[-6e-15,-6e-15],'k:')
    plt.xlim(3600,9800)
    plt.legend()
    plt.tight_layout()
    plt.savefig('%s_sep.png' % (outroot))

    
    
    
    # plt.figure(3,(6,3))
    # plt.clf()
    # plt.plot(10000*fit['lambda'],fit['cflux'],label='Sci-cont')
    # plt.plot(10000*fit['lambda'],fit['mcflux'],label='Sky-cont')
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig('%s_cont.png' % outroot)
    
    # plt.figure(4,(6,8))
    # plt.clf()
    # plt.subplot(311)
    # plt.plot(10000*fit['lambda'],fit['lflux'],label='Sci-lines')
    # plt.legend()
    # plt.tight_layout()
    # plt.subplot(312)
    # plt.plot(10000*fit['lambda'],fit['mlflux'],label='Sky-lines')
    # plt.legend()
    # plt.tight_layout()
    # plt.subplot(313)
    # plt.plot(10000*fit['lambda'],fit['lflux']-fit['mlflux'],label='Diff-lines')
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(outroot+'.png')
    
    
def do_one(xsci='Sci_4852',xsky='SkyE_4852'):
    outroot=run_sky(sci=xsci,sky=xsky,outroot='',print_output=False)
    if outroot==None:
        return
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
