#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This routine is intended for running skycorr.  Based using
lvm_base.par.  Run from the command line it takes two 
skycorr ready files, and then runs sky corr.  


Command line usage (if any):

    usage: RunSkyCorr.py source sky

Description:  

Primary routines:

    doit

Notes:
                                       
History:

220402 ksl Coding begun

'''


# # Run and then evaluate the results of SkyCorr on any two files



import os

from astropy.io import fits,ascii
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import subprocess


xxx='''
# ----------------------------------------------------------------------------
# -------------------- INPUT PARAMETER FILE FOR SKYCORR ----------------------
# ----------------------------------------------------------------------------

# ---------------------------DIRECTORIES + FILES------------------------------

# Absolute path of skycorr installation directory
INST_DIR=/Users/long/Projects/sdss/skycorr/skycorr-1.1.2/

# Absolute or relative (with respect to INST_DIR) path and filename of input
# object spectrum
INPUT_OBJECT_SPECTRUM=/Users/long/Projects/sdss/skycorrdata21/XXOBJECT.fits

# Absolute or relative (with respect to INST_DIR) path and filename of input
# sky spectrum
INPUT_SKY_SPECTRUM=/Users/long/Projects/sdss/skycorr/data21/XXSKY.fits

# Absolute or relative (with respect to INST_DIR) path and filename of output
# directory (will be created if not present; default: <INST_DIR>/output/)
OUTPUT_DIR=/Users/long/Projects/sdss/skycorr/test2204/output

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
COL_NAMES=NONE NONE NONE NONE

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


def update_par(sci='sctest_57454_11_52',sky='sctest_57458_10_53',outroot=''):
    '''
    Read a parameter file with everything but the scince and sky images
    defined, and write out a file that can be run by SkyCorr
    '''
    try:
        x=open('lvm_base.par')
        lines=x.readlines()
        x.close()
    except:
        x=open('lvm_base.par','w')
        x.write(xxx)
        x.close()
        print('lvm_bas.par did not exist, so writing it out for future work')
    
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
    The routione below assumes skycorr is in one's Path
    '''
    
    parfile=update_par(sci,sky,outroot)
    print(parfile)
    result=subprocess.run(['skycorr', parfile],capture_output=True,text=True)
    if print_output:
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)
    outroot=parfile.replace('.par','')
    return outroot
    

def summary_plots(outroot):
    
    x=fits.open('./output/%s_sci.fits' % outroot)
    sci=Table(x[1].data)
    x.close()
    
    x=fits.open('./output/%s_sky.fits' % outroot)
    sky=Table(x[1].data)
    x.close()
    
    x=fits.open('./output/%s_fit.fits' % outroot )
    fit=Table(x[1].data)
    x.close()
    
    plt.figure(1,(6,8))
    plt.clf()
    plt.subplot(311)
    plt.plot(sci['LAMBDA'],sci['FLUX'],label='Sci')
    plt.legend()
    plt.tight_layout()
    plt.subplot(312)
    plt.plot(sky['LAMBDA'],sky['FLUX'],label='Sky')
    plt.legend()
    plt.tight_layout()
    plt.subplot(313)
    plt.plot(sci['LAMBDA'],sci['FLUX']-sky['FLUX'],label='Sci-Sky')
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
    
    

def doit(sci,sky):
    outroot=run_sky(sci,sky)
    print('OK ',outroot)
    summary_plots(outroot)
    


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)==3:
        doit(sys.argv[1],sys.argv[2])
    else:
        print(__doc__)
