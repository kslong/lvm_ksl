#!/usr/bin/env python
# coding: utf-8


'''
                    Space Telescope Science Institute

Synopsis:  

Create a standard plot of a spectrum that has been
extracted from the LVM data. The spectrum should
have been sky subtracted.  


Command line usage (if any):

    usage: PlotSpec.py [-h] [-frac 0.1] [-min ymin] [-max ymax] [-med] [-delta 1e-15]
                       [-mode sep_back] file [files ...]

    Plots are written to Overview_Plot/<basename>.overview.png.

    Default mode (no -mode flag):
        One or more spectrum files are given.  Each is plotted independently.
        If a file contains a BACK_FLUX column (as produced by GetRegSpec.py),
        the background is overlaid in black (alpha 0.4); no further subtraction
        is done because FLUX in that case is already the background-subtracted
        spectrum.

    -mode sep_back:
        Exactly two files are expected: a source spectrum and a separate
        background spectrum.  The background FLUX is subtracted from the
        source FLUX before plotting.

    Scaling options (mutually exclusive; last one wins if combined):
    -frac   autoscale upper limit to frac * max(FLUX) per panel (default 0.1)
    -max    fix the upper y-limit in all panels (also selects fixed-scale mode)
    -min    fix the lower y-limit in all panels
    -med    centre each panel on the median FLUX in that wavelength range,
            with limits median +/- delta
    -delta  half-range for -med mode (default 3e-15)

Description:  

Primary routines:

    doall

Notes:

    The various panels are currently autoscaled.  The
    maximum values plotted are controlled by the 
    frac parameter.
    
History:

240607 ksl Coding begun

'''



import os
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.table import Table,vstack, hstack
from astropy.io import fits
import numpy as np

def do_one_region(spectab,wmin=3600,wmax=4100,frac=0.1):
    extra=10
    xx=spectab[spectab['WAVE']>wmin-extra]
    xx=xx[xx['WAVE']<wmax+extra]
    mask=np.isfinite(xx['FLUX'])
    xx=xx[mask]
    if 'BACK_FLUX' in xx.colnames:
        plt.plot(xx['WAVE'],xx['BACK_FLUX'],'k',alpha=0.4)
    plt.plot(xx['WAVE'],xx['FLUX'])
    plt.xlim(wmin-extra,wmax+extra)
    if frac<1.0 and len(xx)>0:
        ymax=np.max(xx['FLUX'])
        ymin=np.median(xx['FLUX'])-0.05*ymax
        plt.ylim(ymin,frac*ymax)
    return


def do_one_region_fixed(spectab,wmin=3600,wmax=4100,ymin=0, ymax=1e-13):
    extra=10
    xx=spectab[spectab['WAVE']>wmin-extra]
    xx=xx[xx['WAVE']<wmax+extra]
    mask=np.isfinite(xx['FLUX'])
    xx=xx[mask]
    if 'SOURCE_FLUX' in xx.colnames:
        plt.plot(xx['WAVE'],xx['SOURCE_FLUX'],'k',alpha=0.2)
    if 'BACK_FLUX' in xx.colnames:
        plt.plot(xx['WAVE'],xx['BACK_FLUX'],'k',alpha=0.4)
    plt.plot(xx['WAVE'],xx['FLUX'])
    plt.xlim(wmin-extra,wmax+extra)
    plt.ylim(ymin,ymax)
    return


def do_one_region_med(spectab,wmin=3600,wmax=4100,med_delta=1e-15):
    extra=10
    xx=spectab[spectab['WAVE']>wmin-extra]
    xx=xx[xx['WAVE']<wmax+extra]
    mask=np.isfinite(xx['FLUX'])
    xx=xx[mask]
    if len(xx)==0:
        return
    ymed=np.median(xx['FLUX'])
    if 'BACK_FLUX' in xx.colnames:
        plt.plot(xx['WAVE'],xx['BACK_FLUX'],'k',alpha=0.4)
    plt.plot(xx['WAVE'],xx['FLUX'])
    plt.xlim(wmin-extra,wmax+extra)
    plt.ylim(ymed-0.5*med_delta,ymed+med_delta)
    return

def xmark(line='[SIII]',w=9069.3,frac=0.90):

    ymin,ymax=plt.ylim()
    plt.text(w,frac*ymax,line,ha='center',va='top',color='red',clip_on=True)

def do_lines():
    '''
    These are air waveleengs or should be
    '''
    xmark('[OII]', 3728)
    xmark('[OIII]',4363.15,.6)
    xmark('[OIII]',5006.843)
    xmark('[OI]',6300.3)
    xmark('[OI]',6363.77)
    xmark('[SIII]',9069.3)
    xmark('[SIII]',9532.3)
    xmark('HI',8862)
    xmark('HI',9014)
    xmark('[ArIII]',7135.67)
    xmark('[CaII]',7291)
    xmark('[CaII]',7323.72)
    xmark('CaII',3933.61,0.6)
    xmark('[SII]',6720)
    xmark(r'H$\alpha$/[NII]',6563)
    xmark(r'H$\beta$',4861)
    xmark('[NeIII]',3869.06)
    xmark('[NeIII]',3967.79)
    xmark(r'H$\delta$',4104.73)
    xmark(r'H$\gamma$',4340.46)
    xmark('HI', 9229.02)
    xmark('HeII',4199.83)
    xmark('HeII',4685.74)
    xmark('HeI',5875.6)
    xmark('HeI', 6678.2)
    xmark('[NII]',5754.5)
    xmark('[SIII]',6310.05,0.6)
    xmark('FeVII',4987.17,0.6)
    xmark('FeVII',5157.45,0.6)
    xmark('[FeII]',8617.0)
    xmark('[FeII]',8891.9)
    xmark('[SII]',4076.27)
    xmark('HeI',3888.6,0.6)
    xmark('HI',3889.05,0.4)
    xmark('HI',3835.38)
    xmark('HI',3770.6)
    xmark('HI',3797.9)
    xmark('[MgI]',4566.78) 
    xmark('[OII]',7320.24,0.6)
    return



def do_all(xtab,ptype='scale',ymin=0.0,ymax=1e-14,frac=0.1,med_delta=3e-15,title=''):
    '''
    Create the figure
    '''
    plt.close(2)
    plt.figure(2,(8,12))

    wmin=3600
    wmax=9559
    delta=750
    nmax=int((wmax-wmin)/delta)+1
    i=0
    while i<nmax:
        plt.subplot(nmax,1,i+1)
        wwmin=wmin+i*delta
        wwmax=wwmin+delta
        if ptype=='scale':
            do_one_region(xtab,wwmin,wwmax,frac)
        elif ptype=='fixed':
            do_one_region_fixed(xtab,wwmin,wwmax,ymin,ymax)
        elif ptype=='med':
            do_one_region_med(xtab,wwmin,wwmax,med_delta)
        else:
            print('Error: Indecipherable plot type: ',ptype)
            return

        do_lines()
        i+=1
    if title:
        plt.suptitle(title,fontsize=12)
    plt.tight_layout()


def steer(argv):
    '''
    Steering routine for creating plots of the spectra.
    '''

    frac=0.1
    ymin=0.
    ymax=0.
    med_delta=3e-15
    itype='scale'
    mode=''
    filenames=[]

    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-frac':
            i+=1
            frac=eval(argv[i])
        elif argv[i]=='-min':
            i+=1
            ymin=eval(argv[i])
        elif argv[i]=='-max':
            i+=1
            ymax=eval(argv[i])
        elif argv[i]=='-med':
            itype='med'
        elif argv[i]=='-delta':
            i+=1
            med_delta=eval(argv[i])
        elif argv[i]=='-mode':
            i+=1
            mode=argv[i]
        elif argv[i][0]=='-':
            print('Error: Unknown switch: ',argv[i])
            return
        else:
            filenames.append(argv[i])
        i+=1

    if not filenames:
        print('Error: No input files specified')
        return

    if itype!='med' and ymax>0.0:
        itype='fixed'

    os.makedirs('Overview_Plot',exist_ok=True)

    if mode=='sep_back':
        # Two positional args: source file and separate background file
        if len(filenames)<2:
            print('Error: -mode sep_back requires a source file and a background file')
            return
        filename,backname=filenames[0],filenames[1]
        try:
            xtab=ascii.read(filename)
        except:
            print('Error: Could not read %s' % filename)
            return
        try:
            btab=ascii.read(backname)
        except:
            print('Error: Could not read %s' % backname)
            return
        xtab['FLUX']-=btab['FLUX']
        outname=os.path.basename(filename).replace('.txt','').replace('.tab','')
        do_all(xtab,ptype=itype,ymin=ymin,ymax=ymax,frac=frac,med_delta=med_delta,title=os.path.basename(filename))
        plt.savefig('Overview_Plot/%s.overview.png' % outname)
    else:
        # Default: each file is plotted independently
        for filename in filenames:
            try:
                xtab=ascii.read(filename)
            except:
                print('Error: Could not read %s, skipping' % filename)
                continue
            outname=os.path.basename(filename).replace('.txt','').replace('.tab','')
            do_all(xtab,ptype=itype,ymin=ymin,ymax=ymax,frac=frac,med_delta=med_delta,title=os.path.basename(filename))
            plt.savefig('Overview_Plot/%s.overview.png' % outname)




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
