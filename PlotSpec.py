#!/usr/bin/env python
# coding: utf-8


'''
                    Space Telescope Science Institute

Synopsis:  

Create a standard plot of a spectrum that has been
extracted from the LVM data. The spectrum should
have been sky subtracted.  


Command line usage (if any):

    usage: PlotSpec.py [-h] [-frac 0.1] -min whatever - max whateverspectrum [backspec]

    where frac controls the autoscaling for the 
        maximum value in each panel and
        backspec is an optional spectrum to 
            be subtracted from the intial spectrum
        -min fixes the lower limit in all panels
        -max fixes the upper limit in all panels
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
    good=np.where(mask)[0]
    plt.plot(xx['WAVE'],xx['FLUX'])
    plt.xlim(wmin-extra,wmax+extra)
    if frac<1.0:
        ymin,ymax=plt.ylim()
        ymin=np.median(xx['FLUX'])-0.05*ymax
        plt.ylim(ymin,frac*ymax)
    return


def do_one_region_fixed(spectab,wmin=3600,wmax=4100,ymin=0, ymax=1e-13):
    extra=10
    xx=spectab[spectab['WAVE']>wmin-extra]
    xx=xx[xx['WAVE']<wmax+extra]
    mask=np.isfinite(xx['FLUX'])
    xx=xx[mask]
    good=np.where(mask)[0]
    plt.plot(xx['WAVE'],xx['FLUX'])
    plt.xlim(wmin-extra,wmax+extra)
    plt.ylim(ymin,ymax)
    return

def xmark(line='[SIII]',w=9069.3,frac=0.8):

    ymin,ymax=plt.ylim()
    # print(w,ymin,ymax)
    plt.text(w,frac *ymax,line,ha='center',color='red',clip_on=True)

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



def do_all(xtab,ptype='scale',ymin=0.0,ymax=1e-14,frac=0.1):
    '''
    Create the figure
    '''
    plt.figure(2,(8,12))
    plt.clf()
    
    wmin=3600
    wmax=9500
    wmax=9559
    delta=750
    nmax=int((wmax-wmin)/delta)+1
    # print(nmax)
    i=0
    while i<nmax:
        plt.subplot(nmax,1,i+1)
        wwmin=wmin+i*delta
        wwmax=wwmin+delta
        if ptype=='scale':
            do_one_region(xtab,wwmin,wwmax,frac)
        elif ptype=='fixed':
            do_one_region_fixed(xtab,wwmin,wwmax,ymin, ymax)
        else:
            print('Error: Indecipheragle type of plot: ',ptype)
            return

        do_lines()
        i+=1
    plt.tight_layout()


def steer(argv):
    '''
    Steering routine for creating plots of
    the spectra
    '''


    frac=0.1
    ymin=0.
    ymax=0.
    itype='scale'
    filename=''
    backname=''

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
        elif argv[i][0]=='-':
            print('Error: Unknown switch: ',argv)
            return
        elif filename=='':
            filename=argv[i]
        elif backname=='':
            backname=argv[i]
        else:
            print('Error: Badly formed command line args :',argv)
        i+=1

    try:
        xtab=ascii.read(filename)
    except:
        print('Error: Could not read %s' % filename)
        return

    if backname!='':

        try:
            btab=ascii.read(backname)
        except:
            print('Error: Could not read file to subract %s ' % backname)
            return
        xtab['FLUX']-=btab['FLUX']


    
    if ymax>0.0:
        itype='fixed'



    do_all(xtab,ptype=itype,ymin=ymin,ymax=ymax,frac=frac)
    words=filename.split('/')
    outname=words[-1].replace('.txt','')
    outname=outname.replace('.tab','')
    plt.savefig('%s.overview.png' % outname)




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
