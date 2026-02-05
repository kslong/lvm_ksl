#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

This is intended to create a standardized
plot for looking at spectra of SNRs and SNR
candidates extracted from LVM Data


Command line usage (if any):

    usage: SNRPlot.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

241111 ksl Coding begun

'''

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from astropy.io import ascii



def limit_spectrum(xtab,w_cen=3728,vlim=1000,v_offset=280):
    ztab=xtab.copy()
    wmin=w_cen*(1-vlim/3e5)
    wmax=w_cen*(1+vlim/3e5)
    delta=w_cen*v_offset/3e5
    # print(wmin,wmax,delta)
    ztab['WAVE']-=delta
    ztab=ztab[ztab['WAVE']>wmin]
    ztab=ztab[ztab['WAVE']<wmax]
    return ztab


def xplot(filename='Spec_09444_test_ave.txt',outroot='',vlim=1200,v_offset=280):
    xtab=ascii.read(filename)
    plt.close(1)
    fig=plt.figure(1,(8,8))
    gs= GridSpec(3, 3, figure=fig)
    ax1 = fig.add_subplot(gs[0, :])
    ax1.set_title(filename)
    ax1.plot(xtab['WAVE'],xtab['FLUX'],label='Spectrum',zorder=2)
    ax1.plot(xtab['WAVE'],xtab['ERROR'],label='Error',zorder=2)
    ax1.set_xlim(3600,9700)
    ymed=np.nanmedian(xtab['FLUX'])
    ymax=np.nanmax(xtab['FLUX'])
    # min=np.min(xtab['FLUX'])
    ax1.set_ylim(-0.1*ymax,1.1*ymax)
    ax1.legend()
    C=3e5

    ax2 = fig.add_subplot(gs[1, 0])
    ztab=limit_spectrum(xtab,3728,vlim,v_offset)
    ax2.plot(ztab['WAVE'],ztab['FLUX'],label='[OII]')
    ax2.plot(ztab['WAVE'],ztab['ERROR'])
    ax2.set_xlim(ztab['WAVE'][0],ztab['WAVE'][-1])
    ymed=np.median(ztab['FLUX'])
    ymax=np.max(ztab['FLUX'])
    ymin=np.min(ztab['FLUX'])
    if ymin<0 and ymed>0:
        ax2.set_ylim(-0.1*ymax,1.1*ymax)
    ax2.legend()


    ax2 = fig.add_subplot(gs[1, 1])
    ztab=limit_spectrum(xtab,5007,vlim,v_offset)
    ax2.plot(ztab['WAVE'],ztab['FLUX'],label='[OIII]')
    ax2.plot(ztab['WAVE'],ztab['ERROR'])
    ax2.set_xlim(ztab['WAVE'][0],ztab['WAVE'][-1])
    ymed=np.median(ztab['FLUX'])
    ymax=np.max(ztab['FLUX'])
    ymin=np.min(ztab['FLUX'])
    if ymin<0 and ymed>0:
        ax2.set_ylim(-0.1*ymax,1.1*ymax)
    ax2.legend()


    ax2 = fig.add_subplot(gs[1, 2])

    ztab=limit_spectrum(xtab,6300,vlim,v_offset)
    ax2.plot(ztab['WAVE'],ztab['FLUX'],label='[OI]')
    ax2.plot(ztab['WAVE'],ztab['ERROR'])
    ax2.set_xlim(ztab['WAVE'][0],ztab['WAVE'][-1])
    ymed=np.median(ztab['FLUX'])
    ymax=np.max(ztab['FLUX'])
    ymin=np.min(ztab['FLUX'])
    if ymin<0 and ymed>0:
        ax2.set_ylim(-0.1*ymax,1.1*ymax)
    ax2.legend()

    ax2 = fig.add_subplot(gs[2, 0])
    ztab=limit_spectrum(xtab,6563,vlim,v_offset)
    ax2.plot(ztab['WAVE'],ztab['FLUX'],label=r'H$\alpha$')
    ax2.plot(ztab['WAVE'],ztab['ERROR'])
    ax2.set_xlim(ztab['WAVE'][0],ztab['WAVE'][-1])
    ymed=np.median(ztab['FLUX'])
    ymax=np.max(ztab['FLUX'])
    ymin=np.min(ztab['FLUX'])
    if ymin<0 and ymed>0:
        ax2.set_ylim(-0.1*ymax,1.1*ymax)
    ax2.legend()

    ax2 = fig.add_subplot(gs[2, 1])
    ztab=limit_spectrum(xtab,6720,vlim,v_offset)
    ax2.plot(ztab['WAVE'],ztab['FLUX'],label=r'[SII]')
    ax2.plot(ztab['WAVE'],ztab['ERROR'])
    ax2.set_xlim(ztab['WAVE'][0],ztab['WAVE'][-1])
    ymed=np.median(ztab['FLUX'])
    ymax=np.max(ztab['FLUX'])
    ymin=np.min(ztab['FLUX'])
    if ymin<0 and ymed>0:
        ax2.set_ylim(-0.1*ymax,1.1*ymax)
    ax2.legend()

    ax2 = fig.add_subplot(gs[2, 2])
    ztab=limit_spectrum(xtab,9530.6,vlim,v_offset)
    ax2.plot(ztab['WAVE'],ztab['FLUX'],label=r'[SIII]')
    ax2.plot(ztab['WAVE'],ztab['ERROR'])
    ax2.set_xlim(ztab['WAVE'][0],ztab['WAVE'][-1])
    ymed=np.median(ztab['FLUX'])
    ymax=np.max(ztab['FLUX'])
    ymin=np.min(ztab['FLUX'])
    if ymin<0 and ymed>0:
        ax2.set_ylim(-0.1*ymax,1.1*ymax)
    ax2.legend()
    plt.tight_layout()

    if outroot.count('/'):
        xdir=os.path.split(outroot)[0]
        os.makedirs(xdir,exist_ok=True)

    if outroot=='':
        outroot=filename.replace('.txt','')
    plt.savefig(outroot+'.png')
    

def steer(argv):
    '''
    Mkae a plot of spectra extracted from a 
    region of regions of LVM data, with
    an emphasis on lines useful for looking
    at SNR data
    '''
    files=[]

    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i][0]=='-':
            print('Unknown switch : ',argv)
            return
        else:
            files.append(argv[i])

        i+=1


    for one in files:
        xplot(filename=one,outroot='',vlim=1200,v_offset=280)



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
