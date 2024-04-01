#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Plot how well the standards that are observed in an LVM exposure
are calibrated.


Command line usage (if any):

    usage: eval_standard.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

240318 ksl Coding begun

'''



from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from lvmdrp.external import ancillary_func



from scipy.signal import boxcar
from scipy.signal import convolve

def xsmooth(flux,smooth=21):
    '''
    boxcar smooth the flux
    '''
    if (smooth)>1:
        q=convolve(flux,boxcar(smooth)/float(smooth),mode='same')
        return(q)
    else:
        return(flux)



def get_standard(xx,slit='P1-2'):
    '''
    get the spectrum of a standard stat in the lvmCframe
    data
    '''

    ztab=Table(xx['SLITMAP'].data)
    xstar=ztab[ztab['targettype']=='standard']
  
    xline=xstar[xstar['orig_ifulabel']==slit]
    # print(xline)
    wave=xx['WAVE'].data
    flux=xx['FLUX'].data[xline['fiberid']-1]
    wave=wave.flatten()
    flux=flux.flatten()
    return wave,flux



def compare_with_gaia(filename='lvmCFrame-00005059.fits',outroot=''):
    try:
        x=fits.open(filename)
    except:
        print('Error: could not open ',filename)
        return
    header=x[0].header
    exposure=header['EXPOSURE']
    mjd=header['MJD']
    i=1
    good=[]
    obsmag=[]
    gaia_id=[]
    fib=[]
    while i<13:
        try:
            key='STD%dBIN' % i
            obsmag.append(header[key])
            key='STD%dID' % i
            gaia_id.append(header[key])
            key='STD%dFIB' % i
            fib.append(header[key])
            good.append(i)
        except:
            pass
        i+=1
    print(good)
    print(obsmag)
    print(gaia_id)
    print(fib)

    plt.figure(1,(8,8))
    j=0
    while j<len(fib):
        print('start %d' % j)
        wave,flux=ancillary_func.retrive_gaia_star(gaia_id[j],'test')
        swave,sflux=get_standard(x,fib[j])
        plt.semilogy(swave,xsmooth(sflux),label=fib[j])
        plt.semilogy(wave,flux,'k')
        j+=1

    plt.xlim(3500,9500)
    plt.title('MJD %d Exposure %d' % (mjd,exposure))
    # plt.legend()
    print(plt.ylim())
    ylm=plt.ylim()
    plt.ylim(1e-13,ylm[1])
    plt.tight_layout()

    if outroot=='':
        word=filename.split('/')
        outroot=word[-1].replace('.fits','')

    outfile='standard_%s.png' % outroot
    plt.savefig(outfile)

    return

                
def steer(argv):

    files=[]

    i=1
    while i<len(argv):
        if argv[i].count('-h'):
            print(__doc__)
            return
        elif argv[i][0]=='-':
            print('Error: could not process command line: ',argv)
        else:
            files.append(argv[i])
        i+=1

    for one in files:
        compare_with_gaia(one)
                          
        



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
