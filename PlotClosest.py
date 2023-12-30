#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Make a plot of a fiber in a RRS file given either a 
fiber number of and RA and dec


Command line usage (if any):

    usage: PlotClosest.py [-h] exp_no ra dec

Description:  

    where exp_no is the exposure
    and ra and dec are the desrired right ascension, written
        either in degrees or in h:m:s d:m:s

Primary routines:

    doit

Notes:
                                       
History:

231212 ksl Coding begun

'''


from astropy.io import fits,ascii
from astropy.wcs import WCS
from astropy.table import Table
import numpy as np
from astropy import wcs
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt

import fib2radec
from ksl.astro import radec2deg




n103b_ra=77.251042
n103b_dec=-68.72674



def get_closest(fiber_pos,xra=n103b_ra,xdec=n103b_dec):
    xdistance=[]
    pos1=SkyCoord(ra=xra * u.deg, dec=xdec * u.deg, frame='icrs')

    
    xsep=[]
    i=0
    for one in fiber_pos:
        pos2= SkyCoord(ra=one['RA']* u.deg, dec=one['Dec'] * u.deg, frame='icrs')
        angular_distance = pos1.separation(pos2).arcsec
        xsep.append(angular_distance)
        i+=1
        
    fiber_pos['Sep']=xsep
    fiber_pos.sort(['Sep'])
    return fiber_pos


def limit_wavelengths(xspec,wmin,wmax):
    xspec=xspec[xspec['WAVE']>wmin]
    xspec=xspec[xspec['WAVE']<wmax]
    xmed=np.median(xspec['FLUX'])
    xmax=np.max(xspec['FLUX'])
    xdelt=xmax-xmed
    xmin=xmed-0.05*xdelt  
    xmax=xmax+0.05*xdelt
    return xmin,xmax,xspec

def do_plot_all(filename='reduced/lvmCFrame-00007824.fits',fiber_id=770):
    x=fits.open(filename)
    wave=x['WAVE'].data
    flux=x['FLUX'].data[fiber_id-1]
    xtab=Table([wave,flux],names=['WAVE','FLUX'])

    
    plt.figure(1,(8,6))
    plt.clf()
    plt.subplot(2,2,1)
    wmin=3650
    wmax=3820
    ymin,ymax,xspec=limit_wavelengths(xtab,wmin,wmax)
    plt.plot(xspec['WAVE'],xspec['FLUX'])
    plt.xlim(wmin,wmax)
    plt.ylim(ymin,ymax)
    plt.subplot(2,2,2)
    wmin=4825
    wmax=5050
    ymin,ymax,xspec=limit_wavelengths(xtab,wmin,wmax)    
    plt.plot(xspec['WAVE'],xspec['FLUX'])
    plt.xlim(wmin,wmax)
    plt.ylim(ymin,ymax)   
    plt.subplot(2,2,3)
    wmin=6200
    wmax=6490
    ymin,ymax,xspec=limit_wavelengths(xtab,wmin,wmax)
    plt.plot(xspec['WAVE'],xspec['FLUX'])
    plt.xlim(wmin,wmax)
    plt.ylim(ymin,ymax)
    plt.subplot(2,2,4)
    wmin=6500
    wmax=6800
    ymin,ymax,xspec=limit_wavelengths(xtab,wmin,wmax)    
    plt.plot(xspec['WAVE'],xspec['FLUX'])
    plt.xlim(wmin,wmax)
    plt.ylim(ymin,ymax)   
    return


def steer(argv):
    '''
    This is just a steering routine
    '''

    exposure=-99
    ra=None
    dec=None

    i=1
    while i<len(argv):
        if argv[i]=='-h' or len(argv)<4:
            print(__doc__)
            return 
        elif exposure<0:
            exposure=eval(argv[i])
        elif ra==None:
            ra=argv[i]
        elif dec==None:
            dec=argv[i]

        i+=1
        

    ra,dec=radec2deg(argv[2],argv[3])

    print(ra,dec)

    calibrated_file='lvmCFrame-%08d.fits' % exposure
    guider_file='lvm.sci.coadd_s%08d.fits' % exposure

    xcal=fib2radec.locate_file(calibrated_file)
    if xcal==None:
        print('Could not locate calibrated file: %s' % calibrated_file)
        return
    
    xguide=fib2radec.locate_file(guider_file)
    if xguide==None:
        print('Could not locate guider file: %s' % guider_file)
        xguide=''

    xguide,racen,decen,pa,xtab=fib2radec.get_ra_dec(xcal,xguide,outname='xfib.txt')


    fib_no=get_closest(fiber_pos=xtab,xra=ra,xdec=dec)

    print('Plotting fiber %d' % fib_no['fiberid'][0])



    do_plot_all(filename=xcal,fiber_id=fib_no['fiberid'][0])
    plt.savefig('spec%05d_%.2f_%.2f.png' % (exposure,ra,dec) )

    


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
