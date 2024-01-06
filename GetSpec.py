#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Get a spectrum from a fiber or fibers in an RSS  file given an 
RA and DEC and a size 



Command line usage (if any):

    usage: PlotClosest.py [-h] exp_no ra dec

Description:  

    where exp_no is the exposure
    and ra and dec are the desired right ascension, written
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


def radec2deg(ra='05:13:06.2',dec='-10:13:14.2'):
    '''

    Convert an ra dec string to degrees.  The string can already
    be in degrees in which case all that happens is a conversion to
    a float

    If what is transferred is a float, the routine assumes it has been
    given ra and dec in degrees and just returns ra,dec

    170914  ksl Fix error associated with small negative declinations

    '''

    # print 'Before',ra,dec
    try:
        r=ra.split(':')
        d=dec.split(':')
    except AttributeError:
        return ra,dec


    # print 'After',ra,dec

    rr=float(r[0])
    if len(r)>1:
        rr=rr+float(r[1])/60.
    if len(r)>2:
        rr=rr+float(r[2])/3600.
    if len(r)>1:
        rr=15.*rr  # Since we assume ra was in hms

    sign=d[0].count('-')
    dd=abs(float(d[0]))
    x=0
    if len(d)>1:
        x=x+float(d[1])/60.
    if len(d)>2:
        x=x+float(d[2])/3600.

    dd=dd+x
    if sign:
        dd= -dd


    return rr,dd


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


def get_spec(filename,xfib,nfib=1):

    x=fits.open(filename)
    wave=x['WAVE'].data
    xxfib=xfib[0:nfib]
    flux=x['FLUX'].data[xxfib['fiberid']-1]
    print(flux.shape)
    xflux=np.average(flux,axis=0)
    print(xflux.shape)
    xspec=Table([wave,xflux],names=['WAVE','FLUX'])
    return xspec


def steer(argv):
    '''
    This is just a steering routine
    '''

    exposure=-99
    ra=None
    dec=None
    size=0
    root='Spec'

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return 
        elif argv[i]=='-size':
            i+=1
            size=eval(argv[i])
        elif argv[i]=='-root':
            i+=1
            root=argv[i]
        elif ra==None and argv[i][0]=='-':
            print('Error: Incorrect Command line: ',argv)
            return
        elif exposure<0:
            exposure=eval(argv[i])
        elif ra==None:
            ra=argv[i]
        elif dec==None:
            dec=argv[i]

        i+=1


        

    ra,dec=radec2deg(ra,dec)

    print(ra,dec)


    xcal_file='XCFrame-%08d.fits' % exposure
    xcal=fib2radec.locate_file(xcal_file)
    if xcal!=None:
        print('Using %s ' % xcal_file)
        x=fits.open(xcal)
        xtab=Table(x['SLITMAP'].data)
        xtab=xtab[xtab['RA']>0.0]
    else:
        print('Could not locate %s' % xcal_file)
    

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
    xfib=fib_no[fib_no['Sep']<=size]
    nfib=len(xfib)
    if nfib==0:
        nfib=1

    print('test ',nfib, size)
    print(xfib)

    xspec=get_spec(filename=xcal,xfib=fib_no,nfib=nfib)

    outname='%s_%05d_%.2f_%.2f_%02d.txt' % (root,exposure,ra,dec,nfib)

    xspec.write(outname,format='ascii.fixed_width_two_line',overwrite=True)












    


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
