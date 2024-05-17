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
240517 ksl Update to reflect fact that ra and decs are now 
    part of standard processing

'''


from astropy.io import fits,ascii
from astropy.wcs import WCS
from astropy.table import Table
import numpy as np
from astropy import wcs
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt



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




def get_closest(xra=81,xdec=-68,filename='unknown.fits'):

    try:
        x=fits.open(filename)
    except:
        print('Error: could not open %s' % filename)
        return []   

    fiber_tab=Table(x['SLITMAP'].data)

    xdistance=[]
    pos1=SkyCoord(ra=xra * u.deg, dec=xdec * u.deg, frame='icrs')


    xsep=[]
    i=0
    for one in fiber_tab:
        pos2= SkyCoord(ra=one['ra']* u.deg, dec=one['dec'] * u.deg, frame='icrs')
        angular_distance = pos1.separation(pos2).arcsec
        xsep.append(angular_distance)
        i+=1

    fiber_tab['Sep']=xsep
    fiber_tab.sort(['Sep'])

    x.close()
    return fiber_tab


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

    ra=None
    dec=None
    filename=[]

    i=1
    while i<len(argv):
        if argv[i]=='-h' or len(argv)<4:
            print(__doc__)
            return 
        elif ra==None:
            ra=argv[i]
        elif dec==None:
            dec=argv[i]
        elif argv[i].count('fit')>0:
            filename.append(argv[i])

        i+=1
        

    ra,dec=radec2deg(ra,dec)

    print(ra,dec)

    for one_file in filename:

        fib_no=get_closest(ra,dec,one_file)
        print('test')
        if len(fib_no):
            print('Plotting fiber %d' % fib_no['fiberid'][0])

            root=one_file.split('/')[-1]
            root=root.replace('.fits','')



            do_plot_all(filename=one_file,fiber_id=fib_no['fiberid'][0])
            plt.savefig('spec_%s_%.2f_%.2f.png' % (root,ra,dec) )

    


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
