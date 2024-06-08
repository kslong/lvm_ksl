#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Get a spectrum from a fiber or fibers in an RSSfile given an 
RA and DEC and a size 



Command line usage (if any):

    usage: Getspec.py  [-h] [-size xx]  filename ra dec



    where 
    -h prints out this help
    -size xx defines a separation in arc seconsds of fibers to
        incluce in the output spectrum.
    exp_no is the exposure
    and ra and dec are the desired right ascension, written
        either in degrees or in h:m:s d:m:s

Description:  

    The routine looks for a file where the RA's and Dec's have
    been assigned to fibers (by xcal.py) and finds the
    fiber or fibers that are closest to this position 
    and extract the spectral information from for this.

    If multiple fibers are selected the resulting spectrum
    the routine generally returns the average flux in
    the fibers, rather than the sum.  The routine
    also return an average or median of various other
    quantities.

    The errors are the errors for one fiber, not something
    that has been reduced by the number of fibers.
Primary routines:

    steer - unlike most of ksl's routines, at present
    the steering routine manages the entire process.

Notes:

    The routine also produces a regiong file that shows
    what fibers were used.
                                       
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
    '''
    Sort a slit table with RA's and DEC's in order
    of the separattion from a given RA and Dec
    '''

    xdistance=[]
    pos1=SkyCoord(ra=xra * u.deg, dec=xdec * u.deg, frame='icrs')

    
    xsep=[]
    i=0
    for one in fiber_pos:
        pos2= SkyCoord(ra=one['ra']* u.deg, dec=one['dec'] * u.deg, frame='icrs')
        angular_distance = pos1.separation(pos2).arcsec
        xsep.append(angular_distance)
        i+=1
        
    fiber_pos['Sep']=xsep
    fiber_pos.sort(['Sep'])
    return fiber_pos


def get_spec(filename,xfib,nfib=1):
    '''
    Retrieve the spectra from the first
    nfib fibers in the xfib table
    '''

    x=fits.open(filename)
    wave=x['WAVE'].data
    xxfib=xfib[0:nfib]
    flux=x['FLUX'].data[xxfib['fiberid']-1]
    ivar=x['IVAR'].data[xxfib['fiberid']-1]
    mask=x['MASK'].data[xxfib['fiberid']-1]
    sky=x['SKY'].data[xxfib['fiberid']-1]
    sky_ivar=x['SKY_IVAR'].data[xxfib['fiberid']-1]
    lsf=x['LSF'].data[xxfib['fiberid']-1]
    # print(flux.shape)
    xflux=np.nanmean(flux,axis=0)
    xsky=np.nanmean(sky,axis=0)
    xmask=np.sum(mask,axis=0)
    
    xlsf=np.nanmedian(lsf,axis=0)

    # print('z',xerr.shape)
    ivar=np.nansum(ivar,axis=0)
    # print(xerr.shape)
    xerr=1/np.sqrt(ivar)
    # print(xflux.shape)
    # print(xsky.shape)
    # print(xmask.shape)
    xsky_error=np.nansum(sky_ivar,axis=0)
    xsky_error=1/np.sqrt(xsky_error)
    xspec=Table([wave,xflux,xerr,xsky,xsky_error,xmask,xlsf],names=['WAVE','FLUX','ERROR','SKY','SKY_ERROR','MASK','LSF'])
    xspec['WAVE'].format='.1f'
    xspec['FLUX'].format='.3e'
    xspec['ERROR'].format='.3e'
    xspec['SKY'].format='.3e'
    xspec['SKY_ERROR'].format='.3e'
    xspec['LSF'].format='.2f'
    return xspec

header='''
# Region file format: DS9 version 4.1
# Filename: lmc_snr.txt.reg
global color=yellow width=3 font="helvetica 14 bold" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
'''
def write_reg(filename,xtab):
    g=open(filename,'w')
    g.write(header)
    for one in xtab:
        g.write('circle(%f,%f,5.")  # text={%d}\n' % (one['ra'],one['dec'],one['fiberid']))
    g.close()



def steer(argv):
    '''
    This routine parses the command line

    and then oversees the rest of the proces


    In future it would be better to split these
    acitvites apart so the routine can be called
    from a Jupyter script, but that is not
    the way things are at present
    '''

    filename=None
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
        elif filename==None:
            filename=argv[i]
        elif ra==None:
            ra=argv[i]
        elif dec==None:
            dec=argv[i]
        else:
            print('Error: Incorrect Command line: (exta args) ',argv)
            return
        i+=1


        

    # Having parsed the command line, do the work

    ra,dec=radec2deg(ra,dec)

    try:
        x=fits.open(filename)
    except:
        print('Error: could not open %s',filename)
        return



    xtab=Table(x['SLITMAP'].data)
    fib_no=get_closest(fiber_pos=xtab,xra=ra,xdec=dec)
    xfib=fib_no[fib_no['Sep']<=size]

    nfib=len(xfib)
    if nfib==0:
        nfib=1


    if nfib==1 and fib_no['Sep'][0]>70:
        print('Error: This RA and Dec (%.5f %.5f) does not appear to be in the field' % (ra,dec))
        return

    print('Taking spectra from\n',fib_no['fiberid'][0:nfib])
    xspec=get_spec(filename=filename,xfib=fib_no,nfib=nfib)

    exposure=x['PRIMARY'].header['EXPOSURE']

    if root=='Spec':
        outname='%s_%05d_%.2f_%.2f_%02d.txt' % (root,exposure,ra,dec,nfib)
    else:
        outname='%s_%05d_%02d.txt' % (root,exposure,nfib)

    regname=outname.replace('.txt','.reg')
    write_reg(regname,xtab=xfib)

    xspec.meta['comments']=['Filename %s' % filename,'RA %.5f' % ra, 'Dec %.5f' % dec]
    print(xspec)

    xspec.write(outname,format='ascii.fixed_width_two_line',overwrite=True)

    print('The output file is %s' % (outname))












    


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
