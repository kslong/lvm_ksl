#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Get a spectrum from a fiber or fibers in an RSSfile given an 
RA and DEC and a size 



Command line usage (if any):

    usage: Getspec.py  [-h] [-root whatever] -median filename ra dec [rmin] [rmax]



    where 
    -h prints out this help
    -median constructs the median instead of the average spectrum
    -root whatever prepends a root to the standard file name
    exp_no is the exposure
    ra and dec are the desired right ascension, written
        either in degrees or in h:m:s d:m:s
    [rmin] is and optional inner radius, and 
    [rmax] is an optinal outer radius

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

    if rmin is not present then the spectrum of the closest
    fiber will be retrieved

    if only rmin is present then all fibers with centers
    inside that radius will be retried

    if rmin and rmax are rpersend the average or medain apctra
    will be regriived

Primary routines:

    steer - unlike most of ksl's routines, at present
    the steering routine manages the entire process.

Notes:

    The routine also produces a regiong file that shows
    what fibers were used.
                                       
History:

231212 ksl Coding begun

'''


import os
from astropy.io import fits,ascii
from astropy.wcs import WCS
from astropy.table import Table
import numpy as np
from astropy import wcs
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt

import fib2radec

from LocateReduced import read_drpall,find_em
XTOP='/uufs/chpc.utah.edu/common/home/sdss51/'
XXTOP='%s/Projects/lvm_data/sas/' % (os.environ['HOME'])


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


def get_spec(filename,xfib,nfib=1,xtype='ave'):
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

    if xtype=='ave':
        xflux=np.nanmean(flux,axis=0)
        xsky=np.nanmean(sky,axis=0)
        xmask=np.sum(mask,axis=0)
        xlsf=np.nanmean(lsf,axis=0)
    elif xtype=='med':
        xflux=np.nanmedian(flux,axis=0)
        xsky=np.nanmedian(sky,axis=0)
        xmask=np.sum(mask,axis=0)
        xlsf=np.nanmedian(lsf,axis=0)
    else:
        print('Error: getspec: only ave or med allowd for xtype')
        return

    # print('z',xerr.shape)
    ivar=np.nansum(ivar,axis=0)
    # print(xerr.shape)
    xerr=np.select([ivar>1],[1/np.sqrt(ivar)],default=np.nan)
    # print(xflux.shape)
    # print(xsky.shape)
    # print(xmask.shape)
    xsky_error=np.nansum(sky_ivar,axis=0)
    xsky_error=np.select([xsky_error>1],[1/np.sqrt(xsky_error)],default=np.nan)
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
global color=%s width=4 font="helvetica 14 bold" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
'''
def write_reg(filename,xtab,color='yellow'):
    g=open(filename,'w')
    g.write(header % (color) )
    for one in xtab:
        g.write('circle(%f,%f,5.")  # text={%d}\n' % (one['ra'],one['dec'],one['fiberid']))
    g.close()


def get_circle(xtab,ra,dec,radius):
    '''
    Get the fibers which lie within a certain distance of
    a center position

    where xtab is slitmap extension,  
    ra, dec are a position in degrees
    and
    radius is a maxium separation in arcsec
    '''

    fib_no=get_closest(fiber_pos=xtab,xra=ra,xdec=dec)
    xfib=fib_no[fib_no['Sep']<=radius]
    if len(xfib)==0:
        xfib=fib_no[0]
    if xfib['Sep'][0]>70:
        print('Error: This RA and Dec (%.5f %.5f) does not appear to be in the field' % (ra,dec))
        return []
    return xfib


def get_annulus(xtab,ra,dec,rmin,rmax):
    '''
    Get the fibers which lie within an annulus centered on a specific ra and dec

    where xtab is slitmap extension,  
    ra, dec are a position in degrees
    and
    rmin and rmax are the minimum and maximum sizes of the annulus
    '''
    fib_no=get_closest(fiber_pos=xtab,xra=ra,xdec=dec)
    xfib=fib_no[fib_no['Sep']>=rmin]
    xfib=xfib[xfib['Sep']<=rmax]

    if len(xfib)==0:
        print('Error: get_annulus: no fibers RA and DEC ( (%.5f %.5f) and rmin and rmax (%.2f %.2f)' % (ra,dec,rmin,rmax))
    
    return xfib

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
    size_max=0
    root='Spec'
    xtype='ave'

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
        elif argv[i][0:4]=='-med':
            xtype='med'
        elif ra==None and argv[i][0]=='-':
            print('Error: Incorrect Command line: ',argv)
            return
        elif filename==None:
            filename=argv[i]
        elif ra==None:
            ra=argv[i]
        elif dec==None:
            dec=argv[i]
        elif size==0:
            size=eval(argv[i])
        elif size_max==0:
            size_max=eval(argv[i])
        else:
            print('Error: Incorrect Command line: (exta args) ',argv)
            return
        i+=1


        

    # Having parsed the command line, do the work

    ra,dec=radec2deg(ra,dec)

    if filename.count('fits')==0:
        exp_no=int(filename)
        drpall=read_drpall()
        xlocate=find_em(drpall,exp_no,exp_no)
        try:
            filename=xlocate['location'][0]
            print('Found location for %d on mjd %d' % (exp_no,xlocate['mjd'][0]))
        except:
            print('Error: Could not find file for exposure %d in drpall.fits' % exp_no)
            return

    if os.path.isfile(filename)==False:
        xfilename='%s/%s' % (XTOP,filename)
        if os.path.isfile(xfilename)==False:
            xfilename='%s/%s' % (XXTOP,filename)
    else:
        xfilename=filename

    try:
        x=fits.open(xfilename)
    except:
        print('Error: could not open %s' % xfilename)
        return



    xtab=Table(x['SLITMAP'].data)
    if size_max==0:
        fibers=get_circle(xtab,ra,dec,radius=size)
    else:
        print('annulus ',size,size_max)
        fibers=get_annulus(xtab,ra,dec,size,size_max)

    if len(fibers)==0:
        return


    print('Taking spectra from\n',fibers['fiberid'])
    xspec=get_spec(filename=xfilename,xfib=fibers,nfib=len(fibers),xtype=xtype)

    exposure=x['PRIMARY'].header['EXPOSURE']


    outname='%s_%05d_%.2f_%.2f' % (root,exposure,ra,dec)

    if size>0:
        outname='%s_%d' % (outname,size)
    if size_max>0:
        outname='%s_%d' % (outname,size_max)
    if xtype=='med':
        outname='%s_med' % (outname)
    else:
        outname='%s_ave' % (outname)

    if size_max>0:
        write_reg('%s.reg' % outname,xtab=fibers,color='green')
    else:
        write_reg('%s.reg' % outname,xtab=fibers,color='yellow')

    xspec.meta['comments']=['Filename %s' % filename,'RA %.5f' % ra, 'Dec %.5f' % dec, 'nfibers %d' % (len(fibers))]

    xspec.write('%s.txt'% outname,format='ascii.fixed_width_two_line',overwrite=True)

    print('The output file is %s.txt' % (outname))












    


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
