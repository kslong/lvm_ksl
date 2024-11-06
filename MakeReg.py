#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Create a region file containing a region for each fiber in and exposure. 
Alloow for the possibility of changing the color for at least one circulare
region


Command line usage (if any):

    usage: MakeReg.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

241106 ksl Coding begun

'''



# from lvm_ksl import GetSpec
import sys
from astropy.io import fits
from astropy.table import Table
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u




header='''
# Region file format: DS9 version 4.1
# Filename: %s
global color=%s width=4 font="helvetica 14 bold" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
'''
def write_reg(filename,xtab,color='yellow',size=5):
    g=open(filename,'w')
    g.write(header % (filename,color) )
    for one in xtab:
        if one['color']==color:
            g.write('circle(%f,%f,%.1f")  # text={%d}\n' % (one['ra'],one['dec'],size,one['fiberid']))
        else:
            g.write('circle(%f,%f,%.1f")  # color=%s text={%d}\n' % (one['ra'],one['dec'],size,one['color'],one['fiberid']))
        
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



def get_closest(fiber_pos,xra,xdec):
    '''
    Sort a slit table with RA's and DEC's in order
    of the separation from a given RA and Dec
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



def do_circular(filename,ra,dec,radius,outname='',color='yellow',target_type='science'):
    '''
    This version allows one to give a circular region a diffrent color from the other
    fibers
    '''
    x=fits.open(filename)
    xtab=Table(x['SLITMAP'].data)
    xtab.info()
    xtab=xtab[xtab['fibstatus']==0]
    xtab=xtab[xtab['targettype']==target_type]
    xtab['color']=color
    xtab['row_no']=np.arange(len(xtab))
    
    fib_no=get_closest(fiber_pos=xtab.copy(),xra=ra,xdec=dec)
    xfib=fib_no[fib_no['Sep']<=radius]
    # xfib is a table with only those fibers that satisfy the conditions
    print(ra,dec)
    print(xfib['row_no','Sep','ra','dec','fiberid'])
    if len(xfib)==0:
        xfib=fib_no[0]
    if xfib['Sep'][0]>70:
        print('Error: This RA and Dec (%.5f %.5f) does not appear to be in the field' % (ra,dec))
        return []

    for one in xfib:
        xtab['color'][one['row_no']]='red'

    if outname=='':
        word=filename.split('/')
        outname=word[-1].replace('.fits','.reg')
    write_reg(outname,xtab,color)
                    
    

def do_simple(filename,outname='',color='yellow',target_type='science'):
    '''
    This version just makes a region file for all of the good science fibers
    '''
    x=fits.open(filename)
    xtab=Table(x['SLITMAP'].data)
    xtab.info()
    xtab=xtab[xtab['fibstatus']==0]
    xtab=xtab[xtab['targettype']==target_type]
    print(np.unique(xtab['targettype']))
    print(np.unique(xtab['fibstatus'],return_counts=True))
    xtab['color']=color
    if outname=='':
        word=filename.split('/')
        outname=word[-1].replace('.fits','.reg')
    write_reg(outname,xtab,color)
                    
    



def steer(argv):
    '''
    MakeReg.py [-h] [-root whatever] 
    '''
    ra=None
    dec=None
    size=None
    filename=[]

    i=1
    while i<len(argv):
        if argv[i][0:2]=='h':
            print(__doc__return)
        elif argv[i][0]=='-' and ra==None:
            print('Error: cannot parse command line :',argv)
            return
        elif argv[i].count('fits'):
            filename.append(argv[i])
        elif ra==None:
            ra=argv[i]
        elif dec==None:
            dec=argv[i]
        elif size==None:
            size=eval(argv[i])
        i+=1

    if dec==None:
        for one_file in filename:
            do_simple(one_file)
    else:
        ra,dec=radec2deg(ra,dec)
        for one_file in filename:
            do_circular(one_file,ra,dec,size,outname='',color='yellow',target_type='science')





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
