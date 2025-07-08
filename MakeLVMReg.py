#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Read an rss file and create a region file containing all the fibers, and
read in a separate region file with one or more sources defined as 
regions with different colors.  Write out a new region file thaqt
contains all of the fiber positions, and use the second region file
to set the colors of the fibers.


Command line usage (if any):

    usage: MakeLVMReg.py -h *.fits whatever.reg or whatever.reg.txt

Description:  

    

Primary routines:

    doit

Notes:
                                       
History:

241106 ksl Coding begun

'''



# from lvm_ksl import GetSpec
import os
import sys
from astropy.io import fits, ascii
from astropy.table import Table
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import Angle

def check_positions_in_rectangle(table, center_ra, center_dec, width, height, theta):
    """
    Check which positions in an astropy table fall within a rotated rectangle on the sky.
    
    Parameters
    ----------
    table : astropy.table.Table
        Table containing positions. Must have columns 'RA' and 'Dec' in degrees
    center_ra : float
        Right ascension of rectangle center in degrees
    center_dec : float
        Declination of rectangle center in degrees
    width : float
        Width of rectangle in arcseconds
    height : float
        Height of rectangle in arcseconds
    theta : float
        Position angle in degrees from North (positive toward East)
        
    Returns
    -------
    table : astropy.table.Table
        Input table with new boolean column 'in_rectangle' added
    """
    # Convert everything to radians for calculations
    try:
        ra_rad = np.radians(table['RA'])
        dec_rad = np.radians(table['Dec'])
    except:
        ra_rad = np.radians(table['ra'])
        dec_rad = np.radians(table['dec'])
    

    center_ra_rad = np.radians(center_ra)
    center_dec_rad = np.radians(center_dec)
    theta_rad = np.radians(theta)
    
    # Convert width/height to radians
    width_rad = np.radians(width/3600.0)  # convert arcsec to degrees then radians
    height_rad = np.radians(height/3600.0)
    
    # Calculate position differences, accounting for cos(dec) in RA
    delta_ra = (ra_rad - center_ra_rad) * np.cos(center_dec_rad)
    delta_dec = dec_rad - center_dec_rad
    
    # Rotate coordinates to rectangle frame (positive theta rotates toward East)
    x = delta_ra * np.cos(-theta_rad) + delta_dec * np.sin(-theta_rad)
    y = -delta_ra * np.sin(-theta_rad) + delta_dec * np.cos(-theta_rad)
    
    # Check if points are within rectangle bounds
    in_rectangle = (np.abs(x) <= width_rad/2) & (np.abs(y) <= height_rad/2)
    
    # Add result to table
    table['in_area'] = in_rectangle
    
    return table

def check_positions_in_ellipse(table, center_ra, center_dec, semi_major, semi_minor, theta):
    """
    Check which positions fall within a rotated ellipse on the sky.
    
    Parameters
    ----------
    table : astropy.table.Table
        Table containing positions. Must have columns 'RA' and 'Dec' in degrees
    center_ra : float
        Right ascension of ellipse center in degrees
    center_dec : float
        Declination of ellipse center in degrees
    semi_major : float
        Semi-major axis in arcseconds
    semi_minor : float
        Semi-minor axis in arcseconds
    theta : float
        Position angle in degrees from North (positive toward East)
        
    Returns
    -------
    table : astropy.table.Table
        Input table with new boolean column 'in_ellipse' added
    """
    # Convert to radians
    # Convert everything to radians for calculations
    try:
        ra_rad = np.radians(table['RA'])
        dec_rad = np.radians(table['Dec'])
    except:
        ra_rad = np.radians(table['ra'])
        dec_rad = np.radians(table['dec'])
    
    center_ra_rad = np.radians(center_ra)
    center_dec_rad = np.radians(center_dec)
    theta_rad = np.radians(theta)
    
    # Convert axes to radians
    a_rad = np.radians(semi_major/3600.0)
    b_rad = np.radians(semi_minor/3600.0)
    
    # Calculate position differences
    delta_ra = (ra_rad - center_ra_rad) * np.cos(center_dec_rad)
    delta_dec = dec_rad - center_dec_rad
    
    # Rotate coordinates (positive theta rotates toward East)
    x = delta_ra * np.cos(-theta_rad) + delta_dec * np.sin(-theta_rad)
    y = -delta_ra * np.sin(-theta_rad) + delta_dec * np.cos(-theta_rad)
    
    # Check if points are within ellipse using normalized equation
    in_ellipse = (x/a_rad)**2 + (y/b_rad)**2 <= 1
    
    # Add result to table
    table['in_area'] = in_ellipse
    
    return table

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


def add_circular_region(xtab,ra,dec,radius,color):

    xtab['row_no']=np.arange(len(xtab))
    
    fib_no=get_closest(fiber_pos=xtab.copy(),xra=ra,xdec=dec)
    xfib=fib_no[fib_no['Sep']<=radius]
    # xfib is a table with only those fibers that satisfy the conditions
    # print(ra,dec)
    # print(xfib['row_no','Sep','ra','dec','fiberid'])
    if len(xfib)==0:
        xfib=fib_no[0]
    if xfib['Sep'][0]>70:
        print('Error: This RA and Dec (%.5f %.5f) does not appear to be in the field' % (ra,dec))
        return []

    for one in xfib:
        xtab['color'][one['row_no']]=color
    return xtab

    


def do_circular(filename,ra,dec,radius,outname='',color='red',target_type='science'):
    '''
    This version allows one to give a circular region a diffrent color from the other
    fibers
    '''
    xtab=get_good_fibers(filename,color='yellow',target_type='science')
    xtab= add_circular_region(xtab,ra,dec,radius,color)


    if outname=='':
        word=filename.split('/')
        outname=word[-1].replace('.fits','.reg')
    write_reg(outname,xtab,color)
                    
def get_good_fibers(filename,color='yellow',target_type='science'):    
    '''
    Just get the good fibers in the field and return the slitmap table
    along with exposure number
    '''
    try:
        x=fits.open(filename)
    except:
        print('Error: Could not find :',filename)

    xtab=Table(x['SLITMAP'].data)
    xtab=xtab[xtab['fibstatus']==0]
    xtab=xtab[xtab['targettype']==target_type]
    # print(np.unique(xtab['targettype']))
    # print(np.unique(xtab['fibstatus'],return_counts=True))
    xtab['color']=color
    exposure=x[0].header['EXPOSURE']
    return xtab,exposure

def do_simple(filename,outroot='',color='yellow',target_type='science'):
    '''
    This version just makes a region file for all of the good science fibers
    '''
    xtab,exposure=get_good_fibers(filename,color='yellow',target_type='science')
    if outroot=='':
        outname='SlitMap_%05d' % exposure
    else:
        outname='%s_%05d' % (outroot,exposure)
    write_reg(outname+'.reg',xtab,color)
    return xtab,exposure
                    
    


    
def do_one(filename,qtab,outroot=''):
    '''
    This routine creates a potential complex region file for extracting spectra. The fibers
    to extract ar current hardwired to be red, while those to be ignored are in colored yellow
    '''

    icolor=['red','green','blue','cyan','magenta','black','white']
    xtab,exposure=get_good_fibers(filename,color='yellow',target_type='science')

    # Create some sort of root for the file names
    word=filename.split('/')
    root=word[-1].replace('.fits','')
    if outroot!='':
        root='%s_%s' % (root,outroot)
    
    for one_row in qtab:
        print(one_row)
        if one_row['RegType']=='box':
            ftab=check_positions_in_rectangle(xtab, center_ra=one_row['RA'], center_dec=one_row['Dec'], width=one_row['Major'], height=one_row['Minor'], theta=one_row['Theta'])
        elif one_row['RegType']=='ellipse':
            ftab=check_positions_in_ellipse(xtab, center_ra=one_row['RA'], center_dec=one_row['Dec'], semi_major=one_row['Major'], semi_minor=one_row['Minor'], theta=one_row['Theta'])
        elif one_row['RegType']=='circle':
            ftab=check_positions_in_ellipse(xtab, center_ra=one_row['RA'], center_dec=one_row['Dec'], semi_major=one_row['Major'], semi_minor=one_row['Major'], theta=0.0)
        elif one_row['RegType']=='annulus':
            # the copy statement below is necessary because without this, qtab and ftab are pointed to the same exact table.
            ftab=check_positions_in_ellipse(xtab, center_ra=one_row['RA'], center_dec=one_row['Dec'], semi_major=one_row['Major'], semi_minor=one_row['Major'], theta=0.0)
            qtab=ftab.copy()
            ftab=check_positions_in_ellipse(xtab, center_ra=one_row['RA'], center_dec=one_row['Dec'], semi_major=one_row['Minor'], semi_minor=one_row['Minor'], theta=0.0)
            ftab['in_area']=np.select([qtab['in_area']==True],[False],default=ftab['in_area'])
        outfile='%s.%s.reg' % (root,one_row['Source_name'])
        print('Writing outfile: ',outfile)
        ftab['color']='yellow'
        ftab['color'][ftab['in_area']==True]='red'
        write_reg(outfile,ftab,color='yellow')


    return outfile



def steer(argv):
    '''
    MakeReg.py [-h] [-root whatever]  file1 file2   ra dec size ra dec size  
    or

    MakeReg.py -h -root whatever -circle ra dec size -circle ra dec size  file1  file 2

    or 

    MakeReg.py  n103b.txt file 1 file 2 
    '''
    filename=[]
    regfile=''
    masterfile=''
    outroot=''
    icolor=['red','green','blue','cyan','magenta','black','white']

    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-root':
            i+=1
            outroot=argv[i]
        elif argv[i][0]=='-':
            print('Error: cannot parse command line :',argv)
            return
        elif argv[i].count('fits'):
            filename.append(argv[i])
        elif regfile=='' and argv[i].count('.reg'):
            regfile=argv[i]
        elif regfile=='' and masterfile=='':
            masterfile=argv[i]
        i+=1

    # We want to work with a masterfile


    if regfile!='' and masterfile=='':
        print('Making masterfile')
        os.system('reg2master.py %s tmp_reg.txt' % regfile)
        masterfile='tmp_reg.txt'
    if masterfile=='':
        print('There is confusion: No masterfile input or created ', argv)
        return
    try: 
        xtab=ascii.read(masterfile)
    except:
        print('Could not read the masterfile :',masterfile)
        return





    for one_file in filename:
        do_one(one_file,xtab,outroot='')

    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
