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

    usage: MakeLVMReg.py [-h] *.fits masterfile  or whatever.reg.txt

Description:  
    At a base level, this creates a region file for one or more rss 
    files that contains a region for each good science fiber.

    Additionally, from the command line, one must provide either
    a master file which defines a source or a source and background
    region for each file, or altenatively one must provide a ds9 region
    file which contains information that can be converted into a 
    masterfile.

    Most of the time one will input a masterfile, but if one enters
    a regionfile a mastefile will be created on the fly.

    The routine looks at the fiber positions and assigns colors
    to the fiber positions based on the infomation in the masterfile

    The routine then writes out this to a new region file.

    

Primary routines:

    do_complex

Notes:
    The routine creates one region file for each rss file and source
    in the input master or region file.  So if you give it a masterfile
    containing 20 lines, but 10 different sources it will create
    20 region files for each rss file.

    Currently the routine understands region types, box, circle, ellipse
    and (circular) annulus.

    The main reason for allowing multiple rss files from the command line
    is to deal with tiling, where different fibers may be returned 
    because of small offsets

    This routine uses reg2master.py (which is not currently part of lvm_ksl)
                                       
History:

241106 ksl Coding begun

'''



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

    '''
    write reg opens a file and cretes a vailid
    region file with colors set to falues defined
    by xtab, so the color combination can be complex

    '''
    g=open(filename,'w')
    g.write(header % (filename,color) )
    for one in xtab:
        if one['color']==color:
            g.write('circle(%f,%f,%.1f")  # text={%d}\n' % (one['ra'],one['dec'],size,one['fiberid']))
        else:
            g.write('circle(%f,%f,%.1f")  # color=%s text={%d}\n' % (one['ra'],one['dec'],size,one['color'],one['fiberid']))
        
    g.close()



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

def get_fibers_in_region(fiber_tab,region_tab,size_min=17.5):
    '''
    Add a column to fiber_tab to indicate which fibers are in the
    in the reg_tab.  reg_tab can contain multiple rows
    '''

    xtab=fiber_tab.copy()
    fiber_tab['in_area']=False

    for one_row in region_tab:

        if one_row['Major']<size_min:
            print('Setting Major to size_min')
            one_row['Major']=size_min
        if one_row['RegType'!='circle'] and one_row['Minor']<size_min:
            print('Setting Minor to size_min')
            one_row['Minor']=size_min
        print(one_row)

        if one_row['RegType']=='box':
            ftab=check_positions_in_rectangle(xtab, center_ra=one_row['RA'], center_dec=one_row['Dec'], width=one_row['Major'], height=one_row['Minor'], theta=one_row['Theta'])
        elif one_row['RegType']=='ellipse':
            ftab=check_positions_in_ellipse(xtab, center_ra=one_row['RA'], center_dec=one_row['Dec'], semi_major=one_row['Major'], semi_minor=one_row['Minor'], theta=one_row['Theta'])
        elif one_row['RegType']=='circle':
            ftab=check_positions_in_ellipse(xtab, center_ra=one_row['RA'], center_dec=one_row['Dec'], semi_major=one_row['Major'], semi_minor=one_row['Major'], theta=0.0)
        elif one_row['RegType']=='annulus':
            # print('XXX - working on annular region')
            # Appparently annulus is circular, and in this case it's only the outher radius that needs to be greater than the minimu siae
            # the copy statement below is necessary because without this, qtab and ftab are pointed to the same exact table.
            ftab=check_positions_in_ellipse(xtab, center_ra=one_row['RA'], center_dec=one_row['Dec'], semi_major=one_row['Major'], semi_minor=one_row['Major'], theta=0.0)

            # print('outer: ',len(ftab[ftab['in_area']==True]))
            # ftab.write('foo1.txt',format='ascii.fixed_width_two_line',overwrite=True)

            qtab=ftab.copy()
            # qtab has all of the fibers in the outer circle marked as in_area
            ftab=check_positions_in_ellipse(xtab, center_ra=one_row['RA'], center_dec=one_row['Dec'], semi_major=one_row['Minor'], semi_minor=one_row['Minor'], theta=0.0)
            # ftab has all of the fibers in the inner circle marked as in_area

            # print('inner: ',len(ftab[ftab['in_area']==True]))
            # ftab.write('foo2.txt',format='ascii.fixed_width_two_line',overwrite=True)

            ftab['in_area']=np.select([ftab['in_area']==True],[False],default=qtab['in_area'])

            # print('final: ',len(ftab[ftab['in_area']==True]))
            # ftab.write('foo3.txt',format='ascii.fixed_width_two_line',overwrite=True)

            # finished

        else:
            print('Error: Uknown region type :', one_row['RegType'])
            raise ValueError
        fiber_tab['in_area']=fiber_tab['in_area'] | ftab['in_area']
    return fiber_tab




def do_complex(filename,qtab,outroot='',size_min=17.5):
    '''
    Process one rss file producing a region file for
    each source in the masterfile.  There can be
    mulitple lines in the masterfile that correspond
    to a single source, which allows one to create
    a complex set of fibers

    do_complex returns the name of the region file 
    that was written
    
    where 
        filename is the name of an rss file
        qtab is an already opened masterfile
        size is a minium size for any region parameter, 
            design so that one will find at least one fiber
            for each region
    '''
    # print('XXX = Starting ', filename)

    icolor=['red','green','blue','cyan','magenta','black','white']
    # xtab,exposure=get_good_fibers(filename,color='yellow',target_type='science')

    # Create some sort of root for the file names
    word=filename.split('/')
    root=word[-1].replace('.fits','')
    if outroot!='':
        root='%s.%s' % (root,outroot)
    xtab,exposure=get_good_fibers(filename,color='yellow',target_type='science')

    ftab=xtab.copy()

    # Create some sort of root for the file names
    word=filename.split('/')
    root=word[-1].replace('.fits','')
    if outroot!='':
        root='%s_%s' % (root,outroot)

    sources=np.unique(qtab['Source_name'])
    #  print('XXX sources:',sources)
    for one_source in sources:
        outfile='%s.%s.reg' % (root,one_source)
        one_object_tab=qtab[qtab['Source_name']==one_source]
        if 'SourceBack' in one_object_tab.colnames:
            source_tab=one_object_tab[one_object_tab['SourceBack']=='Source']
            back_tab=one_object_tab[one_object_tab['SourceBack']=='Back']
            ftab['color']='yellow'
            if len(back_tab)>0:
                xback_tab=get_fibers_in_region(xtab,back_tab,size_min=17.5)
                ftab['color'][xback_tab['in_area']==True]='green'
                # print(np.unique(ftab['color'],return_counts=True))
            if len(source_tab)>0:
                xsource_tab=get_fibers_in_region(xtab,source_tab,size_min=17.5)
                ftab['color'][xsource_tab['in_area']==True]='red'
                # print(np.unique(ftab['color'],return_counts=True))
        else:
           xsource_tab=get_fibers_in_region(xtab,one_object_tab,size_min=17.5)
           ftab['color'][xsource_tab['in_area']==True]='red'

        write_reg(outfile,ftab,color='yellow')
        # print('XXX Knox ',outfile)
        return outfile





    
def do_one(filename,qtab,outroot='',size_min=17.5):
    '''
    This routine creates a potential complex region file for extracting spectra. The fibers
    to extract ar current hardwired to be red, while those to be ignored are in colored yellow

    The routine reads an rss file, and initally returns a list of 
    all of the good science fibers.  It then produces a region file
    for each row in qtab, where the name is given by the source name
    in qtab, and outroot.  So that at least one fiber will be slected if
    the object is in the field a minium size for the region is set.
    '''

    print('OK Knox')
    icolor=['red','green','blue','cyan','magenta','black','white']
    xtab,exposure=get_good_fibers(filename,color='yellow',target_type='science')

    # Create some sort of root for the file names
    word=filename.split('/')
    root=word[-1].replace('.fits','')
    if outroot!='':
        root='%s.%s' % (root,outroot)
    
    for i in range(len(qtab)):
        one_row=qtab[i:i+1]

        ftab=get_fibers_in_region(xtab,one_row,size_min=17.5)
        outfile='%s.%s.reg' % (root,one_row['Source_name'][0])
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
        print('Making masterfile from a region file')
        os.system('reg2master.py %s tmp_reg.txt' % regfile)
        masterfile='tmp_reg.txt'

    if masterfile=='':
        print('There is confusion: No masterfile input or created, so quitting ', argv)
        return

    try: 
        xtab=ascii.read(masterfile)
    except:
        print('Could not read the masterfile :',masterfile)
        return

    # Process each rss file in turn.
    for one_file in filename:
        do_complex(one_file,xtab,outroot='')

    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
