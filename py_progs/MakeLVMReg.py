#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Read an rss file and create a region file containing all the fibers, and
read in a separate region file with one or more sources defined as
regions with different colors. Write out a new region file that
contains all of the fiber positions, and use the second region file
to set the colors of the fibers.


Command line usage (if any):

    usage: MakeLVMReg.py [-h] filename.fits masterfile

    or

    usage: MakeLVMReg.py [-h] filename.fits whatever.reg.txt

    or (batch mode, one region file per source)

    usage: MakeLVMReg.py [-h] [-all] [-snap snap_dir] [-ctype ave|med]
                         [-regdir reg_dir] [-buffer arcsec] [-imgdir img_dir]
                         masterfile

Description:

At a base level, this creates a region file for one or more rss files that
contains a region for each good science fiber.

Additionally, from the command line, one must provide either a master file
which defines a source or a source and background region for each file, or
alternatively one must provide a ds9 region file which contains information
that can be converted into a masterfile.

Most of the time one will input a masterfile, but if one enters a regionfile
a masterfile will be created on the fly.

The routine looks at the fiber positions and assigns colors to the fiber
positions based on the information in the masterfile.

The routine then writes out this to a new region file.

In batch mode (-all) the routine iterates over every unique Source_name in
the masterfile, locates the corresponding snapshot FITS file in snap_dir,
calls do_complex to assign fiber colors, and optionally calls
LSnap.make_one_image to overlay the fiber assignments on broadband images
found in img_dir.

Primary routines:

    do_complex  - assign fiber colors for one snapshot and one source set
    do_all      - batch-process all sources in a masterfile

Notes:

The routine creates one region file for each rss file and source in the input
master or region file. So if you give it a masterfile containing 20 lines,
but 10 different sources it will create 20 region files for each rss file.

Currently the routine understands region types: box, circle, ellipse, and
circular annulus.

The main reason for allowing multiple rss files from the command line is to
deal with tiling, where different fibers may be returned because of small
offsets.

This routine uses reg2master.py (which is not currently part of lvm_ksl).

History:

241106 ksl Coding begun
260301 ksl Added do_all for batch processing of source catalogs from snapshots

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
    
    Parameters:
        table (astropy.table.Table): Table containing positions with RA and Dec columns in degrees.
        center_ra (float): Right ascension of rectangle center in degrees.
        center_dec (float): Declination of rectangle center in degrees.
        width (float): Width of rectangle in arcseconds.
        height (float): Height of rectangle in arcseconds.
        theta (float): Position angle in degrees from North (positive toward East).

    Returns:
        astropy.table.Table: Input table with new boolean column in_rectangle added.
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
    
    Parameters:
        table (astropy.table.Table): Table containing positions with RA and Dec columns in degrees.
        center_ra (float): Right ascension of ellipse center in degrees.
        center_dec (float): Declination of ellipse center in degrees.
        semi_major (float): Semi-major axis in arcseconds.
        semi_minor (float): Semi-minor axis in arcseconds.
        theta (float): Position angle in degrees from North (positive toward East).

    Returns:
        astropy.table.Table: Input table with new boolean column in_ellipse added.
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
    Sort a slit table with RA and Dec values in order
    of the separation from a given RA and Dec.
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

def get_fibers_in_region(fiber_tab, region_tab, buffer=17.5):
    '''
    Add a column to fiber_tab to indicate which fibers overlap the regions
    in region_tab.

    Each region boundary is expanded by `buffer` arcseconds before testing
    whether fiber centers fall inside it.  With LVM fiber geometry (35 arcsec
    center-to-center spacing, 17.5 arcsec radius), the default buffer of 17.5
    arcsec captures every fiber with any overlap with the nominal region: a
    fiber center exactly at the expanded boundary has its edge just touching
    the original boundary.

    For annular (background) regions both the outer and inner radii are
    expanded by buffer.  This ensures that the expanded source and expanded
    background share a common boundary with no overlap and no gap, as long as
    the background inner radius equals the source outer radius.

    Parameters:
        fiber_tab (astropy.table.Table): Table of fiber positions.
        region_tab (astropy.table.Table): Table of region definitions; may
            contain multiple rows to build a composite selection.
        buffer (float): Arcseconds to add to each region boundary before
            testing containment.  Default 17.5 (one LVM fiber radius).
    '''

    fiber_tab['in_area']=False

    for one_row in region_tab:

        # Use a fresh copy for every region row so that check_positions_*
        # (which modifies its input table in place) never contaminates results
        # from a previous call.
        tmp=fiber_tab.copy()

        if one_row['RegType']=='box':
            # Major/Minor are full widths; expand each side by buffer
            width  = one_row['Major'] + 2 * buffer
            height = one_row['Minor'] + 2 * buffer
            check_positions_in_rectangle(tmp, center_ra=one_row['RA'], center_dec=one_row['Dec'], width=width, height=height, theta=one_row['Theta'])
            new_in_area=np.array(tmp['in_area'])

        elif one_row['RegType']=='ellipse':
            # Major/Minor are semi-axes
            semi_major = one_row['Major'] + buffer
            semi_minor = one_row['Minor'] + buffer
            check_positions_in_ellipse(tmp, center_ra=one_row['RA'], center_dec=one_row['Dec'], semi_major=semi_major, semi_minor=semi_minor, theta=one_row['Theta'])
            new_in_area=np.array(tmp['in_area'])

        elif one_row['RegType']=='circle':
            # Major is the radius
            radius = one_row['Major'] + buffer
            check_positions_in_ellipse(tmp, center_ra=one_row['RA'], center_dec=one_row['Dec'], semi_major=radius, semi_minor=radius, theta=0.0)
            new_in_area=np.array(tmp['in_area'])

        elif one_row['RegType']=='annulus':
            # Major is outer radius, Minor is inner radius.  Expand both by
            # buffer so the background inner boundary aligns with the expanded
            # source outer boundary, guaranteeing no overlap between source
            # and background fiber selections.
            # Use separate copies for outer and inner to avoid in-place
            # overwriting between the two check_positions_in_ellipse calls.
            outer = one_row['Major'] + buffer
            inner = one_row['Minor'] + buffer
            tmp_inner=fiber_tab.copy()
            check_positions_in_ellipse(tmp,       center_ra=one_row['RA'], center_dec=one_row['Dec'], semi_major=outer, semi_minor=outer, theta=0.0)
            check_positions_in_ellipse(tmp_inner, center_ra=one_row['RA'], center_dec=one_row['Dec'], semi_major=inner, semi_minor=inner, theta=0.0)
            new_in_area=np.array(tmp['in_area']) & ~np.array(tmp_inner['in_area'])

        else:
            print('Error: Unknown region type :', one_row['RegType'])
            raise ValueError

        fiber_tab['in_area']=fiber_tab['in_area'] | new_in_area

    return fiber_tab




def do_complex(filename,qtab,outroot='',buffer=17.5,reg_dir=''):
    '''
    Process one rss file producing a region file for each source in
    the masterfile.  There can be multiple lines in the masterfile
    that correspond to a single source, which allows one to create
    a complex set of fibers. Returns the name of the region file
    that was written.

    Parameters:
        filename (str): The name of an rss file.
        qtab (astropy.table.Table): An already opened masterfile.
        outroot (str): Optional root name for output files.
        buffer (float): Arcseconds to expand each region boundary before
            testing fiber containment.  Default 17.5 (one LVM fiber radius),
            which captures every fiber with any overlap with the nominal region.
        reg_dir (str): Directory in which to write the region file.
            Created if it does not exist.  Default is '' (current directory).

    Returns:
        str: The name of the region file that was written.
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

    if reg_dir:
        os.makedirs(reg_dir,exist_ok=True)

    sources=np.unique(qtab['Source_name'])
    #  print('XXX sources:',sources)
    for one_source in sources:
        basename='%s.%s.reg' % (root,one_source)
        outfile=os.path.join(reg_dir,basename) if reg_dir else basename
        one_object_tab=qtab[qtab['Source_name']==one_source]
        if 'SourceBack' in one_object_tab.colnames:
            source_tab=one_object_tab[one_object_tab['SourceBack']=='Source']
            back_tab=one_object_tab[one_object_tab['SourceBack']=='Back']
            ftab['color']='yellow'
            if len(back_tab)>0:
                xback_tab=get_fibers_in_region(xtab,back_tab,buffer=buffer)
                ftab['color'][xback_tab['in_area']==True]='green'
                print('After back assignment:', np.unique(ftab['color'],return_counts=True))
            if len(source_tab)>0:
                xsource_tab=get_fibers_in_region(xtab,source_tab,buffer=buffer)
                ftab['color'][xsource_tab['in_area']==True]='red'
                print('After source assignment:', np.unique(ftab['color'],return_counts=True))
        else:
           xsource_tab=get_fibers_in_region(xtab,one_object_tab,buffer=buffer)
           ftab['color'][xsource_tab['in_area']==True]='red'

        write_reg(outfile,ftab,color='yellow')
        # print('XXX Knox ',outfile)
        return outfile





    
def do_one(filename,qtab,outroot='',buffer=17.5):
    '''
    This routine creates a potential complex region file for extracting spectra. The fibers
    to extract ar current hardwired to be red, while those to be ignored are in colored yellow

    The routine reads an rss file, and initially returns a list of
    all of the good science fibers.  It then produces a region file
    for each row in qtab, where the name is given by the source name
    in qtab, and outroot.  Each region boundary is expanded by buffer
    arcseconds before testing containment, so that edge fibers with any
    overlap are captured.
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

        ftab=get_fibers_in_region(xtab,one_row,buffer=buffer)
        outfile='%s.%s.reg' % (root,one_row['Source_name'][0])
        print('Writing outfile: ',outfile)
        ftab['color']='yellow'
        ftab['color'][ftab['in_area']==True]='red'
    write_reg(outfile,ftab,color='yellow')
    return outfile



def do_all(masterfile, snap_dir='Snap', c_type='ave', reg_dir='FiberReg', buffer=17.5, img_dir=''):
    '''
    Batch-process all sources in a masterfile, creating one region file per
    source from the corresponding snapshot FITS file.

    For every unique Source_name in masterfile the routine:

    1. Locates the snapshot file snap_dir/<source_name>.<c_type>.fits
       (falls back to the .gz variant if the plain file is absent).
    2. Calls do_complex to assign fiber colors (red=source, green=background
       annulus) using the region geometry defined in masterfile.
    3. If img_dir is provided, calls LSnap.make_one_image for every FITS file
       in img_dir whose name begins with the source name, overlaying the fiber
       color assignments on each broadband image.  Output images go to zimage/.

    Parameters:
        masterfile (str): Source catalog with SourceBack, Source_name, RA, Dec,
            RegType, Major, Minor, and Theta columns, such as that produced by
            GenAnnularBackground.py.
        snap_dir (str): Directory containing snapshot FITS files (default: 'Snap').
        c_type (str): Combination-type suffix used in snapshot filenames,
            e.g. 'ave' or 'med' (default: 'ave').
        reg_dir (str): Directory in which region files are written.  Created
            automatically if it does not exist (default: 'FiberReg').
        buffer (float): Arcseconds to expand each region boundary before testing
            fiber containment (default: 17.5, one LVM fiber radius).
        img_dir (str): Directory containing per-source broadband image cutouts.
            If set, LSnap.make_one_image is called for each matching file after
            region-file creation.  Leave empty to skip (default: '').
    '''
    import glob

    if img_dir:
        from lvm_ksl import LSnap

    xtab = ascii.read(masterfile)
    sources = np.unique(xtab['Source_name'])

    n = len(sources)
    for i, source_name in enumerate(sources):
        print('\n!! Processing %s (%d of %d)' % (source_name, i + 1, n))

        filename = '%s/%s.%s.fits' % (snap_dir, source_name, c_type)
        if not os.path.isfile(filename):
            if os.path.isfile(filename + '.gz'):
                filename = filename + '.gz'
            else:
                print('Skipping %s: snapshot file not found' % filename)
                continue

        source_tab = xtab[xtab['Source_name'] == source_name]
        reg_file = do_complex(filename=filename, qtab=source_tab,
                              outroot='source', buffer=buffer, reg_dir=reg_dir)
        if reg_file is None:
            print('Skipping %s: could not create region file' % source_name)
            continue

        if img_dir:
            img_files = glob.glob('%s/%s*.fits*' % (img_dir, source_name))
            if img_files:
                for img_file in sorted(img_files):
                    LSnap.make_one_image(img_file, reg_file, None, None)
            else:
                print('No broadband images found for %s in %s' % (source_name, img_dir))


def steer(argv):
    '''
    MakeReg.py [-h] [-root whatever]  file1 file2   ra dec size ra dec size
    or

    MakeReg.py -h -root whatever -circle ra dec size -circle ra dec size  file1  file 2

    or

    MakeReg.py  n103b.txt file 1 file 2

    or (batch mode)

    MakeReg.py -all [-snap snap_dir] [-ctype ave|med] [-regdir reg_dir]
               [-buffer arcsec] [-imgdir img_dir] masterfile
    '''
    filename=[]
    regfile=''
    masterfile=''
    outroot=''
    xall=False
    snap_dir='Snap'
    c_type='ave'
    reg_dir='FiberReg'
    buffer=17.5
    img_dir=''
    icolor=['red','green','blue','cyan','magenta','black','white']

    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-root':
            i+=1
            outroot=argv[i]
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-snap':
            i+=1
            snap_dir=argv[i]
        elif argv[i]=='-ctype':
            i+=1
            c_type=argv[i]
        elif argv[i]=='-regdir':
            i+=1
            reg_dir=argv[i]
        elif argv[i]=='-buffer':
            i+=1
            buffer=float(argv[i])
        elif argv[i]=='-imgdir':
            i+=1
            img_dir=argv[i]
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

    if xall:
        if masterfile=='':
            print('Error: -all requires a masterfile argument')
            return
        do_all(masterfile, snap_dir=snap_dir, c_type=c_type, reg_dir=reg_dir,
               buffer=buffer, img_dir=img_dir)
        return

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
