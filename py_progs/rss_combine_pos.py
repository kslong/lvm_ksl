#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Create a combined RSS file from multiple SFrame images, centered on a
user-specified position with a user-specified size. This routine is
designed for making small "snapshot" RSS files for examining specific
sources such as supernova remnants.

Command line usage::

    rss_combine_pos.py [-sum] [-med] [-size arcmin] [-out name] [-keep] ra dec filenames

Arguments: ra and dec are the center position in degrees (required).
filenames are the input SFrame FITS files to combine.

Options:

-size arcmin
    Size of output region in arcminutes (default: 20).

-out name
    Output filename root (default: 'test').

-sum
    Use get_nearest for apportionment (no fractional splitting; each input
    fiber's full flux goes to the single nearest output fiber). By default,
    frac_calc2 is used (flux split based on overlap).

-med
    Use median for combining remapped images. By default, mean is used.
    Note that -med and -sum are independent; -med controls combination,
    -sum controls apportionment.

-keep
    Keep temporary files in xtmp/ directory.

Description:

This routine differs from rss_combine.py in that it allows the user
to specify the exact center position and size of the output region,
rather than automatically calculating these from the input files.
This is useful for extracting small regions around specific targets
from a larger set of observations.

The routine creates a WCS centered on the specified RA/Dec with the
given size, generates output fiber positions on a regular grid (35
arcsec spacing), filters input fibers to only those within the
specified region, apportions flux from input fibers to output fibers
using frac_calc2 (based on circular overlap area) or get_nearest
(with -sum flag), combines the remapped images using mean (default)
or median (-med), and writes the output to outroot.fits and outroot.tab.

The -sum and -med options work the same as in rss_combine.py, where
-sum controls apportionment method (frac_calc2 vs get_nearest) and
-med controls combination method (mean vs median).

Unlike rss_combine.py, this routine always uses a regular grid for
output fiber positions (equivalent to fib_type='xy'), so the -orig
option is not available.

Output FITS extensions are the same as rss_combine.py:
PRIMARY, FLUX, IVAR, MASK, WAVE, LSF, SLITMAP, WCS_INFO, EXPOSURE.

Primary routines:

    do_fixed

Notes:

    The goal is to create RSS files containing data centered around a
    particular position, nominally a SNR. This can be executed on a
    large number of positions spread over the LMC or SMC.

    Example positions for N49 and N49B:
    ellipse(81.341994,-65.996487,79.03",77.48",90.0) # text={SNR_N49B}
    ellipse(81.501085,-66.082249,42.35",40.00",110.0) # text={SNR_N49}


History:

250102 ksl Coding begun
250203 ksl Added separate fib_type parameter to align with rss_combine.py

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil



ra_n49b=81.341994
dec_n49b=-65.996487

ra_n49=81.501085
dec_n49=-66.082249








import numpy as np
from astropy.io import ascii,fits
import psutil

from lvm_ksl import rss_combine








def do_fixed(filenames, ra, dec, pa, size, fib_type='xy', c_type='ave', outroot='', keep_tmp=False):
    '''
    Create a combined RSS file centered on a fixed position.

    This is the main processing routine that creates a WCS at the specified
    position, generates output fiber positions on a regular grid, and combines
    flux from multiple input SFrame files.

    Unlike rss_combine.do_combine, this routine:
    * Uses a user-specified center position rather than calculating from files
    * Uses a user-specified size rather than calculating bounding box
    * Filters input fibers to only those within the specified region
    * Always uses a regular grid (no 'orig' option)

    The fib_type and c_type parameters work the same as in rss_combine.py,
    where fib_type controls apportionment method and c_type controls
    combination method.

    Parameters
    ----------
    filenames : list
        List of input SFrame FITS filenames to combine.
    ra : float
        Right Ascension of the center in degrees.
    dec : float
        Declination of the center in degrees.
    pa : float
        Position angle in degrees (rotation of the WCS).
    size : float
        Size of the output region in degrees.
    fib_type : str
        Method for apportioning flux from input to output fibers.
        Options are 'xy' (use frac_calc2, flux split based on overlap)
        or 'sum' (use get_nearest, no fractional splitting).
    c_type : str
        Method for combining the remapped images from multiple input files.
        Options are 'ave' for mean (default) or 'med' for median.
    outroot : str
        Root name for output files. Default creates 'test_square'.
    keep_tmp : bool
        If True, keep the temporary xtmp/ directory. Default is False.

    Returns
    -------
    str
        The output root name used for the files.
    '''
    wcs=rss_combine.create_wcs(ra_deg=ra, dec_deg=dec,pos_ang=pa, size_deg=size)
    new_slitmap_table=rss_combine.generate_grid(wcs,35.)
    # print('new: ',len(new_slitmap_table))
    xmin=np.min(new_slitmap_table['X'])
    xmax=np.max(new_slitmap_table['X'])
    ymin=np.min(new_slitmap_table['Y'])
    ymax=np.max(new_slitmap_table['Y'])
    # print(xmin,xmax,ymin,ymax)
    xtab=rss_combine.xcheck(filenames)
    qtab=xtab[xtab['Good']=='No']
    if len(qtab)>0:
        print('do_fixed: %d in a total of %d were deemed bad' % (len(qtab),len(xtab)))
        print(qtab)
        gtab=xtab[xtab['Good']=='Yes']
        filenames=gtab['Filename']
    slit=rss_combine.prep_tables_square(wcs, filenames)
    # print(slit)
    # this list has many more fibers than we want
    reduced=[]
    for one_tab in slit:
        one_tab=one_tab[one_tab['X']>xmin-1]
        one_tab=one_tab[one_tab['X']<xmax+1]
        one_tab=one_tab[one_tab['Y']>ymin-1]
        one_tab=one_tab[one_tab['Y']<ymax+1]    
        print(len(one_tab))
        reduced.append(one_tab)
    slit=reduced
    print('\nApportioning fractional contributions from individual to final virtual fiber positions')

    zslit=[]
    for one_tab in slit:
        if fib_type=='sum':
            one_tab=rss_combine.get_nearest(new_slitmap_table,one_tab)
        else:
            one_tab=rss_combine.frac_calc2(new_slitmap_table,one_tab)
        # print('toasty\n',one_tab)
        zslit.append(one_tab)

    # zslit is a list of tables, now with added columns indicating how to apportion the data from individual files into the final image

    print('The number of files being combined:',len(slit))

    foo=fits.open(filenames[0])
    final = fits.HDUList()
    final.append(foo['PRIMARY'])
    final.append(foo['FLUX'])
    final.append(foo['IVAR'])
    final.append(foo['MASK'])
    final.append(foo['WAVE'])
    final.append(foo['LSF'])
    final.append(foo['SLITMAP'])
    print('Adding WCS')
    # print(wcs)

    xheader = wcs.to_header()
    data = np.random.rand(6, 6)
    image_hdu = fits.ImageHDU(data=data, header=xheader, name="WCS_INFO")
    final.append(image_hdu)

    shape=[len(new_slitmap_table),12401]
    xzero_array = np.zeros(shape, dtype=np.float32)

    final['FLUX'].data=xzero_array.copy()
    final['IVAR'].data=xzero_array.copy()
    final['LSF'].data=xzero_array.copy()

    final=rss_combine.add_extension_with_same_shape(final, 'FLUX','EXPOSURE')

    zero_array=np.zeros(shape,dtype=np.int64)
    final['MASK'].data=zero_array



    eflux=np.zeros_like(final['FLUX'].data)
    new_slitmap_table['EXPOSURE']=0.0

    # now we want to create various arrays that indicate how each input file would be split into the new fibers system
    VERY_BIG = 1e50
    VERY_SMALL = 1e-50

    print('Begin accumulating data')

    i=0
    while i < len(filenames):
        print('Starting to remap %s' % (filenames[i]),flush=True)
        q=zslit[i]
        for one_row in q:
            new_slitmap_table['EXPOSURE'][one_row['fib_master']-1]+=one_row['frac']
        rss_combine.remap_one(filenames[i],q,new_slitmap_table,shape)
        print('Finished remapping %s\n' % (filenames[i]),flush=True)
        i+=1


    # Now read them all back
    value=rss_combine.check_memory()
    mem=psutil.virtual_memory()
    print(f'Total Memory: {mem.total/1e9:.2f} GB')
    print(f' Used Memory: {mem.used/1e9:.2f} GB')
    print('This is %.2f percent of the available memory' % value)

    print('\nNow read the data back, and construct the final rss spectra')

    xfiles=[]
    for one_file in filenames:
        xfile=one_file.split('/')[-1]
        xfile='xtmp/%s' % xfile
        xfiles.append(xfile)

    xfl=rss_combine.process_remapped_images(file_list=xfiles, extension='FLUX', xproc=c_type,memory_limit=1_000_000_000)
    xexp=rss_combine.process_remapped_images(file_list=xfiles, extension='EXPOSURE', xproc='sum',memory_limit=1_000_000_000)
    xvar=rss_combine.process_remapped_images(file_list=xfiles, extension='IVAR', xproc='sum',memory_limit=1_000_000_000)
    xlsf=rss_combine.process_remapped_images(file_list=xfiles, extension='LSF', xproc=c_type,memory_limit=1_000_000_000)
    final['FLUX'].data=xfl
    final['EXPOSURE'].data=xexp
    final['IVAR'].data=xvar
    final['LSF'].data=xlsf


    print('Finished creating output arrays.\n')
    value=rss_combine.check_memory()
    mem=psutil.virtual_memory()
    print(f'Total Memory: {mem.total/1e9:.2f} GB')
    print(f' Used Memory: {mem.used/1e9:.2f} GB')
    print('This is %.2f percent of the available memory' % value)


    new_slitmap_table['fibstatus']=0
    new_slitmap_table['targettype']='science'

    nbad=0
    for row in new_slitmap_table:
        xflux=final['FLUX'].data[row['fiberid']-1]
        if rss_combine.check_for_nan(xflux):
            row['fibstatus']=95
            nbad+=1
    ngood=len(new_slitmap_table)-nbad
    print('Of %d virtual fibers, %d were filled, and the %d were unfilled  ' % (len(new_slitmap_table),ngood,nbad))


    final['MASK'].data=np.select([final['EXPOSURE'].data==0],[1],default=0)
    # Explicitly set BITPIX to a positive integer value
    final['MASK'].header['BITPIX'] = 64
    # print(final['EXPOSURE'].data[5500])
    # print(final['MASK'].data[5500])

    hdu = fits.BinTableHDU(new_slitmap_table.as_array(), name='SLITMAP')
    final['SLITMAP']=hdu

    # Before writing everything out, we need to add the size to the header

    nxmax=np.max(new_slitmap_table['X'])
    nymax=np.max(new_slitmap_table['Y'])
    final['WCS_INFO'].header['XNAXIS1']=int(nxmax)
    final['WCS_INFO'].header['XNAXIS2']=int(nymax)

    if outroot=='':
        outroot='test_square'

    # Append combination type to outroot
    outroot = outroot + '.' + c_type

    # Add keywords to primary header indicating options used
    final['PRIMARY'].header['RSSCOMB'] = ('rss_combine_pos', 'RSS combination script used')
    final['PRIMARY'].header['FIBTYPE'] = (fib_type, 'Fiber apportionment type (xy/sum)')
    final['PRIMARY'].header['COMBTYPE'] = (c_type, 'Combination method (ave/med)')
    final['PRIMARY'].header['NFILES'] = (len(filenames), 'Number of input files combined')
    final['PRIMARY'].header['CENTRERA'] = (ra, 'Center RA in degrees')
    final['PRIMARY'].header['CENTERDE'] = (dec, 'Center Dec in degrees')
    final['PRIMARY'].header['SIZEDEG'] = (size, 'Region size in degrees')

    print ('\n Final Stats for output image : %s.fits' % outroot)
    print('Combined  FLUX  %10.3e %10.3e %10.e %10.3e'  % (np.nanmedian(final['FLUX'].data),np.nanstd(final['FLUX'].data),np.nanmin(final['FLUX'].data),np.nanmax(final['FLUX'].data)))
    zmax=np.nanmax(final['IVAR'].data)
    print('Combined  IVAR  %10.3e %10.3e %10.e %10.3e\n'  % (np.nanmedian(final['IVAR'].data),zmax*np.nanstd(final['IVAR'].data/zmax),np.nanmin(final['IVAR'].data),np.nanmax(final['IVAR'].data)))


    
    new_slitmap_table.write('%s.tab' % outroot,format='ascii.fixed_width_two_line',overwrite=True)
    final.writeto(outroot+'.fits',overwrite=True)
    print('Wrote new RSS file: %s.fits' % (outroot))

    if keep_tmp==False:
        shutil.rmtree('./xtmp')

    return outroot






def steer(argv):
    '''
    Parse command line arguments and run the RSS combining process.

    Usage:
        rss_combine_pos.py [-sum] [-med] [-size arcmin] [-out name]
                           [-keep] ra dec filenames

    Parameters:
        argv (list): Command line arguments (sys.argv).

    Returns:
        None: Calls do_fixed with parsed arguments.
    '''
    ra=-99.
    dec=-99.
    size=20.
    outroot='test'
    fib_type='xy'
    c_type='ave'
    filenames=[]
    keep_tmp=False

    i=1
    while i < len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][:4]=='-out':
            i+=1
            outroot=argv[i]
        elif argv[i]=='-size':
            i+=1
            size=eval(argv[i])
        elif argv[i]=='-keep':
            keep_tmp=True
        elif argv[i]=='-sum':
            fib_type='sum'
        elif argv[i][0:4]=='-med':
            c_type='med'
        elif ra<0 and argv[i][0]=='-':
            print('Error: Cannot parse command line :',argv)
            return
        elif ra<0:
            ra=eval(argv[i])
        elif dec<-90:
            dec=eval(argv[i])
        elif argv[i].count('fit'):
            filenames.append(argv[i])
        else:
            print('Error: Cannot parse command line :',argv)
            return
        i+=1

    # Validate outroot does not start with '-'
    if outroot.startswith('-'):
        print('Error: outroot cannot start with "-": %s' % outroot)
        print('This may indicate a missing argument after -out')
        return

    do_fixed(filenames, ra, dec, pa=0, size=size/60., fib_type=fib_type, c_type=c_type, outroot=outroot, keep_tmp=keep_tmp)




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
