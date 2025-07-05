#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Make an rss_combined file given a set of images, a position and
a size.  The purpose of this rooutine is to be able to make snap_shot 
sized RSS files for examining a number of sources at one time.


Command line usage (if any):

    usage: rss_combine_pos.py  -c_type med -size 20  -out test ra dec filenames

Description:  

Primary routines:

    doit

Notes:
 
The goal here is to develop a procedure that will allow me to create an rss_file that just contains data centered around a particular position, nominally a SNR.  My plan is to be able to execute this on a large number of postions spread over the LMC or SMC

Here are position for n49 and n49b

ellipse(81.341994,-65.996487,79.03",77.48",90.0) # text={SNR_N49B}
ellipse(81.501085,-66.082249,42.35",40.00",110.0) # text={SNR_N49}


                                       
History:

250702 ksl Coding begun

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt
import os



ra_n49b=81.341994
dec_n49b=-65.996487

ra_n49=81.501085
dec_n49=-66.082249








import numpy as np
from astropy.io import ascii,fits
import psutil

from lvm_ksl import rss_combine








def do_fixed(filenames,ra, dec, pa, size,c_type='ave',outroot=''):
    wcs=rss_combine.create_wcs(ra_deg=ra, dec_deg=dec,pos_ang=pa, size_deg=size)
    new_slitmap_table=rss_combine.generate_grid(wcs,35.)
    print('new: ',len(new_slitmap_table))
    xmin=np.min(new_slitmap_table['X'])
    xmax=np.max(new_slitmap_table['X'])
    ymin=np.min(new_slitmap_table['Y'])
    ymax=np.max(new_slitmap_table['Y'])
    print(xmin,xmax,ymin,ymax)
    slit=rss_combine.prep_tables_square(wcs, filenames)
    # print(slit)
    # this list has many more fibers than we wnat
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
        if c_type=='sum':
            one_tab=rss_combine.get_nearest(new_slitmap_table,one_tab)
        else:
            one_tab=rss_combine.frac_calc2(new_slitmap_table,one_tab)
        # print('toasty\n',one_tab)
        zslit.append(one_tab)

    # zslit is a list of tables, now with added columns indicating how to approtion the data from indvidual files into the final image

    print('The number of files being combined:',len(slit))

    foo=fits.open(filenames[0])
    final = fits.HDUList()
    final.append(foo['PRIMARY'])
    final.append(foo['FLUX'])
    final.append(foo['IVAR'])
    final.append(foo['MASK'])
    final.append(foo['WAVE'])
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

    final=rss_combine.add_extension_with_same_shape(final, 'FLUX','EXPOSURE')

    zero_array=np.zeros(shape,dtype=np.int64)
    final['MASK'].data=zero_array



    eflux=np.zeros_like(final['FLUX'].data)
    new_slitmap_table['EXPOSURE']=0.0

    # now we want to crate various arrays that indicate how each input file would be split into the new filters system
    VERY_BIG = 1e50
    VERY_SMALL = 1e-50

    print('Begin accumulating data')

    i=0
    while i < len(filenames):
        print('starting to rmap  %s' % (filenames[i]),flush=True)
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
    final['FLUX'].data=xfl
    final['EXPOSURE'].data=xexp
    final['IVAR'].data=xvar


    print('Finished creating output arrays.\n')
    value=rss_combine.check_memory()
    mem=psutil.virtual_memory()
    print(f'Total Memory: {mem.total/1e9:.2f} GB')
    print(f' Used Memory: {mem.used/1e9:.2f} GB')
    print('This is %.2f pecent of the available memory' % value)


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

   # Before writing everythin out, we need to add the size to the hdr

    nxmax=np.max(new_slitmap_table['X'])
    nymax=np.max(new_slitmap_table['Y'])
    final['WCS_INFO'].header['XNAXIS1']=int(nxmax)
    final['WCS_INFO'].header['XNAXIS2']=int(nymax)

    if outroot=='':
        outroot='test_square'

    print ('\n Final Stats for output image : %s.fits' % outroot)
    print('Combined  FLUX  %10.3e %10.3e %10.e %10.3e'  % (np.nanmedian(final['FLUX'].data),np.nanstd(final['FLUX'].data),np.nanmin(final['FLUX'].data),np.nanmax(final['FLUX'].data)))
    zmax=np.nanmax(final['IVAR'].data)
    print('Combined  IVAR  %10.3e %10.3e %10.e %10.3e\n'  % (np.nanmedian(final['IVAR'].data),zmax*np.nanstd(final['IVAR'].data/zmax),np.nanmin(final['IVAR'].data),np.nanmax(final['IVAR'].data)))


    
    new_slitmap_table.write('%s.tab' % outroot,format='ascii.fixed_width_two_line',overwrite=True)
    final.writeto(outroot+'.fits',overwrite=True)
    print('Wrote new RSS file: %s.fits' % (outroot))

    return outroot






def steer(argv):
    '''
    rss_combine_pos.py  -c_type med -size 20  -out test ra dec filenames
    '''
    ra=-99.
    dec=-99.
    size=20.
    outroot='test'
    c_type='med'
    filenames=[]

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
        elif argv[i]=='-c_type':
            i+=1
            c_type=argv[i]
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



    do_fixed(filenames,ra, dec, pa=0, size=size/60.,c_type=c_type,outroot=outroot)




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
