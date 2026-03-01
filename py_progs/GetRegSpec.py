#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:

Get a spectrum from a fiber or fibers in an RSS file given an
RA and Dec based on a region file that contains one region per fiber.


Command line usage::

    GetRegSpec.py [-h] [-root whatever] [-med] [-sum] filename[s] source_reg [color] [back_reg] [back_color]

    GetRegSpec.py [-h] [-all] [-snap snap_dir] [-out out_dir] [-med] [-ctype ave|med] masterfile

Options:

-h
    Print this help.

-root whatever
    Prepend a root string to the output filename.

-med
    Construct the median instead of the average source spectrum.

-sum
    Construct the sum instead of the average source spectrum.

-all
    Batch mode: process all sources listed in masterfile. Calls
    MakeLVMReg.do_complex to assign fiber colors for each source,
    then extracts source (red) and background (green) spectra.

-snap snap_dir
    Directory containing snapshot FITS files (default: Snap).
    Used only with -all.

-out out_dir
    Output directory for extracted spectra (default: Snap_spec).
    Used only with -all.

-ctype ave|med
    Combination type used when naming snapshot files, i.e. the
    extension before .fits (default: ave). Used only with -all.

-regdir reg_dir
    Directory in which fiber region files are written by MakeLVMReg
    (default: FiberReg).  When reading region files, the current
    directory is checked first and FiberReg is used as a fallback.
    Used only with -all.

-imgdir img_dir
    Directory containing per-source broadband image cutouts (e.g.
    MCELS H-alpha and [SII] FITS files).  If provided, LSnap is called
    after each spectrum extraction to overlay the fiber color assignments
    on every FITS file whose name begins with the source name.  Output
    images are written to zimage/.  Used only with -all.

Arguments (single-file mode): filename[s] is one or more RSS spectra
files; source_reg is the name of a region file (must have extension
.reg); color defines the color of fibers to extract; back_reg is an
optional background region file; back_color is the corresponding
background color.  The positional arguments source_reg, color,
back_reg, and back_color must appear in order.

Arguments (batch mode): masterfile is a source catalog with SourceBack,
Source_name, RA, Dec, RegType, Major, Minor, and Theta columns, such
as that produced by GenAnnularBackground.py.

Description:

The routine reads a region file produced by MakeLVMReg.py in which
every science fiber has been assigned a color, and extracts the
average (or median) spectrum of fibers with the requested color.
If a background region is provided the median background spectrum
is subtracted.

If multiple fibers are selected the routine returns the average flux
across fibers rather than the sum.  Errors are per-fiber, not reduced
by the number of fibers.

Primary routines:

    do_one   - extract spectrum for a single RSS file
    do_all   - batch-process all sources in a masterfile

Notes:

The routine also writes a region file showing which fibers were used.

History:

231212 ksl Coding begun
260223 ksl Added do_all for batch processing of source catalogs
260223 ksl Added LSnap broadband image overlay support to do_all

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



def read_reg(xfile,color=None):
    '''
    This reads a region file that contains
    all a region for each of the fibers.
    It either selects those fibers of a certain
    color or if color is None, it assumes that
    one wants the fibers which have the fewest
    fibers of a given color.

    If xfile is not found in the current directory the routine
    looks for it in the FiberReg subdirectory before giving up.
    '''
    if not os.path.isfile(xfile):
        fallback=os.path.join('FiberReg',os.path.basename(xfile))
        if os.path.isfile(fallback):
            xfile=fallback
    try:
        f=open(xfile)
        lines=f.readlines()
    except:
        print('Error: GetRegSpec: Could not read regionfile:',xfile)
        return

    fibers=[]
    colors=[]
    for line in lines:
        if line.count('global'):
            word=line.split()
            for one_word in word:
                if one_word.count('color='):
                    xcolor=one_word.replace('color=','')
                    
        if line.count('text='):
            word=line.split()
            zcolor=xcolor
            for one_word in word:
                if one_word.count('text'):
                    foo=one_word.replace('text={','')
                    foo=foo.replace('}','') 
                    xfib=int(foo)
                if one_word.count('color='):
                     zcolor=one_word.replace('color=','')
            fibers.append(xfib)
            colors.append(zcolor)

    xtab=Table([fibers,colors],names=['fiberid','color'])

    if color==None:
        print('There was no defined color, extracting the fibers with the least number of fibers of a color')
        # Return the fibers that have the fewest of a given color
        qcolor,qcount=np.unique(xtab['color'],return_counts=True)
        qtab=Table([qcolor,qcount],names=['color','ncount'])
        qtab.sort('ncount')
        color=qtab['color'][0]
        print(qtab)


    # print('Selecting :', color)
    xtab=xtab[xtab['color']==color]

    return xtab



    





def get_spec(filename,xfib,xtype='ave'):
    '''
    Retrieve the spectra from the first
    nfib fibers in the xfib table

    Returns (xspec, n_good) where n_good is the number of fibers that
    passed the fibstatus filter and were actually included in the average.
    This count is needed by do_one to correctly scale the sum spectrum.
    '''

    x=fits.open(filename)
    # determine what extensiosn exist

    extensions=[]
    sky_exists=False
    i=1
    while i <len(x):
        one_name=x[i].header.get('EXTNAME')
        # print(one_name)
        if one_name.count('SKY'):
            sky_exists=True
            print('Sky Exists')
        i+=1



    try:
        xtab=Table(x['SLITMAP'].data)
        foo_tab=xtab[xfib['fiberid']-1]
        xfib['fibstatus']=foo_tab['fibstatus']
        xfib=xfib[xfib['fibstatus']==0]
        if len(foo_tab)-len(xfib)>0:
            print('Of %d possible fibers, %d were rejected for fibstatus leaving %d' % (len(foo_tab),len(foo_tab)-len(xfib),len(xfib)))
    except:
        print('No SLITMAP')
    wave=x['WAVE'].data
    xxfib=xfib.copy()
    flux=x['FLUX'].data[xxfib['fiberid']-1]
    ivar=x['IVAR'].data[xxfib['fiberid']-1]
    mask=x['MASK'].data[xxfib['fiberid']-1]
    if sky_exists:
        sky=x['SKY'].data[xxfib['fiberid']-1]
        sky_ivar=x['SKY_IVAR'].data[xxfib['fiberid']-1]
        lsf=x['LSF'].data[xxfib['fiberid']-1]


    if xtype=='ave':
        xflux=np.nanmean(flux,axis=0)
        xmask=np.sum(mask,axis=0)
        if sky_exists:
            xsky=np.nanmean(sky,axis=0)
            xlsf=np.nanmean(lsf,axis=0)
    elif xtype=='med':
        xflux=np.nanmedian(flux,axis=0)
        xmask=np.sum(mask,axis=0)
        if sky_exists:
            xsky=np.nanmedian(sky,axis=0)
            xlsf=np.nanmedian(lsf,axis=0)
    else:
        print('Error: getspec: only ave or med allowd for xtype')
        return

    ivar=np.nansum(ivar,axis=0)
    xerr = np.full(ivar.shape, np.nan)
    valid = ivar > 1
    xerr[valid] = 1/np.sqrt(ivar[valid])
    if sky_exists:
        xsky_error=np.nansum(sky_ivar,axis=0)
        valid_sky = xsky_error > 1
        xsky_err = np.full(xsky_error.shape, np.nan)
        xsky_err[valid_sky] = 1/np.sqrt(xsky_error[valid_sky])
        xsky_error = xsky_err
        xspec=Table([wave,xflux,xerr,xsky,xsky_error,xmask,xlsf],names=['WAVE','FLUX','ERROR','SKY','SKY_ERROR','MASK','LSF'])
        xspec['SKY'].format='.3e'
        xspec['SKY_ERROR'].format='.3e'
        xspec['LSF'].format='.2f'
    else:
        xspec=Table([wave,xflux,xerr,xmask],names=['WAVE','FLUX','ERROR','MASK'])
    xspec['WAVE'].format='.1f'
    xspec['FLUX'].format='.3e'
    xspec['ERROR'].format='.3e'
    return xspec, len(xxfib)


def do_one(filename,source_reg,source_reg_color,back_reg=None, back_reg_color=None, xtype='ave',root='Spec'):
    '''
    Extract one spectrum and write it to a txt file.

    Parameters:
        filename (str): The rss spectra file.
        source_reg (str): Region file with all of the fibers.
        source_reg_color (str): Color defining what fibers to extract for the source.
        back_reg (str): Optional region file for background (possibly the same as source_reg).
        back_reg_color (str): Optional color defining what fibers to extract for background.
        xtype (str): Type of spectrum to return (ave, med, or sum). At present the median background spectrum is always extracted.
        root (str): Root name for the output file.

    Returns:
        str: The name of the output file.
    '''

    try:
        x=fits.open(filename)
    except:
        print('Error: could not open %s' % filename)
        return

    print('\n Starting :', filename)

    source_fibers=read_reg(source_reg,source_reg_color)
    print('Source from %d fibers with color: %s' % (len(source_fibers),source_reg_color))

    if len(source_fibers)==0:
        print('Error: No good science fibers for source')
        return

    if xtype=='sum':
        xspec, n_source_good=get_spec(filename=filename,xfib=source_fibers,xtype='ave')
    else:
        xspec, _=get_spec(filename=filename,xfib=source_fibers,xtype=xtype)

    # print('Taking spectra from: ',list(source_fibers['fiberid']))
    if back_reg!=None:
        bfibers=read_reg(back_reg,back_reg_color)
        print('Getting Background from %d fibers with color: %s ' % (len(bfibers),back_reg_color))
        # print('Taking background from: ',list(bfibers['fiberid']))
        bspec, _=get_spec(filename=filename,xfib=bfibers,xtype='med')
        xspec['SOURCE_FLUX']=xspec['FLUX']
        xspec['SOURCE_ERROR']=xspec['ERROR']
        xspec['FLUX']-=bspec['FLUX']
        xspec['ERROR']=np.sqrt(xspec['ERROR']*xspec['ERROR']+bspec['ERROR']*bspec['ERROR'])
        xspec['BACK_FLUX']=bspec['FLUX']
        xspec['BACK_ERROR']=bspec['ERROR']

    if xtype=='sum':
        n_total=len(source_fibers)
        print('Scaling sum spectrum by %d total fibers in region (%d good, %d bad assumed equal to mean)' % (n_total, n_source_good, n_total-n_source_good))
        xspec['FLUX']*=n_total
        xspec['ERROR']*=n_total
        if 'SKY' in xspec.colnames:
            xspec['SKY']*=n_total
            xspec['SKY_ERROR']*=n_total
        if back_reg!=None:
            xspec['BACK_FLUX']*=n_total
            xspec['BACK_ERROR']*=n_total
            xspec['SOURCE_FLUX']*=n_total
            xspec['SOURCE_ERROR']*=n_total

    # Add code do make a directory for storing the results if the root name contains a directory

    if root.count('/'):
        xdir=os.path.split(root)[0]
        os.makedirs(xdir,exist_ok=True)

    if filename.count('SFra'):
        exposure=x['PRIMARY'].header['EXPOSURE']
        outname='%s_%05d' % (root,exposure)
    else:
        word=filename.split('/')
        xroot=word[-1]
        xroot=xroot.replace('.fits','')
        xroot=xroot.replace('.gz','')
        outname='%s_%s' % (root,xroot)

    if xtype=='med':
        outname='%s_med' % (outname)
    elif xtype=='sum':
        outname='%s_sum' % (outname)
    else:
        outname='%s_ave' % (outname)

    if back_reg!=None:
        outname='%s_back' % (outname)




    outname='%s.txt'% outname


    xspec.write(outname,format='ascii.fixed_width_two_line',overwrite=True)

    print('The output file is %s' % (outname))
    return outname




def do_all(masterfile, snap_dir='Snap', out_dir='Snap_spec', c_type='ave', xtype='ave', img_dir='', reg_dir='FiberReg'):
    '''
    Process all sources in a masterfile, creating a background-subtracted
    spectrum for each one.

    For every unique Source_name in masterfile the routine:

    1. Locates the snapshot file snap_dir/<source_name>.<c_type>.fits
    2. Calls MakeLVMReg.do_complex to assign fiber colors (red=source,
       green=background annulus) using the region geometry in masterfile.
    3. Calls do_one to extract the average source spectrum (red fibers)
       and subtract the median background (green fibers).
    4. If img_dir is provided, calls LSnap.make_one_image for every FITS
       file in img_dir whose name begins with the source name, overlaying
       the fiber color assignments on each broadband image.

    Parameters:
        masterfile (str): Source catalog with SourceBack, Source_name,
            RA, Dec, RegType, Major, Minor, and Theta columns, such as
            that produced by GenAnnularBackground.py.
        snap_dir (str): Directory containing snapshot FITS files
            (default: 'Snap').
        out_dir (str): Output directory for extracted spectra
            (default: 'Snap_spec').
        c_type (str): Combination type suffix used in snapshot filenames,
            e.g. 'ave' or 'med' (default: 'ave').
        xtype (str): How to combine source fibers: 'ave' or 'med'
            (default: 'ave').
        reg_dir (str): Directory in which region files are written by
            MakeLVMReg.do_complex.  Created automatically if it does not
            exist (default: 'FiberReg').
        img_dir (str): Directory containing per-source broadband image
            cutouts (e.g. MCELS H-alpha and [SII] FITS files).  If set,
            LSnap.make_one_image is called for each matching file after
            spectrum extraction.  Output images are written to zimage/.
            (default: '', meaning LSnap is not called).
    '''
    import glob
    from lvm_ksl import MakeLVMReg
    from astropy.io import ascii as asc

    if img_dir:
        from lvm_ksl import LSnap

    xtab = asc.read(masterfile)
    sources = np.unique(xtab['Source_name'])
    has_back = 'SourceBack' in xtab.colnames

    os.makedirs(out_dir, exist_ok=True)

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
        reg_file = MakeLVMReg.do_complex(filename=filename, qtab=source_tab,
                                         outroot='source', buffer=17.5, reg_dir=reg_dir)
        if reg_file is None:
            print('Skipping %s: could not create region file' % source_name)
            continue

        back_reg   = reg_file if has_back else None
        back_color = 'green'  if has_back else None

        do_one(filename=filename,
               source_reg=reg_file, source_reg_color='red',
               back_reg=back_reg, back_reg_color=back_color,
               xtype=xtype, root='%s/Spec' % out_dir)

        if img_dir:
            img_files = glob.glob('%s/%s*.fits*' % (img_dir, source_name))
            if img_files:
                for img_file in sorted(img_files):
                    LSnap.make_one_image(img_file, reg_file, None, None)
            else:
                print('No broadband images found for %s in %s' % (source_name, img_dir))


def steer(argv):
    '''
    This routine parses the command line

    and then oversees the rest of the process


    In future it would be better to split these
    acitvites apart so the routine can be called
    from a Jupyter script, but that is not
    the way things are at present
    usage: GetRegSpec.py [-h] [-root whatever] [-med] [-sum] filename[s] source_reg [color] [back_reg] [back_color]
    usage: GetRegSpec.py [-h] [-all] [-snap snap_dir] [-out out_dir] [-med] [-ctype ave|med] masterfile
    '''

    source_reg=None
    source_reg_color=None
    back_reg=None
    back_reg_color=None
    root='Spec'
    xtype='ave'
    filenames=[]
    xall=False
    masterfile=''
    snap_dir='Snap'
    out_dir='Snap_spec'
    c_type='ave'
    reg_dir='FiberReg'
    img_dir=''

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-root':
            i+=1
            root=argv[i]
        elif argv[i][0:4]=='-med':
            xtype='med'
        elif argv[i][0:4]=='-sum':
            xtype='sum'
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-snap':
            i+=1
            snap_dir=argv[i]
        elif argv[i]=='-out':
            i+=1
            out_dir=argv[i]
        elif argv[i]=='-ctype':
            i+=1
            c_type=argv[i]
        elif argv[i]=='-regdir':
            i+=1
            reg_dir=argv[i]
        elif argv[i]=='-imgdir':
            i+=1
            img_dir=argv[i]
        elif argv[i][0]=='-':
            print('Error: Incorrect Command line: ',argv)
            return
        elif xall and masterfile=='':
            masterfile=argv[i]
        elif argv[i].count('.fits')>0:
            filenames.append(argv[i])
        elif source_reg==None and argv[i].count('reg')>0:
            source_reg=argv[i]
        elif back_reg==None and source_reg_color==None and argv[i].count('reg')==0:
            source_reg_color=argv[i]
        elif back_reg==None and argv[i].count('reg')>0:
            back_reg=argv[i]
        elif back_reg!=None and back_reg_color==None:
            back_reg_color=argv[i]
        else:
            print('Error: Incorrect Command line: (extra args) ',argv)
            return
        i+=1

    if xall:
        if masterfile=='':
            print('Error: -all requires a masterfile argument')
            return
        do_all(masterfile, snap_dir=snap_dir, out_dir=out_dir, c_type=c_type, xtype=xtype, img_dir=img_dir, reg_dir=reg_dir)
        return

    if len(filenames)==0:
        print('No files to process: ',argv)
        return

    for filename in filenames:
        print('Starting %s' % filename)
        do_one(filename,source_reg,source_reg_color,back_reg, back_reg_color, xtype,root)


    return








    


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
