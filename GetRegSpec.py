#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Get a spectrum from a fiber or fibers in an RSSfile given an 
RA and DEC based on a region file that conntains the one 
region for fiber.  



Command line usage (if any):

    usage: GetRegSpec.py  [-h] [-root whatever] [-median]-sum filename[s] source_reg  [color] [back_reg] [back_color}



    where 
    -h prints out this help
    -median constructs the median instead of the average spectrum
    -sum constructs the sum instead of the median or average
    -root whatever prepends a root to the standard file name

    filename[s] one or more rss spectra fiiles
    source_reg is the name of a region file (that must have an extension .reg)
    [color] is to define the color of fibers to be extracted
    back_reg is an optional background region file
    back_color is an optinal back_coler

    The lines source_reg color back_reg and back_color need to be in order.
    

Description:  

    The routine looks for a file where the RA's and Dec's have
    been assigned to fibers (by MakeLVMReg.py) and finds the
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

    The routine also produces a region file that shows
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



def read_reg(xfile,color=None):
    '''
    This reads a region file that contains
    all a region for each of the fibers.  
    It either selects those fibers of a certain
    color or if color is None, it assumes that
    one wants the fibers which have the fewest
    fibers of a given color.
    '''
    try:
        f=open(xfile)
        lines=f.readlines()
    except:
        print('Error: Could not read :',xfile)
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
    xerr=np.select([ivar>1],[1/np.sqrt(ivar)],default=np.nan)
    if sky_exists:
        xsky_error=np.nansum(sky_ivar,axis=0)
        xsky_error=np.select([xsky_error>1],[1/np.sqrt(xsky_error)],default=np.nan)
        xspec=Table([wave,xflux,xerr,xsky,xsky_error,xmask,xlsf],names=['WAVE','FLUX','ERROR','SKY','SKY_ERROR','MASK','LSF'])
        xspec['SKY'].format='.3e'
        xspec['SKY_ERROR'].format='.3e'
        xspec['LSF'].format='.2f'
    else:
        xspec=Table([wave,xflux,xerr,xmask],names=['WAVE','FLUX','ERROR','MASK'])
    xspec['WAVE'].format='.1f'
    xspec['FLUX'].format='.3e'
    xspec['ERROR'].format='.3e'
    return xspec


def do_one(filename,source_reg,source_reg_color,back_reg=None, back_reg_color=None, xtype='ave',root='Spec'):

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
        xspec=get_spec(filename=filename,xfib=source_fibers,xtype='ave')
    else:
        xspec=get_spec(filename=filename,xfib=source_fibers,xtype=xtype)

    # print('Taking spectra from: ',list(source_fibers['fiberid']))
    if back_reg!=None:
        bfibers=read_reg(back_reg,back_reg_color)
        print('Getting Background from %d fibers with color: %s ' % (len(bfibers),back_reg_color))
        # print('Taking background from: ',list(bfibers['fiberid']))
        bspec=get_spec(filename=filename,xfib=bfibers,xtype='med')
        xspec['SOURCE_FLUX']=xspec['FLUX']
        xspec['SOURCE_ERROR']=xspec['ERROR']
        xspec['FLUX']-=bspec['FLUX']
        xspec['ERROR']=np.sqrt(xspec['ERROR']*xspec['ERROR']+bspec['ERROR']*bspec['ERROR'])
        xspec['BACK_FLUX']=bspec['FLUX']
        xspec['BACK_ERROR']=bspec['ERROR']

    if xtype=='sum':
        xspec['FLUX']*=len(source_fibers)
        xspec['ERROR']*=len(source_fibers)
        xspec['SKY']*=len(source_fibers)
        xspec['SKY_ERROR']*=len(source_fibers)
        if back_reg!=None:
            xspec['BACK_FLUX']*=len(source_fibers)
            xspec['BACK_ERROR']*=len(source_fibers)
            xspec['SOURCE_FLUX']*=len(source_fibers)
            xspec['SOURCE_ERROR']*=len(source_fibers)

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






def steer(argv):
    '''
    This routine parses the command line

    and then oversees the rest of the process


    In future it would be better to split these
    acitvites apart so the routine can be called
    from a Jupyter script, but that is not
    the way things are at present
    usage: GetRegSpec.py  [-h] [-root whatever] -median -sum filename source_reg color [-back back_reg back_color]
    '''

    source_reg=None
    source_reg_color=None
    back_reg=None
    back_reg_color=None
    root='Spec'
    xtype='ave'
    back=False
    filenames=[]


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
        elif argv[i][0]=='-':
            print('Error: Incorrect Command line: ',argv)
            return
        elif argv[i].count('.fits')>0:
            filenames.append(argv[i])
        elif source_reg==None and argv[i].count('reg')>0:
            source_reg=argv[i]
        elif back_reg==None and source_reg_color==None and argv[i].count('reg')==0:
            source_reg_color=argv[i]
        elif back_reg==None and argv[i].count('reg')>0 :
            back_reg=argv[i]
        elif back_reg !=None and back_reg_color==None:
            back_reg_color=argv[i]
        else:
            print('Error: Incorrect Command line: (exta args) ',argv)
            return
        i+=1

    if len(filenames)==0:
        print('No files to proces: ',argv)
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
