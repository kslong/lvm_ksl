#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Get a spectrum from a fiber or fibers in an RSSfile given an 
RA and DEC based on set of selected region files



Command line usage (if any):

    usage: GetRegSpec.py  [-h] [-root whatever] -median -sum filename source_reg  [-back back_reg]



    where 
    -h prints out this help
    -median constructs the median instead of the average spectrum
    -sum constructs the sum instead of the median or average
    -root whatever prepends a root to the standard file name
    exp_no is the exposure
    ra and dec are the desired right ascension, written
        either in degrees or in h:m:s d:m:s
    [rmin] is and optional inner radius, and 
    [rmax] is an optinal outer radius
    -back a1 a2 subtract a background

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



def read_reg(xfile,color=None):
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
                     print('gotcha', zcolor)
            fibers.append(xfib)
            colors.append(zcolor)

    xtab=Table([fibers,colors],names=['fiberid','color'])

    if color==None:
        # Return the fibers that have the fewest of a given color
        qcolor,qcount=np.unique(xtab['color'],return_counts=True)
        qtab=Table([qcolor,qcount],names=['color','ncount'])
        qtab.sort('ncount')
        color=qtab['color'][0]
        print('Hi knox')
        print(qtab)


    print('Selecting :', color)
    xtab=xtab[xtab['color']==color]




    return xtab



    





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






def steer(argv):
    '''
    This routine parses the command line

    and then oversees the rest of the proces


    In future it would be better to split these
    acitvites apart so the routine can be called
    from a Jupyter script, but that is not
    the way things are at present
    usage: GetRegSpec.py  [-h] [-root whatever] -median -sum filename source_reg color [-back back_reg back_color]
    '''

    filename=None
    source_reg=None
    source_reg_color=None
    back_reg=None
    back_reg_color=None
    root='Spec'
    xtype='ave'
    back=False

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
        elif filename==None:
            filename=argv[i]
        elif source_reg==None:
            source_reg=argv[i]
        elif source_reg_color==None and argv[i].count('reg')==0:
            source_reg_color=argv[i]
        elif back_reg==None:
            back_reg=argv[i]
        elif back_reg_color==None:
            back_reg_color=argv[i]
        else:
            print('Error: Incorrect Command line: (exta args) ',argv)
            return
        i+=1


    print('source.reg',source_reg,source_reg_color)
    print('back.reg',back_reg,back_reg_color)


    source_fibers=read_reg(source_reg,source_reg_color)
    print(source_fibers)


    try:
        x=fits.open(filename)
    except:
        print('Error: could not open %s' % filename)
        return




    if len(source_fibers)==0:
        return

    if xtype=='sum':
        xspec=get_spec(filename=filename,xfib=source_fibers,nfib=len(source_fibers),xtype='ave')
    else:
        xspec=get_spec(filename=filename,xfib=source_fibers,nfib=len(source_fibers),xtype=xtype)

    if back_reg!=None:
        bfibers=read_reg(back_reg,back_reg_color)
        bspec=get_spec(filename=filename,xfib=bfibers,nfib=len(bfibers),xtype='med')
        xspec['SOURCE_FLUX']=xspec['FLUX']
        xspec['SOURCE_ERROR']=xspec['ERROR']
        xspec['FLUX']-=bspec['FLUX']
        xspec['ERROR']=np.sqrt(xspec['ERROR']*xspec['ERROR']+bspec['ERROR']*bspec['ERROR'])
        xspec['BACK_FLUX']=bspec['FLUX']
        xspec['BACK_ERROR']=bspec['ERROR']

    print('Taking spectra from\n',source_fibers['fiberid'])
    if xtype=='sum':
        xspec['FLUX']*=len(source_fibers)
        xspec['ERROR']*=len(source_fibers)
        xspec['SKY']*=len(source_fibers)
        xspec['SKY_ERROR']*=len(source_fibers)
        if back:
            xspec['BACK_FLUX']*=len(source_fibers)
            xspec['BACK_ERROR']*=len(source_fibers)
            xspec['SOURCE_FLUX']*=len(source_fibers)
            xspec['SOURCE_ERROR']*=len(source_fibers)

    exposure=x['PRIMARY'].header['EXPOSURE']


    outname='%s_%05d_test' % (root,exposure)

    if xtype=='med':
        outname='%s_med' % (outname)
    elif xtype=='sum':
        outname='%s_sum' % (outname)
    else:
        outname='%s_ave' % (outname)

    if back_reg!=None:
        outname='%s_back' % (outname)



    # xspec.meta['comments']=['Filename %s' % filename,'RA %.5f' % ra, 'Dec %.5f' % dec, 'nfibers %d' % (len(fibers))]


    xspec.write('%s.txt'% outname,format='ascii.fixed_width_two_line',overwrite=True)

    print('The output file is %s.txt' % (outname))












    


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
