#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Summmarize the median spectrum in Cframe data for a large
number of files


Command line usage (if any):

    usage: SumCFrame [-h] [-emin 900] [-out whatever] [-ver 1.1.3] [-drp_all drp_file] exp_start exp_stop [delta]

Description:  

    where -h prints out this help and quits
            =emin 900           resets the minimum exposure to another value
            -out whatever       sets the root for the output fits file to soemthing other than the default
            -ver 1.1.3          resets the location to lookd for the drp_all file
            -drp_all drp_file   causes the program to read a speccial drp_all file or table named 'drp_file'
                                as opposed to the drp_all file associated with the version of processing
            exp_start           The starting exposure to consider
            exp_stop            The exposure snumber to stop on
            [delta]             The number of 'good' exposures to skip


Primary routines:

    doit

Notes:
                                       
History:

240726 ksl Coding begun

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.table import join, Table
import shutil
from datetime import datetime
from astropy.wcs import WCS


def read_drpall(filename='',drp_ver='1.1.0'):
    '''
    This reads a drp_all fits file or
    an ascii table that contains the table
    and returns the table
    '''

    if filename.count('txt') or filename.count('.tab'):
        try:
            drp_tab=ascii.read(filename)
            return drp_tab
        except:
            print('Error: Could not locate : ', filename)
            return []
    
    elif filename=='':
        DRPFILE='drpall-%s.fits' % (drp_ver)
    else:
        DRPFILE=filename

    # First try to locate the DRP file locally, otherwise
    if os.path.isfile(DRPFILE):
        xfile=DRPFILE
    else:
        BASEDIR='/uufs/chpc.utah.edu/common/home/sdss51/sdsswork/lvm/spectro/redux/%s/' % (drp_ver)
        xfile='%s/%s' % (BASEDIR,DRPFILE)
        if os.path.isfile(xfile)==False:
            print('Error: Could not locate : ', xfile)
            return []

    try:
        drpall=fits.open(xfile)
    except:
        print('Error: Locared but could not read  : ', xfile)
        return []

    drp_tab=Table(drpall[1].data)
    return  drp_tab


def xselect(ztab,exp_start=4000,exp_stop=8000,delta=5,exp_min=900.):
    '''
    select every nth exposure between two start and stop times
    in a drpall table
    '''
    xtab=ztab[ztab['expnum']>=exp_start]
    xtab=xtab[xtab['expnum']<=exp_stop]
    # print(len(xtab))
    # print(np.unique(xtab['exptime'],return_counts=True))
    if exp_min>0:
        xtab=xtab[xtab['exptime']>=exp_min]
    # print(len(xtab))
    # print(np.unique(xtab['exptime'],return_counts=True))
    if delta>1:
        xtab=xtab[::delta]
    return xtab


XTOP='/uufs/chpc.utah.edu/common/home/sdss51/'
XRAINBOW='/Users/long/Projects/lvm_data/sas'
XMUSKIE='/home/long/Projects/lvm_data/sas'

def find_top():
    if os.path.isdir(XTOP):
        loc='Utah'
        topdir=XTOP
    elif os.path.isdir(XRAINBOW):
        loc='Rainbow'
        topdir=XRAINBOW
    elif ox.path.isdir(XMUSKIE):
        loc='Muskie'
        topdir=XMUSKIE
    else:
        print('Error: I donot know where I am:', os.pwd())
        return ''

    print('We am on : ', loc)
    return topdir

        

def scifib(xtab,select='all',telescope=''):
    '''
    Select good fibers from a telescope, of a spefic
    type or all from a telescope from the slitmap table
    of a calbrated file
    '''
    # print(np.unique(xtab['fibstatus']))
    # print(np.unique(xtab['targettype']))
    ztab=xtab[xtab['fibstatus']==0]
    if select=='all':
        ztab=ztab[ztab['targettype']!='standard']
    else:
        ztab=ztab[ztab['targettype']==select]

    if telescope!='' and telescope!='all':
        ztab=ztab[ztab['telescope']==telescope]


    # print('Found %d fibers' % len(ztab))
    return ztab


def sum_frames(xfiles):
    '''
    Given a bunch of rss files containing a FLUX extension, this
    returns the average spectrum for each fiber
    '''
    print('There are %d files to process' % len(xfiles))
    i=0
    for one in xfiles:
        x=fits.open(one)
        if i==0:
            xsum=np.zeros_like(x['FLUX'].data)
            xnorm=np.zeros_like(x['FLUX'].data,)
        q=np.ma.array(x['FLUX'].data,mask=x['MASK'].data,fill_value=np.nan)
        q = q.filled(np.nan)
        xnorm=np.select([x['MASK'].data==0],[xnorm+1],default=0)
        xsum=np.nansum([xsum,q],axis=0)
        i+=1
        if i%10==0:
            print('%d of %d files have been processed' % (i,len(xfiles)))

    xsum/=xnorm
    print('Finished xsum',xsum.shape)
    return xsum

def files_select(filename='',exp_start=4000,exp_stop=8000,delta=5,exp_min=900.,drp_ver='1.1.0'):
    '''
    Read the drp_all table and select files to process from this table
    '''
    data_dir=find_top()
    xtab=read_drpall(filename,drp_ver)
    print('got drp all')
    xtab=xselect(xtab,exp_start,exp_stop,delta)
    return xtab


def process_files(xtab,out_name=''):
    '''
    Having decided what needs processing do the work
    '''
    data_dir=find_top()
    i=0
    select=[]
    xfiles=[]
    while i < len(xtab):
        xfile='%s/%s' % (data_dir,xtab['location'][i])
        if xfile.count('SFrame'):
            xfile=xfile.replace('SFrame','CFrame')
        # print(xfile)
        if os.path.isfile(xfile):
            select.append(i)
            xfiles.append(xfile)
        i+=1

    print(xfiles)
    results=sum_frames(xfiles)

    # At this point we have the results and we just need to build an ouput file

    xdummy=fits.open(xfiles[0])  # so we have something to work with

    hdu1 = fits.PrimaryHDU(data=None)
    hdu1.header['Title'] = 'CFrame_Summmary'
    hdu2=xdummy['FLUX']
    hdu2.data=results
    hdu3= fits.ImageHDU(data=xdummy['WAVE'].data,name='WAVE')
    hdu4=xdummy['SLITMAP']
    hdu5 = fits.BinTableHDU(xtab, name='drp_all')

    hdul = fits.HDUList([hdu1, hdu2, hdu3,hdu4,hdu5])

    hdul.writeto(out_name,overwrite=True)
    return


def doit(filename='',exp_start=4000,exp_stop=8000,delta=5,exp_min=900.,out_name='',drp_ver='1.1.0'):
    '''
    The main routine allowing one to process the data
    '''
    xtab=files_select(filename,exp_start,exp_stop,delta,exp_min,drp_ver)

    if out_name=='':
        out_name='XCframeSum_%d_%d_%d.fits' % (exp_start,exp_stop,delta)
    
    if out_name.count('.fits')==0:
        out_name='%s.fits' % (out_name)


    if len(xtab)>0:
        process_files(xtab,out_name)
    return


def steer(argv):
    '''
    SumCFrame -h -emin 1000 -[-out whatever] [-drp_all drp_file] exp_start expstop [delta]
    '''
    exp_start=-1
    exp_stop=-1
    delta=-1

    exp_min=900
    input_file=''
    out_name=''
    xver='1.1.0'

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-emin':
            i+=1
            exp_min=int(argv[i])
        elif argv[i]=='-out':
            i+=1
            out_name=(argv[i])
        elif argv[i]=='-drp_all':
            i+=1
            input_file=argv[i]
        elif argv[i]=='-ver':
            i+=1
            xver=argv[i]
        elif argv[i][0]=='-':
            print('Unknown option : ',argv)
        elif exp_start<0:
            exp_start=int(argv[i])
        elif exp_stop<0:
            exp_stop=int(argv[i])
        elif delta<0:
            delta=int(argv[i])
        i+=1

    if delta<0:
        delta=1
                

    doit(input_file,exp_start,exp_stop,delta,exp_min,out_name,xver)




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
