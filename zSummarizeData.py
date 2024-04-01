#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Summarize the unprocessed LVM data that exists, creating
and astropy table with information sufficient to allow one
to identify all data associated with a specific object or
position on the sky


Command line usage (if any):

    usage: SummarizeData.py [-h] [-all] [-redo]  [mjd1 mjd2 ...]

    where:

        -h will pring this summary and exit
        -all means to summarize all of the data in various subdirecories
        -redo means to recreate the summary tables, if and only if -all has been
            selected
        [mjd1 mjd2]  run on for specific mjds


Description:  

    The routine reads the headers of the various fits files that are in the various mjd direcories
    and produces and astopy table for each.  It then stacks these tables so that one has all of the
    the files summarized as sell

Primary routines:

    do_summay produces the table for a single day
    stack stacks all fo the tables to produce a single file
    steer controls the overall flow

Notes:

    The output table names are all hardwired.

    The summed result is in All_data.txt

    As written the program needs to be run from the lvm_data/sas/sdsswork/data/lvm/lco directory
    but this would be fairly easy to change

    It would also be straighforward to parallelize.
                                       
History:

231003 ksl Coding begun

'''

import sys
import os
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from glob import glob
from astropy.table import vstack
import numpy as np
import timeit


LVM_DATA_S=os.environ.get('LVM_DATA_S')


def do_summary(directory='60202'):
    print('hello world')
    search_string='%s/%s/*gz' % (LVM_DATA_S,directory)
    print(search_string)
    xfiles=glob(search_string)
    # print(xfiles)
    name=[]
    xobject=[]
    mjd=[]
    xtime=[]
    exptime=[]
    exposure=[]
    ra=[]
    dec=[]
    ra_e=[]
    dec_e=[]
    ra_w=[]
    dec_w=[]
    xtype=[]
    for one_file in xfiles:
        x=fits.open(one_file)
        head=x[0].header
        try:
            name.append(head['FILENAME'])
        except:
            name.append('Unknown')
            
        try:
            xobject.append(head['OBJECT'])
        except:
            xobject.append('Unknown')
            
        try:
            xtype.append(head['IMAGETYP'])
        except:
            xtype.append('Unknown')
            
        try:
            mjd.append(head['MJD'])
        except:
            mjd.append(-99)
            
        try:
            xtime.append(head['OBSTIME'])
        except:
            xtime.append('Unknown')
            
        try:
            exptime.append(head['EXPTIME'])
        except:
            exptime.append(-99)
        
        try:
            exposure.append(head['EXPOSURE'])
        except:
            exposure.append('Unknown')
            
        try:
            ra.append(head['TESCIRA'])
        except:
            try:
                ra.append(head['SCIRA'])
            except:
                ra.append(-99.0)
            
        try:
            dec.append(head['TESCIDE'])
        except:
            try:
                dec.append(head['SCIDEC'])
            except:
                dec.append(-99.0)  
            
        try:
            ra_e.append(head['TESKYERA'])
            dec_e.append(head['TESKYEDE'])
        except:
            try:
                ra_e.append(head['SKYERA'])
                dec_e.append(head['SKYEDEC'])
            except:
                ra_e.append(-99.0)
                dec_e.append(-99.0)
            
        try:
            ra_w.append(head['TESKYWRA'])
            dec_w.append(head['TESKYWDE'])
        except:
            try:
                ra_w.append(head['SKYWRA'])
                dec_w.append(head['SKYWDEC'])
            except:
                ra_w.append(-99.0)
                dec_w.append(-99.0)        
        
        
    
    xtab=Table([exposure,xtype,name,mjd,xtime,ra,dec,ra_e,dec_e,ra_w,dec_w, exptime,xobject],
               names=['Exposure','Type','File','MJD','Time','RA','Dec','RA_E', 'Dec_E','RA_W','Dec_W','Exptime','Object'])
    
    xtab.sort(['Exposure','File'])
    
    words=directory.split('/')
    
    xtab.write('%s/Sum_%s.txt' % (os.environ.get('LVM_KSL_OUT'),words[-1]),format='ascii.fixed_width_two_line',overwrite=True)

    
    return xtab
                    
def stack():
    '''
    Stack all of the Sum files into a single table 
    '''

    files=glob('Sum*.txt')

    x=ascii.read(files[0])
    i=1
    while i<len(files):
        try:
            y=ascii.read(files[i])
            x=vstack([x,y])
        except:
            print('Could not stack %s' % files[i])
        i+=1

    x.sort('Exposure')
    x.write('All_Data.txt',format='ascii.fixed_width_two_line',overwrite=True)
    return

def steer(argv):
    '''
    This is just a steering routine
    '''

    print('xtest ',argv)
    xdays=[]
    xall=False
    redo=False
    xxrange=False


    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            xall=True
        elif argv[i]=='-redo':
            redo=True
        elif argv[i]=='-range':
            xxrange=True
        elif argv[i][0]=='-':
            print('Error: Unknown switch %s ' % argv[i])
            return
        else:
            xdays.append(argv[i])
        i+=1

    if xxrange:
        # convert to integers
        xxint=[]
        for one in xdays:
            xxint.append(int(one))
        xday_min=np.min(xxint)
        xday_max=np.max(xxint)
        print('Work on days %d to %d' % (xday_min,xday_max))
        xdays=[]
        qday=xday_min
        while qday<=xday_max:
            xdays.append('%d' % qday)
            qday+=1
        print(xdays)

    if xall:

        xx=glob('%s/*/*fits*' % (LVM_DIR_S))
        xdir=[]
        for one in xx:
            word=one.split('/')
            xdir.append(word[-2])
        xxdays=np.unique(xdir)
        print(xxdays)
        return

        if redo==False:
            for one in xxdays:
                if os.path.isfile('Sum_%s.txt' % one):
                    print('Omitting %s because Sum_%s.txt exists' % (one,one))
                else:
                    xdays.append(one)
        else:
            xdays=xxdays

    good_days=[]
    for one in xdays:
        xdir='%s/%s' % (LVM_DATA_S,one)
        if os.path.isdir(xdir):
            good_days.append(one)
        else:
            print('Sorry: %s does not appear to be a valid directory' % one)
    xdays=good_days




    if len(xdays)==0:
        print('Sorry: there seems to be nothing to do')
        print('-all not set and no xdays to set up provided' )
        return

    i=1
    for one in xdays:
        time_start=timeit.default_timer()
        print('Starting mjd %s (%d of %d)' % (one,i,len(xdays)))
        do_summary(directory=one)
        elapsed = timeit.default_timer() - time_start
        print('Finished mjd %s in %.1f s' % (one,elapsed))
        i+=1


    stack()


    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
       steer(sys.argv)
    else:
        print(__doc__)

