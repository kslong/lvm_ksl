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
240406 ksl Reworked to add more info, and hopefully be more robust

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
from lvm_ksl.xcal import get_moon_info_las_campanas

LVM_DATA_S=os.environ.get('LVM_DATA_S')

def get_files2use(directory='60202'):
    if LVM_DATA_S==None:
        xfiles=glob('%s/*gz' % directory)
    else:
        xfiles=glob('%s/%s/*gz' % (LVM_DATA_S,directory))

    exposures=[]
    for one_file in xfiles:
        words=one_file.split('-')
        exposures.append(words[-1].replace('.fits.gz',''))
    xtab=Table([xfiles,exposures],names=['Filename','ExpNo'])
    xtab.sort(['ExpNo','Filename'])

    
    exp_no,no_spec=np.unique(xtab['ExpNo'],return_counts=True)
    ytab=Table([exp_no,no_spec],names=['ExpNo','NSpec'])
    
    file2use=[]
    for one_row in ytab:
        qtab=xtab[xtab['ExpNo']==one_row['ExpNo']]
        if len(qtab)==0:
            file2use.append('Unknown')
        else:
            file2use.append(qtab['Filename'][0])
            
    ytab['File2Use']=file2use
    return ytab# 

def get_header_value(header, key, default_value=None, verbose=False):
    '''
    Robust way to get a header value if it exists
    '''

    try:
        value = header[key]
        if isinstance(value, str):
            try:
                value = float(value)  # or int(value) if it's an integer
            except ValueError as e:
                if verbose:
                    print(f"Failed to convert '{value}' to a number for key '{key}': {e}")
                value = default_value
    except KeyError as e:
        if verbose:
            print(f"Key '{key}' not found in header: {e}")
        value = default_value
    return value

def do_summary(directory='60202'):

    ytab=get_files2use(directory)

    xfiles=ytab['File2Use']
    nspec=ytab['NSpec']
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

    pa_sci=[]
    pa_e=[]
    pa_w=[]
    nstandards=[]

    moon_ra=[]
    moon_dec=[]
    moon_alt=[]
    moon_phase=[]
    moon_ill=[]


    for one_file in xfiles:
        # print(one_file)
        x=fits.open(one_file)
        head=x[0].header
        try:
            name.append(head['FILENAME'])
        except:
            name.append('Unknown')
            
        try:
            obj=head['OBJECT']
            if obj=='':
                obj='Unknown'
            xobject.append(obj)
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
            obstime=head['OBSTIME']
            xtime.append(obstime)
        except:
            xtime.append('Unknown')
            
        try:
            exptime.append(head['EXPTIME'])
        except:
            exptime.append(-99)
        
        try:
            exposure.append(head['EXPOSURE'])
        except:
            # exposure.append('Unknown')
            exposure.append(-99.0 )
            

        value=get_header_value(head, 'TESCIRA',-99.0)
        if value==-99.:
            value=get_header_value(head,'SCIRA',-99.0)
        ra.append(value)
            
        value=get_header_value(head,'TESCIDE',-99.)
        if value==-99.:
            value=get_header_value(head,'SCIDEC',-99.)
        dec.append(value)
            
        value=get_header_value(head,'TESKYERA',-99.)
        if value==-99.:
            value=get_header_value(head,'SKYERA',-99.)
        ra_e.append(value)
            
        value=get_header_value(head,'TESKYEDE',-99.)
        if value==-99.:
            value=get_header_value(head,'SKYEDEC',-99.)
        dec_e.append(value)

            
        value=get_header_value(head,'TESKYWRA',-99.)
        if value==-99.:
            value=get_header_value(head,'SKYWRA',-99.)
        ra_w.append(value)
            
        value=get_header_value(head,'TESKYWDE',-99.)
        if value==-99.:
            value=get_header_value(head,'SKYWDEC',-99.)
        dec_w.append(value)

        
        value=get_header_value(head,'POSCIPA',-999.)
        pa_sci.append(value)


        
        value=get_header_value(head,'POSKYEPA',-999.)
        pa_e.append(value)


        value=get_header_value(head,'POSKYWPA',-999.)
        pa_w.append(value)


        
        keys=list(head.keys())
        n=0
        for one_key in keys:
          if one_key.count('STD') and one_key.count('EXP') and  head[one_key]>0.0:
               # if one_key.count('STD') and one_key.count['EXP'] and head[one_key]>0.0:
               # print('gotcha')
               n+=1
        nstandards.append(n)



        try: 
            moon_info=get_moon_info_las_campanas(obstime)
            moon_alt.append(moon_info['MoonAlt'])
            moon_phase.append(moon_info['MoonPhas'])
            moon_ill.append(moon_info['MoonIll'])
            moon_ra.append(moon_info['MoonRA'])
            moon_dec.append(moon_info['MoonDec'])
        except Exception as e:
            print(e)
            moon_alt.append(-99.0)
            moon_phase.append(-99.0)
            moon_ill.append(-99.)
            moon_ra.append(-99.0)
            moon_dec.append(-99.0)

        
    # xtab=Table([exposure,nspec,xtype,name,mjd,xtime,ra,dec,pa_sci,ra_e,dec_e,pa_e,ra_w,dec_w, pa_w,nstandards,moon_alt,moon_phase,exptime,xobject],
    #            names=['Exposure','NSpec','Type','FileUsed','MJD','Time','RA','Dec','PA','RA_E', 'Dec_E','PA_E','RA_W','Dec_W','PA_W','NStandards','MoonAlt','MoonPhas','Exptime','Object'])
    
    xtab=Table([exposure,mjd,nspec,xtype,name,xtime,ra,dec,pa_sci,ra_e,dec_e,pa_e,ra_w,dec_w, pa_w,nstandards,moon_ra,moon_dec,moon_alt,moon_phase,moon_ill,exptime,xobject],
                names=['Exposure','MJD','NSpec','Type','FileUsed','Time','RA','Dec','PA','RA_E', 'Dec_E','PA_E','RA_W','Dec_W','PA_W','NStandards','MoonRA','MoonDec','MoonAlt','MoonPhas','MoonIll','Exptime','Object'])
    
    xtab.sort(['Exposure'])
    
    words=directory.split('/')

    if not os.path.isdir('./SumObs'):
        os.mkdir('./SumObs')
    
    xtab.write('SumObs/Sum_%s.txt' % words[-1],format='ascii.fixed_width_two_line',overwrite=True)
    
    return xtab
                    
def stack():
    '''
    Stack all of the Sum files into a single table 
    '''

    files=glob('SumObs/Sum*.txt')

    x=ascii.read(files[0])
    i=1
    while i<len(files):
        try:
            y=ascii.read(files[i])
            x=vstack([x,y])
        except:
            print('Could not stack %s' % files[i])
        i+=1

    x['RA'].format='.6f'
    x['Dec'].format='.6f'
    x['RA_E'].format='.6f'
    x['Dec_E'].format='.6f'
    x['RA_W'].format='.6f'
    x['Dec_W'].format='.6f'
    x['MoonRA'].format='.6f'
    x['MoonDec'].format='.6f'
    x['MoonAlt'].format='.2f'
    x['MoonPhas'].format='.2f'
    x['MoonIll'].format='.2f'
    x.sort('Exposure')
    x.write('All_Data.txt',format='ascii.fixed_width_two_line',overwrite=True)
    return

def steer(argv):
    '''
    This is just a steering routine
    '''
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
        if LVM_DATA_S==None:
            xx=glob('*/*fits*')
        else:
            xx=glob('%s/*/*fits*' % LVM_DATA_S)

        xdir=[]
        for one in xx:
            word=one.split('/')
            xdir.append(word[-2])
        xxdays=np.unique(xdir)
        print(xxdays)

        if redo==False:
            for one in xxdays:
                if os.path.isfile('SumObs/Sum_%s.txt' % one):
                    print('Omitting %s because Sum_%s.txt exists' % (one,one))
                else:
                    xdays.append(one)
        else:
            xdays=xxdays

    good_days=[]
    for one in xdays:
        if LVM_DATA_S==None and os.path.isdir(one):
            good_days.append(one)
        elif os.path.isdir('%s/%s' % (LVM_DATA_S,one)): 
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

