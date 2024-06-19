#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Locate and if desired retrieve calibrated LMV data
from series of LVM exposures


Command line usage (if any):

    usage: LocatDate.py [-h][-cp] -dir whatever exp_min exp_max

    where exp_min and exp_max are LVM exposure numbers that one wishes to locate  and

    -h prints out this help file
    -cp means not only to locate the files but to copy them to local data directory
    -dir whatever gives an alternative place to copy the data.  Not that -dir implies -cp
         even if it is not given
    -CFrame located CFrame files instead of SFrame files, which is the default

Description:

Primary routines:


Notes:

History:

231223 ksl Coding begun

'''

# # Find the calibrated data
# 
# Read my "standard" processing file and see if I can locate the associated data"
# 
# (This version is adapted from one that I had written previously for W28, but hopefully is will be useful more generally)



from glob import glob
from astropy.table import Table, join
import os
import numpy as np
import shutil
from datetime import datetime

def get_creation_time(file_path):
    try:
        return os.path.getctime(file_path)
    except FileNotFoundError:
        return None

def get_latest_file(files):
    # Get the creation times for each file
    creation_times = [(file, get_creation_time(file)) for file in files]
    
    # Filter out files where creation time couldn't be obtained
    creation_times = [(file, creation_time) for file, creation_time in creation_times if creation_time is not None]

    if not creation_times:
        return None  # No valid creation times found

    # Find the file with the maximum creation time
    latest_file, latest_creation_time = max(creation_times, key=lambda x: x[1])

    return latest_file, datetime.fromtimestamp(latest_creation_time).strftime('%Y-%m-%d %H:%M:%S')



def find_em(exp_start=3596,exp_stop=3599,file_type='C'):

    print('gotcha',exp_start,exp_stop)

    i=exp_start
    expno=[]
    xfiles=[]
    location=[]
    creation_date=[]
    nfiles=[]
    while i<=exp_stop:
        print(i)
        expno.append(i)
        xfile='lvm%cFrame-%08d.fits' % (file_type[0],i)
        xfiles.append(xfile)
        search_string='%s/**/%s' % (os.environ['SAS_BASE_DIR'],xfile)
        print('test :', search_string)
        files=glob('%s/**/%s' % (os.environ['SAS_BASE_DIR'],xfile),recursive=True)
        if len(files)==0:
            location.append('Unknown')
            creation_date.append('Unknown')
            nfiles.append(0)
            print('Could not locate %s' % xfile)
        else:
            result = get_latest_file(files)
            if result:
                latest_file, creation_time = result
                print('Located %s file created %s' % (xfile,creation_time))
                # print(f"The latest file is: {latest_file}")
                # print(f"It was created on: {creation_time}")
                location.append(latest_file)
                creation_date.append(creation_time)
                nfiles.append(len(files))
            else:
                location.append('Unknown')
                creation_date.append('Unknown')
                nfiles.append(0)
        i+=1

    # print(expno)
    xtab=Table([expno,xfiles,nfiles,creation_date,location],
        names=['Exposure','Filename','nfiles','Creation_date','Location'])

    # print(xtab)
    xtab.write('xfound_%05d_%05d.txt' % (exp_start,exp_stop),format='ascii.fixed_width_two_line',overwrite=True)

    return xtab

def get_em(xtab,destination=''):
    '''
    Move data that has been found either into the local directory 
    or a data directory
    '''

    if destination=='':
        destination='data'

    if os.path.isdir(destination)==False:
        os.mkdir(destination)
        
            
    for one in xtab:
        if one['Location']!='Unknown':        
            print('Copying %s to %s' % (one['Location'],destination))
            try:
                shutil.copy(one['Location'], destination)
            except:
                print('Error: Could not copy %s' % one['Location'])
        else:
            print('Note that %s was not located' % (one['Filename']))
    return
               
        
def steer(argv):
    exp_start=0
    exp_stop=0
    xcp=False
    destination=''
    ftype='SFrame'

    i=1
    while i<len(argv):
        if argv[i]=='-h':
               print(__doc__)
               return
        elif argv[i]=='-cp':
            xcp=True
        elif argv[i]=='-dir':
            xcp=True
            i+=1
            destination=argv[i]
        elif argv[i][0:2]=='-C':
            ftype='CFrame'
        elif argv[i][0]=='-':
            print('Error: Unknown switch',argv)
            return
        elif exp_start==0:
            exp_start=eval(argv[i])
        elif exp_stop==0:
            exp_stop=eval(argv[i])
        else:
            print(__doc__)
            print('Error:Improper command line', argv)
            return
        i+=1


    if (exp_start>exp_stop):
        print('Error: exp_stop %d less than exp_start %d, exiting' % (exp_start,exp_stop))
        return 

    if exp_stop==0:
        exp_stop=exp_start

    locate=find_em(exp_start,exp_stop)
    print(locate)
    
    if xcp:
        get_em(locate,destination)    

    locate=find_em(exp_start,exp_stop,ftype)
    if xcp:
        get_em(locate,destination)    

    return
        





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
