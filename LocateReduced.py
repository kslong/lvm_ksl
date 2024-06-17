#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Locate and if desired retrieve calibrated LMV data
from series of LVM exposures on the Utah cluster
using the drpall .fits file


Command line usage (if any):

    usage: LocatReduced.py [-h][-cp] [-dir whatever] [-file_type CFrame]exp_min exp_max

    where exp_min and exp_max are LVM exposure numbers that one wishes to locate  and

    -h prints out this help file
    -cp means not only to locate the files but to copy them to local data directory
    -dir whatever gives an alternative place to copy the data.  Not that -dir implies -cp
         even if it is not given

Description:

Primary routines:


Notes:

    This routine is reads the drp_all file and finds the location of the file 
    within the filestructure associated with the reduced data.  It is 
    not clear that the routine will be used to retrieve files from that location
    a lot, but the routine was really written to show that one could find
    the file location from the drpall file

    The S/W release version is currently hardwired.

History:

240617 ksl Coding begun

'''

# # Find the calibrated data
# 
# Read my "standard" processing file and see if I can locate the associated data"
# 
# (This version is adapted from one that I had written previously for W28, but hopefully is will be useful more generally)



from glob import glob
from astropy.table import Table, join
from astropy.io import fits
import os
import numpy as np
import shutil
from datetime import datetime



def read_drpall(drp_ver='1.0.3'):
    BASEDIR='/uufs/chpc.utah.edu/common/home/sdss51/sdsswork/lvm/spectro/redux/%s/' % (drp_ver)
    DRPFILE='drpall-%s.fits' % (drp_ver)
    xfile='%s/%s' % (BASEDIR,DRPFILE)
    try:
        drpall=fits.open(xfile)
    except:
        print('Error: Could not locate : ', xfile)
        return []

    drp_tab=Table(drpall[1].data)
    return  drp_tab

def find_em(ztab,exp_start,exp_stop):
    xtab=ztab[ztab['expnum']>=exp_start]
    xtab=xtab[xtab['expnum']<=exp_stop]
    return xtab

    
        
XTOP='/uufs/chpc.utah.edu/common/home/sdss51/'

def get_em(xtab,destination='',file_type='SFrame'):
    '''
    Move data that has been found either into the local directory 
    or a data directory
    '''

    nfiles_max=100
    if len(xtab)>nfiles_max:
        print('Error: Trying to retrieve %d files, which is more than %d, and likely a mistake' % (len(xtab),nfiles_max))
        return

    if destination=='':
        destination='data'

    if os.path.isdir(destination)==False:
        os.mkdir(destination)
        
            
    for one in xtab:
        name='%s/%s' %  (XTOP,one['location'])
        if file_type != 'SFrame':
            name=name.replace('SFrame',file_type)
        print('Copying %s to %s' % (name,destination))
        try:
            shutil.copy(name, destination)
        except:
            print('Error: Could not copy %s' % name)
    return
               
        
def steer(argv):
    exp_start=0
    exp_stop=0
    xcp=False
    destination=''
    file_type='SFrame'

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
        elif argv[i]=='-type':
            i+=1
            file_type=argv[i]
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


    drp_tab=read_drpall()

    locate=find_em(drp_tab,exp_start,exp_stop)
    print(locate)
    
    if xcp:
        get_em(locate,destination,file_type)    

    return
        





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
