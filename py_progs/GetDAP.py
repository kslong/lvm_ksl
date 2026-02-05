#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Retrieve DAP files from Utah


Command line usage (if any):

    usage: GetDAP.py [-h] filename.txt

    where filename is an astropy table containg suffiencient information to locate the DAP file

Description:  

Primary routines:

    doit

Notes:
                                       
History:

250204 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy as np
import os
import subprocess




DAP_TOP='/sdsswork/lvm/spectro/analysis/'
def get_dap(drpver, tileid, mjd, expnum):
    '''
    qmjd is a string
    '''

    os.environ["RSYNC_PASSWORD"] = "panoPtic-5"
    xtile='%07d' % tileid
    xtile='%sXX' % xtile[:4]
    print(xtile)
    xfile='%s/%s/%s/%s/%d/%08d/*fits.gz' % (DAP_TOP,drpver,xtile,tileid,mjd,expnum)
    print(xfile)

    if os.path.isdir('DAP')==False:
        os.makedirs('DAP')

    # Get the raw frames
    raw_frames_process = subprocess.run(["rsync", "-av", "--no-motd", f"rsync://sdss5@dtn.sdss.org%s" % xfile, f"./DAP/"])
    if raw_frames_process.returncode == 0:
        print(f"%s successfully downloaded." % xfile)
    else:
        print(f"Failed to download %s." %  xfile)


def steer(argv):

    drpver='1.1.1'
    
    filename=''

    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif  argv[i][0]=='-':
            print('Unknown option :', argv)
            return
        else:
            filename=argv[i]
        i+=1

    xtab=ascii.read(filename)
    for one_row in xtab:
        get_dap(drpver,one_row['tileid'],one_row['mjd'],one_row['expnum'])

    return

                            






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
