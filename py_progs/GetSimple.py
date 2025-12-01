#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Retrieve spefic files from Utah based on their exact location

Command line usage (if any):

    usage: GetSimple.py [-h] [-outdir whatever] filename.txt

    where filename is an astropy table containingg the exact name of  a file to retrive in the column location

Description:  

Primary routines:

    steer
    get_file

Notes:
                                       
History:

250204 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy as np
import os
import subprocess




def get_file(xfile,path='data'):
    '''
    xfile is the name of a spectrum
    '''

    os.environ["RSYNC_PASSWORD"] = "panoPtic-5"
    print(xfile)
    if xfile[0]!='/':
        xfile='/%s' % xfile

    if os.path.isdir(path)==False:
        os.makedirs(path)

    # Get the raw frames
    raw_frames_process = subprocess.run(["rsync", "-av", "--no-motd", f"rsync://sdss5@dtn.sdss.org%s" % xfile, "%s/" % path])
    if raw_frames_process.returncode == 0:
        print(f"%s successfully downloaded." % xfile)
    else:
        print(f"Failed to download %s." %  xfile)


def steer(argv):


    
    path='data/'
    filename=''

    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-outdir':
            i+=1
            path=argv[i]
        elif  argv[i][0]=='-':
            print('Unknown option :', argv)
            return
        else:
            filename=argv[i]
        i+=1
    xtab=ascii.read(filename)
    for one_row in xtab:
        get_file(one_row['location'],path=path)

    return

                            






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
