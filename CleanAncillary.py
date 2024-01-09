#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Remove the Ancillary File directories that are created in
the process of reducing lvm data


Command line usage (if any):

    usage: CleanAncillary.py [-rm]      

Description:  

    Without the -rm the routine similary lists the 
    directories names ancillary.

    With -rm, the directories are actually removed

Primary routines:

    clean

Notes:

    The ancillary data files contain about
    3/4 of the data volume from the lvmdrp
                                       
History:

240109 ksl Coding begun

'''

from glob import glob
import os
import shutil

def remove_directory(directory_path):
    try:
        shutil.rmtree(directory_path)
        print(f"Directory '{directory_path}' successfully removed.")
    except OSError as e:
        print(f"Error removing directory '{directory_path}': {e}")



def clean(dir_name='ancillary',remove=False):
    xxdir='%s/**/%s' % (os.getcwd(),dir_name)
    names=glob(xxdir,recursive=True)

    if len(names)==0:
        print('There are no directories to remove')
        return

    for one in names:
        if remove:
            remove_directory(one)
        else:
            print(one)


def steer(argv):
    xremove=False

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
        elif argv[i]=='-rm':
            xremove=True
        else:
            print('Error: Incorrect command line: ',argv)
            return
        i+=1

    clean(dir_name='ancillary',remove=xremove)

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>0:
        steer(sys.argv)
    else:
        print (__doc__)



