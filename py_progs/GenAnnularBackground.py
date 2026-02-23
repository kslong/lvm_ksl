#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Given a master source table containing source positions and angular
sizes, create a paired annular background region for each source.

Command line usage::

    GenAnnularBackground.py [-h] masterfile

Arguments:

masterfile
    ASCII table containing source names, positions, and angular sizes.
    Must include at least Major, Minor, and Color columns.

Options:

-h
    Prints this documentation.

Description:

For each source in the master table, the script enforces a minimum
angular size (xmin, default 35 arcsec) and then constructs a background
annulus whose inner radius is source_major + xspace and outer radius is
2*source_major + xspace, where xspace (default 35 arcsec) is the gap
between the source and its background annulus.  The source rows and
background annulus rows are interleaved and written to an output table.

Primary routines:

    doit

Notes:

The output file is named <masterfile_root>.ann_reg.txt.

The background annulus rows are marked with SourceBack='Back',
RegType='annulus', and Color='green'.

History:

260223 ksl Coding begun

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt
import os


from astropy.io import ascii
from astropy.table import Table,vstack
import numpy as np


def doit(masterfile='smc_snr_cotton24.txt',outfile='',xmin=35,xspace=35):

    xtab=ascii.read(masterfile)

    xtab['SourceBack']='Source'
    num=np.arange(len(xtab))+1
    xtab['No.']=num
    xtab['Color']=np.array(xtab['Color'],dtype='S10')



    xtab['Major']=np.select([xtab['Major']>xmin],[xtab['Major']],default=xmin)
    xtab['Minor']=np.select([xtab['Minor']>xmin],[xtab['Minor']],default=xmin)




    btab=xtab.copy()
    btab['RegType']='annulus'
    btab['Minor']=btab['Major']+xspace
    btab['Major']=2*btab['Major']+xspace
    btab['SourceBack']='Back'
    btab['Color']='green'



    ftab=vstack([xtab,btab])
    ftab.sort(['No.','Major'])


    ftab['Major'].format='.2f'
    ftab['Minor'].format='.2f'

    if outfile=='':
        outroot=masterfile.replace('.txt','')
        outfile='%s.ann_reg.txt' % (outroot)


    ftab.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)





def steer(argv):
    '''
    This is generally just a steering routine
    '''

    masterfile=''

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][0]=='-':
            print('Error: Could not intepret commands: ',argv)
            return
        elif masterfile=='':
            masterfile=argv[i]
        else:
            print('Error: Could not parse command line :',argv)
        i+=1

    if masterfile!='':
        doit(masterfile)
    else:
        print('Nothing to do with comand line: ',argv)







# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)        
    else:
        print (__doc__)
