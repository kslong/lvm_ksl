#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Get some basic info about SF or CF files that
have been downloaded so one can quickly see 
if the version of the data that one has may be 
out of date.


Command line usage (if any):

    usage: CheckData.py filename[s]

Description:  

Primary routines:

    doit

Notes:

    This prints out a table of some pertinent facts about each file,
    including info about the processing version
                                       
History:

250102 ksl Coding begun

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt

#!/usr/bin/env python
# coding: utf-8

from glob import glob
from astropy.io import fits
from astropy.table import Table
import numpy as np




def doit(xfiles):
    drp=[]
    commit=[]
    fluxcal=[]
    ra=[]
    dec=[]
    tile=[]
    obj=[]
    mjd=[]
    for one in xfiles:
        try:
            f=fits.open(one)
            xhead=f['PRIMARY'].header
            # print('Got ;',one)
        except:
            print('Could not open ',one)
            pass

        try:
            drp.append(xhead['drpver'])
        except:
            drp.append('Unknown')
        try:
            commit.append(xhead['COMMIT'])
        except:
            commit.append('Unknown')

        try:
            fluxcal.append(xhead['FLUXCAL'])
        except:
            fluxcal.append('Unknown')

        try:
            mjd.append(xhead['MJD'])
        except:
            mjd.append(-99)

        try:
            ra.append(xhead['POSCIRA'])
        except:
            ra.append(-37.0)

        try:
            dec.append(xhead['POSCIDE'])
        except:
            dec.append(-99.0)


        try:
            tile.append(xhead['TILE_ID'])
        except:
            tile.append('Unknown')

        try:
            obj.append(xhead['OBJECT'])
        except:
            obj.append('Unknown')

    # print(len(xfiles),len(drp),len(fluxcal),len(ra),len(dec),len(tile),len(obj))

    xtab=Table([xfiles,mjd,drp,commit,fluxcal,ra,dec,tile,obj],names=['Filename','MJD','DRP','Commit','FluxCal','RA','Dec','Tile_ID','Object'])
    xtab.sort((['Filename']))
    print(xtab)
    xtab.write('DataSum.txt',format='ascii.fixed_width_two_line',overwrite=True)
    print('Full summary in DataSum.txt')



def steer(argv):
    
    filenames=[]
    i=1
    while i<len(argv):
        if argv[i][:2]=='h':
            print(__doc__)
            return
        elif argv[i][0]=='-':
            print('Error: unknwn switch :',argv)
            return
        elif argv[i].count('.fits'):
            filenames.append(argv[i])
        i+=1

    if len(filenames)>0:
        doit(filenames)
    else:
        print('Weird - No fits files to examine: ',argv)






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    steer(sys.argv)
