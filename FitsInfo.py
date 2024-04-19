#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Just look at the structure of a fits file


Command line usage (if any):

    usage: FitsInfo filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

240418 ksl Coding begun

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt



def doit(filename='foo.fits'):
    '''
    Do something magnificent

    Description:

    Notes:

    History:


    '''

    try:
        x=fits.open(filename)
    except:
        print('Error: could not locate ',filename)

    x.info()








# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        doit(sys.argv[1])
    else:
        print (__doc__)
