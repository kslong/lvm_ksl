#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

A simple routine to generate a table with RA's and Dec's
of LVM science fiber


Command line usage (if any):

    usage: fib2radec.py [-g guider_file] [-cal calibrated data_file] expno

    where

    -g guider_file is an optional parameter to point to a specific guider file
    -cal calibrated file is an optional parameter to point to a specific calibrated data file
    expno is an (optional actually) paremeter to pint to an specific lvm explsoure

Description:  

    The routine searches for a calibrated rss spectrum and the associate guider file for
    an LVM exposure and generates a table that contains the RA's and Dec's of each fiber.

    The routine searches for these files in the follow order, the current directory of 
    any of its subdirectories and then the standard locations for these files based
    of the defined environment variables LVM_MASTER_DIR and LVMAGCAM_DIR.



Primary routines:


Notes:

    This routine was derived from Kathryn K's routine visualize_genearl.py 
                                       
History:

231128 ksl Coding begun

'''

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
import yaml
from astropy.wcs import wcs
from astropy.table import Table
import os
from astropy.wcs import WCS
from glob import glob


# platescale to convert fiber positions from mm to arcsec
PLATESCALE = 112.36748321030637

# directory where reduced cFrames are stored
DIR_redux = '/data/LVM/sdsswork/lvm/spectro/redux/0.1.1.dev0/1111/'
DIR_agcam = '/data/LVM/sdsswork/data/agcam/lco/'

def rotate(xx,yy,angle):
    # rotate x and y cartesian coordinates by angle (in degrees)
    # about the point (0,0)
    theta = -1.*angle * np.pi / 180. # in radians
    xx1 = np.cos(theta) * xx - np.sin(theta) * yy
    yy1 = np.sin(theta) * xx + np.cos(theta) * yy
    return xx1, yy1

def make_radec(xx0,yy0,ra,dec,pa):
    xx, yy = rotate(xx0,yy0,-pa)
    ra_fib = ra + xx*PLATESCALE/3600./np.cos(dec*np.pi/180.) 
    dec_fib = dec - yy*PLATESCALE/3600. 
    return ra_fib, dec_fib



def get_ra_dec(filename,guider_file,outname=''):
    
    # read in the rss file
    try:
        rsshdu = fits.open(filename)
    except:
        print('Error: Could not find %s' % filename)
        return

    hdr = rsshdu[0].header

    r1 = rsshdu[1].data
    r1_hdr = rsshdu[1].header
    r1_err = rsshdu[2]

    wave=rsshdu[4].data 

    tab = Table(rsshdu['SLITMAP'].data)
    tab=tab[tab['targettype']=='science']
    rsshdu.close()


    # get ra/dec measured from coadd guiders?
    agcam_coadd =guider_file 
    if os.path.isfile(agcam_coadd):
        agcam_hdu = fits.open(agcam_coadd)
        agcam_hdr = agcam_hdu[1].header
        w = WCS(agcam_hdr)
        cen = w.pixel_to_world(2500,1000)
        racen = cen.ra.deg  #agcam_hdr['RAMEAS']
        deccen = cen.dec.deg #agcam_hdr['DECMEAS']
        print('hello ', agcam_hdr['PAMEAS'])
        pa = 180 - agcam_hdr['PAMEAS'] 
#        print(pa)
        agcam_hdu.close()
        xguide=True
    else:
        if guider_file=='':
            print('No guider file provided, using science header')
        else:
            print('No guider file present at %s, using science header' % guider_file)
        racen = hdr['TESCIRA']
        deccen = hdr['TESCIDE']
        pa = hdr['POSCIPA']# *(-1.)
        xguide=False
    
    ra_fib, dec_fib = make_radec(tab['xpmm'], tab['ypmm'], racen, deccen, pa)

    outtab=Table([tab['fiberid'],ra_fib,dec_fib], names=['fiberid','RA','Dec'])

    outtab['RA'].format='.6f'
    outtab['Dec'].format='.6f'

    if outname=='':
        words=filename.split('/')
        outname=words[-1].replace('.fits','.txt')
        outname=outname.replace('lvmCFrame','Fib2radec')
    outtab.write(outname,format='ascii.fixed_width_two_line',overwrite=True)
    return xguide,racen,deccen,pa,outtab



def locate_file(file_name):

    if len(file_name)==0:
        print('Error: cannot look for file without a name')
        return None

    print('Looking for file_name %s' % file_name)
    # Check if the file exists in the current directory
    current_directory = os.getcwd()
    current_file_path = os.path.join(current_directory, file_name)
    if os.path.isfile(current_file_path):
        return current_file_path

    # Check in subdirectories of the current directory
    for root, dirs, files in os.walk(current_directory):
        if file_name in files:
            return os.path.join(root, file_name)

    # Check in a directory specified by an environment variable and its subdirectories
    if file_name.count('CFrame'):
        env_variable_name = 'LVMDATA_DIR'
    elif file_name.count('lvm.sci.coadd'):
        env_variable_name = 'LVMAGCAM_DIR'


    if env_variable_name in os.environ:
        env_variable_value = os.environ[env_variable_name]
        env_variable_path = os.path.join(env_variable_value, file_name)
        if os.path.isfile(env_variable_path):
            return env_variable_path

        for root, dirs, files in os.walk(env_variable_value):
            if file_name in files:
                return os.path.join(root, file_name)

    # If the file is not found in any of the locations, return None
    return None


def steer(argv):
    '''
    This is just a steering routine
    '''
    calibrated_file=''
    guider_file=''
    exposure=-99


    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        if argv[i]=='-g':
            i+=1
            guider_file=argv[i]
        if argv[i]=='-cal':
            i+=1
            calibrated_file=argv[i]
        elif argv[i][0]=='-':
            print('Error: Unknown switch %s ' % argv[i])
            return
        elif exposure==-99:
            exposure=int(argv[i])
        i+=1

    if exposure!=-99 and calibrated_file=='':
        # Construct the name of the calibrated file
        calibrated_file='lvmCFrame-%08d.fits' % exposure
        guider_file='lvm.sci.coadd_s%08d.fits' % exposure

    # Now find the files.  Our hiearchy is to look for them first
    # in the current directory, then in subdirectories of the current
    # directory, and finally using enviroment variables

    print('test cal %s guider %s exposure %d' % (calibrated_file,guider_file,exposure))

    xcal=locate_file(calibrated_file)
    if xcal==None:
        print('Could not locate calibrated file: %s' % calibrated_file)
        return

    if guider_file!='':
        xguide=locate_file(guider_file)
        if xguide==None:
            print('Could not locate guider file: %s' % guider_file)
            xguide=''


    get_ra_dec(xcal,xguide,outname='')


    return
        



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
