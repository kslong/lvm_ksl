#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Reduce one or more consecutive LVM datasets in a single MJD


Command line usage (if any):

    usage: Reduce.py [-h] [-no_sky]  mjd exposure_start [expousure_stop]

    where

    -h  prints this help message
    -no_sky runs the quick-reduction package without sky_subtraction
    -cp causes the routine to copy the reduced frames to a directory .data 
        relative to where the program is being run
    -keep causes the routine to keep the ancilary fits files that are made
        The default is to delete these files to save disk space
    mjd  is the mjd of the exposure or exposures
    exposure_start is the first exposure to process
    exposure_stop is the last exposure to process (numbered consecuatively)

Description:  

    This routine downloadesraw  data to one's local verion of the lvm reduction
    data structure and process the data with the quick-reduction pipleline.
    One can optionally skip sky subtraction, and optionally copy the
    reduced data to a local .data directory.  The normal log fils
    produced by the pipeline can be found in a local ./xlog directory.


Primary routines:

    doit

Notes:
                                       
History:

240404 ksl Coding begun

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os


def format_number(number):
    return "{:08d}".format(number)

def doit(mjd,first_exp,last_exp,sky_sub=False,xcopy=False,clean=True):
    '''
    Process a consecutive sequence of exposures from the same MJD
    with our without sky subtraction.
    '''

    mjd='%s' % mjd
    scani=first_exp
    scanf=last_exp

    if os.path.isdir('./xlog')==False:
        os.mkdir('./xlog')


    os.environ["RSYNC_PASSWORD"] = "panoPtic-5"
    os.environ["LVMAGCAM_DIR"] = os.path.join(os.environ["SAS_BASE_DIR"], "sdsswork/data/agcam/lco/")

    # Get the necessary FITS files
    for i in range(scani, scanf + 1):
        xnumb = format_number(i)
        # Get the raw frames
        raw_frames_process = subprocess.run(["rsync", "-av", "--no-motd", f"rsync://sdss5@dtn.sdss.org/sdsswork/data/lvm/lco/{mjd}/sdR-s-*-{xnumb}.fits.gz", f"{os.environ['SAS_BASE_DIR']}/sdsswork/data/lvm/lco/{mjd}/"])
        if raw_frames_process.returncode == 0:
            print(f"Raw frames for {xnumb} successfully downloaded.")
        else:
            print(f"Failed to download raw frames for {xnumb}.")

        # Get the coadd with astrometry
        os.makedirs(f"{os.environ['SAS_BASE_DIR']}/sdsswork/data/agcam/lco/{mjd}/coadds/", exist_ok=True)
        coadd_process = subprocess.run(["rsync", "-av", "--no-motd", f"rsync://sdss5@dtn.sdss.org/sdsswork/data/agcam/lco/{mjd}/coadds/lvm.sci.coadd_s{xnumb}.fits", f"{os.environ['SAS_BASE_DIR']}/sdsswork/data/agcam/lco/{mjd}/coadds/"])
        if coadd_process.returncode == 0:
            print(f"Coadd for {xnumb} successfully downloaded.")
        else:
            print(f"Failed to download coadd for {xnumb}.")

        # Regenerate metadata
        metadata_process = subprocess.run(["drp", "metadata", "regenerate", "-m", mjd])
        if metadata_process.returncode == 0:
            print("Metadata regenerated successfully.")
        else:
            print("Failed to regenerate metadata.")

        # Quick reduction
        print("Starting reduction with log in xlog/log_{}.txt".format(i))
        with open("xlog/log_{}.txt".format(i), "w") as logfile:
            if sky_sub:
                print("Processing WITH sky subtraction")
                reduction_process = subprocess.run(["drp", "quick-reduction", "-fe", str(i)], stdout=logfile, stderr=subprocess.STDOUT)
            else:
                print("Processing WITHOUT  sky subtraction")
                reduction_process = subprocess.run(["drp", "quick-reduction", "-s","-fe", str(i)], stdout=logfile, stderr=subprocess.STDOUT)
        if reduction_process.returncode == 0:
            print(f"Reduction for {i} completed successfully.")
        else:
            print(f"Failed to complete reduction for {i}.")

        print("Finished reduction of", i)

        # Clean up
        # os.remove(os.path.join(os.environ["LVMDATA_DIR"], mjd, "ancillary", "*.fits"))

    
    if xcopy:
        print("Running LocateData and copying files")
        locate_data_process = subprocess.run(["LocateData.py", "-cp", str(scani), str(scanf)])
    else:
        print("Running LocateData w/o copying files")
        locate_data_process = subprocess.run(["LocateData.py", str(scani), str(scanf)])

    if locate_data_process.returncode == 0:
        print("LocateData.py executed successfully.")
    else:
        print("Failed to execute LocateData.py.")


def steer(argv):
    '''
    Control the flow of the program
    '''

    i=1
    mjd=-1
    first_exp=-1
    last_exp=-1
    sky_sub=True
    xcopy=False
    clean=True

    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-no_sky':
            sky_sub=False
        elif argv[i]=='-cp':
            xcopy=True 
        elif argv[i]=='-keep':
            clean=False
        elif argv[i][0]=='-':
            print('Error: Could not parse command line: ', argv)
            return
        elif mjd<0:
            mjd=int(argv[i])
        elif first_exp<0:
            first_exp=int(argv[i])
        elif last_exp<0:
            last_exp=int(argv[i])
        i+=1

    if first_exp < 0:
        print('Error: Could not parse command line: ', argv)
        return
        
    if first_exp>0 and last_exp < 0:
        last_exp=first_exp
    scani = int(sys.argv[2])

    doit(mjd,first_exp,last_exp,sky_sub,xcopy,clean)

if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
