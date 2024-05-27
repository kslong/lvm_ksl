#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Reduce one or more consecutive LVM datasets in a single MJD


Command line usage (if any):

    usage: Reduce.py [-h] [-keep] [-cp}mjd exposure_start [expousure_stop]

    where

    -h  prints this help message
    -keep retains the ancillary fiels which are otherwise deleteed
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
240526 ksl Adapt to new version of the DRP

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import timeit
import time
import multiprocessing
multiprocessing.set_start_method("spawn",force=True)
import os
import sky_plot
import re
from collections import Counter


def clean_and_count_lines_with_keywords(filename):
    # Regular expression to remove non-ASCII characters
    non_ascii_regex = re.compile(r'[^\x00-\x7F]+')

    # Dictionary to store cleaned lines and their counts
    line_counts = Counter()

    # Read the file line by line
    with open(filename, 'r', encoding='utf-8') as file:
        for line in file:
            # Check if the line contains 'WARNING' or 'ERROR'
            if 'WARNING' in line or 'ERROR' in line:
                # Remove non-ASCII characters
                cleaned_line = non_ascii_regex.sub('', line)

                # Count the cleaned line
                line_counts[cleaned_line.strip()] += 1

    # Display results sorted by frequency (most common first)
    print('### WARNING/ERROR Summary ###')
    for line, count in line_counts.most_common():
        print(f"{count}: {line}")


def format_number(number):
    return "{:08d}".format(number)

def process_one(mjd,i,clean):

        # Regenerate metadata
        # metadata_process = subprocess.run(["drp", "metadata", "regenerate", "-m", mjd])
        # if metadata_process.returncode == 0:
        #     print("Metadata regenerated successfully.")
        # else:
        #     print("Failed to regenerate metadata.")

        print("Starting reduction with log in xlog/log_{}.txt".format(i))
        with open("xlog/log_{}.txt".format(i), "w") as logfile:
            xcommand=["drp", "run", "-m",str(mjd),"-e", str(i)]
            if clean==True:
                xcommand=["drp", "run", "-c","-m",str(mjd),"-e", str(i)]

            print("Begin processing of ",i,"with command: ",xcommand)
            reduction_process = subprocess.run(xcommand, stdout=logfile, stderr=subprocess.STDOUT)
        if reduction_process.returncode == 0:
            print(f"Reduction for {i} completed successfully.")
        else:
            print(f"FAILED to complete reduction for {i} on {mjd}, check log for errors.")

        print('Type for logfile', type(logfile))
        print('Path ',os.path.isfile(logfile.name))
        print('name', logfile.name)

        clean_and_count_lines_with_keywords(logfile.name)

        print("Finished reduction of", i)

        return  reduction_process.returncode


# Function to run rsync command with ignore-existing and dry-run options
def run_forced_dry_rsync(mjd,xnumb,verbose=False):
    # Run rsync command with dry-run and ignore-existing options
    rsync_process = subprocess.Popen(
        ["rsync", "-avn", "--no-motd", "--ignore-existing",
         f"rsync://sdss5@dtn.sdss.org/sdsswork/data/lvm/lco/{mjd}/sdR-s-*-{xnumb}.fits.gz",
         f"{os.environ['SAS_BASE_DIR']}/sdsswork/data/lvm/lco/{mjd}/"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
    )

    # Capture stdout and stderr
    stdout, stderr = rsync_process.communicate()

    # Check if rsync command was successful
    if verbose:
        if rsync_process.returncode == 0:
            print(f"Files for {xnumb} would be successfully downloaded if forced.")
        else:
            print("Rsync forced dry-run failed.")

        # Print stdout and stderr
        print("Standard Output:")
        print(stdout)
        print("Standard Error:")
        print(stderr)

    return rsync_process.returncode


def get_data(mjd,i):
    '''
    qmjd is a string
    '''

    os.environ["RSYNC_PASSWORD"] = "panoPtic-5"
    os.environ["LVMAGCAM_DIR"] = os.path.join(os.environ["SAS_BASE_DIR"], "sdsswork/data/agcam/lco/")
    mjd='%s' % mjd
    xmjd='%d' % (int(mjd)+1)
    xnumb = format_number(i)

    if run_forced_dry_rsync(mjd,xnumb)==0:
        print('All is OK with %s so proceeding' % mjd)
        qmjd=mjd
    elif run_forced_dry_rsync(xmjd,xnumb)==0:
        print('Failed on orginal %s, but succeeded with  %s' % (mjd,xmjd))
        qmjd=xmjd
    else:
        print('Failed with both %s and %s so returning' % (mjd,xmjd))
        return 'Failed'



    # Get the raw frames
    raw_frames_process = subprocess.run(["rsync", "-av", "--no-motd", f"rsync://sdss5@dtn.sdss.org/sdsswork/data/lvm/lco/{qmjd}/sdR-s-*-{xnumb}.fits.gz", f"{os.environ['SAS_BASE_DIR']}/sdsswork/data/lvm/lco/{qmjd}/"])
    if raw_frames_process.returncode == 0:
        print(f"Raw frames for {xnumb} successfully downloaded.")
    else:
        print(f"Failed to download raw frames for {xnumb}.")

    # Get the coadd with astrometry
    os.makedirs(f"{os.environ['SAS_BASE_DIR']}/sdsswork/data/agcam/lco/{qmjd}/coadds/", exist_ok=True)
    coadd_process = subprocess.run(["rsync", "-av", "--no-motd", f"rsync://sdss5@dtn.sdss.org/sdsswork/data/agcam/lco/{qmjd}/coadds/lvm.sci.coadd_s{xnumb}.fits", f"{os.environ['SAS_BASE_DIR']}/sdsswork/data/agcam/lco/{mjd}/coadds/"])
    if coadd_process.returncode == 0:
        print(f"Coadd for {xnumb} successfully downloaded.")
    else:
        print(f"Failed to download coadd for {xnumb}.")

    return qmjd


def do_one(mjd,exp,clean=True):
    '''
    mjd is a string, as is qmjd
    '''

    qmjd=get_data(mjd,exp)

    if qmjd=='Failed':
        return 1

    if mjd!=qmjd:
        print('Although these data were taken on %s, the are in SMJD %s' % (mjd,qmjd))
        mjd=qmjd

    process_one(mjd,exp,clean)

    return 0
    
def get_no_jobs(jobs):
    '''
    Check how many jobs are running
    '''
    njobs=0
    for one in jobs:
        if one.is_alive():
            njobs+=1
    return njobs

def do_many(mjd,first_exp,last_exp,clean=True,xcopy=True,nproc=8):    
    '''
    A routine to run the lvmdrp in parallel
    '''

    if int(mjd) < 60177:
        print('WARNING: THESE DATA ARE UNLIKELY BE CALIBRATABLE WITHOUT SPECIAL EFFORT, AS THEY WERE OBTAINED BEFORE MJD 60177')

    mjd='%s' % mjd
    scani=first_exp
    scanf=last_exp

    if os.path.isdir('./xlog')==False:
        os.mkdir('./xlog')


    start_time = timeit.default_timer()

    jobs=[]
    for i in range(scani, scanf + 1):

        p=multiprocessing.Process(target=do_one,args=[mjd,i,clean])
        jobs.append(p)

    i=0
    while i<nproc and i<len(jobs):
        t = time.localtime()
        one=jobs[i]
        one.start()
        i+=1

    njobs=get_no_jobs(jobs)

    while i<len(jobs):
        time.sleep(2)
        njobs=get_no_jobs(jobs)

        while njobs<nproc and i<len(jobs):
            t = time.localtime()
            one=jobs[i]
            one.start()
            njobs+=1
            i+=1

    p.join()
    p.close()

    elapsed = timeit.default_timer() - start_time
    print('Completed multiprocessing of  %s exposures  %d - %d ' % (mjd,first_exp,last_exp))

    return



def doit(mjd,first_exp,last_exp,clean=True,xcopy=False):
    '''
    Process a consecutive sequence of exposures from the same MJD
    with our without sky subtraction.
    '''

    if int(mjd) < 60177:
        print('WARNING: THESE DATA ARE UNLIKELY BE CALIBRATABLE WITHOUT SPECIAL EFFORT, AS THEY WERE OBTAINED BEFORE MJD 60177')

    mjd='%s' % mjd
    scani=first_exp
    scanf=last_exp

    if os.path.isdir('./xlog')==False:
        os.mkdir('./xlog')



    # Get the necessary FITS files
    for i in range(scani, scanf + 1):
        do_one(mjd,i,clean)
    
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
    clean=True
    xcopy=False
    clean=True
    nproc=1

    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        if argv[i]=='-np':
            i+=1
            nproc=int(argv[i])
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

    if nproc==1:
        doit(mjd,first_exp,last_exp,clean,xcopy)
    else:
        do_many(mjd,first_exp,last_exp,clean,xcopy,nproc)

if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
