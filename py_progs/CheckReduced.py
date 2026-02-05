#!/usr/bin/env python
# coding: utf-8


'''
                    Space Telescope Science Institute

Synopsis:  

Check whether a set of reductions succeeded and whether 
all of the reduced files are from the same commit


Command line usage (if any):

    usage: CheckReduced.py [-h] [-d whatever] [-l whatever]

Description:  

    The routine is intended to be used in conjunction with Reduce.py
    which normally copies reduced files to a subdirectory data,
    and a copy of the drp run log file to a subdirectory xlog.

    With none of the optional switches the routine assumes processed
    data files are in the data directory, and the log files are
    in the xlog. Options: -d whatever specifies where the SFrame files
    to be looked at exist; -l whatever specifies an alternative directory
    for the log files created with Reduce.py; -h prints this documentation
    and quits.

    The results are written to the screen but also stored for reference
    in commits.txt for the commit information and problems.txt for the errors.

Primary routines:

    doit

Notes:

    Unlike most other routines, this one is normally run witout arguments.
                                       
History:

251224 ksl Coding begun

'''

import os
import sys
from astropy.io import ascii,fits
import numpy as np
from glob import glob
from astropy.table import Table
import re



def check_commits(filenames,xout='commits.txt'):

    commit=[]
    exposure=[]
    drpver=[]
    for one_file in filenames:
        x=fits.open(one_file)
        xcom=x[0].header['COMMIT']
        if xcom==None:
            xcom='Unknown'
        commit.append(xcom)
        drpver.append(x[0].header['DRPVER'])
        exposure.append(x[0].header['EXPOSURE'])

    xtab=Table([filenames,exposure,drpver,commit],names=['Filename','Exposure','DRP_ver','Commit'])
    xtab.sort('Exposure')

    # print(commit)

    xver,counts=np.unique(xtab['Commit'],return_counts=True)

    i=0
    while i<len(xver):
        print('There were %3d files with commit %s' % (counts[i],xver[i]))
        i+=1

    print('Writing commit info to ',xout)
    xtab.write(xout,format='ascii.fixed_width_two_line',overwrite=True)



def check_logs(logs,xout='problems.txt'):


    succeeded=[]
    failed=[]
    all_errors=[]

    for one_log in logs:
        xx=open(one_log)
        lines=xx.readlines()
        errors=[]
        for line in lines:
            if line.count('ERROR'):
                errors.append(line)
        if len(errors)>0:
            failed.append(one_log)
            all_errors.append(errors)
        else:
            succeeded.append(one_log)


    print('Of %d logs considered, %d had no errors, and %d had errors' % (len(logs),len(succeeded),len(failed)))

    print('The following had no errors')
    for one in succeeded:
        print(one)


    print('\nBut the following had issues')
    for one in failed:
        print(one)

    print('\nHere are the details on the errors')
    print('Also writing them to ',xout)

    ansi_escape = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')
    f=open(xout,'w')
    f.write('The number of logs with reported problems is %d' % len(failed))
    i=0
    while i<len(failed):

        print(failed[i])
        f.write(failed[i])
        f.write('\n')
        xerr=all_errors[i]
        for one_line in xerr:
            print(one_line.strip())
            clean=ansi_escape.sub('', one_line)
            f.write(clean)
        print('\n')
        f.write('\n')
        i+=1
    f.close()





def doit(data_dir='data',log_dir='xlog'):
    '''
    Do something magnificent

    Description:

    Notes:

    History:


    '''
    os.getcwd()
    filenames=glob('%s/*SF*.fits' % data_dir)
    check_commits(filenames)
    logs=glob('%s/log*txt' % log_dir)
    check_logs(logs)
    
    return


def steer(argv):
    '''
    This is generally just a steering routine
    '''

    data_dir='data'
    log_dir='xlog'

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc_)
            return
        elif argv[i][:2]=='-d':
            i+=1
            data_dir=argv[i]
        elif argv[i][:2]=='-l':
            i+=1
            log_dir_dir=argv[i]
        elif argv[i][0]=='-':
            print('Error: Could not intepret commands: ',argv)
            return
        i+=1


    doit(data_dir,log_dir)




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    steer(sys.argv)
