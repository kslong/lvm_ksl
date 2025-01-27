#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Retrieve drp output data from Utah to a local computer


Command line usage (if any):

    usage: GetFromUtah.py [-h] [-cp] [-CFrame]  -drp 1.1.1 mjd expstart [expstop]

    where:

    -h prints out this help message
    -cp copies the retrieved data to a local directory
    -CFrame causes the CFrame data to be retrieved instead of the 
        SFrame data
    -drp whatever wistchedst the default location to a different processing version
    mjd is the mjd of the observations one wants to retrieve
    expstart is the first exposure to retrieve
    expstart, if given means to retrieve exposures from expstart to
        expstop

Description:  

    The routine retrieves data from Utah and stores it locally
    in the local redux directory.


Primary routines:


Notes:

    At present the routine looks for data in the 1.0.3 directories
    at Utah, which is what is used for the standard processing.

    When it becomes desirable to specify the version of the
    drp pipeline to use, this should be straightforward to
    modify.
                                       
History:

240607 ksl Coding begun; adapted from a routine provided
    by Alfredo

'''
#!/usr/bin/env python
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt

from lvmdrp.core import rss
from lvmdrp.core import image


from lvmdrp import path, log, __version__ as drpver
from lvmdrp.core import rss
import lvmdrp.utils.metadata as md



from sdss_access  import Access

from lvm_ksl import LocateData 


def download_drp_product(drpver, tileid, mjd, expnum, channel=None, kind="SFrame"):
    """Download LVM DRP products: lvmFrame, lvmFFrame, lvmCFrame, lvmSFrame

    Parameters
    ----------
    drpver : str
        DRP version (e.g., '1.0.3')
    tileid : int, str
        Tile ID (e.g., 11111)
    mjd : int, str
        MJD of the observation (e.g., 60275)
    expnum : int, str
        exposure number (e.g., 8432)
    channel : str, optional
        spectrograph channel (e.g., 'b', 'r', 'z')
    kind : str, optional
        LVM DRP product kind/species ('CFrame', 'SFrame'), by default 'SFrame'
    """
    if kind in ["Frame", "FFrame"]:
        kind = f"{kind}-{channel}" if channel in "brz" else f"{kind}-?"
    
    a = Access(release='sdsswork')
    try:
        a.remote()
        a.add('lvm_frame', drpver=drpver, mjd=mjd, tileid=tileid, expnum=expnum, kind=kind)
        a.set_stream()
        a.commit()
        print(f"Downloaded product of {kind = } for {mjd = } - {expnum = }")
    except Exception as e:
        print(f"Error: failed downloading product of {kind = } for {mjd = } - {expnum = }: {e}")


def steer(argv):
    '''
    Run the routine

    GetFromUtah -cp mjd exp_start exp_stop 

    '''
    drp_ver="1.1.0"
    copy=False
    xtile=1028683
    xtile='*'
    xmjd=60281
    xexp=8700
    ftype='SFrame'
    # find_type='S'
    xdest=''

    xmjd=''
    exp_start=-1
    exp_stop=-1

    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-CFrame':
            ftype='CFrame'
        elif argv[i]=='-cp':
            copy=True
        elif argv[i]=='-drp':
            i+=1
            drp_ver=argv[i]
        elif argv[i][0]=='-':
            print('Unknown switch: ',argv)
            return
        elif xmjd=='':
            xmjd=int(argv[i])
        elif exp_start==-1:
            exp_start=int(argv[i])
        elif exp_stop==-1:
            exp_stop=int(argv[i])
        else:
            print('Error: Too many arguments: ',argv)
            return
        i+=1

    if exp_stop==-1:
        exp_stop=exp_start


    exp_now=exp_start
    while exp_now <= exp_stop:
        download_drp_product(drpver=drp_ver, tileid=xtile, mjd=xmjd, expnum=exp_now, kind=ftype)
        exp_now+=1

    xtab=LocateData.find_em(exp_start,exp_stop,ftype)
    print(xtab)
    if copy==True:
        LocateData.get_em(xtab,destination=xdest)



    # download_drp_product(drpver="1.0.3", tileid=1028683, mjd=60281, expnum=8700, kind="SFrame")




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
