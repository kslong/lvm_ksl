#!/usr/bin/env python 

"""GetFromUtah - Retrieve DRP Data from Utah

Space Telescope Science Institute

Synopsis
--------

Retrieve drp output data from Utah to a local computer.

Command Line Usage
------------------

::

    GetFromUtah.py [-h] [-cp] [-CFrame] -drp 1.1.1 mjd expstart [expstop]
    GetFromUtah.py [-h] [-cp] [-CFrame] -drp 1.1.1 filename.txt

    The second form reads an ascii table with ``mjd`` and ``expnum``
    columns (one row per exposure, not required to share a single MJD
    or be contiguous); a ``tileid`` column, if present, is ignored,
    since downloads always use a wildcard tileid.

**Options:**

-h
    Print help message

-cp
    Copy the retrieved data to a local directory

-link
    Create symbolic links in the local directory instead of copying

-CFrame
    Retrieve CFrame data instead of SFrame data

-drp VERSION
    Specify DRP processing version (default: 1.2.0)

-out OUTDIR
    Specify output directory

**Arguments:**

mjd
    MJD of the observations to retrieve

expstart
    First exposure to retrieve

expstop
    Last exposure to retrieve (optional, defaults to expstart)

filename.txt
    Alternative to ``mjd expstart [expstop]``: an ascii table with
    ``mjd``/``expnum`` columns, one row per exposure.

Description
-----------

The routine retrieves data from Utah and stores it locally
in the local redux directory.

Notes
-----

At present the routine looks for data in the 1.0.3 directories
at Utah, which is what is used for the standard processing.

When it becomes desirable to specify the version of the
drp pipeline to use, this should be straightforward to modify.

Version History
---------------

240607 ksl
    Coding begun; adapted from a routine provided by Alfredo

260714 ksl
    Removed the dead ``from lvmdrp...`` imports (rss, image, path, log,
    drpver, md) -- never actually used anywhere in this file, just
    leftover from the routine this was adapted from. This script only
    ever needed sdss_access and lvm_ksl.LocateData, so it now runs in
    any env with sdss_access (e.g. ``ksl``) instead of requiring the
    full ``lvmdrp26`` DRP environment.

260723 ksl
    sdss_access is now imported inside download_drp_product() instead
    of at module level, and steer() checks SAS_BASE_DIR exists before
    doing anything else -- sdss_access.Access eagerly os.makedirs()'s
    SAS_BASE_DIR on import, which previously crashed even -h with a raw
    traceback if e.g. an external drive backing SAS_BASE_DIR wasn't
    mounted (encountered while travelling).

260723 ksl
    Added a table-input mode (a single positional filename instead of
    mjd/expstart/expstop) alongside the existing range-input mode, so
    non-contiguous exposures or exposures spanning several MJDs can be
    requested in one call -- matching GetDAP.py's table convention.
    read_exposures_table() and resolve_tileid() are defined here and
    imported by GetDAP.py, which needs an exact (non-wildcard) tileid
    to build its download path and previously required tileid to
    already be known via its own input table; resolve_tileid() finds
    it with a read-only wildcard listing instead. check_sas_base_dir()
    (formerly private, ``_check_sas_base_dir``) is exported the same
    way, for the same reason.

"""
#!/usr/bin/env python
# coding: utf-8


import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import vstack

from lvm_ksl import LocateData


import re


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before a line consisting of "History:" (or
    "History::"/"Version History" -- whitespace/colon-insensitive), so
    -h stays short even as that section grows -- without hand-
    duplicating the Synopsis/Options text in a second string.  Anchored
    to a whole line (not a bare substring search) so it can't misfire on
    "History:" appearing mid-sentence, and returns doc unchanged if no
    such line is present.
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


def check_sas_base_dir():
    '''
    SAS_BASE_DIR must point to an existing, reachable directory --
    sdss_access's Tree eagerly os.makedirs()'s it as a side effect of
    just importing sdss_access, which raises an unhelpful traceback if
    e.g. the external drive it lives on isn't mounted (as happens while
    travelling). Caught here, up front, with a clear message instead.

    Shared with GetDAP.py (``from lvm_ksl import GetFromUtah``) rather
    than duplicated, since both scripts need the same check.
    '''
    sas_base_dir = os.environ.get('SAS_BASE_DIR')
    if not sas_base_dir:
        print('Error: SAS_BASE_DIR is not set')
        return False
    if not os.path.isdir(sas_base_dir):
        print('Error: SAS_BASE_DIR (%s) does not exist' % sas_base_dir)
        print('If it lives on an external/network drive, check that the drive is mounted.')
        print('Otherwise, point SAS_BASE_DIR at a location that exists, e.g.:')
        print('    export SAS_BASE_DIR=/some/other/path')
        return False
    return True


def read_exposures_table(filename):
    '''
    Read an ascii table with ``mjd`` and ``expnum`` columns and an
    optional ``tileid`` column, the shared table-input format for both
    this script's and GetDAP.py's table-input mode (GetFromUtah.py
    ignores tileid -- its downloads use a wildcard tileid regardless;
    GetDAP.py needs a real one, resolved via resolve_tileid() below if
    the column isn't present).

    Returns a list of (tileid_or_None, mjd, expnum) tuples, or None
    (after printing why) if the file is missing required columns.
    '''
    xtab = ascii.read(filename)
    if 'mjd' not in xtab.colnames or 'expnum' not in xtab.colnames:
        print('Error: %s must have mjd and expnum columns' % filename)
        return None

    has_tileid = 'tileid' in xtab.colnames
    exposures = []
    for row in xtab:
        tileid = int(row['tileid']) if has_tileid else None
        exposures.append((tileid, int(row['mjd']), int(row['expnum'])))
    return exposures


def resolve_tileid(drpver, mjd, expnum, kind='SFrame'):
    '''
    Read-only lookup of the real tileid for a single exposure, via a
    wildcard remote listing against the lvm_frame product (tileid='*',
    Access.set_stream() without .commit() -- nothing is downloaded,
    the same read-only-listing technique lvm_cal_sync.py uses).

    Needed by scripts like GetDAP.py that have to build an exact
    remote path themselves and can't wildcard the tileid the way
    download_drp_product() below does.

    Returns the tileid as an int, or None (after printing why) if it
    could not be resolved to exactly one value.
    '''
    from sdss_access import Access

    a = Access(release='sdsswork')
    a.remote()
    a.add('lvm_frame', drpver=drpver, tileid='*', mjd=mjd, expnum=expnum, kind=kind)
    a.set_stream()

    tileids = set()
    for task in a.stream.task:
        m = re.search(r'/(\d+)/%d/lvm' % mjd, task['location'])
        if m:
            tileids.add(int(m.group(1)))

    if len(tileids) == 1:
        return tileids.pop()
    if len(tileids) == 0:
        print('Error: could not resolve tileid for mjd=%s expnum=%s (drpver=%s) -- '
              'no matching remote files' % (mjd, expnum, drpver))
    else:
        print('Error: ambiguous tileid for mjd=%s expnum=%s -- found %s' %
              (mjd, expnum, sorted(tileids)))
    return None


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
    from sdss_access import Access

    print('what',drpver,tileid,mjd,expnum,channel,kind)
    if kind in ["Frame", "FFrame"]:
        kind = f"{kind}-{channel}" if channel in "brz" else f"{kind}-?"

    a=Access(release='sdsswork')
    q=open('Failed.txt','a')
    
    try:
        a.remote()
        a.add('lvm_frame', drpver=drpver, mjd=mjd, tileid=tileid, expnum=expnum, kind=kind)
        a.set_stream()
        a.commit()
        print(f"Downloaded product of {kind = } for {mjd = } - {expnum = }")
    except Exception as e:
        print(f"Error: failed downloading product of {kind = } for {mjd = } - {expnum = }: {e}")
        q.write(f"Error: failed downloading product of {kind = } for {mjd = } - {expnum = }: {e}")

    q.close()


def steer(argv):
    """Run the GetFromUtah routine.

    Parameters
    ----------
    argv : list
        Command line arguments including the script name

    Examples
    --------
    ::

        GetFromUtah.py -cp [-CFrame] -drp 1.2.0 [-out outdir] mjd exp_start exp_stop

    """
    drp_ver="1.2.0"
    copy=False
    xlink=False
    ftype='SFrame'
    xdest=''

    words=[]

    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(_usage_from_doc(__doc__))
            return
        elif argv[i]=='-CFrame':
            ftype='CFrame'
        elif argv[i]=='-cp':
            copy=True
        elif argv[i]=='-link':
            xlink=True
            copy=True
        elif argv[i]=='-drp':
            i+=1
            drp_ver=argv[i]
        elif argv[i]=='-out':
            i+=1
            xdest=argv[i]
        elif argv[i][0]=='-':
            print('Unknown switch: ',argv)
            return
        else:
            words.append(argv[i])
        i+=1

    if not check_sas_base_dir():
        return

    if len(words)==1:
        exposures=read_exposures_table(words[0])
        if exposures is None:
            return
        exposures=[(mjd,expnum) for (_,mjd,expnum) in exposures]
    elif len(words)>=2:
        xmjd=int(words[0])
        exp_start=int(words[1])
        exp_stop=int(words[2]) if len(words)>2 else exp_start
        exposures=[(xmjd,expnum) for expnum in range(exp_start,exp_stop+1)]
    else:
        print('Error: expected either "mjd expstart [expstop]" or a table filename with mjd/expnum columns, got:', words)
        print(_usage_from_doc(__doc__))
        return

    for xmjd,exp_now in exposures:
        download_drp_product(drpver=drp_ver, tileid='*', mjd=xmjd, expnum=exp_now, kind=ftype)

    tabs=[LocateData.find_em(expnum,expnum,ftype) for _,expnum in exposures]
    xtab=vstack(tabs) if len(tabs)>1 else tabs[0]
    print(xtab)
    if copy==True:
        LocateData.get_em(xtab,destination=xdest,link=xlink)



    # download_drp_product(drpver="1.0.3", tileid=1028683, mjd=60281, expnum=8700, kind="SFrame")




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
