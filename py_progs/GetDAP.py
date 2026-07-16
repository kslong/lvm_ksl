#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Retrieve DAP files from Utah


Command line usage (if any):

    usage: GetDAP.py [-h] filename.txt

    where filename is an astropy table containg suffiencient information to locate the DAP file

Description:  

Primary routines:

    doit

Notes:
                                       
History::

    250204 ksl Coding begun
    260714 ksl get_dap() rewritten to use sdss_access (Access.add_file)
    instead of a hand-rolled rsync subprocess call against
    ~/.sdss_rsync_password -- dtn.sdss.org now requires 2FA for that
    rsync auth path, which sdss_access sidesteps via .netrc. sdss_access
    has no wildcard/glob download, so this now fetches the DAP fits.gz
    by its known dap-<config>-<expnum>.dap.fits.gz name (DAP_CONFIG =
    'rsp108-sn20', matching the convention already assumed by
    DAP2tab.py/DAPGauss2tab.py) instead of the old `*fits.gz` glob --
    confirmed via a real download this is ~5s vs ~4min for listing the
    whole per-exposure directory. Dropped unused sys/np imports.

'''

from astropy.io import ascii
import os
import shutil
from sdss_access import Access


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


DAP_TOP='sdsswork/lvm/spectro/analysis'

# DAP output filenames follow this fixed dap-<config>-<expnum>.dap.fits.gz
# convention (also assumed by DAP2tab.py/DAPGauss2tab.py) rather than being
# discovered via a directory listing, since sdss_access has no wildcard/glob
# download and listing a whole exposure's DAP directory is much slower.
DAP_CONFIG='rsp108-sn20'

def get_dap(drpver, tileid, mjd, expnum):
    '''
    Uses sdss_access (HTTPS/rsync + .netrc) rather than a raw rsync
    subprocess, since dtn.sdss.org now requires 2FA for interactive
    rsync password auth.
    '''

    xtile='%07d' % tileid
    xtile='%sXX' % xtile[:4]
    print(xtile)
    remote_name = 'dap-%s-%08d.dap.fits.gz' % (DAP_CONFIG, expnum)
    local_full = os.path.join(os.environ['SAS_BASE_DIR'], DAP_TOP, str(drpver), xtile,
                               str(tileid), str(mjd), '%08d' % expnum, remote_name)
    print(local_full)

    if os.path.isdir('DAP')==False:
        os.makedirs('DAP')

    a = Access(release='sdsswork')
    try:
        a.remote()
        a.add_file(local_full)
        a.set_stream()
        a.commit()
        shutil.copy(local_full, 'DAP/')
        print(f"%s successfully downloaded." % remote_name)
    except Exception as e:
        print(f"Failed to download %s: %s" % (remote_name, e))


def steer(argv):

    drpver='1.1.1'
    
    filename=''

    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(_usage_from_doc(__doc__))
            return
        elif  argv[i][0]=='-':
            print('Unknown option :', argv)
            return
        else:
            filename=argv[i]
        i+=1

    xtab=ascii.read(filename)
    for one_row in xtab:
        get_dap(drpver,one_row['tileid'],one_row['mjd'],one_row['expnum'])

    return

                            






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
