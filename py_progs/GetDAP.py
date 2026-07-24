#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Retrieve DAP files from Utah


Command line usage (if any):

    usage: GetDAP.py [-h] [-no_cp] [-drp VERSION] filename.txt
    usage: GetDAP.py [-h] [-no_cp] [-drp VERSION] mjd expstart [expstop]

    The first form reads an ascii table with mjd/expnum columns (an optional
    tileid column, if present, is used directly and skips the lookup below).
    The second form retrieves a range of exposures from a single MJD directly
    from the command line, matching GetFromUtah.py's convention; since this
    script needs an exact tileid (unlike GetFromUtah.py, which can wildcard
    it), one is looked up automatically for each exposure via a read-only
    remote listing (see GetFromUtah.resolve_tileid) -- this is an extra
    network round-trip per exposure, so the table form is faster if you
    already know the tileids.
    -drp VERSION sets the DRP version used to locate the DAP files (default: 1.2.1).
    -no_cp skips copying each retrieved DAP file into a local ./DAP directory,
    leaving it only in the SAS_BASE_DIR-mirrored location that sdss_access
    downloads it to (by default it is copied into ./DAP as well).

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
    260723 ksl sdss_access is now imported inside get_dap() instead of
    at module level, and steer() checks SAS_BASE_DIR exists before
    doing anything else -- sdss_access.Access eagerly os.makedirs()'s
    SAS_BASE_DIR on import, which previously crashed even -h with a raw
    traceback if e.g. an external drive backing SAS_BASE_DIR wasn't
    mounted (encountered while travelling).
    260723 ksl Default drpver changed from 1.1.1 to 1.2.1, and it is
    now settable from the command line with -drp, matching
    GetFromUtah.py's convention. Briefly made the ./DAP copy opt-in
    via -cp (matching GetFromUtah.py's -cp), then reverted: copying
    into local ./DAP stays the default behavior, since other tools
    (e.g. DAP2tab.py) assume it is there; -no_cp skips it instead.
    260723 ksl Now also accepts "mjd expstart [expstop]" on the command
    line as an alternative to the table-input mode, matching
    GetFromUtah.py. Since this script (unlike GetFromUtah.py) needs an
    exact tileid to build its download path, tileid is looked up via
    the new GetFromUtah.resolve_tileid() (a read-only remote listing)
    whenever it isn't already supplied by a table's tileid column.
    check_sas_base_dir(), resolve_tileid(), and read_exposures_table()
    are now imported from GetFromUtah.py rather than duplicated.

'''

from astropy.io import ascii
import os
import shutil

from lvm_ksl import GetFromUtah

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

def get_dap(drpver, tileid, mjd, expnum, copy=True):
    '''
    Uses sdss_access (HTTPS/rsync + .netrc) rather than a raw rsync
    subprocess, since dtn.sdss.org now requires 2FA for interactive
    rsync password auth. sdss_access is imported here rather than at
    module level so that just importing/parsing-args-for this script
    (e.g. -h) doesn't touch it -- see GetFromUtah.check_sas_base_dir().

    a.commit() always leaves the file in the SAS_BASE_DIR-mirrored
    location (local_full); copy=True (the default) additionally copies
    it into a local ./DAP directory, which other tools (e.g.
    DAP2tab.py) expect by default -- pass copy=False (-no_cp) to skip.
    '''
    from sdss_access import Access

    xtile='%07d' % tileid
    xtile='%sXX' % xtile[:4]
    print(xtile)
    remote_name = 'dap-%s-%08d.dap.fits.gz' % (DAP_CONFIG, expnum)
    local_full = os.path.join(os.environ['SAS_BASE_DIR'], DAP_TOP, str(drpver), xtile,
                               str(tileid), str(mjd), '%08d' % expnum, remote_name)
    print(local_full)

    a = Access(release='sdsswork')
    try:
        a.remote()
        a.add_file(local_full)
        a.set_stream()
        a.commit()
        print(f"%s successfully downloaded." % remote_name)
        if copy:
            if os.path.isdir('DAP')==False:
                os.makedirs('DAP')
            shutil.copy(local_full, 'DAP/')
    except Exception as e:
        print(f"Failed to download %s: %s" % (remote_name, e))


def _resolve_exposures(words, drpver):
    '''
    Accepts either a single table filename (mjd/expnum columns, and an
    optional tileid column) or "mjd expstart [expstop]" -- matching
    GetFromUtah.py's two input modes. Returns a list of
    (tileid, mjd, expnum) tuples, resolving tileid via
    GetFromUtah.resolve_tileid() for any exposure that doesn't already
    have one (skipping that exposure, with an error already printed by
    resolve_tileid, if it can't be resolved). Returns None (after
    printing why) on a usage error.
    '''
    if len(words)==1:
        exposures=GetFromUtah.read_exposures_table(words[0])
        if exposures is None:
            return None
    elif len(words)>=2:
        mjd=int(words[0])
        exp_start=int(words[1])
        exp_stop=int(words[2]) if len(words)>2 else exp_start
        exposures=[(None,mjd,expnum) for expnum in range(exp_start,exp_stop+1)]
    else:
        print('Error: expected either a table filename (mjd/expnum[/tileid] columns) '
              'or "mjd expstart [expstop]", got:', words)
        return None

    resolved=[]
    for tileid,mjd,expnum in exposures:
        if tileid is None:
            tileid=GetFromUtah.resolve_tileid(drpver,mjd,expnum)
            if tileid is None:
                continue
        resolved.append((tileid,mjd,expnum))
    return resolved


def steer(argv):

    drpver='1.2.1'
    copy=True

    words=[]

    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(_usage_from_doc(__doc__))
            return
        elif argv[i]=='-no_cp':
            copy=False
        elif argv[i]=='-drp':
            i+=1
            drpver=argv[i]
        elif  argv[i][0]=='-':
            print('Unknown option :', argv)
            return
        else:
            words.append(argv[i])
        i+=1

    if not GetFromUtah.check_sas_base_dir():
        return

    exposures=_resolve_exposures(words,drpver)
    if not exposures:
        return

    for tileid,mjd,expnum in exposures:
        get_dap(drpver,tileid,mjd,expnum,copy=copy)

    return

                            






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
