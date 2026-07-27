#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Plot how well the standards that are observed in an LVM exposure
are calibrated.


Command line usage (if any):

    usage: eval_standard.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History::

    240318 ksl Coding begun
    260727 ksl Rewrote compare_with_gaia() -- it was silently failing
        on every current SFrame file for two independent reasons: (1)
        it called ancillary_func.retrive_gaia_star(), which no longer
        exists in the installed lvmdrp (an AttributeError swallowed by
        a bare except); (2) it read STD#BIN/STD#ID/STD#FIB header
        keywords matched against an orig_ifulabel slit label in a
        dedicated 'standard'-targettype fiber set, none of which exist
        any more -- the DRP's current flux calibration
        (fluxCalMethod.py's science_sensitivity) instead identifies
        Gaia-matched field stars among the ordinary science-telescope
        fibers and writes SCI#ID/SCI#FIB (Gaia source id / raw
        fiberid) header keywords, up to 15 slots with gaps. get_standard()
        now indexes FLUX/WAVE directly by that fiberid. GAIA XP spectra
        are now fetched via lvmdrp.core.fluxcal.GaiaXPSpectra, cached
        under $LVM_MASTER_DIR/gaia_cache (the same directory the DRP's
        own flux calibration populates during reduction, so spectra it
        already downloaded are reused instead of re-queried). Retrieval
        is now per-star and failure-tolerant rather than all-or-
        nothing. compare_with_gaia()/qual_eval() now return
        (outfile_or_None, message), so a caller (QuickLook.py) can show
        *why* the comparison failed or partially failed instead of a
        generic could-not-do message.

'''



from astropy.io import fits
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from lvmdrp.core.fluxcal import GaiaXPSpectra



from scipy.signal.windows import boxcar
from scipy.signal import convolve


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


def xsmooth(flux,smooth=21):
    '''
    boxcar smooth the flux
    '''
    if (smooth)>1:
        q=convolve(flux,boxcar(smooth)/float(smooth),mode='same')
        return(q)
    else:
        return(flux)



NSCI_MAX=15

def get_gaia_cache_dir():
    '''
    The Gaia XP spectra cache directory the lvmdrp flux-calibration
    itself uses (LVM_MASTER_DIR/gaia_cache), so spectra it already
    downloaded during reduction are reused instead of re-querying the
    archive. LVM_MASTER_DIR is normally set as soon as lvmdrp is
    imported (via the sdss_access/tree setup in lvmdrp/__init__.py);
    the local './gaia_cache' fallback only applies if that failed.
    '''
    master_dir=os.getenv('LVM_MASTER_DIR')
    if master_dir:
        return os.path.join(master_dir,'gaia_cache')
    return './gaia_cache'


def get_standard(xx,fiberid):
    '''
    Get the spectrum of a single fiber (by fiberid) from an lvmSFrame/
    lvmCFrame file
    '''
    wave=xx['WAVE'].data
    flux=xx['FLUX'].data[fiberid-1]
    return wave,flux


def get_header_stars(header):
    '''
    Retrieve the (fiberid, gaia_id) pairs for the science-telescope
    fibers that land on a Gaia-matched field star, from the SCI#ID/
    SCI#FIB header keywords written by lvmdrp's flux calibration
    (science_sensitivity in fluxCalMethod.py). Not every slot 1..15 is
    populated -- stars that failed acquisition or matching leave gaps.
    '''
    fibers=[]
    gaia_ids=[]
    for i in range(1,NSCI_MAX+1):
        try:
            gaia_id=header['SCI%dID' % i]
            fiber=header['SCI%dFIB' % i]
        except KeyError:
            continue
        fibers.append(fiber)
        gaia_ids.append(gaia_id)
    return fibers,gaia_ids



def compare_with_gaia(filename='lvmSFrame-00005059.fits',outroot=''):
    '''
    Compare the flux-calibrated spectra of the Gaia-matched field
    stars in filename to their Gaia XP spectra.

    Returns (outfile, message): outfile is the plot filename on
    success and None on failure; message explains why on failure, and
    is a non-empty warning (but still returns an outfile) if only some
    of the stars could be retrieved/plotted.
    '''
    try:
        x=fits.open(filename)
    except Exception as e:
        return None,'Could not open %s (%s)' % (filename,e)

    header=x[0].header
    exposure=header['EXPOSURE']
    mjd=header['MJD']

    fibers,gaia_ids=get_header_stars(header)
    if len(fibers)==0:
        return None,'No standard-star header keywords (SCI#ID/SCI#FIB) were found in this file'

    gaia=GaiaXPSpectra(cache_dir=get_gaia_cache_dir())

    plt.figure(1,(8,8))
    plt.clf()
    nfailed=0
    for fiber,gaia_id in zip(fibers,gaia_ids):
        try:
            gaia.fetch_xp_spectra([gaia_id])
            gwave,gflux=gaia.load_xp_spectra(gaia_id)
            swave,sflux=get_standard(x,fiber)
            plt.semilogy(swave,xsmooth(sflux),label=fiber)
            plt.semilogy(gwave,gflux[0],'k')
        except Exception as e:
            nfailed+=1
            print('Error: Failed on GAIA object %s (fiber %s): %s' % (gaia_id,fiber,e))

    if nfailed==len(fibers):
        plt.close(1)
        return None,('Failed to retrieve/plot any of the %d GAIA-matched standard stars '
                      '(no network access to the GAIA archive, and nothing cached locally)' % len(fibers))

    plt.xlim(3500,9500)
    plt.title('MJD %d Exposure %d' % (mjd,exposure))
    ylm=plt.ylim()
    plt.ylim(1e-13,ylm[1])
    plt.tight_layout()

    if outroot=='':
        word=filename.split('/')
        outroot=word[-1].replace('.fits','')
        outfile='standard_%s.png' % outroot
    else:
        outfile=outroot

    message='' if nfailed==0 else '%d of %d standard stars could not be retrieved/plotted' % (nfailed,len(fibers))
    return outfile,message


def qual_eval(filename,outname):
    '''
    This is an extra call so this routine can be run from the qual
    evaluation routines.

    Returns (status, message), see compare_with_gaia.
    '''
    outfile,message=compare_with_gaia(filename,outname)
    if outfile==None:
        return False,message
    plt.savefig(outname)
    plt.close()
    return True,message

                
def steer(argv):

    files=[]

    i=1
    while i<len(argv):
        if argv[i].count('-h'):
            print(_usage_from_doc(__doc__))
            return
        elif argv[i][0]=='-':
            print('Error: could not process command line: ',argv)
        else:
            files.append(argv[i])
        i+=1

    for one in files:
        outfile,message=compare_with_gaia(one)
        if outfile!=None:
            plt.savefig(outfile)
            plt.close()
            if message:
                print('Warning: %s' % message)
        else:
            print('Error: %s' % message)
                          
        



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
