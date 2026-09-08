#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Rerun the DRP's own sky-subtraction routine
    (skyMethod.quick_sky_subtraction) on a single lvmCFrame file,
    producing a corrected lvmSFrame plus diagnostic plots --
    typically used on the output of SubstituteSky.py.

Command line usage (if any)::

    usage: RunSky.py [-h] filename

    where
        -h        prints this documentation and exits
        filename  the lvmCFrame file to run sky subtraction on

Description:

    Loads filename, derives its exposure/MJD/tile from its own header,
    and calls quick_sky_subtraction on it with the 'farlines_nearcont'
    method, writing the resulting lvmSFrame (filename's basename with
    'CFrame' swapped for 'SFrame') and a skyQA PDF into the current
    directory. It then reads back the ancillary skytable
    quick_sky_subtraction just wrote and plots the SCI/SkyE/SkyW/
    SkyE_super/SkyW_super mean spectra for inspection.

    Calling quick_sky_subtraction outside the full science_reduction
    pipeline exposes three missing-directory bugs in lvmdrp -- that
    pipeline happens to pre-create these directories as a side effect
    of an earlier step, so they never surface there. This routine
    works around all three locally rather than patching the vendored
    lvmdrp package: (1) quick_sky_subtraction's own ancillary skytable
    write has no os.makedirs, unlike every other lvmdrp write site --
    created here from the same lvm_anc/kind='sky' path it uses
    internally (via lvmdrp's own path/drpver, so it resolves under
    whatever SAS_BASE_DIR/LVMDRP_VERSION is active); (2)
    writeFitsData's unconditional
    os.makedirs(os.path.dirname(out_file)) raises FileNotFoundError
    when out_file is a bare filename (dirname is '') -- out_sframe is
    always made an absolute path; (3) run_qa's skyQA PDF write into
    <dirname(out_sframe)>/qa/ has the same gap -- pre-created too.

Primary routines:

    steer is the entry point; get_sky_tab/plot_one are the plotting
    helpers

Notes:

    The ancillary skytable is read back from the same lvm_anc path
    just computed for the mkdir fix, rather than searched for on disk
    (find_skytab()'s old rglob-based lookup, removed) -- this both
    avoids an IndexError for an exposure with no prior local
    reduction, and guarantees the plots reflect this run's own sky
    model rather than a stale skytable from some other exposure or
    DRP version that happened to match the filename search.

History::

    240419 ksl Coding begun
    260908 ksl Fixed steer() to accept a CFrame filename directly (as
        the usage line always said), instead of eval()-ing the
        argument as an exposure number and reconstructing a hardcoded
        path; exposure/MJD/tile now come from the file's own header.
        out_sframe is derived from the input filename instead of a
        hardcoded 'foo.fits', and the skip_subtraction kwarg
        (unsupported by the current lvmdrp signature) was removed.
        Worked around three missing-directory bugs in standalone
        quick_sky_subtraction calls (see Description), and replaced
        the rglob-based skytable lookup with a direct read of the
        path quick_sky_subtraction itself just wrote.

'''

import sys
from astropy.io import ascii,fits
import matplotlib.pyplot as plt
from astropy.io import fits,ascii
from astropy.table import Table
import numpy as np
from glob import glob
from lvmdrp.functions.skyMethod import quick_sky_subtraction
from lvmdrp import path, __version__ as drpver
import os
from lvm_ksl import sky_plot


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


def get_sky_tab(filename):
    '''
    Retrive a sky table and split it into 
    tables
    '''
    x=fits.open(filename)
    xsci=Table(x['SCI'].data)
    xsky_e=Table(x['SKYE'].data)
    xsky_w=Table(x['SKYW'].data)
    xsky_e_super=Table(x['SKYE_SUPER'].data)
    xsky_w_super=Table(x['SKYW_SUPER'].data)
    return xsci,xsky_e,xsky_w,xsky_e_super,xsky_w_super


def plot_one(xtab):
    '''
    Plot a single skytab table
    '''


    plt.figure(1,(6,6))
    plt.clf()
    plt.subplot(3,1,1)
    plt.semilogy(xtab['WAVE'],xtab['FLUX'])
    plt.tight_layout()
    try:
        plt.semilogy(xtab['WAVE'],xtab['CONT'])
    except:
        print('No Cont')
        return
    ymin,ymax=plt.ylim()
    plt.xlim(3500,9800)
    plt.subplot(3,1,2)
    plt.semilogy(xtab['WAVE'],xtab['LINES'],label='Lines')
    plt.legend()
    plt.ylim(ymin,ymax)
    plt.xlim(3500,9800)
    plt.tight_layout()
    plt.subplot(3,1,3)
    plt.plot(xtab['WAVE'],xtab['FLUX']-xtab['CONT'])
    plt.plot([3600,9800],[0,0],'r')
    plt.ylim(-.1e-13,0.1e-13)
    plt.xlim(3500,9800)
    plt.tight_layout()
                       


def steer(argv):
    '''
    This is mainly a steering routine
    '''

    i=1
    in_cframe=None
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(_usage_from_doc(__doc__))
            return
        elif argv[i][0]=='-':
            print('Error: Poorly formated command line: ',argv)
            return
        elif in_cframe is None:
            in_cframe=argv[i]
        else:
            print('Error: Poorly formated command line: ',argv)
            return
        i+=1


    if in_cframe is None:
        return

    hdr = fits.getheader(in_cframe, 0)
    exposure = hdr['EXPOSURE']
    base = os.path.basename(in_cframe)
    out_name = base.replace('CFrame', 'SFrame') if 'CFrame' in base else 'sframe_' + base
    # lvmdrp's writeFitsData() does os.makedirs(os.path.dirname(out_file),
    # exist_ok=True) unconditionally; for a bare filename dirname() is ''
    # and os.makedirs('') raises FileNotFoundError. Always pass a path
    # with an explicit directory (here, the current one) to avoid that.
    out_sframe = os.path.join(os.getcwd(), out_name)
    sky_method='farlines_nearcont'

    print('Working on exposure %d (%s)' % (exposure, in_cframe))

    # quick_sky_subtraction writes its own ancillary skytable under the
    # lvm_anc/.../<drpver>/.../ancillary/ directory but, unlike every
    # other write site in lvmdrp, never creates that directory itself --
    # it only exists already when called from within the full
    # science_reduction pipeline, which creates it as a side effect of
    # an earlier step. Recreate that step here so a standalone call
    # doesn't crash with FileNotFoundError. This is also the exact path
    # quick_sky_subtraction will (re)write below, so using it directly
    # -- rather than find_skytab()'s rglob search for a pre-existing
    # file -- both works for an exposure with no prior local reduction
    # and guarantees we read back the fresh skytable this call just
    # produced, not a stale one from some other drpver found by chance.
    skytable_out = path.full('lvm_anc', mjd=hdr['SMJD'], tileid=hdr['TILE_ID'], drpver=drpver,
                              kind='sky', camera='brz', imagetype='table', expnum=exposure)
    os.makedirs(os.path.dirname(skytable_out), exist_ok=True)

    # same gap for the skyQA plots quick_sky_subtraction writes afterward
    # (via run_qa) into <dirname(out_sframe)>/qa/ -- that subdirectory is
    # never created for us either.
    os.makedirs(os.path.join(os.path.dirname(out_sframe), 'qa'), exist_ok=True)

    quick_sky_subtraction(in_cframe, out_sframe, skymethod=sky_method)

    x=fits.open(skytable_out)
    x.info()

    sci,sky_e,sky_w,sky_e_super,sky_w_super=get_sky_tab(skytable_out)

    plot_one(sci)
    plot_one(sky_e)
    plot_one(sky_w)
    plot_one(sky_e_super)
    plot_one(sky_w_super)

    sky_plot.eval_qual(out_sframe)
    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
