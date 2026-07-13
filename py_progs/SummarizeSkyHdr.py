#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Summarize a fixed list of raw acquisition/astrometry PRIMARY-header
    keywords (reported/commanded/adopted telescope positions and sky
    field names -- see py_progs/test.txt) across many exposures selected
    from a drpall table, for diagnosing where in the DRP a sky
    telescope's position and name can end up disagreeing (see
    check_sky_positions.py in .../duplicates, which found this
    disagreement from the drpall side only; this script pulls in the
    per-exposure header values needed to trace where it happens).

    Exposure selection follows the same convention as SummarizeCframe.py/
    SummarizeSciSky.py (exposure-number range against a drpall table),
    but instead of accumulating spectra, this reads a handful of PRIMARY
    header keywords from each exposure's own CFrame file.

Command line usage (if any):

    usage: SummarizeSkyHdr.py [-h] [-emin 900] [-ver 1.2.1] [-drp_all FILE]
                              [-keywords FILE] [-data_dir DIR] [-out ROOT]
                              exp_start exp_stop [delta]

    Arguments::

        exp_start   starting exposure number
        exp_stop    stopping exposure number
        delta       process every delta-th exposure in range (default 1)

    Options::

        -emin N        minimum exposure time to include (default 900)
        -ver VER       DRP version, used to locate drpall-VER.fits (default 1.2.1)
        -drp_all FILE  explicit drpall table to read instead of drpall-VER.fits
        -keywords FILE keyword/definition table (default: this script's own
                       test.txt, in py_progs/ alongside it) -- see Notes
        -data_dir DIR  look for CFrame files directly in DIR by basename
                       first (a flat local cache), before falling back to
                       the standard xtop/location tree layout
        -out ROOT      output filename root; default is
                       SummarizeSkyHdr_<ver>_<exp_start>_<exp_stop>_<delta>

Description:

    For each selected exposure, opens its CFrame file (from the drpall
    table's own `location` column, SFrame renamed to CFrame -- same
    convention as SummarizeCframe.py/SummarizeSciSky.py) and reads the
    PRIMARY header keywords listed in -keywords (default test.txt) --
    NOT spectra, just header values.

    Output FITS structure::

        PRIMARY    header records the calling parameters (DRPVER, EMIN,
                   EXPSTART, EXPSTOP, DELTA, DRPALL, KWFILE, N_PROC) --
                   same convention as SummarizeSciSky.py
        SKY_HDR    one row per successfully-read exposure: EXPNUM (the
                   join key back to DRP_ALL) plus one column per keyword
                   in -keywords.  Each column's FITS TTYPEn comment card
                   is set to that keyword's definition text from
                   -keywords, so the file is self-documenting without
                   needing test.txt alongside it.
        DRP_ALL    the drpall rows for the exposures actually processed
                   (same rows SKY_HDR was built from -- join on expnum
                   if you need columns from both).

Primary routines:

    process_drpall

Notes:

    test.txt (py_progs/test.txt) is a fixed_width_two_line ascii table
    (keyword, definition) -- read directly with astropy.io.ascii, not
    hardcoded here, so keeping the keyword list current only means
    editing that one file.  A keyword's definition ending in "[deg]" is
    treated as numeric (missing/undefined -> NaN); anything else is
    treated as a string (missing/undefined -> the literal string
    "None", matching drpall's own placeholder convention).

    Some of these keywords (SCIASRC/SKYEASRC/SKYWASRC, the astrometry-
    source quality flags) are absent from older/some exposures' headers
    -- this is expected, not an error; they just come back as "None".

    Exposures whose CFrame file cannot be located (checked in -data_dir
    first if given, then the standard xtop/location tree) are silently
    skipped, same as SummarizeCframe.py/SummarizeSciSky.py.

    Deliberately a separate, standalone script for now (per user
    request) even though it overlaps with SummarizeSciSky.py/
    SummarizeCframe.py's exposure-selection machinery -- may get merged
    with one of those later once the sky position/name diagnosis is
    further along.

History:

260713 ksl Coding begun, to diagnose the skye/skyw position-vs-name
    mismatches found by check_sky_positions.py
    (.../lvm_sky2607/duplicates/check_sky_positions.py) from the drpall
    side -- test.txt's reported/commanded/adopted-position keywords let
    that mismatch be traced to a specific stage of DRP processing.

'''

import sys
import os
import numpy as np
from astropy.io import fits, ascii as apy_ascii
from astropy.table import Table


_USAGE = '''Usage:
  SummarizeSkyHdr.py [-emin 900] [-ver 1.2.1] [-drp_all FILE]
                     [-keywords FILE] [-data_dir DIR] [-out ROOT]
                     exp_start exp_stop [delta]

Arguments:
  exp_start   starting exposure number
  exp_stop    stopping exposure number
  delta       process every delta-th exposure in range (default 1)

Options:
  -emin N        minimum exposure time to include (default 900)
  -ver VER       DRP version, used to locate drpall-VER.fits (default 1.2.1)
  -drp_all FILE  explicit drpall table to read instead of drpall-VER.fits
  -keywords FILE keyword/definition table (default: py_progs/test.txt)
  -data_dir DIR  flat local cache to check for CFrame files before the
                 standard xtop/location tree layout
  -out ROOT      output filename root (default:
                 SummarizeSkyHdr_<ver>_<exp_start>_<exp_stop>_<delta>)
'''


# ──────────────────────────────────────────────────────────────
# drpall reading / exposure selection -- local copies, following
# SummarizeSciSky.py's more complete version of the pattern (see that
# module's own docstring for why this is a local copy, not an import)
# ──────────────────────────────────────────────────────────────

def read_drpall(filename='', drp_ver='1.2.1'):
    '''Read a drpall FITS file, or an ascii table, and return the table.'''
    if filename.count('txt') or filename.count('.tab'):
        try:
            return apy_ascii.read(filename)
        except Exception:
            print('Error: Could not locate : ', filename)
            return []

    DRPFILE = filename if filename else 'drpall-%s.fits' % drp_ver

    if os.path.isfile(DRPFILE):
        xfile = DRPFILE
    else:
        BASEDIR = '/uufs/chpc.utah.edu/common/home/sdss51/sdsswork/lvm/spectro/redux/%s/' % drp_ver
        xfile = '%s/%s' % (BASEDIR, DRPFILE)
        if not os.path.isfile(xfile):
            print('Error: Could not locate : ', xfile)
            return []

    try:
        drpall = fits.open(xfile)
        print('Successfully opened ', xfile)
    except Exception:
        print('Error: Located but could not read : ', xfile)
        return []

    return Table(drpall[1].data)


def select_exps(ztab, exp_start=4000, exp_stop=8000, delta=5, exp_min=900.):
    '''Select every delta-th exposure between exp_start and exp_stop.'''
    xtab = ztab[ztab['expnum'] >= exp_start]
    xtab = xtab[xtab['expnum'] <= exp_stop]
    if exp_min > 0:
        xtab = xtab[xtab['exptime'] >= exp_min]
    if delta > 1:
        xtab = xtab[::delta]
    return xtab


_XTOP     = '/uufs/chpc.utah.edu/common/home/sdss51/'
_XRAINBOW = '/Users/long/Projects/lvm_data/sas'
_XMUSKIE  = '/home/long/Projects/lvm_data/sas'


def find_top():
    '''Locate the top of the local redux data tree (Utah / Rainbow / Muskie).'''
    if os.path.isdir(_XTOP):
        loc, topdir = 'Utah', _XTOP
    elif os.path.isdir(_XRAINBOW):
        loc, topdir = 'Rainbow', _XRAINBOW
    elif os.path.isdir(_XMUSKIE):
        loc, topdir = 'Muskie', _XMUSKIE
    else:
        print('Error: I do not know where I am:', os.getcwd())
        return ''
    print('We are on : ', loc)
    return topdir


def resolve_filename(location, xtop, data_dir=''):
    '''
    Build the CFrame path for one exposure.  If data_dir is given, look
    there directly by basename first (a flat local cache, e.g. a
    manually curated test directory) before falling back to the
    standard xtop/location tree layout.
    '''
    basename = os.path.basename(location).replace('SFrame', 'CFrame')
    if data_dir:
        candidate = os.path.join(data_dir, basename)
        if os.path.isfile(candidate):
            return candidate

    xfile = '%s/%s' % (xtop, location)
    if xfile.count('SFrame'):
        xfile = xfile.replace('SFrame', 'CFrame')
    return xfile


# ──────────────────────────────────────────────────────────────
# Keyword table / per-exposure header reading
# ──────────────────────────────────────────────────────────────

_DEFAULT_KEYWORD_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test.txt')


def load_keyword_defs(keyword_file=_DEFAULT_KEYWORD_FILE):
    '''
    Read a fixed_width_two_line (keyword, definition) table -- default
    py_progs/test.txt.  Returns a list of (keyword, definition,
    is_numeric) tuples, in file order.  is_numeric is True when the
    definition ends in "[deg]" (the RA/Dec keywords), False otherwise
    (the sky-field-name and astrometry-source-quality keywords).
    '''
    tab = apy_ascii.read(keyword_file, format='fixed_width_two_line')
    defs = []
    for row in tab:
        kw = str(row['keyword']).strip()
        definition = str(row['definition']).strip()
        is_numeric = definition.endswith('[deg]')
        defs.append((kw, definition, is_numeric))
    return defs


def read_header_keywords(filename, keyword_defs):
    '''
    Read keyword_defs' keywords from one exposure's PRIMARY header.

    Returns a dict: keyword -> value.  A numeric keyword that is
    missing or FITS-undefined comes back as np.nan; a string keyword
    comes back as the literal string "None" in that case (matching
    drpall's own placeholder convention) -- this keeps every column's
    dtype consistent across all rows.  Returns None if the file itself
    could not be opened.
    '''
    try:
        hdr = fits.getheader(filename, 0)
    except Exception as e:
        print('Warning: could not open %s (%s)' % (filename, e))
        return None

    result = {}
    for kw, _definition, is_numeric in keyword_defs:
        val = hdr.get(kw, None)
        if val is None:
            result[kw] = np.nan if is_numeric else 'None'
        else:
            result[kw] = float(val) if is_numeric else str(val)
    return result


# ──────────────────────────────────────────────────────────────
# Batch processing
# ──────────────────────────────────────────────────────────────

def process_drpall(exp_start, exp_stop, delta=1, exp_min=900., drp_ver='1.2.1',
                   drp_all_file='', keyword_file=_DEFAULT_KEYWORD_FILE,
                   data_dir='', outroot=''):
    '''
    Select exposures from a drpall table and write a combined
    SummarizeSkyHdr FITS file (SKY_HDR + DRP_ALL extensions).
    '''
    xtop = find_top()
    xtab = read_drpall(drp_all_file, drp_ver)
    if len(xtab) == 0:
        print('Error: could not read a drpall table; nothing to do.')
        return
    ztab = select_exps(xtab, exp_start, exp_stop, delta, exp_min)
    print('Selected %d exposures from the drpall table' % len(ztab))

    keyword_defs = load_keyword_defs(keyword_file)
    print('Reading %d header keywords from %s' % (len(keyword_defs), keyword_file))

    hdr_rows  = []
    good_rows = []
    for i in range(len(ztab)):
        xfile = resolve_filename(ztab['location'][i], xtop, data_dir)
        if not os.path.isfile(xfile):
            continue
        values = read_header_keywords(xfile, keyword_defs)
        if values is None:
            continue
        values['EXPNUM'] = int(ztab['expnum'][i])
        hdr_rows.append(values)
        good_rows.append(i)
        if len(good_rows) % 10 == 0:
            print('Finished %d of %d' % (len(good_rows), len(ztab)))

    if not hdr_rows:
        print('Error: no exposures were successfully read; nothing written.')
        return
    print('Successfully read %d of %d selected exposures' % (len(hdr_rows), len(ztab)))

    col_order = ['EXPNUM'] + [kw for kw, _d, _n in keyword_defs]
    hdr_tab = Table(rows=[[row[c] for c in col_order] for row in hdr_rows],
                    names=col_order)

    drp_all_tab = ztab[good_rows]

    if outroot == '':
        outroot = 'SummarizeSkyHdr_%s_%d_%d_%d' % (drp_ver, exp_start, exp_stop, delta)
    outfile = outroot if outroot.endswith('.fits') else outroot + '.fits'

    hdr = fits.Header()
    hdr['Title']    = 'SummarizeSkyHdr'
    hdr['DRPVER']   = drp_ver
    hdr['EMIN']     = exp_min
    hdr['EXPSTART'] = exp_start
    hdr['EXPSTOP']  = exp_stop
    hdr['DELTA']    = delta
    hdr['DRPALL']   = drp_all_file if drp_all_file else ('drpall-%s.fits' % drp_ver)
    hdr['KWFILE']   = os.path.basename(keyword_file)
    hdr['N_PROC']   = len(hdr_rows)

    hdr_hdu = fits.BinTableHDU(hdr_tab, name='SKY_HDR')
    definitions = {kw: definition for kw, definition, _n in keyword_defs}
    for i, col in enumerate(hdr_tab.colnames, start=1):
        if col in definitions:
            hdr_hdu.header.comments['TTYPE%d' % i] = definitions[col]

    hdul = fits.HDUList([
        fits.PrimaryHDU(header=hdr),
        hdr_hdu,
        fits.BinTableHDU(drp_all_tab, name='DRP_ALL'),
    ])
    hdul.writeto(outfile, overwrite=True)
    print('Wrote results to %s' % outfile)


def steer(argv):
    exp_start = -1
    exp_stop  = -1
    delta     = -1
    exp_min   = 900.
    drp_ver   = '1.2.1'
    drp_all_file = ''
    keyword_file = _DEFAULT_KEYWORD_FILE
    data_dir  = ''
    outroot   = ''

    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg in ('-h', '--help'):
            print(_USAGE)
            return
        elif arg == '-emin':
            i += 1
            exp_min = float(argv[i])
        elif arg == '-ver':
            i += 1
            drp_ver = argv[i]
        elif arg == '-drp_all':
            i += 1
            drp_all_file = argv[i]
        elif arg == '-keywords':
            i += 1
            keyword_file = argv[i]
        elif arg == '-data_dir':
            i += 1
            data_dir = argv[i]
        elif arg == '-out':
            i += 1
            outroot = argv[i]
        elif arg.startswith('-'):
            print('Error: unknown option "%s"' % arg)
            print(_USAGE)
            return
        elif exp_start < 0:
            exp_start = int(arg)
        elif exp_stop < 0:
            exp_stop = int(arg)
        elif delta < 0:
            delta = int(arg)
        else:
            print('Error: unexpected argument "%s"' % arg)
            print(_USAGE)
            return
        i += 1

    if exp_start < 0 or exp_stop < 0:
        print(_USAGE)
        return
    if delta < 0:
        delta = 1

    process_drpall(exp_start, exp_stop, delta=delta, exp_min=exp_min, drp_ver=drp_ver,
                   drp_all_file=drp_all_file, keyword_file=keyword_file,
                   data_dir=data_dir, outroot=outroot)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
