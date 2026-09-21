#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Summarize the raw-header fix files kept in lvmcore's hdrfix
    directory (one lvmHdrFix-<mjd>.yaml per MJD) into a single flat
    table of (mjd, exposure, keyword, value) rows -- one row per
    distinct (exposure, keyword) pair, holding whichever fix entry's
    value would actually survive in the header (see Notes).

Command line usage (if any):

    usage: GetHdrfix.py [-h] [-hdr_dir DIR] [-exp_dir DIR] [-out ROOT]

    No arguments are required -- with $LVMCORE_DIR set, running with no
    switches at all reads $LVMCORE_DIR/hdrfix and
    $LVMCORE_DIR/exposure_list directly.

    Options::

        -h            print this help and exit
        -hdr_dir DIR  top-level hdrfix directory to read (one
                      subdirectory per MJD, each containing
                      lvmHdrFix-<mjd>.yaml); default is
                      $LVMCORE_DIR/hdrfix
        -exp_dir DIR  directory holding lvmcore's per-MJD
                      exposure_list_<mjd>.parquet files, used to expand
                      fileroot patterns that wildcard the exposure
                      number (see Description); default is
                      $LVMCORE_DIR/exposure_list
        -out ROOT     output filename root; writes ROOT.txt (default
                      GetHdrfix)

Description:

    Each lvmHdrFix-<mjd>.yaml holds a list of fix entries of the form
    (fileroot, keyword, value), where fileroot is an fnmatch pattern
    (e.g. "sdR-*-*-00003449" or "sdR-*-*-0000068[2-4]") matched by
    lvmdrp.utils.hdrfix.apply_hdrfix() against the raw frame's own
    "sdR-<hemi>-<camera>-<expnum>" string to decide whether that fix
    applies to a given raw frame.

    Most fileroot entries pin down a single exact 8-digit exposure
    number and are read directly with no lookup.  A minority use
    fnmatch wildcard/class syntax (a bare "*", or a bracket expression
    such as "[2-4]" or "[33|79]", possibly spanning the camera field
    too) that can match a whole set of exposures -- what set depends on
    which exposures actually occurred that night.  For these,
    get_real_exposures() reads lvmcore's own
    exposure_list/exposure_list_<mjd>.parquet (one row per exposure
    really taken) and expand_fileroot() replays the exact same
    fnmatch.fnmatch() test apply_hdrfix() itself uses -- against every
    camera in lvmdrp's CAMERAS list and hemi='s' (LVM is LCO-only) --
    so the expansion matches what the DRP would really apply, not a
    guessed reinterpretation of the pattern.

    A fileroot pattern that cannot be resolved (its MJD has no
    exposure_list_<mjd>.parquet on disk, or pandas/pyarrow are not
    importable) is skipped with a printed warning rather than raising;
    the entry is simply absent from the output table.

Primary routines:

    get_hdrfix, expand_fileroot, get_real_exposures

Notes:

    Output is a single ascii.fixed_width_two_line table (ROOT.txt),
    sorted on exposure then keyword, with columns mjd, exposure,
    keyword, value -- exactly one row per distinct (exposure, keyword)
    pair.  When more than one fix entry in a file touches the same
    exposure and keyword (whether from one wildcard pattern matching
    several exposures, or several overlapping patterns), only the value
    from whichever entry appears LAST in that file's fixes list is kept
    -- this reproduces apply_hdrfix()'s own behavior of applying fixes
    in file order and letting a later hdr[keyword] = value overwrite an
    earlier one, so the output reflects what actually ends up in the
    header, not just a raw listing of every entry that mentions the
    pair.  (lvmHdrFix-61163.yaml/lvmHdrFix-61164.yaml were found to
    repeat one block of 6 sky-position fixes verbatim 6 times each for
    exposures 58194/58262 -- harmless here since every repeat carries
    the same value, but worth flagging upstream since it is not obvious
    the duplication is a no-op without checking as this script does.)

    value is always written as a string; a yaml null (used by some
    CALIBFIB entries to mean "clear this keyword") is written as the
    literal string "None", matching drpall's own placeholder
    convention.

    Requires pandas and pyarrow to expand wildcard fileroot patterns
    (only used when a pattern's exposure field is not a plain 8-digit
    number); exact-exposure entries -- the large majority -- need
    neither.

History::

    260921 ksl Coding begun.  Reads lvmcore's hdrfix yaml files, expands
        wildcard fileroot patterns against real per-MJD exposure lists
        (replaying apply_hdrfix()'s own fnmatch test), keeps only the
        last-in-file value per (exposure, keyword) pair to match
        apply_hdrfix()'s sequential-overwrite behavior, and writes a
        table sorted on exposure then keyword.  hdr_dir/exp_dir both
        default independently to $LVMCORE_DIR, so no arguments are
        required for normal use.  Found lvmHdrFix-61163.yaml/
        lvmHdrFix-61164.yaml each repeat one block of 6 sky-position
        fixes verbatim 6 times for exposures 58194/58262 -- harmless
        since every repeat carries the same value, but flagged upstream
        as a likely authoring artifact.
'''

import sys
import os
import re
import fnmatch
import functools
from glob import glob

import yaml
from astropy.table import Table


def _usage_from_doc(doc):
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


_USAGE = _usage_from_doc(__doc__)


_CAMERAS = ['b1', 'b2', 'b3', 'r1', 'r2', 'r3', 'z1', 'z2', 'z3']
_HEMI = 's'   # LVM is exclusively at LCO -- see lvmdrp.utils.metadata

_MJD_RE = re.compile(r'lvmHdrFix-(\d+)\.yaml$')
_EXACT_EXPOSURE_RE = re.compile(r'^\d+$')


def find_hdrfix_files(directory):
    '''Return the sorted list of lvmHdrFix-<mjd>.yaml files under directory.'''
    return sorted(glob('%s/*/lvmHdrFix-*.yaml' % directory))


def read_hdrfix_yaml(yaml_file):
    '''Read one lvmHdrFix-<mjd>.yaml file and return its list of fix dicts (or [] on failure).'''
    try:
        with open(yaml_file) as f:
            data = yaml.safe_load(f)
    except Exception as e:
        print('Error: could not read %s: %s' % (yaml_file, e))
        return []
    fixes = data.get('fixes') if data else None
    return fixes if fixes else []


@functools.lru_cache(maxsize=None)
def get_real_exposures(mjd, exp_dir):
    '''
    Return the sorted tuple of real exposure numbers taken on mjd, read
    from exp_dir/exposure_list_<mjd>.parquet.  Returns None (not an
    empty tuple) if the file is missing or pandas/pyarrow are not
    available, so callers can distinguish "expansion not possible" from
    "genuinely no exposures that night".
    '''
    xfile = os.path.join(exp_dir, 'exposure_list_%d.parquet' % mjd)
    if not os.path.isfile(xfile):
        print('Error: could not locate exposure list for mjd %d: %s' % (mjd, xfile))
        return None

    try:
        import pandas as pd
    except ImportError:
        print('Error: pandas/pyarrow not available -- cannot expand wildcard fileroot patterns')
        return None

    try:
        xtab = pd.read_parquet(xfile)
    except Exception as e:
        print('Error: could not read exposure list %s: %s' % (xfile, e))
        return None

    return tuple(sorted(int(x) for x in xtab['exposure_no']))


def expand_fileroot(fileroot, mjd, exp_dir):
    '''
    Return the list of exposure numbers that fileroot (an fnmatch
    pattern as stored in lvmHdrFix-<mjd>.yaml) actually matches.

    An exact 8-digit exposure field is returned directly with no
    lookup.  Anything else (a bare "*", a bracket expression, or a
    wildcarded camera field) is resolved by testing fnmatch.fnmatch()
    against "sdR-<hemi>-<camera>-<expnum>" for every real exposure that
    night (from get_real_exposures) and every camera in _CAMERAS --
    exactly the test lvmdrp.utils.hdrfix.apply_hdrfix() itself performs.

    Returns None (not []) if the pattern cannot even be checked because
    the MJD's exposure list is unavailable (see get_real_exposures) --
    this is distinct from a pattern that was checked but genuinely
    matched no real exposure that night, which returns [].
    '''
    field = fileroot.rsplit('-', 1)[-1]
    if _EXACT_EXPOSURE_RE.match(field):
        return [int(field)]

    real_exposures = get_real_exposures(mjd, exp_dir)
    if real_exposures is None:
        return None

    matched = []
    for expnum in real_exposures:
        tail = '%08d' % expnum
        for camera in _CAMERAS:
            current_file = 'sdR-%s-%s-%s' % (_HEMI, camera, tail)
            if fnmatch.fnmatch(current_file, fileroot):
                matched.append(expnum)
                break
    return matched


def get_hdrfix(hdr_dir='', exp_dir='', outroot='GetHdrfix'):
    '''
    Read every lvmHdrFix-<mjd>.yaml file under hdr_dir and write a
    single ascii.fixed_width_two_line table (outroot.txt) of
    (mjd, exposure, keyword, value) rows, one row per distinct
    (exposure, keyword) pair actually resolved.  See the module
    docstring for the fileroot expansion rules and the last-entry-wins
    rule used when more than one fix entry touches the same exposure
    and keyword.

    hdr_dir and exp_dir each default to $LVMCORE_DIR/hdrfix and
    $LVMCORE_DIR/exposure_list respectively when not given.

    Returns the table, or None if nothing could be read or resolved.
    '''
    lvmcore = os.environ.get('LVMCORE_DIR', '')

    if hdr_dir == '':
        if lvmcore == '':
            print('Error: -hdr_dir not given and LVMCORE_DIR is not set')
            return None
        hdr_dir = os.path.join(lvmcore, 'hdrfix')

    if exp_dir == '':
        if lvmcore == '':
            print('Error: -exp_dir not given and LVMCORE_DIR is not set')
            return None
        exp_dir = os.path.join(lvmcore, 'exposure_list')

    if not os.path.isdir(hdr_dir):
        print('Error: could not locate directory %s' % hdr_dir)
        return None

    yaml_files = find_hdrfix_files(hdr_dir)
    if len(yaml_files) == 0:
        print('Error: no lvmHdrFix-*.yaml files found under %s' % hdr_dir)
        return None

    # (mjd, exposure, keyword) -> value.  Fix entries are applied in the
    # order they appear in each file, so a later entry overwriting an
    # earlier one here reproduces apply_hdrfix()'s own sequential
    # hdr[keyword] = value overwrite -- last entry in the file wins.
    results = {}

    n_files = 0
    n_entries = 0
    n_unavailable = 0
    n_no_match = 0

    for yaml_file in yaml_files:
        m = _MJD_RE.search(yaml_file)
        if m is None:
            print('Error: could not parse mjd from %s -- skipping' % yaml_file)
            continue
        mjd = int(m.group(1))
        n_files += 1

        for fix in read_hdrfix_yaml(yaml_file):
            n_entries += 1
            fileroot = fix.get('fileroot', '')
            keyword = fix.get('keyword', '')
            value = fix.get('value', None)
            value = 'None' if value is None else str(value)

            exposures = expand_fileroot(fileroot, mjd, exp_dir)
            if exposures is None:
                n_unavailable += 1
                continue
            if len(exposures) == 0:
                n_no_match += 1
                continue

            for expnum in exposures:
                results[(mjd, expnum, keyword)] = value

    print('Read %d fix entries from %d files' % (n_entries, n_files))
    if n_unavailable:
        print('Warning: %d fix entries could not be checked because their exposure list was unavailable and were omitted (see errors above)' % n_unavailable)
    if n_no_match:
        print('Note: %d fix entries were checked but matched no real exposure and were omitted' % n_no_match)

    if len(results) == 0:
        print('Error: no fix entries could be resolved to exposure numbers')
        return None

    mjd_col, exp_col, key_col, val_col = [], [], [], []
    for (mjd, expnum, keyword), value in results.items():
        mjd_col.append(mjd)
        exp_col.append(expnum)
        key_col.append(keyword)
        val_col.append(value)

    xtab = Table([mjd_col, exp_col, key_col, val_col],
                 names=['mjd', 'exposure', 'keyword', 'value'])
    xtab.sort(['exposure', 'keyword'])

    outname = '%s.txt' % outroot
    xtab.write(outname, format='ascii.fixed_width_two_line', overwrite=True)
    print('Wrote %d rows to %s' % (len(xtab), outname))

    return xtab


def steer(argv):
    hdr_dir = ''
    exp_dir = ''
    outroot = 'GetHdrfix'

    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg in ('-h', '--help'):
            print(_USAGE)
            return
        elif arg == '-hdr_dir':
            i += 1
            hdr_dir = argv[i]
        elif arg == '-exp_dir':
            i += 1
            exp_dir = argv[i]
        elif arg == '-out':
            i += 1
            outroot = argv[i]
        else:
            print('Error: unknown argument "%s"' % arg)
            print(_USAGE)
            return
        i += 1

    get_hdrfix(hdr_dir=hdr_dir, exp_dir=exp_dir, outroot=outroot)


if __name__ == '__main__':
    # No arguments are required (see Command line usage) -- unlike most
    # scripts here, an empty argv is a normal invocation, not a request
    # for help, so steer() is always called (same convention as
    # CheckData.py/CheckReduced.py).
    steer(sys.argv)
