#!/usr/bin/env python
# coding: utf-8

'''
Check how well the sky positions recorded in a drpall file (or a
SummarizeSkyHdr.py output file) agree with the nominal positions of the
same named fields in final_sky_tiles.csv, one row per problem exposure.

Usage::

    check_sky_positions.py [-csv FILE] [-tol DEG] [-out FILE]
                          [-postype reported|commanded|adopted] [drpall_file]

Arguments::

    drpall_file    drpall FITS file (default: drpall-1.2.1.fits), or a
                   SummarizeSkyHdr.py output file (has a SKY_HDR
                   extension alongside DRP_ALL) if -postype is anything
                   other than 'adopted'

Options::

    -csv FILE      sky tile position catalog (default: final_sky_tiles.csv)
    -tol DEG       agreement tolerance in degrees (default: 0.1)
    -out FILE      output table, ascii.fixed_width_two_line format
                   (default: sky_problem_check_<drpall stem>[_<postype>].tab
                   -- the postype suffix is omitted for the default
                   'adopted'; includes the drpall root name so outputs
                   from different drpall versions don't overwrite each
                   other)
    -postype T     which of the three recorded sky-telescope positions to
                   check against -- 'reported' (what the telescope itself
                   reported back), 'commanded' (the target position sent
                   to it), or 'adopted' (default; the final
                   astrometry-refined position drpall itself uses).  Only
                   'adopted' works on a plain drpall-*.fits; 'reported'/
                   'commanded' require a SummarizeSkyHdr.py output file's
                   SKY_HDR extension, since drpall alone never recorded
                   those.  The sky field NAME (skye_name/skyw_name) is
                   the same regardless of -postype -- SKY_HDR only ever
                   records it once, at the commanded stage -- so this
                   only changes which POSITION the name is compared
                   against.

For every exposure with both skye_name and skyw_name set (non-blank,
non-"None"), each side's recorded (ra, dec) is matched to its single
nearest CSV catalog position (great-circle separation; a match only
counts if that nearest separation is <= tol degrees), independent of
whatever name label it happens to carry.  A row is a problem if skye or
skyw (or both) doesn't match its OWN name's catalog position.

Each problem exposure gets exactly one output row, classified by
problem_type::

    'Swapped'    skye_name is really skyw's position AND skyw_name is
                 really skye's, in the same row -- the specific,
                 already-understood bug pattern
    'HalfMatch'  only ONE of those two relationships holds -- e.g.
                 skyw_name happens to equal skye's true position, but
                 skye_name does NOT equal skyw's true position (or vice
                 versa).  Not a clean pairwise swap, but not random
                 either -- e.g. consistent with a name shifted by one
                 slot across a sequence of pointings, rather than
                 crossed between exactly two.  Worth a look before
                 writing off as unexplained.
    'Unknown'    neither relationship holds -- no known pattern yet

Output table columns: tileid, mjd, expnum, filename, skye_name,
skyw_name, east_true, west_true (each side's real nearest-catalog
position -- the literal string "None" if nothing in the catalog is
within tol), problem_type.  Sorted by problem_type (Swapped, then
HalfMatch, then Unknown) then mjd.

The summary (exposures checked; problem count by type) is printed both
at the start and again at the end, since the table itself can be
hundreds of rows -- easy to lose the first copy off a terminal's visible
scrollback.

History::

    260712 ksl Coding begun, to check a real name/position mismatch
        found via EsoSkyFit.py testing on WHAM_south_08.
    260713 ksl Fixed a real bug: this script opened the drpall table via
        a hardcoded fits.open(...)[1] (extension index 1) -- correct for
        drpall-*.fits (PRIMARY + one unnamed table) but wrong for a
        SummarizeSkyHdr.py output file, which has PRIMARY, SKY_HDR, then
        DRP_ALL -- index 1 there is SKY_HDR, so it was silently reading
        the wrong table.  Added _drpall_hdu(), which looks up 'DRP_ALL'
        by name first and only falls back to index 1 when the table is
        unnamed.
    260713 ksl Added -postype (reported/commanded/adopted), since
        SummarizeSkyHdr.py's SKY_HDR extension carries three independent
        position pairs per sky telescope sharing one commanded name --
        this traces which DRP stage a mismatch first appears at.  Made
        nearest_catalog_match() NaN-safe for real missing SKY_HDR values.
    260713 ksl Added 'HalfMatch' as a third problem_type alongside
        Swapped/Unknown: a real example (tileid 1056437, mjd 61061)
        showed a one-directional name/position match that the strict
        two-way swap test was lumping in with genuinely unexplained
        cases -- turned out to affect 22% of former-Unknown rows.
    260713 ksl Added sci_asrc/skye_asrc/skyw_asrc (astrometry-quality
        flags) to the per-exposure table -- 'CMD position' (no real
        astrometry fit) is a candidate explanation for some HalfMatch/
        Unknown rows, though Swapped rows turned out to universally use
        real 'GDR coadd' astrometry, ruling that out for those.
    260713 ksl Both the per-exposure problem listing and the per-name
        aggregate are always produced together on every run.
'''

import sys
import os
import re
import numpy as np
from astropy.io import fits
from astropy.table import Table


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before a line consisting of "History:"
    (whitespace-insensitive), so -h stays short even as that section
    grows -- without hand-duplicating the Synopsis/Options text in a
    second string.  Anchored to a whole line (not a bare substring
    search) so it can't misfire on "History:" appearing mid-sentence,
    and returns doc unchanged if no such line is present.
    '''
    m = re.search(r'^\s*History:\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


_USAGE = _usage_from_doc(__doc__)


def angsep_deg(ra1, dec1, ra2, dec2):
    '''Vectorized great-circle separation (haversine formula), in degrees.'''
    ra1, dec1, ra2, dec2 = (np.radians(a) for a in (ra1, dec1, ra2, dec2))
    dra = ra2 - ra1
    ddec = dec2 - dec1
    a = np.sin(ddec / 2.0)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(dra / 2.0)**2
    return np.degrees(2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0))))


def load_csv_positions(csv_file):
    '''Return (names, ra, dec) arrays of the CSV catalog, name stripped.'''
    data = np.genfromtxt(csv_file, delimiter=',', names=True, dtype=None,
                        encoding='utf-8', autostrip=True)
    names = np.char.strip(data['name'].astype(str))
    return names, np.array(data['ra'], dtype=float), np.array(data['dec'], dtype=float)


def nearest_catalog_match(ra, dec, cat_names, cat_ra, cat_dec, tol):
    '''
    For each (ra, dec) observation, find its nearest CSV catalog position
    and return the catalog name if that separation is <= tol, else ''.

    Vectorized over the full (n_obs, n_cat) separation matrix -- fine at
    this scale (tens of thousands of observations x hundreds of catalog
    entries is a few hundred MB, not a problem).

    NaN-safe: an observation with a NaN ra or dec (real for the
    reported/commanded SKY_HDR position columns -- e.g. SummarizeSkyHdr.py
    output has 62 NaN POSKYERA values in one real file, for exposures
    where that keyword was missing from the header) is treated as "no
    match" ('') directly, without computing separations for it -- an
    all-NaN row would otherwise make np.argmin's result undefined.
    '''
    matched = np.full(len(ra), '', dtype=cat_names.dtype)
    valid = np.isfinite(ra) & np.isfinite(dec)
    if not valid.any():
        return matched

    sep = angsep_deg(ra[valid, None], dec[valid, None], cat_ra[None, :], cat_dec[None, :])
    best_idx = np.argmin(sep, axis=1)
    best_sep = sep[np.arange(len(best_idx)), best_idx]
    matched[valid] = np.where(best_sep <= tol, cat_names[best_idx], '')
    return matched


def _drpall_hdu(hdul):
    '''
    Return the HDU holding the drpall table, whichever file convention
    this is: the original ``drpall-*.fits`` (a single unnamed table at
    extension 1) or a SummarizeSkyHdr.py-style file (an explicitly
    named 'DRP_ALL' extension alongside others, e.g. SKY_HDR -- extension
    index alone is not reliable there).
    '''
    try:
        return hdul['DRP_ALL']
    except KeyError:
        return hdul[1]


# SKY_HDR (SummarizeSkyHdr.py output) carries three independent position
# pairs per telescope, all keyed off the SAME sky field name (SKYENAME/
# SKYWNAME -- the name is only ever recorded once, at the commanded
# stage) -- see SkyPosKeywords.py/test.txt for the exact definitions:
#   reported  : what the telescope itself reported back (TES*)
#   commanded : the target position sent to the telescope (POS*)
#   adopted   : the final astrometry-refined position (also DRP_ALL's
#               own skye_ra/skye_dec/skyw_ra/skyw_dec -- the default,
#               and the only postype that works without a SKY_HDR
#               extension, since it comes straight from DRP_ALL)
_POSTYPES = {
    'reported':  ('TESKYERA', 'TESKYEDE', 'TESKYWRA', 'TESKYWDE'),
    'commanded': ('POSKYERA', 'POSKYEDE', 'POSKYWRA', 'POSKYWDE'),
    'adopted':   ('SKYERA',   'SKYEDEC',  'SKYWRA',   'SKYWDEC'),
}


def load_drpall_rows(drpall_file, postype='adopted'):
    '''
    Return a Table with one row per drpall exposure, keeping skye and skyw
    row-aligned -- needed to check whether a single row's own east/west
    pair got crossed.  Columns: tileid, mjd, expnum, filename, skye_name,
    skye_ra, skye_dec, skyw_name, skyw_ra, skyw_dec.

    postype : 'reported', 'commanded', or 'adopted' (default) -- see
        _POSTYPES above.  'adopted' with no SKY_HDR extension present
        reads skye_ra/skye_dec/skyw_ra/skyw_dec straight from the drpall
        table itself (the original, pre-SKY_HDR behavior).  Any other
        postype -- or 'adopted' when a SKY_HDR extension IS present --
        is read from SKY_HDR instead, joined back onto the drpall rows
        by EXPNUM; raises ValueError if the file has no SKY_HDR
        extension in that case.
    '''
    if postype not in _POSTYPES:
        raise ValueError('postype must be one of %s' % sorted(_POSTYPES))

    hdul = fits.open(drpall_file)
    tab = Table(_drpall_hdu(hdul).data)
    out = Table()
    out['tileid']   = tab['tileid']
    out['mjd']      = tab['mjd']
    out['expnum']   = tab['expnum']
    out['filename'] = tab['filename']

    if postype == 'adopted' and 'SKY_HDR' not in hdul:
        out['skye_name'] = np.char.strip(np.array(tab['skye_name'], dtype=str))
        out['skye_ra']   = np.array(tab['skye_ra'], dtype=float)
        out['skye_dec']  = np.array(tab['skye_dec'], dtype=float)
        out['skyw_name'] = np.char.strip(np.array(tab['skyw_name'], dtype=str))
        out['skyw_ra']   = np.array(tab['skyw_ra'], dtype=float)
        out['skyw_dec']  = np.array(tab['skyw_dec'], dtype=float)
        return out

    if 'SKY_HDR' not in hdul:
        raise ValueError("postype=%r requires a SKY_HDR extension (from "
                         "SummarizeSkyHdr.py); %s does not have one"
                         % (postype, drpall_file))

    e_ra_col, e_dec_col, w_ra_col, w_dec_col = _POSTYPES[postype]
    sky_hdr = Table(hdul['SKY_HDR'].data)

    # SKY_HDR and DRP_ALL come from the same selected exposures in a
    # SummarizeSkyHdr.py file, but don't assume row order matches --
    # align explicitly on EXPNUM
    idx = {int(e): i for i, e in enumerate(sky_hdr['EXPNUM'])}
    sel = [idx[int(e)] for e in out['expnum']]
    sky_hdr = sky_hdr[sel]

    out['skye_name'] = np.char.strip(np.array(sky_hdr['SKYENAME'], dtype=str))
    out['skye_ra']   = np.array(sky_hdr[e_ra_col], dtype=float)
    out['skye_dec']  = np.array(sky_hdr[e_dec_col], dtype=float)
    out['skyw_name'] = np.char.strip(np.array(sky_hdr['SKYWNAME'], dtype=str))
    out['skyw_ra']   = np.array(sky_hdr[w_ra_col], dtype=float)
    out['skyw_dec']  = np.array(sky_hdr[w_dec_col], dtype=float)

    # astrometry-quality flags ('GDR coadd' = real guider fit, 'CMD
    # position' = fell back to the commanded position) -- only exist in
    # SKY_HDR, and only added if actually present there (e.g. a custom
    # -keywords file passed to SummarizeSkyHdr.py could omit them)
    for src_col, out_col in (('SCIASRC', 'sci_asrc'), ('SKYEASRC', 'skye_asrc'),
                             ('SKYWASRC', 'skyw_asrc')):
        if src_col in sky_hdr.colnames:
            out[out_col] = np.char.strip(np.array(sky_hdr[src_col], dtype=str))

    return out


def _load_and_match_rows(drpall_file, csv_file, tol=0.1, postype='adopted'):
    '''
    Shared helper for check_problems: load drpall keeping skye/skyw
    row-aligned, restrict to rows where both are real
    (non-blank/non-"None") labels, and find each side's own "true"
    (nearest-catalog-position) match independently.  Returns (rows,
    east_true, west_true, skye_name, skyw_name).

    postype: see load_drpall_rows/_POSTYPES.
    '''
    cat_names, cat_ra, cat_dec = load_csv_positions(csv_file)
    rows = load_drpall_rows(drpall_file, postype=postype)

    keep = ((rows['skye_name'] != '') & (rows['skye_name'] != 'None') &
            (rows['skyw_name'] != '') & (rows['skyw_name'] != 'None'))
    rows = rows[keep]

    east_true = nearest_catalog_match(np.array(rows['skye_ra']), np.array(rows['skye_dec']),
                                      cat_names, cat_ra, cat_dec, tol)
    west_true = nearest_catalog_match(np.array(rows['skyw_ra']), np.array(rows['skyw_dec']),
                                      cat_names, cat_ra, cat_dec, tol)

    skye_name = np.array(rows['skye_name'])
    skyw_name = np.array(rows['skyw_name'])
    return rows, east_true, west_true, skye_name, skyw_name


def _swap_mask(east_true, west_true, skye_name, skyw_name):
    '''
    True for rows where skye_name is really skyw's position and vice
    versa (and the two labels differ).
    '''
    return ((east_true == skyw_name) & (west_true == skye_name)
           & (skye_name != skyw_name))


def check_positions(drpall_file, csv_file, tol=0.1, postype='adopted'):
    '''
    Build a per-NAME aggregate: for every catalog sky field, how many
    times it was labeled, how often that label was right, and (crossing
    the other way) how many times its real position was actually
    observed regardless of what name was on it at the time.

    postype selects which of the three SKY_HDR position pairs to check
    against (see load_drpall_rows/_POSTYPES) -- same convention as
    check_problems.

    skye and skyw are pooled together (each drpall row contributes one
    "observation" per side).  For classifying a row's own east/west pair
    as Swapped/HalfMatch, only rows with both skye_name and skyw_name
    real (non-blank/non-"None") are eligible -- same definitions as
    check_problems's _swap_mask/half1^half2 logic -- and both sides of
    such a row share that row's classification.

    Returns (result_table, unmatched_names).  result_table columns::

        name                        catalog field name
        n_observed                  times this name was used as a label
                                     (either side)
        n_correct                   of those, times the recorded position
                                     really was this field
        n_incorrect                 n_observed - n_correct
        n_swapped                   of the incorrect ones, times explained
                                     by a clean two-way swap with the
                                     other side of the same row
        n_halfmatch                 of the incorrect ones, times explained
                                     by a one-directional match (see
                                     check_problems)
        n_unexplained               n_incorrect - n_swapped - n_halfmatch
        mean_offset_unexplained_deg mean separation (deg) between the
                                     recorded position and this name's
                                     catalog position, over only the
                                     n_unexplained rows (NaN if none)
        n_position_observed         times this field's real position was
                                     actually recorded (by either side,
                                     under any label), excluding swap/
                                     halfmatch rows (those are "explained"
                                     elsewhere, not really this field
                                     being freshly observed)
        n_position_mislabeled       of n_position_observed, times the
                                     label on it was NOT this name

    Sorted by name.  Only rows with n_observed > 0 or
    n_position_observed > 0 are included (a field nobody ever pointed at
    or landed on isn't interesting here).

    unmatched_names: labels actually used in drpall (either side) that
    don't match any name in the CSV catalog at all -- these can never
    have n_correct > 0 by construction, since nearest_catalog_match only
    ever returns a catalog name.
    '''
    cat_names, cat_ra, cat_dec = load_csv_positions(csv_file)
    rows = load_drpall_rows(drpall_file, postype=postype)

    east_name = np.char.strip(np.array(rows['skye_name'], dtype=str))
    east_ra   = np.array(rows['skye_ra'], dtype=float)
    east_dec  = np.array(rows['skye_dec'], dtype=float)
    west_name = np.char.strip(np.array(rows['skyw_name'], dtype=str))
    west_ra   = np.array(rows['skyw_ra'], dtype=float)
    west_dec  = np.array(rows['skyw_dec'], dtype=float)

    east_true = nearest_catalog_match(east_ra, east_dec, cat_names, cat_ra, cat_dec, tol)
    west_true = nearest_catalog_match(west_ra, west_dec, cat_names, cat_ra, cat_dec, tol)

    both_real = ((east_name != '') & (east_name != 'None') &
                (west_name != '') & (west_name != 'None'))

    swap_row = np.zeros(len(rows), dtype=bool)
    half_row = np.zeros(len(rows), dtype=bool)
    if both_real.any():
        half1 = (east_true[both_real] == west_name[both_real])
        half2 = (west_true[both_real] == east_name[both_real])
        swap = _swap_mask(east_true[both_real], west_true[both_real],
                          east_name[both_real], west_name[both_real])
        swap_row[both_real] = swap
        half_row[both_real] = (half1 ^ half2) & ~swap

    # pool east+west into flat parallel arrays, one entry per side per
    # row, keeping only real (non-blank) labels
    names   = np.concatenate([east_name, west_name])
    ra      = np.concatenate([east_ra, west_ra])
    dec     = np.concatenate([east_dec, west_dec])
    matched = np.concatenate([east_true, west_true])
    is_swap = np.concatenate([swap_row, swap_row])
    is_half = np.concatenate([half_row, half_row])

    real = (names != '') & (names != 'None')
    names, ra, dec, matched, is_swap, is_half = (
        a[real] for a in (names, ra, dec, matched, is_swap, is_half))

    unmatched_names = sorted(set(np.unique(names)) - set(cat_names))

    cat_ra_of  = dict(zip(cat_names, cat_ra))
    cat_dec_of = dict(zip(cat_names, cat_dec))

    result_rows = []
    for name in cat_names:
        is_label    = (names == name)
        n_observed  = int(is_label.sum())

        correct     = is_label & (matched == name)
        n_correct   = int(correct.sum())
        n_incorrect = n_observed - n_correct

        n_swapped   = int((is_label & is_swap).sum())
        n_halfmatch = int((is_label & is_half & ~is_swap).sum())
        n_unexplained = n_incorrect - n_swapped - n_halfmatch

        unexplained = is_label & ~correct & ~is_swap & ~is_half
        if unexplained.any():
            offsets = angsep_deg(ra[unexplained], dec[unexplained],
                                 cat_ra_of[name], cat_dec_of[name])
            mean_offset = float(np.mean(offsets))
        else:
            mean_offset = np.nan

        is_position = (matched == name) & ~is_swap & ~is_half
        n_position_observed   = int(is_position.sum())
        n_position_mislabeled = int((is_position & (names != name)).sum())

        if n_observed == 0 and n_position_observed == 0:
            continue

        result_rows.append((name, n_observed, n_correct, n_incorrect,
                            n_swapped, n_halfmatch, n_unexplained, mean_offset,
                            n_position_observed, n_position_mislabeled))

    result_table = Table(rows=result_rows,
                         names=['name', 'n_observed', 'n_correct', 'n_incorrect',
                               'n_swapped', 'n_halfmatch', 'n_unexplained',
                               'mean_offset_unexplained_deg',
                               'n_position_observed', 'n_position_mislabeled'])
    result_table.sort('name')

    return result_table, unmatched_names


def check_problems(drpall_file, csv_file, tol=0.1, postype='adopted'):
    '''
    Build a one-row-per-problem-exposure listing for case-by-case
    investigation.

    postype selects which of the three SKY_HDR position pairs to check
    against (see load_drpall_rows/_POSTYPES) -- 'reported', 'commanded',
    or 'adopted' (default; the only one that works without a SKY_HDR
    extension).

    A row is included if skye or skyw (or both) is mislabeled -- i.e.
    its recorded position doesn't match its own name's catalog position.
    Each such exposure gets exactly one row, classified by
    problem_type::

        'Swapped'    skye_name is really skyw's position AND skyw_name is
                     really skye's, in the same row (see _swap_mask) --
                     the specific, already-understood bug pattern
        'HalfMatch'  only ONE of those two relationships holds -- e.g.
                     skyw_name happens to equal skye's true position, but
                     skye_name does NOT equal skyw's true position (or
                     vice versa).  Not a clean pairwise swap, but not
                     random either -- e.g. consistent with a name shifted
                     by one slot across a sequence of pointings, rather
                     than crossed between exactly two.  Worth a look
                     before writing off as unexplained.
        'Unknown'    neither relationship holds -- no known pattern yet

    Returns (problem_table, n_checked).  problem_table columns: tileid,
    mjd, expnum, filename, skye_name, skyw_name, east_true, west_true
    (each side's real nearest-catalog position -- the literal string
    "None" if nothing in the catalog is within tol; this can only happen
    when the position doesn't match any known field, since
    nearest_catalog_match never returns a name unless it's within tol in
    the first place), problem_type, and -- only when reading a
    SummarizeSkyHdr.py file that has them -- sci_asrc/skye_asrc/skyw_asrc
    (the astrometry-quality flags: 'GDR coadd' means a real guider
    astrometry fit was used; 'CMD position' means it fell back to the
    commanded position, i.e. no real astrometry at all -- a likely
    explanation for some HalfMatch/Unknown rows).  Sorted by
    problem_type (Swapped, then HalfMatch, then Unknown) then mjd.
    n_checked is the number of rows where both skye_name and skyw_name
    are real (non-blank/non-"None") labels, i.e. the total this was
    checked against.
    '''
    rows, east_true, west_true, skye_name, skyw_name = _load_and_match_rows(
        drpall_file, csv_file, tol, postype=postype)

    bad = (east_true != skye_name) | (west_true != skyw_name)

    half1 = (east_true == skyw_name)   # west's (wrong) label matches east's truth
    half2 = (west_true == skye_name)   # east's (wrong) label matches west's truth
    swap  = half1 & half2 & (skye_name != skyw_name)
    half_match = (half1 ^ half2) & ~swap

    problem_type = np.full(len(skye_name), 'Unknown', dtype='<U9')
    problem_type[half_match] = 'HalfMatch'
    problem_type[swap] = 'Swapped'

    # nearest_catalog_match already never returns a name beyond tol -- ''
    # here always means "nothing in the catalog is close enough", so it's
    # safe to relabel as the literal string "None" for clarity (matching
    # the same placeholder convention drpall itself uses for skye_name/
    # skyw_name) rather than leaving it ambiguously blank.
    east_display = np.where(east_true == '', 'None', east_true)
    west_display = np.where(west_true == '', 'None', west_true)

    out_cols = ['tileid', 'mjd', 'expnum', 'filename', 'skye_name', 'skyw_name']
    problem_table = rows[bad][out_cols]
    problem_table['east_true']    = east_display[bad]
    problem_table['west_true']    = west_display[bad]
    problem_table['problem_type'] = problem_type[bad]

    # astrometry-quality flags -- only present when reading a
    # SummarizeSkyHdr.py file (see load_drpall_rows)
    for col in ('sci_asrc', 'skye_asrc', 'skyw_asrc'):
        if col in rows.colnames:
            problem_table[col] = np.array(rows[col])[bad]

    _type_order = {'Swapped': 0, 'HalfMatch': 1, 'Unknown': 2}
    problem_table['_order'] = [_type_order[t] for t in problem_table['problem_type']]
    problem_table.sort(['_order', 'mjd'])
    problem_table.remove_column('_order')

    return problem_table, len(rows)


def _print_summary(n_checked, tol, postype, problem_table, n_swapped, n_halfmatch, n_unknown):
    print('Checked %d exposures with both skye/skyw labeled '
          '(tolerance %.3f deg, postype=%s)' % (n_checked, tol, postype))
    print('%d problem exposures: %d Swapped, %d HalfMatch, %d Unknown'
          % (len(problem_table), n_swapped, n_halfmatch, n_unknown))


def steer(argv):
    drpall  = 'drpall-1.2.1.fits'
    csv     = 'final_sky_tiles.csv'
    tol     = 0.1
    outfile = ''
    postype = 'adopted'

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-csv':
            i += 1
            csv = argv[i]
        elif argv[i] == '-tol':
            i += 1
            tol = float(argv[i])
        elif argv[i] == '-out':
            i += 1
            outfile = argv[i]
        elif argv[i] == '-postype':
            i += 1
            postype = argv[i]
        elif argv[i].startswith('-'):
            print('Error: unknown option "%s"' % argv[i])
            print(_USAGE)
            return
        else:
            drpall = argv[i]
        i += 1

    if postype not in _POSTYPES:
        print('Error: -postype must be one of %s' % sorted(_POSTYPES))
        return

    drpall_stem = os.path.splitext(os.path.basename(drpall))[0]
    postype_suffix = '' if postype == 'adopted' else '_%s' % postype

    problem_outfile  = outfile if outfile else (
        'sky_problem_check_%s%s.tab' % (drpall_stem, postype_suffix))
    position_outfile = 'sky_position_check_%s%s.tab' % (drpall_stem, postype_suffix)

    # --- per-exposure problem listing ---------------------------------
    problem_table, n_checked = check_problems(drpall, csv, tol=tol, postype=postype)
    n_swapped   = int((problem_table['problem_type'] == 'Swapped').sum())
    n_halfmatch = int((problem_table['problem_type'] == 'HalfMatch').sum())
    n_unknown   = int((problem_table['problem_type'] == 'Unknown').sum())

    _print_summary(n_checked, tol, postype, problem_table, n_swapped, n_halfmatch, n_unknown)
    print()

    if len(problem_table) > 0:
        print(problem_table)
        problem_table.write(problem_outfile, format='ascii.fixed_width_two_line', overwrite=True)
        print('\nWrote %s' % problem_outfile)
    else:
        print('No problem exposures found; nothing written.')

    print('\n--- Summary ---')
    _print_summary(n_checked, tol, postype, problem_table, n_swapped, n_halfmatch, n_unknown)

    # --- per-name aggregate ---------------------------------------------
    print('\n' + '=' * 70)
    position_table, unmatched = check_positions(drpall, csv, tol=tol, postype=postype)

    print(position_table)
    position_table.write(position_outfile, format='ascii.fixed_width_two_line',
                         overwrite=True, formats={'mean_offset_unexplained_deg': '%.2f'})
    print('\nWrote %s' % position_outfile)

    if unmatched:
        print('\nWarning: %d label(s) used in %s do not match any name in %s:'
              % (len(unmatched), drpall, csv))
        for name in unmatched:
            print('    %s' % name)

    observed    = position_table[position_table['n_observed'] > 0]
    n_clean     = int((observed['n_incorrect'] == 0).sum())
    n_has_swap  = int((observed['n_swapped'] > 0).sum())
    n_has_half  = int((observed['n_halfmatch'] > 0).sum())
    n_has_unk   = int((observed['n_unexplained'] > 0).sum())

    print('\n--- Per-name summary ---')
    print('%d fields observed: %d clean, %d with a Swapped mismatch, '
          '%d with a HalfMatch mismatch, %d with an Unexplained mismatch'
          % (len(observed), n_clean, n_has_swap, n_has_half, n_has_unk))
    print('(these can overlap -- a field can have more than one kind '
          'across its observations)')


if __name__ == '__main__':
    if '-h' in sys.argv or '--help' in sys.argv:
        print(_USAGE)
    else:
        steer(sys.argv)
