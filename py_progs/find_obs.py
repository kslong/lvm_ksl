#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Cross-match LVM observations from a drpall file against a source
catalog to identify which exposures overlap each source position.


Command line usage::

    find_obs.py [-root outroot] [-max sep_arcsec] drpall_file source_catalog

The ``drpall_file`` argument is the LVM drpall FITS file containing one
row per exposure with pointing coordinates, exposure times, and tile IDs.

The ``source_catalog`` argument is an ASCII table with at least three
columns: ``RA`` (degrees), ``Dec`` (degrees), and ``Source_name``.

The ``-root outroot`` option sets the root name for output files.
The default is ``results``, giving ``results.matched.txt`` and
``results.sum.txt``.

The ``-max sep_arcsec`` option sets the cross-match radius in arcseconds.
The default is 900 arcseconds (15 arcminutes), which corresponds roughly
to the radius of one LVM pointing.


Description:

The routine reads the drpall FITS file and extracts the columns needed
for the cross-match: exposure number, MJD, exposure time, pointing RA
and Dec (``sci_ra``, ``sci_dec``), tile ID, and location.  It then reads
the source catalog, which must contain ``RA``, ``Dec``, and
``Source_name`` columns.

A spherical cross-match is performed using astropy SkyCoord, returning
all observation-source pairs within the specified match radius.  Each
matched pair records the observation metadata alongside the source
catalog columns and the angular separation in arcseconds.

A per-source summary table is also produced with the number of matching
exposures, the total integration time, and the median separation between
the source and the overlapping pointings.

Two output files are written in fixed-width ASCII format:

- ``<root>.matched.txt``: All matched pairs, one row per
  observation-source match, including all drpall columns and all
  source-catalog columns together with the angular separation.

- ``<root>.sum.txt``: One row per unique source, reporting the source
  name, the number of matching exposures (``Nobs``), the total
  accumulated exposure time in seconds (``TotExp``), and the median
  angular separation in arcseconds (``MedSep``).


Primary routines:

The main driver is ``steer``, which parses the command line and calls
``load_observations``, ``find_matches``, and ``summarize`` in sequence.


Notes:

The source catalog must be readable by ``astropy.io.ascii.read`` and
must contain columns named ``RA``, ``Dec``, and ``Source_name``.

The drpall file must be a FITS file with the observation table in
extension 1 and at least the columns ``expnum``, ``mjd``, ``exptime``,
``sci_ra``, ``sci_dec``, ``tileid``, and ``location``.

The routine will refuse to overwrite existing output files.  Use a
different ``-root`` name if you want to rerun with different parameters.

The cross-match uses ``search_around_sky``, which considers all pairs
within the match radius, not just the closest match.  A single exposure
can therefore appear multiple times if several sources lie within the
match radius of the same pointing.


History:

260312 ksl  Coding begun

'''

import sys
import os
import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord
from astropy import units as u


def load_observations(drpall_file):
    '''
    Read an LVM drpall FITS file and return a table of observations.

    Opens the FITS file, reads extension 1, and extracts the columns
    needed for position cross-matching: exposure number, MJD, exposure
    time, pointing RA and Dec, tile ID, and location.  The RA and Dec
    column names are normalised to ``RA`` and ``Dec`` respectively.
    A ``Source_name`` column is added containing the zero-padded
    exposure number, which is used later to group matched rows.

    Parameters:
        drpall_file (str): Path to the LVM drpall FITS file.

    Returns:
        astropy.table.Table: Table with columns expnum, mjd, exptime,
        RA, Dec, tileid, location, and Source_name.
    '''
    print('Reading drpall file: %s' % drpall_file)
    x = fits.open(drpall_file)
    xtab = Table(x[1].data)
    x.close()

    obs = xtab['expnum', 'mjd', 'exptime', 'sci_ra', 'sci_dec', 'tileid', 'location']
    obs.rename_column('sci_ra', 'RA')
    obs.rename_column('sci_dec', 'Dec')
    obs['Source_name'] = ['%05d' % row['expnum'] for row in obs]

    print('  Loaded %d observations' % len(obs))
    return obs


def find_matches(source_file, obs, max_sep_arcsec):
    '''
    Cross-match a source catalog against LVM observations.

    Reads the source catalog from ``source_file`` using
    ``astropy.io.ascii.read`` and performs an all-pairs spherical
    cross-match against the observation table ``obs`` using
    ``SkyCoord.search_around_sky``.  All pairs within ``max_sep_arcsec``
    arcseconds are returned.  If observation and source catalog share
    column names, the observation copy is suffixed with ``_obs`` to
    avoid clashes before the two tables are horizontally stacked.

    Parameters:
        source_file (str): Path to the ASCII source catalog.  Must
            contain columns ``RA`` (degrees), ``Dec`` (degrees), and
            ``Source_name``.
        obs (astropy.table.Table): Observation table from
            ``load_observations``.
        max_sep_arcsec (float): Cross-match radius in arcseconds.

    Returns:
        astropy.table.Table: Matched table with one row per
        observation-source pair.  Includes all columns from both
        the source catalog and the observation table, plus a
        ``separation`` column in arcseconds.  Returns an empty
        Table if no matches are found.
    '''
    print('Reading source catalog: %s' % source_file)
    sources = ascii.read(source_file)
    print('  Loaded %d sources' % len(sources))
    print('  Searching within %.1f arcsec ...' % max_sep_arcsec)

    coords_src = SkyCoord(ra=sources['RA'] * u.deg, dec=sources['Dec'] * u.deg)
    coords_obs = SkyCoord(ra=obs['RA'] * u.deg, dec=obs['Dec'] * u.deg)

    max_sep = max_sep_arcsec * u.arcsec
    idx_obs, idx_src, d2d, _ = coords_src.search_around_sky(coords_obs, max_sep)

    if len(idx_src) == 0:
        print('  No matches found.')
        return Table()

    matched_src = sources[idx_src]
    matched_obs = obs[idx_obs]

    src_cols = set(matched_src.colnames)
    for col in list(matched_obs.colnames):
        if col in src_cols:
            matched_obs.rename_column(col, col + '_obs')

    matches = hstack([matched_src, matched_obs])
    matches['separation'] = d2d.to(u.arcsec)
    matches['separation'].format = '.2f'

    print('  Found %d matches' % len(matches))
    return matches


def summarize(matches):
    '''
    Compute a per-source summary from a table of matched pairs.

    Groups the matched pairs by ``Source_name`` and for each unique
    source computes the number of overlapping exposures, the total
    accumulated exposure time, and the median angular separation between
    the source and the matched pointing centres.

    Parameters:
        matches (astropy.table.Table): Matched pair table returned by
            ``find_matches``.  Must contain columns ``Source_name``,
            ``exptime``, and ``separation``.

    Returns:
        astropy.table.Table: Summary table with columns ``Source_name``,
        ``Nobs`` (int), ``TotExp`` (float, seconds), and
        ``MedSep`` (float, arcseconds), one row per unique source.
    '''
    sources = np.unique(matches['Source_name'])
    nobs, tot_exp, med_sep = [], [], []
    for name in sources:
        rows = matches[matches['Source_name'] == name]
        nobs.append(len(rows))
        tot_exp.append(float(np.sum(rows['exptime'])))
        med_sep.append(float(np.median(rows['separation'])))

    summary = Table([sources, nobs, tot_exp, med_sep],
                    names=['Source_name', 'Nobs', 'TotExp', 'MedSep'])
    summary['MedSep'].format = '.1f'
    return summary


def steer(argv):
    '''
    Parse command line arguments and run the observation cross-match.

    Accepts the drpall FITS file and source catalog as positional
    arguments, with optional switches for the output root name and
    the cross-match radius.  Calls ``load_observations``,
    ``find_matches``, and ``summarize`` in sequence and writes the
    results to two fixed-width ASCII files.

    Parameters:
        argv (list): Command line argument list (typically ``sys.argv``).

    Returns:
        None
    '''
    drpall_file = ''
    source_file = ''
    root = 'results'
    max_sep = 900.0
    positional = []

    i = 1
    while i < len(argv):
        if argv[i][0:2] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-root':
            i += 1
            root = argv[i]
        elif argv[i] == '-max':
            i += 1
            max_sep = float(argv[i])
        elif argv[i][0] == '-':
            print('Unknown switch: %s' % argv[i])
            return
        else:
            positional.append(argv[i])
        i += 1

    if len(positional) < 2:
        print('Usage: find_obs.py [-root outroot] [-max sep_arcsec] drpall_file source_catalog')
        return

    drpall_file = positional[0]
    source_file = positional[1]

    obs_output = root + '.matched.txt'
    sum_output = root + '.sum.txt'

    for path in (obs_output, sum_output):
        if os.path.exists(path):
            print('Error: output file already exists: %s' % path)
            print('Use -root to choose a different root name.')
            return

    if not os.path.exists(drpall_file):
        print('Error: drpall file not found: %s' % drpall_file)
        return
    if not os.path.exists(source_file):
        print('Error: source catalog not found: %s' % source_file)
        return

    obs = load_observations(drpall_file)
    matches = find_matches(source_file, obs, max_sep)

    if len(matches) == 0:
        print('No matches found -- no output written.')
        return

    matches.write(obs_output, format='ascii.fixed_width_two_line', overwrite=False)
    print('Matches written to:  %s' % obs_output)

    summary = summarize(matches)
    summary.write(sum_output, format='ascii.fixed_width_two_line', overwrite=False)
    print('Summary written to:  %s  (%d unique sources)' % (sum_output, len(summary)))


# Next lines permit one to run the routine from the command line
if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
