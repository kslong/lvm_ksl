#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Combine per-exposure emission-line fit tables (e.g. the output of
lvm_gaussfit.py run on individual SFrame exposures) into a single
table of line fluxes covering the full region, weighted by both
spatial fiber-overlap and per-fit measurement uncertainty.

This is the fit-then-combine counterpart to rss_combine.py, which
combines the raw spectra first and fits the combined spectrum
afterwards. Because each exposure is fit independently before this
script runs, every measurement already carries its own error bar, so
the combine step can down-weight noisy or unreliable exposures
instead of averaging them in blind.


Command line usage::

    gauss_combine.py [-h] [-outroot xxxx] [-exclude f1,f2,...] filenames

    -outroot xxxx      Root name for the output table (default
                        'gauss_combine'); output is written to
                        '<outroot>.txt'.

    -exclude f1,f2,...  Comma-separated basenames of input files to
                        drop before combining (e.g. an exposure
                        already confirmed bad by OverlapFlux.py).

    filenames           One or more per-exposure line-fit tables
                        (ascii.fixed_width_two_line, one row per
                        science fiber), such as lvm_gaussfit.py's
                        *.gauss.txt output. Wildcards are expanded.


Description::

    The output fiber grid and the spatial fiber-overlap weights are
    built the same way rss_combine.py builds them for raw spectra:
    a regular grid is laid out over the sky area covered by all input
    fibers (create_wcs/generate_grid), and each input fiber's
    fractional area overlap with each output fiber (frac_calc2) is
    calculated from fiber positions alone -- both are imported
    unmodified from rss_combine.py rather than re-derived here.

    Unlike rss_combine.py, this script does not need a -orig/-sum
    fiber-placement option, disk-space pre-flight checks, or a
    memory-limited chunked read-back -- it is combining a couple of
    hundred scalar columns per fiber (already-fit line parameters),
    not a ~12000-wavelength-point spectrum, so everything fits in
    memory as astropy Tables.

    Line keys (e.g. 'ha', 'oii_a') are discovered from the input
    tables' flux_<key> columns rather than hardcoded, so this works
    against any fit table following that naming convention (e.g.
    lvm_gaussfit.py or lvm_snrfit.py output).

    For each line and each output fiber, every (exposure, input
    fiber) contributor whose aperture overlaps the output fiber is
    weighted by frac / eflux**2, where frac is the spatial overlap
    fraction and eflux is that contributor's own fitted flux
    uncertainty -- so a noisy contributor is naturally down-weighted
    (never hard-excluded by significance; see combine_family), on top
    of the spatial apportionment rss_combine.py already does for raw
    spectra. The same per-contributor weight is reused to combine
    wave_<key>, fwhm_<key>, and back_<key> (when present) for that
    line, since a bad flux fit generally means a bad line-center/
    width/background fit from the same contributor.

    A shared background column that isn't literally named back_<key>
    (e.g. oii_a/oii_b's shared back_oii_ab, per lvm_gaussfit.py) is
    not combined in this version -- only exact <base>_<key> columns
    are recognized.


Notes::

    This script does not attempt to correct the exposure-to-exposure
    systematic flux offsets found by OverlapFlux.py (real ~50-100%
    Hb/OII offsets between exposure groups in Vela, cause still
    unresolved per SenSwapCheck.py's negative result). Inverse-
    variance weighting will down-weight a noisy-bad exposure, but not
    a confident-but-wrong one. Instead, each output fiber/line reports
    n_<key> (contributor count with finite flux/eflux) and spread_<key>
    (unweighted std of those contributors' flux), so a bin built from
    disagreeing exposures is visible after the fact; a specific
    exposure confirmed bad can be dropped with -exclude.

    Combined fluxes and errors are not pre-filtered by significance
    (see combine_family), so, like rss_combine.py's own output, a
    combined flux_<key> can legitimately be small, zero, or negative
    for an undetected line -- use eflux_<key> (or flux_<key>/eflux_<key>)
    downstream to judge significance rather than assuming every
    non-NaN value is a detection.

    Line fluxes are combined in whatever units the input tables use
    (e.g. lvm_gaussfit.py's 1e-16 erg/s/cm**2 scaled units) -- this
    script does no unit conversion of its own.

    Only a regular sky grid (rss_combine.py's default, no -orig
    equivalent) is supported in this version.


Primary routines::

    get_size_from_tables  - bounding sky region from input fiber ra/dec
    combine_family         - weighted combine of one line for one output fiber
    do_combine             - main driver
    steer                  - command-line driver


History::

    260826 ksl Coding begun. Reuses rss_combine.py's create_wcs/
        generate_grid/frac_calc2 unmodified for the output grid and
        spatial fiber-overlap geometry; auto-discovers line keys from
        flux_<key> columns so it works against lvm_gaussfit.py or
        lvm_snrfit.py output unchanged. Verified end-to-end against
        real Vela data (21 exposures each from gauss_fit_individual/
        and sky_fit_individual/), cross-checked against rss_combine.py
        + lvm_gaussfit.py's combine-then-fit product on the same
        exposures (median flux agreement within ~1.5%).
    260826 ksl Fixed get_size_from_tables's padding: it computes the
        output grid's bounding box from actual per-fiber ra/dec (which
        already reaches the true footprint edges), but was padding by
        the same 0.25 deg rss_combine.py uses for its own get_size --
        appropriate there because that bounding box comes from one
        telescope-pointing coordinate per exposure and needs the pad to
        reach the fiber-bundle edges, but redundant here, where it
        double-counted that radius and produced an output grid padded
        with a large ring of empty (no-contributor) fibers (confirmed:
        1.744 deg actual vs 2.2 deg produced). Reduced rad to a small
        edge-clipping margin (0.02 deg); confirmed grid size now
        matches rss_combine.py's (1.75 deg) on the same 21 exposures.
    260826 ksl Removed combine_family's flux/eflux >= minsnr and
        flux > 0 requirements. Gating contributors on significance
        before combining introduces Eddington/selection bias: for a
        line near a single exposure's detection threshold, a
        significance cut only keeps the upward noise fluctuations,
        biasing the combined value high exactly where it matters most
        (confirmed against real data: at output fibers with only one
        surviving contributor, Hb/[OII]/[OIII] read 6-13% high relative
        to rss_combine.py's uncut combine-then-fit result, while Ha,
        never near threshold, matched to <0.6%; the cut was also
        dropping the affected bins outright -- Hb coverage was 78% of
        Ha's). Pure frac/eflux**2 weighting has no such bias and needs
        no separate protection against a catastrophic fit failure
        (e.g. an undetected line's fit returning eflux ~1e26) -- its
        weight is already ~0; only a finite-value/eflux>0 sanity check
        on the arithmetic remains. Matches rss_combine.py's own
        approach of never filtering individual pixels/exposures by
        significance.

'''

import os
import re
import sys
from glob import glob

import numpy as np
from astropy.io import ascii
from astropy.table import vstack

from rss_combine import create_wcs, generate_grid, frac_calc2


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before a line consisting of "History:" (or
    "History::"/"Version History" -- whitespace/colon-insensitive), so
    -h stays short even as that section grows.
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


_EXPNUM_RE = re.compile(r'(\d+)')


def parse_exposure(filename):
    '''
    Pull the exposure number out of a filename, e.g.
    lvmSFrame-00009083.gauss.txt -> 9083. Falls back to the bare
    filename (extension stripped) if no digit run is found.

    Parameters:
        filename (str): Path to a per-exposure line-fit table.

    Returns:
        int or str: The exposure number, or the basename if none could
            be parsed.
    '''
    base = os.path.basename(filename)
    m = _EXPNUM_RE.search(base)
    if m is None:
        return os.path.splitext(base)[0]
    return int(m.group(1))


def get_size_from_tables(tables, rad=0.02):
    '''
    Calculate the center and size of a sky region that would encompass
    all fibers in a list of per-exposure line-fit tables.

    Analogous to rss_combine.py's get_size, but reads fiber ra/dec
    columns directly from the tables instead of telescope-pointing
    header keywords from SFrame files. rss_combine.py's get_size
    computes its extent from one telescope-pointing coordinate per
    exposure and needs a large rad (0.25 deg default there) to pad out
    to the actual fiber-bundle edges. Here the extent is already
    computed from every individual fiber's own ra/dec, so it already
    includes the full bundle footprint -- rad only needs to be a small
    margin against edge-fiber clipping, not the bundle radius itself.
    Using rss_combine.py's 0.25 deg default here would double-count
    that radius and produce an output grid padded with a large ring of
    empty (no-contributor) fibers around the real data.

    Parameters:
        tables (list): List of astropy.table.Table, each with ra/dec columns in degrees.
        rad (float): Additional radius in degrees to add as padding around the computed extent. Default 0.02.

    Returns:
        tuple: (ra_center, dec_center, size) in degrees.
    '''
    ra = np.concatenate([np.asarray(t['ra'], dtype=float) for t in tables])
    dec = np.concatenate([np.asarray(t['dec'], dtype=float) for t in tables])

    ra_max = np.max(ra)
    ra_min = np.min(ra)
    dec_max = np.max(dec)
    dec_min = np.min(dec)

    ra_cen = 0.5 * (ra_max + ra_min)
    dec_cen = 0.5 * (dec_max + dec_min)

    delta_dec = dec_max - dec_min
    delta_ra = (ra_max - ra_min) * np.cos(dec_cen / 57.29578)
    size = np.max([delta_ra, delta_dec]) + 2 * rad

    print('RA: %10.6f Dec: %10.6f size: %.1f' % (ra_cen, dec_cen, size))
    return ra_cen, dec_cen, size


def combine_family(sub, key):
    '''
    Combine one line's per-exposure, per-fiber fit values for the
    contributors that map into a single output fiber.

    Contributors are combined with weight frac / eflux**2, where frac
    is the spatial fiber-overlap fraction from frac_calc2 and eflux is
    each contributor's own fitted flux uncertainty -- no significance
    (flux/eflux) or sign (flux > 0) cut is applied first. Gating
    contributors on significance before combining would bias the
    result: for a line near a single exposure's detection threshold,
    noise scatters roughly symmetrically around the true (low) flux,
    but a flux/eflux >= minsnr cut only keeps the upward fluctuations,
    so the combined value would be systematically biased high exactly
    where it matters most. Pure inverse-variance weighting has no such
    bias and needs no separate protection against a catastrophic fit
    failure (e.g. an undetected line's fit returning a huge eflux) --
    its weight is already ~0. This mirrors rss_combine.py's own
    approach of never filtering individual pixels/exposures by
    significance, just co-adding everything and letting the combined
    result's own noise level speak for itself. The same weight and
    surviving-contributor mask are reused for wave_<key>, fwhm_<key>,
    and back_<key> when those columns are present.

    Parameters:
        sub (astropy.table.Table): Rows of the vstacked fractional-contribution table that map to one output fiber.
        key (str): Line key, e.g. 'ha' or 'oii_a'.

    Returns:
        dict: Column name -> combined value for this key, for one output fiber.
    '''
    vcol = 'flux_' + key
    ecol = 'eflux_' + key
    out = {}

    flux = np.asarray(sub[vcol], dtype=float)
    eflux = np.asarray(sub[ecol], dtype=float)
    frac = np.asarray(sub['frac'], dtype=float)

    good = np.isfinite(flux) & np.isfinite(eflux) & (eflux > 0)

    n = int(np.sum(good))
    out['n_' + key] = n

    if n == 0:
        out[vcol] = np.nan
        out[ecol] = np.nan
        out['spread_' + key] = np.nan
        for base in ('wave', 'fwhm', 'back'):
            bcol = base + '_' + key
            if bcol in sub.colnames:
                out[bcol] = np.nan
                out['e' + bcol] = np.nan
        return out

    w = frac[good] / eflux[good] ** 2
    wsum = np.sum(w)

    out[vcol] = np.sum(w * flux[good]) / wsum
    out[ecol] = 1.0 / np.sqrt(wsum)
    out['spread_' + key] = float(np.std(flux[good])) if n > 1 else 0.0

    for base in ('wave', 'fwhm', 'back'):
        bcol = base + '_' + key
        ebcol = 'e' + bcol
        if bcol in sub.colnames and ebcol in sub.colnames:
            vals = np.asarray(sub[bcol], dtype=float)[good]
            valid = np.isfinite(vals)
            if np.any(valid):
                wv = w[valid]
                out[bcol] = np.sum(wv * vals[valid]) / np.sum(wv)
                out[ebcol] = 1.0 / np.sqrt(np.sum(wv))
            else:
                out[bcol] = np.nan
                out[ebcol] = np.nan

    return out


def do_combine(filenames, outroot='', spacing_arcsec=35., dd=35.):
    '''
    Main routine that combines per-exposure line-fit tables into one
    output table covering the full region.

    Steps: (1) read all input tables, (2) build a regular output fiber
    grid covering their combined sky footprint (get_size_from_tables,
    create_wcs, generate_grid, reusing rss_combine.py's geometry),
    (3) compute each input fiber's fractional overlap with the output
    grid (frac_calc2), (4) for each line and each output fiber,
    combine all contributors with frac/eflux**2 weighting, no
    significance cut (combine_family), (5) write the result.

    Parameters:
        filenames (list): Per-exposure line-fit table filenames.
        outroot (str): Root name for the output table. Default writes 'gauss_combine.txt'.
        spacing_arcsec (float): Output grid spacing in arcsec. Default 35 (LVM fiber diameter).
        dd (float): Fiber diameter in pixels used for the overlap-area calculation. Default 35.

    Returns:
        astropy.table.Table or None: The combined output table, or None if nothing could be combined.
    '''
    tables = []
    for fname in filenames:
        tab = ascii.read(fname, format='fixed_width_two_line')
        print('Read %-45s %5d fibers (exposure %s)' % (fname, len(tab), parse_exposure(fname)))
        tables.append(tab)

    ra_center, dec_center, size = get_size_from_tables(tables)
    wcs = create_wcs(ra_center, dec_center, 0., size)
    new_slitmap_table = generate_grid(wcs, spacing_arcsec)

    print('\nApportioning fractional contributions from individual fibers to final virtual fiber positions')
    zslit = []
    for tab in tables:
        tab = tab.copy()
        x_pixels, y_pixels = wcs.world_to_pixel_values(tab['ra'], tab['dec'])
        tab['X'] = x_pixels
        tab['Y'] = y_pixels
        one_tab = frac_calc2(new_slitmap_table, tab, dd=dd)
        if len(one_tab) > 0:
            zslit.append(one_tab)

    if len(zslit) == 0:
        print('Error: do_combine: no fibers could be apportioned onto the output grid')
        return None

    all_frac = vstack(zslit)

    lines = sorted({c[len('flux_'):] for c in all_frac.colnames if c.startswith('flux_')})
    print('\nDiscovered %d lines: %s' % (len(lines), ', '.join(lines)))

    results = {}
    groups = all_frac.group_by('fib_master')
    for sub, fid in zip(groups.groups, groups.groups.keys['fib_master']):
        row = {}
        for key in lines:
            row.update(combine_family(sub, key))
        results[int(fid)] = row

    out = new_slitmap_table['fiberid', 'X', 'Y', 'ra', 'dec'].copy()

    ncols = {}
    for key in lines:
        ncols['flux_' + key] = np.full(len(out), np.nan)
        ncols['eflux_' + key] = np.full(len(out), np.nan)
        ncols['n_' + key] = np.zeros(len(out), dtype=int)
        ncols['spread_' + key] = np.full(len(out), np.nan)
        for base in ('wave', 'fwhm', 'back'):
            bcol = base + '_' + key
            if bcol in all_frac.colnames:
                ncols[bcol] = np.full(len(out), np.nan)
                ncols['e' + bcol] = np.full(len(out), np.nan)

    for i, fid in enumerate(out['fiberid']):
        row = results.get(int(fid))
        if row is None:
            continue
        for colname, val in row.items():
            ncols[colname][i] = val

    for colname, arr in ncols.items():
        out[colname] = arr

    n_filled = int(np.sum(out['n_' + lines[0]] > 0)) if lines else 0
    print('\nOf %d output fibers, %d have at least one surviving contributor for %s' %
          (len(out), n_filled, lines[0] if lines else '(no lines)'))

    if outroot == '':
        outroot = 'gauss_combine'
    outname = '%s.txt' % outroot
    out.write(outname, format='ascii.fixed_width_two_line', overwrite=True)
    print('Wrote combined line-flux table: %s' % outname)

    return out


def steer(argv):
    '''
    Parse command line arguments and run the combining process.

    Parameters:
        argv (list): Command line arguments (sys.argv).

    Returns:
        None: Calls do_combine with parsed arguments.
    '''
    xfiles = []
    outroot = ''
    exclude = []

    i = 1
    while i < len(argv):
        if argv[i][0:2] == '-h':
            print(_usage_from_doc(__doc__))
            return
        elif argv[i] == '-outroot':
            i += 1
            outroot = argv[i]
        elif argv[i] == '-exclude':
            i += 1
            exclude = argv[i].split(',')
        elif argv[i][0] == '-':
            print('Could not parse command line:', argv)
            return
        else:
            xfiles.append(argv[i])
        i += 1

    if outroot.startswith('-'):
        print('Error: outroot cannot start with "-": %s' % outroot)
        print('This may indicate a missing argument after -outroot')
        return

    files = []
    for one in xfiles:
        if one.count('*'):
            files.extend(glob(one))
        else:
            files.append(one)

    if exclude:
        before = len(files)
        files = [f for f in files if os.path.basename(f) not in exclude]
        print('Excluded %d of %d files per -exclude' % (before - len(files), before))

    if len(files) == 0:
        print('There are no files to combine:', argv)
        return

    files.sort()
    print('Combining %d per-exposure line-fit tables' % len(files))
    for one_file in files:
        print(' ', one_file)
    print('')

    do_combine(files, outroot=outroot)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
