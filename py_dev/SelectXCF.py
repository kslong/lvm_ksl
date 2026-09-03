#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Select a criteria-filtered, moon-geometry-stratified random subset of
    rows from an LVM XCframe summary FITS file and write them out in the
    identical XCframe layout (same WAVE/FLUX/SKY_EAST/SKY_WEST/LSF/DRP_ALL
    extensions, just fewer rows).

    Intended as the first step of building a small validation corpus for
    retraining the mlp_ensemble_split_zodi sky-prediction model: this script
    only selects and repackages rows, it does not reformat anything for
    lvmsky's decompose_parallel.py (that is a separate conversion step).

Command line usage (if any):

    usage: SelectXCF.py [-h] [-n N] [-seed SEED] [-nbins NBINS]
                        [-exptime EXPTIME] [-fluxcal FLUXCAL]
                        [-min_glat MIN_GLAT]
                        [-lmc_ra LMC_RA] [-lmc_dec LMC_DEC] [-lmc_radius LMC_RADIUS]
                        [-smc_ra SMC_RA] [-smc_dec SMC_DEC] [-smc_radius SMC_RADIUS]
                        [-orion_ra ORION_RA] [-orion_dec ORION_DEC] [-orion_radius ORION_RADIUS]
                        fits_file [output_file]

    where

    fits_file       is the path to an LVM XCframe summary FITS file (any
                    percentile/by-fiber-or-pixel variant produced by
                    SummarizeCframe.py -- all share the same
                    WAVE/FLUX/SKY_EAST/SKY_WEST/LSF/DRP_ALL layout).

    output_file     output FITS path (default: <stem>_sel<N>.fits in the
                    current working directory, regardless of where
                    fits_file itself lives).

    -n N            number of rows to select (default: 100).

    -seed SEED      random seed for reproducible sampling (default: 42).

    -nbins NBINS    number of quantile bins per axis (moon_alt, moon_fli)
                    for the stratified draw (default: 4, giving up to
                    NBINS**2 strata).

    -exptime EXPTIME
                    required DRP_ALL 'exptime' value in seconds (default:
                    900.0).

    -fluxcal FLUXCAL
                    required DRP_ALL 'fluxcal' value (default: 'MOD').

    -min_glat MIN_GLAT
                    minimum |galactic latitude| in degrees for the Sci
                    pointing (default: 10.0).

    -lmc_ra/-lmc_dec/-lmc_radius
                    LMC exclusion center + radius in degrees (default:
                    80.8942, -69.7561, 6.0 -- py_progs/rss2image.py's
                    values).

    -smc_ra/-smc_dec/-smc_radius
                    SMC exclusion center + radius in degrees (default:
                    13.1583, -72.8003, 3.0 -- py_progs/rss2image.py's
                    values).

    -orion_ra/-orion_dec/-orion_radius
                    Orion Nebula exclusion center + radius in degrees
                    (default: 83.8221, -5.3911, 15.0).

Description:

    1. Reads DRP_ALL from the input file and applies the hard cuts (exptime,
       fluxcal, galactic latitude, LMC/SMC/Orion exclusion) using the Sci
       telescope's pointing (sci_ra/sci_dec).
    2. From the surviving candidate rows, draws -n rows via a stratified
       random sample over a (moon_alt x moon_fli) quantile grid, so the
       small output corpus spans the moon-geometry range the sky-prediction
       model conditions on, rather than clustering at whatever moon state
       happened to be common in the input file.
    3. Writes the selected rows to a new FITS file with the same extension
       structure as the input (PRIMARY header copied, WAVE unchanged,
       FLUX/SKY_EAST/SKY_WEST/LSF and DRP_ALL sliced to the selected rows,
       row order sorted by original index for reproducible diffing).

Notes::

    This script only selects rows -- it does not rename/remap extensions
    for lvmsky's decompose_parallel.py (FLUX -> FLUX_SCI, per-row
    SKY_EAST/SKY_WEST -> FLUX_SKY_NEAR/FLUX_SKY_FAR via the Near/Far label,
    etc). That reformatting is a separate downstream step.

History::

    260901  ksl  Coding begun.
    260903  ksl  Default output_file now lands in the current working
        directory, not next to fits_file -- found running from a fresh
        directory (fits_file pointing at a shared source corpus
        elsewhere) that the old default silently wrote the selection
        back into that source directory instead of where the user was
        actually working.
    260903  ksl  Switched every option from double-dash (--min-glat) to
        single-dash (-min_glat), matching py_progs/'s convention -- see
        BatchPredictSkyESO.py's History for the fuller note.

'''

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u


# ---------------------------------------------------------------------------
# Defaults (py_progs/rss2image.py's LMC/SMC centers+radii; Orion Nebula
# coordinates and a 15 deg exclusion radius per this project's convention)
# ---------------------------------------------------------------------------

DEFAULT_N = 100
DEFAULT_SEED = 42
DEFAULT_NBINS = 4
DEFAULT_EXPTIME = 900.0
DEFAULT_FLUXCAL = 'MOD'
DEFAULT_MIN_GLAT = 10.0
DEFAULT_LMC = dict(ra=80.8942, dec=-69.7561, radius=6.0)
DEFAULT_SMC = dict(ra=13.1583, dec=-72.8003, radius=3.0)
DEFAULT_ORION = dict(ra=83.8221, dec=-5.3911, radius=15.0)


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------

def apply_hard_cuts(drp, args):
    '''
    Apply the exptime/fluxcal/galactic-latitude/LMC/SMC/Orion cuts to a
    DRP_ALL table.

    Parameters
    ----------
    drp : astropy.table.Table
        The DRP_ALL table (or a compatible structured array).
    args : argparse.Namespace
        Parsed command-line arguments carrying the cut parameters.

    Returns
    -------
    numpy.ndarray
        Boolean mask, True for rows passing every cut.
    '''
    sci = SkyCoord(ra=np.asarray(drp['sci_ra'], float) * u.deg,
                    dec=np.asarray(drp['sci_dec'], float) * u.deg,
                    frame='icrs')

    mask = np.isclose(drp['exptime'], args.exptime)
    mask &= (np.char.strip(np.asarray(drp['fluxcal'], str)) == args.fluxcal)
    mask &= np.abs(sci.galactic.b.deg) >= args.min_glat

    lmc = SkyCoord(ra=args.lmc_ra * u.deg, dec=args.lmc_dec * u.deg)
    smc = SkyCoord(ra=args.smc_ra * u.deg, dec=args.smc_dec * u.deg)
    orion = SkyCoord(ra=args.orion_ra * u.deg, dec=args.orion_dec * u.deg)

    mask &= sci.separation(lmc).deg >= args.lmc_radius
    mask &= sci.separation(smc).deg >= args.smc_radius
    mask &= sci.separation(orion).deg >= args.orion_radius

    return mask


def stratified_sample(moon_alt, moon_fli, n, nbins, seed):
    '''
    Draw n indices from range(len(moon_alt)), stratified over a
    (moon_alt x moon_fli) quantile grid.

    Parameters
    ----------
    moon_alt : numpy.ndarray
        Moon altitude (degrees) for each candidate row.
    moon_fli : numpy.ndarray
        Moon fractional illumination for each candidate row.
    n : int
        Number of indices to draw (capped at len(moon_alt)).
    nbins : int
        Quantile bins per axis; the grid has up to nbins**2 occupied cells.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    numpy.ndarray
        Sorted array of selected indices into the input arrays.
    '''
    rng = np.random.default_rng(seed)
    n_candidates = len(moon_alt)
    n = min(int(n), n_candidates)

    alt_edges = np.quantile(moon_alt, np.linspace(0, 1, nbins + 1))
    fli_edges = np.quantile(moon_fli, np.linspace(0, 1, nbins + 1))
    # searchsorted with these edges can put the max value one bin past the
    # last -- clip both axes into [0, nbins - 1].
    alt_bin = np.clip(np.searchsorted(alt_edges, moon_alt, side='right') - 1, 0, nbins - 1)
    fli_bin = np.clip(np.searchsorted(fli_edges, moon_fli, side='right') - 1, 0, nbins - 1)
    cell = alt_bin * nbins + fli_bin

    cells = {}
    for idx, c in enumerate(cell):
        cells.setdefault(c, []).append(idx)
    for c in cells:
        rng.shuffle(cells[c])

    order = list(cells.keys())
    rng.shuffle(order)

    selected = []
    while len(selected) < n:
        progressed = False
        for c in order:
            if cells[c]:
                selected.append(cells[c].pop())
                progressed = True
                if len(selected) == n:
                    break
        if not progressed:
            break

    return np.sort(np.array(selected, dtype=int))


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def select(fits_file, args):
    '''
    Apply the hard cuts and stratified draw to fits_file's DRP_ALL table.

    Returns
    -------
    numpy.ndarray
        Sorted array of selected row indices into the original file.
    '''
    with fits.open(fits_file, memmap=True) as hdul:
        drp = Table(hdul['DRP_ALL'].data)

    mask = apply_hard_cuts(drp, args)
    n_pass = int(mask.sum())
    print(f'{n_pass} / {len(drp)} rows pass the hard cuts '
          f'(exptime={args.exptime}, fluxcal={args.fluxcal!r}, '
          f'|b|>={args.min_glat} deg, LMC/SMC/Orion excluded)')
    if n_pass == 0:
        raise ValueError('No rows survive the selection criteria.')

    candidate_idx = np.flatnonzero(mask)
    sel = stratified_sample(
        np.asarray(drp['moon_alt'], float)[candidate_idx],
        np.asarray(drp['moon_fli'], float)[candidate_idx],
        n=args.n, nbins=args.nbins, seed=args.seed,
    )
    selected_idx = candidate_idx[sel]
    print(f'Selected {len(selected_idx)} rows '
          f'(stratified over moon_alt x moon_fli, {args.nbins}x{args.nbins} grid, '
          f'seed={args.seed})')
    return selected_idx


def write_subset(fits_file, selected_idx, outpath):
    '''
    Write the rows at selected_idx to outpath, preserving the input file's
    extension structure (PRIMARY header, WAVE unchanged, FLUX/SKY_EAST/
    SKY_WEST/LSF and DRP_ALL sliced to selected_idx).
    '''
    with fits.open(fits_file, memmap=True) as hdul:
        primary_hdr = hdul['PRIMARY'].header.copy()
        wave = hdul['WAVE'].data
        wave_hdr = hdul['WAVE'].header.copy()
        flux = hdul['FLUX'].data[selected_idx]
        sky_e = hdul['SKY_EAST'].data[selected_idx]
        sky_w = hdul['SKY_WEST'].data[selected_idx]
        lsf = hdul['LSF'].data[selected_idx]
        drp = Table(hdul['DRP_ALL'].data)[selected_idx]

    primary_hdr['SRCFILE'] = (str(Path(fits_file).name), 'Input XCframe file')
    primary_hdr['NSEL'] = (len(selected_idx), 'Rows selected by SelectXCF.py')

    hdul_out = fits.HDUList([
        fits.PrimaryHDU(header=primary_hdr),
        fits.ImageHDU(data=wave, header=wave_hdr, name='WAVE'),
        fits.ImageHDU(data=flux, name='FLUX'),
        fits.ImageHDU(data=sky_e, name='SKY_EAST'),
        fits.ImageHDU(data=sky_w, name='SKY_WEST'),
        fits.ImageHDU(data=lsf, name='LSF'),
        fits.BinTableHDU(drp, name='DRP_ALL'),
    ])
    hdul_out.writeto(outpath, overwrite=True)
    print(f'Wrote {len(selected_idx)} rows to {outpath}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=('Select a criteria-filtered, moon-geometry-stratified '
                     'random subset of rows from an XCframe file.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('fits_file', help='LVM XCframe summary FITS file')
    p.add_argument('output_file', nargs='?', default=None,
                   help='Output FITS path (default: <stem>_sel<N>.fits)')
    p.add_argument('-n', type=int, default=DEFAULT_N,
                   help='Number of rows to select')
    p.add_argument('-seed', type=int, default=DEFAULT_SEED,
                   help='Random seed for the stratified draw')
    p.add_argument('-nbins', type=int, default=DEFAULT_NBINS,
                   help='Quantile bins per axis for the moon_alt x moon_fli stratification')
    p.add_argument('-exptime', type=float, default=DEFAULT_EXPTIME,
                   help='Required DRP_ALL exptime value (seconds)')
    p.add_argument('-fluxcal', default=DEFAULT_FLUXCAL,
                   help='Required DRP_ALL fluxcal value')
    p.add_argument('-min_glat', type=float, default=DEFAULT_MIN_GLAT,
                   dest='min_glat',
                   help='Minimum |galactic latitude| in degrees for sci_ra/sci_dec')
    p.add_argument('-lmc_ra', type=float, default=DEFAULT_LMC['ra'], dest='lmc_ra')
    p.add_argument('-lmc_dec', type=float, default=DEFAULT_LMC['dec'], dest='lmc_dec')
    p.add_argument('-lmc_radius', type=float, default=DEFAULT_LMC['radius'], dest='lmc_radius')
    p.add_argument('-smc_ra', type=float, default=DEFAULT_SMC['ra'], dest='smc_ra')
    p.add_argument('-smc_dec', type=float, default=DEFAULT_SMC['dec'], dest='smc_dec')
    p.add_argument('-smc_radius', type=float, default=DEFAULT_SMC['radius'], dest='smc_radius')
    p.add_argument('-orion_ra', type=float, default=DEFAULT_ORION['ra'], dest='orion_ra')
    p.add_argument('-orion_dec', type=float, default=DEFAULT_ORION['dec'], dest='orion_dec')
    p.add_argument('-orion_radius', type=float, default=DEFAULT_ORION['radius'], dest='orion_radius')
    args = p.parse_args()

    outpath = args.output_file
    if outpath is None:
        stem = Path(args.fits_file).stem
        outpath = f'{stem}_sel{args.n}.fits'

    selected_idx = select(args.fits_file, args)
    write_subset(args.fits_file, selected_idx, outpath)


if __name__ == '__main__':
    main()
