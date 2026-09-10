#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Decompose one or more of FLUX (science fiber), SKY_EAST, and SKY_WEST
    from an LVM XCframe summary FITS file into a nebula-free "clean sky"
    model, using sky_decomp.lsf_surface_iterative.
    SkyDecompLSFSurfaceIterative with the known nebular emission lines
    excluded from the fit (mask-and-wrap).

Command line usage (if any):

    usage: DecomposeCleanSky.py [-h] [-lvmsky_skysub PATH] [-py_progs_dir PATH]
                                [-row N | -expnum N] [-ext LIST]
                                [-v VEL] [-lmc] [-smc] [-output PATH]
                                fits_file

    where

    fits_file       LVM XCframe summary FITS file (WAVE, FLUX, SKY_EAST,
                    SKY_WEST extensions and a DRP_ALL table, as produced by
                    SummarizeCframe.py -by fiber).

    -row N          selects the exposure by row index (default: 0).

    -expnum N       selects the exposure by DRP_ALL 'expnum' instead of a
                    row index (overrides -row).

    -ext LIST       comma-separated extensions to decompose (default:
                    FLUX,SKY_EAST,SKY_WEST -- SKY_EAST/SKY_WEST are the
                    DRP's actual sky-subtraction inputs; FLUX (the science
                    telescope) is included by default too so its nebular
                    lines are available as a same-exposure comparison
                    point for how much flux the sky telescopes see).

    -v VEL          nebular systemic velocity (km/s) used to Doppler-shift
                    the exclusion windows (default: 0).

    -lmc / -smc     shortcuts for the LMC (~262 km/s) / SMC (~146 km/s)
                    systemic velocity, same convention as sky_gaussfit.py.
                    -v overrides these if both are given.

    -output PATH    output FITS path (default: CleanSky_<expnum>.fits).

    -lvmsky_skysub PATH
                    path to the lvmsky repo's skysub/ directory, which
                    supplies the sky_decomp package this script imports
                    (default: ~/SDSS/lvmsky/skysub).

    -py_progs_dir PATH
                    path to the lvm_ksl repo's py_progs/ directory, which
                    supplies sky_gaussfit.py's NEBULAR_LINES list (default:
                    ~/SDSS/lvm_ksl/py_progs).

Description:

    For each requested extension, this script:

      1. Builds the nebular exclusion mask from sky_gaussfit.NEBULAR_LINES
         (9 lines: Halpha, [OIII] a/b, [OI] a/b, [NII] a/b, [SII] a/b),
         each window Doppler-shifted by the same zz = 1 + vel/3e5
         convention sky_gaussfit.py itself uses, so a given -v/-lmc/-smc
         produces identical windows to that script's own nebular-line fits.
      2. Zeroes IVAR inside those windows (and at any non-finite flux
         pixel) before calling SkyDecompLSFSurfaceIterative.fit() -- the
         fit engine's own good = isfinite(flux) & isfinite(ivar) & (ivar>0)
         then simply never sees the nebular pixels, so none of its
         continuum/sky-line families can absorb nebular flux.
      3. Reconstructs the nebula-free "clean sky" model
         (result.bestfit_lsf) and the full-array residual
         (observed - clean sky, including *inside* the excluded windows --
         that residual there is the nebular-leak signal a companion tool,
         py_progs/sky_nebular_leak_eval.py, measures).

    This is the "mask-and-wrap" half of the continuum/sky-line/nebular-line
    separation effort. A later, more ambitious approach -- a native
    nebular family solved jointly inside SkyDecomp's own design matrix,
    instead of excluding nebular pixels beforehand -- would slot in as a
    different model-building step, without changing the leak-detection
    tooling downstream (sky_nebular_leak_eval.py only ever needs a
    (wave, flux, model) triple).

Primary routines:

    nebular_mask     build the Doppler-shifted exclusion mask.
    build_decomp     construct the SkyDecompLSFSurfaceIterative instance.
    decompose_clean  mask, fit, and reconstruct one extension's clean-sky model.
    write_output     write observed/model/residual/mask/components to FITS.

Notes:

    - The input file has no IVAR extension (same as PredictSky.py), so the
      pre-mask baseline is a uniform ivar = 1, not true photon-noise
      weighting.
    - Uses the plain (non-split-zodi, non-Moon/Zodi-physical) fit-model
      branch of SkyDecompLSFSurfaceIterative -- the same "lsf-surface-
      iterative" mode decompose_parallel.py runs at corpus scale -- not
      the Moon+Zodi geometry-driven variant, which needs per-observation
      ephemeris inputs this script doesn't otherwise use.
    - Imports sky_decomp from an external lvmsky checkout (-lvmsky_skysub)
      and sky_gaussfit.NEBULAR_LINES from this repo's own py_progs/
      (-py_progs_dir) via sys.path; both must resolve for this script to
      run.

History::

    260909  ksl  Coding begun.

'''

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

# ---------------------------------------------------------------------------
# External package setup -- sky_decomp lives in the lvmsky repo, NEBULAR_LINES
# in this repo's py_progs/sky_gaussfit.py.  Inserted before the imports below
# so -lvmsky_skysub / -py_progs_dir can override the defaults first.
# ---------------------------------------------------------------------------

DEFAULT_LVMSKY_SKYSUB = Path('~/SDSS/lvmsky/skysub').expanduser()
DEFAULT_PY_PROGS_DIR = Path('~/SDSS/lvm_ksl/py_progs').expanduser()

_pre = argparse.ArgumentParser(add_help=False)
_pre.add_argument('-lvmsky_skysub', default=str(DEFAULT_LVMSKY_SKYSUB))
_pre.add_argument('-py_progs_dir', default=str(DEFAULT_PY_PROGS_DIR))
_pre_args, _ = _pre.parse_known_args()

# Insert py_progs_dir first, lvmsky_skysub second, so lvmsky_skysub ends up
# at sys.path[0] -- otherwise py_progs/sky_decomp/ (a stale, partial vendor
# copy kept only for PALACE scripts, see py_progs/palace_make_mask.py) would
# shadow the current lvmsky sky_decomp package and break this import with a
# confusing "no submodule lsf_surface_iterative" error.
sys.path.insert(0, _pre_args.py_progs_dir)
sys.path.insert(0, _pre_args.lvmsky_skysub)

from sky_decomp.lsf_surface_iterative import (  # noqa: E402
    SkyDecompLSFSurfaceIterative, LSFSurfaceIterativeConfig,
)
from sky_gaussfit import resolve_nebular_lines  # noqa: E402

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

FACTOR = 1e14  # XCframe flux is erg/s/cm^2/A; the QP fit works in O(1) counts.
MOON_N_KNOTS_DEFAULT = 25  # matches SkyDecomp.__init__'s own default (fit.py).
LMC_VEL = 262.
SMC_VEL = 146.


def resolve_velocity(vel=None, lmc=False, smc=False):
    '''
    Resolve the nebular systemic velocity from -v/-lmc/-smc.

    Same precedence as sky_gaussfit.py's steer(): an explicit -v wins,
    then -lmc, then -smc, else 0.
    '''
    if vel is not None:
        return float(vel)
    if lmc:
        return LMC_VEL
    if smc:
        return SMC_VEL
    return 0.0


def nebular_mask(wave, vel=0.0):
    '''
    Boolean mask, True where a pixel is clean of every genuinely-nebular
    line's Doppler-shifted window (sky_gaussfit.resolve_nebular_lines).

    Lines whose shifted window still overlaps a sky_gaussfit.SKY_LINES
    window (oi_a/oi_b at low velocity -- the same [OI] 6300/6364
    transition as sky6300/sky6363, predominantly sky airglow rather than
    nebular unless the systemic velocity is large enough to separate
    them) are left unmasked, so SkyDecomp's own atom/oh families fit that
    sky-line flux normally instead of losing it to an unconstrained gap.

    Parameters
    ----------
    wave : ndarray
        Wavelength grid (Angstrom).
    vel : float
        Nebular systemic velocity (km/s); windows are shifted by
        zz = 1 + vel/3e5, matching sky_gaussfit.do_one's convention.

    Returns
    -------
    ndarray of bool
        True = clean/keep, False = inside a nebular-line window.
    '''
    wave = np.asarray(wave, dtype=float)
    zz = 1.0 + vel / 3e5
    resolved, dropped = resolve_nebular_lines(vel)
    if dropped:
        print(f"  nebular_mask: treating as sky, not nebular (overlaps a "
              f"SKY_LINES window at vel={vel:.1f}): {', '.join(dropped)}")
    clean = np.ones(wave.shape, dtype=bool)
    for name, center, wmin, wmax in resolved:
        clean &= ~((wave >= zz * wmin) & (wave <= zz * wmax))
    return clean


def read_row(fits_file, row=None, expnum=None, exts=('SKY_EAST', 'SKY_WEST')):
    '''
    Read one exposure's wavelength grid, requested flux extensions, and
    identifying metadata from an XCframe summary file.

    Parameters
    ----------
    fits_file : str or Path
    row : int, optional
    expnum : int, optional
    exts : sequence of str
        Extension names to read (e.g. 'FLUX', 'SKY_EAST', 'SKY_WEST').

    Returns
    -------
    dict
        Keys: row, expnum, mjd, wave, flux (dict ext -> ndarray).
    '''
    with fits.open(fits_file, memmap=True) as hdul:
        drp = Table(hdul['DRP_ALL'].data)
        if expnum is not None:
            matches = np.flatnonzero(np.asarray(drp['expnum']) == int(expnum))
            if matches.size == 0:
                raise ValueError(f"expnum {expnum} not found in {fits_file}")
            i = int(matches[0])
        else:
            i = int(row) if row is not None else 0

        wave = np.asarray(hdul['WAVE'].data, dtype=np.float64)
        flux = {ext: np.asarray(hdul[ext].data[i], dtype=np.float64) for ext in exts}
        drp_row = drp[i]

    mjd = Time(str(drp_row['obstime']), format='isot', scale='utc').mjd
    return dict(row=i, expnum=int(drp_row['expnum']), mjd=mjd, wave=wave, flux=flux)


def build_decomp(wave, lvmsky_skysub):
    '''
    Construct the SkyDecompLSFSurfaceIterative instance decompose_clean
    needs, following decompose_parallel.py's "lsf-surface-iterative"
    fit-model settings (validated at corpus scale) rather than inventing
    new knob values.
    '''
    base_dir = Path(lvmsky_skysub) / 'sky_decomp' / 'data'
    return SkyDecompLSFSurfaceIterative(
        wave, lsf_sigma=1.0, base_dir=base_dir,
        moon_smooth_lambda=0.1, moon_interline_boost=0.0,
        n_spline_knots=MOON_N_KNOTS_DEFAULT,
        config=LSFSurfaceIterativeConfig(n_refinement_cycles=5),
    )


def decompose_clean(decomp, flux_row, clean_mask, label):
    '''
    Fit one spectrum with the nebular windows excluded.

    Parameters
    ----------
    decomp : SkyDecompLSFSurfaceIterative
    flux_row : ndarray
        Observed flux, native XCframe units.
    clean_mask : ndarray of bool
        True = keep in the fit (see nebular_mask).
    label : str
        Extension name, used only for the printed summary.

    Returns
    -------
    dict
        Keys: fit (SkyDecompResult), bestfit (ndarray, physical units),
        resid (ndarray, physical units, full array including the masked
        nebular windows), components (dict of ndarray, physical units).
    '''
    flux = flux_row * FACTOR
    ivar = np.where(clean_mask & np.isfinite(flux), 1.0, 0.0)
    fit = decomp.fit(flux, ivar, verbose=False)
    print(f"  {label} decomp: status={fit.fit_status}, "
          f"chi2_red={fit.reduced_chi2:.2f}, R^2={fit.r2:.4f}, "
          f"n_masked={int((~clean_mask).sum())}")

    bestfit = fit.bestfit_lsf / FACTOR
    resid = flux_row - bestfit
    components = {k: np.asarray(v, dtype=np.float64) / FACTOR
                 for k, v in fit.components.items()}
    return dict(fit=fit, bestfit=bestfit, resid=resid, components=components)


def write_output(row_data, results, clean_mask, fits_file, vel, outpath=None):
    '''
    Write observed/model/residual/mask/component extensions for every
    decomposed extension to one output FITS file.

    Parameters
    ----------
    row_data : dict
        Output of read_row.
    results : dict
        Extension name -> output of decompose_clean.
    clean_mask : ndarray of bool
        The mask used (shared across extensions -- built once from wave).
    fits_file : str or Path
        Source file, recorded for provenance.
    vel : float
        Nebular velocity used, recorded for provenance.
    outpath : str or Path, optional
        Default: CleanSky_<expnum>.fits.

    Returns
    -------
    Path
        The path written.
    '''
    if outpath is None:
        outpath = f"CleanSky_{row_data['expnum']}.fits"
    outpath = Path(outpath)

    hdr = fits.Header()
    hdr['TITLE'] = 'CleanSky'
    hdr['INPUT'] = (str(fits_file), 'source XCframe summary file')
    hdr['ROW'] = (row_data['row'], 'row index in source file')
    hdr['EXPNUM'] = (row_data['expnum'], 'exposure number')
    hdr['MJD'] = (row_data['mjd'], 'UT MJD of exposure')
    hdr['NEBVEL'] = (vel, 'nebular systemic velocity used (km/s)')

    hdus = [fits.PrimaryHDU(header=hdr),
            fits.ImageHDU(data=row_data['wave'].astype(np.float32), name='WAVE'),
            fits.ImageHDU(data=clean_mask.astype(np.uint8), name='NEBMASK')]
    for ext, res in results.items():
        hdus.append(fits.ImageHDU(data=row_data['flux'][ext].astype(np.float32), name=ext))
        hdus.append(fits.ImageHDU(data=res['bestfit'].astype(np.float32), name=f'{ext}_BESTFIT'))
        hdus.append(fits.ImageHDU(data=res['resid'].astype(np.float32), name=f'{ext}_RESID'))
        for key, comp in res['components'].items():
            hdus.append(fits.ImageHDU(data=comp.astype(np.float32),
                                      name=f'{ext}_COMP_{key.upper()}'))

    fits.HDUList(hdus).writeto(outpath, overwrite=True)
    print(f"\nOutput written to {outpath}")
    return outpath


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        parents=[_pre],
        description=('Decompose FLUX/SKY_EAST/SKY_WEST into a '
                     'nebula-free clean-sky model, masking known nebular '
                     'lines out of the SkyDecomp fit.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('-row', type=int, default=None,
                   help='Row index to select (default: 0)')
    p.add_argument('-expnum', type=int, default=None,
                   help='Select by DRP_ALL expnum instead of row index')
    p.add_argument('-ext', default='FLUX,SKY_EAST,SKY_WEST',
                   help='Comma-separated extensions to decompose')
    p.add_argument('-v', dest='vel', type=float, default=None,
                   help='Nebular systemic velocity (km/s)')
    p.add_argument('-lmc', action='store_true',
                   help=f'Use the LMC velocity (~{LMC_VEL:.0f} km/s)')
    p.add_argument('-smc', action='store_true',
                   help=f'Use the SMC velocity (~{SMC_VEL:.0f} km/s)')
    p.add_argument('-output', default=None,
                   help='Output FITS path (default: CleanSky_<expnum>.fits)')
    p.add_argument('fits_file', help='LVM XCframe summary FITS file')
    args = p.parse_args()

    vel = resolve_velocity(args.vel, args.lmc, args.smc)
    exts = [e.strip() for e in args.ext.split(',') if e.strip()]

    row_data = read_row(args.fits_file, row=args.row, expnum=args.expnum, exts=exts)
    print(f"Source: {args.fits_file}")
    print(f"  row={row_data['row']}  expnum={row_data['expnum']}  "
          f"mjd={row_data['mjd']:.5f}  vel={vel:.1f} km/s")

    clean_mask = nebular_mask(row_data['wave'], vel=vel)
    print(f"  nebular mask: {int((~clean_mask).sum())} / {clean_mask.size} "
          f"pixels excluded")

    decomp = build_decomp(row_data['wave'], _pre_args.lvmsky_skysub)
    results = {}
    for ext in exts:
        results[ext] = decompose_clean(decomp, row_data['flux'][ext], clean_mask, ext)

    write_output(row_data, results, clean_mask, args.fits_file, vel, outpath=args.output)


if __name__ == '__main__':
    main()
