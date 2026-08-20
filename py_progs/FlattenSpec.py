#!/usr/bin/env python
# coding: utf-8
'''
                    Space Telescope Science Institute

Synopsis:

    Separate the emission lines of an emission-line region from any
    underlying continuum, which usually arises from stellar contributions
    to the overall spectrum rather than from the nebular gas itself.  Fits
    and subtracts a smooth continuum from a single ascii spectrum table
    (WAVE, FLUX[, ERROR, MASK] columns, e.g. the *_ave_sum.txt files written
    by region-averaging scripts), writing the continuum and residual back
    into the same table.

Command line usage (if any):

    usage: FlattenSpec.py filename.txt [-mask mask.fits] [-kstep N]
                           [-niter N] [-nsigma N] [-out outfile]

    where

    filename.txt   ascii spectrum table with WAVE and FLUX columns; ERROR
                    and MASK columns are used if present.

    Options::

        -h            print this help and exit
        -mask FILE    palace_make_mask.py-format FITS mask (WAVE + MASK
                      extensions, MASK==1 clean / 0 sky-line-affected).
                      Default: the vendored data/sky_mask.fits (same
                      default PlotSpec.py/PlotSpecI.py use for their -mask
                      recolor option).
        -kstep N      B-spline knot spacing in Angstroms (default 100).
        -niter N      number of sigma-clip continuum-refit iterations,
                      rejecting bright emission-line pixels (default 3;
                      0 disables clipping).
        -nsigma N     rejection threshold in robust sigma above the local
                      residual median (default 4.0).
        -out FILE     output filename (default: <stem>_flat.txt).

Description:

    Pixels are excluded from the continuum fit if any of:

    * the mask file flags them sky-line-affected (interpolated onto the
      spectrum's own wavelength grid), since sky subtraction is known to
      leave residuals right where the airglow lines are;
    * the input table's own MASK column already flags them bad (MASK != 0,
      the usual DRP convention used throughout this codebase); or
    * an iterative sigma-clip (see below) flags them as bright
      emission-line pixels.

    The surviving pixels are fit with GetSkyCont.py's plain B-spline design
    matrix (build_design_matrix with no solar/MOON term -- this is a
    science spectrum, not a sky spectrum) via non-negative least squares
    (GetSkyCont.fit_continuum).  The continuum is evaluated at every
    wavelength, including masked regions, which the spline interpolates
    across.

    Because the B-spline coefficients are constrained non-negative, an
    unmasked bright science emission line (e.g. Halpha, [OIII] 5007) can
    only pull the local knot up, never down to compensate -- so the fit
    tracks the line as a continuum bump instead of ignoring it.  To avoid
    this without maintaining a fixed line list (this is a shock spectrum;
    lines can be Doppler-shifted/broadened off their rest wavelengths),
    the continuum is refit iteratively: after each fit, pixels whose
    residual (FLUX - CONT) exceeds ``nsigma`` robust-sigma (MAD-based)
    above the local residual median are excluded, and the continuum is
    refit without them.  Only positive outliers are rejected, since
    science emission adds flux rather than removing it.  Iteration stops
    early if a pass rejects no new pixels.

    Output columns (existing columns are kept; these are added/updated)::

        CONT   fitted continuum
        RESID  FLUX - CONT
        MASK   original MASK bits, OR'd with FLATTEN_SKYLINE_BIT wherever
               the sky-line mask excluded the pixel and/or
               FLATTEN_LINE_BIT wherever the sigma-clip excluded it as a
               bright emission-line pixel.  Pixels already bad in the
               input keep their original bits -- existing flag semantics
               are preserved, not overwritten.

Primary routines:

    flatten_spectrum   fit continuum/residual for one in-memory spectrum
    process_file       read a file, flatten it, write the result

Notes:

    This reuses GetSkyCont.py's fitting machinery directly
    (build_design_matrix/fit_continuum) rather than duplicating it, and
    reuses PlotSpec.py's get_sky_mask/DEFAULT_MASK_FILE for the default
    mask file, so all three scripts share one definition of "the" sky
    mask.  No solar/MOON component is added here, since a science
    spectrum's continuum has no reason to follow the solar Fraunhofer
    spectrum the way scattered moonlight does.

    MASK convention note: this script's output MASK column follows this
    codebase's usual DRP-style convention (0=good; see GetSpec.py/
    lvm_flux.py, which filter on MASK==0) -- the OPPOSITE sense from
    GetSkyCont.load_mask()'s own return value (True/1=clean), which
    describes the palace_mask FITS file itself, not a spectrum table.

History::

    260820 ksl Coding begun
    260820 ksl Added iterative sigma-clip rejection of bright
        emission-line pixels (FLATTEN_LINE_BIT), after finding the
        NNLS continuum tracking unmasked science lines as bumps
'''

import sys
from pathlib import Path

import numpy as np
from astropy.io import ascii

sys.path.insert(0, str(Path(__file__).resolve().parent))
from GetSkyCont import build_design_matrix, fit_continuum, _interp_mask_to_wave
from PlotSpec import get_sky_mask, DEFAULT_MASK_FILE, _usage_from_doc


_USAGE = _usage_from_doc(__doc__)

FLATTEN_SKYLINE_BIT = 1 << 20
FLATTEN_LINE_BIT    = 1 << 21


def flatten_spectrum(wave, flux, mask=None, sky_mask=None, knot_step=100.0,
                      niter=3, nsigma=4.0):
    '''
    Fit and subtract a continuum from one spectrum.

    Parameters
    ----------
    wave : ndarray
    flux : ndarray
    mask : ndarray of int or None
        Input MASK column (0=good); None if the spectrum has none.
    sky_mask : (ndarray, ndarray of bool) or None
        (mask_wave, mask_clean) as returned by PlotSpec.get_sky_mask /
        GetSkyCont.load_mask (True=clean); None to skip sky-line masking.
    knot_step : float
        B-spline knot spacing in Angstroms.
    niter : int
        Number of sigma-clip continuum-refit iterations; 0 disables
        clipping.
    nsigma : float
        Rejection threshold in robust (MAD-based) sigma above the local
        residual median.  Only positive outliers (emission) are rejected.

    Returns
    -------
    cont : ndarray
        Fitted continuum, evaluated at every wavelength.
    resid : ndarray
        flux - cont.
    fit_mask : ndarray of bool
        True where the pixel was actually used in the final fit.
    mask_out : ndarray of int
        Combined MASK column: mask (or 0) OR'd with FLATTEN_SKYLINE_BIT
        wherever the sky-line mask excluded the pixel and/or
        FLATTEN_LINE_BIT wherever the sigma-clip excluded it.
    '''
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    n = len(wave)

    good_input = (np.asarray(mask) == 0) if mask is not None else np.ones(n, dtype=bool)

    if sky_mask is not None:
        mask_wave, mask_clean = sky_mask
        sky_clean = _interp_mask_to_wave(mask_wave, mask_clean, wave)
    else:
        sky_clean = np.ones(n, dtype=bool)

    fit_mask = good_input & sky_clean
    line_excluded = np.zeros(n, dtype=bool)

    A, n_b = build_design_matrix(wave, knot_step=knot_step, solar=None)
    cont, _, _, _, n_clean = fit_continuum(flux, fit_mask, A, n_b)

    for it in range(niter):
        resid = flux - cont
        good_resid = resid[fit_mask]
        med   = np.median(good_resid)
        sigma = 1.4826 * np.median(np.abs(good_resid - med))
        if sigma <= 0:
            break
        new_bad = fit_mask & (resid - med > nsigma * sigma)
        if not new_bad.any():
            break
        line_excluded |= new_bad
        fit_mask = fit_mask & ~new_bad
        cont, _, _, _, n_clean = fit_continuum(flux, fit_mask, A, n_b)
        print(f'  sigma-clip iter {it + 1}: rejected {new_bad.sum()} pixels '
              f'(sigma={sigma:.3g})')

    resid = flux - cont

    mask_out = np.array(mask, dtype=np.int64) if mask is not None else np.zeros(n, dtype=np.int64)
    mask_out = mask_out | np.where(~sky_clean, FLATTEN_SKYLINE_BIT, 0)
    mask_out = mask_out | np.where(line_excluded, FLATTEN_LINE_BIT, 0)

    print(f'Continuum fit used {fit_mask.sum()} / {n} pixels '
          f'({good_input.sum()} pass input MASK, {sky_clean.sum()} pass sky mask, '
          f'{line_excluded.sum()} rejected by sigma-clip)')

    return cont, resid, fit_mask, mask_out


def process_file(filename, mask_file=None, knot_step=100.0, niter=3, nsigma=4.0, outfile=''):
    '''
    Read an ascii spectrum table, flatten it, and write the result.

    Parameters
    ----------
    filename : str or Path
    mask_file : str or Path or None
        Passed to PlotSpec.get_sky_mask; None uses DEFAULT_MASK_FILE.
    knot_step : float
    niter : int
        Sigma-clip refit iterations; 0 disables clipping.
    nsigma : float
        Sigma-clip rejection threshold.
    outfile : str
        Output filename; default '<stem>_flat.txt' if empty.

    Returns
    -------
    tab : astropy.table.Table
        The output table (also written to outfile).
    '''
    tab = ascii.read(filename)
    wave = tab['WAVE']
    flux = tab['FLUX']
    mask_col = tab['MASK'] if 'MASK' in tab.colnames else None

    sky_mask = get_sky_mask(mask_file)
    if sky_mask is None:
        print(f'Warning: no sky-line mask found ({mask_file or DEFAULT_MASK_FILE}); '
              'fitting continuum with sky-line regions included.')

    cont, resid, fit_mask, mask_out = flatten_spectrum(
        wave, flux, mask=mask_col, sky_mask=sky_mask, knot_step=knot_step,
        niter=niter, nsigma=nsigma)

    tab['CONT']  = cont
    tab['RESID'] = resid
    tab['MASK']  = mask_out

    if not outfile:
        outfile = f'{Path(filename).stem}_flat.txt'

    tab.write(outfile, format='ascii.fixed_width_two_line', overwrite=True)
    print(f'Wrote {outfile}')

    return tab


def steer(argv):
    filename  = ''
    mask_file = None
    knot_step = 100.0
    niter     = 3
    nsigma    = 4.0
    outfile   = ''

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-mask':
            i += 1
            mask_file = argv[i]
        elif argv[i] == '-kstep':
            i += 1
            knot_step = float(argv[i])
        elif argv[i] == '-niter':
            i += 1
            niter = int(argv[i])
        elif argv[i] == '-nsigma':
            i += 1
            nsigma = float(argv[i])
        elif argv[i] == '-out':
            i += 1
            outfile = argv[i]
        elif argv[i][0] == '-':
            print('Error: cannot parse command line:', argv)
            return
        elif filename == '':
            filename = argv[i]
        else:
            print('Error: cannot parse command line:', argv)
            return
        i += 1

    if not filename:
        print(_USAGE)
        return

    process_file(filename, mask_file=mask_file, knot_step=knot_step,
                 niter=niter, nsigma=nsigma, outfile=outfile)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
