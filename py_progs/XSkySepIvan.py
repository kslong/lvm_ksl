#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Separate one or more sky spectra into physical components (OH bands,
    atomic airglow, O2, Moon/Zodi continuum, diffuse continuum) using the
    PALACE line model and a B-spline Moon continuum.  This script is a
    wrapper for the SkyDecomp class written by Ivan Katkov (lvmsky
    repository, skysub/sky_decomp/fit.py).

Command line usage (if any):

    # Sky-file mode (auto-detected from file content):
    usage: XSkySepIvan.py sky_file.fits [row_no ...]

    # XCframe / XSFrame summary file mode:
    usage: XSkySepIvan.py xframe.fits ext [row_no ...] [-delta N]

    Arguments:

    sky_file.fits   Sky_<name>.fits from GetSky_from_CFrame_sum.py; all spectra processed if no row numbers are given.
    xframe.fits     an XCframe or XSFrame summary FITS file.
    ext             FITS extension (FLUX, SKY_EAST, SKY_WEST, ...).
    row_no          zero or more 0-based row indices; if omitted all rows are processed (subject to -delta).

    Options:

    -delta N        process every N-th row (0, N, 2N, ...) instead of all; ignored when explicit row numbers are given.
    -out outroot    set the output filename root.
    -lsf FWHM       fixed LSF FWHM in Angstroms (default 1.3 A).
    -refits N       number of iterative LSF kernel refits (default 0).

Description:

    The night sky over LVM is dominated by OH emission bands, atomic airglow
    lines (NaI, KI, [NI], OI), OI recombination multiplets, and the O2
    A-band at ~8650 A, superimposed on a smooth Moon/zodiacal continuum and
    fainter diffuse components (HO2, FeO, O2Ac).  Accurate sky subtraction
    requires decomposing these contributions so they can be removed or
    scaled independently for each exposure.

    The decomposition is performed by Ivan Katkov's SkyDecomp class
    (lvmsky/skysub/sky_decomp/fit.py), which builds a design matrix from
    PALACE (Paranal Airglow Line And Continuum Emission, Noll et al. 2024)
    reference data and solves a non-negative quadratic programme with the
    Clarabel solver.  LSF broadening is initially a fixed Gaussian; with
    -refits N the per-channel LSF kernel is iteratively refined against the
    residuals.  O2 rotational temperature is fitted from the A-band before
    the main QP solve.

    Components fitted::

      OH      402 groups of OH vibrational-rotational lines (HITRAN)
      Moon    B-spline envelope x rebinned solar spectrum (Moon + Zodi)
      Diffuse PALACE diffuse continuum: HO2, FeO, O2Ac
      Atom    atomic airglow: NaI, KI, [NI], OI (green and red)
      ORC     OI recombination multiplets (7774 and 8446 A)
      O2      molecular oxygen A-band

    Two input modes are supported:

    Sky-file mode (files produced by GetSky_from_CFrame_sum.py):
        Each row of the FLUX extension is a sky-telescope spectrum of the
        same sky field from a different exposure.  There is no IVAR
        extension; noise is estimated from pixel-to-pixel differences via
        estimate_ivar().  If no row indices are given on the command line,
        every spectrum is processed.

    XCframe mode (XCframe or XSFrame summary files):
        Standard LVM multi-extension file.  IVAR is read from the file
        when present, otherwise estimated.  Single row or strided batch
        (-many) processing.

    Output FITS structure (both modes)::

        PRIMARY    header with run parameters
        WAVE       float32 (Npix,)          wavelength array (A)
        FLUX       float32 (Nobs, Npix)     input spectra
        LINES      float32 (Nobs, Npix)     total emission (OH+ATOM+ORC+O2)
        CONT       float32 (Nobs, Npix)     total continuum (MOON+DIFFUSE)
        OH         float32 (Nobs, Npix)     OH band component
        ATOM       float32 (Nobs, Npix)     atomic airglow
        ORC        float32 (Nobs, Npix)     OI recombination lines
        O2         float32 (Nobs, Npix)     O2 A-band component
        MOON       float32 (Nobs, Npix)     Moon/zodi continuum
        DIFFUSE    float32 (Nobs, Npix)     diffuse continuum
        RESID      float32 (Nobs, Npix)     fit residuals (FLUX-LINES-CONT)
        COEF       BinTable (Nobs rows)    one column per design-matrix entry
                                            (442 columns) in physical flux units
                                            (coef * flux_scale); column names
                                            match the name column of COEF_META:
                                            OH_v{v}_N{nn}_F{f}, Moon_bs{nn},
                                            HO2, FeO, O2Ac, NaI, KI, NI_forb,
                                            OI_5577, OI_6300, OI_7774, OI_8446,
                                            O2_Aband
        COEF_META  BinTable (Ncoef rows)   one row per coefficient with columns:
                       name         design name (OH_000, ATOM_Na, Moon_bs00, ...)
                       component    OH / Moon / Diffuse / Atom / ORC / O2
                       v_upper      upper vibrational quantum number (-1 non-OH)
                       N_upper      upper rotational quantum number  (-1 non-OH)
                       F_upper      upper spin-rotation component    (-1 non-OH)
                       wave_peak    wavelength (A) of the dominant line (OH:
                                    highest Aij*gi) or peak of the convolved
                                    template (Moon B-spline peak, atomic line
                                    center, ORC feature center, O2 band peak)
                       coef_median  median physical coefficient over spectra
                       coef_mean    mean physical coefficient over spectra
                       coef_nmad    NMAD (robust scatter) over spectra
                       coef_min     minimum coefficient over spectra
                       coef_max     maximum coefficient over spectra
                       coef_skew    skewness over spectra (positive = tail of
                                    bright outliers; all in same flux units as FLUX)
        DRP_ALL    BinTable (Nobs rows)    metadata plus compact coefficient
                                            summary columns using the same names
                                            as COEF_META (all physical flux units):
                       OH_v{3..10}  summed OH amplitude over all N,F at that v
                       NaI          NaI doublet
                       KI           KI doublet
                       NI_forb      [NI] forbidden doublet
                       OI_5577      OI green forbidden line
                       OI_6300      OI red forbidden doublet
                       OI_7774      OI recombination triplet
                       OI_8446      OI recombination doublet
                       O2_Aband     O2 A-band
                       HO2          HO2 diffuse continuum
                       FeO          FeO diffuse continuum
                       O2Ac         O2Ac diffuse continuum
                       Moon_med     median Moon envelope (flux/median-solar)

Primary routines:

    decompose_one          fit one spectrum from an XCframe/XSFrame file
    decompose_one_skyfile  fit one spectrum from a Sky_<name>.fits file
    process_many           batch processing of an XCframe/XSFrame file
    process_skyfile        batch processing of a Sky_<name>.fits file

Notes:

    SkyDecomp is part of the lvmsky package written by Ivan Katkov
    (/Users/long/SDSS/lvmsky/skysub/sky_decomp/fit.py).
    PALACE data (palace/PMD/) and solar spectrum
    (Spectre_HR_LATMOS_Meftah_V1_350_1000nm.txt) are expected under
    base_dir (default /Users/long/Projects/lvm_sky2606/skysub_ivan).

    Output columns are in the same flux units as the input FITS extension
    (no FACTOR scaling is applied).

History:

    260619  ksl  XSkySepPalace.py coding begun
    260628  ksl  Renamed to XSkySepIvan.py; added Sky file mode and -lsf/-refits
    260629  ksl  Added COEF image, COEF_META table, and DRP_ALL coefficient summary
'''

import sys
import os
from pathlib import Path

from astropy.table import Table, join
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d

# SkyDecomp lives in the lvmsky repository
LVMSKY_SKYSUB = Path('/Users/long/SDSS/lvmsky/skysub')
if str(LVMSKY_SKYSUB) not in sys.path:
    sys.path.insert(0, str(LVMSKY_SKYSUB))

from sky_decomp.fit import SkyDecomp, vac_to_air, decode_hitran_id, OH_GROUP_KEYS

# Directory containing palace/PMD/ and the solar spectrum file
DEFAULT_BASE_DIR = Path('/Users/long/Projects/lvm_sky2606/skysub_ivan')

_CAP_WAVE = 5.0   # Å extra margin on each side when selecting PMD lines

_USAGE = '''Usage:
  XSkySepIvan.py sky_file.fits [row_no ...]               (Sky-file mode)
  XSkySepIvan.py xframe.fits ext [row_no ...] [-delta N]  (XCframe mode)

Arguments:
  sky_file.fits  Sky_<name>.fits from GetSky_from_CFrame_sum.py
  xframe.fits    XCframe or XSFrame summary FITS file
  ext            FITS extension (FLUX, SKY_EAST, SKY_WEST, ...)
  row_no         0-based row indices (default: all rows)

Options:
  -delta N    process every N-th row instead of all
  -lsf FWHM   LSF FWHM in Angstroms (default 1.3)
  -refits N   number of iterative LSF kernel refits (default 0)
  -out ROOT   output filename root (default: <stem>_ivan or <stem>_<ext>_ivan)
'''

# Maps PALACE internal atom/ORC names to spectroscopic FITS column names
_ATOM_COEF_COL = {
    'ATOM_K':  'KI',
    'ATOM_N':  'NI_forb',
    'ATOM_Na': 'NaI',
    'ATOM_Og': 'OI_5577',
    'ATOM_Or': 'OI_6300',
}
_ORC_COEF_COL = {
    'ATOM_Orc_OI0777': 'OI_7774',
    'ATOM_Orc_OI0845': 'OI_8446',
}


# ──────────────────────────────────────────────────────────────
# Cached SkyDecomp instance (rebuilt only if wave or LSF changes)
# ──────────────────────────────────────────────────────────────

_decomposer     = None
_decomposer_key = None   # (wave bytes, lsf_sigma bytes, base_dir str, n_spline_knots)


def _get_decomposer(wave, lsf_sigma, base_dir=DEFAULT_BASE_DIR, n_spline_knots=25):
    global _decomposer, _decomposer_key

    wave      = np.asarray(wave, float)
    lsf_sigma = np.asarray(lsf_sigma, float)
    key = (wave.tobytes(), lsf_sigma.tobytes(), str(base_dir), n_spline_knots)

    if _decomposer is None or key != _decomposer_key:
        print('Building PALACE decomposer (this may take ~10-30 s)...')
        _decomposer = SkyDecomp(
            wave,
            lsf_sigma=lsf_sigma if lsf_sigma.ndim > 0 else float(lsf_sigma),
            n_spline_knots=n_spline_knots,
            base_dir=base_dir,
        )
        _decomposer_key = key
        print('Done.')

    return _decomposer


# ──────────────────────────────────────────────────────────────
# File I/O — XCframe / XSFrame
# ──────────────────────────────────────────────────────────────

xfits = None   # cached handle for the current XCframe/XSFrame file


def get_one_spec(row=300, ext='FLUX', filename='XCframe_1.2.1_7325_48860_1_50.fits'):
    global xfits
    if xfits is None:
        xfits = fits.open(filename)

    wave = xfits['WAVE'].data.astype(np.float64)
    flux = xfits[ext].data[row].astype(np.float64)
    xtab = Table(xfits['DRP_ALL'].data)

    if ext == 'FLUX':
        ra, dec = xtab['sci_ra'][row], xtab['sci_dec'][row]
    elif ext == 'SKY_EAST':
        ra, dec = xtab['skye_ra'][row], xtab['skye_dec'][row]
    elif ext == 'SKY_WEST':
        ra, dec = xtab['skyw_ra'][row], xtab['skyw_dec'][row]
    else:
        ra, dec = 0.0, 0.0

    obstime = xtab['obstime'][row]

    ext_names = [h.name for h in xfits]
    ivar = None
    for candidate in ('IVAR', f'{ext}_IVAR'):
        if candidate in ext_names:
            ivar = xfits[candidate].data[row].astype(np.float64)
            break
    if ivar is None:
        ivar = np.ones_like(flux)

    return ra, dec, obstime, Table([wave, flux, ivar], names=['WAVE', 'FLUX', 'IVAR'])


# ──────────────────────────────────────────────────────────────
# File I/O — Sky_<name>.fits (from GetSky_from_CFrame_sum.py)
# ──────────────────────────────────────────────────────────────

_sky_fits_cache = {}   # dict: filename -> opened HDUList


def get_one_spec_skyfile(row=0, filename='Sky_WHAM_south_08.fits'):
    '''Read one spectrum from a Sky_<name>.fits file.

    Parameters
    ----------
    row : int
        Row index (0-based) into the FLUX array.
    filename : str or Path
        Path to the Sky file produced by GetSky_from_CFrame_sum.py.

    Returns
    -------
    ra, dec : float
        Sky position for this observation (from the merged DRP_ALL column).
    obstime : str
        Observation time string.
    xspec : astropy Table
        Columns WAVE, FLUX, IVAR.  IVAR is set to ones here; the calling
        function replaces it with the estimate_ivar() result.
    '''
    global _sky_fits_cache
    fn = str(filename)
    if fn not in _sky_fits_cache:
        _sky_fits_cache[fn] = fits.open(filename)
    hdul = _sky_fits_cache[fn]

    wave = hdul['WAVE'].data.astype(np.float64)
    flux = hdul['FLUX'].data[row].astype(np.float64)
    xtab = Table(hdul['DRP_ALL'].data)

    ra      = float(xtab['ra'][row])      if 'ra'      in xtab.colnames else 0.0
    dec     = float(xtab['dec'][row])     if 'dec'     in xtab.colnames else 0.0
    obstime = str(xtab['obstime'][row])   if 'obstime' in xtab.colnames else ''

    ivar = np.ones_like(flux)   # placeholder; replaced by estimate_ivar in decompose_one_skyfile
    return ra, dec, obstime, Table([wave, flux, ivar], names=['WAVE', 'FLUX', 'IVAR'])


def _is_sky_file(fits_file):
    '''Return True if fits_file is a Sky_*.fits from GetSky_from_CFrame_sum.py.'''
    try:
        with fits.open(fits_file) as hdul:
            ext_names = [h.name for h in hdul]
            if 'DRP_ALL' not in ext_names:
                return False
            tab = Table(hdul['DRP_ALL'].data)
            return 'tel' in tab.colnames
    except Exception:
        return False


# ──────────────────────────────────────────────────────────────
# Noise estimation
# ──────────────────────────────────────────────────────────────

def estimate_ivar(flux, wave, window_a=25.0, sigma_clip=5.0):
    '''Estimate wavelength-dependent IVAR from pixel-to-pixel variations.

    Finite differences of adjacent pixels cancel slow continuum gradients,
    leaving mainly noise.  Bright sky lines produce large differences and
    are excluded by sigma-clipping before the rolling RMS is computed.
    Clipped pixels are replaced by the global noise expectation so that
    line-dense regions do not produce artificially low noise estimates.

    Parameters
    ----------
    flux : array-like
        Raw spectrum (any flux units).
    wave : array-like
        Wavelength array (same length as flux).
    window_a : float
        Width of rolling window in Angstroms.
    sigma_clip : float
        Differences larger than sigma_clip * global_sigma are treated as
        line-affected and replaced by the global noise expectation.

    Returns
    -------
    ivar_est : ndarray
        Estimated IVAR in units of flux^{-2}.
    '''
    flux = np.asarray(flux, float)
    wave = np.asarray(wave, float)

    dwave     = float(np.median(np.diff(wave)))
    window_px = max(int(round(window_a / dwave)), 10)
    if window_px % 2 == 0:
        window_px += 1

    diff    = np.diff(flux)
    diff_sq = diff ** 2

    mad_diff     = float(np.nanmedian(np.abs(diff)))
    sigma_global = mad_diff * 1.4826 / np.sqrt(2.0)
    sigma_global = max(sigma_global, 1e-60)

    clip_thresh  = (sigma_clip * sigma_global) ** 2 * 2.0
    diff_sq_work = np.where(diff_sq > clip_thresh, sigma_global ** 2 * 2.0, diff_sq)

    diff_sq_pad      = np.empty_like(flux)
    diff_sq_pad[:-1] = diff_sq_work
    diff_sq_pad[-1]  = diff_sq_work[-1]

    sigma2_local = uniform_filter1d(diff_sq_pad, size=window_px) / 2.0
    sigma2_est   = np.maximum(sigma2_local, 0.1 * sigma_global ** 2)

    return 1.0 / sigma2_est


# ──────────────────────────────────────────────────────────────
# Core decomposition — XCframe / XSFrame
# ──────────────────────────────────────────────────────────────

def decompose_one(filename='XCframe_1.2.1_7325_48860_1_50.fits', row=300, ext='FLUX',
                  xplot=False, fwhm_lsf=None, n_lsf_refits=0,
                  base_dir=DEFAULT_BASE_DIR):
    '''Decompose one fibre spectrum from an XCframe/XSFrame file.

    Parameters
    ----------
    filename : str
    row : int
    ext : str
        FITS extension ('FLUX', 'SKY_EAST', 'SKY_WEST', etc.).
    xplot : bool
    fwhm_lsf : float, array, or None
        LSF FWHM in Angstroms.  None uses 1.3 A (sigma = 0.552 A).
    n_lsf_refits : int
        Number of iterative LSF kernel refits.
    base_dir : Path or str

    Returns
    -------
    xspec : astropy Table
        Columns WAVE, FLUX, IVAR, IVAR_EST plus fitted components and RESID.
    coef : ndarray, shape (n_components,)
    '''
    ra, dec, obstime, xspec = get_one_spec(row, ext, filename)

    wave = np.asarray(xspec['WAVE'], float)
    flux = np.asarray(xspec['FLUX'], float)

    if fwhm_lsf is None:
        lsf_sigma = 1.3 / 2.355
    elif np.ndim(fwhm_lsf) == 0:
        lsf_sigma = float(fwhm_lsf) / 2.355
    else:
        lsf_sigma = np.asarray(fwhm_lsf, float) / 2.355

    decomposer = _get_decomposer(wave, lsf_sigma, base_dir=base_dir)

    ivar_est          = estimate_ivar(flux, wave)
    xspec['IVAR_EST'] = ivar_est

    finite     = np.isfinite(flux)
    flux_scale = float(np.sqrt(np.nanmean(flux[finite] ** 2))) if finite.any() else 1.0
    flux_scale = max(flux_scale, 1e-30)
    flux_norm  = flux / flux_scale
    ivar_norm  = ivar_est * flux_scale ** 2

    result = decomposer.fit(flux_norm, ivar_norm, n_lsf_refits=n_lsf_refits)

    c = result.components
    xspec['LINES']   = (c['oh'] + c['atom'] + c['orc'] + c['o2']) * flux_scale
    xspec['CONT']    = (c['moon'] + c['diffuse']) * flux_scale
    xspec['OH']      = c['oh']      * flux_scale
    xspec['ATOM']    = c['atom']    * flux_scale
    xspec['ORC']     = c['orc']     * flux_scale
    xspec['O2']      = c['o2']      * flux_scale
    xspec['MOON']    = c['moon']    * flux_scale
    xspec['DIFFUSE'] = c['diffuse'] * flux_scale
    xspec['RESID']   = result.resid * flux_scale

    if xplot:
        resid_plot  = xspec['RESID']
        rms         = float(np.nanstd(resid_plot))
        resid_level = -5 * rms
        bestfit     = flux - resid_plot
        plt.figure(figsize=(12, 5))
        plt.plot(wave, flux,                         color='k',    lw=0.8, label='Obs.')
        plt.plot(wave, bestfit,                      color='red',  lw=1.2, label='Model')
        plt.plot(wave, xspec['MOON'],                color='C2',   lw=1,   ls='--', label='Moon')
        plt.plot(wave, xspec['DIFFUSE'],             color='C3',   lw=1,   ls='--', label='Diffuse')
        plt.plot(wave, resid_level + resid_plot,     color='gray', lw=0.5, label='Residuals')
        plt.axhline(resid_level, color='gray', ls='--')
        plt.legend()

    return xspec, np.asarray(result.coef, float)


# ──────────────────────────────────────────────────────────────
# Core decomposition — Sky_<name>.fits
# ──────────────────────────────────────────────────────────────

def decompose_one_skyfile(filename='Sky_WHAM_south_08.fits', row=0,
                          xplot=False, fwhm_lsf=None, n_lsf_refits=0,
                          base_dir=DEFAULT_BASE_DIR):
    '''Decompose one spectrum from a Sky_<name>.fits file.

    Parameters
    ----------
    filename : str or Path
        Sky file produced by GetSky_from_CFrame_sum.py.
    row : int
        Row index (0-based) into the FLUX array.
    xplot : bool
    fwhm_lsf : float, array, or None
        LSF FWHM in Angstroms.  None uses 1.3 A.
    n_lsf_refits : int
    base_dir : Path or str

    Returns
    -------
    xspec : astropy Table
        Same columns as decompose_one output.
    coef : ndarray
    metrics : dict
        Fit quality: chi2_red, r2, rms_resid, t_o2, t_o2_err, fit_status.
    '''
    ra, dec, obstime, xspec = get_one_spec_skyfile(row, filename)

    wave = np.asarray(xspec['WAVE'], float)
    flux = np.asarray(xspec['FLUX'], float)

    if fwhm_lsf is None:
        lsf_sigma = 1.3 / 2.355
    elif np.ndim(fwhm_lsf) == 0:
        lsf_sigma = float(fwhm_lsf) / 2.355
    else:
        lsf_sigma = np.asarray(fwhm_lsf, float) / 2.355

    decomposer = _get_decomposer(wave, lsf_sigma, base_dir=base_dir)

    ivar_est          = estimate_ivar(flux, wave)
    xspec['IVAR_EST'] = ivar_est

    finite     = np.isfinite(flux)
    flux_scale = float(np.sqrt(np.nanmean(flux[finite] ** 2))) if finite.any() else 1.0
    flux_scale = max(flux_scale, 1e-30)
    flux_norm  = flux / flux_scale
    ivar_norm  = ivar_est * flux_scale ** 2

    result = decomposer.fit(flux_norm, ivar_norm, n_lsf_refits=n_lsf_refits)

    c = result.components
    xspec['LINES']   = (c['oh'] + c['atom'] + c['orc'] + c['o2']) * flux_scale
    xspec['CONT']    = (c['moon'] + c['diffuse']) * flux_scale
    xspec['OH']      = c['oh']      * flux_scale
    xspec['ATOM']    = c['atom']    * flux_scale
    xspec['ORC']     = c['orc']     * flux_scale
    xspec['O2']      = c['o2']      * flux_scale
    xspec['MOON']    = c['moon']    * flux_scale
    xspec['DIFFUSE'] = c['diffuse'] * flux_scale
    xspec['RESID']   = result.resid * flux_scale

    metrics = dict(
        chi2_red   = result.reduced_chi2,
        r2         = result.r2,
        rms_resid  = result.rms_resid,
        t_o2       = result.t_o2,
        t_o2_err   = result.t_o2_err,
        fit_status = result.fit_status,
    )

    if xplot:
        resid_plot  = xspec['RESID']
        rms         = float(np.nanstd(resid_plot))
        resid_level = -5 * rms
        bestfit     = flux - resid_plot
        plt.figure(figsize=(12, 5))
        plt.plot(wave, flux,                         color='k',    lw=0.8, label='Obs.')
        plt.plot(wave, bestfit,                      color='red',  lw=1.2, label='Model')
        plt.plot(wave, xspec['MOON'],                color='C2',   lw=1,   ls='--', label='Moon')
        plt.plot(wave, xspec['DIFFUSE'],             color='C3',   lw=1,   ls='--', label='Diffuse')
        plt.plot(wave, resid_level + resid_plot,     color='gray', lw=0.5, label='Residuals')
        plt.axhline(resid_level, color='gray', ls='--')
        plt.legend()
        plt.title(f'{Path(filename).name}  row={row}  '
                  f'chi2_red={result.reduced_chi2:.3g}  T_O2={result.t_o2:.1f} K')

    return xspec, np.asarray(result.coef, float), metrics


# ──────────────────────────────────────────────────────────────
# Fit quality evaluation
# ──────────────────────────────────────────────────────────────

def eval_fit(xspec):
    '''Summarise fit quality from a decompose_one result.

    Uses IVAR_EST if present in xspec, otherwise estimates it on the fly
    from the FLUX column via estimate_ivar().

    Parameters
    ----------
    xspec : astropy Table
        Output from decompose_one or decompose_one_skyfile.

    Returns
    -------
    overall : dict
        r2, chi2_red, rms, frac_rms.
    comp_table : astropy Table
        One row per component with int_flux, sigma_int, snr, frac_total,
        peak_flux columns.
    '''
    wave  = np.asarray(xspec['WAVE'],  float)
    flux  = np.asarray(xspec['FLUX'],  float)
    resid = np.asarray(xspec['RESID'], float)
    dwave = float(np.median(np.diff(wave)))

    if 'IVAR_EST' in xspec.colnames:
        ivar = np.asarray(xspec['IVAR_EST'], float)
    else:
        ivar = estimate_ivar(flux, wave)

    sigma = np.where((ivar > 0) & np.isfinite(ivar), 1.0 / np.sqrt(ivar), np.nan)
    good  = np.isfinite(flux) & np.isfinite(resid) & (ivar > 0)

    chi2      = float(np.sum(resid[good] ** 2 * ivar[good]))
    chi2_red  = chi2 / max(int(np.sum(good)) - 1, 1)
    r2        = 1.0 - float(np.nanvar(resid[good])) / float(np.nanvar(flux[good]))
    rms       = float(np.nanstd(resid[good]))
    frac_rms  = rms / float(np.nanmedian(np.abs(flux[good])))

    overall = dict(r2=r2, chi2_red=chi2_red, rms=rms, frac_rms=frac_rms)

    comp_names = [n for n in ('OH', 'ATOM', 'ORC', 'O2', 'MOON', 'DIFFUSE', 'LINES', 'CONT')
                  if n in xspec.colnames]

    total_int = 0.0
    for name in ('LINES', 'CONT'):
        if name in xspec.colnames:
            total_int += float(np.nansum(xspec[name])) * dwave

    rows = []
    for name in comp_names:
        comp     = np.asarray(xspec[name], float)
        int_flux = float(np.nansum(comp)) * dwave
        peak     = float(np.nanmax(np.abs(comp)))
        mask     = (np.abs(comp) > 0.01 * peak) & np.isfinite(sigma) if peak > 0 \
                   else np.zeros(len(comp), bool)

        sigma_int = float(np.sqrt(np.nansum(sigma[mask] ** 2))) * dwave if mask.any() else np.nan
        snr       = abs(int_flux) / sigma_int if (np.isfinite(sigma_int) and sigma_int > 0) else np.nan
        frac      = int_flux / total_int       if total_int != 0 else np.nan

        rows.append((name, int_flux, sigma_int, snr, frac, peak))

    comp_table = Table(
        rows=rows,
        names=['component', 'int_flux', 'sigma_int', 'snr', 'frac_total', 'peak_flux'],
    )
    comp_table['int_flux'].format   = '.4e'
    comp_table['sigma_int'].format  = '.4e'
    comp_table['snr'].format        = '.1f'
    comp_table['frac_total'].format = '.3f'
    comp_table['peak_flux'].format  = '.4e'

    return overall, comp_table


def apply_residual_correction(xspec):
    '''Redistribute fit residuals proportionally among sky components.

    At each pixel, every component is multiplied by the same factor
    flux[i] / model[i], so their fractional contributions are preserved
    and their sum equals FLUX exactly.  Pixels where the model is near
    zero (< 1% of its median) are left uncorrected (factor = 1).

    Parameters
    ----------
    xspec : astropy Table
        Output from decompose_one or decompose_one_skyfile.

    Returns
    -------
    xspec_corr : astropy Table
        Copy of xspec with corrected component columns, RESID -> 0,
        RESID_ORIG (original residual), and CORR_FACTOR added.
    '''
    xspec_corr = xspec.copy()
    flux       = np.asarray(xspec['FLUX'], float)

    ind_names = [n for n in ('OH', 'ATOM', 'ORC', 'O2', 'MOON', 'DIFFUSE')
                 if n in xspec.colnames]

    model = np.zeros(len(flux))
    for name in ind_names:
        model += np.asarray(xspec[name], float)

    med_model = float(np.nanmedian(np.abs(model[model > 0]))) if np.any(model > 0) else 1.0
    valid      = (model > 0.01 * med_model) & np.isfinite(flux) & np.isfinite(model)
    factor     = np.where(valid, flux / model, 1.0)

    for name in ind_names:
        xspec_corr[name] = np.asarray(xspec[name], float) * factor

    line_names = [n for n in ('OH', 'ATOM', 'ORC', 'O2') if n in ind_names]
    cont_names = [n for n in ('MOON', 'DIFFUSE')          if n in ind_names]

    if line_names and 'LINES' in xspec.colnames:
        xspec_corr['LINES'] = sum(np.asarray(xspec_corr[n], float) for n in line_names)
    if cont_names and 'CONT' in xspec.colnames:
        xspec_corr['CONT']  = sum(np.asarray(xspec_corr[n], float) for n in cont_names)

    xspec_corr['RESID_ORIG'] = np.asarray(xspec['RESID'], float)
    total_corr               = (np.asarray(xspec_corr['LINES'], float)
                                + np.asarray(xspec_corr['CONT'],  float))
    xspec_corr['RESID']       = flux - total_corr
    xspec_corr['CORR_FACTOR'] = factor

    return xspec_corr


def coef_summary(xspec, coef):
    '''Summarise fit coefficients per component in physical units.

    Parameters
    ----------
    xspec : astropy Table
        Output from decompose_one or decompose_one_skyfile.
    coef : ndarray
        Coefficient vector (second return value from decompose_one*).

    Returns
    -------
    astropy Table with columns component, n_templ, coef_mean, coef_std,
    coef_min, coef_max.
    '''
    if _decomposer is None:
        raise RuntimeError('No decomposer cached; call decompose_one first.')

    flux       = np.asarray(xspec['FLUX'], float)
    finite     = np.isfinite(flux)
    flux_scale = float(np.sqrt(np.nanmean(flux[finite] ** 2))) if finite.any() else 1.0
    flux_scale = max(flux_scale, 1e-30)
    coef_phys  = np.asarray(coef, float) * flux_scale

    d         = _decomposer
    n_oh      = d.matrix_oh.shape[0]
    n_moon    = d.matrix_moon.shape[0]
    n_diffuse = d.matrix_diffuse.shape[0]
    n_atom    = d.matrix_atom.shape[0]
    n_orc     = d.matrix_orc.shape[0]
    n_o2      = d.matrix_o2.shape[0]

    i0 = 0
    slices = {}
    for name, n in [('OH', n_oh), ('MOON', n_moon), ('DIFFUSE', n_diffuse),
                    ('ATOM', n_atom), ('ORC', n_orc), ('O2', n_o2)]:
        slices[name] = slice(i0, i0 + n)
        i0 += n

    moon_sl  = slices['MOON']
    moon_c   = coef_phys[moon_sl]
    bspl_mat = _decomposer.matrix_moon
    solar_rb = _decomposer.vector_moon
    fitted_moon          = moon_c @ bspl_mat
    safe                 = solar_rb > 0.01 * float(np.nanmedian(solar_rb))
    envelope             = np.where(safe, fitted_moon / solar_rb, np.nan)
    moon_envelope_median = float(np.nanmedian(envelope))

    rows = []
    for name, sl in slices.items():
        c = coef_phys[sl]
        if name == 'MOON':
            rows.append((name, int(c.size), moon_envelope_median,
                         float(np.std(c)), float(np.min(c)), float(np.max(c))))
        else:
            rows.append((name, int(c.size), float(np.mean(c)),
                         float(np.std(c)), float(np.min(c)), float(np.max(c))))

    t = Table(rows=rows,
              names=['component', 'n_templ', 'coef_mean', 'coef_std', 'coef_min', 'coef_max'])
    for col in ('coef_mean', 'coef_std', 'coef_min', 'coef_max'):
        t[col].format = '.4e'
    return t


# ──────────────────────────────────────────────────────────────
# Coefficient metadata and DRP summary
# ──────────────────────────────────────────────────────────────

def build_coef_meta(base_dir=DEFAULT_BASE_DIR):
    '''Build a Table mapping every design-matrix column to its physical identity.

    Must be called after at least one decomposition so _decomposer is cached.
    The returned Table has one row per coefficient in the same order as the
    coef vector returned by decompose_one / decompose_one_skyfile.

    Six distribution statistics (coef_median, coef_mean, coef_nmad, coef_min,
    coef_max, coef_skew) are added by process_skyfile / process_many after the
    processing loop, before writing COEF_META to disk.

    Returns
    -------
    astropy.table.Table
        One row per design-matrix entry with columns: ``name`` (descriptive
        string matching the COEF BinTable column names, e.g.
        ``OH_v4_N12_F1``, ``Moon_bs07``, ``NaI``, ``O2_Aband``);
        ``component`` (OH / Moon / Diffuse / Atom / ORC / O2);
        ``v_upper``, ``N_upper``, ``F_upper`` (int16 OH quantum numbers,
        -1 for non-OH); ``wave_peak`` (float32, dominant line wavelength in
        Angstroms for OH groups, or template peak wavelength for all other
        components).
    '''
    d = _decomposer
    if d is None:
        raise RuntimeError('build_coef_meta: run a decomposition first.')

    wave = d.wave
    oh_file = Path(base_dir) / 'palace' / 'PMD' / 'pmd_popmodel_OH.dat'
    oh = Table.read(str(oh_file), format='ascii.basic', guess=False,
                    comment='#', fast_reader=False)
    oh['wave'] = vac_to_air(np.asarray(oh['lam'], float) * 1e4)
    mask = ((oh['wave'] >= wave.min() - _CAP_WAVE) &
            (oh['wave'] <= wave.max() + _CAP_WAVE))
    oh = decode_hitran_id(oh[mask])
    groups = oh.group_by(OH_GROUP_KEYS).groups

    rows = []
    for i, grp in enumerate(groups):
        v   = int(grp['v_upper'][0])
        n   = int(grp['N_upper'][0])
        f   = int(grp['F_upper'][0])
        amp = np.asarray(grp['Aij'] * grp['gi'], float)
        wpeak = float(np.asarray(grp['wave'])[np.argmax(amp)])
        rows.append((f'OH_v{v}_N{n:02d}_F{f}', 'OH', v, n, f, wpeak))

    # For non-OH components use the peak of the convolved template on wave
    solar  = d.vector_moon
    safe   = solar > 0.01 * float(np.nanmedian(solar))

    for i, name in enumerate(d.moon_names):
        envelope = np.where(safe, d.matrix_moon[i] / solar, 0.0)
        wpeak = float(wave[np.argmax(envelope)])
        rows.append((name, 'Moon', -1, -1, -1, wpeak))

    vecs = [d.vector_ho2, d.vector_feo, d.vector_o2ac]
    for name, vec in zip(d.diffuse_names, vecs):
        wpeak = float(wave[np.argmax(vec)])
        rows.append((name, 'Diffuse', -1, -1, -1, wpeak))

    for i, name in enumerate(d.atom_names):
        wpeak = float(wave[np.argmax(d.matrix_atom[i])])
        rows.append((_ATOM_COEF_COL[name], 'Atom', -1, -1, -1, wpeak))

    for i, name in enumerate(d.orc_names):
        wpeak = float(wave[np.argmax(d.matrix_orc[i])])
        rows.append((_ORC_COEF_COL[name], 'ORC', -1, -1, -1, wpeak))

    for i, name in enumerate(d.o2_names):
        wpeak = float(wave[np.argmax(d.matrix_o2[i])])
        rows.append(('O2_Aband', 'O2', -1, -1, -1, wpeak))

    meta = Table(rows=rows,
                 names=['name', 'component', 'v_upper', 'N_upper', 'F_upper',
                        'wave_peak'])
    for col in ('v_upper', 'N_upper', 'F_upper'):
        meta[col] = meta[col].astype(np.int16)
    meta['wave_peak'] = meta['wave_peak'].astype(np.float32)
    return meta


def _coef_drp_summary(coef_phys, meta):
    '''Compute compact per-spectrum DRP summary from physical-unit coefficients.

    Parameters
    ----------
    coef_phys : ndarray
        Coefficient vector scaled to physical flux units (coef * flux_scale).
    meta : astropy Table
        Output of build_coef_meta().

    Returns
    -------
    dict mapping column name -> float, with 20 entries using the same naming
    convention as COEF_META:
        OH_v{3..10}   summed OH amplitude over all N, F groups at that v_upper
        NaI           NaI doublet amplitude
        KI            KI doublet amplitude
        NI_forb       [NI] forbidden doublet amplitude
        OI_5577       OI green forbidden line amplitude
        OI_6300       OI red forbidden doublet amplitude
        OI_7774       OI recombination triplet amplitude
        OI_8446       OI recombination doublet amplitude
        O2_Aband      O2 A-band amplitude
        HO2           HO2 diffuse continuum amplitude
        FeO           FeO diffuse continuum amplitude
        O2Ac          O2Ac diffuse continuum amplitude
        Moon_med      median of fitted Moon envelope (flux / median-solar)
    '''
    d = _decomposer
    coef_phys = np.asarray(coef_phys, float)
    comp = np.asarray(meta['component'])
    v_col = np.asarray(meta['v_upper'])

    summary = {}

    # Per-vibrational-level OH sums
    oh_sel = comp == 'OH'
    for v in range(3, 11):
        summary[f'OH_v{v}'] = float(np.sum(coef_phys[oh_sel & (v_col == v)]))

    # Atomic species — use same names as COEF_META via _ATOM_COEF_COL
    atom_idx = np.where(comp == 'Atom')[0]
    for j, name in enumerate(d.atom_names):
        summary[_ATOM_COEF_COL[name]] = float(coef_phys[atom_idx[j]])

    # ORC groups — use same names as COEF_META via _ORC_COEF_COL
    orc_idx = np.where(comp == 'ORC')[0]
    for j, name in enumerate(d.orc_names):
        summary[_ORC_COEF_COL[name]] = float(coef_phys[orc_idx[j]])

    # O2
    o2_idx = np.where(comp == 'O2')[0]
    summary['O2_Aband'] = float(coef_phys[o2_idx[0]]) if o2_idx.size else np.nan

    # Diffuse continua — use PALACE names directly (HO2, FeO, O2Ac)
    diff_idx = np.where(comp == 'Diffuse')[0]
    for j, name in enumerate(d.diffuse_names):
        summary[name] = float(coef_phys[diff_idx[j]])

    # Moon: median of fitted envelope relative to normalised solar
    moon_idx    = np.where(comp == 'Moon')[0]
    moon_c      = coef_phys[moon_idx]
    fitted_moon = moon_c @ d.matrix_moon
    solar_rb    = d.vector_moon
    safe        = solar_rb > 0.01 * float(np.nanmedian(solar_rb))
    envelope    = np.where(safe, fitted_moon / solar_rb, np.nan)
    summary['Moon_med'] = float(np.nanmedian(envelope))

    return summary


# ──────────────────────────────────────────────────────────────
# Batch processing helpers
# ──────────────────────────────────────────────────────────────

def array_to_table(arr):
    n_rows, n_cols = arr.shape
    col_names = [f'C{i+1:02d}' for i in range(n_cols)]
    table = Table(arr, names=col_names)
    for col in col_names:
        table[col].format = '.4f'
    return table


def table_col_to_array(spectra_list, colname='FLUX'):
    return np.vstack([s[colname] for s in spectra_list])


# ──────────────────────────────────────────────────────────────
# Batch processing — XCframe / XSFrame
# ──────────────────────────────────────────────────────────────

def process_many(filename, ext, rows=None, delta=None, outroot='',
                 fwhm_lsf=None, n_lsf_refits=0, base_dir=DEFAULT_BASE_DIR):
    '''Process rows of an XCframe/XSFrame file and write a FITS output.

    Parameters
    ----------
    filename : str
        Path to the XCframe/XSFrame FITS file.
    ext : str
        FITS extension name for the spectra (e.g. 'SKY_EAST').
    rows : list of int or None
        Explicit row indices to process.  If None, all rows are processed
        (or every delta-th row if delta is set).
    delta : int or None
        Row stride: process rows 0, delta, 2*delta, ...  Ignored if rows
        is not None.
    outroot : str
        Output filename root.  Default: <stem>_<ext>_ivan.
    fwhm_lsf : float or None
    n_lsf_refits : int
    base_dir : Path or str
    '''
    global xfits
    if xfits is None:
        xfits = fits.open(filename)

    n_total = xfits[ext].data.shape[0]

    if rows is not None:
        rows_to_process = list(rows)
    elif delta is not None:
        rows_to_process = list(range(0, n_total, int(delta)))
    else:
        rows_to_process = list(range(n_total))

    xtab     = Table(xfits['DRP_ALL'].data)
    spectra, coefficients = [], []

    for i, row in enumerate(rows_to_process):
        xspec, coef = decompose_one(filename, row, ext,
                                     fwhm_lsf=fwhm_lsf,
                                     n_lsf_refits=n_lsf_refits,
                                     base_dir=base_dir)
        spectra.append(xspec)
        coefficients.append(coef)
        if (i + 1) % 10 == 0 or (i + 1) == len(rows_to_process):
            print(f'Finished {i + 1}/{len(rows_to_process)} (row {row})')

    print(f'Number of spectra processed: {len(spectra)}')

    # Physical-unit coefficients: coef (normalised) * flux_scale per spectrum
    coeff = np.array(coefficients)
    flux_scales = []
    for xspec in spectra:
        flux = np.asarray(xspec['FLUX'], float)
        finite = np.isfinite(flux)
        fs = float(np.sqrt(np.nanmean(flux[finite] ** 2))) if finite.any() else 1.0
        flux_scales.append(max(fs, 1e-30))
    coeff_phys = coeff * np.array(flux_scales)[:, None]

    meta = build_coef_meta(base_dir=base_dir)
    meta['coef_median'] = np.median(coeff_phys, axis=0).astype(np.float32)
    meta['coef_mean']   = np.mean(coeff_phys,   axis=0).astype(np.float32)
    _dev = coeff_phys - np.median(coeff_phys, axis=0)
    meta['coef_nmad']   = (1.4826 * np.median(np.abs(_dev), axis=0)).astype(np.float32)
    meta['coef_min']    = np.min(coeff_phys,    axis=0).astype(np.float32)
    meta['coef_max']    = np.max(coeff_phys,    axis=0).astype(np.float32)
    _std = np.std(coeff_phys, axis=0)
    _std_safe = np.where(_std > 0, _std, np.nan)
    meta['coef_skew']   = (np.mean((_dev / _std_safe) ** 3, axis=0)).astype(np.float32)

    qtab = xtab[rows_to_process]
    qtab['row'] = rows_to_process
    summary_rows = [_coef_drp_summary(coeff_phys[i], meta) for i in range(len(spectra))]
    for col in summary_rows[0]:
        qtab[col] = np.array([r[col] for r in summary_rows], dtype=np.float32)
    xtab_final = qtab

    if outroot == '':
        outroot = f'{Path(filename).stem}_{ext}_ivan'

    spec_cols = ['FLUX', 'LINES', 'CONT', 'OH', 'ATOM', 'ORC', 'O2',
                 'MOON', 'DIFFUSE', 'RESID']
    hdr = fits.Header()
    hdr['INPUT']   = str(filename)
    hdr['EXT']     = ext
    hdr['N_PROC']  = len(rows_to_process)
    hdr['FWHMLSF'] = fwhm_lsf if fwhm_lsf is not None else 1.3
    hdr['NREFITS'] = n_lsf_refits
    if delta is not None:
        hdr['DELTA'] = delta
    hdu_list = [
        fits.PrimaryHDU(header=hdr),
        fits.ImageHDU(data=xfits['WAVE'].data, name='WAVE'),
    ]
    for col in spec_cols:
        hdu_list.append(fits.ImageHDU(
            data=table_col_to_array(spectra, col).astype(np.float32),
            name=col,
        ))
    coef_cols = [fits.Column(name=str(meta['name'][i]), format='E',
                             array=coeff_phys[:, i].astype(np.float32))
                 for i in range(len(meta))]
    coef_hdu = fits.BinTableHDU.from_columns(coef_cols)
    coef_hdu.name = 'COEF'
    hdu_list.append(coef_hdu)
    hdu_list.append(fits.BinTableHDU(meta, name='COEF_META'))
    hdu_list.append(fits.BinTableHDU(xtab_final, name='DRP_ALL'))

    outfile = f'{outroot}.fits'
    fits.HDUList(hdu_list).writeto(outfile, overwrite=True)
    print(f'Wrote results to {outfile}')


# ──────────────────────────────────────────────────────────────
# Batch processing — Sky_<name>.fits
# ──────────────────────────────────────────────────────────────

def process_skyfile(filename, rows=None, delta=None, outroot='',
                    fwhm_lsf=None, n_lsf_refits=0, base_dir=DEFAULT_BASE_DIR):
    '''Process some or all spectra in a Sky_<name>.fits file.

    Parameters
    ----------
    filename : str or Path
        Sky file produced by GetSky_from_CFrame_sum.py.
    rows : list of int or None
        Explicit row indices (0-based) to process.  If None, all rows are
        processed (or every delta-th row if delta is set).
    delta : int or None
        Row stride: process rows 0, delta, 2*delta, ...  Ignored if rows
        is not None.
    outroot : str
        Output filename root.  Defaults to <stem>_ivan.
    fwhm_lsf : float or None
    n_lsf_refits : int
    base_dir : Path or str
    '''
    global _sky_fits_cache
    fn = str(filename)
    if fn not in _sky_fits_cache:
        _sky_fits_cache[fn] = fits.open(filename)
    hdul = _sky_fits_cache[fn]

    n_total = hdul['FLUX'].data.shape[0]
    if rows is not None:
        pass  # use as-is
    elif delta is not None:
        rows = list(range(0, n_total, int(delta)))
    else:
        rows = list(range(n_total))

    xtab = Table(hdul['DRP_ALL'].data)
    spectra, coefficients, metrics_list, processed_rows = [], [], [], []

    for i, row in enumerate(rows):
        print(f'Row {row} ({i + 1}/{len(rows)})')
        xspec, coef, metrics = decompose_one_skyfile(
            filename, row,
            fwhm_lsf=fwhm_lsf,
            n_lsf_refits=n_lsf_refits,
            base_dir=base_dir,
        )
        spectra.append(xspec)
        coefficients.append(coef)
        metrics_list.append(metrics)
        processed_rows.append(row)

    print(f'Number of spectra processed: {len(spectra)}')

    # Physical-unit coefficients: coef (normalised) * flux_scale per spectrum
    coeff = np.array(coefficients)
    flux_scales = []
    for xspec in spectra:
        flux = np.asarray(xspec['FLUX'], float)
        finite = np.isfinite(flux)
        fs = float(np.sqrt(np.nanmean(flux[finite] ** 2))) if finite.any() else 1.0
        flux_scales.append(max(fs, 1e-30))
    coeff_phys = coeff * np.array(flux_scales)[:, None]

    meta = build_coef_meta(base_dir=base_dir)
    meta['coef_median'] = np.median(coeff_phys, axis=0).astype(np.float32)
    meta['coef_mean']   = np.mean(coeff_phys,   axis=0).astype(np.float32)
    _dev = coeff_phys - np.median(coeff_phys, axis=0)
    meta['coef_nmad']   = (1.4826 * np.median(np.abs(_dev), axis=0)).astype(np.float32)
    meta['coef_min']    = np.min(coeff_phys,    axis=0).astype(np.float32)
    meta['coef_max']    = np.max(coeff_phys,    axis=0).astype(np.float32)
    _std = np.std(coeff_phys, axis=0)
    _std_safe = np.where(_std > 0, _std, np.nan)
    meta['coef_skew']   = (np.mean((_dev / _std_safe) ** 3, axis=0)).astype(np.float32)

    qtab = xtab[processed_rows]
    qtab['row'] = processed_rows

    # Fit quality metrics
    qtab['chi2_red']   = [m['chi2_red']   for m in metrics_list]
    qtab['r2']         = [m['r2']         for m in metrics_list]
    qtab['rms_resid']  = [m['rms_resid']  for m in metrics_list]
    qtab['t_o2']       = [m['t_o2']       for m in metrics_list]
    qtab['t_o2_err']   = [m['t_o2_err']   for m in metrics_list]
    qtab['fit_status'] = [m['fit_status'] for m in metrics_list]

    # Compact per-component coefficient summary
    summary_rows = [_coef_drp_summary(coeff_phys[i], meta) for i in range(len(spectra))]
    for col in summary_rows[0]:
        qtab[col] = np.array([r[col] for r in summary_rows], dtype=np.float32)
    xtab_final = qtab

    if outroot == '':
        outroot = f'{Path(filename).stem}_ivan'

    spec_cols = ['FLUX', 'LINES', 'CONT', 'OH', 'ATOM', 'ORC', 'O2',
                 'MOON', 'DIFFUSE', 'RESID']
    hdr = fits.Header()
    hdr['INPUT']    = str(filename)
    hdr['N_PROC']   = len(processed_rows)
    hdr['FWHMLSF']  = fwhm_lsf if fwhm_lsf is not None else 1.3
    hdr['NREFITS']  = n_lsf_refits
    if delta is not None:
        hdr['DELTA'] = delta
    hdu_list = [
        fits.PrimaryHDU(header=hdr),
        fits.ImageHDU(data=hdul['WAVE'].data, name='WAVE'),
    ]
    for col in spec_cols:
        hdu_list.append(fits.ImageHDU(
            data=table_col_to_array(spectra, col).astype(np.float32),
            name=col,
        ))
    coef_cols = [fits.Column(name=str(meta['name'][i]), format='E',
                             array=coeff_phys[:, i].astype(np.float32))
                 for i in range(len(meta))]
    coef_hdu = fits.BinTableHDU.from_columns(coef_cols)
    coef_hdu.name = 'COEF'
    hdu_list.append(coef_hdu)
    hdu_list.append(fits.BinTableHDU(meta, name='COEF_META'))
    hdu_list.append(fits.BinTableHDU(xtab_final, name='DRP_ALL'))

    outfile = '%s.fits' % outroot
    fits.HDUList(hdu_list).writeto(outfile, overwrite=True)
    print('Wrote results to %s' % outfile)


# ──────────────────────────────────────────────────────────────
# Command-line interface
# ──────────────────────────────────────────────────────────────

def steer(argv):
    filename     = ''
    ext          = ''
    delta        = None
    outroot      = ''
    fwhm_lsf     = None
    n_lsf_refits = 0
    positional_ints = []   # explicit row numbers after filename [ext]

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-out':
            i += 1
            outroot = argv[i]
        elif argv[i] == '-lsf':
            i += 1
            fwhm_lsf = float(argv[i])
        elif argv[i] == '-refits':
            i += 1
            n_lsf_refits = int(argv[i])
        elif argv[i] == '-delta':
            i += 1
            delta = int(argv[i])
        elif argv[i][0] == '-':
            print('Error: cannot parse command line:', argv)
            return
        elif filename == '':
            filename = argv[i]
        elif ext == '' and not argv[i].lstrip('-').isdigit():
            ext = argv[i]
        elif argv[i].lstrip('-').isdigit():
            positional_ints.append(int(argv[i]))
        else:
            print('Error: cannot parse command line:', argv)
            return
        i += 1

    if not filename:
        print(_USAGE)
        return

    # ── Sky-file mode ─────────────────────────────────────────
    if ext == '':
        if _is_sky_file(filename):
            rows = positional_ints if positional_ints else None
            process_skyfile(filename, rows=rows, delta=delta, outroot=outroot,
                            fwhm_lsf=fwhm_lsf, n_lsf_refits=n_lsf_refits)
        else:
            print('Error: specify a FITS extension name (e.g. SKY_EAST) for XCframe files.')
            print(_USAGE)
        return

    # ── XCframe mode ──────────────────────────────────────────
    # Explicit row numbers take precedence over -delta.
    rows = positional_ints if positional_ints else None
    process_many(filename, ext, rows=rows, delta=delta, outroot=outroot,
                 fwhm_lsf=fwhm_lsf, n_lsf_refits=n_lsf_refits)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
