#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Fit both a Gaussian and a Moffat profile (each with a constant
    background) to a fixed set of airglow sky lines, on raw
    (pre-subtraction) sky spectra, and compare which profile shape
    better describes the true instrumental line shape -- independent
    of any sky-subtraction algorithm or the PALACE decomposer.

    Fully standalone by design: imports nothing from any other py_progs
    script (only external packages -- numpy, astropy, matplotlib, scipy,
    lmfit).  SKY_LINES and the IVAR estimator are self-contained copies
    (originally from sky_gaussfit.py and XSkySepIvan.py respectively;
    see History) rather than imports, so this script has no dependency
    on lvm_gaussfit.py or anything else that may be reused/changed for
    other, unrelated analyses later.

Command line usage (if any)::

    lvm_line_profile.py [-delta N] [-lines name1,name2,...]
                        [-ext SKY_EAST|SKY_WEST] [-out ROOT] filename

    Arguments::

        filename     XCframe FITS file (SKY_EAST/SKY_WEST extensions) or
                     a Sky_<name>.fits file (FLUX extension) containing
                     raw, pre-subtraction sky spectra

    Options::

        -delta N     process every N-th row (default: 50)
        -lines LIST  comma-separated subset of SKY_LINES names to fit
                     (default: all 18)
        -ext NAME    which raw sky column to use for an XCframe file:
                     SKY_EAST (default) or SKY_WEST; ignored for
                     Sky_<name>.fits files, which use FLUX directly
        -out ROOT    output filename root (default: <stem>_profile)

Description:

    For each requested line:

    1. A median spectrum is formed across all selected rows and fit
       with both profiles -- the primary, highest-S/N comparison,
       plotted as data + both fits + residuals (one PNG per line).
    2. Each individual selected row is also fit with both profiles,
       building a per-row table of fit parameters, AIC/BIC, and
       core/wing residual RMS -- this tests whether any preference for
       one profile is consistent across many independent exposures,
       not just a property of the deep stack.

    Both profiles are parameterized so their fitted 'flux' equals the
    analytic integral of the profile (matching lvm_gaussfit.py's
    convention), making the two directly comparable.  The Moffat profile
    reduces to a Gaussian as beta -> infinity, so a finite,
    well-constrained beta across many lines/rows would indicate genuine
    non-Gaussian wings, not just noise; beta drifting to the fit's upper
    bound would indicate no real preference for Moffat.

    The line center is a free parameter in both fits (bounded only to
    the fit window, not fixed at the nominal catalog wavelength), so a
    real wavelength shift is not assumed away -- the fitted center and
    its shift from the nominal wavelength are recorded for both profiles.

    If the input file has an LSF extension (per-row, per-wavelength FWHM
    in Angstroms -- same one used by SkySubDev2.py), the fitted Gaussian
    FWHM is also compared directly against that file's own stated LSF
    FWHM at each line's wavelength (median-across-rows for the median
    fit; per-row for the per-row fits) -- an independent check of
    whether the LSF extension's values are themselves accurate, entirely
    outside any subtraction algorithm.

    Comparison metrics -- Moffat has one more free parameter (beta) than
    Gaussian, so raw chi2/RMSE always favours it; AIC/BIC penalize that
    extra freedom::

        delta_aic = aic_gaussian - aic_moffat   (positive: Moffat preferred)
        delta_bic = bic_gaussian - bic_moffat
        rms_core_<profile>, rms_wing_<profile>  -- a FIXED, profile-
            independent split of the fit window (core = center +/- 1.5 A,
            wing = the rest), so the two profiles are judged against the
            same pixels regardless of their own fitted width.
        fwhm_diff_g = g_fwhm - lsf_fwhm   (only if an LSF extension is
            present; positive means the independent Gaussian fit found a
            wider line than the file's own stated LSF)
        g_shift, m_shift = fitted center - nominal catalog wavelength

Output:

    <ROOT>_median_<line>.png   one file per line: median spectrum, both
                               fits (top, with fitted FWHM/shift and the
                               LSF extension's FWHM if available), and
                               both fits' residuals (bottom)
    <ROOT>_perrow.txt          ascii table, one row per (input row, line):
                               both fits' parameters (including fitted
                               center/shift), AIC/BIC, core/wing RMS,
                               delta_aic, delta_bic, and (if available)
                               lsf_fwhm/fwhm_diff_g
    <ROOT>_summary.png         per-line distributions of delta_aic and
                               beta across all fitted rows; also fit
                               FWHM - LSF ext FWHM and center shift, if
                               an LSF extension is available

    A summary table (per line: median delta_aic, fraction of rows
    favouring Moffat by delta_aic > 2, median beta and its IQR, and --
    if available -- median FWHM-vs-LSF-ext difference and median center
    shift) is also printed to the screen.

Notes:

    Requires lmfit.  estimate_ivar() (a self-contained copy of
    XSkySepIvan.estimate_ivar()) gives a wavelength-dependent noise
    estimate from pixel-to-pixel scatter, already robust against bright
    sky lines, used as the fit weights for both profiles.

    A row/line combination is skipped (not an error) if there are too
    few finite pixels in its fit window, or if either fit raises an
    exception (e.g. non-convergence); it is simply absent from
    <ROOT>_perrow.txt.

History::

    260709 ksl Coding begun.
    260709 ksl Added comparison against the input file's own LSF
               extension (per-row/median fitted Gaussian FWHM vs that
               row's/median's stated LSF FWHM at each line's
               wavelength), and recorded the fitted line center (and
               its shift from the nominal catalog wavelength) for both
               profiles -- the line center was already a free parameter
               in both fits, just not previously reported anywhere.
    260709 ksl Made fully standalone: SKY_LINES and estimate_ivar() are
               now self-contained copies (from sky_gaussfit.py and
               XSkySepIvan.py respectively) instead of imports, so this
               script has no dependency on any other py_progs module --
               kept separate from lvm_gaussfit.py specifically so that
               script remains free for other, unrelated analyses.  Also
               carries sky_gaussfit.py's sky6553 wavelength correction
               (6553.0 -> 6553.617 A, window widened to 6549-6558 A;
               confirmed via this script's own fit, which converges to
               6553.614 A regardless of which nominal value is supplied).

'''

import sys
import os

import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import gammaln
from scipy.ndimage import uniform_filter1d
from lmfit import Model

# Airglow lines: (name, center_wavelength, window_min, window_max) in
# Angstroms.  Fitted at fixed, unshifted wavelengths (ESO UVES sky
# spectrum atlas).  Self-contained copy of sky_gaussfit.py's SKY_LINES
# (not imported, so this script has no dependency on that module --
# see Notes/History) -- includes the sky6553 wavelength correction
# (6553.0 -> 6553.617 A) confirmed by this script itself.
SKY_LINES = [
    ('sky5577', 5577.34668,  5572.,  5582.),
    ('sky6300', 6300.308594, 6295.,  6305.),
    ('sky6363', 6363.782715, 6358.,  6368.),
    ('sky6533', 6533.04,     6528.,  6538.),
    ('sky6553', 6553.617,    6549.,  6558.),
    ('sky6577', 6577.2,      6572.,  6582.),
    ('sky6912', 6912.623,    6907.,  6917.),
    ('sky6923', 6923.220,    6918.,  6928.),
    ('sky6939', 6939.521,    6934.,  6944.),
    ('sky7358', 7358.680176, 7353.,  7363.),
    ('sky7392', 7392.209961, 7387.,  7397.),
    ('sky7914', 7913.717773, 7908.,  7918.),
    ('sky8344', 8344.613281, 8339.,  8349.),
    ('sky8399', 8399.175781, 8394.,  8404.),
    ('sky8827', 8827.112305, 8822.,  8832.),
    ('sky8988', 8988.383789, 8983.,  8993.),
    ('sky9552', 9552.546875, 9547.,  9557.),
    ('sky9719', 9719.838867, 9714.,  9724.),
]


def estimate_ivar(flux, wave, window_a=25.0, sigma_clip=5.0):
    '''
    Estimate wavelength-dependent IVAR from pixel-to-pixel variations.
    Self-contained copy of XSkySepIvan.estimate_ivar() (not imported, so
    this script has no dependency on that module -- see Notes/History).

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


_USAGE = '''Usage:
  lvm_line_profile.py [-delta N] [-lines name1,name2,...]
                      [-ext SKY_EAST|SKY_WEST] [-out ROOT] filename

Arguments:
  filename     XCframe FITS file or Sky_<name>.fits file with raw,
               pre-subtraction sky spectra

Options:
  -delta N     process every N-th row (default: 50)
  -lines LIST  comma-separated subset of line names (default: all)
  -ext NAME    SKY_EAST (default) | SKY_WEST -- XCframe input only
  -out ROOT    output filename root (default: <stem>_profile)
'''


# ──────────────────────────────────────────────────────────────
# Profile functions -- flux-parameterized: the fitted 'flux' equals the
# analytic integral of the profile, for direct Gaussian/Moffat comparison
# ──────────────────────────────────────────────────────────────

def gaussian_profile(x, flux, center, fwhm, background):
    '''Gaussian + constant background, parameterized by integrated flux.'''
    sigma = fwhm / 2.355
    amp = flux / (sigma * np.sqrt(2 * np.pi))
    return amp * np.exp(-0.5 * ((x - center) / sigma) ** 2) + background


def moffat_profile(x, flux, center, alpha, beta, background):
    '''
    Moffat + constant background: I(x) = A*(1 + ((x-center)/alpha)^2)^(-beta)
    + background.  The analytic integral of the un-normalized Moffat over
    -inf..inf is A * alpha * sqrt(pi) * Gamma(beta-0.5) / Gamma(beta)
    (valid for beta > 0.5), so amp is chosen to make the fitted 'flux'
    equal that integral, matching gaussian_profile's convention.
    '''
    log_ratio = gammaln(beta - 0.5) - gammaln(beta)
    norm = alpha * np.sqrt(np.pi) * np.exp(log_ratio)
    amp = flux / norm
    return amp * (1.0 + ((x - center) / alpha) ** 2) ** (-beta) + background


def moffat_fwhm(alpha, beta):
    '''FWHM of a Moffat profile with the given alpha (width) and beta (shape).'''
    return 2.0 * alpha * np.sqrt(2.0 ** (1.0 / beta) - 1.0)


# ──────────────────────────────────────────────────────────────
# Fitting
# ──────────────────────────────────────────────────────────────

def fit_gaussian(wave, flux, ivar, center, wmin, wmax, init_fwhm=1.6):
    '''Fit gaussian_profile to one spectrum in [wmin, wmax]; None if too few good pixels.'''
    mask = ((wave >= wmin) & (wave <= wmax) & np.isfinite(flux)
           & np.isfinite(ivar) & (ivar > 0))
    x = wave[mask]; y = flux[mask]; w = np.sqrt(ivar[mask])
    if len(x) < 5:
        return None
    dx = float(np.median(np.diff(x)))
    background0 = float(np.median(y))
    init_flux = float(np.sum((y - background0) * dx))

    gmodel = Model(gaussian_profile)
    params = gmodel.make_params(flux=init_flux, center=center, fwhm=init_fwhm,
                                background=background0)
    params['fwhm'].min = 0.3
    params['fwhm'].max = 8.0
    params['center'].min = wmin
    params['center'].max = wmax
    try:
        result = gmodel.fit(y, params, x=x, weights=w)
    except Exception:
        return None
    return dict(x=x, y=y, result=result)


def fit_moffat(wave, flux, ivar, center, wmin, wmax, init_fwhm=1.6, init_beta=3.5):
    '''Fit moffat_profile to one spectrum in [wmin, wmax]; None if too few good pixels.'''
    mask = ((wave >= wmin) & (wave <= wmax) & np.isfinite(flux)
           & np.isfinite(ivar) & (ivar > 0))
    x = wave[mask]; y = flux[mask]; w = np.sqrt(ivar[mask])
    if len(x) < 6:
        return None
    dx = float(np.median(np.diff(x)))
    background0 = float(np.median(y))
    init_flux = float(np.sum((y - background0) * dx))
    init_alpha = (init_fwhm / 2.0) / np.sqrt(2.0 ** (1.0 / init_beta) - 1.0)

    mmodel = Model(moffat_profile)
    params = mmodel.make_params(flux=init_flux, center=center, alpha=init_alpha,
                                beta=init_beta, background=background0)
    params['alpha'].min = 0.1
    params['alpha'].max = 8.0
    params['beta'].min = 0.6
    params['beta'].max = 50.0
    params['center'].min = wmin
    params['center'].max = wmax
    try:
        result = mmodel.fit(y, params, x=x, weights=w)
    except Exception:
        return None
    return dict(x=x, y=y, result=result)


def core_wing_rms(x, y, model_y, center, half_core=1.5):
    '''
    Split residuals (y - model_y) into a fixed, profile-independent core
    (``abs(x-center) <= half_core``) and wing (the rest of the fit
    window), and return the RMS of each.
    '''
    resid = y - model_y
    core = np.abs(x - center) <= half_core
    wing = ~core
    rms_core = float(np.sqrt(np.mean(resid[core] ** 2))) if core.any() else np.nan
    rms_wing = float(np.sqrt(np.mean(resid[wing] ** 2))) if wing.any() else np.nan
    return rms_core, rms_wing


# ──────────────────────────────────────────────────────────────
# Data loading
# ──────────────────────────────────────────────────────────────

def load_raw_sky(filename, ext='SKY_EAST'):
    '''
    Load a raw (pre-subtraction) sky flux array from either an XCframe
    file (SKY_EAST/SKY_WEST extensions) or a Sky_<name>.fits file (FLUX
    extension).  Returns (wave, flux2d, lsf2d), flux2d/lsf2d shape
    (n_rows, n_wave).  lsf2d is None if the file has no LSF extension
    (per-row, per-wavelength FWHM in Angstroms, same convention as
    SkySubDev2.py -- written by SummarizeCframe.py's make_med_spec from
    the source CFrame's own LSF).
    '''
    x = fits.open(filename)
    names = [h.name for h in x]
    wave = np.array(x['WAVE'].data, dtype=float)
    if 'SKY_EAST' in names or 'SKY_WEST' in names:
        if ext not in names:
            raise ValueError('Extension %s not found in %s' % (ext, filename))
        flux2d = np.array(x[ext].data, dtype=float)
    elif 'FLUX' in names:
        flux2d = np.array(x['FLUX'].data, dtype=float)
    else:
        raise ValueError('Could not find a raw sky flux extension '
                         '(SKY_EAST/SKY_WEST or FLUX) in %s' % filename)
    lsf2d = np.array(x['LSF'].data, dtype=float) if 'LSF' in names else None
    x.close()
    if flux2d.ndim == 1:
        flux2d = flux2d[None, :]
    if lsf2d is not None and lsf2d.ndim == 1:
        lsf2d = lsf2d[None, :]
    return wave, flux2d, lsf2d


# ──────────────────────────────────────────────────────────────
# Per-line analysis
# ──────────────────────────────────────────────────────────────

def analyze_median(wave, flux2d, select, line_name, center, wmin, wmax, outroot,
                   lsf2d=None):
    '''
    Fit both profiles to the median spectrum (across select) for one
    line, and save a two-panel (fits + residuals) diagnostic PNG.

    Both profiles fit the line center freely (bounded only to
    [wmin, wmax]), so a real wavelength shift relative to the nominal
    line center is not assumed away -- the fitted centers are reported
    and annotated on the plot.

    lsf2d : optional (n_rows, n_wave) array (see load_raw_sky).  If
        given, the median (across select) of the file's own LSF
        extension at this line's wavelength is annotated on the plot
        and returned, for direct comparison against the two independent
        fits' own FWHM.

    Returns a dict of summary values, or None if either fit failed.
    '''
    med = np.nanmedian(flux2d[select], axis=0)
    ivar = estimate_ivar(med, wave)

    gfit = fit_gaussian(wave, med, ivar, center, wmin, wmax)
    mfit = fit_moffat(wave, med, ivar, center, wmin, wmax)
    if gfit is None or mfit is None:
        print('  %-10s median fit failed (gaussian=%s, moffat=%s)'
              % (line_name, gfit is not None, mfit is not None))
        return None

    gres, mres = gfit['result'], mfit['result']
    xg, yg = gfit['x'], gfit['y']
    xm, ym = mfit['x'], mfit['y']
    gy_fit = gaussian_profile(xg, **gres.best_values)
    my_fit = moffat_profile(xm, **mres.best_values)

    g_rms_core, g_rms_wing = core_wing_rms(xg, yg, gy_fit, center)
    m_rms_core, m_rms_wing = core_wing_rms(xm, ym, my_fit, center)
    m_fwhm = moffat_fwhm(mres.params['alpha'].value, mres.params['beta'].value)
    g_center = gres.params['center'].value
    m_center = mres.params['center'].value

    lsf_fwhm = None
    if lsf2d is not None:
        lsf_med = np.nanmedian(lsf2d[select], axis=0)
        lsf_fwhm = float(np.interp(center, wave, lsf_med))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6), sharex=True,
                                   gridspec_kw=dict(height_ratios=[3, 1]))
    ax1.plot(xg, yg, 'o', color='black', ms=3, label='data (median)')
    xx = np.linspace(wmin, wmax, 300)
    ax1.plot(xx, gaussian_profile(xx, **gres.best_values), '--', color='tab:blue',
             label='Gaussian (FWHM=%.2f A, shift=%+.3f A)'
             % (gres.params['fwhm'].value, g_center - center))
    ax1.plot(xx, moffat_profile(xx, **mres.best_values), '-', color='tab:red',
             label='Moffat (FWHM=%.2f A, beta=%.2f, shift=%+.3f A)'
             % (m_fwhm, mres.params['beta'].value, m_center - center))
    if lsf_fwhm is not None:
        ax1.axvline(center, color='grey', lw=0.6, ls=':')
        ax1.text(0.02, 0.95, 'LSF ext FWHM = %.2f A' % lsf_fwhm,
                 transform=ax1.transAxes, fontsize=8, va='top')
    ax1.set_ylabel('Flux')
    ax1.set_title('%s  (delta_AIC=%.1f, delta_BIC=%.1f)'
                  % (line_name, gres.aic - mres.aic, gres.bic - mres.bic))
    ax1.legend(fontsize=8)

    ax2.axhline(0, color='grey', lw=0.8)
    ax2.plot(xg, yg - gy_fit, '--', color='tab:blue', label='Gaussian resid')
    ax2.plot(xm, ym - my_fit, '-', color='tab:red', label='Moffat resid')
    ax2.set_xlabel('Wavelength (A)')
    ax2.set_ylabel('Residual')
    ax2.legend(fontsize=8)

    plt.tight_layout()
    outfile = '%s_median_%s.png' % (outroot, line_name)
    plt.savefig(outfile, dpi=120)
    plt.close(fig)

    return dict(
        line=line_name, n_rows=len(select),
        g_center=g_center, g_shift=g_center - center,
        g_fwhm=gres.params['fwhm'].value, g_aic=gres.aic, g_bic=gres.bic,
        g_rms_core=g_rms_core, g_rms_wing=g_rms_wing,
        m_center=m_center, m_shift=m_center - center,
        m_fwhm=m_fwhm, m_beta=mres.params['beta'].value,
        m_aic=mres.aic, m_bic=mres.bic,
        m_rms_core=m_rms_core, m_rms_wing=m_rms_wing,
        delta_aic=gres.aic - mres.aic, delta_bic=gres.bic - mres.bic,
        lsf_fwhm=lsf_fwhm,
        plot=outfile,
    )


def analyze_rows(wave, flux2d, select, line_name, center, wmin, wmax, lsf2d=None):
    '''
    Fit both profiles to each individual row in select for one line.
    Both profiles fit the line center freely (bounded to [wmin, wmax]);
    the fitted centers (and their shift from the nominal line center)
    are recorded per row.

    lsf2d : optional (n_rows, n_wave) array (see load_raw_sky).  If
        given, this row's own LSF FWHM at the line's wavelength is
        recorded as 'lsf_fwhm', for direct per-row comparison against
        the two independent fits' own FWHM.

    Returns a list of per-row result dicts (rows where either fit fails
    are simply omitted).
    '''
    records = []
    for i in select:
        row_flux = flux2d[i]
        if not np.any(np.isfinite(row_flux)):
            continue
        ivar = estimate_ivar(row_flux, wave)
        gfit = fit_gaussian(wave, row_flux, ivar, center, wmin, wmax)
        mfit = fit_moffat(wave, row_flux, ivar, center, wmin, wmax)
        if gfit is None or mfit is None:
            continue
        gres, mres = gfit['result'], mfit['result']
        gy_fit = gaussian_profile(gfit['x'], **gres.best_values)
        my_fit = moffat_profile(mfit['x'], **mres.best_values)
        g_rms_core, g_rms_wing = core_wing_rms(gfit['x'], gfit['y'], gy_fit, center)
        m_rms_core, m_rms_wing = core_wing_rms(mfit['x'], mfit['y'], my_fit, center)
        g_center = gres.params['center'].value
        m_center = mres.params['center'].value

        lsf_fwhm = (float(np.interp(center, wave, lsf2d[i]))
                   if lsf2d is not None else np.nan)

        records.append(dict(
            row=i, line=line_name,
            g_flux=gres.params['flux'].value,
            g_center=g_center, g_shift=g_center - center,
            g_fwhm=gres.params['fwhm'].value,
            g_aic=gres.aic, g_bic=gres.bic,
            g_rms_core=g_rms_core, g_rms_wing=g_rms_wing,
            m_flux=mres.params['flux'].value,
            m_center=m_center, m_shift=m_center - center,
            m_alpha=mres.params['alpha'].value, m_beta=mres.params['beta'].value,
            m_fwhm=moffat_fwhm(mres.params['alpha'].value, mres.params['beta'].value),
            m_aic=mres.aic, m_bic=mres.bic,
            m_rms_core=m_rms_core, m_rms_wing=m_rms_wing,
            delta_aic=gres.aic - mres.aic, delta_bic=gres.bic - mres.bic,
            lsf_fwhm=lsf_fwhm, fwhm_diff_g=gres.params['fwhm'].value - lsf_fwhm,
        ))
    return records


# ──────────────────────────────────────────────────────────────
# Batch processing
# ──────────────────────────────────────────────────────────────

def do_all(filename, lines=None, idelta=50, ext='SKY_EAST', outroot=''):
    '''
    Fit both profiles to the requested lines, for the median spectrum and
    for each individual selected row, and write the plots/table/summary
    described in the module docstring.
    '''
    if outroot == '':
        stem = os.path.splitext(os.path.basename(filename))[0]
        outroot = '%s_profile' % stem

    wave, flux2d, lsf2d = load_raw_sky(filename, ext=ext)
    n_rows = flux2d.shape[0]
    select = np.arange(0, n_rows, idelta)
    print('Using %d of %d rows (every %d-th)' % (len(select), n_rows, idelta))
    if lsf2d is None:
        print('Note: no LSF extension in this file; fitted FWHM will not be '
              'compared against a reference LSF value.')

    all_lines = SKY_LINES if lines is None else [l for l in SKY_LINES if l[0] in lines]
    if lines is not None:
        missing = set(lines) - {l[0] for l in all_lines}
        if missing:
            print('Warning: unknown line name(s) ignored: %s' % ', '.join(sorted(missing)))

    median_summaries = []
    all_row_records = []
    for name, center, wmin, wmax in all_lines:
        print('Fitting %s (%.2f A, window %.1f-%.1f)...' % (name, center, wmin, wmax))
        msum = analyze_median(wave, flux2d, select, name, center, wmin, wmax, outroot,
                              lsf2d=lsf2d)
        if msum is not None:
            median_summaries.append(msum)
        rows = analyze_rows(wave, flux2d, select, name, center, wmin, wmax, lsf2d=lsf2d)
        all_row_records.extend(rows)

    # per-row ascii table
    if all_row_records:
        rtab = Table(rows=all_row_records, names=list(all_row_records[0].keys()))
        for col in rtab.colnames:
            if col in ('row',):
                rtab[col].format = 'd'
            elif col != 'line':
                rtab[col].format = '.4g'
        perrow_file = '%s_perrow.txt' % outroot
        rtab.write(perrow_file, format='ascii.fixed_width_two_line', overwrite=True)
        print('Wrote %s' % perrow_file)
    else:
        rtab = None
        print('Warning: no per-row fits succeeded; %s_perrow.txt not written' % outroot)

    # aggregate summary plot + printed table
    if rtab is not None:
        line_names = [l[0] for l in all_lines if l[0] in set(rtab['line'])]
        nlines = len(line_names)
        have_lsf = lsf2d is not None

        if have_lsf:
            # one row per line: delta_AIC, Moffat beta, fit_fwhm - lsf_fwhm, center shift
            fig, axes = plt.subplots(nlines, 4, figsize=(16, 2.0 * nlines), squeeze=False)
            print('\n%-10s %6s %10s %8s %8s %8s %10s %10s' %
                 ('Line', 'N', 'med_dAIC', 'frac>2', 'med_beta', 'IQR_beta',
                  'med_dFWHM', 'med_shift'))
            print('-' * 92)
        else:
            ncols = 2
            nrows_fig = max(1, int(np.ceil(nlines / ncols)))
            fig, axes = plt.subplots(nrows_fig, ncols * 2, figsize=(4 * ncols * 2, 2.2 * nrows_fig))
            axes = np.atleast_2d(axes)
            print('\n%-10s %6s %10s %8s %8s %8s' %
                 ('Line', 'N', 'med_dAIC', 'frac>2', 'med_beta', 'IQR_beta'))
            print('-' * 60)

        for k, name in enumerate(line_names):
            sub = rtab[rtab['line'] == name]
            med_daic = float(np.median(sub['delta_aic']))
            frac_moffat = float(np.mean(sub['delta_aic'] > 2))
            med_beta = float(np.median(sub['m_beta']))
            iqr_beta = float(np.percentile(sub['m_beta'], 75) - np.percentile(sub['m_beta'], 25))

            if have_lsf:
                ax_aic, ax_beta, ax_fwhm, ax_shift = axes[k]
                ax_aic.hist(sub['delta_aic'], bins=20, color='tab:purple')
                ax_aic.axvline(0, color='k', lw=0.8)
                ax_aic.set_title('%s: delta_AIC' % name, fontsize=9)
                ax_beta.hist(sub['m_beta'], bins=20, color='tab:green')
                ax_beta.set_title('%s: Moffat beta' % name, fontsize=9)
                ax_fwhm.hist(sub['fwhm_diff_g'], bins=20, color='tab:orange')
                ax_fwhm.axvline(0, color='k', lw=0.8)
                ax_fwhm.set_title('%s: fit FWHM - LSF ext' % name, fontsize=9)
                ax_shift.hist(sub['g_shift'], bins=20, color='tab:cyan')
                ax_shift.axvline(0, color='k', lw=0.8)
                ax_shift.set_title('%s: center shift (A)' % name, fontsize=9)

                med_dfwhm = float(np.median(sub['fwhm_diff_g']))
                med_shift = float(np.median(sub['g_shift']))
                print('%-10s %6d %10.2f %8.2f %8.2f %8.2f %10.3f %10.4f' %
                     (name, len(sub), med_daic, frac_moffat, med_beta, iqr_beta,
                      med_dfwhm, med_shift))
            else:
                r, c = divmod(k, ncols)
                ax_aic = axes[r, 2 * c]
                ax_beta = axes[r, 2 * c + 1]
                ax_aic.hist(sub['delta_aic'], bins=20, color='tab:purple')
                ax_aic.axvline(0, color='k', lw=0.8)
                ax_aic.set_title('%s: delta_AIC' % name, fontsize=9)
                ax_beta.hist(sub['m_beta'], bins=20, color='tab:green')
                ax_beta.set_title('%s: Moffat beta' % name, fontsize=9)
                print('%-10s %6d %10.2f %8.2f %8.2f %8.2f' %
                     (name, len(sub), med_daic, frac_moffat, med_beta, iqr_beta))

        plt.tight_layout()
        summary_file = '%s_summary.png' % outroot
        plt.savefig(summary_file, dpi=110)
        plt.close(fig)
        print('\nWrote %s' % summary_file)

    return median_summaries, rtab


# ──────────────────────────────────────────────────────────────
# Command-line entry point
# ──────────────────────────────────────────────────────────────

if __name__ == '__main__':
    argv = sys.argv[1:]
    if not argv or '-h' in argv or '--help' in argv:
        print(_USAGE)
        sys.exit(0)

    idelta   = 50
    lines    = None
    ext      = 'SKY_EAST'
    outroot  = ''
    filename = None

    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == '-delta':
            i += 1
            idelta = int(argv[i])
        elif arg == '-lines':
            i += 1
            lines = [s.strip() for s in argv[i].split(',') if s.strip()]
        elif arg == '-ext':
            i += 1
            ext = argv[i]
        elif arg == '-out':
            i += 1
            outroot = argv[i]
        elif arg.startswith('-'):
            print('Error: unknown option "%s"' % arg)
            print(_USAGE)
            sys.exit(1)
        else:
            if filename is not None:
                print('Error: unexpected argument "%s" (filename already set to "%s")'
                      % (arg, filename))
                print(_USAGE)
                sys.exit(1)
            filename = arg
        i += 1

    if filename is None:
        print('Error: no filename supplied')
        print(_USAGE)
        sys.exit(1)

    if not os.path.exists(filename):
        print('Error: file not found: %s' % filename)
        sys.exit(1)

    do_all(filename, lines=lines, idelta=idelta, ext=ext, outroot=outroot)
