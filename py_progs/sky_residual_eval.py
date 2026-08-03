#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Method-agnostic evaluation of a sky (or continuum) subtraction: given
    an observed spectrum and whatever model was subtracted from it --
    a full ESO/PALACE/SkyDecomp decomposition, or something as simple as
    a scaled sky-fiber spectrum -- quantify the residual separately in
    continuum bands and at individual airglow lines, and decompose each
    line's residual into an amplitude error, a wavelength-registration
    error, and an LSF-width error.

    Works on a single (wave, flux, model) triple or, more usefully, on
    arrays of them (one row per exposure or per fiber), so the same
    numbers can be produced directly from an XCframe/XSFrame summary
    table or from a single SFrame image's FLUX/SKY extensions.

    Three companion plotting routines (plot_frac_summary,
    plot_continuum_summary, plot_lines_summary) turn a batch
    summary_table into one-page PNG grids -- quality-fraction
    distributions, continuum diagnostics, and line shape-fit
    diagnostics respectively, each with one row per metric and one
    column per spectrograph arm -- written to a plots_sky_resid/
    subdirectory by default.

Command line usage (if any)::

    sky_residual_eval.py filename [-mask mask.fits] [-nproc N] [-out ROOT]
                        [-plotdir DIR]

    Arguments::

        filename    FITS file with WAVE, FLUX, and SKY extensions (the
                    convention used by SkySubOrig/Drp/Dev1/Dev2.py and
                    read by SkySub_eval.py).  WAVE is 1-D; FLUX/SKY are
                    1-D (single spectrum) or 2-D, n_rows x n_wave.  An
                    IVAR extension is used if present.  NOTE: in this
                    convention FLUX is already sky-subtracted (per
                    SkySub_eval.py: "original flux (FLUX+SKY);
                    sky-subtracted (FLUX)") -- the "observed" spectrum
                    passed to analyze_sky_residuals is reconstructed as
                    FLUX+SKY, with SKY as the model, so the internal
                    flux-model subtraction reduces back to FLUX exactly
                    rather than double-subtracting.

    Options::

        -mask mask.fits  palace_make_mask.py output (WAVE/MASK
                         extensions).  Default: data/sky_mask.fits (the
                         repo's own, same as SkyObsESOCompare.py's
                         default).  Only falls back to deriving clean
                         pixels from the line list's own windows if
                         that file is also missing.
        -nproc N         worker processes for batch rows (default 1).
        -out ROOT        output filename root; writes
                         <ROOT>_summary.fits and <ROOT>_lines.fits
                         (default: <stem of filename>).
        -plotdir DIR     directory for summary plots (default:
                         plots_sky_resid); see plot_frac_summary.

Description:

    ``analyze_sky_residual`` (single spectrum) and ``analyze_sky_residuals``
    (batch) are the reusable entry points; the command-line form is a thin
    smoke-test / convenience wrapper around the batch call.

    residual = flux - model, decomposed two ways:

    1. Continuum bands (default: GetSkyCont.ARM_EVAL_RANGES, i.e. B/R/Z
       with the arm-overlap zones excluded).  Over clean pixels in each
       band: CONT_OFFSET and CONT_NMAD/CONT_RMS (via
       GetSkyCont.arm_continuum_stats); CONT_ALPHA, a dimensionless
       power-law-index mismatch (NOT a flux/Angstrom slope) -- treating
       the true continuum as the model times a small shape correction,
       ``true(w) = model(w)*(w/w_ref)**da``, to first order
       ``resid/model ~= da*ln(w/w_ref)``, so a weighted linear fit of
       (resid/model) against ln(wavelength/w_ref) gives ``da`` directly,
       comparable across bands/exposures of very different absolute
       brightness; CONT_FIT_QUALITY = NMAD/NOISE_PROXY (NOT a formal
       ivar chi-square -- see analyze_sky_residual's Notes for why that
       saturates uselessly for real data); and the fraction of
       individual clean pixels with
       ``|residual|`` AT OR BELOW two fixed absolute flux levels (1e-14,
       1e-15 erg/s/cm^2/Ang) -- CONT_FRAC_1E14/1E15 -- a quality fraction
       (higher = more of the continuum is well subtracted), not a failure
       count.

    2. Individual lines (default: _DEFAULT_LINES -- 16 lines built from
       the DRP's own REF_SKYLINES (see lvmdrp.core.constants) minus
       9719.84 A (outside ARM_EVAL_RANGES' Z band, so it never
       contributed to a per-band aggregate anyway) plus 9 of
       sky_gaussfit.py's SKY_LINES not already in REF_SKYLINES -- all
       except sky6300, excluded by choice, not for any technical reason
       -- ksl, 260803).  For each line, a local Gaussian is fit to the
       MODEL (not the data) to get its own (A0, lam0, sigma0) -- "what
       the model claims this line looks like."  The residual in that
       window is then fit against the analytic basis
       ``[G0, A0*dG0/dlam0, A0*dG0/dsigma, 1, (wave-lam0)]``, a first-order
       Taylor expansion of the line profile in amplitude, center, and
       width.  The three basis coefficients ARE, to first order, the
       physical corrections needed to make the model match the data::

           DELTA_A      flux units   -- amplitude/flux-scale error
           DELTA_LAM    Angstroms    -- wavelength-registration error
           DELTA_SIGMA  Angstroms    -- LSF-width error (the number that
                                        answers "is a small LSF mismatch
                                        contributing to the residual")

       AMP_RATIO = A0_MODEL / (A0_MODEL + DELTA_A) -- predicted/measured
       amplitude ratio (1.0 = exact match; <1 = model under-predicts;
       >1 = over-predicts) -- is derived from DELTA_A and used instead
       of it for aggregation, since DELTA_A is in absolute flux units
       and real line brightness varies so much line to line that a
       flux-unit average is dominated by whichever line happens to be
       brightest (ksl, 260803).  Aggregated per band as
       LINE_AMPRATIO/DLAM/DSIG_MED/SCATTER (DELTA_A itself is still kept
       per-line in line_detail).

    3. LINE_FRAC_10SIG/1E15/1E14: the same fixed-threshold quality
       fraction as the continuum side (CONT_FRAC_1E14/1E15) but over
       line-affected pixels instead of clean ones -- the fraction of
       individual pixels flagged non-clean by ``mask`` whose ``|residual|``
       is AT OR BELOW each threshold.  This deliberately uses the MASK, NOT
       ``line_list`` -- line_list is a handful of named lines
       (wavelength-cal/LSF questions, strong-line stats; see (2) above
       and Notes), too coarse a set for a fraction statistic to be
       anything but a step function.  (An earlier version used the
       fraction of contiguous line REGIONS instead, one vote per
       region regardless of width, matching SkyObsESOCompare.py's own
       convention; in practice this still showed visible digitization
       at the region counts real masks produce, so per-pixel is used
       instead.)  LINE_FRAC_NPIX_<band> records how many pixels the
       fraction was computed over.  Without a real ``mask`` (clean
       pixels derived from line_list's own windows instead), this
       necessarily degrades to line_list's own handful of windows.

Notes:

    The ``fwhm`` argument, if given, is a callable ``fwhm(wave_array) ->
    FWHM_array`` in Angstroms.  It sets each line's fit-window half-width
    (adaptively, since LVM's LSF is not constant across the range) and
    seeds the local Gaussian fit to the model.  Without it, a fixed 5 A
    half-window and a generic sigma guess are used instead, and affected
    line_detail rows are flagged ``no_fwhm_prior`` -- DELTA_SIGMA is
    still returned but with lower confidence on weak lines.

    ``mask``, if given, is either a (wave, bool) pair as returned by
    GetSkyCont.load_mask, or a filename to pass to it -- i.e. a
    palace_make_mask.py output.  This is deliberately the SAME mask
    GetSkyCont.py/SkyObsESOCompare.py already use, not a new one.  If
    omitted (the normal case), DEFAULT_MASK_FILE (data/sky_mask.fits,
    same default as SkyObsESOCompare.py/SkyObsESOStack.py) is used
    automatically; the line-list-derived fallback only applies if that
    file is also missing.

    The per-line derivative-basis fit uses the MODEL's own realized line
    shape as the reference, not an assumed universal instrumental
    profile -- so it is agnostic to which decomposition engine (ESO,
    PALACE, SkyDecomp, or a naive scaled-fiber subtraction) produced
    ``model``.

    CONT_FIT_QUALITY is NMAD/NOISE_PROXY, not a formal ivar chi-square
    (checked directly against a real ivar-bearing file, 260803: median
    chi2/N came out 667-3219, because real sky-subtraction residual --
    both its systematic offset and its scatter -- routinely exceeds the
    formal photon-noise floor implied by ivar by close to an order of
    magnitude, so a formal chi2 saturates at huge values for essentially
    every real exposure rather than discriminating good from bad).
    NOISE_PROXY (from consecutive-pixel differences in flux, same
    formula as SkyObsESOCompare.py's own) needs no ivar at all, so this
    metric is always computable, unlike the chi2 it replaced.

History::

    260803  ksl  First version.  Method-agnostic (wave, flux, model) ->
                 (summary, line_detail) sky-subtraction evaluator, one
                 spectrum or batched (multiprocessing.Pool).  Continuum
                 bands: offset, NMAD/RMS, dimensionless power-law-index
                 mismatch (ALPHA), NMAD/NOISE_PROXY fit quality (not a
                 formal ivar chi-square, which saturates uselessly for
                 real sky-subtraction residual), and fixed-threshold
                 quality fractions.  Lines: per-line derivative-basis
                 decomposition (model vs data) into DELTA_A/AMP_RATIO,
                 DELTA_LAM, DELTA_SIGMA, against a 16-line default list
                 (REF_SKYLINES plus 9 of sky_gaussfit.py's SKY_LINES);
                 line-region quality fractions from the PALACE mask
                 (data/sky_mask.fits default), independent of the line
                 list.  Three summary plots (plot_frac_summary,
                 plot_continuum_summary, plot_lines_summary) written to
                 plots_sky_resid/.
'''

import sys
import re
import multiprocessing
from pathlib import Path

# Force a non-interactive backend BEFORE importing pyplot, but only when
# run as a script -- so PNGs can be written with no display available,
# without breaking a notebook/interactive import (same convention as
# lvm_line_profile.py/SkyObsESO_analysis.py).
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table, vstack
from scipy.optimize import curve_fit

from GetSkyCont import load_mask, _interp_mask_to_wave, ARM_EVAL_RANGES, arm_continuum_stats

# Same convention as SkyObsESOCompare.py/SkyObsESOStack.py: default to the
# repo's own clean-pixel mask rather than requiring -mask every time.
DEFAULT_MASK_FILE = Path(__file__).resolve().parent.parent / 'data' / 'sky_mask.fits'


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before a line consisting of "History:" (or
    "History::"/"Version History" -- whitespace/colon-insensitive), so
    -h stays short even as that section grows.
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


# ──────────────────────────────────────────────────────────────
# Self-contained copies (deliberately not imported -- see lvm_line_profile.py
# for the same choice/rationale): SkyObsESOCompare.py pulls in SkySepESO's
# ESO-sky-model machinery at import time, which this module has no other
# need for.  These three are stable, tiny, and duplicated with attribution
# rather than dragging in that dependency chain.
# ──────────────────────────────────────────────────────────────

def _noise_proxy(flux, clean):
    '''
    Robust per-spectrum noise estimate from consecutive-pixel differences
    within clean pixels (copy of SkyObsESOCompare._noise_proxy)::

        sigma = 1.4826 * median(|diff|) / sqrt(2)
    '''
    vals = np.asarray(flux)[clean]
    vals = vals[np.isfinite(vals)]
    if len(vals) < 3:
        return np.nan
    diffs = np.diff(vals)
    return float(1.4826 * np.median(np.abs(diffs)) / np.sqrt(2))


# Copy of SkyObsESOCompare._LINE_FRAC_SIGMA / _LINE_FRAC_ABS_THRESH.
_LINE_FRAC_SIGMA = 10.0
_LINE_FRAC_ABS_THRESH = (1e-15, 1e-14)


# ──────────────────────────────────────────────────────────────
# Default line list.  REF_SKYLINES is the DRP's own list
# (lvmdrp.core.constants), copied here (not imported -- lvmdrp needs its
# own conda env) so this module has no lvmdrp dependency: "GB hand
# picked isolated bright lines across each channel which are not
# doublets in UVES atlas; true wavelengths taken from UVES sky line
# atlas" -- used there for SKYLINES_FIBERFLAT and wavelength-
# registration checks.  Kept here for reference/attribution even though
# it's no longer what default_line_list() builds from (see below).
# ──────────────────────────────────────────────────────────────

REF_SKYLINES = {
    'b': [5577.346680],
    'r': [6363.782715, 7358.680176, 7392.209961],
    'z': [8399.175781, 8988.383789, 9552.546875, 9719.838867],
}

# The actual default: REF_SKYLINES minus 9719.838867 (falls outside
# ARM_EVAL_RANGES' Z-band evaluation range, so it was never
# contributing to any per-band aggregate anyway -- see History), plus 9
# of sky_gaussfit.py's own SKY_LINES not already in REF_SKYLINES (all
# except sky6300 -- excluded at ksl's request, 260803, no reason given
# beyond preference).  sky6553's wavelength here (6553.617) is
# sky_gaussfit.py's own corrected value (see that file's History), not
# the naive 6553.0.
_DEFAULT_LINES = [
    ('sky5577', 5577.346680),
    ('sky6363', 6363.782715),
    ('sky6533', 6533.04),
    ('sky6553', 6553.617),
    ('sky6577', 6577.2),
    ('sky6912', 6912.623),
    ('sky6923', 6923.220),
    ('sky6939', 6939.521),
    ('sky7358', 7358.680176),
    ('sky7392', 7392.209961),
    ('sky7914', 7913.717773),
    ('sky8344', 8344.613281),
    ('sky8399', 8399.175781),
    ('sky8827', 8827.112305),
    ('sky8988', 8988.383789),
    ('sky9552', 9552.546875),
]


def default_line_list(fwhm=None):
    '''
    Build the default line list from _DEFAULT_LINES.

    Parameters
    ----------
    fwhm : callable or None
        fwhm(wave_array) -> FWHM_array (Angstroms).  If given, each
        line's half-width is max(5, 4*FWHM); otherwise a fixed 5 A.

    Returns
    -------
    list of dict, each with keys name, wave0, half_width.
    '''
    lines = []
    for name, w0 in _DEFAULT_LINES:
        half_width = 5.0
        if fwhm is not None:
            try:
                half_width = max(5.0, 4.0 * float(np.atleast_1d(fwhm(np.array([w0])))[0]))
            except Exception:
                half_width = 5.0
        lines.append({'name': name, 'wave0': float(w0), 'half_width': float(half_width)})
    lines.sort(key=lambda d: d['wave0'])
    return lines


def _gaussian(wave, amp, center, sigma, bkg):
    return amp * np.exp(-0.5 * ((wave - center) / sigma) ** 2) + bkg


def _fit_line_to_model(wave, model, wave0, half_width, fwhm=None):
    '''
    Fit a Gaussian + constant background to the MODEL spectrum in a
    window around wave0.  Returns a dict with sel (boolean selection),
    amp, lam0, sigma0, or None if the window has too few pixels or the
    fit fails.
    '''
    sel = (wave >= wave0 - half_width) & (wave <= wave0 + half_width)
    if sel.sum() < 5:
        return None
    w, m = wave[sel], model[sel]
    finite = np.isfinite(w) & np.isfinite(m)
    if finite.sum() < 5:
        return None
    w, m = w[finite], m[finite]

    bkg_guess = float(np.nanmedian(m))
    amp_guess = float(np.nanmax(m) - bkg_guess)
    if not np.isfinite(amp_guess) or amp_guess <= 0:
        return None
    sigma_guess = 0.6
    if fwhm is not None:
        try:
            sigma_guess = max(0.1, float(np.atleast_1d(fwhm(np.array([wave0])))[0]) / 2.355)
        except Exception:
            sigma_guess = 0.6

    # Real LVM flux is ~1e-14 to 1e-16 -- scipy's bounded 'trf' solver uses a
    # finite-difference step relative to each parameter's own magnitude, and
    # at that scale the step underflows and the amplitude/background never
    # move from p0 (confirmed: unscaled fit silently returns p0 unchanged).
    # Rescaling flux to order-unity before the fit avoids this entirely.
    scale = amp_guess
    p0 = [1.0, wave0, sigma_guess, bkg_guess / scale]
    bounds = ([0.0, wave0 - half_width, 0.05, -np.inf],
              [np.inf, wave0 + half_width, 5.0, np.inf])
    try:
        popt, _ = curve_fit(_gaussian, w, m / scale, p0=p0, bounds=bounds, maxfev=2000)
    except Exception:
        return None
    amp, lam0, sigma0, bkg = popt
    amp, bkg = amp * scale, bkg * scale
    if not (np.isfinite(amp) and np.isfinite(lam0) and np.isfinite(sigma0)) or amp <= 0:
        return None
    return dict(sel=sel & finite_mask_full(wave, finite, sel), amp=float(amp),
                lam0=float(lam0), sigma0=float(sigma0))


def finite_mask_full(wave, finite_within_sel, sel):
    '''Expand a finite-mask computed on wave[sel] back to full-length boolean array.'''
    out = np.zeros(len(wave), dtype=bool)
    idx = np.flatnonzero(sel)
    out[idx[finite_within_sel]] = True
    return out


def _line_shape_fit(wave, resid, ivar, sel, amp, lam0, sigma0):
    '''
    Weighted linear least-squares fit of resid[sel] against the
    derivative basis ``[G0, amp*dG0/dlam0, amp*dG0/dsigma, 1, (wave-lam0)]``.

    Returns dict with delta_a/_err, delta_lam/_err, delta_sigma/_err,
    chi2 (reduced, NaN if no ivar), npix -- or None if underdetermined.
    '''
    w, r = wave[sel], resid[sel]
    if ivar is not None:
        wt = ivar[sel]
        wt = np.where(np.isfinite(wt) & (wt > 0), wt, 0.0)
    else:
        wt = np.ones_like(r)

    g0 = np.exp(-0.5 * ((w - lam0) / sigma0) ** 2)
    dg_dlam = g0 * (w - lam0) / sigma0 ** 2
    dg_dsig = g0 * (w - lam0) ** 2 / sigma0 ** 3
    X = np.column_stack([g0, amp * dg_dlam, amp * dg_dsig,
                         np.ones_like(w), (w - lam0)])

    good = np.isfinite(r) & np.isfinite(wt) & np.all(np.isfinite(X), axis=1) & (wt > 0)
    if good.sum() < X.shape[1] + 2:
        return None
    Xg, rg, wtg = X[good], r[good], wt[good]

    # Column equilibration: amp*dG/dlam0 and amp*dG/dsigma sit ~1e-14 in
    # real flux units while G0/const/linear are O(1-10) -- inverting
    # X.T@X directly on that dynamic range risks garbage covariance even
    # though the coefficient solve itself (SVD-based lstsq) is robust to
    # it.  Scale every column to O(1), solve, then unscale.
    col_scale = np.max(np.abs(Xg), axis=0)
    col_scale = np.where(col_scale > 0, col_scale, 1.0)
    Xn = Xg / col_scale

    sw = np.sqrt(wtg)
    Xw, rw = Xn * sw[:, None], rg * sw

    try:
        coef_n, _, rank, _ = np.linalg.lstsq(Xw, rw, rcond=None)
    except np.linalg.LinAlgError:
        return None
    if rank < X.shape[1]:
        return None

    dof = len(rw) - X.shape[1]
    if dof <= 0:
        return None
    resid_w = rw - Xw @ coef_n
    ssq = float(np.sum(resid_w ** 2))
    try:
        XtX_inv = np.linalg.inv(Xw.T @ Xw)
    except np.linalg.LinAlgError:
        return None

    if ivar is not None:
        chi2 = ssq / dof
        cov_n = XtX_inv
    else:
        chi2 = np.nan
        cov_n = (ssq / dof) * XtX_inv
    errs_n = np.sqrt(np.clip(np.diag(cov_n), 0, None))

    coef = coef_n / col_scale
    errs = errs_n / col_scale

    return dict(delta_a=float(coef[0]), delta_a_err=float(errs[0]),
                delta_lam=float(coef[1]), delta_lam_err=float(errs[1]),
                delta_sigma=float(coef[2]), delta_sigma_err=float(errs[2]),
                chi2=float(chi2), npix=int(good.sum()))


def _band_for(wave0, bands):
    for label, lo, hi in bands:
        if lo <= wave0 <= hi:
            return label
    return 'NONE'


_LINE_DETAIL_COLUMNS = ['LINE_NAME', 'BAND', 'WAVE0', 'A0_MODEL', 'SIGMA0_MODEL',
                        'DELTA_A', 'DELTA_A_ERR', 'AMP_RATIO', 'AMP_RATIO_ERR',
                        'DELTA_LAM', 'DELTA_LAM_ERR',
                        'DELTA_SIGMA', 'DELTA_SIGMA_ERR', 'CHI2', 'N_PIX_WINDOW',
                        'MED_ABS_RESID', 'FLAG']


def analyze_sky_residual(wave, flux, model, ivar=None, fwhm=None,
                          line_list=None, bands=None, mask=None):
    '''
    Evaluate one (wave, flux, model) triple.  See module docstring.

    Returns
    -------
    summary : dict
        Per-band continuum and line-aggregate quantities.
    line_detail : astropy.table.Table
        One row per requested line (see _LINE_DETAIL_COLUMNS).
    '''
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    model = np.asarray(model, dtype=float)
    ivar = None if ivar is None else np.asarray(ivar, dtype=float)
    resid = flux - model

    if bands is None:
        bands = ARM_EVAL_RANGES
    if line_list is None:
        line_list = default_line_list(fwhm)

    if mask is None and DEFAULT_MASK_FILE.exists():
        mask = DEFAULT_MASK_FILE

    if mask is None:
        # No real mask available (DEFAULT_MASK_FILE missing and none
        # given) -- fall back to line_list's own windows.  Necessarily
        # coarse (see LINE_FRAC note in module docstring); a real mask
        # should be used whenever possible.
        clean = np.ones(len(wave), dtype=bool)
        for ln in line_list:
            clean &= ~((wave >= ln['wave0'] - ln['half_width']) &
                      (wave <= ln['wave0'] + ln['half_width']))
    else:
        if isinstance(mask, (str, Path)):
            mask_wave, mask_bool = load_mask(mask)
        else:
            mask_wave, mask_bool = mask
        clean = _interp_mask_to_wave(mask_wave, mask_bool, wave)

    summary = {}

    # ---- continuum bands ----
    cstats = arm_continuum_stats(wave, resid, clean=clean, arm_ranges=bands)
    noise_by_band = {}
    for label, lo, hi in bands:
        st = cstats[label]
        summary['CONT_NPIX_%s' % label] = st['n']
        summary['CONT_OFFSET_%s' % label] = st['med']
        summary['CONT_NMAD_%s' % label] = st['nmad']
        summary['CONT_RMS_%s' % label] = st['rms']

        sel = (wave >= lo) & (wave <= hi) & clean
        noise_by_band[label] = _noise_proxy(flux, sel)
        vals, wv, mv = resid[sel], wave[sel], model[sel]
        finite = np.isfinite(vals) & np.isfinite(wv)
        vals, wv, mv = vals[finite], wv[finite], mv[finite]
        wt = (ivar[sel][finite] if ivar is not None else np.ones_like(vals))
        wt = np.where(np.isfinite(wt) & (wt > 0), wt, 0.0)

        # CONT_ALPHA: dimensionless power-law-index mismatch, not a
        # flux/Angstrom slope.  Picture the true continuum as the model
        # times a small shape correction, true(w) = model(w)*(w/w_ref)^da;
        # to first order resid/model = true/model - 1 ~= da*ln(w/w_ref),
        # so a weighted linear fit of (resid/model) against ln(w/w_ref)
        # gives da directly -- comparable across bands/exposures of very
        # different absolute brightness, unlike a flux/Ang slope (ksl,
        # 260803).
        good_m = np.isfinite(mv) & (mv > 0)
        if np.sum(good_m) >= 3:
            wv_ref = wv[good_m].mean()
            xvar = np.log(wv[good_m] / wv_ref)
            yvar = vals[good_m] / mv[good_m]
            # Var(resid/model) ~= Var(resid)/model**2 (error propagation for
            # a ratio at fixed model), so the correct weight is wt*model**2,
            # not wt alone -- without the model**2 factor, a handful of
            # near-zero-continuum "clean" pixels (real gaps between sky
            # lines, not outliers) get equal weight to bright ones despite
            # their resid/model ratio being enormously noisier, and can
            # single-handedly swing the whole fit to an unphysical alpha.
            w_here = wt[good_m] * mv[good_m] ** 2
            try:
                alpha = float(np.polyfit(xvar, yvar, 1, w=np.sqrt(np.clip(w_here, 0, None)))[0])
            except Exception:
                alpha = np.nan
        else:
            alpha = np.nan
        summary['CONT_ALPHA_%s' % label] = alpha

        # CONT_FIT_QUALITY: NMAD/NOISE_PROXY, not a formal ivar chi-square.
        # A real ivar-based chi2 (sum(ivar*resid**2)/N) saturates at huge
        # values (100s-1000s) for essentially every real exposure, because
        # actual sky-subtraction residual (systematic model error) routinely
        # exceeds the formal photon-noise floor implied by ivar by close to
        # an order of magnitude -- it doesn't discriminate good exposures
        # from bad ones, just confirms "residual != 0" every time (ksl,
        # 260803).  NMAD/NOISE_PROXY is the same fix SkyObsESOCompare.py's
        # own FIT_QUALITY already uses for this exact saturation problem:
        # ~1 means the residual looks like pixel-to-pixel noise; higher
        # means real excess structure.  Works even without ivar (NOISE_PROXY
        # is derived from flux alone), so this metric no longer requires an
        # IVAR extension to be non-NaN.
        nz = noise_by_band[label]
        st_nmad = st['nmad']
        summary['CONT_FIT_QUALITY_%s' % label] = (
            float(st_nmad / nz) if (np.isfinite(st_nmad) and np.isfinite(nz) and nz > 0) else np.nan)

        if len(vals):
            summary['CONT_FRAC_1E14_%s' % label] = float(np.mean(np.abs(vals) <= 1e-14))
            summary['CONT_FRAC_1E15_%s' % label] = float(np.mean(np.abs(vals) <= 1e-15))
        else:
            summary['CONT_FRAC_1E14_%s' % label] = np.nan
            summary['CONT_FRAC_1E15_%s' % label] = np.nan

    # ---- lines ----
    line_rows = []
    for ln in line_list:
        name, wave0, hw = ln['name'], ln['wave0'], ln['half_width']
        band = _band_for(wave0, bands)
        row = dict(LINE_NAME=name, BAND=band, WAVE0=wave0,
                  A0_MODEL=np.nan, SIGMA0_MODEL=np.nan,
                  DELTA_A=np.nan, DELTA_A_ERR=np.nan,
                  AMP_RATIO=np.nan, AMP_RATIO_ERR=np.nan,
                  DELTA_LAM=np.nan, DELTA_LAM_ERR=np.nan,
                  DELTA_SIGMA=np.nan, DELTA_SIGMA_ERR=np.nan,
                  CHI2=np.nan, N_PIX_WINDOW=0, MED_ABS_RESID=np.nan, FLAG='')

        fitres = _fit_line_to_model(wave, model, wave0, hw, fwhm)
        if fitres is None:
            row['FLAG'] = 'model_fit_failed'
            line_rows.append(row)
            continue

        sel = fitres['sel']
        row['A0_MODEL'] = fitres['amp']
        row['SIGMA0_MODEL'] = fitres['sigma0']
        row['N_PIX_WINDOW'] = int(sel.sum())
        if sel.sum():
            row['MED_ABS_RESID'] = float(np.nanmedian(np.abs(resid[sel])))

        shp = _line_shape_fit(wave, resid, ivar, sel, fitres['amp'], fitres['lam0'], fitres['sigma0'])
        flags = []
        if shp is None:
            flags.append('shape_fit_failed')
        else:
            row.update(DELTA_A=shp['delta_a'], DELTA_A_ERR=shp['delta_a_err'],
                      DELTA_LAM=shp['delta_lam'], DELTA_LAM_ERR=shp['delta_lam_err'],
                      DELTA_SIGMA=shp['delta_sigma'], DELTA_SIGMA_ERR=shp['delta_sigma_err'],
                      CHI2=shp['chi2'])
            # AMP_RATIO = predicted/measured = A0_MODEL / (A0_MODEL + DELTA_A),
            # since (A0_MODEL + DELTA_A) is the first-order estimate of the
            # TRUE (observed) line amplitude -- a ratio is comparable across
            # lines of very different brightness, unlike DELTA_A itself
            # (flux units), which is dominated by whichever line happens to
            # be brightest when aggregated (ksl, 260803).  1.0 = perfect
            # match; <1 = model under-predicts; >1 = model over-predicts.
            measured_amp = fitres['amp'] + shp['delta_a']
            if measured_amp > 0:
                row['AMP_RATIO'] = fitres['amp'] / measured_amp
                row['AMP_RATIO_ERR'] = (fitres['amp'] / measured_amp ** 2) * shp['delta_a_err']
            else:
                flags.append('measured_amp_nonpositive')
        if fwhm is None:
            flags.append('no_fwhm_prior')
        row['FLAG'] = ';'.join(flags)
        line_rows.append(row)

    line_detail = (Table(rows=line_rows, names=_LINE_DETAIL_COLUMNS)
                  if line_rows else Table(names=_LINE_DETAIL_COLUMNS,
                                          dtype=[object, object, float, float, float,
                                                float, float, float, float, float, float,
                                                float, float, float, int, float, object]))

    # ---- per-band line aggregates ----
    for label, lo, hi in bands:
        in_band = line_detail[line_detail['BAND'] == label] if len(line_detail) else line_detail
        good = in_band[np.isfinite(in_band['DELTA_A'])] if len(in_band) else in_band
        summary['LINE_N_%s' % label] = int(len(good))

        for key, col in (('AMPRATIO', 'AMP_RATIO'), ('DLAM', 'DELTA_LAM'), ('DSIG', 'DELTA_SIGMA')):
            vals = np.asarray(good[col], dtype=float) if len(good) else np.array([])
            vals = vals[np.isfinite(vals)]
            if len(vals):
                med = float(np.median(vals))
                summary['LINE_%s_MED_%s' % (key, label)] = med
                summary['LINE_%s_SCATTER_%s' % (key, label)] = float(
                    1.4826 * np.median(np.abs(vals - med)))
            else:
                summary['LINE_%s_MED_%s' % (key, label)] = np.nan
                summary['LINE_%s_SCATTER_%s' % (key, label)] = np.nan

        chi2vals = np.asarray(good['CHI2'], dtype=float) if len(good) else np.array([])
        chi2vals = chi2vals[np.isfinite(chi2vals)]
        summary['LINE_CHI2_%s' % label] = float(np.mean(chi2vals)) if len(chi2vals) else np.nan

        # LINE_FRAC_* is the fraction of individual line-affected PIXELS
        # (from the MASK, not line_list -- see module Notes) with
        # |residual| exceeding each threshold -- the same definition as
        # CONT_FRAC_*, just over (~clean) instead of clean.  An earlier
        # version used the fraction of contiguous line REGIONS instead
        # (matching SkyObsESOCompare.py's convention); with this file's
        # actual region counts (order 100/band) that produced visible
        # digitization in the distribution, so per-pixel is used here.
        band_sel = (wave >= lo) & (wave <= hi)
        line_pix = (~clean) & band_sel
        vals = np.abs(resid[line_pix])
        vals = vals[np.isfinite(vals)]
        summary['LINE_FRAC_NPIX_%s' % label] = int(len(vals))
        if len(vals):
            nz = noise_by_band.get(label, np.nan)
            sig_thresh = _LINE_FRAC_SIGMA * nz if np.isfinite(nz) else np.nan
            summary['LINE_FRAC_10SIG_%s' % label] = (
                float(np.mean(vals <= sig_thresh)) if np.isfinite(sig_thresh) else np.nan)
            abs1, abs2 = _LINE_FRAC_ABS_THRESH
            summary['LINE_FRAC_1E15_%s' % label] = float(np.mean(vals <= abs1))
            summary['LINE_FRAC_1E14_%s' % label] = float(np.mean(vals <= abs2))
        else:
            summary['LINE_FRAC_10SIG_%s' % label] = np.nan
            summary['LINE_FRAC_1E15_%s' % label] = np.nan
            summary['LINE_FRAC_1E14_%s' % label] = np.nan

    if ivar is not None:
        finite = np.isfinite(resid) & np.isfinite(ivar) & (ivar > 0)
        summary['CHI2_TOTAL'] = float(np.sum(ivar[finite] * resid[finite] ** 2))
        summary['DOF'] = int(finite.sum())
    else:
        summary['CHI2_TOTAL'] = np.nan
        summary['DOF'] = np.nan

    return summary, line_detail


def _analyze_one_row(wave, flux, model, ivar, fwhm, line_list, bands, mask):
    return analyze_sky_residual(wave, flux, model, ivar=ivar, fwhm=fwhm,
                                line_list=line_list, bands=bands, mask=mask)


def analyze_sky_residuals(wave, flux, model, ivar=None, fwhm=None,
                          line_list=None, bands=None, mask=None, nproc=1):
    '''
    Batch form of analyze_sky_residual.

    Parameters
    ----------
    wave : (n_wave,) or (n_rows, n_wave)
        Shared (1-D) or per-row (2-D) wavelength grid.
    flux, model : (n_wave,) or (n_rows, n_wave)
    ivar : same shape as flux, or None
    nproc : int
        Worker processes for the row loop (default 1 = sequential).
    Other parameters as analyze_sky_residual; line_list/bands are
    resolved once (not per row) for consistency and speed.

    Returns
    -------
    summary_table : astropy.table.Table
        One row per spectrum (plus ROW if batch).
    line_table : astropy.table.Table
        One row per (spectrum, line) (plus ROW if batch).
    '''
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    model = np.asarray(model, dtype=float)

    if flux.ndim == 1:
        summary, line_detail = analyze_sky_residual(wave, flux, model, ivar=ivar, fwhm=fwhm,
                                                     line_list=line_list, bands=bands, mask=mask)
        return Table(rows=[summary]), line_detail

    n = flux.shape[0]
    wave_rows = wave if wave.ndim == 2 else [wave] * n
    if ivar is None:
        ivar_rows = [None] * n
    else:
        ivar = np.asarray(ivar, dtype=float)
        ivar_rows = ivar if ivar.ndim == 2 else [ivar] * n

    if bands is None:
        bands = ARM_EVAL_RANGES
    if line_list is None:
        line_list = default_line_list(fwhm)

    args = [(wave_rows[i], flux[i], model[i], ivar_rows[i], fwhm, line_list, bands, mask)
           for i in range(n)]

    if nproc and nproc > 1:
        with multiprocessing.Pool(nproc) as pool:
            results = pool.starmap(_analyze_one_row, args)
    else:
        results = [_analyze_one_row(*a) for a in args]

    summaries, line_tables = zip(*results)
    summaries = list(summaries)
    for i, s in enumerate(summaries):
        s['ROW'] = i

    line_tables_tagged = []
    for i, lt in enumerate(line_tables):
        if len(lt):
            lt = lt.copy()
            lt['ROW'] = i
            line_tables_tagged.append(lt)

    summary_table = Table(rows=summaries)
    line_table = vstack(line_tables_tagged) if line_tables_tagged else Table()
    return summary_table, line_table


# ──────────────────────────────────────────────────────────────
# Summary plots (written to a subdirectory, plots_sky_resid by default).
# ──────────────────────────────────────────────────────────────

_RESIDUAL_SIGN_NOTE = 'Observations - Sky'


def plot_frac_summary(tab, outdir='plots_sky_resid', bands=None, outname='frac_summary.png'):
    '''
    One-page 2x3 grid of the FRAC_1E15/1E14 REVERSE cumulative
    distributions: rows = continuum/lines, columns = bands.  Each panel
    overlays both thresholds as reverse-cumulative, density, step
    histograms (100 bins, range 0-1), e.g.::

        plt.hist((tab['LINE_FRAC_1E15_R'], tab['LINE_FRAC_1E14_R']), 100,
                 range=(0, 1), cumulative=-1, histtype='step', density=True)

    ``cumulative=-1`` (P(value >= x), not the forward P(value <= x)) is
    deliberate: since the 1e-14 threshold is looser, its per-row FRAC
    value is >= the 1e-15 one for every row, and a REVERSE cumulative is
    guaranteed to place the 1e-14 curve at or above the 1e-15 curve at
    every x for any dataset -- a forward cumulative doesn't have that
    guarantee (see the inline comment in the implementation).

    Parameters
    ----------
    tab : astropy.table.Table or str/Path
        A summary_table from analyze_sky_residuals (or a FITS file of
        one).
    outdir : str or Path
        Directory the PNG is written into (created if missing).
    bands : list of (label, lo, hi) or None
        Defaults to ARM_EVAL_RANGES; only the labels are used, to line
        up with the CONT_FRAC_*/LINE_FRAC_* column suffixes.
    outname : str
        Output filename within outdir.

    Returns
    -------
    Path to the written PNG.
    '''
    if isinstance(tab, (str, Path)):
        tab = Table.read(tab)
    if bands is None:
        bands = ARM_EVAL_RANGES
    labels = [b[0] for b in bands]

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, len(labels), figsize=(4 * len(labels), 8), squeeze=False)

    for i, (prefix, row_title) in enumerate((('CONT_FRAC', 'Continuum'), ('LINE_FRAC', 'Lines'))):
        for j, label in enumerate(labels):
            ax = axes[i][j]
            col15, col14 = '%s_1E15_%s' % (prefix, label), '%s_1E14_%s' % (prefix, label)
            v15 = (np.asarray(tab[col15], dtype=float) if col15 in tab.colnames else np.array([]))
            v14 = (np.asarray(tab[col14], dtype=float) if col14 in tab.colnames else np.array([]))
            v15, v14 = v15[np.isfinite(v15)], v14[np.isfinite(v14)]
            if len(v15) or len(v14):
                # cumulative=-1 (reverse: P(value >= x), not P(value <= x))
                # is deliberate, not cosmetic: since F14 (looser threshold)
                # is >= F15 for every row by construction, a REVERSE
                # cumulative is guaranteed to place the 1e-14 curve at or
                # above the 1e-15 curve at every x, for any dataset -- i.e.
                # "more good pixels at the looser threshold" is a visual
                # invariant, not something that depends on how the two
                # distributions happen to be shaped.  A forward cumulative
                # (P(value <= x)) doesn't have that guarantee: whichever
                # distribution's mass sits further right reads as the
                # LOWER curve for most of the range, which for F14
                # (concentrated near 1) is backwards from the intuitive
                # reading (ksl, 260803).
                ax.hist((v15, v14), bins=100, range=(0, 1), cumulative=-1,
                       histtype='step', density=True, label=['1e-15', '1e-14'])
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1.02)
            ax.set_box_aspect(1)
            ax.set_xlabel('Fraction of pixels below threshold (%s, %s)' % (row_title, label))
            if j == 0:
                ax.set_ylabel('Fraction of rows >= this value')
            ax.legend(fontsize=8, loc='lower right')

    fig.suptitle('Quality Fraction Summary (%s)' % _RESIDUAL_SIGN_NOTE)
    fig.tight_layout()
    outpath = outdir / outname
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print('Wrote %s' % outpath)
    return outpath


def _percentile_xrange(vals, lo=1.0, hi=99.0, pad_frac=0.05):
    '''
    Robust histogram x-range from percentiles rather than min/max, so a
    handful of outlier exposures/fibers don't compress the bulk of a
    real-flux-unit distribution into a sliver of the panel.
    '''
    if len(vals) == 0:
        return (0.0, 1.0)
    p_lo, p_hi = np.percentile(vals, [lo, hi])
    if p_lo == p_hi:
        pad = 1.0 if p_lo == 0 else abs(p_lo) * pad_frac
        return (p_lo - pad, p_hi + pad)
    pad = (p_hi - p_lo) * pad_frac
    return (p_lo - pad, p_hi + pad)



def _plot_metric_grid(tab, prefix, metrics, bands, outdir, outname, row_labels=None,
                      suptitle=None, ref_lines=None):
    '''
    Shared grid-of-histograms machinery for plot_continuum_summary and
    plot_lines_summary: rows = metrics (columns named
    <prefix>_<metric>_<band>), columns = bands.  Each panel is a
    histogram; the x-range is shared across all three bands in a row
    (set from the 1st-99th percentile of the POOLED values across
    bands, not min/max), so the relative width of the distribution in
    each arm is directly visually comparable rather than each panel
    auto-zooming to its own scale -- a few outlier exposures/fibers
    don't compress the bulk of the pooled distribution (points outside
    the shared range are simply excluded from that panel's histogram).
    No panel titles -- the metric's human-readable name (from
    ``row_labels``) plus the band are combined into the x-axis label on
    every panel, so each one reads standalone even if cropped out of
    the page.  ``ref_lines``, if given, is a dict {metric: x_value}
    drawing a dashed vertical reference line in that row (e.g. AMPRATIO
    at 1.0, the "model matches data exactly" point).  A metric with no
    finite values in ANY band is dropped from the grid entirely rather
    than left as a blank row.
    '''
    if isinstance(tab, (str, Path)):
        tab = Table.read(tab)
    if bands is None:
        bands = ARM_EVAL_RANGES
    labels = [b[0] for b in bands]
    row_labels = row_labels or {}
    ref_lines = ref_lines or {}

    all_band_vals = {}
    for metric in metrics:
        band_vals = {}
        for label in labels:
            col = '%s_%s_%s' % (prefix, metric, label)
            vals = np.asarray(tab[col], dtype=float) if col in tab.colnames else np.array([])
            band_vals[label] = vals[np.isfinite(vals)]
        all_band_vals[metric] = band_vals
    metrics = [m for m in metrics if any(len(v) for v in all_band_vals[m].values())]
    if not metrics:
        print('Nothing to plot for %s (%s prefix): every metric was empty.' % (outname, prefix))
        return None

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    nrows = len(metrics)
    fig, axes = plt.subplots(nrows, len(labels), figsize=(4 * len(labels), 4 * nrows),
                             squeeze=False)

    for i, metric in enumerate(metrics):
        band_vals = all_band_vals[metric]
        pooled = (np.concatenate(list(band_vals.values())) if any(len(v) for v in band_vals.values())
                 else np.array([]))
        xlo, xhi = _percentile_xrange(pooled)

        for j, label in enumerate(labels):
            ax = axes[i][j]
            vals = band_vals[label]
            if len(vals):
                ax.hist(vals, bins=40, range=(xlo, xhi))
            if metric in ref_lines:
                ax.axvline(ref_lines[metric], color='k', linestyle='--', linewidth=1)
            ax.set_xlim(xlo, xhi)
            ax.set_box_aspect(1)
            ax.set_xlabel('%s (%s)' % (row_labels.get(metric, metric.title()), label))
            if j == 0:
                ax.set_ylabel('N')

    if suptitle:
        fig.suptitle('%s (%s)' % (suptitle, _RESIDUAL_SIGN_NOTE))
    fig.tight_layout()
    outpath = outdir / outname
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print('Wrote %s' % outpath)
    return outpath


def plot_continuum_summary(tab, outdir='plots_sky_resid', bands=None,
                           outname='continuum_summary.png',
                           metrics=('OFFSET', 'NMAD', 'ALPHA', 'FIT_QUALITY')):
    '''
    One-page grid of continuum-quality histograms: rows = metrics
    (default OFFSET/NMAD/ALPHA/FIT_QUALITY -- bias, scatter,
    dimensionless power-law-index mismatch, NMAD/NOISE_PROXY -- see
    analyze_sky_residual's Notes for why this replaces a formal ivar
    chi-square), columns = bands.  Each panel is a histogram of
    CONT_<METRIC>_<band> across rows of ``tab``.  A metric with no
    finite values in ANY band is dropped from the grid entirely rather
    than plotted as an empty row -- see _plot_metric_grid.

    Parameters
    ----------
    tab : astropy.table.Table or str/Path
        A summary_table from analyze_sky_residuals (or a FITS file of
        one).
    outdir : str or Path
        Directory the PNG is written into (created if missing).
    bands : list of (label, lo, hi) or None
        Defaults to ARM_EVAL_RANGES.
    outname : str
    metrics : tuple of str
        Suffixes of CONT_<METRIC>_<band> columns to plot, one row each.

    Returns
    -------
    Path to the written PNG, or None if every metric was empty.
    '''
    row_labels = {'OFFSET': 'Continuum Offset', 'NMAD': 'Continuum Scatter (NMAD)',
                 'ALPHA': r'Continuum Power-Law Index Mismatch ($\alpha$)',
                 'FIT_QUALITY': 'Continuum Fit Quality (NMAD/NOISE_PROXY)'}
    return _plot_metric_grid(tab, 'CONT', metrics, bands, outdir, outname,
                             row_labels=row_labels, suptitle='Continuum Residual Summary',
                             ref_lines={'OFFSET': 0.0, 'NMAD': 0.0, 'ALPHA': 0.0,
                                       'FIT_QUALITY': 1.0})


def plot_lines_summary(tab, outdir='plots_sky_resid', bands=None,
                       outname='lines_summary.png',
                       metrics=('AMPRATIO_MED', 'DLAM_MED', 'DSIG_MED')):
    '''
    One-page grid of line shape-fit histograms: rows = LINE_AMPRATIO_MED
    (predicted/measured amplitude ratio -- see analyze_sky_residual's
    AMP_RATIO for why a ratio is used instead of the flux-unit DELTA_A),
    LINE_DLAM_MED (registration bias), LINE_DSIG_MED (LSF-width bias --
    the number that answers "is a small LSF mismatch contributing"),
    columns = bands.  Each panel is a histogram of LINE_<METRIC>_<band>
    across rows of ``tab``, with a dashed reference line at 0.0 for
    DLAM_MED/DSIG_MED and 1.0 for AMPRATIO_MED (model exactly matches
    data).  See _plot_metric_grid for the shared-x-range convention.

    Parameters
    ----------
    tab : astropy.table.Table or str/Path
        A summary_table from analyze_sky_residuals (or a FITS file of
        one).
    outdir : str or Path
        Directory the PNG is written into (created if missing).
    bands : list of (label, lo, hi) or None
        Defaults to ARM_EVAL_RANGES.
    outname : str
    metrics : tuple of str
        Suffixes of LINE_<METRIC>_<band> columns to plot, one row each.

    Returns
    -------
    Path to the written PNG.
    '''
    row_labels = {'AMPRATIO_MED': 'Amplitude Ratio (Sky/Obs)',
                 'DLAM_MED': r'Wavelength Registration Bias ($\AA$)',
                 'DSIG_MED': r'LSF-Width Bias ($\AA$)'}
    return _plot_metric_grid(tab, 'LINE', metrics, bands, outdir, outname,
                             row_labels=row_labels, suptitle='Line Residual Summary',
                             ref_lines={'AMPRATIO_MED': 1.0, 'DLAM_MED': 0.0, 'DSIG_MED': 0.0})


# ──────────────────────────────────────────────────────────────
# Minimal CLI: WAVE/FLUX/SKY[/IVAR] convention (SkySub_eval.py's file
# format), primarily as a smoke test / convenience wrapper.
# ──────────────────────────────────────────────────────────────

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description=_usage_from_doc(__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename')
    parser.add_argument('-mask', default=None)
    parser.add_argument('-nproc', type=int, default=1)
    parser.add_argument('-out', default=None)
    parser.add_argument('-plotdir', default='plots_sky_resid')
    args = parser.parse_args()

    hdul = fits.open(args.filename)
    wave = hdul['WAVE'].data.astype(float)
    model = hdul['SKY'].data.astype(float)
    # FLUX in this file convention is already sky-subtracted (see module
    # docstring) -- reconstruct the observed spectrum as FLUX+SKY so the
    # internal flux-model subtraction reduces back to FLUX exactly,
    # rather than double-subtracting SKY.
    flux = hdul['FLUX'].data.astype(float) + model
    ivar = hdul['IVAR'].data.astype(float) if 'IVAR' in hdul else None
    hdul.close()

    summary_table, line_table = analyze_sky_residuals(
        wave, flux, model, ivar=ivar, mask=args.mask, nproc=args.nproc)

    outroot = args.out or Path(args.filename).stem
    summary_table.write('%s_summary.fits' % outroot, overwrite=True)
    line_table.write('%s_lines.fits' % outroot, overwrite=True)
    print(summary_table)
    print('Wrote %s_summary.fits and %s_lines.fits' % (outroot, outroot))

    plot_frac_summary(summary_table, outdir=args.plotdir,
                      outname='%s_frac_summary.png' % outroot)
    plot_continuum_summary(summary_table, outdir=args.plotdir,
                           outname='%s_continuum_summary.png' % outroot)
    plot_lines_summary(summary_table, outdir=args.plotdir,
                       outname='%s_lines_summary.png' % outroot)


if __name__ == '__main__':
    main()
