#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Measure how much flux leaks into each of sky_gaussfit.py's nebular
    emission-line windows in a DecomposeCleanSky.py output's residual
    (observed - clean-sky model) -- a direct, per-line, per-extension test
    of whether the lvmdrp assumption that SKY_EAST/SKY_WEST contain no
    nebular-line flux actually holds.

Command line usage (if any):

    usage: sky_nebular_leak_eval.py [-h] [-ext LIST] [-v VEL] [-lmc] [-smc]
                                    [-sigma S] [-thresh T] [-out ROOT]
                                    [-template_ext EXT] [-template_band LO,HI]
                                    fits_file [fits_file ...]

    where

    fits_file    one or more DecomposeCleanSky.py output FITS files
                (WAVE, <EXT>, <EXT>_BESTFIT[, <EXT>_RESID] extensions).

    -ext LIST    comma-separated extensions to evaluate (default:
                SKY_EAST,SKY_WEST -- must match what DecomposeCleanSky.py
                was run with).

    -v VEL       nebular systemic velocity (km/s) used to Doppler-shift
                the line centers (default: 0). Should match the -v/-lmc/
                -smc DecomposeCleanSky.py was run with, or the fit window
                will be centered on the wrong wavelength.

    -lmc / -smc  shortcuts for the LMC (~262 km/s) / SMC (~146 km/s)
                systemic velocity, same convention as sky_gaussfit.py.
                -v overrides these if both are given.

    -sigma S     initial Gaussian sigma guess in Angstrom (default 1.0).

    -thresh T    ``amplitude / amplitude_error`` above which a line is
                flagged as a candidate leak in the printed summary
                (default 3.0). Does not affect what's written to the
                output table -- every line is written regardless.

    -out ROOT    output table filename root (default: 'nebular_leak').

    -template_ext EXT
                extension (e.g. SKY_EAST) whose RESID is used as a
                shared-systematic template to correct every OTHER
                requested extension's residual, inside -template_band,
                before fitting lines there (see fit_template_scale).
                Confirmed at r~0.95-0.98 between FLUX/SKY_EAST/SKY_WEST
                *within one exposure*, but only r~0.5-0.6 *across
                different exposures* -- always give a fits_file whose
                own extensions include EXT; a template does not transfer
                across files/exposures. Default: no correction applied.

    -template_band LO,HI
                wavelength range (Angstrom) the correction is trusted in
                and applied to (default: 9000,9600 -- the range this was
                empirically validated in). Lines outside this range are
                always left uncorrected even when -template_ext is given.

Description:

    For each requested extension in each input file, fits a Gaussian plus
    constant background directly to the RESID array (observed - model) in
    a window around each of sky_gaussfit.resolve_nebular_lines()'
    Doppler-shifted rest wavelengths -- NEBULAR_LINES with any line
    dropped whose window still overlaps a SKY_LINES window at the given
    velocity (oi_a/oi_b at low velocity: predominantly sky airglow, not
    nebular, since they're the same [OI] 6300/6364 transition as
    sky6300/sky6363 -- see resolve_nebular_lines' docstring).

    This deliberately does NOT reuse sky_residual_eval.py's own line-shape
    machinery (analyze_sky_residual/_fit_line_to_model/_line_shape_fit):
    that machinery first fits a Gaussian to the MODEL to anchor a line's
    expected center/width, then measures the residual's mismatch against
    that model-derived shape -- appropriate for airglow lines the model is
    expected to already contain. A nebular line should have ~no flux in a
    DecomposeCleanSky.py model by construction (its window was excluded
    from the fit), so there is no model peak to anchor on; fitting the
    residual directly with the line's own catalog wavelength as the prior
    is the correct approach here. The plain _gaussian(wave, amp, center,
    sigma, bkg) functional form is reused from sky_residual_eval.py
    (identical shape, no reason to redefine it); the fit path is new.

    A significant nonzero fitted amplitude (LEAK_FLUX / LEAK_FLUX_ERR
    beyond -thresh) at a nebular line's position, in an extension the DRP
    treats as pure sky (SKY_EAST/SKY_WEST), is direct evidence of nebular
    contamination the mask-and-wrap decomposition -- and by extension the
    DRP's own sky subtraction, which makes the same no-nebular-lines
    assumption -- did not remove.

    This script only ever needs (wave, flux_observed, flux_model) per
    extension (RESID is recomputed from FLUX/BESTFIT if not present in
    the input file), so it works unchanged against a future native-
    nebular-family SkyDecomp output too, not just DecomposeCleanSky.py's
    mask-and-wrap models.

    Template correction (-template_ext, off by default): in real data
    (LVM field, exposure with strong nebular emission), the 9000-9600 A
    band residual is dominated by a shared systematic (most likely an
    OH-line/LSF template mismatch in SkyDecomp's Z channel) rather than
    photon noise -- confirmed by cross-correlating FLUX/SKY_EAST/
    SKY_WEST's residuals against each other within one exposure
    (Pearson r = 0.947-0.984) versus the SAME extension's residual
    across two different exposures (r = 0.541-0.629 only). Because the
    systematic is common across fibers of one exposure but real nebular
    flux is not, an extension confirmed to carry little real signal
    (e.g. one that fails a known-fixed-ratio consistency check, such as
    [SIII] 9531/9069 = 2.44 from the Einstein A coefficients) can serve
    as that exposure's own correction template for the others via
    fit_template_scale's least-squares scale fit, restricted to
    -template_band and excluding all NEBULAR_LINES windows from the
    scale fit itself so real signal there can't bias it.

Primary routines:

    fit_leak_line        Gaussian-on-residual fit for one line window.
    fit_template_scale   least-squares scale fit of one residual to another.
    evaluate_file        loop over extensions/lines for one input file.

History::

    260909  ksl  Coding begun.

'''

import argparse
import re
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from scipy.optimize import curve_fit

from sky_gaussfit import resolve_nebular_lines
from sky_residual_eval import _gaussian


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before a line consisting of "History:" (or
    "History::"/"Version History"), so -h stays short.
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


LMC_VEL = 262.
SMC_VEL = 146.

_LEAK_COLUMNS = ['FILE', 'EXT', 'LINE_NAME', 'WAVE0', 'LEAK_FLUX',
                'LEAK_FLUX_ERR', 'SNR', 'LAM0_FIT', 'SIGMA0_FIT', 'N_PIX', 'FLAG',
                'LEAK_FLUX_RAW', 'LEAK_FLUX_ERR_RAW', 'SNR_RAW',
                'TEMPLATE_EXT', 'TEMPLATE_K', 'TEMPLATE_R']


def resolve_velocity(vel=None, lmc=False, smc=False, drp_redshift=None):
    '''
    Resolve the nebular systemic velocity from -v/-lmc/-smc, falling back
    to a per-row DRP_ALL['Redshift'] value when none of those were given.

    Same precedence as sky_gaussfit.py's steer() / DecomposeCleanSky.py's
    resolve_velocity: an explicit -v wins, then -lmc, then -smc. Only
    below all three does drp_redshift apply -- DRP_ALL['Redshift'] is
    SummarizeCframe.py's own RA/Dec-based classification (LMC/SMC/
    Plane/HighLat -> 262./146./0./0. km/s, see that file's Survey/
    Redshift construction), so when the caller hasn't overridden it,
    this is the physically correct per-row default -- NOT 0. Falls back
    to 0 only if drp_redshift itself is also None/non-finite (column
    missing, or a non-nebular field).

    Without this, every nebular-line fit window silently defaulted to
    rest wavelength regardless of target -- fine for Galactic ('Plane'/
    'HighLat') pointings, but for a real LMC/SMC exposure ~262/146 km/s
    (~4/2 A here) is comparable to or larger than the fit window itself,
    so the line being measured is mostly or entirely outside it. Found
    260910 on a real LMC exposure (12226): every doublet SNR-gated out
    as apparent non-detections, which were actually mislocated fits, not
    a genuinely quiescent field (see this session's own investigation).
    '''
    if vel is not None:
        return float(vel)
    if lmc:
        return LMC_VEL
    if smc:
        return SMC_VEL
    if drp_redshift is not None and np.isfinite(drp_redshift):
        return float(drp_redshift)
    return 0.0


def fit_leak_line(wave, resid, wave0, wlo, whi, sigma_guess=1.0):
    '''
    Fit a Gaussian + constant background directly to resid in the window
    [wlo, whi] (no model-derived shape is assumed -- see module
    docstring).

    Parameters
    ----------
    wave, resid : ndarray
        Wavelength grid and residual (observed - model) spectrum.
    wave0 : float
        Doppler-shifted nebular line center (initial center guess/bound
        anchor); need not be the window's midpoint.
    wlo, whi : float
        Fit window bounds (Angstrom). Pass the line's own (Doppler-
        shifted) wmin/wmax directly -- do NOT symmetrize these into a
        single half-width first (that inflates the shorter side to match
        the longer one, which for asymmetric windows like siii_a/siii_b
        sweeps unrelated nearby features into the local-background
        estimate; found via PlotSkySubNebEval.py's per-line panels on a
        real exposure, 260910).
    sigma_guess : float
        Initial Gaussian sigma guess (Angstrom).

    Returns
    -------
    dict or None
        Keys: amp, amp_err, lam0, sigma0, npix, bkg. None if the window
        has too few finite pixels or the fit fails to converge.
    '''
    sel = (wave >= wlo) & (wave <= whi)
    w, r = wave[sel], resid[sel]
    finite = np.isfinite(w) & np.isfinite(r)
    w, r = w[finite], r[finite]
    if len(w) < 5:
        return None

    bkg_guess = float(np.nanmedian(r))
    amp_guess = float(np.nanmax(r) - bkg_guess)
    # Real LVM flux is ~1e-14 to 1e-16 -- rescale to O(1) before the fit,
    # same reasoning as sky_residual_eval._fit_line_to_model.
    scale = max(abs(amp_guess), float(np.nanstd(r)), 1e-30)
    wave0 = min(max(wave0, wlo), whi)  # p0 must lie within its own bounds
    p0 = [amp_guess / scale, wave0, sigma_guess, bkg_guess / scale]
    bounds = ([-np.inf, wlo, 0.1, -np.inf],
             [np.inf, whi, 5.0, np.inf])
    try:
        popt, pcov = curve_fit(_gaussian, w, r / scale, p0=p0, bounds=bounds, maxfev=4000)
    except Exception:
        return None
    if not np.all(np.isfinite(popt)):
        return None
    perr = np.sqrt(np.clip(np.diag(pcov), 0, None))
    amp, lam0, sigma0, bkg = popt
    return dict(amp=float(amp * scale), amp_err=float(perr[0] * scale),
               lam0=float(lam0), sigma0=float(sigma0), npix=int(len(w)),
               bkg=float(bkg * scale))


def line_exclusion_mask(wave, resolved_lines, vel=0.0):
    '''
    Boolean mask, True inside any resolved_lines entry's Doppler-shifted
    window. Used to keep real nebular-line flux from biasing a template
    scale fit -- the inverse role of DecomposeCleanSky.nebular_mask
    (which is True *outside* these windows), reimplemented locally here
    rather than imported, since this module has no lvmsky/py_dev
    dependency and the computation is one line.

    Parameters
    ----------
    wave : ndarray
    resolved_lines : list of (name, center, wmin, wmax)
        Typically resolve_nebular_lines(vel)'s first return value.
    vel : float

    Returns
    -------
    ndarray of bool
    '''
    zz = 1.0 + vel / 3e5
    excl = np.zeros(wave.shape, dtype=bool)
    for _name, _center, wmin, wmax in resolved_lines:
        excl |= (wave >= zz * wmin) & (wave <= zz * wmax)
    return excl


def fit_template_scale(template, target, exclude=None):
    '''
    Least-squares scale k minimizing ||target - k*template||^2, plus the
    Pearson correlation over the same pixels as a fit-quality diagnostic.

    Used to correct a shared systematic (e.g. an OH-line/LSF template
    mismatch common to every fiber of one exposure -- confirmed present
    at r~0.95-0.98 between FLUX/SKY_EAST/SKY_WEST within a single
    exposure, but NOT stable exposure-to-exposure at r~0.5-0.6, so a
    template must come from the same exposure/file as the target, never
    a different one) out of a target extension's residual before fitting
    a line in it.

    Parameters
    ----------
    template, target : ndarray
        Same-length residual arrays (e.g. one extension's RESID used as
        the template for another extension's RESID in the same file).
    exclude : ndarray of bool, optional
        True at pixels to leave out of the fit (e.g. line_exclusion_mask
        for the lines being measured, so real signal there can't bias
        the scale estimate).

    Returns
    -------
    k, r : float
        Fitted scale and Pearson correlation, or (nan, nan) if fewer
        than 10 usable pixels remain.
    '''
    good = np.isfinite(template) & np.isfinite(target)
    if exclude is not None:
        good &= ~exclude
    t, x = template[good], target[good]
    if len(t) < 10:
        return np.nan, np.nan
    denom = np.sum(t * t)
    if denom <= 0:
        return np.nan, np.nan
    k = float(np.sum(t * x) / denom)
    r = float(np.corrcoef(t, x)[0, 1])
    return k, r


def evaluate_file(fits_file, exts, vel=0.0, sigma_guess=1.0, thresh=3.0,
                  template_ext=None, template_band=(9000., 9600.)):
    '''
    Fit every NEBULAR_LINES window against each requested extension's
    residual in one DecomposeCleanSky.py output file.

    Parameters
    ----------
    fits_file : str or Path
    exts : sequence of str
    vel : float
        Nebular systemic velocity (km/s); must match what the input file
        was produced with, or the fit window is centered wrong.
    sigma_guess : float
    thresh : float
        ``SNR`` above which a row's FLAG is set for the printed summary.
    template_ext : str, optional
        Extension whose RESID is used as a template (fit_template_scale)
        to remove a shared systematic from every OTHER requested
        extension's residual before fitting lines inside template_band.
        Confirmed valid only within one exposure/file (FLUX/SKY_EAST/
        SKY_WEST residuals correlate at r~0.95-0.98 within one exposure,
        but only r~0.5-0.6 across different exposures) -- always use a
        template_ext from the SAME fits_file, never a different one.
        template_ext's own rows are still fit normally (uncorrected --
        it has no other extension to draw a template from).
    template_band : (float, float)
        Wavelength range (Angstrom) the correction is trusted in and
        applied to; lines outside this range are always left uncorrected
        even when template_ext is given, since the shared-systematic
        behavior has only been confirmed empirically in this range.

    Returns
    -------
    astropy.table.Table
        One row per (extension, line); see _LEAK_COLUMNS. LEAK_FLUX/
        LEAK_FLUX_ERR/SNR are the template-corrected values when a
        correction was applied (line in-band, template_ext given,
        ext != template_ext, template fit succeeded); otherwise they
        equal the _RAW columns. TEMPLATE_K/TEMPLATE_R are the fitted
        scale/correlation for that (file, ext) pair (NaN if no
        correction was attempted for that ext).
    '''
    zz = 1.0 + vel / 3e5
    resolved, dropped = resolve_nebular_lines(vel)
    if dropped:
        print(f"  {fits_file}: skipping (sky, not nebular, at vel={vel:.1f}): "
              f"{', '.join(dropped)}")

    def load_resid(hdul, ext):
        if f'{ext}_RESID' in hdul:
            return np.asarray(hdul[f'{ext}_RESID'].data, dtype=float)
        if ext in hdul and f'{ext}_BESTFIT' in hdul:
            return (np.asarray(hdul[ext].data, dtype=float)
                   - np.asarray(hdul[f'{ext}_BESTFIT'].data, dtype=float))
        return None

    rows = []
    with fits.open(fits_file) as hdul:
        wave = np.asarray(hdul['WAVE'].data, dtype=float)

        template_resid = None
        if template_ext is not None:
            template_resid = load_resid(hdul, template_ext)
            if template_resid is None:
                print(f"  {fits_file}: -template_ext {template_ext} not found, "
                      f"no correction applied")

        band_sel = (wave >= template_band[0]) & (wave <= template_band[1])
        excl = line_exclusion_mask(wave, resolved, vel)

        for ext in exts:
            resid = load_resid(hdul, ext)
            if resid is None:
                print(f"  {fits_file}: no {ext}/{ext}_BESTFIT extensions, skipping")
                continue

            k, r = np.nan, np.nan
            corrected = resid
            apply_correction = template_resid is not None and ext != template_ext
            if apply_correction:
                k, r = fit_template_scale(template_resid[band_sel], resid[band_sel],
                                          exclude=excl[band_sel])
                if np.isfinite(k):
                    corrected = resid.copy()
                    corrected[band_sel] = resid[band_sel] - k * template_resid[band_sel]
                else:
                    apply_correction = False

            for name, center, wmin, wmax in resolved:
                wave0 = zz * center
                wlo, whi = zz * wmin, zz * wmax
                fit_raw = fit_leak_line(wave, resid, wave0, wlo, whi, sigma_guess=sigma_guess)
                in_band = template_band[0] <= wave0 <= template_band[1]
                use_corrected = apply_correction and in_band
                fit_use = (fit_leak_line(wave, corrected, wave0, wlo, whi, sigma_guess=sigma_guess)
                          if use_corrected else fit_raw)

                if fit_raw is None:
                    amp_raw = amp_err_raw = snr_raw = np.nan
                else:
                    amp_raw, amp_err_raw = fit_raw['amp'], fit_raw['amp_err']
                    snr_raw = amp_raw / amp_err_raw if amp_err_raw > 0 else np.nan

                if fit_use is None:
                    rows.append([str(fits_file), ext, name, wave0, np.nan, np.nan,
                                np.nan, np.nan, np.nan, 0, 'NOFIT',
                                amp_raw, amp_err_raw, snr_raw,
                                template_ext or '', k, r])
                    continue

                snr = fit_use['amp'] / fit_use['amp_err'] if fit_use['amp_err'] > 0 else np.nan
                flag = 'LEAK' if np.isfinite(snr) and abs(snr) >= thresh else ''
                rows.append([str(fits_file), ext, name, wave0, fit_use['amp'], fit_use['amp_err'],
                            snr, fit_use['lam0'], fit_use['sigma0'], fit_use['npix'], flag,
                            amp_raw, amp_err_raw, snr_raw,
                            template_ext or '', k, r])

    return Table(rows=rows, names=_LEAK_COLUMNS)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=_usage_from_doc(__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fits_file', nargs='+',
                        help='DecomposeCleanSky.py output FITS file(s)')
    parser.add_argument('-ext', default='SKY_EAST,SKY_WEST',
                        help='Comma-separated extensions to evaluate')
    parser.add_argument('-v', dest='vel', type=float, default=None,
                        help='Nebular systemic velocity (km/s)')
    parser.add_argument('-lmc', action='store_true',
                        help=f'Use the LMC velocity (~{LMC_VEL:.0f} km/s)')
    parser.add_argument('-smc', action='store_true',
                        help=f'Use the SMC velocity (~{SMC_VEL:.0f} km/s)')
    parser.add_argument('-sigma', dest='sigma_guess', type=float, default=1.0,
                        help='Initial Gaussian sigma guess (Angstrom)')
    parser.add_argument('-thresh', type=float, default=3.0,
                        help='|SNR| threshold for the printed leak summary')
    parser.add_argument('-out', default='nebular_leak',
                        help='Output table filename root')
    parser.add_argument('-template_ext', default=None,
                        help='Extension (e.g. SKY_EAST) whose RESID is used as a '
                        'shared-systematic template to correct every other '
                        'requested extension inside -template_band, before '
                        'fitting lines there. Must come from the same file as '
                        'the extension it corrects (confirmed valid within one '
                        'exposure, not across exposures) -- do not mix files.')
    parser.add_argument('-template_band', default='9000,9600',
                        help='LO,HI (Angstrom) the template correction is '
                        'trusted in and applied to; lines outside this range '
                        'are always left uncorrected')
    args = parser.parse_args()

    vel = resolve_velocity(args.vel, args.lmc, args.smc)
    exts = [e.strip() for e in args.ext.split(',') if e.strip()]
    template_band = tuple(float(x) for x in args.template_band.split(','))
    print(f"Evaluating {len(args.fits_file)} file(s), ext={exts}, vel={vel:.1f} km/s"
          + (f", template_ext={args.template_ext}, template_band={template_band}"
             if args.template_ext else ""))

    tables = [evaluate_file(f, exts, vel=vel, sigma_guess=args.sigma_guess,
                            thresh=args.thresh, template_ext=args.template_ext,
                            template_band=template_band) for f in args.fits_file]
    leak_table = vstack(tables) if tables else Table(names=_LEAK_COLUMNS)

    outpath = f'{args.out}_lines.fits'
    leak_table.write(outpath, overwrite=True)
    print(f"Wrote {outpath} ({len(leak_table)} rows)")

    flagged = leak_table[leak_table['FLAG'] == 'LEAK']
    if len(flagged):
        print(f"\n{len(flagged)} candidate leak(s) (|SNR| >= {args.thresh}):")
        flagged['FILE', 'EXT', 'LINE_NAME', 'WAVE0', 'LEAK_FLUX', 'SNR'].pprint(
            max_lines=-1, max_width=-1)
    else:
        print(f"\nNo candidate leaks at |SNR| >= {args.thresh}")


if __name__ == '__main__':
    main()
