#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Science-specific evaluator for the SkySub* method family: fits
    sky_gaussfit.NEBULAR_LINES directly on each method's sky-subtracted
    FLUX (not a generic sky-diagnostic window), checks fixed-ratio
    doublets against their known atomic-physics value, and measures
    flux/ratio scatter across repeated observations of the same tile --
    the questions SkySub_eval.py's generic sky-residual windows don't
    address, since none of them sit at these wavelengths.

Command line usage (if any)::

        usage: SkySubNebEval.py [-h] [-v VEL] [-lmc] [-smc] [-sigma S]
                                [-mjd_close DAYS] [-snr_min S] [-out ROOT]
                                fits_file [fits_file ...]

        where

        fits_file     one or more SkySub*.py output FITS files (WAVE, FLUX,
                     SKY, DRP_ALL -- SkySubDrp/Orig/Dev1/Dev2/Dev3.py all
                     write this layout). Each file's primary header ROUTINE/
                     METHOD (Title/METHOD keywords) labels its rows in the
                     output.

        -v VEL        nebular systemic velocity (km/s) used to Doppler-shift
                     the line centers, applied to EVERY row regardless of
                     target. Default: none -- instead, each row's velocity
                     is looked up from that row's own DRP_ALL['Redshift']
                     (SummarizeCframe.py's RA/Dec-based LMC/SMC/Plane/HighLat
                     classification -> 262/146/0/0 km/s). Only pass -v/-lmc/
                     -smc to force one velocity across an entire file/run --
                     correct only if every row is known to share one target.

        -lmc / -smc   shortcuts for the LMC (~262 km/s) / SMC (~146 km/s)
                     systemic velocity, same convention as sky_gaussfit.py --
                     same override semantics as -v above.

        -sigma S      initial Gaussian sigma guess in Angstrom (default 1.0).

        -mjd_close DAYS
                     a tileid group's repeat exposures are tagged "closely
                     spaced" when their MJD span is below this (default 7.0
                     days) -- see Description.

        -snr_min S    minimum per-line SNR (both lines of a doublet) required
                     to report that row's ratio (default 5.0) -- most fibers
                     in an arbitrary sample are not pointed at a real
                     emission-line target, so without this gate the doublet-
                     ratio statistics are dominated by noise-fit garbage, not
                     a meaningful measurement (verified: unfiltered median
                     SIII ratio across a real 140-row sample was ~2.98-3.00,
                     nowhere near the true 2.44; restricting to SNR>5 on both
                     lines dropped it to 2.41).

        -out ROOT     output filename root (default: 'nebeval').

Description::

        Two-stage evaluation, one row-level table and one group-level table:

        1. fit_nebular_row: every sky_gaussfit.NEBULAR_LINES entry (all of
           them -- unlike DecomposeCleanSky.py/sky_nebular_leak_eval.py,
           there is no reason to drop oi_a/oi_b here: this fits the
           *sky-subtracted* FLUX, not a pre-subtraction residual being
           screened for sky-line contamination, so oi_a/oi_b's coincidence
           with sky6300/sky6363 is itself part of what's being tested --
           residual sky-subtraction error at that wavelength would bias
           oi_a/oi_b differently from each other and break their fixed
           ratio, exactly the kind of thing this check is for) is fit
           directly on FLUX with sky_nebular_leak_eval.fit_leak_line
           (Gaussian + local constant background, no model-anchoring --
           identical reasoning to that module, just applied to a real
           spectrum instead of a residual). DOUBLETS below then computes
           each fixed-ratio pair's measured ratio, its propagated error,
           and its deviation from the true ratio in units of that error.

        2. repeat_scatter: DRP_ALL['tileid'] groups rows into repeat
           observations of the same pointing, EXCLUDING tileid 11111 (a
           grab-bag placeholder used across many unrelated targets in real
           data, confirmed by inspection -- not a genuine repeat pointing;
           see this session's own investigation). Each real tileid's RA/Dec
           is checked to be tightly clustered before being trusted as a
           repeat group. For every group with >=2 rows, the flux/ratio
           scatter (robust MAD) across its rows is computed per line/
           doublet, and the group is tagged CLOSE (MJD span < -mjd_close)
           or not -- closely-spaced groups are the cleaner test of a
           method's intrinsic precision, since the true astrophysical flux
           should not have changed much between them.

        DOUBLETS has three kinds of ground truth (see its own comment for the
        full rationale):

            fixed          [SIII] 9531/9069 = 2.44 and [NII] 6583/6548 = 3.0
                           were both explicitly confirmed correct for this
                           project (the SIII value matches lvm_gaussfit.py,
                           overriding data/dap_lines.txt's flagged-incorrect
                           wavelengths -- see sky_gaussfit.py's 260909 History
                           entry). [OIII] 5007/4959 uses a standard literature
                           atomic-physics value but has NOT been independently
                           verified against this project's own preferred
                           source the way SIII/NII were -- treat it as lower-
                           confidence until checked. [OI] 6364/6300 is
                           deliberately NOT included here at all -- see the
                           DOUBLETS comment.
            bounded_above  Hb/Ha: reddening can only ever push this ratio
                           DOWN from the intrinsic Case B (no-reddening)
                           value of 1/2.86 = 0.350, never up -- a measured
                           ratio above 0.350 is a real red flag regardless of
                           nebular conditions. There is no single correct
                           value below that ceiling (real reddening varies by
                           target), only the ceiling itself.
            free           [SII] 6731/6716 and [OII] 3729/3726 are both
                           density-dependent with no fixed value at all --
                           useful only through repeat_scatter's across-repeat
                           consistency (they should still be roughly stable
                           for the same field even without a universal
                           truth), never as a per-row deviation check
                           (DEV_SIGMA is always NaN for these two).

Primary routines::

        fit_nebular_row   fit every NEBULAR_LINES entry on one spectrum.
        fit_file          loop fit_nebular_row over every row of one file.
        repeat_scatter     group one file's row table by tileid, compute
                          per-group flux/ratio scatter.

History::

    260909  ksl  Coding begun.

'''

import argparse
import re

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack

from sky_gaussfit import NEBULAR_LINES, fit_double_gaussian_to_spectrum
from sky_nebular_leak_eval import fit_leak_line, resolve_velocity, LMC_VEL, SMC_VEL


def _usage_from_doc(doc):
    '''__doc__ truncated just before "History:", so -h stays short.'''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


# (name, line_a, line_b, kind, value, confidence) -- ratio always reported
# as line_b/line_a. Three kinds of ground truth:
#
#   'fixed'         value is the true ratio (atomic physics, independent of
#                   nebular conditions); any deviation either direction is
#                   suspicious.
#   'bounded_above' value is the maximum physically possible ratio; a
#                   measured ratio ABOVE value is suspicious, below is not
#                   (just means real reddening/conditions, not an error).
#   'free'          no known truth value at all (value=None) -- only usable
#                   via repeat_scatter's across-repeat consistency, never as
#                   a per-row deviation check.
# Ordered by increasing (shorter-wavelength-member) rest wavelength, so any
# table built from this list reads left-to-right blue-to-red -- purely a
# display convention, order has no effect on any computation here.
DOUBLETS = [
    # SII/OII: both density-dependent doublets with no single fixed value --
    # only meaningful via repeat_scatter (should still be roughly constant
    # for repeat visits to the same field, even without a universal truth).
    ('oii_ratio',   'oii_a',  'oii_b',  'free',         None,  'density-dependent; also only 3.8 A '
                                                               'line separation -- see NEBULAR_LINES '
                                                               'comment in sky_gaussfit.py'),
    # Hb/Ha: reddening only ever REDUCES Hb relative to Ha (dust extinguishes
    # the bluer line more), so the observed ratio can never exceed the
    # intrinsic no-reddening (Case B, Te~1e4K) value -- 1/2.86 = 0.350. A
    # measured ratio above that is unphysical regardless of reddening/
    # nebular conditions, and a real diagnostic of a subtraction/fitting
    # problem, not an astrophysical one.
    ('hb_ha_ratio', 'ha',     'hb',     'bounded_above', 0.350, 'Case B no-reddening max; literature'),
    ('oiii_ratio',  'oiii_a', 'oiii_b', 'fixed',        2.98,  'literature'),
    ('nii_ratio',   'nii_a',  'nii_b',  'fixed',        3.0,   'confirmed'),
    # oi_ratio deliberately omitted: [OI] 6300/6364 is dominated by residual
    # sky-subtraction error (it's the same airglow doublet as sky6300/
    # sky6363, see sky_gaussfit.resolve_nebular_lines), not real nebular
    # flux, in most fibers -- a "ground truth" check here mostly measures
    # airglow-subtraction leftovers, not nebular-line recovery.
    ('sii_ratio',   'sii_a',  'sii_b',  'free',         None,  'density-dependent, no fixed value'),
    ('siii_ratio',  'siii_a', 'siii_b', 'fixed',        2.44,  'confirmed'),
]

# Placeholder tileid confirmed (by inspection) to span many unrelated
# targets over huge RA/Dec/MJD ranges -- not a genuine repeat pointing.
PLACEHOLDER_TILEID = 11111
# A genuine repeat group's pointing should not wander more than this
# (degrees) -- guards against any other non-tile-like grouping artifact.
MAX_GROUP_POS_SCATTER_DEG = 0.05


def fit_oii_doublet(wave, flux, wave0_a, wave0_b, half_width, sigma_guess=1.0):
    '''
    Joint double-Gaussian fit for the blended oii_a/oii_b pair, via
    lvm_gaussfit.fit_double_gaussian_to_spectrum (shared FWHM, independent
    amplitudes/centers, one constant background) -- reused directly rather
    than duplicated, since it already does exactly this for this same
    line pair.

    At 3.8 A rest-frame separation, two independent fit_leak_line singlet
    fits cannot deblend this pair (see the NEBULAR_LINES comment in
    sky_gaussfit.py); fit_double_gaussian_to_spectrum fits both
    simultaneously in one window instead. It accepts an input table
    without an ERROR column (falls back to unweighted fitting), which is
    what's built here, since this family's WAVE/FLUX/SKY files don't
    carry per-pixel errors.

    Parameters
    ----------
    wave, flux : ndarray
    wave0_a, wave0_b : float
        Doppler-shifted center wavelengths of the two components.
    half_width : float
        Half-width of the shared fit window around the pair's midpoint.
    sigma_guess : float
        Initial (shared) Gaussian sigma guess -- converted to FWHM
        (x2.3548) for fit_double_gaussian_to_spectrum's own convention.

    Returns
    -------
    (dict, dict) or (None, None)
        Two fit_leak_line-shaped dicts (amp, amp_err, lam0, sigma0, npix)
        -- amp/amp_err here are fit_double_gaussian_to_spectrum's
        *integrated flux* (not peak amplitude), but since both components
        share one fitted FWHM the a/b ratio is identical either way, and
        this pair is only ever used via that ratio (DOUBLETS' oii_ratio
        is 'free', never compared to an absolute value). One component
        per component in (a, b) order, or (None, None) if the fit fails.
    '''
    mid = 0.5 * (wave0_a + wave0_b)
    spectrum_table = Table([wave, flux], names=['WAVE', 'FLUX'])
    try:
        qtab, _fit_table = fit_double_gaussian_to_spectrum(
            spectrum_table, line='oii',
            init_wavelength1=wave0_a, init_wavelength2=wave0_b,
            init_fwhm=sigma_guess * 2.3548,
            wavelength_min=mid - half_width, wavelength_max=mid + half_width)
    except Exception:
        return None, None
    row = qtab[0]
    npix = int(np.sum((wave >= mid - half_width) & (wave <= mid + half_width) & np.isfinite(flux)))
    fields = ('flux_oii_a', 'eflux_oii_a', 'flux_oii_b', 'eflux_oii_b',
             'wave_oii_a', 'wave_oii_b', 'fwhm_oii_a')
    if not all(np.isfinite(row[f]) for f in fields):
        return None, None
    sigma0 = row['fwhm_oii_a'] / 2.3548
    # back_oii_ab is a flux-DENSITY background (same units as the input FLUX
    # column), unlike amp above (integrated flux) -- consistent with how
    # PlotSkySubNebEval.py converts amp back to peak amplitude (also flux
    # density) before adding this background for display.
    bkg = float(row['back_oii_ab']) if np.isfinite(row['back_oii_ab']) else 0.0
    fit_a = dict(amp=float(row['flux_oii_a']), amp_err=float(row['eflux_oii_a']),
                lam0=float(row['wave_oii_a']), sigma0=float(sigma0), npix=npix, bkg=bkg)
    fit_b = dict(amp=float(row['flux_oii_b']), amp_err=float(row['eflux_oii_b']),
                lam0=float(row['wave_oii_b']), sigma0=float(sigma0), npix=npix, bkg=bkg)
    return fit_a, fit_b


def fit_nebular_row(wave, flux, vel=0.0, sigma_guess=1.0):
    '''
    Fit every NEBULAR_LINES entry on one sky-subtracted spectrum.

    oii_a/oii_b are fit jointly (fit_oii_doublet) rather than as two
    independent singlets, since they're too close together (3.8 A) for
    independent fits to deblend -- see fit_oii_doublet. Every other line
    is fit independently (fit_leak_line), unchanged.

    Parameters
    ----------
    wave, flux : ndarray
    vel : float
        Nebular systemic velocity (km/s).
    sigma_guess : float

    Returns
    -------
    dict
        line name -> fit_leak_line-shaped result dict, or None if that
        line's fit failed.
    '''
    zz = 1.0 + vel / 3e5
    out = {}
    oii = {}
    for name, center, wmin, wmax in NEBULAR_LINES:
        wave0 = zz * center
        wlo, whi = zz * wmin, zz * wmax
        if name in ('oii_a', 'oii_b'):
            # fit_oii_doublet still wants a single symmetric half-width
            # (its window is already symmetric by construction -- both
            # oii_a/oii_b share the same [3717,3737] window, see
            # sky_gaussfit.NEBULAR_LINES), unlike fit_leak_line below.
            oii[name] = (wave0, max(wave0 - wlo, whi - wave0))
            continue
        out[name] = fit_leak_line(wave, flux, wave0, wlo, whi, sigma_guess=sigma_guess)

    if 'oii_a' in oii and 'oii_b' in oii:
        wave0_a, hw_a = oii['oii_a']
        wave0_b, hw_b = oii['oii_b']
        out['oii_a'], out['oii_b'] = fit_oii_doublet(
            wave, flux, wave0_a, wave0_b, max(hw_a, hw_b), sigma_guess=sigma_guess)

    return out


def _doublet_ratio(fits_dict, name_a, name_b, kind, value, snr_min=5.0):
    '''(ratio, ratio_err, dev_sigma) from two fit_leak_line results, or
    all-NaN if either line's fit failed, has a non-positive amplitude, or
    either line's own SNR is below snr_min.

    The SNR gate matters: most fibers in an arbitrary sample are not
    pointed at a real emission-line target, so both lines are noise-
    dominated there, and the ratio of two near-zero noisy amplitudes is
    unstable and NOT centered on the true ratio (verified on real data:
    unfiltered median SIII ratio across 140 mixed-target rows was
    ~2.98-3.00, nowhere near the true 2.44 -- restricting to rows with
    both lines' SNR>5 dropped it to 2.41, right on target). Without this
    gate the aggregate statistics are dominated by fit-to-noise garbage,
    not a meaningful measurement.

    dev_sigma's meaning depends on kind (see DOUBLETS): for 'fixed' it is
    signed deviation from the true ratio in either direction; for
    'bounded_above' it is deviation above the physical maximum (positive
    = a real violation, negative = fine, just reddened/different
    conditions); for 'free' it is always NaN (no truth value exists).'''
    fa, fb = fits_dict.get(name_a), fits_dict.get(name_b)
    if fa is None or fb is None or fa['amp'] <= 0 or fb['amp'] <= 0:
        return np.nan, np.nan, np.nan
    snr_a = fa['amp'] / fa['amp_err'] if fa['amp_err'] > 0 else 0.0
    snr_b = fb['amp'] / fb['amp_err'] if fb['amp_err'] > 0 else 0.0
    if snr_a < snr_min or snr_b < snr_min:
        return np.nan, np.nan, np.nan
    ratio = fb['amp'] / fa['amp']
    ratio_err = ratio * np.sqrt((fb['amp_err'] / fb['amp']) ** 2 + (fa['amp_err'] / fa['amp']) ** 2)
    if kind == 'free' or value is None or ratio_err <= 0:
        dev_sigma = np.nan
    else:
        dev_sigma = (ratio - value) / ratio_err
    return ratio, ratio_err, dev_sigma


def fit_file(fits_file, vel=None, lmc=False, smc=False, sigma_guess=1.0, snr_min=5.0):
    '''
    Fit every row of one SkySub*.py output file.

    Parameters
    ----------
    fits_file : str or Path
    vel : float, optional
        Explicit nebular systemic velocity (km/s), overriding lmc/smc AND
        the per-row DRP_ALL['Redshift'] auto-lookup below. Only use this
        for a file where every row shares one real target velocity --
        most files here are per-tile/per-survey mixes (see Redshift
        auto-lookup), so leaving this None is usually correct even for a
        known-LMC run.
    lmc, smc : bool
        Shortcuts for LMC_VEL/SMC_VEL, same precedence as vel above.
    sigma_guess : float
    snr_min : float
        Minimum per-line SNR (both lines) required to report a doublet
        ratio for a row -- see _doublet_ratio.

    Returns
    -------
    astropy.table.Table
        One row per spectrum: FILE, ROUTINE, VARIANT, ROW, TILEID, MJD,
        EXPNUM, SCI_RA, SCI_DEC, SURVEY, VEL (the per-row velocity
        actually used -- see resolve_velocity), then <line>_FLUX/
        _FLUX_ERR/_SNR for every NEBULAR_LINES entry, then <doublet>_
        RATIO/_RATIO_ERR/_DEV_SIGMA for every DOUBLETS entry.

    Notes
    -----
    Unless vel/lmc/smc is given, the velocity used for each row is looked
    up from that row's own DRP_ALL['Redshift'] (SummarizeCframe.py's RA/
    Dec-based LMC/SMC/Plane/HighLat classification -- see
    resolve_velocity) rather than a single value assumed for the whole
    file. A file can and often does mix rows from different surveys
    (e.g. the repeat-observation test sets span HighLat and SMC
    pointings together) -- a single global vel would silently mis-fit
    every row from a different survey than the one assumed.
    '''
    with fits.open(fits_file) as hdul:
        hdr = hdul[0].header
        routine = hdr.get('TITLE', '')
        variant = hdr.get('METHOD', '')
        wave = np.asarray(hdul['WAVE'].data, dtype=float)
        flux = np.asarray(hdul['FLUX'].data, dtype=float)
        drp = Table(hdul['DRP_ALL'].data)

    rows = []
    for i in range(flux.shape[0]):
        row_redshift = (float(drp['Redshift'][i])
                        if 'Redshift' in drp.colnames and np.isfinite(drp['Redshift'][i])
                        else None)
        row_vel = resolve_velocity(vel, lmc, smc, drp_redshift=row_redshift)
        fits_dict = fit_nebular_row(wave, flux[i], vel=row_vel, sigma_guess=sigma_guess)
        row = dict(FILE=str(fits_file), ROUTINE=routine, VARIANT=variant, ROW=i,
                  TILEID=int(drp['tileid'][i]) if 'tileid' in drp.colnames else -1,
                  MJD=float(drp['mjd'][i]) if 'mjd' in drp.colnames else np.nan,
                  EXPNUM=int(drp['expnum'][i]) if 'expnum' in drp.colnames else -1,
                  SCI_RA=float(drp['sci_ra'][i]) if 'sci_ra' in drp.colnames else np.nan,
                  SCI_DEC=float(drp['sci_dec'][i]) if 'sci_dec' in drp.colnames else np.nan,
                  SURVEY=str(drp['Survey'][i]) if 'Survey' in drp.colnames else '',
                  VEL=row_vel)
        for name, _center, _wmin, _wmax in NEBULAR_LINES:
            f = fits_dict[name]
            if f is None:
                row[f'{name}_FLUX'] = np.nan
                row[f'{name}_FLUX_ERR'] = np.nan
                row[f'{name}_SNR'] = np.nan
            else:
                row[f'{name}_FLUX'] = f['amp']
                row[f'{name}_FLUX_ERR'] = f['amp_err']
                row[f'{name}_SNR'] = f['amp'] / f['amp_err'] if f['amp_err'] > 0 else np.nan
        for dname, name_a, name_b, kind, value, _conf in DOUBLETS:
            ratio, ratio_err, dev_sigma = _doublet_ratio(fits_dict, name_a, name_b, kind, value, snr_min=snr_min)
            row[f'{dname.upper()}'] = ratio
            row[f'{dname.upper()}_ERR'] = ratio_err
            row[f'{dname.upper()}_DEV_SIGMA'] = dev_sigma
        rows.append(row)

    return Table(rows=rows)


def repeat_scatter(row_table, mjd_close=7.0):
    '''
    Group row_table by TILEID (excluding PLACEHOLDER_TILEID) and compute
    per-group flux/ratio scatter for every line and doublet.

    Parameters
    ----------
    row_table : astropy.table.Table
        Output of fit_file (or several stacked together, still one
        ROUTINE/VARIANT per group -- groups spanning multiple
        routines are not meaningful and are silently split apart, since
        grouping is done separately per ROUTINE/VARIANT/TILEID triple).
    mjd_close : float
        MJD span (days) below which a group is tagged CLOSE.

    Returns
    -------
    astropy.table.Table
        One row per (ROUTINE, VARIANT, TILEID) group with >=2 rows:
        N_REPEATS, MJD_SPAN, CLOSE, POS_SCATTER_DEG, then <line>_
        FLUX_MED/_FLUX_MAD for every NEBULAR_LINES entry, and
        <doublet>_MED/_MAD for every DOUBLETS entry.
    '''
    line_names = [n for n, *_ in NEBULAR_LINES]
    doublet_names = [d.upper() for d, *_ in DOUBLETS]

    groups = {}
    for row in row_table:
        if row['TILEID'] == PLACEHOLDER_TILEID:
            continue
        key = (row['ROUTINE'], row['VARIANT'], row['TILEID'])
        groups.setdefault(key, []).append(row)

    out_rows = []
    for (routine, variant, tileid), rows in groups.items():
        if len(rows) < 2:
            continue
        ra = np.array([r['SCI_RA'] for r in rows], dtype=float)
        dec = np.array([r['SCI_DEC'] for r in rows], dtype=float)
        pos_scatter = float(np.nanmax([np.nanstd(ra), np.nanstd(dec)])) if len(rows) > 1 else 0.0
        if pos_scatter > MAX_GROUP_POS_SCATTER_DEG:
            continue

        mjd = np.array([r['MJD'] for r in rows], dtype=float)
        mjd_span = float(np.nanmax(mjd) - np.nanmin(mjd)) if np.isfinite(mjd).any() else np.nan

        out = dict(ROUTINE=routine, VARIANT=variant, TILEID=tileid, N_REPEATS=len(rows),
                  MJD_SPAN=mjd_span, CLOSE=bool(np.isfinite(mjd_span) and mjd_span < mjd_close),
                  POS_SCATTER_DEG=pos_scatter)
        for name in line_names:
            vals = np.array([r[f'{name}_FLUX'] for r in rows], dtype=float)
            finite = vals[np.isfinite(vals)]
            out[f'{name}_FLUX_MED'] = float(np.median(finite)) if len(finite) else np.nan
            out[f'{name}_FLUX_MAD'] = (float(1.4826 * np.median(np.abs(finite - np.median(finite))))
                                       if len(finite) > 1 else np.nan)
        for dname in doublet_names:
            vals = np.array([r[dname] for r in rows], dtype=float)
            finite = vals[np.isfinite(vals)]
            out[f'{dname}_MED'] = float(np.median(finite)) if len(finite) else np.nan
            out[f'{dname}_MAD'] = (float(1.4826 * np.median(np.abs(finite - np.median(finite))))
                                   if len(finite) > 1 else np.nan)
        out_rows.append(out)

    return Table(rows=out_rows) if out_rows else Table()


def print_comparison(repeats_table, close_only=True):
    '''Print one summary row per ROUTINE/VARIANT: median (across groups)
    scatter for each doublet, restricted to CLOSE groups by default --
    the head-to-head "which method is more precise" table.'''
    if len(repeats_table) == 0:
        print('No repeat groups found (after excluding the placeholder '
             f'tileid {PLACEHOLDER_TILEID} and position-scatter outliers).')
        return
    t = repeats_table
    if close_only:
        t = t[t['CLOSE']]
    if len(t) == 0:
        print('No closely-spaced repeat groups found; try -mjd_close with a larger value, '
             'or pass close_only=False.')
        return

    doublet_names = [d.upper() for d, *_ in DOUBLETS]
    routines = sorted(set(zip(t['ROUTINE'], t['VARIANT'])))
    print(f"\nDoublet-ratio scatter across {'closely-spaced ' if close_only else ''}"
         f"repeat groups (MAD, smaller = more precise):")
    header = f"{'Routine':<12}{'Variant':<20}{'N groups':>9}"
    for dname in doublet_names:
        header += f"{dname:>16}"
    print(header)
    for routine, variant in routines:
        sel = (t['ROUTINE'] == routine) & (t['VARIANT'] == variant)
        line = f"{routine:<12}{variant:<20}{int(sel.sum()):>9}"
        for dname in doublet_names:
            mad = t[f'{dname}_MAD'][sel]
            finite = np.asarray(mad)[np.isfinite(mad)]
            med_mad = np.median(finite) if len(finite) else np.nan
            line += f"{med_mad:>16.4f}" if np.isfinite(med_mad) else f"{'--':>16}"
        print(line)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=_usage_from_doc(__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fits_file', nargs='+', help='SkySub*.py output FITS file(s)')
    parser.add_argument('-v', dest='vel', type=float, default=None, help='Nebular systemic velocity (km/s)')
    parser.add_argument('-lmc', action='store_true', help=f'Use the LMC velocity (~{LMC_VEL:.0f} km/s)')
    parser.add_argument('-smc', action='store_true', help=f'Use the SMC velocity (~{SMC_VEL:.0f} km/s)')
    parser.add_argument('-sigma', dest='sigma_guess', type=float, default=1.0,
                        help='Initial Gaussian sigma guess (Angstrom)')
    parser.add_argument('-mjd_close', type=float, default=7.0,
                        help='MJD span (days) below which a repeat group is "closely spaced"')
    parser.add_argument('-snr_min', type=float, default=5.0,
                        help='Minimum per-line SNR (both lines) to report a doublet ratio')
    parser.add_argument('-out', default='nebeval', help='Output filename root')
    args = parser.parse_args()

    if args.vel is not None or args.lmc or args.smc:
        vel_msg = f'{resolve_velocity(args.vel, args.lmc, args.smc):.1f} km/s (explicit override, all rows)'
    else:
        vel_msg = "per-row, from each file's own DRP_ALL['Redshift'] (LMC/SMC/Plane/HighLat)"
    print(f'Evaluating {len(args.fits_file)} file(s), vel={vel_msg}, snr_min={args.snr_min:.1f}')

    tables = []
    for f in args.fits_file:
        print(f'  fitting {f} ...')
        tables.append(fit_file(f, vel=args.vel, lmc=args.lmc, smc=args.smc,
                               sigma_guess=args.sigma_guess, snr_min=args.snr_min))
    row_table = vstack(tables) if tables else Table()
    row_out = f'{args.out}_lines.fits'
    row_table.write(row_out, overwrite=True)
    print(f'Wrote {row_out} ({len(row_table)} rows)')

    repeats_table = repeat_scatter(row_table, mjd_close=args.mjd_close)
    repeats_out = f'{args.out}_repeats.fits'
    repeats_table.write(repeats_out, overwrite=True)
    print(f'Wrote {repeats_out} ({len(repeats_table)} groups)')

    print_comparison(repeats_table, close_only=True)


if __name__ == '__main__':
    main()
