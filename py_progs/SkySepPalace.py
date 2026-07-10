#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Perform sky subtraction on an XCframe file using the PALACE airglow
    model instead of the ESO sky model (SkySepESO.py's approach).  A
    first, deliberately scoped-down cut:

        SKY = sum_species  a_species * PALACE_species
            + a_moon * MOON                              (all a >= 0)

    fit jointly by non-negative least squares (all 9 PALACE species
    amplitudes plus one MOON amplitude at once, over the whole spectrum
    -- no line masking is needed since these templates already include
    real predicted line shapes) against the nearest sky fiber's own
    spectrum.  MOON is a single fixed, precomputed spectral-shape
    template (data/moon_base_spectrum.dat, see MakeMoonBase.py) with
    exactly one free amplitude -- deliberately not the flexible multi-
    knot B-spline envelope SkySubDev2.py/SkyDecomp uses for its own
    Moon component (see Notes): a single scalar is a much more
    constrained, more physical assumption than a shape with enough free
    parameters to fit essentially any smooth continuum.  Unlike a bare
    solar spectrum, this template already accounts for both the Moon's
    own wavelength-dependent albedo AND atmospheric (Rayleigh/aerosol)
    scattering, which empirically dominates and reverses the albedo's
    reddening (see MakeMoonBase.py) -- at the cost of using one fixed
    reference phase-angle/airmass rather than this observation's real
    geometry (see Notes).

Command line usage (if any):

    usage: SkySepPalace.py [-delta N] [-out ROOT] filename

    Arguments::

        filename    XCframe FITS file to process

    Options::

        -delta N    process every N-th row; useful for quick tests
                    (default: 1 = all rows)
        -out ROOT   output filename root; default is <stem>_palace

Description:

    For each row (fiber):

    1. RA/Dec/obstime for the science fiber and both sky telescopes are
       read from DRP_ALL; the nearer sky telescope is identified from
       its angular separation from the science fiber (as in
       SkySepESO.py/SkySubOrig.py).

    2. PALACE (see PalaceObs.py) is queried once per unique (near sky
       position, obstime) pair -- cached across rows, since real
       full-IFU XCframe files share only two sky-telescope positions
       across ~1700+ science fibers -- for all 9 species (OH, O2, HO2,
       FeO, Na, K, O, N, H) at high native resolution (resol=20000,
       dlam~0.2 A) over the file's wavelength range.  PALACE's own
       resolution far exceeds LVM's, so these native spectra can be
       rebinned/convolved onto any real per-row LSF without needing to
       account for PALACE's own native width (see Notes).

    3. Each species' native spectrum is flux-conserving rebinned and
       LSF-convolved onto the instrument wavelength grid using this
       row's real LSF (from the file's own LSF extension if present,
       else a reference curve, else a flat default -- same 3-tier
       priority as SkySubDev2.py), via
       sky_decomp.fit.rebin_and_convolve.

    4. The resulting 9 species templates plus the MOON template are fit
       to the near sky fiber's flux by non-negative least squares
       (scipy.optimize.nnls) -- one joint fit, all 10 components, full
       spectrum, no masking.

    5. SKY = the NNLS fit; FLUX (output) = science flux - SKY.

    6. If sky_mask.fits is found, per-arm continuum-fit-quality stats
       (GetSkyCont.arm_continuum_stats/flatten_arm_stats) are also
       computed -- SCI_MED_<arm>/SKY_MED_<arm> (and NMAD/RMS/skew),
       directly comparable to SkySepESO.py's own columns of the same
       name, evaluated against the raw pre-subtraction science/sky
       spectra rather than the final sky-subtracted result.  See
       one_drp()'s Notes for exactly what "continuum" means here (only
       the HO2/FeO/MOON components of the fit, not the full 10).

    DRP_ALL gains per-row columns PALACE_<SPECIES> (the 9 species
    amplitudes), PALACE_MOON (the MOON amplitude), QA_FLAGS,
    ERROR_MSG, and -- if sky_mask.fits is found -- the SCI_*/SKY_*
    continuum-fit-quality columns from step 6.

    QA flag bits stored in DRP_ALL['QA_FLAGS']::

        0x01  NANDATA    NaN/inf found in input flux or sky data
        0x04  MODELFAIL  PALACE prediction failed for this row's
                         geometry (e.g. sky position below the horizon)
        0x08  FAILED     row raised an exception; spectrum filled with
                         NaN

Notes:

    On MOON: an early version of this template used the bare vendored
    solar spectrum (moonlight is reflected sunlight), then a version
    that multiplied it by the Kieffer & Stone (2005) ROLO lunar-albedo
    curve (data/moon_rolo_albedo.dat).  Both were wrong in the same
    direction: validated (260709) against the real ESO Sky Model's own
    MOON output for a fixed real geometry, the ROLO-corrected spectrum
    got *redder* toward longer wavelengths while the real, atmospherically
    -scattered MOON spectrum gets *bluer* -- Rayleigh/aerosol scattering
    (favouring blue) dominates over and reverses the albedo reddening,
    since the quantity being fit is post-atmosphere sky brightness, not
    top-of-atmosphere lunar reflectance.  The current MOON template
    (data/moon_base_spectrum.dat, MakeMoonBase.py) is instead the real
    ESO Sky Model's own MOON column, for one fixed reference geometry,
    used as a static shape -- see MakeMoonBase.py's module docstring for
    the full validation and reasoning.  This keeps SkySepPalace.py free
    of any *runtime* ESO dependency (the ESO model is only needed once,
    offline, to (re)generate the static file) while still using a
    physically complete shape.  The real per-observation dependence on
    lunar phase angle and airmass is not captured by this single fixed
    curve -- see the first bullet below.

    SkySubDev2.py (via Ivan Katkov's SkyDecomp) builds its own "Moon"
    component from the bare solar spectrum shape, multiplied by a
    flexible multi-knot B-spline basis (see sky_decomp/fit.py's
    _build_moon), giving the fit many free parameters to bend that
    shape into whatever smooth continuum is present.  That is close to
    unconstrained curve-fitting, not a physical moon model -- it never
    uses moon phase, moon altitude, or moon-source separation at all,
    and (per the validation above) the underlying shape it's bending is
    itself the wrong color for scattered moonlight.  The single-
    amplitude version here is a deliberately more conservative choice:
    it can only scale a physically-motivated shape, not reshape it, so
    it can't silently absorb continuum that isn't actually moon-colored.
    Zodiacal light is not included in either version.

    Deliberately NOT yet included, in order of what's needed next:

      - Any dependence of the MOON template's *shape* on this
        observation's real lunar phase angle or airmass (currently one
        fixed reference geometry for every row), or any physical
        prediction of its *amplitude* (currently entirely free, with no
        tie to how bright the actual Moon is).  See SkySepESO.py for
        how the ESO model predicts MOON/ZODI from real geometry, once
        that dependency -- or a per-observation call to MakeMoonBase.py's
        approach -- is worth taking on.  Zodiacal light continuum is
        also still missing entirely.

      - Any correction for PALACE's isatm=True output being an
        at-the-telescope (reddened) prediction, whereas LVM's flux
        calibration likely refers fluxes to above the atmosphere.
        Airglow originates near the mesopause (~90 km), so this is a
        real, wavelength- and airmass-dependent mismatch, not a
        rounding error -- flagged for follow-up, not yet addressed.

      - A sky-to-science line-scale correction (SkySepESO.py's
        ksl_bisection step).  Not clearly motivated yet without the
        above two pieces in place; the NNLS amplitudes already give
        each species its own free scaling, unlike ESO's single lumped
        "LINES" template.

    The rebin_and_convolve treatment of MOON (build the template at
    native resolution far exceeding the target, so no quadrature
    correction of its own width is needed) follows the same pattern
    Ivan Katkov's SkyDecomp uses for its solar template -- see
    py_progs/sky_decomp/fit.py around line 668.

History:

    260709 ksl Coding begun.  First cut: PALACE species only, no
        Moon/Zodi, no line-scale correction -- see Notes.
    260709 ksl Added a single-amplitude MOON component (bare solar-
        spectrum shape, no B-spline reshaping -- see Notes).
    260709 ksl Added SCI_MED_<arm>/SKY_MED_<arm> (etc) continuum-fit-
        quality columns, gated on sky_mask.fits, so SkySub_eval.py's
        Figure 5/6 continuum-separation plots -- previously silently
        empty for this script's output -- populate the same as
        SkySepESO.py's, for a direct comparison between the two.
    260709 ksl MOON now uses data/moon_base_spectrum.dat (see
        MakeMoonBase.py) instead of the bare solar spectrum: validated
        against the real ESO Sky Model that a lunar-albedo-only
        correction (data/moon_rolo_albedo.dat, Kieffer & Stone 2005)
        gets the color trend backwards once atmospheric scattering is
        accounted for -- see Notes.

'''

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
from scipy.optimize import nnls
from pathlib import Path

import PalaceObs
from sky_decomp.fit import rebin_and_convolve
from SkySubOrig import obstime_to_mjd
from SkySubDev2 import _load_ref_lsf, _interp_lsf_to_wave

try:
    from GetSkyCont import (load_mask, _interp_mask_to_wave,
                            arm_continuum_stats, flatten_arm_stats)
    _HAVE_MASK = True
except ImportError:
    _HAVE_MASK = False

# ──────────────────────────────────────────────────────────────
# QA flag bits
# ──────────────────────────────────────────────────────────────
QA_NANDATA   = 1
QA_MODELFAIL = 4
QA_FAILED    = 8

_QA_FLAG_NAMES = {
    QA_NANDATA:   'NANDATA',
    QA_MODELFAIL: 'MODELFAIL',
    QA_FAILED:    'FAILED',
}

_USAGE = '''Usage:
  SkySepPalace.py [-delta N] [-out ROOT] filename

Arguments:
  filename         XCframe FITS file to process

Options:
  -delta N         step size through rows for quick tests (default: 1)
  -out ROOT        output filename root (default: <stem>_palace)
'''

PALACE_SPECIES = PalaceObs.SPECIES   # ['OH','O2','HO2','FeO','Na','K','O','N','H']
COMPONENTS     = PALACE_SPECIES + ['MOON']

# The subset of COMPONENTS that are smooth/continuum-like rather than
# line-dominated -- used to reconstruct a "continuum-only" spectrum from an
# already-fit amplitude vector, for the SCI_MED/SKY_MED continuum-fit-quality
# columns (see one_drp() Notes and GetSkyCont.arm_continuum_stats).
CONTINUUM_COMPONENTS = ['HO2', 'FeO', 'MOON']
_CONT_IDX = [COMPONENTS.index(c) for c in CONTINUUM_COMPONENTS]

_NATIVE_RESOL = 20000
_NATIVE_DLAM  = 2e-5   # micron, ~0.2 A -- fine relative to LVM's ~1.3-2 A LSF
_WAVE_PAD_ANG = 5.0    # padding at each end of the native PALACE request

DEFAULT_LSF_FWHM = 1.3   # Angstroms, last-resort fallback (matches XSkySepIvan.py)


# ──────────────────────────────────────────────────────────────
# MOON: a static, precomputed moonlight spectral shape -- the vendored
# solar spectrum reddened by the Moon's own albedo AND then further
# reshaped by atmospheric (Rayleigh/aerosol) scattering, which
# empirically dominates and reverses the albedo reddening (see
# MakeMoonBase.py's module docstring for the validation against the
# real ESO Sky Model that established this).  One free amplitude, no
# B-spline reshaping (see module docstring Notes).
# ──────────────────────────────────────────────────────────────

DEFAULT_MOON_BASE_FILE = (Path(__file__).resolve().parent.parent / 'data' /
                          'moon_base_spectrum.dat')

_moon_base_native = None


def _load_solar_native():
    '''
    Native-resolution MOON_BASE spectrum (Angstroms, already in the
    same air-wavelength convention as the instrument grid -- no
    vac_to_air conversion, unlike a raw vendored solar/atomic-line
    source), normalised to unit median flux.  Despite the name (kept
    for backward compatibility with build_design_matrix's MOON case),
    this is data/moon_base_spectrum.dat's MOON_BASE column -- see
    MakeMoonBase.py -- not the bare solar spectrum.  Loaded once and
    cached at module scope: this is a single fixed reference shape (see
    MakeMoonBase.py Notes), not a function of this observation's
    ra/dec/obstime.
    '''
    global _moon_base_native
    if _moon_base_native is None:
        data = np.loadtxt(DEFAULT_MOON_BASE_FILE)
        wave_native_ang = np.asarray(data[:, 0], dtype=float)
        flux_native = np.asarray(data[:, 2], dtype=float)
        flux_native = flux_native / np.nanmedian(flux_native)
        _moon_base_native = (wave_native_ang, flux_native)
    return _moon_base_native


# ──────────────────────────────────────────────────────────────
# PALACE template fetch (cached per unique sky position/time)
# ──────────────────────────────────────────────────────────────

_template_cache = {}


def _native_wave_range(wave_ang):
    lammin = max(0.30, (float(wave_ang.min()) - _WAVE_PAD_ANG) / 1e4)
    lammax = min(2.50, (float(wave_ang.max()) + _WAVE_PAD_ANG) / 1e4)
    return lammin, lammax


def get_palace_templates(ra, dec, obstime, wave_ang):
    '''
    Native-resolution PALACE species spectra (fine wavelength grid, high
    resol) spanning wave_ang, for one sky position/time.  Cached across
    calls with the same (rounded) ra/dec/obstime, since real XCframe
    files share only two sky-telescope positions across many science
    fibers.

    Returns dict: species name -> (wave_native_ang, flux_native) or None
    for species with no contribution in this wavelength range.
    '''
    key = (round(float(ra), 4), round(float(dec), 4), str(obstime))
    if key in _template_cache:
        return _template_cache[key]

    lammin, lammax = _native_wave_range(wave_ang)
    species_table, _combined = PalaceObs.predict(
        ra, dec, obstime, species_list=PALACE_SPECIES, want_combined=False,
        lammin=lammin, lammax=lammax, dlam=_NATIVE_DLAM, resol=_NATIVE_RESOL)

    out = {}
    for sp in PALACE_SPECIES:
        sub = species_table[species_table['species'] == sp] if len(species_table) else species_table
        if len(sub) == 0:
            out[sp] = None
        else:
            out[sp] = (np.array(sub['lam'], dtype=float) * 1e4,
                      np.array(sub['flux'], dtype=float))

    _template_cache[key] = out
    return out


def build_design_matrix(templates, wave_ang, row_fwhm_ang):
    '''
    Rebin + LSF-convolve each component's native template onto
    wave_ang: the 9 PALACE species plus MOON (see _load_solar_native).
    Species with no native contribution become all-zero columns (kept,
    not dropped, so the amplitude vector always has a fixed, named
    length/order).

    Returns (A, component_names): A is (n_wave, n_components).
    '''
    cols = []
    for comp in COMPONENTS:
        if comp == 'MOON':
            wave_native, flux_native = _load_solar_native()
        else:
            native = templates.get(comp)
            if native is None:
                cols.append(np.zeros_like(wave_ang, dtype=float))
                continue
            wave_native, flux_native = native
        col = rebin_and_convolve(wave_ang, wave_native, flux_native,
                                 row_fwhm_ang, lsf_in_wavelength=True)
        col = np.nan_to_num(np.asarray(col, dtype=float),
                            nan=0.0, posinf=0.0, neginf=0.0)
        cols.append(col)
    return np.column_stack(cols), list(COMPONENTS)


# ──────────────────────────────────────────────────────────────
# Per-row sky subtraction
# ──────────────────────────────────────────────────────────────

def one_drp(x, drp_all, row, wave, row_fwhm_ang, clean_mask=None):
    '''
    Sky-subtract a single spectrum (row) using a PALACE species + MOON
    fit against the nearest sky fiber.

    Returns (scitab, qa_flags, coeffs, err_msg).

    scitab has columns WAVE, FLUX (sky-subtracted), SKY.

    coeffs is a dict with one PALACE_<SPECIES> amplitude per species
    plus PALACE_MOON (from the joint NNLS fit against the near sky
    fiber), and -- if clean_mask is given -- per-arm continuum-fit-
    quality stats (GetSkyCont.arm_continuum_stats/flatten_arm_stats),
    keys 'sci_<stat>_<arm>'/'sky_<stat>_<arm>', directly comparable to
    SkySepESO.py's SCI_MED/SKY_MED columns:

      sky side: reuses the already-fit near-sky amplitudes, just
      reconstructed from only the continuum-like components
      (CONTINUUM_COMPONENTS = HO2, FeO, MOON) instead of all 10 --
      "how well does the continuum-only part of the fit we actually
      used track the real sky continuum".

      sci side: an independent continuum-only fit (same 3 components,
      but PALACE-evaluated at the *science* fiber's own ra/dec/obstime,
      matching SkySepESO's convention of fetching its sky model at the
      science position) restricted to clean_mask pixels only -- since,
      unlike SkySepESO's iteratively-downweighted continuum fit, this
      is a plain NNLS with no other way to keep real emission lines
      (science source or airglow) from biasing it.  Used only to derive
      sci_lines for this stat, exactly as in SkySepESO.py; not part of
      the subtracted SKY.

    err_msg is '' on success, else a description of what failed.
    '''
    qa_flags = 0
    try:
        flux     = np.array(x['FLUX'].data[row],     dtype=float)
        skye     = np.array(x['SKY_EAST'].data[row], dtype=float)
        skyw     = np.array(x['SKY_WEST'].data[row], dtype=float)
        obstime  = drp_all['obstime'][row]
        sci_ra   = drp_all['sci_ra'][row]
        sci_dec  = drp_all['sci_dec'][row]
        skye_ra  = drp_all['skye_ra'][row]
        skye_dec = drp_all['skye_dec'][row]
        skyw_ra  = drp_all['skyw_ra'][row]
        skyw_dec = drp_all['skyw_dec'][row]
    except Exception as e:
        msg = 'could not read row data: %s' % e
        print('Row %d: %s' % (row, msg))
        return None, QA_FAILED, {}, msg

    if not (np.all(np.isfinite(flux)) and np.all(np.isfinite(skye))
            and np.all(np.isfinite(skyw))):
        qa_flags |= QA_NANDATA
        flux = np.nan_to_num(flux, nan=0.0, posinf=0.0, neginf=0.0)
        skye = np.nan_to_num(skye, nan=0.0, posinf=0.0, neginf=0.0)
        skyw = np.nan_to_num(skyw, nan=0.0, posinf=0.0, neginf=0.0)

    sci_coord  = SkyCoord(ra=sci_ra  * u.degree, dec=sci_dec  * u.degree)
    skye_coord = SkyCoord(ra=skye_ra * u.degree, dec=skye_dec * u.degree)
    skyw_coord = SkyCoord(ra=skyw_ra * u.degree, dec=skyw_dec * u.degree)
    de = sci_coord.separation(skye_coord)
    dw = sci_coord.separation(skyw_coord)
    if de < dw:
        sky_near, near_ra, near_dec = skye, skye_ra, skye_dec
    else:
        sky_near, near_ra, near_dec = skyw, skyw_ra, skyw_dec

    try:
        templates = get_palace_templates(near_ra, near_dec, obstime, wave)
    except Exception as e:
        msg = 'PALACE prediction failed: %s' % e
        return None, QA_MODELFAIL, {}, msg

    A, component_names = build_design_matrix(templates, wave, row_fwhm_ang)
    amp, _resid = nnls(A, sky_near)
    sky      = A @ amp
    sci_flux = flux - sky

    scitab = Table([wave, sci_flux, sky], names=['WAVE', 'FLUX', 'SKY'])

    coeffs = {('palace_%s' % comp.lower()): float(a)
             for comp, a in zip(component_names, amp)}

    if clean_mask is not None:
        sky_cont = A[:, _CONT_IDX] @ amp[_CONT_IDX]
        sky_cont_resid = sky_near - sky_cont

        try:
            sci_templates = get_palace_templates(sci_ra, sci_dec, obstime, wave)
            A_sci, _sci_names = build_design_matrix(sci_templates, wave, row_fwhm_ang)
            A_sci_cont = A_sci[:, _CONT_IDX]
            fit_px = clean_mask & np.isfinite(flux)
            if fit_px.sum() > len(_CONT_IDX):
                amp_sci_cont, _r = nnls(A_sci_cont[fit_px], flux[fit_px])
                sci_cont  = A_sci_cont @ amp_sci_cont
                sci_lines = flux - sci_cont
                coeffs.update(flatten_arm_stats(
                    'sci', arm_continuum_stats(wave, sci_lines, clean_mask)))
        except Exception as e:
            print('Row %d: science-side continuum-quality fit failed (%s)' % (row, e))

        coeffs.update(flatten_arm_stats(
            'sky', arm_continuum_stats(wave, sky_cont_resid, clean_mask)))

    return scitab, qa_flags, coeffs, ''


# ──────────────────────────────────────────────────────────────
# Batch processing
# ──────────────────────────────────────────────────────────────

def do_all(filename, idelta=1, outroot=''):
    '''
    Process every N-th row of an XCframe file and write PALACE-only
    sky-subtracted output.

    Parameters
    ----------
    filename : str
        Path to the input XCframe FITS file.
    idelta : int
        Row step size (1 = all rows).
    outroot : str
        Output filename root; defaults to <stem>_palace.
    '''
    x = fits.open(filename)
    drp_all = Table(x['DRP_ALL'].data)
    wave    = np.array(x['WAVE'].data, dtype=float)

    # LSF source, in priority order (same convention as SkySubDev2.py):
    #  1. this file's own per-row, per-wavelength LSF extension
    #  2. a representative wavelength-dependent reference curve
    #     (data/lsf.fits)
    #  3. the flat DEFAULT_LSF_FWHM constant
    have_lsf_ext = 'LSF' in [h.name for h in x]
    if have_lsf_ext:
        lsf_fwhm_ext = np.array(x['LSF'].data, dtype=float)
        lsf_source = 'LSF extension (per-row)'
    else:
        ref_fwhm = None
        _data_dir = Path(__file__).parent.parent / 'data'
        for _candidate in [Path('lsf.fits'), _data_dir / 'lsf.fits']:
            if _candidate.exists():
                try:
                    _ref_wave, _ref_fwhm_raw = _load_ref_lsf(str(_candidate))
                    ref_fwhm = _interp_lsf_to_wave(_ref_wave, _ref_fwhm_raw, wave)
                except Exception as _e:
                    print('Warning: could not load reference LSF %s (%s)' % (_candidate, _e))
                break
        if ref_fwhm is not None:
            lsf_fwhm_ext = np.tile(ref_fwhm, (len(drp_all), 1))
            have_lsf_ext = True
            lsf_source   = 'reference curve (data/lsf.fits)'
        else:
            lsf_fwhm_ext = None
            lsf_source = 'constant default (%.2f A)' % DEFAULT_LSF_FWHM
            print('Warning: no LSF extension or reference LSF curve found; '
                  'using the constant %.2f A for every row' % DEFAULT_LSF_FWHM)
    print('LSF source: %s' % lsf_source)

    # Load sky-line mask once (same convention as SkySepESO.py) purely for
    # the continuum-fit-quality evaluation below; the PALACE fit itself is
    # unaffected if the mask is unavailable.
    clean_mask = None
    if _HAVE_MASK:
        _data_dir = Path(__file__).parent.parent / 'data'
        for _candidate in [Path('sky_mask.fits'), _data_dir / 'sky_mask.fits']:
            if _candidate.exists():
                try:
                    _mask_wave, _mask_arr = load_mask(str(_candidate))
                    clean_mask = _interp_mask_to_wave(_mask_wave, _mask_arr, wave)
                    print('Loaded sky mask for continuum-quality evaluation: %s' % _candidate)
                except Exception as _e:
                    print('Warning: could not load mask %s (%s)' % (_candidate, _e))
                break
    if clean_mask is None:
        print('Warning: sky_mask.fits not found; skipping continuum-quality columns')

    nan_spectrum = np.full(len(wave), np.nan)

    final_flux    = []
    final_sky     = []
    select        = []
    qa_flags_list = []
    coeffs_list   = []
    err_msgs      = []

    i = 0
    while i < len(drp_all):
        if have_lsf_ext:
            row_fwhm = lsf_fwhm_ext[i].copy()
            bad = ~np.isfinite(row_fwhm) | (row_fwhm <= 0)
            if bad.any():
                row_fwhm[bad] = DEFAULT_LSF_FWHM
        else:
            row_fwhm = np.full(len(wave), DEFAULT_LSF_FWHM)

        err_msg = ''
        try:
            ftab, row_flags, coeffs, err_msg = one_drp(x, drp_all, i, wave, row_fwhm,
                                                       clean_mask=clean_mask)
        except Exception as e:
            print('Row %d: exception (%s)' % (i, e))
            ftab = None
            row_flags = 0
            coeffs = {}
            err_msg = 'exception: %s' % e

        if ftab is None:
            row_flags |= QA_FAILED
            final_flux.append(nan_spectrum.copy())
            final_sky.append(nan_spectrum.copy())
        else:
            final_flux.append(np.array(ftab['FLUX']))
            final_sky.append(np.array(ftab['SKY']))

        select.append(i)
        qa_flags_list.append(row_flags)
        coeffs_list.append(coeffs)
        err_msgs.append(err_msg)
        i += idelta
        if i % 100 == 0:
            print('Completed %6d of %d in steps of %d' % (i, len(drp_all), idelta))

    n_failed    = sum(1 for f in qa_flags_list if f & QA_FAILED)
    n_modelfail = sum(1 for f in qa_flags_list if f & QA_MODELFAIL)
    n_warned    = sum(1 for f in qa_flags_list
                      if f != 0 and not (f & (QA_FAILED | QA_MODELFAIL)))
    print('\nProcessed %d rows: %d failed (NaN fill), %d model-fetch '
          'failures, %d with warnings'
          % (len(select), n_failed, n_modelfail, n_warned))

    if outroot == '':
        stem = os.path.splitext(os.path.basename(filename))[0]
        outroot = '%s_palace' % stem

    flagged = [(select[j], qa_flags_list[j], err_msgs[j])
               for j in range(len(select)) if qa_flags_list[j] != 0]
    if flagged:
        report_lines = ['Rows with QA flags:',
                        '  %-8s %-22s %s' % ('Row', 'Flags', 'Reason')]
        for orig_row, flags, msg in flagged:
            active = [name for bit, name in _QA_FLAG_NAMES.items() if flags & bit]
            report_lines.append('  %-8d %-22s %s'
                                % (orig_row, ','.join(active), msg or '(no message)'))
        report = '\n'.join(report_lines)
        print(report)

        errfile = '%s_errors.txt' % outroot
        with open(errfile, 'w') as fh:
            fh.write(report + '\n')
        print('Wrote error summary to %s' % errfile)

    hdu1 = fits.PrimaryHDU(data=None)
    hdu1.header['Title'] = 'SkySepPalace'
    hdu2 = fits.ImageHDU(data=wave,                name='WAVE')
    hdu3 = fits.ImageHDU(data=np.array(final_flux), name='FLUX')
    hdu4 = fits.ImageHDU(data=np.array(final_sky),  name='SKY')

    xtab = Table(drp_all[select])
    xtab['QA_FLAGS'] = np.array(qa_flags_list, dtype=np.int32)
    _known_keys = set()
    for _comp in COMPONENTS:
        _key = 'palace_%s' % _comp.lower()
        _known_keys.add(_key)
        xtab['PALACE_%s' % _comp.upper()] = np.array(
            [c.get(_key, np.nan) for c in coeffs_list], dtype=np.float32)
    # Continuum-fit-quality columns (SCI_MED_<arm>/SKY_MED_<arm>/etc, see
    # one_drp() Notes) -- present only if a sky mask was found.
    _cont_keys = sorted({k for d in coeffs_list for k in d} - _known_keys)
    for _key in _cont_keys:
        xtab[_key.upper()] = np.array(
            [d.get(_key, np.nan) for d in coeffs_list], dtype=np.float32)
    xtab['ERROR_MSG'] = np.array([m[:200] for m in err_msgs])
    if 'obstime' in xtab.colnames:
        xtab['mjd'] = obstime_to_mjd(xtab['obstime'])
    hdu5 = fits.BinTableHDU(xtab, name='DRP_ALL')

    dwave = wave[1] - wave[0] if len(wave) > 1 else 0.5
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [1, 1]
    wcs.wcs.crval = [float(wave[0]), 0]
    wcs.wcs.cdelt = [float(dwave), 1]
    wcs.wcs.ctype = ['WAVE', 'LINE']
    for hdu in (hdu3, hdu4):
        hdu.header.update(wcs.to_header())

    hdul = fits.HDUList([hdu1, hdu2, hdu3, hdu4, hdu5])

    outfile = '%s.fits' % outroot
    hdul.writeto(outfile, overwrite=True)
    print('Wrote results to %s' % outfile)
    x.close()


# ──────────────────────────────────────────────────────────────
# Command-line entry point
# ──────────────────────────────────────────────────────────────

if __name__ == '__main__':
    argv = sys.argv[1:]
    if not argv or '-h' in argv or '--help' in argv:
        print(_USAGE)
        sys.exit(0)

    idelta  = 1
    outroot = ''
    filename = None

    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == '-delta':
            i += 1
            idelta = int(argv[i])
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

    do_all(filename=filename, idelta=idelta, outroot=outroot)
