#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Interface between an LVM pointing (RA, Dec, time) and the PALACE
    airglow model, returning predicted line and continuum spectra split
    out by species (OH, O2, HO2, FeO, Na, K, O, N, H) as well as the
    combined spectrum PALACE would return on its own.

Command line usage (if any)::

    PalaceObs.py [-h] [-srf VALUE] [-species S1,S2,...] [-out ROOT] ra dec obstime

Arguments::

    ra          right ascension in degrees
    dec         declination in degrees
    obstime     UTC observation time (any format accepted by
                SkyModelObs.convert_time, e.g. '2023-08-29T03:20:43.668')

Options::

    -h                print this help and exit
    -srf VALUE        solar radio flux in sfu; default is looked up from
                      data/solar.txt via GetSolar.get_flux
    -species S1,S2    comma-separated species to predict individually;
                      default is all nine (OH, O2, HO2, FeO, Na, K, O, N, H)
    -out ROOT         output FITS filename root; default is PalaceObs_<ra>_<dec>_<mjd>

Output FITS structure::

    WAVE      wavelength [Angstrom], PALACE's own native grid truncated to
              the LVM range (see DEFAULT_LAMMIN/DEFAULT_LAMMAX) and
              generated at fine native resolution (see DEFAULT_DLAM/
              DEFAULT_RESOL)
    FLUX      total airglow flux, above the atmosphere [erg/s/cm^2/Angstrom]
    FLUX_LCO  total airglow flux, as observed on the ground at LCO
              [erg/s/cm^2/Angstrom]
    LINES, LINES_LCO
              PALACE's own total line emission (all species, lines only --
              see _call_model_linecont), above the atmosphere / at LCO
              [erg/s/cm^2/Angstrom] -- matches ESO's flux_ael
    DIFFUSE, DIFFUSE_LCO
              PALACE's own total continuum emission (all species,
              continuum only -- 3 components in palace_cont.fits: HO2,
              FeO, and O2's own separate continuum, VARID='O2Ac'), above
              the atmosphere / at LCO [erg/s/cm^2/Angstrom] -- matches
              ESO's flux_arc.  FLUX = LINES + DIFFUSE (to interpolation
              precision -- FLUX itself comes from PALACE's own separate
              'all' call, not by summing these two)
    OH, O2, HO2, FeO, NA, K, O, N, H
              per-species contribution to FLUX (above the atmosphere only)
              [erg/s/cm^2/Angstrom]

    The primary header records RA, DEC, OBSTIME, MSOLFLUX (the srf
    actually used), and ENGINE='palace' -- the same header convention as
    EsoSkyObs.py.  See Notes for the Rayleigh/nm -> erg/s/cm^2/Angstrom
    conversion and why there is no MOON/ZODI/CONT counterpart (PALACE
    models airglow only -- LINES+DIFFUSE, no scattered moonlight,
    zodiacal light, or telescope/instrument thermal emission).

Description:

    PALACE (Noll et al. 2025, "PALACE v1.0: Paranal Airglow Line And
    Continuum Emission model", Geoscientific Model Development, 18, 4353)
    predicts nighttime airglow line and continuum emission given zenith
    angle, local solar time, month, and solar radio flux.  It is not part
    of this repository -- see Readme.md for how to install it -- and its
    own top-level ``palace.model()`` function only ever returns a single
    combined spectrum: PALACE's internal tables keep every species
    separate (a ``chem``/``CHEMn`` tag per line and per continuum
    component), but ``calclinspec``/``calccontspec`` sum everything into
    one flux array before the table is ever returned.

    This module (a) computes the geometry/time parameters PALACE needs
    from an LVM observation (ra, dec, obstime), reusing
    ``SkyModelObs.get_info_las_campanas`` for the site geometry and
    ``GetSolar.get_flux`` for the solar radio flux, and (b) calls
    ``palace.model()`` once per species (plus once more for the combined
    spectrum) to recover the per-species breakdown that PALACE's own data
    supports but its public API collapses.

    This is a first cut built on PALACE's public, documented API only
    (``palace.parlist`` / ``palace.model``, looped over ``species=``), not
    its internal per-component tables -- simpler and more robust to future
    PALACE versions, at the cost of re-reading/re-scaling PALACE's (small)
    FITS tables once per species instead of once total.  If per-row
    runtime ever becomes a bottleneck, this can be replaced by a version
    that calls PALACE's internal ``readdata``/``calcscalfac``/etc.
    directly and only re-implements the final per-component summing step,
    without changing this module's return shape (see History).

Primary routines:

    predict     main entry point; returns (species_table, combined_table),
                PALACE's own native units (Rayleighs/nm) -- unchanged by
                the homogenization below, since SkySepPalace.py's per-row
                NNLS fit uses these as free-amplitude shapes and doesn't
                care about absolute units
    do_one      CLI-friendly wrapper; writes a homogenized FITS file (see
                Output FITS structure) in physical units, truncated to the
                LVM wavelength range by default

Notes:

    PALACE's climatology is calibrated for Cerro Paranal, not Las
    Campanas; the two sites are close (LCO -29.01,-70.69 vs Paranal
    -24.63,-70.40, both high-altitude Atacama sites) so this is a
    reasonable approximation for now, not an exact match.

    Local solar time (PALACE's ``tbin``) is approximated as
    UTC + longitude/15h, ignoring the equation of time (<= 16 min) --
    negligible next to PALACE's hour-wide time bins.

    predict()/get_palace_templates() (used internally by SkySepPalace.py's
    per-row NNLS fit) still return PALACE's raw Rayleighs/nm -- that fit
    uses each template as a free-amplitude shape, so absolute units don't
    matter there and are deliberately left alone (also avoids doubling the
    palace.model() call count in that per-row/per-fiber hot path).  Only
    do_one()'s output FITS (see below) is converted to physical units.

    do_one()'s output is homogenized with EsoSkyObs.py's convention:
    WAVE in Angstrom, FLUX/FLUX_LCO and per-species columns in
    erg/s/cm^2/Angstrom for one LVM fiber (FIBER_AREA_ARCSEC2, same value
    EsoSkyObs.py uses).  The conversion from PALACE's native Rayleighs/nm
    (see rayleigh_per_nm_to_flux) uses the standard aeronomy definition
    1 Rayleigh = 1e6/(4*pi) photons/cm^2/s/sr (apparent/observed airglow
    surface brightness -- see e.g. Hunten et al. 1956), the same photon
    energy constant as EsoSkyObs.finalize_table, and FIBER_AREA_ARCSEC2 to
    go from per-steradian to per-fiber.  Validated (260710) against the
    real ESO Sky Model's own LINES column at the same position/time:
    PALACE's converted airglow flux agrees with ESO's to within a factor
    of ~1.2-1.5 at the strong OI 5577/6300 lines (reasonable given these
    are independent climatological models for different sites); a few
    individual narrow lines (Na D, OH band heads) differ by an order of
    magnitude or more at a single wavelength bin, consistent with real
    airglow variability and/or the two models' native grids sampling a
    narrow, poorly-resolved line at slightly different points -- not a
    systematic unit error, since the bulk of lines and the overall median
    flux level agree in order of magnitude.

    PALACE's own isatm parameter (default True) already distinguishes the
    two atmospheric states EsoSkyObs.py's FLUX/FLUX_LCO need: isatm=False
    is the above-the-atmosphere emission (no extinction correction
    applied), isatm=True is the as-observed-at-the-telescope (ground)
    value.  So FLUX/FLUX_LCO here come from two separate palace.model()
    calls (isatm=False/True) rather than from one call divided by a
    transmission curve, unlike EsoSkyObs.py's own FLUX/FLUX_LCO split
    (which has to divide by ESO's own trans column, since ESO's API
    doesn't expose an isatm-style toggle per component).

    do_one()'s wavelength range now defaults to the LVM range (0.36-0.98
    micron), truncated from PALACE's own default (0.3-2.5 micron), via
    DEFAULT_LAMMIN/DEFAULT_LAMMAX -- pass lammin/lammax explicitly to
    override.  The wavelength grid points themselves are still PALACE's
    own native grid (dlam/resol -- see next paragraph); it is not
    resampled onto EsoSkyObs.py's exact grid points.

    do_one()'s dlam/resol default to DEFAULT_DLAM/DEFAULT_RESOL (fine
    native resolution: R~20000), not PALACE's own out-of-the-box defaults
    (dlam=1e-4 micron, resol=1e3 -> R~1000, FWHM~5.6A at 5577A).  Without
    this override, PalaceObs.py's output visibly looks lower-resolution
    than EsoSkyObs.py's when plotted side by side -- not a real physical
    difference, just PALACE's coarse default versus ESO's engines both
    natively running at R~5500-14000 (SkyModelObs's local engine:
    wgauss=0.8 px * 0.5 A/px = 0.4A FWHM, R~14000 at 5577A; SkyCalcObs's
    remote engine: lsf_gauss_fwhm=2.0 px * 0.5 A/px = 1.0A FWHM, R~5580).
    DEFAULT_RESOL=20000 matches the native resolution SkySepPalace.py
    already uses internally for its own PALACE templates (_NATIVE_RESOL);
    see that constant's comment for the same "generate fine, convolve to
    whatever's needed downstream" rationale.

    FLUX_LCO can come out very slightly *brighter* than FLUX (a few
    percent, seen at ~45% of wavelength points in one test case) -- this
    looks backwards for extinction but is not a bug: PALACE's own
    isatm=True scattering correction (calcscat, Noll et al. 2012) includes
    a van-Rhijn/in-scattering term for the extended airglow layer whose
    coefficient can go negative, letting the net scattering "transmission"
    slightly exceed 1 for some zenith angles/wavelengths -- i.e. light
    scattered into the line of sight from the rest of the sky can outweigh
    light scattered out of it, for an extended (not point) source.

    LINES vs LINES_LCO (and DIFFUSE vs DIFFUSE_LCO), and which one to
    compare against ESO's LINES/DIFFUSE: as of this writing,
    EsoSkyObs.finalize_table divides FLUX/MOON/ZODI/DIFFUSE/CONT by the
    transmission (above-the-atmosphere convention) but does NOT divide
    LINES -- inherited unchanged from the original SkyModelObs.py/
    SkyCalcObs.py code, and never actually re-examined against the same
    "airglow originates at ~90km and should see the same extinction as
    everything else on the way down" reasoning that FLUX got (see
    EsoSkyObs.py History).  That makes ESO's current LINES a ground-level
    quantity while its DIFFUSE is already above-the-atmosphere -- an
    inconsistency between the two ESO columns, not just between ESO and
    PALACE.  So today: PALACE's LINES_LCO matches ESO's LINES, but
    PALACE's DIFFUSE (not DIFFUSE_LCO) matches ESO's DIFFUSE.  If
    EsoSkyObs.py's LINES is ever also divided by trans for consistency
    with its own DIFFUSE, PALACE's LINES (isatm=False) would become the
    correct comparison instead.

    CONTINUUM_SPECIES/LINE_SPECIES (HO2, FeO vs the other 7) are NOT what
    LINES/DIFFUSE are built from -- an earlier version of this module did
    exactly that (sum whole species by category) and got it wrong: per
    palace_cont.fits' own header (NCONT=3, CHEM1=HO2, CHEM2=FeO, CHEM3=O2,
    VARID3='O2Ac'), O2 has its own separate continuum component in
    addition to its line emission, so classifying all of species='O2' as
    "line" silently folded that continuum into LINES.  LINES/DIFFUSE now
    come from _call_model_linecont, which calls PALACE's internal
    line-table/continuum-table pipeline directly and keeps them separate
    before any species-level regrouping.  CONTINUUM_SPECIES/LINE_SPECIES
    still exist because SkySepPalace.py needs a per-species (not per-line-
    vs-continuum) classification for its own NNLS-amplitude continuum-fit-
    quality columns, where each of the 9 species is one fitted template
    with one amplitude, not separable into a line part and a continuum
    part; SkySepPalace.py imports CONTINUUM_SPECIES from here instead of
    keeping its own separate copy of the HO2/FeO part of that list.

History:

    260709 ksl Coding begun.  First cut: loop the public palace.model()
        over each of PALACE's 9 species plus once for the combined
        spectrum.  Geometry/time reuses SkyModelObs.get_info_las_campanas;
        solar flux reuses GetSolar.get_flux.
    260710 ksl do_one() rewritten to homogenize the output FITS with
        EsoSkyObs.py: single table (WAVE, FLUX, FLUX_LCO, per-species
        columns) in physical units (erg/s/cm^2/Angstrom) instead of the
        native-unit SPECIES/COMBINED extensions, default wavelength range
        truncated to LVM's.  predict()/get_palace_templates() (used by
        SkySepPalace.py) are unchanged.
    260710 ksl Added LINES/LINES_LCO (aliases of FLUX/FLUX_LCO) for direct
        column-name comparison against EsoSkyObs.py's LINES.  Raised
        do_one()'s default dlam/resol from PALACE's own coarse defaults to
        DEFAULT_DLAM/DEFAULT_RESOL (matching SkySepPalace.py's own native
        resolution) -- without this, do_one()'s output was visibly lower
        resolution than EsoSkyObs.py's, not for any physical reason.
    260710 ksl LINES/LINES_LCO were wrongly aliasing the FULL FLUX/FLUX_LCO
        total, including HO2/FeO -- ESO's own DIFFUSE (flux_arc) is
        exactly HO2/FeO-like continuum, kept separate from its LINES
        (flux_ael).  Added CONTINUUM_SPECIES/LINE_SPECIES and redefined
        LINES/LINES_LCO as the LINE_SPECIES-only subset sum, with new
        DIFFUSE/DIFFUSE_LCO columns for the CONTINUUM_SPECIES-only subset
        sum, matching ESO's split.  Requires a full per-species breakdown
        at both isatm states now (previously the ground/isatm=True call
        only fetched the combined total via species_list=[]).
    260710 ksl That per-species split was still wrong: palace_cont.fits'
        own header shows O2 has its own separate continuum component
        (CHEM3='O2', VARID3='O2Ac') in addition to its lines, so
        classifying all of species='O2' as "line" was silently folding
        that continuum into LINES.  Added _call_model_linecont, which
        calls PALACE's internal line-table/continuum-table pipeline
        directly (readdata/calcscalfac/scalelines/scalecont/corratmlines/
        corratmcont/calclinspec/calccontspec/convolvelsf) and keeps the
        two separate instead of summing them (skipping addlinescont) --
        the true split, not a per-species approximation of it.
        LINES/DIFFUSE now come from this; the ground/isatm=True call went
        back to species_list=[] since the per-species ground breakdown it
        was fetching is no longer used for anything.

'''

import os
import sys
import io
import contextlib

import numpy as np
from astropy.table import Table, vstack
from astropy.io import fits

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from palace import palace

from SkyModelObs import get_info_las_campanas, convert_time
from GetSolar import get_flux as get_solar_flux
from EsoSkyObs import FIBER_AREA_ARCSEC2


_USAGE = '''Usage: PalaceObs.py [-h] [-srf VALUE] [-species S1,S2,...] [-out ROOT] ra dec obstime

Predict PALACE airglow spectra (per-species and combined) for an LVM
pointing and time.

Arguments:
    ra          right ascension in degrees
    dec         declination in degrees
    obstime     UTC observation time, e.g. 2023-08-29T03:20:43.668

Options:
    -h                print this help and exit
    -srf VALUE        solar radio flux in sfu (default: looked up from
                      data/solar.txt)
    -species S1,S2    comma-separated species (default: all nine)
    -out ROOT         output FITS filename root
'''


SPECIES = ['OH', 'O2', 'HO2', 'FeO', 'Na', 'K', 'O', 'N', 'H']

# Split of PALACE's 9 species into line-dominated vs continuum-like, for
# do_one()'s LINES/DIFFUSE columns -- matches ESO's own flux_ael (airglow
# emission Lines) / flux_arc (airglow/residual emission) split.  Also the
# classification SkySepPalace.py already uses for its own continuum-fit-
# quality columns (see that module's CONTINUUM_COMPONENTS, which reuses
# this list plus MOON).
CONTINUUM_SPECIES = ['HO2', 'FeO']
LINE_SPECIES = [sp for sp in SPECIES if sp not in CONTINUUM_SPECIES]

# Las Campanas Observatory longitude, matching SkyModelObs.get_info_las_campanas
LCO_LON_DEG = -70.6920

# do_one()'s default output wavelength range, truncated to the LVM range
# (unlike PALACE's own default of 0.3-2.5 micron); micron, matching PALACE's
# own lammin/lammax convention.
DEFAULT_LAMMIN = 0.36
DEFAULT_LAMMAX = 0.98

# do_one()'s default grid step / resolving power.  PALACE's own defaults
# (dlam=1e-4 micron, resol=1e3 -> FWHM~5.6A / R~1000 at 5577A) are much
# coarser than EsoSkyObs.py's engines (both natively R~5500-14000), so
# without overriding these the two outputs look very different in
# resolution when plotted side by side even though the underlying physics
# is comparable.  These defaults instead match the fine native resolution
# SkySepPalace.py already uses internally for its own PALACE templates
# (_NATIVE_RESOL/_NATIVE_DLAM) -- deliberately far exceeding any real
# instrumental LSF, so the output is effectively unconvolved/native and can
# be rebinned/convolved onto whatever comparison or instrumental LSF is
# needed downstream, same rationale as SkySepPalace.py's own choice.
DEFAULT_DLAM = 2e-5    # micron, ~0.2 Angstrom
DEFAULT_RESOL = 20000

# Rayleigh/nm -> erg/s/cm^2/Angstrom, for one LVM fiber (see module Notes).
# 1 Rayleigh = 1e6/(4*pi) photons/cm^2/s/sr (standard aeronomy definition of
# apparent/observed airglow surface brightness, e.g. Hunten et al. 1956).
_RAYLEIGH_TO_PH_CM2_S_SR = 1e6 / (4 * np.pi)
_ARCSEC2_TO_SR = (1.0 / 206264.806) ** 2
_HC_ERG_ANGSTROM = 1.98644586e-8   # hc, for E(lambda)=hc/lambda[Angstrom]
_RAYLEIGH_PER_NM_TO_ESO_CONST = _RAYLEIGH_TO_PH_CM2_S_SR * _ARCSEC2_TO_SR / 10.0 * _HC_ERG_ANGSTROM


def rayleigh_per_nm_to_flux(flux_r_per_nm, wave_ang):
    '''
    Convert a PALACE flux (Rayleighs/nm) to erg/s/cm^2/Angstrom for one LVM
    fiber (FIBER_AREA_ARCSEC2), the same physical convention as
    EsoSkyObs.finalize_table.
    '''
    return flux_r_per_nm * FIBER_AREA_ARCSEC2 * _RAYLEIGH_PER_NM_TO_ESO_CONST / wave_ang


def local_solar_time(obstime):
    '''
    Approximate local mean solar time in hours relative to local midnight
    (range -12 to 12), from UTC and the Las Campanas longitude.  Ignores
    the equation of time (<= 16 min), negligible next to PALACE's
    hour-wide tbin bins.
    '''

    dt = convert_time(obstime, 'datetime')
    ut_hours = dt.hour + dt.minute / 60. + dt.second / 3600.
    lst = (ut_hours + LCO_LON_DEG / 15.) % 24.
    if lst >= 12.:
        lst -= 24.
    return lst


def get_geometry(ra, dec, obstime):
    '''
    Zenith angle (deg), month (1-12), and approximate local solar time
    (hours, see local_solar_time) for an LVM pointing, for use as
    PALACE's z/mbin/tbin parameters.
    '''

    info = get_info_las_campanas(obstime, ra=ra, dec=dec)
    z = 90. - info['SourceAlt']
    month = convert_time(obstime, 'datetime').month
    tbin = local_solar_time(obstime)
    return z, month, tbin


def build_cline(ra, dec, obstime, srf=None, species='all', **overrides):
    '''
    Build a PALACE command-line-style parameter list (see palace.parlist)
    for an LVM pointing/time.
    '''

    z, mbin, tbin = get_geometry(ra, dec, obstime)
    if srf is None:
        srf = float(get_solar_flux(obstime))
    cline = ['z=%.6f' % z, 'mbin=%d' % mbin, 'tbin=%.6f' % tbin,
             'srf=%.4f' % srf, 'species=%s' % species]
    for key, val in overrides.items():
        cline.append('%s=%s' % (key, val))
    return cline


def _call_model(cline):
    '''
    Call palace.model(), working around a real upstream bug: PALACE's
    calclinspec() and calccontspec() use different floating-point
    tolerances (1e-7 vs 1e-8) when building their output wavelength
    grids from the same lammin/lammax/dlam, so for species with both a
    line and a continuum component (only 'O2' among PALACE's nine),
    the two grids can differ by one point for particular (lammin,
    lammax, dlam) combinations.  When that happens, palace.addlinescont
    prints "Line and continuum spectra: different wavelength grids!"
    and silently returns an all-zero spectrum instead of raising an
    error.

    This is detected here (by capturing stdout) and worked around by
    retrying once with dlam nudged by a factor of (1 + 1e-5) -- enough
    to move off the exact floating-point coincidence, negligible for
    the physics.  If the retry also triggers it, the (all-zero) result
    is returned as-is and a warning is printed, since falling back
    further would risk masking a real problem instead.
    '''

    parlist = palace.parlist(cline)
    if len(parlist) == 0:
        return None

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        spec = palace.model(**parlist)
    if 'different wavelength grids' in buf.getvalue():
        retry_cline = list(cline) + ['dlam=%.10g' % (parlist['dlam'] * (1 + 1e-5))]
        retry_parlist = palace.parlist(retry_cline)
        if len(retry_parlist):
            buf2 = io.StringIO()
            with contextlib.redirect_stdout(buf2):
                spec = palace.model(**retry_parlist)
            if 'different wavelength grids' in buf2.getvalue():
                print('Warning: PALACE line/continuum grid mismatch persisted '
                      'after dlam retry for species=%s; result may be all-zero. '
                      'cline=%s' % (parlist.get('species'), cline))
    return spec


def _call_model_linecont(cline):
    '''
    Like _call_model, but returns PALACE's line-only and continuum-only
    spectra SEPARATELY instead of already summed (palace.model()'s
    addlinescont step).  Needed for a correct LINES/DIFFUSE split: PALACE's
    public species= selector operates on whole species, and (per
    palace_cont.fits' own header: NCONT=3, CHEM1=HO2, CHEM2=FeO, CHEM3=O2,
    VARID3='O2Ac') 'O2' has its own separate continuum component in
    addition to its line emission -- there is no public species= value
    that means "all line emission" or "all continuum emission" alone, so
    per-species classification (SPECIES minus CONTINUUM_SPECIES) cannot
    correctly separate O2's continuum from its lines.  This instead
    replicates palace.model()'s own pipeline (readdata -> calcscalfac ->
    scalelines/scalecont -> corratmlines/corratmcont -> calclinspec/
    calccontspec -> convolvelsf) using PALACE's internal functions
    directly, stopping short of the final addlinescont sum -- exactly the
    "internal per-component tables" route PalaceObs.py's own module
    docstring anticipated as a possible future direction.  Only ever
    called with species='all', since do_one() wants PALACE's total line
    emission and total continuum emission, not a per-species breakdown of
    each.

    Returns (linspec, contspec), each an astropy table with lam/flux/dflux
    (empty Table() if that component had no data), or (None, None) if the
    parameter list itself was invalid.
    '''

    parlist = palace.parlist(cline)
    if len(parlist) == 0:
        return None, None

    lintab, conttab, vartab = palace.readdata(**parlist)
    if len(vartab) == 0:
        return None, None

    scalfac = palace.calcscalfac(vartab, **parlist)
    slintab = palace.scalelines(lintab, scalfac, **parlist)
    sconttab = palace.scalecont(conttab, scalfac, **parlist)
    cslintab = palace.corratmlines(slintab, **parlist)
    csconttab = palace.corratmcont(sconttab, **parlist)
    linspec = palace.convolvelsf(palace.calclinspec(cslintab, **parlist), **parlist)
    contspec = palace.convolvelsf(palace.calccontspec(csconttab, **parlist), **parlist)
    return linspec, contspec


def predict(ra, dec, obstime, srf=None, species_list=SPECIES,
            want_combined=True, **overrides):
    '''
    Predict PALACE airglow spectra for an LVM pointing/time, split by
    species.

    Parameters:
        ra, dec       degrees
        obstime       UTC time, any SkyModelObs.convert_time format
        srf           solar radio flux in sfu; looked up from
                      data/solar.txt via GetSolar.get_flux if not given
        species_list  which of PALACE's nine species to predict
                      individually (default: all nine)
        want_combined also compute the single combined ('all') spectrum
        overrides     additional PALACE parameters (e.g. isatm=False,
                      lammin=0.55, resol=4000) passed straight through to
                      every palace.model() call

    Returns:
        species_table  astropy Table with columns lam, flux, dflux,
                       species (species with no contribution in the
                       requested wavelength range are simply omitted)
        combined       astropy Table as returned by
                       palace.model(species='all', ...), or None if
                       want_combined is False
    '''

    if srf is None:
        srf = float(get_solar_flux(obstime))

    rows = []
    for sp in species_list:
        cline = build_cline(ra, dec, obstime, srf=srf, species=sp,
                            **overrides)
        spec = _call_model(cline)
        if spec is None:
            raise ValueError(
                'Invalid PALACE parameter list for species %s '
                '(ra=%s dec=%s obstime=%s)' % (sp, ra, dec, obstime))
        if len(spec) == 0:
            continue
        spec['species'] = sp
        rows.append(spec)

    species_table = vstack(rows) if rows else Table()

    combined = None
    if want_combined:
        cline = build_cline(ra, dec, obstime, srf=srf, species='all',
                            **overrides)
        combined = _call_model(cline)
        if combined is None:
            raise ValueError(
                'Invalid PALACE parameter list '
                '(ra=%s dec=%s obstime=%s)' % (ra, dec, obstime))

    return species_table, combined


def do_one(ra, dec, obstime, srf=None, species_list=SPECIES, outroot='',
          lammin=None, lammax=None, dlam=None, resol=None, **overrides):
    '''
    CLI-friendly wrapper around predict(): writes a FITS file homogenized
    with EsoSkyObs.py's convention (WAVE in Angstrom; FLUX, FLUX_LCO, and
    one column per species in erg/s/cm^2/Angstrom -- see module docstring
    Output FITS structure and Notes) and returns the output filename.

    lammin/lammax (PALACE's own micron units) default to the LVM range
    (DEFAULT_LAMMIN/DEFAULT_LAMMAX), truncated from PALACE's own default
    of 0.3-2.5 micron.  dlam/resol default to DEFAULT_DLAM/DEFAULT_RESOL
    (fine native resolution), not PALACE's own coarse defaults -- see
    those constants' comments for why.
    '''

    if lammin is None:
        lammin = DEFAULT_LAMMIN
    if lammax is None:
        lammax = DEFAULT_LAMMAX
    if dlam is None:
        dlam = DEFAULT_DLAM
    if resol is None:
        resol = DEFAULT_RESOL

    if srf is None:
        srf = float(get_solar_flux(obstime))

    # Above the atmosphere (isatm=False): the primary convention, matching
    # EsoSkyObs's FLUX/MOON/ZODI/DIFFUSE/CONT.  Full per-species breakdown
    # is needed here (for the OH/O2/... columns).
    species_above, combined_above = predict(
        ra, dec, obstime, srf=srf, species_list=species_list,
        want_combined=True, lammin=lammin, lammax=lammax, dlam=dlam,
        resol=resol, isatm=False, **overrides)

    # As actually observed on the ground at LCO (isatm=True): only the
    # combined total is needed here (for FLUX_LCO) -- LINES_LCO/DIFFUSE_LCO
    # come from linecont_columns() below instead of a per-species
    # breakdown, so species_list=[] skips the per-species calls here.
    _unused, combined_ground = predict(
        ra, dec, obstime, srf=srf, species_list=[],
        want_combined=True, lammin=lammin, lammax=lammax, dlam=dlam,
        resol=resol, isatm=True, **overrides)

    wave_ang = np.array(combined_above['lam'], dtype=float) * 1e4
    flux = rayleigh_per_nm_to_flux(np.array(combined_above['flux'], dtype=float), wave_ang)
    flux_lco = rayleigh_per_nm_to_flux(np.array(combined_ground['flux'], dtype=float), wave_ang)

    def species_columns(species_table):
        cols = {}
        for sp in species_list:
            sub = species_table[species_table['species'] == sp] if len(species_table) else species_table
            if len(sub) == 0:
                cols[sp] = np.zeros_like(wave_ang)
            else:
                sub_wave_ang = np.array(sub['lam'], dtype=float) * 1e4
                sub_flux = rayleigh_per_nm_to_flux(np.array(sub['flux'], dtype=float), sub_wave_ang)
                # Interpolated rather than assigned directly: PALACE's
                # grid-mismatch workaround (see _call_model) can nudge
                # dlam very slightly for an individual species, so its own
                # grid isn't always byte-identical to combined_above's.
                cols[sp] = np.interp(wave_ang, sub_wave_ang, sub_flux, left=0.0, right=0.0)
        return cols

    above_cols = species_columns(species_above)

    def linecont_columns(isatm):
        '''
        PALACE's true total line-only and continuum-only spectra (species=
        'all'), converted and interpolated onto wave_ang -- see
        _call_model_linecont for why this can't be built from a per-species
        sum (O2 has its own separate continuum component alongside its
        lines, which per-species classification can't pull apart).
        '''
        cline = build_cline(ra, dec, obstime, srf=srf, species='all',
                            lammin=lammin, lammax=lammax, dlam=dlam,
                            resol=resol, isatm=isatm, **overrides)
        linspec, contspec = _call_model_linecont(cline)
        if linspec is None:
            raise ValueError(
                'Invalid PALACE parameter list '
                '(ra=%s dec=%s obstime=%s)' % (ra, dec, obstime))

        def to_flux(spec):
            if len(spec) == 0:
                return np.zeros_like(wave_ang)
            spec_wave_ang = np.array(spec['lam'], dtype=float) * 1e4
            spec_flux = rayleigh_per_nm_to_flux(np.array(spec['flux'], dtype=float), spec_wave_ang)
            return np.interp(wave_ang, spec_wave_ang, spec_flux, left=0.0, right=0.0)

        return to_flux(linspec), to_flux(contspec)

    out = Table()
    out['WAVE'] = wave_ang
    out['FLUX'] = flux
    out['FLUX_LCO'] = flux_lco
    # LINES/DIFFUSE come from PALACE's own line-table/continuum-table split
    # (see linecont_columns/_call_model_linecont), not from summing whole
    # species -- O2's own separate continuum component would otherwise end
    # up inside LINES.  NOTE: EsoSkyObs.py's own LINES column is currently
    # NOT divided by the transmission (unlike its FLUX/MOON/ZODI/DIFFUSE/
    # CONT), i.e. it is still in the as-observed-at-ground convention --
    # so today, LINES_LCO (not LINES) is the one that matches ESO's LINES
    # like-for-like, while DIFFUSE (not DIFFUSE_LCO) matches ESO's DIFFUSE.
    # See History.
    out['LINES'], out['DIFFUSE'] = linecont_columns(isatm=False)
    out['LINES_LCO'], out['DIFFUSE_LCO'] = linecont_columns(isatm=True)

    for sp in species_list:
        out[sp.upper()] = above_cols[sp]

    if outroot == '':
        mjd = convert_time(obstime, 'mjd')
        outroot = 'PalaceObs_%.5f_%.5f_%08.2f' % (ra, dec, mjd)
    outname = outroot if outroot.endswith('.fits') else outroot + '.fits'

    primary = fits.PrimaryHDU()
    primary.header['RA'] = ra
    primary.header['DEC'] = dec
    primary.header['OBSTIME'] = str(obstime)
    primary.header['MSOLFLUX'] = srf
    primary.header['ENGINE'] = 'palace'
    hdul = fits.HDUList([primary, fits.BinTableHDU(data=out)])
    hdul.writeto(outname, overwrite=True)
    return outname


def is_number(x):
    try:
        eval(x)
        return True
    except Exception:
        return False


def steer(argv):
    '''
    PalaceObs.py [-h] [-srf VALUE] [-species S1,S2,...] [-out ROOT] ra dec obstime
    '''

    ra = None
    dec = None
    obstime = ''
    srf = None
    species_list = SPECIES
    outroot = ''

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-srf':
            i += 1
            srf = float(argv[i])
        elif argv[i] == '-species':
            i += 1
            species_list = argv[i].split(',')
        elif argv[i][:4] == '-out':
            i += 1
            outroot = argv[i]
        elif argv[i][0] == '-' and not is_number(argv[i]):
            print('Error: unknown option:', argv[i])
            print(_USAGE)
            return
        elif ra is None:
            ra = eval(argv[i])
        elif dec is None:
            dec = eval(argv[i])
        elif obstime == '':
            obstime = argv[i]
        else:
            print('Error: too many arguments:', argv)
            print(_USAGE)
            return
        i += 1

    if ra is None or dec is None or obstime == '':
        print(_USAGE)
        return

    outname = do_one(ra, dec, obstime, srf=srf, species_list=species_list,
                     outroot=outroot)
    print('Wrote', outname)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
