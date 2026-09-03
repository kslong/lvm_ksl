#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Generate a predicted sky spectrum for a given RA, Dec, and time using the
real ESO Sky Model.

This routine unifies two previously separate approaches:

- the local ESO Sky Model calculator (calcskymodel), previously run from
  SkyModelObs.py (retired 260711 -- fully absorbed into this module)
- the ESO web-service SkyCalc CLI (skycalc_cli), previously run from
  SkyCalcObs.py (retired 260711 -- fully absorbed into this module)

By default it tries the local install first (it is the version we maintain
and trust) and falls back to the ESO web service only if the local model
is not set up.

Command line usage::

    EsoSkyObs.py [-h] [-engine local|remote|auto] [-msol flux] [-out root] [-site lco|paranal] [-pres hPa] [-keep_workdir] ra dec time

Arguments: ra and dec are positions in the sky in degrees; time is in one
of several formats (date_time string, MJD, or JD).

Options::

    -h                print this help and exit
    -engine E         'local'  use calcskymodel only; error if not available
                      'remote' always use the ESO SkyCalc web service
                      'auto'   try local first, fall back to remote if the
                               local model is not available (default)
    -msol flux        force a specific 10.7 cm solar radio flux instead of
                      looking one up from the historical flux table
    -out root         root name for the output file (default: derived from
                      the position and time)
    -site S           'lco' (default) or 'paranal' -- see Notes on why this
                      is an approximation, not a full site swap
    -pres hPa         override the site's default pressure (lco: 765,
                      paranal: 744) -- see SITE_PRESSURE_HPA; local engine
                      only, ignored by -engine remote
    -keep_workdir     debugging switch: by default, calcskymodel's inputs/
                      outputs (config/*.par, output/*.fits for the local
                      engine; *.json/*.fits scratch files for the remote
                      engine) live in a fresh, auto-deleted temp directory
                      per call -- concurrency-safe, but nothing survives
                      to inspect.  -keep_workdir instead writes them
                      directly into the current directory (./config,
                      ./data as a symlink, ./output) and leaves them there
                      -- the pre-260730 behavior, and where calcskymodel
                      itself looks if you cd there and run it by hand.
                      NOT concurrency-safe -- only use this for one call
                      at a time.

Description:

    Both the local ESO Sky Model and the ESO SkyCalc web service predict the
    sky spectrum as it would actually be observed on the ground, i.e. after
    atmospheric extinction.  LVM spectra are normally compared against the
    above-the-atmosphere equivalent, so FLUX (and the MOON/ZODI/DIFFUSE/CONT
    components) here are corrected for extinction by dividing by the
    transmission.  The value as actually observed on the ground at LCO is
    preserved separately in FLUX_LCO.

Output FITS structure::

    WAVE      wavelength [Angstrom]
    FLUX      total sky flux, corrected for atmospheric extinction
    FLUX_LCO  total sky flux, as actually observed on the ground at LCO
    MOON      scattered moonlight component (extinction-corrected)
    ZODI      zodiacal light component (extinction-corrected)
    LINES     airglow emission line component (as modeled, uncorrected)
    DIFFUSE   diffuse/thermal component (extinction-corrected)
    CONT      MOON + ZODI + DIFFUSE
    trans     atmospheric transmission

    The primary header records RA, DEC, OBSTIME, MSOLFLUX, and ENGINE (which
    of the two engines actually produced the file).

Primary routines::

    run_sky_obs   top-level routine; tries local, falls back to remote
    run_local     run the local calcskymodel binary
    run_remote    run the ESO SkyCalc web service via skycalc_cli

Notes:

    The local engine requires the environment variable ESO_SKY_MODEL to
    point to a working installation of the real ESO Sky Model.  The remote
    engine requires skycalc_cli to be pip-installed, a working
    setuptools/pkg_resources in the environment, and network access to
    eso.org.

    site='paranal' (default 'lco'): for comparing against PALACE, whose
    own atmospheric-physics constants are hardcoded to Cerro Paranal (see
    PalaceObs.py's Notes -- h=2.64 km, p=744 hPa, not overridable via any
    public PALACE parameter), not adjustable to LCO or any other site.
    The two engines handle this differently::

      local engine: only SITE_HEIGHT_KM's sm_h (observatory height/
        pressure physics) changes; the real observing geometry (alt/az,
        moon phase/separation) still comes from get_info_las_campanas's
        actual LCO coordinates -- the same mixed real-geometry/Paranal-
        physics approach PALACE itself uses, so this is the more
        comparable of the two engines.
      remote engine: skycalc_cli's own "observatory" parameter drives
        BOTH the atmosphere physics AND its internal moon/sun almanac
        geometry calculation (see REMOTE_SITE_NAME) -- there is no way to
        keep LCO's real geometry while asking for Paranal's atmosphere
        through the public API, so site='paranal' here also shifts the
        modeled sky to Paranal's real geographic location, not just its
        altitude.  A real difference from what PALACE (and this module's
        own local engine) do, not merely a units/precision detail.

History::

    260710 ksl Coding begun; unifies SkyCalcObs.py and SkyModelObs.py
    260710 ksl Added site='lco'|'paranal' (SITE_HEIGHT_KM/REMOTE_SITE_NAME) to
    let a comparison run approximate PALACE's Paranal-fixed atmosphere
    physics; default behavior (site='lco') unchanged.  Not a full site
    swap for the remote engine -- see Notes.
    260711 ksl Retired SkyModelObs.py and SkyCalcObs.py (nothing else used
    SkyCalcObs.py; SkyModelObs.py was still imported by this module for
    get_info_las_campanas/setup, and directly by SkySepESO.py and
    MakeMoonBase.py).  get_info_las_campanas/setup/safe_remove moved here
    verbatim; SkySepESO.py and MakeMoonBase.py migrated to call
    run_sky_obs()/this module instead of the retired scripts directly (the
    former also gains a real fix along the way: its local engine now
    resolves the historical solar flux via resolve_solar_flux() instead of
    SkyModelObs.py's old hardcoded msolflux=101).
    260711 ksl Real incident: a model fetch run with this repo's own root as
    the working directory (instead of a scratch directory) triggered
    setup()'s old unconditional safe_remove('data') -> shutil.rmtree(),
    destroying this repo's real data/ directory (528MB of vendored PALACE
    reference data, sky_mask.fits, etc.) before symlinking over it.  Fully
    recovered via git checkout (everything was committed), but that was
    luck, not safety.  Fixed: safe_remove() no longer ever removes a real
    (non-symlink) directory, only symlinks; setup() now raises RuntimeError
    if 'data' already exists as a real directory rather than silently
    destroying it.
    260801 ksl Added SITE_PRESSURE_HPA and -pres: the local engine's config
    pres value was previously hardcoded to 744 hPa (the model's own
    built-in default, and PALACE's own fixed value) regardless of site.
    site='lco' now defaults pres to 765 hPa (nominal LCO barometric
    pressure) instead, via the new SITE_PRESSURE_HPA dict (paranal
    unchanged at 744); -pres/create_local_inputs's pressure argument
    overrides either default.  Local engine only -- skycalc_cli exposes no
    separate pressure parameter.
    260801 ksl Added -keep_workdir as a debugging switch: run_local()/
    run_remote() normally write calcskymodel's/skycalc_cli's inputs and
    outputs into a fresh, auto-deleted tempfile.TemporaryDirectory() per
    call (via the new _scratch_dir() helper) -- concurrency-safe, but
    nothing survives a run to inspect.  -keep_workdir instead restores
    the pre-260730 behavior: config/*.par, the data symlink, and
    output/*.fits are written directly into the CALLER's current
    directory and left there, matching exactly where calcskymodel itself
    looks if you cd there and run it by hand -- confirmed by the binary
    itself, which errors with the literal message "File/dir does not
    exist: config/" when that directory isn't present relative to its
    cwd (calcskymodel hardcodes 'config', 'data', and 'output' internally
    and offers no way to rename any of them).  NOT concurrency-safe --
    only use it for one call at a time.
    260903 ksl Added lsf_gauss_fwhm_pix (create_local_inputs/run_local/
    run_sky_obs, local engine only), exposing inst_base's previously
    hardcoded wgauss=0.8 (calcskymodel's own internal LSF convolution,
    ~0.4 A FWHM at the default 0.5 A/pixel grid) as an explicit
    parameter.  Default unchanged, so every existing caller (including
    this script's own CLI) is unaffected.  Motivated by
    BatchPredictSkyESO.py applying its own additional per-row LSF
    convolution downstream -- without this, that stacked with
    calcskymodel's fixed internal kernel rather than replacing it,
    biasing the effective output FWHM slightly too broad.  See
    create_local_inputs()'s docstring for why this isn't just set to a
    much smaller value instead (aliasing risk against the coarse output
    grid, not merely a conservative default).

'''

import os
import sys
import json
import time
import subprocess
import tempfile
import contextlib
import warnings

import numpy as np
from astropy.io import fits
from astropy.table import Table, join
from astropy.time import Time
from astropy.coordinates import (get_body, solar_system_ephemeris, AltAz,
                                 EarthLocation, SkyCoord, GeocentricTrueEcliptic)
import astropy.units as u

from lvm_ksl.GetSolar import convert_time, get_flux

# get_info_las_campanas's cross-frame separation/transform_to calls trigger
# this astropy warning; NonRotationTransformationWarning was removed from
# recent astropy versions but drp still uses it, hence the fallback -- moved
# here (260711) from the now-retired SkyModelObs.py, unchanged.
try:
    from astropy.coordinates.baseframe import NonRotationTransformationWarning
except ImportError:
    from astropy.utils.exceptions import AstropyWarning as NonRotationTransformationWarning
warnings.simplefilter('ignore', NonRotationTransformationWarning)


# Angular area of one LVM fiber, in arcsec^2 -- used to convert the sky
# model's per-solid-angle native flux into the flux collected by one fiber.
# Shared with PalaceObs.py so PALACE-derived fluxes use the same convention.
FIBER_AREA_ARCSEC2 = np.pi * (37 / 2) ** 2


@contextlib.contextmanager
def _scratch_dir(prefix, workdir=None):
    '''
    Yield a per-call scratch directory for run_local()/run_remote().

    workdir=None (default): a fresh tempfile.TemporaryDirectory(),
    removed automatically on exit even on failure/exception -- the
    module's usual concurrency-safe behavior (see module Notes): each
    call gets its own randomly-named directory, so there's nothing
    shared for concurrent calls to race on.

    workdir=<path> (-keep_workdir uses '.'; -workdir uses whatever path
    was given): that directory is used directly instead (created if it
    does not exist) and is NEVER removed -- for a single interactive/
    manual run where you want config/*.par, the data symlink, and
    output/*.fits somewhere predictable you can inspect, hand-edit, and
    rerun calcskymodel in yourself.  NOT concurrency-safe -- do not point
    two simultaneous calls at the same path.
    '''
    if workdir:
        os.makedirs(workdir, exist_ok=True)
        print('Using working directory: %s' % workdir)
        yield workdir
    else:
        with tempfile.TemporaryDirectory(prefix=prefix) as tmp:
            yield tmp


def safe_remove(path):
    '''
    Remove a symlink at path, if one exists.  Deliberately does NOT
    remove a real (non-symlink) directory, even if one is found there --
    see setup()'s docstring for why.
    '''
    if os.path.islink(path):
        os.unlink(path)


def setup(eso_sky_dir='', workdir='.', config=True):
    '''
    Set up the directories (config/, output/, a data/ symlink) that the
    local calcskymodel binary needs, inside workdir (default '.', the
    current working directory -- but run_local() always passes an
    isolated per-call temporary directory instead, unless -keep_workdir;
    see its own docstring and module Notes on why).  Moved here (260711)
    from the now-retired SkyModelObs.py.

    config/, data, and output/ are always named exactly that -- not
    renameable -- because calcskymodel itself hardcodes those names
    (confirmed from the binary: it errors with the literal message "File/
    dir does not exist: config/" when run without them, and filepath=data,
    kernelfile=output/kernel.dat in inst_base hardcode the other two) and
    looks for them relative to its own cwd, wherever calcskymodel is run
    from.

    Refuses to touch '<workdir>/data' if it already exists as a real
    directory (raises RuntimeError) rather than deleting it -- this used
    to call shutil.rmtree() unconditionally on 'data' (via safe_remove),
    which on 260711 destroyed this repo's own real data/ directory
    (528MB of vendored PALACE reference data, sky_mask.fits, etc.) when
    a model fetch was accidentally run with the repo root as the working
    directory instead of a scratch directory -- fully recovered via git
    checkout since everything was committed, but that was luck, not
    safety.  'data' should only ever be a symlink this function itself
    created; a real directory there means setup() is running somewhere
    it shouldn't.  Kept as an unconditional guard even though run_local's
    fresh-per-call temporary directory can no longer trigger it at all
    -- cheap, and defends any future/direct caller that passes a
    workdir of its own.
    '''
    data_path   = os.path.join(workdir, 'data')
    output_path = os.path.join(workdir, 'output')
    config_path = os.path.join(workdir, 'config')

    if config == False:
        icheck = True
        if os.path.isdir(config_path) == False:
            icheck = False
        if os.path.isdir(output_path) == False:
            icheck = False
        if os.path.isdir(data_path) == False and os.path.islink(data_path) == False:
            icheck = False
        if icheck == True:
            return

    if os.path.isdir(data_path) and not os.path.islink(data_path):
        raise RuntimeError(
            "setup(): '%s' already exists as a real directory, not "
            "a symlink -- refusing to remove it. This usually means the "
            "local ESO Sky Model engine is being run from the wrong working "
            "directory (it needs a scratch directory of its own, not one "
            "with real data already in it). cd to a scratch directory, or "
            "remove/rename this data/ directory yourself if you are certain "
            "it is not needed." % data_path)

    xdir = os.getenv('ESO_SKY_MODEL')
    if eso_sky_dir == '':
        eso_sky_dir = xdir

    data_dir = '%s/sm-01_mod2/data' % eso_sky_dir
    if os.path.isdir(data_dir) == False:
        print('Error: %s really does not appear to exist' % data_dir)
        return
    safe_remove(data_path)
    os.symlink(data_dir, data_path)
    os.makedirs(output_path, exist_ok=True)
    os.makedirs(config_path, exist_ok=True)
    return


def get_info_las_campanas(datetime_utc, ra, dec, verbose=False):
    '''
    Get information about the sun, moon, and a source at given RA and Dec
    as a function of UT, from Las Campanas Observatory.  Moved here
    (260711) from the now-retired SkyModelObs.py, unchanged.

    Parameters:
    -----------
    datetime_utc : str or datetime
        UTC time for the observation
    ra : float
        Right ascension of the source in degrees
    dec : float
        Declination of the source in degrees
    verbose : bool, optional
        If True, print the information

    Returns:
    --------
    dict
        Dictionary containing information about the sun, moon, and source
    '''
    if verbose:
        print('get_info_las_campanas,Start: ', datetime_utc, ra, dec)
    # Las Campanas Observatory coordinates
    observatory_location = EarthLocation(lat=-29.0089*u.deg, lon=-70.6920*u.deg, height=2281*u.m)

    obs_time = Time(datetime_utc)

    source_coords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='fk5')

    with solar_system_ephemeris.set('builtin'):
        moon_coords = get_body('moon', obs_time, location=observatory_location)
        sun_coords = get_body('sun', obs_time, location=observatory_location)

    phase_angle = moon_coords.separation(sun_coords).radian
    illumination_fraction = (1 - np.cos(phase_angle))/2

    moon_sun_longitude_diff = (moon_coords.ra - sun_coords.ra).wrap_at(360 * u.deg).value
    if moon_sun_longitude_diff > 0:
        moon_phase = illumination_fraction/2.
    else:
        moon_phase = 1-illumination_fraction/2.
    illumination_fraction *= 100.

    altaz_frame = AltAz(obstime=obs_time, location=observatory_location)
    moon_altaz = moon_coords.transform_to(altaz_frame)
    sun_altaz = sun_coords.transform_to(altaz_frame)
    source_altaz = source_coords.transform_to(altaz_frame)
    if source_altaz.alt.deg < 0:
        print('Error: Source altitude is negative :', source_altaz.alt.deg, ra, dec, datetime_utc)

    moon_ecliptic = moon_coords.transform_to(GeocentricTrueEcliptic(equinox=obs_time))
    sun_ecliptic = sun_coords.transform_to(GeocentricTrueEcliptic(equinox=obs_time))
    source_ecliptic = source_coords.transform_to(GeocentricTrueEcliptic(equinox=obs_time))

    moon_eclip_lon = moon_ecliptic.lon.deg
    if moon_eclip_lon > 180:
        moon_eclip_lon -= 360

    if verbose:
        print('XXX %.1f %.1f -> %.1f ' % (source_ecliptic.lon.deg, sun_ecliptic.lon.deg,
                                          source_ecliptic.lon.deg-sun_ecliptic.lon.deg))

    sun_eclip_lon = sun_ecliptic.lon.deg
    source_eclip_lon = source_ecliptic.lon.deg

    mean_moon_distance = 384400 * u.km
    moon_distance = moon_coords.distance.to(u.km)
    moon_distance_in_mean = moon_distance / mean_moon_distance

    moon_sun_separation = moon_coords.separation(sun_coords).deg
    moon_source_separation = moon_coords.separation(source_coords).deg
    sun_source_separation = sun_coords.separation(source_coords).deg

    xreturn = {
        'SunRA': sun_coords.ra.deg,
        'SunDec': sun_coords.dec.deg,
        'SunAlt': sun_altaz.alt.deg,
        'SunAz': sun_altaz.az.deg,
        'SunEclipLon': sun_eclip_lon,
        'SunEclipLat': sun_ecliptic.lat.deg,
        'MoonRA': moon_coords.ra.deg,
        'MoonDec': moon_coords.dec.deg,
        'MoonAlt': moon_altaz.alt.deg,
        'MoonAz': moon_altaz.az.deg,
        'MoonEclipLon': moon_eclip_lon,
        'MoonEclipLat': moon_ecliptic.lat.deg,
        'MoonPhas': moon_phase,
        'MoonIll': illumination_fraction,
        'MoonDistance': moon_distance.value,
        'MoonDistanceInMeanUnits': moon_distance_in_mean.value,
        'SourceRA': source_coords.ra.deg,
        'SourceDec': source_coords.dec.deg,
        'SourceAlt': source_altaz.alt.deg,
        'SourceAz': source_altaz.az.deg,
        'SourceEclipLon': source_eclip_lon,
        'SourceEclipLat': source_ecliptic.lat.deg,
        'Moon-Sun_Separation': moon_sun_separation,
        'Moon-Source_Separation': moon_source_separation,
        'Sun-Source_Separation': sun_source_separation
    }

    if verbose:
        for key, value in xreturn.items():
            print(f'{key}: {value}')

    return xreturn


_USAGE = '''
Usage:  EsoSkyObs.py [-h] [-engine local|remote|auto] [-msol flux] [-out root] [-site lco|paranal] [-pres hPa] [-keep_workdir] ra dec time

    ra, dec    position in decimal degrees
    time       observation time: ISO string, MJD, or JD

    -engine local|remote|auto   which sky model engine to use (default: auto)
    -msol flux                  force a solar radio flux instead of looking one up
    -out root                   output file root name
    -site lco|paranal           observatory height/pressure physics (default: lco);
                                 see module Notes -- not a full site swap for -engine remote
    -pres hPa                    override the site's default pressure (lco: 765, paranal: 744);
                                  local engine only, see SITE_PRESSURE_HPA
    -keep_workdir                debugging: write directly into ./config, ./data (symlink),
                                  ./output in the current directory -- where calcskymodel itself
                                  looks if you cd there and run it by hand -- instead of the
                                  default per-call temp directory, and leave them there so the
                                  inputs can be inspected/edited and calcskymodel rerun directly.
                                  NOT concurrency-safe: only use this for one call at a time
'''


# ---------------------------------------------------------------------------
# Local engine (calcskymodel) inputs
# ---------------------------------------------------------------------------

inst_base = '''
# Wavelength grid:

# minimum and maximum wavelength [mum]
limlam     = 0.36 0.98

# step size [mum]
dlam       = 0.00005


# Line-spread function:

# radius of convolution kernel [pixels] (N_pixel = 2 x kernrad + 1)
kernrad    = 3

# FWHM of boxcar kernel [pixels]
wbox       = 0.8

# FWHM of Gaussian kernel [pixels] -- substituted by create_local_inputs's
# lsf_gauss_fwhm_pix (default 0.8, i.e. this literal value, so every
# existing caller that doesn't pass it sees unchanged behavior). This is
# calcskymodel's own internal LSF convolution, applied before the FITS
# output is ever seen by anything in lvm_ksl -- unrelated to (and, for
# callers doing their own additional LSF convolution downstream, stacks
# with) any convolution applied later. See BatchPredictSkyESO.py for a
# caller that both overrides this and corrects for it.
wgauss     = %.3f

# FWHM of Lorentzian kernel [pixels]
wlorentz   = 0.8

# variable kernel (width proportional to wavelength)? -> 1 = yes; 0 = no
# if varkern = 1: kernel radius and FWHM for central wavelength
varkern    = 1

# output file for kernel ("stdout": screen; "null": no output)
kernelfile = output/kernel.dat
'''


obs_base = '''
# observatory height in km [2.4, 3.06] (default: 2.64, Cerro Paranal --
# also PALACE's own hardcoded reference height, see EsoSkyObs.py Notes)
sm_h = %.2f

# lower height limit in km (default: 2.0)
sm_hmin = 2.0

# altitude of object above horizon [0,90]
alt      = %.1f

# separation of Sun and Moon as seen from Earth [0,360]
# (> 180 for waning Moon)
alpha    = %.1f

# separation of Moon and object [0,180]
rho      = %.1f

# altitude of Moon above horizon [-90,90]
altmoon  = %.1f

# distance to Moon (mean distance = 1; [0.91,1.08])
moondist = %.2f

# pressure at observer altitude in hPa (default: 744; site defaults below
# -- see SITE_PRESSURE_HPA / -pres)
pres     = %.1f

# single scattering albedo for aerosols [0,1] (default: 0.97)
ssa      = 0.97

# calculation of double scattering of moonlight ('Y' or 'N')
calcds   = N

# relative UV/optical ozone column density (1 -> 258 DU)
o3column = 1.

# scaling factor for scattered moonlight (default: 1.0)
moonscal = 1.0

# heliocentric ecliptic longitude of object [-180,180]
lon_ecl  = %.1f

# ecliptic latitude of object [-90,90]
lat_ecl  = %.1f

# grey-body emissivity (comma-separated list)
emis_str = 0.2

# grey-body temperature in K (comma-separated list)
temp_str = 290.

# monthly-averaged solar radio flux [sfu]
msolflux = %.1f

# bimonthly period (1: Dec/Jan, ..., 6: Oct/Nov; 0: entire year)
season   = 0

# period of the night (x/3 of night, x = 1,2,3; 0: entire night)
time     = 0

# vac[uum] or air wavelengths
vac_air  = air

# precipitable water vapour in mm (-1: bimonthly mean)
pwv      = 3.5

# radiative transfer code L(BLRTM) or R(FM) for molecular spectra
rtcode   = L

# resolution of molecular spectra in library (crucial for run time)
# resol    = 1e6
resol    = 6e4

# path to file sm_filenames.dat for data paths and file names
filepath = data

# inclusion of sky model components
# format: "xxxxxxx" where x = "Y" (yes) or x = "N" (no)
# pos. 1: scattered moonlight
#      2: scattered starlight
#      3: zodiacal light
#      4: thermal emission by telescope/instrument
#      5: molecular emission of lower atmosphere
#      6: sky emission lines of upper atmosphere
#      7: airglow con
incl     = YYYYYYY
'''


# Observatory height [km] used for the local engine's sm_h parameter --
# atmospheric-physics altitude only.  Geometry (alt/moon separation/etc.)
# always comes from get_info_las_campanas's real LCO coordinates regardless
# of this setting -- see run_sky_obs's site parameter and Notes for why
# 'paranal' only swaps the height/pressure physics, not the geometry (the
# same mixed real-geometry / Paranal-physics approach PALACE itself uses,
# see PalaceObs.py's own Notes on its hardcoded Cerro Paranal constants).
SITE_HEIGHT_KM = {
    'lco': 2.5,        # historical default; not an exact match to LCO's
                       # real altitude (2.281 km) or any named ESO preset,
                       # just this module's longstanding value
    'paranal': 2.64,   # Cerro Paranal; matches PALACE's own hardcoded h
}


# Pressure [hPa] used for the local engine's pres parameter, keyed by site --
# same role as SITE_HEIGHT_KM but for pressure rather than height.  Only
# applies to the local engine; skycalc_cli (remote) does not expose pressure
# as a separate parameter, only via its 'observatory' name.  Overridable
# per call via -pres / create_local_inputs's pressure argument.
SITE_PRESSURE_HPA = {
    'lco': 765.0,       # nominal LCO barometric pressure
    'paranal': 744.0,   # the model's own built-in default; matches PALACE's
                        # own hardcoded p, see SITE_HEIGHT_KM
}


def create_local_inputs(ra=296.242608, dec=-14.811007, obstime='2023-08-29T03:20:43.668', msol=0,
                        site='lco', pressure=None, workdir='.', verbose=False,
                        lsf_gauss_fwhm_pix=0.8):
    '''
    Write <workdir>/config/skymodel_etc.par and
    <workdir>/config/instrument_etc.par for calcskymodel, given the
    observing geometry and a resolved solar flux.  If msol<=0, the
    model's own long-term-average default (101 sfu) is used.  site
    selects SITE_HEIGHT_KM's sm_h value ('lco' default, or 'paranal').
    pressure overrides SITE_PRESSURE_HPA's site-keyed default (hPa) if
    given and positive.
    workdir defaults to '.' but run_local() always passes an isolated
    per-call temporary directory instead -- see module Notes.

    lsf_gauss_fwhm_pix: FWHM (pixels) of calcskymodel's own internal
    Gaussian LSF convolution (inst_base's wgauss). Default 0.8 matches
    the value hardcoded here before this became a parameter -- every
    existing caller sees unchanged output. Do not set this much lower
    to try to get an "unsmoothed" spectrum: at the default 0.5 A/pixel
    output grid (dlam), this convolution is also calcskymodel's
    anti-aliasing filter going from the model's native library
    resolution (resol, ~R=60000) down to that grid, so an unphysically
    narrow value risks aliasing artifacts rather than a genuinely finer
    spectrum. See BatchPredictSkyESO.py for a caller that overrides
    this to a known value and corrects for it analytically instead.
    '''
    info = get_info_las_campanas(obstime, ra=ra, dec=dec, verbose=verbose)

    longitude = info['SourceEclipLon'] - info['SunEclipLon']
    if longitude > 180:
        longitude -= 360.
    if longitude < -180.:
        longitude += 360.

    xmsol = msol if msol > 0 else 101.
    sm_h = SITE_HEIGHT_KM[site]
    xpres = pressure if pressure and pressure > 0 else SITE_PRESSURE_HPA[site]

    with open(os.path.join(workdir, 'config', 'skymodel_etc.par'), 'w') as xout:
        xout.write(obs_base % (sm_h, info['SourceAlt'], info['Moon-Sun_Separation'], info['Moon-Source_Separation'],
                                info['MoonAlt'], info['MoonDistanceInMeanUnits'], xpres,
                                longitude, info['SourceEclipLat'], xmsol))

    with open(os.path.join(workdir, 'config', 'instrument_etc.par'), 'w') as xinst:
        xinst.write(inst_base % lsf_gauss_fwhm_pix)


def local_engine_available(eso_sky_dir=''):
    '''
    Check whether a working local ESO Sky Model install can be found, and
    if so return (True, path_to_calcskymodel); otherwise (False, '').
    '''
    xdir = eso_sky_dir or os.getenv('ESO_SKY_MODEL')
    if not xdir:
        return False, ''
    binary = os.path.join(xdir, 'sm-01_mod2', 'bin', 'calcskymodel')
    data_dir = os.path.join(xdir, 'sm-01_mod2', 'data')
    if os.path.isfile(binary) and os.access(binary, os.X_OK) and os.path.isdir(data_dir):
        return True, binary
    return False, ''


def run_local(ra, dec, xtime_iso, msol=0, outroot='', eso_sky_dir='', site='lco', pressure=None,
             keep_workdir=False, lsf_gauss_fwhm_pix=0.8):
    '''
    Run the local ESO Sky Model (calcskymodel) for ra, dec, xtime_iso and
    write outroot.fits (in the CALLER's current directory -- unchanged
    from before).  Returns outroot on success, '' on failure.
    site: 'lco' (default) or 'paranal' -- see SITE_HEIGHT_KM.
    pressure overrides SITE_PRESSURE_HPA's site-keyed default (hPa) if
    given and positive -- see create_local_inputs().
    lsf_gauss_fwhm_pix: passed through to create_local_inputs() -- see
    its docstring; default 0.8 is unchanged behavior.

    By default, everything calcskymodel itself reads/writes (config/*.par,
    output/*.fits/.dat, the data symlink) lives in a fresh
    tempfile.TemporaryDirectory() private to this one call, not the
    caller's working directory -- see module Notes.  This is what makes
    concurrent calls (separate processes OR threads of one process) safe:
    each gets its own directory with a globally-unique random name, so
    there's nothing shared to race on.  The directory and everything in
    it is removed automatically when this function returns, even on
    failure/exception.

    keep_workdir=True (-keep_workdir) instead writes directly into
    './config', './data' (symlink), './output' in the CALLER's current
    directory -- the pre-260730 behavior -- and leaves them there, so the
    inputs can be inspected, hand-edited, and calcskymodel rerun directly
    from that same directory.  NOT concurrency-safe -- do not use this
    from more than one call at a time; see _scratch_dir() and module
    Notes.
    '''
    available, binary = local_engine_available(eso_sky_dir)
    if not available:
        print('Error: local ESO Sky Model is not available (check ESO_SKY_MODEL)')
        return ''

    with _scratch_dir('EsoSkyObs_local_', workdir='.' if keep_workdir else None) as workdir:
        setup(eso_sky_dir, workdir=workdir)
        create_local_inputs(ra=ra, dec=dec, obstime=xtime_iso, msol=msol, site=site, pressure=pressure,
                            workdir=workdir, lsf_gauss_fwhm_pix=lsf_gauss_fwhm_pix)

        result = subprocess.run([binary], capture_output=True, text=True, cwd=workdir)
        if len(result.stderr):
            print('stderr:', result.stderr)
            print('Could not create local sky model: ra %f dec %f time %s' % (ra, dec, xtime_iso))
            return ''

        rad = fits.open(os.path.join(workdir, 'output', 'radspec.fits'))
        trans = fits.open(os.path.join(workdir, 'output', 'transspec.fits'))
        rtab = Table(rad[1].data)
        ttab = Table(trans[1].data)
        ztab = join(rtab, ttab, join_type='left')
        ztab['lam'] *= 1000.   # microns -> nm, to match the remote engine's native units
        header = rad[0].header
        rad.close()
        trans.close()

        ztab = finalize_table(ztab)
        write_output_fits(ztab, header, outroot, engine='local', ra=ra, dec=dec, xtime_iso=xtime_iso, msol=msol,
                          site=site)
    return outroot


# ---------------------------------------------------------------------------
# Remote engine (skycalc_cli / ESO web service)
# ---------------------------------------------------------------------------

default = '''
{
    "ra": 121.75,
    "dec": -29.7,
    "date": "2012-07-17T21:12:14",
    "observatory": "lasilla"
}
'''

xdefaults = '''
{
    "vacair": "air",
    "wmin": 360.0,
    "wmax": 980.0,
    "wgrid_mode": "fixed_wavelength_step",
    "wdelta": 0.05,
    "lsf_type": "Gaussian",
    "lsf_gauss_fwhm": 2.0,
    "lsf_boxcar_fwhm": 2.0,
    "observatory": "lasilla"
}
'''


# skycalc_cli's own accepted observatory names for our site keyword -- see
# skycalc_cli.py's fixObservatory (also accepts '2400'/'2640'/'3060'/
# 'armazones' directly, not needed here).  NOTE unlike the local engine's
# SITE_HEIGHT_KM, switching this ALSO changes skycalc_cli's own internal
# moon/sun almanac geometry to the named site's real location, not just its
# height/pressure physics -- see run_sky_obs's site parameter and Notes.
REMOTE_SITE_NAME = {
    'lco': 'lasilla',
    'paranal': 'paranal',
}


def write_remote_inputs(ra, dec, xtime_iso, msol=0, outroot='test', site='lco', workdir='.'):
    '''
    Write <workdir>/<outroot>.json, the per-observation input file for
    skycalc_cli.  site: 'lco' (default) or 'paranal' -- see
    REMOTE_SITE_NAME.  workdir defaults to '.' but run_remote() always
    passes an isolated per-call temporary directory instead -- see
    module Notes.
    '''
    xdict = json.loads(default)
    xdict.update({'ra': ra})
    xdict.update({'dec': dec})
    xdict.update({'observatory': REMOTE_SITE_NAME[site]})

    xtime_obj = Time(xtime_iso, format='isot', scale='utc')
    isot_whole_seconds = xtime_obj.datetime.strftime('%Y-%m-%dT%H:%M:%S')
    xdict.update({'date': isot_whole_seconds})

    if msol > 0:
        xdict['msolflux'] = msol

    jsonString = json.dumps(xdict, indent=4, sort_keys=False)
    with open(os.path.join(workdir, '%s.json' % outroot), 'w') as jsonFile:
        jsonFile.write(jsonString)


def is_recent(filepath, minutes=10):
    '''Return True if file was modified in the last `minutes`.'''
    if not os.path.isfile(filepath):
        return False
    age_seconds = time.time() - os.path.getmtime(filepath)
    return age_seconds < (minutes * 60)


def run_remote(ra, dec, xtime_iso, msol=0, outroot='', print_output=False, site='lco', keep_workdir=False):
    '''
    Run the ESO SkyCalc web service via skycalc_cli for ra, dec, xtime_iso
    and write outroot.fits (in the CALLER's current directory --
    unchanged from before).  Returns outroot on success, '' on failure.
    site: 'lco' (default) or 'paranal' -- see REMOTE_SITE_NAME.

    By default, skycalc_cli's own scratch inputs/output (<outroot>.json,
    xsky.json, and its raw <outroot>.fits before finalize_table/
    write_output_fits convert it) live in a fresh
    tempfile.TemporaryDirectory() private to this one call, same
    reasoning as run_local() -- see module Notes.

    keep_workdir=True (-keep_workdir) instead writes those scratch files
    directly into the CALLER's current directory and leaves them there --
    the pre-260730 behavior.  NOT concurrency-safe; see _scratch_dir().
    '''
    with _scratch_dir('EsoSkyObs_remote_', workdir='.' if keep_workdir else None) as workdir:
        write_remote_inputs(ra, dec, xtime_iso, msol=msol, outroot=outroot, site=site, workdir=workdir)

        xsky_dict = json.loads(xdefaults)
        xsky_dict['observatory'] = REMOTE_SITE_NAME[site]
        with open(os.path.join(workdir, 'xsky.json'), 'w') as g:
            json.dump(xsky_dict, g, indent=4, sort_keys=False)
        time.sleep(1)

        command_line = ('skycalc_cli -i xsky.json -a %s.json -o %s.fits' % (outroot, outroot)).split()

        try:
            result = subprocess.run(command_line, capture_output=True, text=True, cwd=workdir)
        except FileNotFoundError:
            print("Error: 'skycalc_cli' not found in your PATH.")
            return ''

        Xerror = result.returncode != 0 or 'Traceback' in result.stderr
        if print_output or Xerror:
            print('stdout:', result.stdout)
            print('stderr:', result.stderr)
        if Xerror:
            print('ERROR: skycalc_cli failed to execute successfully.')
            return ''

        raw_fits = os.path.join(workdir, '%s.fits' % outroot)
        if is_recent(raw_fits) == False:
            print('Error: cannot verify that %s.fits was recently created' % outroot)
            return ''

        x = fits.open(raw_fits)
        ztab = Table(x[1].data)
        header = x[0].header
        x.close()

        ztab = finalize_table(ztab)
        write_output_fits(ztab, header, outroot, engine='remote', ra=ra, dec=dec, xtime_iso=xtime_iso, msol=msol,
                          site=site)
    return outroot


# ---------------------------------------------------------------------------
# Shared physics: raw model table -> LVM sky-model FITS convention
# ---------------------------------------------------------------------------

def finalize_table(ztab):
    '''
    Convert a raw ESO Sky Model table (WAVE already in nm; flux components
    in the model's native photon units) into the LVM sky-model convention.

    Both engines compute the sky as it would actually be observed on the
    ground (i.e. after atmospheric extinction).  LVM spectra are normally
    compared against the above-the-atmosphere equivalent, so FLUX and the
    MOON/ZODI/DIFFUSE/CONT components are corrected for extinction here by
    dividing by the transmission; the as-observed value is kept separately
    in FLUX_LCO.
    '''
    ztab.rename_column('lam', 'WAVE')
    ztab.rename_column('flux', 'FLUX')
    ztab.rename_column('flux_sml', 'MOON')
    ztab.rename_column('flux_zl', 'ZODI')
    ztab.rename_column('flux_ael', 'LINES')
    ztab.rename_column('flux_arc', 'DIFFUSE')

    q = 1.98644586e-17 * FIBER_AREA_ARCSEC2 / ztab['WAVE']
    for col in ('FLUX', 'MOON', 'LINES', 'ZODI', 'DIFFUSE'):
        ztab[col] *= q

    ztab['WAVE'] *= 10.   # nm -> Angstrom

    ztab['CONT'] = ztab['MOON'] + ztab['ZODI'] + ztab['DIFFUSE']

    ztab['FLUX_LCO'] = ztab['FLUX'].copy()   # as actually observed on the ground at LCO

    for col in ('FLUX', 'CONT', 'MOON', 'ZODI', 'DIFFUSE'):
        ztab[col] /= ztab['trans']

    return ztab


def write_output_fits(ztab, header, outroot, engine, ra, dec, xtime_iso, msol, site='lco'):
    header = header.copy()
    header['RA'] = ra
    header['DEC'] = dec
    header['OBSTIME'] = xtime_iso
    header['MSOLFLUX'] = msol
    header['ENGINE'] = engine
    header['SITE'] = site
    primary_hdu = fits.PrimaryHDU(header=header)
    table_hdu = fits.BinTableHDU(data=ztab)
    fits.HDUList([primary_hdu, table_hdu]).writeto('%s.fits' % outroot, overwrite=True)


# ---------------------------------------------------------------------------
# Top-level orchestration
# ---------------------------------------------------------------------------

def resolve_solar_flux(xtime_iso, msol=-1):
    '''
    Resolve the 10.7 cm solar radio flux to use.  If msol is already
    positive, use it as given.  Otherwise look one up from the historical
    flux table -- unless the observation is later than 2024, which falls
    outside the range we trust that table for, in which case we fall back
    to the model's own built-in default (signalled by returning 0).
    '''
    year = int(xtime_iso[:4])
    if year > 2024:
        return 0
    if msol > 0:
        return msol
    msol = get_flux(xtime_iso)
    print('Got Solar flux of ', msol)
    return msol


def run_sky_obs(ra, dec, xtime, outroot='', msol=-1, engine='auto', print_output=False, eso_sky_dir='',
                site='lco', pressure=None, keep_workdir=False, lsf_gauss_fwhm_pix=0.8):
    '''
    Generate a predicted sky spectrum for ra, dec, xtime using the real ESO
    Sky Model.  By default (engine='auto') the local calcskymodel install is
    tried first, falling back to the ESO SkyCalc web service if the local
    model is not available.  Returns the output root name (without .fits),
    or '' on failure.

    site: 'lco' (default) or 'paranal' -- see module Notes; primarily for
    comparison runs against PALACE, whose own atmosphere physics is fixed
    to Cerro Paranal.

    pressure overrides SITE_PRESSURE_HPA's site-keyed default (hPa) for the
    local engine only -- see create_local_inputs()/-pres in the CLI help;
    has no effect on the remote engine (skycalc_cli exposes no separate
    pressure parameter).

    lsf_gauss_fwhm_pix: local engine only (no effect on the remote engine,
    which has its own separate, unrelated lsf_gauss_fwhm default in
    xdefaults) -- passed through to run_local()/create_local_inputs();
    default 0.8 is unchanged behavior. See create_local_inputs()'s
    docstring before changing this -- there is a real numerical floor
    below which this stops being "less smoothing" and starts being
    aliasing.

    keep_workdir=True (-keep_workdir) trades the default concurrency-safe
    behavior (a fresh, auto-deleted temp directory per call) for the pre-
    260730 behavior: writes directly into './config', './data' (symlink),
    './output' in the CALLER's current directory and leaves them there,
    for debugging -- the inputs can be inspected, hand-edited, and
    calcskymodel rerun directly from that same directory.  NOT
    concurrency-safe; see run_local()/run_remote()/_scratch_dir().
    '''
    xtime_iso = convert_time(xtime, 'iso_ms')
    msol_resolved = resolve_solar_flux(xtime_iso, msol)

    if outroot == '':
        mjd = convert_time(xtime_iso, 'mjd')
        if dec > 0:
            outroot = 'SkyE_%8.2f_%05.1f_+%04.1f' % (mjd, ra, dec)
        else:
            outroot = 'SkyE_%8.2f_%05.1f_%.1f' % (mjd, ra, dec)

    used_engine = ''
    xroot = ''

    if engine in ('local', 'auto'):
        xroot = run_local(ra, dec, xtime_iso, msol=msol_resolved, outroot=outroot, eso_sky_dir=eso_sky_dir,
                          site=site, pressure=pressure, keep_workdir=keep_workdir,
                          lsf_gauss_fwhm_pix=lsf_gauss_fwhm_pix)
        if xroot:
            used_engine = 'local'
        elif engine == 'local':
            return ''

    if xroot == '' and engine in ('remote', 'auto'):
        if engine == 'auto':
            print('Local ESO Sky Model unavailable or failed; falling back to the ESO SkyCalc web service')
        xroot = run_remote(ra, dec, xtime_iso, msol=msol_resolved, outroot=outroot, print_output=print_output,
                           site=site, keep_workdir=keep_workdir)
        if xroot:
            used_engine = 'remote'

    if xroot == '':
        print('Error: could not generate a sky model with either engine')
        return ''

    print('Used the %s engine; wrote %s.fits' % (used_engine, xroot))
    return xroot


def steer(argv):
    '''
    Usage:  EsoSkyObs.py [-h] [-engine local|remote|auto] [-msol flux] [-out root] ra dec time
    '''

    ra = -1.0
    dec = -100
    xtime = -1
    outroot = ''
    msol = -1
    engine = 'auto'
    site = 'lco'
    pres = -1
    keep_workdir = False

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(_USAGE)
            return
        elif argv[i][:4] == '-out':
            i += 1
            outroot = argv[i]
        elif argv[i][:5] == '-msol':
            i += 1
            msol = eval(argv[i])
        elif argv[i][:7] == '-engine':
            i += 1
            engine = argv[i]
        elif argv[i][:5] == '-site':
            i += 1
            site = argv[i]
        elif argv[i][:5] == '-pres':
            i += 1
            pres = eval(argv[i])
        elif argv[i][:13] == '-keep_workdir':
            keep_workdir = True
        elif argv[i][0] == '-' and ra < 0.0:
            print(_USAGE)
            print('Error: unknown option: ', argv)
            return
        elif ra < 0.0:
            ra = eval(argv[i])
        elif dec < -90:
            dec = eval(argv[i])
        elif xtime == -1:
            xtime = argv[i]
        i += 1

    if engine not in ('local', 'remote', 'auto'):
        print('Error: -engine must be one of local, remote, auto')
        return

    if site not in ('lco', 'paranal'):
        print('Error: -site must be one of lco, paranal')
        return

    run_sky_obs(ra=ra, dec=dec, xtime=xtime, outroot=outroot, msol=msol, engine=engine, site=site,
               pressure=pres if pres > 0 else None, keep_workdir=keep_workdir)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
