#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Generate data/moon_base_spectrum.dat: a static, precomputed
    "reddened and atmospherically-scattered" moonlight spectral shape,
    for use as the MOON template in SkySepPalace.py (and, eventually,
    any other PALACE-based sky model here) without a runtime dependency
    on the ESO Sky Model.

Command line usage (if any):

    MakeMoonBase.py [-ra RA] [-dec DEC] [-obstime TIME] [-out FILE]

Options::

    -ra RA         reference right ascension in degrees (default: a
                   fixed, moderate-airmass, moderate-phase test case)
    -dec DEC       reference declination in degrees
    -obstime TIME  reference UTC time
    -out FILE      output file path (default: data/moon_base_spectrum.dat)

Description:

    A raw solar spectrum reddened only by the Moon's own wavelength-
    dependent albedo is NOT what SkySepPalace.py's MOON template needs
    to match -- that template is fit directly against the *observed*
    sky spectrum, i.e. moonlight *after* it has been scattered by the
    Earth's atmosphere into the line of sight, and Rayleigh/aerosol
    scattering (favouring blue) empirically dominates over and reverses
    the lunar-albedo reddening (verified 260709 by comparing a Kieffer &
    Stone (2005) ROLO-albedo-corrected solar spectrum against the real
    ESO Sky Model's own MOON output for the same geometry: the ROLO-only
    version got redder toward longer wavelengths; the real, fully
    scattered ESO MOON spectrum gets bluer -- the opposite trend).

    Rather than reimplement the ESO Sky Model's full 3D Rayleigh/aerosol
    radiative transfer calculation (see "The Cerro Paranal Advanced Sky
    Model" manual, Sect. 6.2.1) or fit a simplified power-law
    approximation to it (tried: even the best-fit ~lambda^-2.5 power law
    left a much larger residual than just using the real ESO output
    directly), this script runs the real, already-installed, already-
    validated ESO Sky Model (EsoSkyObs.run_sky_obs, engine='local') ONCE,
    for one fixed reference geometry, and stores its MOON column as a static spectral
    shape.  This keeps SkySepPalace.py itself free of any runtime ESO
    dependency (see Readme.md) while still using a physically complete
    albedo+scattering shape rather than an uncorrected or wrongly-
    corrected solar spectrum.

    Output columns (Angstroms, on the ESO model's own 3600-9800 A grid):
    WAVE, SOLAR (vendored solar reference spectrum, interpolated onto
    this grid), MOON_BASE (the ESO MOON column, normalised to unit
    median flux over the same range as SOLAR's normalisation, so the
    two are directly comparable in shape).

Notes:

    This uses ONE fixed reference geometry (a specific lunar phase
    angle and airmass for both target and Moon), not the actual phase
    angle/airmass of each observation being sky-subtracted -- a
    deliberate, documented simplification ("for now"; see
    SkySepPalace.py Notes).  The real per-observation phase-angle
    dependence (steeper reddening near new Moon, flatter near full) and
    airmass dependence of the scattering are both real effects this
    single static curve does not capture.

    EsoSkyObs.setup() (called by run_local(), via run_sky_obs()) writes
    config/, output/, and a data/ symlink into the *current working
    directory* -- run this script from a scratch directory, not from
    py_progs/, or clean those up afterward (they are not meant to be
    tracked).

    data/moon_rolo_albedo.dat (the raw Kieffer & Stone 2005 ROLO
    coefficients) is not used by this script -- it was the first,
    insufficient attempt (albedo alone, wrong direction once atmospheric
    scattering is included) -- but is kept in the repository since the
    physics it documents is still correct and may be useful again if a
    proper per-observation phase-dependent correction is built later
    (e.g. by scaling this base spectrum's *shape* while independently
    tracking the ROLO albedo's phase dependence for the *overall*
    amplitude).

History::

    260709 ksl Coding begun.
    260711 ksl Migrated from SkyModelObs.do_one (retired) to
        EsoSkyObs.run_sky_obs(engine='local') -- forced to 'local' rather
        than 'auto' since this script's whole point is running one
        specific, known model install, not silently falling back to a
        different engine if the local one isn't set up.

'''

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from astropy.io import fits
from astropy.table import Table

import EsoSkyObs

_USAGE = '''Usage:
  MakeMoonBase.py [-ra RA] [-dec DEC] [-obstime TIME] [-out FILE]

Options:
  -ra RA         reference right ascension in degrees
  -dec DEC       reference declination in degrees
  -obstime TIME  reference UTC time
  -out FILE      output file path (default: data/moon_base_spectrum.dat)
'''

# A fixed reference geometry with moderate airmass for both target and
# Moon (source alt ~67 deg, Moon alt ~83 deg) and a moderate phase angle
# (~28 deg from full) -- not an edge case, but see Notes above: this is
# a single fixed stand-in, not the real per-observation geometry.
DEFAULT_RA      = 296.242608
DEFAULT_DEC     = -14.811007
DEFAULT_OBSTIME = '2023-08-29T03:20:43.668'

DEFAULT_SOLAR_FILE = (os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                                   'data', 'palace_ref',
                                   'Spectre_HR_LATMOS_Meftah_V1_350_1000nm.txt'))
DEFAULT_OUTFILE = (os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                                'data', 'moon_base_spectrum.dat'))


def make_moon_base(ra=DEFAULT_RA, dec=DEFAULT_DEC, obstime=DEFAULT_OBSTIME,
                   outfile=DEFAULT_OUTFILE):
    '''
    Run the ESO Sky Model once for (ra, dec, obstime) and write
    data/moon_base_spectrum.dat with WAVE/SOLAR/MOON_BASE columns.
    '''
    outroot = EsoSkyObs.run_sky_obs(ra=ra, dec=dec, xtime=obstime,
                                    outroot='/tmp/MakeMoonBase_tmp', engine='local')
    if outroot == '':
        print('Error: ESO Sky Model call failed; no output written')
        return None

    modelfile = '%s.fits' % outroot
    hdul = fits.open(modelfile)
    model_tab = Table(hdul[1].data)
    hdul.close()
    if os.path.exists(modelfile):
        os.remove(modelfile)

    wave = np.array(model_tab['WAVE'], dtype=float)   # Angstroms
    moon = np.array(model_tab['MOON'], dtype=float)

    sol = np.loadtxt(DEFAULT_SOLAR_FILE, comments=';')
    sol_wave_ang = sol[:, 0] * 10.0   # nm -> Angstroms
    sol_flux = sol[:, 1]
    solar_on_grid = np.interp(wave, sol_wave_ang, sol_flux)

    sel = np.isfinite(moon) & (moon > 0)
    solar_med = np.nanmedian(solar_on_grid[sel])
    moon_med  = np.nanmedian(moon[sel])
    solar_norm = solar_on_grid / solar_med
    moon_norm  = moon / moon_med

    out_tab = Table([wave, solar_norm, moon_norm],
                    names=['WAVE', 'SOLAR', 'MOON_BASE'])

    header = (
        '# Static moonlight base spectrum for SkySepPalace.py MOON template.\n'
        '# WAVE: Angstroms. SOLAR: vendored solar reference spectrum\n'
        '# (data/palace_ref/Spectre_HR_LATMOS_Meftah...), interpolated onto\n'
        '# this grid, normalised to unit median. MOON_BASE: the real ESO\n'
        '# Sky Model MOON column (Rayleigh/aerosol-scattered, ROLO-albedo-\n'
        '# reddened moonlight) for one fixed reference geometry, normalised\n'
        '# to unit median -- see MakeMoonBase.py module docstring for why\n'
        '# this is not simply SOLAR x lunar albedo.\n'
        '# Reference geometry: ra=%.6f dec=%.6f obstime=%s\n'
        '# Generated by MakeMoonBase.py\n' % (ra, dec, obstime)
    )
    with open(outfile, 'w') as fh:
        fh.write(header)
        for row in out_tab:
            fh.write('%10.3f %14.6e %14.6e\n' % (row['WAVE'], row['SOLAR'], row['MOON_BASE']))

    print('Wrote %s (%d rows)' % (outfile, len(out_tab)))
    return outfile


def steer(argv):
    ra = DEFAULT_RA
    dec = DEFAULT_DEC
    obstime = DEFAULT_OBSTIME
    outfile = DEFAULT_OUTFILE

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-ra':
            i += 1
            ra = float(argv[i])
        elif argv[i] == '-dec':
            i += 1
            dec = float(argv[i])
        elif argv[i] == '-obstime':
            i += 1
            obstime = argv[i]
        elif argv[i] == '-out':
            i += 1
            outfile = argv[i]
        else:
            print('Error: unknown argument %s' % argv[i])
            print(_USAGE)
            return
        i += 1

    make_moon_base(ra=ra, dec=dec, obstime=obstime, outfile=outfile)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        make_moon_base()
