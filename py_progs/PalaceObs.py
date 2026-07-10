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

    predict     main entry point; returns (species_table, combined_table)
    do_one      CLI-friendly wrapper; writes a FITS file with SPECIES and
                COMBINED extensions

Notes:

    PALACE's climatology is calibrated for Cerro Paranal, not Las
    Campanas; the two sites are close (LCO -29.01,-70.69 vs Paranal
    -24.63,-70.40, both high-altitude Atacama sites) so this is a
    reasonable approximation for now, not an exact match.

    Local solar time (PALACE's ``tbin``) is approximated as
    UTC + longitude/15h, ignoring the equation of time (<= 16 min) --
    negligible next to PALACE's hour-wide time bins.

    PALACE's spectra are in Rayleighs/nm on its own wavelength grid; no
    unit conversion to LVM's flux-calibrated grid is done here yet.

History:

    260709 ksl Coding begun.  First cut: loop the public palace.model()
        over each of PALACE's 9 species plus once for the combined
        spectrum.  Geometry/time reuses SkyModelObs.get_info_las_campanas;
        solar flux reuses GetSolar.get_flux.

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

# Las Campanas Observatory longitude, matching SkyModelObs.get_info_las_campanas
LCO_LON_DEG = -70.6920


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
          **overrides):
    '''
    CLI-friendly wrapper around predict(): writes a FITS file with
    SPECIES and COMBINED extensions and returns the output filename.
    '''

    species_table, combined = predict(ra, dec, obstime, srf=srf,
                                      species_list=species_list,
                                      want_combined=True, **overrides)

    if outroot == '':
        mjd = convert_time(obstime, 'mjd')
        outroot = 'PalaceObs_%.5f_%.5f_%08.2f' % (ra, dec, mjd)
    outname = outroot if outroot.endswith('.fits') else outroot + '.fits'

    primary = fits.PrimaryHDU()
    primary.header['RA'] = ra
    primary.header['DEC'] = dec
    primary.header['OBSTIME'] = str(obstime)
    hdul = fits.HDUList([primary,
                         fits.BinTableHDU(data=species_table, name='SPECIES'),
                         fits.BinTableHDU(data=combined, name='COMBINED')])
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
