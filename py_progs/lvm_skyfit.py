#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Fit the same set of emission lines lvm_gaussfit.py does, on the
*un-subtracted* FLUX of an lvmCFrame exposure, using a local continuum
model of ``const + scale * SKY(wavelength)`` (the raw spectrum of one
of the two sky telescopes) instead of a flat constant, to correct for
un/mis-subtracted moonlit-sky structure (Fraunhofer features etc.)
biasing these lines' fluxes.

Command line usage::

    usage: lvm_skyfit.py [-h] [-out root] [-v vel] [-near] [-tel EAST|WEST]
                          [-wing] filename ...

    filename is one or more lvmCFrame-*.fits files.  -out root sets the
    output rootname (defaults to each input filename's basename).  -v
    vel applies a velocity offset (km/s) to the initial line-center
    guesses, same convention as lvm_gaussfit.py.

    By default the FARTHER of SKY_EAST/SKY_WEST (larger
    SCI_SKYE_SEP/SCI_SKYW_SEP) is used as the sky template -- this
    avoids a large, diffuse target like Vela's own emission leaking
    into a sky telescope that happens to sit close to the field.  -near
    selects the closer telescope instead.  -tel EAST|WEST forces a
    specific telescope regardless of separation (highest precedence;
    mainly for testing).

    -wing fits const/scale first from the line-excluded wings only,
    then holds them fixed while fitting the line (see
    fit_singlet_with_sky/fit_doublet_with_sky's fix_background) --
    useful when a coincident feature in the sky template would
    otherwise bias a jointly-fit scale.  Default (without -wing) is a
    single joint fit of scale*SKY + const + line, all free together,
    from the same pixels -- the same convention lvm_gaussfit.py uses
    for its flat background, just with a sky-shaped one.

Description::

    lvm_gaussfit.py fits a Gaussian plus a single flat constant
    background, jointly, from the same narrow wavelength window as the
    line -- background and line flux are degenerate free parameters
    with no independent handle on what the true local continuum looks
    like.  This is fine when the continuum really is flat across the
    window, but real moonlit-sky spectra have structure (Fraunhofer
    absorption, molecular bands) near some of these lines, and if that
    structure isn't perfectly removed it can bias the flat-background
    fit in either direction.

    This script instead fits ``const + scale * SKY(wavelength)`` as the
    background, where SKY is the *raw* spectrum of one of the two sky
    telescopes (SKY_EAST or SKY_WEST) -- a single, file-wide choice
    (see -near/-tel above), made from the SCI_SKYE_SEP/SCI_SKYW_SEP
    header keywords.  The fit is done on lvmCFrame's FLUX extension,
    which is *not* sky-subtracted (unlike lvmSFrame), so this SKY
    spectrum genuinely is (approximately) present in the data being
    fit.  A best-fit ``scale`` far from whatever value would mean "no
    local correction needed" is direct, quantitative, per-line evidence
    that the local sky component isn't simply the raw telescope
    spectrum at face value -- neither proof nor disproof of the
    pipeline's own (different, more elaborate) sky subtraction, just an
    independent check using the same raw data.

    Every line lvm_gaussfit.py fits is fit here too, with the same
    windows/initial wavelengths and the same singlet-vs-doublet
    grouping: [OII] is the only joint doublet (two Gaussians sharing one
    FWHM, independent -- freely-ratioed -- fluxes, one shared
    const/scale); every other line (Hbeta, [OIII]a/b, [OI]a/b, Halpha,
    [NII]a/b, [SII]a/b, HeII, HeI, [OIII]4363, [SIII]a/b) is an
    independent singlet with its own const/scale, matching
    lvm_gaussfit.py's do_one exactly -- including [SII]a/[SII]b, which
    are two separate singlets there, not a joint doublet.

    Unlike lvm_gaussfit.py, [OI]a/[OI]b are *not* pre-filtered to
    exclude the 6298-6303/6362-6365 Ang sky-line region -- that
    exclusion existed there specifically because a flat background
    couldn't handle the sky-subtraction residual in that window; the
    scale*SKY term here is meant to model exactly that kind of
    structure directly, so excluding it would prevent ever testing
    whether it does.

Units::

    lvmCFrame stores flux density in erg/s/cm**2/Ang.  As in
    lvm_gaussfit.py, FLUX/SKY/ERROR are all multiplied by 1e16 before
    fitting for numerical convenience.  Because SKY is scaled by the
    same factor as FLUX, the fitted ``scale`` parameter is dimensionless
    and unaffected by this choice.  flux_<line> columns are in units of
    erg/s/cm**2 * 1e16 (divide by 1e16 for physical flux); back_<line>
    (the ``const`` term) is in erg/s/cm**2/Ang * 1e16, matching
    lvm_gaussfit.py's background column exactly.

Primary routines::

    fit_singlet_with_sky   - one singlet line: Gaussian + const + scale*SKY
    fit_doublet_with_sky   - [OII] doublet: two Gaussians (shared fwhm) + const + scale*SKY
    do_one                 - every line group for one fiber's spectrum
    do_all                 - one lvmCFrame file: fiber loop, results table
                              + diagnostic FITS
    do_individual           - loop do_all over multiple lvmCFrame files
    steer                   - command-line driver

Notes::

    The diagnostic FITS output (``<root>.skyfit_model.fits``) is NOT
    shaped like the full input spectrum -- the 16 fit windows together
    total only ~750 of the input's 12401 wavelength pixels (~6%), so a
    full-spectrum-shaped output would be far larger than necessary and
    mostly empty. Instead it holds one small, already-zoomed-in cutout
    per line group: for each line, extensions WAVE_<LINE> (1-D
    wavelength grid for that window), FLUX_<LINE> (original data cutout,
    science fibers only), SKY_<LINE> (fitted const + scale*SKY --
    everything in the model that is not the line, for comparing against
    the shape of FLUX away from the line core), LINE_<LINE> (fitted line
    component alone, zero baseline), and MOD_<LINE> (SKY_<LINE> +
    LINE_<LINE>, the full model). Rows are science fibers only, in the
    same order as the results table, so row i in any extension is the
    same fiber as row i of the results table -- no id-matching needed.
    A FIBERID/RA/DEC table extension is included for identification
    without cross-referencing the original CFrame.

    Uses only the raw SKY_EAST/SKY_WEST telescope spectra, not the
    pipeline's own combined per-fiber sky estimate -- deliberately, to
    test something independent of however that combination was built.

History::

    260825 ksl Coding begun
    260826 ksl Expanded from 4 line groups (oii, hb, ha, sii-as-doublet)
    to the full lvm_gaussfit.py set (16 line groups, [SII]a/[SII]b now
    independent singlets matching lvm_gaussfit.py exactly); dropped
    [OI]'s sky-line pre-exclusion (scale*SKY is meant to model that
    structure, not avoid it); joint fit (scale*SKY+const+line, all free
    together) is now the default, wing-only background moved behind
    -wing; sky-telescope choice now defaults to the FARTHER telescope
    (was nearer), with -near to select the closer one instead of -tel's
    hard override.

'''

import os
import re
import sys

import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table, vstack, hstack
from lmfit import Model

from lvm_gaussfit import scifib, check_for_nan, patch_stderr_fraction


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before a line consisting of "History:" (or
    "History::"/"Version History" -- whitespace/colon-insensitive), so
    -h stays short even as that section grows -- without hand-
    duplicating the Synopsis/Options text in a second string.
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


# (name, init_wavelength, init_fwhm, wavelength_min, wavelength_max) -- same
# wavelengths/windows lvm_gaussfit.py's do_one uses for these lines (in the
# same order it fits them, skipping oii which is DOUBLETS below), so results
# are directly comparable. Unlike lvm_gaussfit.py, oi_a/oi_b are NOT
# pre-filtered to exclude the 6298-6303/6362-6365 Ang sky-line region here
# (see module Description).
SINGLETS = [
    ('hb', 4861.325, 1., 4855, 4870),
    ('oiii_a', 4958.91, 1., 4949, 4969),
    ('oiii_b', 5006.843, 1., 4997, 5017),
    ('oi_a', 6300.304, 1., 6290, 6310),
    ('oi_b', 6363.777, 1., 6353, 6373),
    ('ha', 6562.8, 1., 6555, 6570),
    ('nii_a', 6548.04, 1., 6538, 6558),
    ('nii_b', 6583.46, 1., 6574, 6594),
    ('sii_a', 6716.44, 1., 6706, 6726),
    ('sii_b', 6730.81, 1., 6721, 6741),
    ('heii', 4685.71, 1., 4666, 4706),
    ('hei', 5875.6, 1., 5856, 5886),
    ('oiii_4363', 4363.21, 1., 4343, 4383),
    ('siii_a', 9068.6, 1., 9055, 9090),
    ('siii_b', 9530.6, 1., 9525, 9545),
]

# (name, init_wavelength1, init_wavelength2, init_fwhm, wavelength_min, wavelength_max)
# -- [OII] is the only line lvm_gaussfit.py fits as a true joint doublet;
# everything else, including [SII]a/[SII]b, is in SINGLETS above.
DOUBLETS = [
    ('oii', 3726.092, 3729.875, 1., 3717, 3737),
]

# Default half-width (Ang) excluded around each expected line center when
# fitting const/scale from the wings only (see _fit_wing_background) --
# a real emission line (or a coincident sky-template feature) should not
# extend past this from its rest wavelength given the ~1-2 Ang FWHMs
# typically seen in these fits.
EXCLUDE_HALFWIDTH = 4.0


def _prep_window(spectrum_table, wavelength_min, wavelength_max):
    '''
    Slice spectrum_table (WAVE/FLUX/ERROR/SKY columns) to
    [wavelength_min, wavelength_max], drop non-finite FLUX/ERROR/SKY
    pixels, and floor ERROR to avoid zero/near-zero weights blowing up
    the fit (same 1e-6-of-peak floor lvm_gaussfit.py uses).

    Returns:
        x, y, sky, errors (arrays) -- errors is None if spectrum_table
        has no ERROR column. errors is None if spectrum_table has no
        ERROR column; used for the fit itself.
        x_full, y_full, sky_full (arrays) -- the same window WITHOUT the
        finite filter, always exactly wavelength_min<=WAVE<=wavelength_max
        long (a fixed length set only by the instrument's WAVE grid, the
        same for every fiber in a file) -- used for the diagnostic FITS
        cutouts, so stacking fibers never hits a ragged-array shape
        mismatch from one fiber having a masked/NaN pixel another
        doesn't.
    '''
    mask = (spectrum_table['WAVE'] >= wavelength_min) & (spectrum_table['WAVE'] <= wavelength_max)
    x_full = np.asarray(spectrum_table['WAVE'][mask])
    y_full = np.asarray(spectrum_table['FLUX'][mask])
    sky_full = np.asarray(spectrum_table['SKY'][mask])

    errors_full = np.asarray(spectrum_table['ERROR'][mask]) if 'ERROR' in spectrum_table.colnames else None

    finite = np.isfinite(y_full) & np.isfinite(sky_full)
    if errors_full is not None:
        finite = finite & np.isfinite(errors_full)

    x, y, sky = x_full[finite], y_full[finite], sky_full[finite]
    errors = errors_full[finite] if errors_full is not None else None
    if errors is not None:
        min_error_value = np.max(np.abs(y)) * 1e-6
        errors = np.maximum(errors, min_error_value)

    return x, y, sky, errors, x_full, y_full, sky_full


def _fit_wing_background(x, y, sky, errors, centers, half_width=EXCLUDE_HALFWIDTH):
    '''
    Fit const + scale*sky by weighted linear least squares, using only
    the "wing" pixels more than half_width from every one of centers --
    i.e. excluding the region(s) where the line(s) being fit are
    expected, so a real emission line (or a coincident feature in the
    sky template itself, e.g. geocoronal Halpha or Vela's own diffuse
    emission leaking into a nearby sky telescope's spectrum) cannot bias
    the background estimate. This is what fit_singlet_with_sky/
    fit_doublet_with_sky use to fix const/scale before fitting the line,
    instead of solving for them jointly with the line from the same
    pixels.

    Parameters:
        x, y, sky, errors: arrays -- the finite-filtered fit window (as
            returned by _prep_window). errors may be None.
        centers: list of float -- expected line center wavelength(s) to
            exclude around.
        half_width: float -- Angstroms excluded on each side of each
            center (default EXCLUDE_HALFWIDTH).

    Returns:
        const, scale, const_err, scale_err (floats). Raises ValueError
        if fewer than 3 wing pixels remain (window too narrow / too much
        excluded relative to its width for this to be well-determined).
    '''
    wing = np.ones(len(x), dtype=bool)
    for c in centers:
        wing &= np.abs(x - c) > half_width

    if np.count_nonzero(wing) < 3:
        raise ValueError('Not enough wing pixels outside %.1f Ang of %s to fit a background'
                          % (half_width, centers))

    xw, yw, skyw = x[wing], y[wing], sky[wing]
    w = 1.0 / errors[wing] ** 2 if errors is not None else np.ones_like(yw)

    A = np.column_stack([np.ones_like(skyw), skyw])
    AtWA = A.T @ (A * w[:, None])
    AtWy = A.T @ (yw * w)
    coef, resid, rank, _ = np.linalg.lstsq(AtWA, AtWy, rcond=None)
    const, scale = coef

    cov = np.linalg.inv(AtWA)
    n, k = len(yw), 2
    if n > k:
        chi2 = np.sum(w * (yw - A @ coef) ** 2)
        cov = cov * (chi2 / (n - k))
    const_err, scale_err = np.sqrt(np.diag(cov))

    return const, scale, const_err, scale_err


def fit_singlet_with_sky(spectrum_table, line, init_wavelength, init_fwhm, wavelength_min, wavelength_max,
                          fix_background=False):
    '''
    Fit a single Gaussian plus const + scale*SKY to one emission line.

    Parameters:
        spectrum_table (astropy.table.Table): WAVE (Ang), FLUX
            (erg/s/cm**2/Ang * 1e16, NOT sky-subtracted), ERROR, and SKY
            (same units as FLUX -- the chosen telescope's raw spectrum
            for this fiber) columns.
        line (str): line key, used to build output column names.
        init_wavelength (float): initial center guess (Ang).
        init_fwhm (float): initial FWHM guess (Ang).
        wavelength_min, wavelength_max (float): fit window (Ang).
        fix_background (bool): if False (default), const/scale are fit
            jointly with the line from the same pixels, as
            lvm_gaussfit.py's flat background is. If True, const/scale
            are instead determined first from the wings only
            (_fit_wing_background, excluding +/-EXCLUDE_HALFWIDTH Ang
            around init_wavelength) and held fixed while only
            flux/center/fwhm are fit to the line; this is what keeps a
            real feature in the sky template coincident with the line
            (geocoronal Halpha, Vela's own diffuse emission in a sky
            telescope, etc.) from being able to bias the line flux, but
            cannot help when that feature sits at the line's own
            wavelength in the template itself (see -wing on the CLI).

    Returns:
        qtab (astropy.table.Table): one row with flux_<line>,
            eflux_<line>, wave_<line>, ewave_<line>, fwhm_<line>,
            efwhm_<line>, back_<line>, eback_<line> (the const term),
            skysc_<line>, eskysc_<line> (the scale term), rmse_<line>,
            chi2_<line>.
        fit_table (astropy.table.Table): WAVE, FLUX, SKY_COMP (=const +
            scale*SKY), LINE_COMP (Gaussian alone), MOD (=SKY_COMP +
            LINE_COMP), for the diagnostic FITS output.
    '''
    def model(x, sky, flux, center, fwhm, const, scale):
        sigma = fwhm / 2.355
        amp = flux / (sigma * np.sqrt(2 * np.pi))
        return amp * np.exp(-0.5 * ((x - center) / sigma) ** 2) + const + scale * sky

    x, y, sky, errors, x_full, y_full, sky_full = _prep_window(spectrum_table, wavelength_min, wavelength_max)

    dx = x[1] - x[0]
    init_flux = np.sum(y * dx)

    if fix_background:
        init_const, init_scale, const_err, scale_err = _fit_wing_background(x, y, sky, errors, [init_wavelength])
    else:
        init_scale = 1.0
        init_const = np.median(y - init_scale * sky)

    gmodel = Model(model, independent_vars=['x', 'sky'])
    params = gmodel.make_params(flux=init_flux, center=init_wavelength, fwhm=init_fwhm,
                                 const=init_const, scale=init_scale)
    params['fwhm'].min = 0.1 * init_fwhm
    params['fwhm'].max = 5.0 * init_fwhm
    params['center'].min = wavelength_min
    params['center'].max = wavelength_max
    if fix_background:
        params['const'].vary = False
        params['scale'].vary = False

    if errors is not None:
        result = gmodel.fit(y, params, x=x, sky=sky, weights=1.0 / errors)
        patch_stderr_fraction(result, min_frac=0.01)
    else:
        result = gmodel.fit(y, params, x=x, sky=sky)

    p = result.params
    flux = p['flux'].value
    flux_err = p['flux'].stderr if p['flux'].stderr is not None else np.nan
    wave = p['center'].value
    wave_err = p['center'].stderr if p['center'].stderr is not None else np.nan
    fwhm = p['fwhm'].value
    fwhm_err = p['fwhm'].stderr if p['fwhm'].stderr is not None else np.nan
    const = p['const'].value
    scale = p['scale'].value
    if fix_background:
        # lmfit reports stderr=None for a fixed (vary=False) parameter --
        # use the wing-fit's own uncertainty instead.
        pass
    else:
        const_err = p['const'].stderr if p['const'].stderr is not None else np.nan
        scale_err = p['scale'].stderr if p['scale'].stderr is not None else np.nan

    y_fit = result.best_fit
    rmse = np.sqrt(np.mean((y - y_fit) ** 2))
    reduced_chi2 = result.redchi

    sigma = fwhm / 2.355
    amp = flux / (sigma * np.sqrt(2 * np.pi))
    # Diagnostic components are evaluated on the FULL (fixed-length)
    # window, not the finite-filtered fit arrays -- every fiber's cutout
    # must be the same length for _write_diagnostic_fits to stack them.
    line_comp = amp * np.exp(-0.5 * ((x_full - wave) / sigma) ** 2)
    sky_comp = const + scale * sky_full

    col_names = [f'{p_}_{line}' for p_ in
                 ['flux', 'eflux', 'wave', 'ewave', 'fwhm', 'efwhm', 'back', 'eback',
                  'skysc', 'eskysc', 'rmse', 'chi2']]
    qtab = Table(names=col_names)
    qtab.add_row([flux, flux_err, wave, wave_err, fwhm, fwhm_err, const, const_err,
                  scale, scale_err, rmse, reduced_chi2])

    fit_table = Table([x_full, y_full, sky_comp, line_comp, sky_comp + line_comp],
                       names=('WAVE', 'FLUX', 'SKY_COMP', 'LINE_COMP', 'MOD'))

    return qtab, fit_table


def fit_doublet_with_sky(spectrum_table, line, init_wavelength1, init_wavelength2, init_fwhm,
                          wavelength_min, wavelength_max, fix_background=False):
    '''
    Fit two Gaussians (shared FWHM, independent -- freely ratioed --
    fluxes) plus const + scale*SKY to a close emission-line doublet.

    Parameters: same convention as fit_singlet_with_sky (including
        fix_background), with two initial center guesses
        (init_wavelength1/2) -- both are excluded from the wing fit.

    Returns:
        qtab (astropy.table.Table): one row with flux_<line>_a/_b,
            eflux_<line>_a/_b, wave_<line>_a/_b, ewave_<line>_a/_b,
            fwhm_<line>_a (shared), efwhm_<line>_a, back_<line>_ab,
            eback_<line>_ab, skysc_<line>_ab, eskysc_<line>_ab,
            rmse_<line>_ab, chi2_<line>_ab -- mirrors
            lvm_gaussfit.fit_double_gaussian_to_spectrum's naming.
        fit_table (astropy.table.Table): WAVE, FLUX, SKY_COMP, LINE_COMP,
            MOD, for the diagnostic FITS output.
    '''
    def model(x, sky, flux1, flux2, center1, center2, fwhm, const, scale):
        sigma = fwhm / 2.355
        amp1 = flux1 / (sigma * np.sqrt(2 * np.pi))
        amp2 = flux2 / (sigma * np.sqrt(2 * np.pi))
        g1 = amp1 * np.exp(-0.5 * ((x - center1) / sigma) ** 2)
        g2 = amp2 * np.exp(-0.5 * ((x - center2) / sigma) ** 2)
        return g1 + g2 + const + scale * sky

    x, y, sky, errors, x_full, y_full, sky_full = _prep_window(spectrum_table, wavelength_min, wavelength_max)

    dx = x[1] - x[0]
    if fix_background:
        init_const, init_scale, const_err, scale_err = _fit_wing_background(
            x, y, sky, errors, [init_wavelength1, init_wavelength2])
    else:
        init_scale = 1.0
        init_const = np.median(y - init_scale * sky)
    total_flux = np.sum((y - init_const - init_scale * sky) * dx)
    init_flux1 = total_flux / 2
    init_flux2 = total_flux / 2

    gmodel = Model(model, independent_vars=['x', 'sky'])
    params = gmodel.make_params(flux1=init_flux1, flux2=init_flux2,
                                 center1=init_wavelength1, center2=init_wavelength2,
                                 fwhm=init_fwhm, const=init_const, scale=init_scale)
    params['fwhm'].min = 0.1 * init_fwhm
    params['fwhm'].max = 5.0 * init_fwhm
    params['center1'].min = wavelength_min
    params['center2'].max = wavelength_max
    if init_wavelength1 < init_wavelength2:
        params['center1'].max = init_wavelength2
        params['center2'].min = init_wavelength1
    else:
        params['center1'].max = init_wavelength1
        params['center2'].min = init_wavelength2
    if fix_background:
        params['const'].vary = False
        params['scale'].vary = False

    if errors is not None:
        result = gmodel.fit(y, params, x=x, sky=sky, weights=1.0 / errors)
        patch_stderr_fraction(result, min_frac=0.01)
    else:
        result = gmodel.fit(y, params, x=x, sky=sky)

    p = result.params
    flux1, flux1_err = p['flux1'].value, p['flux1'].stderr or np.nan
    flux2, flux2_err = p['flux2'].value, p['flux2'].stderr or np.nan
    wave1, wave1_err = p['center1'].value, p['center1'].stderr or np.nan
    wave2, wave2_err = p['center2'].value, p['center2'].stderr or np.nan
    fwhm, fwhm_err = p['fwhm'].value, p['fwhm'].stderr or np.nan
    const = p['const'].value
    scale = p['scale'].value
    if not fix_background:
        # lmfit reports stderr=None for a fixed (vary=False) parameter --
        # when fix_background, const_err/scale_err already came from
        # _fit_wing_background above and should be left alone.
        const_err = p['const'].stderr or np.nan
        scale_err = p['scale'].stderr or np.nan

    y_fit = result.best_fit
    rmse = np.sqrt(np.mean((y - y_fit) ** 2))
    reduced_chi2 = result.redchi

    sigma = fwhm / 2.355
    amp1 = flux1 / (sigma * np.sqrt(2 * np.pi))
    amp2 = flux2 / (sigma * np.sqrt(2 * np.pi))
    # Diagnostic components are evaluated on the FULL (fixed-length)
    # window, not the finite-filtered fit arrays -- every fiber's cutout
    # must be the same length for _write_diagnostic_fits to stack them.
    line_comp = (amp1 * np.exp(-0.5 * ((x_full - wave1) / sigma) ** 2)
                 + amp2 * np.exp(-0.5 * ((x_full - wave2) / sigma) ** 2))
    sky_comp = const + scale * sky_full

    col_names = [f'flux_{line}_a', f'eflux_{line}_a', f'flux_{line}_b', f'eflux_{line}_b',
                 f'wave_{line}_a', f'ewave_{line}_a', f'wave_{line}_b', f'ewave_{line}_b',
                 f'fwhm_{line}_a', f'efwhm_{line}_a',
                 f'back_{line}_ab', f'eback_{line}_ab', f'skysc_{line}_ab', f'eskysc_{line}_ab',
                 f'rmse_{line}_ab', f'chi2_{line}_ab']
    qtab = Table(names=col_names)
    qtab.add_row([flux1, flux1_err, flux2, flux2_err, wave1, wave1_err, wave2, wave2_err,
                  fwhm, fwhm_err, const, const_err, scale, scale_err, rmse, reduced_chi2])

    fit_table = Table([x_full, y_full, sky_comp, line_comp, sky_comp + line_comp],
                       names=('WAVE', 'FLUX', 'SKY_COMP', 'LINE_COMP', 'MOD'))

    return qtab, fit_table


def do_one(spectrum_table, vel=0.0, fix_background=False):
    '''
    Fit every line group in DOUBLETS+SINGLETS (matching lvm_gaussfit.py's
    full line set) for one fiber's spectrum_table (WAVE/FLUX/ERROR/SKY).

    Parameters:
        spectrum_table (astropy.table.Table): one fiber's WAVE, FLUX,
            ERROR, SKY columns (all scaled by 1e16, FLUX not sky-subtracted).
        vel (float): velocity offset (km/s) applied to initial line
            centers, same convention as lvm_gaussfit.py.
        fix_background (bool): passed through to fit_singlet_with_sky/
            fit_doublet_with_sky -- default False (joint fit); True
            selects the wing-only background instead (-wing on the CLI).

    Returns:
        qtab (astropy.table.Table): one row, every line group's columns
            hstacked. Empty Table if every line group failed.
        fit_tables (dict): {line_group_name: fit_table}, for whichever
            groups fit successfully (used to build the diagnostic FITS
            cutouts in do_all).
    '''
    zz = 1. + vel / 3e5
    records = []
    fit_tables = {}

    for line, w1, w2, fwhm, wmin, wmax in DOUBLETS:
        try:
            qtab, fit_table = fit_doublet_with_sky(spectrum_table, line, zz * w1, zz * w2, fwhm, wmin, wmax,
                                                    fix_background=fix_background)
            records.append(qtab)
            fit_tables[line] = fit_table
        except Exception as e:
            print(f'Fitting {line}: An exception occurred: {e}')

    for line, w0, fwhm, wmin, wmax in SINGLETS:
        try:
            qtab, fit_table = fit_singlet_with_sky(spectrum_table, line, zz * w0, fwhm, wmin, wmax,
                                                    fix_background=fix_background)
            records.append(qtab)
            fit_tables[line] = fit_table
        except Exception as e:
            print(f'Fitting {line}: An exception occurred: {e}')

    try:
        ztab = hstack(records)
    except Exception:
        print('Nothing fit for this spectrum')
        return Table(), {}

    return ztab, fit_tables


def _pick_sky_telescope(header, mode='far', override=None):
    '''
    Decide which sky telescope's raw spectrum to use for a whole
    exposure, from the SCI_SKYE_SEP/SCI_SKYW_SEP header keywords,
    unless override is given.

    Parameters:
        header (astropy.io.fits.Header): CFrame primary header.
        mode (str): 'far' (default) picks the LARGER separation --
            avoids a large, diffuse target like Vela leaking into a sky
            telescope that happens to sit close to the field. 'near'
            picks the smaller separation instead.
        override (str or None): 'EAST' or 'WEST' to force a choice,
            regardless of mode (highest precedence).

    Returns:
        tel (str): 'EAST' or 'WEST'.
        sep_e, sep_w (float): the two separations (deg), for provenance.
    '''
    sep_e = header['SKY SCI_SKYE_SEP']
    sep_w = header['SKY SCI_SKYW_SEP']
    if override is not None:
        tel = override.upper()
    elif mode == 'near':
        tel = 'EAST' if sep_e <= sep_w else 'WEST'
    elif mode == 'far':
        tel = 'EAST' if sep_e >= sep_w else 'WEST'
    else:
        raise ValueError("mode must be 'far' or 'near', got %r" % mode)
    return tel, sep_e, sep_w


def do_all(filename, vel=0.0, outname='', tel_mode='far', tel_override=None, fix_background=False):
    '''
    Fit all science fibers in one lvmCFrame file and write the results
    table and diagnostic FITS.

    Reads FLUX (erg/s/cm**2/Ang, NOT sky-subtracted), one of
    SKY_EAST/SKY_WEST (per _pick_sky_telescope), and IVAR from the
    file; multiplies FLUX/SKY/ERROR by 1e16 for numerical convenience
    (scale is unaffected, since both FLUX and SKY get the same factor).
    Science fibers are selected via lvm_gaussfit.scifib on the file's
    SLITMAP (same convention as lvmSFrame).

    Parameters:
        filename (str): path to an ``lvmCFrame-*.fits`` file.
        vel (float): velocity offset (km/s), passed to do_one.
        outname (str): output rootname; defaults to filename's basename
            with the .fits suffix stripped.
        tel_mode (str): 'far' (default) or 'near', passed to
            _pick_sky_telescope; ignored if tel_override is given.
        tel_override (str or None): force 'EAST' or 'WEST' instead of
            the automatic telescope choice.
        fix_background (bool): passed through to do_one (default False:
            joint fit; True selects the wing-only background instead).

    Returns:
        results (astropy.table.Table): the per-fiber results table (also
            written to <outname>.skyfit.txt).  Also writes
            <outname>.skyfit_model.fits (see module Notes).
    '''
    try:
        x = fits.open(filename)
    except Exception:
        print('Error: Could not open ', filename)
        return

    wave = x['WAVE'].data
    flux = x['FLUX'].data * 1e16
    with np.errstate(divide='ignore', invalid='ignore'):
        error = 1. / np.sqrt(x['IVAR'].data) * 1e16
        error = np.where(np.isfinite(error), error, 1e30)

    tel, sep_e, sep_w = _pick_sky_telescope(x['PRIMARY'].header, mode=tel_mode, override=tel_override)
    sky = x[f'SKY_{tel}'].data * 1e16
    print('Using SKY_%s as the sky template (SCI_SKYE_SEP=%.2f, SCI_SKYW_SEP=%.2f)' % (tel, sep_e, sep_w))

    slittab = Table(x['SLITMAP'].data)
    good = scifib(slittab, 'science')

    records = []
    fit_tables_by_line = {name: [] for name, *_ in DOUBLETS + SINGLETS}
    fiberids, ras, decs = [], [], []

    for i in range(len(good)):
        j = good['fiberid'][i] - 1
        one_spec = Table([wave, flux[j], error[j], sky[j]], names=['WAVE', 'FLUX', 'ERROR', 'SKY'])
        if check_for_nan(flux[j]) or check_for_nan(sky[j]):
            print('Too many nans for fiber %d at %.2f %.2f' % (good['fiberid'][i], good['ra'][i], good['dec'][i]))
            continue

        rtab, fit_tables = do_one(spectrum_table=one_spec, vel=vel, fix_background=fix_background)
        if len(rtab) == 0:
            print('Nothing fit for fiber %d at %.2f %.2f' % (good['fiberid'][i], good['ra'][i], good['dec'][i]))
            continue

        rtab['fiberid'] = good['fiberid'][i]
        rtab['ra'] = good['ra'][i]
        rtab['dec'] = good['dec'][i]
        records.append(rtab)
        fiberids.append(good['fiberid'][i])
        ras.append(good['ra'][i])
        decs.append(good['dec'][i])
        for name in fit_tables_by_line:
            fit_tables_by_line[name].append(fit_tables.get(name))

    if not records:
        print('Nothing fit in this file')
        return

    results = vstack(records)

    if outname == '':
        outname = os.path.basename(filename).replace('.fits', '')

    for one in results.colnames:
        if one.count('flux'):
            results[one].format = '.3e'
        elif one.count('wave'):
            results[one].format = '.3f'
        elif one.count('fwhm'):
            results[one].format = '.3f'
        elif one.count('rmse'):
            results[one].format = '.3e'
        elif one.count('skysc'):
            results[one].format = '.4f'
        elif one.count('back'):
            results[one].format = '.3e'
        elif one.count('chi'):
            results[one].format = '.3f'
        elif one.count('ra'):
            results[one].format = '.5f'
        elif one.count('dec'):
            results[one].format = '.5f'
        elif one.count('fiberid'):
            results[one].format = 'd'

    table_out = outname + '.skyfit.txt'
    results.write(table_out, format='ascii.fixed_width_two_line', overwrite=True)
    print('Wrote', table_out)

    _write_diagnostic_fits(outname + '.skyfit_model.fits', fit_tables_by_line,
                            fiberids, ras, decs, filename, tel, sep_e, sep_w, x['PRIMARY'].header)

    return results


def _write_diagnostic_fits(outfile, fit_tables_by_line, fiberids, ras, decs,
                            source_filename, tel, sep_e, sep_w, src_header):
    '''
    Assemble and write the compact per-line-window diagnostic FITS
    described in the module Notes: WAVE_<LINE>/FLUX_<LINE>/SKY_<LINE>/
    LINE_<LINE>/MOD_<LINE> for each line group in DOUBLETS+SINGLETS,
    plus a FIBERID/RA/DEC table, rows in the same order as the results
    table.
    '''
    hdus = [fits.PrimaryHDU()]
    hdus[0].header['SRCFILE'] = os.path.basename(source_filename)
    hdus[0].header['SKYTEL'] = ('SKY_%s' % tel, 'sky telescope used as template')
    hdus[0].header['SEP_E'] = (sep_e, 'science-SKY_EAST separation (deg)')
    hdus[0].header['SEP_W'] = (sep_w, 'science-SKY_WEST separation (deg)')
    if 'DRPVER' in src_header:
        hdus[0].header['DRPVER'] = src_header['DRPVER']

    for name in fit_tables_by_line:
        tables = [t for t in fit_tables_by_line[name] if t is not None]
        if not tables:
            continue
        wave_grid = np.asarray(tables[0]['WAVE'])
        flux_stack = np.array([np.asarray(t['FLUX']) for t in fit_tables_by_line[name]
                                if t is not None])
        sky_stack = np.array([np.asarray(t['SKY_COMP']) for t in fit_tables_by_line[name]
                               if t is not None])
        line_stack = np.array([np.asarray(t['LINE_COMP']) for t in fit_tables_by_line[name]
                                if t is not None])
        mod_stack = np.array([np.asarray(t['MOD']) for t in fit_tables_by_line[name]
                               if t is not None])

        upper = name.upper()
        hdus.append(fits.ImageHDU(wave_grid.astype('f4'), name=f'WAVE_{upper}'))
        hdus.append(fits.ImageHDU(flux_stack.astype('f4'), name=f'FLUX_{upper}'))
        hdus.append(fits.ImageHDU(sky_stack.astype('f4'), name=f'SKY_{upper}'))
        hdus.append(fits.ImageHDU(line_stack.astype('f4'), name=f'LINE_{upper}'))
        hdus.append(fits.ImageHDU(mod_stack.astype('f4'), name=f'MOD_{upper}'))

    id_table = Table([fiberids, ras, decs], names=('fiberid', 'ra', 'dec'))
    hdus.append(fits.BinTableHDU(id_table.as_array(), name='FIBERID'))

    fits.HDUList(hdus).writeto(outfile, overwrite=True)
    print('Wrote', outfile)


def do_individual(filenames, vel=0.0, outname='', tel_mode='far', tel_override=None,
                   fix_background=False, do_all_func=do_all):
    '''
    Loop do_all over multiple lvmCFrame files (a batch driver over whole
    FITS files, not lvm_gaussfit.py's ascii-table do_individual).

    Parameters:
        filenames (list of str): ``lvmCFrame-*.fits`` files.
        vel (float): velocity offset, passed through to do_all.
        outname (str): passed through to do_all (each file still gets
            its own output, so this only matters if you want a shared
            rootname override -- normally leave '').
        tel_mode (str): passed through to do_all.
        tel_override (str or None): passed through to do_all.
        fix_background (bool): passed through to do_all.
        do_all_func: override for testing.

    Returns:
        None
    '''
    for f in filenames:
        do_all_func(f, vel=vel, outname=outname, tel_mode=tel_mode, tel_override=tel_override,
                    fix_background=fix_background)


def steer(argv):
    '''
    Parse the command line and run do_all/do_individual.
    '''
    outname = ''
    vel = 0.0
    tel_mode = 'far'
    tel_override = None
    fix_background = False
    files = []

    i = 1
    while i < len(argv):
        if argv[i][0:2] == '-h':
            print(_usage_from_doc(__doc__))
            return
        elif argv[i][0:4] == '-out':
            i += 1
            outname = argv[i]
        elif argv[i] == '-v':
            i += 1
            vel = float(argv[i])
        elif argv[i] == '-near':
            tel_mode = 'near'
        elif argv[i] == '-tel':
            i += 1
            tel_override = argv[i]
        elif argv[i] == '-wing':
            fix_background = True
        elif argv[i][0] == '-':
            print('Unknown option :', argv)
            return
        elif argv[i].count('.fits'):
            files.append(argv[i])
        else:
            print('Unknown option :', argv)
            return
        i += 1

    if not files:
        print('Error: no lvmCFrame file(s) given')
        print(_usage_from_doc(__doc__))
        return

    if len(files) == 1:
        do_all(files[0], vel=vel, outname=outname, tel_mode=tel_mode, tel_override=tel_override,
               fix_background=fix_background)
    else:
        do_individual(files, vel=vel, outname=outname, tel_mode=tel_mode, tel_override=tel_override,
                       fix_background=fix_background)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_usage_from_doc(__doc__))
