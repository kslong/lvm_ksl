#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Fit the Mappings-model emission-line set expected in SNRs (see
data/mappings_snr_lines.txt) to a single sky-subtracted LVM spectrum or
RSS file, the way lvm_gaussfit.py does for its smaller, hardcoded line
list -- but with lines that are too close together to fit independently
(given plausible line broadening) fit jointly instead, sharing one local
background.

Command line usage (if any)::

    usage: lvm_snrfit.py [-h] [-lmc] [-smc] [-v vel] [-stype SOURCE]
                          [-out root] [-plot] [-lines file.txt]
                          [-min_sig N] filename ...

    where

    -h              print this documentation and exit
    -lmc or -smc    apply a fixed velocity offset for fitting (Simbad)
    -v vel          apply a velocity offset that the user specifies
    -stype SOURCE or BACK -- used only for spectra given as ascii text
                    files with SOURCE_FLUX/BACK_FLUX columns
    -out root       root name for the output file
    -plot           save a per-line fit-quality plot for each spectrum
    -lines file.txt reference line list (default data/mappings_snr_lines.txt)
    -min_sig N      significance threshold for keeping a blend component
                    (default 2.0) -- see fit_blend_to_spectrum
    filename        an SFrame-compatible FITS file, or one or more ascii
                    tables with WAVE/FLUX[/ERROR] columns

Description:

    data/mappings_snr_lines.txt gives every line's rest wavelength as
    predicted by a Mappings shock model -- more accurate than the DAP's
    own (approximate) wavelengths -- plus a `group` column tagging lines
    that must be fit together.  Four groups need this (see dap.rst-style
    reasoning in the design discussion that produced this file):

        hei_hI          HeI 3888.64 / HI(H8) 3889.06  (0.4 A apart)
        neiii_hepsilon  [NeIII] 3967.47 / Hepsilon 3970.08 (2.6 A apart)
        ni              [NI] 5197.90 / [NI] 5200.26  (2.4 A apart)
        oii7320_caii    [OII] 7319.99 / [CaII] 7323.89 (5.0 A apart)

    plus the pre-existing `oii` doublet (3726.03/3728.82), fit exactly as
    lvm_gaussfit.py already does via fit_double_gaussian_to_spectrum.  All
    other lines are independent singlets, fit with
    lvm_gaussfit.fit_gaussian_to_spectrum unchanged.

    Why a joint fit: each independent single-Gaussian fit estimates its
    own local background from a window that, for these close pairs, would
    include the neighbour's wing -- biasing both.  fit_blend_to_spectrum
    below fits both lines at once with one shared background, independent
    centers (bounded to each line's own Mappings wavelength by a
    velocity-derived window), and either a shared or independent FWHM per
    group depending on how blended the pair actually is once broadening
    is accounted for (see WIDTH_SHARING below).

    Line widths are bounded by physical velocity rather than an arbitrary
    Angstrom scaling: FWHM_min from V_INSTR (typical instrumental width),
    FWHM_max from V_MAX (the broadest line judged plausible, e.g. shocked
    gas).  A component of a blend fit that comes back statistically
    insignificant (flux/eflux below -min_sig) is dropped and the group is
    refit with the remaining line(s) as an ordinary singlet, so a
    non-detection can't corrupt the shared background/width for a real
    line sitting next to it.

Units:

    Same convention as lvm_gaussfit.py: FLUX is multiplied by 1e16 before
    fitting, so flux_<name>/back_<name> columns are in erg/s/cm**2 * 1e16
    (divide by 1e16 for physical units); wave_<name>/fwhm_<name> are in
    Angstroms.

Primary routines:

    do_one, fit_blend_to_spectrum

Notes::

    File I/O (SFrame/ascii reading, fiber selection, batch output,
    plotting) is not duplicated here -- do_all/do_individual/analyze/
    plot_one/plot_all/save_fit/clean/scifib/check_for_nan are imported
    from lvm_gaussfit.py, which now accepts a do_one_func argument so
    this module's do_one can be run through that same machinery.

History::

    260823 ksl Coding begun

'''

import os
import sys

import numpy as np
from astropy.io import ascii
from astropy.table import Table, hstack
from lmfit import Model

from lvm_gaussfit import (
    _usage_from_doc, fit_gaussian_to_spectrum, fit_double_gaussian_to_spectrum,
    save_fit, clean, scifib, check_for_nan, do_all, do_individual, analyze,
    plot_one, plot_all,
)


C_KMS = 299792.458

# Physically-motivated FWHM bounds: typical instrumental width and the
# broadest line judged plausible (e.g. a shocked SNR knot).  Overridable
# per-call, not currently exposed on the command line.
V_INSTR = 80.0
V_MAX = 200.0

DEFAULT_LINES_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'mappings_snr_lines.txt')

# Groups fit with a single shared FWHM (physically similar / heavily
# blended pairs); everything else uses independent FWHM per component
# (well enough separated, different species).  'oii' is handled
# separately -- it already has a dedicated, working implementation in
# lvm_gaussfit.py.
SHARED_FWHM_GROUPS = {'hei_hI', 'neiii_hepsilon', 'ni'}

# Significance threshold (flux/eflux) for keeping a blend component --
# see fit_blend_to_spectrum.  Overridable with -min_sig.
MIN_SIGNIFICANCE = 2.0

# Pairs closer than roughly their own FWHM: two free centers are
# degenerate with two free fluxes (see fit_blend_to_spectrum's
# tie_centers).  Currently just HeI 3888.65/HI 3889.05 (0.4 A apart).
TIE_CENTER_GROUPS = {'hei_hI'}

# Fixed (un-redshifted -- these are telluric/sky, not object) wavelength
# windows to exclude before fitting specific singlets, matching
# lvm_gaussfit.do_one: [OI] 6300/6364 sit on sky-subtraction-residual-
# prone airglow lines, "seldom subtracted correctly".
SKY_EXCLUDE = {'oi_a': (6298.0, 6303.0), 'oi_b': (6362.0, 6365.0)}


def fwhm_bounds(wave, v_min=V_INSTR, v_max=V_MAX):
    '''
    Velocity-derived (min, max) FWHM bounds in Angstroms at wavelength wave.
    '''
    return wave * v_min / C_KMS, wave * v_max / C_KMS


def window_halfwidth(wave, n_vmax=4.0, floor=8.0, cap=None):
    '''
    Half-width (Angstroms) of a fitting window around wave: n_vmax times
    the V_MAX-broadened FWHM, with a floor so narrow-window numerical
    issues don't show up for lines near the blue end.  If cap is given,
    the result is also capped there -- see line_window_caps(), which
    prevents one line's fitting window from reaching into a neighbour's
    rest wavelength (that neighbour's wing would then bias this line's
    local-background estimate, the same problem the blend groups solve
    for, just for lines close enough that their independent-singlet
    windows overlap without the lines themselves being blended).
    '''
    _, fmax = fwhm_bounds(wave)
    half = max(floor, n_vmax * fmax)
    if cap is not None:
        half = min(half, cap)
    return half


def line_window_caps(tab):
    '''
    For every row in the line table, the max window half-width (Angstroms)
    that keeps its fitting window from reaching a line outside its own
    blend group: half the distance to the nearest such neighbour, or None
    if none is close enough to matter.  Computed on rest wavelengths
    (velocity shifts move all lines together, so relative spacing is
    unaffected).  Keyed by gauss_name.
    '''
    waves = np.asarray(tab['Wave'], dtype=float)
    groups = np.asarray(tab['group'])
    caps = {}
    for i, row in enumerate(tab):
        own_group = row['group']
        other = np.array([j for j in range(len(tab))
                           if j != i and not (own_group != '-' and groups[j] == own_group)])
        if len(other) == 0:
            caps[row['gauss_name']] = None
            continue
        min_dist = np.min(np.abs(waves[other] - waves[i]))
        caps[row['gauss_name']] = min_dist / 2.0
    return caps


def load_lines(filename=None):
    '''
    Read the SNR line-list reference table (default
    data/mappings_snr_lines.txt).  Returns the astropy Table.
    '''
    path = filename if filename else DEFAULT_LINES_FILE
    return ascii.read(str(path))


def get_groups(tab):
    '''
    Split a line table (as read by load_lines) into a dict of
    group_name -> list of rows, preserving the '-' tag for independent
    singlets as its own key.
    '''
    groups = {}
    for row in tab:
        groups.setdefault(row['group'], []).append(row)
    return groups


def _significant(param, min_significance):
    '''
    True if an lmfit flux Parameter's fitted value is significant relative
    to its uncertainty.  False (not significant) if stderr is missing,
    zero, non-finite, or the significance ratio is below threshold.
    '''
    if param.stderr is None or param.stderr == 0 or not np.isfinite(param.stderr):
        return False
    return abs(param.value) / param.stderr >= min_significance


def fit_blend_to_spectrum(spectrum_table, group_name, name1, wave1, name2, wave2,
                           wavelength_min, wavelength_max, shared_fwhm=True,
                           tie_centers=False, min_significance=2.0):
    '''
    Jointly fit two nearby emission lines with one shared local background.

    name1/wave1, name2/wave2 are the (already velocity-shifted) component
    names and initial center wavelengths.  FWHM is a single shared
    parameter if shared_fwhm else independent per component, each bounded
    by fwhm_bounds().

    By default centers are independently fit, each bounded to its own
    initial wavelength by a velocity-derived window (V_MAX).  This is
    only well-conditioned if the two lines are separated by at least
    roughly their own FWHM -- for a pair closer than that (e.g. HeI
    3888.65/HI 3889.05, 0.4 A apart), two free centers plus two free
    fluxes are degenerate: the fit can trade flux between components
    while barely changing chi-square, producing unphysical (even
    negative) individual fluxes.  Set tie_centers=True for such a pair:
    center2 is then locked to center1 + (wave2-wave1), so only one
    systemic shift is fit and both lines move together, leaving flux1/
    flux2 as the only free split between them.

    If either component's fitted flux is not significant (see
    _significant), that component is dropped and the group is refit as a
    single Gaussian (via lvm_gaussfit.fit_gaussian_to_spectrum) using only
    the surviving line -- the dropped line's columns are returned as NaN
    rather than a spurious joint-fit value, so a non-detection can't bias
    the shared background/width used for the real line.

    Returns (qtab, fit_table), qtab a single-row Table with
    flux_<name>/eflux_<name>/wave_<name>/ewave_<name>/fwhm_<name>/
    efwhm_<name> for each of name1/name2, and back_<group_name>/
    eback_<group_name>/rmse_<group_name>/chi2_<group_name> for the shared
    fit; fit_table has the WAVE/FLUX/Fit[/ERROR] arrays actually fit.
    '''

    def blend_gaussian(x, flux1, flux2, center1, center2, fwhm1, fwhm2, background):
        sigma1 = fwhm1 / 2.355
        amp1 = flux1 / (sigma1 * np.sqrt(2 * np.pi))
        sigma2 = fwhm2 / 2.355
        amp2 = flux2 / (sigma2 * np.sqrt(2 * np.pi))
        g1 = amp1 * np.exp(-0.5 * ((x - center1) / sigma1) ** 2)
        g2 = amp2 * np.exp(-0.5 * ((x - center2) / sigma2) ** 2)
        return g1 + g2 + background

    mask = (spectrum_table['WAVE'] >= wavelength_min) & (spectrum_table['WAVE'] <= wavelength_max)
    x = spectrum_table['WAVE'][mask]
    y = spectrum_table['FLUX'][mask]

    if 'ERROR' in spectrum_table.colnames:
        errors = spectrum_table['ERROR'][mask]
    else:
        errors = None

    finite_mask = np.isfinite(y)
    if errors is not None:
        finite_mask = finite_mask & np.isfinite(errors)
    x = x[finite_mask]
    y = y[finite_mask]
    if errors is not None:
        errors = errors[finite_mask]
        min_error_value = np.max(np.abs(y)) * 1e-6
        errors = np.maximum(errors, min_error_value)

    dx = x[1] - x[0]
    total_flux = np.sum(y * dx)
    init_flux1 = total_flux / 2
    init_flux2 = total_flux / 2
    init_background = np.median(y)

    fwhm1_min, fwhm1_max = fwhm_bounds(wave1)
    fwhm2_min, fwhm2_max = fwhm_bounds(wave2)
    shift1 = wave1 * V_MAX / C_KMS
    shift2 = wave2 * V_MAX / C_KMS

    gmodel = Model(blend_gaussian)
    params = gmodel.make_params(
        flux1=init_flux1, flux2=init_flux2,
        center1=wave1, center2=wave2,
        fwhm1=fwhm1_min, fwhm2=fwhm2_min,
        background=init_background,
    )
    params['fwhm1'].min = fwhm1_min
    params['fwhm1'].max = fwhm1_max
    params['fwhm2'].min = fwhm2_min
    params['fwhm2'].max = fwhm2_max
    if shared_fwhm:
        params['fwhm2'].set(expr='fwhm1')
    params['center1'].min = max(wavelength_min, wave1 - shift1)
    params['center1'].max = min(wavelength_max, wave1 + shift1)
    if tie_centers:
        params['center2'].set(expr='center1 + %r' % (wave2 - wave1))
    else:
        params['center2'].min = max(wavelength_min, wave2 - shift2)
        params['center2'].max = min(wavelength_max, wave2 + shift2)

    if errors is not None:
        result = gmodel.fit(y, params, x=x, weights=1.0 / errors)
    else:
        result = gmodel.fit(y, params, x=x)

    sig1 = _significant(result.params['flux1'], min_significance)
    sig2 = _significant(result.params['flux2'], min_significance)

    # Even when both fluxes individually pass the significance test, a
    # heavily overlapping pair (e.g. tie_centers cases) can still have the
    # fit trading flux between components almost freely -- caught by a
    # near-total anti-correlation between flux1 and flux2 that per-
    # parameter stderr alone doesn't reveal.  Keep only the larger (more
    # trustworthy) component in that case.
    correl = result.params['flux1'].correl
    degenerate = bool(correl) and abs(correl.get('flux2', 0.0)) >= 0.95

    if degenerate or not (sig1 and sig2):
        # One (or both) components are not significant, or the split
        # between them is degenerate; fall back to an independent singlet
        # fit for whichever line(s) survive rather than trust a joint fit
        # that can't reliably separate the two.
        if degenerate:
            if abs(result.params['flux1'].value) >= abs(result.params['flux2'].value):
                survivors = [(name1, wave1)]
            else:
                survivors = [(name2, wave2)]
        else:
            survivors = []
            if sig1:
                survivors.append((name1, wave1))
            if sig2:
                survivors.append((name2, wave2))

        if len(survivors) == 1:
            sname, swave = survivors[0]
            half = window_halfwidth(swave)
            init_fwhm = fwhm_bounds(swave)[0]
            sqtab, sfit = fit_gaussian_to_spectrum(
                spectrum_table, line=sname, init_wavelength=swave, init_fwhm=init_fwhm,
                wavelength_min=max(wavelength_min, swave - half),
                wavelength_max=min(wavelength_max, swave + half))
            other_name = name2 if sname == name1 else name1
            for suffix in ['flux', 'eflux', 'wave', 'ewave', 'fwhm', 'efwhm']:
                sqtab['%s_%s' % (suffix, other_name)] = np.nan
            return sqtab, sfit
        # Neither component is significant (or both survived after all,
        # which only happens if min_significance<=0): fall through and
        # report the joint fit as-is -- large uncertainties already flag
        # it as a non-detection to downstream code.

    flux1 = result.params['flux1'].value
    flux1_err = result.params['flux1'].stderr if result.params['flux1'].stderr is not None else np.nan
    flux2 = result.params['flux2'].value
    flux2_err = result.params['flux2'].stderr if result.params['flux2'].stderr is not None else np.nan

    wave1_fit = result.params['center1'].value
    wave1_err = result.params['center1'].stderr if result.params['center1'].stderr is not None else np.nan
    wave2_fit = result.params['center2'].value
    wave2_err = result.params['center2'].stderr if result.params['center2'].stderr is not None else np.nan

    fwhm1_fit = result.params['fwhm1'].value
    fwhm1_err = result.params['fwhm1'].stderr if result.params['fwhm1'].stderr is not None else np.nan
    fwhm2_fit = result.params['fwhm2'].value
    fwhm2_err = result.params['fwhm2'].stderr if result.params['fwhm2'].stderr is not None else np.nan

    background = result.params['background'].value
    back_err = result.params['background'].stderr if result.params['background'].stderr is not None else np.nan

    y_fit = result.best_fit
    rmse = np.sqrt(np.mean((y - y_fit) ** 2))
    reduced_chi2 = result.redchi

    fit_table = Table([x, y, y_fit], names=('WAVE', 'FLUX', 'Fit'))
    if errors is not None:
        fit_table['ERROR'] = errors

    col_names = (
        ['flux_%s' % name1, 'eflux_%s' % name1, 'wave_%s' % name1, 'ewave_%s' % name1,
         'fwhm_%s' % name1, 'efwhm_%s' % name1,
         'flux_%s' % name2, 'eflux_%s' % name2, 'wave_%s' % name2, 'ewave_%s' % name2,
         'fwhm_%s' % name2, 'efwhm_%s' % name2,
         'back_%s' % group_name, 'eback_%s' % group_name,
         'rmse_%s' % group_name, 'chi2_%s' % group_name])
    qtab = Table(names=col_names)
    qtab.add_row([flux1, flux1_err, wave1_fit, wave1_err, fwhm1_fit, fwhm1_err,
                  flux2, flux2_err, wave2_fit, wave2_err, fwhm2_fit, fwhm2_err,
                  background, back_err, rmse, reduced_chi2])

    return qtab, fit_table


def do_one(spectrum_table, vel=0., xplot=False, outroot=''):
    '''
    Fit every line in data/mappings_snr_lines.txt to a single spectrum:
    singlets via lvm_gaussfit.fit_gaussian_to_spectrum, blend groups via
    fit_blend_to_spectrum (or fit_double_gaussian_to_spectrum for the
    pre-existing 'oii' group).  Mirrors lvm_gaussfit.do_one's structure.
    '''
    clean()

    zz = 1. + (vel / 3e5)

    tab = load_lines()
    groups = get_groups(tab)
    caps = line_window_caps(tab)

    records = []
    for row in groups.get('-', []):
        wave = zz * row['Wave']
        name = row['gauss_name']
        fit_table = spectrum_table
        if name in SKY_EXCLUDE:
            # These lines already need to give up a chunk of their window
            # to a known sky-residual gap; don't also let the neighbour-
            # distance cap narrow the window further, or too little data
            # (in some cases none near the line centre at all) survives
            # both cuts for the fit to be well-constrained.
            half = window_halfwidth(wave)
            xmin, xmax = SKY_EXCLUDE[name]
            fit_table = spectrum_table[(spectrum_table['WAVE'] < xmin) | (spectrum_table['WAVE'] > xmax)]
        else:
            half = window_halfwidth(wave, cap=caps[name])
        init_fwhm = fwhm_bounds(wave)[0]
        try:
            results, xspec = fit_gaussian_to_spectrum(
                fit_table, line=name, init_wavelength=wave, init_fwhm=init_fwhm,
                wavelength_min=wave - half, wavelength_max=wave + half)
            records.append(results)
            if xplot:
                save_fit(name, xspec)
        except Exception as e:
            print(f"Fitting {name}: An exception occurred: {e}")

    for group_name, rows in groups.items():
        if group_name == '-':
            continue
        r1, r2 = sorted(rows, key=lambda r: r['Wave'])
        wave1, wave2 = zz * r1['Wave'], zz * r2['Wave']
        wmin = min(wave1, wave2) - window_halfwidth(wave1, cap=caps[r1['gauss_name']])
        wmax = max(wave1, wave2) + window_halfwidth(wave2, cap=caps[r2['gauss_name']])
        try:
            if group_name == 'oii':
                results, xspec = fit_double_gaussian_to_spectrum(
                    spectrum_table, line='oii', init_wavelength1=wave1, init_wavelength2=wave2,
                    init_fwhm=fwhm_bounds(wave1)[0], wavelength_min=wmin, wavelength_max=wmax)
            else:
                shared = group_name in SHARED_FWHM_GROUPS
                tie = group_name in TIE_CENTER_GROUPS
                results, xspec = fit_blend_to_spectrum(
                    spectrum_table, group_name, r1['gauss_name'], wave1, r2['gauss_name'], wave2,
                    wmin, wmax, shared_fwhm=shared, tie_centers=tie,
                    min_significance=MIN_SIGNIFICANCE)
            records.append(results)
            if xplot:
                save_fit(group_name, xspec)
        except Exception as e:
            print(f"Fitting {group_name}: An exception occurred: {e}")

    try:
        ztab = hstack(records)
    except Exception:
        print('Nothing fit for this spectrum')
        return []

    if xplot and outroot != '':
        plot_all(title=outroot)
        import matplotlib.pyplot as plt
        plt.savefig('Gauss_dir/%s.png' % outroot)

    return ztab


def steer(argv):
    '''
    usage: lvm_snrfit.py [-h] [-lmc] [-smc] [-v vel] [-stype SOURCE]
                          [-out root] [-plot] [-lines file.txt] filename ...
    '''
    lmc = 262.
    smc = 146.
    outname = ''
    fitsfiles = []
    specfiles = []
    xplot = False
    stype = ''
    vel = 0

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(_usage_from_doc(__doc__))
            return
        elif argv[i] == '-lmc':
            vel = lmc
        elif argv[i] == '-smc':
            vel = smc
        elif argv[i] == '-stype':
            i += 1
            stype = argv[i]
        elif argv[i] == '-plot':
            xplot = True
        elif argv[i][0:4] == '-out':
            i += 1
            outname = argv[i]
        elif argv[i] == '-v':
            i += 1
            vel = eval(argv[i])
        elif argv[i] == '-lines':
            i += 1
            global DEFAULT_LINES_FILE
            DEFAULT_LINES_FILE = argv[i]
        elif argv[i] == '-min_sig':
            i += 1
            global MIN_SIGNIFICANCE
            MIN_SIGNIFICANCE = eval(argv[i])
        elif argv[i][0] == '-':
            print('Unknown options :', argv)
            return
        elif argv[i].count('.fits'):
            fitsfiles.append(argv[i])
        elif argv[i].count('.txt'):
            specfiles.append(argv[i])
        else:
            print('Unknown options :', argv)
            return
        i += 1

    for one_file in fitsfiles:
        results = do_all(one_file, vel, outname, xplot, do_one_func=do_one)
        analyze(results)

    if len(specfiles) > 0:
        do_individual(specfiles, vel, stype, outname, do_one_func=do_one)

    return


if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
