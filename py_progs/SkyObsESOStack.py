#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Median-stack SkyObsESOCompare.py's observed spectra, binned by moon
    altitude, against the SAME rows' reconstructed scaled continuum and
    line models -- a visual check of WHERE in wavelength (and under what
    moon geometry) the ESO sky model's continuum shape actually diverges
    from real sky, rather than a single blended scalar per exposure.

    Motivated by SkyObsESO_analysis.py's plot_fit_quality_correlation
    (strong FIT_QUALITY_CONT/FIT_QUALITY_LINE correlation, r~0.66-0.98,
    for a first real dataset) and the LINE_SCALE_MEDIAN experiment that
    followed it (no evidence the blended LINE_SCALE is itself distorted
    by a handful of bright regions) -- together pointing at the continuum
    fit, not the line fit, as the thing to actually improve next.  This
    script is the follow-up: look at the continuum mismatch directly.

Command line usage (if any):

    usage: SkyObsESOStack.py [-telescope FAR|NEAR|SCI|SKYE|SKYW]
                             [-engine local|remote|auto] [-delta N]
                             [-mask FILE] [-airmass] [-min-moon-alt DEG]
                             [-min-moon-fli FRAC]
                             [-poly-order N | -spline-knots N]
                             [-mask-grow ANGSTROMS] [-out OUTDIR] filename

    Arguments::

        filename    A SkyObsESOCompare.py XCframe/XSFrame-mode output
                    FITS file (FLUX/SKY_EAST/SKY_WEST + DRP_ALL, with
                    <PREFIX>_CONT_SCALE/<PREFIX>_LINE_SCALE already fit)
                    -- NOT the original XCframe/XSFrame file itself, and
                    not Sky_<name>.fits mode (no Near/Far columns there).

    Options::

        -telescope T     FAR (default) | NEAR | SCI | SKYE | SKYW --
                         which population to stack.  FAR/NEAR are
                         resolved per row from filename's own Near/Far
                         columns (which physical telescope, SKY_EAST or
                         SKY_WEST, was nearer the science field for that
                         exposure); default FAR since that's the
                         population SkyObsESO_analysis.py's
                         FIT_QUALITY_CONT/FIT_QUALITY_LINE correlation
                         was cleanest for (pure sky, no source
                         contamination, r=0.98).
        -engine E        local (default) | remote | auto -- passed to
                         SkySepESO._get_sky_model, re-fetching
                         MOON/ZODI/DIFFUSE/LINES separately per row (see
                         Notes for why this is needed at all).
        -delta N         process every N-th row; default 1 (all rows).
                         Each row costs one fresh model fetch (no
                         re-fitting), so this is the way to keep a
                         first look cheap on a large file.
        -mask FILE       clean-pixel mask FITS file, used only to set
                         plot y-scales from the continuum-region pixels
                         (default: data/sky_mask.fits)
        -airmass         run do_airmass_check (see its own docstring)
                         INSTEAD of the default moon_alt stacking -- the
                         same median-stacked two-panel-PNG-per-bin
                         format, but binned by the Moon's actual airmass
                         (AIRMASS_BINS) rather than moon_alt, a more
                         direct test of whether the assumed extinction/
                         transmission correction is wrong than the
                         moon_alt bins alone.  Not combined with the
                         default run in one invocation -- see Notes.
        -min-moon-alt DEG    -airmass only; Moon-above-horizon cutoff
                             (default 0)
        -min-moon-fli FRAC   -airmass only; minimum Moon illuminated
                             fraction (default 0.5) -- see
                             do_airmass_check's own docstring for why
                             this matters.
        -poly-order N    reconstruct the continuum with a fresh degree-N
                         polynomial multiplicative fit (see
                         _fit_continuum_poly) instead of each row's
                         already-fit CONT_SCALE -- see Description.
                         Mutually exclusive with -spline-knots.
        -spline-knots N  reconstruct the continuum with a fresh cubic-
                         B-spline multiplicative fit with N interior
                         knots instead (_fit_continuum_spline) -- the
                         matched-degrees-of-freedom, locally- rather
                         than globally-supported alternative to
                         -poly-order -- see Description.  Mutually
                         exclusive with -poly-order.
        Default for both -poly-order/-spline-knots (omitted): use the
        stored CONT_SCALE, today's original behavior.
        -mask-grow ANGSTROMS   dilate the clean-pixel mask by this many
                             Angstroms around every line region before
                             fitting (_grow_mask) -- default 0 (mask
                             used as-is).  A test of whether the mask
                             under-masks real line wings, particularly
                             in the IR -- see Description/History.
                             Independent of -poly-order/-spline-knots.
        -out OUTDIR      output directory for PNGs; default <stem>_stack

Description:

    filename's DRP_ALL already carries, per row, the FINAL fit
    coefficients <PREFIX>_CONT_SCALE and <PREFIX>_LINE_SCALE from a
    SkyObsESOCompare.py run -- but not the CONT1/LINES curves themselves
    (SkyObsESOCompare.py never writes those; see its own Description).
    So for each row selected, this script:

    1. Resolves which physical extension (FLUX/SKY_EAST/SKY_WEST) and
       ra/dec to use for the requested -telescope (fixed for SCI/SKYE/
       SKYW; per-row via the Near/Far columns for NEAR/FAR).
    2. Re-fetches MOON/ZODI/DIFFUSE/LINES for that row's ra/dec/obstime
       (SkySepESO._get_sky_model, the same call one_row itself makes --
       reused rather than duplicated) and interpolates onto filename's
       WAVE grid.
    3. Reconstructs the FINAL scaled line curve using the row's OWN
       already-fit <PREFIX>_LINE_SCALE (no refitting)::

           LINES_SCALED  = <PREFIX>_LINE_SCALE * LINES

       The continuum, CONT1_SCALED, is either the row's own already-fit
       <PREFIX>_CONT_SCALE applied to the raw model sum (no refitting,
       the original behavior, still the default)::

           CONT1_SCALED  = <PREFIX>_CONT_SCALE * (MOON + ZODI + DIFFUSE)

       or, if -poly-order N is given, a FRESH degree-N polynomial
       multiplicative fit against the same MOON+ZODI+DIFFUSE sum,
       ignoring the stored CONT_SCALE entirely (_fit_continuum_poly) --
       a prototype for testing whether a low-order wavelength-dependent
       correction (N=1: an overall scale plus a color/slope term) can
       absorb a systematic color-type mismatch that a single scale
       factor (N=0-equivalent) cannot; or, if -spline-knots N is given
       instead, a FRESH cubic-B-spline multiplicative fit with N
       interior knots against the same sum (_fit_continuum_spline) --
       a locally- rather than globally-supported alternative, for when
       the polynomial needs a high order to remove a real mismatch but
       starts overshooting in sparsely-anchored regions (e.g. the line-
       crowded IR, where few clean pixels remain to anchor the fit --
       a global polynomial has no way to be well-behaved there
       independently of how well-constrained it is elsewhere, a spline
       does).  Either way::

           RESID = FLUX - CONT1_SCALED - LINES_SCALED

    Whichever reconstruction is used, the clean-pixel mask feeding it can
    optionally be grown by -mask-grow ANGSTROMS around every line region
    first (_grow_mask) -- real testing found the polynomial and spline
    approaches behave THE SAME (both overpredict the continuum in the
    line-crowded IR) despite their different flexibility, which points
    at biased INPUT (mask under-masking real line wings) rather than a
    fitting-method problem: a locally-supported spline protects OTHER
    regions from a given region's problems, it does not protect a region
    from problems in its OWN input data -- see History.  Direct
    measurement on real data found median clean-pixel gaps of only
    ~1.5 A in the IR (Z arm), versus a median line-region width of 10 A
    there -- gaps easily narrow enough to still be catching unresolved
    line-wing flux.

    Rows are grouped into fixed moon_alt bins (MOON_ALT_BINS below) --
    fixed physical bins rather than quantile bins, since a moon_alt
    regime change (horizon crossing) is what's being isolated, not just
    a statistically convenient split.  The Moon still scatters real
    light into the sky spectrum somewhat below the geometric horizon, so
    -20 to 0 is its own bin rather than folded into a single
    below-horizon catch-all.

    Within each bin, FLUX/CONT1_SCALED/LINES_SCALED/RESID are combined
    with a per-pixel MEDIAN across rows, not a mean -- a mean would let
    whichever few exposures happen to be brightest in that moon_alt bin
    (moon phase and sun altitude both vary independently of moon
    altitude) dominate the stack, the same domination problem this
    codebase already routes around elsewhere (LINE_FRAC_*, the
    LINE_SCALE_MEDIAN experiment).  RESID is the median of each row's OWN
    residual, not (median FLUX) - (median CONT1_SCALED) -
    (median LINES_SCALED) -- the two aren't identical in general, and the
    per-row residual is the physically meaningful "how far off was the
    fit for a typical exposure" question.

    One two-panel PNG is written per non-empty bin: the top panel
    overlays the median observed spectrum against the median scaled
    continuum and median scaled continuum+lines; the bottom panel is the
    median residual alone, at the same wavelength scale.

Notes::

    Each row costs one ESO-model fetch (same per-row cost as
    SkyObsESOCompare.py itself, minus the NNLS/least-squares fitting) --
    processing every row of a large file a second time this way roughly
    doubles the total time already spent producing filename in the first
    place.  Use -delta for an initial look.

    Rows whose <PREFIX>_CONT_SCALE/<PREFIX>_LINE_SCALE/moon_alt aren't
    finite (failed rows from the original SkyObsESOCompare.py run) are
    skipped before any model fetch is attempted, not just NaN-filled
    afterward.  Rows whose model fetch fails here (independently
    possible even though the original run already fetched it once, since
    -engine can differ) are counted and skipped per bin, reported in
    that bin's PNG title.

    do_all and do_airmass_check are deliberately NOT chained in one CLI
    invocation (-airmass replaces the default run rather than adding to
    it): running do_all() then do_airmass_check() back-to-back in the
    same process reproducibly hit a FileExistsError in EsoSkyObs.setup()'s
    'data' symlink handling (pre-existing local-engine setup code, not
    part of this script) on do_airmass_check's very first row, every
    time tested -- even after manually clearing the stale symlink first.
    do_airmass_check() run on its own (a fresh process) worked reliably
    (167/167 rows, no failures) -- so the issue is specific to chaining
    many hundreds of _get_sky_model calls across the two routines within
    one process, not either routine alone.  Root cause not fully
    diagnosed (not this script's own code); if calling both from a
    notebook rather than the CLI, be aware the same chaining may recur.

Primary routines::

    do_all             bin every row by moon_alt, median-stack each bin,
                       one PNG each
    do_airmass_check   bin every qualifying row by Moon airmass, median-
                       stack each bin, one PNG each (same format as
                       do_all, via the shared _stack_bin_rows/_plot_bin)

History::

    260730 ksl  Coding begun -- a companion to SkyObsESOCompare.py/
        SkyObsESO_analysis.py, prompted by the FIT_QUALITY_CONT/
        FIT_QUALITY_LINE correlation and the LINE_SCALE_MEDIAN
        experiment both pointing at the continuum fit as the thing to
        actually improve.  Reuses SkySepESO._get_sky_model (the model
        fetch) and SkyObsESO_analysis.add_near_far_columns (per-row
        Near/Far coefficient resolution) rather than duplicating either.
    260730 ksl  Added do_airmass_check (-airmass): once real moon_alt-
        binned stacks showed a blue-end continuum mismatch that flips
        sign between low and high moon altitude, the natural next test
        is whether that tracks the Moon's actual AIRMASS (Kasten & Young
        1989 -- diverges near the horizon, unlike moon_alt) rather than
        altitude itself -- evidence for/against a wrong extinction/
        transmission correction specifically.  Refactored do_all's
        input-loading/row-filtering into a shared _load_common, reused
        by both routines; do_airmass_check itself reuses _row_arrays
        unchanged (no new model-fetch code) and
        SkyObsESO_analysis._binned_trend for its trend line.
    260730 ksl  Reworked do_airmass_check from a per-row scatter/trend
        plot into the same binned-median-stack format as do_all (one
        two-panel PNG per AIRMASS_BINS bin) -- a real user test on 650
        exposures found the scatter format "difficult to interpret" and
        expected the same per-bin-spectrum layout as the moon_alt run.
        Generalized _plot_bin to take an explicit title/filename from
        the caller instead of deriving them from MOON_ALT_BINS
        internally, and factored the fetch+median-stack loop shared by
        both routines into a new _stack_bin_rows, so do_all/
        do_airmass_check now differ only in which bins they loop over
        and how they build each bin's title/filename -- removed
        BLUE_WINDOW/_blue_fractional_residual/_plot_airmass_dependence
        (the old scatter-plot machinery) as unused once nothing called
        them.  Also added a moon_fli (illuminated fraction) cut
        (-min-moon-fli, default 0.5) alongside the existing moon_alt>0
        one: the same user feedback pointed out that a Moon above the
        horizon but faint (low illumination) doesn't actually dominate
        the continuum, so including those rows dilutes the airmass
        trend this check exists to isolate rather than testing it.
    260730 ksl  Added -poly-order (_fit_continuum_poly): real airmass-
        binned stacks showed a genuine blue excess that a follow-up
        check confirmed is smooth across both instrument arm boundaries
        (not a splice artifact), consistent with a broadband color-type
        mismatch rather than a local one -- the natural next test is
        whether a low-order polynomial multiplicative correction (a
        color/slope term at order 1, curvature at order 2) removes it.
        Deliberately prototyped here rather than in
        SkyObsESOCompare.py/SkyObsPalaceCompare.py (the batch-fitting
        scripts that write the permanent DRP_ALL CONT_SCALE column):
        this script already re-fetches the model per row and already
        has the bin/stack/plot machinery to check the result visually,
        without touching the tested batch scripts or a full per-
        exposure re-run.  Threaded clean/poly_order through
        _row_arrays/_stack_bin_rows/do_all/do_airmass_check; poly_order
        defaults to None everywhere, reproducing the original stored-
        CONT_SCALE reconstruction unchanged.  If this looks real,
        promoting it into the batch scripts as permanent columns
        (enabling population-level analysis of the fitted coefficients,
        not just stacked-median plots) is a separate, later decision.
    260730 ksl  Added -spline-knots (_fit_continuum_spline): real user
        testing found -poly-order needed order 5+ before the correction
        was obviously helping, and at that order the fit visibly
        OVERSHOOTS (model > data) specifically in the line-crowded IR,
        where clean-pixel anchor points are sparse -- classic Runge's-
        phenomenon behavior for a high-order global polynomial weakly
        constrained in one region.  A cubic B-spline's locally-supported
        basis (each basis function nonzero over only a few knot
        intervals) doesn't have this failure mode: a knot in a sparse
        region only has to explain that region, and can't be dragged
        around by how well-anchored the fit is elsewhere.  Reuses
        GetSkyCont.build_design_matrix's cubic-B-spline construction
        pattern (repeated boundary knots, extrapolate=False +
        nan_to_num) but multiplicatively (each basis column x cont0)
        rather than GetSkyCont's additive from-scratch fit, so it's a
        direct matched-DOF alternative to _fit_continuum_poly (N=2
        interior knots ~ order-5 polynomial's 6 coefficients).  Factored
        the poly/spline title-tag/filename-suffix logic (previously
        duplicated in do_all/do_airmass_check) into a shared
        _fit_mode_tag, which also enforces that poly_order and
        spline_knots aren't both given at once (ambiguous -- only one
        continuum reconstruction can be active per run).
    260730 ksl  Added -mask-grow (_grow_mask): real user testing found
        -poly-order and -spline-knots behave THE SAME (both overpredict
        the continuum in the line-crowded IR, no meaningful difference
        between them) -- since the two have genuinely different failure
        modes (global-overshoot vs none), behaving identically points at
        a shared INPUT problem, not a fitting-method one, i.e. the
        clean-pixel mask itself.  Measured this directly on real data:
        median line-region width in the Z arm (IR) is 10 A versus a
        median CLEAN GAP of only ~1.5 A there (versus 4.5 A in the B
        arm) -- gaps easily narrow enough to still be catching
        unresolved line-wing flux, biasing every "clean" anchor point in
        the IR high regardless of what's fit to them.  Also clarified
        (user question) why the spline's local support doesn't already
        fix this: locality protects OTHER regions from a given region's
        problems, it doesn't protect a region from problems in its OWN
        input data -- if anything a global polynomial has a weak
        accidental defense here (forced to also satisfy the well-
        behaved rest of the spectrum) that a spline's genuine
        independence removes.  _grow_mask dilates the clean mask by a
        specified number of Angstroms around every already-masked line
        region (per-segment wavelength range, not a fixed pixel count,
        so it's robust to non-uniform pixel spacing) before either fit;
        wired into _load_common (so both do_all/do_airmass_check and the
        plotting y-scales use the SAME grown mask, not just the fit) and
        folded into _fit_mode_tag so the title/filename tag reports it
        alongside whichever poly/spline setting is active.  Deliberately
        edits nothing outside this script -- data/sky_mask.fits itself
        (used throughout the codebase for many other purposes) is left
        untouched; this is a same reversible-prototype-only spirit as
        -poly-order/-spline-knots.

'''

import sys
import os
from pathlib import Path

# When run as a script (not imported), force a non-interactive backend
# BEFORE importing pyplot, so PNGs can be written with no display
# available -- same convention as SkyObsESO_analysis.py.
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

# ensure py_progs siblings are importable when running directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from SkySepESO import _get_sky_model
from SkyObsESO_analysis import add_near_far_columns
from GetSkyCont import load_mask, _interp_mask_to_wave
from SkyObsESOCompare import _find_segments

DEFAULT_MASK_FILE = Path(__file__).resolve().parent.parent / 'data' / 'sky_mask.fits'

_USAGE = '''Usage:
  SkyObsESOStack.py [-telescope FAR|NEAR|SCI|SKYE|SKYW]
                    [-engine local|remote|auto] [-delta N] [-mask FILE]
                    [-airmass] [-min-moon-alt DEG] [-min-moon-fli FRAC]
                    [-poly-order N | -spline-knots N] [-mask-grow ANGSTROMS]
                    [-out OUTDIR] filename

Arguments:
  filename         SkyObsESOCompare.py XCframe/XSFrame-mode output FITS
                   file

Options:
  -telescope T     FAR (default) | NEAR | SCI | SKYE | SKYW
  -engine E        local (default) | remote | auto
  -delta N         process every N-th row (default: 1 = all rows)
  -mask FILE       clean-pixel mask FITS file, used only to set plot
                   y-scales from the continuum-region pixels (default:
                   data/sky_mask.fits); also the fit mask when
                   -poly-order/-spline-knots is given
  -airmass         run do_airmass_check INSTEAD of the default moon_alt
                   stacking: same median-stacked two-panel-PNG-per-bin
                   format as the default run, but binned by the Moon's
                   actual airmass (AIRMASS_BINS) rather than moon_alt --
                   a direct test of whether the assumed extinction/
                   transmission correction is wrong, rather than
                   moon_alt alone.  Run as a separate invocation, not
                   combined with the default run -- see Notes.
  -min-moon-alt DEG    -airmass only; Moon-above-horizon cutoff
                       (default: 0)
  -min-moon-fli FRAC   -airmass only; minimum Moon illuminated fraction
                       (default: 0.5) -- restricts to conditions where
                       the Moon is actually expected to be bright,
                       since this check is only meaningful where the
                       Moon dominates the continuum budget
  -poly-order N    reconstruct the continuum with a FRESH degree-N
                   polynomial multiplicative fit (see
                   _fit_continuum_poly) instead of each row's already-
                   fit CONT_SCALE -- N=0 is a single overall scale
                   (same quantity CONT_SCALE represents, re-fit here);
                   N=1 adds a color/slope term; N=2 adds curvature.
                   Mutually exclusive with -spline-knots.
  -spline-knots N  reconstruct the continuum with a FRESH cubic-B-spline
                   multiplicative fit with N interior knots instead (see
                   _fit_continuum_spline) -- the matched-DOF, locally-
                   rather than globally-supported alternative to
                   -poly-order (N=2 knots ~ order-5 polynomial's 6
                   coefficients).  Mutually exclusive with -poly-order.
  Default for both (omitted): use the stored CONT_SCALE, i.e. today's
  original behavior.  Both work with either the default moon_alt run or
  -airmass.
  -mask-grow ANGSTROMS   dilate the clean-pixel mask by this many
                   Angstroms around every line region before fitting
                   (default: 0, mask used as-is) -- a test of whether
                   the mask under-masks real line wings, particularly in
                   the IR (median clean gap there measured at only
                   ~1.5 A -- see History).  Independent of -poly-order/
                   -spline-knots -- combine with either, or with neither
                   (still changes the stored-CONT_SCALE reconstruction's
                   y-scale/residual, though not the stored CONT_SCALE
                   value itself, which was fit before this option existed).
  -out OUTDIR      output directory for PNGs (default: <stem>_stack)
'''

# Moon-altitude bin edges in degrees: [lo, hi), plus a filename-safe
# label -- in ascending order, so a numeric index prefixed onto each
# output filename (see _plot_bin) sorts alphabetically from most
# negative to most positive.  -20 to 0 is its own bin (not lumped into
# "below horizon") since the Moon still scatters real light into the sky
# spectrum somewhat below the geometric horizon.
MOON_ALT_BINS = [
    (-90.0, -20.0, 'below_m20'),
    (-20.0,   0.0, 'm20_to_0'),
    (  0.0,  20.0, '0_to_20'),
    ( 20.0,  40.0, '20_to_40'),
    ( 40.0,  60.0, '40_to_60'),
    ( 60.0,  90.0, '60_to_90'),
]

# Airmass bin edges: [lo, hi), ascending, roughly log-spaced (airmass is
# a much more strongly nonlinear function of altitude near the horizon
# than moon_alt itself is) -- for do_airmass_check.  Same zero-padded-
# index-prefix-for-sorting convention as MOON_ALT_BINS.
AIRMASS_BINS = [
    ( 1.0,  1.2, '1p0_to_1p2'),
    ( 1.2,  1.5, '1p2_to_1p5'),
    ( 1.5,  2.0, '1p5_to_2p0'),
    ( 2.0,  3.0, '2p0_to_3p0'),
    ( 3.0,  5.0, '3p0_to_5p0'),
    ( 5.0, 40.0, '5p0_to_40p0'),
]

# Physical image extension + DRP_ALL ra/dec column pair for telescopes
# that don't vary row to row.  NEAR/FAR are resolved per row instead (see
# _row_ext) since which physical telescope is nearer varies by exposure.
_FIXED_EXT = {
    'SCI':  ('FLUX',     'sci_ra',  'sci_dec'),
    'SKYE': ('SKY_EAST', 'skye_ra', 'skye_dec'),
    'SKYW': ('SKY_WEST', 'skyw_ra', 'skyw_dec'),
}


def _row_ext(telescope, drp_all, i):
    '''
    Return (image_ext_name, ra, dec) for DRP_ALL row i under the
    requested -telescope.  Fixed lookup for SCI/SKYE/SKYW; for NEAR/FAR,
    resolved per row from drp_all's own Near/Far columns (which physical
    telescope, SKY_EAST or SKY_WEST, was nearer the science field for
    that exposure).
    '''
    if telescope in _FIXED_EXT:
        ext, ra_col, dec_col = _FIXED_EXT[telescope]
        return ext, float(drp_all[ra_col][i]), float(drp_all[dec_col][i])

    which = drp_all['Near'][i] if telescope == 'NEAR' else drp_all['Far'][i]
    if which == 'SKY_EAST':
        return 'SKY_EAST', float(drp_all['skye_ra'][i]), float(drp_all['skye_dec'][i])
    return 'SKY_WEST', float(drp_all['skyw_ra'][i]), float(drp_all['skyw_dec'][i])


def _fit_continuum_poly(flux, clean, wave, cont0, order):
    '''
    Fit flux[clean] ~ P(wave)*cont0[clean], P a degree-`order` polynomial
    in a normalized wavelength coordinate x = (wave-wave_mid)/wave_half
    (mapped to roughly [-1, 1] over wave[clean]'s own range -- raw
    Angstrom values would badly condition the higher powers).  order=0
    is a single overall scale (the same quantity the stored CONT_SCALE
    represents, freshly re-fit here rather than reusing the stored
    coefficient); order=1 adds a color/slope term (what would absorb a
    smooth blue-excess-type mismatch); order=2 adds curvature.

    Ordinary (unconstrained) linear least squares -- np.linalg.lstsq --
    since the design matrix columns (x^k * cont0) make this linear in
    the coefficients; no NNLS/iteration needed.

    Returns (cont_fit, coeffs): cont_fit = P(wave)*cont0 evaluated at
    every wavelength (not just clean); coeffs is the length-(order+1)
    polynomial coefficient array in the normalized coordinate.
    '''
    wave_clean = wave[clean]
    wave_mid = 0.5 * (wave_clean.min() + wave_clean.max())
    wave_half = 0.5 * (wave_clean.max() - wave_clean.min())
    x = (wave - wave_mid) / wave_half

    powers = np.arange(order + 1)
    A_full = (x[:, None] ** powers[None, :]) * cont0[:, None]
    A_clean = A_full[clean]

    coeffs, _resid, _rank, _sv = np.linalg.lstsq(A_clean, flux[clean], rcond=None)
    cont_fit = A_full @ coeffs
    return cont_fit, coeffs


def _fit_continuum_spline(flux, clean, wave, cont0, n_knots):
    '''
    Fit flux[clean] ~ S(wave)*cont0[clean], S a smooth multiplicative
    correction represented by a cubic B-spline with n_knots interior
    knots (evenly spaced across wave[clean]'s own range) -- the direct
    matched-degrees-of-freedom alternative to _fit_continuum_poly: a
    B-spline basis has LOCAL support (each basis function nonzero over
    only a few knot intervals) instead of a global polynomial's every-
    term-nonzero-everywhere support, so it can flex to match a real
    mismatch in a sparsely-anchored region (e.g. the IR, where dense
    airglow-line contamination leaves few clean pixels) without that
    flexing destabilizing the fit elsewhere -- unlike a polynomial, which
    is prone to exactly this (Runge's-phenomenon-style overshoot) once
    pushed to a similarly high degree of freedom.  n_knots=2 gives 6
    basis functions (n_knots + k + 1, k=3 cubic) -- the same DOF as
    _fit_continuum_poly's order=5.

    Same cubic B-spline construction as GetSkyCont.build_design_matrix
    (repeated boundary knots so the spline is fully determined at the
    domain edges, extrapolate=False + nan_to_num for wavelengths outside
    the knot range) -- reused by pattern, not imported, since here each
    basis column is multiplied by cont0 (a multiplicative correction)
    rather than used as-is in GetSkyCont's additive from-scratch fit.

    Ordinary (unconstrained) linear least squares -- np.linalg.lstsq --
    same reasoning as _fit_continuum_poly: linear in the coefficients
    regardless of basis choice.

    Returns (cont_fit, coeffs) -- same shape as _fit_continuum_poly;
    coeffs has length n_knots + 4 (k + 1, k=3).
    '''
    from scipy.interpolate import BSpline

    wave_clean = wave[clean]
    wmin, wmax = wave_clean.min(), wave_clean.max()
    k = 3
    interior = np.linspace(wmin, wmax, n_knots + 2)[1:-1] if n_knots > 0 else np.array([])
    t = np.r_[[wmin] * (k + 1), interior, [wmax] * (k + 1)]
    n_b = len(t) - k - 1

    A_full = np.zeros((len(wave), n_b), dtype=float)
    for j in range(n_b):
        c = np.zeros(n_b)
        c[j] = 1.0
        A_full[:, j] = np.nan_to_num(BSpline(t, c, k, extrapolate=False)(wave))
    A_full = A_full * cont0[:, None]
    A_clean = A_full[clean]

    coeffs, _resid, _rank, _sv = np.linalg.lstsq(A_clean, flux[clean], rcond=None)
    cont_fit = A_full @ coeffs
    return cont_fit, coeffs


def _row_arrays(x, drp_all, wave, telescope, i, engine, clean=None, poly_order=None,
                spline_knots=None, verbose=False):
    '''
    For DRP_ALL row i, reconstruct and return (flux, cont1_scaled,
    lines_scaled, resid) -- the observed spectrum actually analyzed, the
    FINAL scaled line curve rebuilt from the row's own already-fit
    <telescope>_LINE_SCALE coefficient, and a freshly refetched ESO model
    (see module Description for why the fetch has to be redone).

    The continuum, cont1_scaled, comes from one of three places (at most
    one of poly_order/spline_knots may be given -- see do_all/
    do_airmass_check, which reject both at once):

    poly_order=None, spline_knots=None (default): the row's own already-
        fit <telescope>_CONT_SCALE coefficient, i.e. today's/the
        original reconstruction --
        CONT1_SCALED = CONT_SCALE * (MOON+ZODI+DIFFUSE).
    poly_order=N (an int): a FRESH degree-N polynomial multiplicative fit
        (_fit_continuum_poly) against MOON+ZODI+DIFFUSE, ignoring the
        stored CONT_SCALE entirely -- requires clean (the clean-pixel
        mask) to be given.  See _fit_continuum_poly's own docstring.
    spline_knots=N (an int): a FRESH cubic-B-spline multiplicative fit
        with N interior knots (_fit_continuum_spline), same idea as
        poly_order but with a locally- rather than globally-supported
        basis -- also requires clean.  See _fit_continuum_spline's own
        docstring.

    Returns None if the model fetch fails for this row.
    '''
    ext, ra, dec = _row_ext(telescope, drp_all, i)
    flux = np.asarray(x[ext].data[i], dtype=float)

    obstime = drp_all['obstime'][i]
    model_tab, err = _get_sky_model(ra, dec, obstime, engine=engine, verbose=verbose)
    if model_tab is None:
        return None

    moon    = np.interp(wave, model_tab['WAVE'], model_tab['MOON'])
    zodi    = np.interp(wave, model_tab['WAVE'], model_tab['ZODI'])
    diffuse = np.interp(wave, model_tab['WAVE'], model_tab['DIFFUSE'])
    lines   = np.interp(wave, model_tab['WAVE'], model_tab['LINES'])

    if poly_order is None and spline_knots is None:
        cont_scale = float(drp_all['%s_CONT_SCALE' % telescope][i])
        cont1_scaled = cont_scale * (moon + zodi + diffuse)
    elif poly_order is not None:
        cont0 = moon + zodi + diffuse
        cont1_scaled, _coeffs = _fit_continuum_poly(flux, clean, wave, cont0, poly_order)
    else:
        cont0 = moon + zodi + diffuse
        cont1_scaled, _coeffs = _fit_continuum_spline(flux, clean, wave, cont0, spline_knots)

    line_scale = float(drp_all['%s_LINE_SCALE' % telescope][i])
    lines_scaled = line_scale * lines
    resid = flux - cont1_scaled - lines_scaled
    return flux, cont1_scaled, lines_scaled, resid


def _continuum_yscale(flux, cont, clean, margin=1.3):
    '''
    Y-limits for the top (spectrum) panel, set from the CONTINUUM-region
    pixels of flux/cont rather than the full spectrum -- airglow line
    cores are ~1-2 orders of magnitude brighter than the continuum, so a
    scale that fits them in too squashes the continuum comparison (the
    actual point of this plot) flat against zero.  Lines still plot, just
    clipped at the top of the frame.  Returns (0, ymax), or None if no
    finite clean-pixel values are available.
    '''
    ref = np.concatenate([flux[clean], cont[clean]])
    ref = ref[np.isfinite(ref)]
    if ref.size == 0:
        return None
    return 0.0, float(np.nanpercentile(ref, 99.5)) * margin


def _residual_yscale(resid, clean, margin=1.3):
    '''
    Y-limits for the residual panel, set from the CONTINUUM-region
    pixels' residual spread rather than the full array -- a handful of
    bright-line-core residuals are otherwise far larger than the
    continuum-region residual actually being diagnosed, and squash it
    flat the same way _continuum_yscale's flux/cont would be.  Returns
    (-half, half), or None if no finite clean-pixel values are available.
    '''
    ref = resid[clean]
    ref = ref[np.isfinite(ref)]
    if ref.size == 0:
        return None
    half = max(abs(np.nanpercentile(ref, 1)), abs(np.nanpercentile(ref, 99))) * margin
    return -half, half


def _plot_bin(wave, flux, cont, lines, resid, clean, title, filename, outdir):
    '''
    One two-panel PNG for a single stacked bin (moon_alt or airmass --
    generic, the caller builds both the title and the filename): top
    panel overlays the median observed spectrum against the median
    scaled continuum and median scaled continuum+lines; bottom panel is
    the median per-row residual alone, same wavelength axis.  Both
    panels' y-limits are set from the continuum-region (clean-pixel)
    values (see _continuum_yscale/_residual_yscale), not autoscaled, so
    the continuum comparison this plot exists for is actually visible.
    '''
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 7), sharex=True, gridspec_kw={'height_ratios': [2, 1]})

    ax1.plot(wave, flux, lw=0.6, color='k', label='Observed (median)')
    ax1.plot(wave, cont, lw=0.6, color='C0', label='Scaled continuum (median)')
    ax1.plot(wave, cont + lines, lw=0.5, color='C1', alpha=0.7,
             label='Scaled continuum + lines (median)')
    ax1.set_ylabel('Flux')
    ax1.set_title(title)
    ax1.legend(fontsize=8)
    top_ylim = _continuum_yscale(flux, cont, clean)
    if top_ylim is not None:
        ax1.set_ylim(*top_ylim)

    ax2.axhline(0, color='gray', lw=0.5)
    ax2.plot(wave, resid, lw=0.6, color='C3')
    ax2.set_ylabel('Median residual\n(flux - continuum - lines)')
    ax2.set_xlabel('Wavelength')
    bottom_ylim = _residual_yscale(resid, clean)
    if bottom_ylim is not None:
        ax2.set_ylim(*bottom_ylim)

    fig.tight_layout()
    outfile = outdir / filename
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print('Wrote %s' % outfile)


def _stack_bin_rows(x, drp_all, wave, telescope, engine, bin_rows, progress_label,
                    clean=None, poly_order=None, spline_knots=None):
    '''
    Shared per-bin worker for do_all/do_airmass_check: fetch+reconstruct
    every row in bin_rows (_row_arrays) and median-combine FLUX/
    CONT1_SCALED/LINES_SCALED/RESID across them -- see module Description
    for why MEDIAN, not mean, and why RESID is the median of each row's
    own residual rather than built from the other three medians.

    clean/poly_order/spline_knots are passed straight through to
    _row_arrays -- see its own docstring; poly_order=spline_knots=None
    (default) reproduces the original stored-CONT_SCALE reconstruction,
    unchanged from before these parameters existed.

    Returns (med_flux, med_cont, med_lines, med_resid, n_ok, n_fail), or
    None if every row in bin_rows failed its model fetch.
    '''
    fluxes, conts, lineses, resids = [], [], [], []
    for n, i in enumerate(bin_rows):
        out = _row_arrays(x, drp_all, wave, telescope, int(i), engine,
                          clean=clean, poly_order=poly_order, spline_knots=spline_knots)
        if out is not None:
            fluxes.append(out[0]); conts.append(out[1])
            lineses.append(out[2]); resids.append(out[3])
        if (n + 1) % 100 == 0 or (n + 1) == len(bin_rows):
            print('  %s: %d of %d' % (progress_label, n + 1, len(bin_rows)))

    n_fail = len(bin_rows) - len(fluxes)
    if not fluxes:
        return None

    med_flux  = np.median(np.array(fluxes), axis=0)
    med_cont  = np.median(np.array(conts), axis=0)
    med_lines = np.median(np.array(lineses), axis=0)
    med_resid = np.median(np.array(resids), axis=0)
    return med_flux, med_cont, med_lines, med_resid, len(fluxes), n_fail


def _grow_mask(clean, wave, grow_ang):
    '''
    Dilate the line-affected (not clean) regions of clean by grow_ang
    Angstroms on each side, using each line segment's own wavelength
    range directly (robust to non-uniform pixel spacing, rather than a
    fixed pixel-count dilation).  A test of the hypothesis that the
    clean-pixel mask doesn't extend far enough into real line wings,
    especially in the IR where OH bandheads are wide (median 10 A) and
    densely packed, leaving median clean gaps of only ~1.5 A (measured
    directly on real data, see module History) -- narrow enough that
    "clean" pixels there plausibly still catch unresolved wing flux,
    biasing every continuum fit (regardless of basis) high in exactly
    the region where overprediction was seen.

    grow_ang=0 (default) returns clean unchanged.  Returns a new boolean
    array; clean itself is not modified in place.
    '''
    if grow_ang <= 0:
        return clean
    grown = clean.copy()
    for s, e in _find_segments(~clean):
        lo, hi = wave[s] - grow_ang, wave[e - 1] + grow_ang
        grown &= ~((wave >= lo) & (wave <= hi))
    return grown


def _load_common(filename, telescope, mask_file, idelta, mask_grow=0.0):
    '''
    Shared setup for do_all/do_airmass_check: open filename, load the
    clean-pixel mask (optionally grown by mask_grow Angstroms around
    every line region -- see _grow_mask), resolve NEAR/FAR columns if
    needed, and return the idelta-strided rows with finite moon_alt/
    <telescope>_CONT_SCALE/<telescope>_LINE_SCALE (failed rows from the
    original SkyObsESOCompare.py run skipped up front, not just NaN-
    filled afterward).

    Returns (x, wave, drp_all, clean, moon_alt, rows) -- x (the open
    fits.HDUList) is the caller's responsibility to close.
    '''
    x = fits.open(filename)
    wave = np.asarray(x['WAVE'].data, dtype=float)
    drp_all = Table(x['DRP_ALL'].data)

    mask_wave, mask_arr = load_mask(mask_file)
    clean = _interp_mask_to_wave(mask_wave, mask_arr, wave)
    clean = _grow_mask(clean, wave, mask_grow)

    if telescope in ('NEAR', 'FAR'):
        if 'Near' not in drp_all.colnames or 'Far' not in drp_all.colnames:
            raise ValueError('telescope=%r needs Near/Far columns -- is filename an '
                             'XCframe/XSFrame-mode SkyObsESOCompare.py output?' % telescope)
        add_near_far_columns(drp_all)

    cont_col, line_col = '%s_CONT_SCALE' % telescope, '%s_LINE_SCALE' % telescope
    if cont_col not in drp_all.colnames or line_col not in drp_all.colnames:
        raise ValueError('%s/%s not in filename -- is telescope=%r right for this file?'
                         % (cont_col, line_col, telescope))

    moon_alt = np.asarray(drp_all['moon_alt'], dtype=float)
    good = (np.isfinite(moon_alt) & np.isfinite(np.asarray(drp_all[cont_col], dtype=float))
           & np.isfinite(np.asarray(drp_all[line_col], dtype=float)))
    rows = np.arange(0, len(drp_all), idelta)
    rows = rows[good[rows]]

    return x, wave, drp_all, clean, moon_alt, rows


def _kasten_young_airmass(alt_deg):
    '''
    Kasten & Young (1989) airmass approximation from altitude (degrees),
    valid down to the horizon -- unlike the plane-parallel sec(zenith)
    approximation, which diverges to infinity as altitude -> 0 rather
    than leveling off near ~38 at the horizon.  Only physically
    meaningful for alt_deg > 0 (a real line-of-sight path through the
    atmosphere to the source) -- callers are expected to have already
    restricted to that regime (see do_airmass_check's min_moon_alt).
    '''
    alt = np.asarray(alt_deg, dtype=float)
    return 1.0 / (np.sin(np.radians(alt)) + 0.50572 * (alt + 6.07995) ** -1.6364)


def _fit_mode_tag(poly_order, spline_knots, mask_grow=0.0):
    '''
    Build the (title_tag, filename_suffix) pair describing which
    continuum reconstruction mode -- and mask-growth setting -- is
    active, shared by do_all/do_airmass_check so the two don't duplicate
    this logic (and so a future mode/setting only needs to be taught
    here).

    Raises ValueError if both poly_order and spline_knots are given --
    only one continuum reconstruction can be active for a given run.
    mask_grow is independent of that choice (applies to the mask before
    either fit) and is folded into the same tag/suffix so callers only
    need one pair of strings.
    '''
    if poly_order is not None and spline_knots is not None:
        raise ValueError('poly_order and spline_knots are mutually exclusive -- give at most one')
    if poly_order is not None:
        tag, suffix = '   poly order %d' % poly_order, '_poly%d' % poly_order
    elif spline_knots is not None:
        tag, suffix = '   spline knots %d' % spline_knots, '_spline%d' % spline_knots
    else:
        tag, suffix = '', ''
    if mask_grow > 0:
        tag += '   mask+%gA' % mask_grow
        suffix += '_grow%g' % mask_grow
    return tag, suffix


def do_airmass_check(filename, telescope='FAR', engine='local', idelta=1,
                     mask_file=DEFAULT_MASK_FILE, outdir=None,
                     min_moon_alt=0.0, min_moon_fli=0.5, poly_order=None, spline_knots=None,
                     mask_grow=0.0):
    '''
    Bin every idelta-th qualifying row of filename by the Moon's actual
    AIRMASS (AIRMASS_BINS, Kasten & Young 1989 from moon_alt) instead of
    moon_alt itself -- airmass is the physically relevant quantity for
    an atmospheric-extinction error (diverges steeply near the horizon,
    unlike the roughly linear moon_alt), so this is a more direct test
    of "is the assumed extinction/transmission correction wrong" than
    the moon_alt-binned stacks in do_all.  Same median-stacking, same
    two-panel PNG-per-bin format as do_all (_stack_bin_rows/_plot_bin,
    reused directly, not reimplemented) -- just binned by a different
    quantity.

    Rows are restricted to moon_alt > min_moon_alt (default 0 -- Moon
    above the horizon, a real line-of-sight path, where airmass is
    physically meaningful at all) AND moon_fli > min_moon_fli (default
    0.5 -- at least half illuminated).  The moon_fli cut matters: this
    check is specifically about the MOON's own predicted continuum
    shape, so it only makes sense where the Moon is actually the
    dominant continuum source -- a faint (low-illumination) Moon, even
    if above the horizon, contributes little enough that its blue-end
    mismatch would be swamped by everything else in the residual
    (airglow continuum, real per-exposure noise), diluting rather than
    testing the airmass trend this check exists to isolate.

    poly_order/spline_knots=None (default): the continuum reconstruction
    uses each row's already-fit CONT_SCALE, exactly as before these
    parameters existed.  poly_order=N: a fresh degree-N polynomial
    multiplicative fit instead (_fit_continuum_poly) -- prototype for
    testing whether a color-term/curvature correction removes the blue-
    excess residual seen in the CONT_SCALE (order-0-equivalent) plots.
    spline_knots=N: a fresh cubic-B-spline multiplicative fit with N
    interior knots instead (_fit_continuum_spline) -- the matched-DOF,
    locally- rather than globally-supported alternative to poly_order,
    for when a high polynomial order is needed but starts overshooting
    in sparsely-anchored regions (e.g. the line-crowded IR).  At most
    one of poly_order/spline_knots may be given (see _fit_mode_tag).
    Output filenames/titles get a tag when either is active, so they
    don't collide with the default-reconstruction plots in the same
    outdir.

    Writes one PNG per non-empty airmass bin
    (stack_<telescope>_airmass_<NN>_<label>[_poly<N>|_spline<N>].png)
    to outdir (created if needed; default <stem>_stack, same convention
    as do_all).

    mask_grow : Angstroms to dilate the clean mask by around every line
        region before fitting (_grow_mask) -- default 0 (mask used as-is).
        A test of whether the clean-pixel mask under-masks real line
        wings, particularly in the IR (see module History).

    Don't call this immediately after do_all() in the same process (the
    CLI's -airmass deliberately runs this INSTEAD of do_all, not after
    it) -- see module Notes for the reproducible FileExistsError that
    causes.
    '''
    x, wave, drp_all, clean, moon_alt, rows = _load_common(filename, telescope, mask_file, idelta,
                                                            mask_grow=mask_grow)

    moon_fli = np.asarray(drp_all['moon_fli'], dtype=float)
    rows = rows[(moon_alt[rows] > min_moon_alt)
               & np.isfinite(moon_fli[rows]) & (moon_fli[rows] > min_moon_fli)]

    if len(rows) == 0:
        print('No rows with moon_alt > %g and moon_fli > %g -- nothing to check.'
             % (min_moon_alt, min_moon_fli))
        x.close()
        return

    if outdir is None:
        outdir = '%s_stack' % Path(filename).stem
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Computed only at rows already restricted to moon_alt > min_moon_alt
    # (>= 0 by default) -- the Kasten & Young formula's (alt+6.08)**-1.6364
    # term produces spurious NaN/warnings for negative altitudes it was
    # never meant to be evaluated at.
    airmass_full = np.full(len(drp_all), np.nan)
    airmass_full[rows] = _kasten_young_airmass(moon_alt[rows])

    fit_tag, fit_suffix = _fit_mode_tag(poly_order, spline_knots, mask_grow=mask_grow)

    for idx, (lo, hi, label) in enumerate(AIRMASS_BINS):
        bin_rows = rows[(airmass_full[rows] >= lo) & (airmass_full[rows] < hi)]
        if len(bin_rows) == 0:
            print('airmass %s: no rows -- skipping' % label)
            continue

        stacked = _stack_bin_rows(x, drp_all, wave, telescope, engine, bin_rows,
                                  progress_label='airmass %s' % label,
                                  clean=clean, poly_order=poly_order, spline_knots=spline_knots)
        if stacked is None:
            print('airmass %s: all %d rows failed model fetch -- skipping'
                 % (label, len(bin_rows)))
            continue
        med_flux, med_cont, med_lines, med_resid, n_ok, n_fail = stacked

        fail_note = ', %d model-fetch failures' % n_fail if n_fail else ''
        title = ('%s   airmass %g to %g   moon_fli>%g   (N=%d%s)%s'
                 % (telescope, lo, hi, min_moon_fli, n_ok, fail_note, fit_tag))
        filename_out = 'stack_%s_airmass_%02d_%s%s.png' % (telescope.lower(), idx, label, fit_suffix)
        _plot_bin(wave, med_flux, med_cont, med_lines, med_resid, clean, title, filename_out, outdir)

    x.close()


def do_all(filename, telescope='FAR', engine='local', idelta=1, mask_file=DEFAULT_MASK_FILE,
          outdir=None, poly_order=None, spline_knots=None, mask_grow=0.0, verbose=False):
    '''
    Bin every idelta-th row of filename by moon_alt (MOON_ALT_BINS),
    median-stack each non-empty bin's observed spectrum and reconstructed
    scaled continuum/line curves, and write one two-panel PNG per bin to
    outdir (created if needed; default <stem>_stack).  See module
    Description for what's actually computed and why.

    mask_file : clean-pixel mask (default data/sky_mask.fits, the same
        mask SkyObsESOCompare.py itself uses) -- used only to set each
        plot's y-scale from the continuum-region pixels (see
        _continuum_yscale/_residual_yscale), not to restrict which
        pixels are stacked or plotted.  Also the fit mask used when
        poly_order/spline_knots is not None (see below).

    poly_order/spline_knots=None (default): the continuum reconstruction
    uses each row's already-fit CONT_SCALE, exactly as before these
    parameters existed.  poly_order=N: a fresh degree-N polynomial
    multiplicative fit instead (_fit_continuum_poly).  spline_knots=N: a
    fresh cubic-B-spline multiplicative fit with N interior knots instead
    (_fit_continuum_spline) -- same prototype/mutual-exclusivity as
    do_airmass_check's own poly_order/spline_knots (see its docstring
    and _fit_mode_tag).  Output filenames/titles get a tag when either
    is active.

    mask_grow : Angstroms to dilate the clean mask by around every line
        region before fitting (_grow_mask) -- default 0 (mask used
        as-is).  See do_airmass_check's own docstring.
    '''
    x, wave, drp_all, clean, moon_alt, rows = _load_common(filename, telescope, mask_file, idelta,
                                                            mask_grow=mask_grow)

    if outdir is None:
        outdir = '%s_stack' % Path(filename).stem
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fit_tag, fit_suffix = _fit_mode_tag(poly_order, spline_knots, mask_grow=mask_grow)

    for idx, (lo, hi, label) in enumerate(MOON_ALT_BINS):
        bin_rows = rows[(moon_alt[rows] >= lo) & (moon_alt[rows] < hi)]
        if len(bin_rows) == 0:
            print('moon_alt %s: no rows -- skipping' % label)
            continue

        stacked = _stack_bin_rows(x, drp_all, wave, telescope, engine, bin_rows,
                                  progress_label='moon_alt %s' % label,
                                  clean=clean, poly_order=poly_order, spline_knots=spline_knots)
        if stacked is None:
            print('moon_alt %s: all %d rows failed model fetch -- skipping'
                 % (label, len(bin_rows)))
            continue
        med_flux, med_cont, med_lines, med_resid, n_ok, n_fail = stacked

        fail_note = ', %d model-fetch failures' % n_fail if n_fail else ''
        title = '%s   moon_alt %g to %g deg   (N=%d%s)%s' % (telescope, lo, hi, n_ok, fail_note, fit_tag)
        filename = 'stack_%s_moonalt_%02d_%s%s.png' % (telescope.lower(), idx, label, fit_suffix)
        _plot_bin(wave, med_flux, med_cont, med_lines, med_resid, clean, title, filename, outdir)

    x.close()


if __name__ == '__main__':
    argv = sys.argv[1:]
    if not argv or '-h' in argv or '--help' in argv:
        print(_USAGE)
        sys.exit(0)

    telescope = 'FAR'
    engine = 'local'
    idelta = 1
    mask_file = DEFAULT_MASK_FILE
    airmass = False
    min_moon_alt = 0.0
    min_moon_fli = 0.5
    poly_order = None
    spline_knots = None
    mask_grow = 0.0
    outdir = None
    filename = None

    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == '-telescope':
            i += 1
            telescope = argv[i]
        elif arg == '-engine':
            i += 1
            engine = argv[i]
        elif arg == '-delta':
            i += 1
            idelta = int(argv[i])
        elif arg == '-mask':
            i += 1
            mask_file = argv[i]
        elif arg == '-airmass':
            airmass = True
        elif arg == '-min-moon-alt':
            i += 1
            min_moon_alt = float(argv[i])
        elif arg == '-min-moon-fli':
            i += 1
            min_moon_fli = float(argv[i])
        elif arg == '-poly-order':
            i += 1
            poly_order = int(argv[i])
        elif arg == '-spline-knots':
            i += 1
            spline_knots = int(argv[i])
        elif arg == '-mask-grow':
            i += 1
            mask_grow = float(argv[i])
        elif arg == '-out':
            i += 1
            outdir = argv[i]
        elif arg.startswith('-'):
            print('Error: unknown option "%s"' % arg)
            print(_USAGE)
            sys.exit(1)
        else:
            if filename is not None:
                print('Error: unexpected argument "%s" (filename already set to "%s")'
                     % (arg, filename))
                sys.exit(1)
            filename = arg
        i += 1

    if filename is None:
        print('Error: filename is required')
        print(_USAGE)
        sys.exit(1)

    if poly_order is not None and spline_knots is not None:
        print('Error: -poly-order and -spline-knots are mutually exclusive -- give at most one')
        sys.exit(1)

    if airmass:
        # Deliberately NOT combined with do_all() in one process -- see
        # module Notes (a reproducible FileExistsError in EsoSkyObs.py's
        # local-engine setup when the two are chained back-to-back).
        do_airmass_check(filename, telescope=telescope, engine=engine, idelta=idelta,
                         mask_file=mask_file, outdir=outdir,
                         min_moon_alt=min_moon_alt, min_moon_fli=min_moon_fli,
                         poly_order=poly_order, spline_knots=spline_knots, mask_grow=mask_grow)
    else:
        do_all(filename, telescope=telescope, engine=engine, idelta=idelta,
              mask_file=mask_file, outdir=outdir, poly_order=poly_order,
              spline_knots=spline_knots, mask_grow=mask_grow)
