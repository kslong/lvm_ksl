#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Reusable analysis routines for SkyObsESOCompare.py output, plus a
    command-line driver that produces a standard set of diagnostic plots
    (one PNG each) and text summary/outlier tables for a single input
    file in one call.  Moved here from a personal working script
    (sky_eso_analysis.py) once the plot set stabilized enough to be worth
    sharing; see History for that script's own development.

Command line usage (if any):

    usage: SkyObsESO_analysis.py [-out OUTDIR] [-outlier-thresh VAL] filename

    Arguments::

        filename         SkyObsESOCompare.py output FITS file (DRP_ALL
                         with SCI_/SKYE_/SKYW_ columns; NEAR_/FAR_ and
                         FLUX_ERROR_* are added on load)

    Options::

        -out OUTDIR      output directory for PNGs and text files;
                         default <stem>_plots
        -outlier-thresh VAL
                         <PREFIX>_LINE_SCALE outlier threshold -- rows
                         below this are flagged (default: 0.4)

Description:

    Run as a script, produces the "standard set" of plots (see
    make_standard_plots) for SCI/NEAR/FAR -- not SCI/SKYE/SKYW, since
    nothing so far has shown SKY_EAST/SKY_WEST behaving differently as
    populations, whereas NEAR/FAR (which telescope was nearer the science
    field for that row) is the split that's actually shown structure
    (e.g. a real near-field continuum excess when the Moon is below the
    horizon).  Imported as a module (e.g. in a notebook), the individual
    load/add_*/plot_*/find_outliers/summarize functions below are
    reusable independently of the standard set -- see each one's
    docstring, and the Quickstart above the History section.

Primary routines::

    load                    read + derive NEAR_/FAR_/FLUX_ERROR_* columns
    add_near_far_columns, add_flux_error_columns, add_fit_quality_columns
    weighted_mean, summarize, find_outliers
    plot_scale_cdf, plot_scale_hist, plot_scale_cdf_sep_vs_final, plot_cross_telescope,
    plot_fit_quality_correlation, plot_vs_condition, plot_vs_brightness
    make_standard_plots     the standard-set driver this script's CLI calls

Notes::

    Plotting functions here are all matplotlib, not the Plotly/HTML
    convention eval_sky.py/GetSkyCont_eval.py use elsewhere in this
    codebase -- kept consistent with how this script was actually
    developed (interactively, in a notebook) rather than rewritten.

    make_standard_plots skips (with a printed note, not an error) any
    plot whose columns aren't present -- e.g. the SEP_/FINAL fit-quality
    and LINE_FRAC plots need a SkyObsESOCompare.py run from 260728 or
    later; an older output file will just produce fewer PNGs.

Quickstart (interactive use, in a notebook)::

    from SkyObsESO_analysis import (load, summarize, plot_scale_cdf,
        plot_cross_telescope, plot_vs_condition, plot_vs_brightness,
        find_outliers, make_standard_plots)

    tab = load('hoo.fits')          # DRP_ALL + NEAR_/FAR_ + FLUX_ERROR_* columns
    summarize(tab)                  # weighted-mean/median table, all telescopes

    plot_scale_cdf(tab, 'CONT_SCALE', telescopes=('SCI', 'NEAR', 'FAR'))
    plot_cross_telescope(tab, 'CONT_SCALE', 'NEAR', 'FAR')

    plot_vs_condition(tab['moon_alt'], tab['NEAR_CONT_SCALE'] - tab['FAR_CONT_SCALE'],
                      xlabel='moon_alt', ylabel='Near - Far Cont Scale', ylims=(-3, 3))

    # ANY quantity vs brightness, ANY telescopes overlaid on one plot --
    # works like plot_scale_cdf's telescopes= selector:
    plot_vs_brightness(tab, 'FLUX_ERROR_CONT', telescopes=('SCI', 'NEAR', 'FAR'))
    plot_vs_brightness(tab, 'SEP_FIT_QUALITY_CONT', telescopes=('NEAR', 'FAR'))

    find_outliers(tab, 'SKYE_LINE_SCALE', low=0.4)

    # Or just get the whole standard set in one call, same as the CLI:
    make_standard_plots(tab, 'my_plots')

History::

    260728 ksl  Coding begun as a personal script
        (/Users/long/Projects/lvm_sky2607/GetSkyData/sky_eso_analysis.py),
        generalizing Today260727.ipynb's ad hoc, copy-pasted cells into
        reusable functions -- load/add_near_far_columns/
        add_flux_error_columns/add_fit_quality_columns/weighted_mean/
        summarize/find_outliers/plot_scale_cdf/plot_cross_telescope/
        plot_vs_condition/plot_vs_brightness -- developed and verified in
        lockstep against real SkyObsESOCompare.py output as that script
        gained its SEP_/final fit-quality split, NMAD/RMS pair, and
        LINE_FRAC region metrics the same day (see its own History).
        Moved into py_progs as SkyObsESO_analysis.py once the plot set
        stabilized, gaining make_standard_plots and a command-line driver
        (-out/-outlier-thresh) so the whole standard set can be produced
        for a new file in one call.  Standard set defaults to SCI/NEAR/
        FAR, not SCI/SKYE/SKYW -- no evidence so far that SKY_EAST/
        SKY_WEST behave differently as populations, whereas NEAR/FAR
        (which telescope was nearer the science field) is the split
        that's actually shown structure.  Absolute paths (this repo's own
        py_progs/data dirs, hardcoded in the personal-script version)
        switched to the standard __file__-relative resolution used
        throughout this codebase.
    260729 ksl  Plotting refinements from real use: every plotting
        function gained an optional title= parameter, and a new
        plot_scale_cdf_sep_vs_final overlays a quantity's SEP (3-parameter
        separation-stage, dashed) and final (unprefixed, solid -- the
        result that actually matters) curves with both a distinguishing
        line style AND legend label per telescope (plot_scale_cdf's own
        label=tel alone gives duplicate, indistinguishable entries when
        the same telescope's SEP and final curves share one axes); a
        show_sep=False option (and make_standard_plots'/-final-only's
        matching final_only) plots just the final curves where the SEP
        diagnostic isn't needed.  All five SEP-vs-final plots (fit-
        quality x2, line-region-fraction x3, the last now including the
        two fixed-absolute-threshold flavors alongside the original
        10-sigma one) use this.  Added plot_scale_hist for a non-
        cumulative view (CONT_SCALE and LINE_SCALE, final only, all
        telescopes at once -- hist_cont_scale.png/hist_line_scale.png in
        the standard set): one ax.hist() call per telescope rather than
        one call given a list of arrays (which places bars side-by-side
        within a bin, unreadable with more than one or two series), each
        unfilled (histtype='step') so the outlines overlay cleanly, the
        same convention plot_scale_cdf already used for its cumulative
        curves.  All three histogram-based plotting functions
        (plot_scale_cdf, plot_scale_hist, plot_scale_cdf_sep_vs_final)
        gained thicker line widths and axis spines/ticks for presentation
        use (_HIST_LINEWIDTH/_AXIS_LINEWIDTH module-level constants, via
        a shared _thicken_axes helper); the scatter/trend functions
        (plot_cross_telescope, plot_vs_condition, plot_vs_brightness)
        were left as they were.

'''

import sys
import os
from pathlib import Path

# When run as a script (not imported), force a non-interactive backend
# BEFORE importing pyplot, so PNGs can be written with no display
# available.  When imported (e.g. into a notebook), __name__ isn't
# '__main__', so this is skipped and whatever backend is already active
# (interactive or not) is left alone.
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from scipy.stats import pearsonr, spearmanr

# ensure py_progs siblings are importable when running directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from GetSkyCont import load_mask, _interp_mask_to_wave

DEFAULT_MASK_FILE = Path(__file__).resolve().parent.parent / 'data' / 'sky_mask.fits'

_USAGE = '''Usage:
  SkyObsESO_analysis.py [-out OUTDIR] [-outlier-thresh VAL] [-final-only] filename

Arguments:
  filename         SkyObsESOCompare.py output FITS file

Options:
  -out OUTDIR          output directory for PNGs/text files
                       (default: <stem>_plots)
  -outlier-thresh VAL  <PREFIX>_LINE_SCALE outlier threshold
                       (default: 0.4)
  -final-only          on the five SEP-vs-final plots, plot only the
                       final (solid) curves, skipping the 3-parameter
                       separation-stage (dashed) diagnostic curves
'''


# ──────────────────────────────────────────────────────────────
# Loading and derived columns
# ──────────────────────────────────────────────────────────────

def load(filename, ext='DRP_ALL'):
    '''
    Read the DRP_ALL table from a SkyObsESOCompare.py output file and add
    the ``NEAR_``/``FAR_`` and FLUX_ERROR_* derived columns, so the returned
    table is ready for every function below.
    '''
    tab = Table(fits.open(filename)[ext].data)
    add_near_far_columns(tab)
    add_flux_error_columns(tab)
    return tab


def add_near_far_columns(tab):
    '''
    For every SKYE_<suffix> column that has a matching SKYW_<suffix>
    column, add NEAR_<suffix> and FAR_<suffix>, picking each row's value
    from whichever of ``SKYE_``/``SKYW_`` is the Near/Far telescope for that row
    (using the Near/Far columns already in DRP_ALL).  Modifies tab in
    place and returns it.
    '''
    if 'Near' not in tab.colnames or 'Far' not in tab.colnames:
        raise ValueError("tab has no Near/Far columns -- is this a "
                         "SkyObsESOCompare.py XCframe-mode output?")

    is_east_near = np.asarray(tab['Near']) == 'SKY_EAST'
    suffixes = [c[len('SKYE_'):] for c in tab.colnames if c.startswith('SKYE_')]
    for suf in suffixes:
        skye_col, skyw_col = 'SKYE_' + suf, 'SKYW_' + suf
        if skyw_col not in tab.colnames:
            continue
        skye_vals = np.asarray(tab[skye_col])
        skyw_vals = np.asarray(tab[skyw_col])
        tab['NEAR_' + suf] = np.where(is_east_near, skye_vals, skyw_vals)
        tab['FAR_' + suf]  = np.where(is_east_near, skyw_vals, skye_vals)
    return tab


# Prefixes to compute FLUX_ERROR_* for; NEAR/FAR only produce anything if
# add_near_far_columns has already been run on this table.
_FLUX_ERROR_PREFIXES = ('SCI', 'SKYE', 'SKYW', 'NEAR', 'FAR')


def add_flux_error_columns(tab):
    '''
    For every prefix in _FLUX_ERROR_PREFIXES with the needed columns, adds the
    following columns::

        <PREFIX>_FLUX_ERROR_CONT = (<PREFIX>_CONT_SCALE - 1) *
            (<PREFIX>_MOON_FLUX + <PREFIX>_ZODI_FLUX + <PREFIX>_DIFFUSE_FLUX)
        <PREFIX>_FLUX_ERROR_LINE = (<PREFIX>_LINE_SCALE - 1) * <PREFIX>_LINE_FLUX

    Modifies tab in place and returns it.
    '''
    for p in _FLUX_ERROR_PREFIXES:
        cont_scale = '%s_CONT_SCALE' % p
        moon, zodi, diffuse = '%s_MOON_FLUX' % p, '%s_ZODI_FLUX' % p, '%s_DIFFUSE_FLUX' % p
        line_scale, line_flux = '%s_LINE_SCALE' % p, '%s_LINE_FLUX' % p

        if all(c in tab.colnames for c in (cont_scale, moon, zodi, diffuse)):
            tab['%s_FLUX_ERROR_CONT' % p] = ((tab[cont_scale] - 1.0)
                                             * (tab[moon] + tab[zodi] + tab[diffuse]))
        if line_scale in tab.colnames and line_flux in tab.colnames:
            tab['%s_FLUX_ERROR_LINE' % p] = (tab[line_scale] - 1.0) * tab[line_flux]
    return tab


# Which SkyObsESOCompare.py output extension holds each telescope's
# observed spectrum.
_FIT_QUALITY_EXT = {'SCI': 'FLUX', 'SKYE': 'SKY_EAST', 'SKYW': 'SKY_WEST'}


def add_fit_quality_columns(tab, filename, mask_file=DEFAULT_MASK_FILE,
                            telescopes=('SCI', 'SKYE', 'SKYW')):
    '''
    LEGACY / retroactive helper -- for SkyObsESOCompare.py output files
    produced before the ``SEP_``/FINAL naming existed at all, which only have
    a single, unprefixed <PREFIX>_CONT_NMAD and <PREFIX>_LINE_NMAD (no
    NOISE_PROXY, no FIT_QUALITY columns, and only one line fit -- against
    a CONT3/separation-stage baseline, same as current SkyObsESOCompare.py's
    SEP line fit, not its FINAL one).  Adds::

        <PREFIX>_NOISE_PROXY          per-row noise estimate (see below)
        <PREFIX>_SEP_FIT_QUALITY_CONT = <PREFIX>_CONT_NMAD  / <PREFIX>_NOISE_PROXY
        <PREFIX>_SEP_FIT_QUALITY_LINE = <PREFIX>_LINE_NMAD  / <PREFIX>_NOISE_PROXY

    NOISE_PROXY is a per-row, model-independent noise estimate from
    consecutive-pixel differences within clean pixels (sigma =
    ``1.4826*median(|diff|)/sqrt(2)``), computed here from the observed
    FLUX/SKY_EAST/SKY_WEST arrays read directly from filename -- DRP_ALL
    alone, as returned by load(), doesn't carry them.  Each ratio is a
    reduced-chi-squared-like number, directly comparable across rows
    regardless of brightness: ~1 means the model's shape is basically a
    perfect match (leftover residual is consistent with noise); values
    much greater than 1 mean real, unmodeled shape mismatch.

    Only reproduces the SEPARATION-stage metrics.  The FINAL/unprefixed
    FIT_QUALITY_CONT and FIT_QUALITY_LINE (current SkyObsESOCompare.py's
    CONT1-baseline versions) are NOT retroactively reproducible here at
    all: they need CONT1's own residual (for FIT_QUALITY_CONT) or a
    second line fit against a CONT1 baseline (for FIT_QUALITY_LINE),
    both of which need the unscaled MOON/ZODI/DIFFUSE model templates,
    not just the observed flux an old output file has.  If you need
    those, rerun with current SkyObsESOCompare.py, which computes them
    natively at fit time.

    For any telescope where tab already has a native
    <PREFIX>_SEP_FIT_QUALITY_CONT column (i.e. a SkyObsESOCompare.py run
    from after the ``SEP_``/FINAL split existed), that telescope is skipped
    with a printed note instead of silently recomputing a redundant
    duplicate.

    Re-runs add_near_far_columns(tab) at the end, so ``NEAR_``/``FAR_`` versions
    of whatever this adds are produced too if Near/Far are present; safe
    to call even though add_near_far_columns already ran once inside
    load() (also already produces ``NEAR_``/``FAR_`` versions of any NATIVE
    ``SEP_``/FINAL fit-quality columns automatically, with no need to call
    this function at all, for current-format files).

    Raises ValueError if filename's FLUX/SKY_EAST/SKY_WEST row count
    doesn't match len(tab) (e.g. tab was filtered/subset after load()).
    '''
    x = fits.open(filename)
    wave = np.asarray(x['WAVE'].data, dtype=float)
    mask_wave, mask_arr = load_mask(mask_file)
    clean = _interp_mask_to_wave(mask_wave, mask_arr, wave)

    for tel in telescopes:
        ext = _FIT_QUALITY_EXT.get(tel)
        if ext is None or ext not in x:
            continue

        if '%s_SEP_FIT_QUALITY_CONT' % tel in tab.colnames:
            print('%s_SEP_FIT_QUALITY_CONT already present (260728+ '
                 'SkyObsESOCompare.py output) -- skipping %s.' % (tel, tel))
            continue

        cont_nmad_col, line_nmad_col = '%s_CONT_NMAD' % tel, '%s_LINE_NMAD' % tel
        if cont_nmad_col not in tab.colnames and line_nmad_col not in tab.colnames:
            continue

        flux_clean = np.asarray(x[ext].data, dtype=float)[:, clean]
        if flux_clean.shape[0] != len(tab):
            raise ValueError('%s in %s has %d rows but tab has %d -- was '
                             'tab filtered/subset after load()?'
                             % (ext, filename, flux_clean.shape[0], len(tab)))

        diffs = np.diff(flux_clean, axis=1)
        with np.errstate(invalid='ignore'):
            sigma = 1.4826 * np.nanmedian(np.abs(diffs), axis=1) / np.sqrt(2)
        tab['%s_NOISE_PROXY' % tel] = sigma

        with np.errstate(divide='ignore', invalid='ignore'):
            if cont_nmad_col in tab.colnames:
                tab['%s_SEP_FIT_QUALITY_CONT' % tel] = np.asarray(tab[cont_nmad_col], dtype=float) / sigma
            if line_nmad_col in tab.colnames:
                # Legacy files' one-and-only line fit used the CONT3
                # (separation-stage) residual as its baseline -- what
                # this script calls the SEP line fit, NOT the unprefixed/
                # FINAL one (which uses CONT1 and needs a rerun to get at
                # all).  Name it accordingly so it isn't mistaken for the
                # FINAL metric it isn't.
                tab['%s_SEP_FIT_QUALITY_LINE' % tel] = np.asarray(tab[line_nmad_col], dtype=float) / sigma

    x.close()
    add_near_far_columns(tab)
    return tab


# ──────────────────────────────────────────────────────────────
# Small stats helpers
# ──────────────────────────────────────────────────────────────

def weighted_mean(values, errors):
    '''
    Inverse-variance-weighted mean and its uncertainty, ignoring
    non-finite values/errors and errors <= 0.  Returns (mean, err, n_used).
    '''
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    good = np.isfinite(values) & np.isfinite(errors) & (errors > 0)
    if not np.any(good):
        return np.nan, np.nan, 0
    w = 1.0 / errors[good] ** 2
    mean = float(np.sum(w * values[good]) / np.sum(w))
    err = float(1.0 / np.sqrt(np.sum(w)))
    return mean, err, int(good.sum())


def summarize(tab, quantities=('CONT_SCALE', 'LINE_SCALE'),
             prefixes=('SCI', 'SKYE', 'SKYW', 'NEAR', 'FAR')):
    '''
    Print the error-weighted mean (using the matching _ERR column if
    present) and plain median of <PREFIX>_<quantity> for every
    prefix/quantity combination present in tab -- replaces re-running
    np.nanmedian by hand on six differently-named columns.  Defaults to
    all five prefixes (including SKYE/SKYW individually, not just
    SCI/NEAR/FAR) since a text table isn't visually cluttered the way
    overlaying many plot lines would be, and it's a cheap way to confirm
    SKY_EAST/SKY_WEST really do stay similar to each other.
    '''
    print('%-6s %-16s %10s %10s %10s %8s' %
         ('Tel', 'Quantity', 'WMean', 'WErr', 'Median', 'N'))
    for q in quantities:
        for p in prefixes:
            col = '%s_%s' % (p, q)
            if col not in tab.colnames:
                continue
            errcol = col + '_ERR'
            vals = np.asarray(tab[col], dtype=float)
            good = np.isfinite(vals)
            median = float(np.median(vals[good])) if np.any(good) else np.nan
            if errcol in tab.colnames:
                wmean, werr, n = weighted_mean(tab[col], tab[errcol])
            else:
                wmean, werr, n = np.nan, np.nan, int(good.sum())
            print('%-6s %-16s %10.4f %10.4f %10.4f %8d' %
                 (p, q, wmean, werr, median, n))


def find_outliers(tab, column, low=None, high=None,
                  id_cols=('tileid', 'mjd', 'expnum', 'object', 'Near', 'Far')):
    '''
    Return the subset of tab where tab[column] < low or > high (either
    bound may be omitted), restricted to id_cols + column, sorted by
    column -- for tracing a scatter-plot outlier back to specific
    exposures instead of leaving it as a dangling subset.
    '''
    vals = np.asarray(tab[column], dtype=float)
    sel = np.zeros(len(tab), dtype=bool)
    if low is not None:
        sel |= vals < low
    if high is not None:
        sel |= vals > high
    cols = [c for c in id_cols if c in tab.colnames] + [column]
    out = tab[cols][sel].copy()
    out.sort(column)
    return out


# ──────────────────────────────────────────────────────────────
# Plots
# ──────────────────────────────────────────────────────────────

# Presentation-friendly styling for the three histogram-based plotting
# functions (plot_scale_cdf, plot_scale_hist, plot_scale_cdf_sep_vs_final):
# thicker line widths and axis spines/ticks so the plots stay legible
# projected or shrunk into a slide, not just on a laptop screen at full
# size.  Adjust these two constants to retune all of them at once.
_HIST_LINEWIDTH = 2.5
_AXIS_LINEWIDTH = 1.5


def _thicken_axes(ax, linewidth=_AXIS_LINEWIDTH):
    '''
    Thicken ax's spines (the plot border) and tick marks to linewidth --
    called by the histogram plotting functions so the axes stay visually
    proportionate to the thicker (_HIST_LINEWIDTH) histogram outlines.
    '''
    for spine in ax.spines.values():
        spine.set_linewidth(linewidth)
    ax.tick_params(width=linewidth)


def plot_scale_cdf(tab, quantity, telescopes=('SCI', 'SKYE', 'SKYW'),
                   xlims=(0, 3), nbins=100, title=None, ax=None):
    '''
    Cumulative distribution of <TEL>_<quantity> for each telescope,
    overlaid.  title, if given, is set via ax.set_title -- otherwise the
    axis labels are left to speak for themselves.
    '''
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    for tel in telescopes:
        col = '%s_%s' % (tel, quantity)
        if col not in tab.colnames:
            continue
        vals = np.asarray(tab[col], dtype=float)
        vals = vals[np.isfinite(vals)]
        ax.hist(vals, nbins, range=xlims, cumulative=True, histtype='step',
               density=True, linewidth=_HIST_LINEWIDTH, label=tel)
    ax.set_xlim(*xlims)
    ax.set_xlabel(quantity)
    ax.set_ylabel('Cumulative fraction')
    if title:
        ax.set_title(title)
    ax.legend()
    _thicken_axes(ax)
    ax.figure.tight_layout()
    return ax


def plot_scale_hist(tab, quantity, telescopes=('SCI', 'SKYE', 'SKYW'),
                    xlims=(0, 3), nbins=50, title=None, ax=None):
    '''
    Non-cumulative distribution of <TEL>_<quantity> for each telescope,
    overlaid -- one separate ax.hist() call per telescope (NOT one call
    given a list of arrays, which would place each telescope's bars
    side-by-side within a bin, "interspersed" and unreadable with more
    than one or two series), each drawn unfilled (histtype='step') so the
    outlines overlay cleanly instead of obscuring each other -- the same
    style plot_scale_cdf already uses for its cumulative curves, just
    without cumulative=True.
    '''
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    for tel in telescopes:
        col = '%s_%s' % (tel, quantity)
        if col not in tab.colnames:
            continue
        vals = np.asarray(tab[col], dtype=float)
        vals = vals[np.isfinite(vals)]
        ax.hist(vals, nbins, range=xlims, histtype='step', density=True,
               linewidth=_HIST_LINEWIDTH, label=tel)
    ax.set_xlim(*xlims)
    ax.set_xlabel(quantity)
    ax.set_ylabel('Probability density')
    if title:
        ax.set_title(title)
    ax.legend()
    _thicken_axes(ax)
    ax.figure.tight_layout()
    return ax


def plot_scale_cdf_sep_vs_final(tab, quantity, telescopes=('SCI', 'SKYE', 'SKYW'),
                                xlims=(0, 3), nbins=100, show_sep=True,
                                title=None, ax=None):
    '''
    CDF of the FINAL (unprefixed) <TEL>_<quantity> -- solid lines, since
    it's the result that actually matters -- with the SEP (3-parameter
    continuum-separation-stage) version overlaid too (dashed lines) unless
    show_sep=False.  E.g. quantity='FIT_QUALITY_CONT' plots FIT_QUALITY_CONT
    (solid) and, by default, SEP_FIT_QUALITY_CONT (dashed) too, one pair of
    curves per telescope.

    show_sep : bool, default True.  SEP is a diagnostic for whether the
        continuum-separation stage is helping, not itself the answer --
        pass show_sep=False for a final-only view (e.g. once SEP has
        already been checked and isn't needed on every plot going
        forward).  With show_sep=False this behaves like plot_scale_cdf
        on the final quantity alone (plain "<TEL>" legend labels, no
        line-style distinction needed since there's only one curve per
        telescope).

    plot_scale_cdf's own label=tel alone would give duplicate,
    indistinguishable legend entries when show_sep=True, since the same
    telescope then appears in both an SEP curve and a final curve on the
    same axes -- so each telescope gets a fixed color (consistent between
    its two curves), and each curve is distinguished BOTH by legend label
    ("<TEL> (3-param)" / "<TEL> (final)") and by line style, not just one
    or the other.
    '''
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    variants = [('', '-', 'final')]
    if show_sep:
        variants.insert(0, ('SEP_', '--', '3-param'))
    for i, tel in enumerate(telescopes):
        color = 'C%d' % (i % 10)
        for prefix, ls, tag in variants:
            col = '%s_%s%s' % (tel, prefix, quantity)
            if col not in tab.colnames:
                continue
            vals = np.asarray(tab[col], dtype=float)
            vals = vals[np.isfinite(vals)]
            label = '%s (%s)' % (tel, tag) if show_sep else tel
            ax.hist(vals, nbins, range=xlims, cumulative=True, histtype='step',
                   density=True, linestyle=ls, color=color,
                   linewidth=_HIST_LINEWIDTH, label=label)
    ax.set_xlim(*xlims)
    ax.set_xlabel(quantity)
    ax.set_ylabel('Cumulative fraction')
    if title:
        ax.set_title(title)
    ax.legend(fontsize=8)
    _thicken_axes(ax)
    ax.figure.tight_layout()
    return ax


def plot_cross_telescope(tab, quantity, tel_x, tel_y, lims=(0, 2), alpha=0.1,
                         title=None, ax=None):
    '''
    Scatter plot of <tel_x>_<quantity> vs <tel_y>_<quantity>, with the
    Pearson and Spearman correlation coefficients (computed over all
    finite points, not just those inside lims) annotated in the title.
    title, if given, is prepended above the correlation-stats line rather
    than replacing it.
    '''
    xcol, ycol = '%s_%s' % (tel_x, quantity), '%s_%s' % (tel_y, quantity)
    x = np.asarray(tab[xcol], dtype=float)
    y = np.asarray(tab[ycol], dtype=float)
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    ax.plot(x, y, '.', alpha=alpha)
    ax.set_xlim(*lims)
    ax.set_ylim(*lims)
    ax.set_xlabel('%s %s' % (tel_x, quantity))
    ax.set_ylabel('%s %s' % (tel_y, quantity))

    stats_line = ''
    if len(x) > 2:
        pear = pearsonr(x, y)
        spear = spearmanr(x, y)
        stats_line = ('Pearson r=%.3f (p=%.1e)   Spearman rho=%.3f (p=%.1e)  N=%d'
                      % (pear[0], pear[1], spear[0], spear[1], len(x)))
    full_title = '%s\n%s' % (title, stats_line) if title else stats_line
    if full_title:
        ax.set_title(full_title, fontsize=9)
    ax.figure.tight_layout()
    return ax


def plot_fit_quality_correlation(tab, telescopes=('SCI', 'NEAR', 'FAR'), lims=(0, 15),
                                 alpha=0.1, title=None, ax=None):
    '''
    Scatter of <TEL>_FIT_QUALITY_CONT (x) vs <TEL>_FIT_QUALITY_LINE (y),
    one color per telescope, Pearson r annotated per telescope in the
    legend (e.g. "SCI (r=0.62)") rather than pooled across telescopes --
    SCI carries real source flux on top of sky while NEAR/FAR don't (see
    SkyObsESOCompare.py Notes), so pooling could conflate genuinely
    different populations' continuum-vs-line relationships into one
    misleading number.

    Motivation: the FINAL line residual is measured against flux - CONT1,
    with CONT1 fit only on clean pixels and then extrapolated (never
    refit) into the line-affected pixels, and the line fit itself is a
    single non-negative amplitude on the LINES template with no freedom
    to absorb a smooth continuum offset/slope (see SkyObsESOCompare.py's
    one_row).  So a continuum SHAPE error can leak straight through into
    FIT_QUALITY_LINE, indistinguishable there from a genuine line-shape
    mismatch.  A strong positive correlation here is evidence that's
    actually happening at scale -- i.e. that solving the continuum
    problem would substantially fix the apparent line problem too, rather
    than the two being independent issues that both need separate work.
    '''
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    for i, tel in enumerate(telescopes):
        xcol, ycol = '%s_FIT_QUALITY_CONT' % tel, '%s_FIT_QUALITY_LINE' % tel
        if xcol not in tab.colnames or ycol not in tab.colnames:
            continue
        x = np.asarray(tab[xcol], dtype=float)
        y = np.asarray(tab[ycol], dtype=float)
        good = np.isfinite(x) & np.isfinite(y)
        x, y = x[good], y[good]
        color = 'C%d' % (i % 10)
        label = tel
        if len(x) > 2:
            r = pearsonr(x, y)[0]
            label = '%s (r=%.2f)' % (tel, r)
        ax.plot(x, y, '.', alpha=alpha, color=color, label=label)

    ax.set_xlim(*lims)
    ax.set_ylim(*lims)
    ax.set_xlabel('FIT_QUALITY_CONT')
    ax.set_ylabel('FIT_QUALITY_LINE')
    if title:
        ax.set_title(title, fontsize=9)
    ax.legend(fontsize=8)
    _thicken_axes(ax)
    ax.figure.tight_layout()
    return ax


def _binned_trend(x, y, yerr=None, nbins=20):
    '''
    Core binning logic shared by plot_vs_condition and plot_vs_brightness:
    split x into nbins equal-COUNT (quantile) bins -- so a skewed
    x-distribution like a flux (mostly near zero with a long tail) still
    gets well-populated bins rather than one giant first bin -- and
    summarize each bin's y as the inverse-variance-weighted mean (using
    yerr as the per-point uncertainty) if yerr is given, else a robust
    median with an NMAD-based standard error.  Returns (xc, yc, yce)
    lists, one entry per non-empty (>=3 points) bin.
    '''
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    good = np.isfinite(x) & np.isfinite(y)
    if yerr is not None:
        yerr = np.asarray(yerr, dtype=float)
        good &= np.isfinite(yerr) & (yerr > 0)
    x, y = x[good], y[good]
    ye = yerr[good] if yerr is not None else None

    edges = np.unique(np.percentile(x, np.linspace(0, 100, nbins + 1)))
    bin_idx = np.digitize(x, edges[1:-1])

    xc, yc, yce = [], [], []
    for i in range(len(edges) - 1):
        sel = bin_idx == i
        if sel.sum() < 3:
            continue
        xc.append(float(np.median(x[sel])))
        if ye is not None:
            m, e, _ = weighted_mean(y[sel], ye[sel])
        else:
            m = float(np.median(y[sel]))
            e = float(1.4826 * np.median(np.abs(y[sel] - m)) / np.sqrt(sel.sum()))
        yc.append(m)
        yce.append(e)
    return xc, yc, yce


def plot_vs_condition(x, y, yerr=None, xlabel='', ylabel='', xlims=None, ylims=None,
                      nbins=20, alpha=0.1, logx=False, title=None, ax=None):
    '''
    Scatter plot of y vs x, with a robust binned trend overlaid (see
    _binned_trend).

    x, y, yerr : array-like, same length (e.g. tab['moon_alt'],
        tab['NEAR_CONT_SCALE'] - tab['FAR_CONT_SCALE']).

    For overlaying several telescopes' trends on one plot (rather than
    one x/y series with a raw scatter cloud), use plot_vs_brightness.
    '''
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    good = np.isfinite(x) & np.isfinite(y)
    if yerr is not None:
        good &= np.isfinite(np.asarray(yerr, dtype=float))
    x_plot, y_plot = x[good], y[good]

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 6))
    ax.plot(x_plot, y_plot, '.', alpha=alpha, color='0.6', zorder=1)

    xc, yc, yce = _binned_trend(x, y, yerr=yerr, nbins=nbins)
    ax.errorbar(xc, yc, yerr=yce, fmt='o-', color='C3', zorder=2,
               label='binned weighted mean' if yerr is not None else 'binned median')
    if logx:
        ax.set_xscale('log')
    if xlims:
        ax.set_xlim(*xlims)
    if ylims:
        ax.set_ylim(*ylims)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    ax.legend()
    ax.figure.tight_layout()
    return ax


def _brightness(tab, tel, component):
    '''Which MEAN-brightness columns a given "component" plots against,
    per telescope prefix -- see plot_vs_brightness.'''
    if component.upper() == 'CONT':
        b = (np.asarray(tab['%s_MOON_FLUX' % tel], dtype=float)
            + np.asarray(tab['%s_ZODI_FLUX' % tel], dtype=float)
            + np.asarray(tab['%s_DIFFUSE_FLUX' % tel], dtype=float))
    else:
        b = np.asarray(tab['%s_LINE_FLUX' % tel], dtype=float)
    return np.clip(b, 1e-20, None)


def plot_vs_brightness(tab, quantity, telescopes=('SCI', 'SKYE', 'SKYW'), component='CONT',
                       nbins=20, ylims=None, logx=True, title=None, ax=None):
    '''
    <TEL>_<quantity> vs that telescope's own predicted brightness (log
    x-axis by default), one binned trend line per telescope in
    telescopes, overlaid on one plot -- works like plot_scale_cdf's
    telescopes= selector (so SCI/SKYE/SKYW or NEAR/FAR can be compared
    directly, e.g. plot_vs_brightness(tab, 'SEP_FIT_QUALITY_CONT',
    telescopes=('NEAR','FAR'))).

    component : 'CONT' (brightness = MOON_FLUX+ZODI_FLUX+DIFFUSE_FLUX,
        the continuum-fit footprint) or 'LINE' (brightness = LINE_FLUX,
        the line-fit footprint) -- which brightness proxy to use as the
        x-axis; independent of what quantity itself measures, so e.g. a
        line-fit quantity can still be plotted against continuum
        brightness if that's the comparison wanted.

    No raw scatter cloud is drawn here (unlike plot_vs_condition) -- with
    several telescopes overlaid, multiple scatter clouds would obscure
    the trend lines that are the point of the comparison.
    '''
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 6))

    for tel in telescopes:
        ycol = '%s_%s' % (tel, quantity)
        if ycol not in tab.colnames:
            continue
        brightness = _brightness(tab, tel, component)
        xc, yc, yce = _binned_trend(brightness, tab[ycol], nbins=nbins)
        ax.errorbar(xc, yc, yerr=yce, fmt='o-', label=tel)

    if logx:
        ax.set_xscale('log')
    if ylims:
        ax.set_ylim(*ylims)
    ax.set_xlabel('Predicted %s Brightness' % component.capitalize())
    ax.set_ylabel(quantity)
    if title:
        ax.set_title(title)
    ax.legend()
    ax.figure.tight_layout()
    return ax


# ──────────────────────────────────────────────────────────────
# Standard plot set
# ──────────────────────────────────────────────────────────────

_STANDARD_TELESCOPES = ('SCI', 'NEAR', 'FAR')


def make_standard_plots(tab, outdir, telescopes=_STANDARD_TELESCOPES, outlier_thresh=0.4,
                        final_only=False):
    '''
    Produce the standard set of diagnostic plots (one PNG each) plus a
    summary text file and a LINE_SCALE outlier table, in outdir (created
    if it doesn't exist).  tab should already be load()'s output (i.e.
    has ``NEAR_``/``FAR_``/FLUX_ERROR_* columns).

    final_only : bool, default False.  The five SEP-vs-final plots
        (fit_quality_cont/line, line_frac_10sig/1e15/1e14) normally show
        both the 3-parameter separation-stage curve (dashed) and the
        final curve (solid, the result that actually matters); pass
        final_only=True to plot just the final curves on all five,
        skipping the SEP diagnostic (see plot_scale_cdf_sep_vs_final's
        show_sep).

    Plots (all using telescopes, default SCI/NEAR/FAR -- see module
    Description for why NEAR/FAR rather than SKYE/SKYW)::

        cdf_cont_scale.png, cdf_line_scale.png, hist_cont_scale.png,
        hist_line_scale.png
        cross_near_far_cont_scale.png, cross_near_far_line_scale.png
            (skipped if 'NEAR' and 'FAR' aren't both in telescopes)
        near_far_vs_moonalt_cont_scale.png, near_far_vs_moonalt_line_scale.png
            (skipped under the same condition, or if moon_alt is absent)
        flux_error_cont_vs_brightness.png, flux_error_line_vs_brightness.png
        fit_quality_cont_sep_vs_final.png, fit_quality_line_sep_vs_final.png,
        fit_quality_cont_vs_line.png
            (skipped if SEP_FIT_QUALITY_CONT isn't present for every
            telescope -- needs a 260728+ SkyObsESOCompare.py run)
        line_frac_10sig_sep_vs_final.png, line_frac_1e15_sep_vs_final.png,
        line_frac_1e14_sep_vs_final.png
            (skipped under the same condition)

    Plus summary.txt (summarize() output, always all five SCI/SKYE/SKYW/
    NEAR/FAR prefixes regardless of telescopes, per summarize's own
    default) and outliers_line_scale.txt (find_outliers on
    <PREFIX>_LINE_SCALE < outlier_thresh, for each telescope in
    telescopes).
    '''
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    def _save(fig, name):
        path = outdir / name
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print('Wrote %s' % path)

    # 1. Scale distributions
    ax = plot_scale_cdf(tab, 'CONT_SCALE', telescopes=telescopes,
                        title='Continuum Scale Factor (CONT_SCALE)')
    _save(ax.figure, 'cdf_cont_scale.png')
    ax = plot_scale_cdf(tab, 'LINE_SCALE', telescopes=telescopes,
                        title='Line Scale Factor (LINE_SCALE)')
    _save(ax.figure, 'cdf_line_scale.png')

    ax = plot_scale_hist(tab, 'CONT_SCALE', telescopes=telescopes,
                         title='Continuum Scale Factor Distribution')
    _save(ax.figure, 'hist_cont_scale.png')
    ax = plot_scale_hist(tab, 'LINE_SCALE', telescopes=telescopes,
                         title='Line Scale Factor Distribution')
    _save(ax.figure, 'hist_line_scale.png')

    have_near_far = 'NEAR' in telescopes and 'FAR' in telescopes

    # 2. NEAR vs FAR consistency
    if have_near_far:
        ax = plot_cross_telescope(tab, 'CONT_SCALE', 'NEAR', 'FAR',
                                  title='NEAR vs FAR Continuum Scale Consistency')
        _save(ax.figure, 'cross_near_far_cont_scale.png')
        ax = plot_cross_telescope(tab, 'LINE_SCALE', 'NEAR', 'FAR',
                                  title='NEAR vs FAR Line Scale Consistency')
        _save(ax.figure, 'cross_near_far_line_scale.png')
    else:
        print('Skipping NEAR/FAR cross-telescope plots: telescopes=%r '
             'does not include both NEAR and FAR.' % (telescopes,))

    # 3. Near/Far geometry vs moon_alt
    if have_near_far and 'moon_alt' in tab.colnames:
        ax = plot_vs_condition(tab['moon_alt'], tab['NEAR_CONT_SCALE'] - tab['FAR_CONT_SCALE'],
                               xlabel='moon_alt', ylabel='NEAR - FAR CONT_SCALE', ylims=(-3, 3),
                               title='Near-Far Continuum Scale Difference vs Moon Altitude')
        _save(ax.figure, 'near_far_vs_moonalt_cont_scale.png')
        ax = plot_vs_condition(tab['moon_alt'], tab['NEAR_LINE_SCALE'] - tab['FAR_LINE_SCALE'],
                               xlabel='moon_alt', ylabel='NEAR - FAR LINE_SCALE', ylims=(-3, 3),
                               title='Near-Far Line Scale Difference vs Moon Altitude')
        _save(ax.figure, 'near_far_vs_moonalt_line_scale.png')
    elif have_near_far:
        print('Skipping NEAR/FAR-vs-moon_alt plots: moon_alt not in tab.')

    # 4. Flux error vs brightness
    ax = plot_vs_brightness(tab, 'FLUX_ERROR_CONT', telescopes=telescopes, component='CONT',
                            title='Continuum Flux Error vs Predicted Continuum Brightness')
    _save(ax.figure, 'flux_error_cont_vs_brightness.png')
    ax = plot_vs_brightness(tab, 'FLUX_ERROR_LINE', telescopes=telescopes, component='LINE',
                            title='Line Flux Error vs Predicted Line Brightness')
    _save(ax.figure, 'flux_error_line_vs_brightness.png')

    # 5. Fit quality, SEP (3-parameter separation stage) vs final
    have_fit_quality = all('%s_SEP_FIT_QUALITY_CONT' % t in tab.colnames for t in telescopes)
    if have_fit_quality:
        ax = plot_scale_cdf_sep_vs_final(
            tab, 'FIT_QUALITY_CONT', telescopes=telescopes, xlims=(0, 15),
            show_sep=not final_only,
            title='Continuum Fit Quality: Excess/Intrinsic Variation')
        _save(ax.figure, 'fit_quality_cont_sep_vs_final.png')

        ax = plot_scale_cdf_sep_vs_final(
            tab, 'FIT_QUALITY_LINE', telescopes=telescopes, xlims=(0, 20),
            show_sep=not final_only,
            title='Line Fit Quality: Excess/Intrinsic Variation')
        _save(ax.figure, 'fit_quality_line_sep_vs_final.png')

        ax = plot_fit_quality_correlation(
            tab, telescopes=telescopes,
            title='FIT_QUALITY_CONT vs FIT_QUALITY_LINE: shared continuum origin?')
        _save(ax.figure, 'fit_quality_cont_vs_line.png')
    else:
        print('Skipping fit-quality plots: SEP_FIT_QUALITY_CONT not present '
             'for every telescope (older SkyObsESOCompare.py output?).')

    # 6. Line region-fraction, SEP vs final -- one plot per threshold
    # flavor (10-sigma, relative to NOISE_PROXY; 1e-15 and 1e-14, fixed
    # absolute flux levels that don't rescale with an exposure's own
    # noise -- see SkyObsESOCompare.py's Description for why both kinds
    # exist).
    have_line_frac = all('%s_SEP_LINE_FRAC_10SIG' % t in tab.colnames for t in telescopes)
    if have_line_frac:
        ax = plot_scale_cdf_sep_vs_final(
            tab, 'LINE_FRAC_10SIG', telescopes=telescopes, xlims=(0, 1),
            show_sep=not final_only,
            title='Fraction of Line Regions Exceeding 10-Sigma')
        _save(ax.figure, 'line_frac_10sig_sep_vs_final.png')

        ax = plot_scale_cdf_sep_vs_final(
            tab, 'LINE_FRAC_1E15', telescopes=telescopes, xlims=(0, 1),
            show_sep=not final_only,
            title='Fraction of Line Regions Exceeding 1e-15 erg')
        _save(ax.figure, 'line_frac_1e15_sep_vs_final.png')

        ax = plot_scale_cdf_sep_vs_final(
            tab, 'LINE_FRAC_1E14', telescopes=telescopes, xlims=(0, 1),
            show_sep=not final_only,
            title='Fraction of Line Regions Exceeding 1e-14 erg')
        _save(ax.figure, 'line_frac_1e14_sep_vs_final.png')
    else:
        print('Skipping line-fraction plots: SEP_LINE_FRAC_10SIG not present '
             'for every telescope (older SkyObsESOCompare.py output?).')

    # Summary text file
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        summarize(tab)
    summary_path = outdir / 'summary.txt'
    summary_path.write_text(buf.getvalue())
    print('Wrote %s' % summary_path)

    # Outlier table
    lines = []
    for tel in telescopes:
        col = '%s_LINE_SCALE' % tel
        if col not in tab.colnames:
            continue
        out = find_outliers(tab, col, low=outlier_thresh)
        lines.append('=== %s (%s < %.3g): %d rows ===' % (tel, col, outlier_thresh, len(out)))
        lines.extend(out.pformat(max_lines=-1, max_width=-1))
        lines.append('')
    outlier_path = outdir / 'outliers_line_scale.txt'
    outlier_path.write_text('\n'.join(lines))
    print('Wrote %s' % outlier_path)


# ──────────────────────────────────────────────────────────────
# Command-line entry point
# ──────────────────────────────────────────────────────────────

if __name__ == '__main__':
    argv = sys.argv[1:]
    if not argv or '-h' in argv or '--help' in argv:
        print(_USAGE)
        sys.exit(0)

    outdir = ''
    outlier_thresh = 0.4
    final_only = False
    filename = None

    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == '-out':
            i += 1
            outdir = argv[i]
        elif arg == '-outlier-thresh':
            i += 1
            outlier_thresh = float(argv[i])
        elif arg == '-final-only':
            final_only = True
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

    if not outdir:
        stem = os.path.splitext(os.path.basename(filename))[0]
        outdir = '%s_plots' % stem

    tab = load(filename)
    make_standard_plots(tab, outdir, outlier_thresh=outlier_thresh, final_only=final_only)
