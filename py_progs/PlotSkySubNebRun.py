#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Run-level companion to PlotSkySubNebEval.py: aggregates
    SkySubNebEval.py's per-row nebular-line fits and per-tileid
    repeat_scatter groups across a whole run (one or more SkySub*.py
    output files, each with many rows -- not just one exposure) into a
    single, real HTML report (actual <h1>/<h2>/<h3> headings around
    several small Plotly figures, not one giant Plotly canvas), so a
    method comparison doesn't require paging through one file per
    exposure.

Command line usage (if any)::

        usage: PlotSkySubNebRun.py [-h] [-lines_file PATH]
                                   [-v VEL] [-lmc] [-smc] [-sigma S]
                                   [-snr_min S] [-mjd_close DAYS]
                                   [-mask PATH]
                                   [-title TITLE] [-outfile PATH]
                                   [fits_file ...]

        where

        fits_file      one or more SkySub*.py output FITS files (WAVE,
                      FLUX, SKY, DRP_ALL), the same contract
                      SkySubNebEval.py's own CLI takes. Ignored if
                      -lines_file is given.

        -lines_file PATH
                      an existing SkySubNebEval.py run's <root>_lines.fits
                      output -- skips the (currently unparallelized)
                      per-row Gaussian refit, useful for iterating on
                      this plot against a large sample without redoing
                      it every time. Exactly one of fits_file/-lines_file
                      is required. repeat_scatter's own grouping/MAD
                      summary is always recomputed from whichever row
                      table is in hand either way -- that step is cheap
                      (aggregation only, no curve fitting), so there is
                      no separate "-repeats_file" fast path.

        -v VEL / -lmc / -smc
                      nebular systemic velocity override, same
                      precedence/default (per-row DRP_ALL['Redshift']
                      lookup) as SkySubNebEval.py. Ignored if
                      -lines_file is given (the velocity used is
                      whatever that file was fit with).

        -sigma S      initial Gaussian sigma guess in Angstrom (default
                      1.0). Ignored if -lines_file is given.

        -snr_min S    minimum per-line SNR required to report a
                      doublet ratio (default 5.0) -- same gate
                      SkySubNebEval.py uses; applied here when
                      building the row table (fits_file mode) or
                      already baked into a loaded -lines_file's own
                      RATIO columns (recomputing it post-hoc from a
                      loaded file is not supported).

        -mjd_close DAYS
                      MJD span below which a repeat group is "closely
                      spaced" (default 7.0), passed to repeat_scatter.

        -mask PATH    clean-pixel mask (WAVE/MASK extensions, a
                      palace_make_mask.py-style file) for the
                      Continuum Residual (B/R/Z) section -- default
                      data/sky_mask.fits, the project-wide standard.
                      Only used in fits_file mode (the continuum
                      section needs the raw WAVE/FLUX/SKY arrays, so
                      it is skipped entirely in -lines_file mode).

        -title TITLE  report title (default: 'SkySubNebRun').

        -outfile PATH output HTML path (default: nebrun.html). Rerun
                      with the same -outfile to update it in place
                      rather than accumulating one file per attempt.

Description::

        Two input modes:

        1. Fresh: one or more raw SkySub*.py output files. Calls
           SkySubNebEval.fit_file() on each (per-row Gaussian fits to
           every NEBULAR_LINES entry, DOUBLETS ratios), stacks the
           results into one row table.
        2. Precomputed: -lines_file points at an existing
           SkySubNebEval.py run's own <root>_lines.fits, read directly
           -- no refitting.

        Either way, SkySubNebEval.repeat_scatter(..., return_rows=True)
        is then run on the row table to get both the per-(ROUTINE,
        VARIANT, TILEID) MAD summary AND the individual surviving rows
        in each group (for a per-row scatter view, not just the summary
        statistic) -- see that function's own docstring for the
        tileid-exclusion/position-clustering filtering it applies.

        The report has one real HTML heading per section, each wrapping
        its OWN small Plotly figure (built by build_figures() below) --
        deliberately NOT one giant multi-row Plotly canvas with
        hand-computed pixel margins standing in for section breaks. That
        approach kept needing another manually-tuned margin/spacer-row
        fix every time a new section was added; real HTML block flow
        (headings with their own margins) reserves space between
        sections automatically and can't overlap:

            <h1>title</h1>
            [run-summary table]
            <h2>Line Ratios</h2>
            <h3>Ratio vs. Line Flux</h3>
            [one row per doublet, one column per method]
            <h3>Repeat-Group Scatter (MAD) vs. Median Flux</h3>
            [same grid, y = that group's ratio MAD]
            <h2>Repeat-Observation Flux Consistency</h2>
            <h3>Repeat-Group Fractional Scatter vs. Median Flux</h3>
            [one row per FLUX_METRICS entry, one column per method;
             x = a repeat group's median flux, y = that group's
             fractional MAD (robust MAD / median) of the same flux
             across its repeat exposures]
            <h3>Summary</h3>
            [interactive MAD/fractional-MAD table]
            <h2>Continuum Residual (B/R/Z)</h2>
            <h3>Post-Subtraction Continuum Level (All Exposures)</h3>
            [one row per spectrograph arm, one column per method; every
             exposure, not just repeat groups; x = that exposure's own
             pre-subtraction continuum median, y = its post-subtraction
             continuum median (linear, unclipped) with a dashed red
             line at zero -- a physical continuum cannot go negative]
            [percent-negative summary table]
            <h3>Post-Subtraction Continuum Level Consistency</h3>
            [repeat groups only; one row per spectrograph arm, one
             column per method; x = a repeat group's PRE-subtraction
             continuum median in that arm (shared across methods), y =
             that group's fractional MAD of the POST-subtraction
             continuum median, normalized against the PRE-subtraction
             brightness (log y) -- NOT against its own median, which
             would be unstable near zero]
            <h3>Post-Subtraction Continuum RMS Consistency</h3>
            [same grid; y = fractional MAD of the POST-subtraction
             continuum RMS/NMAD, normalized against its own median
             (safe -- a noise floor is never near zero)]
            <h3>Summary</h3>
            [interactive fractional-MAD table, level and RMS]

        Each section/subsection is preceded by a short <p> caption
        explaining what is plotted -- added 260911 so the report reads
        standalone without needing this docstring open alongside it.

        This script produces exactly one HTML file, never one file per
        exposure -- the per-exposure picture is PlotSkySubNebEval.py's
        job; this one is the aggregate/distribution view across an
        entire run.

Primary routines:

    build_run_tables   the two input modes -> (row_table, repeats_table,
                       group_rows), reusing SkySubNebEval.py throughout.
    build_figures      -> dict of the small per-section Plotly figures
                       described above (values are None where there
                       were no repeat groups to plot).
    build_continuum_table, build_continuum_run_tables, build_continuum_figures
                       fits_files (raw only, no -lines_file equivalent)
                       -> per-arm B/R/Z pre-/post-subtraction continuum
                       table -> repeat-tileid groups (reusing
                       SkySubNebEval.group_repeat_exposures) -> the
                       Continuum Residual section's figures (the all-
                       exposures negative-continuum check plus the two
                       repeat-only consistency panels). See
                       build_continuum_table's own docstring for why
                       one function (GetSkyCont.arm_continuum_stats on
                       raw flux, no local fit) works for both the
                       pre- and post-subtraction measurement.

History::

        260910  ksl  Coding begun: build_run_tables() (two input modes,
                     reusing SkySubNebEval.fit_file/repeat_scatter) and
                     build_figures()/write_report() (the real-HTML,
                     multi-figure report described above).

                     Two earlier designs were tried and replaced entirely
                     rather than incrementally patched, each for a
                     concrete reason:

                     1. Per-method-category box plots (one box per method
                        within a shared panel) were replaced by the
                        scatter-vs-flux design: pooling every row/group
                        into one box hid the dependence of ratio/flux
                        scatter on how bright the line actually was, and
                        a column-per-repeat-group layout does not scale
                        past a handful of groups.
                     2. A single giant multi-row Plotly canvas (one
                        figure holding every table and panel, with
                        in-plot "spacer" annotation rows standing in for
                        section headings) was replaced by build_figures()
                        -> several focused per-section figures assembled
                        into a real HTML document (actual <h1>/<h2>/<h3>
                        tags, CSS margins between sections) by
                        write_report(). The recurring overlap/crowding
                        bugs along the way (title vs. spectrum, legend
                        vs. title, spacer-row vs. content) were all
                        symptoms of one root cause -- reserving vertical
                        space for non-plot text with Plotly margins/
                        annotations inside a single canvas instead of
                        normal HTML block flow, which reserves space
                        between elements automatically and cannot
                        overlap.

                     Also found and fixed a real bug in
                     SkySubNebEval.repeat_scatter() while building the
                     per-repeat-group panels: its new return_rows option
                     (Table(rows=<Row objects>)) was silently replacing
                     masked/NaN entries with 0.0 (see that file's own
                     History) -- caught by cross-checking this report's
                     MAD-summary numbers against
                     SkySubNebEval.print_comparison()'s printed values.

                     Design conventions settled on: column order is
                     input-file order (_ordered_routines), not
                     alphabetical; one row per metric, one column per
                     method (not all methods overlaid in one panel,
                     which stops being legible past a handful of points);
                     shared x/y range per metric row (1st/99th
                     percentile, not literal min/max); axis/table labels
                     are bare line names for flux totals ('OII') and the
                     actual wavelength pair for ratios ('OII 3730:3726'),
                     computed live from NEBULAR_LINES; low-density-limit
                     reference lines for OII/SII (dashed, not a truth/
                     ceiling -- density-sensitive ratios do vary, this is
                     just the value observed at most average sky
                     positions).

        260911  ksl  Nebular-line repeat-consistency panels reworked, all
                     found by direct inspection of the rendered report:

                     Repeat-Observation Flux Consistency replotted --
                     it used to show the raw per-exposure flux cloud
                     (every exposure's own value against its group's
                     median), the only panel in the report with no
                     collapsed dispersion number. Now one point per
                     repeat group, same design as the ratio-MAD panel,
                     plotting fractional MAD (robust MAD / group
                     median) -- flux spans many decades across lines/
                     methods, so only the fractional error is
                     comparable panel to panel (same reasoning
                     _mad_summary_figure's flux columns already used).

                     x-axis brightness proxy unified: the ratio panels
                     used to plot the MEAN of a doublet's two member
                     fluxes while the flux panel used a differently-
                     defined per-metric total (both members summed for
                     OII/SII, only the dominant member for NII/OIII/
                     SIII) -- the same line sat at a different x-scale
                     depending which panel you looked at (2x for OII/
                     SII, ~1.2-1.3x for NII/SIII/OIII). _build_metrics()
                     now carries proxy_names (always the line's full
                     DOUBLETS pair, used for x everywhere), decoupled
                     from flux_names (what's actually measured on y --
                     unchanged, still just the dominant member for
                     NII/OIII/SIII, to avoid adding the weak member's
                     noise to that measurement). x is the TOTAL (not
                     mean) of both members throughout, so 'Ratio vs.
                     Line Flux' relabeled 'Mean line flux' ->
                     'Total line flux'. Beyond the formula, the actual
                     RANGE shown could still differ panel to panel
                     (~8e-15 vs ~3e-14 for the same line) since the
                     repeat-only panels only have ~10-30 pooled points
                     to compute a 1st/99th percentile from -- too few
                     to reliably trim outliers. Fixed by having
                     _ratio_vs_flux_figure (thousands of rows, the
                     stable estimate) return the range it computed per
                     line and having the repeat-only panels force that
                     same range (_shared_row_range's forced_xr) instead
                     of recomputing their own.

                     Every section/subsection also gained a short <p>
                     caption stating what's plotted (self-explanatory
                     without this docstring open alongside it) and,
                     since the first section covers ~20x more rows than
                     the repeat-only ones, an explicit repeat-tileid-
                     group/exposure count where the report switches
                     from the full sample to CLOSE groups only -- with
                     its own divider line (h3.subsection-break), so a
                     line with a low per-row detection rate (e.g. SIII,
                     in the z-arm's densest OH forest) showing very few
                     points from that point on is not a mystery.
                     x-axis tick labels also render at a fixed 45
                     degrees everywhere now, not Plotly's row-dependent
                     'auto' (which picked a different angle per row
                     depending on that row's own tick-label lengths).

        260911  ksl  New Continuum Residual (B/R/Z) section: tests the
                     continuum (not emission-line) flux left behind in
                     each spectrograph arm after sky subtraction --
                     planned at length in this session, as a set of
                     conversation turns, before any code was written.
                     Three subsections, in order:

                     Post-Subtraction Continuum Level (All Exposures) --
                     a physical continuum flux cannot be negative (a
                     negative value means too much sky was removed), a
                     check that needs no repeat observation at all, so
                     it runs on every exposure rather than just the
                     repeat-tileid groups the next two subsections use.
                     One row per arm, a dashed red line at y=0, and a
                     percent-negative summary table.

                     Post-Subtraction Continuum Level Consistency and
                     RMS Consistency -- repeat-tileid groups only (via
                     SkySubNebEval.group_repeat_exposures(), repeat_
                     scatter's own grouping logic factored out so this
                     differently-shaped per-arm table can reuse it
                     without needing NEBULAR_LINES/DOUBLETS columns it
                     doesn't have). Both plot fractional MAD of a post-
                     subtraction quantity (median for Level, NMAD for
                     RMS) against the group's PRE-subtraction continuum
                     brightness (FLUX+SKY -- SkySub*.py's own
                     convention that FLUX is already sky-subtracted and
                     SKY is the model removed, so FLUX+SKY reconstructs
                     the original CFrame spectrum with no need to
                     re-open the input file; identical across all 5
                     methods since they share one input CFrame).
                     Deliberately NOT SkySubOrig.py's existing pre-
                     subtraction SCI_MED_B/etc. columns, which are a
                     continuum-FIT-quality diagnostic (raw flux minus a
                     locally-fit polynomial), not a brightness
                     measurement, and never computed post-subtraction
                     by any of the 5 methods.

                     The two panels normalize their fractional MAD
                     differently, on purpose: RMS divides by the
                     group's own median NMAD (safe -- a noise floor is
                     never near zero). Level instead divides by the
                     group's PRE-subtraction brightness, NOT by its own
                     post-subtraction median -- a post-subtraction
                     level is expected to sit near zero for a good
                     method, so self-normalizing is unstable and can
                     invert the ranking entirely (confirmed on real
                     data: the single worst-looking point under the
                     naive metric had a median of -1.4e-17, essentially
                     perfect, while the best-looking point had the
                     LARGEST systematic residual in the set). Once
                     normalized against pre-subtraction brightness
                     instead, the fractional values span a genuine ~2
                     decades across groups, so the Level panel's y-axis
                     is log-scale (a log axis alone, tried first, does
                     NOT fix the self-normalization instability -- it
                     only compresses the same wrong ordering).
                     _continuum_repeat_figure takes normalize_by
                     ('own'/'pre') and log_y for this;
                     _shared_row_range gained matching log_y support.

                     Underlying arm definition is GetSkyCont.
                     ARM_EVAL_RANGES (B/R/Z with overlap zones and outer
                     edges excluded) -- the canonical definition
                     SkySubOrig.py already uses, distinct from two
                     other, unrelated arm-range definitions elsewhere in
                     the repo. GetSkyCont.arm_continuum_stats is the
                     one function used for every measurement (pre- and
                     post-subtraction alike, no local continuum fit):
                     data/sky_mask.fits already excludes every sky-
                     line-dominated pixel (built from an actual PALACE
                     sky-emission model, not a crude window list), so
                     the median/NMAD of an arm's clean pixels already
                     are the brightness/scatter estimates needed.

        260911  ksl  Two rendering bugs in _grid_figure (shared by every
                     scatter-grid figure, nebular-line ones included),
                     found against the rendered continuum panels:
                     vertical_spacing was a hardcoded FRACTION (0.05)
                     of total figure height, so the actual pixel gap
                     between rows shrank on a short figure -- the 3-row
                     continuum panels got only ~50px (not enough room
                     for a tickangle=45 tick-label stack plus the next
                     row's subplot-title annotation without overlap)
                     while the 6-row nebular panels got ~92px from the
                     same fraction and looked fine. Now computed from a
                     fixed pixel target (new ROW_GAP_PX=90) so every
                     grid gets the same absolute gap regardless of row
                     count. Separately, subplot_titles used to repeat
                     the method-label row above EVERY row instead of
                     just the first -- redundant (column order is fixed
                     for the whole figure) and ate into that same
                     vertical space; now only row 1 carries labels.
                     Also shortened the continuum panels' x-axis title
                     ('Pre-subtraction continuum (median)', 35
                     characters -- long enough to collide with its own
                     copy in the next ~260px-wide column) to 'Pre-sub
                     continuum' (18 characters, matching every other
                     x-axis title's length in this module).
        260913  ksl  chmod +x -- had a shebang and __main__ block but
                     was missing the executable bit.

'''

import argparse
import re
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from SkySubNebEval import (fit_file, repeat_scatter, group_repeat_exposures,
                           DOUBLETS, NEBULAR_LINES, PLACEHOLDER_TILEID)
from sky_nebular_leak_eval import LMC_VEL, SMC_VEL
from GetSkyCont import (load_mask, _interp_mask_to_wave, arm_continuum_stats,
                        ARM_EVAL_RANGES)

METHOD_COLORS = ['black', 'firebrick', 'seagreen', 'steelblue', 'darkorange', 'purple']
SUBPLOT_TITLE_SIZE = 14
TABLE_ROW_HEIGHT = 300
GROUP_ROW_HEIGHT = 280
MARGIN_T = 90
MARGIN_B = 60
# Target ABSOLUTE pixel gap between grid rows (260911) -- vertical_spacing
# is a fraction of the whole figure height, so a fixed fraction (the
# previous hardcoded 0.05) gives very different actual pixel gaps on a
# short 3-row figure (e.g. the continuum panels) vs. a tall 6-row one
# (the nebular-line panels): 0.05 on a ~990px 3-row figure is only ~50px,
# not enough room for a rotated (tickangle=45) tick label stack plus the
# next row's subplot-title annotation without them overlapping -- found
# directly from the rendered report. Computing vertical_spacing from this
# fixed pixel target instead keeps the gap visually consistent regardless
# of how many rows a given figure has.
ROW_GAP_PX = 90

# Spectrograph arms, in the same order/labels as GetSkyCont.ARM_EVAL_RANGES
# (the project-wide canonical B/R/Z definition -- see build_continuum_table).
ARM_LABELS = [label for label, _lo, _hi in ARM_EVAL_RANGES]

# Default clean-pixel mask for the continuum-residual repeat test, same
# file/convention SkySubOrig.py and sky_residual_eval.py already default to.
_DEFAULT_MASK_FILE = Path(__file__).parent.parent / 'data' / 'sky_mask.fits'

# (name, [NEBULAR_LINES member names to sum]) -- total line flux for
# repeat-observation consistency checks: ratio consistency alone can't
# catch a method that shifts both lines' absolute flux together (e.g.
# via continuum over/under-subtraction right at the line). NII/OIII/
# SIII use only their stronger member (nii_b/oiii_b/siii_b) -- the
# weaker member is a fixed fraction of it (DOUBLETS' fixed ratio), so
# summing would add noise, not information. OII/SII sum both members
# since neither doublet has one member dominating the way the fixed-
# ratio pairs do.
FLUX_METRICS = [
    ('OII_FLUX',  ['oii_a', 'oii_b']),
    ('OIII_FLUX', ['oiii_b']),
    ('NII_FLUX',  ['nii_b']),
    ('SII_FLUX',  ['sii_a', 'sii_b']),
    ('SIII_FLUX', ['siii_b']),
]


def _ordered_routines(row_table):
    '''Distinct (ROUTINE, VARIANT) pairs in first-occurrence order, i.e.
    the order the input files were actually given/created in -- NOT
    alphabetically sorted. Columns/rows built from this should read in
    that same order (the order the user ran the routines), not whatever
    order string-sorting happens to produce (found confusing: "Dev1,
    Dev2, Dev3, Drp, Orig" from sorting -- arbitrary relative to how the
    run was actually done). Every figure in this module is built from
    the SAME (labels, routines) pair returned by _common_setup, so
    column order is identical across every section by construction.'''
    if len(row_table) == 0:
        return []
    return list(dict.fromkeys(zip(row_table['ROUTINE'], row_table['VARIANT'])))


def _total_flux(row_table, names, snr_min):
    '''
    Row-by-row total flux summed over `names` NEBULAR_LINES entries,
    gated on EVERY listed line's own SNR (same per-row reliability gate
    SkySubNebEval._doublet_ratio uses for DOUBLETS) -- NaN, not a partial
    sum, if any component fails the gate, so "total" here always means
    every component was individually usable, not just whichever ones
    happened to pass.
    '''
    n = len(row_table)
    total = np.zeros(n)
    ok = np.ones(n, dtype=bool)
    for name in names:
        flux = np.asarray(row_table[f'{name}_FLUX'], dtype=float)
        snr = np.asarray(row_table[f'{name}_SNR'], dtype=float)
        ok &= np.isfinite(flux) & np.isfinite(snr) & (snr >= snr_min)
        total += np.nan_to_num(flux)
    total = total.copy()
    total[~ok] = np.nan
    return total


def _usage_from_doc(doc):
    '''__doc__ truncated just before "History:", so -h stays short.'''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


def build_run_tables(fits_files=None, lines_file=None, vel=None, lmc=False, smc=False,
                     sigma_guess=1.0, snr_min=5.0, mjd_close=7.0):
    '''
    Assemble the row-level and repeat-group tables for one run, from
    either raw SkySub*.py output files or a precomputed
    SkySubNebEval.py _lines.fits file.

    Parameters
    ----------
    fits_files : sequence of str or Path, optional
        Raw SkySub*.py output files. Mutually exclusive with lines_file.
    lines_file : str or Path, optional
        An existing SkySubNebEval.py run's <root>_lines.fits. Mutually
        exclusive with fits_files.
    vel, lmc, smc, sigma_guess, snr_min : see SkySubNebEval.fit_file.
        Ignored (a ValueError is not raised, but they have no effect) if
        lines_file is given instead of fits_files.
    mjd_close : float
        Passed to repeat_scatter.

    Returns
    -------
    row_table : astropy.table.Table
        One row per spectrum (see SkySubNebEval.fit_file).
    repeats_table : astropy.table.Table
        One row per surviving (ROUTINE, VARIANT, TILEID) repeat group.
    group_rows : dict
        {(routine, variant, tileid): astropy.table.Table} of the
        individual rows in each surviving group.
    '''
    if (fits_files is None) == (lines_file is None):
        raise ValueError('exactly one of fits_files or lines_file is required')

    if lines_file is not None:
        row_table = Table.read(lines_file)
    else:
        tables = [fit_file(f, vel=vel, lmc=lmc, smc=smc,
                           sigma_guess=sigma_guess, snr_min=snr_min)
                 for f in fits_files]
        row_table = vstack(tables) if tables else Table()

    repeats_table, group_rows = repeat_scatter(row_table, mjd_close=mjd_close, return_rows=True)
    return row_table, repeats_table, group_rows


def build_continuum_table(fits_files, mask_file=None):
    '''
    Per-exposure pre- and post-subtraction continuum brightness/scatter
    in each spectrograph arm (B/R/Z) -- the data side of the continuum-
    residual repeat-consistency test (see module docstring's 260911
    History entry).

    One measurement, reused verbatim for both quantities:
    GetSkyCont.arm_continuum_stats(wave, flux, clean_mask) called
    directly on RAW flux, with no local continuum fit/detrending first.
    That works here specifically because clean_mask (from
    data/sky_mask.fits, a palace_make_mask.py-built mask driven by an
    actual PALACE sky-emission model) already excludes every sky-line-
    dominated pixel -- the OH forest, atomic lines, the O2 A-band -- so
    the median of whatever flux remains in an arm's clean pixels IS
    already a continuum-brightness estimate, and its NMAD IS already a
    scatter estimate, without needing a separate polynomial fit the way
    SkySubOrig.py's existing (pre-subtraction, fit-QUALITY, not
    brightness) SCI_MED_B/etc. columns do -- see this session's own
    investigation into why those columns aren't reusable for this.

    The two quantities::

        PRE  : FLUX + SKY. SkySub*.py's own convention is that FLUX is
               already sky-subtracted and SKY is the model subtracted
               from it, so FLUX+SKY reconstructs the ORIGINAL, pre-
               subtraction CFrame science spectrum -- no need to
               re-open the input CFrame file at all. This is identical
               across all 5 methods (they share one input CFrame), so
               it's computed here for every row of every file, but
               downstream (_continuum_pre_x_by_group) only the first
               routine's copy is actually used, the same "one
               representative method" pattern build_figures already
               uses for n_close_groups.
        POST : FLUX directly -- what a given method actually left
               behind.

    Parameters
    ----------
    fits_files : sequence of str or Path
        SkySub*.py output FITS files (WAVE, FLUX, SKY, DRP_ALL) -- the
        same contract build_run_tables' fits_files takes. Unlike that
        path, there is no -lines_file equivalent here: the full WAVE/
        FLUX/SKY arrays are required, not just per-row scalar fit
        results, so this always needs the raw files.
    mask_file : str or Path, optional
        Defaults to _DEFAULT_MASK_FILE (data/sky_mask.fits), the
        project-wide standard.

    Returns
    -------
    astropy.table.Table
        One row per spectrum: FILE, ROUTINE, VARIANT, ROW, TILEID, MJD,
        SCI_RA, SCI_DEC, then PRE_MED_<ARM>/PRE_NMAD_<ARM>/
        POST_MED_<ARM>/POST_NMAD_<ARM> for ARM in ARM_LABELS.
    '''
    mask_file = mask_file or _DEFAULT_MASK_FILE
    mask_wave, mask_arr = load_mask(str(mask_file))

    rows = []
    for f in fits_files:
        with fits.open(f) as hdul:
            hdr = hdul[0].header
            routine = hdr.get('TITLE', '')
            variant = hdr.get('METHOD', '')
            wave = np.asarray(hdul['WAVE'].data, dtype=float)
            flux = np.asarray(hdul['FLUX'].data, dtype=float)
            sky = np.asarray(hdul['SKY'].data, dtype=float)
            drp = Table(hdul['DRP_ALL'].data)
        clean = _interp_mask_to_wave(mask_wave, mask_arr, wave)

        for i in range(flux.shape[0]):
            post_stats = arm_continuum_stats(wave, flux[i], clean)
            pre_stats = arm_continuum_stats(wave, flux[i] + sky[i], clean)
            row = dict(FILE=str(f), ROUTINE=routine, VARIANT=variant, ROW=i,
                      TILEID=int(drp['tileid'][i]) if 'tileid' in drp.colnames else -1,
                      MJD=float(drp['mjd'][i]) if 'mjd' in drp.colnames else np.nan,
                      SCI_RA=float(drp['sci_ra'][i]) if 'sci_ra' in drp.colnames else np.nan,
                      SCI_DEC=float(drp['sci_dec'][i]) if 'sci_dec' in drp.colnames else np.nan)
            for arm in ARM_LABELS:
                row[f'PRE_MED_{arm}']  = pre_stats[arm]['med']
                row[f'PRE_NMAD_{arm}'] = pre_stats[arm]['nmad']
                row[f'POST_MED_{arm}']  = post_stats[arm]['med']
                row[f'POST_NMAD_{arm}'] = post_stats[arm]['nmad']
            rows.append(row)

    return Table(rows=rows)


def build_continuum_run_tables(fits_files, mjd_close=7.0, mask_file=None):
    '''fits_files -> (continuum_table, cont_meta_rows, cont_group_rows),
    the continuum analog of build_run_tables -- build_continuum_table()
    for the data, then SkySubNebEval.group_repeat_exposures() for the
    SAME tileid-exclusion/position-clustering/CLOSE-flag grouping the
    nebular-line panels use (reused, not reimplemented -- see that
    function's 260911 History entry). cont_meta_rows is the list-of-
    dict form group_repeat_exposures returns directly (ROUTINE/VARIANT/
    TILEID/N_REPEATS/MJD_SPAN/CLOSE/POS_SCATTER_DEG per group) -- there
    is no line-specific MED/MAD to add on top here the way
    repeat_scatter adds for the nebular row_table, so this table is
    used as-is rather than wrapped in an astropy Table. No -lines_file
    fast path here (see build_continuum_table); always requires the
    raw fits_files.'''
    cont_table = build_continuum_table(fits_files, mask_file=mask_file)
    cont_meta_rows, cont_group_rows = group_repeat_exposures(cont_table, mjd_close=mjd_close)
    return cont_table, cont_meta_rows, cont_group_rows


def _summarize(row_table, repeats_table, group_rows):
    '''Print what would be plotted -- a quick sanity check on the data
    side, independent of anything Plotly-related.'''
    print(f'\nrow_table: {len(row_table)} rows, {len(row_table.colnames)} columns')
    if len(row_table) == 0:
        return
    routines = _ordered_routines(row_table)
    print(f'  methods ({len(routines)}): {routines}')

    print(f'\nDoublet SNR-pass counts (out of {len(row_table)} rows), per method:')
    header = f"  {'Method':<32}" + ''.join(f'{d.upper():>14}' for d, *_ in DOUBLETS)
    print(header)
    for routine, variant in routines:
        sel = (row_table['ROUTINE'] == routine) & (row_table['VARIANT'] == variant)
        line = f"  {routine + '/' + variant:<32}"
        for dname, *_ in DOUBLETS:
            vals = np.asarray(row_table[dname.upper()])[sel]
            npass = int(np.isfinite(vals).sum())
            line += f'{npass:>14}'
        print(line)

    n_real_groups = int(np.sum(repeats_table['TILEID'] != PLACEHOLDER_TILEID)) if len(repeats_table) else 0
    print(f'\nrepeats_table: {len(repeats_table)} groups ({n_real_groups} real tileid groups)')
    print(f'group_rows: {len(group_rows)} entries')
    for key, sub in list(group_rows.items())[:3]:
        print(f'  {key}: {len(sub)} rows')
    if len(group_rows) > 3:
        print(f'  ... and {len(group_rows) - 3} more')


def _hline(fig, row_i, col_j, xrange, value, color='grey'):
    '''Plain go.Scatter reference line spanning xrange, not
    fig.add_hline/add_shape's row=/col= convenience: that resolution
    walks every subplot to compute axis anchors and raises a
    PlotlyKeyError as soon as it reaches a go.Table in the SAME figure.
    None of the figures built here mix a Table with scatter subplots any
    more (each section is its own figure), so this is no longer strictly
    required for that reason, but it's kept for consistency/robustness.'''
    if value is None or xrange is None:
        return
    fig.add_trace(go.Scatter(x=list(xrange), y=[value, value], mode='lines',
                             line=dict(color=color, width=2.5, dash='dot'),
                             showlegend=False, hoverinfo='skip'),
                  row=row_i, col=col_j)


def _shared_row_range(fig, n_cols, row_i, all_x, all_y, forced_xr=None, log_y=False):
    '''Pool x/y across every method-column in one metric row and apply
    one shared range to all of them (1st/99th percentile, not literal
    min/max -- a single noisy method's outlier point would otherwise set
    a shared range wide enough to flatten every other column's real
    signal, the same lesson PlotSkySubNebEval.py's per-line-group rows
    already learned). x is always log-flux; y is linear by default,
    log if log_y=True (260911, for the continuum-level-consistency
    panel -- a fractional-MAD quantity that genuinely spans ~2 decades
    once normalized against a stable denominator, see
    _continuum_repeat_figure's normalize_by). Returns xr (linear-space
    tuple) for _hline.

    The two extra parameters::

        forced_xr : (lo, hi) or None
            When given, use this range verbatim instead of computing
            one from all_x -- added 260911 so the ratio-MAD/flux-MAD
            panels can reuse the SAME x-range _ratio_vs_flux_figure
            already computed from the full sample for this line,
            rather than recomputing their own 1st/99th percentile
            from just the handful of repeat-group points they have.
            With as few as ~10-30 pooled group-level points for a
            line like SIII, a 1st/99th percentile barely trims
            anything (that's already close to min/max at low n), so
            the SAME line's "shared" x-range could come out very
            different panel to panel purely from which few repeat
            groups happened to survive -- not a real brightness
            difference (found by direct comparison this session:
            e.g. an ~8e-15 vs ~3e-14 upper limit for the same line,
            panel to panel). Reusing one authoritative range (from
            the thousands-of-rows full sample, the most stable
            estimate available) makes a line's x-axis window
            identical everywhere it appears, by construction rather
            than by coincidence. y still gets its own per-panel
            percentile below regardless -- MAD/fractional-MAD are
            genuinely different quantities panel to panel and are
            not part of this fix.
        log_y : bool
            True to make y a log axis instead of linear -- see this
            function's own log_y description above for when/why.'''
    if forced_xr is not None:
        xr = forced_xr
    else:
        xr = tuple(np.percentile(all_x, [1, 99])) if all_x.size else None
        if xr and xr[0] >= xr[1] and xr[0] > 0:
            xr = (xr[0] / 3.0, xr[1] * 3.0)  # degenerate (identical/single point)
    pos_y = all_y[all_y > 0] if (log_y and all_y.size) else None
    for col_j in range(1, n_cols + 1):
        if xr and xr[0] > 0:
            fig.update_xaxes(range=[np.log10(xr[0]), np.log10(xr[1])], row=row_i, col=col_j)
        if log_y:
            if pos_y is not None and pos_y.size:
                ylo, yhi = np.percentile(pos_y, [1, 99])
                if ylo >= yhi:
                    ylo, yhi = ylo / 3.0, yhi * 3.0
                fig.update_yaxes(type='log', range=[np.log10(ylo), np.log10(yhi)], row=row_i, col=col_j)
            else:
                fig.update_yaxes(type='log', row=row_i, col=col_j)
        elif all_y.size:
            ylo, yhi = np.percentile(all_y, [1, 99])
            pad = 0.08 * (yhi - ylo) if yhi > ylo else max(abs(ylo), 1e-30) * 0.1
            fig.update_yaxes(range=[ylo - pad, yhi + pad], row=row_i, col=col_j)
    return xr


_LINE_WAVE = {name: center for name, center, _wmin, _wmax in NEBULAR_LINES}


def _wave_key(m):
    '''Shortest-member rest wavelength of a metric's constituent lines --
    standing convention in this project (also how SkySubNebEval.DOUBLETS
    itself is ordered): panels/rows/table columns are always presented
    in increasing wavelength, not declaration order, so this is computed
    from NEBULAR_LINES directly rather than relying on DOUBLETS/
    FLUX_METRICS happening to already be listed that way.'''
    names = m['flux_names']
    return min(_LINE_WAVE[n] for n in names)


def _build_metrics():
    '''
    One combined list describing every panel this figure plots, sorted
    by increasing wavelength (see _wave_key) -- ratio and flux metrics
    are interleaved by their own line wavelength, not grouped by kind,
    so e.g. OII_FLUX sits next to OII_RATIO rather than after every
    ratio metric regardless of its own (much longer) wavelength.

    Each entry: dict(label, kind, ratio_col, flux_names, proxy_names,
    ref_kind, ref_value).

    proxy_names=(name_a, name_b) is ALWAYS that line's full DOUBLETS
    pair -- the x-axis brightness proxy used by every panel this metric
    appears in (see _metric_arrays), so a given line sits at the same
    x-scale (total of both members) in the ratio-vs-flux, ratio-MAD,
    and flux-MAD sections alike. flux_names is the (possibly different,
    for kind='flux') set of members that actually go into the plotted
    y-VALUE.

    kind='ratio' (one per DOUBLETS entry): flux_names == proxy_names
    (the doublet's own two members); ratio_col is both the ratio-vs-
    flux y-value and the repeat-group MAD's dispersion source.

    kind='flux' (one per FLUX_METRICS entry): ratio_col is None --
    there is no "ratio vs its own flux" panel (that would just be a
    trivial 1:1 relation), so these only appear in the flux-MAD
    section. flux_names there is deliberately NOT always the same as
    proxy_names -- NII/OIII/SIII's plotted VALUE uses only the
    dominant member (see FLUX_METRICS comment: adding the fixed-ratio
    weak member would add noise, not information, to that
    measurement), but proxy_names still pulls in both members so the
    x-axis brightness scale matches the other two panels regardless
    (260911: replaces flux_names doing double duty for both x and y,
    which is what caused the same line to sit at different x-scales in
    different panels -- see this session's own investigation).
    '''
    doublet_pair = {d.rsplit('_ratio', 1)[0]: (name_a, name_b)
                    for d, name_a, name_b, _kind, _value, _c in DOUBLETS}
    metrics = [dict(label=d.upper(), kind='ratio', ratio_col=d.upper(),
                    flux_names=(name_a, name_b), proxy_names=(name_a, name_b),
                    ref_kind=kind, ref_value=value)
              for d, name_a, name_b, kind, value, _c in DOUBLETS]
    metrics += [dict(label=name, kind='flux', ratio_col=None, flux_names=names,
                     proxy_names=doublet_pair[name.rsplit('_FLUX', 1)[0].lower()],
                     ref_kind=None, ref_value=None)
               for name, names in FLUX_METRICS]
    metrics.sort(key=_wave_key)
    return metrics


def _short_label(m):
    '''Bare display name for a metric, e.g. 'OII', 'SIII' -- no '_RATIO'/
    '_FLUX' suffix (requested explicitly: the suffix was noise, not
    information, once the section heading already says which one it is).'''
    return m['label'].replace('_RATIO', '').replace('_FLUX', '')


def _ratio_wave_label(m):
    '''Display name for a RATIO metric as its actual wavelength pair,
    e.g. 'OII 3730:3726' -- unambiguous about which line is the
    numerator (unlike a bare name, which doesn't say), computed from
    NEBULAR_LINES so it stays correct if a doublet's orientation is ever
    changed (as SII's was, to match sii_a/sii_b convention) rather than
    a hand-maintained string.'''
    name_a, name_b = m['flux_names']
    return f'{_short_label(m)} {round(_LINE_WAVE[name_b])}:{round(_LINE_WAVE[name_a])}'


def _row_label(m):
    '''The one function used everywhere a metric needs a display name
    (axis titles, table headers) -- wavelength-pair notation for a
    ratio, bare line name for a flux total.'''
    return _ratio_wave_label(m) if m['kind'] == 'ratio' else _short_label(m)


def _metric_arrays(sub, m, snr_min):
    '''
    Row-aligned (flux_proxy, value) arrays for one metric on one table
    (either the whole row_table, filtered to one method, or one repeat-
    group's sub-table) -- gated to rows where both are finite and
    flux_proxy > 0 (needed for the log x-axis either way).

    flux_proxy is the x-axis brightness measure: ALWAYS the total flux
    (sum, not mean) of m['proxy_names'] -- that line's full DOUBLETS
    pair -- so the same line shares one x-scale across every panel it
    appears in (260911, replacing a design where the ratio panels used
    the pair's MEAN while the flux-MAD panel used a different,
    metric-specific total: see _build_metrics' docstring).

    value is the y-axis quantity: the ratio itself for a 'ratio'
    metric, or flux_names' own SNR-gated total (_total_flux) for a
    'flux' metric -- unchanged, and deliberately independent of
    flux_proxy's proxy_names above (a 'flux' metric's y can still be
    just the dominant member while its x uses both).
    '''
    name_a, name_b = m['proxy_names']
    fa = np.asarray(sub[f'{name_a}_FLUX'], dtype=float)
    fb = np.asarray(sub[f'{name_b}_FLUX'], dtype=float)
    flux_proxy = fa + fb
    if m['kind'] == 'ratio':
        value = np.asarray(sub[m['ratio_col']], dtype=float)
    else:
        value = _total_flux(sub, m['flux_names'], snr_min)
    ok = np.isfinite(flux_proxy) & (flux_proxy > 0) & np.isfinite(value)
    return flux_proxy[ok], value[ok]


def _common_setup(row_table, repeats_table, group_rows):
    '''Everything every figure-builder below needs, computed once so
    column order/colors/metric lists are identical across every section
    by construction (rather than each builder re-deriving them and
    risking drift between sections).'''
    routines = _ordered_routines(row_table)
    # Short label (routine with the common 'SkySub' prefix stripped) --
    # method names in every run seen so far (Orig/Drp/Dev1/Dev2/Dev3)
    # collapse to short, distinct labels this way; variant is dropped
    # since one run normally has only one variant per routine.
    labels = [routine.replace('SkySub', '') for routine, _variant in routines]
    colors = {label: METHOD_COLORS[i % len(METHOD_COLORS)] for i, label in enumerate(labels)}
    tileids = sorted(set(tileid for (_r, _v, tileid) in group_rows))
    close_lookup = {(r['ROUTINE'], r['VARIANT'], int(r['TILEID'])): bool(r['CLOSE'])
                    for r in repeats_table} if len(repeats_table) else {}
    metrics = _build_metrics()
    return dict(routines=routines, labels=labels, colors=colors, tileids=tileids,
               close_lookup=close_lookup, metrics=metrics,
               ratio_metrics=[m for m in metrics if m['kind'] == 'ratio'],
               flux_metrics=[m for m in metrics if m['kind'] == 'flux'],
               have_groups=bool(tileids), n_cols=max(len(labels), 1))


def _table_figure(header, cols, width=None):
    '''A standalone go.Table figure -- no make_subplots needed for a
    single table, and no risk of the add_hline/go.Table interaction
    issue since nothing else shares this figure.'''
    fig = go.Figure(data=[go.Table(
        header=dict(values=header, font=dict(size=SUBPLOT_TITLE_SIZE), align='center'),
        cells=dict(values=cols, font=dict(size=SUBPLOT_TITLE_SIZE), align='center', height=28))])
    fig.update_layout(height=TABLE_ROW_HEIGHT, width=width or max(900, 130 * len(header)),
                      margin=dict(t=10, b=10, l=10, r=10))
    return fig


def _run_summary_figure(row_table, setup):
    '''Row 1 (of the whole report): rows fit and SNR-pass count/fraction
    per doublet, one row per method.'''
    labels, routines, ratio_metrics = setup['labels'], setup['routines'], setup['ratio_metrics']
    header = ['Method', 'N rows'] + [_row_label(m) for m in ratio_metrics]
    table_cols = [labels,
                 [int(np.sum((row_table['ROUTINE'] == r) & (row_table['VARIANT'] == v)))
                  for r, v in routines]]
    for m in ratio_metrics:
        col = []
        for r, v in routines:
            sel = (row_table['ROUTINE'] == r) & (row_table['VARIANT'] == v)
            n = int(sel.sum())
            npass = int(np.isfinite(np.asarray(row_table[m['ratio_col']], dtype=float)[sel]).sum())
            col.append(f'{npass} ({100 * npass / n:.0f}%)' if n else '--')
        table_cols.append(col)
    return _table_figure(header, table_cols, width=max(1000, 150 * len(header)))


def _grid_figure(n_rows, n_cols, labels):
    '''Shared make_subplots scaffold for every scatter-grid section in
    this module -- one row per metric, one column per method (see
    module docstring for why: overlaying every method in one panel
    stops being legible past a handful of points; a single-method
    panel stays readable regardless of how large the run is).

    subplot_titles carries the method-name labels ONLY for row 1, not
    repeated on every row (260911) -- they're the same 5 names every
    time, since column order is fixed for the whole figure, so
    repeating them above every row was pure redundancy that also ate
    into the vertical space available between rows (see ROW_GAP_PX).
    One label row at the very top is enough.

    vertical_spacing is computed from ROW_GAP_PX rather than hardcoded
    as a fraction, so the actual pixel gap between rows stays constant
    regardless of n_rows (a fixed fraction gives a much smaller ABSOLUTE
    gap on a short figure than a tall one -- see ROW_GAP_PX's own
    comment).'''
    specs = [[{}] * n_cols for _ in range(n_rows)]
    subplot_titles = labels + [''] * (n_cols * (n_rows - 1))
    row_heights = [GROUP_ROW_HEIGHT] * n_rows
    total_height = sum(row_heights) + MARGIN_T + MARGIN_B
    vertical_spacing = min(0.3, ROW_GAP_PX / total_height) if n_rows > 1 else 0.05
    fig = make_subplots(rows=n_rows, cols=n_cols, specs=specs, row_heights=row_heights,
                        subplot_titles=subplot_titles, vertical_spacing=vertical_spacing,
                        horizontal_spacing=0.04)
    fig.update_annotations(font_size=SUBPLOT_TITLE_SIZE)
    return fig, row_heights


def _finish_grid_figure(fig, row_heights, n_cols):
    fig.update_layout(height=sum(row_heights) + MARGIN_T + MARGIN_B, width=max(1200, 260 * n_cols),
                      margin=dict(t=MARGIN_T, b=MARGIN_B), showlegend=False)
    return fig


def _ratio_vs_flux_figure(row_table, setup, snr_min):
    '''One row per doublet, one column per method: x = total flux of its
    two lines (log), y = the ratio, one point per SNR-passing row,
    pooled across the whole run (not restricted to repeat groups).
    Returns (fig, x_ranges) -- x_ranges maps each metric's proxy_names
    to the (lo, hi) linear-flux range computed here (the full-sample
    percentile, by far the most stable of the three panels since it
    pools thousands of rows rather than a few dozen repeat-group
    points), so _ratio_mad_figure/_flux_mad_figure can force the exact
    same x-window for a given line instead of recomputing their own
    from a much smaller, less representative sample (see
    _shared_row_range's forced_xr docstring, 260911).'''
    labels, routines, colors = setup['labels'], setup['routines'], setup['colors']
    ratio_metrics, n_cols = setup['ratio_metrics'], setup['n_cols']
    fig, row_heights = _grid_figure(len(ratio_metrics), n_cols, labels)
    x_ranges = {}

    for i, m in enumerate(ratio_metrics):
        row_i = i + 1
        all_x, all_y = [], []
        for col_j, (label, (r, v)) in enumerate(zip(labels, routines), start=1):
            sel = (row_table['ROUTINE'] == r) & (row_table['VARIANT'] == v)
            flux_proxy, ratio = _metric_arrays(row_table[sel], m, snr_min)
            fig.add_trace(go.Scatter(x=flux_proxy, y=ratio, mode='markers',
                                     marker=dict(color=colors[label], size=5, opacity=0.6),
                                     showlegend=False),
                          row=row_i, col=col_j)
            # exponentformat='e': flux values here are ~1e-14 to 1e-17 --
            # Plotly's default SI-prefix abbreviation renders these as
            # "10f"/"100p" (femto/pico) instead of scientific notation,
            # a recurring problem in this project -- see the
            # feedback_plotly_log_axes memory (Rule 5): apply this on
            # sight to any flux-valued axis, do not wait to be told.
            # tickangle=45, fixed rather than Plotly's default 'auto':
            # auto picks a rotation per subplot based on that
            # subplot's own tick-label lengths/spacing, which differs
            # row to row (different flux magnitudes -> different
            # exponents/tick counts) even at identical column width --
            # left unset, the same axis style read horizontal in one
            # row, 45 degrees in another, and vertical in a third
            # (found by direct inspection this session).
            fig.update_xaxes(type='log', exponentformat='e', tickangle=45, row=row_i, col=col_j)
            if flux_proxy.size:
                all_x.append(flux_proxy)
                all_y.append(ratio)
        all_x = np.concatenate(all_x) if all_x else np.array([])
        all_y = np.concatenate(all_y) if all_y else np.array([])
        xr = _shared_row_range(fig, n_cols, row_i, all_x, all_y)
        x_ranges[m['proxy_names']] = xr
        if m['ref_value'] is not None:
            for col_j in range(1, n_cols + 1):
                _hline(fig, row_i, col_j, xr, m['ref_value'])
        fig.update_yaxes(title_text=_row_label(m), row=row_i, col=1)
        if i == len(ratio_metrics) - 1:
            for col_j in range(1, n_cols + 1):
                fig.update_xaxes(title_text='Total line flux', row=row_i, col=col_j)

    return _finish_grid_figure(fig, row_heights, n_cols), x_ranges


def _ratio_mad_figure(group_rows, setup, snr_min, x_ranges):
    '''One row per doublet, one column per method: x = a repeat group's
    median total line flux (log), y = that group's ratio MAD, one
    point per repeat group, CLOSE groups only (same restriction
    SkySubNebEval.print_comparison() uses). x_ranges (from
    _ratio_vs_flux_figure) supplies the SAME x-window this metric used
    there, forced via _shared_row_range's forced_xr rather than
    recomputed from this panel's own much smaller repeat-group sample
    (260911 -- see that function's docstring). Returns (fig, mad_rows)
    -- mad_rows feeds the summary table so it's built from exactly the
    points plotted here, not recomputed a second way.'''
    labels, routines, colors = setup['labels'], setup['routines'], setup['colors']
    ratio_metrics, tileids = setup['ratio_metrics'], setup['tileids']
    close_lookup, n_cols = setup['close_lookup'], setup['n_cols']
    fig, row_heights = _grid_figure(len(ratio_metrics), n_cols, labels)
    mad_rows = {label: {} for label in labels}

    for i, m in enumerate(ratio_metrics):
        row_i = i + 1
        all_x, all_y = [], []
        for col_j, (label, (r, v)) in enumerate(zip(labels, routines), start=1):
            pts_x, pts_y = [], []
            for tileid in tileids:
                sub = group_rows.get((r, v, tileid))
                if sub is None or not close_lookup.get((r, v, tileid), False):
                    continue
                flux_proxy, value = _metric_arrays(sub, m, snr_min)
                if flux_proxy.size < 2:
                    continue
                med = float(np.median(value))
                mad = float(1.4826 * np.median(np.abs(value - med)))
                pts_x.append(float(np.median(flux_proxy)))
                pts_y.append(mad)
            fig.add_trace(go.Scatter(x=pts_x, y=pts_y, mode='markers',
                                     marker=dict(color=colors[label], size=8, opacity=0.8),
                                     showlegend=False),
                          row=row_i, col=col_j)
            fig.update_xaxes(type='log', exponentformat='e', tickangle=45, row=row_i, col=col_j)
            mad_rows[label][m['label']] = pts_y
            if pts_x:
                all_x.extend(pts_x)
                all_y.extend(pts_y)
        all_x, all_y = np.asarray(all_x), np.asarray(all_y)
        _shared_row_range(fig, n_cols, row_i, all_x, all_y, forced_xr=x_ranges.get(m['proxy_names']))
        fig.update_yaxes(title_text=_row_label(m), row=row_i, col=1)
        if i == len(ratio_metrics) - 1:
            for col_j in range(1, n_cols + 1):
                fig.update_xaxes(title_text='Median flux', row=row_i, col=col_j)

    return _finish_grid_figure(fig, row_heights, n_cols), mad_rows


def _flux_mad_figure(group_rows, setup, snr_min, x_ranges):
    '''One row per FLUX_METRICS entry, one column per method: x = a
    repeat group's median total-doublet flux (log, m['proxy_names'] --
    same brightness scale _ratio_vs_flux_figure/_ratio_mad_figure use
    for this line, NOT necessarily the same members as the y-value
    below -- see _metric_arrays; x_ranges forces the exact same
    x-window those two panels used for this line, 260911, rather than
    recomputing one from this panel's own small repeat-group sample),
    y = that group's fractional MAD
    (1.4826*MAD/median) of the line's own plotted flux total
    (m['flux_names'], still just the dominant member for NII/OIII/SIII)
    across the group's repeat exposures -- same one-point-per-group
    design as _ratio_mad_figure, but fractional rather than absolute
    (flux spans many decades across lines/methods, so only the
    relative error is comparable across panels; ratios are already
    O(1) so absolute MAD is fine there). Replaced an earlier design
    (260910) that plotted every individual exposure's raw flux against
    its group's median -- see module docstring's 260911 History entry.
    Returns (fig, mad_rows) -- mad_rows are the exact fractional-MAD
    values plotted, reused directly by the summary table.'''
    labels, routines, colors = setup['labels'], setup['routines'], setup['colors']
    flux_metrics, tileids = setup['flux_metrics'], setup['tileids']
    close_lookup, n_cols = setup['close_lookup'], setup['n_cols']
    fig, row_heights = _grid_figure(len(flux_metrics), n_cols, labels)
    mad_rows = {label: {} for label in labels}

    for i, m in enumerate(flux_metrics):
        row_i = i + 1
        all_x, all_y = [], []
        for col_j, (label, (r, v)) in enumerate(zip(labels, routines), start=1):
            pts_x, pts_y = [], []
            for tileid in tileids:
                sub = group_rows.get((r, v, tileid))
                if sub is None or not close_lookup.get((r, v, tileid), False):
                    continue
                flux_proxy, value = _metric_arrays(sub, m, snr_min)
                if flux_proxy.size < 2:
                    continue
                med = float(np.median(value))
                if not med:
                    continue
                mad = float(1.4826 * np.median(np.abs(value - med)))
                pts_x.append(float(np.median(flux_proxy)))
                pts_y.append(mad / med)
            fig.add_trace(go.Scatter(x=pts_x, y=pts_y, mode='markers',
                                     marker=dict(color=colors[label], size=8, opacity=0.8),
                                     showlegend=False),
                          row=row_i, col=col_j)
            fig.update_xaxes(type='log', exponentformat='e', tickangle=45, row=row_i, col=col_j)
            mad_rows[label][m['label']] = pts_y
            if pts_x:
                all_x.extend(pts_x)
                all_y.extend(pts_y)
        all_x, all_y = np.asarray(all_x), np.asarray(all_y)
        _shared_row_range(fig, n_cols, row_i, all_x, all_y, forced_xr=x_ranges.get(m['proxy_names']))
        fig.update_yaxes(title_text=_row_label(m) + ' frac MAD', row=row_i, col=1)
        if i == len(flux_metrics) - 1:
            for col_j in range(1, n_cols + 1):
                fig.update_xaxes(title_text='Median flux', row=row_i, col=col_j)

    return _finish_grid_figure(fig, row_heights, n_cols), mad_rows


def _mad_summary_figure(setup, mad_rows):
    '''Interactive summary table: median-across-groups absolute MAD for
    ratios (same numbers/units as SkySubNebEval.print_comparison()),
    median-across-groups fractional MAD (MAD/median) for flux totals --
    fractional because flux is ~1e-14 while ratio MAD is O(1), so
    sharing one '.4f' format would silently round every flux entry to
    0.0000 (found via direct inspection).'''
    labels, metrics = setup['labels'], setup['metrics']
    header = ['Method'] + [f'{_row_label(m)} (frac MAD)' if m['kind'] == 'flux' else f'{_row_label(m)} MAD'
                           for m in metrics]
    cols = [labels]
    for m in metrics:
        col = []
        for label in labels:
            vals = mad_rows[label].get(m['label'], [])
            col.append(f'{np.median(vals):.4f}' if vals else '--')
        cols.append(col)
    return _table_figure(header, cols, width=max(1000, 150 * len(header)))


def _continuum_pre_x_by_group(cont_group_rows, setup):
    '''One shared pre-subtraction brightness value per (tileid, arm),
    used as the x-axis in BOTH continuum repeat-consistency panels so a
    given repeat group sits at the identical x position in every
    column -- PRE_MED_<arm> is method-invariant by construction (see
    build_continuum_table: all 5 methods reconstruct the same input
    CFrame via FLUX+SKY), but each method computed its own copy, so
    this takes ONE (from the first routine, same "representative
    method" convention build_figures already uses) rather than 5
    near-identical copies that could drift a pixel apart from
    floating-point noise between methods.'''
    rep_r, rep_v = setup['routines'][0]
    out = {}
    for (r, v, tileid), sub in cont_group_rows.items():
        if r != rep_r or v != rep_v:
            continue
        for arm in ARM_LABELS:
            vals = np.asarray(sub[f'PRE_MED_{arm}'], dtype=float)
            finite = vals[np.isfinite(vals)]
            if finite.size:
                out[(tileid, arm)] = float(np.median(finite))
    return out


def _continuum_repeat_figure(cont_group_rows, cont_tileids, cont_close_lookup, setup,
                             stat_col, pre_x, normalize_by='own', log_y=False):
    '''One row per spectrograph arm (B/R/Z), one column per method: x =
    the repeat group's PRE-subtraction continuum median in that arm
    (log, shared across every column -- see _continuum_pre_x_by_group),
    y = that group's fractional MAD of stat_col (POST_MED or
    POST_NMAD) for that arm across the group's repeat exposures. CLOSE
    groups only, same restriction the nebular-line panels use. Used
    for BOTH the post-subtraction LEVEL panel and the post-subtraction
    RMS panel -- identical structure, only the column name/normalizer
    differs, so one function serves both rather than two near-duplicate
    ~40-line bodies (unlike _ratio_mad_figure/_flux_mad_figure, which
    differ in their actual gating/value logic, not just a column name).

    The two extra parameters::

        normalize_by : 'own' or 'pre'
            'own' (RMS panel): divide MAD by abs(median(stat_col)) --
            safe here, since POST_NMAD (a noise floor) is never near
            zero. 'pre' (LEVEL panel): divide MAD by the group's PRE-
            subtraction brightness (the same value plotted on x)
            instead -- POST_MED (the level) is EXPECTED to be near
            zero for a good method, so dividing by its own magnitude
            is unstable/inverted (a near-perfect group can score
            worse than a badly-biased one purely from a near-zero
            denominator; found and fixed 260911 by direct inspection
            of real per-group numbers -- see this session's own
            investigation). Dividing by the group's pre-subtraction
            brightness instead is always well-behaved (that
            denominator is never near zero) and physically
            meaningful: leftover scatter as a fraction of how bright
            the sky was to begin with.
        log_y : bool
            True for the LEVEL panel: once normalized by 'pre'
            instead of 'own', the fractional values span a genuine
            ~2 decades across groups -- worth a log axis, same
            reasoning as every flux axis elsewhere in this module
            (see _shared_row_range's log_y). False for the RMS panel
            (still normalized 'own', unaffected by this change).

    Returns (fig, mad_rows) -- mad_rows feeds the summary table, same
    pattern as the nebular-line panels.'''
    labels, routines, colors = setup['labels'], setup['routines'], setup['colors']
    n_cols = setup['n_cols']
    fig, row_heights = _grid_figure(len(ARM_LABELS), n_cols, labels)
    mad_rows = {label: {} for label in labels}

    for i, arm in enumerate(ARM_LABELS):
        row_i = i + 1
        all_x, all_y = [], []
        for col_j, (label, (r, v)) in enumerate(zip(labels, routines), start=1):
            pts_x, pts_y = [], []
            for tileid in cont_tileids:
                sub = cont_group_rows.get((r, v, tileid))
                if sub is None or not cont_close_lookup.get((r, v, tileid), False):
                    continue
                x = pre_x.get((tileid, arm))
                vals = np.asarray(sub[f'{stat_col}_{arm}'], dtype=float)
                finite = vals[np.isfinite(vals)]
                if x is None or x <= 0 or finite.size < 2:
                    continue
                med = float(np.median(finite))
                mad = float(1.4826 * np.median(np.abs(finite - med)))
                denom = x if normalize_by == 'pre' else abs(med)
                if not denom:
                    continue
                pts_x.append(x)
                pts_y.append(mad / denom)
            fig.add_trace(go.Scatter(x=pts_x, y=pts_y, mode='markers',
                                     marker=dict(color=colors[label], size=8, opacity=0.8),
                                     showlegend=False),
                          row=row_i, col=col_j)
            fig.update_xaxes(type='log', exponentformat='e', tickangle=45, row=row_i, col=col_j)
            mad_rows[label][arm] = pts_y
            if pts_x:
                all_x.extend(pts_x)
                all_y.extend(pts_y)
        all_x, all_y = np.asarray(all_x), np.asarray(all_y)
        _shared_row_range(fig, n_cols, row_i, all_x, all_y, log_y=log_y)
        fig.update_yaxes(title_text=f'{arm} frac MAD', row=row_i, col=1)
        if i == len(ARM_LABELS) - 1:
            for col_j in range(1, n_cols + 1):
                fig.update_xaxes(title_text='Pre-sub continuum', row=row_i, col=col_j)

    return _finish_grid_figure(fig, row_heights, n_cols), mad_rows


def _continuum_level_all_figure(cont_table, setup):
    '''One row per spectrograph arm (B/R/Z), one column per method: x =
    that exposure's own PRE-subtraction continuum median (log), y =
    its POST-subtraction continuum median -- LINEAR, and deliberately
    NOT clipped to positive values the way every flux axis elsewhere
    in this module is, because the point of this panel is exactly
    whether y goes negative. A negative post-subtraction continuum is
    not physically possible for a correct subtraction (continuum flux
    cannot be negative; a negative value means too much sky was
    removed) -- a dashed red line at y=0 marks the boundary directly.

    Unlike every other panel in this section, this uses EVERY exposure
    in cont_table, not just repeat groups: over-subtraction is visible
    on a single exposure with no repeat needed at all, so restricting
    to the ~34 CLOSE groups would throw away the vast majority of the
    check for no reason (260911, added directly at the user's request
    -- this session's own investigation into the LEVEL-consistency
    panel surfaced the "should never be negative" question along the
    way, but it deserves its own full-sample check, not a repeat-only
    one).

    Returns (fig, frac_neg) -- frac_neg[label][arm] is the fraction of
    that method's exposures with a negative post-subtraction median in
    that arm, feeding _continuum_negative_summary_figure.'''
    labels, routines, colors = setup['labels'], setup['routines'], setup['colors']
    n_cols = setup['n_cols']
    fig, row_heights = _grid_figure(len(ARM_LABELS), n_cols, labels)
    frac_neg = {label: {} for label in labels}

    for i, arm in enumerate(ARM_LABELS):
        row_i = i + 1
        all_x, all_y = [], []
        for col_j, (label, (r, v)) in enumerate(zip(labels, routines), start=1):
            sel = (cont_table['ROUTINE'] == r) & (cont_table['VARIANT'] == v)
            x = np.asarray(cont_table[f'PRE_MED_{arm}'], dtype=float)[sel]
            y = np.asarray(cont_table[f'POST_MED_{arm}'], dtype=float)[sel]
            ok = np.isfinite(x) & (x > 0) & np.isfinite(y)
            x, y = x[ok], y[ok]
            fig.add_trace(go.Scatter(x=x, y=y, mode='markers',
                                     marker=dict(color=colors[label], size=4, opacity=0.35),
                                     showlegend=False),
                          row=row_i, col=col_j)
            fig.update_xaxes(type='log', exponentformat='e', tickangle=45, row=row_i, col=col_j)
            fig.update_yaxes(exponentformat='e', row=row_i, col=col_j)
            frac_neg[label][arm] = float(np.mean(y < 0)) if y.size else np.nan
            if x.size:
                all_x.append(x)
                all_y.append(y)
        all_x = np.concatenate(all_x) if all_x else np.array([])
        all_y = np.concatenate(all_y) if all_y else np.array([])
        xr = _shared_row_range(fig, n_cols, row_i, all_x, all_y)
        for col_j in range(1, n_cols + 1):
            _hline(fig, row_i, col_j, xr, 0.0, color='crimson')
        fig.update_yaxes(title_text=f'{arm} post-sub level', row=row_i, col=1)
        if i == len(ARM_LABELS) - 1:
            for col_j in range(1, n_cols + 1):
                fig.update_xaxes(title_text='Pre-sub continuum', row=row_i, col=col_j)

    return _finish_grid_figure(fig, row_heights, n_cols), frac_neg


def _continuum_negative_summary_figure(setup, frac_neg):
    '''Interactive summary table: percentage of ALL exposures (not just
    repeat groups) with a negative post-subtraction continuum median,
    one column per arm -- the scalar readout for _continuum_level_all_figure.'''
    labels = setup['labels']
    header = ['Method'] + [f'{arm} % negative' for arm in ARM_LABELS]
    cols = [labels]
    for arm in ARM_LABELS:
        cols.append([f'{100 * frac_neg[l][arm]:.1f}%' if arm in frac_neg[l] and np.isfinite(frac_neg[l][arm])
                    else '--' for l in labels])
    return _table_figure(header, cols, width=max(900, 150 * len(header)))


def _continuum_summary_figure(setup, level_mad_rows, rms_mad_rows):
    '''Interactive summary table: median-across-groups fractional MAD
    for the post-subtraction continuum LEVEL and, separately, for its
    RMS/NMAD -- one column per arm per quantity, same '.4f' fractional-
    MAD convention _mad_summary_figure uses for flux totals.'''
    labels = setup['labels']
    header = (['Method'] + [f'{arm} level (frac MAD)' for arm in ARM_LABELS]
             + [f'{arm} RMS (frac MAD)' for arm in ARM_LABELS])
    cols = [labels]
    for arm in ARM_LABELS:
        cols.append([f'{np.median(level_mad_rows[l][arm]):.4f}' if level_mad_rows[l].get(arm) else '--'
                    for l in labels])
    for arm in ARM_LABELS:
        cols.append([f'{np.median(rms_mad_rows[l][arm]):.4f}' if rms_mad_rows[l].get(arm) else '--'
                    for l in labels])
    return _table_figure(header, cols, width=max(1000, 150 * len(header)))


def build_continuum_figures(cont_table, cont_meta_rows, cont_group_rows, setup):
    '''
    Build the continuum-residual section's figures.

    Parameters
    ----------
    cont_table : astropy.table.Table
        Full per-exposure table from build_continuum_table -- used
        directly by the all-exposures negative-continuum check
        (_continuum_level_all_figure), which needs every row, not just
        the repeat groups the other two panels are restricted to.
    cont_meta_rows, cont_group_rows : see build_continuum_run_tables.
    setup : dict
        From _common_setup(row_table, ...) on the NEBULAR-line
        row_table -- reused here (routines/labels/colors/n_cols only)
        rather than recomputed, since main() always builds both tables
        from the identical fits_files list in the identical order, so
        the (ROUTINE, VARIANT) columns necessarily match.

    Returns
    -------
    dict with keys 'level_all', 'neg_summary' (always present -- these
    need no repeat groups at all), 'level', 'rms', 'summary' (None for
    all three if there are no CLOSE groups), plus 'n_close_groups'/
    'n_close_exposures' (0 if none) -- same counts/convention as
    build_figures' nebular-line ones, reported in write_report's
    caption for this section.
    '''
    level_all_fig, frac_neg = _continuum_level_all_figure(cont_table, setup)
    figs = {
        'level_all': level_all_fig,
        'neg_summary': _continuum_negative_summary_figure(setup, frac_neg),
        'level': None, 'rms': None, 'summary': None,
        'n_close_groups': 0, 'n_close_exposures': 0,
    }

    tileids = sorted(set(tileid for (_r, _v, tileid) in cont_group_rows))
    if not tileids:
        return figs
    close_lookup = {(m['ROUTINE'], m['VARIANT'], m['TILEID']): bool(m['CLOSE']) for m in cont_meta_rows}

    pre_x = _continuum_pre_x_by_group(cont_group_rows, setup)
    level_fig, level_mad_rows = _continuum_repeat_figure(
        cont_group_rows, tileids, close_lookup, setup, 'POST_MED', pre_x,
        normalize_by='pre', log_y=True)
    rms_fig, rms_mad_rows = _continuum_repeat_figure(
        cont_group_rows, tileids, close_lookup, setup, 'POST_NMAD', pre_x,
        normalize_by='own', log_y=False)
    figs['level'] = level_fig
    figs['rms'] = rms_fig
    figs['summary'] = _continuum_summary_figure(setup, level_mad_rows, rms_mad_rows)

    rep_r, rep_v = setup['routines'][0]
    close_tileids = [t for t in tileids if close_lookup.get((rep_r, rep_v, t), False)]
    figs['n_close_groups'] = len(close_tileids)
    figs['n_close_exposures'] = sum(len(cont_group_rows[(rep_r, rep_v, t)]) for t in close_tileids)
    return figs


def build_figures(row_table, repeats_table, group_rows, snr_min=5.0):
    '''
    Build every figure the HTML report needs.

    Parameters
    ----------
    row_table, repeats_table, group_rows : see build_run_tables.
    snr_min : float
        Must match whatever build_run_tables/fit_file used -- needed
        here again only for FLUX_METRICS' own SNR gate (DOUBLETS ratios
        are already gated in row_table/group_rows).

    Returns
    -------
    dict with keys 'summary', 'ratio_vs_flux' (always present),
    'ratio_mad', 'flux_mad', 'mad_summary' (None if there were no
    repeat groups to plot), 'n_cols' (method-column count, for
    write_report's captions), and 'n_close_groups'/'n_close_exposures'
    (0 if none) -- the counts write_report reports just below the
    Repeat-Group Scatter heading, so a reader isn't left guessing how
    much data that panel (and everything after it) is actually built
    from, given the first section covers something like 20x more rows
    (see this session's own investigation into why e.g. SIII has so
    few points there). All 5 methods share the same underlying repeat-
    tileid groups/exposures (they're re-fitting the same input
    observations), so one representative (routine, variant) is enough
    -- taking setup['routines'][0] rather than summing across methods.
    '''
    setup = _common_setup(row_table, repeats_table, group_rows)
    ratio_vs_flux_fig, x_ranges = _ratio_vs_flux_figure(row_table, setup, snr_min)
    figs = {
        'summary': _run_summary_figure(row_table, setup),
        'ratio_vs_flux': ratio_vs_flux_fig,
        'ratio_mad': None, 'flux_mad': None, 'mad_summary': None,
        'n_close_groups': 0, 'n_close_exposures': 0,
    }
    if setup['have_groups']:
        ratio_mad_fig, ratio_mad_rows = _ratio_mad_figure(group_rows, setup, snr_min, x_ranges)
        flux_fig, flux_mad_rows = _flux_mad_figure(group_rows, setup, snr_min, x_ranges)
        combined_mad_rows = {label: {**ratio_mad_rows[label], **flux_mad_rows[label]}
                             for label in setup['labels']}
        figs['ratio_mad'] = ratio_mad_fig
        figs['flux_mad'] = flux_fig
        figs['mad_summary'] = _mad_summary_figure(setup, combined_mad_rows)

        rep_r, rep_v = setup['routines'][0]
        close_tileids = [t for t in setup['tileids']
                         if setup['close_lookup'].get((rep_r, rep_v, t), False)]
        figs['n_close_groups'] = len(close_tileids)
        figs['n_close_exposures'] = sum(len(group_rows[(rep_r, rep_v, t)]) for t in close_tileids)
    figs['n_cols'] = setup['n_cols']
    return figs


_REPORT_CSS = '''
body { font-family: Helvetica, Arial, sans-serif; margin: 20px 40px 100px; color: #222; }
h1 { margin: 0 0 30px; }
h2 { margin: 70px 0 12px; padding-top: 24px; border-top: 3px solid #444; }
h3 { margin: 40px 0 8px; color: #555; }
h3.subsection-break { margin-top: 70px; padding-top: 24px; border-top: 3px solid #444; }
.section { margin-bottom: 50px; }
p.caption { margin: 0 0 16px; max-width: 900px; color: #555; line-height: 1.4; }
'''


def write_report(figs, title, outfile):
    '''
    Assemble build_figures()'s output into one real HTML document --
    actual <h1>/<h2>/<h3> headings with CSS margins reserving space
    between sections, not Plotly annotations/margins standing in for
    them inside a single canvas (see module docstring's History for why
    that kept causing overlap bugs). Only the first embedded figure
    loads the Plotly JS library (via CDN); the rest reuse it. A short
    <p class="caption"> precedes every section/subsection (added
    260911) stating what is plotted, so the report is readable without
    this docstring open alongside it.
    '''
    parts = ['<!DOCTYPE html><html><head><meta charset="utf-8">',
            f'<title>{title}</title>', f'<style>{_REPORT_CSS}</style>',
            '</head><body>', f'<h1>{title}</h1>']

    def _add(fig, first=False):
        parts.append('<div class="section">')
        parts.append(fig.to_html(full_html=False, include_plotlyjs='cdn' if first else False))
        parts.append('</div>')

    def _caption(text):
        parts.append(f'<p class="caption">{text}</p>')

    _caption('Rows successfully fit per method, and how many of those rows passed the '
             'SNR gate for each doublet ratio (N and % of that method\'s own row count).')
    _add(figs['summary'], first=True)

    parts.append('<h2>Line Ratios</h2>')
    parts.append('<h3>Ratio vs. Line Flux</h3>')
    _caption('One point per spectrum, pooled across the whole run: x is the total flux of the '
             'doublet\'s two lines (log scale), y is the measured ratio. The dashed line, '
             'where shown, is the reference value -- a fixed atomic-physics ratio, a '
             'reddening ceiling, or a low-density-limit value (see SkySubNebEval.DOUBLETS). '
             'A method with a real ratio bias sits systematically off that line; scatter '
             'should shrink toward higher flux regardless.')
    _add(figs['ratio_vs_flux'])

    if figs['ratio_mad'] is not None:
        parts.append('<h3 class="subsection-break">Repeat-Group Scatter (MAD) vs. Median Flux</h3>')
        _caption('One point per repeat group -- closely-spaced repeat observations of the '
                 'same tile: x is the group\'s median total line flux, y is the robust scatter '
                 '(MAD) of the ratio across that group\'s exposures. This is a method\'s own '
                 'reproducibility (lower is better), independent of whether it is biased '
                 'relative to the reference line above.')
        _caption(f'From here through the rest of the report, every panel draws on only '
                 f'{figs["n_close_groups"]} distinct repeat-tileid groups with closely-spaced '
                 f'repeats, comprising {figs["n_close_exposures"]} individual exposures in '
                 f'total (the same {figs["n_close_groups"]} groups are evaluated separately '
                 f'for each of the {figs["n_cols"]} methods above) -- a much smaller pool than '
                 f'the full sample in the section above, which is why a line with a lower '
                 f'per-row detection rate can show very few points from here on.')
        _add(figs['ratio_mad'])

        parts.append('<h2>Repeat-Observation Flux Consistency</h2>')
        parts.append('<h3>Repeat-Group Fractional Scatter vs. Median Flux</h3>')
        _caption('One point per repeat group: x is the group\'s median total line flux, y is '
                 'the fractional MAD (robust MAD / median) of that same flux across the '
                 'group\'s exposures. Fractional rather than absolute because flux spans many '
                 'orders of magnitude across lines and methods, so only the relative error is '
                 'comparable panel to panel. Lower is better.')
        _add(figs['flux_mad'])

        parts.append('<h3>Summary</h3>')
        _caption('Median, across all repeat groups, of the scatter plotted in the two '
                 'sections above: MAD for line ratios (same units as the ratio itself), '
                 'fractional MAD for flux totals (dimensionless, comparable across lines).')
        _add(figs['mad_summary'])

    if figs.get('cont_level_all') is not None:
        parts.append('<h2>Continuum Residual (B/R/Z)</h2>')
        _caption('Independent of the emission-line sections above: for each spectrograph arm '
                 '(B/R/Z), GetSkyCont.arm_continuum_stats is run directly on raw flux, restricted '
                 'to the sky-line-free pixels of data/sky_mask.fits -- no local continuum fit -- '
                 'giving a median (brightness) and NMAD (scatter) per arm. The SAME function is '
                 'used both on the pre-subtraction CFrame flux (reconstructed as FLUX+SKY, '
                 'identical across all 5 methods) and on each method\'s own post-subtraction '
                 'FLUX, so the two are directly comparable.')

        parts.append('<h3>Post-Subtraction Continuum Level (All Exposures)</h3>')
        _caption('One point per EXPOSURE (not just repeat groups -- every row in the run), per '
                 'arm: x is that exposure\'s own pre-subtraction continuum median (log scale), '
                 'y is its post-subtraction continuum median, linear and unclipped. A physical '
                 'continuum flux cannot be negative -- a point below the dashed red line means '
                 'that exposure had more sky removed from it than was actually there, in that '
                 'arm. This needs no repeat observation at all, so it runs on the full sample.')
        _add(figs['cont_level_all'])
        _caption('Percentage of all exposures (not just repeat groups) with a negative '
                 'post-subtraction continuum median, per arm -- the scalar readout for the '
                 'plot above.')
        _add(figs['cont_neg_summary'])

        if figs.get('cont_level') is not None:
            parts.append('<h3 class="subsection-break">Post-Subtraction Continuum Level Consistency</h3>')
            _caption('One point per repeat group, per arm: x is the group\'s PRE-subtraction '
                     'continuum median in that arm (log scale, shared across every method -- the '
                     'original CFrame brightness, for context), y is the fractional MAD (robust '
                     'MAD / that group\'s pre-subtraction brightness) of the POST-subtraction '
                     'continuum median across the group\'s repeat exposures -- log scale. '
                     'Normalized against the PRE-subtraction brightness rather than the post-'
                     'subtraction level\'s own median: the leftover level is expected to sit '
                     'near zero for a good method, so dividing by its own magnitude is unstable '
                     '(a near-perfect, near-zero group can score worse than a badly-biased one '
                     'purely from a near-zero denominator -- found by direct inspection of real '
                     'per-group numbers). Tests whether a method leaves a stable, repeatable '
                     'continuum-level bias, as a fraction of the original sky brightness, rather '
                     'than one that jumps around from exposure to exposure. Lower is better.')
            _caption(f'From here through the rest of this section, every panel draws on only '
                     f'{figs["cont_n_close_groups"]} distinct repeat-tileid groups with closely-'
                     f'spaced repeats, comprising {figs["cont_n_close_exposures"]} individual '
                     f'exposures in total -- computed independently of the emission-line grouping '
                     f'above (same underlying logic, but a different, arm-based row table), so the '
                     f'count need not match {figs["n_close_groups"]}.')
            _add(figs['cont_level'])

            parts.append('<h3>Post-Subtraction Continuum RMS Consistency</h3>')
            _caption('Same x and same repeat groups as above; y is instead the fractional MAD of '
                     'the POST-subtraction continuum RMS/NMAD (the leftover scatter within clean '
                     'pixels, not its central level), normalized against its OWN median as before '
                     '(safe here -- a noise floor is never near zero the way a level can be) -- '
                     'tests whether the leftover noise level itself is repeatable, separately from '
                     'any systematic offset.')
            _add(figs['cont_rms'])

            parts.append('<h3>Summary</h3>')
            _caption('Median, across all repeat groups, of the fractional MAD plotted in the two '
                     'sections above -- level (relative to pre-subtraction brightness) and RMS '
                     '(relative to its own median), one column per arm.')
            _add(figs['cont_summary'])

    parts.append('</body></html>')
    Path(outfile).write_text('\n'.join(parts))


def main():
    parser = argparse.ArgumentParser(
        description=_usage_from_doc(__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fits_file', nargs='*', help='SkySub*.py output FITS file(s)')
    parser.add_argument('-lines_file', default=None,
                        help='Precomputed SkySubNebEval.py _lines.fits (skips refitting)')
    parser.add_argument('-v', dest='vel', type=float, default=None,
                        help='Nebular systemic velocity (km/s) override')
    parser.add_argument('-lmc', action='store_true', help=f'Use the LMC velocity (~{LMC_VEL:.0f} km/s)')
    parser.add_argument('-smc', action='store_true', help=f'Use the SMC velocity (~{SMC_VEL:.0f} km/s)')
    parser.add_argument('-sigma', dest='sigma_guess', type=float, default=1.0,
                        help='Initial Gaussian sigma guess (Angstrom)')
    parser.add_argument('-snr_min', type=float, default=5.0,
                        help='Minimum per-line SNR to report a doublet ratio')
    parser.add_argument('-mjd_close', type=float, default=7.0,
                        help='MJD span (days) below which a repeat group is "closely spaced"')
    parser.add_argument('-mask', default=None,
                        help='Clean-pixel mask FITS for the continuum-residual section '
                             '(default: data/sky_mask.fits)')
    parser.add_argument('-title', default='SkySubNebRun')
    parser.add_argument('-outfile', default='nebrun.html', help='Output HTML path')
    args = parser.parse_args()

    if bool(args.fits_file) == bool(args.lines_file):
        parser.error('give either fits_file(s) or -lines_file, not both/neither')

    row_table, repeats_table, group_rows = build_run_tables(
        fits_files=args.fits_file or None, lines_file=args.lines_file,
        vel=args.vel, lmc=args.lmc, smc=args.smc, sigma_guess=args.sigma_guess,
        snr_min=args.snr_min, mjd_close=args.mjd_close)

    _summarize(row_table, repeats_table, group_rows)

    figs = build_figures(row_table, repeats_table, group_rows, snr_min=args.snr_min)

    if args.fits_file:
        setup = _common_setup(row_table, repeats_table, group_rows)
        cont_table, cont_meta_rows, cont_group_rows = build_continuum_run_tables(
            args.fits_file, mjd_close=args.mjd_close, mask_file=args.mask)
        print(f'\ncontinuum_table: {len(cont_table)} rows '
             f'({len(cont_meta_rows)} repeat-tileid groups)')
        cont_figs = build_continuum_figures(cont_table, cont_meta_rows, cont_group_rows, setup)
        figs.update({f'cont_{k}': v for k, v in cont_figs.items()})
    else:
        print('\n-lines_file given: skipping the continuum-residual section (it needs the raw '
             'WAVE/FLUX/SKY arrays, not just the cached per-row line fits).')
        figs.update(cont_level_all=None, cont_neg_summary=None,
                    cont_level=None, cont_rms=None, cont_summary=None,
                    cont_n_close_groups=0, cont_n_close_exposures=0)

    write_report(figs, args.title, args.outfile)
    print(f'\nWrote {args.outfile}')


if __name__ == '__main__':
    main()
