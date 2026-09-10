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
            [one row per FLUX_METRICS entry, one column per method;
             x = a repeat group's median flux, y = EACH individual
             exposure's own measured flux in that group -- not a
             collapsed dispersion number]
            <h3>Summary</h3>
            [interactive MAD/fractional-STD table]

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

'''

import argparse
import re
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from SkySubNebEval import (fit_file, repeat_scatter, DOUBLETS, NEBULAR_LINES,
                           PLACEHOLDER_TILEID)
from sky_nebular_leak_eval import LMC_VEL, SMC_VEL

METHOD_COLORS = ['black', 'firebrick', 'seagreen', 'steelblue', 'darkorange', 'purple']
SUBPLOT_TITLE_SIZE = 14
TABLE_ROW_HEIGHT = 300
GROUP_ROW_HEIGHT = 280
MARGIN_T = 90
MARGIN_B = 60

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


def _shared_row_range(fig, n_cols, row_i, all_x, all_y):
    '''Pool x/y across every method-column in one metric row and apply
    one shared range to all of them (1st/99th percentile, not literal
    min/max -- a single noisy method's outlier point would otherwise set
    a shared range wide enough to flatten every other column's real
    signal, the same lesson PlotSkySubNebEval.py's per-line-group rows
    already learned). x is always log-flux; y is linear throughout this
    module. Returns xr (linear-space tuple) for _hline.'''
    xr = tuple(np.percentile(all_x, [1, 99])) if all_x.size else None
    if xr and xr[0] >= xr[1] and xr[0] > 0:
        xr = (xr[0] / 3.0, xr[1] * 3.0)  # degenerate (identical/single point)
    for col_j in range(1, n_cols + 1):
        if xr and xr[0] > 0:
            fig.update_xaxes(range=[np.log10(xr[0]), np.log10(xr[1])], row=row_i, col=col_j)
        if all_y.size:
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

    Each entry: dict(label, kind, ratio_col, flux_names, ref_kind, ref_value).

    kind='ratio' (one per DOUBLETS entry): flux_names=(name_a, name_b)
    is the doublet's own two NEBULAR_LINES members, used as the x-axis
    brightness proxy (their mean flux); ratio_col is both the ratio-vs-
    flux y-value and the repeat-group MAD's dispersion source.

    kind='flux' (one per FLUX_METRICS entry): ratio_col is None --
    there is no "ratio vs its own flux" panel (that would just be a
    trivial 1:1 relation), so these only appear in the flux-cloud
    section, where flux_names is used for both the x-axis (group
    median) and the y-axis (each individual exposure's own value).
    '''
    metrics = [dict(label=d.upper(), kind='ratio', ratio_col=d.upper(),
                    flux_names=(name_a, name_b), ref_kind=kind, ref_value=value)
              for d, name_a, name_b, kind, value, _c in DOUBLETS]
    metrics += [dict(label=name, kind='flux', ratio_col=None, flux_names=names,
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

    flux_proxy is the x-axis brightness measure: the doublet's own two
    lines' mean flux for a 'ratio' metric, or the metric's own SNR-
    gated total (_total_flux) for a 'flux' metric. value is the ratio
    itself for a 'ratio' metric, or the same array as flux_proxy for a
    'flux' metric (both x and y come from the same per-row totals there
    -- x is their median, y is each one directly).
    '''
    if m['kind'] == 'ratio':
        name_a, name_b = m['flux_names']
        fa = np.asarray(sub[f'{name_a}_FLUX'], dtype=float)
        fb = np.asarray(sub[f'{name_b}_FLUX'], dtype=float)
        flux_proxy = (fa + fb) / 2.0
        value = np.asarray(sub[m['ratio_col']], dtype=float)
    else:
        flux_proxy = _total_flux(sub, m['flux_names'], snr_min)
        value = flux_proxy
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
    '''Shared make_subplots scaffold for the three scatter sections
    below -- one row per metric, one column per method (see module
    docstring for why: overlaying every method in one panel stops being
    legible past a handful of points; a single-method panel stays
    readable regardless of how large the run is).'''
    specs = [[{}] * n_cols for _ in range(n_rows)]
    subplot_titles = labels * n_rows
    row_heights = [GROUP_ROW_HEIGHT] * n_rows
    fig = make_subplots(rows=n_rows, cols=n_cols, specs=specs, row_heights=row_heights,
                        subplot_titles=subplot_titles, vertical_spacing=0.05,
                        horizontal_spacing=0.04)
    fig.update_annotations(font_size=SUBPLOT_TITLE_SIZE)
    return fig, row_heights


def _finish_grid_figure(fig, row_heights, n_cols):
    fig.update_layout(height=sum(row_heights) + MARGIN_T + MARGIN_B, width=max(1200, 260 * n_cols),
                      margin=dict(t=MARGIN_T, b=MARGIN_B), showlegend=False)
    return fig


def _ratio_vs_flux_figure(row_table, setup, snr_min):
    '''One row per doublet, one column per method: x = mean flux of its
    two lines (log), y = the ratio, one point per SNR-passing row,
    pooled across the whole run (not restricted to repeat groups).'''
    labels, routines, colors = setup['labels'], setup['routines'], setup['colors']
    ratio_metrics, n_cols = setup['ratio_metrics'], setup['n_cols']
    fig, row_heights = _grid_figure(len(ratio_metrics), n_cols, labels)

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
            fig.update_xaxes(type='log', exponentformat='e', row=row_i, col=col_j)
            if flux_proxy.size:
                all_x.append(flux_proxy)
                all_y.append(ratio)
        all_x = np.concatenate(all_x) if all_x else np.array([])
        all_y = np.concatenate(all_y) if all_y else np.array([])
        xr = _shared_row_range(fig, n_cols, row_i, all_x, all_y)
        if m['ref_value'] is not None:
            for col_j in range(1, n_cols + 1):
                _hline(fig, row_i, col_j, xr, m['ref_value'])
        fig.update_yaxes(title_text=_row_label(m), row=row_i, col=1)
        if i == len(ratio_metrics) - 1:
            for col_j in range(1, n_cols + 1):
                fig.update_xaxes(title_text='Mean line flux', row=row_i, col=col_j)

    return _finish_grid_figure(fig, row_heights, n_cols)


def _ratio_mad_figure(group_rows, setup, snr_min):
    '''One row per doublet, one column per method: x = a repeat group's
    median flux (log), y = that group's ratio MAD, one point per
    repeat group, CLOSE groups only (same restriction
    SkySubNebEval.print_comparison() uses). Returns (fig, mad_rows) --
    mad_rows feeds the summary table so it's built from exactly the
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
            fig.update_xaxes(type='log', exponentformat='e', row=row_i, col=col_j)
            mad_rows[label][m['label']] = pts_y
            if pts_x:
                all_x.extend(pts_x)
                all_y.extend(pts_y)
        all_x, all_y = np.asarray(all_x), np.asarray(all_y)
        _shared_row_range(fig, n_cols, row_i, all_x, all_y)
        fig.update_yaxes(title_text=_row_label(m), row=row_i, col=1)
        if i == len(ratio_metrics) - 1:
            for col_j in range(1, n_cols + 1):
                fig.update_xaxes(title_text='Median flux', row=row_i, col=col_j)

    return _finish_grid_figure(fig, row_heights, n_cols), mad_rows


def _flux_cloud_figure(group_rows, setup, snr_min):
    '''One row per FLUX_METRICS entry, one column per method: x = a
    repeat group's median flux (log), y = EACH individual exposure's
    own measured flux in that group (all sharing that group's x
    position) -- the actual measurements, not a collapsed dispersion
    number, per explicit request. Returns (fig, mad_rows) -- mad_rows
    here holds each group's fractional STD (std/median), used only by
    the summary table, not by this plot.'''
    labels, routines, colors = setup['labels'], setup['routines'], setup['colors']
    flux_metrics, tileids = setup['flux_metrics'], setup['tileids']
    close_lookup, n_cols = setup['close_lookup'], setup['n_cols']
    fig, row_heights = _grid_figure(len(flux_metrics), n_cols, labels)
    mad_rows = {label: {} for label in labels}

    for i, m in enumerate(flux_metrics):
        row_i = i + 1
        all_x, all_y = [], []
        for col_j, (label, (r, v)) in enumerate(zip(labels, routines), start=1):
            pts_x, pts_y, frac_std = [], [], []
            for tileid in tileids:
                sub = group_rows.get((r, v, tileid))
                if sub is None or not close_lookup.get((r, v, tileid), False):
                    continue
                flux_proxy, value = _metric_arrays(sub, m, snr_min)
                if flux_proxy.size < 2:
                    continue
                med = float(np.median(value))
                pts_x.extend([med] * value.size)
                pts_y.extend(value.tolist())
                if med:
                    frac_std.append(float(np.std(value, ddof=1)) / med)
            fig.add_trace(go.Scatter(x=pts_x, y=pts_y, mode='markers',
                                     marker=dict(color=colors[label], size=6, opacity=0.6),
                                     showlegend=False),
                          row=row_i, col=col_j)
            fig.update_xaxes(type='log', exponentformat='e', row=row_i, col=col_j)
            fig.update_yaxes(exponentformat='e', row=row_i, col=col_j)
            mad_rows[label][m['label']] = frac_std
            if pts_x:
                all_x.extend(pts_x)
                all_y.extend(pts_y)
        all_x, all_y = np.asarray(all_x), np.asarray(all_y)
        _shared_row_range(fig, n_cols, row_i, all_x, all_y)
        fig.update_yaxes(title_text=_row_label(m), row=row_i, col=1)
        if i == len(flux_metrics) - 1:
            for col_j in range(1, n_cols + 1):
                fig.update_xaxes(title_text='Median flux', row=row_i, col=col_j)

    return _finish_grid_figure(fig, row_heights, n_cols), mad_rows


def _mad_summary_figure(setup, mad_rows):
    '''Interactive summary table: median-across-groups absolute MAD for
    ratios (same numbers/units as SkySubNebEval.print_comparison()),
    median-across-groups fractional STD (std/median) for flux totals --
    fractional because flux is ~1e-14 while ratio MAD is O(1), so
    sharing one '.4f' format would silently round every flux entry to
    0.0000 (found via direct inspection).'''
    labels, metrics = setup['labels'], setup['metrics']
    header = ['Method'] + [f'{_row_label(m)} (frac STD)' if m['kind'] == 'flux' else f'{_row_label(m)} MAD'
                           for m in metrics]
    cols = [labels]
    for m in metrics:
        col = []
        for label in labels:
            vals = mad_rows[label].get(m['label'], [])
            col.append(f'{np.median(vals):.4f}' if vals else '--')
        cols.append(col)
    return _table_figure(header, cols, width=max(1000, 150 * len(header)))


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
    dict with keys 'summary', 'ratio_vs_flux' (always present) and
    'ratio_mad', 'flux_cloud', 'mad_summary' (None if there were no
    repeat groups to plot).
    '''
    setup = _common_setup(row_table, repeats_table, group_rows)
    figs = {
        'summary': _run_summary_figure(row_table, setup),
        'ratio_vs_flux': _ratio_vs_flux_figure(row_table, setup, snr_min),
        'ratio_mad': None, 'flux_cloud': None, 'mad_summary': None,
    }
    if setup['have_groups']:
        ratio_mad_fig, ratio_mad_rows = _ratio_mad_figure(group_rows, setup, snr_min)
        flux_fig, flux_mad_rows = _flux_cloud_figure(group_rows, setup, snr_min)
        combined_mad_rows = {label: {**ratio_mad_rows[label], **flux_mad_rows[label]}
                             for label in setup['labels']}
        figs['ratio_mad'] = ratio_mad_fig
        figs['flux_cloud'] = flux_fig
        figs['mad_summary'] = _mad_summary_figure(setup, combined_mad_rows)
    return figs


_REPORT_CSS = '''
body { font-family: Helvetica, Arial, sans-serif; margin: 20px 40px 100px; color: #222; }
h1 { margin: 0 0 30px; }
h2 { margin: 70px 0 12px; padding-top: 24px; border-top: 3px solid #444; }
h3 { margin: 40px 0 8px; color: #555; }
.section { margin-bottom: 50px; }
'''


def write_report(figs, title, outfile):
    '''
    Assemble build_figures()'s output into one real HTML document --
    actual <h1>/<h2>/<h3> headings with CSS margins reserving space
    between sections, not Plotly annotations/margins standing in for
    them inside a single canvas (see module docstring's History for why
    that kept causing overlap bugs). Only the first embedded figure
    loads the Plotly JS library (via CDN); the rest reuse it.
    '''
    parts = ['<!DOCTYPE html><html><head><meta charset="utf-8">',
            f'<title>{title}</title>', f'<style>{_REPORT_CSS}</style>',
            '</head><body>', f'<h1>{title}</h1>']

    def _add(fig, first=False):
        parts.append('<div class="section">')
        parts.append(fig.to_html(full_html=False, include_plotlyjs='cdn' if first else False))
        parts.append('</div>')

    _add(figs['summary'], first=True)
    parts.append('<h2>Line Ratios</h2>')
    parts.append('<h3>Ratio vs. Line Flux</h3>')
    _add(figs['ratio_vs_flux'])
    if figs['ratio_mad'] is not None:
        parts.append('<h3>Repeat-Group Scatter (MAD) vs. Median Flux</h3>')
        _add(figs['ratio_mad'])
        parts.append('<h2>Repeat-Observation Flux Consistency</h2>')
        _add(figs['flux_cloud'])
        parts.append('<h3>Summary</h3>')
        _add(figs['mad_summary'])
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
    write_report(figs, args.title, args.outfile)
    print(f'\nWrote {args.outfile}')


if __name__ == '__main__':
    main()
