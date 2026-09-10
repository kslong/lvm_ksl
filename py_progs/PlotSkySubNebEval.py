#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Interactive Plotly visualization of SkySubNebEval.py's per-line/
    doublet fits for ONE exposure across one or more SkySub*.py output
    files (methods): the actual sky-subtracted FLUX spectrum around each
    line group, with the fitted Gaussian(s) overlaid and SNR/ratio
    annotated -- the picture behind SkySubNebEval.py's numbers, one
    exposure at a time.

Command line usage (if any)::

        usage: PlotSkySubNebEval.py [-h] [-row N | -expnum N] [-v VEL]
                                    [-lmc] [-smc] [-sigma S] [-title TITLE]
                                    [-outfile PATH]
                                    fits_file [fits_file ...]

        where

        fits_file     one or more SkySub*.py output FITS files (WAVE, FLUX,
                     SKY, DRP_ALL), all for the SAME exposure -- e.g. the
                     same expnum run through -routine drp/orig/dev1/dev2/
                     dev3. Each file's primary header ROUTINE/METHOD (Title/
                     METHOD keywords) labels its column.

        -row N        row index to select from each file (default: 0 --
                     correct as-is for a single-row file, as produced by
                     extracting one exposure; for a multi-row file, use
                     -expnum instead).

        -expnum N     select by DRP_ALL 'expnum' instead of row index
                     (overrides -row).

        -v VEL        nebular systemic velocity (km/s), overriding the
                     default per-row DRP_ALL['Redshift'] auto-lookup (see
                     SkySubNebEval.py's fit_file docstring / resolve_velocity)
                     -- same convention as SkySubNebEval.py.

        -lmc / -smc   shortcuts for the LMC (~262 km/s) / SMC (~146 km/s)
                     systemic velocity -- same override semantics as -v.

        -sigma S      initial Gaussian sigma guess in Angstrom (default 1.0).

        -title TITLE  plot title (default: the exposure's expnum).

        -outfile PATH output HTML path (default:
                     Overview_Plot/nebeval_<expnum>.html).

Description::

        Row 1: a pointing/Moon/Sun geometry table -- Target (Science/SkyE/
        SkyW/Moon/Sun), RA, Dec., PA, Ang. dist. (from the science field),
        Alt., Illum. (%), Astrometry Src, Shadow Ht (km). Same column set as
        QualSFrame.py's create_overview() pointing table (see
        _POINTING_TABLE_HEADER); that version recomputes these from SFrame
        header keywords, this one just reads the equivalent DRP_ALL columns
        (sci_skye_sep, sci_moon_sep, moon_fli, etc.), already computed once
        by SummarizeCframe.py. Pulled from the first input file's DRP_ALL
        row, since all columns are the same exposure by construction.

        Row 2: the full sky-subtracted spectrum (log y-axis, fixed range
        SPECTRUM_YMIN-SPECTRUM_YMAX regardless of the data -- so every
        exposure's overview panel is directly comparable at a glance, all
        input files/methods overlaid in one panel) -- the big-picture view,
        not the zoomed line windows below. Each trace is median-filtered
        (SPECTRUM_SMOOTH_PIX) before plotting, not raw per-pixel flux, so
        inter-method differences aren't swamped by per-pixel noise on this
        log-y, full-range panel; a color key (method name/color, inside the
        plot's own domain, not the figure margin -- see below) identifies
        each trace, alongside a compact exposure-identification line (expnum,
        tileid, MJD, DRP_ALL Survey classification, the nebular velocity
        actually used for the fits below) -- RA/Dec/Moon geometry are NOT
        repeated here, see Row 1's fuller table instead.

        Row 3: a summary table, one row per input file (method), one column
        per SkySubNebEval.DOUBLETS entry -- the measured ratio (SNR-gated,
        -snr_min), with the true/bound value in the column header and a
        warning marker on any 'fixed' ratio deviating >=3 sigma or any
        'bounded_above' ratio exceeding its ceiling. Fixed pixel height
        (TABLE_ROW_HEIGHT), not a fraction of the figure -- Plotly tables
        clip/scroll silently rather than growing to fit, so a fraction that
        shrinks as more line groups are added will clip real rows out of
        view even though the underlying data is complete (hit this at 8
        line groups: only ~2 of 5 methods' rows were visible).

        Below that, one row of panels per line group, one column per input
        file, each row sharing one y-axis range across all its columns (so
        amplitude/shape are directly comparable method to method, not each
        silently auto-scaled to its own data):

            OII      oii_a/oii_b, joint double-Gaussian fit (fit_oii_doublet)
            Hb       hb
            OIII_a   oiii_a
            OIII_b   oiii_b
            NII+Ha   nii_a, ha, nii_b together (close enough to share one
                     window)
            SII      sii_a/sii_b
            SIII_a   siii_a
            SIII_b   siii_b

        Each panel: observed FLUX (black, raw/unsmoothed -- unlike Row 2),
        the fitted Gaussian(s) (colored, reconstructed from SkySubNebEval.
        fit_nebular_row's own fit results -- same numbers as that tool's
        tables, not recomputed independently), and a title annotated with
        each line's SNR (and, for the DOUBLETS pair in that group, the
        measured ratio and its status against the known/bounded/free truth --
        see SkySubNebEval.DOUBLETS). Each row's shared y-range is bounded by
        the 2nd/98th percentile of the observed FLUX pooled across all
        columns (not literal min/max) plus the fitted curves' true min/max --
        a single noisy method's occasional extreme pixel would otherwise set
        a shared range wide enough to flatten every other (well-behaved)
        column's fit into a near-flat line, which reads as "this method's fit
        is bad" when it's actually the shared axis being dominated by a
        different column's outlier pixels.

        OI is deliberately not shown, matching SkySubNebEval.py's DOUBLETS
        (dominated by sky-subtraction residual, not real nebular flux).

Primary routines:

    plot_row   build the full figure for one exposure across N files.

History::

    260910  ksl  Coding begun.

'''

import argparse
from pathlib import Path

import numpy as np
from scipy.ndimage import median_filter
from astropy.io import fits
from astropy.table import Table
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from sky_nebular_leak_eval import _gaussian, resolve_velocity, LMC_VEL, SMC_VEL
from SkySubNebEval import fit_nebular_row, _doublet_ratio, DOUBLETS, NEBULAR_LINES

TRACE_WIDTH = 2.0
FIT_WIDTH = 2.5
SUBPLOT_TITLE_SIZE = 14
ROW_HEIGHT = 300
SPECTRUM_ROW_HEIGHT = 380 * 4
# Fixed pixel height for the summary table row, not a fraction of the
# figure -- a fraction shrinks as more line-group rows are added, which
# clipped the table down to ~2 visible rows out of 6 (header + 5 methods)
# once GROUPS grew; Plotly tables clip/scroll silently rather than
# growing to fit, so this must stay an absolute size, not scale with
# n_line_rows. 220 still clipped one row (4/5 visible) -- table chrome
# (header padding, borders) needs more headroom than a bare
# n_rows*cell_height sum suggests; 300 is a more generous margin.
TABLE_ROW_HEIGHT = 300
# Same reasoning/value as TABLE_ROW_HEIGHT -- the pointing table (header +
# Science/SkyE/SkyW/Moon/Sun = 6 rows) is the same shape as the doublet
# summary table (header + 5 methods = 6 rows), so it needs the same
# headroom to avoid the identical clipping bug.
POINTING_TABLE_ROW_HEIGHT = 300
# Fixed axis range (not data-derived) -- user-specified so every exposure's
# overview panel is directly comparable at a glance.
SPECTRUM_YMIN = 1e-17
SPECTRUM_YMAX = 1e-12
# Median-filter window (pixels) applied to the overview panel only -- at
# native pixel sampling, per-pixel noise on the log-y overview dominates
# the plot and hides genuine inter-method continuum/residual differences.
# A median filter (not a mean/boxcar) is used so it doesn't get dragged
# around by single-pixel cosmic-ray-like spikes -- consistent with using
# percentile, not max(), for this plot's y-ceiling (see
# feedback_plotly_log_axes memory).
SPECTRUM_SMOOTH_PIX = 11
# Explicit top/bottom margins (rather than relying on Plotly defaults for
# b) so the figure's total height exactly equals sum(row_heights) +
# MARGIN_T + MARGIN_B -- otherwise the domain area is a bit shorter than
# row_heights implies, and the title/legend (positioned in paper
# fractions of the domain) land closer to the plot than intended. This
# was the root cause of the spectrum panel overwriting the title text:
# with only margin.t set, the actual domain top crept up under the title.
MARGIN_T = 190
MARGIN_B = 80

# (group_name, [line names], pad_angstrom) -- window is the union of the
# member lines' own Doppler-shifted windows, padded by pad_angstrom.
GROUPS = [
    ('OII',    ['oii_a', 'oii_b'],          3.0),
    ('Hb',     ['hb'],                       3.0),
    ('OIII_a', ['oiii_a'],                   3.0),
    ('OIII_b', ['oiii_b'],                   3.0),
    ('NII+Ha', ['nii_a', 'ha', 'nii_b'],     3.0),
    ('SII',    ['sii_a', 'sii_b'],           3.0),
    ('SIII_a', ['siii_a'],                   3.0),
    ('SIII_b', ['siii_b'],                   3.0),
]

_LINE_LOOKUP = {name: (center, wmin, wmax) for name, center, wmin, wmax in NEBULAR_LINES}
_DOUBLET_FOR_LINE = {}
for _dname, _la, _lb, _kind, _val, _conf in DOUBLETS:
    _DOUBLET_FOR_LINE[_la] = (_dname, _la, _lb, _kind, _val)
    _DOUBLET_FOR_LINE[_lb] = (_dname, _la, _lb, _kind, _val)


def _group_window(names, vel, pad):
    zz = 1.0 + vel / 3e5
    los, his = [], []
    for name in names:
        center, wmin, wmax = _LINE_LOOKUP[name]
        los.append(zz * wmin)
        his.append(zz * wmax)
    return min(los) - pad, max(his) + pad


def _snr(fit):
    if fit is None or fit['amp_err'] <= 0:
        return np.nan
    return fit['amp'] / fit['amp_err']


def _doublet_header(dname, kind, value):
    if kind == 'fixed' and value is not None:
        return f'{dname}<br>(={value:.3g})'
    if kind == 'bounded_above' and value is not None:
        return f'{dname}<br>(≤{value:.3g})'
    return f'{dname}<br>(free)'


def _exposure_info_text(drp, i, vel):
    '''One-line exposure identification string for the overview panel --
    expnum/tileid/MJD/survey/velocity, pulled from the selected row's own
    DRP_ALL entry. vel is the per-row systemic velocity actually used for
    the nebular-line fits (see resolve_velocity) -- shown explicitly
    since silently using the wrong one (e.g. 0 for an LMC/SMC field) is
    exactly the failure mode this is guarding against (see
    resolve_velocity's docstring, 260910). RA/Dec and Moon geometry are
    NOT repeated here -- see _pointing_table_rows' fuller table instead.'''
    def g(name, fmt='{}'):
        if name not in drp.colnames:
            return None
        val = drp[name][i]
        return fmt.format(val)

    parts = []
    expnum = g('expnum')
    if expnum is not None:
        parts.append(f'expnum {expnum}')
    tileid = g('tileid')
    if tileid is not None:
        parts.append(f'tileid {tileid}')
    mjd = g('mjd', '{:.2f}')
    if mjd is not None:
        parts.append(f'MJD {mjd}')
    survey = g('Survey')
    if survey is not None:
        parts.append(f'survey {str(survey).strip()}')
    parts.append(f'vel {vel:.0f} km/s')
    return ' &nbsp;|&nbsp; '.join(parts)


# Column set matches QualSFrame.py's create_overview() pointing table
# (Target, RA, Dec., PA, Ang. dist., Alt., Illum. (%), Astrometry Src,
# Shadow Ht (km)) -- same table shown "at the top" of that tool's plots,
# reused here for consistency. That version recomputes Ang. dist./Moon/
# Sun geometry from SFrame header keywords via distance()/
# get_moon_info_las_campanas(); here the same quantities are already
# columns on DRP_ALL (sci_skye_sep, sci_skyw_sep, sci_moon_sep, moon_fli,
# etc. -- computed once by SummarizeCframe.py), so this just reads them
# rather than recomputing.
_POINTING_TABLE_HEADER = ['Target', 'RA', 'Dec.', 'PA', 'Ang. dist.', 'Alt.',
                          'Illum. (%)', 'Astrometry Src', 'Shadow Ht (km)']


def _pointing_table_rows(drp, i):
    '''Science/SkyE/SkyW/Moon/Sun pointing-geometry rows for the selected
    DRP_ALL row -- see _POINTING_TABLE_HEADER comment. Blank cells where
    a quantity doesn't apply (e.g. PA/Illum. for Science; RA/Dec for Sun,
    not carried in DRP_ALL), same convention as QualSFrame.py's table.'''
    def g(name, fmt='{:.2f}'):
        if drp is None or name not in drp.colnames:
            return ''
        val = drp[name][i]
        try:
            if not np.isfinite(val):
                return ''
        except TypeError:
            pass
        return fmt.format(val)

    def gs(name):
        if drp is None or name not in drp.colnames:
            return ''
        return str(drp[name][i]).strip()

    illum = g('moon_fli', '{:.2f}')
    illum_pct = f'{100 * float(illum):.1f}' if illum else ''

    return [
        list(_POINTING_TABLE_HEADER),
        ['Science', g('sci_ra'), g('sci_dec'), g('sci_pa'), '', g('sci_alt'),
         '', gs('sci_astsrc'), g('sci_sh_hght', '{:.1f}')],
        ['SkyE', g('skye_ra'), g('skye_dec'), g('skye_pa'), g('sci_skye_sep'), g('skye_alt'),
         '', gs('skye_astsrc'), g('skye_sh_hght', '{:.1f}')],
        ['SkyW', g('skyw_ra'), g('skyw_dec'), g('skyw_pa'), g('sci_skyw_sep'), g('skyw_alt'),
         '', gs('skyw_astsrc'), g('skyw_sh_hght', '{:.1f}')],
        ['Moon', g('moon_ra'), g('moon_dec'), '', g('sci_moon_sep'), g('moon_alt'),
         illum_pct, '', ''],
        ['Sun', '', '', '', '', g('sun_alt'), '', '', ''],
    ]


def plot_row(fits_files, row=None, expnum=None, vel=None, lmc=False, smc=False,
            sigma_guess=1.0, snr_min=5.0, title=None):
    '''
    Build the summary-table + per-line-group x per-file grid for one
    exposure.

    Parameters
    ----------
    fits_files : sequence of str or Path
        SkySub*.py output files, all for the same exposure.
    row : int, optional
    expnum : int, optional
    vel : float, optional
        Explicit nebular systemic velocity (km/s) override. Default:
        None -- the velocity is instead looked up per column from that
        column's own DRP_ALL['Redshift'] (see SkySubNebEval.fit_file's
        docstring and resolve_velocity) -- correct for the normal case
        of all columns being the same exposure/target; only pass vel/
        lmc/smc to force a velocity the DRP_ALL lookup would get wrong.
    lmc, smc : bool
        Shortcuts for LMC_VEL/SMC_VEL, same precedence as vel above.
    sigma_guess : float
    snr_min : float
        Passed to _doublet_ratio for the summary table (see
        SkySubNebEval.py -- same SNR gate that tool's own ratio columns use).
    title : str, optional

    Returns
    -------
    plotly.graph_objects.Figure
    '''
    cols = []
    exp_info_text = None
    pointing_rows = None
    for f in fits_files:
        with fits.open(f) as hdul:
            hdr = hdul[0].header
            routine = hdr.get('TITLE', Path(f).stem)
            wave = np.asarray(hdul['WAVE'].data, dtype=float)
            flux2d = np.atleast_2d(np.asarray(hdul['FLUX'].data, dtype=float))
            drp = Table(hdul['DRP_ALL'].data) if 'DRP_ALL' in hdul else None
            if expnum is not None and drp is not None:
                matches = np.flatnonzero(np.asarray(drp['expnum']) == int(expnum))
                if matches.size == 0:
                    raise ValueError(f'expnum {expnum} not found in {f}')
                i = int(matches[0])
            else:
                i = int(row) if row is not None else 0
            flux = flux2d[i]
            row_redshift = None
            if drp is not None and 'Redshift' in drp.colnames and np.isfinite(drp['Redshift'][i]):
                row_redshift = float(drp['Redshift'][i])
            row_vel = resolve_velocity(vel, lmc, smc, drp_redshift=row_redshift)
            if exp_info_text is None and drp is not None:
                exp_info_text = _exposure_info_text(drp, i, row_vel)
                pointing_rows = _pointing_table_rows(drp, i)
        fits_dict = fit_nebular_row(wave, flux, vel=row_vel, sigma_guess=sigma_guess)
        cols.append(dict(label=routine, wave=wave, flux=flux, fits=fits_dict, vel=row_vel))

    n_line_rows, n_cols = len(GROUPS), len(cols)
    # +1 pointing table, +1 overall-spectrum row, +1 doublet-ratio summary table
    n_rows = n_line_rows + 3

    subplot_titles = ['', 'Overall spectrum by method', '']  # '' = table rows' unused titles
    for gname, names, _pad in GROUPS:
        for col in cols:
            parts = []
            for name in names:
                snr = _snr(col['fits'].get(name))
                parts.append(f'{name} SNR={snr:.1f}' if np.isfinite(snr) else f'{name} --')
            subplot_titles.append(f"{col['label']}: {gname} ({', '.join(parts)})")

    specs = [[{'type': 'table', 'colspan': n_cols}] + [None] * (n_cols - 1),
            [{'colspan': n_cols}] + [None] * (n_cols - 1),
            [{'type': 'table', 'colspan': n_cols}] + [None] * (n_cols - 1)]
    specs += [[{}] * n_cols for _ in GROUPS]
    # Absolute pixel weights (SPECTRUM_ROW_HEIGHT/TABLE_ROW_HEIGHT/ROW_HEIGHT),
    # matched below by setting the figure's total height to their sum -- see
    # TABLE_ROW_HEIGHT's comment for why this can't be a fraction of n_rows.
    row_heights = ([POINTING_TABLE_ROW_HEIGHT, SPECTRUM_ROW_HEIGHT, TABLE_ROW_HEIGHT]
                  + [ROW_HEIGHT] * n_line_rows)
    vspace = min(0.03, 0.9 / max(n_rows - 1, 1))
    fig = make_subplots(rows=n_rows, cols=n_cols, specs=specs, row_heights=row_heights,
                        subplot_titles=subplot_titles, vertical_spacing=vspace,
                        horizontal_spacing=0.04)
    # Applied here, before the color-key annotation is added below, so it
    # only touches the subplot-title annotations (not the color key, which
    # sets its own smaller size).
    fig.update_annotations(font_size=SUBPLOT_TITLE_SIZE)

    # --- Row 1: pointing table (Science/SkyE/SkyW/Moon/Sun geometry) ---
    if pointing_rows is not None:
        pt_header, *pt_body = pointing_rows
        pt_cols = list(zip(*pt_body)) if pt_body else [[] for _ in pt_header]
        fig.add_trace(go.Table(
            header=dict(values=pt_header, font=dict(size=SUBPLOT_TITLE_SIZE - 1), align='center'),
            cells=dict(values=pt_cols, font=dict(size=SUBPLOT_TITLE_SIZE - 1), align='center',
                      height=26)),
            row=1, col=1)

    # --- Row 2: overall spectrum, all methods overlaid, full spectral range ---
    # Median-filtered (SPECTRUM_SMOOTH_PIX), not raw per-pixel flux -- at
    # native sampling the per-pixel noise on this log-y, full-range panel
    # swamps the (generally much smaller) inter-method differences we
    # actually want to see here. The per-line-group panels below still
    # show the raw, unsmoothed flux.
    spec_colors = ['black', 'firebrick', 'seagreen', 'steelblue', 'darkorange', 'purple']
    for j, col in enumerate(cols):
        smoothed = median_filter(col['flux'], size=SPECTRUM_SMOOTH_PIX, mode='nearest')
        yclip = np.clip(smoothed, SPECTRUM_YMIN, None)
        color = spec_colors[j % len(spec_colors)]
        fig.add_trace(go.Scatter(x=col['wave'], y=yclip, mode='lines', name=col['label'],
                                 line=dict(color=color, width=1.5)),
                      row=2, col=1)
    fig.update_yaxes(type='log', range=[np.log10(SPECTRUM_YMIN), np.log10(SPECTRUM_YMAX)],
                     title_text='Flux (erg s⁻¹ cm⁻² Å⁻¹)', row=2, col=1, exponentformat='e')
    fig.update_xaxes(title_text='Wavelength (Å)', row=2, col=1)
    # Color key lives INSIDE the row-1 plot area (x/y domain of that
    # subplot, not figure 'paper' coordinates) so it can never collide
    # with the figure title or subplot-title annotations that live in the
    # top margin -- that collision (spectrum panel overwriting the title
    # text) was the actual bug being fixed here.
    color_key = ' &nbsp; '.join(
        f'<span style="color:{spec_colors[j % len(spec_colors)]}">&#9632;</span> {col["label"]}'
        for j, col in enumerate(cols))
    # No row=/col= kwargs here: add_annotation's row/col resolution walks
    # every subplot in the grid to compute axis anchors and raises a
    # PlotlyKeyError as soon as it reaches the go.Table subplot (same
    # failure mode documented below for add_hline) -- passing explicit
    # 'x domain'/'y domain' refs (which default to the first xaxis/yaxis,
    # i.e. row 1 col 1) sidesteps that resolution entirely.
    info_lines = [f'{SPECTRUM_SMOOTH_PIX}-pixel median filter: {color_key}']
    if exp_info_text:
        info_lines.insert(0, exp_info_text)
    fig.add_annotation(text='<br>'.join(info_lines),
                       xref='x domain', yref='y domain', x=0.5, y=0.98, xanchor='center', yanchor='top',
                       showarrow=False, bgcolor='rgba(255,255,255,0.75)',
                       font=dict(size=SUBPLOT_TITLE_SIZE - 1))

    # --- Summary table: one row per method, one column per doublet ratio ---
    header = ['Method'] + [_doublet_header(d.upper(), kind, val) for d, _a, _b, kind, val, _c in DOUBLETS]
    table_cols = [[c['label'] for c in cols]]
    for _dname, name_a, name_b, kind, value, _conf in DOUBLETS:
        cells = []
        for c in cols:
            ratio, _err, dev = _doublet_ratio(c['fits'], name_a, name_b, kind, value, snr_min=snr_min)
            if not np.isfinite(ratio):
                cells.append('--')
            elif kind == 'fixed' and np.isfinite(dev) and abs(dev) >= 3:
                cells.append(f'{ratio:.3f} ⚠')
            elif kind == 'bounded_above' and np.isfinite(dev) and dev > 0:
                cells.append(f'{ratio:.3f} ⚠')
            else:
                cells.append(f'{ratio:.3f}')
        table_cols.append(cells)
    fig.add_trace(go.Table(
        header=dict(values=header, font=dict(size=SUBPLOT_TITLE_SIZE - 1), align='center'),
        cells=dict(values=table_cols, font=dict(size=SUBPLOT_TITLE_SIZE - 1), align='center',
                  height=26)),
        row=3, col=1)

    # --- Per-line-group panels, one shared y-range per row across all methods ---
    colors = ['firebrick', 'seagreen', 'steelblue']
    for i, (gname, names, pad) in enumerate(GROUPS):
        row_i = i + 4
        raw_vals = []   # observed FLUX samples, all columns -- percentile-bounded below
        fit_vals = []   # fitted Gaussian curve samples, all columns -- already smooth, use min/max
        for j, col in enumerate(cols):
            col_j = j + 1
            wlo, whi = _group_window(names, col['vel'], pad)
            sel = (col['wave'] >= wlo) & (col['wave'] <= whi)
            ww, ff = col['wave'][sel], col['flux'][sel]
            fig.add_trace(go.Scatter(x=ww, y=ff, mode='lines', name='FLUX',
                                     line=dict(color='black', width=TRACE_WIDTH),
                                     showlegend=False),
                          row=row_i, col=col_j)
            finite_ff = ff[np.isfinite(ff)]
            if finite_ff.size:
                raw_vals.append(finite_ff)

            for k, name in enumerate(names):
                fit = col['fits'].get(name)
                if fit is None:
                    continue
                amp = fit['amp']
                if name in ('oii_a', 'oii_b'):
                    # fit_oii_doublet's amp is integrated flux, not peak
                    # amplitude (see its docstring) -- convert for plotting.
                    amp = amp / (fit['sigma0'] * np.sqrt(2 * np.pi))
                yfit = _gaussian(ww, amp, fit['lam0'], fit['sigma0'], fit.get('bkg', 0.0))
                fig.add_trace(go.Scatter(x=ww, y=yfit, mode='lines', name=name,
                                         line=dict(color=colors[k % len(colors)], width=FIT_WIDTH),
                                         showlegend=False),
                              row=row_i, col=col_j)
                finite_yfit = yfit[np.isfinite(yfit)]
                if finite_yfit.size:
                    fit_vals.append(finite_yfit)

            # Plain zero-reference trace, not fig.add_hline(): add_hline/
            # add_shape's row/col resolution scans every subplot in the
            # figure to work out axis references, and raises a
            # PlotlyKeyError as soon as it reaches the row-1 go.Table
            # (which has no xaxis/yaxis at all) -- even when targeting an
            # unrelated row. A real (invisible-in-legend) trace sidesteps
            # that code path entirely.
            fig.add_trace(go.Scatter(x=[wlo, whi], y=[0, 0], mode='lines',
                                     line=dict(color='grey', width=1, dash='dot'),
                                     showlegend=False, hoverinfo='skip'),
                          row=row_i, col=col_j)
            fig.update_xaxes(title_text='Wavelength (Å)', row=row_i, col=col_j)
            fig.update_yaxes(title_text='Flux (erg s⁻¹ cm⁻² Å⁻¹)', row=row_i, col=col_j,
                             exponentformat='e')

        # Apply one shared y-range across every column in this line-group row,
        # so amplitude/shape are directly comparable method-to-method instead
        # of each panel silently auto-scaling to its own data. The raw FLUX
        # trace is bounded by its 2nd/98th percentile ACROSS ALL COLUMNS
        # (pooled), not its literal min/max: a single noisy column (e.g. one
        # method with a much noisier residual) otherwise sets a shared range
        # so wide it flattens every OTHER column's real, well-behaved fit
        # into a near-flat line -- visually reading as "this method's fit is
        # bad" when the fit itself is fine and the shared axis is simply
        # dominated by a different column's outlier pixels. The fitted
        # Gaussian curves (already smooth, no per-pixel noise) still use
        # true min/max.
        if raw_vals or fit_vals:
            lo_candidates, hi_candidates = [], []
            if raw_vals:
                pooled = np.concatenate(raw_vals)
                lo, hi = np.nanpercentile(pooled, [2, 98])
                lo_candidates.append(lo)
                hi_candidates.append(hi)
            if fit_vals:
                pooled_fit = np.concatenate(fit_vals)
                lo_candidates.append(float(np.nanmin(pooled_fit)))
                hi_candidates.append(float(np.nanmax(pooled_fit)))
            row_ymin, row_ymax = min(lo_candidates), max(hi_candidates)
            if row_ymax > row_ymin:
                pad_y = 0.08 * (row_ymax - row_ymin)
                yr = [row_ymin - pad_y, row_ymax + pad_y]
                for j in range(n_cols):
                    fig.update_yaxes(range=yr, row=row_i, col=j + 1)

    fig.update_layout(
        title=dict(text=title or 'SkySubNebEval', y=1, yanchor='top', font=dict(size=20)),
        height=sum(row_heights) + MARGIN_T + MARGIN_B, width=380 * n_cols,
        margin=dict(t=MARGIN_T, b=MARGIN_B),
        showlegend=False)
    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Visualize SkySubNebEval.py\'s per-line/doublet fits for one exposure.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fits_file', nargs='+', help='SkySub*.py output FITS file(s)')
    parser.add_argument('-row', type=int, default=None, help='Row index to select (default: 0)')
    parser.add_argument('-expnum', type=int, default=None,
                        help='Select by DRP_ALL expnum instead of row index')
    parser.add_argument('-v', dest='vel', type=float, default=None,
                        help='Nebular systemic velocity (km/s), overriding the default per-row '
                             "DRP_ALL['Redshift'] auto-lookup")
    parser.add_argument('-lmc', action='store_true', help=f'Use the LMC velocity (~{LMC_VEL:.0f} km/s)')
    parser.add_argument('-smc', action='store_true', help=f'Use the SMC velocity (~{SMC_VEL:.0f} km/s)')
    parser.add_argument('-sigma', dest='sigma_guess', type=float, default=1.0,
                        help='Initial Gaussian sigma guess (Angstrom)')
    parser.add_argument('-snr_min', type=float, default=5.0,
                        help='Minimum per-line SNR (both lines) for the summary table\'s ratios')
    parser.add_argument('-title', default=None)
    parser.add_argument('-outfile', default=None,
                        help='Output HTML path (default: Overview_Plot/nebeval_<expnum>.html)')
    args = parser.parse_args()

    fig = plot_row(args.fits_file, row=args.row, expnum=args.expnum,
                   vel=args.vel, lmc=args.lmc, smc=args.smc,
                   sigma_guess=args.sigma_guess, snr_min=args.snr_min, title=args.title)

    label = args.expnum if args.expnum is not None else (args.row if args.row is not None else 0)
    outpath = Path(args.outfile) if args.outfile else Path('Overview_Plot') / f'nebeval_{label}.html'
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(outpath))
    print(f'Wrote {outpath}')


if __name__ == '__main__':
    main()
