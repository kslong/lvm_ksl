#!/usr/bin/env python
# coding: utf-8


'''
                    Space Telescope Science Institute

Synopsis:

Create an interactive (Plotly) plot of a spectrum extracted from the LVM
data, laid out as stacked wavelength panels, with emission lines from a
reference line list overlaid.  The spectrum should have been sky subtracted.


Command line usage (if any):

    usage: PlotSpecI.py [-h] [-wmin w] [-wmax w] [-width w] [-npanel n]
                        [-frac 0.1] [-min ymin] [-max ymax] [-med] [-delta 1e-15]
                        [-mask] [-mask_file file.fits]
                        [-lines file.txt] [-no_lines]
                        [-sky_lines] [-sky_lines_file file.txt]
                        [-mode sep_back] file [files ...]

    Plots are written to Overview_Plot/<basename>.overview.html.

    Default mode (no -mode flag):
        One or more spectrum files are given.  Each is plotted independently.
        If a file contains a BACK_FLUX column (as produced by GetRegSpec.py),
        the background is overlaid in black (alpha 0.4); no further subtraction
        is done because FLUX in that case is already the background-subtracted
        spectrum.

    -mode sep_back:
        Exactly two files are expected: a source spectrum and a separate
        background spectrum.  The background FLUX is subtracted from the
        source FLUX before plotting.

    Panel layout (defaults match PlotSpec.py: 3600-9559 AA in 750 AA panels,
    i.e. 8 panels):
    -wmin    overall lower wavelength bound (default 3600)
    -wmax    overall upper wavelength bound (default 9559)
    -width   panel width in AA (default 750)
    -npanel  number of panels; overrides -width by setting width=(wmax-wmin)/npanel

    -mask
        Overlay pixels flagged as sky-line-contaminated in data/sky_mask.fits
        (see palace_make_mask.py) as a light-grey trace on top of the black
        spectrum.  The overlay is a legend entry ("Sky-line masked") that can
        be toggled on/off by clicking it; toggling one panel's legend entry
        toggles all panels at once.  Off by default.

    -mask_file file.fits
        Use this mask FITS file (WAVE/MASK extensions, as produced by
        palace_make_mask.py) instead of the data/sky_mask.fits default.
        Implies -mask.

    -lines file.txt
        Overlay emission lines from this reference table.  Defaults to
        data/dap_lines.txt (the 215-line DAP reference table).  The
        wavelength column is taken as the first of Wave, WAVE, Wave_air
        present; the label column is taken as the first of LineID, Name,
        name, Ion, DAP_name present.  Every matched line gets a vertical
        tick mark and
        a hover tooltip; a subset also get an always-attempted static text
        label, decluttered live in the browser as you zoom/pan so a label
        hidden at the default view becomes visible text once you zoom in on
        it (not just hover-only).

    -no_lines
        Disable the line overlay entirely.

    -sky_lines
        Overlay a second, independent line list of strong sky lines as blue
        tick marks (position only, no text label -- unlike the scientific
        -lines overlay in red), drawn alongside it, from data/sky_lines.txt
        as produced by palace_make_mask.py's --line-output.  Hovering over a
        tick still shows the line ID and wavelength.  Intended to flag
        candidate sky-subtraction residuals rather than to identify features.
        Off by default.

    -sky_lines_file file.txt
        Use this sky line list instead of the data/sky_lines.txt default.
        Same table format as -lines (Wave_air/LineID or equivalent
        columns).  Implies -sky_lines.

    Scaling options (mutually exclusive; last one wins if combined):
    -frac   autoscale upper limit to frac * max(FLUX) per panel (default 0.1)
    -max    fix the upper y-limit in all panels (also selects fixed-scale mode)
    -min    fix the lower y-limit in all panels
    -med    centre each panel on the median FLUX in that wavelength range, with limits median +/- delta
    -delta  half-range for -med mode (default 3e-15)

Description:

    Interactive counterpart to PlotSpec.py.  Panel geometry and y-scaling
    logic are ported directly from PlotSpec.py; mask loading is imported
    from PlotSpec.py rather than duplicated.  The line-overlay declutter
    logic runs client-side in the exported HTML (a small script embedded via
    plotly's post_script), so the output stays a single self-contained file
    with no server required.

Primary routines:

    do_all

Notes:

    The label-declutter algorithm (see the embedded JS built by
    _build_post_script) is: for the lines visible in the current x-range of
    a panel, sorted by wavelength, greedily keep a label only if it does not
    collide (in estimated pixel width) with the last kept label; otherwise
    that line keeps its tick mark and hover tooltip but no static text this
    view.  This is re-run on every plotly_relayout event (zoom, pan, double-
    click reset, window resize), so it adapts as the user explores the plot.

History::

    260814 ksl Coding begun
    260815 ksl Added -sky_lines/-sky_lines_file: a second, independent line
                list (default data/sky_lines.txt, as produced by
                palace_make_mask.py's --line-output) overlaid alongside the
                scientific -lines list. Drawn as blue tick marks with a hover
                tooltip but, unlike the scientific list, no static text
                label -- these mark candidate sky-subtraction-residual
                positions rather than identified features, and the strong-
                line list can run to several hundred entries (mostly OH), so
                labeling would defeat the declutter logic. The tick+label
                rendering block in add_panel was factored out into
                _add_line_overlay(color, show_labels=...) so both lists
                share one implementation instead of two near-duplicate
                blocks.

'''


import os
import json
from pathlib import Path
import re

import numpy as np
from astropy.io import ascii

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from PlotSpec import get_sky_mask, DEFAULT_MASK_FILE, _usage_from_doc
from GetSkyCont import _interp_mask_to_wave


DEFAULT_LINES_FILE = Path(__file__).resolve().parent.parent / 'data' / 'dap_lines.txt'
DEFAULT_SKY_LINES_FILE = Path(__file__).resolve().parent.parent / 'data' / 'sky_lines.txt'

WAVE_COLS = ('Wave', 'WAVE', 'Wave_air')
NAME_COLS = ('LineID', 'Name', 'name', 'Ion', 'DAP_name')

FIG_WIDTH = 1000
PANEL_HEIGHT = 130
MARGIN = dict(l=70, r=30, t=40, b=60)


def load_lines(filename=None):
    '''
    Read a line-list reference table (defaulting to data/dap_lines.txt) and
    return (names, waves) arrays.  The wavelength column is the first of
    Wave/WAVE/Wave_air present; the label column is the first of
    LineID/Name/name/Ion/DAP_name present.  Raises ValueError if neither can
    be found, rather than silently skipping the overlay.
    '''
    path = filename if filename else DEFAULT_LINES_FILE
    tab = ascii.read(str(path))
    wave_col = next((c for c in WAVE_COLS if c in tab.colnames), None)
    name_col = next((c for c in NAME_COLS if c in tab.colnames), None)
    if wave_col is None or name_col is None:
        raise ValueError(
            'Could not find a wavelength/name column in %s (have columns %s; '
            'looked for wavelength in %s and name in %s)'
            % (path, tab.colnames, WAVE_COLS, NAME_COLS))
    names = [str(x) for x in tab[name_col]]
    waves = np.asarray(tab[wave_col], dtype=float)
    return names, waves


def _slice_region(spectab, wmin, wmax, extra=10):
    xx = spectab[spectab['WAVE'] > wmin - extra]
    xx = xx[xx['WAVE'] < wmax + extra]
    finite = np.isfinite(xx['FLUX'])
    return xx[finite]


def _panel_ylim(xx, ptype, ymin, ymax, frac, med_delta):
    '''
    Same scaling math as PlotSpec.py's do_one_region/do_one_region_fixed/
    do_one_region_med, returning (ylo,yhi) or None (meaning: let Plotly
    autorange).
    '''
    if ptype == 'scale':
        if frac < 1.0 and len(xx) > 0:
            yy_max = np.max(xx['FLUX'])
            yy_min = np.median(xx['FLUX']) - 0.05 * yy_max
            return yy_min, frac * yy_max
        return None
    elif ptype == 'fixed':
        return ymin, ymax
    elif ptype == 'med':
        if len(xx) == 0:
            return None
        ymed = np.median(xx['FLUX'])
        return ymed - 0.25 * med_delta, ymed + med_delta
    else:
        raise ValueError('Indecipherable plot type: %s' % ptype)


def _axis_suffix(row):
    return '' if row == 1 else str(row)


def _add_line_overlay(fig, row, wmin, wmax, extra, top_y, line_names, line_waves,
                      color, meta, show_labels=True):
    '''
    Add tick marks and an always-hoverable marker trace for one line list to
    panel `row`, in `color`.  If show_labels is True, also add placeholder
    text labels (decluttered client-side, see _build_post_script), appending
    one dict per label to `meta` in the same order fig.add_annotation is
    called, so the JS can match gd.layout.annotations[i] to meta[i].  If
    show_labels is False, only the tick marks and hover tooltip are drawn --
    useful for a dense reference list (e.g. sky lines) where only the
    position matters, not a readable label.  No-op if line_names is None or
    no line in this list falls in [wmin-extra, wmax+extra].
    '''
    if line_names is None:
        return

    panel_idx = [i for i, w in enumerate(line_waves) if wmin - extra <= w <= wmax + extra]
    if not panel_idx:
        return

    panel_waves = [line_waves[i] for i in panel_idx]
    panel_names = [line_names[i] for i in panel_idx]

    for w in panel_waves:
        fig.add_shape(type='line', x0=w, x1=w, y0=0, y1=1,
                      xref='x%s' % _axis_suffix(row), yref='y%s domain' % _axis_suffix(row),
                      line=dict(color=color, width=1, dash='dot'), opacity=0.35,
                      row=row, col=1)

    fig.add_trace(go.Scatter(x=panel_waves, y=[top_y] * len(panel_waves), mode='markers',
                              marker=dict(size=10, opacity=0), showlegend=False,
                              hoverinfo='text',
                              hovertext=['%s  %.2f Å' % (nm, w) for nm, w in zip(panel_names, panel_waves)]),
                  row=row, col=1)

    if not show_labels:
        return

    for nm, w in zip(panel_names, panel_waves):
        fig.add_annotation(x=w, y=top_y, text=nm, showarrow=False, visible=True,
                           textangle=-90, font=dict(size=9, color=color),
                           xanchor='center', yanchor='bottom', row=row, col=1)
        meta.append({'row': row, 'wave': float(w), 'name': str(nm)})


def add_panel(fig, row, spectab, wmin, wmax, ptype, ymin, ymax, frac, med_delta,
              mask, mask_file, line_names, line_waves,
              sky_line_names, sky_line_waves, show_mask_legend, meta):
    '''
    Add one panel (row) to fig: the flux trace (plus SOURCE_FLUX/BACK_FLUX
    overlays and the sky-mask overlay if requested), the scientific line
    overlay (red) from line_names/line_waves, and the sky line overlay
    (blue) from sky_line_names/sky_line_waves.  Each overlay is tick
    shapes, an always-hoverable marker trace, and placeholder annotations
    that the embedded JS declutters; see _add_line_overlay.
    '''
    extra = 10
    xx = _slice_region(spectab, wmin, wmax, extra)
    wave = np.asarray(xx['WAVE'], dtype=float)
    flux = np.asarray(xx['FLUX'], dtype=float)

    if 'SOURCE_FLUX' in xx.colnames:
        fig.add_trace(go.Scatter(x=wave, y=np.asarray(xx['SOURCE_FLUX'], dtype=float),
                                  mode='lines', line=dict(color='black', width=1),
                                  opacity=0.2, showlegend=False, hoverinfo='skip'),
                      row=row, col=1)
    if 'BACK_FLUX' in xx.colnames:
        fig.add_trace(go.Scatter(x=wave, y=np.asarray(xx['BACK_FLUX'], dtype=float),
                                  mode='lines', line=dict(color='black', width=1),
                                  opacity=0.4, showlegend=False, hoverinfo='skip'),
                      row=row, col=1)

    fig.add_trace(go.Scatter(x=wave, y=flux, mode='lines',
                              line=dict(color='black', width=1),
                              name='Flux', showlegend=False, hoverinfo='skip'),
                  row=row, col=1)

    sky_mask = get_sky_mask(mask_file) if mask else None
    if sky_mask is not None:
        mask_wave, mask_bool = sky_mask
        clean = _interp_mask_to_wave(mask_wave, mask_bool, wave)
        masked_flux = np.where(~clean, flux, np.nan)
        fig.add_trace(go.Scatter(x=wave, y=masked_flux, mode='lines',
                                  line=dict(color='lightgrey', width=1.5),
                                  name='Sky-line masked', legendgroup='mask',
                                  showlegend=show_mask_legend, hoverinfo='skip'),
                      row=row, col=1)

    ylim = _panel_ylim(xx, ptype, ymin, ymax, frac, med_delta)
    fig.update_xaxes(range=[wmin - extra, wmax + extra], row=row, col=1)
    if ylim is not None:
        fig.update_yaxes(range=[float(ylim[0]), float(ylim[1])], row=row, col=1)
        top_y = ylim[0] + 0.92 * (ylim[1] - ylim[0])
    else:
        finite_flux = flux[np.isfinite(flux)]
        top_y = float(np.max(finite_flux) * 0.92) if len(finite_flux) else 1.0

    _add_line_overlay(fig, row, wmin, wmax, extra, top_y, line_names, line_waves,
                      'firebrick', meta)
    _add_line_overlay(fig, row, wmin, wmax, extra, top_y, sky_line_names, sky_line_waves,
                      'steelblue', meta, show_labels=False)


def _build_post_script(meta):
    '''
    JS injected into the exported HTML (via plotly's post_script) that
    declutters the line-overlay annotations client-side: on load and on
    every plotly_relayout (zoom/pan/reset/resize), it recomputes -- per
    panel, from that panel's currently visible x-range and the actual
    rendered plot width -- which labels have room for static text, using
    the same greedy left-to-right collision pass described in this script's
    Notes.  A reentrancy flag stops the programmatic Plotly.relayout call
    used to apply the result from retriggering itself.
    '''
    meta_json = json.dumps(meta)
    return (
        "var gd = document.getElementById('{plot_id}');\n"
        "var lsMeta = " + meta_json + ";\n"
        "var lsUpdating = false;\n"
        "var lsPxPerChar = 6.5, lsPad = 6;\n"
        "function lsAxisKey(row) { return row === 1 ? 'xaxis' : ('xaxis' + row); }\n"
        "function lsRecompute() {\n"
        "  if (lsUpdating || !gd._fullLayout) return;\n"
        "  var fl = gd._fullLayout;\n"
        "  var plotWidthPx = fl.width - fl.margin.l - fl.margin.r;\n"
        "  var byRow = {};\n"
        "  for (var i = 0; i < lsMeta.length; i++) {\n"
        "    var r = lsMeta[i].row;\n"
        "    if (!byRow[r]) byRow[r] = [];\n"
        "    byRow[r].push(i);\n"
        "  }\n"
        "  var newVisible = new Array(lsMeta.length).fill(false);\n"
        "  Object.keys(byRow).forEach(function(rowStr) {\n"
        "    var row = parseInt(rowStr, 10);\n"
        "    var ax = fl[lsAxisKey(row)];\n"
        "    if (!ax || !ax.range) return;\n"
        "    var x0 = ax.range[0], x1 = ax.range[1];\n"
        "    if (x1 === x0) return;\n"
        "    var pxPerUnit = plotWidthPx / (x1 - x0);\n"
        "    var idxs = byRow[row].slice();\n"
        "    idxs.sort(function(a, b) { return lsMeta[a].wave - lsMeta[b].wave; });\n"
        "    var lastEdge = -Infinity;\n"
        "    idxs.forEach(function(i) {\n"
        "      var w = lsMeta[i].wave;\n"
        "      if (w < x0 || w > x1) return;\n"
        "      var xPx = (w - x0) * pxPerUnit;\n"
        "      var halfW = (lsMeta[i].name.length * lsPxPerChar + lsPad) / 2.0;\n"
        "      if (xPx - halfW >= lastEdge) {\n"
        "        newVisible[i] = true;\n"
        "        lastEdge = xPx + halfW;\n"
        "      }\n"
        "    });\n"
        "  });\n"
        "  var update = {};\n"
        "  for (var i = 0; i < lsMeta.length; i++) {\n"
        "    update['annotations[' + i + '].visible'] = newVisible[i];\n"
        "  }\n"
        "  lsUpdating = true;\n"
        "  Plotly.relayout(gd, update).then(function() { lsUpdating = false; });\n"
        "}\n"
        "gd.on('plotly_relayout', lsRecompute);\n"
        "gd.on('plotly_autosize', lsRecompute);\n"
        "window.addEventListener('resize', lsRecompute);\n"
        "lsRecompute();\n"
    )


def do_all(xtab, wmin=3600, wmax=9559, width=750, npanel=None,
          ptype='scale', ymin=0.0, ymax=1e-14, frac=0.1, med_delta=3e-15,
          title='', mask=False, mask_file=None, lines_file=None, no_lines=False,
          sky_lines=False, sky_lines_file=None):
    '''
    Build the interactive figure: stacked wavelength panels (default
    geometry matches PlotSpec.py) with optional sky-mask overlay, a
    scientific line-list overlay (red), and a sky line-list overlay (blue).
    Returns (fig, post_script).
    '''
    if npanel:
        width = (wmax - wmin) / float(npanel)
        nmax = int(npanel)
    else:
        nmax = int((wmax - wmin) / width) + 1

    line_names = line_waves = None
    if not no_lines:
        line_names, line_waves = load_lines(lines_file)

    sky_line_names = sky_line_waves = None
    if sky_lines:
        sky_line_names, sky_line_waves = load_lines(sky_lines_file or DEFAULT_SKY_LINES_FILE)

    fig = make_subplots(rows=nmax, cols=1, vertical_spacing=0.4 / nmax)

    meta = []
    show_mask_legend = True
    for i in range(nmax):
        row = i + 1
        wwmin = wmin + i * width
        wwmax = wwmin + width
        add_panel(fig, row, xtab, wwmin, wwmax, ptype, ymin, ymax, frac, med_delta,
                  mask, mask_file, line_names, line_waves,
                  sky_line_names, sky_line_waves, show_mask_legend, meta)
        if mask:
            show_mask_legend = False

    fig.update_xaxes(title_text='Wavelength (Å)', row=nmax, col=1)
    fig.update_yaxes(title_text='Flux (erg s⁻¹ cm⁻² Å⁻¹)', row=(nmax + 1) // 2, col=1)
    fig.update_yaxes(exponentformat='e', showexponent='all')
    fig.update_layout(width=FIG_WIDTH, height=PANEL_HEIGHT * nmax + MARGIN['t'] + MARGIN['b'],
                      margin=MARGIN, title=title, showlegend=mask,
                      legend=dict(orientation='h', yanchor='bottom', y=1.0))

    post_script = _build_post_script(meta) if meta else None
    return fig, post_script


def steer(argv):
    '''
    Steering routine for creating interactive plots of the spectra.
    '''

    wmin = 3600.
    wmax = 9559.
    width = 750.
    npanel = None
    frac = 0.1
    ymin = 0.
    ymax = 0.
    med_delta = 3e-15
    itype = 'scale'
    mode = ''
    mask = False
    mask_file = None
    lines_file = None
    no_lines = False
    sky_lines = False
    sky_lines_file = None
    filenames = []

    i = 1
    while i < len(argv):
        if argv[i][0:2] == '-h':
            print(_usage_from_doc(__doc__))
            return
        elif argv[i] == '-wmin':
            i += 1
            wmin = eval(argv[i])
        elif argv[i] == '-wmax':
            i += 1
            wmax = eval(argv[i])
        elif argv[i] == '-width':
            i += 1
            width = eval(argv[i])
        elif argv[i] == '-npanel':
            i += 1
            npanel = int(argv[i])
        elif argv[i] == '-frac':
            i += 1
            frac = eval(argv[i])
        elif argv[i] == '-min':
            i += 1
            ymin = eval(argv[i])
        elif argv[i] == '-max':
            i += 1
            ymax = eval(argv[i])
        elif argv[i] == '-med':
            itype = 'med'
        elif argv[i] == '-delta':
            i += 1
            med_delta = eval(argv[i])
        elif argv[i] == '-mask':
            mask = True
        elif argv[i] == '-mask_file':
            i += 1
            mask_file = argv[i]
            mask = True
        elif argv[i] == '-lines':
            i += 1
            lines_file = argv[i]
        elif argv[i] == '-no_lines':
            no_lines = True
        elif argv[i] == '-sky_lines':
            sky_lines = True
        elif argv[i] == '-sky_lines_file':
            i += 1
            sky_lines_file = argv[i]
            sky_lines = True
        elif argv[i] == '-mode':
            i += 1
            mode = argv[i]
        elif argv[i][0] == '-':
            print('Error: Unknown switch: ', argv[i])
            return
        else:
            filenames.append(argv[i])
        i += 1

    if not filenames:
        print('Error: No input files specified')
        return

    if itype != 'med' and ymax > 0.0:
        itype = 'fixed'

    os.makedirs('Overview_Plot', exist_ok=True)

    def _plot_and_save(xtab, outname, plot_title):
        fig, post_script = do_all(xtab, wmin=wmin, wmax=wmax, width=width, npanel=npanel,
                                  ptype=itype, ymin=ymin, ymax=ymax, frac=frac, med_delta=med_delta,
                                  title=plot_title, mask=mask, mask_file=mask_file,
                                  lines_file=lines_file, no_lines=no_lines,
                                  sky_lines=sky_lines, sky_lines_file=sky_lines_file)
        outfile = 'Overview_Plot/%s.overview.html' % outname
        fig.write_html(outfile, include_plotlyjs=True, post_script=post_script,
                       default_width=FIG_WIDTH, default_height=fig.layout.height)
        print('Wrote', outfile)

    if mode == 'sep_back':
        if len(filenames) < 2:
            print('Error: -mode sep_back requires a source file and a background file')
            return
        filename, backname = filenames[0], filenames[1]
        try:
            xtab = ascii.read(filename)
        except Exception:
            print('Error: Could not read %s' % filename)
            return
        try:
            btab = ascii.read(backname)
        except Exception:
            print('Error: Could not read %s' % backname)
            return
        xtab['FLUX'] -= btab['FLUX']
        outname = os.path.basename(filename).replace('.txt', '').replace('.tab', '')
        _plot_and_save(xtab, outname, os.path.basename(filename))
    else:
        for filename in filenames:
            try:
                xtab = ascii.read(filename)
            except Exception:
                print('Error: Could not read %s, skipping' % filename)
                continue
            outname = os.path.basename(filename).replace('.txt', '').replace('.tab', '')
            _plot_and_save(xtab, outname, os.path.basename(filename))


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
