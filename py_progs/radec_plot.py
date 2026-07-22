#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    General-purpose routines for plotting a quantity from an astropy
    Table against RA and Dec, either as a per-row scatter plot or as
    an interpolated image.  Not tied to any fixed set of column names
    or a particular upstream script -- works with any table that has
    RA/Dec columns and a value column, e.g. the per-fiber/per-exposure
    tables produced by lvm_gaussfit.py, sky_gaussfit.py,
    gauss_offset.py, or rss_snap.py.

    Generalizes rss_snap.py's plot_one()/plot_one_interpolated() (RA/Dec
    column names, and an irregular-footprint mask for the interpolated
    version -- see Description).

Command line usage::

    radec_plot.py [-h] [-scatter | -interp] [-ra col] [-dec col]
                  [-ext N] [-ymin val] [-ymax val] [-label text]
                  [-grid N] [-mask_radius deg] [-out file.png]
                  filename var

    Arguments::

        filename   FITS table (or anything astropy.table.Table.read
                   understands) containing RA/Dec columns and var
        var        name of the column to plot

    Options::

        -scatter        one point per row (default)
        -interp         interpolated image
        -ra col         RA column name (default: auto-detect ra/RA/Ra)
        -dec col        Dec column name (default: auto-detect dec/Dec/DEC)
        -ext N          FITS extension/HDU to read (default: 1)
        -ymin val       color scale lower limit (default: auto, 20th pct)
        -ymax val       color scale upper limit (default: auto, 80th pct)
        -label text     colorbar label (default: var)
        -grid N         interpolation grid resolution (default: 100)
        -mask_radius deg  see Description (default: auto)
        -out file.png   output filename (default: <stem>_<var>.png)

Description:

    Primary routines (meant to be imported into notebooks; the CLI
    above is a thin convenience wrapper around the same two calls):

    plot_scatter(table, var, ...)
        One scatter point per row, colored by var.  Always shows
        exactly the rows present in the table -- no assumption about
        the footprint shape.

    plot_interpolated(table, var, ...)
        Interpolates var onto a regular RA/Dec grid (scipy.griddata)
        and renders it as an image -- smoother than the scatter plot,
        but griddata only knows to blank grid points outside the
        data's convex hull.  LVM footprints are frequently concave or
        have interior gaps (multiple pointings, missing fibers), so a
        convex-hull mask alone lets the interpolation smear values
        across regions with no real data.  plot_interpolated() instead
        masks any grid point whose distance to the nearest real data
        point exceeds mask_radius, which follows the true (possibly
        concave) boundary and interior gaps.  If mask_radius is None
        (the default), it is set automatically to twice the median
        nearest-neighbor spacing among the real data points.

    Both routines accept an ax= argument (default: current axis) so
    they can be used to build a mosaic, e.g. with plt.subplots()::

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        plot_interpolated(tab, 'flux_ha',  ax=axes[0, 0])
        plot_interpolated(tab, 'flux_sii', ax=axes[0, 1])

Notes:

    RA is plotted with a cos(dec) correction so angular scale is
    preserved, and the RA axis is inverted (increasing to the left),
    matching astronomical convention.

    Requires scipy (griddata, cKDTree) for plot_interpolated(); no
    matplotlib backend is forced at import time, so interactive display
    works normally when imported into a notebook -- the CLI path saves
    a PNG explicitly instead of relying on a non-interactive backend.

    The CLI's table loader (_read_table) tries Table.read() first, then
    falls back to astropy.io.ascii.read() for formats Table.read cannot
    auto-identify -- notably lvm_gaussfit.py's/sky_gaussfit.py's
    <root>.gauss.txt output (ascii.fixed_width_two_line), which has no
    auto-identify support in Table.read.

History::

    260722 ksl Coding begun.

'''

import os
import re
import sys

import numpy as np
from astropy.table import Table
from astropy.io import ascii as astropy_ascii
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata
from scipy.spatial import cKDTree


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before a line consisting of "History:" (or
    "History::"/"Version History" -- whitespace/colon-insensitive), so
    -h stays short even as that section grows -- without hand-
    duplicating the Synopsis/Options text in a second string.  Anchored
    to a whole line (not a bare substring search) so it can't misfire on
    "History:" appearing mid-sentence, and returns doc unchanged if no
    such line is present.
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


_USAGE = _usage_from_doc(__doc__)


def _read_table(filename, ext=1):
    '''
    Read filename into an astropy Table, e.g. a FITS table or one of
    lvm_gaussfit.py's/sky_gaussfit.py's <root>.gauss.txt output files
    (ascii.fixed_width_two_line -- not auto-identified by Table.read,
    so the more permissive ascii.read guesser is tried as a fallback).
    '''
    if filename.lower().endswith(('.fits', '.fits.gz', '.fit')):
        return Table.read(filename, hdu=ext)
    try:
        return Table.read(filename)
    except Exception:
        return astropy_ascii.read(filename)


def _find_col(table, col, candidates):
    '''
    Resolve a column name in table: use col if given (must exist),
    otherwise return the first of candidates that exists.
    '''
    if col is not None:
        if col in table.colnames:
            return col
        raise KeyError('Column %r not found in table (columns: %s)'
                        % (col, table.colnames))
    for name in candidates:
        if name in table.colnames:
            return name
    raise KeyError('None of %s found in table (columns: %s)'
                    % (candidates, table.colnames))


def _resolve_radec(table, ra_col, dec_col):
    ra_col = _find_col(table, ra_col, ['ra', 'RA', 'Ra'])
    dec_col = _find_col(table, dec_col, ['dec', 'Dec', 'DEC'])
    return ra_col, dec_col


def _resolve_value(table, var):
    if isinstance(var, str):
        return np.asarray(table[var], dtype=float)
    return np.asarray(var, dtype=float)


def _auto_marker_size(ax, x, y, fill_factor=1.3):
    '''
    Estimate a matplotlib scatter `s` (marker area, points^2) so markers
    approximate the median nearest-neighbor spacing of (x, y) at the
    axes' actual rendered size on the page/screen -- avoids a fixed
    marker size bloating dense tables (overlapping circles smear the
    apparent footprint past its true boundary) or leaving gaps in
    sparse ones.  Reuses the same nearest-neighbor spacing calculation
    plot_interpolated() uses for its default mask_radius.  fill_factor
    slightly over-sizes markers so adjacent points touch/overlap a
    little rather than leaving hairline gaps.

    Requires the axes' final data limits and divider/colorbar layout to
    already be set, and forces one draw so the rendered box size is
    known -- see plot_scatter() call order.  Falls back to a fixed
    marker area if there are too few points or the renderer/window
    extent cannot be obtained (e.g. an unusual backend).
    '''
    if len(x) < 2:
        return 240.0

    points = np.column_stack((x, y))
    tree = cKDTree(points)
    nn_dist, _ = tree.query(points, k=2)
    spacing = np.median(nn_dist[:, 1])

    fig = ax.figure
    fig.canvas.draw()
    try:
        bbox = ax.get_window_extent(renderer=fig.canvas.get_renderer())
    except Exception:
        return 240.0

    data_span = abs(ax.get_xlim()[1] - ax.get_xlim()[0])
    if data_span == 0 or bbox.width == 0:
        return 240.0

    points_per_data_unit = (bbox.width / data_span) * 72.0 / fig.dpi
    diameter_points = spacing * points_per_data_unit * fill_factor
    area = np.pi * (diameter_points / 2.0) ** 2
    return float(np.clip(area, 1.0, 2000.0))


def plot_scatter(table, var, ra_col=None, dec_col=None, ymin=0, ymax=0,
                  label='', ax=None, marker_size=None, cmap='viridis'):
    '''
    Scatter-plot a quantity from table against RA/Dec, one point per row.

    Parameters
    ----------
    table : astropy.table.Table
        Must contain RA/Dec columns (see ra_col/dec_col) and var.
    var : str or array-like
        Column name in table, or an array of values (same length as
        table) to plot instead.
    ra_col, dec_col : str or None
        Column names to use; if None, auto-detected from
        ra/RA/Ra and dec/Dec/DEC respectively.
    ymin, ymax : float
        Color scale limits; if ymax == 0 (default) the limits are the
        20th/80th percentile of var.
    label : str
        Colorbar label (default: var if a column name, else 'Value').
    ax : matplotlib.axes.Axes or None
        Axis to plot into (default: current axis) -- pass explicitly
        when building a mosaic with plt.subplots().
    marker_size : float or None
        Marker area in points^2 (matplotlib scatter's `s`).  If None
        (default), auto-computed from the data's nearest-neighbor
        spacing and the axes' actual rendered size -- see
        _auto_marker_size().  Pass a number to override.

    Returns
    -------
    matplotlib.collections.PathCollection
        The scatter artist, e.g. for further customization.
    '''
    ra_col, dec_col = _resolve_radec(table, ra_col, dec_col)
    ra = np.asarray(table[ra_col], dtype=float)
    dec = np.asarray(table[dec_col], dtype=float)
    value = _resolve_value(table, var)

    if ymin < ymax:
        value = np.clip(value, ymin, ymax)

    cos_dec_factor = np.cos(np.radians(np.nanmean(dec)))
    ra_scaled = ra * cos_dec_factor

    if ax is None:
        ax = plt.gca()

    finite_pos = np.isfinite(ra_scaled) & np.isfinite(dec)
    ra_min, ra_max = np.nanmin(ra_scaled), np.nanmax(ra_scaled)
    dec_min, dec_max = np.nanmin(dec), np.nanmax(dec)
    ra_pad = (ra_max - ra_min) * 0.05
    dec_pad = (dec_max - dec_min) * 0.05

    # Set limits/aspect/colorbar layout up front (matching
    # plot_interpolated's convention) so the axes' rendered box is
    # final before _auto_marker_size() measures it.
    ax.set_xlim(ra_max + ra_pad, ra_min - ra_pad)
    ax.set_ylim(dec_min - dec_pad, dec_max + dec_pad)
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    ax.set_aspect('equal')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)

    if marker_size is None:
        marker_size = _auto_marker_size(ax, ra_scaled[finite_pos], dec[finite_pos])

    sc = ax.scatter(ra_scaled, dec, s=marker_size, c=value, cmap=cmap, alpha=0.8)

    if ymax == 0.0:
        sc.set_clim(np.nanpercentile(value, 20), np.nanpercentile(value, 80))
    else:
        sc.set_clim(ymin, ymax)

    cbar = plt.colorbar(sc, cax=cax)
    cbar.set_label(label if label else (var if isinstance(var, str) else 'Value'))

    xlim = ax.get_xlim()
    tick_locs = ax.get_xticks()
    ax.set_xticks(tick_locs)
    ax.set_xticklabels(['%.1f' % (t / cos_dec_factor) for t in tick_locs])
    ax.set_xlim(xlim)

    return sc


def plot_interpolated(table, var, ra_col=None, dec_col=None, ymin=0, ymax=0,
                       label='', ax=None, grid_resolution=100,
                       interpolation_method='linear', mask_radius=None,
                       cmap='viridis'):
    '''
    Interpolate a quantity from table onto a regular RA/Dec grid and
    display it as an image, masked to the true (possibly concave or
    gappy) footprint of the data -- see module Description.

    Parameters
    ----------
    table : astropy.table.Table
        Must contain RA/Dec columns (see ra_col/dec_col) and var.
    var : str or array-like
        Column name in table, or an array of values (same length as
        table) to plot instead.
    ra_col, dec_col : str or None
        Column names to use; if None, auto-detected from
        ra/RA/Ra and dec/Dec/DEC respectively.
    ymin, ymax : float
        Color scale limits; if ymax == 0 (default) the limits are the
        20th/80th percentile of the interpolated grid.
    label : str
        Colorbar label (default: var if a column name, else 'Value').
    ax : matplotlib.axes.Axes or None
        Axis to plot into (default: current axis) -- pass explicitly
        when building a mosaic with plt.subplots().
    grid_resolution : int
        Number of grid points along each axis.
    interpolation_method : str
        Passed to scipy.interpolate.griddata ('linear', 'nearest',
        or 'cubic').
    mask_radius : float or None
        Maximum allowed distance (degrees, in the same cos(dec)-scaled
        RA frame used for interpolation) from a grid point to the
        nearest real data point; farther grid points are blanked
        (NaN).  If None (default), set automatically to twice the
        median nearest-neighbor spacing among the real data points.

    Returns
    -------
    matplotlib.image.AxesImage
        The image artist, e.g. for further customization.
    '''
    ra_col, dec_col = _resolve_radec(table, ra_col, dec_col)
    ra = np.asarray(table[ra_col], dtype=float)
    dec = np.asarray(table[dec_col], dtype=float)
    value = _resolve_value(table, var)

    finite = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(value)
    ra, dec, value = ra[finite], dec[finite], value[finite]

    if ymin < ymax:
        value = np.clip(value, ymin, ymax)

    cos_dec_factor = np.cos(np.radians(np.mean(dec)))
    ra_scaled = ra * cos_dec_factor
    points = np.column_stack((ra_scaled, dec))

    if ax is None:
        ax = plt.gca()

    ra_min, ra_max = ra_scaled.min(), ra_scaled.max()
    dec_min, dec_max = dec.min(), dec.max()
    ra_pad = (ra_max - ra_min) * 0.05
    dec_pad = (dec_max - dec_min) * 0.05

    grid_ra = np.linspace(ra_min - ra_pad, ra_max + ra_pad, grid_resolution)
    grid_dec = np.linspace(dec_min - dec_pad, dec_max + dec_pad, grid_resolution)
    grid_ra_mesh, grid_dec_mesh = np.meshgrid(grid_ra, grid_dec)

    grid_values = griddata(points, value, (grid_ra_mesh, grid_dec_mesh),
                            method=interpolation_method, fill_value=np.nan)

    tree = cKDTree(points)
    if mask_radius is None:
        nn_dist, _ = tree.query(points, k=2)
        mask_radius = 2.0 * np.median(nn_dist[:, 1])
    grid_points = np.column_stack((grid_ra_mesh.ravel(), grid_dec_mesh.ravel()))
    grid_dist, _ = tree.query(grid_points)
    grid_dist = grid_dist.reshape(grid_ra_mesh.shape)
    grid_values = np.where(grid_dist <= mask_radius, grid_values, np.nan)

    extent = [ra_max + ra_pad, ra_min - ra_pad, dec_min - dec_pad, dec_max + dec_pad]
    grid_values = np.fliplr(grid_values)
    im = ax.imshow(grid_values, extent=extent, origin='lower', cmap=cmap,
                    alpha=0.8, aspect='equal')

    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')

    if ymax == 0.0:
        im.set_clim(np.nanpercentile(grid_values, 20), np.nanpercentile(grid_values, 80))
    else:
        im.set_clim(ymin, ymax)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label(label if label else (var if isinstance(var, str) else 'Value'))

    tick_locs = ax.get_xticks()
    valid_ticks = tick_locs[(tick_locs >= extent[1]) & (tick_locs <= extent[0])]
    ax.set_xticks(valid_ticks)
    ax.set_xticklabels(['%.1f' % (t / cos_dec_factor) for t in valid_ticks])

    return im


def main(argv):
    '''
    Command-line entry point -- see _USAGE / module docstring.
    '''
    if not argv or '-h' in argv or '--help' in argv:
        print(_USAGE)
        return

    mode = 'scatter'
    ra_col = None
    dec_col = None
    ext = 1
    ymin = 0
    ymax = 0
    label = ''
    grid_resolution = 100
    mask_radius = None
    out = ''
    positional = []

    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == '-scatter':
            mode = 'scatter'
        elif arg == '-interp':
            mode = 'interp'
        elif arg == '-ra':
            i += 1
            ra_col = argv[i]
        elif arg == '-dec':
            i += 1
            dec_col = argv[i]
        elif arg == '-ext':
            i += 1
            ext = int(argv[i])
        elif arg == '-ymin':
            i += 1
            ymin = float(argv[i])
        elif arg == '-ymax':
            i += 1
            ymax = float(argv[i])
        elif arg == '-label':
            i += 1
            label = argv[i]
        elif arg == '-grid':
            i += 1
            grid_resolution = int(argv[i])
        elif arg == '-mask_radius':
            i += 1
            mask_radius = float(argv[i])
        elif arg == '-out':
            i += 1
            out = argv[i]
        elif arg.startswith('-'):
            print('Error: unknown option "%s"' % arg)
            print(_USAGE)
            return
        else:
            positional.append(arg)
        i += 1

    if len(positional) != 2:
        print('Error: expected filename and var, got %r' % positional)
        print(_USAGE)
        return
    filename, var = positional

    if not os.path.exists(filename):
        print('Error: file not found: %s' % filename)
        return

    table = _read_table(filename, ext=ext)

    fig, ax = plt.subplots(figsize=(8, 7))
    if mode == 'scatter':
        plot_scatter(table, var, ra_col=ra_col, dec_col=dec_col, ymin=ymin,
                     ymax=ymax, label=label, ax=ax)
    else:
        plot_interpolated(table, var, ra_col=ra_col, dec_col=dec_col, ymin=ymin,
                           ymax=ymax, label=label, ax=ax,
                           grid_resolution=grid_resolution, mask_radius=mask_radius)

    ax.set_title(label if label else var)
    plt.tight_layout()

    if out == '':
        stem = os.path.splitext(os.path.basename(filename))[0]
        out = '%s_%s.png' % (stem, var)
    plt.savefig(out, dpi=110)
    print('Wrote %s' % out)


if __name__ == '__main__':
    main(sys.argv[1:])
