Plotting Outputs
=================

The lvm_ksl package is moving toward a common, reusable layer for
plotting analysis results, rather than each script carrying its own
copy of RA/Dec plotting logic. ``radec_plot.py`` is the first piece of
this: general-purpose routines for plotting a quantity from an astropy
Table against RA and Dec. It replaces duplicated plotting code that
previously lived separately in ``rss_snap.py`` and
``plot_sky_gaussfit.py``, and both now delegate to it (see Integration
below).

This is the beginning, not the end, of that effort -- see Future Plans.


Overview
--------

Available tools:

- ``radec_plot.py`` - Plot a table quantity against RA/Dec, as a
  per-row scatter plot or an interpolated image

Unlike the ``visualization.rst`` tools, which build a WCS image
directly from an RSS/CFrame file, ``radec_plot.py`` plots quantities
that already exist as columns in a table -- typically the per-fiber
fit-result tables produced by ``lvm_gaussfit.py``, ``sky_gaussfit.py``,
``gauss_offset.py``, or ``rss_snap.py``. It is not tied to any fixed
set of column names: RA/Dec columns are auto-detected (``ra``/``RA``/
``Ra`` and ``dec``/``Dec``/``DEC``), and any other table column (or a
plain array of the same length as the table) can be plotted.

The module is meant primarily to be imported into notebooks -- both
routines accept an ``ax=`` argument so they can be composed into
mosaics -- with a command-line wrapper included for quick, one-off
looks at a table.


radec_plot.py
--------------

plot_scatter() — Per-Row Scatter Plot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One point per table row, colored by the requested quantity. Always
shows exactly the rows present in the table, with no assumption about
the footprint's shape.

**Python usage**::

    from lvm_ksl import radec_plot

    radec_plot.plot_scatter(table, var, ra_col=None, dec_col=None,
                             ymin=0, ymax=0, label='', ax=None,
                             marker_size=None, cmap='viridis')

**Parameters:**

table
    An astropy Table containing RA/Dec columns and ``var``.

var
    Column name in ``table`` to plot, or an array of values (same
    length as ``table``).

ra_col, dec_col
    Column names to use; if ``None`` (default), auto-detected from
    ``ra``/``RA``/``Ra`` and ``dec``/``Dec``/``DEC`` respectively.

ymin, ymax
    Color scale limits. If ``ymax == 0`` (default), the limits are the
    20th/80th percentile of ``var``.

label
    Colorbar label (default: ``var`` if it's a column name, else
    ``'Value'``).

ax
    Axis to plot into (default: current axis) -- pass explicitly when
    building a mosaic with ``plt.subplots()``.

marker_size
    Marker area in points² (matplotlib scatter's ``s``). If ``None``
    (default), auto-computed from the data's own nearest-neighbor
    spacing and the axes' actual rendered size -- see Auto Marker
    Sizing below. Pass a number to override.


plot_interpolated() — Interpolated Image
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Interpolates the requested quantity onto a regular RA/Dec grid
(``scipy.interpolate.griddata``) and renders it as an image --
smoother than the scatter plot, but with a masking step needed to
avoid misrepresenting the data's true footprint (see Irregular
Footprint Masking below).

**Python usage**::

    radec_plot.plot_interpolated(table, var, ra_col=None, dec_col=None,
                                  ymin=0, ymax=0, label='', ax=None,
                                  grid_resolution=100,
                                  interpolation_method='linear',
                                  mask_radius=None, cmap='viridis')

**Parameters:**

table, var, ra_col, dec_col, ymin, ymax, label, ax
    Same as ``plot_scatter()`` above.

grid_resolution
    Number of grid points along each axis (default 100).

interpolation_method
    Passed to ``scipy.interpolate.griddata``: ``'linear'`` (default),
    ``'nearest'``, or ``'cubic'``.

mask_radius
    Maximum allowed distance (degrees, in the same cos(dec)-scaled RA
    frame used for interpolation) from a grid point to the nearest
    real data point; farther grid points are blanked (``NaN``). If
    ``None`` (default), set automatically to twice the median
    nearest-neighbor spacing among the real data points.


Irregular Footprint Masking
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``griddata`` only knows to blank grid points outside the data's
*convex hull*. LVM footprints are frequently concave (multiple
pointings stitched together) or have interior gaps (missing/failed
fibers), so a convex-hull mask alone lets the interpolation smear
values across sky positions with no real data -- filling in concave
notches as if they were observed.

``plot_interpolated()`` instead masks any grid point whose distance to
the nearest real data point exceeds ``mask_radius``, which follows the
true, possibly concave, boundary and any interior gaps. Verified
against both a synthetic L-shaped (concave) footprint and a regular
hex fiber grid with a wedge cut out of it -- in both cases the masked
region matches the true missing footprint rather than the convex hull.


Auto Marker Sizing
^^^^^^^^^^^^^^^^^^^

Early versions of the RA/Dec scatter plotting in this package (and its
predecessor in ``rss_snap.py``) used a fixed marker size regardless of
how many points were being plotted or how much sky they covered. At
low fiber counts this looked fine, but on a dense table (for example,
~19,000 fibers stacked from many exposures of one source) a fixed
marker size bloats each fiber into an oversized circle -- overlapping
circles visually smear the plotted footprint's edges into rounded
blobs well beyond the true coverage, and no longer agree with the
sharper, mask-derived footprint that ``plot_interpolated()`` shows for
the same table.

When ``marker_size=None`` (the default), ``plot_scatter()`` instead
computes the median nearest-neighbor spacing between points and the
axes' actual rendered size (forcing one draw pass to measure it), and
sizes markers to approximately fill that spacing -- avoiding both
bloat on dense tables and gaps on sparse ones. Pass an explicit number
to ``marker_size`` to override this.


Building Mosaics
^^^^^^^^^^^^^^^^^

Both routines accept an ``ax=`` argument (default: current axis), so
they can be composed into a multi-panel figure with
``plt.subplots()``::

    import matplotlib.pyplot as plt
    from lvm_ksl import radec_plot

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    radec_plot.plot_interpolated(tab, 'flux_ha',  ax=axes[0, 0], label='Hα Flux')
    radec_plot.plot_interpolated(tab, 'flux_sii', ax=axes[0, 1], label='[SII] Flux')
    radec_plot.plot_interpolated(tab, 's2:ha',    ax=axes[1, 0], label='[SII]:Hα')
    radec_plot.plot_interpolated(tab, 'fwhm_ha',  ax=axes[1, 1], label='Hα FWHM')

Without an explicit ``ax=``, both routines fall back to
``plt.gca()``, which also works with the older
``plt.subplot(2, 2, N)`` style used by ``rss_snap.py``'s ``fig1()``.


Command-Line Usage
^^^^^^^^^^^^^^^^^^^

A thin CLI wrapper is included for a quick, one-off look at a table
without writing any notebook code::

    radec_plot.py [-h] [-scatter | -interp] [-ra col] [-dec col]
                  [-ext N] [-ymin val] [-ymax val] [-label text]
                  [-grid N] [-mask_radius deg] [-out file.png]
                  filename var

**Options:**

-scatter
    One point per row (default).

-interp
    Interpolated image.

-ra col, -dec col
    RA/Dec column names (default: auto-detect).

-ext N
    FITS extension/HDU to read (default: 1).

-ymin val, -ymax val
    Color scale limits (default: auto, 20th/80th percentile).

-label text
    Colorbar and title label (default: ``var``).

-grid N
    Interpolation grid resolution (default: 100).

-mask_radius deg
    See Irregular Footprint Masking above (default: auto).

-out file.png
    Output filename (default: ``<stem>_<var>.png``).

**Arguments:**

filename
    A FITS table, or an ASCII table such as one of
    ``lvm_gaussfit.py``'s/``sky_gaussfit.py``'s ``<root>.gauss.txt``
    output files (see Reading Non-FITS Tables below).

var
    Name of the column to plot.

**Example**::

    radec_plot.py -label "Hα Flux" W28.ave.gauss.txt flux_ha


Reading Non-FITS Tables
^^^^^^^^^^^^^^^^^^^^^^^^

The CLI's table loader tries ``astropy.table.Table.read()`` first,
then falls back to ``astropy.io.ascii.read()`` for formats
``Table.read()`` cannot auto-identify. This matters specifically for
``lvm_gaussfit.py``'s and ``sky_gaussfit.py``'s
``ascii.fixed_width_two_line`` output tables (``<root>.gauss.txt``),
which ``Table.read()`` alone cannot open without an explicit
``format=`` argument.


Integration
-----------

``radec_plot.py``'s two plotting routines are now the single
implementation behind the RA/Dec plots produced elsewhere in the
package:

- ``rss_snap.py``'s ``plot_one()``/``plot_one_interpolated()`` (used
  by ``fig1()`` -- see :doc:`snapshots`) are thin wrappers around
  ``radec_plot.plot_scatter()``/``plot_interpolated()``.
- ``plot_sky_gaussfit.py``'s ``plot_panel()`` (see :doc:`data_quality`)
  delegates its spatial rendering to ``radec_plot.plot_scatter()``,
  keeping only its own median-subtraction/MAD-std residual computation
  and per-quantity color range local.

Both call sites gained the cos(dec) correction, RA-axis inversion,
and (where relevant) auto marker sizing described above simply by
switching to the shared implementation.


Typical Workflow
-----------------

::

    # Fit emission lines across an RSS file
    lvm_gaussfit.py lvmSFrame-00012345.fits

    # Quick look at the result from the command line
    radec_plot.py -label "Hα Flux" lvmSFrame-00012345.gauss.txt flux_ha

    # Or, from a notebook, build a mosaic of several lines
    from astropy.table import Table
    from lvm_ksl import radec_plot
    import matplotlib.pyplot as plt

    tab = Table.read('lvmSFrame-00012345.gauss.txt', format='ascii.fixed_width_two_line')
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    radec_plot.plot_scatter(tab, 'flux_ha', ax=axes[0], label='Hα Flux')
    radec_plot.plot_interpolated(tab, 'flux_ha', ax=axes[1], label='Hα Flux')


Future Plans
------------

This module currently covers table-column-vs-RA/Dec plots only.
Planned next steps include bringing image-based outputs (currently
handled independently by the ``visualization.rst`` tools -- ``kslmap.py``,
``quick_map.py``, ``line_map.py``, ``rss2image.py``) into the same
integrated approach, so that spatial maps built from fitted-line
tables and spatial maps built directly from RSS images share a common,
consistent plotting layer.


See Also
--------

- :doc:`visualization` - Image-based spatial maps built directly from RSS/CFrame files
- :doc:`snapshots` - ``rss_snap.py``, whose diagnostic plots now use ``radec_plot.py``
- :doc:`data_quality` - ``plot_sky_gaussfit.py``, whose residual maps now use ``radec_plot.py``
- :doc:`spectral_fitting` - ``lvm_gaussfit.py``/``sky_gaussfit.py`` produce the tables typically plotted here
- :doc:`api/radec_plot/index` - API documentation
