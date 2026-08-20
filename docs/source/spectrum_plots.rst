Spectrum Overview Plots
========================

Three related tools turn an extracted, sky-subtracted 1D spectrum (a
``WAVE``/``FLUX`` ascii table, as produced by ``GetSpec.py`` or
``GetRegSpec.py``) into a multi-panel overview plot spanning the full LVM
wavelength range, with common emission lines marked:

- ``PlotSpec.py`` - static (matplotlib) plot, 8 stacked panels, portrait
  layout
- ``PlotSpec3.py`` - static (matplotlib) plot, 3 stacked panels, 16:9
  landscape layout sized for presentation slides
- ``PlotSpecI.py`` - interactive (Plotly) plot, panel geometry adjustable
  from the command line, with a file-driven line overlay

This is a distinct kind of plot from the spatial maps covered in
:doc:`visualization` and :doc:`plotting_outputs` -- those plot a quantity
against RA/Dec across the IFU footprint; the tools on this page plot flux
against wavelength for a single already-extracted spectrum.


Common Concepts
----------------

All three tools share the same input format, sky-masking option, and
y-axis scaling modes -- ``PlotSpecI.py`` imports ``get_sky_mask()``
directly from ``PlotSpec.py`` rather than re-implementing it, and the
scaling math is ported line-for-line.

**Input table:** an ascii table with ``WAVE`` and ``FLUX`` columns.
``SOURCE_FLUX``/``BACK_FLUX`` columns, if present (as produced by
``GetRegSpec.py``), are overlaid automatically.

**Sky-line masking (** ``-mask`` **/** ``-mask_file`` **):** highlights
pixels flagged as sky-line-contaminated in ``data/sky_mask.fits`` (see
``palace_make_mask.py``, :doc:`sky_subtraction`) in light grey. Off by
default.

**Y-axis scaling** (mutually exclusive; last one given wins):

-frac frac
    Autoscale the upper limit to ``frac * max(FLUX)`` per panel (default
    ``0.1``).

-max ymax / -min ymin
    Fix the y-limits in all panels.

-med
    Centre each panel on the median ``FLUX`` in that panel's wavelength
    range, with limits ``median +/- delta``.

-delta delta
    Half-range for ``-med`` mode (default ``3e-15``).


PlotSpec.py
------------

The default, general-purpose overview plot: 8 panels of 750 Å each,
covering 3600-9559 Å, stacked in an 8x12-inch portrait figure.

**Usage**::

    PlotSpec.py [-h] [-frac 0.1] [-min ymin] [-max ymax] [-med] [-delta 1e-15]
               [-mask] [-mask_file file.fits] [-mode sep_back] file [files ...]

**Output:** ``Overview_Plot/<basename>.overview.png``

Emission lines are marked from a hardcoded list of ~39 lines
(``do_lines()``/``xmark()`` in the script), each drawn as red text at a
fixed position near the top of its panel. This list was hand-picked
rather than read from a file, and a handful of its wavelengths are
several Å off from the corresponding DAP-derived values in
``data/dap_lines.txt`` (see :doc:`dap`) -- a few (e.g. ``[SII]`` at
6720 Å, ``[OII]`` at 3728 Å) look like blended-doublet midpoints rather
than individual line wavelengths, and ``HeII`` at 4199.83 Å has no
counterpart in ``data/dap_lines.txt`` at all. This hardcoded list is not
currently kept in sync with ``data/dap_lines.txt``; for a machine-
generated, DAP-sourced line list instead, use ``PlotSpecI.py`` below.

**-mode sep_back:** exactly two files (source, background); the
background ``FLUX`` is subtracted from the source before plotting.


PlotSpec3.py
-------------

Presentation variant of ``PlotSpec.py``: 3 panels of 2000 Å each,
stacked in a 16:9 landscape figure sized for full-screen display.

**Usage**::

    PlotSpec3.py [-h] [-frac 0.1] [-min ymin] [-max ymax] [-med] [-delta 1e-15]
                [-mask] [-mask_file file.fits] [-mode sep_back] file [files ...]

**Output:** ``Overview_Plot/<basename>.overview.png``

Same options, same hardcoded line list, and the same masking/scaling
logic as ``PlotSpec.py`` -- the only difference is panel count/geometry.


PlotSpecI.py
-------------

Interactive (Plotly) counterpart to ``PlotSpec.py``. Panel geometry
defaults to the same layout as ``PlotSpec.py`` (8 panels, 750 Å,
3600-9559 Å) but is adjustable from the command line, and the line
overlay is read from a file rather than hardcoded.

**Usage**::

    PlotSpecI.py [-h] [-wmin w] [-wmax w] [-width w] [-npanel n]
                [-frac 0.1] [-min ymin] [-max ymax] [-med] [-delta 1e-15]
                [-mask] [-mask_file file.fits]
                [-lines file.txt] [-no_lines]
                [-sky_lines] [-sky_lines_file file.txt]
                [-mode sep_back] file [files ...]

**Output:** ``Overview_Plot/<basename>.overview.html`` -- a single
self-contained file (Plotly.js is embedded); no server is needed to
view it, just a browser.

**Panel geometry options:**

-wmin w, -wmax w
    Overall wavelength range (default 3600/9559 Å, matching
    ``PlotSpec.py``).

-width w
    Panel width in Å (default 750, giving 8 panels by default).

-npanel n
    Set the panel count directly (overrides ``-width`` by computing
    ``width = (wmax-wmin)/npanel``).

**Line overlay options:**

-lines file.txt
    Reference table of lines to mark (default ``data/dap_lines.txt``,
    the 215-line DAP-derived table -- see :doc:`dap`). The wavelength
    column is the first of ``Wave``, ``WAVE``, ``Wave_air`` present in
    the file; the label column is the first of ``LineID``, ``Name``,
    ``name``, ``Ion``, ``DAP_name`` present. Any table with columns
    matching one name from each list can be used in place of the
    default.

-no_lines
    Disable the line overlay entirely.

-sky_lines
    Overlay a second, independent line list of strong sky lines (default
    ``data/sky_lines.txt``, as produced by ``palace_make_mask.py``'s
    ``--line-output`` -- see :doc:`sky_subtraction`), drawn as blue tick
    marks alongside the scientific ``-lines`` overlay (red). Unlike the
    scientific list, sky lines get no static text label -- only a tick
    mark and a hover tooltip (name + wavelength) -- since the point is to
    see *where* a strong sky line falls (e.g. to judge whether a feature
    in the spectrum might be a sky-subtraction residual), not to identify
    it by name, and the list can run to several hundred entries (mostly
    OH) in the Z arm. Off by default.

-sky_lines_file file.txt
    Use this sky line list instead of the ``data/sky_lines.txt`` default.
    Same table format as ``-lines``. Implies ``-sky_lines``.

Sky masking (``-mask``) is a legend-toggleable trace here rather than a
baked-in recolor: clicking the "Sky-line masked" legend entry once
shows/hides the masked overlay across all panels at once.


Zoom-Reactive Line Labels
^^^^^^^^^^^^^^^^^^^^^^^^^^

The default line list has 215 entries, far too many to label as static
text all at once without overlapping -- some wavelength regions (e.g.
the ``[FeII]``/``[FeIII]`` forest around 4400-5400 Å) have a dozen or
more lines within a few Å of each other.

Every matched line always gets a thin vertical tick mark and a hover
tooltip (name + wavelength), regardless of local crowding. A subset of
lines also get an always-attempted static text label, and *which*
subset is recomputed live in the browser on every zoom, pan, double-
click reset, or window resize:

1. For each panel, take the lines currently inside its visible x-range,
   sorted by wavelength.
2. Sweep left to right, estimating each label's rendered pixel width
   from its character count; keep a label only if it doesn't collide
   with the last label kept, otherwise leave it as tick-plus-hover-only
   for this view.
3. Apply the result via a targeted Plotly relayout call (one boolean
   per annotation), guarded against re-triggering itself.

The practical effect: a line that's hover-only at the default zoomed-out
view (because its neighbors are too close together) gets its name
rendered as real text as soon as you zoom in far enough to give it room
-- you don't have to already know which tick to hover over to identify
it. The whole mechanism is a small script embedded in the exported HTML
via Plotly's ``post_script`` hook; it runs entirely client-side, so the
output stays one self-contained file.


Choosing Between Them
-----------------------

.. list-table::
   :header-rows: 1

   * - Tool
     - Output
     - Default panels
     - Line source
     - Best for
   * - ``PlotSpec.py``
     - PNG
     - 8, portrait
     - hardcoded ~39
     - quick default look, printable
   * - ``PlotSpec3.py``
     - PNG
     - 3, 16:9 landscape
     - hardcoded ~39
     - presentation slides
   * - ``PlotSpecI.py``
     - HTML
     - 8, portrait (adjustable)
     - file-driven, default 215
     - interactive exploration, identifying crowded/faint lines


See Also
--------

- :doc:`dap` - ``data/dap_lines.txt``, the default line-overlay
  reference table for ``PlotSpecI.py``
- :doc:`spectral_fitting_local` - ``GetSpec.py``/emission-line fitting tools
  that produce the spectra plotted here
- :doc:`sky_subtraction` - ``palace_make_mask.py``, which produces
  ``data/sky_mask.fits`` (the ``-mask`` overlay) and
  ``data/sky_lines.txt`` (the ``-sky_lines`` overlay)
- :doc:`visualization` - Spatial (RA/Dec) image maps, a different kind
  of plot from the wavelength panels on this page
- :doc:`plotting_outputs` - Spatial (RA/Dec) table-column plots, also
  distinct from the wavelength panels on this page
- :doc:`api/PlotSpec/index` - API documentation
- :doc:`api/PlotSpec3/index` - API documentation
- :doc:`api/PlotSpecI/index` - API documentation
