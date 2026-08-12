Extracting Spectra from Custom DS9 Regions
===========================================

This page describes an ad hoc workflow for extracting a single combined
spectrum from an arbitrary region (or set of regions) that you have drawn
by eye in DS9, or generated with your own program, on top of a single LVM
RSS file.

Use this workflow when you have one region (or a handful) to extract from
one file. If instead you have a catalog of many sources spread across many
exposures, use the batch pipeline in :doc:`mc_snr_analysis` instead, which
automates the same underlying steps (cross-matching, per-source snapshot
combination, and batch region/spectrum generation).

The input RSS file can be any of:

- a single ``lvmSFrame`` file
- a single ``lvmCFrame`` file
- the output of :doc:`rss_combining`'s ``rss_combine.py`` or
  ``rss_combine_pos.py``

All three share the same ``FLUX``/``IVAR``/``MASK``/``WAVE``/``SLITMAP``
extension layout, so the tools below work identically on any of them. If
you have several exposures of the same region, combine them first with
``rss_combine.py``/``rss_combine_pos.py`` (see :doc:`rss_combining`) so
the steps below only need to run once, on the combined file.

The workflow has four steps:

1. :ref:`region-step1` — draw the region(s) in DS9
2. :ref:`region-step2` — convert the region file into a "masterfile" table
3. :ref:`region-step3` — assign each science fiber to a region (``MakeLVMReg.py``)
4. :ref:`region-step4` — extract the combined spectrum (``GetRegSpec.py``)


.. _region-step1:

Step 1 — Draw Regions in DS9
------------------------------

Open the LVM data (or a broadband image of the field, e.g. from
``LSnap.py``) in DS9 and draw one or more regions marking the area(s) you
want to extract. Supported shapes are ``circle``, ``ellipse``, ``box``,
and ``annulus``.

For each region:

- Give it a name via the region's **text** field (DS9: double-click the
  region → Text). This becomes ``Source_name`` in the masterfile below.
  Multiple regions may share the same name to build a composite selection
  (e.g. two overlapping ellipses that together define one source).
- Optionally set the region's **color**. Color is not interpreted by the
  extraction tools directly, but it is useful for keeping source and
  background regions visually distinct while you work, and can be carried
  through to the masterfile's ``Color`` column.

Save the regions as a DS9 region file (**Region → Save Regions**, format
``ds9``, coordinate system ``fk5``).


.. _region-step2:

Step 2 — Build a Masterfile Table
------------------------------------

The region-assignment tool (``MakeLVMReg.py``, step 3) does not read the
DS9 region file directly. It reads an ASCII table — referred to as a
**masterfile** — readable by ``astropy.io.ascii.read``, with one row per
region:

==============  =====================================================================
Column          Description
==============  =====================================================================
Source_name     Name shared by all rows belonging to one source (see Step 1)
RA              Region center, degrees
Dec             Region center, degrees
RegType         ``circle``, ``ellipse``, ``box``, or ``annulus``
Major           Semi-major axis / radius / full width, arcsec (see below)
Minor           Semi-minor axis / inner radius / full height, arcsec (see below)
Theta           Position angle, degrees East of North (ignored for ``circle``)
==============  =====================================================================

Per-``RegType`` meaning of ``Major``/``Minor``:

- ``circle``: ``Major`` is the radius; ``Minor`` is unused.
- ``ellipse``: ``Major``/``Minor`` are the semi-major/semi-minor axes.
- ``box``: ``Major``/``Minor`` are the full width/height.
- ``annulus``: ``Major`` is the outer radius, ``Minor`` is the inner radius.

An optional ``SourceBack`` column (values ``'Source'`` or ``'Back'``)
marks which rows of a given ``Source_name`` are the source region(s) and
which are the background region(s). This is what enables background
subtraction in Steps 3–4 — see :doc:`region_files` for an automated way to
generate a background annulus for every source (``GenAnnularBackground.py``),
if that fits your case better than drawing it by hand. Without this
column, all rows for a source are treated as source fibers and no
background is subtracted.

``reg2master.py`` builds this table directly from the DS9 region file
saved in Step 1, reading each region's shape, position, size, ``text=``
name, and ``color=``:

**Command**::

    reg2master.py regionfile.reg [masterfile]

If ``masterfile`` is omitted, the output is named ``<regionfile>.txt``.
The output table has the columns above (with ``Color`` taken from each
region's DS9 color) but **no** ``SourceBack`` column — add that column by
hand afterward if you want background subtraction, e.g. by giving the
source and background regions in Step 1 distinct colors and setting
``SourceBack`` based on ``Color`` before moving on to Step 3. Multiple
rows sharing one ``Source_name`` (e.g. a source region plus a background
annulus, both named the same in DS9) are preserved as separate rows, as
needed for the composite selections described above.

The inverse conversion, ``master2reg.py masterfile [regionfile]``, turns a
masterfile back into a DS9 region file — useful for visually checking a
hand-edited masterfile (e.g. after adding a ``SourceBack`` column) before
running Step 3.

**Example masterfile** (source + background annulus for one SNR, written
as ``ascii.fixed_width_two_line``, with ``SourceBack`` added by hand)::

    Source_name    RA          Dec        RegType    Major    Minor   Theta   SourceBack
    ------------   ---------   --------   --------   ------   -----   -----   ----------
    J0056-7209     14.150000   -72.150    ellipse     79.03    77.48   90.0   Source
    J0056-7209     14.150000   -72.150    annulus    193.51   114.03    0.0   Back

See :doc:`api/reg2master/index` and :doc:`api/master2reg/index` for API
documentation.


.. _region-step3:

Step 3 — Assign Fibers to Regions (``MakeLVMReg.py``)
---------------------------------------------------------

``MakeLVMReg.py`` reads the RSS file's ``SLITMAP`` extension, tests every
good science fiber against each region in the masterfile (expanding each
region boundary by a fixed buffer so that fibers with any overlap are
included), and writes a new DS9 region file with **one small circle per
fiber**, colored to show the assignment: red = source, green =
background, yellow = everything else.

**Command**::

    MakeLVMReg.py rss_file.fits masterfile.txt

**Output**

- ``<rss_root>.<source_name>.reg`` — one per unique ``Source_name`` in the
  masterfile, written to the current directory.

**Current behavior/limitations** (single-file mode, i.e. without ``-all``):

- The fiber-containment buffer is fixed at 17.5 arcsec (one LVM fiber
  radius) and is not currently adjustable from the command line in this
  mode.
- Output region files are always written to the current directory, named
  from the input RSS filename and the source name; there is currently no
  way to redirect them elsewhere or rename them in this mode.

**Verification**

- Load the output ``.reg`` file in DS9 on top of the RSS pointing (or the
  broadband cutout used in Step 1) and confirm red fibers cover the
  intended source area and green fibers (if any) form the intended
  background region, with no unexpected gaps or overlaps.

See :doc:`api/MakeLVMReg/index` for API documentation.


.. _region-step4:

Step 4 — Extract the Spectrum (``GetRegSpec.py``)
------------------------------------------------------

``GetRegSpec.py`` reads the per-fiber region file from Step 3, selects the
fibers of the requested color, and combines them into a single spectrum.
If a background color/region is also given, its median spectrum is
subtracted.

**Command (source only, sum of fibers)**::

    GetRegSpec.py rss_file.fits rss_root.source_name.reg red

**Command (with background subtraction)**::

    GetRegSpec.py rss_file.fits rss_root.source_name.reg red rss_root.source_name.reg green

The background region file is often the same file as the source region
file, since ``MakeLVMReg.py`` writes red/green/yellow assignments into one
file per source.

**Options:**

-root name
    Prepend a root name to the output file (default: ``Spec``).

-ave
    Return the per-fiber average spectrum instead of the sum.

-med
    Return the per-fiber median spectrum instead of the sum.

With no ``-ave``/``-med`` flag, the default is the **sum**: the per-fiber
average is computed and multiplied by the total number of fibers in the
region (bad fibers are assumed to equal the mean of the good ones). The
background spectrum, when subtracted, is always the median of the
background fibers regardless of this setting.

**Output**

- ``<root>_<rss_root>_<ave|med|sum>[_back].txt`` — ASCII spectrum table.

Output spectrum columns:

==================  ================================================
Column              Description
==================  ================================================
WAVE                Wavelength in Angstroms
FLUX                Source flux (background-subtracted, if requested)
ERROR                Combined source (and background) uncertainty
SOURCE_FLUX         Raw source spectrum before subtraction (only if background given)
SOURCE_ERROR        Uncertainty on source spectrum (only if background given)
BACK_FLUX           Background spectrum (only if background given)
BACK_ERROR           Uncertainty on background spectrum (only if background given)
SKY                 Sky spectrum, if present in the input file
SKY_ERROR            Sky spectrum uncertainty, if present in the input file
MASK                Sum of mask flags across selected fibers
LSF                 Mean line spread function, if present in the input file
==================  ================================================

**Verification**

- Plot ``FLUX`` vs ``WAVE`` and confirm expected emission/absorption
  features are present.
- If background subtraction was used, compare ``SOURCE_FLUX`` and
  ``BACK_FLUX`` to confirm the subtraction removed continuum without
  over-subtracting line emission.

See :doc:`api/GetRegSpec/index` for API documentation.


Notes
-----

- If you only need a simple circular (or annular) aperture centered on a
  known position — no hand-drawn shapes, no masterfile — ``GetSpec.py``
  (see :doc:`visualization`) does the position-to-spectrum extraction
  directly from an RA/Dec and radius, without any of the steps above.
- This page documents the ad hoc, single-file path through
  ``MakeLVMReg.py``/``GetRegSpec.py``. Both scripts also have a batch
  (``-all``) mode, driven by an annular-background catalog and a
  directory of per-source snapshot files, documented in
  :doc:`mc_snr_analysis`.


See Also
--------

- :doc:`rss_combining` — combine multiple exposures before extraction
- :doc:`region_files` — generate paired source/background annulus rows automatically
- :doc:`mc_snr_analysis` — the batch/catalog-driven version of this workflow
- :doc:`visualization` — ``GetSpec.py``, a simpler alternative for plain circular apertures
- :doc:`api/MakeLVMReg/index` — API documentation for ``MakeLVMReg.py``
- :doc:`api/GetRegSpec/index` — API documentation for ``GetRegSpec.py``
- :doc:`api/reg2master/index` — API documentation for ``reg2master.py``
- :doc:`api/master2reg/index` — API documentation for ``master2reg.py``
