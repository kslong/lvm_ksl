Magellanic Cloud SNR Analysis
==============================

This page describes an end-to-end procedure for analysing individual
sources — primarily supernova remnants (SNRs) — observed with LVM in
the Magellanic Clouds.  The procedure verifies that snapshot combination,
background region assignment, fiber selection, and spectral extraction
all work as expected.

The four steps are:

1. :ref:`mc-step1` — combine dithered exposures into per-source RSS snapshots
2. :ref:`mc-step2` — pair each source with an annular background region
3. :ref:`mc-step3` — map source and background regions onto individual fibers
4. :ref:`mc-step4` — extract background-subtracted spectra


Starting Point: the Input Files
---------------------------------

Two files drive the entire procedure.

**Source catalog** (e.g. ``smc_snr_cotton24.txt``)
    One row per source.  Contains the source name, sky position, region
    geometry, and physical/radio properties of each SNR.  Key columns:

    ============  =========================================================
    Column        Description
    ============  =========================================================
    Source_name   Unique identifier (e.g. ``J0056-7209``)
    RA, Dec       Source center in degrees
    RegType       Region shape: ``ellipse``, ``circle``, ``box``, ``annulus``
    Major         Semi-major axis in arcseconds
    Minor         Semi-minor axis in arcseconds
    Theta         Position angle in degrees
    Color         Display color (``yellow`` = unassigned initially)
    ============  =========================================================

**Observed file** (e.g. ``smc_snr_observed.txt``)
    All columns from the source catalog, plus one row per exposure per
    source (multiple rows for sources covered by more than one pointing).
    Additional columns include:

    ============  =========================================================
    Column        Description
    ============  =========================================================
    expnum        Exposure number
    mjd           Modified Julian Date of the observation
    exptime       Exposure time in seconds
    tileid        LVM tile identifier
    location      Relative path to the SFrame FITS file
    ============  =========================================================

    This file is typically produced by cross-matching the source catalog
    against the LVM DRP-All file.


.. _mc-step1:

Step 1 — Create Spectral Snapshots (``rss_snap.py``)
------------------------------------------------------

``rss_snap.py`` reads the observed file, groups rows by ``Source_name``,
combines all exposures for each source, runs Gaussian emission-line
fitting, and produces a 4-panel diagnostic plot.

**Command**::

    rss_snap.py -all smc_snr_observed.txt

To reprocess sources that have already been run, add ``-redo``.  To use
median rather than mean combination, add ``-med``.  To process a single
source::

    rss_snap.py smc_snr_observed.txt J0056-7209

**Outputs**

- ``Snap/<source_name>.ave.fits`` — combined RSS FITS snapshot
- ``Snap/<source_name>.ave.tab`` — ASCII fiber position table
- ``Snap_gauss/<source_name>.gauss.txt`` — Gaussian fit parameters per fiber
- ``Snap_fig/<source_name>.png`` — 4-panel H-alpha / [SII] / ratio / FWHM map

**Verification**

- Confirm that ``Snap/`` contains one ``.ave.fits`` file per source.
- Open the PNG files in ``Snap_fig/`` and check that the H-alpha and
  [SII] maps show coherent spatial structure consistent with the known
  source morphology.
- Inspect ``Snap_gauss/<source_name>.gauss.txt`` and verify that
  H-alpha fluxes and FWHM values are physically reasonable.

See :doc:`snapshots` for full documentation of this step.


.. _mc-step2:

Step 2 — Generate Annular Background Regions (``GenAnnularBackground.py``)
---------------------------------------------------------------------------

``GenAnnularBackground.py`` reads the source catalog (one row per
source) and creates a paired background annulus for each entry.  It is
important to use the source catalog here rather than the observed file,
because the observed file contains multiple rows per source and would
produce duplicate region entries.

**Command**::

    GenAnnularBackground.py smc_snr_cotton24.txt

**Output**

- ``smc_snr_cotton24.ann_reg.txt`` — interleaved table with one source
  row (``SourceBack='Source'``) and one background annulus row
  (``SourceBack='Back'``) per source, sorted by source number.

The annulus geometry applied to each source:

- Source aperture: ``max(Major, 35")`` × ``max(Minor, 35")``
- Annulus inner radius: ``source_major + 35"``
- Annulus outer radius: ``2 × source_major + 35"``

**Verification**

- Open the output file and confirm that each source row is immediately
  followed by a background annulus row with the same ``Source_name``.
- Verify that annulus radii are larger than the source aperture by the
  expected gap.
- Confirm that background rows have ``Color='green'`` and
  ``RegType='annulus'``.

See :doc:`region_files` for full documentation of this step.


.. _mc-step3:

Step 3 — Assign Fibers to Regions (``MakeLVMReg.py``)
------------------------------------------------------

``MakeLVMReg.py`` reads the snapshot FITS file and the annotated region
table from step 2.  For each source it produces a DS9-format ``.reg``
file in which every good science fiber is colour-coded:

- **red** — fiber falls within the source aperture
- **green** — fiber falls within the background annulus
- **yellow** — fiber is outside both regions

**Command**::

    MakeLVMReg.py Snap/<source_name>.ave.fits smc_snr_cotton24.ann_reg.txt

To process all sources in a single call, pass all snapshot files at once::

    MakeLVMReg.py Snap/*.ave.fits smc_snr_cotton24.ann_reg.txt

An optional ``-root`` label can be appended to output filenames to
distinguish different region runs::

    MakeLVMReg.py -root source Snap/<source_name>.ave.fits smc_snr_cotton24.ann_reg.txt

**Output**

One ``.reg`` file per (snapshot, source) combination, named::

    <snapshot_root>.<source_name>.reg
    # e.g.  J0056-7209.ave.J0056-7209.reg
    # or    J0056-7209.ave.source.J0056-7209.reg  (with -root source)

**Verification**

- Load a snapshot FITS file and the corresponding ``.reg`` file in DS9.
- Confirm that red fibers lie on the source and green fibers form a
  surrounding ring.
- Verify that the number of red and green fibers is consistent with
  the source angular size and the 35 arcsec fiber spacing.

See :doc:`api/MakeLVMReg/index` for API documentation.


.. _mc-step4:

Step 4 — Extract Spectra (``GetRegSpec.py``)
--------------------------------------------

``GetRegSpec.py`` reads the snapshot RSS file and the ``.reg`` file from
step 3.  It averages the spectra of the red (source) fibers, subtracts
the median of the green (background) fibers, and writes the result to
an ASCII table.  The same ``.reg`` file is used for both source and
background because both fiber sets are encoded in a single file using
color labels.

**Command**::

    GetRegSpec.py Snap/<source_name>.ave.fits \
        <source_name>.ave.<source_name>.reg red \
        <source_name>.ave.<source_name>.reg green

To store outputs in a subdirectory::

    GetRegSpec.py -root Snap_spec/<source_name> \
        Snap/<source_name>.ave.fits \
        <source_name>.ave.<source_name>.reg red \
        <source_name>.ave.<source_name>.reg green

**Output**

An ASCII fixed-width table with a name constructed from the root, the
snapshot filename, the combination type (``ave``), and a ``back`` suffix
indicating background subtraction was applied, e.g.
``Snap_spec/J0056-7209_J0056-7209.ave_ave_back.txt``.

Columns:

==================  ================================================
Column              Description
==================  ================================================
WAVE                Wavelength in Angstroms
FLUX                Background-subtracted flux
ERROR               Combined source and background uncertainty
SOURCE_FLUX         Raw source spectrum before subtraction
SOURCE_ERROR        Uncertainty on source spectrum
BACK_FLUX           Background spectrum
BACK_ERROR          Uncertainty on background spectrum
SKY                 Sky spectrum (if present in input)
SKY_ERROR           Sky spectrum uncertainty
MASK                Sum of mask flags across selected fibers
LSF                 Mean line spread function
==================  ================================================

**Verification**

- Plot ``FLUX`` vs ``WAVE`` and confirm that emission lines (H-alpha,
  [NII] 6548/6583, [SII] 6717/6731) are present at the expected velocity.
- Compare ``SOURCE_FLUX`` and ``BACK_FLUX`` to confirm that the background
  subtraction has removed continuum without over-subtracting line emission.
- Check that ``ERROR`` values are consistent with the signal-to-noise
  expected for the source.

See :doc:`api/GetRegSpec/index` for API documentation.


Complete Worked Example
-----------------------

The following sequence processes a single SMC SNR from start to finish::

    # Step 1: create snapshot, fit lines, make diagnostic plot
    rss_snap.py smc_snr_observed.txt J0056-7209

    # Step 2: generate source + background annulus table (run once for all sources)
    GenAnnularBackground.py smc_snr_cotton24.txt

    # Step 3: assign fiber colors for this source
    MakeLVMReg.py Snap/J0056-7209.ave.fits smc_snr_cotton24.ann_reg.txt

    # Step 4: extract background-subtracted spectrum
    GetRegSpec.py -root Snap_spec/J0056-7209 \
        Snap/J0056-7209.ave.fits \
        J0056-7209.ave.J0056-7209.reg red \
        J0056-7209.ave.J0056-7209.reg green

To run the full catalog through steps 1 and 2::

    rss_snap.py -all smc_snr_observed.txt
    GenAnnularBackground.py smc_snr_cotton24.txt


See Also
--------

- :doc:`snapshots` — detailed documentation for ``rss_snap.py``
- :doc:`region_files` — detailed documentation for ``GenAnnularBackground.py``
- :doc:`api/MakeLVMReg/index` — API documentation for ``MakeLVMReg.py``
- :doc:`api/GetRegSpec/index` — API documentation for ``GetRegSpec.py``
