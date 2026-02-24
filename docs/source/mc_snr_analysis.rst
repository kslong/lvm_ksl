Magellanic Cloud SNR Analysis
==============================

This page describes an end-to-end procedure for analysing individual
sources — primarily supernova remnants (SNRs) — observed with LVM in
the Magellanic Clouds.  The procedure verifies that snapshot combination,
background region assignment, spectral extraction, and broadband image
context all work as expected.

The main workflow has five steps:

1. :ref:`mc-step1` — combine dithered exposures into per-source RSS snapshots
2. :ref:`mc-step2` — pair each source with an annular background region
3. :ref:`mc-step3` — extract background-subtracted spectra and overlay fiber maps
4. :ref:`mc-step4` — inspect extracted spectra with overview plots
5. :ref:`mc-step5` — fit emission lines in the extracted spectra

An independent preparation step, :ref:`mc-prep`, creates per-source
broadband image cutouts from survey mosaics (e.g. MCELS).  This can be
done at any time and does not depend on the main workflow steps.


Starting Point: the Input Files
---------------------------------

Two files drive the main workflow.

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
    source.  Additional columns include:

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


.. _mc-prep:

Preparation — Broadband Image Cutouts (``LSnap.py``)
-----------------------------------------------------

This step is independent of the main workflow and can be run at any time.
``LSnap.py`` reads a large broadband survey mosaic (e.g. an MCELS H-alpha
or [SII] image) and produces per-source cutout images and FITS files using
the source catalog.  The cutout FITS files in ``xdata/`` are later used by
step 3 to overlay LVM fiber assignments on the broadband images.

**Commands**::

    LSnap.py -size 10 -type ha  mcels_ha.fits  smc_snr_cotton24.txt
    LSnap.py -size 10 -type sii mcels_sii.fits smc_snr_cotton24.txt

The ``-size`` value sets the cutout side length in arcminutes.  ``-type``
labels the output filenames.

**Outputs**

- ``ximage/<source_name>.ha.png`` — H-alpha cutout plot with source
  regions overlaid
- ``ximage/<source_name>.sii.png`` — [SII] cutout plot with source
  regions overlaid
- ``xdata/<source_name>_<image_file>`` — FITS cutout files used in step 3

**Verification**

- Confirm that ``ximage/`` contains ``.ha.png`` and ``.sii.png`` for
  each source.
- Check that the source region ellipse is centered correctly on the
  known SNR position in each image.

See :doc:`api/LSnap/index` for API documentation.


.. _mc-step1:

Step 1 — Create Spectral Snapshots (``rss_snap.py``)
------------------------------------------------------

``rss_snap.py`` reads the observed file, groups rows by ``Source_name``,
combines all exposures for each source, runs Gaussian emission-line
fitting, and produces a 4-panel diagnostic plot.

**Command**::

    rss_snap.py -all smc_snr_observed.txt

To reprocess existing outputs add ``-redo``; to use median combination
add ``-med``.  To process a single source::

    rss_snap.py smc_snr_observed.txt J0056-7209

**Outputs**

- ``Snap/<source_name>.ave.fits`` — combined RSS FITS snapshot
- ``Snap/<source_name>.ave.tab`` — ASCII fiber position table
- ``Snap_gauss/<source_name>.gauss.txt`` — Gaussian fit parameters per fiber
- ``Snap_fig/<source_name>.png`` — 4-panel H-alpha / [SII] / ratio / FWHM map

**Verification**

- Confirm that ``Snap/`` contains one ``.ave.fits`` per source.
- Check the PNG files in ``Snap_fig/`` for coherent spatial structure.
- Inspect ``Snap_gauss/<source_name>.gauss.txt`` and verify that H-alpha
  fluxes and FWHM values are physically reasonable.

See :doc:`snapshots` for full documentation.


.. _mc-step2:

Step 2 — Generate Annular Background Regions (``GenAnnularBackground.py``)
---------------------------------------------------------------------------

``GenAnnularBackground.py`` reads the source catalog (one row per source)
and creates a paired background annulus for each entry.  Use the source
catalog here, not the observed file, to avoid duplicate entries.

**Command**::

    GenAnnularBackground.py smc_snr_cotton24.txt

**Output**

- ``smc_snr_cotton24.ann_reg.txt`` — interleaved table with one source
  row (``SourceBack='Source'``) and one background annulus row
  (``SourceBack='Back'``) per source.

The annulus geometry enforces a minimum source size of 35 arcsec, places
the annulus inner radius at ``source_major + 35"`` and the outer radius
at ``2 × source_major + 35"``.

**Verification**

- Confirm each source row is immediately followed by a background annulus
  row with the same ``Source_name``.
- Verify background rows have ``Color='green'`` and ``RegType='annulus'``.

See :doc:`region_files` for full documentation.


.. _mc-step3:

Step 3 — Extract Spectra and Overlay Fiber Maps (``GetRegSpec.py -all``)
------------------------------------------------------------------------

``GetRegSpec.py -all`` processes every source in the region table.  For
each source it:

1. Calls ``MakeLVMReg.do_complex()`` to read the snapshot and assign fiber
   colors: **red** for source fibers, **green** for background annulus
   fibers, **yellow** for everything else.  The ``.reg`` file is written
   to disk for inspection in DS9.
2. Calls ``do_one()`` to average the red fibers, subtract the median
   green fibers, and write the spectrum.
3. If ``-imgdir`` is provided, calls ``LSnap.make_one_image()`` for every
   FITS file in that directory whose name begins with the source name.
   This overlays the fiber color assignments on the broadband cutouts
   produced in the preparation step, writing output plots to ``zimage/``.

**Command without broadband overlay**::

    GetRegSpec.py -all smc_snr_cotton24.ann_reg.txt

**Command with broadband fiber overlay** (requires preparation step)::

    GetRegSpec.py -all -imgdir xdata smc_snr_cotton24.ann_reg.txt

**Outputs**

- ``FiberReg/<source_name>.ave.source.<source_name>.reg`` — DS9 region
  file with fiber color assignments (source=red, background=green,
  other=yellow)
- ``Snap_spec/Spec_<source_name>.ave_ave_back.txt`` — background-subtracted
  spectrum
- ``zimage/<image_file>.<reg_root>.png`` — broadband cutout with LVM fiber
  overlay (only when ``-imgdir`` is used)

Output spectrum columns:

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

- Plot ``FLUX`` vs ``WAVE`` and confirm emission lines (H-alpha,
  [NII] 6548/6583, [SII] 6717/6731) at the expected velocity.
- Compare ``SOURCE_FLUX`` and ``BACK_FLUX`` to confirm the background
  subtraction removed continuum without over-subtracting line emission.
- Load the ``.reg`` file from ``FiberReg/`` in DS9 to confirm red fibers
  lie on the source and green fibers form a surrounding ring.
- Check the ``zimage/`` overlay plots to confirm the fiber assignments
  are consistent with the known SNR morphology in the broadband images.

See :doc:`api/GetRegSpec/index` for API documentation.


.. _mc-step4:

Step 4 — Inspect Spectra (``PlotSpec.py`` and ``SNRPlot.py``)
--------------------------------------------------------------

Two complementary tools are run after step 3 to visually verify the
extracted spectra before committing to full Gaussian fitting.

**PlotSpec.py** — full-wavelength panel view
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``PlotSpec.py`` divides the full LVM wavelength range (3600–9560 Å) into
panels with common emission lines labelled.  When ``BACK_FLUX`` is present
(as it always is in spectra produced by step 3), the background is overlaid
in black so source and background levels can be compared in every panel.

**Command**::

    PlotSpec.py -med Snap_spec/*back*txt

The ``-med`` flag scales each panel around the local median flux, which works
well for SNR spectra that are faint in the continuum but bright in lines.
The default half-range is ``3e-15`` erg/s/cm²/Å; adjust with ``-delta`` if
needed::

    PlotSpec.py -med -delta 1e-15 Snap_spec/*back*txt

**Outputs**

- ``Overview_Plot/<spectrum_name>.overview.png`` — one file per input spectrum

**SNRPlot.py** — zoomed line panels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SNRPlot.py`` produces a complementary layout: a top panel with the full
spectrum and six sub-panels zoomed to ±1200 km/s around the key diagnostic
lines ([OII] 3728, [OIII] 5007, [OI] 6300, Hα 6563, [SII] 6720,
[SIII] 9531).  ``FLUX`` and ``ERROR`` are overlaid in each panel.  Use
``-smc`` or ``-lmc`` to set the systemic velocity offset.

**Command**::

    SNRPlot.py -smc Snap_spec/*back*txt

**Outputs**

- ``SNRPlot_fig/<spectrum_name>.png`` — one file per input spectrum
- ``SNRPlot_fig/README.txt`` — records the date, command, and velocity
  offset used to generate the plots

**Verification**

- In the ``PlotSpec.py`` panels, confirm that H-alpha, [NII] 6548/6583,
  and [SII] 6717/6731 are positive peaks in the ``FLUX`` (blue) trace and
  that the ``BACK_FLUX`` (black) trace is smooth with no sharp emission
  features.
- In the ``SNRPlot.py`` panels, check that the line centres are at the
  expected wavelength for the systemic velocity and that the ERROR is
  well below the line flux.
- If continuum levels look unreasonably high or negative, re-examine the
  fiber assignments from step 3.


.. _mc-step5:

Step 5 — Fit Emission Lines (``lvm_gaussfit.py``)
--------------------------------------------------

``lvm_gaussfit.py`` fits Gaussian profiles to the prominent emission lines
in the extracted spectra.  The spectra produced by step 3 each contain
three spectral columns: the background-subtracted flux (``FLUX``), the raw
source flux before subtraction (``SOURCE_FLUX``), and the background flux
(``BACK_FLUX``).  The ``-stype`` option selects which column to fit;
omitting it fits the background-subtracted spectrum.

All three variants are normally run to allow comparison of line fluxes
across the source, background, and net spectra.  The ``-smc`` flag sets
the systemic velocity to 146 km/s (SMC); use ``-lmc`` for LMC targets
(262 km/s).

**Commands**::

    lvm_gaussfit.py -smc -stype BACK   Snap_spec/*back*txt
    lvm_gaussfit.py -smc -stype SOURCE Snap_spec/*back*txt
    lvm_gaussfit.py -smc               Snap_spec/*back*txt

**What each command fits**

===================  ========================================================
Command              Spectrum fitted
===================  ========================================================
``-stype BACK``      Background annulus (``BACK_FLUX`` / ``BACK_ERROR``)
``-stype SOURCE``    Raw source before subtraction (``SOURCE_FLUX`` / ``SOURCE_ERROR``)
no ``-stype``        Background-subtracted source (``FLUX`` / ``ERROR``)
===================  ========================================================

**Outputs**

Output files are named ``Gauss_<date>.<stype>.txt`` for the ``BACK`` and
``SOURCE`` runs and ``Gauss_<date>.txt`` for the background-subtracted run,
where ``<date>`` is today's date in YYMMDD format.  To use a fixed name
rather than a date, add ``-out <name>``::

    lvm_gaussfit.py -smc -stype BACK   -out smc_snr Snap_spec/*back*txt
    lvm_gaussfit.py -smc -stype SOURCE -out smc_snr Snap_spec/*back*txt
    lvm_gaussfit.py -smc               -out smc_snr Snap_spec/*back*txt

This produces ``Gauss_smc_snr.BACK.txt``, ``Gauss_smc_snr.SOURCE.txt``,
and ``Gauss_smc_snr.txt``.

**Verification**

- Open the output files and confirm that H-alpha fluxes are positive and
  physically plausible for SNRs (line widths, velocity offsets).
- Compare ``SOURCE`` and ``BACK`` fits to confirm that the background
  contribution is being properly characterised.
- Check that the background-subtracted fits (no ``-stype``) show enhanced
  [SII]/H-alpha ratios characteristic of shock-ionised SNR gas.

See :doc:`api/lvm_gaussfit/index` for API documentation.




Complete Worked Example
-----------------------

The preparation step and main workflow can be run in either order::

    # Preparation (independent): create per-source broadband image cutouts
    LSnap.py -size 10 -type ha  mcels_ha.fits  smc_snr_cotton24.txt
    LSnap.py -size 10 -type sii mcels_sii.fits smc_snr_cotton24.txt

    # Step 1: create snapshots, fit lines, make diagnostic plots
    rss_snap.py -all smc_snr_observed.txt

    # Step 2: generate source + background annulus table
    GenAnnularBackground.py smc_snr_cotton24.txt

    # Step 3: extract spectra and overlay fiber maps on broadband images
    GetRegSpec.py -all -imgdir xdata smc_snr_cotton24.ann_reg.txt

    # Step 4: inspect extracted spectra
    PlotSpec.py -med Snap_spec/*back*txt
    SNRPlot.py -smc  Snap_spec/*back*txt

    # Step 5: fit emission lines in all three spectral variants
    lvm_gaussfit.py -smc -stype BACK   -out smc_snr Snap_spec/*back*txt
    lvm_gaussfit.py -smc -stype SOURCE -out smc_snr Snap_spec/*back*txt
    lvm_gaussfit.py -smc               -out smc_snr Snap_spec/*back*txt


Output Summary
--------------

After a complete run the working directory should contain the following
subdirectories.  All are created automatically by the relevant script.
File counts are per source unless otherwise noted.

``Snap/`` — RSS snapshots (Step 1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Created by ``rss_snap.py``.

- ``<source>.ave.fits`` — combined RSS FITS file; one per source.
  Extensions: ``FLUX``, ``ERROR``, ``MASK``, ``SKY``, ``LSF``, and a
  fibre-position table.  This is the primary input to steps 3 onward.
- ``<source>.ave.tab`` — ASCII fibre position table (RA, Dec, fibre index).

``Snap_gauss/`` — per-fibre Gaussian fits (Step 1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Created by ``rss_snap.py``.

- ``<source>.gauss.txt`` — ASCII table of Gaussian fit parameters
  (centre, amplitude, FWHM) for H-alpha and [SII] in every fibre.
  Used to generate the Snap_fig diagnostic maps.

``Snap_fig/`` — snapshot diagnostic plots (Step 1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Created by ``rss_snap.py``.

- ``<source>.png`` — 4-panel image: H-alpha flux map, [SII] flux map,
  [SII]/H-alpha ratio map, and H-alpha FWHM map.  Coloured circles show
  individual fibre values.  Used to confirm that the snapshot combination
  has captured coherent spatial structure.

``ximage/`` — broadband cutout plots (Preparation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Created by ``LSnap.py -size``.

- ``<source>.ha.png`` — H-alpha image cutout with source region overlaid.
- ``<source>.sii.png`` — [SII] image cutout with source region overlaid.

Inspect these to verify that the source ellipse is correctly centred on
the known SNR position.

``xdata/`` — broadband FITS cutouts (Preparation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Created by ``LSnap.py -size``.

- ``<source>_<image_file>`` — FITS cutout used by step 3 to overlay LVM
  fibre assignments on the broadband image.  One file per source per
  input mosaic (e.g. one for H-alpha, one for [SII]).

``FiberReg/`` — DS9 region files (Step 3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Created by ``GetRegSpec.py -all`` (via ``MakeLVMReg.do_complex()``).

- ``<source>.ave.source.<source>.reg`` — DS9 region file assigning each
  LVM fibre a colour: **red** for source fibres inside the SNR ellipse,
  **green** for background annulus fibres, **yellow** for all others.

Load these in DS9 alongside the snapshot FITS file to confirm that the
red fibres cover the SNR shell and the green fibres form a clean
surrounding ring free of source emission.

``Snap_spec/`` — extracted spectra (Step 3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Created by ``GetRegSpec.py -all``.

- ``Spec_<source>.ave_ave_back.txt`` — background-subtracted spectrum;
  one per source.  Columns:

  ================  ==============================================
  Column            Contents
  ================  ==============================================
  WAVE              Wavelength in Angstroms
  FLUX              Background-subtracted flux (source − background)
  ERROR             Combined source and background uncertainty
  SOURCE_FLUX       Mean spectrum of red (source) fibres
  SOURCE_ERROR      Uncertainty on source spectrum
  BACK_FLUX         Median spectrum of green (background) fibres
  BACK_ERROR        Uncertainty on background spectrum
  SKY               Sky spectrum (if present in SFrame)
  SKY_ERROR         Sky uncertainty
  MASK              Sum of mask flags across selected fibres
  LSF               Mean line spread function
  ================  ==============================================

Use ``FLUX`` for science measurements.  ``SOURCE_FLUX`` and ``BACK_FLUX``
are retained so the subtraction can be inspected and, if necessary, redone
with different weights.

``zimage/`` — broadband fibre overlay plots (Step 3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Created by ``GetRegSpec.py -all -imgdir xdata``.

- ``<image_file>.<reg_root>.png`` — broadband cutout with the DS9 fibre
  colour assignments overlaid.  One file per source per input broadband
  image.  Confirms that the fibre selection matches the known SNR
  morphology.

``Overview_Plot/`` — full-wavelength overview plots (Step 4)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Created by ``PlotSpec.py``.

- ``<spectrum_name>.overview.png`` — full-wavelength plot (3600–9560 Å)
  divided into 750 Å panels, each autoscaled to the local median.
  Common emission lines are labelled in red.  Where ``BACK_FLUX`` is
  present it is overlaid in black.

Check that ``FLUX`` shows positive peaks at H-alpha, [NII] 6548/6583,
and [SII] 6717/6731, and that ``BACK_FLUX`` is smooth with no sharp line
features.

``SNRPlot_fig/`` — zoomed line-panel plots (Step 4)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Created by ``SNRPlot.py``.

- ``<spectrum_name>.png`` — one plot per spectrum.  Top panel: full
  spectrum with ``FLUX`` and ``ERROR``.  Lower six panels: ±1200 km/s
  zoom around [OII] 3728, [OIII] 5007, [OI] 6300, Hα 6563, [SII] 6720,
  and [SIII] 9531, each with ``FLUX`` and ``ERROR``.
- ``README.txt`` — records the date, command line, and velocity offset
  used to generate the plots.

Check that line centres fall at the expected wavelength for the systemic
velocity and that ``ERROR`` is well below the line peak flux.

Gaussian fit tables — emission-line measurements (Step 5)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Written to the **current working directory** by ``lvm_gaussfit.py``.

- ``Gauss_<name>.BACK.txt`` — Gaussian fits to the background spectrum
  (``BACK_FLUX``).
- ``Gauss_<name>.SOURCE.txt`` — Gaussian fits to the raw source spectrum
  (``SOURCE_FLUX``).
- ``Gauss_<name>.txt`` — Gaussian fits to the background-subtracted
  spectrum (``FLUX``).

Each file has one row per source.  Check that H-alpha fluxes are positive,
line widths are physically reasonable (FWHM ≲ a few hundred km/s for
thermal/shock-ionised gas), and that ``Gauss_<name>.txt`` shows elevated
[SII]/H-alpha ratios (> 0.4) characteristic of shock-ionised SNR gas.


See Also
--------

- :doc:`snapshots` — detailed documentation for ``rss_snap.py``
- :doc:`region_files` — detailed documentation for ``GenAnnularBackground.py``
- :doc:`api/LSnap/index` — API documentation for ``LSnap.py``
- :doc:`api/GetRegSpec/index` — API documentation for ``GetRegSpec.py``
- :doc:`api/MakeLVMReg/index` — API documentation for ``MakeLVMReg.py``
- :doc:`api/PlotSpec/index` — API documentation for ``PlotSpec.py``
- :doc:`api/SNRPlot/index` — API documentation for ``SNRPlot.py``
- :doc:`api/lvm_gaussfit/index` — API documentation for ``lvm_gaussfit.py``
