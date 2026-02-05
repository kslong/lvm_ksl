Creating Snapshots from Multiple Exposures
==========================================

The ``rss_snap.py`` script provides an automated workflow for creating
stacked RSS spectral files centered on specific astronomical sources.
This is particularly useful for analyzing extended objects like supernova
remnants (SNRs), planetary nebulae, or HII regions that have been observed
in multiple LVM exposures.  This routine not only creates row-stacked-spectra
of the shanpshots, but also carries out simple gaussian fits to some
of the strongs lines.  The RSS spectra are suitable for more detailed
fitting with the DAP. Though not strictly limited to this it is 
primarily intended for creating RSS spectra of dithered exposures,
such as those that have been obtained for the LMC and SMC, where the region
of interest is relatively small..

Overview
--------

When studying individual sources in LVM data, you often want to:

1. Combine multiple exposures that overlap a source position
2. Extract a small region centered on the source
3. Fit emission lines across all fibers in that region
4. Create diagnostic maps of the source

The ``rss_snap.py`` script automates this entire workflow, making it easy
to process catalogs of sources in batch mode.


Input File Format
-----------------

The primary input is an ASCII table (xfile) that associates source names
with their LVM exposures. This file is typically created by cross-matching
a source catalog with the LVM observation database, that is to say
the drp\_all file. 

Required columns:

============  ===========================================================
Column        Description
============  ===========================================================
Source_name   Unique identifier for each source (e.g., "SNR_N49")
RA            Right Ascension of the source center in degrees
Dec           Declination of the source center in degrees
location      Path to the SFrame FITS file for each exposure
============  ===========================================================

Each row represents one exposure for a given source. Sources observed
multiple times will have multiple rows with the same Source_name, RA,
and Dec but different location values.

Example input file::

    Source_name     RA          Dec         location
    SNR_N49         81.501085   -66.082249  dr1/spectro/lvm/lvm1/0012/12345/lvmSFrame-00012345.fits
    SNR_N49         81.501085   -66.082249  dr1/spectro/lvm/lvm1/0012/12346/lvmSFrame-00012346.fits
    SNR_N49         81.501085   -66.082249  dr1/spectro/lvm/lvm1/0012/12400/lvmSFrame-00012400.fits
    SNR_N49B        81.341994   -65.996487  dr1/spectro/lvm/lvm1/0012/12345/lvmSFrame-00012345.fits
    SNR_N49B        81.341994   -65.996487  dr1/spectro/lvm/lvm1/0012/12350/lvmSFrame-00012350.fits


The input file is read using the routine astropy.io.ascii.read, and may contain addtional columns.

Processing Steps
----------------

For each source, ``rss_snap.py`` performs the following steps:

1. **File Location**

   Reads the input table and locates all SFrame FITS files associated
   with the requested source. The script automatically detects whether
   it is running at Utah/CHPC, or on a local machine with data mirrors.

2. **RSS Combination**

   Calls ``rss_combine_pos.do_fixed()`` to combine the exposures:

   - Creates a WCS centered on the source RA/Dec
   - Generates output fiber positions on a regular 35 arcsec grid
   - Apportions flux from input fibers to output virtual fibers based
     on circular aperture overlap
   - Combines multiple exposures using mean (default) or median (-med)
   - Writes the combined RSS FITS file

3. **Velocity Determination**

   Looks up the expected radial velocity at the source position using
   a velocity model (based on position in LMC, SMC, or Galaxy).

4. **Emission Line Fitting**

   Runs ``lvm_gaussfit.do_all()`` to fit Gaussian profiles to emission
   lines (H-alpha, [NII], [SII], etc.) for every fiber in the combined
   RSS file.

5. **Visualization**

   Creates a 4-panel diagnostic plot showing spatial maps of:

   - H-alpha flux
   - [SII] flux (sum of 6717 and 6731 lines)
   - [SII]/H-alpha ratio
   - H-alpha FWHM


Output Files and Directories
----------------------------

The script creates three output directories:

**Snap/**
    Combined RSS FITS files and associated tables.

    - ``<source_name>.ave.fits`` - Combined RSS using mean
    - ``<source_name>.med.fits`` - Combined RSS using median (if -med used)
    - ``<source_name>.ave.tab`` - ASCII table of fiber positions

**Snap_gauss/**
    Gaussian fit results for each source.

    - ``<source_name>.fits`` - FITS table with fit parameters per fiber
    - ``<source_name>.tab`` - ASCII version of fit results

**Snap_fig/**
    Diagnostic plots.

    - ``<source_name>.png`` - 4-panel flux/ratio/FWHM map

**xtmp/** (temporary)
    Intermediate remapped files, deleted after processing unless
    the ``-keep`` option is used.


Command Line Usage
------------------

Basic syntax::

    rss_snap.py [-h] [-keep] [-redo] [-all] [-med] [-size arcmin] xfile source_name

**Arguments:**

xfile
    Input file containing source names and exposure locations.

source_name
    One or more source names to process (must match Source_name column).

**Options:**

-h
    Print help message.

-keep
    Retain temporary files in the xtmp/ directory.

-redo
    Recreate output files even if they already exist.

-all
    Process all unique sources in xfile. This is the typical mode for
    batch processing a source catalog.

-med
    Use median instead of mean when combining exposures. The median is
    more robust to outliers and cosmic rays but may reduce S/N for
    faint sources.

-size arcmin
    Size of the output region in arcminutes (default: 10).


Examples
--------

Process a single source::

    rss_snap.py lmc_snr.out SNR_N49

Process multiple named sources::

    rss_snap.py lmc_snr.out SNR_N49 SNR_N49B SNR_N63A

Process all sources in a catalog::

    rss_snap.py -all lmc_snr.out

Reprocess all sources using median combination::

    rss_snap.py -all -redo -med lmc_snr.out

Process with a larger extraction region (15 arcmin)::

    rss_snap.py -all -size 15 lmc_snr.out


Notes
-----

- The output FITS files use the same format as ``rss_combine_pos``, with
  extensions: PRIMARY, FLUX, IVAR, MASK, WAVE, LSF, SLITMAP, WCS_INFO,
  EXPOSURE.

- The combination type (ave or med) is appended to the output filename,
  so you can compare mean and median results for the same source.

- When using ``-redo``, only the RSS combination step is repeated. The
  Gaussian fitting and plotting are always performed.

- The script determines whether files exist locally before attempting
  to process them. Missing files are skipped with a warning.


See Also
--------

- :doc:`api/rss_combine_pos/index` - Lower-level RSS combination for fixed positions
- :doc:`api/rss_combine/index` - General RSS combination utilities
- :doc:`api/lvm_gaussfit/index` - Emission line fitting
