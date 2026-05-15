Summarizing LVM Data
====================

The lvm_ksl package includes several scripts for summarizing LVM data
at different stages of processing. These scripts are useful for:

- Creating catalogs of available observations
- Evaluating data quality across many exposures
- Comparing sky levels and sky subtraction residuals over time
- Identifying trends or problems in the data

Five summarization scripts are provided, each operating on different
data products:

- ``SummarizeData.py`` - Catalog raw/unprocessed data files
- ``SummarizeCframe.py`` - Summarize CFrame spectra (before sky subtraction)
- ``SummarizeSframe.py`` - Summarize SFrame spectra (after sky subtraction)
- ``SummarizeRings.py`` - Summarize SFrame spectra by fiber ring position
- ``SummarizeSpec.py`` - Summarize CFrame or SFrame spectra by spectrograph


SummarizeData.py - Cataloging Raw Data
--------------------------------------

This script reads FITS headers from raw LVM data files and creates an
astropy table with metadata for each exposure. This is useful for
identifying all observations of a particular object or sky position.

**Command line usage**::

    SummarizeData.py [-h] [-all] [-redo] [-range] [mjd1 mjd2 ...]

**Options:**

-h
    Print help and exit.

-all
    Summarize all data in subdirectories.

-redo
    Recreate summary tables (only with -all).

-range
    Process all MJDs between mjd1 and mjd2 (inclusive).

mjd1 mjd2 ...
    Specific MJD values to process.

**Output:**

Individual summary tables are created for each MJD, then stacked into
a combined file ``All_data.txt``.

**Notes:**

This script should be run from the ``lvm_data/sas/sdsswork/data/lvm/lco``
directory (or equivalent data location).


SummarizeCframe.py - CFrame Spectral Summary
--------------------------------------------

This script processes CFrame files (wavelength-calibrated, flux-calibrated,
but before sky subtraction) and computes the median spectrum across all
science fibers for each exposure. This creates a compact representation
of how the raw spectra (including sky) vary over time.

**Command line usage**::

    SummarizeCframe.py [-h] [-out file_out] [-emin 900] [-ver drp_ver] exp_start exp_stop delta

**Options:**

-h
    Print help and exit.

-out file_out
    Set output filename (default: auto-generated).

-emin N
    Minimum exposure time in seconds to include (default: 900).

-ver drp_ver
    DRP version to use (default: 1.2.0).

**Arguments:**

exp_start
    Starting exposure number.

exp_stop
    Ending exposure number.

delta
    Process every Nth exposure (use 1 for all).

**Output:**

A FITS file with the following extensions:

=========  =============  =====================================================
Extension  Type           Contents
=========  =============  =====================================================
PRIMARY    ImageHDU       No data; header records Title
WAVE       ImageHDU       1-D wavelength array
FLUX       ImageHDU       2-D array (n_exposures × n_wavelengths), median science fiber flux
SKY_EAST   ImageHDU       2-D array, SKY_EAST telescope spectrum
SKY_WEST   ImageHDU       2-D array, SKY_WEST telescope spectrum
LSF        ImageHDU       2-D array, line spread function
drp_all    BinTableHDU    drpall metadata for the included exposures
=========  =============  =====================================================

The default output filename follows the pattern
``XCframe_<ver>_<exp_start>_<exp_stop>_<delta>_<percent>.fits``,
e.g. ``XCframe_1.2.0_10000_20000_10_50.fits``.

**Use cases:**

- Monitoring sky brightness variations over a night or survey
- Identifying exposures with unusual sky conditions
- Comparing sky telescope spectra with science fiber sky


SummarizeSframe.py - SFrame Spectral Summary
--------------------------------------------

This script processes SFrame files (after sky subtraction) and computes
the median (or other percentile) spectrum across all science fibers for
each exposure. This is useful for evaluating sky subtraction quality.

**Command line usage**::

    SummarizeSframe.py [-h] [-ver drp_ver] [-percent 50] [-emin 900] [-out whatever] exp_start exp_stop delta

**Options:**

-h
    Print help and exit.

-ver drp_ver
    DRP version to use (default: 1.2.0).

-percent N
    Percentile to compute (default: 50 = median).

-emin N
    Minimum exposure time in seconds to include (default: 900).

-out name
    Output filename root.

**Arguments:**

exp_start
    Starting exposure number.

exp_stop
    Ending exposure number.

delta
    Process every Nth exposure.

**Output:**

A FITS file with the following extensions:

=========  =============  =====================================================
Extension  Type           Contents
=========  =============  =====================================================
PRIMARY    ImageHDU       No data; header records Title
WAVE       ImageHDU       1-D wavelength array
FLUX       ImageHDU       2-D array (n_exposures × n_wavelengths), percentile sky-subtracted flux
SKY        ImageHDU       2-D array, sky spectrum
IVAR       ImageHDU       2-D array, inverse variance
LSF        ImageHDU       2-D array, line spread function
drp_all    BinTableHDU    drpall metadata for the included exposures
=========  =============  =====================================================

The default output filename follows the pattern
``XSFrame_<ver>_<exp_start>_<exp_stop>_<delta>_<percent>.fits``,
e.g. ``XSFrame_1.2.0_10000_20000_10_50.fits``.

**Use cases:**

- Evaluating sky subtraction residuals across many exposures
- Identifying systematic patterns in sky subtraction
- Comparing different DRP versions


SummarizeRings.py - Ring-Based Spectral Summary
-----------------------------------------------

This script is similar to SummarizeSframe.py but computes separate
summaries for three radial ring sets (inner, middle, outer) based on
the fiber ring number in the IFU. This is useful for detecting radial
variations in sky subtraction quality.

**Command line usage**::

    SummarizeRings.py [-h] [-ver drp_ver] [-percent 50] [-emin 900] [-out whatever]
                      [-inner 1 9] [-middle 10 19] [-outer 20 25]
                      exp_start exp_stop delta

**Options:**

-h
    Print help and exit.

-ver drp_ver
    DRP version to use (default: 1.2.0).

-percent N
    Percentile to compute (default: 50 = median).

-emin N
    Minimum exposure time in seconds to include (default: 900).

-out name
    Output filename root.

-inner min max
    Ring number range for inner set (default: 1-9).

-middle min max
    Ring number range for middle set (default: 10-19).

-outer min max
    Ring number range for outer set (default: 20-25).

**Arguments:**

exp_start
    Starting exposure number.

exp_stop
    Ending exposure number.

delta
    Process every Nth exposure.

**Output:**

A FITS file with the following extensions:

===========  =============  =====================================================
Extension    Type           Contents
===========  =============  =====================================================
PRIMARY      ImageHDU       No data; header records PERCENT, INNER, MIDDLE, OUTER
WAVE         ImageHDU       1-D wavelength array
FLUX_INNER   ImageHDU       2-D array (n_exposures × n_wavelengths), percentile flux for inner ring fibers
FLUX_MIDDLE  ImageHDU       2-D array, percentile flux for middle ring fibers
FLUX_OUTER   ImageHDU       2-D array, percentile flux for outer ring fibers
SKY          ImageHDU       2-D array, median sky spectrum
drp_all      BinTableHDU    drpall metadata for the included exposures
===========  =============  =====================================================

The default output filename follows the pattern
``XRings_<ver>_<exp_start>_<exp_stop>_<delta>_<percent>.fits``,
e.g. ``XRings_1.2.0_10000_20000_10_50.fits``.

**Use cases:**

- Detecting radial gradients in sky subtraction residuals
- Evaluating whether sky subtraction performs differently at
  different distances from the IFU center
- Diagnosing fiber-dependent calibration issues


SummarizeSpec.py - Per-Spectrograph Spectral Summary
-----------------------------------------------------

This script processes CFrame files (the default) or SFrame files and
computes the percentile spectrum separately for the science fibers
belonging to each of the three LVM spectrographs.  Only fibers with
``telescope == 'Sci'`` and ``fibstatus == 0`` are used, and the
spectrograph is identified from the ``spectrographid`` column of the
SLITMAP extension.

An exposure is included in the output only if all three spectrographs
have valid science fibers.  Exposures missing one or more spectrographs
are skipped and recorded in a companion ASCII table.

**Command line usage**::

    SummarizeSpec.py [-sf] [-ver drp_ver] [-percent 50] [-emin 900] [-out name]
                     exp_start exp_stop delta

**Options:**

-h
    Print help and exit.

-sf
    Read SFrame (sky-subtracted) files instead of the default CFrame
    (flux-calibrated, before sky subtraction) files.

-ver drp_ver
    DRP version to use (default: 1.2.0).

-percent N
    Percentile to compute across fibers (default: 50 = median).

-emin N
    Minimum exposure time in seconds to include (default: 900).

-out name
    Output filename root.  If omitted, a name is generated automatically
    from the file type, DRP version, and exposure range.

**Arguments:**

exp_start
    Starting exposure number.

exp_stop
    Ending exposure number.

delta
    Process every Nth exposure (use 1 for all).

**Output files:**

A FITS file with the following extensions:

=========  =============  =====================================================
Extension  Type           Contents
=========  =============  =====================================================
PRIMARY    ImageHDU       No data; header records PERCENT, SP1/SP2/SP3, Title
WAVE       ImageHDU       1-D wavelength array
FLUX1      ImageHDU       2-D array (n_exposures × n_wavelengths), spectrograph 1
FLUX2      ImageHDU       2-D array (n_exposures × n_wavelengths), spectrograph 2
FLUX3      ImageHDU       2-D array (n_exposures × n_wavelengths), spectrograph 3
drp_all    BinTableHDU    drpall metadata for the accepted exposures only
=========  =============  =====================================================

A fixed-width ASCII table ``<out>.skipped.txt`` is also written when any
exposures are rejected, listing each one with columns ``expnum``, ``mjd``,
``tileid``, ``SP1``, ``SP2``, ``SP3`` (``True`` = present, ``False`` =
absent).

**Notes:**

The drpall ``location`` column records SFrame paths.  For CFrame files the
script replaces ``'SFrame'`` with ``'CFrame'`` in the path automatically.

The default output filename follows the pattern
``XSpec_<type>_<ver>_<exp_start>_<exp_stop>_<delta>_<percent>.fits``,
where ``<type>`` is ``CFrame`` or ``SFrame``, e.g.
``XSpec_CFrame_1.2.0_10000_20000_10_50.fits``.  The skipped-exposure file
takes the same root with a ``.skipped.txt`` suffix, e.g.
``XSpec_CFrame_1.2.0_10000_20000_10_50.skipped.txt``.  Using the file type
in the name means CFrame and SFrame runs do not overwrite each other.

**Use cases:**

- Detecting spectrograph-to-spectrograph offsets in flux calibration
- Identifying whether sky subtraction residuals are confined to one
  spectrograph
- Comparing CFrame and SFrame runs to isolate sky subtraction artefacts
  by spectrograph


Typical Workflow
----------------

A typical workflow for evaluating data quality might be:

1. **Catalog the data**::

       SummarizeData.py -all

   This creates ``All_data.txt`` with metadata for all exposures.

2. **Summarize CFrame spectra** (to see raw sky levels)::

       SummarizeCframe.py -out cframe_summary 10000 20000 10

   This samples every 10th exposure from 10000 to 20000.

3. **Summarize SFrame spectra** (to evaluate sky subtraction)::

       SummarizeSframe.py -out sframe_summary 10000 20000 10

4. **Check for radial variations**::

       SummarizeRings.py -out rings_summary 10000 20000 10

5. **Check for spectrograph-to-spectrograph variations**::

       SummarizeSpec.py -out spec_summary 10000 20000 10

   Add ``-sf`` to compare the same exposures after sky subtraction::

       SummarizeSpec.py -sf -out spec_summary_sf 10000 20000 10


gauss_offset.py — Airglow Line Fitting
---------------------------------------

Fits single Gaussians to a fixed set of airglow lines in each exposure row
of a SummarizeCframe FITS file, producing one output row per exposure.
This is useful for tracking wavelength-calibration drifts or changes in sky
brightness over a night or survey.

**Command line usage**::

    gauss_offset.py [-ext FLUX] [-out root] [-ver drp_ver] [-file xfile]
                    [exp_start [exp_stop [delta]]]
                    filename [filename ...]

**Options:**

-ext FLUX|SKY_EAST|SKY_WEST
    Extension to fit (default: FLUX).

-out root
    Root name for the output file.  If omitted, the output filename is
    derived automatically (see Output below).

-ver drp_ver
    DRP version for drpall lookup when resolving file paths (default 1.2.1).

-file xfile
    Read lvmCFrame filenames from a table that has a column named
    ``filename`` or ``Filename``.  SFrame paths are converted to CFrame
    automatically.  These are combined with any files from
    ``exp_start``/``exp_stop``.

**Arguments:**

exp_start
    First exposure number.  Combined with ``exp_stop``, queries drpall
    to build the list of CFrame files to process.

exp_stop
    Last exposure number (inclusive).

delta
    Use every Nth exposure in [exp_start, exp_stop] (default 1 = all).

filename
    One or more SummarizeCframe FITS files (``XCframe_*.fits``).

**Two input modes:**

*CFrame mode* (``-file`` or ``exp_start``/``exp_stop``): each CFrame file
is opened individually, the median spectrum is computed across all good
science fibers (``telescope == 'Sci'``, ``fibstatus == 0``) using the
MASK extension, and Gaussians are fit to that median.  SFrame paths are
converted to CFrame automatically; drpall is used as a fallback if the
literal path does not exist.

*SummarizeCframe mode* (direct ``filename`` arguments): each row of the
chosen extension in an ``XCframe_*.fits`` file is treated as a
pre-computed median spectrum and fit directly.

**Airglow lines fitted:**

=========  ===========
Line name  Wavelength
=========  ===========
sky5577    5577.34 A
sky6300    6300.31 A
sky6363    6363.78 A
sky7358    7358.68 A
sky7392    7392.21 A
sky7914    7913.72 A
sky8344    8344.61 A
sky8399    8399.18 A
sky8827    8827.11 A
sky8988    8988.38 A
sky9552    9552.55 A
sky9719    9719.84 A
=========  ===========

Wavelengths are from the ESO UVES sky spectrum atlas.

**Output:**

A FITS binary table with one row per exposure.  Columns cover the fit
parameters (flux, wave, fwhm, back, rmse) for each line, plus ``expnum``
and ``mjd``.

- CFrame mode: ``Gauss_<ext>.<YYMMDD>.fits`` by default.
- SummarizeCframe mode: ``<input_stem>_gauss.<ext>.fits`` by default.

Supplying ``-out root`` overrides both defaults.

**Examples**::

    # CFrame mode: read file list from a drpall-derived table
    gauss_offset.py -file drpall_dr20.fits

    # CFrame mode: use drpall to select exposures 11000–12000, every 5th
    gauss_offset.py -ver 1.2.0 11000 12000 5

    # CFrame mode: fit the east sky telescope
    gauss_offset.py -ext SKY_EAST -file drpall_dr20.fits

    # SummarizeCframe mode: fit a pre-built summary file
    gauss_offset.py XCframe_1.2.0_10000_20000_10_50.fits

    # SummarizeCframe mode: fit sky east extension
    gauss_offset.py -ext SKY_EAST XCframe_1.2.0_10000_20000_10_50.fits


See Also
--------

- :doc:`api/SummarizeData/index` - API documentation
- :doc:`api/SummarizeCframe/index` - API documentation
- :doc:`api/SummarizeSframe/index` - API documentation
- :doc:`api/SummarizeRings/index` - API documentation
- :doc:`api/SummarizeSpec/index` - API documentation
