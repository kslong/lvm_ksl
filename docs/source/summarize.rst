Summarizing LVM Data
====================

The lvm_ksl package includes several scripts for summarizing LVM data
at different stages of processing. These scripts are useful for:

- Creating catalogs of available observations
- Evaluating data quality across many exposures
- Comparing sky levels and sky subtraction residuals over time
- Identifying trends or problems in the data

Four summarization scripts are provided, each operating on different
data products:

- ``SummarizeData.py`` - Catalog raw/unprocessed data files
- ``SummarizeCframe.py`` - Summarize CFrame spectra (before sky subtraction)
- ``SummarizeSframe.py`` - Summarize SFrame spectra (after sky subtraction)
- ``SummarizeRings.py`` - Summarize SFrame spectra by fiber ring position


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
    DRP version to use (default: 1.1.1).

**Arguments:**

exp_start
    Starting exposure number.

exp_stop
    Ending exposure number.

delta
    Process every Nth exposure (use 1 for all).

**Output:**

A FITS file with one row per exposure, containing:

- Median science fiber spectrum
- SKY_EAST telescope spectrum
- SKY_WEST telescope spectrum
- Exposure metadata (MJD, coordinates, etc.)

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

A FITS file with one row per exposure, containing:

- Percentile sky-subtracted flux spectrum
- Sky spectrum
- Inverse variance spectrum
- Exposure metadata

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

A FITS file with one row per exposure, containing:

- Percentile spectrum for inner ring fibers
- Percentile spectrum for middle ring fibers
- Percentile spectrum for outer ring fibers
- Median sky spectrum
- Exposure metadata

**Use cases:**

- Detecting radial gradients in sky subtraction residuals
- Evaluating whether sky subtraction performs differently at
  different distances from the IFU center
- Diagnosing fiber-dependent calibration issues


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


See Also
--------

- :doc:`api/SummarizeData/index` - API documentation
- :doc:`api/SummarizeCframe/index` - API documentation
- :doc:`api/SummarizeSframe/index` - API documentation
- :doc:`api/SummarizeRings/index` - API documentation
