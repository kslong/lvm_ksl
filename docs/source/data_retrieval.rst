Data Retrieval and Reduction
============================

The lvm_ksl package provides tools for retrieving LVM data from the SDSS
data archive at Utah and for running the LVM Data Reduction Pipeline (DRP)
on local machines.


Overview
--------

The typical workflow for obtaining reduced LVM data is:

1. **Retrieve data from Utah** using ``GetFromUtah.py`` or ``LocateReduced.py``
2. **Run the DRP locally** using ``Reduce.py`` (if processing raw data)
3. **Locate processed files** using ``LocateData.py`` or ``LocateReduced.py``

These tools simplify working with the complex LVM data directory structure
and allow batch processing of multiple exposures.


Retrieving Data from Utah
-------------------------

GetFromUtah.py
^^^^^^^^^^^^^^

Retrieves DRP-processed data (SFrame or CFrame files) from the Utah
data archive to a local machine.

**Usage**::

    GetFromUtah.py [-h] [-cp] [-CFrame] [-drp VERSION] [-out OUTDIR] mjd expstart [expstop]

**Options:**

-h
    Print help and exit.

-cp
    Copy retrieved data to a local directory.

-CFrame
    Retrieve CFrame files instead of SFrame files.

-drp VERSION
    Specify DRP version (default: 1.2.0).

-out OUTDIR
    Specify output directory.

**Arguments:**

mjd
    MJD of the observations.

expstart
    First exposure number to retrieve.

expstop
    Last exposure number (optional, defaults to expstart).

**Example**::

    # Retrieve exposures 4155-4160 from MJD 60188
    GetFromUtah.py -drp 1.2.0 60188 4155 4160

LocateReduced.py
^^^^^^^^^^^^^^^^

Locates reduced data files using the drpall FITS file, which contains
metadata for all DRP-processed exposures.

**Usage**::

    LocateReduced.py [-h] [-cp] [-dir whatever] [-file_type CFrame] exp_min exp_max

**Options:**

-h
    Print help and exit.

-cp
    Copy files to a local directory.

-dir path
    Alternative directory for copied files.

-file_type
    File type to locate (default: SFrame).

**Arguments:**

exp_min, exp_max
    Range of exposure numbers to locate.

**Notes:**

This routine reads the drpall file to find file locations within the
reduced data directory structure. It works locally and does not retrieve
data from Utah.


Running the DRP Locally
-----------------------

Reduce.py
^^^^^^^^^

Downloads raw LVM data and processes it through the DRP pipeline locally.
This is useful for reprocessing data with different parameters or for
processing data that is not yet available in the Utah archive.

**Usage**::

    Reduce.py [-h] [-keep] [-cp] [-np N] exposures_to_process

**Options:**

-h
    Print help and exit.

-keep
    Retain ancillary files (normally deleted to save disk space).

-cp
    Copy reduced frames to a local ``./data`` directory.

-np N
    Number of parallel threads for processing.

**Arguments:**

The ``exposures_to_process`` argument uses a flexible format:

- Numbers > 50000 are treated as MJDs
- Numbers < 500 are exposure numbers
- Ranges like ``500-510`` process exposures 500 through 510
- Comma-separated lists like ``500,505,513`` process specific exposures

**Examples**::

    # Process exposures 4155-4157 from MJD 60188
    Reduce.py 60188 4155-4157

    # Process multiple MJDs with various exposures
    Reduce.py 60188 4155-4157 60189 4321 60190 5011,5012,5013

    # Process with 4 parallel threads
    Reduce.py -np 4 60188 4155-4160

**Output:**

- Reduced data is stored in the local LVM data structure
- Log files are written to the ``./xlog`` directory
- With ``-cp``, reduced frames are copied to ``./data``


Locating Processed Data
-----------------------

LocateData.py
^^^^^^^^^^^^^

Searches local directories to find calibrated LVM data files for a
range of exposures.

**Usage**::

    LocateData.py [-h] [-cp] [-dir whatever] [-drp VERSION] [-CFrame] exp_start [exp_stop]

**Options:**

-h
    Print help and exit.

-cp
    Copy files to a local data directory.

-dir path
    Alternative destination for copied files (implies -cp).

-drp VERSION
    Select data from a specific DRP version (e.g., 1.2.0).

-CFrame
    Locate CFrame files instead of SFrame files.

**Arguments:**

exp_start
    First exposure number to locate.

exp_stop
    Last exposure number (optional).

**Output:**

File locations are saved to a file in the ``xlog`` directory. If multiple
versions of a file exist, only the most recent is reported.

**Example**::

    # Locate and copy SFrame files for exposures 4155-4160
    LocateData.py -cp -drp 1.2.0 4155 4160


Typical Workflows
-----------------

Retrieving Existing Reduced Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the data has already been processed by the DRP at Utah::

    # Option 1: Direct retrieval
    GetFromUtah.py -cp -drp 1.2.0 60188 4155 4160

    # Option 2: Locate using drpall, then copy
    LocateReduced.py -cp 4155 4160

Reprocessing Data Locally
^^^^^^^^^^^^^^^^^^^^^^^^^

To reprocess data with the local DRP installation::

    # Process and copy to local directory
    Reduce.py -cp 60188 4155-4160

    # Check the log files for any issues
    ls xlog/

Finding Data Across Multiple MJDs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # Locate all exposures in a range
    LocateData.py -drp 1.2.0 10000 20000

    # The output file in xlog/ lists all found files


Notes
-----

- The scripts assume a standard LVM data directory structure
- Environment variables may need to be set to point to data locations
- The DRP must be installed locally to use ``Reduce.py``
- Processing large numbers of exposures can be parallelized with ``-np``


See Also
--------

- :doc:`summarize` - Tools for cataloging and summarizing data
- :doc:`api/GetFromUtah/index` - API documentation
- :doc:`api/Reduce/index` - API documentation
- :doc:`api/LocateData/index` - API documentation
- :doc:`api/LocateReduced/index` - API documentation
