Locating the Data For Multiple Objects
======================================

A common task is to determine which LVM exposures cover a set of
objects — for example, a catalog of supernova remnants, HII regions,
or galaxies — and how much total integration time has been accumulated
on each one.  The ``find_obs.py`` script does this by cross-matching
the pointing positions recorded in an LVM drpall file against a
user-supplied source catalog.

The drpall file is a FITS table maintained by the LVM Data Reduction
Pipeline that contains one row per reduced exposure, including the
telescope pointing position (``sci_ra``, ``sci_dec``), the MJD, the
exposure time, the tile ID, and the file location on disk.


Overview
--------

``find_obs.py`` performs a spherical cross-match between all pointings
in the drpall file and all sources in the catalog.  Every exposure
whose pointing centre falls within a specified angular radius of a
source is considered a match.  The match radius defaults to 900
arcseconds (15 arcminutes), which is approximately the radius of the
LVM science IFU.

The script produces two output files:

``<root>.matched.txt``
    One row for every observation-source pair that satisfied the match
    criterion.  Contains all columns from the source catalog and the
    key drpall columns (exposure number, MJD, exposure time, RA, Dec,
    tile ID, location), plus the angular separation in arcseconds
    between the source position and the pointing centre.

``<root>.sum.txt``
    One row per unique source in the catalog, summarising how many
    exposures matched (``Nobs``), the total accumulated exposure time
    in seconds (``TotExp``), and the median angular separation in
    arcseconds (``MedSep``).


Source Catalog Format
---------------------

The source catalog must be an ASCII file readable by
``astropy.io.ascii.read`` (whitespace-separated, CSV, ECSV, etc.) and
must contain at least the following three columns:

``RA``
    Right ascension in decimal degrees.

``Dec``
    Declination in decimal degrees.

``Source_name``
    A string identifier for each source.  This is used to group
    matched rows in the output and to label rows in the summary table.

Additional columns in the catalog are preserved and included in the
matched output file.  A minimal catalog looks like::

    Source_name   RA        Dec
    SNR_001       83.8221   -5.3911
    SNR_002       84.1055   -5.7812
    NGC_2070     83.8221   -69.0756


Command line usage::

    find_obs.py [-root outroot] [-max sep_arcsec] drpall_file source_catalog

**Positional arguments:**

drpall_file
    LVM drpall FITS file.

source_catalog
    ASCII source catalog with ``RA``, ``Dec``, and ``Source_name``
    columns.

**Options:**

-h
    Print the full documentation and exit.

-root outroot
    Root name for output files.  Default is ``results``, which gives
    ``results.matched.txt`` and ``results.sum.txt``.  The script
    refuses to overwrite an existing file; use a different root name
    to rerun with changed parameters.

-max sep_arcsec
    Cross-match radius in arcseconds.  Default is 900.  Use a smaller
    value if the catalog objects are densely packed and you want to
    restrict matches to pointings that are centred close to each source.


Typical Workflows
-----------------

Basic cross-match
^^^^^^^^^^^^^^^^^

Cross-match a catalog of LMC objects against the drpall file, using
the default 900 arcsecond radius, and write results with the root name
``lmc``::

    find_obs.py -root lmc drpall-1.2.0.fits lmc_catalog.txt

This writes ``lmc.matched.txt`` (all matched pairs) and
``lmc.sum.txt`` (one row per source with exposure counts and total
integration time).

Tighter match radius
^^^^^^^^^^^^^^^^^^^^

For a catalog of compact objects where you only want pointings that
are centred within 5 arcminutes::

    find_obs.py -root snr_300 -max 300 drpall-1.2.0.fits snr_catalog.txt

Checking coverage before downloading data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A useful pattern is to run ``find_obs.py`` first to identify which
objects have been observed and how much data exists, and then use
``LocateReduced.py`` or ``GetFromUtah.py`` to fetch the relevant
files::

    # Step 1: find which exposures cover your sources
    find_obs.py -root my_targets drpall-1.2.0.fits my_catalog.txt

    # Step 2: inspect the summary
    more my_targets.sum.txt

    # Step 3: retrieve the exposures for a source with good coverage
    #         (exposure numbers are listed in my_targets.matched.txt)
    LocateReduced.py -cp 12345 12350


Understanding the Output
------------------------

The matched file contains all the information needed to identify and
retrieve individual exposures.  The ``expnum`` column gives the LVM
exposure number, which can be passed directly to ``LocateReduced.py``
or ``GetFromUtah.py``.  The ``location`` column gives the path to the
reduced file in the standard LVM data tree.  The ``separation`` column
records how far the pointing centre was from the source position, which
can be used to filter matches by proximity.

The summary file gives a quick overview of coverage.  Sources with
``Nobs = 0`` do not appear (they had no matches at all).  A source
with a large ``TotExp`` and small ``MedSep`` is well-centred in many
exposures; a large ``MedSep`` indicates the source is near the edge
of the IFU or being covered by an adjacent pointing rather than a
dedicated one.


Notes
-----

Because ``search_around_sky`` returns all pairs within the match
radius, a single exposure can appear more than once in the matched
file if several catalog sources lie within the match radius of the
same pointing.  Similarly, a single source can be matched by many
overlapping exposures, including dithered or repeated visits.

The drpall file must be obtained separately.  Current drpall files are
available from the SDSS-V data archive at Utah.


See Also
--------

- :doc:`api/find_obs/index` - API documentation for ``find_obs.py``
- :doc:`data_retrieval` - Retrieving LVM data from Utah
- :doc:`summarize` - Summarizing LVM data across many exposures
