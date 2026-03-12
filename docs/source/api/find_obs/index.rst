find_obs
========

.. py:module:: find_obs

.. autoapi-nested-parse::

                       Space Telescope Science Institute

   Synopsis:

   Cross-match LVM observations from a drpall file against a source
   catalog to identify which exposures overlap each source position.


   Command line usage::

       find_obs.py [-root outroot] [-max sep_arcsec] drpall_file source_catalog

   The ``drpall_file`` argument is the LVM drpall FITS file containing one
   row per exposure with pointing coordinates, exposure times, and tile IDs.

   The ``source_catalog`` argument is an ASCII table with at least three
   columns: ``RA`` (degrees), ``Dec`` (degrees), and ``Source_name``.

   The ``-root outroot`` option sets the root name for output files.
   The default is ``results``, giving ``results.matched.txt`` and
   ``results.sum.txt``.

   The ``-max sep_arcsec`` option sets the cross-match radius in arcseconds.
   The default is 900 arcseconds (15 arcminutes), which corresponds roughly
   to the radius of one LVM pointing.


   Description:

   The routine reads the drpall FITS file and extracts the columns needed
   for the cross-match: exposure number, MJD, exposure time, pointing RA
   and Dec (``sci_ra``, ``sci_dec``), tile ID, and location.  It then reads
   the source catalog, which must contain ``RA``, ``Dec``, and
   ``Source_name`` columns.

   A spherical cross-match is performed using astropy SkyCoord, returning
   all observation-source pairs within the specified match radius.  Each
   matched pair records the observation metadata alongside the source
   catalog columns and the angular separation in arcseconds.

   A per-source summary table is also produced with the number of matching
   exposures, the total integration time, and the median separation between
   the source and the overlapping pointings.

   Two output files are written in fixed-width ASCII format:

   - ``<root>.matched.txt``: All matched pairs, one row per
     observation-source match, including all drpall columns and all
     source-catalog columns together with the angular separation.

   - ``<root>.sum.txt``: One row per unique source, reporting the source
     name, the number of matching exposures (``Nobs``), the total
     accumulated exposure time in seconds (``TotExp``), and the median
     angular separation in arcseconds (``MedSep``).


   Primary routines:

   The main driver is ``steer``, which parses the command line and calls
   ``load_observations``, ``find_matches``, and ``summarize`` in sequence.


   Notes:

   The source catalog must be readable by ``astropy.io.ascii.read`` and
   must contain columns named ``RA``, ``Dec``, and ``Source_name``.

   The drpall file must be a FITS file with the observation table in
   extension 1 and at least the columns ``expnum``, ``mjd``, ``exptime``,
   ``sci_ra``, ``sci_dec``, ``tileid``, and ``location``.

   The routine will refuse to overwrite existing output files.  Use a
   different ``-root`` name if you want to rerun with different parameters.

   The cross-match uses ``search_around_sky``, which considers all pairs
   within the match radius, not just the closest match.  A single exposure
   can therefore appear multiple times if several sources lie within the
   match radius of the same pointing.


   History:

   260312 ksl  Coding begun



Functions
---------

.. autoapisummary::

   find_obs.find_matches
   find_obs.load_observations
   find_obs.steer
   find_obs.summarize


Module Contents
---------------

.. py:function:: find_matches(source_file, obs, max_sep_arcsec)

   Cross-match a source catalog against LVM observations.

   Reads the source catalog from ``source_file`` using
   ``astropy.io.ascii.read`` and performs an all-pairs spherical
   cross-match against the observation table ``obs`` using
   ``SkyCoord.search_around_sky``.  All pairs within ``max_sep_arcsec``
   arcseconds are returned.  If observation and source catalog share
   column names, the observation copy is suffixed with ``_obs`` to
   avoid clashes before the two tables are horizontally stacked.

   Parameters:
       source_file (str): Path to the ASCII source catalog.  Must
           contain columns ``RA`` (degrees), ``Dec`` (degrees), and
           ``Source_name``.
       obs (astropy.table.Table): Observation table from
           ``load_observations``.
       max_sep_arcsec (float): Cross-match radius in arcseconds.

   Returns:
       astropy.table.Table: Matched table with one row per
       observation-source pair.  Includes all columns from both
       the source catalog and the observation table, plus a
       ``separation`` column in arcseconds.  Returns an empty
       Table if no matches are found.


.. py:function:: load_observations(drpall_file)

   Read an LVM drpall FITS file and return a table of observations.

   Opens the FITS file, reads extension 1, and extracts the columns
   needed for position cross-matching: exposure number, MJD, exposure
   time, pointing RA and Dec, tile ID, and location.  The RA and Dec
   column names are normalised to ``RA`` and ``Dec`` respectively.
   A ``Source_name`` column is added containing the zero-padded
   exposure number, which is used later to group matched rows.

   Parameters:
       drpall_file (str): Path to the LVM drpall FITS file.

   Returns:
       astropy.table.Table: Table with columns expnum, mjd, exptime,
       RA, Dec, tileid, location, and Source_name.


.. py:function:: steer(argv)

   Parse command line arguments and run the observation cross-match.

   Accepts the drpall FITS file and source catalog as positional
   arguments, with optional switches for the output root name and
   the cross-match radius.  Calls ``load_observations``,
   ``find_matches``, and ``summarize`` in sequence and writes the
   results to two fixed-width ASCII files.

   Parameters:
       argv (list): Command line argument list (typically ``sys.argv``).

   Returns:
       None


.. py:function:: summarize(matches)

   Compute a per-source summary from a table of matched pairs.

   Groups the matched pairs by ``Source_name`` and for each unique
   source computes the number of overlapping exposures, the total
   accumulated exposure time, and the median angular separation between
   the source and the matched pointing centres.

   Parameters:
       matches (astropy.table.Table): Matched pair table returned by
           ``find_matches``.  Must contain columns ``Source_name``,
           ``exptime``, and ``separation``.

   Returns:
       astropy.table.Table: Summary table with columns ``Source_name``,
       ``Nobs`` (int), ``TotExp`` (float, seconds), and
       ``MedSep`` (float, arcseconds), one row per unique source.


