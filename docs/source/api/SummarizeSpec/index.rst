SummarizeSpec
=============

.. py:module:: SummarizeSpec

.. autoapi-nested-parse::

                       Space Telescope Science Institute

   Synopsis:

   Create a summary of CFrame (flux-calibrated, before sky subtraction) data
   where each row contains the percentile spectrum for each of the three LVM
   spectrographs from a single exposure.  Useful for evaluating
   spectrograph-to-spectrograph variations before and after sky subtraction.

   Command line usage::

       SummarizeSpec.py [-h] [-sf] [-ver drp_ver] [-percent 50] [-emin 900]
                        [-out name] exp_start exp_stop delta

   Description:

   This script processes multiple CFrame files (the default) and computes the
   percentile spectrum for science fibers belonging to each of the three LVM
   spectrographs (spectrographid = 1, 2, 3) as recorded in the SLITMAP
   extension.  Only fibers with telescope == 'Sci' and fibstatus == 0 are
   included.  The -sf switch causes SFrame files to be read instead.

   File paths are taken from the location column of the drpall file, which
   records SFrame paths.  For CFrame files the path is adjusted automatically
   by replacing 'SFrame' with 'CFrame' in the location string.

   An exposure is included in the output only if science fibers are found
   for all three spectrographs.  Exposures where one or more spectrographs
   are missing are skipped and recorded in a separate ASCII table.

   Two output files are produced.  The main FITS file contains one row per
   valid exposure in each flux extension.  The skipped-exposure table lists
   every exposure that was rejected, with columns indicating which
   spectrographs were present and which were absent.

   Options: -h prints this documentation; -sf reads SFrame instead of CFrame
   files; -ver drp_ver sets the DRP version (default 1.2.0); -percent N sets
   the percentile to compute (default 50 = median); -emin N sets the minimum
   exposure time in seconds to include (default 900); -out name sets the root
   name of the output files.

   Positional arguments: exp_start is the first exposure number to consider;
   exp_stop is the last exposure number; delta processes every Nth exposure.

   Primary routines:

       doit

   Notes:

   The three flux extensions are named FLUX1, FLUX2, and FLUX3, corresponding
   to spectrographid values 1, 2, and 3 respectively.  Each is a 2-D array
   with one row per valid exposure and one column per wavelength pixel.

   The skipped-exposure table is always written, even if no exposures were
   skipped, so there is always a record that the completeness check was
   performed.

   History:

   260312 ksl Coding begun, based on SummarizeRings.py



Attributes
----------

.. autoapisummary::

   SummarizeSpec.XMUSKIE
   SummarizeSpec.XRAINBOW
   SummarizeSpec.XTOP


Functions
---------

.. autoapisummary::

   SummarizeSpec.augment_drp_all
   SummarizeSpec.doit
   SummarizeSpec.find_top
   SummarizeSpec.get_all_spec_specs
   SummarizeSpec.get_spec
   SummarizeSpec.make_spec_specs
   SummarizeSpec.read_drpall
   SummarizeSpec.select
   SummarizeSpec.steer


Module Contents
---------------

.. py:data:: XMUSKIE
   :value: '/home/long/Projects/lvm_data/sas'


.. py:data:: XRAINBOW
   :value: '/Users/long/Projects/lvm_data/sas'


.. py:data:: XTOP
   :value: '/uufs/chpc.utah.edu/common/home/sdss51/'


.. py:function:: augment_drp_all(xtab)

   Add survey classification, near/far sky telescope, and redshift columns
   to a drpall table.

   Parameters:
       xtab (astropy.table.Table): drpall table with sci_ra, sci_dec, skye_ra, skye_dec, skyw_ra, skyw_dec columns.

   Returns:
       astropy.table.Table: Input table with Survey, Near, Far, and
       Redshift columns added.


.. py:function:: doit(exp_start=4000, exp_stop=8000, delta=5, exp_min=900.0, out_name='', drp_ver='1.2.0', percentile=50, file_type='CFrame')

   Top-level driver: read drpall, select exposures, compute spectrograph spectra.

   Parameters:
       exp_start (int): First exposure number to process.
       exp_stop (int): Last exposure number to process.
       delta (int): Process every Nth exposure.
       exp_min (float): Minimum exposure time in seconds.
       out_name (str): Output file root name.
       drp_ver (str): DRP version string.
       percentile (int): Percentile to compute across fibers.
       file_type (str): 'CFrame' (default) or 'SFrame'.

   Returns:
       None


.. py:function:: find_top()

   Determine the top-level data directory for the current machine.

   Returns:
       str: Path to the top-level data directory, or an empty string
       if the location cannot be determined.


.. py:function:: get_all_spec_specs(filename, percentile=50)

   Compute percentile flux spectra for each spectrograph from one SFrame file.

   Opens the SFrame FITS file and checks whether science fibers are
   present for each of the three spectrographs (spectrographid 1, 2, 3).
   If any spectrograph has no valid science fibers the function returns
   None for all arrays together with a status dictionary that records
   which spectrographs were found.

   Parameters:
       filename (str): Path to the SFrame FITS file.
       percentile (int): Percentile to compute across fibers (default 50).

   Returns:
       tuple: (wav, flux1, flux2, flux3, status) where wav is the
       wavelength array, flux1/flux2/flux3 are 1-D percentile spectra
       for spectrographs 1/2/3, and status is a dict with keys 'SP1',
       'SP2', 'SP3' set to True if fibers were found and False otherwise.
       If the file cannot be opened, returns (None, None, None, None, None).
       If any spectrograph is missing, returns (None, None, None, None, status).


.. py:function:: get_spec(xtab, spectrograph_id)

   Select good science fibers for a given spectrograph from a SLITMAP table.

   Filters to rows where telescope == 'Sci', fibstatus == 0, and
   spectrographid == spectrograph_id.

   Parameters:
       xtab (astropy.table.Table): SLITMAP table from an SFrame file.
       spectrograph_id (int): Spectrograph identifier (1, 2, or 3).

   Returns:
       astropy.table.Table: Filtered table of matching fibers.


.. py:function:: make_spec_specs(xtab, data_dir, outfile='', percentile=50, file_type='CFrame')

   Process multiple CFrame or SFrame files and write a FITS summary by spectrograph.

   Iterates over the exposures in xtab, calling get_all_spec_specs for
   each file that exists on disk.  Exposures for which not all three
   spectrographs have valid science fibers are skipped and recorded in a
   separate ASCII table written alongside the main output file.

   The drpall location column records SFrame paths.  When file_type is
   'CFrame' (the default), 'SFrame' is replaced with 'CFrame' in the path
   before checking for the file.

   The main FITS file contains extensions WAVE, FLUX1, FLUX2, FLUX3, and
   drp_all.  FLUX1/FLUX2/FLUX3 are 2-D arrays (n_exposures x n_wavelengths)
   containing the percentile spectrum for spectrographs 1, 2, and 3.

   The skipped-exposure table is always written and contains one row per
   rejected exposure with columns expnum, mjd, tileid, SP1, SP2, SP3
   showing which spectrographs were present (True) or absent (False).

   Parameters:
       xtab (astropy.table.Table): drpall table of selected exposures.
       data_dir (str): Top-level data directory prepended to the location column of xtab to form full file paths.
       outfile (str): Output filename root (default 'test_spec.fits').
       percentile (int): Percentile to compute across fibers (default 50).
       file_type (str): 'CFrame' (default) or 'SFrame'.

   Returns:
       None


.. py:function:: read_drpall(drp_ver='1.2.0')

   Read the drpall FITS file and return an augmented metadata table.

   Looks for the file locally first, then falls back to the Utah path.

   Parameters:
       drp_ver (str): DRP version string (default '1.2.0').

   Returns:
       astropy.table.Table: drpall table with survey classification
       columns added, or an empty list if the file cannot be found
       or read.


.. py:function:: select(ztab, exp_start=4000, exp_stop=8000, delta=5, exp_min=900.0)

   Select every Nth exposure between exp_start and exp_stop.

   Parameters:
       ztab (astropy.table.Table): drpall table.
       exp_start (int): First exposure number to include.
       exp_stop (int): Last exposure number to include.
       delta (int): Step size; selects every Nth row.
       exp_min (float): Minimum exposure time in seconds.

   Returns:
       astropy.table.Table: Filtered table.


.. py:function:: steer(argv)

   Parse command line arguments and run the spectrograph summary.

   Parameters:
       argv (list): Command line argument list (typically sys.argv).

   Returns:
       None


