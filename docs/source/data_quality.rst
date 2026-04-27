Data Quality
============

The lvm_ksl package includes tools for assessing the quality of reduced
LVM data.  Currently the focus is on wavelength calibration stability,
characterised by measuring per-fiber wavelength offsets across an exposure
using Fourier cross-correlation.


fourier_offset.py — Wavelength Offset Analysis
-----------------------------------------------

Measures per-fiber wavelength offsets in a CFrame or SFrame file by
cross-correlating each fiber spectrum against a median reference template
within a user-specified wavelength window.  The offsets quantify how much
the wavelength solution drifts from fiber to fiber, which is a useful
diagnostic for wavelength calibration quality.

**Command line usage**::

    fourier_offset.py [-h] [-wmin 3900] [-wmax 4000] filename [filename ...]

**Options:**

-h
    Print help and exit.

-wmin WAVE
    Short-wavelength edge of the cross-correlation window in Angstroms
    (default: 3900).

-wmax WAVE
    Long-wavelength edge of the cross-correlation window in Angstroms
    (default: 4000).

**Arguments:**

filename
    One or more lvmCFrame or lvmSFrame FITS files to analyse.

**Method:**

Each fiber spectrum is cross-correlated with the median of all 1754
science fibers in the wavelength window using FFTs.  Sub-pixel precision
is obtained by parabolic interpolation of the cross-correlation peak.
All offsets are therefore measured relative to the global median template.

**Output:**

*Diagnostic figures* — one PNG per input file, saved to ``Fig_Qual/``.
Each figure has four panels:

1. Flux percentiles (5th, 50th, 95th) vs wavelength
2. Histogram of per-fiber wavelength offsets
3. Spectra at the 5th- and 95th-percentile offsets
4. Spatial map of offsets (RA/Dec), coloured by :math:`\delta\lambda`

*Summary table* — ``Fourier_offsets_<wmin>_<wmax>.txt``, one row per
exposure, with the following columns:

- ``Exposure`` — exposure number from the primary FITS header
- ``dw_5pct``, ``dw_95pct`` — 5th and 95th percentile offsets (Å)
- ``dw_med`` — median offset across all fibers (Å)
- ``dw_med_sp1``, ``dw_med_sp2``, ``dw_med_sp3`` — median offset per
  spectrograph
- ``dw_r01_04_sp1/2/3`` — median offset for rings 1–4 combined, per
  spectrograph (rings 1–4 are combined because they contain too few
  fibers to measure reliably on their own)
- ``dw_r05_sp1/2/3`` through ``dw_r25_sp1/2/3`` — median offset per
  ring (5–25) per spectrograph

The table therefore has 73 columns in total: 7 summary columns plus
66 ring/spectrograph columns (22 ring groups × 3 spectrographs).

**Example**::

    # Analyse a single SFrame file using the default Ca H&K window
    fourier_offset.py lvmSFrame-00012345.fits

    # Analyse multiple files with a custom wavelength window
    fourier_offset.py -wmin 6540 -wmax 6590 lvmSFrame-*.fits


CheckData.py — Inspect Downloaded Files
----------------------------------------

Prints and saves a summary table of key header keywords for one or more
downloaded CFrame or SFrame files.  Useful for quickly checking whether
the files on disk are up to date with the current DRP version.

**Command line usage**::

    CheckData.py filename [filename ...]

**Arguments:**

filename
    One or more lvmCFrame or lvmSFrame FITS files to inspect.

**Output:**

Prints a table to the screen and writes it to ``DataSum.txt``.  Columns
include: ``Filename``, ``MJD``, ``DRP``, ``Commit``, ``FluxCal``,
``BUNIT``, ``RA``, ``Dec``, ``HelioV``, ``Tile_ID``, ``EXPTIME``,
``Source_name``.

**Example**::

    CheckData.py data/lvmSFrame-*.fits


CheckReduced.py — Verify Local DRP Reductions
----------------------------------------------

Checks whether a set of locally reduced files all succeeded and whether
they were all produced from the same DRP git commit.  Intended to be run
after ``Reduce.py``, which stores reduced SFrame files in ``data/`` and
DRP log files in ``xlog/``.

**Command line usage**::

    CheckReduced.py [-h] [-d data_dir] [-l log_dir]

**Options:**

-h
    Print help and exit.

-d data_dir
    Directory containing the reduced SFrame files (default: ``data``).

-l log_dir
    Directory containing the DRP log files written by ``Reduce.py``
    (default: ``xlog``).

**Output:**

- Prints a per-commit count to the screen.
- Writes per-file commit information to ``commits.txt``.
- Writes any DRP ERROR lines found in the log files to ``problems.txt``.

**Example**::

    # Run with defaults after Reduce.py has finished
    CheckReduced.py

    # Specify alternative directories
    CheckReduced.py -d mydata -l mylogs


eval_sky.py — Sky Subtraction Quality Plot
-------------------------------------------

Produces a four-panel plot showing median flux, sky, and total
(flux + sky) spectra for the science fibers, the east sky telescope
(SkyE), and the west sky telescope (SkyW) of a single SFrame file.
Provides a quick visual check of how well sky subtraction has worked.

**Command line usage**::

    eval_sky.py filename [filename ...]

**Arguments:**

filename
    One or more lvmSFrame FITS files to evaluate.

**Output:**

One PNG file per input file, named ``sky_<basename>.png``, written to
the current directory.

**Notes:**

This routine operates on individual SFrame files only — it does not work
on median-combined spectra.

**Example**::

    eval_sky.py data/lvmSFrame-00012345.fits


eval_standard.py — Flux Standard Calibration Plot
--------------------------------------------------

Compares the observed standard star spectra in a CFrame file against
their Gaia BP/RP reference spectra, providing a visual check of the
flux calibration quality.

**Command line usage**::

    eval_standard.py filename [filename ...]

**Arguments:**

filename
    One or more lvmCFrame FITS files to evaluate.

**Output:**

One PNG file per input file, named ``standard_<basename>.png``, written
to the current directory.  If no Gaia spectra can be retrieved for any
standard in the file, no plot is produced and a warning is printed.

**Notes:**

Requires ``lvmdrp`` to be installed (uses ``ancillary_func.retrive_gaia_star``
to fetch Gaia reference spectra).

**Example**::

    eval_standard.py data/lvmCFrame-00012345.fits


See Also
--------

- :doc:`api/fourier_offset/index` - API documentation
- :doc:`api/CheckData/index` - API documentation
- :doc:`api/CheckReduced/index` - API documentation
- :doc:`api/eval_sky/index` - API documentation
- :doc:`api/eval_standard/index` - API documentation
- :doc:`summarize` - Tools for cataloging and summarizing exposures
