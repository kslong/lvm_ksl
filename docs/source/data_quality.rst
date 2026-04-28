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
All offsets are measured relative to the global median template.  A
scale-invariant quality value in [0, 1] is also returned for each fiber;
values below ~0.5 indicate unreliable measurements (featureless window,
bad pixels, or low S/N).

**Output:**

*Diagnostic figures* — one PNG per input file, saved to ``Fig_Qual/``.
The supertitle shows the exposure number and wavelength window.
Each figure has four panels:

1. Flux percentiles (5th, 50th, 95th) vs wavelength, y-scale set from
   the region ``[wmin, wmax]`` only.
2. Histogram of per-fiber wavelength offsets (range ±0.2 Å), annotated
   with the 5th and 95th percentile values.
3. Spatial map (RA/Dec) of the per-fiber maximum flux in ``[wmin, wmax]``,
   coloured by :math:`F_{\rm max}` with vmin/vmax at the 5th/95th
   percentile of all fiber values.  RA increases to the left.
4. Spatial map (RA/Dec) of the per-fiber wavelength offset, coloured by
   :math:`\delta\lambda`.  RA increases to the left.

*Summary table* — ``Fourier_<wmin>_<wmax>.<YYMMDD>.txt``, one row per
exposure.  If the file already exists, rows for exposures already present
are replaced with the newly computed values; rows for other exposures are
preserved.  The table is sorted by exposure number before writing.

Columns:

- ``Exposure`` — exposure number from the primary FITS header
- ``mjd`` — observation time converted to MJD (float)
- ``dw_5pct``, ``dw_95pct`` — 5th and 95th percentile offsets (Å)
- ``dw_med`` — median offset across all fibers (Å)
- ``dw_med_sp1``, ``dw_med_sp2``, ``dw_med_sp3`` — median offset per
  spectrograph
- ``med_quality`` — median cross-correlation quality across all fibers
- ``dw_r01_04_sp1/2/3`` — median offset for rings 1–4 combined, per
  spectrograph (rings 1–4 are combined because they contain too few
  fibers to measure reliably on their own)
- ``dw_r05_sp1/2/3`` through ``dw_r25_sp1/2/3`` — median offset per
  ring (5–25) per spectrograph

The table has 75 columns in total: 9 summary columns plus 66
ring/spectrograph columns (22 ring groups × 3 spectrographs).

After each exposure is processed, the 95th-percentile offset and the
median quality are printed to the terminal::

    95th pct offset: 0.0312 AA   median quality: 0.847

**Example**::

    # Analyse a single SFrame file using the default Ca H&K window
    fourier_offset.py lvmSFrame-00012345.fits

    # Analyse multiple files with a custom wavelength window
    fourier_offset.py -wmin 6540 -wmax 6590 lvmSFrame-*.fits


fourier_offset_check.py — Synthetic-Line Injection Test
--------------------------------------------------------

Adds a synthetic Gaussian emission line to the FLUX extension of an
lvmCFrame file and writes the result to a new file.  Intended for
verifying that ``fourier_offset.py`` recovers known per-spectrograph
wavelength offsets.

**Command line usage**::

    fourier_offset_check.py [-h] [-wave 4200] [-flux 1e-11]
                            [-off1 0.0] [-off2 0.1] [-off3 -0.1] filename

**Options:**

-h
    Print help and exit.

-wave WAVE
    Central wavelength of the synthetic line in Angstroms (default: 4200).

-flux FLUX
    Integrated line flux in erg/s/cm² (default: 1e-11).

-off1 OFF, -off2 OFF, -off3 OFF
    Wavelength offset in Angstroms applied to all fibers in spectrographs
    1, 2, and 3 respectively (defaults: 0.0, +0.1, −0.1).

**Arguments:**

filename
    Path to the lvmCFrame FITS file to modify.

**Method:**

The FWHM at the requested wavelength is read from the LSF extension for
each fiber and converted to a Gaussian sigma.  The peak amplitude is
derived from the integrated flux as :math:`A = F / (\sigma \sqrt{2\pi})`.
The line centre is shifted by the per-spectrograph offset before the
Gaussian is evaluated and added to the FLUX array.  Fibers whose LSF
FWHM is zero are skipped.

**Output:**

A new FITS file named ``test_<basename>.fits`` in the current directory.
A one-line summary is printed to the terminal showing the output filename,
line wavelength, flux, and per-spectrograph offsets.

**Example**::

    # Inject a line at 4200 Å with default offsets (0, +0.1, -0.1 AA)
    fourier_offset_check.py lvmCFrame-00004171.fits

    # Custom line and offsets
    fourier_offset_check.py -wave 5007 -flux 1e-12 -off2 0.2 lvmCFrame-*.fits


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
- :doc:`api/fourier_offset_check/index` - API documentation
- :doc:`api/CheckData/index` - API documentation
- :doc:`api/CheckReduced/index` - API documentation
- :doc:`api/eval_sky/index` - API documentation
- :doc:`api/eval_standard/index` - API documentation
- :doc:`summarize` - Tools for cataloging and summarizing exposures
