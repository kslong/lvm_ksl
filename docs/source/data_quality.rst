Data Quality
============

The lvm_ksl package includes tools for assessing the quality of reduced
LVM data.  They fall into a few distinct categories, grouped below by
what they actually check rather than listed alphabetically:

- **Wavelength calibration** — how much the wavelength solution drifts
  from fiber to fiber within an exposure (``fourier_offset.py``,
  ``fourier_offset_check.py``).
- **Reduction bookkeeping** — whether a batch of downloaded or locally
  reduced files is complete, current, and free of DRP errors
  (``CheckData.py``, ``CheckReduced.py``).
- **Sky subtraction** — how well the subtracted sky matches what was
  actually there, both spectroscopically and fiber-by-fiber
  (``eval_sky.py``, ``plot_sky_gaussfit.py``).
- **Flux calibration** — how well calibrated standard-star spectra agree
  with their Gaia reference spectra (``eval_standard.py``).
- **Sky telescope pointing** — whether a sky exposure's recorded
  position actually agrees with the sky field name it was labelled with
  (``SummarizeSkyHdr.py``, ``check_sky_positions.py``).
- **Combined quality report** — a single-exposure HTML report combining
  header overview, sky subtraction, and flux calibration checks
  (``Quicklook.py``).


Wavelength Calibration
-----------------------

Checks whether the wavelength solution is stable across fibers within a
single exposure, using a synthetic-line injection test to validate the
measurement itself.

fourier_offset.py — Wavelength Offset Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


Reduction Bookkeeping
----------------------

Checks whether a batch of downloaded or locally reduced files is
complete, up to date with the current DRP version, and free of DRP
processing errors — not a check on the science content of the data
itself.

CheckData.py — Inspect Downloaded Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


Sky Subtraction
-----------------

Checks how well the subtracted sky spectrum matches what was actually
present, both as a per-exposure spectroscopic comparison and as
fiber-by-fiber spatial maps of individual airglow-line fit residuals.

eval_sky.py — Sky Subtraction Quality Plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


plot_sky_gaussfit.py — Sky Subtraction Residual Maps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Visualises fiber-by-fiber sky Gaussian fit results produced by
``sky_gaussfit.py`` as spatial scatter maps (fiber RA vs Dec coloured by
the residual), one page per fitted quantity, via
:doc:`plotting_outputs`'s ``radec_plot.plot_scatter()``.  Three PNG files
are written per input file:

- ``<root>_wave.png`` — centroid wavelength residuals (Å)
- ``<root>_flux.png`` — flux residuals (fractional)
- ``<root>_fwhm.png`` — FWHM residuals (km/s)

Each page is a 6-row × 3-column grid, one panel per airglow line (18 lines
total).  Each panel subtracts the per-panel median before plotting so the
colorbar shows the residual; the suptitle records the colour range applied.
Panel titles show the line name, median value, and MAD-based standard
deviation.

**Command line usage**::

    plot_sky_gaussfit.py [-out root] [-s size] filename [filename ...]

**Options:**

-out root
    Root name for the output PNGs when a single input file is given.
    Ignored when multiple files are provided (each uses its own stem).

-s size
    Scatter marker size in points² (default: auto -- sized from the
    fiber spacing and the axes' rendered size; see
    :doc:`plotting_outputs`).

**Arguments:**

filename
    One or more ASCII fixed-width tables written by ``sky_gaussfit.py``.
    Each file is processed independently; no stacking is performed.

**Colour limits:**

The limits are controlled by ``VRANGE_FRAC`` (for wave and flux) and
``FWHM_VRANGE_KMS`` (for FWHM) near the top of the script.

.. list-table::
   :header-rows: 1
   :widths: 10 50 20

   * - Quantity
     - Scale
     - Default range
   * - wave
     - ± (median_wave × 5 km/s / c) Å
     - ±5 km/s equiv.
   * - flux
     - Fractional variation of median flux
     - ±2 %
   * - fwhm
     - Absolute, in km/s (FWHM_Ang / lambda × c)
     - ±30 km/s

FWHM is converted to velocity space using each line's central wavelength
so that all 18 panels share a common km/s scale, consistent with the
quadrature broadening formula FWHM_obs² = FWHM_inst² + σ².

**Output:**

Three PNG files saved to ``Figs_gaussfit_sky/``:

- ``Figs_gaussfit_sky/<root>_wave.png``
- ``Figs_gaussfit_sky/<root>_flux.png``
- ``Figs_gaussfit_sky/<root>_fwhm.png``

**Typical workflow**::

    # 1. Fit all science fibers in an SFrame file
    sky_gaussfit.py -lmc lvmSFrame-00012345.fits

    # 2. Plot spatial residual maps from the fit results
    plot_sky_gaussfit.py lvmSFrame-00012345.txt

    # 3. Process many files at once (one PNG set per file)
    plot_sky_gaussfit.py lvmSFrame-*.txt


Flux Calibration
-----------------

Checks how well flux-calibrated standard-star spectra agree with their
Gaia BP/RP reference spectra.

eval_standard.py — Flux Standard Calibration Plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


Sky Telescope Pointing
------------------------

Checks whether a sky exposure's *recorded* telescope position actually
agrees with the sky field name it was *labelled* with, at each stage of
the DRP's astrometry processing.  This grew out of a real, confirmed
data-quality bug: some SkyE/SkyW exposures carry a field name that
doesn't match where the telescope was actually pointed — most commonly a
clean east/west name swap, but also partial and unexplained variants
(see ``check_sky_positions.py``'s module docstring history for the full
diagnostic trail).

SummarizeSkyHdr.py — Per-Exposure Astrometry Keyword Extraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pulls a fixed, hardwired list of raw acquisition/astrometry PRIMARY
header keywords out of each exposure's CFrame file, for a range of
exposures selected from a drpall table.  These keywords carry three
independent position pairs per sky telescope — what it *reported* back,
what was *commanded*, and the final *adopted* (astrometry-refined)
position — plus the astrometry-source quality flags (SCIASRC/SKYEASRC/
SKYWASRC).  drpall itself only ever stores the adopted position, so this
is the tool that makes the earlier pipeline stages visible for
diagnosis; ``check_sky_positions.py`` is what actually checks them.

**Command line usage**::

    SummarizeSkyHdr.py [-h] [-emin 900] [-ver 1.2.1] [-drp_all FILE]
                       [-keywords FILE] [-data_dir DIR] [-out ROOT]
                       exp_start exp_stop [delta]

**Arguments:**

exp_start, exp_stop
    Exposure number range to select from the drpall table.

delta
    Process every delta-th exposure in the range (default 1).

**Options:**

-h
    Print help and exit.

-emin N
    Minimum exposure time to include (default 900).

-ver VER
    DRP version, used to locate ``drpall-VER.fits`` (default 1.2.1).

-drp_all FILE
    Explicit drpall table to read instead of ``drpall-VER.fits``.

-keywords FILE
    Optional (keyword, definition) table, overriding the hardwired
    keyword list built into the script.  Not needed for normal use — the
    script has no external file dependency by default, so it can be
    copied to and run standalone at Utah.

-data_dir DIR
    Look for CFrame files directly in this flat local cache by basename
    first, before falling back to the standard xtop/location tree.

-out ROOT
    Output filename root (default:
    ``SummarizeSkyHdr_<ver>_<exp_start>_<exp_stop>_<delta>``).

**Output:**

A FITS file with three extensions:

- ``PRIMARY`` — records the calling parameters (DRPVER, EMIN, EXPSTART,
  EXPSTOP, DELTA, DRPALL, KWFILE, N_PROC).
- ``SKY_HDR`` — one row per successfully-read exposure: EXPNUM (the join
  key back to DRP_ALL) plus one column per keyword.  Each column's FITS
  ``TTYPEn`` comment card carries that keyword's definition, so the file
  is self-documenting.
- ``DRP_ALL`` — the drpall rows for the exposures actually processed.

**Example**::

    # Summarize astrometry keywords for a range of exposures
    SummarizeSkyHdr.py 7325 48860 1


check_sky_positions.py — Sky Field Name/Position Cross-Check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cross-references the recorded skye/skyw positions in a drpall file (or a
``SummarizeSkyHdr.py`` output file) against the nominal catalog positions
in ``final_sky_tiles.csv``, to find exposures where a sky field's
recorded position doesn't match its own name.  Every problem exposure is
classified as ``Swapped`` (skye and skyw's names are cleanly crossed with
each other), ``HalfMatch`` (only one direction of that crossing holds —
suggestive of a shifted-by-one-slot pattern rather than a clean swap), or
``Unknown`` (no recognized pattern).

**Command line usage**::

    check_sky_positions.py [-csv FILE] [-tol DEG] [-out FILE]
                           [-postype reported|commanded|adopted] [drpall_file]

**Arguments:**

drpall_file
    A drpall FITS file (default: ``drpall-1.2.1.fits``), or a
    ``SummarizeSkyHdr.py`` output file if ``-postype`` is anything other
    than ``adopted``.

**Options:**

-h
    Print help and exit.

-csv FILE
    Sky tile position catalog (default: ``final_sky_tiles.csv``).

-tol DEG
    Agreement tolerance in degrees (default: 0.1).

-out FILE
    Override the per-exposure problem table's output filename (default:
    ``sky_problem_check_<drpall stem>[_<postype>].tab``).

-postype T
    Which of the three recorded sky-telescope positions to check against
    — ``reported``, ``commanded``, or ``adopted`` (default; the only one
    that works on a plain ``drpall-*.fits``).  ``reported``/``commanded``
    require a ``SummarizeSkyHdr.py`` output file's SKY_HDR extension.
    Comparing across postype values traces which DRP stage a mismatch
    first appears at.

**Output:**

Two ASCII fixed-width tables are always written together:

- ``sky_problem_check_<drpall stem>[_<postype>].tab`` — one row per
  problem exposure: tileid, mjd, expnum, filename, skye_name, skyw_name,
  each side's real nearest-catalog position, problem_type, and (when
  available) the SCIASRC/SKYEASRC/SKYWASRC astrometry-quality flags.
- ``sky_position_check_<drpall stem>[_<postype>].tab`` — one row per
  catalog field: how many times it was labelled, how often that label
  was correct, and a breakdown of the incorrect ones into Swapped/
  HalfMatch/Unexplained, plus how many times the field's real position
  was observed under any label (and how many of those were mislabeled).

Both are also printed to the terminal, each followed by a summary
(problem counts by type for the per-exposure table; clean/Swapped/
HalfMatch/Unexplained field counts for the per-name table).

**Example**::

    # Check the adopted (final) positions in a plain drpall file
    check_sky_positions.py drpall-1.2.1.fits

    # Trace the same check back to the commanded stage, loosening the
    # tolerance to 1 degree
    check_sky_positions.py -postype commanded -tol 1 SummarizeSkyHdr_1.2.1_7325_48860_1.fits


Combined Quality Report
-------------------------

Produces a single self-contained HTML report for one exposure, combining
a header overview with the sky-subtraction and flux-calibration checks
above, so that an exposure can be assessed at a glance without running
several scripts separately.

Quicklook.py — Per-Exposure HTML Quality Report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reads an lvmSFrame file and builds an HTML file containing header
information plus the science/sky spectral comparison, Hα/[SII]/continuum
images, and (if possible) the flux-calibrated standard-star comparison.

**Command line usage**::

    Quicklook.py [-h] SFrame1 SFrame2 ...

**Options:**

-h
    Print help and exit.

**Arguments:**

filename
    One or more lvmSFrame FITS files to analyse. Each produces its own
    HTML report.

**Output:**

- ``<root>.html`` — the report, written to the current working directory
  (``<root>`` is the SFrame filename with its directory and ``.fits``
  extension stripped). Image links in the file are relative, so the
  report and the ``figs_qual/`` directory below must be kept together.
- ``figs_qual/`` — subdirectory holding all PNGs referenced by the
  report (science/sky spectra from ``eval_qual_sframe``, line and
  continuum images from ``make_images``, and the standard-star
  comparison from ``eval_standard.qual_eval``).

**Overview section:**

The top of the report lists, from the SFrame's PRIMARY header: exposure
number, MJD, observation time, object name, DRP version (``DRPVER``) and
commit hash (``COMMIT``), the science and sky-telescope RA/Dec/PA (with
angular distance from the science pointing), and the Moon/Sun RA, Dec,
altitude, and (for the Moon) illumination at Las Campanas. Any of these
header keywords that are missing falls back to a placeholder (``Unknown``
for strings, ``-999.0`` for numbers) rather than raising an error, since
not every keyword is present in every DRP version's headers.

**Notes:**

The standard-star comparison panel requires ``lvmdrp`` (via
``eval_standard.py``); if it isn't available, the report notes that the
comparison could not be done and continues without it.

**Example**::

    Quicklook.py data/lvmSFrame-00012345.fits


See Also
--------

- :doc:`api/fourier_offset/index` - API documentation
- :doc:`api/fourier_offset_check/index` - API documentation
- :doc:`api/CheckData/index` - API documentation
- :doc:`api/CheckReduced/index` - API documentation
- :doc:`api/eval_sky/index` - API documentation
- :doc:`api/eval_standard/index` - API documentation
- :doc:`api/SummarizeSkyHdr/index` - API documentation
- :doc:`api/check_sky_positions/index` - API documentation
- :doc:`api/QuickLook/index` - API documentation
- :doc:`summarize` - Tools for cataloging and summarizing exposures
- :doc:`spectral_fitting` - ``sky_gaussfit.py`` produces the input tables for ``plot_sky_gaussfit.py``
- :doc:`plotting_outputs` - ``radec_plot.py``, which ``plot_sky_gaussfit.py`` now uses for its spatial rendering
