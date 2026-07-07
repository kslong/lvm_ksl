Sky Subtraction
===============

Sky subtraction is one of the most critical steps in LVM data reduction.
The lvm_ksl package provides tools for both evaluating the standard DRP
sky subtraction and for experimenting with alternative approaches using
the ESO SkyCorr tool.

This page describes the available sky subtraction and sky modeling tools.


Overview
--------

LVM uses dedicated sky telescopes (SKY_EAST and SKY_WEST) to measure the
sky spectrum simultaneously with science observations. The DRP uses these
sky spectra to subtract the sky from the science fibers.

The tools in lvm_ksl allow you to:

- Evaluate sky subtraction quality in DRP-processed data
- Run alternative sky subtraction using ESO's SkyCorr tool, a polynomial or
  B-spline continuum fit (SkySubOrig/SkySubDev1), a PALACE decomposition
  (SkySubDev2), the lvmdrp routine directly (SkySubDrp), or the ESO Sky
  Model itself (SkySepESO)
- Generate theoretical sky models for comparison
- Visualize sky residuals and identify problems, including separately
  evaluating sky *line* subtraction quality and continuum *separation*
  quality (SkySub_eval.py)


Evaluating Sky Subtraction
--------------------------

eval_sky.py
^^^^^^^^^^^

Creates diagnostic plots to evaluate sky subtraction quality in SFrame
files (sky-subtracted data from the DRP).

**Usage**::

    eval_sky.py filename1 filename2 ...

Plots median, min, max, and average spectra across science fibers to
identify residual sky features.

sky_plot.py
^^^^^^^^^^^

Similar to eval_sky.py, creates plots to evaluate sky subtraction
quality with additional visualization options.

**Usage**::

    sky_plot.py filename1 filename2 ...


Alternative Sky Subtraction with SkyCorr
----------------------------------------

ESO's SkyCorr is a tool for sky subtraction that can handle cases where
the sky varies spatially or temporally. The following scripts provide
an interface to SkyCorr for LVM data.

Prep4SkyCorr.py
^^^^^^^^^^^^^^^

Prepares LVM CFrame data for use with SkyCorr by extracting mean spectra
from science and sky fibers into a format SkyCorr can read.

**Usage**::

    Prep4SkyCorr.py [-h] [-dir whatever] [-all] filename1 filename2 ...

**Options:**

-h
    Print help and exit.

-dir path
    Search for files in the specified directory.

-all
    Process all lvmCFrame files in the directory.

**Output:**

Creates FITS files with mean spectra suitable for SkyCorr input.

RunSkyCorr.py
^^^^^^^^^^^^^

Runs the ESO SkyCorr tool on prepared input files.

**Usage**::

    RunSkyCorr.py sci.fits sky.fits

where ``sci.fits`` is the science spectrum and ``sky.fits`` is the sky
spectrum (both prepared by Prep4SkyCorr.py).

**Requirements:**

1. SkyCorr must be installed and the ``skycorr`` executable in your PATH
2. A parameter file ``lvm_base.par`` must exist in the working directory

**Setup:**

On first run, if ``lvm_base.par`` does not exist, the script will create
a template. You must edit this file to set the correct paths::

    INST_DIR=/path/to/skycorr/installation/
    INPUT_OBJECT_SPECTRUM=/path/to/data/XXOBJECT.fits
    INPUT_SKY_SPECTRUM=/path/to/data/XXSKY.fits
    OUTPUT_DIR=/path/to/output/

Note: SkyCorr requires absolute paths. The placeholders XXOBJECT and
XXSKY are replaced automatically by the script.

SkySub.py
^^^^^^^^^

Performs simple sky subtraction by scaling the sky spectrum to match
the science spectrum and subtracting.

**Usage**::

    SkySub.py filename

This provides a quick alternative to the DRP sky subtraction for testing
purposes. The scaling factor is determined by minimizing residuals in
regions dominated by sky emission.


Sky Modeling
------------

These tools generate theoretical sky spectra using ESO models, which can
be compared to observed sky spectra for validation.

SkyCalcObs.py
^^^^^^^^^^^^^

Uses ESO's SkyCalc web service to generate a theoretical sky spectrum
for a given position and time.

**Usage**::

    SkyCalcObs.py [-h] [-out name] ra dec time

**Arguments:**

ra, dec
    Sky position in degrees.

time
    Observation time as a date string, JD, or MJD.

**Options:**

-h
    Print help and exit.

-out name
    Set output filename root.

**Requirements:**

The ``skycalc_cli`` package must be installed::

    pip install skycalc_cli

**Output:**

A FITS file containing the theoretical sky spectrum.

SkyModelObs.py
^^^^^^^^^^^^^^

Uses the ESO Sky Model (local installation) to generate theoretical
sky spectra. This is faster than SkyCalcObs for batch processing.

**Usage**::

    SkyModelObs.py [-h] [-config] [-data data_dir] [-out name] ra dec time

**Arguments:**

ra, dec
    Sky position in degrees.

time
    Observation time as a date string, MJD, or JD.

**Options:**

-h
    Print help and exit.

-config
    Force reconfiguration of directories.

-data data_dir
    Set data directory (overrides ESO_SKY_MODEL environment variable).

-out name
    Set output filename root.

**Requirements:**

The ESO Sky Model package must be installed locally, with the
``calcskymodel`` executable available.


PALACE-based Sky Line Masking and Sky Spectrum Collection
----------------------------------------------------------

These tools support the development of improved sky subtraction by identifying
sky-line-free wavelength windows for continuum fitting and by assembling
stacked sky spectra from the LVM sky telescopes.  Together they are intended
to characterise the sky background well enough to constrain physical models
of the airglow emission.

XSkySepIvan.py
^^^^^^^^^^^^^^

Decomposes LVM sky spectra into physical emission components using the
PALACE (Paranal Airglow Line And Continuum Emission, Noll et al. 2024) line
model combined with a B-spline Moon/zodiacal continuum.  This script is a
wrapper for the ``SkyDecomp`` class written by Ivan Katkov (lvmsky repository,
``skysub/sky_decomp/fit.py``).

The decomposition solves a non-negative quadratic programme (Clarabel solver)
for the amplitudes of six component families:

- **OH** — 402 groups of hydroxyl vibrational-rotational lines (HITRAN data)
- **Moon** — B-spline envelope multiplied by a rebinned solar spectrum (captures Moon reflected light and zodiacal continuum)
- **Diffuse** — PALACE diffuse continuum: HO\ :sub:`2`, FeO, O\ :sub:`2`\Ac
- **Atom** — atomic airglow: NaI, KI, [NI], OI (green and red)
- **ORC** — OI recombination multiplets at 7774 and 8446 Å
- **O2** — molecular oxygen A-band (~8650 Å); rotational temperature fitted

The script supports two input modes:

*Sky-file mode* (files produced by ``GetSky_from_CFrame_sum.py``): each row
of the FLUX array is a sky-telescope spectrum of the same field from a
different exposure.  Noise is estimated from pixel-to-pixel differences since
no IVAR extension is present.

*XCframe mode*: standard LVM summary file with FLUX/SKY_EAST/SKY_WEST
extensions.  IVAR is read from the file if present, otherwise estimated.

**Usage**::

    # Process all spectra in a Sky file
    XSkySepIvan.py Sky_WHAM_south_08.fits

    # Process specific rows of a Sky file
    XSkySepIvan.py Sky_WHAM_south_08.fits 0 5 10

    # Process all rows of an XCframe extension
    XSkySepIvan.py XCframe_file.fits SKY_EAST

    # Process specific rows of an XCframe extension
    XSkySepIvan.py XCframe_file.fits SKY_EAST 0 300 600

    # Process every 1000th row of an XCframe extension
    XSkySepIvan.py XCframe_file.fits SKY_EAST -delta 1000

**Key options:**

-delta N
    Process every N-th row (0, N, 2N, ...) instead of all rows.
    Ignored if explicit row numbers are given.

-lsf FWHM
    Fixed LSF FWHM in Angstroms (default 1.3 Å).

-refits N
    Number of iterative per-channel LSF kernel refits (default 0).

-out outroot
    Set output filename root.

**Output:**

A FITS file named ``<stem>_<ext>_ivan.fits`` (XCframe mode) or
``<stem>_ivan.fits`` (Sky-file mode) with the following extensions:

*Spectral image extensions* (float32, N_obs × N_pix):

- ``WAVE`` — wavelength array (Å)
- ``FLUX`` — input sky spectra
- ``LINES`` — total emission (OH + ATOM + ORC + O2)
- ``CONT`` — total continuum (MOON + DIFFUSE)
- ``OH`` — OH band component
- ``ATOM`` — atomic airglow (NaI, KI, [NI], OI)
- ``ORC`` — OI recombination multiplets
- ``O2`` — O2 A-band component
- ``MOON`` — Moon/zodiacal continuum
- ``DIFFUSE`` — diffuse airglow continuum
- ``RESID`` — fit residuals (FLUX − LINES − CONT)

*Coefficient tables:*

- ``COEF`` — BinTable with 442 named columns (one per design-matrix
  entry), storing the physical-unit coefficient for each spectrum.
  Column names match the ``name`` column of COEF_META:
  ``OH_v{v}_N{nn}_F{f}``, ``Moon_bs{nn}``, ``HO2``, ``FeO``,
  ``O2Ac``, ``NaI``, ``KI``, ``NI_forb``, ``OI_5577``,
  ``OI_6300``, ``OI_7774``, ``OI_8446``, ``O2_Aband``.

- ``COEF_META`` — BinTable with 442 rows, one per design-matrix
  entry, carrying physical identity (``name``, ``component``,
  ``v_upper``, ``N_upper``, ``F_upper``, ``wave_peak``) and
  cross-spectrum distribution statistics (``coef_median``,
  ``coef_mean``, ``coef_nmad``, ``coef_min``, ``coef_max``,
  ``coef_skew``).

- ``DRP_ALL`` — observation metadata plus 20 compact coefficient
  summary columns (physical flux units, same names as COEF_META):
  ``OH_v3``…``OH_v10`` (total OH flux per vibrational band),
  ``NaI``, ``KI``, ``NI_forb``, ``OI_5577``, ``OI_6300``,
  ``OI_7774``, ``OI_8446``, ``O2_Aband``, ``HO2``, ``FeO``,
  ``O2Ac``, ``Moon_med``.

XSkySepIvan_eval.py
^^^^^^^^^^^^^^^^^^^

Evaluates the quality of a PALACE sky decomposition produced by
``XSkySepIvan.py`` by creating an interactive three-panel HTML plot over a
chosen wavelength window.

Each panel shows the median spectrum (coloured line) surrounded by a shaded
10th/90th percentile band, with a random sample of individual spectra overlaid
in grey so outliers and systematic trends are immediately visible:

- **Panel 1 (Flux)** — input sky spectra
- **Panel 2 (Residual)** — fit residuals (FLUX − model)
- **Panel 3 (Lines)** — total fitted emission (OH + Atom + ORC + O2)
- **Panel 4 (Continuum)** — smooth Moon/zodiacal + diffuse continuum, log scale

Each panel carries its own legend.  The individual-spectrum alpha is set
automatically as ``min(0.7, 3/√N)`` so traces remain legible regardless of
how many spectra are in the file.

**Usage**::

    # Full wavelength range
    XSkySepIvan_eval.py Sky_WHAM_south_08_ivan.fits

    # Specific window
    XSkySepIvan_eval.py Sky_WHAM_south_08_ivan.fits 9300 9600

    # Limit random sample overlay
    XSkySepIvan_eval.py Sky_WHAM_south_08_ivan.fits 6000 7000 -num 10

**Output:**

An HTML file named ``<stem>_<wmin>_<wmax>.html`` (e.g.
``Sky_WHAM_south_08_ivan_9300_9600.html``) that can be opened in any
browser for interactive zoom, pan, and hover inspection.

palace_make_mask.py
^^^^^^^^^^^^^^^^^^^

Builds a sky-line contamination mask across the full LVM wavelength range
(3600-9800 Å) using the PALACE (Paranal Airglow Line And Continuum Emission,
Noll et al. 2024) sky emission model.  Four components are rendered onto the
LVM wavelength grid: OH vibrational-rotational bands, OI recombination lines,
atomic forbidden/permitted lines (NaI, KI, [NI], OI), and the O2 A-band.
Each component is normalised to its own peak before summing so that no single
family dominates the mask.  The combined model is scaled to the observed sky
spectrum via a least-squares fit to bright OH pixels in the Z arm, and pixels
where the predicted contamination exceeds a user-specified threshold are
flagged as unusable for continuum fitting.

**Usage**::

    palace_make_mask.py fits_file palace_dir [--threshold T] [--plot] ...

**Arguments:**

fits_file
    LVM XCframe FITS file (provides WAVE, sky spectrum, and LSF).

palace_dir
    Path to the ``palace/PMD`` directory containing the PALACE data files.

**Key options:**

--threshold T
    Contamination threshold in FACTOR-scaled flux units.  Lower values give a
    stricter mask.  Default is 0.01 (= 1×10⁻¹⁶ erg s⁻¹ cm⁻² Å⁻¹ with the
    default FACTOR of 10¹⁴).  The tradeoff between mask strictness and the
    number of clean pixels available for continuum fitting is the primary
    tuning parameter.

--plot
    Display the diagnostic plot interactively (it is always saved as a PNG).

**Output:**

A FITS file (``<stem>_mask.fits``) containing:

- ``WAVE`` — wavelength array (Å)
- ``SKY`` — median observed sky spectrum (FACTOR-scaled)
- ``CONTINUUM`` — sky spectrum with contaminated pixels set to NaN
- ``MASK`` — boolean mask (1 = clean, 0 = contaminated)

A PNG diagnostic plot (``<stem>_mask.png``) showing all three
spectrograph arms on a log flux scale with the PALACE model, threshold line,
the full sky spectrum, and the clean continuum pixels highlighted.

GetSky_from_CFrame_sum.py
^^^^^^^^^^^^^^^^^^^^^^^^^

Inventories which sky fields have been observed in an XCframe summary file
and, optionally, extracts all sky spectra for a chosen field into a single
FITS file for stacking or modelling.

LVM records both an east (SKY_EAST) and west (SKY_WEST) sky telescope
pointing for each science exposure.  The DRP_ALL table inside an XCframe
summary file lists the field names in the ``skye_name`` and ``skyw_name``
columns.  This program has two operating modes:

*Summary mode* (no source name given): counts how many times each sky field
appears across both sky telescopes, prints the top N entries sorted by
observation count, and writes the full count table to
``sky_summary_<stem>.fits``.

*Extraction mode* (source name given): selects all rows where the field name
matches, extracts the corresponding spectra from SKY_EAST and SKY_WEST,
merges paired metadata columns (``skye_ra``/``skyw_ra`` → ``ra``, etc.) into
single columns for the relevant telescope, removes science-pointing columns,
and writes the result to ``Sky_<source_name>.fits``.

**Usage**::

    # Inventory sky fields
    GetSky_from_CFrame_sum.py fits_file

    # Extract spectra for one field
    GetSky_from_CFrame_sum.py fits_file source_name [--output PATH]

**Output (extraction mode):**

A FITS file containing:

- ``WAVE`` — wavelength array (Å)
- ``FLUX`` — sky spectra, shape (N_obs, N_pix), one row per observation
- ``DRP_ALL`` — metadata table with merged sky columns and a ``tel`` column
  indicating which sky telescope (SKY_EAST or SKY_WEST) each row came from

GetSkyCont.py
^^^^^^^^^^^^^

Fits a smooth two-component B-spline continuum to LVM sky spectra, separating
the continuum into a MOON/zodiacal component (B-splines modulated by a solar
spectrum) and a DIFFUSE component (plain B-splines for airglow continuum).
The fit is performed only on pixels flagged as line-free by a
``palace_make_mask.py`` mask; the result is evaluated at all wavelengths so
the continuum interpolates smoothly across emission-line regions.

The solar spectrum is read from a high-resolution reference file
(Meftah et al., LATMOS); its Fraunhofer absorption structure is kept fixed
while the B-spline envelope adjusts the colour and amplitude of the
Moon-reflected and zodiacal-light contribution.

**Usage**::

    # Sky-file mode (output of GetSky_from_CFrame_sum.py):
    GetSkyCont.py sky_file.fits -mask mask.fits [row_no ...] [-delta N]

    # XCframe / XSFrame summary file mode:
    GetSkyCont.py xframe.fits ext -mask mask.fits [row_no ...] [-delta N]

**Arguments:**

sky_file.fits
    Sky_<name>.fits produced by ``GetSky_from_CFrame_sum.py``.

xframe.fits, ext
    XCframe or XSFrame summary FITS file and the extension to read
    (e.g. ``SKY_EAST``, ``SKY_WEST``, ``FLUX``).

row_no
    Zero or more 0-based row indices.  If omitted all rows are processed
    (subject to ``-delta``).

**Key options:**

-mask file
    (required) mask FITS file from ``palace_make_mask.py``
    (MASK extension: 1 = clean, 0 = line-affected).

-delta N
    Process every N-th row (0, N, 2N, ...) instead of all rows.

-kstep N
    B-spline knot spacing in Angstroms (default 100 Å, giving ~66 basis
    functions per component).

-out outroot
    Set output filename root.

**Output:**

A FITS file (``<stem>_cont.fits``) containing:

- ``WAVE`` — wavelength array (Å)
- ``FLUX`` — input sky spectra (N_obs × N_pix)
- ``CONT`` — total continuum = MOON + DIFFUSE
- ``MOON`` — B-spline × solar component (Moon/zodiacal light)
- ``DIFFUSE`` — plain B-spline component (diffuse airglow continuum)
- ``RESID`` — residual FLUX − CONT (line emission isolated from continuum)
- ``MASK`` — boolean clean-pixel mask from the input palace_mask file
- ``DRP_ALL`` — observation metadata table

GetSkyCont_eval.py
^^^^^^^^^^^^^^^^^^

Evaluates the continuum fit produced by ``GetSkyCont.py`` by writing a single
HTML file containing three interactive Plotly figures.  Also updates the
DRP_ALL table in the input FITS file with per-spectrum fit-quality statistics.

**Figure 1 — four-panel spectral overview:**

- **Panel 1 (Flux + Continuum, log)** — observed sky spectra with the fitted
  total continuum median overlaid in red.  A green line near the bottom of
  the panel marks the wavelengths included in the fit (gaps where sky lines
  were masked); grey vertical bands shade the excluded regions.
- **Panel 2 (Total Continuum, log)** — the CONT band (10th–90th percentile
  and median) showing the overall level and smoothness of the B-spline model.
- **Panel 3 (Components, log)** — MOON (red) and DIFFUSE (orange) components
  on the same y-axis as Panel 2, with the total CONT median in purple.  The
  relative amplitude of the components shows how much of the continuum is
  Moon/zodiacal versus diffuse airglow.
- **Panel 4 (Residual, linear)** — FLUX − CONT; should be near zero in clean
  regions and show sky-line emission in the masked regions.

**Figure 2 — per-arm residual histograms:**

Three panels (Blue 3600–5900 Å, Red 5900–7600 Å, NIR 7600–9800 Å) each show
the distribution of all residual flux values — across every spectrum and every
clean (unmasked) pixel in that arm.  The x-axis is residual flux; the y-axis
is count N.  A Gaussian with center = median and σ = NMAD is overlaid in
black.  An annotation box reports N, the median, NMAD, skewness, and the
10th/90th percentiles.  The histogram range is clipped to median ± 5·NMAD to
suppress extreme outliers.

**Figure 3 — per-spectrum fit quality:**

- **Row 1** — three scatter plots (Blue, Red, NIR) of per-spectrum median
  residual (x-axis) vs NMAD (y-axis).  Each point is one spectrum; hovering
  shows the spectrum row index and its median and NMAD values.  Axis limits
  are clipped to the 2nd–98th percentile range so extreme outliers do not
  compress the scale for the bulk of the spectra; any off-scale spectra are
  noted in an annotation box.  A dashed vertical line marks x = 0 (ideal
  median); a dotted horizontal line marks the ensemble NMAD for reference.
- **Row 2** — NMAD vs original spectrum number (log y-scale) for all three
  arms on one panel, so that clusters of temporally adjacent poor fits are
  immediately visible.  Dotted horizontal lines mark the ensemble NMAD per arm.

**Per-spectrum statistics written to DRP_ALL:**

For each arm (``blue``, ``red``, ``nir``) and for all arms combined
(``all``), four float32 columns are added or updated in the DRP_ALL table:

+---------------------+---------------------------------------------+
| Column              | Description                                 |
+=====================+=============================================+
| ``resid_med_<arm>`` | median residual in clean pixels             |
+---------------------+---------------------------------------------+
| ``resid_nmad_<arm>``| NMAD of residuals in clean pixels           |
+---------------------+---------------------------------------------+
| ``resid_rms_<arm>`` | RMS of residuals in clean pixels            |
+---------------------+---------------------------------------------+
| ``resid_skew_<arm>``| skewness of residuals in clean pixels       |
+---------------------+---------------------------------------------+

**Usage**::

    # Full wavelength range
    GetSkyCont_eval.py Sky_WHAM_south_08_cont.fits

    # Specific window
    GetSkyCont_eval.py Sky_WHAM_south_08_cont.fits 6000 7000

    # Limit the random sample overlay
    GetSkyCont_eval.py Sky_WHAM_south_08_cont.fits 6000 7000 -num 10

**Output:**

An HTML file named ``<stem>_<wmin>_<wmax>.html`` containing all three
figures, viewable in any browser with interactive zoom, pan, and hover
inspection.  The input FITS file is updated in place with the DRP_ALL
statistics columns.


XCframe Sky Subtraction — Alternative Methods and Evaluation
-------------------------------------------------------------

This group of scripts performs per-spectrum sky subtraction directly on
XCframe summary FITS files, providing five independent algorithms that can be
run and compared side by side using ``SkySub_eval.py``.  All five subtraction
scripts produce output FITS files with WAVE, FLUX (sky-subtracted), SKY, and
DRP_ALL extensions (SkySepESO adds MOON/ZODI/DIFFUSE; see below).  DRP_ALL
carries a QA_FLAGS column recording per-row quality issues.

**Common QA flag bits**

======  ==========  ===========================================================
0x01    NANDATA     NaN or inf found in the input flux or sky data.
0x02    ZEROSKY     Sky line vector is all-zero; scale factor is unreliable.
0x04    POORFIT     Continuum fit is poorly conditioned (SkySubOrig, SkySubDev1).
0x04    MODELFAIL   ESO sky model (local SM-01 binary + SkyCalc web-service
                    fallback) both failed for this row (SkySepESO only).
0x08    FAILED      Row raised an exception; spectrum filled with NaN.
======  ==========  ===========================================================

**Common DRP_ALL columns beyond QA_FLAGS**

- ``mjd`` — precise MJD recomputed from the ``obstime`` ISO timestamp string
  via ``SkySubOrig.obstime_to_mjd()``, replacing the truncated-integer MJD
  carried through from the input file (all five scripts, plus SkySubSci.py/
  SummarizeSciSky.py below).
- ``LINE_SCALE`` — the bisection factor *r* used to scale the sky lines
  (SkySubOrig, SkySubDev1, SkySepESO only; SkySubDev2 applies no scale factor,
  and SkySubDrp wraps an external lvmdrp routine that does not expose its
  internal factor).
- ``SCI_MED_<arm>``/``SCI_NMAD_<arm>``/``SCI_RMS_<arm>``/``SCI_SKEW_<arm>`` and
  the equivalent ``SKY_*`` columns (arm = ``B``, ``R``, ``Z``) — per-arm
  continuum-fit-quality statistics (SkySubOrig, SkySubDev1, SkySepESO only;
  requires ``sky_mask.fits``, searched for automatically).  These test the
  continuum fit *itself* against the raw, pre-subtraction science and sky
  spectra, independent of the line-scaling step — see
  ``GetSkyCont.arm_continuum_stats()`` and the ``SkySub_eval.py``
  "Continuum Separation" figures below for how they are used.

SkySubOrig.py
^^^^^^^^^^^^^

Perform per-spectrum sky subtraction on an XCframe file using a degree-4
polynomial continuum and a custom bisection search for the sky-line scale
factor.

**Usage**::

    SkySubOrig.py [-method METHOD] [-delta N] [-out ROOT] filename

**Arguments:**

filename
    XCframe FITS file to process.

**Options:**

-method METHOD
    Sky subtraction method:

    - ``nearest`` — subtract the nearest sky telescope spectrum (continuum
      and lines from the same telescope).
    - ``farthest`` — subtract the farthest sky telescope spectrum.
    - ``farlines_nearcont`` — scale sky emission lines from the far sky
      telescope and continuum from the near sky telescope (default).

-delta N
    Process every N-th row; useful for quick tests (default: 1 = all rows).

-out ROOT
    Output filename root.  Default: ``<stem>_<method>``.

**Description:**

Each spectrum is decomposed into a polynomial continuum (degree 4, fitted
with sigma clipping on clean pixels) and line residuals.  Near and far sky
telescopes are identified from RA/Dec separations in DRP_ALL.  A global
line scale factor *r* is found by bisection minimising::

    sum |sci_lines × (sci_lines − r × sky_lines)| / ‖sky_lines‖²

For ``farlines_nearcont`` the sky model is::

    sky = cont_near + r × lines_far

The sky is then subtracted: ``flux_out = flux_sci − sky``.  Only the sky
lines are scaled by *r*; the continuum (``cont_near``) is used exactly as
fitted, with no additional scaling.

**Output:**

A FITS file ``<ROOT>.fits`` with extensions WAVE, FLUX (sky-subtracted),
SKY, and DRP_ALL (with QA_FLAGS, ``mjd``, ``LINE_SCALE``, and per-arm
continuum-fit-quality columns — see "Common DRP_ALL columns" above).

**See Also:** :doc:`api/SkySubOrig/index`


SkySubDrp.py
^^^^^^^^^^^^

Perform per-spectrum sky subtraction using the lvmdrp
``create_skysub_spectrum`` routine.

**Usage**::

    SkySubDrp.py [-method METHOD] [-delta N] [-out ROOT] filename

**Arguments:**

filename
    XCframe FITS file to process.

**Options:**

-method METHOD
    Sky subtraction method: ``nearest`` | ``farlines_nearcont`` (default).

-delta N
    Process every N-th row (default: 1).

-out ROOT
    Output filename root.  Default: ``<stem>_drp_<method>``.

**Description:**

For each row a minimal SCI/SKYE/SKYW BinTable HDUList is constructed from
the per-row FLUX, SKY_EAST, and SKY_WEST spectra together with RA/Dec
information from DRP_ALL::

    PRIMARY  (empty)
    SCI      BinTable  WAVE, FLUX, ERROR   RA/DEC from sci_ra/sci_dec
    SKYE     BinTable  WAVE, FLUX, ERROR   RA/DEC from skye_ra/skye_dec
    SKYW     BinTable  WAVE, FLUX, ERROR   RA/DEC from skyw_ra/skyw_dec

Flux errors are estimated as ``sqrt(|flux|)`` since XCframe files do not
carry IVAR.  These errors affect only the internally propagated sky error;
the sky model itself does not depend on them.  ``create_skysub_spectrum``
then selects the sky telescope and computes the sky model according to the
chosen method.

Requires the ``lvmdrp`` conda environment (``lvmdrp26``).

**Output:**

A FITS file ``<ROOT>.fits`` with extensions WAVE, FLUX (sky-subtracted),
SKY, and DRP_ALL (with QA_FLAGS and precise ``mjd``).  ``create_skysub_spectrum``
only returns the sky spectrum and its error, not an internal scale factor, so
no ``LINE_SCALE`` or continuum-fit-quality columns are written here.

**See Also:** :doc:`api/SkySubDrp/index`


SkySubDev1.py
^^^^^^^^^^^^^

Perform per-spectrum sky subtraction using the B-spline continuum
separation from ``GetSkyCont.py`` instead of the polynomial fit used by
SkySubOrig.  The bisection scale search and sky assembly are otherwise
identical to SkySubOrig.

**Usage**::

    SkySubDev1.py [-method METHOD] [-delta N] [-kstep N]
                  [-out ROOT] [-mask mask.fits] filename

**Arguments:**

filename
    XCframe FITS file to process.

**Options:**

-mask file
    palace_mask FITS file produced by ``palace_make_mask.py``
    (MASK extension: 1 = clean, 0 = sky-line affected).
    If omitted, ``sky_mask.fits`` is searched for in the current directory
    and then in the ``lvm_ksl data/`` directory.

-method METHOD
    ``nearest`` | ``farlines_nearcont`` (default).

-delta N
    Process every N-th row (default: 1).

-kstep N
    B-spline knot spacing in Angstroms (default: 100).

-out ROOT
    Output filename root.  Default: ``<stem>_dev1_<method>``.

**Description:**

A two-component B-spline design matrix (DIFFUSE + MOON, same as
``GetSkyCont.py``) is built once from the full wavelength grid and the
palace mask, then reused for every row.  For each row:

1. Science and sky spectra are decomposed into continuum (non-negative
   least squares on clean pixels) and line residuals.
2. Near/far sky telescopes are identified from RA/Dec separations.
3. A global line scale factor *r* is found by bisection (same objective
   function as SkySubOrig).
4. The sky model is assembled and subtracted::

       farlines_nearcont:  sky = cont_near + r × lines_far
       nearest:            sky = cont_near + r × lines_near

   ``cont_near`` is always the continuum used in the final SKY (both
   methods), used exactly as fitted with no additional scaling.
5. If ``sky_mask.fits`` is available, per-arm continuum-fit-quality stats
   (``SCI_*``/``SKY_*`` columns) are computed from the raw pre-subtraction
   science spectrum (``lines_sci``) and the near-sky spectrum
   (``lines_near``), via ``GetSkyCont.arm_continuum_stats()``.

**Output:**

A FITS file ``<ROOT>.fits`` with extensions WAVE, FLUX (sky-subtracted),
SKY, and DRP_ALL (with QA_FLAGS, precise ``mjd``, ``LINE_SCALE``, and
per-arm continuum-fit-quality columns — see "Common DRP_ALL columns" above).

**See Also:** :doc:`api/SkySubDev1/index`


SkySubDev2.py
^^^^^^^^^^^^^

Perform per-spectrum sky subtraction using the PALACE spectral
decomposition (``XSkySepIvan.py``), without any additional scaling.

**Usage**::

    SkySubDev2.py [-method METHOD] [-delta N] [-lsf FWHM]
                  [-out ROOT] filename

**Arguments:**

filename
    XCframe FITS file to process.

**Options:**

-method METHOD
    ``nearest`` | ``farlines_nearcont`` (default).

-delta N
    Process every N-th row (default: 1).

-lsf FWHM
    LSF FWHM in Angstroms (default: 1.3).

-out ROOT
    Output filename root.  Default: ``<stem>_dev2_<method>``.

**Description:**

A PALACE decomposer (``SkyDecomp``) is built once from the full wavelength
grid and LSF, then reused for every row.  For each row:

1. Near and far sky telescopes are identified from RA/Dec separations.
2. Each sky spectrum is decomposed by PALACE into emission-line and
   continuum components::

       LINES = oh + atom + orc + o2
       CONT  = moon + diffuse

3. The sky model is assembled without any scale factor::

       farlines_nearcont:  sky = CONT_near + LINES_far
       nearest:            sky = CONT_near + LINES_near

4. ``flux_out = flux_sci − sky``.

Requires the PALACE library (``lvmsky/skysub/sky_decomp``) and the
PALACE data files; paths are taken from ``XSkySepIvan.py``.

**Output:**

A FITS file ``<ROOT>.fits`` with extensions WAVE, FLUX (sky-subtracted),
SKY, and DRP_ALL (with QA_FLAGS and precise ``mjd``).  No scale factor is
applied here (see docstring: "no scaling"), so there is no ``LINE_SCALE``
column or continuum-fit-quality columns to record.

**See Also:** :doc:`api/SkySubDev2/index`


SkySepESO.py
^^^^^^^^^^^^

Perform per-spectrum sky subtraction using the **ESO Sky Model** to separate
the sky continuum into its physical MOON, ZODI, and DIFFUSE components,
using the same bisection line-scale search as SkySubOrig.  Ported from a
prototype developed in ``lvm_sky2506`` (``SkySepMod.py`` +
``SkySubModDev250624.ipynb``).

**Usage**::

    SkySepESO.py [-method METHOD] [-delta N] [-out ROOT] filename

**Arguments:**

filename
    XCframe FITS file to process.

**Options:**

-method METHOD
    ``nearest`` | ``farthest`` | ``farlines_nearcont`` (default) — same
    meaning as SkySubOrig.

-delta N
    Process every N-th row (default: 1).  Each row requires 2–3 ESO sky
    model fetches (science fiber, plus one or two sky fibers depending on
    method), so this script is far slower per row than the other four —
    measured at ~2.5s/row for ``farlines_nearcont`` (3 model fetches + 3
    continuum fits), i.e. roughly an hour for a full ~1750-fiber exposure.
    Use ``-delta`` liberally for quick tests.

-out ROOT
    Output filename root.  Default: ``<stem>_eso_<method>``.

**Description:**

For each row, the ESO sky model is fetched for the relevant fiber's
coordinates and observation time (``SkyModelObs.do_one``, using the local
ESO SM-01 binary; falls back automatically to the SkyCalc web service via
``SkyCalcObs.py`` if the local model call fails — e.g. because the model
rejects the target/Moon geometry).  The model gives MOON/ZODI/DIFFUSE
continuum templates, which are fit to the observed flux as a non-negative,
iteratively-downweighted 3-component linear combination
(``CONT = a·MOON + b·ZODI + c·DIFFUSE``) so that sky/airglow line pixels do
not bias the continuum fit upward.  The resulting line residual is then
scaled against the science spectrum's own line residual by the same
bisection search used in SkySubOrig, and the sky is subtracted::

    sky = sky_CONT + r × sky_LINES
    flux_out = flux_sci − sky

Only the sky lines are scaled by *r*; ``sky_CONT`` (whichever sky fiber's
ESO-model fit produced it — near for ``nearest``/``farlines_nearcont``, far
for ``farthest``) is used exactly as fitted.

Each ESO-model fetch writes a small per-call FITS file
(``SkyM_*.fits``/``SkyC_*.fits``) to the current working directory; this
file is read and deleted immediately, so nothing accumulates on disk across
a run.

**Output:**

A FITS file ``<ROOT>.fits`` with extensions WAVE, FLUX (sky-subtracted),
SKY, **MOON**, **ZODI**, **DIFFUSE** (the *scaled* ESO-model components of
whichever sky fit produced SKY's continuum term — they sum to SKY's
continuum; the remaining line term is ``SKY − (MOON+ZODI+DIFFUSE)``), and
DRP_ALL.  DRP_ALL carries QA_FLAGS (with the additional ``MODELFAIL`` bit
0x04 — both the local ESO model and the SkyCalc fallback failed for this
row), precise ``mjd``, ``LINE_SCALE``, the fitted ESO-model coefficients
(``SCI_MOON``/``SCI_ZODI``/``SCI_DIFFUSE``, ``SKY_MOON``/``SKY_ZODI``/
``SKY_DIFFUSE``), per-arm continuum-fit-quality columns (see "Common
DRP_ALL columns" above), and ``ERROR_MSG`` — a per-row failure-reason
string (truncated to 200 characters) for any row with a non-zero QA_FLAGS
value.  A summary of flagged rows is printed at the end of the run and also
written to ``<ROOT>_errors.txt``.

**Requirements:**

The local ESO SM-01 sky model binary, gated by the ``ESO_SKY_MODEL``
environment variable (see ``SkyModelObs.py`` above); without it, every row
falls back to the SkyCalc web service (slower, needs network access).

**See Also:** :doc:`api/SkySepESO/index`


SkySub_eval.py
^^^^^^^^^^^^^^

Evaluate sky subtraction quality for one or more output FITS files
produced by SkySubOrig, SkySubDrp, SkySubDev1, SkySubDev2, or SkySepESO.
When multiple files are given they are overlaid in the same figures for
direct method comparison.

**Usage**::

    SkySub_eval.py [wmin wmax] [-num N] [-out outroot] filename [filename ...]

**Arguments:**

filename
    One or more SkySub output FITS files to evaluate.

wmin, wmax
    Wavelength range in Angstroms for the spectral overview panels
    (defaults: 3600 and 9800).

**Options:**

-num N
    Overlay N randomly selected individual spectra on the band panels
    (default: 20; 0 = band only).

-out outroot
    Combine all files into one HTML file with this root.  Without ``-out``,
    each file produces its own ``<stem>_eval.html``.

**Description — HTML output:**

The output HTML file is titled "Sky Subtraction Quality Check" and is
organized into two named sections: **Sky Line Subtraction** (Figures 1–4
plus a statistics table) and **Continuum Separation** (Figures 5–6 plus a
second statistics table).  In total there are six interactive Plotly
figures and two inline statistics tables.

*Figure 1 — spectral overview (3 panels, linear scale):*

All three panels share a y-range of −1×10⁻¹⁴ to 1×10⁻¹³
erg s⁻¹ cm⁻² Å⁻¹, chosen to make sky-subtracted residuals visible.
Panel 1 shows FLUX + SKY (before subtraction), Panel 2 shows FLUX
(after subtraction), and Panel 3 shows the SKY model.  Each panel shows
the median and 10th/90th percentile band together with up to N individual
spectra in light grey.  One colour per input file.  This figure sits above
the "Sky Line Subtraction" section header (it isn't part of either named
section).

*Figure 2 — residual histograms:*

For each of three diagnostic windows ([OI] 5577 Å: 5560–5594 Å,
[OI] 6300 Å: 6280–6320 Å, IR OH: 9300–9500 Å) the distribution of the
per-pixel **HF (high-frequency, continuum-subtracted) residual** — not
raw FLUX — is plotted with a Gaussian overlay.  Using the HF residual
means leftover continuum in the window does not bias the reported
median/NMAD away from genuine sky-line residuals.  A statistics box (N,
median, NMAD, skewness) appears inside each panel; skewness is computed
on the same median ± 5·NMAD clipped range shown in the plot, so it isn't
dominated by a handful of outliers invisible in the display.

*Figure 3 — diagnostic window median spectra:*

For each diagnostic window a wider search region (5400–5750,
6100–6500, 9000–9800 Å) is plotted showing the median and 10th/90th
percentile band for the original (dotted) and sky-subtracted (solid)
spectra.  Red shading marks the diagnostic (signal) window; green shading
marks mask-selected sky-line-free pixels used to estimate the noise floor.

*Statistics table:*

An HTML table between Figure 3 and Figure 4 reports per-file per-window
the noise floor, sky-line RMS before and after subtraction, and
percentiles of the HF RMS ratio (see below).  The same information is
printed to the terminal.

*Figure 4 — HF RMS ratio per spectrum:*

For each diagnostic window the noise-corrected HF RMS ratio is plotted
against spectrum index, or against MJD (recomputed precisely from
``obstime`` rather than the truncated-integer ``mjd`` column) when there
are more than 100 spectra and every overlaid file has that information —
otherwise all files fall back to spectrum index so every trace stays on a
common axis.  Points are plotted as markers (not connected lines).  The
three panels share a linked x-axis so zooming one panel pans all three.
See the algorithm description below.


**HF RMS quality metric — algorithm:**

The central diagnostic is a per-spectrum, per-window noise-corrected
high-frequency (HF) RMS ratio.  The following steps are applied
independently to both the original spectra (FLUX + SKY) and the
sky-subtracted spectra (FLUX).

*Step 1 — high-frequency residuals.*

A smooth estimate of the underlying continuum is subtracted from each
spectrum::

    hf(λ) = flux(λ) − smooth(λ)

The smooth is a Gaussian-weighted running average with σ = 50 pixels
(≈ 25 Å at the LVM pixel scale of 0.5 Å pixel⁻¹).  Sky-line pixels
(mask = 0 from the palace mask) are excluded by a weighted convolution::

    smooth(λ) = Σ_λ' [ flux(λ') · w(λ') · G_σ(λ−λ') ]
              / Σ_λ' [ w(λ') · G_σ(λ−λ') ]

where w = 1 for mask-clean pixels and w = 0 for sky-line pixels.
This prevents bright sky lines from leaking into the smooth estimate and
artificially reducing the HF amplitude in the diagnostic window.

After this step, ``hf`` retains only structure narrower than ≈25 Å —
the scale of individual sky emission lines — while broad continuum
mismatches from polynomial or B-spline fitting are suppressed.

*Step 2 — why RMS, not NMAD.*

Sky lines are spatially sparse.  [OI] 5577 spans roughly 3–4 pixels in
the 34-pixel diagnostic window; individual OH lines in the IR cover a
similarly small fraction.  NMAD (the median absolute deviation) is
dominated by the majority of clean, noise-floor pixels and is nearly
insensitive to a small number of bright outliers.  RMS responds to the
squared amplitude, so even 2–4 bright sky-line pixels contribute
proportionally to their intensity — exactly the signal we want to track.

*Step 3 — noise floor subtraction.*

Even after perfect sky subtraction the HF residuals carry photon and
detector noise.  This noise contributes to the RMS measured in the
diagnostic window.  A nearby background region — mask-selected
line-free pixels (mask = 1) within the broader search range but outside
the diagnostic window — provides an independent noise estimate
``rms_bg``.  The sky-line contribution is then isolated in quadrature::

    sky_rms = sqrt( max(0, rms_diag² − rms_bg²) )

For a perfect subtraction ``sky_rms`` approaches zero regardless of the
noise level.  The assumption is that the noise is approximately stationary
across the search region, which is valid as long as the continuum level
does not change drastically within a few hundred ångströms.

*Step 4 — the ratio and its interpretation.*

For each spectrum::

    ratio = sky_rms_sub / sky_rms_orig

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - ratio
     - Meaning
   * - 0
     - Sky lines completely removed.
   * - 0 – 0.5
     - Substantial improvement; more than half the sky-line power removed.
   * - ≈ 1
     - No improvement; subtraction left sky-line power unchanged.
   * - > 1
     - Sky lines made worse; subtraction introduced new residuals.
   * - NaN
     - ``sky_rms_orig`` ≈ 0; no detectable sky-line power in the original
       spectrum — excluded from statistics.

*Step 5 — summary statistics.*

The statistics table reports the following per file and per window:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Description
   * - ``noise_med``
     - Median of ``rms_bg`` across all spectra — the noise floor
       (erg s⁻¹ cm⁻² Å⁻¹).  Tells you how sensitive the measurement is.
   * - ``sky_orig_med``
     - Median ``sky_rms`` in the original spectra — the typical sky-line
       amplitude before subtraction.
   * - ``sky_sub_med``
     - Median ``sky_rms`` after subtraction — the residual sky-line amplitude.
       Should be much smaller than ``sky_orig_med`` for a good method.
   * - ``p50 ratio``
     - Median ratio across all spectra; characterises typical performance.
   * - ``p90 ratio``
     - 90th-percentile ratio; indicates how the worst 10% of spectra behave.
   * - ``p95 ratio``
     - 95th-percentile ratio; the tail of poor performance.
   * - ``frac < 0.5``
     - Fraction of spectra with ratio < 0.5, i.e. where sky-line power was
       more than halved.  A simple pass/fail rate at this threshold.


**Continuum Separation section — Figures 5 & 6 and the second statistics table:**

Unlike the HF RMS metric above (which targets sky *line* residuals), this
section evaluates *continuum* separation quality, using three
spectrograph-arm bands with the B/R/Z overlap zones and outer edges
excluded — B (3650–5775 Å), R (5800–7520 Å), Z (7570–9600 Å).  Both
figures share the same three-row layout and row order:

- **Row 1 (SCI)** — DRP_ALL['SCI_MED_<arm>']: the per-spectrum median of
  (raw science flux − sci_cont) in clean (sky-line-free) pixels of the
  *raw, pre-subtraction* science spectrum.  This is the actual
  science-side continuum-fit-quality test.
- **Row 2 (SKY)** — DRP_ALL['SKY_MED_<arm>']: the same, but for the raw
  sky-telescope spectrum whose continuum fit was actually used in SKY.
  This is the actual sky-side continuum-fit-quality test, and is the one
  that matters directly for the final result's continuum level.
- **Row 3 (Subtracted)** — the per-spectrum median of the *raw*
  sky-subtracted FLUX (not the HF residual used by Figure 2) in clean
  pixels.  This is the final, post-subtraction leftover signal.  It does
  **not** test continuum-fit quality — the science-side continuum fit is
  only ever used to derive the bisection line-scale target and never
  enters the subtracted result — so this row instead reflects real source
  (stellar) continuum entangled with any net error in the sky-side
  continuum estimate.  Comparing this row across overlaid methods on the
  same input isolates sky-continuum quality specifically, since real
  source continuum is identical across methods.

Rows 1/2 require the ``SCI_MED_<arm>``/``SKY_MED_<arm>`` DRP_ALL columns
already present in the input file (written by SkySubOrig.py, SkySubDev1.py,
or SkySepESO.py); a file without them (SkySubDev2.py/SkySubDrp.py output,
or an older file predating this feature) simply leaves those panels empty.

*Figure 5* histograms all three rows (pooled/distribution across spectra,
same N/median/NMAD/skew annotation style as Figure 2).  *Figure 6* plots
the same three quantities per spectrum against spectrum index or MJD
(same x-axis logic as Figure 4), so specific bad exposures/fibers are
visible rather than just an aggregate number.  Both figures use a
colour-coded suptitle above the plot grid (rather than a floating legend,
which would otherwise sit on top of a subplot in the 3×3 grid).

The second statistics table (screen + HTML) reports, per file/arm/kind
(SCI, SKY, Subtracted), the N/median/NMAD/skew of the row's distribution.

**DRP_ALL write-back:** the per-spectrum Subtracted-row statistics
(``resid_med_<arm>``, ``resid_nmad_<arm>``, ``resid_rms_<arm>``,
``resid_skew_<arm>``) are written back into each evaluated file's own
DRP_ALL table, in place — the same convention used by
``GetSkyCont_eval.py`` for its own fit-quality columns.  (The SCI_MED/
SKY_MED columns themselves are not written by this script — they already
exist in the input file, written by whichever SkySub* script produced it.)

**Output:**

An HTML file named ``<stem>_eval.html`` (one per input file without
``-out``; a single combined file when ``-out outroot`` is given).  Any
evaluated file with a DRP_ALL table is updated in place with the
Subtracted-row statistics above; no backup is created.

**See Also:** :doc:`api/SkySub_eval/index`


Comparing Sky Subtraction Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A typical workflow for running all five methods on a single XCframe file
and comparing the results::

    # 1. Build the sky-line mask (if not already present)
    palace_make_mask.py XCframe_file.fits /path/to/palace/PMD

    # 2. Run the five subtraction methods
    SkySubOrig.py  XCframe_file.fits
    SkySubDrp.py   XCframe_file.fits
    SkySubDev1.py  XCframe_file.fits -mask sky_mask.fits
    SkySubDev2.py  XCframe_file.fits
    SkySepESO.py   XCframe_file.fits -delta 50   # slow; use -delta for a quick look

    # 3. Evaluate and compare in a single HTML file
    SkySub_eval.py -out compare \
        XCframe_file_farlines_nearcont.fits \
        XCframe_file_drp_farlines_nearcont.fits \
        XCframe_file_dev1_farlines_nearcont.fits \
        XCframe_file_dev2_farlines_nearcont.fits \
        XCframe_file_eso_farlines_nearcont.fits

Open ``compare_eval.html`` in a browser.  Figure 1 overlays the median
spectra for all five methods; Figure 4 (Sky Line Subtraction section) shows
the per-spectrum HF RMS ratio for each diagnostic window; Figures 5/6
(Continuum Separation section) show per-arm continuum-fit quality for the
three methods that record it (SkySubOrig, SkySubDev1, SkySepESO) — together
these make it straightforward to identify which method best suppresses
sky lines *and* which best separates continuum from lines for the
observation.


Science-Fiber-Based Sky Estimation
------------------------------------

Unlike the XCframe methods above (which use the dedicated SKY_EAST/SKY_WEST
sky telescopes), these two scripts estimate a sky spectrum directly from the
science IFU itself.  At a typical LVM pointing most of the 1801 science
fibers see mostly sky rather than an astronomical source; ranking fibers by
sky-line-free continuum flux and averaging the faintest fibers gives a sky
proxy without needing a sky telescope pointing at all.  No scale factor is
applied to emission lines here, so the method is best suited to fields
where a genuinely sky-dominated tail of faint fibers exists (e.g. diffuse
or extended sources, not compact point sources filling the IFU).

Both scripts write output with WAVE/SCI/SKY/FLUX/DRP_ALL extensions
compatible with ``SkySub_eval.py`` (FLUX = sky-subtracted, SKY = sky
model), so they can be evaluated and compared alongside the four XCframe
methods above.

SkySubSci.py
^^^^^^^^^^^^

Estimate a sky spectrum for one or more LVM CFrame exposures directly from
the science fibers, without using the dedicated sky telescopes.

**Usage**::

    SkySubSci.py [-low PCT] [-high PCT] [-navg N] [-sigma S] [-maxiters K]
                [-mask FILE] [-stat median|mean] [-out ROOT]
                filename [filename ...]

**Arguments:**

filename
    One or more lvmCFrame FITS files.  One output row is written per file,
    in the order given.

**Options:**

-low PCT
    Percentile rank (0-100) of the faint/sky-like fiber (default 10).

-high PCT
    Percentile rank (0-100) of the bright/science-like fiber (default 90).

-navg N
    Number of fibers, ranked closest to ``-low``/``-high``, combined with a
    sigma-clipped (robust) mean (default 10; use 1 to reproduce the
    original single-fiber behaviour).

-sigma S
    Sigma-clipping threshold for the robust mean (default 3.0).

-maxiters K
    Sigma-clipping iteration limit (default 5).

-mask FILE
    palace_mask FITS file from ``palace_make_mask.py`` (WAVE/MASK
    extensions, MASK=1 means clean/sky-line-free).  If omitted,
    ``sky_mask.fits`` is searched for in the current directory, then in
    the ``lvm_ksl data/`` directory.

-stat STAT
    Statistic used to rank fibers by continuum flux: ``median`` (default)
    or ``mean``.

-out ROOT
    Output filename root.  Default:
    ``SkySubSci_<first_expnum>_<last_expnum>`` (or ``SkySubSci_<expnum>``
    for a single file).

**Description:**

For each exposure, science fibers are selected from SLITMAP (``scifib``,
imported from ``SummarizeCframe.py``).  Each fiber's continuum flux is
measured as the median (or mean) FLUX over sky-line-free pixels (from the
palace mask, resampled onto the file's own wavelength grid).  Fibers are
sorted by that continuum level; around each of the ``-low``/``-high``
percentile ranks, a window of the ``-navg`` fibers whose rank is closest to
that target is combined pixel-by-pixel with a sigma-clipped mean
(``astropy.stats.sigma_clipped_stats``).  The ``-low`` window is the sky
estimate; the ``-high`` window is a bright/science-like reference; their
difference is the sky-subtracted result.  Averaging several fibers per
percentile trades spatial locality in rank-space for lower per-pixel noise
(~1/√navg for clean pixels); sigma clipping keeps a fiber with a faint
source or a defect from dominating the average.

**Output:**

A FITS file with WAVE, SCI (robust mean of the high-percentile window,
reference only), SKY (robust mean of the low-percentile window, the sky
model), FLUX (SCI − SKY, sky-subtracted), and DRP_ALL (one row per
exposure: filename, expnum, exptime, obstime, mjd, fiber counts and IDs
used, RA/Dec, spectrograph IDs, and continuum flux for both windows)
extensions.

**See Also:** :doc:`api/SkySubSci/index`


SummarizeSciSky.py
^^^^^^^^^^^^^^^^^^

Drpall-driven, remote-friendly version of ``SkySubSci.py``: runs the same
science-fiber sky estimate over a range of exposure numbers selected from a
drpall table (as ``SummarizeCframe.py`` does) instead of an explicit file
list, so it can be run unattended over many exposures (e.g. at Utah).

**Usage**::

    SummarizeSciSky.py [-emin 900] [-ver 1.2.1] [-drp_all FILE]
                       [-low 10] [-high 90] [-navg 10] [-sigma 3.0]
                       [-maxiters 5] [-mask FILE] [-stat median|mean]
                       [-out ROOT] exp_start exp_stop [delta]

**Arguments:**

exp_start
    Starting exposure number.

exp_stop
    Stopping exposure number.

delta
    Process every delta-th exposure in range (default 1).

**Options:**

-emin N
    Minimum exposure time to include (default 900).

-ver VER
    DRP version, used to locate ``drpall-VER.fits`` (default 1.2.1).

-drp_all FILE
    Explicit drpall table to read instead of ``drpall-VER.fits`` (FITS, or
    ascii if the name contains ``txt``/``.tab``).

-low PCT, -high PCT, -navg N, -sigma S, -maxiters K, -mask FILE, -stat STAT
    Same meaning as in ``SkySubSci.py``.

-out ROOT
    Output filename root.  Default:
    ``SummarizeSciSky_<ver>_<exp_start>_<exp_stop>_<delta>``.

**Description:**

Exposure selection and file-path resolution follow the same pattern as
``SummarizeCframe.py``'s ``read_drpall``/``select_exps``/``find_top``, but
are re-implemented locally here (with an added ``-drp_all`` override)
rather than imported from ``SumCframe.py``, so this script has no
dependency on the optional ``dask`` package that ``SumCframe.py`` requires
for an unrelated function.  The per-exposure algorithm is identical to
``SkySubSci.py`` (also re-implemented locally rather than imported, so the
script is standalone).  Rather than building a fresh metadata table, the
calculated values (continuum flux, fiber IDs used, positions, spectrograph
IDs) are added as new columns directly onto the selected drpall rows, which
become the DRP_ALL extension of the output.

**Output:**

Same extension structure as ``SkySubSci.py`` (WAVE, SCI, SKY, FLUX,
DRP_ALL), with DRP_ALL being the selected drpall rows plus the added
columns.

**See Also:** :doc:`api/SummarizeSciSky/index`


Typical Workflows
-----------------

Evaluating DRP Sky Subtraction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Obtain SFrame files from the DRP
2. Run eval_sky.py to visualize residuals::

       eval_sky.py lvmSFrame-00012345.fits

3. Look for systematic residuals at sky line wavelengths

Testing Alternative Sky Subtraction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Start with CFrame files (before DRP sky subtraction)
2. Prepare files for SkyCorr::

       Prep4SkyCorr.py lvmCFrame-00012345.fits

3. Run SkyCorr::

       RunSkyCorr.py sci_mean.fits sky_mean.fits

4. Compare results with DRP sky subtraction

Comparing with Sky Models
^^^^^^^^^^^^^^^^^^^^^^^^^

1. Generate a theoretical sky for the observation::

       SkyCalcObs.py 81.5 -66.0 60000.5

2. Compare with observed sky from SKY_EAST or SKY_WEST telescopes
3. Identify discrepancies that may indicate calibration issues

Fitting and Evaluating the Sky Continuum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Build a palace line mask for the field::

       palace_make_mask.py XCframe_file.fits /path/to/palace/PMD

2. Collect sky spectra from repeated observations of the field::

       GetSky_from_CFrame_sum.py XCframe_file.fits Sky_WHAM_south_08

3. Fit the two-component B-spline continuum::

       GetSkyCont.py Sky_WHAM_south_08.fits -mask XCframe_file_mask.fits

4. Evaluate the fit interactively::

       GetSkyCont_eval.py skycont_Sky_WHAM_south_08.fits


Notes
-----

- Sky subtraction quality depends strongly on observing conditions
- The SKY_EAST and SKY_WEST telescopes point at different positions
  than the science telescope, which can lead to spatial sky variations
- SkyCorr can model and remove airglow emission lines more accurately
  than simple scaling methods
- Theoretical sky models are useful for identifying instrumental
  artifacts vs. real sky features


See Also
--------

- :doc:`summarize` - Tools for evaluating sky subtraction across many exposures
- :doc:`api/eval_sky/index` - API documentation
- :doc:`api/Prep4SkyCorr/index` - API documentation
- :doc:`api/RunSkyCorr/index` - API documentation
- :doc:`api/SkySub/index` - API documentation
- :doc:`api/SkyCalcObs/index` - API documentation
- :doc:`api/SkyModelObs/index` - API documentation
- :doc:`api/palace_make_mask/index` - API documentation
- :doc:`api/GetSky_from_CFrame_sum/index` - API documentation
- :doc:`api/XSkySepIvan/index` - API documentation
- :doc:`api/XSkySepIvan_eval/index` - API documentation
- :doc:`api/GetSkyCont/index` - API documentation
- :doc:`api/GetSkyCont_eval/index` - API documentation
- :doc:`api/SkySubOrig/index` - API documentation
- :doc:`api/SkySubDrp/index` - API documentation
- :doc:`api/SkySubDev1/index` - API documentation
- :doc:`api/SkySubDev2/index` - API documentation
- :doc:`api/SkySepESO/index` - API documentation
- :doc:`api/SkySub_eval/index` - API documentation
- :doc:`api/SkySubSci/index` - API documentation
- :doc:`api/SummarizeSciSky/index` - API documentation
- :doc:`summarize` - SummarizeCframe.py, whose drpall selection logic SummarizeSciSky.py mirrors
