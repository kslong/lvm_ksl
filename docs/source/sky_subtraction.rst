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
- Run alternative sky subtraction using ESO's SkyCorr tool
- Generate theoretical sky models for comparison
- Visualize sky residuals and identify problems


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
XCframe summary FITS files, providing four independent algorithms that can be
run and compared side by side using ``SkySub_eval.py``.  All four subtraction
scripts produce output FITS files with WAVE, FLUX (sky-subtracted), SKY, and
DRP_ALL extensions.  DRP_ALL carries a QA_FLAGS column recording per-row
quality issues.

**Common QA flag bits**

======  ==========  ===========================================================
0x01    NANDATA     NaN or inf found in the input flux or sky data.
0x02    ZEROSKY     Sky line vector is all-zero; scale factor is unreliable.
0x04    POORFIT     Continuum fit is poorly conditioned (SkySubOrig only).
0x08    FAILED      Row raised an exception; spectrum filled with NaN.
======  ==========  ===========================================================

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

The sky is then subtracted: ``flux_out = flux_sci − sky``.

**Output:**

A FITS file ``<ROOT>.fits`` with extensions WAVE, FLUX (sky-subtracted),
SKY, and DRP_ALL (with QA_FLAGS column).

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
SKY, and DRP_ALL (with QA_FLAGS column).

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

**Output:**

A FITS file ``<ROOT>.fits`` with extensions WAVE, FLUX (sky-subtracted),
SKY, and DRP_ALL (with QA_FLAGS column).

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
SKY, and DRP_ALL (with QA_FLAGS column).

**See Also:** :doc:`api/SkySubDev2/index`


SkySub_eval.py
^^^^^^^^^^^^^^

Evaluate sky subtraction quality for one or more output FITS files
produced by SkySubOrig, SkySubDrp, SkySubDev1, or SkySubDev2.  When
multiple files are given they are overlaid in the same figures for direct
method comparison.

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

The output HTML file contains four interactive Plotly figures and an
inline statistics table.

*Figure 1 — spectral overview (3 panels, linear scale):*

All three panels share a y-range of −1×10⁻¹⁴ to 1×10⁻¹³
erg s⁻¹ cm⁻² Å⁻¹, chosen to make sky-subtracted residuals visible.
Panel 1 shows FLUX + SKY (before subtraction), Panel 2 shows FLUX
(after subtraction), and Panel 3 shows the SKY model.  Each panel shows
the median and 10th/90th percentile band together with up to N individual
spectra in light grey.  One colour per input file.

*Figure 2 — residual histograms:*

For each of three diagnostic windows ([OI] 5577 Å: 5560–5594 Å,
[OI] 6300 Å: 6280–6320 Å, IR OH: 9300–9500 Å) the distribution of all
per-pixel sky-subtracted values is plotted with a Gaussian overlay.
A statistics box (N, median, NMAD, skewness) appears inside each panel.

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
against spectrum index.  The three panels share a linked x-axis so
zooming one panel pans all three.  See the algorithm description below.


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

**Output:**

An HTML file named ``<stem>_eval.html`` (one per input file without
``-out``; a single combined file when ``-out outroot`` is given).

**See Also:** :doc:`api/SkySub_eval/index`


Comparing Sky Subtraction Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A typical workflow for running all four methods on a single XCframe file
and comparing the results::

    # 1. Build the sky-line mask (if not already present)
    palace_make_mask.py XCframe_file.fits /path/to/palace/PMD

    # 2. Run the four subtraction methods
    SkySubOrig.py  XCframe_file.fits
    SkySubDrp.py   XCframe_file.fits
    SkySubDev1.py  XCframe_file.fits -mask sky_mask.fits
    SkySubDev2.py  XCframe_file.fits

    # 3. Evaluate and compare in a single HTML file
    SkySub_eval.py -out compare \
        XCframe_file_farlines_nearcont.fits \
        XCframe_file_drp_farlines_nearcont.fits \
        XCframe_file_dev1_farlines_nearcont.fits \
        XCframe_file_dev2_farlines_nearcont.fits

Open ``compare_eval.html`` in a browser.  Figure 1 overlays the
median spectra for all four methods; Figure 4 shows the per-spectrum HF
RMS ratio for each diagnostic window, making it straightforward to
identify which method best suppresses sky lines for the observation.


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
- :doc:`api/SkySub_eval/index` - API documentation
