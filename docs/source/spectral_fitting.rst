Spectral Fitting: local routines
================================

The lvm_ksl package provides tools for fitting emission lines in LVM
spectra. These tools can fit Gaussian profiles to standard emission
lines across all fibers in an RSS file, or perform more detailed fits
to individual spectra.


Overview
--------

The spectral fitting tools include:

- ``lvm_gaussfit.py`` - Fit standard emission lines across an RSS file
- ``sky_gaussfit.py`` - Fit nebular and airglow lines fiber-by-fiber in SFrame files
- ``lvm_line_profile.py`` - Compare Gaussian vs Moffat airglow line profiles on raw sky spectra
- ``lvm_double.py`` - Fit single or double Gaussian profiles to a line
- ``lvm_triple.py`` - Fit up to triple Gaussian profiles
- ``lvm_flux.py`` - Calculate fluxes using parameters from other fits

These tools use the ``lmfit`` package for robust non-linear fitting
with proper error estimation.


Standard Emission Line Fitting
------------------------------

lvm_gaussfit.py
^^^^^^^^^^^^^^^

The primary tool for emission line fitting. It fits Gaussian profiles
to a standard set of emission lines (H-alpha, [NII], [SII], etc.) for
every fiber in an RSS file.

**Usage**::

    lvm_gaussfit.py [-h] [-lmc] [-smc] [-v vel] [-stype SOURCE] [-out root] filename ...

**Options:**

-h
    Print help and exit.

-lmc
    Apply LMC velocity offset for fitting (~262 km/s).

-smc
    Apply SMC velocity offset for fitting (~146 km/s).

-v vel
    Apply a custom velocity offset in km/s.

-stype SOURCE|BACK
    For extracted text spectra, fit the source or background column.

-out root
    Set the root name for output files.

**Arguments:**

filename
    One or more SFrame-compatible FITS files, or text files containing
    extracted spectra with WAVE, FLUX, ERROR columns.

**Output:**

A table (FITS and ASCII) with one row per fiber containing:

- Fitted line fluxes (flux_ha, flux_nii, flux_sii_a, flux_sii_b, etc.)
- Line widths (fwhm_ha, fwhm_nii, etc.)
- Velocity offsets
- Fit quality metrics
- Fiber positions (ra, dec)

**Units:**

The LVM DRP stores flux density in units of erg/s/cm²/Å.  Before fitting,
all spectra are multiplied by 10¹⁶ to bring values into a numerically
convenient range.  Because the Gaussian model is parameterized by its
integrated flux (i.e. the analytic integral over wavelength), the output
``flux_*`` columns have units of **erg/s/cm² × 10¹⁶**.  To convert to
physical integrated line fluxes divide by 10¹⁶::

    flux_physical [erg/s/cm²]  =  flux_col  ×  1e-16

Wavelengths (line centers, FWHM) are in Ångströms.  The ``back_*`` columns
retain the scaled flux-density units (erg/s/cm²/Å × 10¹⁶).

**Example**::

    # Fit lines in an LMC observation
    lvm_gaussfit.py -lmc -out snr_n49 lvmSFrame-00012345.fits

    # Fit with custom velocity
    lvm_gaussfit.py -v 280 -out my_source data.fits

**Lines Fitted:**

The standard line list includes:

- H-alpha (6563 A)
- [NII] 6548, 6584 A
- [SII] 6717, 6731 A
- [OIII] 5007 A
- H-beta (4861 A)
- And others depending on wavelength coverage


sky_gaussfit.py — Fiber-by-Fiber Nebular and Airglow Fitting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fits single Gaussians to a fixed set of nebular emission lines and airglow
lines for every good science fiber in a sky-subtracted LVM SFrame file,
producing one output row per fiber.  This is complementary to
``lvm_gaussfit.py``: whereas ``lvm_gaussfit.py`` is a general-purpose
fitter for individual spectra or RSS files, ``sky_gaussfit.py`` is tuned
for survey-scale fiber-by-fiber analysis of SFrame data and includes a
comprehensive airglow line set for monitoring sky-subtraction residuals.

**Usage**::

    sky_gaussfit.py [-lmc] [-smc] [-v vel] [-out root] [-np nproc] filename [filename ...]

**Options:**

-lmc
    Apply the LMC radial velocity (~262 km/s) to nebular lines.

-smc
    Apply the SMC radial velocity (~146 km/s) to nebular lines.

-v vel
    Apply an arbitrary radial velocity in km/s to nebular lines.

-out root
    Root name for the output file (default: derived from input filename).

-np nproc
    Number of parallel processes for fiber fitting (default: 8).

**Arguments:**

filename
    One or more SFrame FITS files.  ASCII spectrum files with WAVE and FLUX
    columns are also accepted (processed via ``do_individual``).

**Lines fitted:**

Nebular lines — wavelengths are shifted by the supplied velocity:

=========  ===========  ==========
Line       Wavelength   Column tag
=========  ===========  ==========
[OIII]     4958.911 A   oiii_a
[OIII]     5006.843 A   oiii_b
[OI]       6300.309 A   oi_a
[OI]       6363.783 A   oi_b
[NII]      6548.04  A   nii_a
Ha         6562.80  A   ha
[NII]      6583.46  A   nii_b
[SII]      6716.440 A   sii_a
[SII]      6730.815 A   sii_b
=========  ===========  ==========

Airglow lines — fitted at fixed, unshifted wavelengths (ESO UVES atlas):

=========  ===========
Line name  Wavelength
=========  ===========
sky5577    5577.34 A
sky6300    6300.31 A
sky6363    6363.78 A
sky6533    6533.04 A
sky6553    6553.617 A
sky6577    6577.2  A
sky6912    6912.62 A
sky6923    6923.22 A
sky6939    6939.52 A
sky7358    7358.68 A
sky7392    7392.21 A
sky7914    7913.72 A
sky8344    8344.61 A
sky8399    8399.18 A
sky8827    8827.11 A
sky8988    8988.38 A
sky9552    9552.55 A
sky9719    9719.84 A
=========  ===========

Note that sky6300 and sky6363 overlap the nebular [OI] doublet.  Both are
fit independently: the sky lines at fixed wavelengths, the nebular lines
velocity-shifted.

**Output:**

When the input is an SFrame FITS file, the output is an ASCII fixed-width
table (one row per successfully fit fiber) with columns covering the fit
parameters (flux, wave, fwhm, back, rmse) for each line together with
fiberid, ra, and dec.  The output filename defaults to the input filename
with ``.fits`` replaced by ``.txt``, or ``<root>.txt`` if ``-out`` is
supplied.

When the input is one or more ASCII spectrum files, output is written to
``Gauss_<stem>.txt`` (single file) or ``Gauss_<root>.txt`` (multiple files).

**Performance:**

Fiber fitting is parallelized using ``multiprocessing.Pool``.  The spectrum
is pre-trimmed to the fitting range and per-line index arrays are
pre-computed once per file, so each worker operates only on the small
wavelength window needed for each line.  Use ``-np 1`` to disable
parallelism for debugging.

**Example**::

    # Fit all fibers in an LMC SFrame file, 12 parallel processes
    sky_gaussfit.py -lmc -np 12 lvmSFrame-00012345.fits

    # Fit with an arbitrary velocity
    sky_gaussfit.py -v 280 -out my_field lvmSFrame-00009088.fits

    # Fit individual ASCII spectra (no velocity shift)
    sky_gaussfit.py spectrum1.txt spectrum2.txt


Spatial maps of the per-fiber fit results (wavelength residuals, flux
residuals, and FWHM residuals for all 18 airglow lines) are produced by
``plot_sky_gaussfit.py``; see :doc:`data_quality` for full documentation.


lvm_line_profile.py — Gaussian vs Moffat Airglow Line Profiles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fits both a Gaussian and a Moffat profile (each with a constant
background) to the same 18 airglow lines used by ``sky_gaussfit.py``,
but on **raw, pre-subtraction** sky spectra rather than sky-subtracted
science fibers — the goal is to characterize the true instrumental line
shape directly, independent of any sky-subtraction algorithm or the
PALACE decomposer used by ``SkySubDev2.py`` (see :doc:`sky_subtraction`).

Fully standalone: it imports nothing from any other ``py_progs`` script
(only external packages — numpy, astropy, matplotlib, scipy, lmfit).
Its own copy of the airglow line list and its own IVAR estimator are
self-contained, deliberately not shared with ``lvm_gaussfit.py`` or
``sky_gaussfit.py``, so those scripts remain free to be reused/changed
for other analyses without affecting this one.

**Usage**::

    lvm_line_profile.py [-delta N] [-lines name1,name2,...]
                        [-ext SKY_EAST|SKY_WEST] [-out ROOT] filename

**Arguments:**

filename
    An XCframe FITS file (``SKY_EAST``/``SKY_WEST`` extensions) or a
    ``Sky_<name>.fits`` file (``FLUX`` extension) containing raw,
    pre-subtraction sky spectra.

**Options:**

-delta N
    Process every N-th row (default: 50).

-lines LIST
    Comma-separated subset of line names to fit (default: all 18).

-ext NAME
    Which raw sky column to use for an XCframe file: ``SKY_EAST``
    (default) or ``SKY_WEST``; ignored for ``Sky_<name>.fits`` files,
    which use ``FLUX`` directly.

-out ROOT
    Output filename root (default: ``<stem>_profile``).

**Description:**

For each requested line, a median spectrum is formed across all
selected rows and fit with both profiles — the primary, highest-S/N
comparison, plotted as data + both fits + residuals. Each individual
selected row is also fit with both profiles, building a per-row table
used to test whether any preference for one profile is consistent
across many independent exposures, not just a property of the deep
stack.

Both profiles are parameterized so their fitted flux equals the
analytic integral of the profile (matching ``lvm_gaussfit.py``'s
convention), and both fit the line **center freely** (bounded only to
the fit window, not fixed at the nominal catalog wavelength) — so a
real wavelength shift is not assumed away. The Moffat profile reduces
to a Gaussian as beta → ∞, so a finite, well-constrained beta across
many lines/rows would indicate genuine non-Gaussian wings; beta
drifting to the fit's upper bound indicates no real preference for
Moffat.

Because Moffat has one more free parameter (beta) than Gaussian, raw
chi², RMSE always favour it — comparisons use AIC/BIC instead, which
penalize that extra freedom::

    delta_aic = aic_gaussian - aic_moffat   (positive: Moffat preferred)
    delta_bic = bic_gaussian - bic_moffat
    rms_core_<profile>, rms_wing_<profile>  -- a fixed, profile-
        independent split of the fit window (core = center +/- 1.5 A,
        wing = the rest), so both profiles are judged against the same
        pixels regardless of their own fitted width

If the input file has an LSF extension (per-row, per-wavelength FWHM in
Ångströms — the same one ``SkySubDev2.py`` uses), the fitted Gaussian
FWHM is also compared directly against that file's own stated LSF FWHM
at each line's wavelength — an independent check of whether the LSF
extension's values are themselves accurate, entirely outside any
subtraction algorithm.

**Output:**

- ``<ROOT>_median_<line>.png`` — one file per line: median spectrum,
  both fits (with fitted FWHM/shift and the LSF extension's FWHM if
  available), and both fits' residuals.
- ``<ROOT>_perrow.txt`` — ascii table, one row per (input row, line):
  both fits' parameters (including fitted center/shift), AIC/BIC,
  core/wing RMS, and (if available) the LSF-extension comparison.
- ``<ROOT>_summary.png`` — per-line distributions of delta_AIC and beta
  across all fitted rows; also fit-FWHM-vs-LSF-extension difference and
  center shift, if an LSF extension is available.

A summary table (per line: median delta_AIC, fraction of rows favouring
Moffat by delta_AIC > 2, median beta and its IQR, and — if available —
median FWHM-vs-LSF-extension difference and median center shift) is
also printed to the screen.

**Example**::

    # Compare profile shapes for all 18 lines, sampling every 50th row
    lvm_line_profile.py -delta 50 XCframe_1.2.1_7325_48860_1_50.fits

    # Just the two brightest optical lines, every row
    lvm_line_profile.py -delta 1 -lines sky5577,sky6300 XCframe_file.fits

**Findings so far (260709):**

Run across the full 18-line set on real XCframe data (280 rows), the
Gaussian profile is preferred (or at least not improved upon) at every
single line — delta_AIC negative and beta pinned near its fit ceiling
everywhere. Separately, the fitted Gaussian FWHM comes out consistently
wider than the file's own LSF extension value at every line (a modest
but real few-percent-to-~12% gap, not a smooth function of wavelength);
that irregular, line-specific pattern is more consistent with several
of these catalog lines being unresolved blends of close airglow
transitions than with a genuine, correctable LSF calibration error —
still under investigation. This same run also directly caught and
confirmed a wavelength error in ``sky_gaussfit.py``'s own line list:
``sky6553`` was catalogued at 6553.0 Å but the data consistently show
its true centroid at 6553.617 Å (now corrected in both
``sky_gaussfit.py`` and this script's own line-list copy).


Multi-Component Fitting
-----------------------

lvm_double.py
^^^^^^^^^^^^^

Experimental routine for fitting a single emission line with either
a single or double Gaussian profile. This is useful for detecting
velocity structure (e.g., expanding shells, multiple velocity components).

**Usage**::

    lvm_double.py [-h] [-wmin wavelength] [-wmax wavelength] filename

**Options:**

-h
    Print help and exit.

-wmin, -wmax
    Wavelength range for fitting.

**Arguments:**

filename
    FITS or text file containing the spectrum.

**Output:**

Fit parameters for both single and double Gaussian models, along with
statistical comparison to determine if the double Gaussian is justified.

lvm_triple.py
^^^^^^^^^^^^^

Similar to lvm_double.py but allows fitting up to three Gaussian
components. Useful for complex line profiles.

**Usage**::

    lvm_triple.py [-h] [-wmin wavelength] [-wmax wavelength] filename

**Options:**

Same as lvm_double.py.


Flux Calculation
----------------

lvm_flux.py
^^^^^^^^^^^

Calculates emission line fluxes using velocity and width parameters
determined from another line (typically H-alpha). This is useful when
fitting weak lines that cannot be fit independently.

**Usage**::

    lvm_flux.py filename

**Description:**

Given fit parameters from a strong line (e.g., velocity and FWHM from
H-alpha), this routine calculates fluxes for other lines by fixing the
velocity and width to the known values and only fitting the amplitude.


Working with Extracted Spectra
------------------------------

The fitting tools can work with either:

1. **RSS FITS files** - Standard SFrame or combined RSS files, where
   fitting is performed for each fiber automatically.

2. **Text files** - Extracted spectra in ASCII format with columns for
   WAVE, FLUX, and ERROR. Additional columns (SOURCE_FLUX, BACK_FLUX)
   can be used with the ``-stype`` option.

For extracted spectra::

    # Fit the flux column (background-subtracted if available)
    lvm_gaussfit.py spectrum.txt

    # Fit the original source spectrum (before background subtraction)
    lvm_gaussfit.py -stype SOURCE spectrum.txt

    # Fit the background spectrum
    lvm_gaussfit.py -stype BACK spectrum.txt


Typical Workflows
-----------------

Basic Emission Line Mapping
^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Obtain an SFrame or combined RSS file
2. Run lvm_gaussfit to fit all fibers::

       lvm_gaussfit.py -lmc -out my_source lvmSFrame-00012345.fits

3. Use the output table for analysis or mapping

Velocity Structure Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Extract a spectrum from a region of interest
2. Use lvm_double or lvm_triple to check for multiple components::

       lvm_double.py -wmin 6550 -wmax 6580 extracted_spectrum.txt

3. Interpret the fit statistics to determine if multiple components
   are present

Fitting Weak Lines
^^^^^^^^^^^^^^^^^^

1. First fit strong lines to get velocity and width::

       lvm_gaussfit.py -lmc -out initial_fit data.fits

2. Use lvm_flux to measure weak lines with fixed parameters::

       lvm_flux.py initial_fit.fits


Notes
-----

- Velocity offsets are important for correct line identification,
  especially for LMC/SMC targets
- The fitting uses lmfit with Levenberg-Marquardt optimization
- Fits include error estimates from the covariance matrix
- Failed fits are flagged in the output table
- For RSS files, the output table includes fiber positions for mapping
- All ``flux_*`` output columns are integrated line fluxes in units of
  erg/s/cm² × 10¹⁶.  Divide by 10¹⁶ to obtain physical values in erg/s/cm².


See Also
--------

- :doc:`dap` - Running LVM-DAP for stellar-population + emission-line
  fits, an alternative to the ``lmfit``-based tools on this page
- :doc:`snapshots` - Batch processing with automatic fitting
- :doc:`summarize` - Summarizing exposures; ``gauss_offset.py`` for airglow monitoring
- :doc:`data_quality` - ``plot_sky_gaussfit.py`` for spatial maps of sky Gaussian fit residuals
- :doc:`api/lvm_gaussfit/index` - API documentation
- :doc:`api/sky_gaussfit/index` - API documentation
- :doc:`api/lvm_line_profile/index` - API documentation
- :doc:`api/lvm_double/index` - API documentation
- :doc:`api/lvm_triple/index` - API documentation
- :doc:`api/lvm_flux/index` - API documentation
