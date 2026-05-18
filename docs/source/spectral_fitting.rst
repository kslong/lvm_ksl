Spectral Fitting
================

The lvm_ksl package provides tools for fitting emission lines in LVM
spectra. These tools can fit Gaussian profiles to standard emission
lines across all fibers in an RSS file, or perform more detailed fits
to individual spectra.


Overview
--------

The spectral fitting tools include:

- ``lvm_gaussfit.py`` - Fit standard emission lines across an RSS file
- ``sky_gaussfit.py`` - Fit nebular and airglow lines fiber-by-fiber in SFrame files
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
sky6553    6553.0  A
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

- :doc:`snapshots` - Batch processing with automatic fitting
- :doc:`summarize` - Summarizing exposures; ``gauss_offset.py`` for airglow monitoring
- :doc:`data_quality` - ``plot_sky_gaussfit.py`` for spatial maps of sky Gaussian fit residuals
- :doc:`api/lvm_gaussfit/index` - API documentation
- :doc:`api/sky_gaussfit/index` - API documentation
- :doc:`api/lvm_double/index` - API documentation
- :doc:`api/lvm_triple/index` - API documentation
- :doc:`api/lvm_flux/index` - API documentation
