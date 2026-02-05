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


See Also
--------

- :doc:`snapshots` - Batch processing with automatic fitting
- :doc:`api/lvm_gaussfit/index` - API documentation
- :doc:`api/lvm_double/index` - API documentation
- :doc:`api/lvm_triple/index` - API documentation
- :doc:`api/lvm_flux/index` - API documentation
