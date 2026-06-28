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

A FITS file (``palace_mask_<stem>.fits``) containing:

- ``WAVE`` — wavelength array (Å)
- ``SKY`` — median observed sky spectrum (FACTOR-scaled)
- ``CONTINUUM`` — sky spectrum with contaminated pixels set to NaN
- ``MASK`` — boolean mask (1 = clean, 0 = contaminated)

A PNG diagnostic plot (``palace_mask_<stem>.png``) showing all three
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
