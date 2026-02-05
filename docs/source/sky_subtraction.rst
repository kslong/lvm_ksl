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
