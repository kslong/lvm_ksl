Visualization and Mapping
=========================

The lvm_ksl package provides tools for creating images and maps from
LVM spectral data. These tools convert RSS (row-stacked spectra) files
into 2D images by interpolating fiber data onto a regular grid.


Overview
--------

LVM data is collected as spectra from individual fibers arranged in a
hexagonal pattern. To create images, the fiber positions must be mapped
to a 2D grid. The visualization tools handle this conversion and provide
various options for creating emission line maps and broadband images.

Available tools:

- ``kslmap.py`` - Create images with predefined or custom filters
- ``quick_map.py`` - Quick image creation with flexible band options
- ``line_map.py`` - Create emission line maps
- ``rss2image.py`` - Convert RSS files to FITS image files
- ``rss_explore.py`` - Interactive image and spectrum explorer
- ``MakeVideo.py`` - Create videos from image sequences


Creating Images from RSS Files
------------------------------

kslmap.py
^^^^^^^^^

Creates images from LVM exposures using predefined emission line
filters or custom wavelength ranges.

**Usage**::

    kslmap.py [-no_back] [-image_type filter] filename

**Options:**

-no_back
    Do not subtract background from the image.

-image_type filter
    Predefined filter to use. Options include:

    - ``ha`` - H-alpha emission
    - ``sii`` - [SII] emission

**Arguments:**

filename
    SFrame or CFrame FITS file.

**Output:**

A FITS image file with proper WCS (World Coordinate System) that can
be displayed in tools like DS9.

**Example**::

    # Create an H-alpha image
    kslmap.py -image_type ha lvmSFrame-00012345.fits

    # Create image without background subtraction
    kslmap.py -no_back -image_type ha lvmSFrame-00012345.fits

quick_map.py
^^^^^^^^^^^^

A more flexible version of kslmap.py with additional band options
and custom wavelength ranges.

**Usage**::

    quick_map.py [-no_back] [-band filter] filename

**Options:**

-no_back
    Do not subtract background from the image.

-band filter
    Predefined filter or custom range. Options:

    - ``ha`` - H-alpha emission
    - ``sii`` - [SII] emission
    - ``x wmin wmax`` - Custom wavelength range (e.g., ``-band x 6700 6800``)

**Arguments:**

filename
    SFrame or CFrame FITS file.

**Output:**

- For SFrame files: output file starts with ``x``
- For CFrame files: output file starts with ``z``

**Examples**::

    # Create H-alpha map
    quick_map.py -band ha lvmSFrame-00012345.fits

    # Create custom band image (6700-6800 A)
    quick_map.py -band x 6700 6800 lvmSFrame-00012345.fits

line_map.py
^^^^^^^^^^^

Creates emission line maps, similar to kslmap.py but with different
options for line selection.

**Usage**::

    line_map.py [-no_back] [-image_type filter] filename

**Options:**

Same as kslmap.py.

**Notes:**

Uses the same predefined filters as kslmap.py. The image orientation
is set for standard DS9 display.


Extracting Spectra
------------------

GetSpec.py
^^^^^^^^^^

Extracts spectra from fibers at a specified position, with options
for averaging over a region.

**Usage**::

    GetSpec.py [-h] [-root name] [-median] [-sum] filename ra dec [rmin] [rmax] [-back r1 r2]

**Options:**

-h
    Print help and exit.

-root name
    Prepend a root name to the output file.

-median
    Use median instead of average when combining fibers.

-sum
    Use sum instead of average when combining fibers.

**Arguments:**

filename
    RSS FITS file.

ra, dec
    Position in degrees or sexagesimal format (h:m:s d:m:s).

rmin
    Inner radius in arcseconds (optional).

rmax
    Outer radius in arcseconds (optional).

-back r1 r2
    Subtract background from annulus between r1 and r2 arcseconds.

**Extraction modes:**

- No radius: Extract single closest fiber
- rmin only: Extract all fibers within rmin
- rmin and rmax: Extract fibers in annulus between rmin and rmax

**Output:**

- ASCII table with WAVE, FLUX, ERROR columns
- DS9 region file showing which fibers were used

**Examples**::

    # Extract single fiber at position
    GetSpec.py lvmSFrame-00012345.fits 81.5 -66.0

    # Extract and average fibers within 30 arcsec
    GetSpec.py lvmSFrame-00012345.fits 81.5 -66.0 30

    # Extract with background subtraction
    GetSpec.py lvmSFrame-00012345.fits 81.5 -66.0 30 -back 60 90


Interactive Exploration
-----------------------

rss_explore.py
^^^^^^^^^^^^^^

Interactive RSS image and spectrum explorer.  Displays a spatial flux
image from an RSS FITS file with two adjustable sliders (center
wavelength and passband width) and a click-to-spectrum panel below.
Can be used from a Jupyter notebook or run directly from the command line.

**Command-line usage**::

    rss_explore.py filename [--center WAVELENGTH] [--width WIDTH] [--ext EXT]

**Options:**

--center WAVELENGTH
    Starting center wavelength in Angstroms (default: 6563, H-alpha).

--width WIDTH
    Starting passband width in Angstroms (default: 20).

--ext EXT
    FITS extension to display (default: FLUX).

**Command-line examples**::

    # Open at H-alpha
    rss_explore.py lvmSFrame-00007373.fits

    # Open at [O III] with a 15 A window
    rss_explore.py lvmSFrame-00007373.fits --center 5007 --width 15

    # Display the IVAR extension
    rss_explore.py lvmSFrame-00007373.fits --ext IVAR

**Jupyter notebook usage**::

    %matplotlib widget          # must be the first line, before imports

    from rss_explore import RSSExplorer

    exp = RSSExplorer('lvmSFrame-00007373.fits')
    exp.show()                       # start at H-alpha
    exp.show(center=5007, width=15)  # start at [O III]

Multiple calls to ``show()`` can be active simultaneously in different
cells; each cell maintains its own independent sliders, image, and
spectrum panel.

**Interactions:**

+------------------------------+-------------------------------------------+
| Action                       | Effect                                    |
+==============================+===========================================+
| Move Center or Width slider  | Redraws the spatial flux image for the    |
|                              | new wavelength passband                   |
+------------------------------+-------------------------------------------+
| Click on the image           | Shows the spectrum of the nearest fibre   |
|                              | in the lower panel; marks the fibre with  |
|                              | a red cross; displays fibre number,       |
|                              | RA and Dec in the panel title             |
+------------------------------+-------------------------------------------+
| Drag a box on the spectrum   | Zooms into the selected wavelength and    |
|                              | flux region.  Switching fibres while      |
|                              | zoomed keeps the same axes so fibres can  |
|                              | be compared on the same scale             |
+------------------------------+-------------------------------------------+
| Z button / Reset zoom button | Restores the full spectrum view           |
+------------------------------+-------------------------------------------+

**Requirements:**

The notebook version requires ``ipympl`` for the interactive canvas::

    pip install ipywidgets ipympl

**Notes:**

The RA/Dec tick labels on the image are approximate: RA values are
evaluated at the middle row of the image and Dec values at the middle
column.  The error is negligible over a typical LVM field (~0.5 deg).

On first load the tool precomputes pixel masks for every fibre
(a few seconds).  All subsequent slider and click interactions are
fast (< 1 s).


Additional Tools
----------------

rss2image.py
^^^^^^^^^^^^

Converts RSS files to FITS image files by summing flux over a
specified wavelength range and projecting fibre positions onto a
regular pixel grid.

MakeVideo.py
^^^^^^^^^^^^

Creates video sequences from a series of images, useful for visualizing
time-series data or spectral scans.


Coordinate Conversion
---------------------

fib2radec.py
^^^^^^^^^^^^

Converts fiber positions to RA/Dec coordinates. This is used internally
by the mapping tools but can also be run standalone.

**Usage**::

    fib2radec.py filename

**Description:**

Reads an LVM FITS file and computes RA/Dec for each fiber based on
the telescope pointing and fiber positions in the IFU.


Typical Workflows
-----------------

Quick Look at an Exposure
^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # Create an H-alpha image for quick inspection
    quick_map.py -band ha lvmSFrame-00012345.fits

    # View in DS9
    ds9 xlvmSFrame-00012345_ha.fits

Creating Multi-Band Images
^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # Create images in multiple bands
    quick_map.py -band ha lvmSFrame-00012345.fits
    quick_map.py -band sii lvmSFrame-00012345.fits
    quick_map.py -band x 5000 5020 lvmSFrame-00012345.fits  # [OIII]

Extracting and Plotting a Spectrum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # Extract spectrum at a position
    GetSpec.py lvmSFrame-00012345.fits 81.5 -66.0 20

    # The output can be plotted with standard tools or
    # used as input to lvm_gaussfit.py


Notes
-----

- Images are created with proper WCS for overlay with other data
- The default orientation matches standard astronomical conventions
- Velocity corrections for LMC/SMC are applied automatically based
  on the pointing position
- Interpolation uses scipy's griddata with linear interpolation
- Background subtraction uses regions away from the emission lines


See Also
--------

- :doc:`rss_combining` - Combining multiple exposures before mapping
- :doc:`spectral_fitting` - Fitting emission lines in extracted spectra
- :doc:`api/kslmap/index` - API documentation
- :doc:`api/quick_map/index` - API documentation
- :doc:`api/GetSpec/index` - API documentation
- :doc:`api/rss_explore/index` - API documentation for rss_explore
