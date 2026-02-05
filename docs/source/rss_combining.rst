Combining RSS Spectra from Multiple Exposures of Extended Regions
=================================================================

The LVM Data Reduction Pipeline produces row-stacked spectra (RSS) in
SFrame files, where each row corresponds to a single fiber. When a region
has been observed multiple times, you may want to combine these exposures
into a single, deeper RSS file. Two scripts are provided for this purpose:

- ``rss_combine.py`` - Combine exposures covering an extended region
- ``rss_combine_pos.py`` - Extract and combine a region centered on a specific position

Both scripts produce RSS FITS files that can be used for further analysis,
including spectral fitting with the DAP or custom analysis tools. Both of
these scripts can be used combine exposures of a single dithered tile, or
alternatively of an extended region that extends over multiple tiles.


Why Combine RSS Files?
----------------------

LVM observations are taken as individual exposures, each producing an
SFrame file. Combining multiple exposures provides several benefits:

- **Increased signal-to-noise**: Stacking exposures improves S/N
- **Outlier rejection**: Using median combination removes outliers
- **Fill coverage gaps**: Dithered observations fill gaps between fibers
- **Simplified analysis**: Work with a single file instead of many

The combination process creates a new set of "virtual fibers" on a regular
grid, and apportions flux from the input fibers to these output positions
based on geometric overlap.


rss_combine.py - General RSS Combination
----------------------------------------

Use ``rss_combine.py`` when you want to combine all exposures of a region
into a single RSS file. The script automatically determines the spatial
extent from the input files.

**Command line usage**::

    rss_combine.py [-orig] [-sum] [-med] [-keep] [-outroot name] filenames

**Arguments:**

filenames
    One or more SFrame FITS files to combine.

**Options:**

-outroot name
    Root name for output files (default: creates name from input).

-orig
    Use average input fiber positions instead of a regular grid.
    Useful for single-tile observations where you want to preserve
    the original fiber geometry.

-sum
    Assign each input fiber's full flux to the nearest output fiber
    (no fractional splitting). Best for non-dithered observations.

-med
    Use median instead of mean for combining exposures. More robust
    to outliers and cosmic rays.

-keep
    Retain temporary files in xtmp/ directory for debugging.

**Examples**::

    # Combine all exposures of a tile using defaults (regular grid, mean)
    rss_combine.py lvmSFrame-*.fits

    # Combine preserving original fiber positions
    rss_combine.py -orig -outroot tile1234 lvmSFrame-*.fits

    # Combine using median (robust to outliers)
    rss_combine.py -med -outroot deep_field lvmSFrame-*.fits


rss_combine_pos.py - Position-Centered Extraction
-------------------------------------------------

Use ``rss_combine_pos.py`` when you want to extract a specific region
centered on a known position (e.g., a supernova remnant, planetary nebula,
or other discrete source).

**Command line usage**::

    rss_combine_pos.py [-sum] [-med] [-size arcmin] [-out name] [-keep] ra dec filenames

**Arguments:**

ra dec
    Center position in degrees (required).

filenames
    One or more SFrame FITS files to combine.

**Options:**

-size arcmin
    Size of output region in arcminutes (default: 20).

-out name
    Output filename root (default: 'test').

-sum
    Assign full flux to nearest output fiber (no fractional splitting).

-med
    Use median instead of mean for combining.

-keep
    Retain temporary files in xtmp/ directory.

**Examples**::

    # Extract 10 arcmin region around SNR N49
    rss_combine_pos.py -size 10 -out N49 81.501 -66.082 lvmSFrame-*.fits

    # Extract with median combination
    rss_combine_pos.py -med -size 15 -out N49B 81.342 -65.996 lvmSFrame-*.fits


How Flux Apportionment Works
----------------------------

When combining exposures, flux must be redistributed from input fibers
to output fiber positions. Two methods are available:

**Fractional splitting (default)**

Each input fiber's flux is split among nearby output fibers based on
the geometric overlap area between circular apertures (35 arcsec diameter).
This preserves flux and handles sub-fiber offsets correctly.

**Nearest neighbor (-sum)**

Each input fiber's full flux goes to the single nearest output fiber.
This is faster and appropriate when input fibers are already well-aligned
(non-dithered observations).


Combination Methods
-------------------

After remapping each input exposure to the output fiber grid, the
exposures are combined:

**Mean (default)**

Simple average of all valid exposures at each fiber/wavelength.
Maximizes S/N but sensitive to outliers.

**Median (-med)**

Median of all valid exposures. More robust to cosmic rays, bad pixels,
and other outliers. May slightly reduce S/N compared to mean.


Output File Format
------------------

Both scripts produce a FITS file with the following extensions:

============  ================================================================
Extension     Description
============  ================================================================
PRIMARY       Primary header from one of the input images
FLUX          Combined flux array (n_fibers x n_wavelengths)
IVAR          Inverse variance array
MASK          Quality mask (0=valid, 1=no data)
WAVE          Wavelength array (same as input)
LSF           Line spread function array
SLITMAP       Table of fiber positions (fiberid, RA, Dec, X, Y)
WCS_INFO      WCS defining the fiber coordinate system
EXPOSURE      Effective exposure count per fiber/wavelength
============  ================================================================

An ASCII table (``.tab``) with fiber positions is also written.

The output filename includes the combination method as a suffix:

- ``name.ave.fits`` - Mean combination
- ``name.med.fits`` - Median combination


Typical Workflows
-----------------

**Combining a single tile (multiple exposures)**::

    # Get all exposures of tile 12345
    ls /path/to/data/lvmSFrame-*12345*.fits > tile_files.txt

    # Combine using original fiber positions (best for single tile)
    rss_combine.py -orig -med -outroot tile12345 $(cat tile_files.txt)

**Creating a mosaic from multiple tiles**::

    # Combine all tiles in a region onto a regular grid
    rss_combine.py -med -outroot lmc_mosaic lvmSFrame-*.fits

**Extracting a source for detailed analysis**::

    # Extract SNR with 10 arcmin field of view
    rss_combine_pos.py -size 10 -med -out snr_n49 81.501 -66.082 lvmSFrame-*.fits


Notes
-----

- Input files should be flux-calibrated SFrame files from the DRP
- The scripts do not weight by exposure time or variance; all inputs
  are treated equally
- Large numbers of input files may require significant memory
- The ``-keep`` option is useful for debugging apportionment issues


See Also
--------

- :doc:`snapshots` - Batch processing of source catalogs with fitting
- :doc:`api/rss_combine/index` - API documentation for rss_combine
- :doc:`api/rss_combine_pos/index` - API documentation for rss_combine_pos
