Combining RSS Spectra from Multiple Exposures of Extended Regions
=================================================================

The LVM Data Reduction Pipeline produces row-stacked spectra (RSS) in
SFrame files, where each row corresponds to a single fiber. When a region
has been observed multiple times, you may want to combine these exposures
into a single, deeper dataset. Three scripts are provided for this,
covering two different combination strategies (see `Two Combination
Strategies`_ below):

- ``rss_combine.py`` - Combine exposures covering an extended region
- ``rss_combine_pos.py`` - Extract and combine a region centered on a specific position
- ``gauss_combine.py`` - Combine per-exposure emission-line fits onto the same kind of output grid

The first two produce RSS FITS files that can be used for further analysis,
including spectral fitting with the DAP or custom analysis tools. Both of
these scripts can be used combine exposures of a single dithered tile, or
alternatively of an extended region that extends over multiple tiles.

Once you have a combined RSS file, a common next step is to collapse it
into a 2D FITS image for a quick look -- see `Viewing the Combined Output`_
below, which uses ``rss2image.py``.


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


Two Combination Strategies
---------------------------

There are two different ways to combine multiple exposures of the same
region into a single, deeper measurement, and the right one depends on
what you are combining:

**Combine-then-fit** (``rss_combine.py``, ``rss_combine_pos.py``)
    Raw spectra are apportioned onto the output fiber grid and combined
    first; any line fitting (e.g. with the DAP or ``lvm_gaussfit.py``)
    happens afterward, on the already-combined spectrum. This is the
    approach documented in the rest of this page up to `gauss_combine.py
    - Fit-Then-Combine`_ below.

**Fit-then-combine** (``gauss_combine.py``)
    Each exposure is fit independently first (e.g. with
    ``lvm_gaussfit.py`` or ``lvm_snrfit.py``), so every per-fiber line
    measurement already carries its own uncertainty. Those per-exposure
    *measurements* -- not the raw spectra -- are then apportioned onto
    the output grid, weighted by both spatial fiber-overlap and
    measurement uncertainty, so a noisy or unreliable exposure is
    naturally down-weighted rather than averaged in blind. See
    `gauss_combine.py - Fit-Then-Combine`_ below.

Combine-then-fit is the more general approach: any line list can be fit
on the combined spectrum after the fact. Fit-then-combine gives you a
principled per-contributor noise weight for free, at the cost of
committing to whatever line list the per-exposure fit already measured.


rss_combine.py - General RSS Combination
----------------------------------------

Use ``rss_combine.py`` when you want to combine all exposures of a region
into a single RSS file. The script automatically determines the spatial
extent from the input files.

**Command line usage**::

    rss_combine.py [-orig] [-sum] [-med] [-keep] [-helio] [-outroot name] filenames

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

-helio
    **Experimental.** Apply a heliocentric whole-pixel wavelength shift
    to each input file before combining, based on its ``WAVE HELIORV_SCI``
    header value (see :doc:`api/rss_combine/index` for
    ``compute_helio_shifts``). Off by default.

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

    rss_combine_pos.py [-sum] [-med] [-size arcmin] [-out name] [-keep] [-helio] ra dec filenames

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

-helio
    **Experimental.** Apply a heliocentric whole-pixel wavelength shift
    to each input file before combining, based on its ``WAVE HELIORV_SCI``
    header value (same mechanism as ``rss_combine.py``'s ``-helio``). Off
    by default.

**Examples**::

    # Extract 10 arcmin region around SNR N49
    rss_combine_pos.py -size 10 -out N49 81.501 -66.082 lvmSFrame-*.fits

    # Extract with median combination
    rss_combine_pos.py -med -size 15 -out N49B 81.342 -65.996 lvmSFrame-*.fits


Viewing the Combined Output
----------------------------

The FITS files produced by ``rss_combine.py`` and ``rss_combine_pos.py``
are row-stacked spectra, not images, so they cannot be displayed directly
in a tool like DS9. ``rss2image.py`` collapses the FLUX extension (or any
other extension of the same shape) over a wavelength range and projects
the fiber positions onto a regular pixel grid, producing a standard 2D
FITS image with WCS.

**Command line usage**::

    rss2image.py [-no_back] [-band filter] [-ext IVAR] filename(s)

**Example**::

    # Combine a tile, then make a quick H-alpha image of the result
    rss_combine.py -orig -med -outroot tile12345 lvmSFrame-*12345*.fits
    rss2image.py -band ha tile12345.med.fits

    # View the image (output is named FLUX_<input_root>_<band>.fits)
    ds9 FLUX_tile12345.med_ha.fits

This is often the fastest way to sanity-check a combination before moving
on to spectral fitting. See :doc:`visualization` for the full set of
options (predefined and custom wavelength bands, extension selection,
background subtraction).


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


gauss_combine.py - Fit-Then-Combine
------------------------------------

Use ``gauss_combine.py`` when you already have per-exposure emission-line
fits (e.g. from ``lvm_gaussfit.py`` or ``lvm_snrfit.py``) and want a
single table of line fluxes covering the full region, rather than a
combined spectrum you still need to fit. Because each exposure's fit
already carries its own uncertainty, the combine step can down-weight a
noisy or unreliable exposure instead of averaging it in blind -- see
`Two Combination Strategies`_ above.

**Command line usage**::

    gauss_combine.py [-h] [-outroot xxxx] [-exclude f1,f2,...] filenames

**Arguments:**

filenames
    One or more per-exposure line-fit tables (``ascii.fixed_width_two_line``,
    one row per science fiber), such as ``lvm_gaussfit.py``'s
    ``*.gauss.txt`` output. Wildcards are expanded.

**Options:**

-outroot xxxx
    Root name for the output table (default ``gauss_combine``); output
    is written to ``<outroot>.txt``.

-exclude f1,f2,...
    Comma-separated basenames of input files to drop before combining
    (e.g. an exposure already confirmed bad by ``OverlapFlux.py``).

**How it works**

The output fiber grid and each input fiber's fractional overlap with it
are built exactly as ``rss_combine.py`` builds them for raw spectra
(``create_wcs``/``generate_grid``/``frac_calc2``, imported unmodified) --
the difference is what gets apportioned. Line keys (e.g. ``ha``,
``oii_a``) are discovered automatically from the input tables'
``flux_<key>`` columns, so this works against any fit table following
that naming convention.

For each line and each output fiber, every contributing (exposure, input
fiber) pair is weighted by ``frac / eflux**2`` -- the spatial overlap
fraction times inverse-variance -- and combined without any significance
or sign cut. This mirrors ``rss_combine.py``'s own philosophy of never
filtering individual contributors before combining: a significance cut
(``flux/eflux >= minsnr``) would bias the combined result high for lines
near a single exposure's detection threshold, by keeping only the upward
noise fluctuations (Eddington bias) -- confirmed on real Vela data (Hb/
[OII]/[OIII] read 6-13% high at single-contributor output fibers when
such a cut was tried; see the script's History).

**Output columns**

The output table (``<outroot>.txt``, same ``ascii.fixed_width_two_line``
format as the input) has one row per output fiber (``fiberid``, ``X``,
``Y``, ``ra``, ``dec``), plus per discovered line key:

============  ================================================================
Column        Description
============  ================================================================
flux_<key>    Combined line flux (inverse-variance weighted)
eflux_<key>   Combined flux uncertainty
n_<key>       Number of surviving contributors
spread_<key>  Unweighted std of contributing fluxes (n>1 only)
wave_<key>    Combined line-center wavelength, when present in the input
fwhm_<key>    Combined FWHM, when present in the input
back_<key>    Combined local background, when present in the input
============  ================================================================

``n_<key>`` and ``spread_<key>`` let you see, after the fact, whether a
bin was built from disagreeing exposures -- useful since the script does
not itself correct exposure-to-exposure systematic offsets (as opposed to
per-exposure noise), only down-weight noisy ones. If a specific exposure
is confirmed bad (e.g. via ``OverlapFlux.py``), drop it with ``-exclude``
rather than relying on the weighting to fix it.

**Example**::

    # Combine per-exposure lvm_gaussfit.py fits into one line-flux table
    gauss_combine.py -outroot vela_lines lvmSFrame-*.gauss.txt

    # Same, but drop one exposure already confirmed bad
    gauss_combine.py -outroot vela_lines -exclude lvmSFrame-00014771.gauss.txt lvmSFrame-*.gauss.txt


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
- :doc:`visualization` - rss2image.py and other imaging/mapping tools
- :doc:`spectral_fitting` - lvm_gaussfit.py/lvm_snrfit.py, whose per-exposure output feeds gauss_combine.py
- :doc:`api/rss_combine/index` - API documentation for rss_combine
- :doc:`api/rss_combine_pos/index` - API documentation for rss_combine_pos
- :doc:`api/rss2image/index` - API documentation for rss2image
- :doc:`api/gauss_combine/index` - API documentation for gauss_combine
