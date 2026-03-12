Generating Region Files
=======================

The ``GenAnnularBackground.py`` script automates the creation of paired
source and background annulus region entries from a master source catalog.
It is primarily intended for use with catalogs of extended objects such as
supernova remnants in the LMC or SMC, where each source requires a local
background region for spectral analysis.


Overview
--------

When extracting spectra of extended sources it is important to sample a
local background that is free of source emission.  This script takes a
master table of source positions and angular sizes and, for each source,
constructs a background annulus that:

- surrounds the source with a configurable gap (``xspace``)
- has an area comparable to the source itself
- is marked distinctly in the output table (green, ``SourceBack='Back'``)

The resulting table interleaves source and background rows and can be
used directly as input to region-based spectral extraction tools.


Input File Format
-----------------

The input is an ASCII table readable by ``astropy.io.ascii.read``.

Required columns:

=========  =====================================================================
Column     Description
=========  =====================================================================
Major      Semi-major axis of the source in arcseconds
Minor      Semi-minor axis of the source in arcseconds
Color      Display color for the source region (string)
=========  =====================================================================

Additional columns (e.g. source name, RA, Dec) are preserved unchanged
in the output.


Processing Steps
----------------

For each source row, the script:

1. Enforces a minimum size: ``Major`` and ``Minor`` are replaced by
   ``xmin`` (default 35 arcsec) if they are smaller.

2. Creates a background annulus row with:

   - Inner radius = ``Major`` + ``xspace``
   - Outer radius = 2 × ``Major`` + ``xspace``
   - ``RegType`` = ``'annulus'``
   - ``SourceBack`` = ``'Back'``
   - ``Color`` = ``'green'``

3. Stacks and sorts source and background rows by source number and
   ascending ``Major``, so each source is immediately followed by its
   background annulus.


Output File
-----------

The output is written as an ASCII fixed-width table named
``<masterfile_root>.ann_reg.txt``, e.g. ``smc_snr_cotton24.ann_reg.txt``.

New columns added to the output table:

============  ==============================================================
Column        Description
============  ==============================================================
No.           Sequential source number (1-based)
SourceBack    ``'Source'`` for source rows, ``'Back'`` for background rows
RegType       Set to ``'annulus'`` for background rows
============  ==============================================================


Command Line Usage
------------------

Basic syntax::

    GenAnnularBackground.py [-h] masterfile

**Arguments:**

masterfile
    Input ASCII table containing source sizes and positions.

**Options:**

-h
    Print help message and exit.


Examples
--------

Generate annular background regions for an SMC SNR catalog::

    GenAnnularBackground.py smc_snr_cotton24.txt

This produces ``smc_snr_cotton24.ann_reg.txt`` with paired source and
background rows for every entry in the catalog.


Notes
-----

- The ``xmin`` (minimum source size, default 35 arcsec) and ``xspace``
  (gap between source and annulus, default 35 arcsec) parameters are
  currently only adjustable by calling ``doit()`` directly from Python,
  not from the command line.

- The annulus width equals ``Major`` (after enforcing ``xmin``), so the
  background aperture area is comparable to the source aperture area.

- The script does not modify positional columns (RA, Dec, etc.).


See Also
--------

- :doc:`api/GenAnnularBackground/index` - API documentation
- :doc:`snapshots` - Creating spectral snapshots of individual sources
