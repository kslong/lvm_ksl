Sky from the Science Field
==========================

This page is part of :doc:`sky_subtraction`.

Unlike the methods in :doc:`sky_methods` (which use the dedicated SKY_EAST/SKY_WEST
sky telescopes), the scripts in this section take the sky from the science
IFU itself.  The sky then comes from the same direction, at the same time
and through the same telescope as the source, so it avoids the sky
telescopes' different pointing (and, e.g., one placed too close to the
Moon).  The price is that whatever source emission the chosen fibers
contain is subtracted too: in an extended nebula the result is the
emission in *excess* of the faintest part of the field, not the absolute
emission.

There are two approaches:

- ``SkySubSci.py`` and ``SummarizeSciSky.py`` rank fibers by sky-line-free
  continuum flux and average the faintest ones into one sky spectrum per
  exposure.  No scale factor is applied to emission lines, so this suits
  fields where a genuinely sky-dominated tail of faint fibers exists (e.g.
  diffuse or extended sources, not compact point sources filling the IFU).
  Both write WAVE/SCI/SKY/FLUX/DRP_ALL extensions compatible with
  ``SkySub_eval.py`` (FLUX = sky-subtracted, SKY = sky model), so they can
  be evaluated alongside the methods in :doc:`sky_methods` (see
  :doc:`sky_method_eval`).
- ``SkySubPatch.py`` builds the sky from a compact *patch* of fibers with
  the least nebular emission, then transforms it separately for every
  fiber to that fiber's own wavelength offset, line-spread function and
  throughput before subtracting it, and writes a full lvmSFrame-layout
  file.  ``lsf_kernel.py`` provides the line-spread-function matching.


SkySubSci.py
------------

Estimate a sky spectrum for one or more LVM CFrame exposures directly from
the science fibers, without using the dedicated sky telescopes.

**Usage**::

    SkySubSci.py [-low PCT] [-high PCT] [-navg N] [-sigma S] [-maxiters K]
                [-mask FILE] [-stat median|mean] [-out ROOT]
                filename [filename ...]

**Arguments:**

filename
    One or more lvmCFrame FITS files.  One output row is written per file,
    in the order given.

**Options:**

-low PCT
    Percentile rank (0-100) of the faint/sky-like fiber (default 10).

-high PCT
    Percentile rank (0-100) of the bright/science-like fiber (default 90).

-navg N
    Number of fibers, ranked closest to ``-low``/``-high``, combined with a
    sigma-clipped (robust) mean (default 10; use 1 to reproduce the
    original single-fiber behaviour).

-sigma S
    Sigma-clipping threshold for the robust mean (default 3.0).

-maxiters K
    Sigma-clipping iteration limit (default 5).

-mask FILE
    palace_mask FITS file from ``palace_make_mask.py`` (WAVE/MASK
    extensions, MASK=1 means clean/sky-line-free).  If omitted,
    ``sky_mask.fits`` is searched for in the current directory, then in
    the ``lvm_ksl data/`` directory.

-stat STAT
    Statistic used to rank fibers by continuum flux: ``median`` (default)
    or ``mean``.

-out ROOT
    Output filename root.  Default:
    ``SkySubSci_<first_expnum>_<last_expnum>`` (or ``SkySubSci_<expnum>``
    for a single file).

**Description:**

For each exposure, science fibers are selected from SLITMAP (``scifib``,
imported from ``SummarizeCframe.py``).  Each fiber's continuum flux is
measured as the median (or mean) FLUX over sky-line-free pixels (from the
palace mask, resampled onto the file's own wavelength grid).  Fibers are
sorted by that continuum level; around each of the ``-low``/``-high``
percentile ranks, a window of the ``-navg`` fibers whose rank is closest to
that target is combined pixel-by-pixel with a sigma-clipped mean
(``astropy.stats.sigma_clipped_stats``).  The ``-low`` window is the sky
estimate; the ``-high`` window is a bright/science-like reference; their
difference is the sky-subtracted result.  Averaging several fibers per
percentile trades spatial locality in rank-space for lower per-pixel noise
(~1/√navg for clean pixels); sigma clipping keeps a fiber with a faint
source or a defect from dominating the average.

**Output:**

A FITS file with WAVE, SCI (robust mean of the high-percentile window,
reference only), SKY (robust mean of the low-percentile window, the sky
model), FLUX (SCI − SKY, sky-subtracted), and DRP_ALL (one row per
exposure: filename, expnum, exptime, obstime, mjd, fiber counts and IDs
used, RA/Dec, spectrograph IDs, and continuum flux for both windows)
extensions.

**See Also:** :doc:`api/SkySubSci/index`


SummarizeSciSky.py
------------------

Drpall-driven, remote-friendly version of ``SkySubSci.py``: runs the same
science-fiber sky estimate over a range of exposure numbers selected from a
drpall table (as ``SummarizeCframe.py`` does) instead of an explicit file
list, so it can be run unattended over many exposures (e.g. at Utah).

**Usage**::

    SummarizeSciSky.py [-emin 900] [-ver 1.3.2] [-drp_all FILE]
                       [-low 10] [-high 90] [-navg 10] [-sigma 3.0]
                       [-maxiters 5] [-mask FILE] [-stat median|mean]
                       [-out ROOT] exp_start exp_stop [delta]

**Arguments:**

exp_start
    Starting exposure number.

exp_stop
    Stopping exposure number.

delta
    Process every delta-th exposure in range (default 1).

**Options:**

-emin N
    Minimum exposure time to include (default 900).

-ver VER
    DRP version, used to locate ``drpall-VER.fits`` (default 1.3.2).

-drp_all FILE
    Explicit drpall table to read instead of ``drpall-VER.fits`` (FITS, or
    ascii if the name contains ``txt``/``.tab``).

-low PCT, -high PCT, -navg N, -sigma S, -maxiters K, -mask FILE, -stat STAT
    Same meaning as in ``SkySubSci.py``.

-out ROOT
    Output filename root.  Default:
    ``SummarizeSciSky_<ver>_<exp_start>_<exp_stop>_<delta>``.

**Description:**

Exposure selection and file-path resolution follow the same pattern as
``SummarizeCframe.py``'s ``read_drpall``/``select_exps``/``find_top``, but
are re-implemented locally here rather than imported from ``SumCframe.py``
(deprecated -- see ``deprecated/SumCframe.py``), so this script has no
dependency on the optional ``dask`` package that ``SumCframe.py`` requires
for an unrelated function.  The per-exposure algorithm is identical to
``SkySubSci.py``; the rank-window/robust-mean helpers themselves are no
longer duplicated per script -- both ``SkySubSci.py`` and
``SummarizeSciSky.py`` import the canonical ``_rank_window``/``_robust_mean``
from ``SummarizeCframe.py``, which also uses them for its own
``-by fiber`` selection mode (see :doc:`summarize`).  Rather than building
a fresh metadata table, the
calculated values (continuum flux, fiber IDs used, positions, spectrograph
IDs) are added as new columns directly onto the selected drpall rows, which
become the DRP_ALL extension of the output.

**Output:**

Same extension structure as ``SkySubSci.py`` (WAVE, SCI, SKY, FLUX,
DRP_ALL), with DRP_ALL being the selected drpall rows plus the added
columns.

**See Also:** :doc:`api/SummarizeSciSky/index`


SkySubPatch.py
--------------

Sky-subtract an lvmCFrame using a sky taken from a patch of the science
field itself, with the sky adjusted separately for every fiber, and write
an lvmSFrame-layout file with each fiber's own sky in its SKY extension.
Still under development.

**Usage**::

    SkySubPatch.py [-h] [-mode groups|center|slit] [-lsf broaden|none]
                   [-nbkg N] [-bkg FILE] [-ngroup N] [-mask FILE]
                   [-np N] [-out ROOT] cframe [cframe ...]

**Arguments:**

cframe
    One or more lvmCFrame FITS files, each processed independently.

**Options:**

-mode M
    Calibration unit (default ``groups``): ``groups`` = non-overlapping
    groups of about ``-ngroup`` neighbouring fibers on the sky, each fiber
    getting an inverse-distance average of its 3 nearest groups'
    calibrations; ``center`` = every fiber calibrated on itself plus its
    ``ngroup``-1 nearest neighbours (slowest); ``slit`` = blocks of
    ``ngroup`` consecutive fiberids within a spectrograph.

-lsf L
    ``broaden`` (default): at each wavelength broaden whichever of the fiber
    and the sky is sharper (the fiber's FLUX, IVAR and LSF are updated
    where the fiber is broadened); ``none``: no LSF matching.

-nbkg N
    Fibers in the background patch (default 20).

-bkg FILE
    Use these fiberids (one per line, or a table with a ``fiberid``
    column) as the background instead of the automatic choice.

-ngroup N
    Fibers per calibration group (default 7: a fiber and its surrounding
    hexagon).

-mask FILE
    Sky-line mask (``palace_make_mask.py`` output; default
    ``data/sky_mask.fits``).

-np N
    Parallel worker processes (default 8).

-out ROOT
    Output root (default ``lvmSFrame-<exposure>.patch`` in the current
    directory; with several inputs ``ROOT_<exp>``).

**Description:**

1. *Background patch.*  Every science fiber's nebular lines (Hα,
   [N II] 6583, [S II] 6716/6731, [O III] 5007, plus Hβ when the Moon is
   below the horizon; never [O II]) and bright sky lines are fitted with
   Gaussians (``lvm_gaussfit.py``'s fitter).  After screening out fibers
   with masked pixels, bright continuum (stars) or anomalous sky-line
   throughput, the patch is the compact group of ``-nbkg`` fibers in one
   spectrograph with the lowest emission that is faint in every scored
   line.  The background is the mean of the patch spectra.

2. *Calibration.*  For each calibration unit (``-mode``), the unit's
   combined spectrum is compared with the background on sky-line pixels
   only -- nebular lines (including [O I] 6300/6364, which can come from
   the source as well as the sky) and the b/r (5750-5810 Å) and r/z
   (7450-7650 Å) arm boundaries, where the flux calibration is
   problematic, are excluded.  Three things are fitted:

   - a wavelength shift (constant in b, quadratic in wavelength in r and z);
   - the relative line-spread function, by two-sided kernel matching with
     ``lsf_kernel.py``;
   - a throughput factor for each arm, after the LSF matching.

3. *Per-fiber sky.*  Each science fiber's sky is throughput × (background
   continuum + background sky lines shifted and, where the background is
   the sharper, broadened).  Where the *fiber* is the sharper, the fiber
   is broadened instead.  Fibers that are not good science fibers (SkyE,
   SkyW, standards, ``fibstatus`` ≠ 0) get the uncorrected background.

The throughput factors are fitted on the sky lines only, but applied to
the whole sky, continuum included -- i.e. the difference is assumed to be
a flat-field-like effect that dims lines and continuum alike.

**Output:**

An lvmSFrame-layout file: FLUX (= the possibly broadened CFrame FLUX
minus SKY), IVAR (including the sky variance), MASK, WAVE, LSF (updated
where a fiber was broadened), SKY, SKY_IVAR, FLUXCAL_*, SLITMAP, plus

- BACKGROUND -- the one-dimensional background spectrum (WAVE, BKG,
  BKG_ERR, BKG_LSF);
- BKGFIBERS -- the patch fibers;
- CALIB -- per fiber: calibration unit, wavelength shift, throughput
  (``thr_b``/``thr_r``/``thr_z``) and the fraction of pixels in which the
  fiber was broadened in each arm.

The PRIMARY header gains ``SSP*`` keywords recording the mode, the
patch's spectrograph, and how the patch was chosen.

**Results so far:** on three Vela exposures (9087, 14964, 16998) the
sky-line residual (rms over the line / line peak) is about 0.010, compared
with about 0.019 for the same patch sky subtracted without per-fiber
adjustments; the photon-noise floor is 0.006-0.008.  An exposure takes
about 2 minutes with 8 processes.  The default ``groups``/``broaden``
setting is the one that has been tested; ``center`` and ``slit`` are
less exercised.

**Caveats:**

- FLUX is emission in *excess* of the patch.  In Vela the patch still
  contains roughly 70-80% of the field-median Hα and 30-60% of the
  field-median [O III] 5007.
- Kernels only broaden, so wherever the patch is broader than a fiber,
  that fiber's spectrum is degraded to the patch's resolution.
- Run ``QualSFrame.py`` on the output and on the DRP SFrame in separate
  directories: its map files are named by exposure number only.

**See Also:** :doc:`api/SkySubPatch/index`


lsf_kernel.py
-------------

Module (no command line) that fits and applies an empirical,
wavelength-dependent *relative* line-spread-function kernel between two
LVM spectra, so that the sharper one can be broadened to match the other
before sky lines are subtracted.  Used by ``SkySubPatch.py``.

It is adapted from Ivan Katkov's ``sky_decomp.lsf_surface_iterative`` in
the lvmsky repository (reference commit e30b730) but is an independent,
simplified re-implementation that needs only numpy and scipy.  Per arm,
the kernel has 11 taps whose weights vary with wavelength as B-splines
(constant in b); it is constrained to be non-negative, to sum to one (so
it conserves flux) and, by default, to be single-peaked, and is solved
as a quadratic program with ``scipy.optimize.minimize`` (SLSQP).
Because a kernel can only broaden, ``two_sided_match()`` fits both
directions and at each wavelength broadens whichever spectrum is sharper.

Primary routines: ``fit_kernel``, ``apply_kernel``, ``kernel_moments``,
``two_sided_match``.

The sky-line list ``SkySubPatch.py`` uses to find sky-line pixels,
``data/lvm_sky_lines_all.dat``, is also vendored from lvmsky (file commit
677c304); its header records the provenance.

**See Also:** :doc:`api/lsf_kernel/index`


Typical Workflow
----------------

1. Sky-subtract one or more CFrames with the default settings::

       SkySubPatch.py lvmCFrame-00016998.fits

   The printed summary gives the patch's spectrograph and, for each
   scored line, the patch's flux as a fraction of the field median -- how
   much source emission is being subtracted from every fiber.

2. Inspect the result: which fibers formed the patch (BKGFIBERS), the
   background spectrum (BACKGROUND), and each fiber's shift, throughput
   and broadening (CALIB).  Throughput maps from CALIB show the smooth,
   flat-field-like pattern the method corrects for.

3. Compare with the DRP's SFrame for the same exposure by running
   ``QualSFrame.py`` on each, in separate directories.  Its sections on
   the line emission in the subtracted sky, the continuum and the
   sky-line residuals use only FLUX, SKY and IVAR, so the two reports
   are directly comparable (see :doc:`data_quality`).


See Also
--------

- :doc:`summarize` - ``SummarizeCframe.py``, whose drpall selection
  logic ``SummarizeSciSky.py`` mirrors
- :doc:`data_quality` - ``QualSFrame.py``, used to compare the result
  with the DRP
- :doc:`sky_method_eval` - ``SkySub_eval.py``, which reads the
  ``SkySubSci.py``/``SummarizeSciSky.py`` output
- :doc:`api/SkySubSci/index` - API documentation
- :doc:`api/SummarizeSciSky/index` - API documentation
- :doc:`api/SkySubPatch/index` - API documentation
- :doc:`api/lsf_kernel/index` - API documentation
