Spectral Fitting -- Local
==========================

The lvm_ksl package provides tools for fitting emission lines in LVM
spectra. These tools can fit Gaussian profiles to standard emission
lines across all fibers in an RSS file, or perform more detailed fits
to individual spectra.


Overview
--------

The spectral fitting tools include:

- ``lvm_gaussfit.py`` - Fit standard emission lines across an RSS file
- ``lvm_snrfit.py`` - Fit the Mappings-model line set expected in SNRs,
  jointly for lines too close together to fit independently
- ``sky_gaussfit.py`` - Fit nebular and airglow lines fiber-by-fiber in SFrame files
- ``lvm_line_profile.py`` - Compare Gaussian vs Moffat airglow line profiles on raw sky spectra
- ``lvm_double.py`` - Fit single or double Gaussian profiles to a line
- ``lvm_triple.py`` - Fit up to triple Gaussian profiles
- ``lvm_flux.py`` - Calculate fluxes using parameters from other fits
- ``FlattenSpec.py`` - Separate emission lines from an underlying
  (typically stellar) continuum in an extracted spectrum

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


lvm_snrfit.py — SNR Emission Line Fitting with Joint Blend Handling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fits the fuller, SNR-relevant emission-line set in
``data/mappings_snr_lines.txt`` (36 lines, wavelengths from a Mappings
shock model rather than the DAP's own approximate values) to a single
spectrum or RSS file, the way ``lvm_gaussfit.py`` does for its smaller,
hardcoded line list. Reuses ``lvm_gaussfit.py``'s file I/O, fiber
selection, batch output, and plotting machinery by import
(``do_all``/``do_individual`` now accept a ``do_one_func`` argument for
exactly this) rather than duplicating it -- only the per-spectrum fitting
routine (``do_one``) and the model itself are new.

**Why a separate tool:** four of the added lines are close enough in
wavelength that, once physically plausible line broadening is allowed
for, an independent per-line local-background fit (``lvm_gaussfit.py``'s
approach) would let each line's neighbour bias its background estimate.
``lvm_snrfit.py`` fits those pairs jointly instead, with one shared local
background.

**Usage**::

    lvm_snrfit.py [-h] [-lmc] [-smc] [-v vel] [-stype SOURCE] [-out root]
                  [-plot] [-lines file.txt] [-min_sig N] filename ...

**Options:**

-h
    Print help and exit.

-lmc, -smc, -v vel
    Same velocity-offset options as ``lvm_gaussfit.py``.

-stype SOURCE|BACK
    Same as ``lvm_gaussfit.py``, for extracted text spectra.

-out root
    Root name for the output file.

-plot
    Save a per-line fit-quality plot for each spectrum.

-lines file.txt
    Reference line list (default ``data/mappings_snr_lines.txt``).
    Passing a file with a subset of rows restricts which lines are fit --
    there is no separate on/off flag per line.

-min_sig N
    Significance threshold (flux/eflux) for keeping a blend component
    before it's dropped as a non-detection (default 2.0).

**Arguments:**

filename
    An SFrame-compatible FITS file, or one or more ascii tables with
    WAVE/FLUX[/ERROR] columns.

**The line-list reference table, data/mappings_snr_lines.txt:**

Vendored from a Mappings v100 shock model run
(``~/Projects/Mappings26/test/Mapping_v_DAP.txt``), matched against the
DAP's own line set where possible (``DAP2tab.py``'s ``get_radec_fluxes()``
line/name table follows the same naming convention -- see :doc:`dap`).
Columns: ``Wave`` (Mappings rest wavelength, authoritative), ``Ion``,
``Kind``/``Accuracy`` (carried through from Mappings, informational only),
``gauss_name`` (the output column-name suffix), and ``group`` (blend-group
tag, ``-`` for an independent singlet).

Five groups need a joint fit (four new, plus the pre-existing ``oii``
doublet, still fit via ``lvm_gaussfit.fit_double_gaussian_to_spectrum``
unchanged):

.. list-table::
   :header-rows: 1
   :widths: 16 34 8 20

   * - Group
     - Members
     - Gap
     - FWHM handling
   * - ``oii``
     - [OII] 3726.03 / 3728.82
     - 2.79 A
     - shared (existing)
   * - ``hei_hI``
     - HeI 3888.64 / HI (H8) 3889.06
     - 0.40 A
     - shared, centers tied
   * - ``neiii_hepsilon``
     - [NeIII] 3967.47 / Hepsilon 3970.08
     - 2.61 A
     - shared
   * - ``ni``
     - [NI] 5197.90 / [NI] 5200.26
     - 2.36 A
     - shared
   * - ``oii7320_caii``
     - [OII] 7319.99 / [CaII] 7323.89
     - 4.96 A
     - independent

Gap sizes look small in isolation, but the relevant comparison is to the
line's own FWHM once broadening is allowed for -- see the next section.

**Line-width bounds and the shared-vs-independent-FWHM choice:**

FWHM bounds come from velocity, not an arbitrary Angstrom scaling:
``V_INSTR`` (80 km/s, typical instrumental width) sets the minimum,
``V_MAX`` (200 km/s, the broadest line judged plausible for a shocked SNR
knot) sets the maximum. Comparing each blend pair's wavelength gap to its
FWHM at ``V_MAX`` (using the already-working ``oii`` doublet, gap/FWHM ~
1.0-1.1, as the "this is known to work" benchmark) is what decided which
groups share a single FWHM parameter versus fit two independent ones:
``hei_hI``, ``neiii_hepsilon``, and ``ni`` all come out at or beyond that
benchmark's blend severity at 200 km/s, so a shared FWHM is doing real
stabilizing work, not just convenience; ``oii7320_caii`` sits at a
comparable margin to ``oii`` itself, different species, so independent
widths were judged numerically tractable.

**HeI/HI (H8) is a special case:** at 0.4 A separation, two independent
centers plus two independent fluxes are degenerate -- the fit can trade
flux between the components while barely changing chi-square, producing
unphysical (even negative) individual fluxes. ``fit_blend_to_spectrum``'s
``tie_centers`` option locks the second line's center to the first's plus
their fixed Mappings wavelength offset, so only one systemic shift is fit
and the only freedom left is how flux splits between the two -- and even
that split is checked (see below) and reported as one combined feature
under the dominant line's name if it's still not reliably separable.

**Weak-line / degenerate-fit safeguard:** a blend component is dropped,
and the group refit as an ordinary singlet with the other line, if either:

- its fitted flux is not significant relative to its uncertainty
  (``-min_sig``), or
- flux1 and flux2 are found to be nearly perfectly anti-correlated
  (\|correlation\| >= 0.95) -- this catches high-S/N cases where each
  flux individually looks "significant" by its own marginal stderr, but
  the fit still can't actually tell the two components apart (their
  covariance reveals what the marginal error doesn't).

The dropped line's columns are reported as NaN rather than a spurious
joint-fit value, so a non-detection (or an irresolvable blend) can't
corrupt the shared background/width used for the line that *is*
measurable.

**Fitting-window sizing:** each line's window is ``n_vmax`` times its
``V_MAX``-broadened FWHM (with a floor), but capped at half the distance
to the nearest *other* line in the table (excluding its own blend-group
partner) -- found necessary during validation: without the cap, some
"safe" (non-blended) singlets still had large enough windows to overlap
a neighbour's and bias each other's local background (e.g. ``[FeX]``
6374.51 sitting only 10.7 A from ``[OI]`` 6363.78). ``[OI]`` 6300.3/6363.78
additionally exclude a small, fixed (un-redshifted) window matching
``lvm_gaussfit.do_one`` -- both sit on sky-subtraction-residual-prone
airglow lines "seldom subtracted correctly" -- and are exempted from the
neighbour cap, since giving up part of the window to that exclusion
already needs the room back.

**Validated (260823):** synthetic spectra with known injected fluxes
recover all lines except the genuinely-irresolvable HeI/HI pair (reports
correctly as one combined feature); a real Vela SNR shock spectrum
(``~/Projects/lvm_science/Vela/Shock_vela.ave_sum.txt``) fits all 36
lines with plausible flux/significance for every line except one --
``[OI]`` 6363.78, whose fit is numerically unstable on this particular
spectrum because of an unusually strong sky-subtraction residual there.
That instability was confirmed to be a **pre-existing** property of
``lvm_gaussfit.fit_gaussian_to_spectrum`` itself (reproduced with
``lvm_gaussfit.py``'s own original window/exclusion on the same file),
not something introduced by this tool.

**Example**::

    # Fit the full SNR line set to a Vela shock spectrum
    lvm_snrfit.py Shock_vela.ave_sum.txt

    # Restrict to a subset of lines
    lvm_snrfit.py -lines my_lines.txt Shock_vela.ave_sum.txt


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


Continuum / Line Separation
----------------------------

FlattenSpec.py
^^^^^^^^^^^^^^

Fits and subtracts a smooth continuum from a single extracted spectrum
(a ``WAVE``/``FLUX`` ascii table, e.g. a region-averaged
``*_ave_sum.txt`` file). The purpose is to separate the emission lines
of an emission-line region from any underlying continuum, which
usually arises from stellar contributions to the overall spectrum
(scattered/foreground starlight, an embedded stellar population, ...)
rather than from the nebular gas itself. Isolating that continuum is
useful both to characterize it on its own and to leave a clean,
continuum-free residual for line measurement.

**Usage**::

    FlattenSpec.py filename.txt [-mask mask.fits] [-kstep N]
                   [-niter N] [-nsigma N] [-out outfile]

**Options:**

-mask FILE
    ``palace_make_mask.py``-format FITS mask flagging sky-line-affected
    pixels (default: the vendored ``data/sky_mask.fits``, the same
    default ``PlotSpec.py``/``PlotSpecI.py`` use for their ``-mask``
    option -- see :doc:`sky_subtraction`).

-kstep N
    B-spline knot spacing in Angstroms (default 100).

-niter N
    Number of sigma-clip continuum-refit iterations (default 3; 0
    disables clipping -- see below).

-nsigma N
    Rejection threshold in robust sigma above the local residual
    median (default 4.0).

-out FILE
    Output filename (default ``<stem>_flat.txt``).

**Method:**

Pixels are excluded from the continuum fit if the sky-line mask flags
them, if the input table's own ``MASK`` column already flags them bad,
or if the sigma-clip (below) flags them as a bright emission line. The
surviving pixels are fit with ``GetSkyCont.py``'s plain B-spline design
matrix via non-negative least squares (no solar/MOON term -- this is a
science spectrum, not a sky spectrum), reusing that machinery directly
rather than duplicating it.

Because the B-spline coefficients are constrained non-negative, an
*unmasked* bright science emission line (H-alpha, [OIII] 5007, ...)
can only pull the local knot up, never down to compensate -- so a
single least-squares pass tracks the line as a spurious continuum
bump rather than ignoring it. Since this tool is meant to work on any
emission-line spectrum without maintaining a fixed line list (and
lines can be Doppler-shifted or broadened off their rest wavelengths,
e.g. in a shock), the fit instead iterates: after each pass, pixels
whose residual (``FLUX - CONT``) exceeds ``-nsigma`` robust-sigma
(MAD-based) above the local residual median are excluded, and the
continuum is refit without them. Only positive outliers are rejected,
since science emission adds flux rather than removing it. Iteration
stops early once a pass rejects no new pixels.

**Output columns** (existing columns are kept; these are added/updated)::

    CONT   fitted continuum
    RESID  FLUX - CONT
    MASK   original MASK bits, OR'd with bit 20 wherever the sky-line
           mask excluded the pixel and/or bit 21 wherever the sigma-clip
           excluded it as a bright emission-line pixel

**Example**::

    # Default: 3 sigma-clip iterations at 4-sigma
    FlattenSpec.py Shock_vela.ave_sum.txt

    # Disable sigma-clipping (single pass, matches a fixed line mask only)
    FlattenSpec.py Shock_vela.ave_sum.txt -niter 0

**Validated (260820):** run on a Vela shock spectrum
(``Shock_vela.ave_sum.txt``), the single-pass fit produced a
spurious continuum bump centered on every bright unmasked line
(H-alpha, [OIII] 4959/5007, [OII] 3727, ...) -- e.g. the continuum
under [OIII] 5007 peaked at ~9x its surrounding baseline. With the
default sigma-clip enabled, those bumps disappear and the continuum
tracks the true underlying baseline through ~90 rejected line-affected
regions (not just the handful of strong nebular lines -- also weaker
lines like [OIII] 4363 and a cluster of narrow features that were
sky-line residuals not fully caught by the sky mask).


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
- :doc:`api/lvm_snrfit/index` - API documentation
- :doc:`api/sky_gaussfit/index` - API documentation
- :doc:`api/lvm_line_profile/index` - API documentation
- :doc:`api/lvm_double/index` - API documentation
- :doc:`api/lvm_triple/index` - API documentation
- :doc:`api/lvm_flux/index` - API documentation
- :doc:`api/FlattenSpec/index` - API documentation
