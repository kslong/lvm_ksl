Evaluating Sky-Subtraction Methods
==================================

This page is part of :doc:`sky_subtraction`.

The tools here measure how well a
sky-subtraction method removes the sky: ``SkySub_eval.py`` compares the
output files of the methods in :doc:`sky_methods` (and the
science-field methods in :doc:`sky_from_science`) side by side, and
``sky_residual_eval.py`` scores any predicted sky spectrum against the
observed one, independent of how the prediction was made.


SkySub_eval.py
--------------

Evaluate sky subtraction quality for one or more output FITS files
produced by SkySubOrig, SkySubDrp, SkySubDev1, SkySubDev2, or SkySepESO.
When multiple files are given they are overlaid in the same figures for
direct method comparison.

**Usage**::

    SkySub_eval.py [wmin wmax] [-num N] [-out outroot] filename [filename ...]

**Arguments:**

filename
    One or more SkySub output FITS files to evaluate.

wmin, wmax
    Wavelength range in Angstroms for the spectral overview panels
    (defaults: 3600 and 9800).

**Options:**

-num N
    Overlay N randomly selected individual spectra on the band panels
    (default: 20; 0 = band only).

-out outroot
    Combine all files into one HTML file with this root.  Without ``-out``,
    each file produces its own ``<stem>_eval.html``.

**Description — HTML output:**

The output HTML file is titled "Sky Subtraction Quality Check" and is
organized into two named sections: **Sky Line Subtraction** (Figures 1–4
plus a statistics table) and **Continuum Separation** (Figures 5–6 plus a
second statistics table).  In total there are six interactive Plotly
figures and two inline statistics tables.

*Figure 1 — spectral overview (3 panels, linear scale):*

All three panels share a y-range of −1×10⁻¹⁴ to 1×10⁻¹³
erg s⁻¹ cm⁻² Å⁻¹, chosen to make sky-subtracted residuals visible.
Panel 1 shows FLUX + SKY (before subtraction), Panel 2 shows FLUX
(after subtraction), and Panel 3 shows the SKY model.  Each panel shows
the median and 10th/90th percentile band together with up to N individual
spectra in light grey.  One colour per input file.  This figure sits above
the "Sky Line Subtraction" section header (it isn't part of either named
section).

*Figure 2 — residual histograms:*

For each of three diagnostic windows ([OI] 5577 Å: 5560–5594 Å,
[OI] 6300 Å: 6280–6320 Å, IR OH: 9300–9500 Å) the distribution of the
per-pixel **HF (high-frequency, continuum-subtracted) residual** — not
raw FLUX — is plotted with a Gaussian overlay.  Using the HF residual
means leftover continuum in the window does not bias the reported
median/NMAD away from genuine sky-line residuals.  A statistics box (N,
median, NMAD, skewness) appears inside each panel; skewness is computed
on the same median ± 5·NMAD clipped range shown in the plot, so it isn't
dominated by a handful of outliers invisible in the display.

*Figure 3 — diagnostic window median spectra:*

For each diagnostic window a wider search region (5400–5750,
6100–6500, 9000–9800 Å) is plotted showing the median and 10th/90th
percentile band for the original (dotted) and sky-subtracted (solid)
spectra.  Red shading marks the diagnostic (signal) window; green shading
marks mask-selected sky-line-free pixels used to estimate the noise floor.

*Statistics table:*

An HTML table between Figure 3 and Figure 4 reports per-file per-window
the noise floor, sky-line RMS before and after subtraction, and
percentiles of the HF RMS ratio (see below).  The same information is
printed to the terminal.

*Figure 4 — HF RMS ratio per spectrum:*

For each diagnostic window the noise-corrected HF RMS ratio is plotted
against spectrum index, or against MJD (recomputed precisely from
``obstime`` rather than the truncated-integer ``mjd`` column) when there
are more than 100 spectra and every overlaid file has that information —
otherwise all files fall back to spectrum index so every trace stays on a
common axis.  Points are plotted as markers (not connected lines).  The
three panels share a linked x-axis so zooming one panel pans all three.
See the algorithm description below.


**HF RMS quality metric — algorithm:**

The central diagnostic is a per-spectrum, per-window noise-corrected
high-frequency (HF) RMS ratio.  The following steps are applied
independently to both the original spectra (FLUX + SKY) and the
sky-subtracted spectra (FLUX).

*Step 1 — high-frequency residuals.*

A smooth estimate of the underlying continuum is subtracted from each
spectrum::

    hf(λ) = flux(λ) − smooth(λ)

The smooth is a Gaussian-weighted running average with σ = 50 pixels
(≈ 25 Å at the LVM pixel scale of 0.5 Å pixel⁻¹).  Sky-line pixels
(mask = 0 from the palace mask) are excluded by a weighted convolution::

    smooth(λ) = Σ_λ' [ flux(λ') · w(λ') · G_σ(λ−λ') ]
              / Σ_λ' [ w(λ') · G_σ(λ−λ') ]

where w = 1 for mask-clean pixels and w = 0 for sky-line pixels.
This prevents bright sky lines from leaking into the smooth estimate and
artificially reducing the HF amplitude in the diagnostic window.

After this step, ``hf`` retains only structure narrower than ≈25 Å —
the scale of individual sky emission lines — while broad continuum
mismatches from polynomial or B-spline fitting are suppressed.

*Step 2 — why RMS, not NMAD.*

Sky lines are spatially sparse.  [OI] 5577 spans roughly 3–4 pixels in
the 34-pixel diagnostic window; individual OH lines in the IR cover a
similarly small fraction.  NMAD (the median absolute deviation) is
dominated by the majority of clean, noise-floor pixels and is nearly
insensitive to a small number of bright outliers.  RMS responds to the
squared amplitude, so even 2–4 bright sky-line pixels contribute
proportionally to their intensity — exactly the signal we want to track.

*Step 3 — noise floor subtraction.*

Even after perfect sky subtraction the HF residuals carry photon and
detector noise.  This noise contributes to the RMS measured in the
diagnostic window.  A nearby background region — mask-selected
line-free pixels (mask = 1) within the broader search range but outside
the diagnostic window — provides an independent noise estimate
``rms_bg``.  The sky-line contribution is then isolated in quadrature::

    sky_rms = sqrt( max(0, rms_diag² − rms_bg²) )

For a perfect subtraction ``sky_rms`` approaches zero regardless of the
noise level.  The assumption is that the noise is approximately stationary
across the search region, which is valid as long as the continuum level
does not change drastically within a few hundred ångströms.

*Step 4 — the ratio and its interpretation.*

For each spectrum::

    ratio = sky_rms_sub / sky_rms_orig

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - ratio
     - Meaning
   * - 0
     - Sky lines completely removed.
   * - 0 – 0.5
     - Substantial improvement; more than half the sky-line power removed.
   * - ≈ 1
     - No improvement; subtraction left sky-line power unchanged.
   * - > 1
     - Sky lines made worse; subtraction introduced new residuals.
   * - NaN
     - ``sky_rms_orig`` ≈ 0; no detectable sky-line power in the original
       spectrum — excluded from statistics.

*Step 5 — summary statistics.*

The statistics table reports the following per file and per window:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Description
   * - ``noise_med``
     - Median of ``rms_bg`` across all spectra — the noise floor
       (erg s⁻¹ cm⁻² Å⁻¹).  Tells you how sensitive the measurement is.
   * - ``sky_orig_med``
     - Median ``sky_rms`` in the original spectra — the typical sky-line
       amplitude before subtraction.
   * - ``sky_sub_med``
     - Median ``sky_rms`` after subtraction — the residual sky-line amplitude.
       Should be much smaller than ``sky_orig_med`` for a good method.
   * - ``p50 ratio``
     - Median ratio across all spectra; characterises typical performance.
   * - ``p90 ratio``
     - 90th-percentile ratio; indicates how the worst 10% of spectra behave.
   * - ``p95 ratio``
     - 95th-percentile ratio; the tail of poor performance.
   * - ``frac < 0.5``
     - Fraction of spectra with ratio < 0.5, i.e. where sky-line power was
       more than halved.  A simple pass/fail rate at this threshold.


**Continuum Separation section — Figures 5 & 6 and the second statistics table:**

Unlike the HF RMS metric above (which targets sky *line* residuals), this
section evaluates *continuum* separation quality, using three
spectrograph-arm bands with the B/R/Z overlap zones and outer edges
excluded — B (3650–5775 Å), R (5800–7520 Å), Z (7570–9600 Å).  Both
figures share the same three-row layout and row order:

- **Row 1 (SCI)** — DRP_ALL['SCI_MED_<arm>']: the per-spectrum median of
  (raw science flux − sci_cont) in clean (sky-line-free) pixels of the
  *raw, pre-subtraction* science spectrum.  This is the actual
  science-side continuum-fit-quality test.
- **Row 2 (SKY)** — DRP_ALL['SKY_MED_<arm>']: the same, but for the raw
  sky-telescope spectrum whose continuum fit was actually used in SKY.
  This is the actual sky-side continuum-fit-quality test, and is the one
  that matters directly for the final result's continuum level.
- **Row 3 (Subtracted)** — the per-spectrum median of the *raw*
  sky-subtracted FLUX (not the HF residual used by Figure 2) in clean
  pixels.  This is the final, post-subtraction leftover signal.  It does
  **not** test continuum-fit quality — the science-side continuum fit is
  only ever used to derive the bisection line-scale target and never
  enters the subtracted result — so this row instead reflects real source
  (stellar) continuum entangled with any net error in the sky-side
  continuum estimate.  Comparing this row across overlaid methods on the
  same input isolates sky-continuum quality specifically, since real
  source continuum is identical across methods.

Rows 1/2 require the ``SCI_MED_<arm>``/``SKY_MED_<arm>`` DRP_ALL columns
already present in the input file (written by SkySubOrig.py, SkySubDev1.py,
or SkySepESO.py); a file without them (SkySubDev2.py/SkySubDrp.py output,
or an older file predating this feature) simply leaves those panels empty.

*Figure 5* histograms all three rows (pooled/distribution across spectra,
same N/median/NMAD/skew annotation style as Figure 2).  *Figure 6* plots
the same three quantities per spectrum against spectrum index or MJD
(same x-axis logic as Figure 4), so specific bad exposures/fibers are
visible rather than just an aggregate number.  Both figures use a
colour-coded suptitle above the plot grid (rather than a floating legend,
which would otherwise sit on top of a subplot in the 3×3 grid).

The second statistics table (screen + HTML) reports, per file/arm/kind
(SCI, SKY, Subtracted), the N/median/NMAD/skew of the row's distribution.

**DRP_ALL write-back:** the per-spectrum Subtracted-row statistics
(``resid_med_<arm>``, ``resid_nmad_<arm>``, ``resid_rms_<arm>``,
``resid_skew_<arm>``) are written back into each evaluated file's own
DRP_ALL table, in place — the same convention used by
``GetSkyCont_eval.py`` for its own fit-quality columns.  (The SCI_MED/
SKY_MED columns themselves are not written by this script — they already
exist in the input file, written by whichever SkySub* script produced it.)

**Output:**

An HTML file named ``<stem>_eval.html`` (one per input file without
``-out``; a single combined file when ``-out outroot`` is given).  Any
evaluated file with a DRP_ALL table is updated in place with the
Subtracted-row statistics above; no backup is created.

**See Also:** :doc:`api/SkySub_eval/index`


Comparing Sky Subtraction Methods
----------------------------------

A typical workflow for running the XCframe methods on a single file and
comparing the results::

    # 1. Build the sky-line mask (if not already present)
    palace_make_mask.py XCframe_file.fits

    # 2. Run the subtraction methods (SkySubRun.py avoids tracking each
    #    routine's own output-naming convention by hand)
    SkySubRun.py -routine orig  XCframe_file.fits
    SkySubRun.py -routine drp   XCframe_file.fits
    SkySubRun.py -routine dev1  XCframe_file.fits -mask sky_mask.fits
    SkySubRun.py -routine dev2  XCframe_file.fits
    SkySubRun.py -routine dev3  XCframe_file.fits
    SkySepESO.py XCframe_file.fits -delta 50   # slow; use -delta for a quick look

    # 3. Evaluate and compare in a single HTML file
    SkySub_eval.py -out compare \
        sky_runs/orig/XCframe_file_orig_farlines_nearcont.fits \
        sky_runs/drp/XCframe_file_drp_farlines_nearcont.fits \
        sky_runs/dev1/XCframe_file_dev1_farlines_nearcont.fits \
        sky_runs/dev2/XCframe_file_dev2_scilines_nearcont.fits \
        sky_runs/dev3/XCframe_file_dev3_farlines_nearcont.fits \
        XCframe_file_eso_farlines_nearcont.fits

Open ``compare_eval.html`` in a browser.  Figure 1 overlays the median
spectra for all methods; Figure 4 (Sky Line Subtraction section) shows
the per-spectrum HF RMS ratio for each diagnostic window; Figures 5/6
(Continuum Separation section) show per-arm continuum-fit quality for the
methods that record it (SkySubOrig, SkySubDev1, SkySubDev2, SkySepESO) —
together these make it straightforward to identify which method best
suppresses sky lines *and* which best separates continuum from lines for
the observation.

This tells you how well a method suppresses *airglow* residuals, but
nothing about whether it also distorts or destroys real *nebular* signal
on the science fiber — a method that oversubtracts is invisible to a
generic sky-residual metric, since the airglow residual it's built from
looks equally good either way.  For that question, see :doc:`sky_nebular_eval`.


sky_residual_eval.py
--------------------

A method-agnostic sky-subtraction evaluator, complementary to
SkySub_eval.py above.  Where SkySub_eval.py overlays several methods'
own output files for a visual side-by-side comparison, sky_residual_eval.py
takes any single (observed, model) pair — whatever produced it, a full
ESO/PALACE/SkyDecomp decomposition or something as simple as a scaled
sky-fiber spectrum — and reduces it to a fixed set of numeric quality
metrics, separately for the continuum and for individual airglow lines,
suitable for batch comparison across many exposures or fibers at once.

**Usage**::

    sky_residual_eval.py filename [-mask mask.fits] [-nproc N]
                         [-out ROOT] [-plotdir DIR]

**Arguments:**

filename
    FITS file with WAVE, FLUX, and SKY extensions (the convention used by
    SkySubOrig/Drp/Dev1/Dev2.py and read by SkySub_eval.py above).  WAVE
    is 1-D; FLUX/SKY are 1-D (single spectrum) or 2-D, ``n_rows x
    n_wave``.  An IVAR extension is used if present.  In this file
    convention FLUX is already sky-subtracted, so the observed spectrum
    is reconstructed internally as FLUX+SKY, with SKY as the model.

**Options:**

-mask mask.fits
    A palace_make_mask.py output (WAVE/MASK extensions).  Defaults to
    ``data/sky_mask.fits``, the same default GetSkyCont.py and
    SkyObsESOCompare.py use.

-nproc N
    Worker processes for batch rows (default 1).

-out ROOT
    Output filename root; writes ``<ROOT>_summary.fits`` and
    ``<ROOT>_lines.fits`` (default: the input file's stem).

-plotdir DIR
    Directory for the three summary plots described below (default
    ``plots_sky_resid``).

**Description:**

For every row, residual = Observations − Sky is decomposed two ways:

*Continuum bands* (B/R/Z, arm-overlap zones excluded) — over clean
pixels identified from the mask: an offset, a dimensionless power-law
index mismatch (``CONT_ALPHA``, treating the residual as a small
``(wavelength/reference)**alpha`` shape correction to the model rather
than a flux-per-Angstrom slope, so it is comparable across bands and
exposures of very different brightness), a fit-quality ratio
(``CONT_FIT_QUALITY`` = NMAD/NOISE_PROXY — not a formal ivar
chi-square, which saturates uselessly large for real sky-subtraction
residual since it routinely exceeds the formal photon-noise floor), and
the fraction of clean pixels with ``|residual|`` at or below two fixed
absolute flux levels (1e-14 and 1e-15 erg/s/cm^2/Angstrom).

*Individual lines* — a fixed default list of 16 airglow lines (the
DRP's own ``REF_SKYLINES`` plus 9 lines from ``sky_gaussfit.py``'s
``SKY_LINES``).  For each line, a local Gaussian is fit to the MODEL
(not the data) to find its own amplitude, center, and width, and the
residual in that window is fit against the analytic first-order
derivative of that Gaussian in amplitude, center, and width.  The three
resulting numbers are the physical corrections needed to make the model
match the data: an amplitude ratio (``AMP_RATIO`` — predicted/measured,
used instead of a flux-unit amplitude difference so it is comparable
across lines of very different brightness), a wavelength-registration
offset (``DELTA_LAM``, Angstroms), and an LSF-width offset
(``DELTA_SIGMA``, Angstroms — the number that answers whether a small
LSF mismatch is contributing to the residual).  The same fixed-threshold
quality fractions as the continuum side are also computed, but over the
mask's own line-affected pixels rather than the line list, since the
line list is deliberately too short a set for a fraction to be more than
a coarse step function.

**Output files:**

``<ROOT>_summary.fits``
    One row per input spectrum with all of the per-band continuum and
    line-aggregate quantities above.

``<ROOT>_lines.fits``
    One row per (spectrum, line), with each line's individual fit
    results.

Three PNG summary plots are also written to ``-plotdir`` (default
``plots_sky_resid/``):

``<ROOT>_frac_summary.png``
    Reverse-cumulative distributions of the quality fractions, one row
    for continuum and one for lines, one column per arm.

``<ROOT>_continuum_summary.png``
    Histograms of the continuum diagnostics (offset, scatter, alpha,
    fit quality), one row per metric, one column per arm; a metric with
    no data in any arm (e.g. fit quality is always available, but a
    metric that genuinely has none) is dropped rather than left blank.

``<ROOT>_lines_summary.png``
    Histograms of the per-band line diagnostics (amplitude ratio,
    registration bias, LSF-width bias), same layout.

In all three plots, a metric's row shares one x-axis range across all
three arms (from the metric's pooled 1st–99th percentile, not the
min/max), so the relative width of the distribution in each arm is
directly comparable, and a dashed reference line marks the "model
matches data exactly" value (0 for offsets/biases, 1 for ratios).

**Example**::

    sky_residual_eval.py XCframe_file_drp_farlines_nearcont.fits -nproc 8

This writes ``XCframe_file_drp_farlines_nearcont_summary.fits``,
``XCframe_file_drp_farlines_nearcont_lines.fits``, and the three PNGs
under ``plots_sky_resid/``, evaluating every row of the input file
against ``data/sky_mask.fits`` in parallel across 8 worker processes.

``analyze_sky_residual``/``analyze_sky_residuals`` are also directly
importable for use outside the command line -- see the API reference
for their full parameter and return-value documentation.


See Also
--------

- :doc:`sky_nebular_eval` - evaluating methods by nebular-line recovery
- :doc:`sky_methods` - the methods being evaluated
- :doc:`api/SkySub_eval/index` - API documentation
- :doc:`api/sky_residual_eval/index` - API documentation
