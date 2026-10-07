Evaluating Methods by Nebular-Line Recovery
===========================================

This page is part of :doc:`sky_subtraction`.

Every generic sky-residual metric in :doc:`sky_method_eval` (``SkySub_eval.py``,
``sky_residual_eval.py``) measures how well *airglow* is suppressed. None
of them can tell you whether a method also distorts real *nebular*
emission on the science fiber — an oversubtracting method looks just as
good by those metrics as one that doesn't, since both are judged purely
on leftover airglow residual. This group of tools instead measures
recovery of known nebular-line physics directly: fixed atomic-physics
line ratios ([SIII] 9531/9069 = 2.44, [NII] 6583/6548 = 3.0), a
reddening-bounded ratio (Hβ/Hα ≤ 0.350, Case B no-reddening ceiling), and
scatter across repeat observations of the same tile — the same
philosophy as ``sky_residual_eval.py``'s method-agnostic design
(``[[feedback_method_agnostic_eval]]``), applied to line ratios instead
of residual shape. It also directly tests the lvmdrp assumption that
SKY_EAST/SKY_WEST themselves carry no nebular-line flux, an assumption
every method in :doc:`sky_methods` (this evaluation family included, for the far/near
sky-line component) implicitly depends on.

The tools split into two lines of investigation sharing a common line
catalog and velocity convention:

- **Does the sky itself leak nebular flux?** (``DecomposeCleanSky.py`` →
  ``sky_nebular_leak_eval.py`` → ``PlotNebularLeak.py``) — fits a
  nebula-free "clean sky" model to FLUX/SKY_EAST/SKY_WEST by masking out
  nebular-line windows before the fit, then measures whether real flux
  still leaked into those windows anyway.
- **Which SkySub* method best recovers real nebular flux?**
  (``SkySubRun.py`` → ``SkySubNebEval.py`` → ``PlotSkySubNebEval.py`` /
  ``PlotSkySubNebRun.py``) — fits the same nebular-line catalog directly
  on each method's sky-*subtracted* output and compares the measured
  ratios/repeat scatter across methods, either one exposure at a time
  (``PlotSkySubNebEval.py``) or aggregated across a whole run
  (``PlotSkySubNebRun.py``).

Both share ``sky_gaussfit.NEBULAR_LINES``/``resolve_nebular_lines()``
(see :doc:`spectral_fitting_local`) for the line catalog and Doppler-shift
convention, and ``sky_nebular_leak_eval.resolve_velocity()``/``LMC_VEL``/
``SMC_VEL`` for turning ``-v``/``-lmc``/``-smc`` (or, in
``SkySubNebEval.py``/``PlotSkySubNebEval.py``, a per-row DRP_ALL
``Redshift`` lookup — see below) into a systemic velocity::

    SummarizeCframe.py -by fiber
            |
            v
    XCframe summary FITS (WAVE, FLUX, SKY_EAST, SKY_WEST, DRP_ALL[, LSF])
            |
            +---------------------------------+
            |                                 |
            v                                 v
    DecomposeCleanSky.py                SkySubRun.py -routine
    (py_dev/ -- mask nebular lines,      {drp,orig,dev1,dev2,dev3}
     fit clean-sky continuum with              |
     SkyDecompLSFSurfaceIterative)              v
            |                            <routine> output FITS
            v                            (WAVE, FLUX, SKY, DRP_ALL)
    CleanSky_<expnum>.fits                     |
            |                                  v
            v                            SkySubNebEval.py
    sky_nebular_leak_eval.py              (fit NEBULAR_LINES on FLUX,
     (Gaussian-on-residual fit             DOUBLETS ratio checks,
      at each nebular-line window)         repeat_scatter across tileid)
            |                                  |
            v                                  v
    PlotNebularLeak.py                  PlotSkySubNebEval.py (one exposure)
     (per-line before/after,             or PlotSkySubNebRun.py (whole run)
      shared-systematic template          (pointing table / run-summary
      correction panels)                  table, ratio-vs-flux and repeat-
                                           group scatter -- one column
                                           per method)

DRP_ALL's own ``Survey``/``Redshift`` columns (``SummarizeCframe.py``'s
RA/Dec-based LMC/SMC/Plane/HighLat classification, 262/146/0/0 km/s) are
the per-row systemic-velocity default for ``SkySubNebEval.py``/
``PlotSkySubNebEval.py`` — important because a single file/run routinely
mixes rows from different surveys (e.g. a repeat-observation test set
spanning HighLat and SMC tiles together), so one global ``-v``/``-lmc``/
``-smc`` assumed for an entire run is wrong for whichever rows don't
match it. Getting this wrong silently mis-locates every nebular-line fit
window by the missed velocity's Doppler shift (comparable to or larger
than the fit window itself for LMC/SMC fields), which reads as "no
nebular signal detected" rather than "fit at the wrong wavelength" —
found on a real LMC exposure where every doublet ratio came back
SNR-gated as an apparent non-detection under the wrong assumed velocity.
Passing ``-v``/``-lmc``/``-smc`` explicitly still overrides the per-row
lookup, for the case where every row in a run is known to share one
real target velocity.


DecomposeCleanSky.py
--------------------

Decomposes one or more of FLUX (science fiber), SKY_EAST, and SKY_WEST
from an XCframe summary file into a nebula-free "clean sky" model, using
``sky_decomp.lsf_surface_iterative.SkyDecompLSFSurfaceIterative`` with the
known nebular emission lines excluded from the fit ("mask-and-wrap").
Lives in ``py_dev/`` (not Sphinx-API-documented — see :doc:`sky_models`
for that tier's conventions) since it live-imports the ``lvmsky`` repo's
``skysub`` package, the same pattern as ``py_dev/PredictSky.py``.

**Usage**::

    DecomposeCleanSky.py [-lvmsky_skysub PATH] [-py_progs_dir PATH]
                         [-row N | -expnum N] [-ext LIST]
                         [-v VEL] [-lmc] [-smc] [-output PATH]
                         fits_file

**Arguments:**

fits_file
    XCframe summary FITS file (WAVE, FLUX, SKY_EAST, SKY_WEST, DRP_ALL).

**Options:**

-row N / -expnum N
    Select the exposure by row index (default 0) or DRP_ALL ``expnum``
    (overrides ``-row``).

-ext LIST
    Comma-separated extensions to decompose (default:
    ``FLUX,SKY_EAST,SKY_WEST`` — FLUX is included by default too, as a
    same-exposure comparison point for how much nebular flux the sky
    telescopes see relative to the science fiber).

-v VEL / -lmc / -smc
    Nebular systemic velocity (km/s), same convention as
    ``sky_gaussfit.py`` (default: 0).

-output PATH
    Output FITS path (default: ``CleanSky_<expnum>.fits``).

-lvmsky_skysub PATH / -py_progs_dir PATH
    Paths to the ``lvmsky`` repo's ``skysub/`` directory and the
    ``lvm_ksl`` repo's ``py_progs/`` directory (defaults
    ``~/SDSS/lvmsky/skysub`` / ``~/SDSS/lvm_ksl/py_progs``).

**Description:**

For each requested extension: builds the nebular exclusion mask from
``sky_gaussfit.NEBULAR_LINES``, each window Doppler-shifted by the same
``zz = 1 + vel/3e5`` convention ``sky_gaussfit.py`` itself uses; zeroes
IVAR inside those windows (and at any non-finite flux pixel) before
fitting, so the fit engine's own
``good = isfinite(flux) & isfinite(ivar) & (ivar>0)`` never sees the
nebular pixels and none of its continuum/sky-line families can absorb
nebular flux; reconstructs the clean-sky model (``bestfit_lsf``) and the
full-array residual (observed − clean sky, including *inside* the
excluded windows — that residual there is exactly the nebular-leak
signal ``sky_nebular_leak_eval.py`` measures).

This is the "mask-and-wrap" half of the continuum/sky-line/nebular-line
separation effort; a later, more ambitious native nebular family solved
jointly inside SkyDecomp's own design matrix would slot in as a different
model-building step without changing the leak-detection tooling
downstream (which only ever needs a ``(wave, flux, model)`` triple).

**Output:**

A FITS file with, per requested extension: observed ``<EXT>``,
``<EXT>_BESTFIT`` (clean-sky model), ``<EXT>_RESID``, ``<EXT>_NEBMASK``
(the boolean mask actually used), and one extension per fit component.


sky_nebular_leak_eval.py
------------------------

Measures how much flux leaks into each nebular emission-line window in a
``DecomposeCleanSky.py`` output's residual — a direct, per-line,
per-extension test of whether SKY_EAST/SKY_WEST actually contain no
nebular-line flux, as the DRP's own sky subtraction assumes.

**Usage**::

    sky_nebular_leak_eval.py [-ext LIST] [-v VEL] [-lmc] [-smc]
                             [-sigma S] [-thresh T] [-out ROOT]
                             [-template_ext EXT] [-template_band LO,HI]
                             fits_file [fits_file ...]

**Arguments:**

fits_file
    One or more ``DecomposeCleanSky.py`` output files.

**Options:**

-ext LIST
    Comma-separated extensions to evaluate (default: ``SKY_EAST,SKY_WEST``
    — must match what ``DecomposeCleanSky.py`` was run with).

-v VEL / -lmc / -smc
    Nebular systemic velocity — must match what ``DecomposeCleanSky.py``
    was run with, or the fit window is centered on the wrong wavelength.

-sigma S
    Initial Gaussian sigma guess in Angstrom (default 1.0).

-thresh T
    ``|amplitude/amplitude_error|`` above which a line is flagged as a
    candidate leak in the printed summary (default 3.0); does not affect
    what's written to the output table.

-out ROOT
    Output table filename root (default: ``nebular_leak``).

-template_ext EXT / -template_band LO,HI
    Optional shared-systematic template correction: an extension whose
    RESID is used to correct every other requested extension's residual
    inside ``-template_band`` (default: off; band default 9000,9600 Å if
    used). Confirmed at r~0.95-0.98 between FLUX/SKY_EAST/SKY_WEST
    *within one exposure*, but only r~0.5-0.6 *across different
    exposures*, so a template never transfers across files. Useful when
    the band is dominated by a shared instrumental systematic (e.g. an
    OH-line/LSF template mismatch in the Z channel) rather than photon
    noise — an extension confirmed to carry little real signal (e.g. one
    that fails a known-fixed-ratio check) can serve as that exposure's own
    correction template for the others.

**Description:**

Fits a Gaussian plus constant background directly to the RESID array
(observed − model) at each of ``sky_gaussfit.resolve_nebular_lines()``'s
Doppler-shifted rest wavelengths — deliberately not anchored to a
model-fit shape the way ``sky_residual_eval.py``'s airglow-line fitting
is, since a nebular line has ~no flux in a ``DecomposeCleanSky.py`` model
by construction (its window was excluded from the fit): there is no
model peak to anchor on, so the line's own catalog wavelength is the
prior instead. A significant nonzero fitted amplitude at a nebular line's
position in SKY_EAST/SKY_WEST is direct evidence of nebular contamination
the mask-and-wrap decomposition — and by extension the DRP's own sky
subtraction — did not remove.

Only ever needs ``(wave, flux_observed, flux_model)`` per extension, so
it works unchanged against a future native-nebular-family SkyDecomp
output too, not just ``DecomposeCleanSky.py``'s mask-and-wrap models.

**Output:**

``<ROOT>_lines.fits`` — one row per (file, extension, line) with the
fitted flux/error/SNR/center/width and (if ``-template_ext`` was used)
both the raw and template-corrected values.

**See Also:** :doc:`api/sky_nebular_leak_eval/index`


PlotNebularLeak.py
------------------

Interactive Plotly visualization of the ``DecomposeCleanSky.py``/
``sky_nebular_leak_eval.py`` workflow for one exposure: the
shared-systematic residual pattern across FLUX/SKY_EAST/SKY_WEST in a
chosen wavelength band, plus zoomed before/after panels at each requested
nebular line showing the raw residual, the scaled template being
subtracted, the corrected residual, and the fitted Gaussian.

**Usage**::

    PlotNebularLeak.py [-lines LIST] [-targets LIST]
                       [-template_ext EXT] [-template_band LO,HI]
                       [-leak_file PATH] [-title TITLE] [-outfile PATH]
                       decomp_file

**Arguments:**

decomp_file
    A ``DecomposeCleanSky.py`` output FITS file.

**Options:**

-lines LIST
    Comma-separated ``NEBULAR_LINES`` names to show zoomed panels for
    (default: ``siii_a,siii_b``).

-targets LIST
    Comma-separated extensions to show zoom panels for (default:
    ``FLUX,SKY_WEST``).

-template_ext EXT / -template_band LO,HI
    Same convention as ``sky_nebular_leak_eval.py`` (default template
    extension: ``SKY_EAST``; band: 9000,9600). Recomputed locally with
    ``fit_template_scale``, not read from ``-leak_file``, so this plot
    works even without one.

-leak_file PATH
    ``sky_nebular_leak_eval.py`` output (``<root>_lines.fits``). If given,
    the fitted Gaussian curves are drawn from its own fitted columns
    (matching that table's numbers exactly); otherwise this script fits
    them itself with the same ``fit_leak_line`` routine.

-title TITLE / -outfile PATH
    Plot title (default: input filename) / output HTML path (default:
    ``Overview_Plot/<stem>.nebleak.html``).

**Description:**

Top panel: ``<EXT>_RESID`` for every extension present, overlaid over the
template band, sharing one y-axis — shows the shared-systematic
correlation directly (same shape, different amplitude, across
extensions within one exposure). One row of zoom panels per ``-lines``
entry, one column per ``-targets`` entry, each showing the raw residual,
scaled template, corrected residual, and fitted Gaussian.

**See Also:** :doc:`api/PlotNebularLeak/index`


SkySubNebEval.py
----------------

Science-specific evaluator for the SkySub* method family: fits
``sky_gaussfit.NEBULAR_LINES`` directly on each method's sky-subtracted
FLUX, checks fixed-ratio doublets against their known atomic-physics
value, and measures flux/ratio scatter across repeated observations of
the same tile.

**Usage**::

    SkySubNebEval.py [-v VEL] [-lmc] [-smc] [-sigma S] [-mjd_close DAYS]
                     [-snr_min S] [-out ROOT] fits_file [fits_file ...]

**Arguments:**

fits_file
    One or more SkySub*.py output FITS files (WAVE, FLUX, SKY, DRP_ALL —
    the layout every SkySubDrp/Orig/Dev1/Dev2/Dev3.py output shares).
    Each file's primary header ``TITLE``/``METHOD`` keywords label its
    rows.

**Options:**

-v VEL / -lmc / -smc
    Nebular systemic velocity, applied to EVERY row regardless of target
    if given. Default: none — each row's velocity is instead looked up
    from that row's own DRP_ALL ``Redshift`` (see the introduction at the
    top of this page). Only pass these to force one velocity
    across an entire run.

-sigma S
    Initial Gaussian sigma guess in Angstrom (default 1.0).

-mjd_close DAYS
    A tileid group's repeat exposures are tagged "closely spaced" when
    their MJD span is below this (default 7.0 days).

-snr_min S
    Minimum per-line SNR (both lines of a doublet) required to report
    that row's ratio (default 5.0) — verified necessary: an unfiltered
    median [SIII] ratio across a real 140-row sample was ~2.98-3.00,
    nowhere near the true 2.44 (most fibers in an arbitrary sample aren't
    pointed at a real emission-line target); restricting to SNR>5 on both
    lines dropped it to 2.41.

-out ROOT
    Output filename root (default: ``nebeval``).

**Description:**

Two-stage evaluation:

1. ``fit_nebular_row`` fits every ``NEBULAR_LINES`` entry (Gaussian +
   local constant background via ``sky_nebular_leak_eval.fit_leak_line``,
   with a joint double-Gaussian fit for the blended [OII] 3726/3729 pair
   via ``lvm_gaussfit.fit_double_gaussian_to_spectrum``) directly on the
   sky-*subtracted* FLUX. ``DOUBLETS`` then computes each ratio pair's
   measured value, propagated error, and deviation from the true ratio.
   Three kinds of ground truth: ``fixed`` ([SIII] 9531/9069 = 2.44, [NII]
   6583/6548 = 3.0, [OIII] 5007/4959 = 2.98 literature-only/lower
   confidence), ``bounded_above`` (Hβ/Hα ≤ 0.350 — reddening can only push
   this down, never up), and ``free`` ([SII] 6731/6716, [OII] 3729/3726 —
   density-dependent, no fixed value, useful only through repeat-scatter
   consistency). [OI] 6364/6300 is deliberately excluded — dominated by
   sky-subtraction residual (the same airglow doublet as
   sky6300/sky6363), not real nebular signal.
2. ``repeat_scatter`` groups rows by DRP_ALL ``tileid`` (excluding
   tileid 11111, a confirmed grab-bag placeholder spanning unrelated
   targets, not a genuine repeat pointing), checks each group's RA/Dec is
   tightly clustered, and computes robust (MAD) flux/ratio scatter per
   group, tagging groups CLOSE when their MJD span is below
   ``-mjd_close``.

**Output:**

``<ROOT>_lines.fits`` — one row per spectrum (FILE, ROUTINE, VARIANT,
ROW, TILEID, MJD, EXPNUM, SCI_RA, SCI_DEC, SURVEY, VEL, then per-line
FLUX/FLUX_ERR/SNR and per-doublet RATIO/RATIO_ERR/DEV_SIGMA).

``<ROOT>_repeats.fits`` — one row per (ROUTINE, VARIANT, TILEID) group
with ≥2 rows, plus a printed head-to-head table of median scatter per
doublet per method, restricted to CLOSE groups by default.

**See Also:** :doc:`api/SkySubNebEval/index`


PlotSkySubNebEval.py
--------------------

Interactive Plotly visualization of ``SkySubNebEval.py``'s per-line/
doublet fits for ONE exposure across one or more SkySub*.py output files
(methods) — the picture behind that tool's numbers, one exposure at a
time; the primary day-to-day diagnostic for "is this method doing
something sensible here."

**Usage**::

    PlotSkySubNebEval.py [-row N | -expnum N] [-v VEL] [-lmc] [-smc]
                         [-sigma S] [-title TITLE] [-outfile PATH]
                         fits_file [fits_file ...]

**Arguments:**

fits_file
    One or more SkySub*.py output files, all for the same exposure (e.g.
    the same expnum run through ``SkySubRun.py -routine`` drp/orig/dev1/
    dev2/dev3). Each file's column is labeled from its primary header
    ``TITLE``/``METHOD``.

**Options:**

-row N / -expnum N
    Row index (default 0) or DRP_ALL ``expnum`` (overrides ``-row``).

-v VEL / -lmc / -smc
    Explicit velocity override — same per-row DRP_ALL ``Redshift``
    default as ``SkySubNebEval.py`` above.

-sigma S / -title TITLE / -outfile PATH
    Initial Gaussian sigma guess (default 1.0); plot title (default: the
    exposure's expnum); output HTML path (default:
    ``Overview_Plot/nebeval_<expnum>.html``).

**Description:**

**Row 1** — a pointing/Moon/Sun geometry table: Target
(Science/SkyE/SkyW/Moon/Sun), RA, Dec., PA, angular distance from the
science field, altitude, Moon illumination (%), astrometry source, shadow
height. Same column layout as :doc:`data_quality`'s ``QualSFrame.py``
pointing table, read from the equivalent DRP_ALL columns (already
computed once by ``SummarizeCframe.py``) rather than recomputed.

**Row 2** — the full sky-subtracted spectrum, all methods overlaid, log
y-axis with a fixed range so every exposure's overview is directly
comparable at a glance. Each trace is median-filtered (11 pixels, not raw
per-pixel flux) so inter-method differences aren't swamped by per-pixel
noise; a color key and a compact exposure-ID line (expnum, tileid, MJD,
Survey, the velocity actually used) sit inside the panel itself, not the
figure margin.

**Row 3** — a summary table, one row per method, one column per
``SkySubNebEval.DOUBLETS`` entry (ordered by increasing rest wavelength),
with the true/bound value in the column header and a warning marker on
any deviant ratio.

**Rows 4+** — one row of panels per line group (OII, Hβ, [OIII] a/b,
[NII]+Hα, [SII], [SIII] a/b), one column per method, each row sharing one
y-axis range across its columns. That shared range is bounded by the
2nd/98th percentile of observed flux pooled across all columns (not
literal min/max) plus the fitted curves' true min/max — a single noisy
method's occasional extreme pixel otherwise sets a shared range wide
enough to flatten every *other*, well-behaved method's panel into a
near-flat line, which reads as "this method's fit is bad" when it's
really just the shared axis being dominated by a different column's
outlier.

**See Also:** :doc:`api/PlotSkySubNebEval/index`


PlotSkySubNebRun.py
-------------------

Run-level companion to ``PlotSkySubNebEval.py``: aggregates
``SkySubNebEval.py``'s per-row nebular-line fits and per-tileid
``repeat_scatter`` groups across a whole run (one or more SkySub*.py
output files, each with many rows — not just one exposure) into a single
HTML report, so a method comparison doesn't require paging through one
file per exposure. This is the "distributions" tool the single-exposure
diagnostic above was always meant to be followed by.

**Usage**::

    PlotSkySubNebRun.py [-lines_file PATH] [-v VEL] [-lmc] [-smc]
                        [-sigma S] [-snr_min S] [-mjd_close DAYS]
                        [-mask PATH] [-title TITLE] [-outfile PATH]
                        [fits_file ...]

**Arguments:**

fits_file
    One or more SkySub*.py output files, the same contract
    ``SkySubNebEval.py``'s own CLI takes. Ignored if ``-lines_file`` is
    given. Required (not ignorable) for the Continuum Residual section
    below, which needs the raw WAVE/FLUX/SKY arrays — that section is
    silently skipped in ``-lines_file`` mode.

**Options:**

-lines_file PATH
    An existing ``SkySubNebEval.py`` run's ``<root>_lines.fits`` output —
    skips the (currently unparallelized) per-row Gaussian refit, useful
    for iterating on the emission-line sections against a large sample
    without redoing it every time. Exactly one of ``fits_file``/
    ``-lines_file`` is required. ``repeat_scatter``'s own grouping/MAD
    summary is always recomputed from whichever row table is in hand
    either way — that step is cheap (aggregation only, no curve
    fitting), so there is no separate ``-repeats_file`` fast path. The
    Continuum Residual section needs the raw spectra, so it is skipped
    entirely in this mode (a printed note says so).

-v VEL / -lmc / -smc
    Nebular systemic velocity override, same precedence/default (per-row
    DRP_ALL ``Redshift`` lookup) as ``SkySubNebEval.py``. Ignored if
    ``-lines_file`` is given.

-sigma S / -snr_min S / -mjd_close DAYS
    Same meaning as the equivalent ``SkySubNebEval.py`` options.

-mask PATH
    Clean-pixel mask (WAVE/MASK extensions, a ``palace_make_mask.py``-
    style file) for the Continuum Residual section — default
    ``data/sky_mask.fits``, the project-wide standard also used by
    ``SkySubOrig.py`` and ``sky_residual_eval.py``.

-title TITLE / -outfile PATH
    Report title (default: ``SkySubNebRun``) / output HTML path
    (default: ``nebrun.html`` — rerun with the same path to update it in
    place rather than accumulating one file per attempt).

**Description:**

The report is a real HTML document — actual ``<h1>``/``<h2>``/``<h3>``
headings with CSS margins around several small, focused Plotly figures
(``build_figures()``/``build_continuum_figures()``), **not** one giant
multi-row Plotly canvas with hand-tuned pixel margins standing in for
section breaks. That approach was tried first and kept needing another
manually-tuned margin/spacer-row fix every time a new section was added;
normal HTML block flow reserves space between sections automatically and
cannot overlap, which a Plotly-internal annotation used as a section
divider cannot guarantee.

Every scatter section uses the same grid convention as
``PlotSkySubNebEval.py``: one row per metric, one **column per method**
(not all methods overlaid in one panel) — overlaying every method's
points in a single panel is fine at the tens-of-points scale of a small
test set, but stops being legible at real survey scale (hundreds to
thousands of rows), whereas a single-method panel stays readable
regardless of how large the run is. Each row's method-name labels are
shown once, above the first row of the whole grid, rather than repeated
above every row, and the vertical gap between rows is held to a fixed
~90 pixels regardless of how many rows a given figure has — both a fixed
fraction of figure height and repeated per-row labels looked fine on a
tall (6-row) grid but visibly collided on a short (3-row) one::

    <h1>title</h1>
    [run-summary table: rows fit, SNR-pass count/fraction per doublet]
    <h2>Line Ratios</h2>
    <h3>Ratio vs. Line Flux</h3>
    [one row per doublet, one column per method; x = total flux of its
     two lines (log), y = the ratio, one point per SNR-passing row,
     pooled across the whole run]
    <h3>Repeat-Group Scatter (MAD) vs. Median Flux</h3>
    [same grid; x = a repeat group's median total line flux (log,
     the SAME x-range _ratio_vs_flux_figure computed from the full
     sample, forced rather than recomputed from this panel's own much
     smaller repeat-group sample), y = that group's ratio MAD, one
     point per repeat group, CLOSE groups only]
    <h2>Repeat-Observation Flux Consistency</h2>
    <h3>Repeat-Group Fractional Scatter vs. Median Flux</h3>
    [one row per FLUX_METRICS entry (total flux for OII, OIII_b, NII_b,
     SII_a+SII_b, SIII_b — chosen to use only lines whose partner is a
     fixed multiple, or sum both when neither dominates), one column per
     method; x = the SAME shared total-line-flux x as above, y = that
     group's fractional MAD (robust MAD / median) of the flux total
     across the group's repeat exposures]
    <h3>Summary</h3>
    [interactive table: median-across-groups absolute MAD for ratios,
     median-across-groups fractional MAD for flux totals]
    <h2>Continuum Residual (B/R/Z)</h2>
    <h3>Post-Subtraction Continuum Level (All Exposures)</h3>
    [one row per spectrograph arm, one column per method; EVERY
     exposure, not just repeat groups; x = that exposure's own
     pre-subtraction continuum median (log), y = its post-subtraction
     continuum median (linear, unclipped) with a dashed red line at
     zero -- a physical continuum flux cannot go negative]
    [percent-negative summary table]
    <h3>Post-Subtraction Continuum Level Consistency</h3>
    [repeat groups only; x = a repeat group's PRE-subtraction continuum
     median in that arm (shared across methods), y = that group's
     fractional MAD of the POST-subtraction continuum median, divided
     by the group's PRE-subtraction brightness (log y)]
    <h3>Post-Subtraction Continuum RMS Consistency</h3>
    [same grid; y = fractional MAD of the POST-subtraction continuum
     RMS/NMAD, divided by its own median instead]
    <h3>Summary</h3>
    [interactive fractional-MAD table, level and RMS, one column per arm]

Both nebular-line scatter designs replace an earlier pooled-box-plot
version: pooling every row (or every repeat group) into one box per
method hid the dependence of ratio/flux scatter on how bright the line
actually was in a given fiber — plotting against flux directly controls
for that confound instead of comparing methods across an uncontrolled
mix of bright and faint measurements.

Axis/table labels are kept short: a flux total is labeled with the bare
line name (``OII``, not ``OII_FLUX``); a ratio is labeled with its actual
wavelength pair (e.g. ``OII 3730:3726``, ``SII 6716:6731``) computed live
from ``sky_gaussfit.NEBULAR_LINES`` rather than hand-typed, so the label
stays correct if a doublet's numerator/denominator convention is ever
changed (as ``SII`` was — see ``SkySubNebEval.DOUBLETS``'s own comment).
Low-density-limit reference lines (dashed, not a truth/ceiling — density-
sensitive ratios can and do vary; the low-density limit is just the value
observed at most average sky positions) are drawn wherever
``SkySubNebEval.DOUBLETS`` gives one for a ``'free'`` entry (OII ≈ 1.42,
SII ≈ 1.5), the same visual treatment as the ``'fixed'``/``'bounded_above'``
truth lines for OIII/NII/SIII/Hβ:Hα. Every method's x-axis brightness
proxy for a given line is the SAME total-doublet-flux value everywhere
that line appears (ratio-vs-flux, repeat-group MAD, and flux-consistency
panels alike) — an earlier version mixed a mean-of-two-members proxy in
the ratio panels with a different per-metric total in the flux panel,
which put the same line at a visibly different x-position/scale
depending which panel you looked at.

The **Continuum Residual (B/R/Z)** section is a different, independent
test: not a per-line measurement but the leftover *continuum* level in
each spectrograph arm (``GetSkyCont.ARM_EVAL_RANGES`` — B: 3650–5775 Å,
R: 5800–7520 Å, Z: 7570–9600 Å, with arm-overlap zones and outer edges
already excluded — the same canonical definition ``SkySubOrig.py`` uses,
distinct from two other, unrelated arm-range definitions elsewhere in the
repo used for different purposes). One function,
``GetSkyCont.arm_continuum_stats``, is run directly on raw flux with no
local continuum fit — ``data/sky_mask.fits`` (built from an actual
PALACE sky-emission model, not a crude window list) already excludes
every sky-line-dominated pixel, so the median of what remains in an arm
already is a continuum-brightness estimate and its NMAD already is a
scatter estimate. The SAME call is used both on the *pre*-subtraction
CFrame flux (reconstructed as ``FLUX+SKY``, identical across all 5
methods since they share one input CFrame — SkySub*.py's own convention
is that ``FLUX`` is already sky-subtracted and ``SKY`` is the model
subtracted from it) and on each method's own *post*-subtraction
``FLUX``, so the two are directly comparable — deliberately NOT the
pre-subtraction ``SCI_MED_B``/etc. columns some of the SkySub*.py scripts
already write, which are a continuum-FIT-quality diagnostic (raw flux
minus a locally-fit polynomial), not a brightness measurement, and are
never computed post-subtraction by any of the 5 methods.

The first subsection (**Post-Subtraction Continuum Level (All
Exposures)**) needs no repeat observation at all: a physical continuum
cannot be negative, so a per-exposure check for that requires nothing
more than one exposure. It therefore runs on every exposure in the run,
not just the repeat-tileid groups the other two continuum subsections
are restricted to (via ``SkySubNebEval.group_repeat_exposures`` —
``repeat_scatter``'s own tileid-exclusion/position-clustering/CLOSE-flag
grouping logic, factored out so this differently-shaped per-arm table
can reuse it without needing NEBULAR_LINES/DOUBLETS columns it doesn't
have).

The two repeat-only continuum subsections normalize their fractional MAD
differently, and deliberately so. The RMS subsection divides by the
group's own median NMAD — safe, since a noise floor is never near zero.
The LEVEL subsection instead divides by the group's *pre*-subtraction
brightness (the same value on the x-axis), not by its own post-
subtraction median: a post-subtraction continuum level is *expected* to
sit near zero for a good method, so dividing by its own magnitude is
unstable and can invert the ranking entirely — a near-perfect, near-zero
group can score far worse than a badly, but consistently, biased one
purely because of a near-zero denominator (confirmed on real per-group
numbers: the single worst-looking point under the naive metric had a
median residual of order 1e-17 — essentially perfect — while the single
best-looking point had the largest systematic residual in the whole
set). The level subsection's y-axis is log-scale for this reason too —
once normalized against a stable denominator, the fractional values
genuinely span about two decades, worth a log axis the same way every
flux axis elsewhere in this report is.

Column order (which method is which column) is the order the input
files/rows were actually given/created in (``_ordered_routines`` — first
occurrence in the row table), not an alphabetical sort, and is identical
across every section by construction (all figures are built from the
same ``_common_setup()`` result; the continuum section reuses that same
routine/label/color ordering rather than recomputing it).

**Output:** one HTML file, never one file per exposure — the
per-exposure picture is ``PlotSkySubNebEval.py``'s job; this one is the
aggregate/distribution view across an entire run.

**See Also:** :doc:`api/PlotSkySubNebRun/index`


Typical Workflow
----------------

1. Run all five methods through the dispatcher::

       SkySubRun.py -routine orig  XCframe_file.fits
       SkySubRun.py -routine drp   XCframe_file.fits
       SkySubRun.py -routine dev1  XCframe_file.fits -mask sky_mask.fits
       SkySubRun.py -routine dev2  XCframe_file.fits
       SkySubRun.py -routine dev3  XCframe_file.fits

2. Look at one exposure across all five methods side by side::

       PlotSkySubNebEval.py -expnum 12226 \
           sky_runs/orig/XCframe_file_orig_farlines_nearcont.fits \
           sky_runs/drp/XCframe_file_drp_farlines_nearcont.fits \
           sky_runs/dev1/XCframe_file_dev1_farlines_nearcont.fits \
           sky_runs/dev2/XCframe_file_dev2_scilines_nearcont.fits \
           sky_runs/dev3/XCframe_file_dev3_farlines_nearcont.fits

3. For a batch view across many exposures/repeat groups instead of one
   exposure, run ``PlotSkySubNebRun.py`` on the same file list (or
   ``-lines_file`` an existing ``SkySubNebEval.py`` run's own
   ``_lines.fits`` to skip refitting)::

       PlotSkySubNebRun.py -outfile nebrun.html \
           sky_runs/orig/XCframe_file_orig_farlines_nearcont.fits \
           sky_runs/drp/XCframe_file_drp_farlines_nearcont.fits \
           sky_runs/dev1/XCframe_file_dev1_farlines_nearcont.fits \
           sky_runs/dev2/XCframe_file_dev2_scilines_nearcont.fits \
           sky_runs/dev3/XCframe_file_dev3_farlines_nearcont.fits

   or run ``SkySubNebEval.py`` directly on the same file list and inspect
   its printed doublet-scatter comparison (or the ``_repeats.fits`` table)
   for the equivalent numbers without the report.

4. To check whether SKY_EAST/SKY_WEST leak nebular flux in the first
   place (an assumption every sky-telescope-based method depends on)::

       DecomposeCleanSky.py -ext SKY_EAST,SKY_WEST XCframe_file.fits
       sky_nebular_leak_eval.py CleanSky_<expnum>.fits
       PlotNebularLeak.py CleanSky_<expnum>.fits


See Also
--------

- :doc:`sky_method_eval` - generic sky-residual evaluation
- :doc:`sky_methods` - the methods being evaluated
- :doc:`api/sky_nebular_leak_eval/index` - API documentation
- :doc:`api/PlotNebularLeak/index` - API documentation
- :doc:`api/SkySubNebEval/index` - API documentation
- :doc:`api/PlotSkySubNebEval/index` - API documentation
- :doc:`api/PlotSkySubNebRun/index` - API documentation
- :doc:`spectral_fitting_local` - ``sky_gaussfit.py``'s NEBULAR_LINES/
  resolve_nebular_lines, shared by the tools on this page
- :doc:`data_quality` - ``QualSFrame.py``'s pointing table, reused by
  ``PlotSkySubNebEval.py``
