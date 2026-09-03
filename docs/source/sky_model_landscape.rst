Sky Model Comparison Landscape
===============================

This page is a map, not a manual: it says *where* each piece of the
sky-model-comparison and LSF-investigation effort lives, *how* the pieces
connect into end-to-end workflows, and what state each piece is in. For
the detailed usage of an individual committed script, follow the links
out to :doc:`sky_subtraction` or :doc:`spectral_fitting_local`, which
already document the pieces that live in this repository in full.

.. note::
   Unlike the rest of this documentation, this page is **hand-maintained
   prose, not autodoc** — Sphinx/AutoAPI only scans this repository's own
   ``py_progs/``, so it has no way to notice changes in the other two
   locations described below. Treat this page as a snapshot, most
   recently updated 2026-09-03, and expect it to drift; it is deliberately
   kept to an index of *what* and *where*, not a copy of details that
   live (and change) elsewhere.


Three Places This Work Lives
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 15 63

   * - Location
     - Git status
     - What's there
   * - ``lvm_ksl/py_progs/``
     - Tracked, Sphinx-documented
     - The production/analysis tier: ``EsoSkyObs.py``/``PalaceObs.py``
       (fetch physical sky models), the ``SkySub*``/``SkySep*`` family
       (subtraction methods, see :doc:`sky_subtraction`), the ESO/PALACE
       continuum-comparison toolchain (``SkyObsESOCompare.py`` etc.,
       see below), ``sky_residual_eval.py`` (the method-agnostic
       evaluator both candidates below use), and a vendored, older
       snapshot of PALACE's ``SkyDecomp`` (``py_progs/sky_decomp/fit.py``
       + ``data/palace_ref/``, commit ``6f877c0``).
   * - ``lvm_ksl/py_dev/``
     - Tracked, **not** Sphinx-documented (outside ``autoapi_dirs``)
     - Second tier: research code that's already parameter-driven and
       reusable across datasets, but not yet stable/documented enough
       for ``py_progs``. First contents (2026-09-03): the MLP
       prediction/evaluation/test-set-curation scripts promoted from
       ``niv/`` — ``PredictSky.py``, ``BatchPredictSky.py``,
       ``BatchPredictSkyESO.py``, ``EvalFluxResiduals.py``,
       ``MasterResidualByMoon.py``, ``SelectXCF.py`` — copied, not
       moved, so the ``niv/`` originals still exist too.
   * - ``lvmsky/skysub/``
     - Separate git repo, not controlled by this documentation
     - Ivan Katkov's ongoing ``SkyDecomp`` research (``sky_decomp/``) and
       the ``mlp_predictor`` package that trains/serializes the
       ``mlp_ensemble_split_zodi`` model. Has grown past what was
       vendored into ``py_progs`` — see `The SkyDecomp Fork`_ below.
   * - ``~/Projects/lvm_sky2609/niv/``
     - **Not a git repo.** Working directory only.
     - The MLP candidate's training pipeline (still here, not promoted —
       see `Known Gaps and Promotion Candidates`_) plus
       ``EvalCoefResiduals.py``/``PlotPredictSky.py`` and the (now also
       copied to ``py_dev/``) prediction/evaluation scripts. Several
       training-pipeline scripts hardcode local paths and are one-off
       drivers, not reusable CLI tools.

Two informal, untracked notes at the ``lvm_ksl`` repo root go one level
deeper than this page on two specific sub-threads: ``ivan.md``
(``SkyDecomp``'s internal LSF handling) and ``mlp_sky_predictor_notes.md``
(the MLP retraining/evaluation effort). They are deliberately left out of
git as fast-moving scratch notes — see `Known Gaps and Promotion
Candidates`_.


Current Workflows
------------------

There are two independent "candidates" for predicting the sky spectrum
at a given time/position, evaluated against real LVM data by the same
downstream tooling. Both ultimately answer the question "how well does
this model match reality," but they get there by different routes.

ESO Sky Model Comparison Workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    EsoSkyObs.py (fetch model)
        -> SkyObsESOCompare.py (per-exposure continuum/line fit vs. real data)
        -> SkyObsESO_analysis.py / SkyObsESOStack.py (aggregate, stack, plot)
    SkyObsPalaceCompare.py -- same interface, PALACE lines + ESO moon/zodi

All four scripts are committed, in ``py_progs/``, with API stubs, but
**have no narrative documentation page yet** (see `Known Gaps and
Promotion Candidates`_) — this
is the toolchain that established the real, still-unresolved findings:
a genuine blue continuum excess at low airmass (not an arm-splice
artifact, checked directly), and an IR continuum overprediction traced
to the clean-pixel mask under-accounting for blended OH-line wings in
the Z arm, not a fitting-method problem (confirmed with polynomial vs.
matched-DOF B-spline continuum reconstructions).

MLP Ensemble Prediction Workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Test-set curation (routine — the usual entry point for testing the
*existing* trained checkpoint against fresh data, as opposed to
retraining)::

    SummarizeCframe.py -by fiber (source XCframe summary, if one doesn't
                                    already exist)
        -> SelectXCF.py [py_dev] (criteria-filtered, moon-geometry-stratified subset)

``SelectXCF.py`` does not check its selection against the training
corpus's own row list, so it does not by itself guarantee a genuinely
held-out test set — worth a manual cross-check if that matters for what
you're testing.

Training (occasional, only needed for a *new* checkpoint)::

    SelectXCF.py [py_dev] (criteria-filtered exposure subset)
        -> ConvertForDecompose.py (reformat for lvmsky's decompose_parallel.py)
        -> [lvmsky: decompose_parallel.py --fit-model lsf-surface-iterative-split-zodi]
        -> stage1_wavecache.py -> stage1_train.py -> stage2_train.py
        -> mlp_ensemble_stage2_production.pt (checkpoint)

Prediction and evaluation (the part that runs against an existing
checkpoint; scripts marked ``[py_dev]`` are promoted and tracked, the
rest are still ``niv/``-only)::

    PredictSky.py [py_dev] (single exposure)  --or--  BatchPredictSky.py [py_dev] (parallel batch)
        -> EvalFluxResiduals.py [py_dev] (method-agnostic, via sky_residual_eval.py)
        -> EvalCoefResiduals.py (model-specific, coefficient-space)
        -> MasterResidualByMoon.py [py_dev] (moon-brightness-stacked residual spectra)
        -> PlotPredictSky.py (per-exposure interactive diagnostic)

``BatchPredictSkyESO.py`` [py_dev] runs the *ESO* candidate through this
same batch/evaluation interface (``EsoSkyObs.run_sky_obs`` in place of
the MLP), so the two candidates land in a directly comparable
``WAVE/FLUX_OBS/FLUX_PRED/LINE_PRED/META`` layout and can be evaluated
by the same downstream tools — the multi-candidate harness design is
confirmed working (ESO continuum found ~25-35% too bright, worse than
the MLP, as expected for a model with no data-driven correction).

The training pipeline above and the two still-``niv/``-only evaluation
scripts remain outside git — see `Known Gaps and Promotion Candidates`_.


How To: Test the MLP Candidate Against New Data
--------------------------------------------------

Concrete recipe for the routine case — you have an existing trained
checkpoint and want to see how well it does against a fresh set of
observations, not retrain anything. All of ``py_dev/`` is written to
take its inputs as CLI arguments (see `Known Gaps and Promotion
Candidates`_ for why that mattered enough to fix), so this is meant to
run unmodified from a pristine directory.

**Prerequisites** (one-time shell setup, not per-run):

- ``conda activate niv`` — the ``niv`` conda env has everything the
  chain needs (torch, astropy, scipy, plotly, clarabel); the project's
  usual ``lvmdrp26`` env does **not** have torch installed, so
  ``PredictSky.py``/``BatchPredictSky.py`` fail there at import time.
- Both ``py_progs/`` and ``py_dev/`` on ``PATH`` and ``PYTHONPATH``
  (``py_dev`` is not on either by default — add it yourself the same
  way ``py_progs`` is added, per the Setup section at the top of this
  documentation).

**The four steps**, from a pristine directory, given an existing source
XCframe summary FITS (``source.fits`` below) and the production
checkpoint::

    mkdir my_test_run && cd my_test_run

    # 1. Curate a test set
    SelectXCF.py source.fits my_test_set.fits \
        --n 200 --seed 42 --exptime 900 --fluxcal MOD --min-glat 10

    # 2. Run the MLP candidate against it
    BatchPredictSky.py my_test_set.fits \
        --model /path/to/mlp_ensemble_stage2_production.pt \
        --n 200 --n-workers 8 --outfile my_test_set_batch.fits

    # 3. Quantify accuracy
    EvalFluxResiduals.py my_test_set_batch.fits \
        -out my_test_set_eval -plotdir plots_my_test_set

    # 4. Produce the comparison plots
    MasterResidualByMoon.py my_test_set_batch.fits \
        --outfile my_test_set_master_resid.fits --html my_test_set_master_resid.html

Step 2 is the expensive one (real per-row PALACE decomposition + LSF
reconstruction per exposure; minutes, not seconds, for a few hundred
rows even with 8 workers) — worth backgrounding for anything past a
handful of test rows.

**Outputs**: ``my_test_set.fits`` (the curated subset), ``my_test_set_batch.fits``
(``WAVE``/``FLUX_OBS``/``FLUX_PRED``/``META``), the eval summary table +
``plots_my_test_set/`` figures, and ``my_test_set_master_resid.fits``/``.html``
(the dark/medium/bright moon-tercile comparison — this is the "relation
plots" step).

**To add the ESO candidate to the same comparison**, run
``BatchPredictSkyESO.py my_test_set.fits --outfile my_test_set_batch_eso.fits``
against the identical test set, then steps 3-4 unchanged on that file —
this is exactly the multi-candidate harness confirmed working under
`MLP Ensemble Prediction Workflow`_ above. As of 2026-09-03 this step
also applies each row's real LSF (from ``my_test_set.fits``'s own
``LSF`` extension) to the ESO prediction, matching the convolution the
MLP candidate already gets — see `Current Open Problem: The
Instrumental LSF`_ below.

**Caveat repeated from above**: ``SelectXCF.py`` draws from whatever
rows pass its hard cuts without checking against any particular training
corpus's row list, so step 1 does not by itself guarantee the result is
held-out relative to a specific checkpoint's training data. Confirmed
concretely on this run: the 200-row selection came from the identical
3,501-row filtered pool the training corpus was itself built from.

**Verified end-to-end 2026-09-03** on a real 200-row test set (source:
``XCframe_1.2.1_7325_48860_1_10_fiber.fits``, the same base corpus the
production checkpoint was trained from): 200/200 rows predicted
successfully in 132.7s (0.66s/row, 8 workers); moon-brightness terciles
split 67/66/67. Headline continuum/line accuracy from
``EvalFluxResiduals.py``'s summary table: per-arm continuum
``CONT_FIT_QUALITY`` medians 1.44 (B) / 1.52 (R) / 2.36 (Z) (units of
noise-normalized residual — Z arm reads worst, consistent with the
Z-band continuum problem already flagged as unresolved during the
original retraining/evaluation effort); per-arm line-amplitude ratio
(predicted/observed) medians 0.87 (B) / 0.94 (R) / 0.92 (Z), consistent
with the ~6-12% line under-prediction already known from the original
training-time evaluation.


The SkyDecomp Fork
--------------------

Both workflows ultimately lean on the same underlying decomposition
machinery (PALACE line templates + a Moon/Zodi continuum model fit by
constrained least squares), but two different vintages of it are in play
at once, in two different places:

- ``py_progs/sky_decomp/fit.py`` — vendored once (commit ``6f877c0``),
  self-contained, no ``lvmdrp``/``lvmsky`` dependency. Used directly by
  ``XSkySepIvan.py``/``SkySubDev2.py``. LSF handling here is the
  two-layer design in ``ivan.md``: a fixed Gaussian convolution at
  construction, optionally refined (``-refits N``) into a free,
  non-negative 11-tap kernel **per spectrograph arm** (B/R/Z — not
  continuous in wavelength).
- ``lvmsky/skysub/sky_decomp/lsf_surface_iterative.py`` — a newer,
  **not vendored**, more general design: the kernel taps are
  B-spline functions of wavelength (smoothly varying within an arm, not
  one fixed kernel per arm), refined jointly with the continuum over
  several iterations. ``py_dev/PredictSky.py`` imports this directly
  from ``lvmsky`` (``sys.path.insert('/Users/long/SDSS/lvmsky/skysub')``)
  to reconstruct the MLP's predicted flux using each exposure's own
  per-row LSF.

So a change to ``lvmsky``'s LSF-surface code reaches the MLP prediction
workflow immediately (live import) but has no effect on
``XSkySepIvan.py``/``SkySubDev2.py`` unless and until it's re-vendored.


Current Open Problem: The Instrumental LSF
--------------------------------------------

This is the live thread motivating this page. Status as of 2026-09-03:

- The DRP writes a real per-row, per-wavelength LSF into CFrame/XCframe
  files (the ``LSF`` extension). ``lvm_line_profile.py`` (see
  :doc:`spectral_fitting_local`) independently fit 18 raw airglow lines
  and found the fitted Gaussian FWHM consistently *wider* than this
  header value, by an amount that is **not smooth in wavelength**
  (1.4-12%, line-specific) — evidence pointing at several catalog lines
  being unresolved blends rather than a genuine, smooth LSF calibration
  error. This has not yet been re-tested per-exposure (only one pooled,
  280-row sample so far), so it's not yet known whether any real residual
  gap is a *stable* correction (safe to apply as one curve on top of the
  header value) or varies exposure-to-exposure.
- ``SkySubDev2.py`` tested a flat multiplicative LSF-width correction
  (``-lsf_boost``) inside the full PALACE pipeline and found it makes
  fit residuals monotonically *worse* with more boost — consistent with
  the header LSF being close to right at the single-component level,
  reinforcing the blend-not-calibration-error read above.
- The **ESO Sky Model's own internal LSF is a fixed, LVM-untuned
  constant** — found 2026-09-03 that ``calcskymodel`` (the local engine)
  convolves with a fixed 0.8-pixel (~0.4 A at the reference wavelength)
  Gaussian, identical for every exposure, with none of the real
  exposure-to-exposure variation the DRP's own per-row LSF carries.

**Fixed 2026-09-03**: ``BatchPredictSkyESO.py`` now applies each row's
real per-row LSF (from the same file's ``LSF`` extension) on top of
ESO's own internal ~0.4 A kernel, using the identical construction
``PredictSky.py`` uses for the MLP candidate
(``sky_decomp.lsf_surface_iterative.build_lsf_operator``) — so the two
candidates are now compared on equal footing. This also required
changing ``BatchPredictSkyESO.py``'s interface from a
``(meta_file, batch_file)`` pair to a single ``fits_file`` argument
(the same XCframe-layout file ``BatchPredictSky.py`` already takes
directly): the old ``meta_file`` (a ``*_meta_only.fits`` artifact of the
training pipeline) never carried spectral extensions at all, and
``batch_file`` (``BatchPredictSky.py``'s own output) drops the ``LSF``
extension that the source XCframe file has — so neither of the old
inputs could have supplied it. Verified numerically on real data
(row with LSF FWHM ~1.57 A): unconvolved ESO peak at 5577 A
4.15e-12, convolved peak 1.85e-12 — broader and lower, as physically
expected, with integrated flux over the line roughly conserved.

Planned next step: re-run ``lvm_line_profile.py``-style fits across
multiple *separate* exposures to determine whether the DRP's own
header-vs-fit LSF gap (first bullet above) is stable or exposure
dependent, which decides whether a single global correction curve on
top of the header LSF is enough or a per-exposure one is needed. That
question is now more directly testable than before, since the ESO
candidate's comparison is no longer confounded by its own missing
convolution.


Known Gaps and Promotion Candidates
--------------------------------------

**Documentation gaps** (code exists and is committed, no narrative page):

- ``SkyObsESOCompare.py`` / ``SkyObsPalaceCompare.py`` /
  ``SkyObsESOStack.py`` / ``SkyObsESO_analysis.py`` — the whole ESO/
  PALACE continuum-comparison toolchain, despite real findings (blue
  continuum excess, IR overprediction) resting on it.
- ``lvm_skyfit.py`` — has an API stub but no narrative section in
  :doc:`spectral_fitting_local` yet.

**Recently promoted** (2026-09-03): ``PredictSky.py``,
``BatchPredictSky.py``, ``BatchPredictSkyESO.py``, ``EvalFluxResiduals.py``,
``MasterResidualByMoon.py``, ``SelectXCF.py`` copied from ``niv/``
(untracked) to ``py_dev/`` (tracked, not yet Sphinx-documented) — the
first test of the ``py_dev`` tier described in the table above (``niv/``
originals left in place, not deleted). Along the way: a placeholder
``DEFAULT_MODEL`` path in ``PredictSky.py``/``BatchPredictSky.py``
(``~/foo/goo/...``, never filled in) was found and, on reflection, not
just fixed but removed — ``--model`` is now a required argument with no
default at all, since this tier is meant to work against different
trained models for different purposes and a silent fallback (even a
correct one) risks going stale the moment that stops being true.
``EvalFluxResiduals.py``/``MasterResidualByMoon.py``/
``BatchPredictSkyESO.py``'s ``py_progs``-location constants were also
made CLI-overridable for the same reason, matching the pattern
``PredictSky.py`` already used for ``--lvmsky-skysub``.

**Remaining promotion candidates**:

- ``EvalCoefResiduals.py`` / ``PlotPredictSky.py`` — still ``niv/``-only;
  reasonably parameter-driven already (take a file path as their main
  argument) and could follow the same path once there's a reason to.
- ``lvmsky/skysub/sky_decomp/lsf_surface_iterative.py`` — candidate for
  a second, updated vendoring pass (into ``py_progs``, following PALACE's
  own precedent — see `The SkyDecomp Fork`_) if the LSF investigation
  above concludes the smoother kernel is worth the extra complexity.

**Not promotion candidates** (tied to one specific exercise, not
reusable tools): the rest of ``niv``'s training pipeline
(``ConvertForDecompose.py``, ``stage1_wavecache.py``, ``stage1_train.py``,
``stage2_train.py``) — ``SelectXCF.py`` itself turned out to be general
enough to promote (see above): it's used for curating a fresh test set
just as much as for building a training corpus, so it wasn't really
training-pipeline-specific after all.


See Also
--------

- :doc:`sky_subtraction` - Full documentation of every committed sky
  subtraction/modeling script, including the ``SkyDecomp``/PALACE
  LSF discussion under ``SkySubDev2.py``
- :doc:`spectral_fitting_local` - ``lvm_line_profile.py``'s Gaussian-vs-Moffat
  airglow line-shape investigation
- :doc:`rss_combining` - The two multi-exposure combination strategies
  (unrelated to this page's model-comparison focus, but shares
  ``gauss_combine.py``/``lvm_skyfit.py`` as the newest additions to
  ``py_progs/``)
