Physically Based Sky Models
===========================

The DRP, and most of the tools in :doc:`sky_subtraction`, estimate the
sky in the science field from the two sky telescopes (SkyE and SkyW)
observed at the same time. This page is about a different idea: predict
the sky spectrum from a *model* that starts from the physics of the
night sky -- airglow emission from the upper atmosphere, sunlight
scattered by the Moon and by interplanetary dust, and the atmosphere
those photons pass through on the way to the telescope.

Three such approaches are being explored for LVM:

1. **The ESO Sky Model** -- a complete physical model of the sky at
   Cerro Paranal, run here for the LVM site and observing geometry.
2. **PALACE** -- a newer, detailed physical model of the airglow alone
   (emission lines and airglow continuum), also built for Paranal.
3. **The semi-empirical machine-learning (ML) approach** -- PALACE's
   line lists and physical moon/zodiacal-light shapes are used as a
   template set, the templates are fitted to real LVM sky spectra, and
   a neural network learns how those fitted amplitudes in the sky
   telescopes map onto the sky in the science field.

They differ mainly in how much they rely on physics versus LVM's own
data: the ESO Sky Model and PALACE predict the sky from the observing
conditions alone, while the semi-empirical approach uses physics only
for the *shapes* of the sky components and learns their *strengths* from
LVM observations.

.. note::
   This page is hand-maintained prose. Much of the code it describes
   lives in ``py_dev/`` (not covered by the API documentation) or in the
   separate ``lvmsky`` repository, so Sphinx cannot notice when it
   changes. Treat it as a snapshot, most recently updated 2026-10-07.


The Three Approaches
--------------------

The ESO Sky Model
^^^^^^^^^^^^^^^^^

**What it is.** The Cerro Paranal Advanced Sky Model developed for ESO
(Noll et al. 2012), available as a local program (``calcskymodel``) and
as the SkyCalc web service. Given a pointing, a time, and the solar
radio flux, it computes every major contribution to the night-sky
spectrum, including scattered moonlight, zodiacal light, airglow
emission lines, the airglow continuum, and thermal emission, together
with the atmospheric transmission.

**What it is calibrated on.** Paranal observations. Nothing in it is
fitted to LVM data, so any difference between it and an LVM sky
spectrum is a real test of the model at Las Campanas. ``EsoSkyObs.py``
can run it with LCO or Paranal site parameters (altitude, pressure).

**Tools in lvm_ksl** (all in ``py_progs/``, documented in
:doc:`sky_model_tools` and :doc:`sky_methods`):

- ``EsoSkyObs.py`` -- predict the sky for one pointing and time, with
  the MOON, ZODI, DIFFUSE and LINES components kept separate.
- ``SkySepESO.py`` -- sky subtraction of an XCframe summary file using
  the ESO model's MOON, ZODI and DIFFUSE components as the continuum.
- ``SkyObsESOCompare.py``, ``SkyObsESOStack.py``,
  ``SkyObsESO_analysis.py`` -- compare the model's continuum and lines
  with observed sky spectra, exposure by exposure and stacked by moon
  altitude and airmass.
- ``BatchPredictSkyESO.py`` [``py_dev/``] -- run the ESO model through
  the same prediction/evaluation harness as the semi-empirical model
  (see `Comparing the Approaches`_), including convolution with each
  exposure's own LSF.

**What we have found.** The ESO continuum is about 25-35% too bright
compared with LVM sky spectra. There is a genuine excess of observed
blue continuum at low airmass (checked not to be an artifact of joining
the spectrograph arms). There is also an apparent over-prediction of
the infrared continuum; checks so far point to the clean-pixel mask not
fully accounting for the wings of blended OH lines in the Z arm, rather
than to the fitting method, but the problem is not resolved.


PALACE
^^^^^^

**What it is.** PALACE, the "Paranal Airglow Line And Continuum
Emission" model (Noll et al. 2025), is a physical model of the airglow
alone. It predicts thousands of individual emission lines -- OH, O2,
and atomic lines of Na, K, O, N and H -- and the airglow continuum (from
HO2, FeO and O2), as a function of time, solar activity and observing
geometry. It does *not* include moonlight or zodiacal light.

**What it is calibrated on.** Paranal observations, like the ESO Sky
Model. PALACE is an external Python package (see ``Readme.md`` for
installation).

PALACE is used in two quite different ways:

1. **As a forward model**: predict the airglow for a pointing and time
   and compare with the data.

   - ``PalaceObs.py`` -- PALACE prediction for one pointing and time,
     split by species, in the same output convention as
     ``EsoSkyObs.py``.
   - ``SkyObsPalaceCompare.py`` -- the same comparison as
     ``SkyObsESOCompare.py``, using PALACE lines plus the ESO model's
     moon and zodiacal light.
   - ``SkySepPalace.py`` -- sky subtraction with one free amplitude per
     PALACE species plus a moon amplitude.

2. **As a template set** fitted to an observed spectrum. PALACE's line
   lists fix *where* each sky line falls and the relative strengths of
   lines that share an upper energy level; a fit then finds the
   amplitude of each group of lines, plus a smooth continuum, directly
   from the data. This "sky decomposition" was developed by Ivan Katkov
   and is the foundation of the semi-empirical approach below.

   - ``XSkySepIvan.py`` / ``SkySubDev2.py`` -- an older, self-contained
     copy of the decomposition vendored into ``py_progs/sky_decomp/``.
   - ``SkySubDev3.py`` / ``DecomposeCleanSky.py`` -- the newer
     decomposition, used directly from the ``lvmsky`` repository.
   - ``palace_make_mask.py`` -- uses the PALACE line lists to build the
     mask of sky-line-free pixels used by many continuum tools.

See `How the Decomposition Works`_ for what the fitted model contains.


The Semi-Empirical Machine-Learning Approach
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**What it is.** A sky predictor developed in the ``lvmsky`` repository
(Niv Drory and Ivan Katkov) that combines physical templates with
learning from LVM data. It works in three steps:

1. **Decompose.** Every sky-telescope and science spectrum in a large
   training set is fitted with the PALACE-based decomposition, so that
   each spectrum is reduced to a few hundred coefficients: amplitudes
   for groups of OH lines, atomic lines and the O2 band, and the
   parameters of a moonlight, zodiacal-light and airglow-continuum
   model.
2. **Learn.** An ensemble of neural networks (multi-layer perceptrons,
   hence ``mlp`` in the file names) is trained to predict the *science*
   field's coefficients from the coefficients of the two sky telescopes
   plus the observing geometry (positions of the Sun and Moon, airmass,
   time, solar activity).
3. **Rebuild.** For a new exposure, the sky-telescope spectra are
   decomposed, the network predicts the science-field coefficients, and
   the decomposition's templates turn those back into a sky spectrum,
   convolved with the science fibers' own LSF.

Physics enters through the line lists and component shapes, so every
predicted spectrum is built from physically meaningful pieces; the data
determine how those pieces scale from the sky telescopes to the
science field. The model does need real SkyE/SkyW spectra for each
exposure -- it is not a pure function of position and time.

**Tools in lvm_ksl** (all in ``py_dev/``; the decomposition and network
code itself is in ``lvmsky``):

- ``SelectXCF.py`` -- choose a training or test set of exposures from
  an XCframe summary file.
- ``TrainSkyModel.py`` -- the full training pipeline in one command
  (uses ``ConvertForDecompose.py``).
- ``PredictSky.py`` / ``BatchPredictSky.py`` -- predict the sky for one
  exposure, or many in parallel.
- ``EvalFluxResiduals.py``, ``EvalCoefResiduals.py``,
  ``PlotSkyResiduals.py``, ``PlotPredictSky.py`` -- evaluate a trained
  model.

**Status.** A model trained on DRP 1.2.1 data (3,501 exposures) exists,
but it has to be retrained: the data are now at DRP 1.3.2, and the
``lvmsky`` branch now in use (``skydecomp-telluric-corrected-lines``)
saves models in a new format that cannot load the old one. See
`How To: Train a New Model`_ for what that involves.


The Three Approaches at a Glance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 19 27 27 27

   * -
     - ESO Sky Model
     - PALACE (forward model)
     - Semi-empirical ML
   * - Components
     - moonlight, zodiacal light, airglow lines and continuum, thermal
     - airglow lines and continuum only
     - airglow lines (OH, O2, atomic) and continuum, moonlight,
       zodiacal light
   * - Inputs for a prediction
     - pointing, time, solar flux
     - pointing, time, solar flux
     - pointing, time **and** the SkyE/SkyW spectra
   * - Fitted to LVM data
     - nothing
     - nothing
     - template amplitudes and the network that predicts them
   * - Calibrated for
     - Paranal
     - Paranal
     - LVM (Las Campanas)
   * - Main lvm_ksl tools
     - ``EsoSkyObs.py``, ``SkySepESO.py``, ``SkyObsESO*.py``
     - ``PalaceObs.py``, ``SkyObsPalaceCompare.py``,
       ``SkySepPalace.py``
     - ``TrainSkyModel.py``, ``PredictSky.py``,
       ``BatchPredictSky.py``


How the Decomposition Works
---------------------------

Because the semi-empirical approach, ``SkySubDev2.py`` and
``SkySubDev3.py`` all rest on the PALACE-based decomposition, it helps to
know what the fitted model contains. A sky spectrum is described as a
sum of:

- **Line groups.** Each OH group is a set of lines sharing an upper
  energy level, whose relative strengths are fixed by PALACE; one
  amplitude scales the whole group. Atomic lines and the O2 band are
  handled the same way. Because the lines are placed at their catalog
  wavelengths and convolved with an LSF, the fitted amplitudes do not
  depend on the pixel grid or on the instrumental resolution.
- **Moonlight.** A solar spectrum multiplied by a smooth B-spline in
  wavelength. In the newest version the solar spectrum is first given
  the color of scattered moonlight (lunar albedo and a Rayleigh
  :math:`\lambda^{-4}` factor), so the spline only has to make a mild
  correction.
- **Zodiacal light** (newer versions only). A slightly reddened solar
  spectrum, attenuated for the target's airmass, again times a smooth
  B-spline with very few knots.
- **Airglow continuum.** Three fixed PALACE templates (HO2, FeO and the
  O2 continuum), each with a single amplitude.

All amplitudes are constrained to be non-negative. The LSF is refined
from the sky lines themselves during the fit.

On its own, the fit cannot tell moonlight from zodiacal light, or the
three airglow continua from each other, very reliably. The newest
version therefore adds constraints, using the geometry of each
exposure: limits on how fast the moonlight and zodiacal-light colors
can change, a prediction of the moonlight fraction from the positions
of the Moon and Sun (so that the moonlight term is held near zero when
the Moon is down), a bracket on the zodiacal-light brightness, and
limits on the ratios of the three airglow continua. These constraints
are applied by ``lvmsky``'s ``decompose_parallel.py`` driver, not by the
decomposition class on its own, so scripts that call the class directly
(``SkySubDev3.py``, ``DecomposeCleanSky.py``) do not get them.

Several versions of the decomposition are in use at once:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Version
     - Used by
   * - ``py_progs/sky_decomp/fit.py`` -- copied into lvm_ksl in July
       2026 (commit ``6f877c0``) from an April 2026 ``lvmsky`` version;
       moonlight + airglow continuum only, no separate zodiacal light
     - ``XSkySepIvan.py``, ``SkySubDev2.py``. Self-contained; changes in
       ``lvmsky`` do not reach it.
   * - ``lvmsky/skysub/sky_decomp/`` -- the live repository, imported
       directly (default ``~/SDSS/lvmsky/skysub``, override with
       ``-lvmsky_skysub``)
     - ``SkySubDev3.py``, ``DecomposeCleanSky.py`` and all of the
       semi-empirical ``py_dev/`` scripts. Whatever branch of ``lvmsky``
       is checked out is what runs.

As of 2026-10-07 the ``lvmsky`` checkout is on the
``skydecomp-telluric-corrected-lines`` branch, which is ahead of
``lvmsky``'s ``main`` and is where development is happening. Its
production decomposition adds the constraints above, telluric
absorption of each sky line, a continuous two-dimensional LSF model,
photon-noise weighting, and masking of bright nebular lines.


Comparing the Approaches
------------------------

Two kinds of comparison are set up:

- **Sky subtraction.** ``SkySepESO.py``, ``SkySepPalace.py``,
  ``SkySubDev2.py`` and ``SkySubDev3.py`` all read an XCframe summary
  file and write sky-subtracted spectra in a common layout, so they can
  be run through ``SkySubRun.py`` and compared with the DRP and the
  other methods using ``SkySub_eval.py`` and the nebular-line tools (see
  :doc:`sky_methods`, :doc:`sky_method_eval` and :doc:`sky_nebular_eval`).
- **Sky prediction.** ``BatchPredictSky.py`` (semi-empirical model) and
  ``BatchPredictSkyESO.py`` (ESO Sky Model) both write the observed and
  predicted sky for the same exposures in one layout (WAVE, FLUX_OBS,
  FLUX_PRED, LINE_PRED, META), which ``EvalFluxResiduals.py`` and
  ``PlotSkyResiduals.py`` evaluate in exactly the same way. The
  evaluation depends only on the observed and predicted spectra, never
  on how a model works internally, so further candidates can be added
  without changing it.

The semi-empirical model does not yet have a sky-subtraction wrapper in
the ``SkySubRun.py`` layout; ``BatchPredictSky.py`` produces a
prediction, not a sky-subtracted file.


How To: Test the Semi-Empirical Model Against New Data
------------------------------------------------------

For the case where a trained model exists and you want to see how well
it does on a fresh set of exposures.

**Prerequisites**:

- ``conda activate niv`` -- the semi-empirical scripts need PyTorch,
  which the ``lvmdrp26`` environment does not have.
- Both ``py_progs/`` and ``py_dev/`` on ``PATH`` and ``PYTHONPATH``
  (``py_dev`` is added the same way as ``py_progs``, see
  :doc:`installation`).
- The model must have been trained against the same ``lvmsky`` version
  that is checked out.

**Steps**, given a source XCframe summary file (``source.fits``, from
``SummarizeCframe.py -by fiber``, see :doc:`summarize`) and a trained
model (``model.pt``)::

    # 1. Choose a test set
    SelectXCF.py source.fits my_test_set.fits -n 200 -seed 42

    # 2. Predict the sky for each exposure (the slow step)
    BatchPredictSky.py model.pt my_test_set.fits -np 8 \
        -outfile my_test_set_batch.fits

    # 3. Quantify accuracy
    EvalFluxResiduals.py my_test_set_batch.fits \
        -out my_test_set_eval -plotdir plots_my_test_set

    # 4. Residual spectra stacked by moon brightness
    PlotSkyResiduals.py my_test_set_batch.fits \
        -outfile my_test_set_master_resid.fits -html my_test_set_master_resid.html

To compare with the ESO Sky Model on the same exposures, run
``BatchPredictSkyESO.py my_test_set.fits`` and then steps 3-4 on its
output.

``SelectXCF.py`` does not check its selection against the exposures a
model was trained on, so a test set drawn from the same source file may
overlap the training set. Keep the two apart by hand if that matters.

Step 2 decomposes two spectra per exposure; 200 exposures took about
two minutes with 8 workers.


How To: Train a New Model
-------------------------

**Steps**::

    # 0. Build the source XCframe summary from the current DRP version
    SummarizeCframe.py -by fiber -ver 1.3.2 exp_start exp_stop delta

    # 1. Choose the training set
    SelectXCF.py source.fits my_train_set.fits -n 3500 -seed 42

    # 2. Run the whole training pipeline
    TrainSkyModel.py my_train_set.fits -np 8

``TrainSkyModel.py`` runs four stages -- ``convert`` (reformat the
XCframe for ``lvmsky``), ``decompose`` (``lvmsky``'s
``decompose_parallel.py``), ``wavecache`` (assemble and filter the
training coefficients) and ``train`` (fit the network ensemble) -- and
writes everything to ``my_train_set_train/``, ending with
``mlp_ensemble.pt``. Each stage is skipped if its output already exists,
so an interrupted run can simply be restarted; ``-start_at`` and
``-stop_after`` rerun selected stages. Evaluate the result with
``EvalCoefResiduals.py my_train_set_train`` and the test recipe above.

The ``train`` stage needs enough rows for every coefficient group: in
practice at least ~150-200 exposures for a smoke test, and of order
1,000-3,500 for a real model.

**Open items before retraining on the current lvmsky branch.**
``TrainSkyModel.py`` and ``ConvertForDecompose.py`` were written against
``lvmsky``'s ``main`` (2026-09) and have not yet been updated:

- the decompose stage requests the older ``lsf-surface-iterative-split-zodi``
  fit model rather than the branch's new default;
- ``ConvertForDecompose.py`` does not carry the per-fiber LSF or the
  per-telescope airmass, PWV and pointing columns that the branch's
  geometry priors and telluric models need, so those are silently off
  (or, for the telluric models, cannot run);
- ``PredictSky.py`` does not yet pass the new per-telescope LSF inputs
  to the branch's prediction routine, nor apply its new
  post-prediction corrections.


Open Questions
--------------

- **The instrumental LSF.** Gaussian fits to isolated airglow lines
  (``lvm_line_profile.py``, see :doc:`spectral_fitting_local`) come out
  1-12% wider than the DRP's LSF, by an amount that varies from line to
  line rather than smoothly with wavelength -- more consistent with
  unresolved line blends than with a calibration error. This has only
  been tested on one pooled sample, so it is not known whether any
  residual difference is stable from exposure to exposure. A flat LSF
  broadening inside ``SkySubDev2.py`` made the fits worse.
- **The continuum.** The blue continuum excess and the Z-arm continuum
  problem seen with the ESO Sky Model are not explained. The
  semi-empirical model also fits the Z arm worst (median continuum
  residual about 2.4 times the noise, against about 1.5 in B and R).
  The newer, more constrained decomposition leaves less freedom to
  absorb such structure into the moonlight term, which should make the
  question easier to test.
- **Line amplitudes.** The DRP 1.2.1 semi-empirical model predicted the
  science-field line strengths about 6-13% too low.
- **The ESO Sky Model's own LSF** is a fixed ~0.4 A Gaussian.
  ``BatchPredictSkyESO.py`` therefore adds each exposure's LSF in
  quadrature on top of it, so the ESO and semi-empirical predictions are
  compared at the same resolution.


Where the Code Lives
--------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Location
     - Contents
   * - ``lvm_ksl/py_progs/``
     - The ESO and PALACE tools (``EsoSkyObs.py``, ``PalaceObs.py``,
       ``SkySep*.py``, ``SkyObs*.py``), the sky-subtraction methods
       including ``SkySubDev2.py`` and ``SkySubDev3.py``, the shared
       evaluator ``sky_residual_eval.py``, and the vendored
       decomposition in ``py_progs/sky_decomp/`` (reference data in
       ``data/palace_ref/``). Documented in :doc:`sky_subtraction`.
   * - ``lvm_ksl/py_dev/``
     - The semi-empirical training, prediction and evaluation scripts,
       plus ``BatchPredictSkyESO.py`` and ``DecomposeCleanSky.py``.
       Tracked in git, but outside the API documentation; scripts here
       depend on ``lvmsky`` and may move or disappear with it.
   * - ``lvmsky/skysub/``
     - Separate repository (sdss/lvmsky). The current decomposition
       (``sky_decomp/``), its batch driver (``decompose_parallel.py``),
       and the network code (``mlp_predictor/``).


See Also
--------

- :doc:`sky_subtraction` -- full documentation of the ESO, PALACE and
  decomposition-based sky tools in ``py_progs/``, and of the sky
  subtraction methods they are compared with
- :doc:`summarize` -- ``SummarizeCframe.py``, which builds the XCframe
  summary files every tool on this page reads
- :doc:`spectral_fitting_local` -- ``lvm_line_profile.py`` and the
  airglow line-shape investigation behind the LSF question
- Noll, S., et al. 2012, A&A, 543, A92 -- the ESO (Cerro Paranal)
  sky model
- Noll, S., et al. 2025 -- PALACE v1.0, Paranal Airglow Line And
  Continuum Emission model
