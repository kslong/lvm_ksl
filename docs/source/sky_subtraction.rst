Sky Subtraction
===============

Sky subtraction is one of the most critical steps in LVM data reduction.
LVM uses dedicated sky telescopes (SKY_EAST and SKY_WEST) to measure the
sky spectrum simultaneously with the science observations, and the DRP
uses these sky spectra to subtract the sky from the science fibers.

The lvm_ksl package provides tools for evaluating the DRP's sky
subtraction, for repairing it when a sky telescope was badly placed, and
for experimenting with alternative approaches -- both ones that still use
the sky telescopes (SkyCorr, continuum fits, PALACE, the ESO Sky Model)
and ones that take the sky from the science IFU itself. They are
described on the pages below.

Measuring how uniform the sky lines are across the IFU (an instrumental
flat-field question rather than a sky-subtraction one) is covered in
:doc:`data_quality`.


Which Page Do I Need?
---------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - If you want to...
     - See
   * - check the DRP's sky subtraction in an SFrame, or recover an
       exposure whose sky telescope was badly placed (e.g. near the Moon)
     - :doc:`sky_drp` (``eval_sky.py``, ``sky_plot.py``,
       ``SubstituteSky.py``, ``RunSky.py``)
   * - subtract the sky with ESO's SkyCorr
     - :doc:`sky_skycorr` (``Prep4SkyCorr.py``, ``RunSkyCorr.py``,
       ``SkySub.py``)
   * - run several sky-subtraction approaches on the same XCframe
       summary file so they can be compared
     - :doc:`sky_methods` (``SkySubOrig.py``, ``SkySubDrp.py``,
       ``SkySubDev1/2/3.py``, ``SkySepESO.py``, ``SkySepPalace.py``,
       ``SkySubRun.py``)
   * - measure and compare how well those approaches remove the sky
     - :doc:`sky_method_eval` (``SkySub_eval.py``,
       ``sky_residual_eval.py``)
   * - compare approaches by how well real nebular emission survives, or
       test whether the sky telescopes see nebular emission themselves
     - :doc:`sky_nebular_eval` (``DecomposeCleanSky.py``,
       ``SkySubNebEval.py``, ``PlotSkySubNebRun.py`` and related tools)
   * - take the sky from the science fibers instead of the sky telescopes
     - :doc:`sky_from_science` (``SkySubSci.py``, ``SummarizeSciSky.py``,
       ``SkySubPatch.py``, ``lsf_kernel.py``)
   * - build a sky-line mask, collect sky spectra, or fit the sky
       continuum
     - :doc:`sky_continuum` (``palace_make_mask.py``, ``XSkySepIvan.py``,
       ``GetSky_from_CFrame_sum.py``, ``GetSkyCont.py``)
   * - predict the sky for a pointing and time with the ESO Sky Model or
       PALACE
     - :doc:`sky_model_tools` (``EsoSkyObs.py``, ``PalaceObs.py``)
   * - understand the physically based sky models and the
       semi-empirical machine-learning approach
     - :doc:`sky_models`

Most of the comparison tools work on XCframe summary files -- one summary
spectrum per exposure for the science fibers and for each sky telescope --
made by ``SummarizeCframe.py`` (see :doc:`summarize`).

.. toctree::
   :maxdepth: 1

   sky_drp
   sky_skycorr
   sky_methods
   sky_method_eval
   sky_nebular_eval
   sky_from_science
   sky_continuum
   sky_model_tools
   sky_models


Notes
-----

- Sky subtraction quality depends strongly on observing conditions
- The SKY_EAST and SKY_WEST telescopes point at different positions
  than the science telescope, which can lead to spatial sky variations
- SkyCorr can model and remove airglow emission lines more accurately
  than simple scaling methods
- Theoretical sky models are useful for identifying instrumental
  artifacts vs. real sky features
- A method's generic sky-residual quality (SkySub_eval.py,
  sky_residual_eval.py) and its nebular-line recovery quality
  (SkySubNebEval.py) are independent axes -- a method can suppress
  airglow well while still distorting real nebular signal, or vice versa;
  check both before trusting a single "which method is better" answer
- A sky taken from the science field itself (SkySubSci.py,
  SkySubPatch.py) avoids the sky telescopes' different pointing, but in an
  extended source it removes the emission present in the chosen fibers
  along with the sky; the result measures emission relative to the
  faintest part of the field
- Fiber-to-fiber differences in throughput (1-3%), line width and
  wavelength (~10 mÅ) limit how well any single sky spectrum subtracts
  bright sky lines from every fiber; SkySubPatch.py corrects for all
  three, using groups of neighbouring fibers, which are much more alike
  than random pairs



See Also
--------

- :doc:`sky_models` - The physically based sky models (ESO Sky Model,
  PALACE, semi-empirical machine learning) and how they are compared
- :doc:`summarize` - ``SummarizeCframe.py``, which builds the XCframe
  summary files, and tools for looking at many exposures at once
- :doc:`data_quality` - ``QualSFrame.py``/``QualCFrame.py``, and sky-line
  flatness across the IFU
- :doc:`spectral_fitting_local` - ``sky_gaussfit.py`` and
  ``lvm_line_profile.py``
