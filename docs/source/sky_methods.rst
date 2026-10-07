Packaged Sky-Subtraction Approaches
===================================

This page is part of :doc:`sky_subtraction`.

The scripts on this page are tools that package different sky-subtraction
approaches -- the DRP's own routine, polynomial and B-spline continuum
fits, the PALACE-based decomposition, and the ESO and PALACE sky models --
so that they can all be run the same way and their results evaluated
side by side.  Each one reads an XCframe summary file (one summary
spectrum per exposure for the science fibers and for SKY_EAST and
SKY_WEST, from ``SummarizeCframe.py``; see :doc:`summarize`), subtracts
the sky row by row, and writes the result in a common layout: WAVE, FLUX
(sky-subtracted), SKY, and DRP_ALL extensions (SkySepESO adds
MOON/ZODI/DIFFUSE; see below).  ``SkySubRun.py`` runs any of them through
one command line with a common output-naming scheme.

The outputs are then compared with the tools in :doc:`sky_method_eval`
(how well the sky is removed) and :doc:`sky_nebular_eval` (how well real
nebular emission survives).

DRP_ALL carries a QA_FLAGS column recording per-row quality issues.

**Common QA flag bits**

======  ==========  ===========================================================
0x01    NANDATA     NaN or inf found in the input flux or sky data.
0x02    ZEROSKY     Sky line vector is all-zero; scale factor is unreliable.
0x04    POORFIT     Continuum fit is poorly conditioned (SkySubOrig, SkySubDev1).
0x04    MODELFAIL   ESO sky model (local SM-01 binary + SkyCalc web-service
                    fallback) both failed for this row (SkySepESO only).
0x04    NOTSOLVED   At least one of the three per-row SkyDecomp fits
                    reported a status other than Solved/AlmostSolved
                    (SkySubDev3 only).
0x08    FAILED      Row raised an exception; spectrum filled with NaN.
======  ==========  ===========================================================

**Common DRP_ALL columns beyond QA_FLAGS**

- ``mjd`` — precise MJD recomputed from the ``obstime`` ISO timestamp string
  via ``SkySubOrig.obstime_to_mjd()``, replacing the truncated-integer MJD
  carried through from the input file (all five scripts, plus SkySubSci.py/
  SummarizeSciSky.py in :doc:`sky_from_science`).
- ``LINE_SCALE`` — the bisection factor *r* used to scale the sky lines
  (SkySubOrig, SkySubDev1, SkySepESO only; SkySubDev2 applies no scale factor,
  and SkySubDrp wraps an external lvmdrp routine that does not expose its
  internal factor).
- ``SCI_MED_<arm>``/``SCI_NMAD_<arm>``/``SCI_RMS_<arm>``/``SCI_SKEW_<arm>`` and
  the equivalent ``SKY_*`` columns (arm = ``B``, ``R``, ``Z``) — per-arm
  continuum-fit-quality statistics (SkySubOrig, SkySubDev1, SkySepESO,
  SkySubDev2 only; requires ``sky_mask.fits``, searched for automatically —
  SkySubDrp does not have these, since it wraps an external lvmdrp routine
  that doesn't expose an internal continuum fit).  These test the
  continuum fit *itself* against the raw, pre-subtraction science and sky
  spectra, independent of the line-scaling step — see
  ``GetSkyCont.arm_continuum_stats()`` and the ``SkySub_eval.py``
  "Continuum Separation" figures (:doc:`sky_method_eval`) for how they
  are used.


SkySubOrig.py
-------------

Perform per-spectrum sky subtraction on an XCframe file using a degree-4
polynomial continuum and a custom bisection search for the sky-line scale
factor.

**Usage**::

    SkySubOrig.py [-method METHOD] [-delta N] [-out ROOT] filename

**Arguments:**

filename
    XCframe FITS file to process.

**Options:**

-method METHOD
    Sky subtraction method:

    - ``nearest`` — subtract the nearest sky telescope spectrum (continuum
      and lines from the same telescope).
    - ``farthest`` — subtract the farthest sky telescope spectrum.
    - ``farlines_nearcont`` — scale sky emission lines from the far sky
      telescope and continuum from the near sky telescope (default).

-delta N
    Process every N-th row; useful for quick tests (default: 1 = all rows).

-out ROOT
    Output filename root.  Default: ``<stem>_orig_<method>``.

**Description:**

Each spectrum is decomposed into a polynomial continuum (degree 4, fitted
with sigma clipping on clean pixels) and line residuals.  Near and far sky
telescopes are identified from RA/Dec separations in DRP_ALL.  A global
line scale factor *r* is found by bisection minimising::

    sum |sci_lines × (sci_lines − r × sky_lines)| / ‖sky_lines‖²

For ``farlines_nearcont`` the sky model is::

    sky = cont_near + r × lines_far

The sky is then subtracted: ``flux_out = flux_sci − sky``.  Only the sky
lines are scaled by *r*; the continuum (``cont_near``) is used exactly as
fitted, with no additional scaling.

**Output:**

A FITS file ``<ROOT>.fits`` with extensions WAVE, FLUX (sky-subtracted),
SKY, and DRP_ALL (with QA_FLAGS, ``mjd``, ``LINE_SCALE``, and per-arm
continuum-fit-quality columns — see "Common DRP_ALL columns" above).

**See Also:** :doc:`api/SkySubOrig/index`


SkySubDrp.py
------------

Perform per-spectrum sky subtraction using the lvmdrp
``create_skysub_spectrum`` routine.

**Usage**::

    SkySubDrp.py [-method METHOD] [-delta N] [-out ROOT] filename

**Arguments:**

filename
    XCframe FITS file to process.

**Options:**

-method METHOD
    Sky subtraction method: ``nearest`` | ``farlines_nearcont`` (default).

-delta N
    Process every N-th row (default: 1).

-out ROOT
    Output filename root.  Default: ``<stem>_drp_<method>``.

**Description:**

For each row a minimal SCI/SKYE/SKYW BinTable HDUList is constructed from
the per-row FLUX, SKY_EAST, and SKY_WEST spectra together with RA/Dec
information from DRP_ALL::

    PRIMARY  (empty)
    SCI      BinTable  WAVE, FLUX, ERROR   RA/DEC from sci_ra/sci_dec
    SKYE     BinTable  WAVE, FLUX, ERROR   RA/DEC from skye_ra/skye_dec
    SKYW     BinTable  WAVE, FLUX, ERROR   RA/DEC from skyw_ra/skyw_dec

Flux errors are estimated as ``sqrt(|flux|)`` since XCframe files do not
carry IVAR.  These errors affect only the internally propagated sky error;
the sky model itself does not depend on them.  ``create_skysub_spectrum``
then selects the sky telescope and computes the sky model according to the
chosen method.

Requires the ``lvmdrp`` conda environment (``lvmdrp26``).

**Output:**

A FITS file ``<ROOT>.fits`` with extensions WAVE, FLUX (sky-subtracted),
SKY, and DRP_ALL (with QA_FLAGS and precise ``mjd``).  ``create_skysub_spectrum``
only returns the sky spectrum and its error, not an internal scale factor, so
no ``LINE_SCALE`` or continuum-fit-quality columns are written here.

**See Also:** :doc:`api/SkySubDrp/index`


SkySubDev1.py
-------------

Perform per-spectrum sky subtraction using the B-spline continuum
separation from ``GetSkyCont.py`` instead of the polynomial fit used by
SkySubOrig.  The bisection scale search and sky assembly are otherwise
identical to SkySubOrig.

**Usage**::

    SkySubDev1.py [-method METHOD] [-delta N] [-kstep N]
                  [-out ROOT] [-mask mask.fits] filename

**Arguments:**

filename
    XCframe FITS file to process.

**Options:**

-mask file
    palace_mask FITS file produced by ``palace_make_mask.py``
    (MASK extension: 1 = clean, 0 = sky-line affected).
    If omitted, ``sky_mask.fits`` is searched for in the current directory
    and then in the ``lvm_ksl data/`` directory.

-method METHOD
    ``nearest`` | ``farlines_nearcont`` (default).

-delta N
    Process every N-th row (default: 1).

-kstep N
    B-spline knot spacing in Angstroms (default: 100).

-out ROOT
    Output filename root.  Default: ``<stem>_dev1_<method>``.

**Description:**

A two-component B-spline design matrix (DIFFUSE + MOON, same as
``GetSkyCont.py``) is built once from the full wavelength grid and the
palace mask, then reused for every row.  For each row:

1. Science and sky spectra are decomposed into continuum (non-negative
   least squares on clean pixels) and line residuals.
2. Near/far sky telescopes are identified from RA/Dec separations.
3. A global line scale factor *r* is found by bisection (same objective
   function as SkySubOrig).
4. The sky model is assembled and subtracted::

       farlines_nearcont:  sky = cont_near + r × lines_far
       nearest:            sky = cont_near + r × lines_near

   ``cont_near`` is always the continuum used in the final SKY (both
   methods), used exactly as fitted with no additional scaling.
5. If ``sky_mask.fits`` is available, per-arm continuum-fit-quality stats
   (``SCI_*``/``SKY_*`` columns) are computed from the raw pre-subtraction
   science spectrum (``lines_sci``) and the near-sky spectrum
   (``lines_near``), via ``GetSkyCont.arm_continuum_stats()``.

**Output:**

A FITS file ``<ROOT>.fits`` with extensions WAVE, FLUX (sky-subtracted),
SKY, and DRP_ALL (with QA_FLAGS, precise ``mjd``, ``LINE_SCALE``, and
per-arm continuum-fit-quality columns — see "Common DRP_ALL columns" above).

**See Also:** :doc:`api/SkySubDev1/index`


SkySubDev2.py
-------------

Perform per-spectrum sky subtraction using the PALACE spectral
decomposition (``XSkySepIvan.py``), without any additional scaling.

**Usage**::

    SkySubDev2.py [-method METHOD] [-delta N] [-lsf FWHM]
                  [-lsf_boost FACTOR] [-out ROOT] filename

**Arguments:**

filename
    XCframe FITS file to process.

**Options:**

-method METHOD
    ``scilines_nearcont`` (default) | ``nearest`` | ``farlines_nearcont``.

-delta N
    Process every N-th row (default: 1).

-lsf FWHM
    Constant LSF FWHM in Angstroms (default: 1.3), used only if the input
    file has no ``LSF`` extension.

-lsf_boost FACTOR
    Multiplicative correction applied to whatever LSF is used, any source
    (default: **1.0**, i.e. no correction).  A single-row peak/area test
    on real data initially suggested a boost around 1.08 would fix
    fitted lines coming out too narrow, but a proper 280-row-median test
    showed the opposite — increasing the boost made the noise-corrected
    HF RMS residual monotonically *worse*, and made the oversubtracted-
    center/adjacent-bump pattern at [OI] 5577 monotonically deeper, not
    smaller. Left available for further investigation, but do not
    re-enable as the default without re-testing against a multi-row
    sample.

-out ROOT
    Output filename root.  Default: ``<stem>_dev2_<method>``.

**Description:**

The PALACE decomposer (``SkyDecomp``) needs an assumed LSF, chosen in
priority order:

1. This file's own per-row, per-wavelength LSF, if it has an ``LSF``
   extension (written by ``SummarizeCframe.py``'s ``make_med_spec`` from
   the source CFrame's own LSF — one row of FWHM-vs-wavelength per
   exposure).  The decomposer is rebuilt whenever a row's LSF differs
   from the previous row's (see "Output" below for the performance cost).
2. Otherwise, a representative wavelength-dependent reference curve
   (``data/lsf.fits``, derived once from real per-row LSF data), if
   found — the same curve every row, so the decomposer is only ever
   built once for the whole file (no extra per-row cost).
3. Otherwise, the constant ``-lsf`` FWHM, built once and reused for
   every row.

A flat 1.3 Å FWHM (case 3's default) is normally too narrow for at least
part of the real spectrum — cases 1 and 2 both capture the true
wavelength dependence instead.  Note, however, that even case 1 (a real,
correct per-row LSF) does not fully eliminate an oversubtracted-center/
adjacent-bump residual pattern seen at bright lines like [OI] 5577; a
flat multiplicative width correction (``-lsf_boost``) was tested and
found to make that pattern worse, not better, so the residual's root
cause is evidently not a pure LSF-width deficit and remains under
investigation (see ``lvm_line_profile.py`` in :doc:`spectral_fitting_local`,
built specifically to investigate this independently of PALACE — it
finds a real but line-specific, non-smooth-in-wavelength gap between
the LSF extension's stated FWHM and an independently-fit Gaussian FWHM
on raw sky lines, more consistent with several catalog lines being
unresolved blends than with a genuine LSF calibration error).  For each
row:

1. Near and far sky telescopes are identified from RA/Dec separations.
2. The near-sky spectrum, and whichever other spectrum the chosen method
   needs (far-sky for ``farlines_nearcont``, the science flux itself for
   ``scilines_nearcont``), are decomposed by PALACE into emission-line
   and continuum components::

       LINES = oh + atom + orc + o2
       CONT  = moon + diffuse

3. The sky model is assembled without any scale factor::

       scilines_nearcont:  sky = CONT_near + LINES_sci
       farlines_nearcont:  sky = CONT_near + LINES_far
       nearest:            sky = CONT_near + LINES_near

   CONT always comes from the near-sky telescope in all three methods —
   fitting a continuum from the science spectrum itself would be
   confounded by real astrophysical continuum, so PALACE's CONT
   component is never taken from there, only its LINES component
   (``scilines_nearcont``).  ``scilines_nearcont`` is the default:
   fitting the sky lines directly from the fiber being corrected avoids
   relying on a sky telescope's lines matching the science fiber's
   actual sky-line amplitude, which ``farlines_nearcont``/``nearest``
   implicitly assume.

4. ``flux_out = flux_sci − sky``.
5. If ``sky_mask.fits`` is available, per-arm continuum-fit-quality
   stats are computed from the raw science and near-sky spectra via
   ``GetSkyCont.arm_continuum_stats()``.  For ``scilines_nearcont`` the
   science-side decomposition needed here is the same one already
   computed for the subtraction (no extra cost); for the other two
   methods it's an extra PALACE decomposition run purely for this check.

Requires the PALACE library, vendored in ``py_progs/sky_decomp/`` (260709;
no longer an external dependency, and no longer requires ``lvmdrp`` to be
installed) and its reference data in ``data/palace_ref/``; paths are taken
from ``XSkySepIvan.py``'s ``DEFAULT_BASE_DIR``.  ``clarabel`` (the QP
solver used for the fits) is still a separate, normal dependency.

Rebuilding the PALACE decomposer for a new LSF costs roughly 2 seconds
(measured), on top of the roughly 1 second already spent per PALACE
decomposition in ``one_drp()`` — real added cost across a full run when
using this file's own per-row LSF extension (not just a one-time setup
cost, as with the constant-LSF or reference-curve paths), since
``_get_decomposer()`` only ever caches one decomposer instance at a time
and rebuilds whenever the row-to-row LSF changes.  The reference curve
(``data/lsf.fits``, searched for in the current directory then in the
lvm_ksl ``data/`` directory, same convention as ``sky_mask.fits``) is the
same array every row, so it costs nothing extra despite going through the
same per-row code path.  Any non-finite or non-positive FWHM value at a
given wavelength in a row's own LSF array falls back to the constant
``-lsf`` default at that wavelength only (not the whole row).

**Output:**

A FITS file ``<ROOT>.fits`` with extensions WAVE, FLUX (sky-subtracted),
SKY, and DRP_ALL (with QA_FLAGS, precise ``mjd``, and per-arm
continuum-fit-quality columns — see "Common DRP_ALL columns" above).  No
scale factor is applied anywhere in this method (see docstring: "no
scaling"), so there is still no ``LINE_SCALE`` column, unlike SkySubOrig/
SkySubDev1/SkySepESO.  When an ``LSF`` extension or reference curve was
used, DRP_ALL also gains ``LSF_FWHM_MED`` (this row's median LSF FWHM
actually used), and the primary header gains ``LSFSRC`` recording which of
the three LSF sources was used for the whole run.

**See Also:** :doc:`api/SkySubDev2/index`


SkySubDev3.py
-------------

Perform per-spectrum sky subtraction using
``sky_decomp.lsf_surface_iterative.SkyDecompLSFSurfaceIterative`` (the
same physically-motivated, per-row-LSF-refined decomposition
``DecomposeCleanSky.py`` uses — see :doc:`sky_nebular_eval`) for the continuum/line separation, with the near/far
recipe and bisection scale search otherwise identical to SkySubDev1.py.
The decomposition is imported from the ``lvmsky`` repository; see
:doc:`sky_models` for what it fits and which ``lvmsky`` version is in
use.

**Usage**::

    SkySubDev3.py [-method METHOD] [-delta N] [-v VEL] [-lmc] [-smc]
                  [-out ROOT] [-lvmsky_skysub PATH] filename

**Arguments:**

filename
    XCframe FITS file to process.

**Options:**

-method METHOD
    ``nearest`` | ``farlines_nearcont`` (default).

-delta N
    Process every N-th row (default: 1).

-v VEL / -lmc / -smc
    Nebular systemic velocity (km/s) used to Doppler-shift the exclusion
    windows masked out of every spectrum's continuum/line fit (default:
    0). Same convention as ``sky_gaussfit.py``/``DecomposeCleanSky.py``.

-out ROOT
    Output filename root. Default: ``<stem>_dev3_<method>``.

-lvmsky_skysub PATH
    Path to the ``lvmsky`` repo's ``skysub/`` directory, which supplies
    the ``sky_decomp`` package this script imports (default:
    ``~/SDSS/lvmsky/skysub``). Unlike SkySubDev1/Dev2, no ``-mask``
    argument exists — ``SkyDecomp`` fits its own OH/atomic-line families
    directly against the data, needing no external palace_mask file.

**Description:**

The ``SkyDecompLSFSurfaceIterative`` instance is built once (from the
full wavelength grid) and reused for every row, via
``DecomposeCleanSky.build_decomp`` (imported directly, not duplicated).
For each row:

1. Science and sky spectra are read; near/far sky is identified from
   RA/Dec separations (identical to SkySubDev1.py).
2. Each of the three spectra (science, near sky, far sky) is fit with
   SkyDecomp, with ``sky_gaussfit.resolve_nebular_lines()``'s windows
   (Doppler-shifted by ``-v``/``-lmc``/``-smc``) excluded from the fit
   via ``ivar=0`` — the same mask-and-wrap approach
   ``DecomposeCleanSky.py`` uses. This matters specifically for the two
   sky-telescope spectra: the nebular-leak validation work in :doc:`sky_nebular_eval` found
   real nebular-line leak in SKY_WEST on at least one tested exposure,
   so masking it out of the *continuum* fit (rather than assuming, as
   SkySubDev1/Dev2/Drp implicitly do, that the sky telescopes are
   nebula-free) keeps that leak from biasing the fitted continuum.
   ``cont = bestfit_lsf - (oh+atom+orc+o2)``; ``lines = spectrum - cont``
   (the same "lines = observed - continuum" definition SkySubDev1.py
   uses, kept identical so the two are comparable apples-to-apples).
3. A global line scale factor *r* is found by ``ksl_bisection`` (from
   SkySubOrig.py), exactly as in SkySubDev1.py.
4. The sky model is assembled and subtracted::

       farlines_nearcont:  sky = cont_near + r × lines_far
       nearest:            sky = cont_near + r × lines_near

**Output:**

A FITS file ``<ROOT>.fits`` with extensions WAVE, FLUX (sky-subtracted),
SKY, and DRP_ALL — the same layout as SkySubOrig/Drp/Dev1/Dev2.py, so
SkySub_eval.py and SkySubNebEval.py both read it exactly like those, with
no changes needed there. QA_FLAGS reuses 0x01/NANDATA and 0x02/ZEROSKY,
adds a new 0x04/NOTSOLVED (at least one of the three per-row SkyDecomp
fits reported a status other than Solved/AlmostSolved), and 0x08/FAILED.

**See Also:** :doc:`api/SkySubDev3/index`


SkySubRun.py
------------

Dispatcher for the SkySub* family (SkySubDrp.py, SkySubOrig.py,
SkySubDev1.py, SkySubDev2.py, SkySubDev3.py): runs one named routine on
an XCframe file by calling its ``do_all()`` directly (no subprocess),
writes its output into a routine-specific subdirectory so different
routines — or repeated runs of the same routine with a different variant
— never collide on a filename, and optionally chains ``SkySub_eval.py``
on the result afterward. This is the recommended entry point for running
and comparing the five methods, rather than calling each script's own
CLI directly: it removes the bookkeeping of tracking which routine
produced which file, and its own naming convention
(``DIR/<routine>/<input_stem>_<routine>_<variant>.fits``) is what
SkySubNebEval.py/PlotSkySubNebEval.py's own multi-file comparisons
(:doc:`sky_nebular_eval`) expect.

**Usage**::

    SkySubRun.py -routine {drp,orig,dev1,dev2,dev3} [-variant NAME]
                 [-delta N] [-mask FILE] [-kstep N]
                 [-fwhm_lsf F] [-lsf_boost F]
                 [-v VEL] [-lmc] [-smc] [-lvmsky_skysub PATH]
                 [-outdir DIR] [-eval] [-eval_out ROOT]
                 filename

**Arguments:**

filename
    XCframe FITS file to process.

**Options:**

-routine NAME
    Which SkySub* routine to run — ``drp`` | ``orig`` | ``dev1`` |
    ``dev2`` | ``dev3`` (required; no default, so a run always names its
    own routine explicitly).

-variant NAME
    The sky-construction variant passed as that routine's own
    ``-method`` (default: that routine's own default variant).

-delta N
    Process every N-th row (default: 1).

-mask FILE / -kstep N
    Passed through to SkySubDev1.py only.

-fwhm_lsf F / -lsf_boost F
    Passed through to SkySubDev2.py only (defaults 1.3 / 1.0).

-v VEL / -lmc / -smc / -lvmsky_skysub PATH
    Passed through to SkySubDev3.py only.

-outdir DIR
    Top-level directory under which each routine gets its own
    subdirectory, ``DIR/<routine>/`` (default: ``sky_runs``).

-eval
    After the routine finishes, also run ``SkySub_eval.plot_eval()`` on
    its output (a convenience — running ``SkySub_eval.py`` by hand on the
    written file afterward gives identical results).

-eval_out ROOT
    Output root for the ``-eval`` HTML (default: the routine's own
    output stem + ``"_eval"``, in the same subdirectory).

**Description:**

Each routine's module is imported lazily (only once ``-routine``
selects it), so running e.g. ``-routine drp`` never pays the cost of
importing dev3's ``lvmsky`` dependency or dev1/dev2's GetSkyCont/PALACE
machinery. Calls each routine's ``do_all()`` directly within this
process (not via subprocess) — faster, and errors surface as normal
Python tracebacks rather than being swallowed into a subprocess return
code.

**Example**::

    SkySubRun.py -routine dev3 -lmc -eval XCframe_file.fits

    # writes sky_runs/dev3/XCframe_file_dev3_farlines_nearcont.fits
    # and    sky_runs/dev3/XCframe_file_dev3_farlines_nearcont_eval.html

**See Also:** :doc:`api/SkySubRun/index`


SkySepESO.py
------------

Perform per-spectrum sky subtraction using the **ESO Sky Model** to separate
the sky continuum into its physical MOON, ZODI, and DIFFUSE components,
using the same bisection line-scale search as SkySubOrig.  Ported from a
prototype developed in ``lvm_sky2506`` (``SkySepMod.py`` +
``SkySubModDev250624.ipynb``).

**Usage**::

    SkySepESO.py [-method METHOD] [-delta N] [-out ROOT] filename

**Arguments:**

filename
    XCframe FITS file to process.

**Options:**

-method METHOD
    ``nearest`` | ``farthest`` | ``farlines_nearcont`` (default) — same
    meaning as SkySubOrig.

-delta N
    Process every N-th row (default: 1).  Each row requires 2–3 ESO sky
    model fetches (science fiber, plus one or two sky fibers depending on
    method), so this script is far slower per row than the other four —
    measured at ~2.5s/row for ``farlines_nearcont`` (3 model fetches + 3
    continuum fits), i.e. roughly an hour for a full ~1750-fiber exposure.
    Use ``-delta`` liberally for quick tests.

-out ROOT
    Output filename root.  Default: ``<stem>_eso_<method>``.

**Description:**

For each row, the ESO sky model is fetched for the relevant fiber's
coordinates and observation time (``EsoSkyObs.run_sky_obs``, ``engine='auto'``:
local ESO SM-01 binary first; falls back automatically to the SkyCalc web
service if the local model call fails — e.g. because the model
rejects the target/Moon geometry).  The model gives MOON/ZODI/DIFFUSE
continuum templates, which are fit to the observed flux as a non-negative,
iteratively-downweighted 3-component linear combination
(``CONT = a·MOON + b·ZODI + c·DIFFUSE``) so that sky/airglow line pixels do
not bias the continuum fit upward.  The resulting line residual is then
scaled against the science spectrum's own line residual by the same
bisection search used in SkySubOrig, and the sky is subtracted::

    sky = sky_CONT + r × sky_LINES
    flux_out = flux_sci − sky

Only the sky lines are scaled by *r*; ``sky_CONT`` (whichever sky fiber's
ESO-model fit produced it — near for ``nearest``/``farlines_nearcont``, far
for ``farthest``) is used exactly as fitted.

Each ESO-model fetch writes a small per-call FITS file
(``EsoSkyModel_<pid>_<uuid>.fits``, a process-and-call-unique name rather
than the content-derived ``SkyE_*.fits`` default, so two concurrent
callers requesting the same or coincidentally same-rounded geometry can't
race on each other's read-then-delete) to the current working directory;
this file is read and deleted immediately, so nothing accumulates on disk
across a run.  ``EsoSkyObs.run_local``/``run_remote``'s own internal
scratch files (``config``/``output``/``data``) are independently isolated
per call in a private temporary directory -- see ``EsoSkyObs.py`` in
:doc:`sky_model_tools`.

**Output:**

A FITS file ``<ROOT>.fits`` with extensions WAVE, FLUX (sky-subtracted),
SKY, **MOON**, **ZODI**, **DIFFUSE** (the *scaled* ESO-model components of
whichever sky fit produced SKY's continuum term — they sum to SKY's
continuum; the remaining line term is ``SKY − (MOON+ZODI+DIFFUSE)``), and
DRP_ALL.  DRP_ALL carries QA_FLAGS (with the additional ``MODELFAIL`` bit
0x04 — both the local ESO model and the SkyCalc fallback failed for this
row), precise ``mjd``, ``LINE_SCALE``, the fitted ESO-model coefficients
(``SCI_MOON``/``SCI_ZODI``/``SCI_DIFFUSE``, ``SKY_MOON``/``SKY_ZODI``/
``SKY_DIFFUSE``), per-arm continuum-fit-quality columns (see "Common
DRP_ALL columns" above), and ``ERROR_MSG`` — a per-row failure-reason
string (truncated to 200 characters) for any row with a non-zero QA_FLAGS
value.  A summary of flagged rows is printed at the end of the run and also
written to ``<ROOT>_errors.txt``.

**Requirements:**

The local ESO SM-01 sky model binary, gated by the ``ESO_SKY_MODEL``
environment variable (see ``EsoSkyObs.py`` in :doc:`sky_model_tools`); without it, every row
falls back to the SkyCalc web service (slower, needs network access).

**See Also:** :doc:`api/SkySepESO/index`


SkySepPalace.py
---------------

Perform per-spectrum sky subtraction using the **PALACE airglow model**
instead of the ESO sky model (SkySepESO.py's approach) -- a deliberately
scoped-down first cut::

    SKY = sum_species  a_species * PALACE_species  +  a_moon * MOON     (all a >= 0)

fit jointly by non-negative least squares (all 9 PALACE species
amplitudes plus one MOON amplitude at once, full spectrum, no line
masking needed) against the nearest sky fiber's own spectrum.

**Usage**::

    SkySepPalace.py [-delta N] [-out ROOT] filename

**Arguments:**

filename
    XCframe FITS file to process.

**Options:**

-delta N
    Process every N-th row; useful for quick tests (default: 1 = all rows).

-out ROOT
    Output filename root; default is ``<stem>_palace``.

**Description:**

For each row, RA/Dec/obstime for the science fiber and both sky
telescopes are read from DRP_ALL, and the nearer sky telescope is
identified (as in SkySepESO.py).  PALACE (via ``PalaceObs.predict()``,
**not** the homogenized ``do_one()`` output -- see Notes) is queried once
per unique (near sky position, obstime) pair, cached across rows since a
real XCframe file shares only two sky-telescope positions across
~1700+ fibers, for all 9 species at high native resolution
(``resol=20000``, ``dlam~0.2 A`` -- matching what ``PalaceObs.py``'s own
``do_one()`` output now defaults to as well).  Each species' native
spectrum is flux-conserving rebinned and LSF-convolved onto the
instrument grid using the row's real LSF (same 3-tier priority as
SkySubDev2.py: this file's own LSF extension, else a reference curve,
else a flat default).  The 9 species templates plus MOON are then fit by
NNLS to the near-sky fiber's flux; ``SKY`` is the fit, ``FLUX = FLUX_sci
- SKY``.

**MOON:** a single, static, precomputed spectral shape
(``data/moon_base_spectrum.dat``, built offline by ``MakeMoonBase.py``)
with exactly one free amplitude -- not the flexible multi-knot B-spline
envelope SkySubDev2.py uses for its own Moon component.  The template is
the real ESO Sky Model's own MOON column, evaluated once for one fixed
reference geometry (an earlier version tried the bare solar spectrum,
then solar x ROLO lunar albedo -- both wrong in the same direction;
validated against the real ESO Sky Model that atmospheric Rayleigh/
aerosol scattering reverses the albedo reddening, so the current template
uses ESO's own already-scattered MOON output instead).  Neither the
template's shape (fixed reference geometry) nor its amplitude (entirely
free in the fit) depends on this observation's real lunar phase, moon
altitude, or moon-target separation -- flagged as the top item still
needed.  Zodiacal light is not included at all.

**Output:**

A FITS file ``<ROOT>.fits`` with extensions WAVE, FLUX (sky-subtracted),
SKY, and DRP_ALL.  DRP_ALL gains per-row ``PALACE_<SPECIES>`` (the 9
species amplitudes), ``PALACE_MOON``, ``QA_FLAGS``, ``ERROR_MSG``, and --
if ``sky_mask.fits`` is found -- per-arm continuum-fit-quality columns
(``SCI_MED_<arm>``/``SKY_MED_<arm>`` etc., directly comparable to
SkySepESO.py's own columns of the same name; see "Common DRP_ALL columns"
above), evaluated only from the HO2/FeO/MOON (continuum-like) subset of
the fit.

**QA flag bits:**

======  ==========  ===========================================================
0x01    NANDATA     NaN/inf found in input flux or sky data.
0x04    MODELFAIL   PALACE prediction failed for this row's geometry (e.g.
                    sky position below the horizon).
0x08    FAILED      Row raised an exception; spectrum filled with NaN.
======  ==========  ===========================================================

**Notes:**

Uses ``PalaceObs.predict()``'s raw PALACE units (Rayleighs/nm), not
``do_one()``'s homogenized physical-unit output (erg/s/cm^2/Angstrom) --
see PalaceObs.py in :doc:`sky_model_tools`.  Since each species/MOON amplitude is entirely
free in the NNLS fit, this is not a problem for the sky subtraction
itself (any constant unit-conversion factor is absorbed into the fitted
amplitude, since the fit only needs the right relative wavelength shape,
not an absolute scale), but it does mean the ``PALACE_<SPECIES>``/
``PALACE_MOON`` columns written to DRP_ALL are **not** physically
meaningful brightnesses -- they can't be compared across species or
across observations in absolute terms, since the missing unit conversion
and the missing at-telescope-vs-above-atmosphere extinction correction
are both silently baked into that one free scalar per component.

Deliberately not yet included (in order of what's needed next): real
lunar-phase/airmass dependence for MOON's shape or a physical prediction
of its amplitude; zodiacal light continuum; a sky-to-science line-scale
correction (SkySepESO.py's bisection step) -- not clearly motivated yet
since the NNLS amplitudes already give each species its own free scaling,
unlike ESO's single lumped LINES template.

**Requirements:**

The PALACE Python package (external dependency; see ``PalaceObs.py``
in :doc:`sky_model_tools` and Readme.md).

**See Also:** :doc:`api/SkySepPalace/index`


See Also
--------

- :doc:`sky_method_eval` - evaluating and comparing the output of
  these methods, including a complete run-and-compare workflow
- :doc:`sky_nebular_eval` - comparing the methods by how well they
  recover real nebular emission
- :doc:`sky_models` - the physical models behind ``SkySepESO.py``,
  ``SkySepPalace.py``, ``SkySubDev2.py`` and ``SkySubDev3.py``
- :doc:`api/SkySubOrig/index` - API documentation
- :doc:`api/SkySubDrp/index` - API documentation
- :doc:`api/SkySubDev1/index` - API documentation
- :doc:`api/SkySubDev2/index` - API documentation
- :doc:`api/SkySubDev3/index` - API documentation
- :doc:`api/SkySubRun/index` - API documentation
- :doc:`api/SkySepESO/index` - API documentation
- :doc:`api/SkySepPalace/index` - API documentation
