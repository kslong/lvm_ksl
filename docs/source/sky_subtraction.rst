Sky Subtraction
===============

Sky subtraction is one of the most critical steps in LVM data reduction.
The lvm_ksl package provides tools for both evaluating the standard DRP
sky subtraction and for experimenting with alternative approaches using
the ESO SkyCorr tool.

This page describes the available sky subtraction and sky modeling tools.


Overview
--------

LVM uses dedicated sky telescopes (SKY_EAST and SKY_WEST) to measure the
sky spectrum simultaneously with science observations. The DRP uses these
sky spectra to subtract the sky from the science fibers.

The tools in lvm_ksl allow you to:

- Evaluate sky subtraction quality in DRP-processed data
- Correct an exposure whose DRP sky subtraction is compromised by a bad
  sky-telescope pointing (e.g. too close to the Moon), by substituting in
  a clean sky-telescope's fiber data from a different exposure and
  rerunning the DRP's own sky-subtraction routine (SubstituteSky.py,
  RunSky.py)
- Run alternative sky subtraction using ESO's SkyCorr tool, a polynomial or
  B-spline continuum fit (SkySubOrig/SkySubDev1), a PALACE decomposition
  (SkySubDev2), a SkyDecomp continuum with nebular-line masking
  (SkySubDev3), the lvmdrp routine directly (SkySubDrp), or the ESO Sky
  Model itself (SkySepESO) — dispatched uniformly through ``SkySubRun.py``
  so a run/evaluate workflow doesn't need to remember which routine's own
  CLI to call
- Generate theoretical sky models for comparison
- Visualize sky residuals and identify problems, including separately
  evaluating sky *line* subtraction quality and continuum *separation*
  quality (SkySub_eval.py)
- Test whether the SKY_EAST/SKY_WEST telescopes themselves leak real
  nebular-line flux (an assumption every method above depends on), and
  compare how well each SkySub* method recovers known nebular-line ratios
  on the science fiber itself (``DecomposeCleanSky.py``,
  ``sky_nebular_leak_eval.py``, ``SkySubNebEval.py`` — see "Nebular-Line-
  Based Method Evaluation" below)


Evaluating Sky Subtraction
--------------------------

eval_sky.py
^^^^^^^^^^^

Creates diagnostic plots to evaluate sky subtraction quality in SFrame
files (sky-subtracted data from the DRP).

**Usage**::

    eval_sky.py filename1 filename2 ...

Plots median, min, max, and average spectra across science fibers to
identify residual sky features.

sky_plot.py
^^^^^^^^^^^

Similar to eval_sky.py, creates plots to evaluate sky subtraction
quality with additional visualization options.

**Usage**::

    sky_plot.py filename1 filename2 ...


Correcting and Rerunning DRP Sky Subtraction
---------------------------------------------

These two scripts work together to recover an exposure whose DRP sky
subtraction is compromised by a bad sky-telescope pointing (e.g. SkyW
too close to the Moon): substitute in a clean sky telescope's fiber
data from a different exposure, then rerun the production DRP's own
sky-subtraction routine on the result.

SubstituteSky.py
^^^^^^^^^^^^^^^^^

Replaces one sky telescope's fiber data (FLUX/IVAR/MASK/LSF) and
associated header metadata in an lvmCFrame with the corresponding data
from a different lvmCFrame.

**Usage**::

    SubstituteSky.py [-h] [-o outfile] target_cframe target_tel
                      source_cframe source_tel

**Arguments:**

target_cframe
    lvmCFrame file to be corrected.

target_tel
    ``SkyE`` or ``SkyW`` (case-insensitive) -- telescope whose fibers
    in ``target_cframe`` are replaced.

source_cframe
    lvmCFrame file to draw the replacement data from.

source_tel
    ``SkyE`` or ``SkyW`` (case-insensitive) -- telescope in
    ``source_cframe`` supplying the replacement data.

**Options:**

-o outfile
    Output filename.  Default: ``target_cframe``'s basename (directory
    stripped) with ``.sky_subst`` inserted before the extension,
    written to the current directory.  An existing outfile is
    overwritten, with a warning printed first.

**Description:**

An lvmCFrame's SLITMAP assigns every fiber to a telescope (Sci, SkyE,
SkyW, Spec) -- fixed by fiber-plugging hardware, identical
fiberid-for-fiberid across all exposures.  The DRP's own sky-
subtraction routine (``skyMethod.quick_sky_subtraction``) builds its
sky spectrum from the raw FLUX/IVAR at the fibers tagged SkyE/SkyW in
the CFrame being reduced -- not the extrapolated SKY_EAST/SKY_WEST
extensions, which the current production method ignores.  So fixing a
contaminated sky telescope means replacing the FLUX, IVAR, MASK, and
LSF rows for that telescope's fibers, matched by fiberid.

Since SkyE/SkyW fiber assignment is fixed hardware, requesting
different telescopes on the two sides (``target_tel`` != ``source_tel``)
fails with a clear error -- there is no physically meaningful
fiber-by-fiber correspondence between them.

Every PRIMARY header keyword tied to the source telescope (pointing,
altitude, airmass, guider frames, sky-field name, heliocentric
velocity, moon/shadow geometry, ecliptic coordinates, etc.) is copied
too, renamed to the target telescope's own keyword names -- except the
SKYEW/SKYWW combination weights, a joint SkyE+SkyW property recomputed
elsewhere.  ``SKY SCI_SKYW_SEP`` (or the SkyE equivalent) is relative
to the *target's* own science pointing, so it is recomputed from the
target's real SCIRA/SCIDEC and the newly-copied sky position
(``lvmdrp.core.sky.ang_distance``) rather than copied as-is.  New
``SKY SUBST_*`` provenance keywords record what was substituted and
from where.  ``target_cframe``/``source_cframe`` are never modified.

**Output:**

An lvmCFrame FITS file, structurally identical to the input, with the
named telescope's FLUX/IVAR/MASK/LSF rows and header block replaced.

**See Also:** :doc:`api/SubstituteSky/index`


RunSky.py
^^^^^^^^^

Reruns the DRP's own ``quick_sky_subtraction`` on a single lvmCFrame
file -- typically the output of ``SubstituteSky.py`` above -- producing
a corrected lvmSFrame.

**Usage**::

    RunSky.py [-h] filename

**Arguments:**

filename
    lvmCFrame file to run sky subtraction on.

**Description:**

Calling ``quick_sky_subtraction`` outside the full ``science_reduction``
pipeline exposes three missing-directory bugs in ``lvmdrp`` (that
pipeline happens to pre-create these directories as a side effect of an
earlier step, so they never surface there): its own ancillary skytable
write, ``writeFitsData``'s output write when given a bare filename, and
``run_qa``'s skyQA PDF write.  ``RunSky.py`` works around all three
locally rather than patching the vendored ``lvmdrp`` package, and reads
back the freshly-written ancillary skytable -- by constructing the same
path ``quick_sky_subtraction`` uses internally, via ``lvmdrp``'s own
``path``/``drpver`` -- rather than searching for a pre-existing one, so
the diagnostic plots reflect this run's own sky model, not a stale
skytable from some other exposure or DRP version.

**Output:**

- ``lvmSFrame-<...>.fits`` (or ``sframe_<name>`` if the input name
  doesn't contain ``CFrame``) in the current directory
- ``qa/skyQA_<expnum>.pdf`` diagnostic plot
- diagnostic matplotlib figures for the SCI/SkyE/SkyW/SkyE_super/
  SkyW_super mean spectra

**See Also:** :doc:`api/RunSky/index`


Alternative Sky Subtraction with SkyCorr
----------------------------------------

ESO's SkyCorr is a tool for sky subtraction that can handle cases where
the sky varies spatially or temporally. The following scripts provide
an interface to SkyCorr for LVM data.

Prep4SkyCorr.py
^^^^^^^^^^^^^^^

Prepares LVM CFrame data for use with SkyCorr by extracting mean spectra
from science and sky fibers into a format SkyCorr can read.

**Usage**::

    Prep4SkyCorr.py [-h] [-dir whatever] [-all] filename1 filename2 ...

**Options:**

-h
    Print help and exit.

-dir path
    Search for files in the specified directory.

-all
    Process all lvmCFrame files in the directory.

**Output:**

Creates FITS files with mean spectra suitable for SkyCorr input.

RunSkyCorr.py
^^^^^^^^^^^^^

Runs the ESO SkyCorr tool on prepared input files.

**Usage**::

    RunSkyCorr.py sci.fits sky.fits

where ``sci.fits`` is the science spectrum and ``sky.fits`` is the sky
spectrum (both prepared by Prep4SkyCorr.py).

**Requirements:**

1. SkyCorr must be installed and the ``skycorr`` executable in your PATH
2. A parameter file ``lvm_base.par`` must exist in the working directory

**Setup:**

On first run, if ``lvm_base.par`` does not exist, the script will create
a template. You must edit this file to set the correct paths::

    INST_DIR=/path/to/skycorr/installation/
    INPUT_OBJECT_SPECTRUM=/path/to/data/XXOBJECT.fits
    INPUT_SKY_SPECTRUM=/path/to/data/XXSKY.fits
    OUTPUT_DIR=/path/to/output/

Note: SkyCorr requires absolute paths. The placeholders XXOBJECT and
XXSKY are replaced automatically by the script.

SkySub.py
^^^^^^^^^

Performs simple sky subtraction by scaling the sky spectrum to match
the science spectrum and subtracting.

**Usage**::

    SkySub.py filename

This provides a quick alternative to the DRP sky subtraction for testing
purposes. The scaling factor is determined by minimizing residuals in
regions dominated by sky emission.


Sky Modeling
------------

These tools generate theoretical sky spectra using the ESO Sky Model and
the PALACE airglow model, which can be compared to observed sky spectra
for validation.

EsoSkyObs.py
^^^^^^^^^^^^

Generates a predicted sky spectrum for a given RA, Dec, and time using the
real ESO Sky Model.  Unifies two previously separate approaches
(SkyCalcObs.py, SkyModelObs.py) into one script and one output convention;
both of those scripts have since been retired and removed (see note below).

**Usage**::

    EsoSkyObs.py [-h] [-engine local|remote|auto] [-msol flux] [-out root] [-site lco|paranal] [-pres hPa] [-keep_workdir] ra dec time

**Arguments:**

ra, dec
    Sky position in degrees.

time
    Observation time as a date string, MJD, or JD.

**Options:**

-h
    Print help and exit.

-engine local\|remote\|auto
    ``local`` uses the local ``calcskymodel`` binary only, erroring if it
    isn't available; ``remote`` always uses the ESO SkyCalc web service;
    ``auto`` (default) tries local first, falling back to remote only if
    the local model isn't set up.

-msol flux
    Force a specific 10.7 cm solar radio flux instead of looking one up
    from the historical flux table (``GetSolar.py``).

-out root
    Output filename root; default is ``SkyE_<mjd>_<ra>_<dec>``, the same
    naming convention adopted by PalaceObs.py's ``SkyP_`` (see below).

-site lco\|paranal
    Observatory height/pressure physics used by the model (default
    ``lco``).  See Notes below — this is an approximation for comparing
    against PALACE, not a full site swap.

-pres hPa
    Override the site's default pressure (``lco``: 765, ``paranal``: 744
    — see ``SITE_PRESSURE_HPA``).  Local engine only; ignored by
    ``-engine remote`` (``skycalc_cli`` exposes no separate pressure
    parameter).

-keep_workdir
    Debugging switch.  By default, ``calcskymodel``'s/``skycalc_cli``'s
    inputs and outputs live in a fresh, automatically-deleted temporary
    directory per call — concurrency-safe, but nothing survives a run to
    inspect.  ``-keep_workdir`` instead writes them directly into
    ``./config``, ``./data`` (a symlink), and ``./output`` in the current
    directory and leaves them there — exactly where ``calcskymodel``
    itself looks if you ``cd`` there and run it by hand (confirmed from
    the binary itself, which hardcodes those three names and cannot be
    told to use anything else).  **Not concurrency-safe** — only use it
    for one call at a time.

**Output FITS structure:**

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Description
   * - ``WAVE``
     - wavelength [Angstrom]
   * - ``FLUX``
     - total sky flux, corrected for atmospheric extinction
   * - ``FLUX_LCO``
     - total sky flux, as actually observed on the ground
   * - ``MOON``
     - scattered moonlight component (extinction-corrected)
   * - ``ZODI``
     - zodiacal light component (extinction-corrected)
   * - ``LINES``
     - airglow emission line component (as modeled, **not**
       extinction-corrected -- see Notes)
   * - ``DIFFUSE``
     - diffuse/residual airglow continuum (extinction-corrected)
   * - ``CONT``
     - MOON + ZODI + DIFFUSE
   * - ``trans``
     - atmospheric transmission

The primary header records ``RA``, ``DEC``, ``OBSTIME``, ``MSOLFLUX``,
``ENGINE`` (``local`` or ``remote`` -- which engine actually produced the
file), and ``SITE`` (``lco`` or ``paranal``).

**Requirements:**

The local engine requires the environment variable ``ESO_SKY_MODEL`` to
point to a working installation of the real ESO Sky Model.  The remote
engine requires ``skycalc_cli`` to be pip-installed, a working
setuptools/pkg_resources in the environment (recent ``setuptools``
releases, roughly >=81, dropped ``pkg_resources`` entirely -- pin
``setuptools<81`` if ``skycalc_cli`` fails with
``ModuleNotFoundError: No module named 'pkg_resources'``), and network
access to eso.org.

**Notes:**

Both engines compute the sky as it would actually be observed on the
ground (after atmospheric extinction).  LVM spectra are normally compared
against the above-the-atmosphere equivalent, so FLUX and the
MOON/ZODI/DIFFUSE/CONT components here are corrected for extinction by
dividing by the transmission; the as-observed value is kept separately in
FLUX_LCO.  LINES is **not** currently divided by the transmission --
inherited unchanged from the original SkyModelObs.py/SkyCalcObs.py code
and never re-examined against the same reasoning FLUX got (airglow
originates at ~90 km and should see the same extinction as everything
else on the way down).  This is a real inconsistency between LINES and
the other components worth revisiting.

``-site paranal`` exists for comparing against PALACE, whose own
atmospheric-physics constants are hardcoded to Cerro Paranal
(``h=2.64 km``, ``p=744 hPa``, not overridable via any public PALACE
parameter -- see PalaceObs.py below).  The two engines handle this
differently: the local engine only changes the observatory height/
pressure physics (``SITE_HEIGHT_KM``/``SITE_PRESSURE_HPA``), keeping the
real LCO observing geometry (alt/az, moon phase/separation) -- the same
mixed real-geometry/Paranal-physics approach PALACE itself uses, so this
is the more directly comparable of the two.  ``-site lco``'s local-engine
default pressure is 765 hPa (nominal LCO barometric pressure, not 744) --
``-pres`` overrides either site's default directly.  The remote engine's
``observatory`` parameter drives both the atmosphere physics *and*
skycalc_cli's own internal moon/sun almanac geometry
(``REMOTE_SITE_NAME``), so ``-site paranal`` there also shifts the
modeled sky to Paranal's real geographic location, not just its altitude
-- a real, different kind of approximation than the local engine's.

**See Also:** :doc:`api/EsoSkyObs/index`


PalaceObs.py
^^^^^^^^^^^^

Generates a predicted airglow spectrum for a given RA, Dec, and time using
the PALACE model (Noll et al. 2025, "PALACE v1.0: Paranal Airglow Line And
Continuum Emission model"), in the same physical-unit, homogenized output
convention as EsoSkyObs.py, for direct comparison.  PALACE is an external
dependency, not vendored in this repository -- see Readme.md for
installation.

**Usage**::

    PalaceObs.py [-h] [-srf VALUE] [-species S1,S2,...] [-out ROOT] ra dec obstime

**Arguments:**

ra, dec
    Sky position in degrees.

obstime
    UTC observation time, e.g. ``2023-08-29T03:20:43.668``.

**Options:**

-h
    Print help and exit.

-srf VALUE
    Solar radio flux in sfu; default is looked up from ``data/solar.txt``
    via ``GetSolar.get_flux``.

-species S1,S2
    Comma-separated species to predict individually; default is all nine
    (OH, O2, HO2, FeO, Na, K, O, N, H).

-out ROOT
    Output FITS filename root; default is ``SkyP_<mjd>_<ra>_<dec>``,
    matching EsoSkyObs.py's ``SkyE_`` naming convention.

**Output FITS structure:**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Column
     - Description
   * - ``WAVE``
     - wavelength [Angstrom], PALACE's own native grid truncated to the
       LVM range (0.36-0.98 micron) and generated at fine native
       resolution (R~20000)
   * - ``FLUX``, ``FLUX_LCO``
     - total airglow flux, above the atmosphere / as observed on the
       ground at LCO
   * - ``LINES``, ``LINES_LCO``
     - PALACE's own total line emission, above the atmosphere / at LCO --
       matches ESO's flux_ael
   * - ``DIFFUSE``, ``DIFFUSE_LCO``
     - PALACE's own total continuum emission (3 components: HO2, FeO,
       and O2's own separate continuum), above the atmosphere / at LCO --
       matches ESO's flux_arc
   * - ``OH``, ``O2``, ``HO2``, ``FeO``, ``NA``, ``K``, ``O``, ``N``, ``H``
     - per-species contribution to FLUX (above the atmosphere only)

All flux columns are in erg/s/cm^2/Angstrom for one LVM fiber.  The
primary header records ``RA``, ``DEC``, ``OBSTIME``, ``MSOLFLUX``, and
``ENGINE='palace'``.

**Requirements:**

The PALACE Python package (external dependency; see Readme.md).

**Notes:**

Both the local ESO Sky Model and PALACE compute the sky as it would
actually be observed on the ground; FLUX/LINES/DIFFUSE here are the
above-the-atmosphere value (PALACE's own ``isatm=False``), with the
as-observed value in the ``_LCO`` columns (``isatm=True``).  Since ESO's
own LINES column is not currently divided by the transmission (see
EsoSkyObs.py Notes), today PALACE's ``LINES_LCO`` (not ``LINES``) matches
ESO's LINES like-for-like, while PALACE's ``DIFFUSE`` (not
``DIFFUSE_LCO``) matches ESO's DIFFUSE (which *is* divided by the
transmission) -- an inconsistency inherited from ESO's side, not PALACE's.

LINES/DIFFUSE are built from PALACE's own internal line-table/continuum-
table split (calling ``readdata``/``calcscalfac``/``scalelines``/
``scalecont``/``corratmlines``/``corratmcont``/``calclinspec``/
``calccontspec``/``convolvelsf`` directly, stopping short of the final
line+continuum sum), not by summing whole species: PALACE's public
``species=`` selector operates on whole species, and O2 has its own
separate continuum component (``palace_cont.fits``' own header:
``NCONT=3``, ``CHEM3='O2'``, ``VARID3='O2Ac'``) in addition to its line
emission, so a per-species line/continuum classification would silently
fold that continuum into LINES.

``FLUX_LCO`` can come out very slightly *brighter* than ``FLUX`` in
places -- not a bug: PALACE's own ``isatm=True`` scattering correction
includes a van-Rhijn/in-scattering term for the extended airglow layer
that can slightly exceed unit transmission at some zenith angles/
wavelengths (light scattered into the line of sight from the rest of the
sky outweighing light scattered out of it).

If your local ``$ESO_SKY_MODEL`` install's airglow continuum data file
has been patched to use PALACE's own continuum shape (check
``sm_filenames.dat``'s ``acontname`` entry), verify its absolute scaling
carefully before comparing ESO's local-engine DIFFUSE against PALACE's --
a stale legacy scale factor left over from the original (much cruder)
ESO continuum template can inflate the local engine's DIFFUSE relative to
PALACE's by a large factor while leaving the wavelength *shape* similar
(scale and shape errors look very different, which is how this was first
noticed). Diagnosed and fixed in this way on 2026-07-10 for this
project's own ``$ESO_SKY_MODEL`` install: the patched ``palace_cont.dat``
carried its predecessor ``airglow_cont.dat``'s ``scale`` header value
unchanged even though the tabulated continuum was no longer normalized to
1.0 at the model's 0.543 micron reference wavelength, and the corrected
absolute level (informed by comparison with real LVM data, not the
ESO/PALACE ratio alone) needed a further empirical adjustment beyond that
units fix -- see the ``scale`` header value and the backup copies kept
alongside ``palace_cont.dat`` in that installation's ``data/`` directory
for the full history.

**See Also:** :doc:`api/PalaceObs/index`


.. note::
   SkyCalcObs.py and SkyModelObs.py (the two previously separate
   approaches EsoSkyObs.py unified above) were retired and removed on
   2026-07-11: everything that used them (EsoSkyObs.py itself,
   PalaceObs.py, SkySepESO.py, MakeMoonBase.py) has been migrated to
   EsoSkyObs.py/its ``get_info_las_campanas``.  Use EsoSkyObs.py's
   ``-engine local``/``-engine remote`` in their place.


PALACE-based Sky Line Masking and Sky Spectrum Collection
----------------------------------------------------------

These tools support the development of improved sky subtraction by identifying
sky-line-free wavelength windows for continuum fitting and by assembling
stacked sky spectra from the LVM sky telescopes.  Together they are intended
to characterise the sky background well enough to constrain physical models
of the airglow emission.

XSkySepIvan.py
^^^^^^^^^^^^^^

Decomposes LVM sky spectra into physical emission components using the
PALACE (Paranal Airglow Line And Continuum Emission, Noll et al. 2024) line
model combined with a B-spline Moon/zodiacal continuum.  This script is a
wrapper for the ``SkyDecomp`` class written by Ivan Katkov, vendored into
``py_progs/sky_decomp/fit.py`` (260709; previously an external dependency on
the ``lvmsky`` repository).  The vendored copy has no dependency on the
``lvmdrp`` package — the one function PALACE used from it (a rebin+convolve
utility, no real DRP logic) is inlined directly; ``clarabel`` (the QP solver
used for the fits themselves) remains a normal, separate dependency that
must still be installed in whichever environment runs this.

The decomposition solves a non-negative quadratic programme (Clarabel solver)
for the amplitudes of six component families:

- **OH** — 402 groups of hydroxyl vibrational-rotational lines (HITRAN data)
- **Moon** — B-spline envelope multiplied by a rebinned solar spectrum (captures Moon reflected light and zodiacal continuum)
- **Diffuse** — PALACE diffuse continuum: HO\ :sub:`2`, FeO, O\ :sub:`2`\Ac
- **Atom** — atomic airglow: NaI, KI, [NI], OI (green and red)
- **ORC** — OI recombination multiplets at 7774 and 8446 Å
- **O2** — molecular oxygen A-band (~8650 Å); rotational temperature fitted

The script supports two input modes:

*Sky-file mode* (files produced by ``GetSky_from_CFrame_sum.py``): each row
of the FLUX array is a sky-telescope spectrum of the same field from a
different exposure.  Noise is estimated from pixel-to-pixel differences since
no IVAR extension is present.

*XCframe mode*: standard LVM summary file with FLUX/SKY_EAST/SKY_WEST
extensions.  IVAR is read from the file if present, otherwise estimated.

**Usage**::

    # Process all spectra in a Sky file
    XSkySepIvan.py Sky_WHAM_south_08.fits

    # Process specific rows of a Sky file
    XSkySepIvan.py Sky_WHAM_south_08.fits 0 5 10

    # Process all rows of an XCframe extension
    XSkySepIvan.py XCframe_file.fits SKY_EAST

    # Process specific rows of an XCframe extension
    XSkySepIvan.py XCframe_file.fits SKY_EAST 0 300 600

    # Process every 1000th row of an XCframe extension
    XSkySepIvan.py XCframe_file.fits SKY_EAST -delta 1000

**Key options:**

-delta N
    Process every N-th row (0, N, 2N, ...) instead of all rows.
    Ignored if explicit row numbers are given.

-lsf FWHM
    Fixed LSF FWHM in Angstroms (default 1.3 Å).

-refits N
    Number of iterative per-channel LSF kernel refits (default 0).

-out outroot
    Set output filename root.

**Output:**

A FITS file named ``<stem>_<ext>_ivan.fits`` (XCframe mode) or
``<stem>_ivan.fits`` (Sky-file mode) with the following extensions:

*Spectral image extensions* (float32, N_obs × N_pix):

- ``WAVE`` — wavelength array (Å)
- ``FLUX`` — input sky spectra
- ``LINES`` — total emission (OH + ATOM + ORC + O2)
- ``CONT`` — total continuum (MOON + DIFFUSE)
- ``OH`` — OH band component
- ``ATOM`` — atomic airglow (NaI, KI, [NI], OI)
- ``ORC`` — OI recombination multiplets
- ``O2`` — O2 A-band component
- ``MOON`` — Moon/zodiacal continuum
- ``DIFFUSE`` — diffuse airglow continuum
- ``RESID`` — fit residuals (FLUX − LINES − CONT)

*Coefficient tables:*

- ``COEF`` — BinTable with 442 named columns (one per design-matrix
  entry), storing the physical-unit coefficient for each spectrum.
  Column names match the ``name`` column of COEF_META:
  ``OH_v{v}_N{nn}_F{f}``, ``Moon_bs{nn}``, ``HO2``, ``FeO``,
  ``O2Ac``, ``NaI``, ``KI``, ``NI_forb``, ``OI_5577``,
  ``OI_6300``, ``OI_7774``, ``OI_8446``, ``O2_Aband``.

- ``COEF_META`` — BinTable with 442 rows, one per design-matrix
  entry, carrying physical identity (``name``, ``component``,
  ``v_upper``, ``N_upper``, ``F_upper``, ``wave_peak``) and
  cross-spectrum distribution statistics (``coef_median``,
  ``coef_mean``, ``coef_nmad``, ``coef_min``, ``coef_max``,
  ``coef_skew``).

- ``DRP_ALL`` — observation metadata plus 20 compact coefficient
  summary columns (physical flux units, same names as COEF_META):
  ``OH_v3``…``OH_v10`` (total OH flux per vibrational band),
  ``NaI``, ``KI``, ``NI_forb``, ``OI_5577``, ``OI_6300``,
  ``OI_7774``, ``OI_8446``, ``O2_Aband``, ``HO2``, ``FeO``,
  ``O2Ac``, ``Moon_med``.

XSkySepIvan_eval.py
^^^^^^^^^^^^^^^^^^^

Evaluates the quality of a PALACE sky decomposition produced by
``XSkySepIvan.py`` by creating an interactive three-panel HTML plot over a
chosen wavelength window.

Each panel shows the median spectrum (coloured line) surrounded by a shaded
10th/90th percentile band, with a random sample of individual spectra overlaid
in grey so outliers and systematic trends are immediately visible:

- **Panel 1 (Flux)** — input sky spectra
- **Panel 2 (Residual)** — fit residuals (FLUX − model)
- **Panel 3 (Lines)** — total fitted emission (OH + Atom + ORC + O2)
- **Panel 4 (Continuum)** — smooth Moon/zodiacal + diffuse continuum, log scale

Each panel carries its own legend.  The individual-spectrum alpha is set
automatically as ``min(0.7, 3/√N)`` so traces remain legible regardless of
how many spectra are in the file.

**Usage**::

    # Full wavelength range
    XSkySepIvan_eval.py Sky_WHAM_south_08_ivan.fits

    # Specific window
    XSkySepIvan_eval.py Sky_WHAM_south_08_ivan.fits 9300 9600

    # Limit random sample overlay
    XSkySepIvan_eval.py Sky_WHAM_south_08_ivan.fits 6000 7000 -num 10

**Output:**

An HTML file named ``<stem>_<wmin>_<wmax>.html`` (e.g.
``Sky_WHAM_south_08_ivan_9300_9600.html``) that can be opened in any
browser for interactive zoom, pan, and hover inspection.

palace_make_mask.py
^^^^^^^^^^^^^^^^^^^

Builds a sky-line contamination mask across the full LVM wavelength range
(3600-9800 Å) using the PALACE (Paranal Airglow Line And Continuum Emission,
Noll et al. 2024) sky emission model.  Four components are rendered onto the
LVM wavelength grid: OH vibrational-rotational bands, OI recombination lines,
atomic forbidden/permitted lines (NaI, KI, [NI], OI), and the O2 A-band.
Each component is normalised to its own peak before summing so that no single
family dominates the mask.  The combined model is scaled to the observed sky
spectrum via a least-squares fit to bright OH pixels in the Z arm, and pixels
where the predicted contamination exceeds a user-specified threshold are
flagged as unusable for continuum fitting.

**Usage**::

    palace_make_mask.py fits_file [palace_dir] [--threshold T] [--plot] ...

**Arguments:**

fits_file
    LVM XCframe FITS file (provides WAVE, sky spectrum, and LSF).

palace_dir
    Path to the ``palace/PMD`` directory containing the PALACE data files.
    Optional; defaults to the vendored copy at
    ``data/palace_ref/palace/PMD`` (computed relative to
    ``palace_make_mask.py``'s own location), so it only needs to be given
    to point at a different PMD installation.

**Key options:**

--threshold T
    Contamination threshold in FACTOR-scaled flux units.  Lower values give a
    stricter mask.  Default is 0.01 (= 1×10⁻¹⁶ erg s⁻¹ cm⁻² Å⁻¹ with the
    default FACTOR of 10¹⁴).  The tradeoff between mask strictness and the
    number of clean pixels available for continuum fitting is the primary
    tuning parameter.  This threshold constrains only the modelled PALACE
    sky *line* flux (OH + OI + atomic + O2, summed and scaled to the
    observed sky) — it says nothing about the real observed continuum
    level, S/N, or local brightness, and applies as one fixed absolute
    flux value across the whole 3600–9800 Å range regardless of arm.

--plot
    Display the diagnostic plot interactively (it is always saved as a PNG).

--line-output PATH
    Output path for the strong-sky-line list (default ``<stem>_lines.txt``).
    See Output below.

**Console report:**

For each of the B, R, and Z arms (and an ``ALL`` row for the full range),
prints total pixel count, clean pixel count and percentage, number of
distinct clean windows, and the threshold expressed as a physical flux
value (``threshold / factor``, erg s⁻¹ cm⁻² Å⁻¹) — the same value in every
row since one global threshold applies everywhere, shown per arm so it's
visible alongside each arm's clean fraction. As a reference point, the Z
(NIR) arm is dominated by the OH airglow forest densely enough that even a
threshold 20× looser than the 0.01 default (0.2, i.e. 2×10⁻¹⁵ erg s⁻¹
cm⁻² Å⁻¹) still leaves under half the Z arm flagged clean.

**Output:**

A FITS file (``<stem>_mask.fits``) containing:

- ``WAVE`` — wavelength array (Å)
- ``SKY`` — median observed sky spectrum (FACTOR-scaled)
- ``CONTINUUM`` — sky spectrum with contaminated pixels set to NaN
- ``MASK`` — boolean mask (1 = clean, 0 = contaminated)

A PNG diagnostic plot (``<stem>_mask.png``) showing all three
spectrograph arms on a log flux scale with the PALACE model, threshold line,
the full sky spectrum, and the clean continuum pixels highlighted.

An ascii table (``<stem>_lines.txt``, ``--line-output`` to change the path)
listing strong sky lines as *named positions* rather than a pixel mask.  For
each PALACE line group already used to build the mask (OH grouped by
(v_upper, N_upper, F_upper), OI recombination grouped by reffeat, atomic
lines grouped by feat -- e.g. NaI0589, OI0558, KI0770 -- and the O2 A-band
as one group), the group's single brightest transition is kept if the
combined, scaled contamination model exceeds the same ``--threshold`` used
for the mask.  Columns: ``Wave_air``, ``LineID``, ``Component``, ``Ampl``,
written as ``ascii.fixed_width_two_line``.  Because it uses the same
threshold as the mask, in the Z arm this list can run to several hundred
entries (mostly OH); it is intended as a per-line reference for judging
whether a spectral feature might be a sky-subtraction residual, not as a
short curated list. See :doc:`spectrum_plots` (PlotSpecI.py's
``-sky_lines``/``-sky_lines_file``) for how this file is meant to be used
-- it overlays these positions as unlabeled tick marks, distinct from the
labeled scientific line list.

GetSky_from_CFrame_sum.py
^^^^^^^^^^^^^^^^^^^^^^^^^

Inventories which sky fields have been observed in an XCframe summary file
and, optionally, extracts all sky spectra for a chosen field into a single
FITS file for stacking or modelling.

LVM records both an east (SKY_EAST) and west (SKY_WEST) sky telescope
pointing for each science exposure.  The DRP_ALL table inside an XCframe
summary file lists the field names in the ``skye_name`` and ``skyw_name``
columns.  This program has two operating modes:

*Summary mode* (no source name given): counts how many times each sky field
appears across both sky telescopes, prints the top N entries sorted by
observation count, and writes the full count table to
``sky_summary_<stem>.fits``.

*Extraction mode* (source name given): selects all rows where the field name
matches, extracts the corresponding spectra from SKY_EAST and SKY_WEST,
merges paired metadata columns (``skye_ra``/``skyw_ra`` → ``ra``, etc.) into
single columns for the relevant telescope, removes science-pointing columns,
and writes the result to ``Sky_<source_name>.fits``.

**Usage**::

    # Inventory sky fields
    GetSky_from_CFrame_sum.py fits_file

    # Extract spectra for one field
    GetSky_from_CFrame_sum.py fits_file source_name [--output PATH]

**Output (extraction mode):**

A FITS file containing:

- ``WAVE`` — wavelength array (Å)
- ``FLUX`` — sky spectra, shape (N_obs, N_pix), one row per observation
- ``DRP_ALL`` — metadata table with merged sky columns and a ``tel`` column
  indicating which sky telescope (SKY_EAST or SKY_WEST) each row came from

GetSkyCont.py
^^^^^^^^^^^^^

Fits a smooth two-component B-spline continuum to LVM sky spectra, separating
the continuum into a MOON/zodiacal component (B-splines modulated by a solar
spectrum) and a DIFFUSE component (plain B-splines for airglow continuum).
The fit is performed only on pixels flagged as line-free by a
``palace_make_mask.py`` mask; the result is evaluated at all wavelengths so
the continuum interpolates smoothly across emission-line regions.

The solar spectrum is read from a high-resolution reference file
(Meftah et al., LATMOS); its Fraunhofer absorption structure is kept fixed
while the B-spline envelope adjusts the colour and amplitude of the
Moon-reflected and zodiacal-light contribution.

**Usage**::

    # Sky-file mode (output of GetSky_from_CFrame_sum.py):
    GetSkyCont.py sky_file.fits -mask mask.fits [row_no ...] [-delta N]

    # XCframe / XSFrame summary file mode:
    GetSkyCont.py xframe.fits ext -mask mask.fits [row_no ...] [-delta N]

**Arguments:**

sky_file.fits
    Sky_<name>.fits produced by ``GetSky_from_CFrame_sum.py``.

xframe.fits, ext
    XCframe or XSFrame summary FITS file and the extension to read
    (e.g. ``SKY_EAST``, ``SKY_WEST``, ``FLUX``).

row_no
    Zero or more 0-based row indices.  If omitted all rows are processed
    (subject to ``-delta``).

**Key options:**

-mask file
    (required) mask FITS file from ``palace_make_mask.py``
    (MASK extension: 1 = clean, 0 = line-affected).

-delta N
    Process every N-th row (0, N, 2N, ...) instead of all rows.

-kstep N
    B-spline knot spacing in Angstroms (default 100 Å, giving ~66 basis
    functions per component).

-out outroot
    Set output filename root.

**Output:**

A FITS file (``<stem>_cont.fits``) containing:

- ``WAVE`` — wavelength array (Å)
- ``FLUX`` — input sky spectra (N_obs × N_pix)
- ``CONT`` — total continuum = MOON + DIFFUSE
- ``MOON`` — B-spline × solar component (Moon/zodiacal light)
- ``DIFFUSE`` — plain B-spline component (diffuse airglow continuum)
- ``RESID`` — residual FLUX − CONT (line emission isolated from continuum)
- ``MASK`` — boolean clean-pixel mask from the input palace_mask file
- ``DRP_ALL`` — observation metadata table

GetSkyCont_eval.py
^^^^^^^^^^^^^^^^^^

Evaluates the continuum fit produced by ``GetSkyCont.py`` by writing a single
HTML file containing three interactive Plotly figures.  Also updates the
DRP_ALL table in the input FITS file with per-spectrum fit-quality statistics.

**Figure 1 — four-panel spectral overview:**

- **Panel 1 (Flux + Continuum, log)** — observed sky spectra with the fitted
  total continuum median overlaid in red.  A green line near the bottom of
  the panel marks the wavelengths included in the fit (gaps where sky lines
  were masked); grey vertical bands shade the excluded regions.
- **Panel 2 (Total Continuum, log)** — the CONT band (10th–90th percentile
  and median) showing the overall level and smoothness of the B-spline model.
- **Panel 3 (Components, log)** — MOON (red) and DIFFUSE (orange) components
  on the same y-axis as Panel 2, with the total CONT median in purple.  The
  relative amplitude of the components shows how much of the continuum is
  Moon/zodiacal versus diffuse airglow.
- **Panel 4 (Residual, linear)** — FLUX − CONT; should be near zero in clean
  regions and show sky-line emission in the masked regions.

**Figure 2 — per-arm residual histograms:**

Three panels (Blue 3600–5900 Å, Red 5900–7600 Å, NIR 7600–9800 Å) each show
the distribution of all residual flux values — across every spectrum and every
clean (unmasked) pixel in that arm.  The x-axis is residual flux; the y-axis
is count N.  A Gaussian with center = median and σ = NMAD is overlaid in
black.  An annotation box reports N, the median, NMAD, skewness, and the
10th/90th percentiles.  The histogram range is clipped to median ± 5·NMAD to
suppress extreme outliers.

**Figure 3 — per-spectrum fit quality:**

- **Row 1** — three scatter plots (Blue, Red, NIR) of per-spectrum median
  residual (x-axis) vs NMAD (y-axis).  Each point is one spectrum; hovering
  shows the spectrum row index and its median and NMAD values.  Axis limits
  are clipped to the 2nd–98th percentile range so extreme outliers do not
  compress the scale for the bulk of the spectra; any off-scale spectra are
  noted in an annotation box.  A dashed vertical line marks x = 0 (ideal
  median); a dotted horizontal line marks the ensemble NMAD for reference.
- **Row 2** — NMAD vs original spectrum number (log y-scale) for all three
  arms on one panel, so that clusters of temporally adjacent poor fits are
  immediately visible.  Dotted horizontal lines mark the ensemble NMAD per arm.

**Per-spectrum statistics written to DRP_ALL:**

For each arm (``blue``, ``red``, ``nir``) and for all arms combined
(``all``), four float32 columns are added or updated in the DRP_ALL table:

+---------------------+---------------------------------------------+
| Column              | Description                                 |
+=====================+=============================================+
| ``resid_med_<arm>`` | median residual in clean pixels             |
+---------------------+---------------------------------------------+
| ``resid_nmad_<arm>``| NMAD of residuals in clean pixels           |
+---------------------+---------------------------------------------+
| ``resid_rms_<arm>`` | RMS of residuals in clean pixels            |
+---------------------+---------------------------------------------+
| ``resid_skew_<arm>``| skewness of residuals in clean pixels       |
+---------------------+---------------------------------------------+

**Usage**::

    # Full wavelength range
    GetSkyCont_eval.py Sky_WHAM_south_08_cont.fits

    # Specific window
    GetSkyCont_eval.py Sky_WHAM_south_08_cont.fits 6000 7000

    # Limit the random sample overlay
    GetSkyCont_eval.py Sky_WHAM_south_08_cont.fits 6000 7000 -num 10

**Output:**

An HTML file named ``<stem>_<wmin>_<wmax>.html`` containing all three
figures, viewable in any browser with interactive zoom, pan, and hover
inspection.  The input FITS file is updated in place with the DRP_ALL
statistics columns.


XCframe Sky Subtraction — Alternative Methods and Evaluation
-------------------------------------------------------------

This group of scripts performs per-spectrum sky subtraction directly on
XCframe summary FITS files, providing five independent algorithms that can be
run and compared side by side using ``SkySub_eval.py``.  All five subtraction
scripts produce output FITS files with WAVE, FLUX (sky-subtracted), SKY, and
DRP_ALL extensions (SkySepESO adds MOON/ZODI/DIFFUSE; see below).  DRP_ALL
carries a QA_FLAGS column recording per-row quality issues.

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
  SummarizeSciSky.py below).
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
  "Continuum Separation" figures below for how they are used.

SkySubOrig.py
^^^^^^^^^^^^^

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
^^^^^^^^^^^^

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
^^^^^^^^^^^^^

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
^^^^^^^^^^^^^

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
^^^^^^^^^^^^^

Perform per-spectrum sky subtraction using
``sky_decomp.lsf_surface_iterative.SkyDecompLSFSurfaceIterative`` (the
same physically-motivated, per-row-LSF-refined decomposition
``DecomposeCleanSky.py`` uses — see "Nebular-Line-Based Method
Evaluation" below) for the continuum/line separation, with the near/far
recipe and bisection scale search otherwise identical to SkySubDev1.py.

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
   sky-telescope spectra: the nebular-leak validation work below found
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
^^^^^^^^^^^^

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
below expect.

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
^^^^^^^^^^^^

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
per call in a private temporary directory -- see ``EsoSkyObs.py`` above.

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
environment variable (see ``EsoSkyObs.py`` above); without it, every row
falls back to the SkyCalc web service (slower, needs network access).

**See Also:** :doc:`api/SkySepESO/index`


SkySepPalace.py
^^^^^^^^^^^^^^^

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
see PalaceObs.py above.  Since each species/MOON amplitude is entirely
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
above and Readme.md).

**See Also:** :doc:`api/SkySepPalace/index`


SkySub_eval.py
^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
looks equally good either way.  For that question, see "Nebular-Line-
Based Method Evaluation" below.

sky_residual_eval.py
^^^^^^^^^^^^^^^^^^^^

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


Nebular-Line-Based Method Evaluation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Every generic sky-residual metric above (``SkySub_eval.py``,
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
every method above (this evaluation family included, for the far/near
sky-line component) implicitly depends on.

The tools split into two lines of investigation sharing a common line
catalog and velocity convention:

- **Does the sky itself leak nebular flux?** (``DecomposeCleanSky.py`` →
  ``sky_nebular_leak_eval.py`` → ``PlotNebularLeak.py``) — fits a
  nebula-free "clean sky" model to FLUX/SKY_EAST/SKY_WEST by masking out
  nebular-line windows before the fit, then measures whether real flux
  still leaked into those windows anyway.
- **Which SkySub* method best recovers real nebular flux?**
  (``SkySubRun.py`` → ``SkySubNebEval.py`` → ``PlotSkySubNebEval.py``) —
  fits the same nebular-line catalog directly on each method's
  sky-*subtracted* output and compares the measured ratios/repeat
  scatter across methods.

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
    PlotNebularLeak.py                  PlotSkySubNebEval.py
     (per-line before/after,             (pointing table, spectrum
      shared-systematic template          overview, ratio summary
      correction panels)                  table, per-line-group panels
                                           -- one column per method)

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
^^^^^^^^^^^^^^^^^^^^

Decomposes one or more of FLUX (science fiber), SKY_EAST, and SKY_WEST
from an XCframe summary file into a nebula-free "clean sky" model, using
``sky_decomp.lsf_surface_iterative.SkyDecompLSFSurfaceIterative`` with the
known nebular emission lines excluded from the fit ("mask-and-wrap").
Lives in ``py_dev/`` (not Sphinx-API-documented — see :doc:`sky_model_landscape`
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
^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^

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
    from that row's own DRP_ALL ``Redshift`` (see "Nebular-Line-Based
    Method Evaluation" above). Only pass these to force one velocity
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
^^^^^^^^^^^^^^^^^^^^

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


Science-Fiber-Based Sky Estimation
------------------------------------

Unlike the XCframe methods above (which use the dedicated SKY_EAST/SKY_WEST
sky telescopes), these two scripts estimate a sky spectrum directly from the
science IFU itself.  At a typical LVM pointing most of the 1801 science
fibers see mostly sky rather than an astronomical source; ranking fibers by
sky-line-free continuum flux and averaging the faintest fibers gives a sky
proxy without needing a sky telescope pointing at all.  No scale factor is
applied to emission lines here, so the method is best suited to fields
where a genuinely sky-dominated tail of faint fibers exists (e.g. diffuse
or extended sources, not compact point sources filling the IFU).

Both scripts write output with WAVE/SCI/SKY/FLUX/DRP_ALL extensions
compatible with ``SkySub_eval.py`` (FLUX = sky-subtracted, SKY = sky
model), so they can be evaluated and compared alongside the four XCframe
methods above.

SkySubSci.py
^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^

Drpall-driven, remote-friendly version of ``SkySubSci.py``: runs the same
science-fiber sky estimate over a range of exposure numbers selected from a
drpall table (as ``SummarizeCframe.py`` does) instead of an explicit file
list, so it can be run unattended over many exposures (e.g. at Utah).

**Usage**::

    SummarizeSciSky.py [-emin 900] [-ver 1.2.1] [-drp_all FILE]
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
    DRP version, used to locate ``drpall-VER.fits`` (default 1.2.1).

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
are re-implemented locally here (with an added ``-drp_all`` override)
rather than imported from ``SumCframe.py``, so this script has no
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


Typical Workflows
-----------------

Evaluating DRP Sky Subtraction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Obtain SFrame files from the DRP
2. Run eval_sky.py to visualize residuals::

       eval_sky.py lvmSFrame-00012345.fits

3. Look for systematic residuals at sky line wavelengths

Recovering an Exposure with a Bad Sky Pointing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Identify a clean exposure whose sky telescope pointed well away from
   the Moon (e.g. by checking ``SKY SKYW_MOON_SEP``/``SKY SKYE_MOON_SEP``
   in each candidate CFrame's header)
2. Substitute that telescope's fiber data into the compromised
   exposure's CFrame::

       SubstituteSky.py lvmCFrame-00014964.fits SkyW \
                         lvmCFrame-00014771.fits SkyW

3. Rerun the DRP's sky subtraction on the result::

       RunSky.py lvmCFrame-00014964.sky_subst.fits

4. Evaluate the corrected lvmSFrame as usual (eval_sky.py, sky_plot.py)

Testing Alternative Sky Subtraction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Start with CFrame files (before DRP sky subtraction)
2. Prepare files for SkyCorr::

       Prep4SkyCorr.py lvmCFrame-00012345.fits

3. Run SkyCorr::

       RunSkyCorr.py sci_mean.fits sky_mean.fits

4. Compare results with DRP sky subtraction

Comparing with Sky Models
^^^^^^^^^^^^^^^^^^^^^^^^^

1. Generate a theoretical sky for the observation::

       EsoSkyObs.py 81.5 -66.0 60000.5

2. Compare with observed sky from SKY_EAST or SKY_WEST telescopes
3. Identify discrepancies that may indicate calibration issues

Fitting and Evaluating the Sky Continuum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Build a palace line mask for the field::

       palace_make_mask.py XCframe_file.fits

2. Collect sky spectra from repeated observations of the field::

       GetSky_from_CFrame_sum.py XCframe_file.fits Sky_WHAM_south_08

3. Fit the two-component B-spline continuum::

       GetSkyCont.py Sky_WHAM_south_08.fits -mask XCframe_file_mask.fits

4. Evaluate the fit interactively::

       GetSkyCont_eval.py skycont_Sky_WHAM_south_08.fits

Comparing Methods by Nebular-Line Recovery
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

3. For a batch statistic across many exposures/repeat groups instead of
   one exposure, run ``SkySubNebEval.py`` on the same file list and
   inspect its printed doublet-scatter comparison (or the ``_repeats.fits``
   table directly)

4. To check whether SKY_EAST/SKY_WEST leak nebular flux in the first
   place (an assumption every method above depends on)::

       DecomposeCleanSky.py -ext SKY_EAST,SKY_WEST XCframe_file.fits
       sky_nebular_leak_eval.py CleanSky_<expnum>.fits
       PlotNebularLeak.py CleanSky_<expnum>.fits


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


See Also
--------

- :doc:`summarize` - Tools for evaluating sky subtraction across many exposures
- :doc:`api/eval_sky/index` - API documentation
- :doc:`api/SubstituteSky/index` - API documentation
- :doc:`api/RunSky/index` - API documentation
- :doc:`api/Prep4SkyCorr/index` - API documentation
- :doc:`api/RunSkyCorr/index` - API documentation
- :doc:`api/SkySub/index` - API documentation
- :doc:`api/EsoSkyObs/index` - API documentation
- :doc:`api/PalaceObs/index` - API documentation
- :doc:`api/palace_make_mask/index` - API documentation
- :doc:`api/GetSky_from_CFrame_sum/index` - API documentation
- :doc:`api/XSkySepIvan/index` - API documentation
- :doc:`api/XSkySepIvan_eval/index` - API documentation
- :doc:`api/GetSkyCont/index` - API documentation
- :doc:`api/GetSkyCont_eval/index` - API documentation
- :doc:`api/SkySubOrig/index` - API documentation
- :doc:`api/SkySubDrp/index` - API documentation
- :doc:`api/SkySubDev1/index` - API documentation
- :doc:`api/SkySubDev2/index` - API documentation
- :doc:`api/SkySubDev3/index` - API documentation
- :doc:`api/SkySubRun/index` - API documentation
- :doc:`api/SkySepESO/index` - API documentation
- :doc:`api/SkySepPalace/index` - API documentation
- :doc:`api/SkySub_eval/index` - API documentation
- :doc:`api/sky_residual_eval/index` - API documentation
- :doc:`api/sky_nebular_leak_eval/index` - API documentation
- :doc:`api/PlotNebularLeak/index` - API documentation
- :doc:`api/SkySubNebEval/index` - API documentation
- :doc:`api/PlotSkySubNebEval/index` - API documentation
- :doc:`api/SkySubSci/index` - API documentation
- :doc:`api/SummarizeSciSky/index` - API documentation
- :doc:`summarize` - SummarizeCframe.py, whose drpall selection logic SummarizeSciSky.py mirrors
- :doc:`spectral_fitting_local` - ``sky_gaussfit.py``'s NEBULAR_LINES/resolve_nebular_lines, shared by the nebular-line evaluation tools above
- :doc:`data_quality` - ``QualSFrame.py``'s pointing table, reused by PlotSkySubNebEval.py
