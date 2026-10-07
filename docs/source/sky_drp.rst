Evaluating and Repairing DRP Sky Subtraction
============================================

This page is part of :doc:`sky_subtraction`.

This page covers checking the quality of the
DRP's own sky subtraction in an SFrame, and recovering an exposure whose
sky subtraction is compromised by a badly placed sky telescope (for
example one too close to the Moon).


eval_sky.py
-----------

Creates diagnostic plots to evaluate sky subtraction quality in SFrame
files (sky-subtracted data from the DRP).

**Usage**::

    eval_sky.py filename1 filename2 ...

Plots median, min, max, and average spectra across science fibers to
identify residual sky features.

sky_plot.py
-----------

Similar to eval_sky.py, creates plots to evaluate sky subtraction
quality with additional visualization options.

**Usage**::

    sky_plot.py filename1 filename2 ...


Correcting and Rerunning DRP Sky Subtraction
--------------------------------------------

These two scripts work together to recover an exposure whose DRP sky
subtraction is compromised by a bad sky-telescope pointing (e.g. SkyW
too close to the Moon): substitute in a clean sky telescope's data --
the other telescope in the same exposure, or the same telescope from a
different exposure -- then rerun the production DRP's own
sky-subtraction routine (or any other method) on the result.


SubstituteSky.py
-----------------

Replaces one sky telescope's fiber data (FLUX/IVAR/MASK/LSF), its
SKY_EAST/SKY_WEST sky-model extension, and associated header metadata
in an lvmCFrame with data from the same or the other sky telescope, in
the same or a different lvmCFrame, producing a new self-consistent
lvmCFrame usable by any sky-subtraction method.

**Example**::

    # 14964's SkyW was 4.5 deg from the Moon; use its own SkyE instead
    SubstituteSky.py lvmCFrame-00014964.fits SkyW lvmCFrame-00014964.fits SkyE

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
extensions, which the current production method ignores.  Many other
tools (``SummarizeCframe.py``, ``GetTelData.py``, ``lvm_skyfit.py``,
``Prep4SkyCorr.py``, ``QualCFrame.py``, ...) do read SKY_EAST/SKY_WEST,
so both are replaced.

For the same telescope on both sides, the FLUX, IVAR, MASK, and LSF
rows are matched by fiberid (fiber assignment is fixed hardware,
identical across exposures).  For different telescopes
(``target_tel`` != ``source_tel``) there is no such correspondence, so
a warning is printed and the target rows are filled in order from the
source telescope's good (``fibstatus==0``) fibers, reused cyclically
if there are fewer; since the DRP and the XCframe summaries average
over a telescope's fibers, the ordering does not matter.  The
target telescope's SKY_EAST/SKY_WEST extension (and ``_IVAR``) is
replaced by the source telescope's.

Every PRIMARY header keyword tied to the source telescope (pointing,
altitude, airmass, guider frames, sky-field name, heliocentric
velocity, moon/shadow geometry, ecliptic coordinates, etc.) is copied
too, renamed to the target telescope's own keyword names.  The
SKYEW/SKYWW combination weights are recomputed from the updated
pointings with ``skyMethod.combine_skies``'s own formula (inverse
angular distance to the science field, normalized).  ``SKY SCI_SKYW_SEP`` (or the SkyE equivalent) is relative
to the *target's* own science pointing, so it is recomputed from the
target's real SCIRA/SCIDEC and the newly-copied sky position
(``lvmdrp.core.sky.ang_distance``) rather than copied as-is.  New
``SKY SUBST_*`` provenance keywords record what was substituted and
from where.  ``target_cframe``/``source_cframe`` are never modified.

**Output:**

An lvmCFrame FITS file, structurally identical to the input, with the
named telescope's FLUX/IVAR/MASK/LSF rows, SKY_EAST/SKY_WEST extension,
and header block replaced.

**See Also:** :doc:`api/SubstituteSky/index`


RunSky.py
---------

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

1. Identify the bad telescope and a clean replacement, by checking
   ``SKY SKYW_MOON_SEP``/``SKY SKYE_MOON_SEP`` (and
   ``SKY SCI_SKYW_SEP``/``SKY SCI_SKYE_SEP``) in the CFrame header.  The
   replacement can be:

   - the *other* sky telescope in the same exposure, if it is clean --
     the simplest option, since it was observed at the same time and
     airmass; or
   - the same telescope from a different exposure taken near in time.

2. Substitute it into the compromised exposure's CFrame.  The result is
   a new, self-consistent CFrame (fiber rows, SKY_EAST/SKY_WEST, header,
   SKYEW/SKYWW all updated) that any sky-subtraction method can use::

       # 14964's SkyW was 4.5 deg from the Moon; its SkyE was clean
       SubstituteSky.py lvmCFrame-00014964.fits SkyW \
                         lvmCFrame-00014964.fits SkyE

       # or: same telescope from another exposure
       SubstituteSky.py lvmCFrame-00014964.fits SkyW \
                         lvmCFrame-00014771.fits SkyW

3. Rerun the DRP's sky subtraction on the result (or run any other
   method, e.g. via ``SummarizeCframe.py`` + ``SkySubRun.py``)::

       RunSky.py lvmCFrame-00014964.sky_subst.fits

   ``RunSky.py`` also writes an ancillary sky table under
   ``$SAS_BASE_DIR`` (``.../ancillary/lvm-skytable-brz-<expnum>.fits``).

4. Compare the corrected exposure with the original::

       QualSFrame.py lvmSFrame-00014964.sky_subst.fits
       QualCFrame.py lvmCFrame-00014964.sky_subst.fits

   For 14964 (SkyE substituted for SkyW), the result was judged
   plausibly better than the original though not perfect.  Note that
   when the substituted telescope is close to the science field (SkyE
   was 2.2 deg from Vela), any emission from the target itself in that
   sky field is subtracted too.

``SubstituteSky.py``, ``RunSky.py``, ``QualSFrame.py`` and
``QualCFrame.py`` all import ``lvmdrp``, so they must be run in an
environment that has it, with ``$SAS_BASE_DIR`` and ``$LVM_MASTER_DIR``
reachable.


See Also
--------

- :doc:`api/eval_sky/index` - API documentation
- :doc:`api/SubstituteSky/index` - API documentation
- :doc:`api/RunSky/index` - API documentation
- :doc:`data_quality` - ``QualSFrame.py`` and ``QualCFrame.py``, used to
  compare a corrected exposure with the original
- :doc:`sky_methods` - running other sky-subtraction approaches on the
  corrected exposure
