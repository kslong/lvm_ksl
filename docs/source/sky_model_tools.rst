ESO and PALACE Sky Model Tools
==============================

This page is part of :doc:`sky_subtraction`.

These tools generate theoretical sky spectra using the ESO Sky Model and
the PALACE airglow model, which can be compared to observed sky spectra
for validation.  See :doc:`sky_models` for an overview of these models
and of the semi-empirical machine-learning approach built on PALACE.


EsoSkyObs.py
------------

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
------------

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


Typical Workflow
----------------

1. Generate a theoretical sky for the observation::

       EsoSkyObs.py 81.5 -66.0 60000.5

2. Compare with observed sky from SKY_EAST or SKY_WEST telescopes
3. Identify discrepancies that may indicate calibration issues


See Also
--------

- :doc:`sky_models` - how the ESO Sky Model, PALACE and the
  semi-empirical machine-learning approach relate to each other
- :doc:`sky_methods` - ``SkySepESO.py`` and ``SkySepPalace.py``, which use
  these models for sky subtraction
- :doc:`api/EsoSkyObs/index` - API documentation
- :doc:`api/PalaceObs/index` - API documentation
