Sky-Line Masks, Sky Spectra and Continuum Fits
==============================================

This page is part of :doc:`sky_subtraction`.

These tools support the development of improved sky subtraction by identifying
sky-line-free wavelength windows for continuum fitting and by assembling
stacked sky spectra from the LVM sky telescopes.  Together they are intended
to characterise the sky background well enough to constrain physical models
of the airglow emission.  What the PALACE-based decomposition used here
fits, and how its versions differ, is described in :doc:`sky_models`.


XSkySepIvan.py
--------------

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
-------------------

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
-------------------

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
-------------------------

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
-------------

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
------------------

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


Typical Workflow
----------------

1. Build a palace line mask for the field::

       palace_make_mask.py XCframe_file.fits

2. Collect sky spectra from repeated observations of the field::

       GetSky_from_CFrame_sum.py XCframe_file.fits Sky_WHAM_south_08

3. Fit the two-component B-spline continuum::

       GetSkyCont.py Sky_WHAM_south_08.fits -mask XCframe_file_mask.fits

4. Evaluate the fit interactively::

       GetSkyCont_eval.py skycont_Sky_WHAM_south_08.fits


See Also
--------

- :doc:`sky_models` - the PALACE-based decomposition used by
  ``XSkySepIvan.py``
- :doc:`sky_methods` - sky-subtraction methods that use the mask and
  continuum fits from this page
- :doc:`api/XSkySepIvan/index` - API documentation
- :doc:`api/XSkySepIvan_eval/index` - API documentation
- :doc:`api/palace_make_mask/index` - API documentation
- :doc:`api/GetSky_from_CFrame_sum/index` - API documentation
- :doc:`api/GetSkyCont/index` - API documentation
- :doc:`api/GetSkyCont_eval/index` - API documentation
