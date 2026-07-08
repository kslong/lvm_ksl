#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Evaluate sky subtraction quality for one or more FITS files produced
    by SkySubOrig.py, SkySubDrp.py, SkySubDev1.py, or SkySubDev2.py.
    All four scripts write files with WAVE, FLUX (sky-subtracted), and
    SKY extensions that this script reads directly.  When multiple files
    are supplied they are overlaid in the same figures, making it easy
    to compare sky subtraction methods on a common footing.

    Writes a single interactive HTML file containing six Plotly figures
    and two inline statistics tables.  Also writes per-arm continuum-
    residual statistics back into each evaluated file's own DRP_ALL
    table, in place, for detailed analysis (see Figure 5/6 below).

Command line usage (if any):

    usage: SkySub_eval.py [wmin wmax] [-num N] [-out outroot]
                          filename [filename ...]

    Arguments:

    filename         one or more SkySub output FITS files to evaluate.
    wmin             minimum wavelength in Angstroms (default 3600).
    wmax             maximum wavelength in Angstroms (default 9800).

    Options:

    -num N           overlay N randomly-selected individual spectra (default 20).
    -out outroot     combine all files into one HTML (default: per-file <stem>_eval.html).

Description:

    **Figure 1 — spectral overview (3 panels, shared x-axis, linear scale):**

    Panel 1 (Original Flux): FLUX + SKY, showing the pre-subtraction
    spectrum.  Median and 10th/90th percentile band, with N individual
    spectra overlaid in grey.  One colour per file when multiple files
    are supplied.

    Panel 2 (Sky-subtracted Flux): FLUX extension on the same linear
    scale.

    Panel 3 (Sky Model): SKY extension on the same scale.

    All three panels share a common y-range of -1e-14 to 1e-13
    erg s⁻¹ cm⁻² Å⁻¹, set to make sky-subtracted residuals visible.
    The x-axis is downsampled to ~2000 display pixels for performance.

    **Figure 2 — residual histograms at three diagnostic sky-line windows:**

    Diagnostic windows: [OI] 5577 (5560–5594 Å), [OI] 6300
    (6280–6320 Å), and IR OH (9300–9500 Å).  For each window the
    distribution of per-pixel HF (high-frequency) residual values across
    every spectrum is plotted as a bar histogram.  The HF residual is
    each spectrum minus a Gaussian-smoothed (clean-pixel-weighted)
    version of itself (see ``_hf_residual``), so leftover continuum in
    the window does not bias the histogram — only sky-line-scale
    structure remains.  A Gaussian with the same median and NMAD is
    overlaid as a dotted curve.  A statistics box (N, median, NMAD,
    skewness) is placed inside each panel; boxes stack vertically when
    multiple files are given.  The histogram range is a fixed
    ``_LINE_HIST_XLIM`` (-2e-15 to 2e-15), the same for every overlaid
    file, rather than each file's own median ± 5·NMAD -- this way
    differences in histogram width/shape between models are real, not
    an artifact of each file being binned over a different range.  The
    reported skewness is computed on that same fixed, clipped range so
    it matches what's plotted rather than being dominated by a handful
    of outlier pixels beyond the display range.

    **Figure 3 (Figure 4 in HTML) — diagnostic window median spectra:**

    For each of the three diagnostic windows a wider search region
    (5400–5750, 6100–6500, 9000–9800 Å) is plotted showing the median
    and 10th/90th percentile band of the original (dotted) and
    sky-subtracted (solid) spectra.  Two types of shading show which
    pixels are used by the HF RMS analysis:

    - Red shading — signal window (diagnostic window, used to measure
      sky-line RMS)
    - Green shading — background pixels (mask-selected sky-line-free
      pixels in the search region, used to estimate noise RMS)

    The sky-line mask is read from sky_mask.fits (palace_make_mask.py
    output) if available; otherwise all pixels outside the diagnostic
    window are treated as background.

    **Statistics table:**

    An HTML table between Figure 3 and Figure 4 reports per-file per-window:
    noise RMS (median background RMS), sky-line RMS before/after
    subtraction, and the 50th, 90th, and 95th percentiles of the
    HF RMS ratio (sky_sub / sky_orig) across all spectra, plus the
    fraction of spectra with ratio < 0.5.  The same information is
    printed to the terminal.

    **Figure 4 (Figure 3 in HTML) — HF RMS ratio per spectrum:**

    For each diagnostic window the noise-corrected high-frequency RMS
    ratio (sky-subtracted / original) is plotted against spectrum index,
    or against MJD if there are more than 100 spectra and a DRP_ALL
    OBSTIME column is available (e.g. the per-exposure rows written by
    SkySubSci.py) — MJD is derived from OBSTIME rather than the table's
    own MJD column, since that is stored as a truncated integer.  If any
    file in the run lacks this information, all files fall back to plain
    spectrum index so every trace stays on a common axis.
    The ratio is computed as follows:

    1. A Gaussian smooth (σ = 50 pixels ≈ 25 Å) using only mask-selected
       clean pixels is subtracted from each spectrum to produce a
       high-frequency (HF) residual array that retains sky-line-scale
       (≤ LSF FWHM) structure.
    2. The RMS of the HF residuals in the diagnostic window gives
       ``rms_diag`` (sky lines + noise).
    3. The RMS of the HF residuals in the background region gives
       ``rms_bg`` (noise only).
    4. The sky-line-only amplitude is ``sky_rms = √max(0, rms_diag² −
       rms_bg²)``.
    5. The ratio ``sky_sub / sky_orig`` is plotted.  A value of 0 means
       the sky lines were completely removed; a value of 1 means no
       change.  The panels share a linked x-axis (spectrum index) so
       zooming one panel pans all three.

    **Figure 5 — per-arm continuum-quality histograms (3 rows x 3 cols):**

    Three columns, one per spectrograph arm (B 3650-5775, R 5800-7520,
    Z 7570-9600 Å — the B/R and R/Z overlap zones and the outer edges
    are excluded).  Same annotation style throughout (N, median, NMAD,
    skew; Gaussian overlay), but each row plots a different quantity,
    and they answer different questions.  All nine panels share the
    same fixed x-axis range (``_CONT_HIST_XLIM``, -5e-15 to 5e-15)
    rather than each file's own median ± 5·NMAD, so overlaid files are
    binned identically and differences in width/shape between models
    are real rather than a range artifact (as with Figure 2's
    ``_LINE_HIST_XLIM``).  As with Figure 2, N/median/NMAD reflect the
    full underlying distribution; only the reported skewness is
    recomputed on the clipped/displayed range.

    Row 1 -- distribution (across spectra) of DRP_ALL['SCI_MED_<arm>'],
    i.e. the per-spectrum median of (raw science flux − sci_cont) in
    clean pixels of the *raw, pre-subtraction* science spectrum.  This
    **is** the science-side continuum-fit-quality test.

    Row 2 -- distribution of DRP_ALL['SKY_MED_<arm>'], the per-spectrum
    median of (raw sky flux − sky_cont) in clean pixels of the *raw,
    pre-subtraction* sky-telescope spectrum whose CONT was actually used
    in the final SKY.  This **is** the sky-side continuum-fit-quality
    test, and is the one that matters directly for the final result's
    continuum level, since sky_cont is literally part of SKY.

    Row 3 -- distribution (across spectra) of the per-spectrum median of
    the *raw* sky-subtracted FLUX (not the HF residual used by Figure 2)
    in clean (sky-line-free) pixels -- one value per spectrum, exactly
    like rows 1/2, so all three rows are directly comparable (not a pooled
    per-pixel distribution).  This is the final, post-subtraction leftover
    signal.  It does **not** test continuum-fit quality: the science
    spectrum's own continuum fit (sci_cont) is only ever used to derive
    the bisection line-scale target and never enters the subtracted
    result, so this row reflects real source (stellar) continuum
    entangled with any net error in the sky-side continuum estimate --
    see Notes.

    Rows 1/2 require SCI_MED_<arm>/SKY_MED_<arm> columns already present
    in the input file's DRP_ALL (written by SkySubOrig.py/SkySubDev1.py/
    SkySepESO.py); a file without them (SkySubDev2.py/SkySubDrp.py output,
    or an older file predating this feature) leaves those panels empty.

    **Continuum-quality stats table**, and **Figure 6 — per-arm,
    per-spectrum continuum residual (same row order: SCI, SKY, Subtracted):**

    A table (screen + HTML) reports per-file, per-arm N/median/NMAD/skew
    of the pooled clean-pixel residual (row 3 above).  Figure 6 plots the
    per-spectrum SCI_MED_<arm>/SKY_MED_<arm> (rows 1/2) and the per-spectrum
    median of the row-3 residual (row 3), each against spectrum index or
    MJD (same x-axis logic as Figure 4), so specific bad exposures/fibers
    are visible rather than just an aggregate number.  For row 3, a
    nonzero residual can reflect real source continuum as well as a sky
    error; comparing this panel across overlaid methods on the same
    input isolates sky-continuum quality, since source continuum is
    identical across methods.

    Every panel shares the same fixed y-axis range (``_CONT_TREND_YLIM``,
    -2e-14 to 2e-14) rather than an adaptive per-panel range, for the same
    cross-file-comparability reason as Figures 2/5's fixed x-axis ranges.
    Since a fixed range means real outliers are simply off-scale and
    invisible, each panel's annotation reports N and the percentage of
    that file's exposures falling below/above the range in each
    direction (``_add_trend_panel``), so an exposure that's off the
    visible plot is still accounted for.

    The same per-spectrum statistics (resid_med_<arm>, resid_nmad_<arm>,
    resid_rms_<arm>, resid_skew_<arm> for arm in b, r, z) are written
    back into each evaluated file's own DRP_ALL table, in place, for
    detailed analysis (matching GetSkyCont_eval.py's convention of
    updating the input file rather than only reporting in the HTML).

Notes:

    Requires plotly (pip install plotly) and scipy.
    If GetSkyCont.py is on the Python path, sky_mask.fits is loaded
    automatically from the current directory or the lvm_ksl data/
    directory; without the mask the Gaussian smooth is unweighted.

    Figures 5/6 and the DRP_ALL update require a DRP_ALL extension in
    the input file; files without one are skipped for this diagnostic
    only (the other figures still work).  Input files are modified in
    place when DRP_ALL is updated; no backup is created.

History::

    260630 ksl  Initial version (matplotlib, mode-based: eval/xplot/bigplot).
    260630 ksl  Rewritten to Plotly HTML with multi-file overlay.
    260630 ksl  Added HF RMS ratio analysis, diagnostic window spectra, stats table.
    260704 ksl  Figure 2 histograms now use the HF (continuum-subtracted)
                residual instead of raw FLUX, so leftover continuum in the
                diagnostic window no longer biases the reported median/NMAD.
    260704 ksl  Figure 2 skewness is now computed on the same median ± 5·NMAD
                clipped range as the histogram, instead of all pixels, so it
                is not dominated by a few outliers invisible in the plot.
    260704 ksl  Figure 3 (HF RMS ratio) plots against MJD (from DRP_ALL's
                OBSTIME) instead of spectrum index when n_spec > 100 and
                the information is available in every file, falling back
                to spectrum index otherwise.
    260706 ksl  Added Figures 5/6 and a stats table for per-arm (B/R/Z)
                continuum-quality: raw (not HF) sky-subtracted flux in
                clean pixels, since only sky lines were scaled by any
                SkySub* method -- the continuum is used exactly as
                fitted, so this is the diagnostic that can actually catch
                a bad continuum estimate.  Per-spectrum stats are written
                back into each evaluated file's DRP_ALL, in place.
    260707 ksl  Correction: Figure 5 row 1 (raw sky-subtracted flux in clean
                pixels) does NOT test continuum-fit quality -- the science
                spectrum's own continuum fit never enters the subtracted
                result, only the sky-side one does, and even that is
                entangled with real source continuum.  Added Figure 5 rows
                2/3: histograms of the DRP_ALL SCI_MED_<arm>/SKY_MED_<arm>
                columns across spectra (written by SkySubOrig.py/
                SkySubDev1.py/SkySepESO.py) -- these are the actual
                continuum-fit-quality metrics, evaluated against the raw
                pre-subtraction science/sky spectra respectively.  Figure 5
                is now 3 rows x 3 cols; row-1 histogram/annotation logic
                factored into a shared _add_hist_panel() helper.
    260707 ksl  Reordered Figure 5's rows to SCI, SKY, Subtracted (was
                Subtracted, SCI, SKY) -- continuum-fit-quality tests first,
                the net post-subtraction diagnostic last.  Figure 6 expanded
                from 1 row to the same 3-row x 3-col layout/order, plotting
                SCI_MED_<arm>/SKY_MED_<arm> per-spectrum (not just the
                Subtracted-row residual) against index/MJD.  Legend
                visibility is now tracked explicitly per file (first panel
                that actually has data), since SCI/SKY rows can be empty for
                files lacking those columns -- no longer tied to a fixed
                row/col.
    260707 ksl  Increased Figure 5/6 vertical_spacing (0.08/0.06 -> 0.14) and
                heights so rows no longer crowd each other.  Added top-level
                "Sky Line Subtraction" / "Continuum Separation" HTML section
                headers.  Replaced Figures 5/6's floating Plotly legend
                (which sat on top of a subplot in the 3x3 grid) with a
                colour-coded suptitle above each figure (_suptitle_text()).
                Extended the continuum-quality stats table/printout to cover
                SCI/SKY (not just the Subtracted row), with a Kind column.
    260707 ksl  Figure 5/6 row 3 (Subtracted) now histograms/plots the
                per-spectrum median of the raw sky-subtracted flux in clean
                pixels (one value per spectrum), not the pooled per-pixel
                distribution -- matches rows 1/2 (SCI_MED/SKY_MED are
                already per-spectrum medians), so all three rows are
                directly comparable.  Removed the now-unused pooled
                _window_stats() call in this loop.
    260707 ksl  Added an overall page title "Sky Subtraction Quality Check"
                (HTML <title> + top-of-page <h1>).  Moved the "Sky Line
                Subtraction" section header to sit directly above Figure 2
                (residual histograms) rather than above Figure 1 (spectral
                overview, which isn't part of either named section).
                Section headers demoted to <h2> under the new page <h1>.
    260708 ksl  Figure 2's per-window histograms now use a fixed x-axis
                range (_LINE_HIST_XLIM, -2e-15 to 2e-15) instead of each
                file's own median +/- 5*NMAD, so overlaid files are binned
                identically and differences in width/shape between models
                are directly comparable rather than range-dependent.
    260708 ksl  Same fix applied to Figure 5's continuum-quality histograms
                (SCI/SKY/Subtracted rows): _add_hist_panel() now accepts an
                optional fixed xlim, used here as _CONT_HIST_XLIM (-5e-15 to
                5e-15). N/median/NMAD still reflect the full distribution;
                only the displayed skewness is clipped, unchanged from
                before.
    260708 ksl  Figure 6 (continuum quantities vs index/MJD) now uses a
                fixed y-axis range too (_CONT_TREND_YLIM, -2e-14 to 2e-14,
                wider than _CONT_HIST_XLIM since per-spectrum values scatter
                more than the histogrammed bulk), replacing the previous
                adaptive 1st/99th-percentile-per-panel range. Since a fixed
                range hides real outliers off-scale, each panel now reports
                the percentage of exposures below/above the range via a new
                shared _add_trend_panel() helper (replacing the inlined
                fig6.add_trace() calls). Removed the now-unused
                all_arm_sci/all_arm_sky/all_arm_resid accumulators that only
                existed to feed the old adaptive range.

'''

import sys
import os
import warnings
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from GetSkyCont import load_mask, _interp_mask_to_wave
    _HAVE_MASK = True
except ImportError:
    _HAVE_MASK = False


# Diagnostic sky-line windows: (label, win_lo, win_hi, search_lo, search_hi)
# The search range is the broader region scanned for clean background pixels
# (mask=1) outside the diagnostic window itself.
_DIAG_WINDOWS = [
    ('OI 5577',  5560.0, 5594.0,  5400.0, 5750.0),
    ('OI 6300',  6280.0, 6320.0,  6100.0, 6500.0),
    ('IR 9300',  9300.0, 9500.0,  9000.0, 9800.0),
]

# Spectrograph arm ranges used for continuum-quality evaluation, trimmed to
# exclude the B/R and R/Z overlap zones (5775-5800 and 7520-7570) and the
# outer edges below 3650 / above 9600.  Distinct from palace_make_mask.py's
# ARM_RANGES, which is a full-coverage (no-gap) split for a different
# purpose (PALACE line-model mask construction).
_ARM_EVAL_RANGES = [
    ('B', 3650.0, 5775.0),
    ('R', 5800.0, 7520.0),
    ('Z', 7570.0, 9600.0),
]

# Fixed x-axis range for Figure 2's per-window HF (line) residual histograms.
# A common, fixed range (rather than each file's own median +/- 5*NMAD) means
# every overlaid file is binned identically, so differences in histogram
# width/shape between models are directly comparable rather than an artifact
# of each file having its own auto-scaled range.
_LINE_HIST_XLIM = (-2e-15, 2e-15)

# Same idea for Figure 5's continuum-separation histograms (SCI_MED/SKY_MED/
# Subtracted); wider than _LINE_HIST_XLIM since these are per-spectrum
# medians rather than per-pixel residuals.
_CONT_HIST_XLIM = (-5e-15, 5e-15)

# Fixed y-axis range for Figure 6 (the same quantities as Figure 5, plotted
# per-spectrum against index/MJD instead of histogrammed).  Wider than
# _CONT_HIST_XLIM: individual per-spectrum values scatter more than the bulk
# histogrammed distribution, so a wider window is needed to keep most
# exposures visible while still excluding genuine outliers.
_CONT_TREND_YLIM = (-2e-14, 2e-14)

# Per-file colour pairs: (solid line, shaded band)
_FILE_COLORS = [
    ('rgb(31,119,180)',   'rgba(31,119,180,0.20)'),
    ('rgb(214,39,40)',    'rgba(214,39,40,0.20)'),
    ('rgb(44,160,44)',    'rgba(44,160,44,0.20)'),
    ('rgb(148,103,189)', 'rgba(148,103,189,0.20)'),
    ('rgb(255,127,14)',  'rgba(255,127,14,0.20)'),
    ('rgb(23,190,207)',  'rgba(23,190,207,0.20)'),
]

_USAGE = '''Usage: SkySub_eval.py [wmin wmax] [-num N] [-out outroot] filename [filename ...]

Arguments:
  filename   one or more SkySub FITS files (SkySubOrig/Drp/Dev1/Dev2 output)
  wmin       spectral overview minimum wavelength in Angstroms (default 3600)
  wmax       spectral overview maximum wavelength in Angstroms (default 9800)

Options:
  -num N        number of individual spectra to overlay (default 20; 0 = band only)
  -out outroot  overlay all files into one HTML with this output root;
                without -out, each file produces its own <stem>_eval.html
'''


# ──────────────────────────────────────────────────────────────
# Analysis helpers
# ──────────────────────────────────────────────────────────────

def _read_mjd(hdul, n_spec):
    '''Precise per-row MJD from the DRP_ALL table's OBSTIME column.

    OBSTIME is an ISO-format timestamp string (e.g. from SkySubSci.py's
    per-exposure DRP_ALL rows); it is used instead of the table's MJD
    column because that MJD is stored as a truncated integer.

    Returns a 1-D float array of length n_spec, or None if DRP_ALL is
    missing, has no OBSTIME-like column, row count doesn't match
    n_spec, or the timestamps can't be parsed.
    '''
    if 'DRP_ALL' not in hdul:
        return None
    data = hdul['DRP_ALL'].data
    names = {n.lower(): n for n in data.columns.names}
    if 'obstime' not in names or len(data) != n_spec:
        return None
    try:
        obstime = np.asarray(data[names['obstime']], dtype=str)
        return Time(obstime, format='isot', scale='utc').mjd.astype(float)
    except Exception:
        return None


def _window_stats(flux, wave, wmin, wmax, clean=None):
    '''Pixel-level residual statistics within a wavelength window.

    Parameters
    ----------
    flux : 2-D array (n_spec, n_wave), sky-subtracted flux
    wave : 1-D array (n_wave,)
    wmin, wmax : float
    clean : 1-D boolean array (n_wave,), optional
        If given, True = sky-line-free pixel; only these pixels within
        the window are used (e.g. for continuum-quality statistics,
        where sky-line pixels must be excluded).

    Returns
    -------
    dict with keys n, med, nmad, skew, p10, p90, vals — or None if no data.
    '''
    sel  = (wave >= wmin) & (wave <= wmax)
    if clean is not None:
        sel = sel & clean
    if not sel.any():
        return None
    vals = flux[:, sel].astype(float).flatten()
    vals = vals[np.isfinite(vals)]
    if len(vals) < 10:
        return None
    med  = float(np.median(vals))
    nmad = float(1.4826 * np.median(np.abs(vals - med)))
    mn   = float(np.mean(vals))
    sd   = float(np.std(vals, ddof=0))
    skew = float(np.mean(((vals - mn) / sd) ** 3)) if sd > 0 else 0.0
    p10, p90 = (float(v) for v in np.percentile(vals, [10, 90]))
    return dict(n=len(vals), med=med, nmad=nmad, skew=skew,
                p10=p10, p90=p90, vals=vals)


def _hf_residual(flux_2d, clean_mask=None, sigma_pix=50):
    '''High-frequency residual: each spectrum minus a Gaussian-smoothed version.

    sigma_pix=50 ≈ 25 Å at LVM sampling, well above the 1.3 Å LSF FWHM,
    so the residual retains only sky-line-scale and finer structure.

    clean_mask : 1-D boolean array (n_wave,), True = sky-line-free pixel.
        Sky-line pixels are zeroed out before smoothing so their flux does
        not leak into the continuum estimate (same masking as GetSkyCont).
        If None, all pixels are treated as clean.

    Uses a weighted Gaussian smooth so masked pixels are reconstructed from
    their neighbours rather than pulled toward zero.
    '''
    filled = np.where(np.isfinite(flux_2d), flux_2d, 0.0)
    if clean_mask is None:
        smooth = gaussian_filter1d(filled, sigma=sigma_pix, axis=1)
    else:
        w = clean_mask.astype(float)               # (n_wave,)
        w2d = np.broadcast_to(w, filled.shape).copy()
        num = gaussian_filter1d(filled * w2d, sigma=sigma_pix, axis=1)
        den = gaussian_filter1d(w2d,           sigma=sigma_pix, axis=1)
        smooth = num / np.where(den > 1e-6, den, 1.0)
    hf = flux_2d - smooth
    return hf


def _per_spec_hf_rms(hf_2d, wave, wmin, wmax):
    '''Per-spectrum RMS of high-frequency residuals within a wavelength window.

    RMS = sqrt(mean(hf^2)) responds to squared amplitude, so even a few
    bright sky-line pixels dominate over the many clean pixels — unlike NMAD
    which is suppressed by the majority of non-line pixels.  The quadrature
    formula rms_sky = sqrt(rms_diag^2 - rms_clean^2) then isolates the
    sky-line contribution.  Equally sensitive to positive and negative
    deviations.

    Returns 1-D array of length n_spec, or None if no pixels fall in window.
    '''
    sel = (wave >= wmin) & (wave <= wmax)
    if not sel.any():
        return None
    r = hf_2d[:, sel].astype(float)
    r[~np.isfinite(r)] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        rms = np.sqrt(np.nanmean(r ** 2, axis=1))
    return rms


def _per_spec_hf_rms_sel(hf_2d, sel):
    '''Per-spectrum RMS of HF residuals for an arbitrary pixel selection.

    sel : 1-D boolean array (n_wave,)
    Returns 1-D array of length n_spec, or None if sel is empty.
    '''
    if not sel.any():
        return None
    r = hf_2d[:, sel].astype(float)
    r[~np.isfinite(r)] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        rms = np.sqrt(np.nanmean(r ** 2, axis=1))
    return rms


def _per_spec_arm_stats(wave, clean, flux, wlo, whi):
    '''Per-spectrum continuum-residual statistics within [wlo, whi].

    Unlike the HF-residual helpers above, this operates on the raw
    (not Gaussian-smoothed) sky-subtracted flux, restricted to
    sky-line-free (clean) pixels — i.e. it measures continuum-level
    residual, not line-level residual.  Mirrors GetSkyCont_eval.py's
    _per_spec_stats.

    Parameters
    ----------
    wave  : 1-D array (n_wave,), full-resolution wavelength
    clean : 1-D bool array (n_wave,) or None, True = sky-line-free pixel
    flux  : 2-D array (n_spec, n_wave), sky-subtracted flux
    wlo, whi : float, wavelength bounds

    Returns
    -------
    dict with 1-D float arrays of length n_spec: med, nmad, rms, skew
    or None if no clean pixels fall within the range.
    '''
    sel = (wave >= wlo) & (wave <= whi)
    if clean is not None:
        sel = sel & clean
    if not sel.any():
        return None
    r = flux[:, sel].astype(float)
    r[~np.isfinite(r)] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        med  = np.nanmedian(r, axis=1)
        nmad = 1.4826 * np.nanmedian(np.abs(r - med[:, np.newaxis]), axis=1)
        rms  = np.sqrt(np.nanmean(r ** 2, axis=1))
        mn   = np.nanmean(r, axis=1)
        sd   = np.nanstd(r, axis=1)
        skew = np.where(sd > 0, np.nanmean(((r - mn[:, np.newaxis])
                                            / sd[:, np.newaxis]) ** 3, axis=1), 0.0)
    return dict(med=med, nmad=nmad, rms=rms, skew=skew)


# ──────────────────────────────────────────────────────────────
# Plotly trace helpers
# ──────────────────────────────────────────────────────────────

def _band_traces(wave, arr, name, color_med, color_band,
                 legend_ref='legend', show_legend=True):
    '''Median line + shaded 10/90-pct band traces.'''
    p10 = np.nanpercentile(arr, 10, axis=0)
    p90 = np.nanpercentile(arr, 90, axis=0)
    med = np.nanmedian(arr, axis=0)
    return [
        go.Scatter(x=wave, y=p10, mode='lines', line=dict(width=0),
                   showlegend=False, hoverinfo='skip'),
        go.Scatter(x=wave, y=p90, mode='lines', line=dict(width=0),
                   fill='tonexty', fillcolor=color_band,
                   name=f'{name} 10–90 pct', showlegend=False,
                   hoverinfo='skip'),
        go.Scatter(x=wave, y=med, mode='lines',
                   line=dict(color=color_med, width=1.5),
                   name=name, showlegend=show_legend, legend=legend_ref),
    ]


def _sample_traces(wave, arr, idx):
    '''Semi-transparent individual spectrum traces.'''
    if len(idx) == 0:
        return []
    alpha = min(0.7, 3.0 / max(1.0, len(idx) ** 0.5))
    return [
        go.Scatter(x=wave, y=arr[i], mode='lines',
                   line=dict(color=f'rgba(120,120,120,{alpha:.2f})', width=0.6),
                   showlegend=False, hoverinfo='skip')
        for i in idx
    ]


def _fmt(v):
    '''Format a flux value compactly in scientific notation.'''
    if not np.isfinite(v) or abs(v) < 1e-30:
        return '0'
    return f'{v:.2e}'


def _suptitle_text(file_legend):
    '''
    Build an HTML colour-swatch legend string ("<sq> label   <sq> label ...")
    from a list of (label, colour) pairs, for use as a figure-level suptitle.

    Used instead of Plotly's floating legend box for figures with a dense
    subplot grid, where the floating legend (positioned in the figure's
    top-right corner regardless of subplot boundaries) ends up sitting on
    top of -- and obscuring -- whichever subplot happens to be there.
    '''
    parts = [f'<span style="color:{c}">&#9632;</span> {lbl}' for lbl, c in file_legend]
    return '&nbsp;&nbsp;&nbsp;'.join(parts)


def _add_hist_panel(fig, row, col, n_cols, vals, label, c_med, c_hist, k,
                    show_legend, stat_label='Med', xlim=None):
    '''
    Histogram + Gaussian overlay + stats annotation in one subplot panel,
    from a 1-D array of values already collected for one file.  Shared by
    Figure 5's three rows (SCI_MED, SKY_MED, pooled pixel residual) so the
    stats/annotation logic isn't tripled.

    row, col : 1-based subplot position.
    n_cols : number of columns in the subplot grid (for annotation axis refs).
    vals : 1-D array; non-finite values are dropped internally.
    k : file index, used to stack annotation boxes vertically for
        multiple overlaid files.
    show_legend : whether this panel's bar trace shows the file's legend
        entry.  Caller's responsibility -- since some rows can be empty
        for a given file (e.g. SCI_MED/SKY_MED absent for SkySubDev2.py/
        SkySubDrp.py output), the legend must be shown on the first panel
        that actually has data, not hardcoded to a fixed row/col.
    stat_label : label for the central-tendency stat in the annotation
        (e.g. 'Med').
    xlim : optional (lo, hi) tuple giving a fixed histogram/display range,
        the same for every overlaid file, so differences in width/shape
        between files are real rather than an artifact of each file being
        binned over its own adaptive range.  If omitted, falls back to
        this file's own median +/- 5*NMAD (or p10/p90 if that's degenerate).

    Returns a stats dict (n, med, nmad, skew), or None if too few
    finite values to plot (caller should not count this as "shown").
    N/med/nmad always reflect the full (unclipped) vals array; only the
    displayed skewness is recomputed on the clipped/displayed range, since
    skewness (unlike median/NMAD) is dominated by a handful of outliers
    that may sit outside the display window.
    '''
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) < 3:
        return None

    med  = float(np.median(vals))
    nmad = float(1.4826 * np.median(np.abs(vals - med)))
    n    = len(vals)

    if xlim is not None:
        lo, hi = xlim
    else:
        lo = med - 5.0 * nmad
        hi = med + 5.0 * nmad
        if lo >= hi:
            p10, p90 = np.percentile(vals, [10, 90])
            lo, hi = p10 - 1e-20, p90 + 1e-20

    vals_disp = vals[(vals >= lo) & (vals <= hi)]
    if len(vals_disp) > 2:
        mn_d, sd_d = np.mean(vals_disp), np.std(vals_disp, ddof=0)
        skew_disp = (float(np.mean(((vals_disp - mn_d) / sd_d) ** 3))
                    if sd_d > 0 else 0.0)
    else:
        mn, sd = np.mean(vals), np.std(vals)
        skew_disp = float(np.mean(((vals - mn) / sd) ** 3)) if sd > 0 else 0.0

    counts, edges = np.histogram(vals, bins=60, range=(lo, hi))
    centers = 0.5 * (edges[:-1] + edges[1:])
    bw      = edges[1] - edges[0]

    fig.add_trace(
        go.Bar(x=centers, y=counts, name=label, legendgroup=label,
              showlegend=show_legend,
              marker_color=c_hist, marker_line_width=0),
        row=row, col=col,
    )
    xg  = np.linspace(lo, hi, 300)
    sig = max(nmad, 1e-30)
    yg  = (n * bw / (np.sqrt(2 * np.pi) * sig)
          * np.exp(-0.5 * ((xg - med) / sig) ** 2))
    fig.add_trace(
        go.Scatter(x=xg, y=yg, mode='lines',
                  line=dict(color=c_med, width=1.5, dash='dot'),
                  showlegend=False, hoverinfo='skip'),
        row=row, col=col,
    )

    subplot_idx = (row - 1) * n_cols + col
    ax_sfx = '' if subplot_idx == 1 else str(subplot_idx)
    y_top  = 0.97 - k * 0.23
    if y_top > 0.08:
        fig.add_annotation(
            x=0.97, y=y_top,
            xref=f'x{ax_sfx} domain', yref=f'y{ax_sfx} domain',
            xanchor='right', yanchor='top',
            showarrow=False, align='right',
            font=dict(size=9, family='monospace'),
            bgcolor='rgba(255,255,255,0.80)',
            bordercolor=c_med, borderwidth=1,
            text=(f'N    = {n}<br>'
                  f'{stat_label}  = {_fmt(med)}<br>'
                  f'NMAD = {_fmt(nmad)}<br>'
                  f'Skew = {skew_disp:.2f}'),
        )
    return dict(n=n, med=med, nmad=nmad, skew=skew_disp)


def _add_trend_panel(fig, row, col, n_cols, x, y, label, c_med, k,
                     show_legend, ylim):
    '''
    Scatter trace + out-of-range annotation in one Figure 6 subplot panel.

    Every file/panel shares the same fixed ylim, so the fraction of
    exposures actually visible on the plot (vs. excluded by the axis
    range) is meaningful when comparing overlaid files -- a point outside
    ylim is otherwise simply invisible with no indication it exists.
    Reports N and the percentage of finite values below/above ylim in an
    annotation box, stacked vertically per file (index k).

    row, col : 1-based subplot position.
    n_cols : number of columns in the subplot grid (for annotation axis refs).
    ylim : fixed (lo, hi) y-axis range shared by every panel/file.
    '''
    y = np.asarray(y, dtype=float)
    fig.add_trace(
        go.Scatter(x=x, y=y, mode='markers',
                  marker=dict(color=c_med, size=4),
                  name=label, legendgroup=label,
                  showlegend=show_legend),
        row=row, col=col,
    )

    good = y[np.isfinite(y)]
    lo, hi = ylim
    n = len(good)
    if n > 0:
        pct_lo = 100.0 * np.count_nonzero(good < lo) / n
        pct_hi = 100.0 * np.count_nonzero(good > hi) / n
    else:
        pct_lo = pct_hi = 0.0

    subplot_idx = (row - 1) * n_cols + col
    ax_sfx = '' if subplot_idx == 1 else str(subplot_idx)
    y_top  = 0.97 - k * 0.10
    if y_top > 0.03:
        fig.add_annotation(
            x=0.97, y=y_top,
            xref=f'x{ax_sfx} domain', yref=f'y{ax_sfx} domain',
            xanchor='right', yanchor='top',
            showarrow=False, align='right',
            font=dict(size=8, family='monospace', color=c_med),
            bgcolor='rgba(255,255,255,0.80)',
            bordercolor=c_med, borderwidth=1,
            text=f'N={n}  <lo:{pct_lo:.0f}%  >hi:{pct_hi:.0f}%',
        )


# ──────────────────────────────────────────────────────────────
# Main evaluation function
# ──────────────────────────────────────────────────────────────

def plot_eval(filenames, wmin=3600.0, wmax=9800.0, n_sample=20, outroot=''):
    '''Build a three-figure Plotly HTML evaluation file.

    Parameters
    ----------
    filenames : str or list of str
        SkySub output FITS file(s).
    wmin, wmax : float
        Wavelength range for the spectral overview panels.
    n_sample : int
        Number of randomly-selected spectra to overlay on band panels.
    outroot : str
        Output filename root.  Default: <first stem>_eval.
    '''
    if isinstance(filenames, str):
        filenames = [filenames]

    # Validate up front: warn and drop missing or unreadable files
    valid = []
    for fn in filenames:
        if not Path(fn).exists():
            print(f'Warning: file not found, skipping: {fn}')
        else:
            valid.append(fn)
    filenames = valid
    if not filenames:
        print('Error: no valid input files.')
        return

    # Load sky-line mask (same as GetSkyCont) to avoid leaking lines into smooth.
    # Raw mask is interpolated to each file's wave grid inside the per-file loop.
    _raw_mask = None
    if _HAVE_MASK:
        _data_dir = Path(__file__).parent.parent / 'data'
        for _candidate in [Path('sky_mask.fits'), _data_dir / 'sky_mask.fits']:
            if _candidate.exists():
                try:
                    _raw_mask = load_mask(str(_candidate))  # (mwave, marr)
                    print(f'Loaded sky mask: {_candidate}')
                except Exception as _e:
                    print(f'Warning: could not load mask {_candidate}: {_e}')
                break
    if _raw_mask is None:
        print('Warning: sky_mask.fits not found; smoothing without sky-line mask')

    rng = np.random.default_rng(42)

    _ax = dict(showline=True, linewidth=1, linecolor='black',
               mirror=True, ticks='outside', ticklen=4, showticklabels=True)
    _leg_style = dict(xanchor='right', yanchor='top', x=0.99,
                      bgcolor='rgba(255,255,255,0.7)')

    # ══════════════════════════════════════════════════════════
    # Figure 1 — spectral overview (orig flux + sky-sub + sky, 3 rows)
    # ══════════════════════════════════════════════════════════
    fig1 = make_subplots(
        rows=3, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.04,
        subplot_titles=[
            f'Original Flux  (sky-sub + sky model)  ({len(filenames)} file(s))',
            'Sky-subtracted Flux',
            'Sky Model',
        ],
    )

    # ══════════════════════════════════════════════════════════
    # Figure 2 — per-window residual histograms (1 row × 3 cols)
    # ══════════════════════════════════════════════════════════
    fig2 = make_subplots(
        rows=1, cols=3,
        horizontal_spacing=0.08,
        subplot_titles=[
            f'{lbl}  ({int(lo)}–{int(hi)} Å)'
            for lbl, lo, hi, *_ in _DIAG_WINDOWS
        ],
    )

    # ══════════════════════════════════════════════════════════
    # Figure 3 — per-spectrum HF RMS ratio at each diagnostic window
    # ══════════════════════════════════════════════════════════
    fig3 = make_subplots(
        rows=1, cols=3,
        horizontal_spacing=0.08,
        subplot_titles=[f'{lbl}  HF RMS ratio (sub/orig)' for lbl, *_ in _DIAG_WINDOWS],
    )

    # ══════════════════════════════════════════════════════════
    # Figure 4 — median spectra per diagnostic window with shading
    # ══════════════════════════════════════════════════════════
    fig4 = make_subplots(
        rows=1, cols=3,
        horizontal_spacing=0.08,
        subplot_titles=[f'{lbl}  ({int(slo)}–{int(shi)} Å)'
                        for lbl, *_, slo, shi in _DIAG_WINDOWS],
    )

    # ══════════════════════════════════════════════════════════
    # Figure 5 — per-arm continuum-quality histograms (3 rows × 3 cols)
    # Row 1: distribution of SCI_MED_<arm> across spectra (science-side
    #        continuum-fit quality, against the raw pre-subtraction spectrum)
    # Row 2: distribution of SKY_MED_<arm> across spectra (sky-side
    #        continuum-fit quality, against the raw pre-subtraction spectrum)
    # Row 3: pooled raw sky-subtracted flux in clean pixels (post-subtraction
    #        leftover signal -- source continuum + any net sky-continuum error)
    # Rows 1/2 read columns written by SkySubOrig.py/SkySubDev1.py/SkySepESO.py;
    # a file without them just leaves those panels empty (see main loop).
    # ══════════════════════════════════════════════════════════
    fig5 = make_subplots(
        rows=3, cols=3,
        horizontal_spacing=0.08, vertical_spacing=0.14,
        subplot_titles=(
            [f'SCI {arm}  continuum-fit quality' for arm, *_ in _ARM_EVAL_RANGES]
            + [f'SKY {arm}  continuum-fit quality' for arm, *_ in _ARM_EVAL_RANGES]
            + [f'{arm}  ({int(lo)}–{int(hi)} Å)' for arm, lo, hi in _ARM_EVAL_RANGES]
        ),
    )

    # ══════════════════════════════════════════════════════════
    # Figure 6 — per-arm, per-spectrum continuum residual vs index/MJD
    # (3 rows × 3 cols, same row order as Figure 5: SCI, SKY, Subtracted)
    # ══════════════════════════════════════════════════════════
    fig6 = make_subplots(
        rows=3, cols=3,
        horizontal_spacing=0.08, vertical_spacing=0.14,
        subplot_titles=(
            [f'SCI {arm}  continuum residual' for arm, *_ in _ARM_EVAL_RANGES]
            + [f'SKY {arm}  continuum residual' for arm, *_ in _ARM_EVAL_RANGES]
            + [f'{arm}  continuum residual (median)' for arm, *_ in _ARM_EVAL_RANGES]
        ),
    )

    all_n_spec   = []
    all_used_mjd = []      # per-file: True if Figure 3 x-axis used MJD
    all_mjd_vals = []      # collected MJD values, for tight axis range
    all_hf_ratio = {lbl: [] for lbl, *_ in _DIAG_WINDOWS}  # for y-range
    _fig4_shaded = False   # add background/signal shading once (first file)
    _stats_rows  = []      # collected stats for screen + HTML output
    _arm_stats_rows = []   # per-file, per-arm continuum stats
    _file_legend = []      # (label, colour) pairs, for Figures 5/6 suptitles

    for k, filename in enumerate(filenames):
        c_med, c_band = _FILE_COLORS[k % len(_FILE_COLORS)]
        # Semi-transparent solid colour for histogram bars
        c_hist = c_med.replace('rgb(', 'rgba(').replace(')', ',0.55)')
        label  = Path(filename).stem

        try:
            hdul = fits.open(filename)
            wave = hdul['WAVE'].data.astype(float)
            flux = hdul['FLUX'].data.astype(float)
            sky  = hdul['SKY'].data.astype(float)
            mjd_col = _read_mjd(hdul, flux.shape[0] if flux.ndim > 1 else 1)
            drp_all = (Table(hdul['DRP_ALL'].data)
                      if 'DRP_ALL' in hdul else None)
            hdul.close()
        except Exception as e:
            print(f'Warning: could not read {filename} ({e}), skipping.')
            continue

        _file_legend.append((label, c_med))

        if flux.ndim == 1:
            flux = flux[np.newaxis, :]
            sky  = sky[np.newaxis, :]

        n_spec = len(flux)
        all_n_spec.append(n_spec)
        use_mjd = n_spec > 100 and mjd_col is not None
        all_used_mjd.append(use_mjd)
        if use_mjd:
            all_mjd_vals.append(mjd_col)

        # Wavelength window + downsample to ~2000 display pixels for Figure 1
        wmask  = (wave >= wmin) & (wave <= wmax)
        ds     = max(1, int(wmask.sum()) // 2000)
        pidx   = np.where(wmask)[0][::ds]
        w      = wave[pidx]
        f_disp = flux[:, pidx]
        s_disp = sky[:,  pidx]

        sample_idx = (rng.choice(n_spec, size=min(n_sample, n_spec), replace=False)
                      if n_sample > 0 else [])

        # ── HF residuals: build clean_mask and Gaussian-smooth-subtract ─
        # (needed by Figure 2 histograms as well as Figures 3 & 4)
        clean_mask_full = (_interp_mask_to_wave(_raw_mask[0], _raw_mask[1], wave)
                           if _raw_mask is not None and _HAVE_MASK else None)
        # Also on the downsampled grid used for Figure 1
        clean_mask = (_interp_mask_to_wave(_raw_mask[0], _raw_mask[1], w)
                      if _raw_mask is not None and _HAVE_MASK else None)

        hf_orig = _hf_residual(f_disp + s_disp, clean_mask=clean_mask)
        hf_sub  = _hf_residual(f_disp,           clean_mask=clean_mask)
        # Figure 3 x-axis: MJD (from OBSTIME) when there are enough spectra
        # to make a time trend meaningful and the information is available;
        # otherwise fall back to plain spectrum index.
        spec_x  = mjd_col if use_mjd else np.arange(n_spec)

        # ── Figures 5 & 6: per-arm continuum-quality diagnostics ──
        # Row order (both figures): SCI, SKY, Subtracted -- science-side and
        # sky-side continuum-fit quality first (the actual "did we separate
        # continuum from lines well" tests), post-subtraction leftover signal
        # last (a different, net diagnostic -- see module docstring).
        # Legend visibility is tracked explicitly per file rather than tied to
        # a fixed row/col, since SCI/SKY panels can be empty for files that
        # lack those DRP_ALL columns (e.g. SkySubDev2.py/SkySubDrp.py output).
        _arm_drp_cols = {}
        _shown5 = False
        _shown6 = False
        for col_idx, (arm, alo, ahi) in enumerate(_ARM_EVAL_RANGES, start=1):
            per_spec = _per_spec_arm_stats(wave, clean_mask_full, flux, alo, ahi)

            if per_spec is not None:
                _arm_drp_cols['resid_med_'  + arm.lower()] = per_spec['med']
                _arm_drp_cols['resid_nmad_' + arm.lower()] = per_spec['nmad']
                _arm_drp_cols['resid_rms_'  + arm.lower()] = per_spec['rms']
                _arm_drp_cols['resid_skew_' + arm.lower()] = per_spec['skew']

            # Figures 5 & 6, rows 1/2: distribution / trend of the per-spectrum
            # SCI_MED_<arm>/SKY_MED_<arm> DRP_ALL columns -- the actual
            # continuum-fit-quality metrics (against the raw pre-subtraction
            # science/sky spectra), written by SkySubOrig.py/SkySubDev1.py/
            # SkySepESO.py.  Absent for other methods, or files predating this
            # feature -- panels just stay empty in that case.
            for row_idx, side in [(1, 'SCI'), (2, 'SKY')]:
                col_name = f'{side}_MED_{arm}'
                if drp_all is not None and col_name in drp_all.colnames:
                    vals = np.asarray(drp_all[col_name], dtype=float)
                    hstats = _add_hist_panel(fig5, row_idx, col_idx, 3, vals,
                                             label, c_med, c_hist, k,
                                             show_legend=not _shown5,
                                             xlim=_CONT_HIST_XLIM)
                    if hstats is not None:
                        _shown5 = True
                        _arm_stats_rows.append(dict(
                            file=label, arm=arm, kind=side, n=hstats['n'],
                            med=hstats['med'], nmad=hstats['nmad'],
                            skew=hstats['skew'],
                        ))
                    _add_trend_panel(fig6, row_idx, col_idx, 3, spec_x, vals,
                                     label, c_med, k, show_legend=not _shown6,
                                     ylim=_CONT_TREND_YLIM)
                    _shown6 = True

            # Figures 5 & 6, row 3: per-spectrum median of the raw
            # sky-subtracted flux in clean pixels (post-subtraction leftover
            # signal).  Histogrammed as per-spectrum medians -- one value per
            # spectrum -- like rows 1/2 (SCI_MED/SKY_MED), not pooled
            # per-pixel values, so all three rows are directly comparable.
            if per_spec is not None:
                row3_stats = _add_hist_panel(fig5, 3, col_idx, 3, per_spec['med'],
                                             label, c_med, c_hist, k,
                                             show_legend=not _shown5,
                                             xlim=_CONT_HIST_XLIM)
                if row3_stats is not None:
                    _shown5 = True
                    _arm_stats_rows.append(dict(
                        file=label, arm=arm, kind='Subtracted', n=row3_stats['n'],
                        med=row3_stats['med'], nmad=row3_stats['nmad'],
                        skew=row3_stats['skew'],
                    ))

                _add_trend_panel(fig6, 3, col_idx, 3, spec_x, per_spec['med'],
                                 label, c_med, k, show_legend=not _shown6,
                                 ylim=_CONT_TREND_YLIM)
                _shown6 = True

        # ── Write continuum-quality stats back into this file's DRP_ALL ──
        if drp_all is not None and _arm_drp_cols:
            for col, vals in _arm_drp_cols.items():
                drp_all[col] = np.asarray(vals, dtype=np.float32)
            try:
                with fits.open(filename, mode='update') as hdul_upd:
                    new_drp = fits.table_to_hdu(drp_all)
                    new_drp.name = 'DRP_ALL'
                    for i, h in enumerate(hdul_upd):
                        if h.name == 'DRP_ALL':
                            hdul_upd[i] = new_drp
                            break
                    hdul_upd.flush()
                print(f'Updated DRP_ALL in {filename} with continuum-residual stats')
            except Exception as e:
                print(f'Warning: could not update DRP_ALL in {filename} ({e})')

        # ── Figure 1, Panel 1: original flux (FLUX+SKY, linear) ─
        orig_disp = f_disp + s_disp
        for tr in _sample_traces(w, orig_disp, sample_idx):
            fig1.add_trace(tr, row=1, col=1)
        for tr in _band_traces(w, orig_disp, label, c_med, c_band,
                               legend_ref='legend', show_legend=True):
            fig1.add_trace(tr, row=1, col=1)

        # ── Figure 1, Panel 2: sky-subtracted FLUX (linear) ───
        for tr in _sample_traces(w, f_disp, sample_idx):
            fig1.add_trace(tr, row=2, col=1)
        for tr in _band_traces(w, f_disp, label, c_med, c_band,
                               legend_ref='legend2', show_legend=True):
            fig1.add_trace(tr, row=2, col=1)

        # ── Figure 1, Panel 3: SKY (linear) ───────────────────
        for tr in _band_traces(w, s_disp, label, c_med, c_band,
                               legend_ref='legend3', show_legend=True):
            fig1.add_trace(tr, row=3, col=1)

        # ── Figure 2: diagnostic histograms ───────────────────
        # Histogram the HF (continuum-subtracted) residual, not raw FLUX,
        # so leftover continuum in the diagnostic window doesn't bias the
        # median/NMAD away from genuine sky-line residuals.
        for col_idx, (win_lbl, win_lo, win_hi, *_) in enumerate(_DIAG_WINDOWS, start=1):
            st = _window_stats(hf_sub, w, win_lo, win_hi)
            if st is None:
                continue

            # Fixed range (not each file's own median +/- 5*NMAD) so that
            # overlaid files are binned identically and directly comparable.
            lo, hi = _LINE_HIST_XLIM

            # Skew of the moment-based estimator is dominated by rare
            # outliers beyond the display clip, so recompute it on the
            # same clipped pixels shown in the histogram — this is what
            # the eye actually judges "symmetric" or "skewed" against.
            vals_disp = st['vals']
            vals_disp = vals_disp[(vals_disp >= lo) & (vals_disp <= hi)]
            if len(vals_disp) > 2:
                mn_disp = np.mean(vals_disp)
                sd_disp = np.std(vals_disp, ddof=0)
                skew_disp = (float(np.mean(((vals_disp - mn_disp) / sd_disp) ** 3))
                             if sd_disp > 0 else 0.0)
            else:
                skew_disp = st['skew']

            counts, edges = np.histogram(st['vals'], bins=60, range=(lo, hi))
            centers = 0.5 * (edges[:-1] + edges[1:])
            bw      = edges[1] - edges[0]

            fig2.add_trace(
                go.Bar(x=centers, y=counts,
                       name=label, legendgroup=label,
                       showlegend=(col_idx == 1),
                       marker_color=c_hist, marker_line_width=0),
                row=1, col=col_idx,
            )

            # Gaussian overlay (dotted, same solid colour)
            xg  = np.linspace(lo, hi, 300)
            sig = max(st['nmad'], 1e-30)
            yg  = (st['n'] * bw / (np.sqrt(2 * np.pi) * sig)
                   * np.exp(-0.5 * ((xg - st['med']) / sig) ** 2))
            fig2.add_trace(
                go.Scatter(x=xg, y=yg, mode='lines',
                           line=dict(color=c_med, width=1.5, dash='dot'),
                           showlegend=False, hoverinfo='skip'),
                row=1, col=col_idx,
            )

            # Stats annotation; stack vertically for multiple files
            ax_sfx = '' if col_idx == 1 else str(col_idx)
            y_top  = 0.97 - k * 0.23
            if y_top > 0.08:
                fig2.add_annotation(
                    x=0.97, y=y_top,
                    xref=f'x{ax_sfx} domain', yref=f'y{ax_sfx} domain',
                    xanchor='right', yanchor='top',
                    showarrow=False, align='right',
                    font=dict(size=9, family='monospace'),
                    bgcolor='rgba(255,255,255,0.80)',
                    bordercolor=c_med, borderwidth=1,
                    text=(f'N    = {st["n"]}<br>'
                          f'Med  = {_fmt(st["med"])}<br>'
                          f'NMAD = {_fmt(st["nmad"])}<br>'
                          f'Skew = {skew_disp:.2f}'),
                )

        # ── Figures 3 & 4: HF RMS ratio and median spectra ────
        for col_idx, (win_lbl, win_lo, win_hi, slo, shi) in \
                enumerate(_DIAG_WINDOWS, start=1):

            # ── pixel selections on the downsampled grid ──────
            diag_sel  = (w >= win_lo) & (w <= win_hi)
            # background: search range outside diagnostic window, mask=clean
            bg_sel = (w >= slo) & (w <= shi) & ~diag_sel
            if clean_mask is not None:
                bg_sel &= clean_mask

            rms_orig_diag = _per_spec_hf_rms_sel(hf_orig, diag_sel)
            rms_sub_diag  = _per_spec_hf_rms_sel(hf_sub,  diag_sel)
            rms_orig_bg   = _per_spec_hf_rms_sel(hf_orig, bg_sel)
            rms_sub_bg    = _per_spec_hf_rms_sel(hf_sub,  bg_sel)
            if rms_orig_diag is None or rms_sub_diag is None:
                continue
            noise_orig = rms_orig_bg if rms_orig_bg is not None else np.zeros(n_spec)
            noise_sub  = rms_sub_bg  if rms_sub_bg  is not None else np.zeros(n_spec)
            sky_orig = np.sqrt(np.maximum(0.0, rms_orig_diag**2 - noise_orig**2))
            sky_sub  = np.sqrt(np.maximum(0.0, rms_sub_diag**2  - noise_sub**2))
            with np.errstate(invalid='ignore', divide='ignore'):
                ratio = np.where(sky_orig > 0, sky_sub / sky_orig, np.nan)
            good = ratio[np.isfinite(ratio) & (ratio >= 0)]
            all_hf_ratio[win_lbl].extend(good.tolist())

            # Figure 3 trace
            hover = [f'spec {i}<br>ratio = {ratio[i]:.3f}<br>'
                     f'sky orig = {_fmt(sky_orig[i])}<br>'
                     f'sky sub  = {_fmt(sky_sub[i])}'
                     for i in range(n_spec)]
            fig3.add_trace(
                go.Scatter(
                    x=spec_x, y=ratio,
                    mode='markers',
                    marker=dict(color=c_med, size=4),
                    name=label, legendgroup=label,
                    showlegend=(col_idx == 1),
                    text=hover,
                    hovertemplate='%{text}<extra></extra>',
                ),
                row=1, col=col_idx,
            )

            # Collect statistics for screen + HTML output
            noise_med    = float(np.nanmedian(noise_sub))
            sky_orig_med = float(np.nanmedian(sky_orig))
            sky_sub_med  = float(np.nanmedian(sky_sub))
            if len(good) > 3:
                p50, p90, p95 = np.nanpercentile(good, [50, 90, 95])
                frac_lt_half  = float(np.mean(good < 0.5)) * 100
            else:
                p50 = p90 = p95 = frac_lt_half = float('nan')
            n_bg = int(bg_sel.sum())
            _stats_rows.append(dict(
                file=label, window=win_lbl,
                n_spec=int(n_spec), n_valid=len(good), n_bg_pix=n_bg,
                noise_med=noise_med,
                sky_orig_med=sky_orig_med, sky_sub_med=sky_sub_med,
                p50=p50, p90=p90, p95=p95, frac_lt_half=frac_lt_half,
            ))

            # ── Figure 4: median spectra in search range ──────
            # Use full-resolution arrays for the search range
            s4_sel = (wave >= slo) & (wave <= shi)
            if not s4_sel.any():
                continue
            w4      = wave[s4_sel]
            orig4   = (flux + sky)[:, s4_sel]
            sub4    = flux[:, s4_sel]

            for arr4, dash in [(orig4, 'dot'), (sub4, 'solid')]:
                p10_4 = np.nanpercentile(arr4, 10, axis=0)
                p90_4 = np.nanpercentile(arr4, 90, axis=0)
                med4  = np.nanmedian(arr4, axis=0)
                fig4.add_trace(go.Scatter(
                    x=w4, y=p10_4, mode='lines', line=dict(width=0),
                    showlegend=False, hoverinfo='skip'), row=1, col=col_idx)
                fig4.add_trace(go.Scatter(
                    x=w4, y=p90_4, mode='lines', line=dict(width=0),
                    fill='tonexty', fillcolor=c_band,
                    showlegend=False, hoverinfo='skip'), row=1, col=col_idx)
                fig4.add_trace(go.Scatter(
                    x=w4, y=med4, mode='lines',
                    line=dict(color=c_med, width=1.2, dash=dash),
                    showlegend=False, hoverinfo='skip'), row=1, col=col_idx)

            # Add shading on first file only (shapes are per-figure)
            if not _fig4_shaded:
                # diagnostic window — red/pink
                fig4.add_vrect(x0=win_lo, x1=win_hi,
                               fillcolor='rgba(220,50,50,0.30)',
                               line_width=1, line_color='rgba(180,0,0,0.5)',
                               row=1, col=col_idx)
                # background: mask-selected clean pixels outside diagnostic window
                if clean_mask_full is not None:
                    bg_full = (wave >= slo) & (wave <= shi) & \
                              ~((wave >= win_lo) & (wave <= win_hi)) & \
                              clean_mask_full
                else:
                    bg_full = (wave >= slo) & (wave <= shi) & \
                              ~((wave >= win_lo) & (wave <= win_hi))
                # shade contiguous clean runs
                idx = np.where(np.diff(np.concatenate([[False], bg_full, [False]])))[0]
                for r0, r1 in zip(idx[::2], idx[1::2]):
                    fig4.add_vrect(x0=wave[r0], x1=wave[r1 - 1],
                                   fillcolor='rgba(30,160,30,0.25)',
                                   line_width=0, row=1, col=col_idx)

        _fig4_shaded = True   # only add shapes once

    # ── Post-loop axis clipping ────────────────────────────────
    # Figure 3: linear y-range 0 to 99th-pct + margin; reference line at 1
    for col_idx, (win_lbl, *_) in enumerate(_DIAG_WINDOWS, start=1):
        vals = all_hf_ratio[win_lbl]
        if vals:
            v = np.array(vals)
            v = v[np.isfinite(v) & (v >= 0)]
            if len(v):
                y_hi = float(np.nanpercentile(v, 99)) + 0.05
                fig3.update_yaxes(range=[-0.05, y_hi], row=1, col=col_idx)
        # reference line at ratio = 1 (no improvement)
        fig3.add_hline(y=1.0, line_dash='dot', line_color='grey',
                       line_width=1, row=1, col=col_idx)

    # Figure 6: fixed y-range (_CONT_TREND_YLIM) on every panel, the same
    # for every row/arm/file, so the fraction of exposures actually visible
    # is comparable across overlaid files (the out-of-range percentage is
    # reported directly by _add_trend_panel's annotation).  Reference line
    # at 0 (no continuum bias).  Row order matches Figure 5: SCI, SKY,
    # Subtracted.
    for row_idx in (1, 2, 3):
        for col_idx in (1, 2, 3):
            fig6.update_yaxes(range=list(_CONT_TREND_YLIM), row=row_idx, col=col_idx)
            fig6.add_hline(y=0.0, line_dash='dot', line_color='grey',
                           line_width=1, row=row_idx, col=col_idx)

    # ── Figure 1 layout ───────────────────────────────────────
    _lin_range = [-1e-14, 1e-13]   # common linear scale for all three panels

    _title_stem = Path(filenames[0]).stem if len(filenames) == 1 \
                  else f'{len(filenames)} files'
    fig1.update_layout(
        title=f'Sky Subtraction Evaluation — {_title_stem}  {int(wmin)}–{int(wmax)} Å',
        height=1575, template='simple_white',
        legend  = dict(**_leg_style, y=0.99),
        legend2 = dict(**_leg_style, y=0.66),
        legend3 = dict(**_leg_style, y=0.33),
    )
    fig1.update_xaxes(**_ax)
    fig1.update_yaxes(**_ax, exponentformat='e', showexponent='all')
    fig1.update_xaxes(range=[wmin, wmax])
    fig1.update_xaxes(showticklabels=False, row=1, col=1)
    fig1.update_xaxes(showticklabels=False, row=2, col=1)
    fig1.update_xaxes(title_text='Wavelength (Å)', row=3, col=1)
    fig1.update_yaxes(title_text='Flux',         range=_lin_range, row=1, col=1)
    fig1.update_yaxes(title_text='Sky-sub Flux', range=_lin_range, row=2, col=1)
    fig1.update_yaxes(title_text='Sky',          range=_lin_range, row=3, col=1)

    # ── Figure 2 layout ───────────────────────────────────────
    fig2.update_layout(
        height=480, template='simple_white',
        barmode='overlay',
        legend=dict(**_leg_style, y=0.99),
    )
    fig2.update_xaxes(**_ax, exponentformat='e', showexponent='all',
                      title_text='HF residual flux', range=list(_LINE_HIST_XLIM))
    fig2.update_yaxes(**_ax, exponentformat='none', title_text='N')

    # ── Figure 3 layout ───────────────────────────────────────
    max_n = max(all_n_spec) if all_n_spec else 1
    # Use MJD only if every file plotted it that way; a mix would put
    # index- and MJD-based traces on the same numeric axis.
    use_mjd_axis = bool(all_used_mjd) and all(all_used_mjd)
    fig3.update_layout(
        height=420, template='simple_white',
        legend=dict(**_leg_style, y=0.99),
    )
    if use_mjd_axis:
        mjd_all = np.concatenate(all_mjd_vals)
        mjd_lo, mjd_hi = float(np.min(mjd_all)), float(np.max(mjd_all))
        pad = 0.01 * (mjd_hi - mjd_lo) if mjd_hi > mjd_lo else 0.5
        fig3.update_xaxes(**_ax, title_text='MJD', matches='x',
                          exponentformat='none', tickformat='.1f',
                          range=[mjd_lo - pad, mjd_hi + pad])
    else:
        fig3.update_xaxes(**_ax, title_text='Spectrum index',
                          range=[-0.5, max_n - 0.5], matches='x',
                          exponentformat='none')
    fig3.update_yaxes(**_ax, title_text='HF RMS ratio (sub/orig)')

    # ── Figure 4 layout ───────────────────────────────────────
    fig4.update_layout(
        height=420, template='simple_white',
        showlegend=False,
    )
    fig4.update_xaxes(**_ax, title_text='Wavelength (Å)')
    fig4.update_yaxes(**_ax, exponentformat='e', showexponent='all',
                      title_text='Flux')

    # ── Figure 5 layout ───────────────────────────────────────
    # No floating legend (it would sit on top of the top-right subplot in a
    # 3x3 grid); a colour-coded suptitle above the whole grid replaces it.
    fig5.update_layout(
        height=1450, template='simple_white',
        barmode='overlay',
        showlegend=False,
        margin=dict(t=110),
        title=dict(text=_suptitle_text(_file_legend), x=0.5, xanchor='center',
                  y=0.99, yanchor='top', font=dict(size=13)),
    )
    fig5.update_xaxes(**_ax, exponentformat='e', showexponent='all',
                     range=list(_CONT_HIST_XLIM))
    fig5.update_yaxes(**_ax, exponentformat='none', title_text='N')
    fig5.update_xaxes(title_text='SCI continuum residual (per-spectrum median)', row=1)
    fig5.update_xaxes(title_text='SKY continuum residual (per-spectrum median)', row=2)
    fig5.update_xaxes(title_text='Sky-subtracted flux (clean pixels)', row=3)

    # ── Figure 6 layout ───────────────────────────────────────
    fig6.update_layout(
        height=1250, template='simple_white',
        showlegend=False,
        margin=dict(t=110),
        title=dict(text=_suptitle_text(_file_legend), x=0.5, xanchor='center',
                  y=0.99, yanchor='top', font=dict(size=13)),
    )
    if use_mjd_axis:
        fig6.update_xaxes(**_ax, title_text='MJD', matches='x',
                          exponentformat='none', tickformat='.1f',
                          range=[mjd_lo - pad, mjd_hi + pad])
    else:
        fig6.update_xaxes(**_ax, title_text='Spectrum index',
                          range=[-0.5, max_n - 0.5], matches='x',
                          exponentformat='none')
    fig6.update_yaxes(**_ax, exponentformat='e', showexponent='all')
    fig6.update_yaxes(title_text='SCI continuum residual', row=1)
    fig6.update_yaxes(title_text='SKY continuum residual', row=2)
    fig6.update_yaxes(title_text='Sky-subtracted flux residual', row=3)

    # ── Print statistics to screen and build HTML stats block ─
    _hdr = (f'{"File":<55} {"Window":<10} {"N":>5} {"Noise":>10} '
            f'{"Sky orig":>10} {"Sky sub":>10} {"Med":>6} '
            f'{"p90":>6} {"p95":>6} {"<50%":>6}')
    print()
    print(_hdr)
    print('-' * len(_hdr))
    _html_rows = []
    for r in _stats_rows:
        line = (f'{r["file"]:<55} {r["window"]:<10} {r["n_valid"]:>5} '
                f'{r["noise_med"]:>10.2e} {r["sky_orig_med"]:>10.2e} '
                f'{r["sky_sub_med"]:>10.2e} {r["p50"]:>6.3f} '
                f'{r["p90"]:>6.3f} {r["p95"]:>6.3f} {r["frac_lt_half"]:>5.0f}%')
        print(line)
        _html_rows.append(
            f'<tr><td>{r["file"]}</td><td>{r["window"]}</td>'
            f'<td>{r["n_valid"]}</td>'
            f'<td>{r["noise_med"]:.2e}</td>'
            f'<td>{r["sky_orig_med"]:.2e}</td>'
            f'<td>{r["sky_sub_med"]:.2e}</td>'
            f'<td>{r["p50"]:.3f}</td><td>{r["p90"]:.3f}</td>'
            f'<td>{r["p95"]:.3f}</td><td>{r["frac_lt_half"]:.0f}%</td></tr>'
        )
    print()

    _stats_html = '''
<div style="font-family:monospace; font-size:13px; margin:10px 20px;">
<h3>Sky subtraction diagnostics — HF RMS ratio (sky-subtracted / original)</h3>
<p>Noise: median per-spectrum background RMS (clean pixels near each window).<br>
Sky orig/sub: median noise-corrected sky-line RMS before/after subtraction.<br>
Ratio columns: 50th (median), 90th, 95th percentile across spectra.<br>
&lt;50%: fraction of spectra where the ratio is below 0.5.</p>
<table border="1" cellpadding="4" cellspacing="0" style="border-collapse:collapse;">
<tr style="background:#ddd;">
  <th>File</th><th>Window</th><th>N valid</th>
  <th>Noise</th><th>Sky orig</th><th>Sky sub</th>
  <th>p50</th><th>p90</th><th>p95</th><th>&lt;50%</th>
</tr>
''' + '\n'.join(_html_rows) + '\n</table></div>\n'

    _kind_order = {'SCI': 0, 'SKY': 1, 'Subtracted': 2}
    _arm_stats_rows.sort(key=lambda r: (r['file'], _kind_order.get(r['kind'], 9), r['arm']))

    _arm_hdr = (f'{"File":<55} {"Kind":<11} {"Arm":<5} {"N":>8} '
               f'{"Med":>10} {"NMAD":>10} {"Skew":>6}')
    print()
    print('Continuum separation quality (SCI/SKY vs raw pre-subtraction spectra; '
         'Subtracted = post-subtraction leftover):')
    print(_arm_hdr)
    print('-' * len(_arm_hdr))
    _arm_html_rows = []
    for r in _arm_stats_rows:
        line = (f'{r["file"]:<55} {r["kind"]:<11} {r["arm"]:<5} {r["n"]:>8} '
                f'{_fmt(r["med"]):>10} {_fmt(r["nmad"]):>10} {r["skew"]:>6.2f}')
        print(line)
        _arm_html_rows.append(
            f'<tr><td>{r["file"]}</td><td>{r["kind"]}</td><td>{r["arm"]}</td>'
            f'<td>{r["n"]}</td>'
            f'<td>{_fmt(r["med"])}</td><td>{_fmt(r["nmad"])}</td>'
            f'<td>{r["skew"]:.2f}</td></tr>'
        )
    print()

    _arm_stats_html = '''
<div style="font-family:monospace; font-size:13px; margin:10px 20px;">
<h3>Continuum separation quality — per file, per arm</h3>
<p><b>Kind = SCI</b>: science-side continuum-fit quality, against the raw
pre-subtraction science spectrum (DRP_ALL['SCI_MED_&lt;arm&gt;']).<br>
<b>Kind = SKY</b>: sky-side continuum-fit quality, against the raw
pre-subtraction sky-telescope spectrum whose fit was used in SKY
(DRP_ALL['SKY_MED_&lt;arm&gt;']).<br>
<b>Kind = Subtracted</b>: per-spectrum median of the raw sky-subtracted flux
in clean (sky-line-free) pixels -- the final post-subtraction leftover
signal, <i>not</i> a continuum-fit-quality test (see Figure 5 caption).</p>
<p>A median near zero for Kind=SCI/SKY indicates that side's continuum fit
tracked the true continuum well in that arm. Per-spectrum values of all
three are also written to this file's own DRP_ALL table for detailed
analysis (resid_med_&lt;arm&gt; etc. for Subtracted; SCI_MED_&lt;arm&gt;/
SKY_MED_&lt;arm&gt; etc. already present in the input, for SCI/SKY).</p>
<table border="1" cellpadding="4" cellspacing="0" style="border-collapse:collapse;">
<tr style="background:#ddd;">
  <th>File</th><th>Kind</th><th>Arm</th><th>N</th><th>Med</th><th>NMAD</th><th>Skew</th>
</tr>
''' + '\n'.join(_arm_html_rows) + '\n</table></div>\n'

    _fig2_html_header = '''
<div style="margin:10px 20px;">
<b>Residual histograms</b> — distribution of per-pixel HF (high-frequency)
residual flux within each diagnostic sky-line window, pooled across all
spectra. Each spectrum has a Gaussian-smoothed, clean-pixel-weighted
continuum estimate subtracted first (see Notes), so the histogram
reflects leftover sky-line-scale structure rather than continuum level.
A distribution centred on zero with small NMAD indicates good sky-line
removal; the dotted curve is a Gaussian with the same median and NMAD
for reference.
</div>
'''

    _fig4_html_header = '''
<div style="margin:10px 20px;">
<b>Diagnostic window spectra</b> — median with 10–90% band.
Dotted = original flux (FLUX+SKY); solid = sky-subtracted (FLUX).
Each file shown with its own colour.<br>
<span style="background:rgba(220,50,50,0.35); padding:2px 8px;">&#9632;</span>
Signal window (used for sky-line RMS)&nbsp;&nbsp;
<span style="background:rgba(30,160,30,0.30); padding:2px 8px;">&#9632;</span>
Background (mask-selected clean pixels, used for noise RMS).
</div>
'''

    _fig5_html_header = '''
<div style="margin:10px 20px;">
<b>Continuum-quality histograms</b> — three rows, one per diagnostic:<br>
<b>Row 1 (SCI)</b> — distribution of DRP_ALL['SCI_MED_&lt;arm&gt;'] across
spectra: the actual science-side continuum-fit-quality metric, evaluated
against the raw pre-subtraction science spectrum.<br>
<b>Row 2 (SKY)</b> — distribution of DRP_ALL['SKY_MED_&lt;arm&gt;'] across
spectra: the actual sky-side continuum-fit-quality metric, evaluated against
the raw pre-subtraction sky-telescope spectrum whose continuum fit was used
in SKY.<br>
<b>Row 3 (Subtracted)</b> — distribution of the per-spectrum median of raw
sky-subtracted FLUX in clean (sky-line-free) pixels, per spectrograph arm
(overlap zones and outer edges excluded) — one value per spectrum, same as
rows 1/2, so all three rows are directly comparable. This is the final
post-subtraction leftover signal — real source continuum entangled with any
net sky-continuum error, <i>not</i> a test of continuum-fit quality (the
science-side continuum fit never enters the subtracted result).<br>
A distribution centred on zero in rows 1/2 indicates that side's continuum
fit tracked the true continuum well in that arm. Rows 1/2 require
SkySubOrig.py/SkySubDev1.py/SkySepESO.py output (blank otherwise).
</div>
'''

    _fig6_html_header = '''
<div style="margin:10px 20px;">
<b>Per-arm continuum residual per spectrum</b> — same three rows and same
quantities as the histograms above (SCI, SKY, Subtracted), but plotted
per-spectrum against spectrum index or MJD instead of pooled into a
histogram, so specific bad exposures/fibers are visible rather than just an
aggregate number. Rows 1/2 require SkySubOrig.py/SkySubDev1.py/SkySepESO.py
output. For row 3 (Subtracted), a nonzero value can also reflect real
source continuum rather than a sky-subtraction defect -- comparing this
row across overlaid methods on the same input isolates sky-continuum
quality specifically, since real source continuum is identical across
methods.
</div>
'''

    # ── Write combined HTML ───────────────────────────────────
    if outroot == '':
        outroot = Path(filenames[0]).stem + '_eval'
    outfile = f'{outroot}.html'

    html1 = fig1.to_html(full_html=False, include_plotlyjs=True)
    html2 = fig2.to_html(full_html=False, include_plotlyjs=False)
    html3 = fig3.to_html(full_html=False, include_plotlyjs=False)
    html4 = fig4.to_html(full_html=False, include_plotlyjs=False)
    html5 = fig5.to_html(full_html=False, include_plotlyjs=False)
    html6 = fig6.to_html(full_html=False, include_plotlyjs=False)

    _page_title = 'Sky Subtraction Quality Check'
    _page_header = f'''
<div style="margin:20px 20px 10px 20px;">
<h1 style="margin-bottom:4px;">{_page_title}</h1>
</div>
'''
    _section_header_sky = '''
<div style="margin:30px 20px 10px 20px; border-bottom:3px solid #333;">
<h2 style="margin-bottom:6px;">Sky Line Subtraction</h2>
</div>
'''
    _section_header_cont = '''
<div style="margin:50px 20px 10px 20px; border-bottom:3px solid #333;">
<h2 style="margin-bottom:6px;">Continuum Separation</h2>
</div>
'''

    with open(outfile, 'w') as fh:
        fh.write(f'<!DOCTYPE html>\n<html>\n<head><title>{_page_title}</title></head>\n<body>\n')
        fh.write(_page_header)
        fh.write(html1)
        fh.write('\n')
        fh.write(_section_header_sky)
        fh.write(_fig2_html_header)
        fh.write(html2)
        fh.write('\n')
        fh.write(_fig4_html_header)
        fh.write(html4)
        fh.write('\n')
        fh.write(_stats_html)
        fh.write(html3)
        fh.write('\n')
        fh.write(_section_header_cont)
        fh.write(_fig5_html_header)
        fh.write(html5)
        fh.write('\n')
        fh.write(_arm_stats_html)
        fh.write(_fig6_html_header)
        fh.write(html6)
        fh.write('\n</body>\n</html>\n')

    print(f'Wrote {outfile}')


# ──────────────────────────────────────────────────────────────
# Command-line entry point
# ──────────────────────────────────────────────────────────────

def steer(argv):
    filenames       = []
    positional_nums = []
    n_sample        = 20
    outroot         = ''

    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg in ('-h', '--help'):
            print(_USAGE)
            return
        elif arg == '-num':
            i += 1
            n_sample = int(argv[i])
        elif arg == '-out':
            i += 1
            outroot = argv[i]
        elif arg.startswith('-'):
            print('Error: unknown option "%s"' % arg)
            print(_USAGE)
            return
        else:
            try:
                positional_nums.append(float(arg))
            except ValueError:
                filenames.append(arg)
        i += 1

    if not filenames:
        print(_USAGE)
        return

    wmin = positional_nums[0] if len(positional_nums) > 0 else 3600.0
    wmax = positional_nums[1] if len(positional_nums) > 1 else 9800.0

    if outroot:
        # Explicit output root: overlay all files in one HTML
        plot_eval(filenames, wmin=wmin, wmax=wmax, n_sample=n_sample, outroot=outroot)
    else:
        # Default: one HTML per file
        for filename in filenames:
            plot_eval(filename, wmin=wmin, wmax=wmax, n_sample=n_sample, outroot='')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
