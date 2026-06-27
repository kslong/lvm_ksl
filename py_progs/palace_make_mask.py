#!/usr/bin/env python
"""
palace_clean_windows.py
=======================

Synopsis
--------
Build a sky-line mask for LVM spectra using PALACE sky emission models, and
identify the wavelength windows that are clean enough to use for continuum
measurement.

Description
-----------
Night-sky emission in the LVM wavelength range (3600–9800 Å) is dominated by
OH vibrational-rotational lines in the near-infrared (Z arm, 7400–9800 Å),
atomic forbidden lines ([NI], OI, NaI, KI) in the blue and red, and the O2
atmospheric band near 8650 Å.  These features contaminate any attempt to fit
a smooth continuum through the spectrum.

This program builds a model of sky line emission using the PALACE (Paranal
Airglow Line And Continuum Emission) reference data set (Noll et al. 2024),
renders each sky line component onto the LVM wavelength grid with the
instrument LSF, scales the combined model to match the observed sky spectrum,
and produces a boolean mask that is True (good) at pixels where the predicted
sky line contamination is below a user-specified flux threshold.

Sky line model
--------------
Four PALACE components are used:

  OH   — hydroxyl vibrational-rotational bands from pmd_popmodel_OH.dat.
          Line amplitudes are proportional to Einstein Aij × degeneracy gi,
          grouped by (v_upper, N_upper, F_upper) quantum numbers.

  OI   — oxygen recombination lines from pmd_intmodel_Orc.dat.
          Two multiplet groups (OI 7774 Å and OI 8446 Å).

  atom — atomic forbidden/permitted lines (NaI, KI, [NI], OI green/red) from
          pmd_intdata_atom.dat.  Intensities from the PALACE climatology in
          Rayleigh units.  Hydrogen recombination and OI recombination lines
          are excluded (handled separately above).

  O2   — molecular oxygen A-band near 8650 Å from pmd_popmodel_O2.dat,
          computed for a rotational temperature of 191.5 K.

Each component is normalised to its own peak before being summed into the
combined model.  This prevents the very bright OH lines from completely
suppressing the atomic lines in the mask.

The combined model is then scaled by a single least-squares factor so that
it best matches the observed sky spectrum at bright OH pixels in the Z arm
(7400–9800 Å), where the scaling is most reliable.  The resulting scaled
model is in the same flux units as the sky spectrum (factor × flux, where
factor = 1e14 by default so that typical sky flux values are of order unity).

Threshold and tradeoff
----------------------
A pixel is flagged as clean (mask = True) when:

    palace_model_scaled[i] < threshold

The threshold is expressed in the same scaled flux units as the sky spectrum.
A lower threshold means stricter masking (fewer clean pixels but lower
contamination), while a higher threshold is more permissive (more clean pixels
but more residual sky line flux).  The default of 0.01 in FACTOR=1e14 units
corresponds to a sky line contamination level of 1e-16 erg s⁻¹ cm⁻² Å⁻¹,
which is well below the typical Z-arm continuum (~7e-15 erg s⁻¹ cm⁻² Å⁻¹).

Output FITS file
----------------
The output FITS file (written by default to palace_mask_<input_stem>.fits)
contains four extensions:

  WAVE       float32  (Npix,)   Wavelength array in Å (air wavelengths)
  SKY        float32  (Npix,)   Input sky spectrum (median over all fibers,
                                scaled by FACTOR)
  CONTINUUM  float32  (Npix,)   Sky spectrum with sky-line-contaminated pixels
                                set to NaN; these are the pixels that can be
                                used for continuum fitting
  MASK       uint8    (Npix,)   Boolean mask: 1 = clean, 0 = contaminated

The primary HDU header records THRESH (threshold), SCALE (model scale factor),
FACTOR (flux scale factor), and SKYEXT (source sky extension name).

Usage
-----
    python palace_clean_windows.py <fits_file> <palace_dir> [options]

    # Basic run with default threshold and output file
    python palace_clean_windows.py XCframe_1.2.1_7325_48860_1_50.fits palace/PMD

    # Adjust threshold and show diagnostic plot
    python palace_clean_windows.py XCframe_1.2.1_7325_48860_1_50.fits palace/PMD \\
        --threshold 0.05 --plot

    # Save to a specific output file
    python palace_clean_windows.py XCframe_1.2.1_7325_48860_1_50.fits palace/PMD \\
        --threshold 0.01 --output my_mask.fits --plot

Arguments
---------
    fits_file    LVM XCframe FITS file.  Must contain WAVE, a sky spectrum
                 extension (default SKY_EAST), and optionally LSF.
    palace_dir   Path to the palace/PMD directory containing the PALACE data
                 files (pmd_popmodel_OH.dat, pmd_intdata_atom.dat, etc.).

Options
-------
    --sky-ext NAME     FITS extension name for the sky spectrum.
                       (default: SKY_EAST)
    --threshold T      Sky line contamination threshold in FACTOR*flux units.
                       Pixels where the scaled PALACE model exceeds this value
                       are masked.  Lower = stricter mask, fewer clean pixels.
                       (default: 0.01)
    --factor F         Multiplicative factor applied to the raw flux values
                       when reading from the FITS file.  Brings raw cgs flux
                       (~7e-15 erg s⁻¹ cm⁻² Å⁻¹) into convenient numeric
                       range (~0.7 in FACTOR units).
                       (default: 1e14)
    --no-lsf           If set, use a fixed Gaussian LSF sigma (--lsf-sigma)
                       instead of reading the per-pixel LSF from the FITS file.
    --lsf-sigma S      Fixed Gaussian LSF sigma in Å.  Only used when
                       --no-lsf is specified.
                       (default: 0.65 Å ≈ FWHM 1.5 Å)
    --output PATH      Write output FITS to PATH.  If not given, defaults to
                       palace_mask_<stem>.fits in the current directory, where
                       <stem> is the input filename without extension.
    --min-window W     Minimum clean window width in Å to include in the
                       printed window table.
                       (default: 5 Å)
    --plot             Display a three-panel diagnostic plot (one per arm)
                       showing the sky spectrum, PALACE model, threshold, and
                       the clean continuum.

Notes
-----
History
    260627  Written (Knox Long / Claude)
"""

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table


FACTOR_DEFAULT   = 1.0e14
LSF_SIGMA_DEFAULT = 0.65       # Å  (≈ FWHM 1.5 Å / 2.355)
THRESHOLD_DEFAULT = 0.01
CAP = 5.0                      # Å padding when loading line catalogues

ARM_RANGES = {
    "B": (3600.0, 5900.0),
    "R": (5900.0, 7400.0),
    "Z": (7400.0, 9800.0),
}

# Z arm wavelength used to anchor model-to-sky scale factor
Z_SCALE_MIN = 7400.0


# ---------------------------------------------------------------------------
# Wavelength conversion
# ---------------------------------------------------------------------------

def vac_to_air(lam_vac_a):
    """Convert vacuum wavelengths (Å) to air wavelengths using the
    Edlén (1966) formula as parameterised by Ciddor (1996).

    Parameters
    ----------
    lam_vac_a : array-like
        Vacuum wavelengths in Å.

    Returns
    -------
    ndarray
        Air wavelengths in Å.
    """
    lam = np.asarray(lam_vac_a, float)
    s2  = (1e4 / lam) ** 2
    n   = 1.0 + 8.34254e-5 + 2.406147e-2 / (130.0 - s2) + 1.5998e-4 / (38.9 - s2)
    return lam / n


# ---------------------------------------------------------------------------
# HITRAN ID decoder (internal helper)
# ---------------------------------------------------------------------------

def _decode_hitran_id(table):
    """Parse PALACE/HITRAN line ID strings and add quantum number columns.

    The ID format encodes (in positions 4–12): v_upper, v_lower, branch_N,
    branch_J, F_upper, F_lower, N_lower, parity.  The symbol 'X' for v_upper
    denotes level 10 in PMD notation.

    Adds columns v_upper, v_lower, N_upper, N_lower, F_upper, F_lower,
    branch_N to the input Table and returns it.
    """
    ids_str  = np.asarray(table["ID"].astype(str), dtype=str)
    if np.min(np.char.str_len(ids_str)) < 13:
        raise ValueError("HITRAN ID strings shorter than expected (need ≥13 chars).")
    v_up_str = np.array([s[4:5]   for s in ids_str])
    v_low    = np.array([s[5:6]   for s in ids_str], dtype=int)
    branch_n = np.array([s[6:7]   for s in ids_str])
    f_up     = np.array([s[8:9]   for s in ids_str], dtype=int)
    f_low    = np.array([s[9:10]  for s in ids_str], dtype=int)
    n_low    = np.array([s[10:12] for s in ids_str], dtype=int)
    v_up_str = np.where(v_up_str == "X", "10", v_up_str)
    v_up     = v_up_str.astype(int)
    delta_map = {"O": -2, "P": -1, "Q": 0, "R": 1, "S": 2}
    d_n      = np.vectorize(delta_map.get)(branch_n)
    n_up     = n_low + d_n
    table["v_upper"]  = v_up
    table["v_lower"]  = v_low
    table["N_upper"]  = n_up
    table["N_lower"]  = n_low
    table["F_upper"]  = f_up
    table["F_lower"]  = f_low
    table["branch_N"] = branch_n
    return table


# ---------------------------------------------------------------------------
# Line rendering
# ---------------------------------------------------------------------------

def _grp2vector(line_wave, line_amp, wave, lsf_sigma):
    """Render a set of emission lines onto the wavelength grid.

    Each line is represented as a Gaussian with sigma given by lsf_sigma.
    Lines are summed into a single 1-D array on the output grid.

    Parameters
    ----------
    line_wave : array-like
        Central wavelengths of the lines (Å).
    line_amp : array-like
        Amplitude of each line (arbitrary units).
    wave : ndarray
        Output wavelength grid (Å).
    lsf_sigma : float or ndarray
        LSF sigma in Å.  If a scalar, all lines use the same sigma.  If an
        array with the same length as `wave`, the sigma at each line centre
        is interpolated from the grid (wavelength-dependent LSF).

    Returns
    -------
    ndarray
        1-D array of shape (len(wave),) containing the rendered line profile.
    """
    cent = np.asarray(line_wave, float)
    amp  = np.asarray(line_amp,  float)
    if np.ndim(lsf_sigma) > 0:
        sig = np.interp(cent, wave, lsf_sigma)
        yy  = (wave[:, None] - cent[None, :]) / sig[None, :]
    else:
        sig = float(lsf_sigma)
        yy  = (wave[:, None] - cent[None, :]) / sig
    return np.sum(amp[None, :] * np.exp(-0.5 * yy ** 2), axis=1)


# ---------------------------------------------------------------------------
# Per-component PALACE model builders (internal)
# ---------------------------------------------------------------------------

def _build_oh(wave, palace_dir, lsf_sigma, cap):
    """Render OH vibrational-rotational line groups from pmd_popmodel_OH.dat."""
    oh = Table.read(str(palace_dir / "pmd_popmodel_OH.dat"),
                    format="ascii.basic", guess=False, comment="#", fast_reader=False)
    oh["wave"] = vac_to_air(np.asarray(oh["lam"], float) * 1e4)
    oh = oh[(oh["wave"] >= wave.min() - cap) & (oh["wave"] <= wave.max() + cap)]
    oh = _decode_hitran_id(oh)
    model = np.zeros(wave.size)
    for grp in oh.group_by(("v_upper", "N_upper", "F_upper")).groups:
        amp = np.asarray(grp["Aij"], float) * np.asarray(grp["gi"], float)
        model += _grp2vector(grp["wave"], amp, wave, lsf_sigma)
    return model


def _build_orc(wave, palace_dir, lsf_sigma, cap):
    """Render OI recombination line multiplets from pmd_intmodel_Orc.dat."""
    orc = Table.read(str(palace_dir / "pmd_intmodel_Orc.dat"),
                     format="ascii.basic", guess=False, comment="#", fast_reader=False)
    orc["wave"] = vac_to_air(np.asarray(orc["lam"], float) * 1e4)
    orc = orc[(orc["wave"] >= wave.min() - cap) & (orc["wave"] <= wave.max() + cap)]
    model = np.zeros(wave.size)
    for grp in orc.group_by("reffeat").groups:
        model += _grp2vector(grp["wave"], np.asarray(grp["I"], float), wave, lsf_sigma)
    return model


def _build_atom(wave, palace_dir, lsf_sigma, cap):
    """Render atomic sky lines (NaI, KI, [NI], OI) from pmd_intdata_atom.dat.

    Hydrogen recombination lines and OI recombination lines (class 'H' and
    'Orc') are excluded; they are either negligible or handled separately.
    Intensities are PALACE reference values in Rayleigh units.
    """
    atom = Table.read(str(palace_dir / "pmd_intdata_atom.dat"),
                      format="ascii.basic", guess=False, comment="#", fast_reader=False)
    atom["wave"] = vac_to_air(np.asarray(atom["lam"], float) * 1e4)
    atom = atom[(atom["wave"] >= wave.min() - cap) & (atom["wave"] <= wave.max() + cap)]
    atom = atom[~np.isin(np.asarray(atom["class"], str), ["H", "Orc"])]
    model = np.zeros(wave.size)
    for grp in atom.group_by("class").groups:
        model += _grp2vector(grp["wave"], np.asarray(grp["I"], float), wave, lsf_sigma)
    return model


def _build_o2(wave, palace_dir, lsf_sigma, t_o2=191.5):
    """Render the O2 A-band near 8650 Å from pmd_popmodel_O2.dat.

    Only the v_i=0 vibrational level is used.  Line amplitudes follow a
    Boltzmann population distribution at temperature t_o2.

    Parameters
    ----------
    t_o2 : float
        O2 rotational temperature in K (default 191.5 K from PALACE).
    """
    o2_min, o2_max = 8600.0, 8715.0
    hc_kB = 1.4387769   # hc / k_B in cm K
    pop = Table.read(str(palace_dir / "pmd_popmodel_O2.dat"),
                     format="ascii.basic", guess=False, comment="#", fast_reader=False)
    pop["wave"] = vac_to_air(np.asarray(pop["lam"], float) * 1e4)
    o2 = pop[
        (pop["wave"] >= o2_min) & (pop["wave"] <= o2_max) &
        (np.asarray(pop["vi"], int) == 0)
    ]
    Ei  = np.asarray(o2["Ei"], float)
    rel = (np.asarray(o2["Aij"], float) * np.asarray(o2["gi"], float) *
           np.exp(-hc_kB * (Ei - Ei.min()) / t_o2))
    rel /= rel.sum()
    return _grp2vector(np.asarray(o2["wave"], float), rel, wave, lsf_sigma)


# ---------------------------------------------------------------------------
# Combined PALACE model
# ---------------------------------------------------------------------------

def build_palace_line_model(wave, palace_dir, lsf_sigma=LSF_SIGMA_DEFAULT,
                             cap=CAP, t_o2=191.5, verbose=True):
    """Build a combined sky line contamination model from all PALACE components.

    Each component (OH, OI recombination, atomic, O2) is normalised to its
    own peak before summing, so that every sky line family contributes
    equally to the mask regardless of its absolute brightness.  The combined
    model therefore represents sky line *presence* across the spectrum rather
    than a physically scaled flux prediction.  Physical scaling is applied
    separately by scale_model_to_sky().

    Parameters
    ----------
    wave : ndarray
        1-D wavelength array (Å, air wavelengths).
    palace_dir : str or Path
        Path to the palace/PMD directory.
    lsf_sigma : float or ndarray
        Gaussian LSF sigma in Å.  A scalar uses a fixed width; an array of
        the same length as `wave` gives a wavelength-dependent LSF.
    cap : float
        Wavelength padding in Å beyond the wave grid edges when loading
        line catalogues (prevents edge truncation of broad lines).
    t_o2 : float
        O2 rotational temperature in K (default 191.5 K).
    verbose : bool
        Print progress messages while loading components.

    Returns
    -------
    combined : ndarray
        Sum of all normalised components, shape (len(wave),).
    parts : dict
        Individual normalised components keyed by 'OH', 'OI', 'atom', 'O2'.
    """
    palace_dir = Path(palace_dir)

    if verbose:
        print("  Loading OH lines ...", end=" ", flush=True)
    raw_oh   = _build_oh(  wave, palace_dir, lsf_sigma, cap)
    if verbose:
        print("OI recombination ...", end=" ", flush=True)
    raw_orc  = _build_orc( wave, palace_dir, lsf_sigma, cap)
    if verbose:
        print("atomic ...", end=" ", flush=True)
    raw_atom = _build_atom(wave, palace_dir, lsf_sigma, cap)
    if verbose:
        print("O2 band ...", flush=True)
    raw_o2   = _build_o2(  wave, palace_dir, lsf_sigma, t_o2=t_o2)

    def _norm(v):
        pk = np.nanmax(v)
        return v / pk if pk > 0 else v.copy()

    parts = {
        "OH":   _norm(raw_oh),
        "OI":   _norm(raw_orc),
        "atom": _norm(raw_atom),
        "O2":   _norm(raw_o2),
    }
    combined = sum(parts.values())
    return combined, parts


# ---------------------------------------------------------------------------
# Scaling model to observed sky
# ---------------------------------------------------------------------------

def scale_model_to_sky(wave, palace_model, sky_flux,
                        z_wave_min=Z_SCALE_MIN, bright_frac=0.10):
    """Fit a single scale factor to match the PALACE model to an observed sky.

    The fit is restricted to pixels in the Z arm (λ ≥ z_wave_min) where the
    PALACE model exceeds bright_frac × its Z-arm peak, i.e., pixels at or
    near the cores of bright OH lines.  These pixels provide the most
    reliable constraint on the overall amplitude because OH dominates the Z
    arm and the model's relative line ratios are well determined by HITRAN
    transition probabilities.

    The scale factor is found by minimising the sum of squared residuals
    between scale × model and sky over the selected pixels, giving:

        scale = (A · b) / (A · A)

    where A = model[bright_pixels] and b = sky[bright_pixels].

    Parameters
    ----------
    wave : ndarray
        Wavelength array (Å).
    palace_model : ndarray
        Combined normalised PALACE model (output of build_palace_line_model).
    sky_flux : ndarray
        Observed sky spectrum in FACTOR-scaled flux units.
    z_wave_min : float
        Blue edge of the Z arm used for fitting (Å).
    bright_frac : float
        Fraction of the Z-arm model peak used as the lower amplitude cutoff
        for selecting bright-line pixels.

    Returns
    -------
    float
        Non-negative scale factor converting PALACE model units to sky units.
    """
    z_mask = (wave >= z_wave_min) & np.isfinite(sky_flux)
    if not z_mask.any():
        return 1.0
    model_z_peak = np.nanmax(palace_model[z_mask])
    bright = z_mask & (palace_model > bright_frac * model_z_peak)
    if bright.sum() < 10:
        bright = z_mask & (palace_model > 0)
    A = palace_model[bright]
    b = sky_flux[bright]
    scale = float(np.dot(A, b) / np.dot(A, A))
    return max(scale, 0.0)


# ---------------------------------------------------------------------------
# Mask construction
# ---------------------------------------------------------------------------

def make_sky_mask(palace_model_scaled, threshold):
    """Return a boolean clean-pixel mask.

    Parameters
    ----------
    palace_model_scaled : ndarray
        PALACE sky line model in the same flux units as the observed sky
        (i.e., after applying the scale factor from scale_model_to_sky).
    threshold : float
        Contamination threshold in the same units.  Pixels where the model
        exceeds this value are flagged as contaminated (mask = False).

    Returns
    -------
    ndarray of bool
        True where predicted sky line flux < threshold (clean pixels).
    """
    return palace_model_scaled < threshold


# ---------------------------------------------------------------------------
# Window finding
# ---------------------------------------------------------------------------

def find_clean_windows(wave, mask):
    """Identify contiguous runs of clean (True) pixels.

    Parameters
    ----------
    wave : ndarray
        Wavelength array (Å).
    mask : ndarray of bool
        Clean-pixel mask (True = good).

    Returns
    -------
    list of (wave_start, wave_end, n_pix, width_A)
        Each element describes one contiguous clean window.  The list is
        sorted by decreasing window width.
    """
    windows = []
    in_window = False
    i_start = 0
    for i in range(len(mask)):
        if mask[i] and not in_window:
            in_window = True
            i_start = i
        elif not mask[i] and in_window:
            in_window = False
            ws, we = wave[i_start], wave[i - 1]
            windows.append((ws, we, i - i_start, float(we - ws)))
    if in_window:
        ws, we = wave[i_start], wave[-1]
        windows.append((ws, we, len(mask) - i_start, float(we - ws)))
    return sorted(windows, key=lambda x: -x[3])


# ---------------------------------------------------------------------------
# Console reporting
# ---------------------------------------------------------------------------

def report_coverage(wave, mask, arm_ranges=ARM_RANGES):
    """Print a per-arm summary of clean pixel counts and window counts."""
    print(f"\n{'Arm':<5}  {'Range (Å)':<16}  {'N total':>8}  "
          f"{'N clean':>8}  {'% clean':>8}  {'N windows':>10}")
    print("-" * 65)
    for arm, (wmin, wmax) in arm_ranges.items():
        sel     = (wave >= wmin) & (wave < wmax)
        n_tot   = int(sel.sum())
        n_clean = int((mask & sel).sum())
        pct     = 100.0 * n_clean / n_tot if n_tot else 0.0
        wins    = find_clean_windows(wave[sel], mask[sel])
        print(f"{arm:<5}  {wmin:.0f}–{wmax:.0f} Å       "
              f"{n_tot:>8d}  {n_clean:>8d}  {pct:>7.1f}%  {len(wins):>10d}")
    n_tot   = len(wave)
    n_clean = int(mask.sum())
    wins    = find_clean_windows(wave, mask)
    print(f"{'ALL':<5}  {wave.min():.0f}–{wave.max():.0f} Å  "
          f"{n_tot:>8d}  {n_clean:>8d}  {100*n_clean/n_tot:>7.1f}%  {len(wins):>10d}")


def report_windows(wave, mask, arm_ranges=ARM_RANGES, min_width=5.0, n_top=25):
    """Print a table of the largest clean windows per arm."""
    for arm, (wmin, wmax) in arm_ranges.items():
        sel  = (wave >= wmin) & (wave < wmax)
        wins = [w for w in find_clean_windows(wave[sel], mask[sel])
                if w[3] >= min_width]
        print(f"\nClean windows ≥ {min_width:.0f} Å in {arm} arm "
              f"({wmin:.0f}–{wmax:.0f} Å) — {len(wins)} found:")
        if not wins:
            print("  (none)")
            continue
        print(f"  {'Start':>8}  {'End':>8}  {'Width':>8}  {'N pix':>6}")
        for ws, we, np_, dw in wins[:n_top]:
            print(f"  {ws:>8.1f}  {we:>8.1f}  {dw:>7.1f} Å  {np_:>6d}")
        if len(wins) > n_top:
            print(f"  ... ({len(wins) - n_top} more not shown)")


# ---------------------------------------------------------------------------
# Output FITS
# ---------------------------------------------------------------------------

def save_output(outpath, wave, sky_flux, mask, palace_model_scaled,
                threshold, scale, factor, sky_ext):
    """Write the mask and associated spectra to a FITS file.

    Extensions
    ----------
    PRIMARY    No data.  Header records the run parameters.
    WAVE       float32  Wavelength array (Å, air).
    SKY        float32  Input sky spectrum (FACTOR-scaled flux), the median
                        over all fibers from the sky extension of the input
                        FITS file.
    CONTINUUM  float32  SKY spectrum with sky-line-contaminated pixels set
                        to NaN.  These are the pixels available for continuum
                        measurement.
    MASK       uint8    Clean-pixel mask: 1 = clean (usable), 0 = masked
                        (sky line contaminated).

    Parameters
    ----------
    outpath : str or Path
        Output file path.
    wave : ndarray
        Wavelength array (Å).
    sky_flux : ndarray
        Input sky spectrum in FACTOR-scaled units.
    mask : ndarray of bool
        Clean-pixel mask (True = clean).
    palace_model_scaled : ndarray
        PALACE model scaled to sky units (stored in SKYMODEL extension).
    threshold : float
        Contamination threshold used to build the mask.
    scale : float
        PALACE-to-sky scale factor.
    factor : float
        Flux scale factor applied when reading the input spectra.
    sky_ext : str
        Name of the sky spectrum extension in the input FITS file.
    """
    continuum = np.where(mask, sky_flux, np.nan)

    hdr = fits.Header()
    hdr["THRESH"]  = (threshold, "sky line contamination threshold (FACTOR units)")
    hdr["SCALE"]   = (scale,     "PALACE model to sky flux scale factor")
    hdr["FACTOR"]  = (factor,    "flux scale factor applied on read")
    hdr["SKYEXT"]  = (sky_ext,   "sky spectrum source extension")

    fits.HDUList([
        fits.PrimaryHDU(header=hdr),
        fits.ImageHDU(data=wave.astype(np.float32),                  name="WAVE"),
        fits.ImageHDU(data=sky_flux.astype(np.float32),              name="SKY"),
        fits.ImageHDU(data=continuum.astype(np.float32),             name="CONTINUUM"),
        fits.ImageHDU(data=mask.astype(np.uint8),                    name="MASK"),
    ]).writeto(outpath, overwrite=True)
    print(f"\nOutput written to {outpath}")


# ---------------------------------------------------------------------------
# Diagnostic plot
# ---------------------------------------------------------------------------

def make_plot(wave, sky_flux, palace_model_scaled, mask, threshold,
              arm_ranges=ARM_RANGES, title="", savepath=None, show=True):
    """Save and optionally display a three-panel diagnostic plot.

    One panel per spectrograph arm.  Each panel shows:
      - Sky spectrum (full, faded) — the median observed sky spectrum
      - Sky continuum (bold black) — sky spectrum with masked pixels blanked
        to NaN, showing only the pixels available for continuum fitting
      - PALACE model (coloured) — the scaled sky line contamination model
      - Threshold (red dashed) — the contamination cutoff level
      - Green shading — wavelength regions where the mask is True (clean)

    The y-axis is logarithmic so that both the faint inter-line continuum and
    the bright OH line peaks are simultaneously visible.

    Parameters
    ----------
    wave : ndarray
        Wavelength array (Å).
    sky_flux : ndarray
        Median sky spectrum in FACTOR-scaled flux units.
    palace_model_scaled : ndarray
        PALACE model in the same units as sky_flux.
    mask : ndarray of bool
        Clean-pixel mask (True = clean).
    threshold : float
        Contamination threshold (plotted as a horizontal red dashed line).
    arm_ranges : dict
        Mapping of arm name to (wmin, wmax) wavelength range.
    title : str
        Figure suptitle string.
    savepath : str or Path or None
        If given, save the figure to this path before displaying.
    show : bool
        If True, call plt.show() to display the figure interactively.
    """
    import matplotlib.pyplot as plt

    arm_colors = {"B": "#4e79a7", "R": "#e15759", "Z": "#59a14f"}
    n_arms = len(arm_ranges)
    fig, axes = plt.subplots(n_arms, 1, figsize=(15, 4 * n_arms))
    if n_arms == 1:
        axes = [axes]

    for ax, (arm, (wmin, wmax)) in zip(axes, arm_ranges.items()):
        sel = (wave >= wmin) & (wave < wmax)
        w   = wave[sel]
        sky = sky_flux[sel]
        mod = palace_model_scaled[sel]
        m   = mask[sel]

        # blank contaminated pixels for the continuum trace
        sky_clean = np.where(m, sky, np.nan)

        # log-scale y limits: positive floor from p1 of sky, top at p99.5
        sky_pos = sky[np.isfinite(sky) & (sky > 0)]
        ymin = float(np.nanpercentile(sky_pos, 1)) * 0.3 if len(sky_pos) else 0.01
        ymax = float(np.nanpercentile(sky, 99.5)) * 2.0

        # green shading for clean regions (drawn first, behind spectra)
        ax.fill_between(w, ymin, ymax, where=m,
                        alpha=0.15, color="limegreen", step="mid",
                        label="Clean (mask=True)")

        ax.step(w, sky,       color="black",         lw=0.5, where="mid",
                alpha=0.30,   label="Sky (full)")
        ax.step(w, sky_clean, color="black",         lw=2.0, where="mid",
                label="Sky (clean only)")
        ax.step(w, mod,       color=arm_colors[arm], lw=0.9, where="mid",
                label="PALACE model (scaled)")
        ax.axhline(threshold, color="red", lw=1.1, ls="--",
                   label=f"Threshold = {threshold:.3g}")

        ax.set_yscale("log")
        ax.set_xlim(wmin, wmax)
        ax.set_ylim(ymin, ymax)
        ax.set_ylabel("Flux × FACTOR", fontsize=9)
        ax.set_title(f"{arm} arm  ({wmin:.0f}–{wmax:.0f} Å)", fontsize=10)
        ax.legend(fontsize=8, loc="upper right", ncol=2)

    axes[-1].set_xlabel("Wavelength (Å)", fontsize=10)
    if title:
        fig.suptitle(title, fontsize=10, y=1.002)
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath, bbox_inches="tight")
        print(f"Plot saved to {savepath}")
    if show:
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Build a sky-line-free mask from PALACE models for LVM continuum fitting.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("fits_file",    help="LVM XCframe FITS file")
    p.add_argument("palace_dir",   help="Path to palace/PMD directory")
    p.add_argument("--sky-ext",    default="SKY_EAST",
                   help="FITS extension name for the sky spectrum")
    p.add_argument("--threshold",  type=float, default=THRESHOLD_DEFAULT,
                   help="Sky line contamination threshold in FACTOR*flux units")
    p.add_argument("--factor",     type=float, default=FACTOR_DEFAULT,
                   help="Flux scale factor applied when reading spectra")
    p.add_argument("--no-lsf",     action="store_true",
                   help="Use fixed --lsf-sigma instead of reading LSF from file")
    p.add_argument("--lsf-sigma",  type=float, default=LSF_SIGMA_DEFAULT,
                   help="Fixed Gaussian LSF sigma in Å (used with --no-lsf)")
    p.add_argument("--output",     default=None,
                   help="Output FITS path (default: palace_mask_<stem>.fits)")
    p.add_argument("--min-window", type=float, default=5.0,
                   help="Minimum window width in Å for the printed window table")
    p.add_argument("--plot",       action="store_true",
                   help="Display the diagnostic plot interactively (always saved as PDF)")
    args = p.parse_args()

    # default output filename
    outpath = args.output or f"palace_mask_{Path(args.fits_file).stem}.fits"

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    print(f"Reading {args.fits_file} ...")
    with fits.open(args.fits_file) as hdul:
        wave    = hdul["WAVE"].data.astype(np.float64)
        sky_all = hdul[args.sky_ext].data.astype(np.float64) * args.factor
        lsf_fwhm_all = (hdul["LSF"].data.astype(np.float64)
                        if (not args.no_lsf and "LSF" in hdul) else None)

    sky_median = np.nanmedian(sky_all, axis=0)

    if lsf_fwhm_all is not None:
        lsf_sigma = np.nanmedian(lsf_fwhm_all, axis=0) / 2.355
        print(f"LSF from file (σ range: {lsf_sigma.min():.2f}–{lsf_sigma.max():.2f} Å)")
    else:
        lsf_sigma = args.lsf_sigma
        print(f"Fixed LSF σ = {lsf_sigma:.2f} Å")

    print(f"Wave: {wave.min():.1f}–{wave.max():.1f} Å  ({len(wave)} pixels)")
    print(f"Sky: {args.sky_ext}  median over {sky_all.shape[0]} fibers")

    # ------------------------------------------------------------------
    # Build PALACE model
    # ------------------------------------------------------------------
    print(f"\nBuilding PALACE sky line model from {args.palace_dir} ...")
    palace_model, _ = build_palace_line_model(
        wave, args.palace_dir, lsf_sigma=lsf_sigma, verbose=True
    )

    # ------------------------------------------------------------------
    # Scale to observed sky and build mask
    # ------------------------------------------------------------------
    scale = scale_model_to_sky(wave, palace_model, sky_median)
    palace_model_scaled = palace_model * scale
    print(f"Scale factor (PALACE → sky units): {scale:.4g}")

    mask = make_sky_mask(palace_model_scaled, args.threshold)
    print(f"\nThreshold: {args.threshold:.4g}  "
          f"(= {args.threshold / args.factor:.2e} erg s⁻¹ cm⁻² Å⁻¹)")

    # ------------------------------------------------------------------
    # Report and save
    # ------------------------------------------------------------------
    report_coverage(wave, mask)
    report_windows(wave, mask, min_width=args.min_window)

    save_output(outpath, wave, sky_median, mask, palace_model_scaled,
                args.threshold, scale, args.factor, args.sky_ext)

    plot_path = Path(outpath).with_suffix(".pdf")
    title = (f"{Path(args.fits_file).name}  |  "
             f"threshold = {args.threshold:.3g}  |  "
             f"scale = {scale:.3g}  |  "
             f"FACTOR = {args.factor:.1e}")
    make_plot(wave, sky_median, palace_model_scaled, mask, args.threshold,
              title=title, savepath=plot_path, show=args.plot)


if __name__ == "__main__":
    main()
