#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Build a sky-line contamination mask for LVM spectra using PALACE sky
    emission models, identifying wavelength windows clean enough for
    continuum measurement.

Command line usage (if any):

    usage: palace_make_mask.py [-h] [--sky-ext NAME] [--threshold T]
                               [--factor F] [--no-lsf] [--lsf-sigma S]
                               [--output PATH] [--min-window W] [--plot]
                               fits_file [palace_dir]

    where

    fits_file     is an LVM XCframe FITS file providing WAVE, a sky spectrum
                  extension (default SKY_EAST), and optionally LSF.

    palace_dir    is the path to the palace/PMD directory containing the
                  PALACE data files (pmd_popmodel_OH.dat, etc.).  Optional;
                  defaults to the vendored copy at
                  data/palace_ref/palace/PMD (relative to this script's own
                  location), so it only needs to be given to point at a
                  different PMD installation.

    --sky-ext NAME    sky spectrum extension name (default: SKY_EAST).

    --threshold T     contamination threshold in FACTOR*flux units; pixels
                      where the model exceeds this are masked.  Lower values
                      give a stricter mask with fewer clean pixels.
                      (default: 0.01)

    --factor F        flux scale factor applied on read; brings raw cgs flux
                      (~7e-15 erg/s/cm2/A) to order unity (~0.7).
                      (default: 1e14)

    --no-lsf          use fixed --lsf-sigma instead of reading LSF from file.

    --lsf-sigma S     fixed Gaussian LSF sigma in Å (default: 0.65).

    --output PATH     output FITS path (default: <stem>_mask.fits).

    --line-output PATH  output strong-sky-line-list path
                      (default: <stem>_lines.txt).

    --min-window W    minimum clean window width in Å for the printed table
                      (default: 5).

    --plot            display the diagnostic plot interactively (it is always
                      saved as a PNG alongside the output FITS).

Description:

    Night-sky emission across the LVM wavelength range (3600-9800 A) is
    dominated by OH vibrational-rotational lines in the Z arm (7400-9800 A),
    atomic forbidden lines ([NI], OI, NaI, KI) in the B and R arms, and the
    O2 A-band near 8650 A.  These features contaminate attempts to fit a
    smooth stellar or nebular continuum through the spectrum.

    This program constructs a model of the sky line emission from the PALACE
    (Paranal Airglow Line And Continuum Emission, Noll et al. 2024) reference
    data set, renders each component onto the LVM wavelength grid convolved
    with the instrument LSF, scales the result to match an observed sky
    spectrum, and produces a boolean mask that is True at pixels where the
    predicted contamination is below the threshold.

    Four PALACE components are used:
      OH    hydroxyl vibrational-rotational bands (pmd_popmodel_OH.dat);
            amplitudes proportional to Einstein Aij x degeneracy gi, grouped
            by (v_upper, N_upper, F_upper) quantum numbers.
      OI    oxygen recombination multiplets at 7774 and 8446 A
            (pmd_intmodel_Orc.dat).
      atom  atomic lines NaI, KI, [NI], OI green/red from the PALACE
            climatology in Rayleigh units (pmd_intdata_atom.dat).
      O2    molecular oxygen A-band near 8650 A at T_rot = 191.5 K
            (pmd_popmodel_O2.dat).

    Each component is normalised to its own peak before summing so that no
    single family dominates the mask.  The combined model is then scaled to
    the observed sky spectrum via a least-squares fit to the bright OH pixels
    in the Z arm.

    The threshold controls the tradeoff between mask strictness and the
    fraction of pixels available for continuum fitting.  The default of 0.01
    (FACTOR=1e14 units) corresponds to 1e-16 erg/s/cm2/A, well below the
    typical Z-arm continuum (~7e-15 erg/s/cm2/A).

    Output FITS extensions:
      WAVE       float32 (Npix,)   wavelength array (A, air)
      SKY        float32 (Npix,)   median sky spectrum (FACTOR-scaled)
      CONTINUUM  float32 (Npix,)   sky with masked pixels set to NaN
      MASK       uint8   (Npix,)   1 = clean, 0 = contaminated

    A three-panel diagnostic PNG is always written alongside the FITS output.

    In addition to the mask, a strong-sky-line list is written to an ascii
    table (--line-output, default <stem>_lines.txt).  This is a *labeled
    line position* view of the same contamination model, rather than a
    per-pixel mask: for each PALACE line group (OH by (v_upper, N_upper,
    F_upper), OI recombination by reffeat, atomic lines by feat, O2 as one
    group) the brightest transition's wavelength is kept if the combined,
    scaled contamination model exceeds --threshold there -- the same
    threshold that builds the mask.  Columns: Wave_air, LineID, Component,
    Ampl.  The column names match what PlotSpecI.py's -lines/load_lines()
    mechanism expects, so this file can be used directly as a second,
    sky-line overlay via PlotSpecI.py's -sky_lines flag, alongside the
    usual scientific line list (data/dap_lines.txt).

Primary routines:

    build_palace_line_model
    scale_model_to_sky
    make_sky_mask
    find_clean_windows
    find_strong_sky_lines

Notes:

    The model scale factor is derived from the Z arm only (OH-dominated),
    so the absolute threshold values are most meaningful in the Z arm.  The
    mask is still useful in the B and R arms but atomic line heights may be
    approximate relative to OH.

History::

    260627  ksl  Coding begun
    260628  ksl  Renamed to palace_make_mask.py; added PNG output; threshold default 0.01
    260810  ksl  palace_dir made optional (nargs='?'), defaulting to the vendored
                  copy at data/palace_ref/palace/PMD (DEFAULT_PALACE_DIR, computed
                  relative to this file's own location, same pattern as
                  XSkySepIvan.py's DEFAULT_BASE_DIR); still overridable by passing
                  an explicit path. report_coverage() gained optional threshold/
                  factor args to print a per-arm "Flux (cgs)" column (the
                  threshold expressed as threshold/factor, erg/s/cm2/A) alongside
                  the existing per-arm % clean.
    260815  ksl  Added find_strong_sky_lines() and --line-output: alongside the
                  pixel mask, write a named strong-sky-line list (one row per
                  PALACE line group -- OH by (v_upper,N_upper,F_upper), OI
                  recombination by reffeat, atomic lines by feat, O2 as one
                  group -- at that group's brightest transition, kept if the
                  combined scaled model exceeds --threshold there). Output
                  columns Wave_air/LineID/Component/Ampl match what
                  PlotSpecI.py's -sky_lines overlay expects. The four
                  _build_* component loaders were refactored to share
                  _load_oh_table/_load_orc_table/_load_atom_table/
                  _load_o2_table so the new finder reuses the same catalogue
                  parsing rather than duplicating it. Output written as
                  ascii.fixed_width_two_line (round-trips cleanly through
                  astropy.io.ascii.read()'s format guesser, unlike plain
                  ascii.fixed_width).
'''

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table


FACTOR_DEFAULT   = 1.0e14
LSF_SIGMA_DEFAULT = 0.65       # Å  (≈ FWHM 1.5 Å / 2.355)
THRESHOLD_DEFAULT = 0.01
CAP = 5.0                      # Å padding when loading line catalogues

# Vendored PMD location (see py_progs/XSkySepIvan.py's DEFAULT_BASE_DIR,
# which points at data/palace_ref -- this is that same directory's palace/PMD
# subdirectory, computed the same way, relative to this file's location).
DEFAULT_PALACE_DIR = (Path(__file__).resolve().parent.parent
                       / 'data' / 'palace_ref' / 'palace' / 'PMD')

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

def _load_oh_table(wave, palace_dir, cap):
    """Load, wavelength-window, and quantum-number-decode pmd_popmodel_OH.dat."""
    oh = Table.read(str(palace_dir / "pmd_popmodel_OH.dat"),
                    format="ascii.basic", guess=False, comment="#", fast_reader=False)
    oh["wave"] = vac_to_air(np.asarray(oh["lam"], float) * 1e4)
    oh = oh[(oh["wave"] >= wave.min() - cap) & (oh["wave"] <= wave.max() + cap)]
    return _decode_hitran_id(oh)


def _build_oh(wave, palace_dir, lsf_sigma, cap):
    """Render OH vibrational-rotational line groups from pmd_popmodel_OH.dat."""
    oh = _load_oh_table(wave, palace_dir, cap)
    model = np.zeros(wave.size)
    for grp in oh.group_by(("v_upper", "N_upper", "F_upper")).groups:
        amp = np.asarray(grp["Aij"], float) * np.asarray(grp["gi"], float)
        model += _grp2vector(grp["wave"], amp, wave, lsf_sigma)
    return model


def _load_orc_table(wave, palace_dir, cap):
    """Load and wavelength-window pmd_intmodel_Orc.dat."""
    orc = Table.read(str(palace_dir / "pmd_intmodel_Orc.dat"),
                     format="ascii.basic", guess=False, comment="#", fast_reader=False)
    orc["wave"] = vac_to_air(np.asarray(orc["lam"], float) * 1e4)
    return orc[(orc["wave"] >= wave.min() - cap) & (orc["wave"] <= wave.max() + cap)]


def _build_orc(wave, palace_dir, lsf_sigma, cap):
    """Render OI recombination line multiplets from pmd_intmodel_Orc.dat."""
    orc = _load_orc_table(wave, palace_dir, cap)
    model = np.zeros(wave.size)
    for grp in orc.group_by("reffeat").groups:
        model += _grp2vector(grp["wave"], np.asarray(grp["I"], float), wave, lsf_sigma)
    return model


def _load_atom_table(wave, palace_dir, cap):
    """Load and wavelength-window pmd_intdata_atom.dat.

    Hydrogen recombination lines and OI recombination lines (class 'H' and
    'Orc') are excluded; they are either negligible or handled separately.
    Intensities are PALACE reference values in Rayleigh units.
    """
    atom = Table.read(str(palace_dir / "pmd_intdata_atom.dat"),
                      format="ascii.basic", guess=False, comment="#", fast_reader=False)
    atom["wave"] = vac_to_air(np.asarray(atom["lam"], float) * 1e4)
    atom = atom[(atom["wave"] >= wave.min() - cap) & (atom["wave"] <= wave.max() + cap)]
    return atom[~np.isin(np.asarray(atom["class"], str), ["H", "Orc"])]


def _build_atom(wave, palace_dir, lsf_sigma, cap):
    """Render atomic sky lines (NaI, KI, [NI], OI) from pmd_intdata_atom.dat."""
    atom = _load_atom_table(wave, palace_dir, cap)
    model = np.zeros(wave.size)
    for grp in atom.group_by("class").groups:
        model += _grp2vector(grp["wave"], np.asarray(grp["I"], float), wave, lsf_sigma)
    return model


def _load_o2_table(palace_dir, t_o2=191.5):
    """Load pmd_popmodel_O2.dat and compute relative line weights.

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
    o2["rel"] = rel / rel.sum()
    return o2


def _build_o2(wave, palace_dir, lsf_sigma, t_o2=191.5):
    """Render the O2 A-band near 8650 Å from pmd_popmodel_O2.dat."""
    o2 = _load_o2_table(palace_dir, t_o2=t_o2)
    return _grp2vector(np.asarray(o2["wave"], float), np.asarray(o2["rel"], float),
                       wave, lsf_sigma)


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
# Strong sky-line list (for overlay on spectrum plots, not masking)
# ---------------------------------------------------------------------------

def find_strong_sky_lines(wave, palace_model_scaled, palace_dir, threshold, cap=CAP):
    """Build a named list of strong sky lines for overlay on spectrum plots.

    For each PALACE line group already used to build the contamination
    model (OH grouped by (v_upper, N_upper, F_upper), OI recombination
    grouped by reffeat, atomic lines grouped by feat -- the named
    line/multiplet, e.g. NaI0589, OI0558, KI0770 -- and the O2 A-band
    treated as one group), this takes that group's single brightest
    transition wavelength and evaluates the combined, scaled contamination
    model (palace_model_scaled -- the same array used to build the mask)
    there.  A group is kept if that value exceeds `threshold`, the same
    threshold used to build the mask, so the mask (pixel coverage) and this
    list (named line positions) describe the same underlying contamination
    level from two different angles.

    Parameters
    ----------
    wave : ndarray
        Wavelength array (Å) matching palace_model_scaled.
    palace_model_scaled : ndarray
        Combined PALACE model scaled to sky flux units (as used for the
        mask).
    palace_dir : str or Path
        Path to the palace/PMD directory.
    threshold : float
        Contamination threshold in the same units as palace_model_scaled.
    cap : float
        Wavelength padding in Å beyond the wave grid edges when loading
        line catalogues.

    Returns
    -------
    astropy.table.Table
        Columns Wave_air, LineID, Component, Ampl, sorted by Wave_air.
        Empty (but correctly typed) if no group exceeds threshold.
    """
    palace_dir = Path(palace_dir)
    rows = []

    oh = _load_oh_table(wave, palace_dir, cap)
    for grp in oh.group_by(("v_upper", "N_upper", "F_upper")).groups:
        amp = np.asarray(grp["Aij"], float) * np.asarray(grp["gi"], float)
        i = int(np.argmax(amp))
        label = "OH_%d-%d" % (int(grp["v_upper"][i]), int(grp["v_lower"][i]))
        rows.append((float(grp["wave"][i]), label, "OH"))

    orc = _load_orc_table(wave, palace_dir, cap)
    for grp in orc.group_by("reffeat").groups:
        i = int(np.argmax(np.asarray(grp["I"], float)))
        rows.append((float(grp["wave"][i]), str(grp["reffeat"][i]), "OI"))

    atom = _load_atom_table(wave, palace_dir, cap)
    for grp in atom.group_by("feat").groups:
        i = int(np.argmax(np.asarray(grp["I"], float)))
        rows.append((float(grp["wave"][i]), str(grp["feat"][i]), "atom"))

    o2 = _load_o2_table(palace_dir)
    o2 = o2[(o2["wave"] >= wave.min() - cap) & (o2["wave"] <= wave.max() + cap)]
    if len(o2):
        i = int(np.argmax(np.asarray(o2["rel"], float)))
        rows.append((float(o2["wave"][i]), "O2", "O2"))

    if not rows:
        out = Table(names=("Wave_air", "LineID", "Component", "Ampl"),
                    dtype=(float, "U16", "U8", float))
        return out

    line_wave = np.array([r[0] for r in rows])
    amp = np.interp(line_wave, wave, palace_model_scaled)
    keep = amp > threshold

    out = Table(rows=[r for r, k in zip(rows, keep) if k],
               names=("Wave_air", "LineID", "Component"))
    out["Wave_air"] = np.round(out["Wave_air"], 2)
    out["Ampl"] = np.round(amp[keep], 4)
    out.sort("Wave_air")
    return out


# ---------------------------------------------------------------------------
# Console reporting
# ---------------------------------------------------------------------------

def report_coverage(wave, mask, arm_ranges=ARM_RANGES, threshold=None, factor=None):
    """Print a per-arm summary of clean pixel counts and window counts.

    If threshold and factor are given, also prints the threshold expressed
    as a physical flux value (threshold / factor, erg s⁻¹ cm⁻² Å⁻¹) in its
    own column -- the same absolute value in every row, since one global
    threshold applies across the whole spectrum, but shown per arm so the
    cutoff is visible alongside each arm's clean fraction without having to
    look back at the threshold line printed above the table.
    """
    show_flux = threshold is not None and factor is not None
    thresh_flux = (threshold / factor) if show_flux else None
    flux_hdr = f"  {'Flux (cgs)':>11}" if show_flux else ""
    print(f"\n{'Arm':<5}  {'Range (Å)':<16}  {'N total':>8}  "
          f"{'N clean':>8}  {'% clean':>8}  {'N windows':>10}{flux_hdr}")
    print("-" * (65 + (13 if show_flux else 0)))
    for arm, (wmin, wmax) in arm_ranges.items():
        sel     = (wave >= wmin) & (wave < wmax)
        n_tot   = int(sel.sum())
        n_clean = int((mask & sel).sum())
        pct     = 100.0 * n_clean / n_tot if n_tot else 0.0
        wins    = find_clean_windows(wave[sel], mask[sel])
        flux_col = f"  {thresh_flux:>11.2e}" if show_flux else ""
        print(f"{arm:<5}  {wmin:.0f}–{wmax:.0f} Å       "
              f"{n_tot:>8d}  {n_clean:>8d}  {pct:>7.1f}%  {len(wins):>10d}{flux_col}")
    n_tot   = len(wave)
    n_clean = int(mask.sum())
    wins    = find_clean_windows(wave, mask)
    flux_col = f"  {thresh_flux:>11.2e}" if show_flux else ""
    print(f"{'ALL':<5}  {wave.min():.0f}–{wave.max():.0f} Å  "
          f"{n_tot:>8d}  {n_clean:>8d}  {100*n_clean/n_tot:>7.1f}%  {len(wins):>10d}{flux_col}")


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

    Extensions: PRIMARY (no data, header only), WAVE (float32 wavelengths),
    SKY (float32 median sky spectrum, FACTOR-scaled), CONTINUUM (float32 sky
    with contaminated pixels set to NaN), MASK (uint8 clean-pixel mask).

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
        fig.savefig(savepath, bbox_inches="tight", dpi=150)
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
    p.add_argument("palace_dir",   nargs="?", default=str(DEFAULT_PALACE_DIR),
                   help="Path to palace/PMD directory")
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
                   help="Output FITS path (default: <stem>_mask.fits)")
    p.add_argument("--line-output", default=None,
                   help="Output strong-sky-line-list path "
                        "(default: <stem>_lines.txt)")
    p.add_argument("--min-window", type=float, default=5.0,
                   help="Minimum window width in Å for the printed window table")
    p.add_argument("--plot",       action="store_true",
                   help="Display the diagnostic plot interactively (always saved as PDF)")
    args = p.parse_args()

    # default output filename
    outpath = args.output or f"{Path(args.fits_file).stem}_mask.fits"

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
    report_coverage(wave, mask, threshold=args.threshold, factor=args.factor)
    report_windows(wave, mask, min_width=args.min_window)

    save_output(outpath, wave, sky_median, mask, palace_model_scaled,
                args.threshold, scale, args.factor, args.sky_ext)

    # ------------------------------------------------------------------
    # Strong sky-line list (for overlay on spectrum plots, e.g. PlotSpecI.py)
    # ------------------------------------------------------------------
    line_table = find_strong_sky_lines(wave, palace_model_scaled, args.palace_dir,
                                       args.threshold)
    line_outpath = args.line_output or f"{Path(args.fits_file).stem}_lines.txt"
    line_table.write(line_outpath, format="ascii.fixed_width_two_line", overwrite=True)
    print(f"\nStrong sky lines: {len(line_table)} found "
          f"(threshold {args.threshold:.4g})  ->  {line_outpath}")

    plot_path = Path(outpath).with_suffix(".png")
    title = (f"{Path(args.fits_file).name}  |  "
             f"threshold = {args.threshold:.3g}  |  "
             f"scale = {scale:.3g}  |  "
             f"FACTOR = {args.factor:.1e}")
    make_plot(wave, sky_median, palace_model_scaled, mask, args.threshold,
              title=title, savepath=plot_path, show=args.plot)


if __name__ == "__main__":
    main()
