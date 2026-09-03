#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Predict the sky spectrum at the Sci telescope's pointing for a single
    exposure, using the trained mlp_ensemble_split_zodi MLP ensemble and the
    SkyE/SkyW spectra recorded in an LVM XCframe summary FITS file.

Command line usage (if any):

    usage: PredictSky.py [-h] [-row N | -expnum N] -model PATH
                         [-output PATH] [-lvmsky_skysub PATH]
                         [-lvmcore_dir PATH]
                         fits_file

    where

    fits_file       is the path to an LVM XCframe summary FITS file
                    containing WAVE, FLUX, SKY_EAST, SKY_WEST, LSF
                    extensions and a DRP_ALL table (as produced by
                    SummarizeCframe.py in -by fiber mode).

    -row N          selects the exposure by row index into the file
                    (default: 0).

    -expnum N       selects the exposure by its DRP_ALL 'expnum' value
                    instead of a row index (overrides -row).

    -model PATH     path to the trained ensemble .pt archive. Required,
                    no default -- this script is meant to work against
                    different trained models for different purposes, so
                    the checkpoint is always named explicitly rather than
                    silently falling back to whichever one was current
                    when the script was last edited.

    -output PATH    output FITS path (default: PredictedSky_<expnum>.fits).

    -lvmsky_skysub PATH
                    path to the lvmsky repo's skysub/ directory, which
                    supplies the mlp_predictor and sky_decomp packages this
                    script imports (default: ~/SDSS/lvmsky/skysub).

    -lvmcore_dir PATH
                    sets LVMCORE_DIR, needed by mlp_predictor.data for the
                    LCO extinction curve (default: ~/SDSS/lvmcore).

Description:

    For the selected exposure, this program:

      1. Physically decomposes the observed SKY_EAST and SKY_WEST spectra
         into the 433-dim coefficient vector the trained ensemble consumes,
         via sky_decomp.lsf_surface_iterative.SkyDecompLSFSurfaceIterative
         (a deterministic QP fit -- no training/model-building involved).
      2. Runs those coefficients, plus the exposure's MJD and the three
         telescopes' RA/Dec, through the frozen, already-trained 10-seed
         ensemble (mlp_predictor.inference.predict_sky_from_minimal_inputs)
         to get the predicted coefficient vector at the Sci pointing, with
         a per-coefficient ensemble-spread uncertainty and a scalar
         confidence score.
      3. Reconstructs a predicted flux spectrum from those coefficients,
         using the exposure's own LSF (which SummarizeCframe.py -by fiber
         draws from the same science-fiber window as FLUX, so it is already
         the Sci-arm LSF -- no separate LSF lookup needed).
      4. Writes the observed and predicted spectra, the raw coefficient
         vectors, and diagnostics to an output FITS file.

    This is an inference-only pipeline: the ensemble's weights are loaded
    and used exactly as trained. Step 1 is unavoidable physical feature
    extraction (the network cannot consume raw flux), not model-building.

Primary routines:

    read_row
    decompose_arm
    predict_row
    write_output

Notes:

    - The input file has no IVAR extension, so both arm decompositions use
      a uniform ivar = 1 fallback (same precedent as decompose_parallel.py)
      rather than true photon-noise weighting.
    - The single LSF extension is the Sci-fiber window's own LSF (per
      SummarizeCframe.py), reused here for the SkyE/SkyW decomposition
      fits' *initial* LSF guess and, more importantly, for the final Sci
      reconstruction step -- an approximation that assumes the three LVM
      telescopes have near-identical internal optical performance.
    - This script imports mlp_predictor/sky_decomp from an external lvmsky
      checkout (-lvmsky_skysub, default ~/SDSS/lvmsky/skysub) via
      sys.path; it only works when that checkout's tree contains those
      packages (confirmed merged into lvmsky's main branch as of this
      writing) and the palace/PMD tables matching PALACE_OH_SUFFIX /
      PALACE_DIFFUSE_SUFFIX below.
    - PALACE_OH_SUFFIX / PALACE_DIFFUSE_SUFFIX are left as None so
      SkyDecomp resolves its own bundled defaults (lvmsky main commit
      4c2a1e4 made this the correct choice -- the ensemble's coef_names
      imply 402 OH groups, which the bundled "_h_family_default_ef_v1"
      table produces; an earlier hardcoded "_joint_v2_updated" override
      here only had 357 groups and made predict_sky_from_minimal_inputs
      raise a coef-count mismatch).

History::

    260901  ksl  Coding begun.
    260903  ksl  Switched every option from double-dash (--lvmsky-skysub)
        to single-dash (-lvmsky_skysub), matching py_progs/'s convention
        -- see BatchPredictSkyESO.py's History for the fuller note; this
        also required fixing BatchPredictSky.py's worker-init sys.argv
        injection, which fakes a command line for this module's own
        module-level _pre parser and had hardcoded the old spelling.

'''

import argparse
import os
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

# ---------------------------------------------------------------------------
# External package setup -- mlp_predictor / sky_decomp live in the lvmsky
# repo, not in this project.  Inserted before argparse runs so -lvmsky_skysub
# / -lvmcore_dir can override the defaults before the imports below fire.
# ---------------------------------------------------------------------------

DEFAULT_LVMSKY_SKYSUB = Path("~/SDSS/lvmsky/skysub").expanduser()
DEFAULT_LVMCORE_DIR = Path("~/SDSS/lvmcore").expanduser()

_pre = argparse.ArgumentParser(add_help=False)
_pre.add_argument("-lvmsky_skysub", default=str(DEFAULT_LVMSKY_SKYSUB))
_pre.add_argument("-lvmcore_dir", default=str(DEFAULT_LVMCORE_DIR))
_pre_args, _ = _pre.parse_known_args()

sys.path.insert(0, _pre_args.lvmsky_skysub)
os.environ.setdefault("LVMCORE_DIR", _pre_args.lvmcore_dir)

from mlp_predictor import serialization, inference, wavelengths as _wave_mod  # noqa: E402
from scipy.stats import norm  # noqa: E402
from sky_decomp.lsf_surface_iterative import (  # noqa: E402
    SkyDecompLSFSurfaceIterative, LSFSurfaceIterativeConfig,
    LSF_TAP_OFFSETS, build_lsf_operator,
)

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

FACTOR = 1e14  # XCframe flux is erg/s/cm^2/A; the QP fit works in O(1) counts.
# None -> SkyDecomp resolves its own bundled PALACE OH/diffuse defaults
# (as of lvmsky main 4c2a1e4, this is the "_h_family_default_ef_v1" OH
# table, which matches the 402-group coef_names the checkpoint expects;
# the previously hardcoded "_joint_v2_updated" table only has 357 groups
# and was never what this ensemble was trained against).
PALACE_OH_SUFFIX = None
PALACE_DIFFUSE_SUFFIX = None
COMPONENT_KEYS = ("oh", "moon", "zodi", "diffuse", "atom", "orc", "o2")
# Discrete/emission-line components vs. smooth/broadband ones, per
# SkyDecomp's own physical basis (see ivan.md): OH Meinel bands, atomic
# lines (Na/K/[NI]/[OI]5577), ORC (OI recombination 7774/8446), and the O2
# A-band are all line-like; moon-scattered light, zodiacal light, and the
# diffuse airglow continua (HO2/FeO/O2Ac) are smooth continuum. Continuum
# is not summed separately here -- get it via FLUX_PRED - LINE_PRED.
LINE_KEYS = ("oh", "atom", "orc", "o2")


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

def read_row(fits_file, row=None, expnum=None):
    """Read one exposure's spectra and metadata from an XCframe summary file.

    Parameters
    ----------
    fits_file : str or Path
        Path to the XCframe summary FITS file (WAVE, FLUX, SKY_EAST,
        SKY_WEST, LSF extensions; DRP_ALL table).
    row : int, optional
        Row index to select (default 0 if expnum is also None).
    expnum : int, optional
        DRP_ALL 'expnum' value to select instead of a row index.

    Returns
    -------
    dict
        Keys: row, wave, flux_sci, flux_e, flux_w, lsf_fwhm, mjd, sci_ra,
        sci_dec, skye_ra, skye_dec, skyw_ra, skyw_dec, expnum, drp_row.
    """
    with fits.open(fits_file, memmap=True) as hdul:
        drp = Table(hdul["DRP_ALL"].data)

        if expnum is not None:
            matches = np.flatnonzero(np.asarray(drp["expnum"]) == int(expnum))
            if matches.size == 0:
                raise ValueError(f"expnum {expnum} not found in {fits_file}")
            i = int(matches[0])
        else:
            i = int(row) if row is not None else 0

        wave = np.asarray(hdul["WAVE"].data, dtype=np.float64)
        flux_sci = np.asarray(hdul["FLUX"].data[i], dtype=np.float64)
        flux_e = np.asarray(hdul["SKY_EAST"].data[i], dtype=np.float64)
        flux_w = np.asarray(hdul["SKY_WEST"].data[i], dtype=np.float64)
        lsf_fwhm = np.asarray(hdul["LSF"].data[i], dtype=np.float64)
        # SummarizeCframe.py's sigma-clipped fiber-window average leaves the
        # first/last pixel NaN (no interior neighbour on one side); fill by
        # interpolation so the LSF-based reconstruction kernel downstream
        # (build_lsf_operator) never sees a non-finite value.
        bad = ~np.isfinite(lsf_fwhm)
        if bad.any():
            good = ~bad
            lsf_fwhm[bad] = np.interp(
                np.flatnonzero(bad), np.flatnonzero(good), lsf_fwhm[good])
        drp_row = drp[i]

    mjd = Time(str(drp_row["obstime"]), format="isot", scale="utc").mjd

    return dict(
        row=i,
        wave=wave,
        flux_sci=flux_sci, flux_e=flux_e, flux_w=flux_w, lsf_fwhm=lsf_fwhm,
        mjd=mjd,
        sci_ra=float(drp_row["sci_ra"]), sci_dec=float(drp_row["sci_dec"]),
        skye_ra=float(drp_row["skye_ra"]), skye_dec=float(drp_row["skye_dec"]),
        skyw_ra=float(drp_row["skyw_ra"]), skyw_dec=float(drp_row["skyw_dec"]),
        expnum=int(drp_row["expnum"]),
        drp_row=drp_row,
    )


# ---------------------------------------------------------------------------
# Decomposition
# ---------------------------------------------------------------------------

def decompose_arm(decomp, flux_row, label):
    """Scale and QP-decompose one sky-arm spectrum; print a one-line summary.

    Parameters
    ----------
    decomp : SkyDecompLSFSurfaceIterative
        Decomposition instance (reused across both arms).
    flux_row : ndarray
        Raw observed flux for this arm, native XCframe units.
    label : str
        Arm label used only for the printed summary (e.g. 'SkyE').

    Returns
    -------
    SkyDecompResult
        The fit result; `.coef` is the 433-dim coefficient vector.
    """
    flux = flux_row * FACTOR
    ivar = np.ones_like(flux)  # no IVAR extension in the source file
    fit = decomp.fit(flux, ivar, verbose=False)
    print(f"  {label} decomp: status={fit.fit_status}, "
          f"chi2_red={fit.reduced_chi2:.2f}, R^2={fit.r2:.4f}")
    return fit


# ---------------------------------------------------------------------------
# Prediction
# ---------------------------------------------------------------------------

def build_decomp(wave, ensemble):
    """Construct the SkyDecompLSFSurfaceIterative instance predict_row needs.

    Split out of predict_row so a batch loop over many exposures (same
    XCframe file, same WAVE grid) can build this once -- loading the
    PALACE tables is a real cost, and decompose_arm() already proves one
    instance can be fit() repeatedly (SkyE then SkyW) safely, so reuse
    across rows is equally safe.
    """
    n_moon_knots, split_zodi, n_zodi_knots = _wave_mod.infer_spline_knots(
        ensemble["coef_names"])
    base_dir = Path(_pre_args.lvmsky_skysub) / "sky_decomp" / "data"
    return SkyDecompLSFSurfaceIterative(
        wave, lsf_sigma=1.0, n_spline_knots=n_moon_knots,
        base_dir=base_dir,
        palace_oh_suffix=PALACE_OH_SUFFIX,
        palace_diffuse_suffix=PALACE_DIFFUSE_SUFFIX,
        split_zodi=split_zodi, n_zodi_spline_knots=n_zodi_knots,
        config=LSFSurfaceIterativeConfig(n_refinement_cycles=5),
    )


def predict_row(row_data, model_path=None, ensemble=None, decomp=None):
    """Decompose SkyE/SkyW and predict the Sci-pointing sky for one exposure.

    Parameters
    ----------
    row_data : dict
        Output of `read_row`.
    model_path : str or Path, optional
        Path to the trained ensemble .pt archive. Ignored if `ensemble` is
        given; required otherwise (no default checkpoint -- see module
        Synopsis).
    ensemble : dict, optional
        An already-loaded ensemble (`serialization.load_ensemble`'s
        return value) -- pass this in a batch loop over many exposures so
        the (56 MB) checkpoint is loaded once, not once per exposure.
    decomp : SkyDecompLSFSurfaceIterative, optional
        An already-built decomposer (see build_decomp) -- pass this in a
        batch loop over many exposures from the same XCframe file (same
        WAVE grid) so the PALACE tables load once, not once per exposure.
        Its LSF state is overwritten per-call, so reuse across rows is
        safe as long as they share the same wave grid.

    Returns
    -------
    dict
        Keys: wave, flux_pred, flux_pred_lo, flux_pred_hi, coef, coef_std,
        confidence, coef_names, fit_east, fit_west, ensemble_config.
    """
    if ensemble is None:
        if model_path is None:
            raise ValueError("predict_row: either model_path or ensemble must be given")
        ensemble = serialization.load_ensemble(str(model_path))
    if decomp is None:
        decomp = build_decomp(row_data["wave"], ensemble)

    fit_east = decompose_arm(decomp, row_data["flux_e"], "SkyE")
    fit_west = decompose_arm(decomp, row_data["flux_w"], "SkyW")

    result = inference.predict_sky_from_minimal_inputs(
        ensemble,
        obstime_mjd=[row_data["mjd"]],
        sci_ra=[row_data["sci_ra"]], sci_dec=[row_data["sci_dec"]],
        sky_e_ra=[row_data["skye_ra"]], sky_e_dec=[row_data["skye_dec"]],
        sky_w_ra=[row_data["skyw_ra"]], sky_w_dec=[row_data["skyw_dec"]],
        coef_e=fit_east.coef.astype(np.float32)[None, :],
        coef_w=fit_west.coef.astype(np.float32)[None, :],
    )
    coef = result["coef"][0].astype(np.float64)
    coef_std = result["coef_std"][0].astype(np.float64)
    confidence = float(result["confidence"][0])

    # Install this exposure's own (Sci-fiber-window) LSF for reconstruction,
    # replacing whatever LSF state SkyW's decomposition fit left behind.
    dlam_pix = float(np.median(np.diff(decomp.wave)))
    sigma_pix = row_data["lsf_fwhm"] / 2.355 / dlam_pix
    taps = np.asarray(LSF_TAP_OFFSETS, dtype=np.float64)
    kernel = norm.pdf(taps[None, :], loc=0.0, scale=sigma_pix[:, None])
    kernel /= kernel.sum(axis=1, keepdims=True)
    decomp.lsf_surface_state = None
    decomp._lsf_surface = kernel
    decomp._lsf_operator = build_lsf_operator(decomp.wave, kernel)
    mats = decomp._assemble_refined_matrices()

    def _sum_components(coef_vec, keys=COMPONENT_KEYS):
        comps = decomp._components_from_coef(coef_vec, mats)
        total = np.zeros_like(row_data["wave"], dtype=np.float64)
        for key in keys:
            arr = comps.get(key)
            if arr is not None:
                total = total + np.asarray(arr, dtype=np.float64)
        return total

    # Predicted coef is in the same FACTOR-scaled fit space decompose_parallel.py
    # trains against (flux_row * FACTOR); divide back out so flux_pred/line_pred
    # land in the same physical erg/s/cm^2/A units as FLUX_OBS. Continuum is
    # FLUX_PRED - LINE_PRED, not computed/stored separately.
    flux_pred = _sum_components(coef) / FACTOR
    flux_pred_lo = _sum_components(np.maximum(coef - coef_std, 0.0)) / FACTOR
    flux_pred_hi = _sum_components(np.maximum(coef + coef_std, 0.0)) / FACTOR
    line_pred = _sum_components(coef, keys=LINE_KEYS) / FACTOR

    return dict(
        wave=row_data["wave"],
        flux_pred=flux_pred, flux_pred_lo=flux_pred_lo, flux_pred_hi=flux_pred_hi,
        line_pred=line_pred,
        coef=coef, coef_std=coef_std, confidence=confidence,
        coef_names=list(result["coef_names"]),
        fit_east=fit_east, fit_west=fit_west,
        near_is_east=bool(result["near_is_east"][0]),
    )


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_output(row_data, prediction, fits_file, outpath=None):
    """Write observed + predicted spectra and coefficients to a FITS file.

    Parameters
    ----------
    row_data : dict
        Output of `read_row`.
    prediction : dict
        Output of `predict_row`.
    fits_file : str or Path
        Source XCframe file (recorded in the output header for provenance).
    outpath : str or Path, optional
        Output path (default: PredictedSky_<expnum>.fits).

    Returns
    -------
    Path
        The path written.
    """
    if outpath is None:
        outpath = f"PredictedSky_{row_data['expnum']}.fits"
    outpath = Path(outpath)

    hdr = fits.Header()
    hdr["TITLE"] = "PredictedSky"
    hdr["INPUT"] = (str(fits_file), "source XCframe summary file")
    hdr["ROW"] = (row_data["row"], "row index in source file")
    hdr["EXPNUM"] = (row_data["expnum"], "exposure number")
    hdr["MJD"] = (row_data["mjd"], "UT MJD of exposure")
    hdr["CONFID"] = (prediction["confidence"], "ensemble prediction confidence (0,1]")
    hdr["NEARE"] = (prediction["near_is_east"], "True if SkyE was the near arm")
    hdr["SCIRA"] = (row_data["sci_ra"], "Sci pointing RA (deg)")
    hdr["SCIDEC"] = (row_data["sci_dec"], "Sci pointing Dec (deg)")

    coef_tab = Table(
        [prediction["coef_names"], prediction["coef"], prediction["coef_std"]],
        names=["Name", "Coef", "Coef_Err"],
    )

    fits.HDUList([
        fits.PrimaryHDU(header=hdr),
        fits.ImageHDU(data=prediction["wave"].astype(np.float32), name="WAVE"),
        fits.ImageHDU(data=row_data["flux_sci"].astype(np.float32), name="FLUX_OBS"),
        fits.ImageHDU(data=prediction["flux_pred"].astype(np.float32), name="FLUX_PRED"),
        fits.ImageHDU(data=prediction["flux_pred_lo"].astype(np.float32), name="FLUX_PRED_LO"),
        fits.ImageHDU(data=prediction["flux_pred_hi"].astype(np.float32), name="FLUX_PRED_HI"),
        fits.ImageHDU(data=prediction["line_pred"].astype(np.float32), name="LINE_PRED"),
        fits.BinTableHDU(coef_tab, name="COEF"),
    ]).writeto(outpath, overwrite=True)

    print(f"\nOutput written to {outpath}")
    return outpath


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        parents=[_pre],
        description=("Predict the Sci-pointing sky spectrum for a single "
                     "exposure from an XCframe summary file."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("fits_file", help="LVM XCframe summary FITS file")
    p.add_argument("-row", type=int, default=None,
                   help="Row index to select (default: 0)")
    p.add_argument("-expnum", type=int, default=None,
                   help="Select by DRP_ALL expnum instead of row index")
    p.add_argument("-model", required=True,
                   help="Trained ensemble .pt archive (required, no default)")
    p.add_argument("-output", default=None,
                   help="Output FITS path (default: PredictedSky_<expnum>.fits)")
    args = p.parse_args()

    row_data = read_row(args.fits_file, row=args.row, expnum=args.expnum)
    print(f"Source: {args.fits_file}")
    print(f"  row={row_data['row']}  expnum={row_data['expnum']}  "
          f"mjd={row_data['mjd']:.5f}")
    print(f"  sci=({row_data['sci_ra']:.4f}, {row_data['sci_dec']:.4f})  "
          f"skyE=({row_data['skye_ra']:.4f}, {row_data['skye_dec']:.4f})  "
          f"skyW=({row_data['skyw_ra']:.4f}, {row_data['skyw_dec']:.4f})")

    prediction = predict_row(row_data, model_path=args.model)

    print(f"  confidence = {prediction['confidence']:.3f}  "
          f"near_arm={'SkyE' if prediction['near_is_east'] else 'SkyW'}")
    rel = prediction["coef_std"] / np.maximum(np.abs(prediction["coef"]), 1e-12)
    order = np.argsort(rel)
    worst, best = order[-3:][::-1], order[:3]
    names = prediction["coef_names"]
    c, std = prediction["coef"], prediction["coef_std"]
    print("  most uncertain coefs: " + ", ".join(
        f"{names[k]}={c[k]:+.3g}±{std[k]:.2g}" for k in worst))
    print("  most certain coefs:   " + ", ".join(
        f"{names[k]}={c[k]:+.3g}±{std[k]:.2g}" for k in best))

    write_output(row_data, prediction, args.fits_file, outpath=args.output)


if __name__ == "__main__":
    main()
