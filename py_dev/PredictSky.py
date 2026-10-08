#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Predict the sky spectrum at the Sci telescope's pointing for a single
    exposure, using a trained semi-empirical machine-learning sky model
    (lvmsky's mlp_ensemble_split_zodi ensemble, from TrainSkyModel.py) and
    the SkyE/SkyW spectra recorded in an LVM XCframe summary FITS file.
    Written for lvmsky branch skydecomp-telluric-corrected-lines.

Command line usage (if any):

    usage: PredictSky.py [-h] [-lvmsky_skysub PATH] [-lvmcore_dir PATH]
                         [-row N | -expnum N] [-output PATH]
                         [-no_arm_correction] [-no_line_scaling]
                         model fits_file

    where

    -row N          selects the exposure by row index into the file
                    (default: 0).

    -expnum N       selects the exposure by its DRP_ALL 'expnum' value
                    instead of a row index (overrides -row).

    -output PATH    output FITS path (default: PredictedSky_<expnum>.fits).

    -no_arm_correction
                    do not apply the sky-arm residual correction (see
                    Description, step 5).

    -no_line_scaling
                    do not rescale the predicted sky lines on the science
                    spectrum (see Description, step 6).

    -lvmsky_skysub PATH
                    path to the lvmsky repo's skysub/ directory, which
                    supplies decompose_parallel.py and the mlp_predictor
                    and sky_decomp packages (default: ~/SDSS/lvmsky/skysub).

    -lvmcore_dir PATH
                    sets LVMCORE_DIR, needed by mlp_predictor.data for the
                    LCO extinction curve (default: ~/SDSS/lvmcore).

    model           path to the trained ensemble .pt archive. Positional
                    and required, no default -- the model is always named
                    explicitly.

    fits_file       an LVM XCframe summary FITS file (WAVE, FLUX, SKY_EAST,
                    SKY_WEST, LSF, DRP_ALL), as from SummarizeCframe.py
                    -by fiber, preferably selected with SelectXCF.py -hdr so
                    that DRP_ALL carries each exposure's PWV.

Description:

    For the selected exposure, following lvmsky's
    notebook_example_predict_sky.ipynb:

      1. Builds a one-row decomposition stack with
         ConvertForDecompose.build_stack -- the same code that builds the
         training stacks.
      2. Decomposes SkyE and SkyW with
         decompose_parallel.decompose_in_process, i.e. with exactly the
         fit model, constraints, telluric transmission (PWV, airmass),
         photon weights and science-line mask the training corpus was
         decomposed with.  The coefficients fed to the network must come
         from the same decomposition it was trained on.
      3. Predicts the Sci-pointing coefficients with
         inference.predict_sky_from_minimal_inputs, giving it the
         wavelength grid, the LSFs and date_obs, which the physical
         moon/zodiacal-light context features need.
      4. Rebuilds the predicted sky on the training basis
         (data.make_reconstruction_decomposer, with this exposure's
         telluric transmission along the science line of sight and the
         science LSF; the O2 band shape is taken from the near arm's fit).
      5. Sky-arm residual correction (sky_arm_correction): both sky arms'
         residuals against their own fits, mostly solar absorption lines
         in the blue that the moon/zodiacal-light basis misses, are
         combined and added to the prediction (below 5000 A, tapered).
      6. Sky-line scaling (sky_line_scaling): one scale per OH band, per
         atomic line family and for O2 is fitted on the science spectrum
         itself after a high-pass, pulled toward the prediction; Na D and
         K I are left as predicted.
      7. Writes the observed and predicted spectra, the corrections and
         the coefficients to a FITS file.

    Output FITS structure::

        WAVE           wavelength [A]
        FLUX_OBS       observed science spectrum [erg/s/cm^2/A]
        FLUX_PRED      predicted sky, including any corrections
        FLUX_PRED_LO/HI
                       FLUX_PRED with the coefficients moved down/up by
                       their ensemble spread
        FLUX_PRED_RAW  predicted sky before the corrections
        LINE_PRED      sky-line part of FLUX_PRED (OH, atomic, O2 lines,
                       including the line-scaling correction)
        SKY_CORR       sky-arm residual correction (zero if not applied)
        LINE_CORR      sky-line scaling correction (zero if not applied)
        COEF           predicted coefficients and their ensemble spread

Primary routines:

    read_row
    build_decomp
    decompose_row
    predict_row
    write_output

Notes::

    - Only the science fibers' LSF is in the XCframe; it is used for the
      sky telescopes too (as in training -- see ConvertForDecompose.py).
    - Photon-noise weights in steps 5-6 assume the summary spectra combine
      n_sci_fibers science fibers (DRP_ALL) and SKY_ARM_FIBRES sky fibers.
    - The science spectrum in a -by fiber XCframe includes the target's
      own light; the line scaling only uses the narrow sky lines after a
      high-pass, so it is not affected by the continuum.
    - This script imports from an external lvmsky checkout
      (-lvmsky_skysub).  The model must have been trained against the same
      lvmsky version (see TrainSkyModel.py).

History::

    260901  ksl  Coding begun.
    260903  ksl  Switched every option from double-dash (--lvmsky-skysub)
        to single-dash (-lvmsky_skysub), matching py_progs/'s convention
        -- see BatchPredictSkyESO.py's History for the fuller note; this
        also required fixing BatchPredictSky.py's worker-init sys.argv
        injection, which fakes a command line for this module's own
        module-level _pre parser and had hardcoded the old spelling.
    260904  ksl  -model switched from a required dashed option to a
        plain positional argument, matching py_progs/'s convention that
        required inputs are positional and only truly optional settings
        get a -flag (does not affect BatchPredictSky.py's worker-init
        sys.argv injection, which only drives the module-level _pre
        parser above, never this main() parser). Positional order is
        model then fits_file.
    261008  ksl  Rewritten for lvmsky branch skydecomp-telluric-corrected-lines,
        following its notebook_example_predict_sky.ipynb: SkyE/SkyW are
        decomposed with decompose_parallel.decompose_in_process on a
        one-row stack (ConvertForDecompose.build_stack), i.e. exactly as
        the training corpus; the prediction is given the wavelength grid,
        LSFs and date_obs; the sky is rebuilt on the telluric basis
        (data.make_reconstruction_decomposer); the sky-arm residual
        correction and sky-line scaling are applied by default
        (-no_arm_correction / -no_line_scaling).  Output gains SKY_CORR,
        LINE_CORR and FLUX_PRED_RAW.


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

from mlp_predictor import serialization, inference, data as mp_data  # noqa: E402
from mlp_predictor import wavelengths as _wave_mod  # noqa: E402
from mlp_predictor.sky_arm_correction import (  # noqa: E402
    SkyArm, arm_photon_variance, sky_arm_residual_correction)
from mlp_predictor.sky_line_scaling import (  # noqa: E402
    line_templates, sky_line_scaling_correction)
from sky_decomp.moon_zodi_model import LSF_FWHM_TO_SIGMA  # noqa: E402
import decompose_parallel as _dp  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ConvertForDecompose  # noqa: E402  -- sibling script in py_dev/

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

FACTOR = 1e14  # XCframe flux is erg/s/cm^2/A; the decomposition works in O(1) units.
# Must match TrainSkyModel.py's DECOMP_VARIANT (the basis the model was trained on).
DECOMP_VARIANT = "telluric-palacecorr"
DECOMP_SUFFIX = mp_data.DECOMP_VARIANTS[DECOMP_VARIANT]["suffix"]
COMPONENT_KEYS = ("oh", "moon", "zodi", "diffuse", "atom", "orc", "o2")
# Line-like components (OH, atomic, O I recombination, O2 band); the rest
# (moonlight, zodiacal light, airglow continua) is the continuum.
LINE_KEYS = ("oh", "atom", "orc", "o2")
# Fibers combined into an XCframe sky-telescope spectrum, for photon weights.
SKY_ARM_FIBRES = 50
# Post-prediction correction settings: the defaults of lvmsky's
# notebook_example_predict_sky.ipynb.
RESIDUAL_SMOOTHING_A = "auto"
RESIDUAL_CORRECTION_MAX_A = 5000.0
LINE_SCALING_HIGHPASS_A = 25.0
LINE_SCALING_PRIOR = 0.1
LINE_SCALING_TILT = True


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
        Keys: row, wave, flux_sci, flux_e, flux_w, lsf_raw, drp (one-row
        Table), mjd, date_obs, sci_ra, sci_dec, skye_ra, skye_dec, skyw_ra,
        skyw_dec, expnum, n_sci_fibers.
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
        lsf_raw = np.asarray(hdul["LSF"].data[i], dtype=np.float64)
        drp_row = drp[i:i + 1]

    date_obs = str(drp_row["obstime"][0]).strip()
    n_sci = (int(drp_row["n_sci_fibers"][0]) if "n_sci_fibers" in drp_row.colnames
             else None)
    return dict(
        row=i,
        wave=wave,
        flux_sci=flux_sci, flux_e=flux_e, flux_w=flux_w, lsf_raw=lsf_raw,
        drp=drp_row,
        mjd=Time(date_obs, format="isot", scale="utc").mjd,
        date_obs=date_obs,
        sci_ra=float(drp_row["sci_ra"][0]), sci_dec=float(drp_row["sci_dec"][0]),
        skye_ra=float(drp_row["skye_ra"][0]), skye_dec=float(drp_row["skye_dec"][0]),
        skyw_ra=float(drp_row["skyw_ra"][0]), skyw_dec=float(drp_row["skyw_dec"][0]),
        expnum=int(drp_row["expnum"][0]),
        n_sci_fibers=n_sci,
    )


# ---------------------------------------------------------------------------
# Decomposition
# ---------------------------------------------------------------------------

def build_decomp(wave, ensemble):
    """Settings for rebuilding spectra from this ensemble's coefficients.

    Returns a dict (spline knot counts and the PALACE data directory)
    that predict_row needs; pass it as predict_row's decomp argument in a
    batch loop so it is worked out once.  The decomposition itself is per
    exposure (its telluric transmission depends on PWV and airmass), so
    nothing heavier can be shared.
    """
    n_moon_knots, split_zodi, n_zodi_knots = _wave_mod.infer_spline_knots(
        ensemble["coef_names"])
    return dict(n_moon_knots=n_moon_knots, split_zodi=split_zodi,
                n_zodi_knots=n_zodi_knots,
                base_dir=mp_data._infer_base_dir_for_reconstruction())


def decompose_row(row_data):
    """Decompose one exposure's SkyE and SkyW exactly as the training corpus.

    Builds a one-row stack (ConvertForDecompose.build_stack) and runs
    decompose_parallel.decompose_in_process on its near and far sky arms.

    Returns
    -------
    dict
        Keys: stack (the in-memory HDUList), lsf (filled science LSF),
        fit_near, fit_far, fit_east, fit_west, near_is_east, flags_near,
        flags_far.
    """
    stack = ConvertForDecompose.build_stack(
        row_data["wave"], row_data["flux_sci"], row_data["flux_e"],
        row_data["flux_w"], row_data["lsf_raw"], row_data["drp"], verbose=True)
    lsf = np.asarray(stack["LSF_SCI"].data[0], dtype=np.float64)
    near_is_east = str(stack["META"].data["sky_near_label"][0]).strip() == "SkyE"

    fits_ = _dp.decompose_in_process(stack, rows=(0,), kinds=("sky1", "sky2"))
    fit_near, flags_near = fits_[("sky1", 0)]
    fit_far, flags_far = fits_[("sky2", 0)]
    for label, fit, flags in (("near", fit_near, flags_near), ("far", fit_far, flags_far)):
        print(f"  {label} ({'SkyE' if (label == 'near') == near_is_east else 'SkyW'}) "
              f"decomp: status={fit.fit_status}, chi2_red={fit.reduced_chi2:.3f} "
              f"(photon-weighted), R^2={fit.r2:.4f}, "
              f"reliability={flags.get('reliability') if flags else None}")
    fit_east, fit_west = (fit_near, fit_far) if near_is_east else (fit_far, fit_near)
    return dict(stack=stack, lsf=lsf, fit_near=fit_near, fit_far=fit_far,
                fit_east=fit_east, fit_west=fit_west, near_is_east=near_is_east,
                flags_near=flags_near, flags_far=flags_far)


# ---------------------------------------------------------------------------
# Prediction
# ---------------------------------------------------------------------------

def _sum_components(comps, keys=COMPONENT_KEYS):
    total = np.zeros_like(np.asarray(comps["oh"], dtype=np.float64))
    for key in keys:
        arr = comps.get(key)
        if arr is not None:
            total = total + np.asarray(arr, dtype=np.float64)
    return total


def predict_row(row_data, model_path=None, ensemble=None, decomp=None,
                arm_correction=True, line_scaling=True):
    """Predict the Sci-pointing sky for one exposure.

    Parameters
    ----------
    row_data : dict
        Output of `read_row`.
    model_path : str or Path, optional
        Path to the trained ensemble .pt archive.  Ignored if `ensemble` is
        given; required otherwise.
    ensemble : dict, optional
        An already-loaded ensemble (`serialization.load_ensemble`), so a
        batch loop loads it once.
    decomp : dict, optional
        build_decomp's output, so a batch loop works it out once.
    arm_correction, line_scaling : bool
        Apply the sky-arm residual correction / the sky-line scaling.

    Returns
    -------
    dict
        Keys: wave, flux_pred, flux_pred_lo, flux_pred_hi, flux_pred_raw,
        line_pred, sky_corr, line_corr, coef, coef_std, confidence,
        reliability, coef_names, fit_east, fit_west, near_is_east,
        line_scales.
    """
    if ensemble is None:
        if model_path is None:
            raise ValueError("predict_row: either model_path or ensemble must be given")
        ensemble = serialization.load_ensemble(str(model_path))
    wave = row_data["wave"]
    if decomp is None:
        decomp = build_decomp(wave, ensemble)

    dec = decompose_row(row_data)
    lsf = dec["lsf"]
    if dec["fit_east"].coef.size != len(ensemble["coef_names"]):
        raise ValueError(
            f"decomposition has {dec['fit_east'].coef.size} coefficients but the "
            f"model expects {len(ensemble['coef_names'])}: the model was trained "
            f"on a different decomposition")

    result = inference.predict_sky_from_minimal_inputs(
        ensemble,
        obstime_mjd=[row_data["mjd"]],
        sci_ra=[row_data["sci_ra"]], sci_dec=[row_data["sci_dec"]],
        sky_e_ra=[row_data["skye_ra"]], sky_e_dec=[row_data["skye_dec"]],
        sky_w_ra=[row_data["skyw_ra"]], sky_w_dec=[row_data["skyw_dec"]],
        coef_e=dec["fit_east"].coef.astype(np.float32)[None, :],
        coef_w=dec["fit_west"].coef.astype(np.float32)[None, :],
        wave=wave, lsf_e=lsf, lsf_w=lsf, lsf_sci=lsf,
        date_obs=row_data["date_obs"], expnum=row_data["expnum"],
    )
    coef = result["coef"][0].astype(np.float64)
    coef_std = result["coef_std"][0].astype(np.float64)

    # Rebuild on the training basis, with this exposure's telluric
    # transmission along the science line of sight and the science LSF.
    meta = Table(dec["stack"]["META"].data)
    tel_sci = mp_data.telluric_row_kwargs(
        meta, 0, "sci", wave, lsf,
        palace_oh_suffix=mp_data.palace_oh_suffix_for(DECOMP_SUFFIX))
    recon = mp_data.make_reconstruction_decomposer(
        wave, n_spline_knots=decomp["n_moon_knots"], base_dir=decomp["base_dir"],
        split_zodi=decomp["split_zodi"], n_zodi_spline_knots=decomp["n_zodi_knots"],
        telluric=tel_sci, lsf_sigma=lsf / LSF_FWHM_TO_SIGMA)
    recon._set_lsf_state(recon._nominal_state("predictsky_sci_nominal_lsf"))
    mats = recon._assemble_refined_matrices()
    # The O2 band template is a per-row fitted shape; use the near arm's.
    mats["o2"] = np.asarray(dec["fit_near"].vector_o2, dtype=np.float64).ravel()[None, :]

    comps = recon._components_from_coef(coef, mats)
    flux_raw = _sum_components(comps) / FACTOR
    flux_lo = _sum_components(recon._components_from_coef(np.maximum(coef - coef_std, 0), mats)) / FACTOR
    flux_hi = _sum_components(recon._components_from_coef(np.maximum(coef + coef_std, 0), mats)) / FACTOR
    line_raw = _sum_components(comps, LINE_KEYS) / FACTOR
    rev_bits, _rev_info = inference.reliability_from_components(comps, recon.wave)
    reliability = int(result["reliability"][0]) | int(rev_bits)

    sky_corr = np.zeros_like(flux_raw)
    if arm_correction:
        def _arm(flux_obs, fit):
            model = np.asarray(fit.bestfit_lsf, dtype=np.float64) / FACTOR
            solar = (np.asarray(fit.components["moon"], dtype=np.float64)
                     + np.asarray(fit.components.get("zodi", 0.0), dtype=np.float64)) / FACTOR
            return SkyArm(observed=flux_obs, model=model, solar=solar,
                          variance=arm_photon_variance(flux_obs, wave, [SKY_ARM_FIBRES]))
        solar_pred = (np.asarray(comps["moon"], dtype=np.float64)
                      + np.asarray(comps.get("zodi", 0.0), dtype=np.float64)) / FACTOR
        sky_corr = sky_arm_residual_correction(
            wave, flux_raw, solar_pred,
            [_arm(row_data["flux_e"], dec["fit_east"]),
             _arm(row_data["flux_w"], dec["fit_west"])],
            max_wavelength=RESIDUAL_CORRECTION_MAX_A, smoothing=RESIDUAL_SMOOTHING_A)

    line_corr = np.zeros_like(flux_raw)
    line_scales = {}
    if line_scaling:
        templates = {k: v / FACTOR for k, v in line_templates(
            recon, mats, coef, tilt=LINE_SCALING_TILT).items()}
        flux_near = row_data["flux_e"] if dec["near_is_east"] else row_data["flux_w"]
        sci_mask = mp_data.science_line_mask_rows(
            dec["stack"], wave, row_data["flux_sci"], flux_near)
        n_sci = row_data["n_sci_fibers"]
        line_corr, info = sky_line_scaling_correction(
            wave, row_data["flux_sci"], flux_raw + sky_corr, templates,
            variance=arm_photon_variance(row_data["flux_sci"], wave,
                                         None if n_sci is None else [n_sci]),
            mask=sci_mask, highpass_A=LINE_SCALING_HIGHPASS_A,
            prior_sigma=LINE_SCALING_PRIOR, return_info=True)
        line_scales = {k: v for k, v in info["scales"].items() if not k.endswith("_tilt")}

    corr = sky_corr + line_corr
    return dict(
        wave=wave,
        flux_pred=flux_raw + corr,
        flux_pred_lo=flux_lo + corr, flux_pred_hi=flux_hi + corr,
        flux_pred_raw=flux_raw,
        line_pred=line_raw + line_corr,
        sky_corr=sky_corr, line_corr=line_corr,
        coef=coef, coef_std=coef_std,
        confidence=float(result["confidence"][0]),
        reliability=reliability,
        coef_names=list(result["coef_names"]),
        fit_east=dec["fit_east"], fit_west=dec["fit_west"],
        near_is_east=bool(result["near_is_east"][0]),
        line_scales=line_scales,
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
    conf = prediction["confidence"]
    hdr["CONFID"] = (float(conf) if np.isfinite(conf) else -1.0,
                     "ensemble confidence (0,1]; -1 if undefined (1 seed)")
    hdr["RELIAB"] = (prediction["reliability"], "reliability bitmask (0 = clean)")
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
        fits.ImageHDU(data=prediction["flux_pred_raw"].astype(np.float32), name="FLUX_PRED_RAW"),
        fits.ImageHDU(data=prediction["line_pred"].astype(np.float32), name="LINE_PRED"),
        fits.ImageHDU(data=prediction["sky_corr"].astype(np.float32), name="SKY_CORR"),
        fits.ImageHDU(data=prediction["line_corr"].astype(np.float32), name="LINE_CORR"),
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
    p.add_argument("-row", type=int, default=None,
                   help="Row index to select (default: 0)")
    p.add_argument("-expnum", type=int, default=None,
                   help="Select by DRP_ALL expnum instead of row index")
    p.add_argument("-output", default=None,
                   help="Output FITS path (default: PredictedSky_<expnum>.fits)")
    p.add_argument("-no_arm_correction", action="store_true",
                   help="Do not apply the sky-arm residual correction")
    p.add_argument("-no_line_scaling", action="store_true",
                   help="Do not rescale the predicted sky lines on the science spectrum")
    p.add_argument("model", help="Trained ensemble .pt archive")
    p.add_argument("fits_file", help="LVM XCframe summary FITS file")
    args = p.parse_args()

    row_data = read_row(args.fits_file, row=args.row, expnum=args.expnum)
    print(f"Source: {args.fits_file}")
    print(f"  row={row_data['row']}  expnum={row_data['expnum']}  "
          f"mjd={row_data['mjd']:.5f}")
    print(f"  sci=({row_data['sci_ra']:.4f}, {row_data['sci_dec']:.4f})  "
          f"skyE=({row_data['skye_ra']:.4f}, {row_data['skye_dec']:.4f})  "
          f"skyW=({row_data['skyw_ra']:.4f}, {row_data['skyw_dec']:.4f})")

    prediction = predict_row(row_data, model_path=args.model,
                             arm_correction=not args.no_arm_correction,
                             line_scaling=not args.no_line_scaling)

    print(f"  confidence = {prediction['confidence']:.3f}  "
          f"reliability = {prediction['reliability']}  "
          f"near_arm={'SkyE' if prediction['near_is_east'] else 'SkyW'}")
    if prediction["line_scales"]:
        print("  line scales: " + ", ".join(
            f"{k}={v:.3f}" for k, v in prediction["line_scales"].items()))
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
