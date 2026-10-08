#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Consolidated, resumable driver that turns a SelectXCF.py-selected
    XCframe corpus into a trained semi-empirical machine-learning sky
    model (lvmsky's mlp_ensemble_split_zodi ensemble), running every
    intermediate step in one call: ConvertForDecompose.py, lvmsky's
    decompose_parallel.py, the triplet build/filter/augment step, and
    the network training.  Written for lvmsky branch
    skydecomp-telluric-corrected-lines (2026-10).

Command line usage (if any):

    usage: TrainSkyModel.py [-h] [-work_dir DIR] [-start_at STAGE]
                            [-stop_after STAGE] [-force] [-np N]
                            [-atom_k_max ATOM_K_MAX] [-epochs EPOCHS]
                            [-seeds SEEDS] [-flux_mse_groups GROUPS]
                            [-train_frac F] [-val_frac F] [-split_seed SEED]
                            [-n_bins N] [-output PATH]
                            [-lvmsky_skysub PATH] [-lvmcore_dir PATH]
                            xcf_fits

    where

    xcf_fits        is an XCframe-layout FITS file (SelectXCF.py output),
                    the corpus of exposures to train on.  Select it with
                    SelectXCF.py -hdr so that DRP_ALL carries each
                    exposure's PWV (see Notes).

    -work_dir DIR   directory for every pipeline artifact -- decomp_input
                    FITS, decompose_parallel.py outputs, filtered_triplet
                    pickle, trained checkpoint (default: <xcf_fits
                    stem>_train/, created if missing).

    -start_at STAGE first stage to (re)run: convert, decompose, wavecache
                    or train (default: convert). Use this to resume a
                    partial run -- e.g. after running the decompose stage
                    separately with a different -np, or after this
                    script's own decompose stage completed but a later
                    stage failed.

    -stop_after STAGE
                    last stage to run, same choices as -start_at (default:
                    train).

    -force          redo every stage from -start_at onward even if its
                    output files already exist (default: skip a stage
                    whose expected output is already present).

    -np N           decompose_parallel.py worker processes (default: 4,
                    matching its own default; matches py_progs/Reduce.py
                    and py_progs/sky_gaussfit.py's -np convention for
                    process count).

    -atom_k_max ATOM_K_MAX
                    upper hard bound for the atom_k coefficient filter,
                    passed to apply_triplet_filters' hard_coef_bounds
                    (default: 10.01, as in lvmsky's training notebook).

    -epochs EPOCHS  override n_epochs (default: 400, from
                    TRAIN_CFG_OVERRIDES; pass a small value for a quick
                    smoke-test run).

    -seeds SEEDS    comma-separated ensemble seeds, e.g. "42,43" (default:
                    unset, uses the config's own 10-seed production
                    ensemble).

    -flux_mse_groups GROUPS
                    comma-separated group names for the flux-MSE loss
                    term (default: moon,zodi,continuum, from
                    TRAIN_CFG_OVERRIDES).

    -train_frac F   moon-phase-stratified split fraction for training
                    (default: 0.8, mlp_predictor's own default).

    -val_frac F     moon-phase-stratified split fraction for validation;
                    the remainder is held out as test (default: 0.1).

    -split_seed SEED
                    RNG seed for the moon-phase-stratified split (default:
                    42).

    -n_bins N       number of moon-phase bins used to stratify the split
                    (default: 10).

    -output PATH    trained ensemble .pt path (default: <work_dir>/
                    mlp_ensemble.pt).

    -lvmsky_skysub PATH
                    path to the lvmsky repo's skysub/ directory, which
                    supplies decompose_parallel.py and the mlp_predictor
                    package this script imports (default:
                    ~/SDSS/lvmsky/skysub).

    -lvmcore_dir PATH
                    sets LVMCORE_DIR, needed by mlp_predictor.data for the
                    LCO extinction curve (default: ~/SDSS/lvmcore).

Description:

    Runs four stages in order, each independently resumable (skipped if
    its expected output files already exist, unless -force is given).
    Stages 3 and 4 follow lvmsky's own training notebook
    (notebook_sky_interpolation_triplet_dual_encoder_group_mlp_split_zodi_
    module.ipynb):

      1. convert    -- ConvertForDecompose.convert(): reformats xcf_fits
                       into decompose_parallel.py's stack layout (FLUX_*,
                       LSF_*, META with airmasses, date_obs and pwv_med).
      2. decompose  -- runs lvmsky's decompose_parallel.py as a subprocess
                       (--fit-model palacecorr-aijc-vnf-split-zodi-lsf-
                       spline2d, the branch default; -np workers),
                       producing the decomp_{sky1,sky2,sci},
                       {sky1,sky2,sci}_meta_coef, meta_only and every10
                       FITS files.
      3. wavecache  -- builds the (near sky, far sky, science) coefficient
                       triplet, applies the notebook's quality filters
                       (coefficient bounds, LMC/SMC exclusion, moon/zodi
                       role reversal, science colour excess, collapsed
                       airglow continuum), adds the ecliptic, physics-prior
                       and physical-moon-model context features (the last
                       builds a cache beside the corpus on first use), and
                       resolves each coefficient's wavelength and
                       extinction; pickles the result as
                       work_dir/filtered_triplet.pkl.
      4. train      -- moon-phase-stratified train/val/test split,
                       per-group compressors, then the seed ensemble
                       (Trainer.run_ensemble) with the notebook's settings
                       (TRAIN_CFG_OVERRIDES), saved with save_ensemble.
                       Also writes work_dir/filtered_triplet_augmented.pkl
                       (the corpus EvalCoefResiduals.py reads).

    The B-spline knot counts (n_moon_knots, split_zodi, n_zodi_knots) that
    run_ensemble needs are inferred from the decomposed coef_names via
    mlp_predictor.wavelengths.infer_spline_knots.

Notes::

    - PWV: the branch's telluric decomposition recomputes the telluric
      transmission the DRP applied, which used PWV_MED from each CFrame
      header.  PWV is not in drpall, so it comes from a
      SummarizeSkyHdr.py file via SelectXCF.py -hdr.  Rows without a
      valid pwv_med fall back to the DRP default of 15 mm, which is
      wrong for most exposures; the convert stage warns how many.
    - LSF: the XCframe carries only the science fibers' LSF, which is
      used for the sky telescopes too, and its NaN pixels at the ends of
      the wavelength range are filled (see ConvertForDecompose.py).
    - TRAIN_CFG_OVERRIDES copies the settings in which the notebook's
      train_cfg cell (lvmsky commit 59687bb) differs from
      mlp_predictor.trainer.default_dual_group_config.  Recheck it when
      lvmsky is updated.
    - The decompose stage runs about 1.5 s per spectrum per worker (three
      spectra per exposure); the moon-model cache about 0.55 s per
      exposure per worker.
    - decompose is the only stage run as a subprocess rather than an
      in-process call: it is lvmsky's own long-running, CPU-thread-
      pinned, parallel script, and keeping it separate preserves its own
      progress bar, thread limits and -np tuning.
    - No evaluation stage is included by design -- run EvalCoefResiduals.py
      on work_dir and/or BatchPredictSky.py + EvalFluxResiduals.py
      against the saved model afterwards.
    - mlp_predictor's filters build Plotly diagnostic histograms and call
      fig.show(), which would open browser tabs.  Patched at import time
      here (not in lvmsky) so the figures are written to
      work_dir/plots/*.html instead.

History::

    260904  ksl  Coding begun -- consolidates ConvertForDecompose.py,
        lvmsky's decompose_parallel.py, and the wavecache/train logic
        previously duplicated between niv/stage1_wavecache.py and
        niv/stage1_train.py + niv/stage2_train.py into one resumable
        driver.
    261007  ksl  Updated for lvmsky branch skydecomp-telluric-corrected-lines:
        branch default fit model and output suffix; the notebook's
        filters, context augments (incl. the physical moon model) and
        training settings (TRAIN_CFG_OVERRIDES); run_ensemble given the
        flux stack and decomposition suffix; PWV taken from DRP_ALL
        (SelectXCF.py -hdr); default split 0.8/0.1.

'''

import argparse
import copy
import os
import pickle
import subprocess
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

# ---------------------------------------------------------------------------
# External package setup -- mlp_predictor lives in the lvmsky repo, not in
# this project.  Inserted before argparse runs so -lvmsky_skysub /
# -lvmcore_dir can override the defaults before the imports below fire.
# ---------------------------------------------------------------------------

DEFAULT_LVMSKY_SKYSUB = Path("~/SDSS/lvmsky/skysub").expanduser()
DEFAULT_LVMCORE_DIR = Path("~/SDSS/lvmcore").expanduser()

_pre = argparse.ArgumentParser(add_help=False)
_pre.add_argument("-lvmsky_skysub", default=str(DEFAULT_LVMSKY_SKYSUB))
_pre.add_argument("-lvmcore_dir", default=str(DEFAULT_LVMCORE_DIR))
_pre_args, _ = _pre.parse_known_args()

sys.path.insert(0, _pre_args.lvmsky_skysub)
os.environ.setdefault("LVMCORE_DIR", _pre_args.lvmcore_dir)

from mlp_predictor.config import PipelineConfig  # noqa: E402
from mlp_predictor import data as mp_data  # noqa: E402
from mlp_predictor import wavelengths as mp_wave  # noqa: E402
from mlp_predictor import compressor as mp_comp  # noqa: E402
from mlp_predictor import trainer as mp_trainer  # noqa: E402
from mlp_predictor import serialization as mp_ser  # noqa: E402
from mlp_predictor import ml_utils as mp_ml  # noqa: E402

import ConvertForDecompose  # noqa: E402  -- sibling script in py_dev/

# ---------------------------------------------------------------------------
# mlp_predictor.data.apply_triplet_filters() (called by stage_wavecache)
# builds a couple of Plotly diagnostic histograms and calls fig.show() on
# each, which pops a new browser tab per Plotly's default non-notebook
# renderer -- undesirable for a batch/CLI driver, and not something we
# touch lvmsky's own source to fix. Patched here, at the BaseFigure level
# (covers both plotly.express and plotly.graph_objects figures), so every
# .show() call anywhere in the imported lvmsky code writes its figure to
# work_dir/plots/ instead of opening a browser. _PLOT_DIR is repointed at
# the real work_dir once main() knows it.
# ---------------------------------------------------------------------------
import re  # noqa: E402
import plotly.basedatatypes as _pbd  # noqa: E402

_PLOT_DIR = Path(".")
_plot_counter = {"n": 0}


def _save_instead_of_show(self, *_args, **_kwargs):
    _plot_counter["n"] += 1
    try:
        title = self.layout.title.text
    except Exception:
        title = None
    slug = re.sub(r"[^a-zA-Z0-9]+", "_", title).strip("_").lower() if title else "figure"
    _PLOT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = _PLOT_DIR / f"{_plot_counter['n']:02d}_{slug}.html"
    self.write_html(str(out_path), auto_open=False)
    print(f"[plot] saved diagnostic figure to {out_path} (browser auto-open disabled)")


_pbd.BaseFigure.show = _save_instead_of_show

STAGES = ["convert", "decompose", "wavecache", "train"]

# lvmsky branch skydecomp-telluric-corrected-lines: decompose_parallel.py's
# default fit model (2026-09-24) and the matching mlp_predictor variant.
FIT_MODEL = "palacecorr-aijc-vnf-split-zodi-lsf-spline2d"
DECOMP_VARIANT = "telluric-palacecorr"
DECOMP_SUFFIX = mp_data.DECOMP_VARIANTS[DECOMP_VARIANT]["suffix"]

# Training settings that differ from mlp_predictor.trainer.default_dual_group_config,
# copied from the train_cfg cell of lvmsky's training notebook
# (notebook_sky_interpolation_triplet_dual_encoder_group_mlp_split_zodi_module.ipynb,
# lvmsky commit 59687bb, 2026-10-07) -- the configuration the lvmsky team
# trains their deployed model with.  Every other setting is the library default.
TRAIN_CFG_OVERRIDES = {
    "n_epochs": 400,
    "weight_decay": 2e-4,
    "zodi_head_extra_dims": (),
    "moon_group_weight": 3.0,
    "flux_mse_groups": ("moon", "zodi", "continuum"),
    "flux_amp_lambda": {"moon": 5.0, "zodi": 0.0},
    "blend_init_alpha": 0.85,
    "alpha_lr_mult": 30.0,
    "ensemble_workers": 4,
    "zodi_ctx_restriction": (
        "airmass", "vanrhijn_285km", "ecl_beta_deg", "ecl_lon_sin", "ecl_lon_cos",
        "zodi_log10_v", "sun_sep", "sun_alt", "alt", "moon_alt", "moon_sep",
        "moon_phase_sin", "moon_phase_cos", "moon_fli", "moon_up_smooth",
        "moon_airmass_up", "moon_signal_proxy", "zodi_po_log10", "moon_frac_po",
    ),
    "moon_zodi_ctx_restriction": (
        "airmass", "vanrhijn_285km", "ecl_beta_deg", "ecl_lon_sin", "ecl_lon_cos",
        "zodi_log10_v", "sun_sep", "moon_alt", "moon_sep", "moon_phase_sin",
        "moon_phase_cos", "moon_fli", "moon_up_smooth", "moon_airmass_up",
        "moon_signal_proxy", "moon_fli_x_phase_cos", "moon_sig_x_lon_cos",
        "moon_sig_x_lon_sin", "zodi_po_log10", "moon_frac_po",
    ),
}


def stage_convert(xcf_fits, decomp_input, force=False):
    if decomp_input.exists() and not force:
        print(f"[convert] exists, skip: {decomp_input}")
        return decomp_input
    print(f"[convert] {xcf_fits} -> {decomp_input}")
    ConvertForDecompose.convert(str(xcf_fits), str(decomp_input))
    return decomp_input


def stage_decompose(decomp_input, work_dir, lvmsky_skysub, n_workers, force=False):
    stem = decomp_input.stem
    outputs = [
        work_dir / f"{stem}_meta_only.fits",
        work_dir / f"{stem}_sky1_meta_coef{DECOMP_SUFFIX}.fits",
        work_dir / f"{stem}_sky2_meta_coef{DECOMP_SUFFIX}.fits",
        work_dir / f"{stem}_sci_meta_coef{DECOMP_SUFFIX}.fits",
        work_dir / f"{stem}_every10.fits",
    ] + [work_dir / f"{stem}_decomp_{arm}{DECOMP_SUFFIX}.fits"
         for arm in ("sky1", "sky2", "sci")]
    if all(p.exists() for p in outputs) and not force:
        print(f"[decompose] outputs exist, skip ({len(outputs)} files)")
        return

    script = Path(lvmsky_skysub) / "decompose_parallel.py"
    cmd = [
        sys.executable, str(script), str(decomp_input),
        "--fit-model", FIT_MODEL,
        "--n-workers", str(n_workers),
        "--output-dir", str(work_dir),
    ]
    print("[decompose] running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    missing = [str(p) for p in outputs if not p.exists()]
    if missing:
        raise RuntimeError(
            f"decompose_parallel.py did not produce expected outputs: {missing}"
        )


def _pipeline_config(work_dir, decomp_stem):
    cfg = PipelineConfig()
    cfg.data.decomp_data_root = str(work_dir)
    cfg.data.decomp_stem = decomp_stem
    cfg.data.decomp_suffix = DECOMP_SUFFIX
    return cfg


def stage_wavecache(work_dir, decomp_stem, atom_k_max, n_workers, force=False):
    """Build, filter and augment the training triplet; resolve coefficient
    wavelengths.  Follows cells 3 and 6 of lvmsky's training notebook."""
    cfg = _pipeline_config(work_dir, decomp_stem)
    pkl_path = work_dir / "filtered_triplet.pkl"

    if pkl_path.exists() and not force:
        print(f"[wavecache] exists, skip: {pkl_path}")
        with open(pkl_path, "rb") as fh:
            filtered = pickle.load(fh)
        return filtered, cfg

    trip = mp_data.build_triplet_coef_dataset(
        input_fits_path=cfg.data.input_fits_meta,
        sky_near_decomp_fits_path=cfg.data.coef_fits("sky1"),
        sky_far_decomp_fits_path=cfg.data.coef_fits("sky2"),
        sci_decomp_fits_path=cfg.data.coef_fits("sci"),
        context_columns=list(cfg.data.context_columns),
        return_chi2=True,
    )
    with fits.open(cfg.data.input_fits_for_basis) as hdul:
        wave = np.asarray(hdul["WAVE"].data, dtype=float)
        wave = wave if wave.ndim == 1 else wave[0]

    filtered = mp_data.apply_triplet_filters(
        trip,
        thin_every_n=1, chi2_qmax=90.0, chi2_min=0.0, chi2_max=10.0,
        hard_coef_bounds={"feo": (0.0, 36.01), "atom_k": (0.0, atom_k_max)},
        kappa=8.0, kappa_iter=3, oh_kappa=6.0, oh_kappa_iter=3,
        exclude_field_regions=[mp_data.LMC_EXCLUSION, mp_data.SMC_EXCLUSION],
        reversal_decomp_fits={
            "near": cfg.data.decomp_fits("sky1"),
            "far": cfg.data.decomp_fits("sky2"),
            "sci": cfg.data.decomp_fits("sci"),
        },
        reversal_wave=wave,
        reversal_min_component_frac=0.05,
        reversal_min_separation=0.0,
        colour_excess_input_fits=cfg.data.input_fits_flux,
        colour_excess_max=mp_data.SCI_COLOUR_EXCESS_MAX,
        diffuse_zeroed_frac=mp_data.DIFFUSE_ZEROED_FRAC,
    )
    print("[wavecache] filtered n_rows:", filtered["coef_sci"].shape[0])

    mp_data._augment_triplet_with_ecliptic(
        filtered, force=True, meta_fits_path=cfg.data.input_fits_meta)
    mp_data._augment_triplet_with_physics_priors(filtered, force=True)
    # Physical moon-model context feature; builds a cache beside the corpus
    # the first time (roughly 0.55 s per row per worker).
    mp_data._augment_triplet_with_moon_model(
        filtered, cfg.data.decomp_prefix, force=True, n_workers=n_workers)
    print(f"[wavecache] n_ctx={len(filtered['ctx_names'])}")

    ext = mp_wave.resolve_wavelengths_and_extinction(
        filtered,
        input_fits_for_basis=cfg.data.input_fits_for_basis,
        use_fitted_extinction=True,
        verbose=True,
        decomp_suffix=cfg.data.decomp_suffix,
    )
    filtered["wavecache_group_indices"] = ext.group_indices

    with open(pkl_path, "wb") as fh:
        pickle.dump(filtered, fh)
    print(f"[wavecache] wrote {pkl_path}")
    return filtered, cfg


def stage_train(filtered, cfg, work_dir, output_path, *, epochs, seeds,
                 flux_mse_groups, train_frac, val_frac, split_seed, n_bins,
                 force=False):
    """Fit the compressors and train the ensemble.  Follows cells 8-11 of
    lvmsky's training notebook, with TRAIN_CFG_OVERRIDES."""
    if output_path.exists() and not force:
        print(f"[train] exists, skip: {output_path}")
        return output_path

    group_indices = filtered["wavecache_group_indices"]
    print("[train] groups:", {g: len(idx) for g, idx in group_indices.items()})

    moon_phase = mp_ml.moon_phase_deg_from_ctx(filtered)
    train_idx, val_idx, test_idx = mp_ml.split_indices_by_moon_phase(
        filtered["obstime_mjd"], moon_phase,
        train_frac=train_frac, val_frac=val_frac, seed=split_seed, n_bins=n_bins,
    )
    print(f"[train] split: train={train_idx.size} val={val_idx.size} test={test_idx.size}")

    compressors, geom_kwargs = mp_comp.fit_all_group_compressors(
        filtered, group_indices, train_idx=train_idx, held_idx=val_idx,
        xarm_threshold=mp_comp.COMPRESSION_XARM_THRESHOLD,
    )
    filtered["compress_train_idx"] = train_idx
    filtered["compress_val_idx"] = val_idx
    filtered["compress_test_idx"] = test_idx

    n_moon_knots, split_zodi, n_zodi_knots = mp_wave.infer_spline_knots(filtered["coef_names"])
    print(f"[train] inferred spline knots: n_moon_knots={n_moon_knots} "
          f"split_zodi={split_zodi} n_zodi_knots={n_zodi_knots}")

    train_cfg = copy.deepcopy(mp_trainer.default_dual_group_config)
    train_cfg.update(copy.deepcopy(TRAIN_CFG_OVERRIDES))
    if epochs is not None:
        train_cfg["n_epochs"] = epochs
    if seeds is not None:
        train_cfg["ensemble_seeds"] = tuple(seeds)
    if flux_mse_groups is not None:
        train_cfg["flux_mse_groups"] = tuple(flux_mse_groups)
    print(f"[train] n_epochs={train_cfg['n_epochs']} "
          f"seeds={list(train_cfg['ensemble_seeds'])} "
          f"flux_mse_groups={train_cfg['flux_mse_groups']}")

    trainer = mp_trainer.Trainer(cfg=train_cfg)
    artifacts = trainer.run_ensemble(
        filtered, compressors, group_indices, geom_kwargs,
        input_fits_for_basis=cfg.data.input_fits_for_basis,
        input_fits_flux=cfg.data.input_fits_flux,
        n_moon_knots=n_moon_knots, split_zodi=split_zodi, n_zodi_knots=n_zodi_knots,
        decomp_suffix=cfg.data.decomp_suffix,
    )

    mp_ser.save_ensemble(artifacts.mlp_artifacts, output_path)
    print(f"[train] saved ensemble to {output_path}")

    aug_path = work_dir / "filtered_triplet_augmented.pkl"
    with open(aug_path, "wb") as fh:
        pickle.dump(filtered, fh)
    print(f"[train] saved augmented triplet to {aug_path}")
    return output_path


def main():
    p = argparse.ArgumentParser(
        parents=[_pre],
        description="Consolidated driver: SelectXCF.py output -> trained "
                    "mlp_ensemble_split_zodi ensemble.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("xcf_fits", help="XCframe-layout FITS file (SelectXCF.py output)")
    p.add_argument("-work_dir", default=None,
                   help="Directory for all pipeline artifacts (default: <xcf_fits stem>_train/)")
    p.add_argument("-start_at", choices=STAGES, default="convert", help="First stage to (re)run")
    p.add_argument("-stop_after", choices=STAGES, default="train", help="Last stage to run")
    p.add_argument("-force", action="store_true",
                   help="Redo stages even if their outputs already exist")
    p.add_argument("-np", dest="n_workers", type=int, default=4,
                   help="decompose_parallel.py worker processes")
    p.add_argument("-atom_k_max", type=float, default=10.01,
                   help="Upper hard bound for the atom_k coefficient filter")
    p.add_argument("-epochs", type=int, default=None,
                   help="Override n_epochs (default: 400, from TRAIN_CFG_OVERRIDES)")
    p.add_argument("-seeds", default=None,
                   help="Comma-separated ensemble seeds (default: the config's 10-seed ensemble)")
    p.add_argument("-flux_mse_groups", default=None,
                   help="Comma-separated groups for the flux-MSE loss term "
                        "(default: moon,zodi,continuum, from TRAIN_CFG_OVERRIDES)")
    p.add_argument("-train_frac", type=float, default=0.8, help="Training split fraction")
    p.add_argument("-val_frac", type=float, default=0.1, help="Validation split fraction")
    p.add_argument("-split_seed", type=int, default=42, help="Moon-phase split RNG seed")
    p.add_argument("-n_bins", type=int, default=10, help="Moon-phase stratification bins")
    p.add_argument("-output", default=None,
                   help="Trained ensemble .pt path (default: <work_dir>/mlp_ensemble.pt)")
    args = p.parse_args()

    if STAGES.index(args.start_at) > STAGES.index(args.stop_after):
        p.error(f"-start_at {args.start_at} comes after -stop_after {args.stop_after}")
    stages_to_run = STAGES[STAGES.index(args.start_at):STAGES.index(args.stop_after) + 1]

    xcf_path = Path(args.xcf_fits)
    stem = xcf_path.stem
    work_dir = Path(args.work_dir) if args.work_dir else Path(f"{stem}_train")
    work_dir.mkdir(parents=True, exist_ok=True)

    global _PLOT_DIR
    _PLOT_DIR = work_dir / "plots"

    decomp_stem = f"{stem}_decomp_input"
    decomp_input = work_dir / f"{decomp_stem}.fits"
    output_path = Path(args.output) if args.output else work_dir / "mlp_ensemble.pt"

    print(f"Work dir: {work_dir}")
    print(f"Stages to run: {stages_to_run}")

    if "convert" in stages_to_run:
        decomp_input = stage_convert(xcf_path, decomp_input, force=args.force)

    if "decompose" in stages_to_run:
        stage_decompose(decomp_input, work_dir, _pre_args.lvmsky_skysub,
                        args.n_workers, force=args.force)

    filtered, cfg = None, None
    if "wavecache" in stages_to_run:
        filtered, cfg = stage_wavecache(work_dir, decomp_stem, args.atom_k_max,
                                        args.n_workers, force=args.force)

    if "train" in stages_to_run:
        if filtered is None:
            cfg = _pipeline_config(work_dir, decomp_stem)
            with open(work_dir / "filtered_triplet.pkl", "rb") as fh:
                filtered = pickle.load(fh)
        seeds = [int(s) for s in args.seeds.split(",")] if args.seeds else None
        flux_groups = (None if args.flux_mse_groups is None
                       else [g for g in args.flux_mse_groups.split(",") if g])
        stage_train(
            filtered, cfg, work_dir, output_path,
            epochs=args.epochs, seeds=seeds, flux_mse_groups=flux_groups,
            train_frac=args.train_frac, val_frac=args.val_frac,
            split_seed=args.split_seed, n_bins=args.n_bins, force=args.force,
        )

    print("Done.")


if __name__ == "__main__":
    main()
