#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Consolidated, resumable driver that turns a SelectXCF.py-selected
    XCframe corpus into a trained mlp_ensemble_split_zodi ensemble
    checkpoint, running every intermediate step in one call instead of
    hand-invoking ConvertForDecompose.py / decompose_parallel.py / a
    wavecache script / a training script in sequence.

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
                    the corpus of exposures to train on.

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
                    (default: 10.0, widened from that function's own
                    default of 1.0 -- see Notes).

    -epochs EPOCHS  override mlp_predictor.trainer.default_dual_group_
                    config's n_epochs (default: unset, uses that config's
                    own default of 50; pass a small value for a quick
                    smoke-test run).

    -seeds SEEDS    comma-separated ensemble seeds, e.g. "42,43" (default:
                    unset, uses the config's own 10-seed production
                    ensemble).

    -flux_mse_groups GROUPS
                    comma-separated group names to add to the flux-MSE
                    loss term, e.g. "moon,zodi" (default: empty/off -- see
                    Notes; this pipeline does not build the per-row LSF
                    path that loss term needs).

    -train_frac F   moon-phase-stratified split fraction for training
                    (default: 0.7).

    -val_frac F     moon-phase-stratified split fraction for validation;
                    the remainder is held out as test (default: 0.15).

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
    its expected output files already exist, unless -force is given):

      1. convert    -- ConvertForDecompose.convert(): reformats xcf_fits
                       into the WAVE/FLUX_SCI/FLUX_SKY_NEAR/FLUX_SKY_FAR/
                       META layout decompose_parallel.py expects.
      2. decompose  -- runs lvmsky's decompose_parallel.py as a subprocess
                       (--fit-model lsf-surface-iterative-split-zodi,
                       -np workers), producing the meta_only,
                       {sky1,sky2,sci}_meta_coef, and every10 (LSF/flux
                       basis) FITS files.
      3. wavecache  -- mlp_predictor.data.build_triplet_coef_dataset +
                       apply_triplet_filters, then mlp_predictor.
                       wavelengths.resolve_wavelengths_and_extinction;
                       pickles the result as work_dir/filtered_triplet.pkl.
      4. train      -- augments the triplet with ecliptic/physics-prior
                       context, builds a moon-phase-stratified train/val/
                       test split, fits per-group compressors
                       (mlp_predictor.compressor.fit_all_group_
                       compressors), then trains the seed ensemble
                       (mlp_predictor.trainer.Trainer.run_ensemble) and
                       saves it with mlp_predictor.serialization.
                       save_ensemble. Also writes work_dir/
                       filtered_triplet_augmented.pkl (the same corpus
                       EvalCoefResiduals.py expects as its triplet_pkl
                       argument).

    The B-spline knot counts (n_moon_knots, split_zodi, n_zodi_knots) that
    run_ensemble needs are inferred from the decomposed coef_names via
    mlp_predictor.wavelengths.infer_spline_knots rather than hand-
    transcribed, removing a class of copy/paste error present in the
    scripts this replaces.

Notes::

    - -atom_k_max's default of 10.0 (vs. apply_triplet_filters' own 1.0)
      carries forward a widened bound adopted during the 260901-260902
      retraining session: the default hard bound rejected ~95-97% of
      rows because ATOM_K ran ~2x high vs. that bound for reasons not
      understood at the time (ruled out by-fiber-vs-by-pixel XCframe
      construction and multi- vs. single-exposure stacking as causes).
      Root cause is still open; this is a working default, not a fix.
    - flux_mse_groups defaults to off for the same reason both prior
      training runs (stage1_train.py/stage2_train.py) left it off: it
      needs a real per-row LSF path into decompose_parallel.py's input,
      which this pipeline does not build (ConvertForDecompose.py
      deliberately drops LSF -- lsf-surface-iterative-split-zodi's
      worker uses a single global --lsf-sigma scalar instead).
    - decompose is the only stage run as a subprocess rather than an
      in-process call: it is lvmsky's own long-running, CPU-thread-
      pinned, ProcessPoolExecutor-parallel script, and keeping it as a
      separate process preserves its own progress bar, thread-limiting
      env vars, and -np tuning independent of this driver.
    - No evaluation stage is included by design -- run EvalCoefResiduals.py
      (triplet_pkl=work_dir/filtered_triplet_augmented.pkl, meta_fits=
      work_dir/<decomp_stem>_meta_only.fits) and/or BatchPredictSky.py +
      EvalFluxResiduals.py against the saved checkpoint afterward.
    - mlp_predictor.data.apply_triplet_filters() (inside the wavecache
      stage) builds two Plotly diagnostic histograms and calls fig.show()
      on each; left alone this pops a browser tab per run under Plotly's
      default non-notebook renderer. Patched at import time here (not in
      lvmsky) by overriding plotly.basedatatypes.BaseFigure.show, so both
      figures are written to work_dir/plots/*.html instead, with no
      browser auto-open -- view them locally at will.

History::

    260904  ksl  Coding begun -- consolidates ConvertForDecompose.py,
        lvmsky's decompose_parallel.py, and the wavecache/train logic
        previously duplicated between niv/stage1_wavecache.py and
        niv/stage1_train.py + niv/stage2_train.py into one resumable
        driver.

'''

import argparse
import copy
import os
import pickle
import subprocess
import sys
from pathlib import Path

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

DECOMP_SUFFIX = "_lsf_surface_iterative_split_zodi"


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
    ]
    if all(p.exists() for p in outputs) and not force:
        print(f"[decompose] outputs exist, skip ({len(outputs)} files)")
        return

    script = Path(lvmsky_skysub) / "decompose_parallel.py"
    cmd = [
        sys.executable, str(script), str(decomp_input),
        "--fit-model", "lsf-surface-iterative-split-zodi",
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
    return cfg


def stage_wavecache(work_dir, decomp_stem, atom_k_max, force=False):
    cfg = _pipeline_config(work_dir, decomp_stem)
    pkl_path = work_dir / "filtered_triplet.pkl"

    if pkl_path.exists() and not force:
        print(f"[wavecache] exists, skip: {pkl_path}")
        with open(pkl_path, "rb") as fh:
            filtered = pickle.load(fh)
        return filtered, cfg

    trip = mp_data.build_triplet_coef_dataset(
        cfg.data.input_fits_meta,
        cfg.data.coef_fits("sky1"),
        cfg.data.coef_fits("sky2"),
        cfg.data.coef_fits("sci"),
        context_columns=cfg.data.context_columns,
    )
    filtered = mp_data.apply_triplet_filters(
        trip, hard_coef_bounds={"feo": (0.0, 1.0), "atom_k": (0.0, atom_k_max)},
    )
    print("[wavecache] filtered n_rows:", filtered["coef_sci"].shape[0])

    mp_wave.resolve_wavelengths_and_extinction(
        filtered,
        cache_path=cfg.data.wavelength_cache_path,
        input_fits_for_basis=cfg.data.input_fits_for_basis,
    )

    with open(pkl_path, "wb") as fh:
        pickle.dump(filtered, fh)
    print(f"[wavecache] wrote {pkl_path}")
    return filtered, cfg


def stage_train(filtered, cfg, work_dir, output_path, *, epochs, seeds,
                 flux_mse_groups, train_frac, val_frac, split_seed, n_bins,
                 force=False):
    if output_path.exists() and not force:
        print(f"[train] exists, skip: {output_path}")
        return output_path

    mp_data._augment_triplet_with_ecliptic(filtered, meta_fits_path=cfg.data.input_fits_meta)
    mp_data._augment_triplet_with_physics_priors(filtered)

    group_indices = mp_data._build_group_indices(filtered["coef_names"])
    print("[train] groups:", {g: len(idx) for g, idx in group_indices.items()})

    moon_phase = mp_ml.moon_phase_deg_from_ctx(filtered, arm="sci")
    train_idx, val_idx, test_idx = mp_ml.split_indices_by_moon_phase(
        filtered["obstime_mjd"], moon_phase,
        train_frac=train_frac, val_frac=val_frac, seed=split_seed, n_bins=n_bins,
    )
    print(f"[train] split: train={train_idx.size} val={val_idx.size} test={test_idx.size}")
    filtered["compress_train_idx"] = train_idx
    filtered["compress_val_idx"] = val_idx
    filtered["compress_test_idx"] = test_idx

    compressors, geom_kwargs = mp_comp.fit_all_group_compressors(
        filtered, group_indices, train_idx=train_idx, held_idx=val_idx,
    )

    n_moon_knots, split_zodi, n_zodi_knots = mp_wave.infer_spline_knots(filtered["coef_names"])
    print(f"[train] inferred spline knots: n_moon_knots={n_moon_knots} "
          f"split_zodi={split_zodi} n_zodi_knots={n_zodi_knots}")

    train_cfg = copy.deepcopy(mp_trainer.default_dual_group_config)
    if epochs is not None:
        train_cfg["n_epochs"] = epochs
    if seeds is not None:
        train_cfg["ensemble_seeds"] = tuple(seeds)
    train_cfg["flux_mse_groups"] = tuple(flux_mse_groups)

    trainer = mp_trainer.Trainer(cfg=train_cfg)
    artifacts = trainer.run_ensemble(
        filtered, compressors, group_indices, geom_kwargs,
        input_fits_for_basis=cfg.data.input_fits_for_basis,
        n_moon_knots=n_moon_knots, split_zodi=split_zodi, n_zodi_knots=n_zodi_knots,
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
    p.add_argument("-atom_k_max", type=float, default=10.0,
                   help="Upper hard bound for the atom_k coefficient filter")
    p.add_argument("-epochs", type=int, default=None,
                   help="Override n_epochs (default: default_dual_group_config's own 50)")
    p.add_argument("-seeds", default=None,
                   help="Comma-separated ensemble seeds (default: the config's 10-seed ensemble)")
    p.add_argument("-flux_mse_groups", default="",
                   help="Comma-separated groups for the flux-MSE loss term (default: off)")
    p.add_argument("-train_frac", type=float, default=0.7, help="Training split fraction")
    p.add_argument("-val_frac", type=float, default=0.15, help="Validation split fraction")
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
        filtered, cfg = stage_wavecache(work_dir, decomp_stem, args.atom_k_max, force=args.force)

    if "train" in stages_to_run:
        if filtered is None:
            cfg = _pipeline_config(work_dir, decomp_stem)
            with open(work_dir / "filtered_triplet.pkl", "rb") as fh:
                filtered = pickle.load(fh)
        seeds = [int(s) for s in args.seeds.split(",")] if args.seeds else None
        flux_groups = [g for g in args.flux_mse_groups.split(",") if g]
        stage_train(
            filtered, cfg, work_dir, output_path,
            epochs=args.epochs, seeds=seeds, flux_mse_groups=flux_groups,
            train_frac=args.train_frac, val_frac=args.val_frac,
            split_seed=args.split_seed, n_bins=args.n_bins, force=args.force,
        )

    print("Done.")


if __name__ == "__main__":
    main()
