#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Layer-1 model-quality evaluation: per-exposure, per-coefficient-group
    weighted RMSE between the trained mlp_ensemble_split_zodi ensemble's
    predicted Sci-arm coefficients and the real (decomposed) ones, across
    an entire corpus at once -- no re-decomposition or flux reconstruction
    needed, since the corpus already carries real coef_near/coef_far/
    coef_sci/context for every row.  Scales to any corpus size that has
    been through decompose_parallel.py; wall-clock is dominated by
    loading the ensemble, not by the corpus size.

Command line usage (if any):

    usage: EvalCoefResiduals.py [-h] [-triplet_pkl PATH] [-meta_fits PATH]
                                [-model PATH] [-top N] [-outfile PATH]
                                [-lvmsky_skysub PATH]
                                work_dir

    where

    work_dir        a TrainSkyModel.py work_dir (e.g. the directory
                    holding its mlp_ensemble.pt). triplet_pkl, meta_fits
                    and model are all auto-derived from this single
                    directory -- see Notes -- since TrainSkyModel.py
                    always writes them together; pass the three -*
                    overrides below only to mix files from different
                    directories (e.g. evaluate a different checkpoint
                    against this work_dir's own triplet).

    -triplet_pkl PATH
                    pickled, augmented filtered-triplet dict (default:
                    work_dir/filtered_triplet_augmented.pkl -- must
                    already carry ctx_names/coef_names/coef_near/
                    coef_far/coef_sci/coef_err_sci/row_index/obstime_mjd
                    and compress_train_idx/compress_val_idx/
                    compress_test_idx).

    -meta_fits PATH the corpus's "*_meta_only.fits" file, for the expnum
                    column (default: the single such file found directly
                    in work_dir; everything else needed is already in
                    triplet_pkl).

    -model PATH     trained ensemble .pt archive (default:
                    work_dir/mlp_ensemble.pt, if that exists -- pass
                    explicitly to evaluate a checkpoint from elsewhere
                    against this work_dir's triplet).

    -top N          print the N worst exposures by overall WRMSE
                    (default: 20).

    -outfile PATH   output table path (default: work_dir/coef_residuals.fits).

    -lvmsky_skysub PATH
                    path to the lvmsky repo's skysub/ directory (default:
                    ~/SDSS/lvmsky/skysub).

Description:

    For every row in triplet_pkl (train + val + test together -- this is
    a corpus-quality audit, not a held-out-only metric):

      1. mlp_predictor.trainer.predict_sci_coefficients_default() predicts
         coef_sci for the whole corpus in one call (ensemble-averaged in
         physical space, same function used for the ensemble's own
         reported test metrics).
      2. mlp_predictor.metrics.weighted_rmse_per_row() gives one overall
         WRMSE per exposure, and is called again per coefficient group
         (moon/zodi/continuum/mesospheric/ionospheric/atomic) restricted
         to that group's own columns, using the same per-group sigma
         floors the trainer itself uses
         (DEFAULT_COEF_ERR_SIGMA_FLOOR_BY_GROUP) -- so a bad exposure can
         be attributed to a specific physical component, not just flagged
         as "bad".

    Output is one row per exposure: expnum, mjd, moon_alt, moon_fli,
    split (train/val/test), WRMSE_ALL, WRMSE_<group> for each group --
    sorted worst-first by WRMSE_ALL and written to -outfile for later
    filtering/thresholding.

Notes::

    work_dir's auto-derivation is a directory-local default, not the
    kind of silent global fallback -model deliberately avoids (see
    260904 History entry below): work_dir/mlp_ensemble.pt is *the
    checkpoint that specific run produced*, not an arbitrary "current"
    model that could go stale as the project moves on to new ones. The
    -triplet_pkl/-meta_fits/-model overrides exist precisely for the
    case that default doesn't cover -- evaluating a checkpoint from one
    run against a triplet from another.

History::

    260902  ksl  Coding begun.
    260904  ksl  Promoted from niv/ to py_dev/. Switched every option
        from double-dash (--model) to single-dash (-model), matching
        py_progs/'s convention. Removed the DEFAULT_MODEL placeholder
        (pointed at one specific niv/ checkpoint) -- -model was made a
        required argument with no default, matching PredictSky.py/
        BatchPredictSky.py's precedent: this tier works against
        different trained models for different purposes, so a silent
        fallback risks going stale.
    260904  ksl  Replaced the triplet_pkl/meta_fits positional pair (and
        required -model) with a single work_dir positional, since
        TrainSkyModel.py always writes all three into the same
        directory -- auto-derived, with -triplet_pkl/-meta_fits/-model
        left as overrides for mixing files across directories.

'''

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table

DEFAULT_LVMSKY_SKYSUB = Path('~/SDSS/lvmsky/skysub').expanduser()

_pre = argparse.ArgumentParser(add_help=False)
_pre.add_argument('-lvmsky_skysub', default=str(DEFAULT_LVMSKY_SKYSUB))
_pre_args, _ = _pre.parse_known_args()
sys.path.insert(0, _pre_args.lvmsky_skysub)

from mlp_predictor import serialization, data as mp_data, trainer as mp_trainer  # noqa: E402
from mlp_predictor.metrics import weighted_rmse_per_row  # noqa: E402


def resolve_inputs(work_dir, triplet_pkl=None, meta_fits=None, model=None):
    '''
    Fill in triplet_pkl/meta_fits/model from work_dir wherever an
    explicit override wasn't given.

    Returns
    -------
    (triplet_pkl, meta_fits, model) : tuple of str
    '''
    work_dir = Path(work_dir)

    if triplet_pkl is None:
        triplet_pkl = work_dir / 'filtered_triplet_augmented.pkl'
        if not triplet_pkl.exists():
            raise FileNotFoundError(
                f'{triplet_pkl} not found -- pass -triplet_pkl explicitly '
                f'if this work_dir uses a different layout.')

    if meta_fits is None:
        matches = sorted(work_dir.glob('*_meta_only.fits'))
        if len(matches) != 1:
            raise FileNotFoundError(
                f'Expected exactly one *_meta_only.fits in {work_dir}, '
                f'found {len(matches)} ({[m.name for m in matches]}) -- '
                f'pass -meta_fits explicitly.')
        meta_fits = matches[0]

    if model is None:
        model = work_dir / 'mlp_ensemble.pt'
        if not model.exists():
            raise FileNotFoundError(
                f'{model} not found -- pass -model explicitly to evaluate '
                f'a checkpoint from elsewhere against this work_dir.')

    return str(triplet_pkl), str(meta_fits), str(model)


def evaluate(triplet, meta_fits, ensemble):
    '''
    Compute per-exposure overall + per-group WRMSE for every row in
    triplet against ensemble's predicted coef_sci.

    Returns
    -------
    astropy.table.Table
    '''
    coef_near = np.asarray(triplet['coef_near'], dtype=np.float32)
    coef_far = np.asarray(triplet['coef_far'], dtype=np.float32)
    coef_sci = np.asarray(triplet['coef_sci'], dtype=np.float32)
    ctx_near = np.asarray(triplet['ctx_near'], dtype=np.float32)
    ctx_far = np.asarray(triplet['ctx_far'], dtype=np.float32)
    ctx_sci = np.asarray(triplet['ctx_sci'], dtype=np.float32)
    coef_names = list(triplet['coef_names'])
    ctx_names = list(triplet['ctx_names'])
    n_rows = coef_sci.shape[0]

    coef_err_sci = np.asarray(
        triplet.get('coef_err_sci', np.full_like(coef_sci, np.nan)), dtype=np.float32)

    print(f'Predicting coef_sci for all {n_rows} rows ...')
    coef_pred = mp_trainer.predict_sci_coefficients_default(
        ensemble,
        coef_near_phys=coef_near, coef_far_phys=coef_far,
        ctx_near_phys=ctx_near, ctx_far_phys=ctx_far, ctx_sci_phys=ctx_sci,
    ).astype(np.float32)

    group_indices = mp_data._build_group_indices(coef_names)
    floor_by_group = dict(mp_trainer.DEFAULT_COEF_ERR_SIGMA_FLOOR_BY_GROUP)

    wrmse_all = weighted_rmse_per_row(
        coef_sci, coef_pred, coef_err_sci, group_indices, floor_by_group)

    per_group_wrmse = {}
    for gname, gidx in group_indices.items():
        gidx = np.asarray(gidx, dtype=int)
        local_group_indices = {gname: np.arange(gidx.size)}
        per_group_wrmse[gname] = weighted_rmse_per_row(
            coef_sci[:, gidx], coef_pred[:, gidx], coef_err_sci[:, gidx],
            local_group_indices, floor_by_group)

    # Per-exposure metadata: expnum from meta_fits (indexed by row_index,
    # the original-corpus row each triplet row came from); mjd/moon_alt/
    # moon_fli already carried by triplet itself.
    row_index = np.asarray(triplet['row_index'], dtype=int)
    with fits.open(meta_fits) as hdul:
        expnum = np.asarray(hdul['META'].data['expnum'])[row_index]

    ctx_idx = {n: i for i, n in enumerate(ctx_names)}
    moon_alt = ctx_sci[:, ctx_idx['moon_alt']] if 'moon_alt' in ctx_idx else np.full(n_rows, np.nan)
    moon_fli = ctx_sci[:, ctx_idx['moon_fli']] if 'moon_fli' in ctx_idx else np.full(n_rows, np.nan)

    split = np.full(n_rows, 'unassigned', dtype='<U11')
    split[np.asarray(triplet['compress_train_idx'], dtype=int)] = 'train'
    split[np.asarray(triplet['compress_val_idx'], dtype=int)] = 'val'
    split[np.asarray(triplet['compress_test_idx'], dtype=int)] = 'test'

    mjd = np.asarray(triplet['obstime_mjd'], dtype=np.float64)
    night_id = np.floor(mjd - 0.5).astype(int)

    cols = dict(
        expnum=expnum,
        mjd=mjd, night_id=night_id,
        moon_alt=moon_alt, moon_fli=moon_fli,
        split=split,
        WRMSE_ALL=wrmse_all,
    )
    for gname in group_indices:
        cols[f'WRMSE_{gname}'] = per_group_wrmse[gname]

    # WRMSE_ALL is dominated by whichever group has the most coefficients
    # (mesospheric = 358/388) -- it's a proxy for OH quality, not a
    # balanced cross-process summary. Instead, normalise each group's
    # WRMSE to its own corpus-wide median (a dimensionless ratio, ~1.0 for
    # a typical row regardless of that group's absolute WRMSE scale), then
    # combine across groups two ways: the mean (flags exposures broadly
    # off across many processes) and the max (flags exposures with one
    # severely bad process even if the rest are fine).
    norm_cols = np.column_stack([
        per_group_wrmse[g] / max(float(np.nanmedian(per_group_wrmse[g])), 1e-30)
        for g in group_indices
    ])
    cols['MEAN_NORM_WRMSE'] = np.nanmean(norm_cols, axis=1)
    cols['MAX_NORM_WRMSE'] = np.nanmax(norm_cols, axis=1)
    worst_group_idx = np.nanargmax(norm_cols, axis=1)
    group_names = list(group_indices.keys())
    cols['WORST_GROUP'] = np.asarray([group_names[i] for i in worst_group_idx])

    tab = Table(cols)
    tab.sort('MEAN_NORM_WRMSE', reverse=True)
    return tab


def main():
    p = argparse.ArgumentParser(
        parents=[_pre],
        description='Per-exposure, per-group coefficient-space model-quality audit.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('work_dir', help='a TrainSkyModel.py work_dir')
    p.add_argument('-triplet_pkl', default=None,
                   help='override: pickled augmented filtered-triplet dict '
                        '(default: work_dir/filtered_triplet_augmented.pkl)')
    p.add_argument('-meta_fits', default=None,
                   help='override: corpus *_meta_only.fits file '
                        '(default: the one found in work_dir)')
    p.add_argument('-model', default=None,
                   help='override: trained ensemble .pt archive '
                        '(default: work_dir/mlp_ensemble.pt)')
    p.add_argument('-top', type=int, default=20, help='print the N worst exposures')
    p.add_argument('-outfile', default=None, help='output table path')
    args = p.parse_args()

    triplet_pkl, meta_fits, model = resolve_inputs(
        args.work_dir, args.triplet_pkl, args.meta_fits, args.model)

    with open(triplet_pkl, 'rb') as fh:
        triplet = pickle.load(fh)

    ensemble = serialization.load_ensemble(model)
    tab = evaluate(triplet, meta_fits, ensemble)

    print(f'\nWorst {args.top} exposures by MEAN_NORM_WRMSE '
          f'(per-group WRMSE normalised to that group\'s own corpus median, '
          f'averaged across groups -- not dominated by coefficient count):')
    tab[:args.top].pprint(max_width=200)

    outfile = args.outfile
    if outfile is None:
        outfile = str(Path(args.work_dir) / 'coef_residuals.fits')
    tab.write(outfile, overwrite=True)
    print(f'\nWrote {len(tab)} rows to {outfile}')


if __name__ == '__main__':
    main()
