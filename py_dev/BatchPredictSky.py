#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Run PredictSky.py's prediction (unmodified logic) across many rows of
    an XCframe file in parallel worker processes, collecting WAVE/
    FLUX_OBS/FLUX_PRED into batch arrays for a method-agnostic flux-space
    evaluator (e.g. sky_residual_eval.py) -- this script only produces the
    (observed, predicted) pairs for this one candidate (the MLP ensemble);
    it has no evaluation logic itself.

Command line usage (if any):

    usage: BatchPredictSky.py [-h] [--model PATH] [--n N] [--rows R [R ...]]
                              [--n-workers N] [--outfile PATH]
                              [--lvmsky-skysub PATH]
                              fits_file

    where

    fits_file       LVM XCframe summary FITS file.

    --model PATH    trained ensemble .pt archive (default: this
                    session's Stage 2 production checkpoint).

    --n N           use the first N rows of fits_file (default: 20).
                    Ignored if --rows is given.

    --rows R [R ...]
                    explicit row indices to use instead of --n.

    --n-workers N   parallel worker processes (default: 8).

    --outfile PATH  output FITS path for the batch WAVE/FLUX_OBS/
                    FLUX_PRED/EXPNUM arrays (default:
                    <stem>_batch_predictsky.fits).

Description:

    Each row involves two real QP decomposition fits (SkyE/SkyW) plus
    ephemeris/reconstruction work -- CPU-bound, independent across rows.
    A first version of this script (single process) measured only a ~19%
    speedup from reusing the ensemble/SkyDecomp instance across rows,
    because that process was already running near 800% CPU via NumPy/
    SciPy's automatic BLAS multithreading -- an inefficient kind of
    parallelism for matrices this size (high overhead relative to actual
    work), not something added deliberately.

    This version follows decompose_parallel.py's proven pattern instead:
    clamp every worker to 1 BLAS thread (_clamp_native_threads, same
    implementation) and get real parallelism from N independent
    single-threaded worker PROCESSES (multiprocessing, spawn context --
    fork would inherit the parent's already-initialised BLAS pools and
    undermine the clamp, same reasoning as decompose_parallel.py's own
    comment on this). Each worker loads the ensemble and builds its own
    SkyDecomp instance once (in init_worker), not once per row.

History::

    260902  ksl  Coding begun.  Rewritten same day to use
                 decompose_parallel.py's clamp-to-1-thread +
                 multiprocess-across-rows pattern after a single-process,
                 BLAS-auto-threaded version only gave a ~19% speedup.

'''

import argparse
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

DEFAULT_LVMSKY_SKYSUB = Path('~/SDSS/lvmsky/skysub').expanduser()
_pre = argparse.ArgumentParser(add_help=False)
_pre.add_argument('--lvmsky-skysub', default=str(DEFAULT_LVMSKY_SKYSUB))
_pre_args, _ = _pre.parse_known_args()

THIS_DIR = str(Path(__file__).resolve().parent)


def _clamp_native_threads(n=1):
    '''Force every loaded thread pool (BLAS/OpenMP/etc.) to `n` threads.

    Same implementation as decompose_parallel.py's own -- must run before
    numpy/scipy/torch are imported in this process to actually take
    effect (env-var based; a library that's already initialised its
    thread pool won't re-read these).
    '''
    for var in (
        'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
        'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS',
    ):
        os.environ[var] = str(n)


_clamp_native_threads(1)

# ---------------------------------------------------------------------------
# Worker-process globals (populated once per worker by init_worker, not
# once per row) -- multiprocessing.Pool with the spawn context gives each
# worker its own fresh process, so these are safely worker-local.
# ---------------------------------------------------------------------------
_WORKER_PS = None
_WORKER_ENSEMBLE = None
_WORKER_DECOMP = None
_WORKER_FITS_FILE = None


def init_worker(fits_file, model_path, lvmsky_skysub):
    global _WORKER_PS, _WORKER_ENSEMBLE, _WORKER_DECOMP, _WORKER_FITS_FILE
    _clamp_native_threads(1)

    sys.path.insert(0, THIS_DIR)
    sys.argv = [sys.argv[0], '--lvmsky-skysub', lvmsky_skysub]  # PredictSky's module-level _pre parses this
    import PredictSky as ps  # noqa: E402 -- deferred so the clamp above lands first
    _WORKER_PS = ps
    _WORKER_FITS_FILE = fits_file

    _WORKER_ENSEMBLE = ps.serialization.load_ensemble(str(model_path))
    with fits.open(fits_file, memmap=True) as hdul:
        wave0 = np.asarray(hdul['WAVE'].data, dtype=np.float64)
    _WORKER_DECOMP = ps.build_decomp(wave0, _WORKER_ENSEMBLE)


def predict_one_row(row):
    ps = _WORKER_PS
    t0 = time.perf_counter()
    try:
        row_data = ps.read_row(_WORKER_FITS_FILE, row=row)
        prediction = ps.predict_row(row_data, ensemble=_WORKER_ENSEMBLE, decomp=_WORKER_DECOMP)
    except Exception as exc:
        return dict(row=row, ok=False, error=f'{type(exc).__name__}: {exc}')
    dt = time.perf_counter() - t0
    return dict(
        row=row, ok=True, dt=dt,
        expnum=row_data['expnum'], confidence=prediction['confidence'],
        wave=prediction['wave'], flux_obs=row_data['flux_sci'],
        flux_pred=prediction['flux_pred'], line_pred=prediction['line_pred'],
    )


def main():
    p = argparse.ArgumentParser(
        parents=[_pre],
        description='Batch-run PredictSky.py predictions in parallel, collecting '
                    'WAVE/FLUX_OBS/FLUX_PRED arrays for a flux-space evaluator.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('fits_file', help='LVM XCframe summary FITS file')
    p.add_argument('--model',
                   default=str(Path(
                       '~/Projects/lvm_sky2609/niv/moon_zodi_stage2/'
                       'mlp_ensemble_stage2_production.pt').expanduser()),
                   help='trained ensemble .pt archive')
    p.add_argument('--n', type=int, default=20, help='use the first N rows')
    p.add_argument('--rows', type=int, nargs='+', default=None,
                   help='explicit row indices (overrides --n)')
    p.add_argument('--n-workers', type=int, default=8, help='parallel worker processes')
    p.add_argument('--outfile', default=None, help='output FITS path')
    args = p.parse_args()

    rows = args.rows if args.rows is not None else list(range(args.n))

    print(f'Predicting {len(rows)} rows with {args.n_workers} workers '
          f'(model={args.model}) ...')
    t0 = time.perf_counter()
    ctx = mp.get_context('spawn')
    with ctx.Pool(
        args.n_workers, initializer=init_worker,
        initargs=(args.fits_file, args.model, args.lvmsky_skysub),
    ) as pool:
        results = pool.map(predict_one_row, rows)
    total = time.perf_counter() - t0

    wave = None
    flux_obs_list, flux_pred_list, line_pred_list, expnum_list, row_list = [], [], [], [], []
    n_fail = 0
    for r in results:
        if not r['ok']:
            print(f'  row {r["row"]}: FAILED ({r["error"]}), skipping')
            n_fail += 1
            continue
        if wave is None:
            wave = r['wave']
        flux_obs_list.append(r['flux_obs'])
        flux_pred_list.append(r['flux_pred'])
        line_pred_list.append(r['line_pred'])
        expnum_list.append(r['expnum'])
        row_list.append(r['row'])

    n_ok = len(row_list)
    print(f'\n{n_ok}/{len(rows)} rows succeeded ({n_fail} failed) in {total:.1f}s '
          f'({total / max(n_ok, 1):.2f}s/row wall-clock, {args.n_workers} workers)')

    outfile = args.outfile
    if outfile is None:
        outfile = str(Path(args.fits_file).with_suffix('').as_posix()) + '_batch_predictsky.fits'

    hdul = fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(data=np.asarray(wave, dtype=np.float32), name='WAVE'),
        fits.ImageHDU(data=np.asarray(flux_obs_list, dtype=np.float32), name='FLUX_OBS'),
        fits.ImageHDU(data=np.asarray(flux_pred_list, dtype=np.float32), name='FLUX_PRED'),
        fits.ImageHDU(data=np.asarray(line_pred_list, dtype=np.float32), name='LINE_PRED'),
        fits.BinTableHDU.from_columns([
            fits.Column(name='row', format='K', array=np.asarray(row_list, dtype=np.int64)),
            fits.Column(name='expnum', format='K', array=np.asarray(expnum_list, dtype=np.int64)),
        ], name='META'),
    ])
    hdul.writeto(outfile, overwrite=True)
    print(f'Wrote {n_ok} rows to {outfile}')


if __name__ == '__main__':
    main()
