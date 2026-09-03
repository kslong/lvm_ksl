#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Run EsoSkyObs.py's run_sky_obs() across many rows of the same
    exposure corpus used by BatchPredictSky.py (the MLP candidate),
    producing a directly comparable WAVE/FLUX_OBS/FLUX_PRED/LINE_PRED/
    META batch file for the same evaluation tools (EvalFluxResiduals.py,
    MasterResidualByMoon.py) -- this is the ESO-sky-model candidate, not
    a replacement for the MLP one.

Command line usage (if any)::

    usage: BatchPredictSkyESO.py [-h] [--n N] [--rows R [R ...]]
                                 [--n-workers N] [--outfile PATH]
                                 [--site {lco,paranal}] [--engine E]
                                 [--lvm-ksl-progs PATH]
                                 meta_file batch_file

    where

    meta_file       the *_meta_only.fits file (has sci_ra/sci_dec/obstime
                    per row -- same row order as batch_file's META).
    batch_file      BatchPredictSky.py's own output (supplies WAVE and
                    the ground-truth FLUX_OBS/row/expnum to reuse
                    unchanged, so the two candidates are compared on
                    the exact same exposures and the exact same
                    observed spectra).

    --n N           use the first N rows of batch_file (default: 20).
    --rows R [R ...]
                    explicit row indices (overrides --n; indices are
                    into batch_file's FLUX_OBS/META, not meta_file).
    --n-workers N   parallel worker processes (default: 8).
    --outfile PATH  output FITS path (default:
                    <batch_file stem>_eso.fits).
    --site S        'lco' (default) or 'paranal' -- passed to
                    EsoSkyObs.run_sky_obs.
    --engine E      'auto' (default), 'local', or 'remote'.
    --lvm-ksl-progs PATH
                    path to the lvm_ksl repo's py_progs/ directory, which
                    supplies EsoSkyObs.py (default: ~/SDSS/lvm_ksl/py_progs).

Description:

    For each row, looks up sci_ra/sci_dec/obstime from meta_file (joined
    to batch_file via expnum) and calls EsoSkyObs.run_sky_obs(ra, dec,
    obstime, engine=..., site=...), which writes a small per-row FITS
    with WAVE (already on the same 3600-9800 A / 0.5 A grid as our
    XCframe corpus -- verified to sub-pixel precision) and FLUX/CONT/
    LINES columns (already extinction-corrected to the same
    above-the-atmosphere convention LVM spectra are compared against).

    FLUX_PRED is the FLUX column as-is. LINE_PRED is FLUX - CONT
    (verified numerically equal to LINES/trans, the model's own
    line-only prediction correctly extinction-corrected) -- ESO's own
    physically-based line/continuum split, not a mask or refit of ours.

    Each per-row scratch FITS is deleted after being read (EsoSkyObs.py
    already isolates every call in its own scratch tempdir internally
    for the calcskymodel inputs/outputs -- see
    project_esoskyobs_concurrency_fix.md -- so this script's own
    per-row output filenames only need to not collide with each other,
    which a row-indexed name guarantees).

History::

    260902  ksl  Coding begun.

'''

import argparse
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

THIS_DIR = str(Path(__file__).resolve().parent)
DEFAULT_LVM_KSL_PY_PROGS = str(Path('~/SDSS/lvm_ksl/py_progs').expanduser())

_WORKER_ESO = None
_WORKER_OUTDIR = None
_WORKER_SITE = None
_WORKER_ENGINE = None


def init_worker(outdir, site, engine, lvm_ksl_progs):
    global _WORKER_ESO, _WORKER_OUTDIR, _WORKER_SITE, _WORKER_ENGINE
    sys.path.insert(0, lvm_ksl_progs)
    import EsoSkyObs as eso  # noqa: E402
    _WORKER_ESO = eso
    _WORKER_OUTDIR = outdir
    _WORKER_SITE = site
    _WORKER_ENGINE = engine


def predict_one_row(task):
    row, expnum, sci_ra, sci_dec, obstime = task
    eso = _WORKER_ESO
    t0 = time.perf_counter()
    outroot = os.path.join(_WORKER_OUTDIR, f'eso_r{row:05d}')
    try:
        xroot = eso.run_sky_obs(
            float(sci_ra), float(sci_dec), obstime,
            outroot=outroot, engine=_WORKER_ENGINE, site=_WORKER_SITE,
        )
        if not xroot:
            return dict(row=row, ok=False, error='run_sky_obs returned empty root')
        with fits.open(xroot + '.fits') as hdul:
            d = hdul[1].data
            wave = np.asarray(d['WAVE'], dtype=np.float64)
            flux_pred = np.asarray(d['FLUX'], dtype=np.float64)
            line_pred = flux_pred - np.asarray(d['CONT'], dtype=np.float64)
    except Exception as exc:
        return dict(row=row, ok=False, error=f'{type(exc).__name__}: {exc}')
    finally:
        try:
            os.remove(outroot + '.fits')
        except OSError:
            pass
    dt = time.perf_counter() - t0
    return dict(
        row=row, ok=True, dt=dt, expnum=expnum,
        wave=wave, flux_pred=flux_pred, line_pred=line_pred,
    )


def main():
    p = argparse.ArgumentParser(
        description="Batch-run EsoSkyObs.py's run_sky_obs across the same corpus "
                    "used for BatchPredictSky.py's MLP output, for a directly "
                    "comparable candidate dataset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('meta_file', help='*_meta_only.fits with sci_ra/sci_dec/obstime/expnum')
    p.add_argument('batch_file', help="BatchPredictSky.py's output (WAVE/FLUX_OBS/row/expnum)")
    p.add_argument('--n', type=int, default=20, help='use the first N rows of batch_file')
    p.add_argument('--rows', type=int, nargs='+', default=None,
                   help='explicit row indices into batch_file (overrides --n)')
    p.add_argument('--n-workers', type=int, default=8, help='parallel worker processes')
    p.add_argument('--outfile', default=None, help='output FITS path')
    p.add_argument('--site', default='lco', choices=['lco', 'paranal'])
    p.add_argument('--engine', default='auto', choices=['auto', 'local', 'remote'])
    p.add_argument('--lvm-ksl-progs', default=DEFAULT_LVM_KSL_PY_PROGS,
                   help='path to the lvm_ksl repo\'s py_progs/ directory, which supplies EsoSkyObs.py')
    args = p.parse_args()

    with fits.open(args.meta_file) as hdul:
        meta = hdul[1].data
        meta_ra = {int(e): float(r) for e, r in zip(meta['expnum'], meta['sci_ra'])}
        meta_dec = {int(e): float(d) for e, d in zip(meta['expnum'], meta['sci_dec'])}
        meta_obstime = {int(e): str(o) for e, o in zip(meta['expnum'], meta['obstime'])}

    with fits.open(args.batch_file) as hdul:
        wave_ref = np.asarray(hdul['WAVE'].data, dtype=np.float64)
        flux_obs_all = np.asarray(hdul['FLUX_OBS'].data, dtype=np.float32)
        batch_row = np.asarray(hdul['META'].data['row'], dtype=np.int64)
        batch_expnum = np.asarray(hdul['META'].data['expnum'], dtype=np.int64)

    rows = args.rows if args.rows is not None else list(range(min(args.n, len(batch_row))))

    tasks = []
    for idx in rows:
        expnum = int(batch_expnum[idx])
        if expnum not in meta_ra:
            print(f'  batch row {idx} (expnum {expnum}): not found in {args.meta_file}, skipping')
            continue
        tasks.append((idx, expnum, meta_ra[expnum], meta_dec[expnum], meta_obstime[expnum]))

    outdir = str(Path(args.batch_file).resolve().parent / '_eso_scratch')
    os.makedirs(outdir, exist_ok=True)

    print(f'Predicting {len(tasks)} rows with {args.n_workers} workers '
          f'(engine={args.engine}, site={args.site}) ...')
    t0 = time.perf_counter()
    ctx = mp.get_context('spawn')
    with ctx.Pool(
        args.n_workers, initializer=init_worker,
        initargs=(outdir, args.site, args.engine, args.lvm_ksl_progs),
    ) as pool:
        results = pool.map(predict_one_row, tasks)
    total = time.perf_counter() - t0

    flux_pred_list, line_pred_list, flux_obs_list, expnum_list, row_list = [], [], [], [], []
    n_fail = 0
    for r in results:
        if not r['ok']:
            print(f'  row {r["row"]}: FAILED ({r["error"]}), skipping')
            n_fail += 1
            continue
        idx = r['row']
        # ESO's own grid matches wave_ref to sub-pixel precision (verified this
        # session) but interpolate anyway for exactness / robustness to future
        # grid changes on either side.
        flux_pred = np.interp(wave_ref, r['wave'], r['flux_pred'])
        line_pred = np.interp(wave_ref, r['wave'], r['line_pred'])
        flux_pred_list.append(flux_pred)
        line_pred_list.append(line_pred)
        flux_obs_list.append(flux_obs_all[idx])
        expnum_list.append(r['expnum'])
        row_list.append(idx)

    n_ok = len(row_list)
    print(f'\n{n_ok}/{len(tasks)} rows succeeded ({n_fail} failed) in {total:.1f}s '
          f'({total / max(n_ok, 1):.2f}s/row wall-clock, {args.n_workers} workers)')

    outfile = args.outfile
    if outfile is None:
        outfile = str(Path(args.batch_file).with_suffix('').as_posix()) + '_eso.fits'

    hdul = fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(data=np.asarray(wave_ref, dtype=np.float32), name='WAVE'),
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

    try:
        os.rmdir(outdir)
    except OSError:
        pass  # leave it if anything's still in there (failed cleanup on some row)


if __name__ == '__main__':
    main()
