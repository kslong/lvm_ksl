#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Run EsoSkyObs.py's run_sky_obs() across many rows of an XCframe-
    layout corpus (e.g. SelectXCF.py's output -- the same file
    BatchPredictSky.py itself takes directly), producing a WAVE/
    FLUX_OBS/FLUX_PRED/LINE_PRED/META batch file for the same
    evaluation tools (EvalFluxResiduals.py, PlotSkyResiduals.py) --
    this is the ESO-sky-model candidate, not a replacement for the MLP
    one. Applies the same per-row instrumental LSF convolution
    PredictSky.py applies to the MLP candidate (via
    sky_decomp.lsf_surface_iterative's build_lsf_operator), so the two
    candidates are compared on equal footing rather than one being
    convolved and the other not.

Command line usage (if any)::

    usage: BatchPredictSkyESO.py [-h] [-n N] [-rows R [R ...]]
                                 [-np N] [-outfile PATH]
                                 [-no_lsf_convolve]
                                 [-site {lco,paranal}] [-engine E]
                                 [-lvm_ksl_progs PATH]
                                 [-lvmsky_skysub PATH]
                                 fits_file

    where

    fits_file       an XCframe-layout FITS file (e.g. SelectXCF.py's
                    output) with WAVE/FLUX/LSF extensions and a DRP_ALL
                    table (sci_ra/sci_dec/obstime/expnum) -- the same
                    file BatchPredictSky.py itself takes directly, so
                    both candidates can be run independently against the
                    exact same test set with no ordering dependency
                    between them.

    -n N            use only the first N rows of fits_file (default: all
                    rows).
    -rows R [R ...]
                    explicit row indices (overrides -n).
    -np N           parallel worker processes (default: 8; matches
                    py_progs/Reduce.py and py_progs/sky_gaussfit.py's
                    process-count convention).
    -outfile PATH   output FITS path (default: <fits_file stem>_eso_lsf.fits
                    with the LSF convolution, <fits_file stem>_eso_nolsf.fits
                    with -no_lsf_convolve -- the suffix always records which
                    was used, so the two are never confusable on disk).
    -no_lsf_convolve
                    skip the per-row LSF convolution described below,
                    writing ESO's raw interpolated-onto-fits_file's-WAVE
                    prediction instead -- for seeing the effect of that
                    convolution directly against the default output for
                    the same fits_file.
    -site S         'lco' (default) or 'paranal' -- passed to
                    EsoSkyObs.run_sky_obs.
    -engine E       'auto' (default), 'local', or 'remote'.
    -lvm_ksl_progs PATH
                    path to the lvm_ksl repo's py_progs/ directory, which
                    supplies EsoSkyObs.py (default: ~/SDSS/lvm_ksl/py_progs).
    -lvmsky_skysub PATH
                    path to the lvmsky repo's skysub/ directory, which
                    supplies sky_decomp.lsf_surface_iterative's LSF
                    convolution operator (default: ~/SDSS/lvmsky/skysub).

Description:

    For each row, reads sci_ra/sci_dec/obstime/expnum straight from
    fits_file's DRP_ALL table and calls EsoSkyObs.run_sky_obs(ra, dec,
    obstime, engine=..., site=...), which writes a small per-row FITS
    with WAVE (already on the same 3600-9800 A / 0.5 A grid as our
    XCframe corpus -- verified to sub-pixel precision) and FLUX/CONT/
    LINES columns (already extinction-corrected to the same
    above-the-atmosphere convention LVM spectra are compared against).

    FLUX_PRED is the FLUX column as-is, LINE_PRED is FLUX - CONT
    (verified numerically equal to LINES/trans, the model's own
    line-only prediction correctly extinction-corrected) -- ESO's own
    physically-based line/continuum split, not a mask or refit of ours
    -- both interpolated onto fits_file's own WAVE grid and then
    convolved with a Gaussian kernel built from that row's own LSF
    extension (FWHM -> sigma_pix = fwhm/2.355/dlam_pix, evaluated at
    LSF_TAP_OFFSETS, each row normalized to sum to 1, applied via
    build_lsf_operator -- identical construction to PredictSky.py's own
    reconstruction step). Convolution is linear, so convolving FLUX_PRED
    and LINE_PRED independently with the same per-row kernel is
    equivalent to convolving FLUX_PRED and CONT independently and
    resumming. Before this change, the ESO candidate was compared
    unconvolved against LVM data while the MLP candidate already had
    this LSF reconstruction applied -- see docs/source/
    sky_model_landscape.rst's "Current Open Problem" section for the
    fuller context (the ESO local engine's own internal convolution is
    a fixed, LVM-untuned ~0.4 A kernel, unrelated to this per-row step).

    Each per-row scratch FITS is deleted after being read (EsoSkyObs.py
    already isolates every call in its own scratch tempdir internally
    for the calcskymodel inputs/outputs -- see
    project_esoskyobs_concurrency_fix.md -- so this script's own
    per-row output filenames only need to not collide with each other,
    which a row-indexed name guarantees).

History::

    260902  ksl  Coding begun.
    260903  ksl  Switched from a (meta_file, batch_file) pair to a single
        fits_file argument (the same XCframe-layout file
        BatchPredictSky.py itself already takes) -- meta_file's
        sci_ra/sci_dec/obstime/expnum are all in fits_file's own
        DRP_ALL table, and batch_file's WAVE/FLUX_OBS are identical to
        fits_file's WAVE/FLUX (FLUX_OBS is FLUX unchanged, verified in
        PredictSky.read_row). This also removes the previous ordering
        dependency on BatchPredictSky.py having already been run.
        Added the LSF convolution step described above, made possible
        by fits_file's own LSF extension being available now that
        meta_file (which never carried spectral extensions at all) is
        no longer in the loop.
    260903  ksl  Added -no_lsf_convolve to see the convolution's effect
        directly (skips it, writing ESO's raw interpolated prediction
        instead), and changed the default -outfile suffix from _eso to
        _eso_lsf/_eso_nolsf so a given fits_file's two possible outputs
        never collide and the filename itself records which was used.
    260903  ksl  Switched every option from double-dash (--n-workers) to
        single-dash (-n_workers), matching py_progs/'s convention --
        py_dev had drifted onto ordinary Python argparse habits instead
        (double-dash, hyphens) without reference to it; only
        EvalFluxResiduals.py had already followed py_progs's style.
    260904  ksl  -n_workers renamed to -np, matching py_progs/Reduce.py/
        sky_gaussfit.py's own process-count spelling more closely (this
        script's -n_workers had matched BatchPredictSky.py instead, a
        second, less-precise convention that had crept in alongside it).
        -n's default changed from 20 to all rows in fits_file.

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
DEFAULT_LVMSKY_SKYSUB = str(Path('~/SDSS/lvmsky/skysub').expanduser())

# calcskymodel's own internal LSF convolution (EsoSkyObs.py's
# lsf_gauss_fwhm_pix, local engine only) -- passed explicitly to
# run_sky_obs() below AND used in main()'s quadrature correction, so the
# two stay in sync even if EsoSkyObs.py's own default ever changes. Do
# not lower this to try to reduce ESO's internal smoothing further --
# see EsoSkyObs.create_local_inputs()'s docstring for the aliasing floor
# that makes 0.8 pixels close to as narrow as it should safely go.
ESO_LSF_GAUSS_FWHM_PIX = 0.8

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
            lsf_gauss_fwhm_pix=ESO_LSF_GAUSS_FWHM_PIX,
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


def _clean_lsf(fwhm):
    '''
    Fill any non-finite value in a per-row LSF FWHM array by linear
    interpolation from its own neighbours (same fallback PredictSky.py's
    read_row uses) so the convolution kernel below never sees a NaN.
    '''
    fwhm = np.asarray(fwhm, dtype=np.float64)
    bad = ~np.isfinite(fwhm)
    if bad.any():
        good = ~bad
        if good.any():
            fwhm = fwhm.copy()
            fwhm[bad] = np.interp(np.flatnonzero(bad), np.flatnonzero(good), fwhm[good])
    return fwhm


def main():
    p = argparse.ArgumentParser(
        description="Batch-run EsoSkyObs.py's run_sky_obs across an XCframe-layout "
                    "corpus, applying the same per-row LSF convolution PredictSky.py "
                    "applies to the MLP candidate, for a directly comparable dataset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('fits_file', help='XCframe-layout FITS file (WAVE/FLUX/LSF/DRP_ALL)')
    p.add_argument('-n', type=int, default=None,
                   help='use only the first N rows of fits_file (default: all rows)')
    p.add_argument('-rows', type=int, nargs='+', default=None,
                   help='explicit row indices into fits_file (overrides -n)')
    p.add_argument('-np', dest='n_workers', type=int, default=8,
                   help='parallel worker processes')
    p.add_argument('-outfile', default=None,
                   help='output FITS path (default: <fits_file stem>_eso_lsf.fits with the '
                        'LSF convolution, <fits_file stem>_eso_nolsf.fits with -no_lsf_convolve)')
    p.add_argument('-no_lsf_convolve', action='store_true',
                   help='skip the per-row LSF convolution (see module Description) and write '
                        'ESO\'s raw interpolated-onto-wave_ref prediction instead -- for seeing '
                        'the effect of that convolution directly, e.g. against the default output '
                        'for the same fits_file')
    p.add_argument('-site', default='lco', choices=['lco', 'paranal'])
    p.add_argument('-engine', default='auto', choices=['auto', 'local', 'remote'])
    p.add_argument('-lvm_ksl_progs', default=DEFAULT_LVM_KSL_PY_PROGS,
                   help='path to the lvm_ksl repo\'s py_progs/ directory, which supplies EsoSkyObs.py')
    p.add_argument('-lvmsky_skysub', default=DEFAULT_LVMSKY_SKYSUB,
                   help='path to the lvmsky repo\'s skysub/ directory, which supplies '
                        'sky_decomp.lsf_surface_iterative')
    args = p.parse_args()

    with fits.open(args.fits_file) as hdul:
        wave_ref = np.asarray(hdul['WAVE'].data, dtype=np.float64)
        flux_obs_all = np.asarray(hdul['FLUX'].data, dtype=np.float32)
        lsf_all = np.asarray(hdul['LSF'].data, dtype=np.float64)
        drp = hdul['DRP_ALL'].data
        sci_ra_all = np.asarray(drp['sci_ra'], dtype=np.float64)
        sci_dec_all = np.asarray(drp['sci_dec'], dtype=np.float64)
        obstime_all = np.asarray(drp['obstime'], dtype=str)
        expnum_all = np.asarray(drp['expnum'], dtype=np.int64)

    n_rows = len(sci_ra_all)
    if args.rows is not None:
        rows = args.rows
    elif args.n is not None:
        rows = list(range(min(args.n, n_rows)))
    else:
        rows = list(range(n_rows))

    tasks = [
        (idx, int(expnum_all[idx]), float(sci_ra_all[idx]), float(sci_dec_all[idx]),
         str(obstime_all[idx]))
        for idx in rows
    ]

    outdir = str(Path(args.fits_file).resolve().parent / '_eso_scratch')
    os.makedirs(outdir, exist_ok=True)

    sys.path.insert(0, args.lvmsky_skysub)
    from scipy.stats import norm  # noqa: E402
    from sky_decomp.lsf_surface_iterative import LSF_TAP_OFFSETS, build_lsf_operator  # noqa: E402

    dlam_pix = float(np.median(np.diff(wave_ref)))
    taps = np.asarray(LSF_TAP_OFFSETS, dtype=np.float64)
    eso_internal_fwhm_ang = ESO_LSF_GAUSS_FWHM_PIX * dlam_pix
    n_floored = [0]
    n_floored_total = [0]

    def _convolve(spec, lsf_fwhm):
        # calcskymodel already convolved this spectrum with a fixed
        # ~eso_internal_fwhm_ang Gaussian (see ESO_LSF_GAUSS_FWHM_PIX
        # above). Quadrature-subtract that from the target LSF so the
        # *combined* effective FWHM matches lsf_fwhm, instead of naively
        # stacking two convolutions and overshooting it.
        target_fwhm = _clean_lsf(lsf_fwhm)
        add_var = target_fwhm ** 2 - eso_internal_fwhm_ang ** 2
        floored = add_var < 0.0
        n_floored[0] += int(np.sum(floored))
        n_floored_total[0] += floored.size
        add_fwhm = np.sqrt(np.clip(add_var, 0.0, None))
        sigma_pix = np.clip(add_fwhm / 2.355 / dlam_pix, 1.0e-6, None)
        kernel = norm.pdf(taps[None, :], loc=0.0, scale=sigma_pix[:, None])
        kernel /= kernel.sum(axis=1, keepdims=True)
        operator = build_lsf_operator(wave_ref, kernel)
        return np.asarray(spec[None, :] @ operator.T).ravel()

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
        if not args.no_lsf_convolve:
            # Apply this row's own real instrumental LSF -- ESO's local engine
            # only applies its own fixed, LVM-untuned ~0.4 A kernel internally
            # (see module Description), so without this step the ESO candidate
            # would be compared unconvolved while the MLP candidate already has
            # its own per-row LSF reconstruction applied.
            flux_pred = _convolve(flux_pred, lsf_all[idx])
            line_pred = _convolve(line_pred, lsf_all[idx])
        flux_pred_list.append(flux_pred)
        line_pred_list.append(line_pred)
        flux_obs_list.append(flux_obs_all[idx])
        expnum_list.append(r['expnum'])
        row_list.append(idx)

    n_ok = len(row_list)
    print(f'\n{n_ok}/{len(tasks)} rows succeeded ({n_fail} failed) in {total:.1f}s '
          f'({total / max(n_ok, 1):.2f}s/row wall-clock, {args.n_workers} workers)')
    if n_floored[0]:
        print(f'  Note: {n_floored[0]}/{n_floored_total[0]} wavelength pixels had a target LSF '
              f'FWHM narrower than ESO\'s own internal {eso_internal_fwhm_ang:.3f} A kernel -- '
              f'the quadrature correction floored at zero there (no additional narrowing possible; '
              f'output at those pixels is ESO\'s own internal resolution, not the true target LSF).')

    outfile = args.outfile
    if outfile is None:
        suffix = '_eso_nolsf.fits' if args.no_lsf_convolve else '_eso_lsf.fits'
        outfile = str(Path(args.fits_file).with_suffix('').as_posix()) + suffix

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
