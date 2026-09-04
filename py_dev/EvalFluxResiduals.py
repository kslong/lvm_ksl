#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Step 1 (general, method-agnostic accuracy quantification): run
    py_progs/sky_residual_eval.py's batch (wave, flux, model) evaluator
    over a BatchPredictSky.py output file's FLUX_OBS/FLUX_PRED pairs.
    This script has no knowledge of how FLUX_PRED was produced -- it
    would work identically on a batch of ESO/PALACE/any-other-candidate
    predicted skies, as long as they come in the same WAVE/FLUX_OBS/
    FLUX_PRED layout.

Command line usage (if any):

    usage: EvalFluxResiduals.py [-h] [-mask MASK] [-np N] [-out ROOT]
                                [-plotdir DIR] [-residuals FITS]
                                batch_fits

    where

    -mask MASK      palace_make_mask.py mask FITS (default:
                    data/sky_mask.fits, same as sky_residual_eval.py's
                    own default).

    -np N           worker processes (default: 8; matches
                    py_progs/Reduce.py and py_progs/sky_gaussfit.py's
                    process-count convention).

    -out ROOT       output filename root (default: <batch_fits stem>).

    -plotdir DIR    directory for summary plots (default: plots_sky_resid).

    -residuals FITS optional EvalCoefResiduals.py output table -- if
                    given, its expnum/split (train/val/test) columns are
                    merged into the summary table by expnum, so
                    in-sample vs. genuinely-held-out rows can be told
                    apart afterward. Purely informational (join key
                    only); sky_residual_eval.py itself never sees it.

    -py_progs_dir PATH
                    path to the lvm_ksl repo's py_progs/ directory, which
                    supplies sky_residual_eval.py (default:
                    ~/SDSS/lvm_ksl/py_progs).

    batch_fits      a BatchPredictSky.py output FITS file (WAVE,
                    FLUX_OBS, FLUX_PRED, META[row,expnum]).

Description:

    Thin wrapper around sky_residual_eval.analyze_sky_residuals() +
    its three summary-plot routines, following the same pattern as that
    module's own CLI (main()) but reading our WAVE/FLUX_OBS/FLUX_PRED
    file layout instead of its WAVE/FLUX/SKY convention.

History::

    260902  ksl  Coding begun.
    260904  ksl  -nproc renamed to -np, matching py_progs/Reduce.py/
        sky_gaussfit.py's process-count spelling. batch_fits (the only
        positional) moved to the end of the "where" list in this
        docstring to match where it already appears in the usage line
        and -h output -- no code change, just doc ordering.

'''

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table, join

DEFAULT_PY_PROGS_DIR = Path('~/SDSS/lvm_ksl/py_progs').expanduser()

_pre = argparse.ArgumentParser(add_help=False)
_pre.add_argument('-py_progs_dir', default=str(DEFAULT_PY_PROGS_DIR))
_pre_args, _ = _pre.parse_known_args()

sys.path.insert(0, _pre_args.py_progs_dir)
from sky_residual_eval import (  # noqa: E402
    analyze_sky_residuals, plot_frac_summary, plot_continuum_summary,
    plot_lines_summary,
)


def main():
    p = argparse.ArgumentParser(
        parents=[_pre],
        description="Method-agnostic flux-space accuracy evaluation of a "
                    "BatchPredictSky.py output file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('-mask', default=None, help='palace_make_mask.py mask FITS')
    p.add_argument('-np', dest='nproc', type=int, default=8, help='worker processes')
    p.add_argument('-out', default=None, help='output filename root')
    p.add_argument('-plotdir', default='plots_sky_resid', help='summary-plot directory')
    p.add_argument('-residuals', default=None,
                   help='EvalCoefResiduals.py output table, for expnum/split merge')
    p.add_argument('batch_fits', help='BatchPredictSky.py output FITS file')
    args = p.parse_args()

    with fits.open(args.batch_fits) as hdul:
        wave = np.asarray(hdul['WAVE'].data, dtype=float)
        flux_obs = np.asarray(hdul['FLUX_OBS'].data, dtype=float)
        flux_pred = np.asarray(hdul['FLUX_PRED'].data, dtype=float)
        meta = Table(hdul['META'].data)

    print(f'Evaluating {flux_obs.shape[0]} rows (nproc={args.nproc}) ...')
    summary_table, line_table = analyze_sky_residuals(
        wave, flux_obs, flux_pred, mask=args.mask, nproc=args.nproc)

    summary_table['expnum'] = meta['expnum']
    line_table['expnum'] = meta['expnum'][np.asarray(line_table['ROW'])]

    if args.residuals:
        res = Table.read(args.residuals)[['expnum', 'split']]
        summary_table = join(summary_table, res, keys='expnum', join_type='left')
        summary_table['split'] = summary_table['split'].filled('unknown') \
            if hasattr(summary_table['split'], 'filled') else summary_table['split']

    outroot = args.out or Path(args.batch_fits).stem
    summary_table.write(f'{outroot}_summary.fits', overwrite=True)
    line_table.write(f'{outroot}_lines.fits', overwrite=True)
    print(f'Wrote {outroot}_summary.fits ({len(summary_table)} rows) '
          f'and {outroot}_lines.fits ({len(line_table)} rows)')

    # plot_*_summary's outname is joined with outdir internally
    # (os.path.join-style) -- must be a bare filename, not outroot's full
    # path, or an absolute outname silently discards outdir.
    outstem = Path(outroot).name
    plot_frac_summary(summary_table, outdir=args.plotdir,
                      outname=f'{outstem}_frac_summary.png')
    plot_continuum_summary(summary_table, outdir=args.plotdir,
                           outname=f'{outstem}_continuum_summary.png')
    plot_lines_summary(summary_table, outdir=args.plotdir,
                       outname=f'{outstem}_lines_summary.png')
    print(f'Wrote summary plots to {args.plotdir}/')


if __name__ == '__main__':
    main()
