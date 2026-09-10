#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Dispatcher for the SkySub* family (SkySubDrp.py, SkySubOrig.py,
    SkySubDev1.py, SkySubDev2.py, SkySubDev3.py): runs one named
    sky-subtraction routine on an XCframe file by calling its do_all()
    directly (no subprocess), writes its output into a routine-specific
    subdirectory so different routines -- or repeated runs of the same
    routine with a different variant -- never collide on a filename, and
    optionally runs SkySub_eval.py on the result afterward.

Command line usage (if any):

    usage: SkySubRun.py -routine {drp,orig,dev1,dev2,dev3} [-variant NAME]
                        [-delta N] [-mask FILE] [-kstep N]
                        [-fwhm_lsf F] [-lsf_boost F]
                        [-v VEL] [-lmc] [-smc] [-lvmsky_skysub PATH]
                        [-outdir DIR] [-eval] [-eval_out ROOT]
                        filename

    where

    filename        XCframe FITS file to process.

    -routine NAME   which SkySub* routine to run. One of::

                        drp    SkySubDrp.py  (calls lvmdrp's own
                               create_skysub_spectrum -- "Method 0")
                        orig   SkySubOrig.py (from-scratch near/far
                               recipe, polynomial continuum)
                        dev1   SkySubDev1.py (near/far recipe,
                               GetSkyCont.py B-spline continuum)
                        dev2   SkySubDev2.py (near/far recipe, PALACE
                               decomposition)
                        dev3   SkySubDev3.py (near/far recipe, SkyDecomp
                               continuum + nebular-line masking --
                               "Method 1")

                    Required -- there is no default, so a run always
                    names its own routine explicitly rather than falling
                    back to whichever happened to be current.

    -variant NAME   the sky-construction variant passed as that
                    routine's own -method (e.g. nearest,
                    farlines_nearcont, farthest, scilines_nearcont --
                    see each routine's own Synopsis for which variants
                    it supports). Default: that routine's own default
                    variant.

    -delta N        process every N-th row (default: 1 = all rows).

    -mask FILE      palace_mask FITS file -- required by -routine dev1
                    only. If omitted, searched for as sky_mask.fits in
                    the current directory then the lvm_ksl data/
                    directory, same fallback SkySubDev1.py's own CLI
                    uses.

    -kstep N        B-spline knot spacing in Angstrom -- dev1 only
                    (default 100).

    -fwhm_lsf F / -lsf_boost F
                    PALACE LSF knobs -- dev2 only (defaults 1.3 / 1.0).

    -v VEL / -lmc / -smc
                    nebular systemic velocity for the exclusion mask --
                    dev3 only (default: 0). Same convention as
                    DecomposeCleanSky.py/sky_nebular_leak_eval.py.

    -lvmsky_skysub PATH
                    path to lvmsky's skysub/ directory -- dev3 only
                    (default: ~/SDSS/lvmsky/skysub).

    -outdir DIR     top-level directory under which each routine gets
                    its own subdirectory, DIR/<routine>/ (default:
                    sky_runs). Created if it doesn't exist.

    -eval           after the routine finishes, also run
                    SkySub_eval.plot_eval() on its output (same function
                    SkySub_eval.py's own CLI calls -- this is a
                    convenience, not a separate code path; running
                    SkySub_eval.py by hand on the written file afterward
                    gives identical results).

    -eval_out ROOT  output root for the -eval HTML (default: the
                    routine's own output stem + "_eval", in the same
                    subdirectory).

Description:

    Output naming: DIR/<routine>/<input_stem>_<routine>_<variant>.fits.
    Passing this explicit, non-empty outroot to every routine's do_all()
    means each one writes exactly there and nowhere else -- none of
    SkySubDrp/Orig/Dev1/Dev2/Dev3.py needed to change; every one of them
    already uses outroot verbatim (appending only ".fits") whenever it's
    given non-empty. The per-routine subdirectory means two different
    routines, or the same routine run with two different variants, can
    never collide even if someone reuses an outdir across many runs.

    Each routine's module is imported lazily (only once -routine
    selects it), so running e.g. -routine drp never pays the cost of
    importing dev3's lvmsky dependency or dev1/dev2's GetSkyCont/PALACE
    machinery.

Primary routines:

    run_one   dispatch to one routine, return the output FITS path.

Notes:

    This script calls each routine's do_all() function directly (an
    ordinary Python function call within this process), not via
    subprocess -- faster, and errors surface as normal Python
    tracebacks instead of being swallowed into a subprocess return code.

History::

    260909  ksl  Coding begun.

'''

import argparse
import os
import sys
import importlib
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ──────────────────────────────────────────────────────────────
# Registry: routine key -> module name, allowed variants, default
# variant, and which extra CLI-level knobs it consumes.
# ──────────────────────────────────────────────────────────────

ROUTINES = {
    'drp':  dict(module='SkySubDrp',
                variants=('nearest', 'farlines_nearcont'),
                default_variant='farlines_nearcont', extra=()),
    'orig': dict(module='SkySubOrig',
                variants=('nearest', 'farthest', 'farlines_nearcont'),
                default_variant='farlines_nearcont', extra=()),
    'dev1': dict(module='SkySubDev1',
                variants=('nearest', 'farlines_nearcont'),
                default_variant='farlines_nearcont', extra=('mask', 'kstep')),
    'dev2': dict(module='SkySubDev2',
                variants=('scilines_nearcont', 'nearest', 'farlines_nearcont'),
                default_variant='scilines_nearcont', extra=('fwhm_lsf', 'lsf_boost')),
    'dev3': dict(module='SkySubDev3',
                variants=('nearest', 'farlines_nearcont'),
                default_variant='farlines_nearcont', extra=('vel', 'lvmsky_skysub')),
}


def _resolve_mask_file(mask_file):
    '''Same cwd-then-data/ fallback SkySubDev1.py's own CLI uses.'''
    if mask_file:
        if not os.path.exists(mask_file):
            raise FileNotFoundError(f'mask file not found: {mask_file}')
        return mask_file
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
    for candidate in (os.path.join(os.getcwd(), 'sky_mask.fits'),
                      os.path.join(data_dir, 'sky_mask.fits')):
        if os.path.exists(candidate):
            print(f'Using default mask: {candidate}')
            return candidate
    raise FileNotFoundError(
        '-mask not given and sky_mask.fits not found in the current '
        'directory or data/ (required for -routine dev1)')


def run_one(routine, filename, variant=None, delta=1, outdir='sky_runs',
           mask=None, kstep=100.0, fwhm_lsf=1.3, lsf_boost=1.0,
           vel=0.0, lvmsky_skysub=None):
    '''
    Run one SkySub* routine on filename, writing into
    outdir/routine/<stem>_<routine>_<variant>.fits.

    Parameters
    ----------
    routine : str
        One of ROUTINES' keys.
    filename : str
        XCframe FITS file to process.
    variant : str, optional
        The routine's own -method variant; defaults to its registry
        default_variant.
    delta : int
        Row step size (1 = all rows).
    outdir : str
        Top-level output directory; routine gets its own subdirectory
        under it.
    mask, kstep, fwhm_lsf, lsf_boost, vel, lvmsky_skysub
        Routine-specific extra options; only the ones the selected
        routine's registry entry lists in "extra" are actually used.

    Returns
    -------
    str
        Path to the FITS file written.
    '''
    if routine not in ROUTINES:
        raise ValueError(f'-routine must be one of: {", ".join(sorted(ROUTINES))}')
    spec = ROUTINES[routine]
    variant = variant or spec['default_variant']
    if variant not in spec['variants']:
        raise ValueError(f'-routine {routine} does not support -variant {variant!r}; '
                         f'valid: {", ".join(spec["variants"])}')

    abs_filename = str(Path(filename).expanduser().resolve())
    if not os.path.exists(abs_filename):
        raise FileNotFoundError(f'file not found: {filename}')

    rundir = os.path.join(outdir, routine)
    os.makedirs(rundir, exist_ok=True)
    stem = Path(filename).stem
    outroot = os.path.join(rundir, f'{stem}_{routine}_{variant}')

    kwargs = dict(filename=abs_filename, method=variant, idelta=delta, outroot=outroot)
    if 'mask' in spec['extra']:
        kwargs['mask_file'] = _resolve_mask_file(mask)
    if 'kstep' in spec['extra']:
        kwargs['knot_step'] = kstep
    if 'fwhm_lsf' in spec['extra']:
        kwargs['fwhm_lsf'] = fwhm_lsf
    if 'lsf_boost' in spec['extra']:
        kwargs['lsf_boost'] = lsf_boost
    if 'vel' in spec['extra']:
        kwargs['vel'] = vel
    if 'lvmsky_skysub' in spec['extra'] and lvmsky_skysub:
        kwargs['lvmsky_skysub'] = lvmsky_skysub

    print(f'Running -routine {routine} (module {spec["module"]}), '
         f'variant={variant}, delta={delta}')
    module = importlib.import_module(spec['module'])
    module.do_all(**kwargs)

    outfile = f'{outroot}.fits'
    print(f'-routine {routine} wrote {outfile}')
    return outfile


# ──────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description='Dispatch one SkySub* routine and optionally evaluate its output.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-routine', required=True, choices=sorted(ROUTINES),
                   help='which SkySub* routine to run')
    p.add_argument('-variant', default=None,
                   help="that routine's own -method variant (default: routine-specific)")
    p.add_argument('-delta', type=int, default=1, help='process every N-th row')
    p.add_argument('-mask', default=None, help='palace_mask FITS file (dev1 only)')
    p.add_argument('-kstep', type=float, default=100.0, help='B-spline knot spacing, A (dev1 only)')
    p.add_argument('-fwhm_lsf', type=float, default=1.3, help='PALACE LSF FWHM, A (dev2 only)')
    p.add_argument('-lsf_boost', type=float, default=1.0, help='PALACE LSF boost (dev2 only)')
    p.add_argument('-v', dest='vel', type=float, default=None, help='nebular velocity, km/s (dev3 only)')
    p.add_argument('-lmc', action='store_true', help='LMC velocity shortcut for -v (dev3 only)')
    p.add_argument('-smc', action='store_true', help='SMC velocity shortcut for -v (dev3 only)')
    p.add_argument('-lvmsky_skysub', default=None, help='lvmsky skysub/ dir (dev3 only)')
    p.add_argument('-outdir', default='sky_runs', help='top-level output directory')
    p.add_argument('-eval', action='store_true', help='also run SkySub_eval.py on the output')
    p.add_argument('-eval_out', default=None, help='output root for -eval (default: <output>_eval)')
    p.add_argument('filename', help='XCframe FITS file to process')
    args = p.parse_args()

    vel = args.vel
    if vel is None:
        vel = 262. if args.lmc else (146. if args.smc else 0.)

    outfile = run_one(args.routine, args.filename, variant=args.variant, delta=args.delta,
                      outdir=args.outdir, mask=args.mask, kstep=args.kstep,
                      fwhm_lsf=args.fwhm_lsf, lsf_boost=args.lsf_boost,
                      vel=vel, lvmsky_skysub=args.lvmsky_skysub)

    if args.eval:
        from SkySub_eval import plot_eval
        eval_out = args.eval_out or (os.path.splitext(outfile)[0] + '_eval')
        print(f'Running SkySub_eval on {outfile} -> {eval_out}.html')
        plot_eval([outfile], outroot=eval_out)


if __name__ == '__main__':
    main()
