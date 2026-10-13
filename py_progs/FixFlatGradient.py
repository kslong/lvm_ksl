#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Apply to DRP 1.3.2 lvmCFrame files the science-IFU gradient that the
    master fiber flats record in their headers but that the DRP never
    applied, removing most of the large-scale (spectrograph-to-
    spectrograph) pattern across the science IFU without rerunning the
    pipeline.

Command line usage (if any):

    usage: FixFlatGradient.py [-h] [-outdir DIR] [-calib DIR] [-epochs FILE]
                              [-epoch MJD] cframe [cframe ...]

    Arguments::

        cframe       one or more lvmCFrame files (DRP 1.3.2)

    Options::

        -h           print this help and exit
        -outdir DIR  directory for the corrected files (default: the
                     current directory); each output is the input's
                     basename with '.gradfix' inserted before '.fits'
        -calib DIR   calibration directory holding <epoch>/
                     lvm-mfiberflat_twilight-{b,r,z}.fits (default
                     $LVM_MASTER_DIR)
        -epochs FILE lvmcore calibration-epochs.yaml (default
                     $LVMCORE_DIR/calibrations/calibration-epochs.yaml)
        -epoch MJD   use this calibration epoch for every file instead of
                     looking it up

Description::

    In DRP 1.3.2 the per-spectrograph fiber-flat factors are measured once
    per "ffactor epoch" from science exposures, jointly with a planar IFU
    gradient, and multiplied into the master twilight flats by
    lvmdrp run_twilights.do_ffactor_correction.  That routine keeps only
    the factors part of the model,

        _, _, flatfield_corr = mflat.eval_ifu_gradient(...)

    (eval_ifu_gradient returns joint, gradient, factors), so the gradient
    is written to the flat header (<ARM> FIBERFLAT GCOEFF1-4) but not
    into the flat.  Because the factors were fitted together with the
    gradient they do not absorb it, and every science exposure keeps the
    gradient -- seen as sp1 high / sp3 low by up to 3-4% (growing with
    time) plus a smooth gradient across the IFU.

    For each file this routine:

    1. finds the calibration epoch: the latest epoch in
       calibration-epochs.yaml with epoch MJD <= the exposure MJD (this
       reproduces drpall's calib_mjd for all 21793 1.3.2 exposures
       checked), unless -epoch is given;
    2. for each arm reads GCOEFF1-4 from that epoch's
       lvm-mfiberflat_twilight-<arm>.fits and evaluates the model exactly
       as lvmdrp IFUGradient.ifu_gradient does: G = c1 + c2*y + c3*x +
       c4*x*y with x, y = SLITMAP xpmm, ypmm, normalized to its mean over
       the Sci fibers (the DRP normalizes over all fibers; the two
       differ by < 0.01% because the IFU is centred on x = y = 0);
    3. builds a per-fiber, per-pixel correction from the three arms,
       using the arm limits of lvmdrp SPEC_CHANNELS (b 3600-5800,
       r 5775-7570, z 7520-9800) and a linear ramp across the two
       overlaps (5775-5800, 7520-7570) -- the DRP combines the arms
       there with inverse-variance weights that the CFrame does not
       keep;
    4. divides FLUX by the correction and multiplies IVAR by its square
       (what a flat multiplied by G would have given) for the Sci fibers
       only.  The gradient was measured on the science IFU; the SkyE,
       SkyW and Spec fibers sit in their own bundles, so G evaluated at
       their xpmm, ypmm has no physical meaning (for Spec it would scale
       individual standard-star fibers by up to +-7%).  Their rows and
       the SKY_EAST/SKY_WEST extensions are left unchanged (G averages to
       1.0000 over the SkyE and SkyW fibers in any case);
    5. records the epoch and coefficients in the primary header
       (FGRADFIX, FGRDEPOC, FGRD<arm><k>) and refuses files already
       carrying FGRADFIX.

Primary routines:

    find_epoch, read_gradient, correction, fix_cframe

Notes:

    The flux calibration is not redone and is not affected for FLUXCAL
    MOD (or STD): those sensitivities come from the standard-star (Spec)
    fibers, which are not changed, and G averages to 1 over the Sci
    fibers, so the mean science flux level is unchanged; only its
    distribution across the IFU changes.  For FLUXCAL SCI the
    sensitivity came from field stars in the science IFU, which carried
    the uncorrected gradient at their positions, so a small zero-point
    error remains; a warning is printed for such files.

    The correction removes the error made in building the flats; it does
    not follow changes within an ffactor epoch, which DRP 1.3.2 cannot
    either (about +-1-1.5% sp1/sp3 late in 2025-2026).  SFrames cannot
    be corrected by this routine -- correct the CFrame and redo the sky
    subtraction.

    Investigation and validation: lvm_moonsky README (Rainbow3),
    scripts ffactor_windows.py, gradient_sector.py,
    apply_header_gradient.py.

History::

    261009 ksl Coding begun.  Corrects Sci fibers only (Spec fibers'
               nominal xpmm, ypmm would scale standard stars by up
               to +-7%); validated on 60 moonsky exposures against
               DRP 1.2.1 and a per-fiber test (agreement 0.01-0.07%)

'''

import os
import re
import sys
import numpy as np
import yaml
from astropy.io import fits

def _usage_from_doc(doc):
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc

_USAGE = _usage_from_doc(__doc__)

# lvmdrp 1.3.2 core/constants.py SPEC_CHANNELS
SPEC_CHANNELS = {'b': (3600., 5800.), 'r': (5775., 7570.), 'z': (7520., 9800.)}


def default_epochs_file():
    core = os.getenv('LVMCORE_DIR', '')
    return os.path.join(core, 'calibrations', 'calibration-epochs.yaml') if core else ''


def load_epochs(epochs_file):
    '''sorted array of calibration-epoch MJDs, or None'''
    try:
        with open(epochs_file) as f:
            return np.array(sorted(int(k) for k in yaml.safe_load(f)['epochs']))
    except Exception as err:
        print('FixFlatGradient: could not read epochs file %s (%s)' % (epochs_file, err))
        return None


def find_epoch(mjd, epochs):
    '''latest calibration epoch <= mjd, or None'''
    i = np.searchsorted(epochs, mjd, side='right') - 1
    return int(epochs[i]) if i >= 0 else None


def read_gradient(calib, epoch, arm):
    '''GCOEFF1-4 for one arm of a master twilight flat, or None'''
    fname = os.path.join(calib, str(epoch), 'lvm-mfiberflat_twilight-%s.fits' % arm)
    try:
        h = fits.getheader(fname)
        return [float(h['%s FIBERFLAT GCOEFF%d' % (arm.upper(), k)]) for k in (1, 2, 3, 4)]
    except Exception as err:
        print('FixFlatGradient: no gradient coefficients in %s (%s)' % (fname, err))
        return None


def arm_weights(wave):
    '''per-pixel weight of each arm: 1 inside an arm, linear ramps across the overlaps'''
    w = {arm: np.zeros(wave.size) for arm in 'brz'}
    br0, br1 = SPEC_CHANNELS['r'][0], SPEC_CHANNELS['b'][1]
    rz0, rz1 = SPEC_CHANNELS['z'][0], SPEC_CHANNELS['r'][1]
    w['b'][wave < br0] = 1.
    m = (wave >= br0) & (wave <= br1)
    w['r'][m] = (wave[m] - br0) / (br1 - br0)
    w['b'][m] = 1. - w['r'][m]
    w['r'][(wave > br1) & (wave < rz0)] = 1.
    m = (wave >= rz0) & (wave <= rz1)
    w['z'][m] = (wave[m] - rz0) / (rz1 - rz0)
    w['r'][m] = 1. - w['z'][m]
    w['z'][wave > rz1] = 1.
    return w


def correction(coeffs, x, y, wave):
    '''(fibers, pixels) correction for the given fibers: sum over arms of arm weight x
    gradient normalized to its mean over those fibers'''
    w = arm_weights(wave)
    corr = np.zeros((x.size, wave.size))
    grads = {}
    for arm in 'brz':
        c = coeffs[arm]
        g = c[0] + c[1] * y + c[2] * x + c[3] * x * y
        grads[arm] = g / np.nanmean(g)
        corr += grads[arm][:, None] * w[arm][None, :]
    return corr, grads


def fix_cframe(filename, outdir='.', calib='', epochs=None, epoch=None):
    '''write a gradient-corrected copy of one CFrame; returns the output name or None'''
    try:
        x = fits.open(filename)
    except Exception as err:
        print('FixFlatGradient: could not open %s (%s)' % (filename, err))
        return None
    hdr = x['PRIMARY'].header
    if hdr.get('FGRADFIX', False):
        print('FixFlatGradient: %s already corrected, skipped' % filename)
        return None
    if hdr.get('DRPVER', '') != '1.3.2':
        print('FixFlatGradient: warning: %s has DRPVER %s; the diagnosis was made for 1.3.2'
              % (filename, hdr.get('DRPVER')))
    if epoch is None:
        epoch = find_epoch(int(hdr['MJD']), epochs)
        if epoch is None:
            print('FixFlatGradient: no calibration epoch for MJD %s (%s)' % (hdr['MJD'], filename))
            return None
    coeffs = {arm: read_gradient(calib, epoch, arm) for arm in 'brz'}
    if any(c is None for c in coeffs.values()):
        return None

    if hdr.get('FLUXCAL', 'MOD') not in ('MOD', 'STD'):
        print('FixFlatGradient: warning: %s has FLUXCAL %s; its sensitivity came from uncorrected '
              'science-IFU data, so a small zero-point error remains' % (filename, hdr.get('FLUXCAL')))

    st = x['SLITMAP'].data
    sci = np.char.strip(np.asarray(st['telescope']).astype(str)) == 'Sci'
    xp, yp = np.asarray(st['xpmm'], float)[sci], np.asarray(st['ypmm'], float)[sci]
    wave = np.asarray(x['WAVE'].data, float)
    corr, grads = correction(coeffs, xp, yp, wave)

    flux, ivar = x['FLUX'].data, x['IVAR'].data
    flux[sci] = (flux[sci] / corr).astype(flux.dtype)
    ivar[sci] = (ivar[sci] * corr ** 2).astype(ivar.dtype)

    hdr['FGRADFIX'] = (True, 'Sci IFU gradient from flat header applied')
    hdr['FGRDEPOC'] = (epoch, 'calibration epoch of the flats used')
    for arm in 'brz':
        for k, c in enumerate(coeffs[arm]):
            hdr['FGRD%s%d' % (arm.upper(), k + 1)] = (c, '%s FIBERFLAT GCOEFF%d applied' % (arm, k + 1))
    hdr.add_history('FixFlatGradient: Sci-fiber FLUX/IVAR corrected for the IFU gradient in the '
                    'flat headers of epoch %d (not applied by lvmdrp do_ffactor_correction)' % epoch)

    base = os.path.basename(filename)
    root = base[:-5] if base.endswith('.fits') else base
    outname = os.path.join(outdir, root + '.gradfix.fits')
    x.writeto(outname, overwrite=True)
    x.close()
    span = {arm: 100 * (grads[arm].max() / grads[arm].min() - 1) for arm in 'brz'}
    print('FixFlatGradient: %s -> %s  epoch %d  Sci gradient max/min b %.1f%% r %.1f%% z %.1f%%'
          % (base, outname, epoch, span['b'], span['r'], span['z']))
    return outname


def steer(argv):
    files, outdir, calib, epochs_file, epoch = [], '.', os.getenv('LVM_MASTER_DIR', ''), default_epochs_file(), None
    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-outdir':
            i += 1
            outdir = argv[i]
        elif argv[i] == '-calib':
            i += 1
            calib = argv[i]
        elif argv[i] == '-epochs':
            i += 1
            epochs_file = argv[i]
        elif argv[i] == '-epoch':
            i += 1
            epoch = int(argv[i])
        elif argv[i][0] == '-':
            print('Error: unknown option %s' % argv[i])
            print(_USAGE)
            return
        else:
            files.append(argv[i])
        i += 1
    if not files:
        print(_USAGE)
        return
    if not calib:
        print('FixFlatGradient: no calibration directory: set $LVM_MASTER_DIR or use -calib')
        return
    epochs = None
    if epoch is None:
        if not epochs_file:
            print('FixFlatGradient: no epochs file: set $LVMCORE_DIR, or use -epochs or -epoch')
            return
        epochs = load_epochs(epochs_file)
        if epochs is None:
            return
    os.makedirs(outdir, exist_ok=True)
    for f in files:
        fix_cframe(f, outdir, calib, epochs, epoch)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
