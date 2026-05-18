#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Add a synthetic Gaussian emission line to the FLUX extension of an
lvmCFrame file, for use in testing fourier_offset.py.

Command line usage (if any):

    usage: check_fourier.py [-h] [-wave 4000] [-flux 1e-13]
                            [-off1 0.0] [-off2 0.1] [-off3 -0.1] filename

    where filename is an lvmCFrame FITS file.
    -wave sets the central wavelength of the line in Angstroms (default 4200).
    -flux sets the integrated line flux in erg/s/cm^2 (default 1e-11).
    -off1, -off2, -off3 set wavelength offsets in Angstroms applied to all
    fibers in spectrographs 1, 2, and 3 respectively (defaults 0, +0.1, -0.1).

Description:

    Reads the WAVE, FLUX, LSF, and SLITMAP extensions of the input file.
    For each fiber the FWHM at the requested wavelength is taken from the
    LSF extension and converted to a Gaussian sigma.  The line centre is
    shifted by the per-spectrograph offset before the Gaussian is computed.
    The peak amplitude is derived from the integrated flux as
    A = flux / (sigma * sqrt(2*pi)), and the Gaussian is added to the FLUX
    array.  The modified file is written to test_<basename>.fits in the
    current directory.

Primary routines:

    do_one

Notes:

    The LSF extension is assumed to contain FWHM values in Angstroms on
    the same wavelength grid as the WAVE extension.

History:

260428 ksl Coding begun

'''

import os
import sys
import numpy as np
from astropy.io import fits
from astropy.table import Table


def do_one(filename, wave=4200.0, flux=1e-11, off1=0.0, off2=0.1, off3=-0.1):
    '''
    Add a synthetic Gaussian emission line to the FLUX extension of a
    CFrame file and write the result to test_<basename>.fits.

    Parameters
    ----------
    filename : str
        Path to the input lvmCFrame FITS file.
    wave : float
        Central wavelength of the line in Angstroms.
    flux : float
        Integrated line flux in erg/s/cm^2.
    off1, off2, off3 : float
        Wavelength offset in Angstroms applied to all fibers in
        spectrographs 1, 2, and 3 respectively.
    '''
    try:
        hdul = fits.open(filename)
    except Exception as e:
        print('Error: could not open %s: %s' % (filename, e))
        return

    wave_arr  = hdul['WAVE'].data.astype(float)    # (N_wave,)
    flux_data = hdul['FLUX'].data.astype(float)    # (N_fibers, N_wave)
    lsf_data  = hdul['LSF'].data.astype(float)     # (N_fibers, N_wave) FWHM in AA
    slit_tab  = Table(hdul['SLITMAP'].data)

    # Per-fiber line centre from spectrograph offsets
    offsets = {1: off1, 2: off2, 3: off3}
    spec_id = np.array(slit_tab['spectrographid'])
    centers = np.full(len(spec_id), wave)
    for sp, off in offsets.items():
        centers[spec_id == sp] += off

    # FWHM at the nearest wavelength pixel, per fiber
    idx   = int(np.argmin(np.abs(wave_arr - wave)))
    fwhm  = lsf_data[:, idx]                              # (N_fibers,)
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))    # (N_fibers,)

    # Only add the Gaussian to fibers with a valid (non-zero) LSF
    good  = sigma > 0.0
    n_good = good.sum()
    n_bad  = (~good).sum()
    print('Adding line to %d fibers (%d skipped — zero LSF)' % (n_good, n_bad))

    # Peak amplitude from integrated flux: F = A * sigma * sqrt(2*pi)
    amplitude = np.zeros_like(sigma)
    amplitude[good] = flux / (sigma[good] * np.sqrt(2.0 * np.pi))

    # Vectorised Gaussian — each fiber has its own centre
    dwave = wave_arr[np.newaxis, :] - centers[:, np.newaxis]   # (N_fibers, N_wave)
    gauss = np.zeros_like(flux_data)
    gauss[good] = amplitude[good, np.newaxis] * np.exp(
        -0.5 * (dwave[good] / sigma[good, np.newaxis]) ** 2
    )

    flux_data += gauss
    hdul['FLUX'].data = flux_data.astype(hdul['FLUX'].data.dtype)

    outname = 'test_' + os.path.basename(filename)
    hdul.writeto(outname, overwrite=True)
    hdul.close()
    print('Wrote %s  (line at %.1f AA, flux %.2e, offsets sp1=%.3f sp2=%.3f sp3=%.3f AA)'
          % (outname, wave, flux, off1, off2, off3))


def steer(argv):
    wave     = 4200.0
    flux     = 1e-11
    off1     =  0.0
    off2     =  0.1
    off3     = -0.1
    filename = ''

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-wave':
            i += 1
            wave = float(argv[i])
        elif argv[i] == '-flux':
            i += 1
            flux = float(argv[i])
        elif argv[i] == '-off1':
            i += 1
            off1 = float(argv[i])
        elif argv[i] == '-off2':
            i += 1
            off2 = float(argv[i])
        elif argv[i] == '-off3':
            i += 1
            off3 = float(argv[i])
        elif argv[i][0] == '-':
            print('Error: unknown switch: %s' % argv[i])
            return
        else:
            filename = argv[i]
        i += 1

    if filename == '':
        print(__doc__)
        return

    do_one(filename, wave, flux, off1, off2, off3)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
