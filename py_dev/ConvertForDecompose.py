#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Reformat an LVM XCframe-layout FITS file (WAVE/FLUX/SKY_EAST/SKY_WEST/
    LSF/DRP_ALL) into the stack layout lvmsky's decompose_parallel.py
    expects (WAVE, FLUX_SCI, FLUX_SKY_NEAR, FLUX_SKY_FAR, LSF_SCI,
    LSF_SKY_NEAR, LSF_SKY_FAR, META) -- the same layout lvmsky's own
    lvm_medians tool writes.

Command line usage (if any):

    usage: ConvertForDecompose.py [-h] fits_file [output_file]

    where

    fits_file       is the path to an XCframe-layout FITS file (e.g. the
                    output of SelectXCF.py, preferably run with -hdr so
                    that DRP_ALL carries pwv_med).

    output_file     output FITS path (default: <stem>_decomp_input.fits).

Description:

    1. Copies WAVE unchanged and renames FLUX -> FLUX_SCI.
    2. Builds FLUX_SKY_NEAR / FLUX_SKY_FAR by selecting, row by row,
       SKY_EAST or SKY_WEST according to that row's DRP_ALL Near/Far label
       ('SKY_EAST'/'SKY_WEST') -- which telescope is nearer varies from
       exposure to exposure.
    3. Writes LSF_SCI, LSF_SKY_NEAR and LSF_SKY_FAR.  The XCframe's LSF
       extension holds only the science fibers' LSF, so it is used for all
       three (see Notes).  Non-finite or non-positive LSF pixels -- in
       practice the NaNs at the ends of the wavelength range -- are first
       filled from the nearest valid pixels (fill_lsf).
    4. Copies DRP_ALL through as META, adding the columns, under the exact
       lower-case names, that decompose_parallel.py's telluric fit models
       and mlp_predictor read:

         sky_near_ra/dec, sky_far_ra/dec   per-row Near/Far selection
         sky_near_label, sky_far_label     'SkyE'/'SkyW'
         sci_airmass, skye_airmass,        from sci_amass/skye_amass/
         skyw_airmass                      skyw_amass
         date_obs                          from obstime
         pwv_med                           from DRP_ALL pwv_med if present
                                           (SelectXCF.py -hdr), else NaN

       decompose_parallel.py's worker checks these names case-sensitively,
       so they must be lower case.

Notes::

    PWV: the telluric fit models recompute the DRP's own telluric
    transmission, which used PWV_MED from the CFrame header.  A row with
    no valid pwv_med falls back to the DRP default of 15 mm, which is
    wrong for most exposures (typically 1-4 mm); a warning reports how
    many rows that affects.

    LSF: the sky telescopes' LSF is approximated by the science fibers'
    LSF until SummarizeCframe.py records per-telescope LSFs.  The
    decomposition still refines the LSF from each spectrum's own sky
    lines; the input LSF is used for the telluric transmission, the
    moon/zodiacal-light geometry prior and the science-line mask.

History::

    260901  ksl  Coding begun.
    261007  ksl  Updated for lvmsky branch skydecomp-telluric-corrected-lines:
        write LSF_SCI/LSF_SKY_NEAR/LSF_SKY_FAR and the lower-case META
        columns its telluric fit models read (airmasses, date_obs,
        pwv_med, Near/Far pointings and labels); fill NaN LSF pixels at
        the ends of the wavelength range, which otherwise make
        decompose_parallel.py reject the row.  The stack is built by
        build_stack(), which PredictSky.py also uses.

'''

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table


def fill_lsf(lsf):
    '''
    Replace non-finite or non-positive LSF pixels by linear interpolation
    between the nearest valid pixels of the same row (held constant past
    the first/last valid pixel).

    SummarizeCframe.py's LSF is NaN at the ends of the wavelength range
    (3600-3604 A, 9800 A) in many rows, inherited from CFrame fibers whose
    own LSF is NaN there; decompose_parallel.py rejects a row with a bad
    LSF pixel at the array edge.  The LSF varies smoothly with wavelength,
    so the filled values are the obvious ones.

    Parameters
    ----------
    lsf : numpy.ndarray
        (n_rows, n_wave) LSF FWHM array.

    Returns
    -------
    tuple
        (filled float32 array, number of rows changed, number of pixels
        filled).  Rows with no valid pixel at all are left unchanged.
    '''
    lsf = np.array(lsf, dtype=np.float32)
    good = np.isfinite(lsf) & (lsf > 0)
    pix = np.arange(lsf.shape[1])
    n_rows = n_pix = 0
    for i in np.flatnonzero(~good.all(axis=1)):
        if not good[i].any():
            continue
        bad = ~good[i]
        lsf[i, bad] = np.interp(pix[bad], pix[good[i]], lsf[i, good[i]])
        n_rows += 1
        n_pix += int(bad.sum())
    return lsf, n_rows, n_pix


def build_stack(wave, flux_sci, sky_east, sky_west, lsf, drp,
                wave_hdr=None, primary_hdr=None, verbose=True):
    '''
    Build decompose_parallel.py's input stack, in memory, from XCframe
    arrays.  Used by convert() for training and by PredictSky.py for one
    exposure, so both decompose exactly the same input.

    Parameters
    ----------
    wave : numpy.ndarray
        Wavelength grid.
    flux_sci, sky_east, sky_west, lsf : numpy.ndarray
        (n_rows, n_wave) XCframe FLUX, SKY_EAST, SKY_WEST and LSF rows.
    drp : astropy.table.Table
        The matching DRP_ALL rows (a copy is modified, not drp itself).
    wave_hdr, primary_hdr : astropy.io.fits.Header, optional
        Headers to carry over.
    verbose : bool
        Print the PWV and LSF-fill messages.

    Returns
    -------
    astropy.io.fits.HDUList
        WAVE, FLUX_SCI, FLUX_SKY_NEAR, FLUX_SKY_FAR, LSF_SCI, LSF_SKY_NEAR,
        LSF_SKY_FAR and META.
    '''
    drp = Table(drp, copy=True)
    flux_sci = np.atleast_2d(flux_sci)
    sky_east = np.atleast_2d(sky_east)
    sky_west = np.atleast_2d(sky_west)

    near = np.char.strip(np.asarray(drp['Near'], str))
    far = np.char.strip(np.asarray(drp['Far'], str))
    valid = {'SKY_EAST', 'SKY_WEST'}
    bad_near = set(np.unique(near)) - valid
    bad_far = set(np.unique(far)) - valid
    if bad_near or bad_far:
        raise ValueError(f'Unexpected Near/Far labels: near={bad_near}, far={bad_far}')

    is_near_east = (near == 'SKY_EAST')
    is_far_east = (far == 'SKY_EAST')
    flux_sky_near = np.where(is_near_east[:, None], sky_east, sky_west)
    flux_sky_far = np.where(is_far_east[:, None], sky_east, sky_west)

    skye_ra = np.asarray(drp['skye_ra'], float)
    skye_dec = np.asarray(drp['skye_dec'], float)
    skyw_ra = np.asarray(drp['skyw_ra'], float)
    skyw_dec = np.asarray(drp['skyw_dec'], float)
    drp['sky_near_ra'] = np.where(is_near_east, skye_ra, skyw_ra)
    drp['sky_near_dec'] = np.where(is_near_east, skye_dec, skyw_dec)
    drp['sky_far_ra'] = np.where(is_far_east, skye_ra, skyw_ra)
    drp['sky_far_dec'] = np.where(is_far_east, skye_dec, skyw_dec)
    drp['sky_near_label'] = np.where(is_near_east, 'SkyE', 'SkyW')
    drp['sky_far_label'] = np.where(is_far_east, 'SkyE', 'SkyW')
    drp['sci_airmass'] = np.asarray(drp['sci_amass'], float)
    drp['skye_airmass'] = np.asarray(drp['skye_amass'], float)
    drp['skyw_airmass'] = np.asarray(drp['skyw_amass'], float)
    drp['date_obs'] = np.char.strip(np.asarray(drp['obstime'], str))

    if 'pwv_med' in drp.colnames:
        pwv = np.asarray(drp['pwv_med'], float)
    else:
        pwv = np.full(len(drp), np.nan)
        drp['pwv_med'] = pwv
    n_bad = int(np.sum(~(np.isfinite(pwv) & (pwv > 0))))
    if n_bad and verbose:
        print(f'Warning: {n_bad} of {len(drp)} rows have no valid pwv_med; '
              f'decompose_parallel.py will use the DRP default of 15 mm for them '
              f'(run SelectXCF.py with -hdr to add PWV)')

    lsf_sci, n_rows, n_pix = fill_lsf(np.atleast_2d(lsf))
    if n_rows and verbose:
        print(f'Filled {n_pix} bad LSF pixels in {n_rows} of {len(drp)} rows '
              f'(mostly the ends of the wavelength range)')

    primary_hdr = fits.Header() if primary_hdr is None else primary_hdr.copy()
    primary_hdr['LSFSKY'] = ('LSF_SCI', 'Sky-telescope LSFs approximated by LSF_SCI')

    return fits.HDUList([
        fits.PrimaryHDU(header=primary_hdr),
        fits.ImageHDU(data=wave, header=wave_hdr, name='WAVE'),
        fits.ImageHDU(data=flux_sci.astype(np.float32), name='FLUX_SCI'),
        fits.ImageHDU(data=flux_sky_near.astype(np.float32), name='FLUX_SKY_NEAR'),
        fits.ImageHDU(data=flux_sky_far.astype(np.float32), name='FLUX_SKY_FAR'),
        fits.ImageHDU(data=lsf_sci, name='LSF_SCI'),
        fits.ImageHDU(data=lsf_sci, name='LSF_SKY_NEAR'),
        fits.ImageHDU(data=lsf_sci, name='LSF_SKY_FAR'),
        fits.BinTableHDU(drp, name='META'),
    ])


def convert(fits_file, outpath):
    '''
    Write fits_file (XCframe layout) as a decompose_parallel.py input
    stack at outpath.  See the module docstring for the layout.
    '''
    with fits.open(fits_file, memmap=True) as hdul:
        primary_hdr = hdul['PRIMARY'].header.copy()
        primary_hdr['SRCFILE'] = (str(Path(fits_file).name), 'Input XCframe-layout file')
        stack = build_stack(
            hdul['WAVE'].data, hdul['FLUX'].data, hdul['SKY_EAST'].data,
            hdul['SKY_WEST'].data, hdul['LSF'].data, Table(hdul['DRP_ALL'].data),
            wave_hdr=hdul['WAVE'].header.copy(), primary_hdr=primary_hdr)
        stack.writeto(outpath, overwrite=True)
    print(f'Wrote {len(stack["META"].data)} rows to {outpath} '
          f'(FLUX_SCI/FLUX_SKY_NEAR/FLUX_SKY_FAR, LSF_*, META)')


def main():
    p = argparse.ArgumentParser(
        description=('Reformat an XCframe-layout FITS file into the layout '
                     'decompose_parallel.py expects.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('fits_file', help='XCframe-layout FITS file (e.g. SelectXCF.py output)')
    p.add_argument('output_file', nargs='?', default=None,
                   help='Output FITS path (default: <stem>_decomp_input.fits)')
    args = p.parse_args()

    outpath = args.output_file
    if outpath is None:
        stem = Path(args.fits_file).stem
        outpath = str(Path(args.fits_file).parent / f'{stem}_decomp_input.fits')

    convert(args.fits_file, outpath)


if __name__ == '__main__':
    main()
