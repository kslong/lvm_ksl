#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Reformat an LVM XCframe-layout FITS file (WAVE/FLUX/SKY_EAST/SKY_WEST/
    LSF/DRP_ALL) into the layout lvmsky's decompose_parallel.py expects for
    --fit-model lsf-surface-iterative-split-zodi (WAVE/FLUX_SCI/
    FLUX_SKY_NEAR/FLUX_SKY_FAR/META).

Command line usage (if any):

    usage: ConvertForDecompose.py [-h] fits_file [output_file]

    where

    fits_file       is the path to an XCframe-layout FITS file (e.g. the
                    output of SelectXCF.py).

    output_file     output FITS path (default: <stem>_decomp_input.fits).

Description:

    decompose_parallel.py's worker (for --fit-model
    lsf-surface-iterative-split-zodi) only reads WAVE, FLUX_SCI,
    FLUX_SKY_NEAR, and FLUX_SKY_FAR from the per-row fit, and copies META
    through verbatim afterwards (extract_meta_and_coef_products) -- it does
    not need any LSF_* extensions for this fit model (that only applies to
    the separate moon-zodi-lsf-surface-iterative fit model, not used here).

    This script therefore:

      1. Copies WAVE unchanged.
      2. Renames FLUX -> FLUX_SCI.
      3. Builds FLUX_SKY_NEAR / FLUX_SKY_FAR by selecting, row by row,
         SKY_EAST or SKY_WEST according to that row's own DRP_ALL Near/Far
         label (confirmed values 'SKY_EAST'/'SKY_WEST', matching the input
         file's own extension names) -- this is a per-row choice, not a
         fixed rename, since which physical telescope (SkyE or SkyW) is
         nearer varies exposure to exposure.
      4. Copies DRP_ALL through as the META extension (lvmsky's
         extract_meta_and_coef_products does a byte-for-byte HDU copy, so
         whatever columns are here pass straight through to the training
         notebook), plus adds the derived columns mlp_predictor.data
         actually requires and DRP_ALL doesn't have under those names:
         SKY_NEAR_RA/DEC, SKY_FAR_RA/DEC (per-row Near/Far selection, same
         idea as FLUX_SKY_NEAR/FAR), and SKY_NEAR_LABEL/SKY_FAR_LABEL
         (values 'SKYE'/'SKYW' -- confirmed exact strings _resolve_context_
         feature's 'ew' branch matches against, not 'SKY_EAST'/'SKY_WEST').
         Read _build_context_matrix/_pointing_ra_dec_columns/
         _resolve_context_feature in mlp_predictor/data.py this session to
         confirm these are the only columns missing under expected names;
         everything else (SCI_RA/DEC, OBSTIME, moon_phase) resolves via
         DRP_ALL's existing lowercase columns through meta_upper's
         case-insensitive lookup, and alt/az/airmass/moon_*/sun_*/vanrhijn_*
         are computed internally from RA/Dec/MJD, not read from META.

Notes::

    LSF is intentionally dropped here -- not required by
    lsf-surface-iterative-split-zodi's worker (it uses a single global
    --lsf-sigma scalar instead of per-row LSF).

History::

    260901  ksl  Coding begun.

'''

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table


def convert(fits_file, outpath):
    with fits.open(fits_file, memmap=True) as hdul:
        primary_hdr = hdul['PRIMARY'].header.copy()
        wave = hdul['WAVE'].data
        wave_hdr = hdul['WAVE'].header.copy()
        flux_sci = hdul['FLUX'].data
        sky_east = hdul['SKY_EAST'].data
        sky_west = hdul['SKY_WEST'].data
        drp = Table(hdul['DRP_ALL'].data)

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
    drp['SKY_NEAR_RA'] = np.where(is_near_east, skye_ra, skyw_ra)
    drp['SKY_NEAR_DEC'] = np.where(is_near_east, skye_dec, skyw_dec)
    drp['SKY_FAR_RA'] = np.where(is_far_east, skye_ra, skyw_ra)
    drp['SKY_FAR_DEC'] = np.where(is_far_east, skye_dec, skyw_dec)
    drp['SKY_NEAR_LABEL'] = np.where(is_near_east, 'SKYE', 'SKYW')
    drp['SKY_FAR_LABEL'] = np.where(is_far_east, 'SKYE', 'SKYW')

    primary_hdr['SRCFILE'] = (str(Path(fits_file).name), 'Input XCframe-layout file')

    hdul_out = fits.HDUList([
        fits.PrimaryHDU(header=primary_hdr),
        fits.ImageHDU(data=wave, header=wave_hdr, name='WAVE'),
        fits.ImageHDU(data=flux_sci.astype(np.float32), name='FLUX_SCI'),
        fits.ImageHDU(data=flux_sky_near.astype(np.float32), name='FLUX_SKY_NEAR'),
        fits.ImageHDU(data=flux_sky_far.astype(np.float32), name='FLUX_SKY_FAR'),
        fits.BinTableHDU(drp, name='META'),
    ])
    hdul_out.writeto(outpath, overwrite=True)
    print(f'Wrote {len(drp)} rows to {outpath} '
          f'(FLUX_SCI/FLUX_SKY_NEAR/FLUX_SKY_FAR/META, no LSF)')


def main():
    p = argparse.ArgumentParser(
        description=('Reformat an XCframe-layout FITS file into the layout '
                     'decompose_parallel.py expects for '
                     '--fit-model lsf-surface-iterative-split-zodi.'),
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
