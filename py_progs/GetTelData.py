#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Extract the FLUX/IVAR/MASK/LSF/SLITMAP data for ONE requested telescope
    (Sci, SkyE, or SkyW -- or whichever of SkyE/SkyW is nearer/farther from
    the science field) out of a single lvmCFrame or lvmSFrame file, and
    either return it in memory (for notebook use) or write it to a
    standalone output FITS file (command-line use).

Command line usage (if any):

    usage: GetTelData.py telescope filename [-out FILE] [-fiberid N]

    where

    telescope   Sci | SkyE | SkyW | Near | Far (case-insensitive). Near/Far
                are resolved per-exposure from the header's
                SCI_SKYE_SEP/SCI_SKYW_SEP separations.

    filename    an lvmCFrame or lvmSFrame FITS file. An SFrame path is
                transparently swapped for the matching CFrame, since the
                per-telescope FLUX/IVAR/MASK/LSF data live there.

    Options::

        -h            print this help and exit
        -out FILE     output FITS filename (default:
                      <stem>_<TELESCOPE>.fits)
        -fiberid N    restrict to a single fiberid (must belong to the
                      requested telescope)
        -sky          also include the SKY_EAST/SKY_WEST broadcast arrays
                      (Sci telescope only; see Notes)

Description:

    A raw CFrame carries one shared FLUX/IVAR/MASK/LSF/SLITMAP set covering
    every fiber -- Sci, SkyE, and SkyW alike, distinguished only by the
    SLITMAP 'telescope' column. Every consumer needing just one telescope's
    data (sky-subtraction methods, quality-assurance plots, ad hoc notebook
    checks) has so far re-opened the file and re-derived that selection
    independently. This script is the single shared place that does it.

    get_tel_data(filename, telescope) returns FLUX/IVAR/MASK/LSF/SLITMAP
    for exactly the telescope requested -- asking for 'Sci' never returns
    SkyE/SkyW data, and only the requested fibers' rows are read out of the
    FITS arrays. Near/Far are resolved from the DRP's own
    SCI_SKYE_SEP/SCI_SKYW_SEP header keywords rather than recomputed.

Primary routines:

    get_tel_data     read one telescope's data from one CFrame/SFrame
    write_tel_data   write that data to a standalone FITS file

Notes:

    This is a standalone script (no imports from lvm_gaussfit.py or
    SummarizeCframe.py) -- the fiber-selection filter and near/far lookup
    are small, self-contained copies of the logic used elsewhere, not a
    shared import, so this can be used independently while those other
    routines are still on their existing, separate implementations.

    The SKY_EAST/SKY_WEST broadcast extensions (the sky-model spectra
    repeated across every Sci fiber row) are NOT returned by default --
    they belong to a Sci-fiber selection conceptually, not to SkyE/SkyW's
    own fibers, so returning them unconditionally would mean asking for
    one telescope could silently hand back another's data. Pass
    include_sky=True (Sci selection only) to opt in.

History::

    260913 ksl Coding begun. get_tel_data()/write_tel_data() extract one
        telescope's (Sci/SkyE/SkyW/Near/Far) FLUX/IVAR/MASK/LSF/SLITMAP
        from a single lvmCFrame/lvmSFrame, replacing the copy-pasted
        per-script fiber selection used across Prep4SkyCorr.py,
        kslmap.py/quick_map.py/line_map.py, and others. Near/Far are
        resolved from the DRP's own SCI_SKYE_SEP/SCI_SKYW_SEP header
        keywords rather than recomputed. Asking for one telescope never
        returns another's data.
    260913 ksl Added include_sky (-sky) to get_tel_data(): opt-in
        SKY_EAST/SKY_WEST/SKY_EAST_IVAR/SKY_WEST_IVAR broadcast arrays
        for the Sci selection, needed by SummarizeCframe.py's
        get_med_spec()/get_fiber_spec() and QualCFrame.py/QualSFrame.py.
'''

import sys
import re
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before a line consisting of "History:" (or
    "History::"/"Version History" -- whitespace/colon-insensitive), so
    -h stays short even as that section grows -- without hand-
    duplicating the Synopsis/Options text in a second string.  Anchored
    to a whole line (not a bare substring search) so it can't misfire on
    "History:" appearing mid-sentence, and returns doc unchanged if no
    such line is present.
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


_USAGE = _usage_from_doc(__doc__)

_TEL_NAME = {'SCI': 'Sci', 'SKYE': 'SkyE', 'SKYW': 'SkyW'}


def _select_fibers(xtab, telescope):
    '''
    Good (fibstatus==0) fibers of one telescope from a SLITMAP table.
    '''
    ztab = xtab[xtab['fibstatus'] == 0]
    return ztab[ztab['telescope'] == telescope]


def get_tel_data(filename, telescope='Sci', fiberid=None, include_sky=False):
    '''
    Read one raw lvmCFrame (or lvmSFrame, transparently swapped to the
    matching CFrame) and return FLUX/IVAR/MASK/LSF/SLITMAP for exactly
    the requested telescope -- nothing else, unless include_sky is set.

    Parameters:
        filename: str
            Path to an lvmCFrame or lvmSFrame file.
        telescope: str
            'Sci' | 'SkyE' | 'SkyW' | 'Near' | 'Far' (case-insensitive).
            'Near'/'Far' resolve to SkyE or SkyW per-exposure from the
            header's SCI_SKYE_SEP/SCI_SKYW_SEP separations.
        fiberid: int or None
            If given, restrict to this one fiberid (must belong to the
            requested telescope, or the result has zero rows).
        include_sky: bool
            If True, also return the SKY_EAST/SKY_WEST broadcast arrays
            (skye_flux, skye_ivar, skyw_flux, skyw_ivar) -- the sky-model
            spectra the DRP repeats across every Sci fiber row. Only
            valid when telescope resolves to 'Sci'; otherwise an error is
            printed and None is returned, since those arrays don't mean
            anything relative to SkyE/SkyW's own fibers.

    Returns::

        dict, or None if the file/telescope/fiberid could not be resolved:
            wave           : 1-D wavelength array
            flux, ivar, mask, lsf : 2-D arrays (n_fibers, n_wave) for the
                             requested fibers only
            slitmap        : Table of the matching SLITMAP rows,
                             row-aligned with the arrays above
            telescope      : 'Sci'/'SkyE'/'SkyW' actually selected
            requested      : the telescope argument as given
            near, far      : 'SkyE'/'SkyW', or None if the separation
                             header keywords are not present
            sep_e, sep_w   : float (deg), or None
            header         : the file's PRIMARY header (for provenance)
            skye_flux, skye_ivar, skyw_flux, skyw_ivar : 2-D arrays,
                             only present if include_sky was True
    '''
    if filename.count('SFrame'):
        filename = filename.replace('SFrame', 'CFrame')

    try:
        x = fits.open(filename)
    except Exception:
        print('get_tel_data: could not open %s' % filename)
        return None

    hdr = x['PRIMARY'].header
    try:
        sep_e = float(hdr['SKY SCI_SKYE_SEP'])
        sep_w = float(hdr['SKY SCI_SKYW_SEP'])
    except KeyError:
        sep_e = sep_w = None

    if sep_e is not None:
        near = 'SkyE' if sep_e <= sep_w else 'SkyW'
        far = 'SkyW' if sep_e <= sep_w else 'SkyE'
    else:
        near = far = None

    tel_key = telescope.upper()
    if tel_key in ('NEAR', 'FAR'):
        resolved = near if tel_key == 'NEAR' else far
        if resolved is None:
            print('get_tel_data: SCI_SKYE_SEP/SCI_SKYW_SEP not found '
                  'in header; cannot resolve %s' % telescope)
            return None
    elif tel_key in _TEL_NAME:
        resolved = _TEL_NAME[tel_key]
    else:
        print('get_tel_data: telescope must be one of '
              'Sci, SkyE, SkyW, Near, Far (got %s)' % telescope)
        return None

    if include_sky and resolved != 'Sci':
        print('get_tel_data: include_sky is only valid for the Sci '
              'selection (got telescope=%s, resolved to %s)'
              % (telescope, resolved))
        return None

    xtab = Table(x['SLITMAP'].data)
    rows = _select_fibers(xtab, resolved)
    if fiberid is not None:
        rows = rows[rows['fiberid'] == fiberid]

    if len(rows) == 0:
        print('get_tel_data: no matching fibers for %s in %s'
              % (telescope, filename))
        return None

    idx = np.asarray(rows['fiberid']) - 1

    result = {
        'wave': x['WAVE'].data,
        'flux': x['FLUX'].data[idx],
        'ivar': x['IVAR'].data[idx],
        'mask': x['MASK'].data[idx],
        'lsf': x['LSF'].data[idx],
        'slitmap': rows,
        'telescope': resolved,
        'requested': telescope,
        'near': near,
        'far': far,
        'sep_e': sep_e,
        'sep_w': sep_w,
        'header': hdr,
    }

    if include_sky:
        result['skye_flux'] = x['SKY_EAST'].data[idx]
        result['skye_ivar'] = x['SKY_EAST_IVAR'].data[idx]
        result['skyw_flux'] = x['SKY_WEST'].data[idx]
        result['skyw_ivar'] = x['SKY_WEST_IVAR'].data[idx]

    return result


def write_tel_data(data, filename, outfile=None):
    '''
    Write get_tel_data()'s output to a standalone FITS file
    (WAVE/FLUX/IVAR/MASK/LSF/SLITMAP extensions).

    Parameters:
        data: dict
            Return value of get_tel_data().
        filename: str
            Original input filename, used for the default output name.
        outfile: str or None
            Output FITS filename (default: <stem>_<TELESCOPE>.fits).

    Returns:
        outfile (str).
    '''
    if not outfile:
        outfile = '%s_%s.fits' % (Path(filename).stem, data['telescope'].upper())

    hdr = fits.Header()
    hdr['SRCFILE'] = (Path(filename).name, 'source CFrame/SFrame file')
    hdr['TELSEL'] = (data['telescope'], 'telescope selection written to this file')
    for key in ('MJD', 'EXPOSURE'):
        if key in data['header']:
            hdr[key] = data['header'][key]
    if data['requested'].upper() in ('NEAR', 'FAR'):
        hdr['NEARFAR'] = (data['requested'].upper(),
                           'Near/Far resolved to %s' % data['telescope'])
    if data['sep_e'] is not None:
        hdr['SEP_E'] = (data['sep_e'], 'science-SkyE separation (deg)')
        hdr['SEP_W'] = (data['sep_w'], 'science-SkyW separation (deg)')

    hdus = fits.HDUList([
        fits.PrimaryHDU(header=hdr),
        fits.ImageHDU(data=data['wave'], name='WAVE'),
        fits.ImageHDU(data=data['flux'], name='FLUX'),
        fits.ImageHDU(data=data['ivar'], name='IVAR'),
        fits.ImageHDU(data=data['mask'], name='MASK'),
        fits.ImageHDU(data=data['lsf'], name='LSF'),
        fits.BinTableHDU(data=data['slitmap'], name='SLITMAP'),
    ])
    if 'skye_flux' in data:
        hdus.append(fits.ImageHDU(data=data['skye_flux'], name='SKY_EAST'))
        hdus.append(fits.ImageHDU(data=data['skye_ivar'], name='SKY_EAST_IVAR'))
        hdus.append(fits.ImageHDU(data=data['skyw_flux'], name='SKY_WEST'))
        hdus.append(fits.ImageHDU(data=data['skyw_ivar'], name='SKY_WEST_IVAR'))
    hdus.writeto(outfile, overwrite=True)
    print('Wrote %s (%d fibers)' % (outfile, len(data['slitmap'])))
    return outfile


def steer(argv):
    telescope = ''
    filename = ''
    outfile = ''
    fiberid = None
    include_sky = False

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-out':
            i += 1
            outfile = argv[i]
        elif argv[i] == '-fiberid':
            i += 1
            fiberid = int(argv[i])
        elif argv[i] == '-sky':
            include_sky = True
        elif argv[i][0] == '-':
            print('Error: cannot parse command line:', argv)
            return
        elif telescope == '':
            telescope = argv[i]
        elif filename == '':
            filename = argv[i]
        else:
            print('Error: cannot parse command line:', argv)
            return
        i += 1

    if not telescope or not filename:
        print(_USAGE)
        return

    data = get_tel_data(filename, telescope=telescope, fiberid=fiberid,
                         include_sky=include_sky)
    if data is None:
        return

    write_tel_data(data, filename, outfile=outfile if outfile else None)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
