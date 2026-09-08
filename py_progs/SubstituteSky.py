#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Replace one sky telescope's fiber data in an lvmCFrame file with
    the corresponding fiber data from a different lvmCFrame, so the
    DRP's sky-subtraction step can be rerun using a better sky
    measurement (e.g. when SkyE or SkyW was pointed too close to a
    bright Moon).

Command line usage (if any)::

    usage: SubstituteSky.py [-h] [-o outfile] target_cframe target_tel
                             source_cframe source_tel

    where
        -h             prints this documentation and exits
        -o outfile     output filename (default: target_cframe's
                       basename -- directory stripped -- with
                       '.sky_subst' inserted before the extension,
                       written to the current directory). An
                       existing outfile is overwritten, with a
                       warning printed first.

        target_cframe  the lvmCFrame file to be corrected
        target_tel     SkyE or SkyW (case-insensitive) -- which
                       telescope's fibers in target_cframe are to
                       be replaced
        source_cframe  the lvmCFrame file to draw the replacement
                       data from
        source_tel     SkyE or SkyW (case-insensitive) -- which
                       telescope's fibers in source_cframe supply
                       the replacement data

    Example (exposure 14964's SkyW pointed 4.5 deg from the Moon;
    substitute in SkyW data from a clean exposure 14960):

        SubstituteSky.py lvmCFrame-00014964.fits SkyW \\
                         lvmCFrame-00014960.fits SkyW

Description:

    An lvmCFrame's SLITMAP assigns every fiber to a telescope (Sci,
    SkyE, SkyW, Spec) -- fixed by fiber-plugging hardware, identical
    fiberid-for-fiberid across all exposures. The DRP's sky-
    subtraction routine (skyMethod.quick_sky_subtraction) builds its
    sky spectrum from the raw FLUX/IVAR at the fibers tagged
    SkyE/SkyW in the CFrame -- not the extrapolated SKY_EAST/SKY_WEST
    extensions, which the current production method ignores. So
    fixing a contaminated sky telescope means replacing the FLUX,
    IVAR, MASK, and LSF rows for that telescope's fibers.

    This routine copies those rows (matched by fiberid) from
    source_cframe/source_tel into a copy of target_cframe/target_tel,
    along with every PRIMARY header keyword tied to the source
    telescope -- pointing, altitude, airmass, guider frames,
    sky-field name, heliocentric velocity, moon/shadow geometry,
    ecliptic coordinates, etc. (but not the SKYEW/SKYWW combination
    weights, a joint SkyE+SkyW property recomputed elsewhere). It
    adds SKY SUBST_* provenance keywords, then writes the result to
    outfile. target_cframe and source_cframe are never modified.

    Since SkyE/SkyW fiber assignment is fixed hardware, requesting
    different telescopes on the two sides (target_tel != source_tel)
    fails with a clear error -- there is no physically meaningful
    fiber-by-fiber correspondence between them.

Notes:

    Only SKYWRA/SKYWDEC/SKYWALT (or the SkyE equivalents) are
    actually read downstream by the DRP; the rest of the header
    block is for provenance/QA.

    SKY SCI_SKYW_SEP (or SkyE) is relative to the *target's* science
    pointing, so it is recomputed from the target's SCIRA/SCIDEC and
    the newly-copied sky position (lvmdrp.core.sky.ang_distance)
    rather than copied as-is. Other geometry keys (moon/shadow
    separation, ecliptic coords) describe the conditions the
    *substituted* data was actually observed under, so those are
    copied unchanged.

    Rerun lvmdrp.functions.skyMethod.quick_sky_subtraction on
    outfile to produce a corrected lvmSFrame.

History::

    260908 ksl Coding begun. Substitutes FLUX/IVAR/MASK/LSF fiber rows
        (matched by fiberid) and the associated PRIMARY header block
        between lvmCFrame files, for recovering an exposure whose sky
        telescope was pointed too close to the Moon. SKY SCI_<tel>_SEP
        is recomputed (via lvmdrp.core.sky.ang_distance) rather than
        copied, since it is relative to the target's own science
        pointing. Telescope names are case-insensitive; default output
        filename strips the input directory and writes to the current
        one; an existing outfile is overwritten with a warning, no -f
        flag needed.

'''

import os
import re
import sys
import numpy as np
from astropy.io import fits
from astropy.table import Table
from lvmdrp.core.sky import ang_distance


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

# fiber/RSS extensions that describe a single fiber's spectrum and must
# travel together when a fiber's data is substituted from another exposure
EXTENSIONS_TO_COPY = ('FLUX', 'IVAR', 'MASK', 'LSF')

# maps a telescope name to the 4-character token used throughout the
# lvmCFrame PRIMARY header to mark keywords specific to that telescope
TEL_TOKENS = {'SkyE': 'SKYE', 'SkyW': 'SKYW'}

# SkyE/SkyW combination weights: a joint property of both telescopes
# (recomputed by skyMethod.combine_skies), not a per-telescope one
WEIGHT_KEYS = {'SKYEW', 'SKYWW'}


def normalize_tel(tel):
    '''
    Map a user-supplied telescope name (case-insensitive) to the
    canonical 'SkyE'/'SkyW' form used in TEL_TOKENS and the SLITMAP
    'telescope' column. Raises ValueError for anything else.
    '''
    key = tel.strip().lower()
    for canonical in TEL_TOKENS:
        if key == canonical.lower():
            return canonical
    raise ValueError(f"invalid telescope '{tel}': must be one of {list(TEL_TOKENS)}")


def substitute_sky(target_cframe, target_tel, source_cframe, source_tel, outfile=None):
    '''
    Replace the FLUX/IVAR/MASK/LSF fiber rows for target_tel in
    target_cframe with the corresponding rows for source_tel from
    source_cframe, copy the associated telescope header block, and
    write the result to outfile.

    Parameters
    ----------
    target_cframe: str
        lvmCFrame file to be corrected
    target_tel: str
        'SkyE' or 'SkyW' -- telescope whose fibers get replaced
    source_cframe: str
        lvmCFrame file to draw replacement data from
    source_tel: str
        'SkyE' or 'SkyW' -- telescope in source_cframe supplying the data
    outfile: str, optional
        output filename; default is target_cframe with '.sky_subst'
        inserted before the extension

    Returns
    -------
    str
        the path written
    '''

    target_tel = normalize_tel(target_tel)
    source_tel = normalize_tel(source_tel)

    with fits.open(target_cframe, memmap=False) as ht, fits.open(source_cframe, memmap=False) as hs:

        wave_t = ht['WAVE'].data
        wave_s = hs['WAVE'].data
        if wave_t.shape != wave_s.shape or not np.allclose(wave_t, wave_s):
            raise ValueError(
                f"'{target_cframe}' and '{source_cframe}' do not share the same WAVE grid; "
                "cannot substitute fiber rows directly"
            )

        slit_t = Table(ht['SLITMAP'].data)
        slit_s = Table(hs['SLITMAP'].data)

        fiberids = np.asarray(slit_t['fiberid'][slit_t['telescope'] == target_tel])
        if fiberids.size == 0:
            raise ValueError(f"no fibers tagged '{target_tel}' found in target SLITMAP")

        rows = fiberids - 1

        src_tel_at_rows = np.asarray(slit_s['telescope'])[rows]
        bad = src_tel_at_rows != source_tel
        if bad.any():
            example = fiberids[bad][0]
            example_tel = src_tel_at_rows[bad][0]
            raise ValueError(
                f"{bad.sum()} of the target's {target_tel} fibers are not tagged '{source_tel}' "
                f"in the source SLITMAP (e.g. fiberid {example} is '{example_tel}' there). "
                "SkyE/SkyW fiber assignment is fixed by hardware and identical across exposures, "
                "so substituting between different telescopes fiber-by-fiber is not physically "
                "meaningful -- did you mean the same telescope on both sides?"
            )

        for ext in EXTENSIONS_TO_COPY:
            ht[ext].data[rows] = hs[ext].data[rows]

        src_token = TEL_TOKENS[source_tel]
        dst_token = TEL_TOKENS[target_tel]
        hdr_t = ht['PRIMARY'].header
        hdr_s = hs['PRIMARY'].header

        for card in hdr_s.cards:
            kw = card.keyword
            if src_token not in kw or kw in WEIGHT_KEYS or kw.startswith('SKY SUBST'):
                continue
            new_kw = kw.replace(src_token, dst_token)
            hdr_t[new_kw] = (card.value, card.comment)

        # SKY SCI_{tel}_SEP is the angular separation between the science
        # field and this sky telescope -- a relation to the target's own
        # science pointing, not an intrinsic property of the source
        # exposure, so the block copy above leaves it wrong (it still holds
        # the source's own SCI-to-sky separation). Recompute it from the
        # target's real SCIRA/SCIDEC and the just-copied sky RA/Dec.
        sep_kw = f'SKY SCI_{dst_token}_SEP'
        if sep_kw in hdr_t:
            new_sep = ang_distance(hdr_t['SCIRA'], hdr_t['SCIDEC'],
                                    hdr_t[f'{dst_token}RA'], hdr_t[f'{dst_token}DEC'])
            hdr_t[sep_kw] = (round(float(new_sep), 4), hdr_t.comments[sep_kw])

        hdr_t['HIERARCH SKY SUBST_TEL'] = (target_tel, 'telescope whose fiber data was substituted')
        hdr_t['HIERARCH SKY SUBST_SRC'] = (os.path.basename(source_cframe), 'source CFrame for substituted sky')
        hdr_t['HIERARCH SKY SUBST_SRCTEL'] = (source_tel, 'telescope in source CFrame supplying the data')
        hdr_t['HIERARCH SKY SUBST_SRCEXP'] = (hdr_s.get('EXPOSURE', -1), 'source exposure number')
        hdr_t['HIERARCH SKY SUBST_NFIB'] = (int(rows.size), 'number of fibers substituted')

        if outfile is None:
            base, ext_ = os.path.splitext(os.path.basename(target_cframe))
            outfile = f"{base}.sky_subst{ext_}"

        if os.path.exists(outfile):
            print(f"Warning: overwriting existing file '{outfile}'")

        ht.writeto(outfile, overwrite=True)

    return outfile


def steer(argv):
    '''
    Just a steering routine
    '''

    outfile = None
    positional = []

    i = 1
    while i < len(argv):
        if argv[i][0:2] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-o':
            i += 1
            outfile = argv[i]
        elif argv[i][0] == '-':
            print('Error: Unknown optional parameter; improperly formatted command line: ', argv)
            return
        else:
            positional.append(argv[i])
        i += 1

    if len(positional) != 4:
        print('Error: expected 4 positional arguments, got', len(positional), positional)
        print(_USAGE)
        return

    target_cframe, target_tel, source_cframe, source_tel = positional

    try:
        outname = substitute_sky(target_cframe, target_tel, source_cframe, source_tel,
                                  outfile=outfile)
    except (ValueError, OSError) as err:
        print('Error:', err)
        return

    print('Wrote', outname)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
