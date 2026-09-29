#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Replace one sky telescope's data in an lvmCFrame file -- its fiber
    rows and its SKY_EAST/SKY_WEST sky-model extension -- with data
    from the same or the other sky telescope, in the same or a
    different lvmCFrame, producing a new, self-consistent lvmCFrame on
    which any sky-subtraction method can be run with a better sky
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

    Examples (exposure 14964's SkyW pointed 4.5 deg from the Moon):

        # use 14964's own (clean) SkyE in place of its SkyW
        SubstituteSky.py lvmCFrame-00014964.fits SkyW \\
                         lvmCFrame-00014964.fits SkyE

        # or substitute SkyW data from a clean exposure 14960
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

    Many non-DRP tools (SummarizeCframe.py, GetTelData.py,
    lvm_skyfit.py, Prep4SkyCorr.py, QualCFrame.py, ...) instead read
    the SKY_EAST/SKY_WEST extensions, so both must be replaced for the
    output to be usable by every method.

    This routine writes a copy of target_cframe in which, for
    target_tel:

    - the FLUX/IVAR/MASK/LSF rows of its fibers come from
      source_tel's fibers in source_cframe. For the same telescope
      they are matched by fiberid (SkyE/SkyW fiber assignment is
      fixed hardware, identical across exposures). For different
      telescopes there is no such correspondence, so a warning is
      printed and the target rows are filled in order from the source
      telescope's good (fibstatus==0) fibers, reusing them cyclically
      if there are fewer; since the DRP and the XCframe summaries
      average over a telescope's fibers, the ordering does not matter.
      Each row's LSF travels with its spectrum.
    - SKY_EAST/SKY_WEST (whichever is target_tel's) and its _IVAR are
      replaced by source_tel's in source_cframe.
    - every PRIMARY header keyword tied to the telescope -- pointing,
      altitude, airmass, guider frames, sky-field name, heliocentric
      velocity, moon/shadow geometry, ecliptic coordinates, etc. --
      is copied from source_tel.
    - the SKYEW/SKYWW combination weights are recomputed from the
      updated pointings with skyMethod.combine_skies's own formula
      (inverse angular distance to the science field, normalized).

    SKY SUBST_* provenance keywords are added. target_cframe and
    source_cframe are never modified. The target's SLITMAP is left
    unchanged.

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

    When SkyE is substituted into SkyW of the same exposure (or vice
    versa), both telescopes then sit at the same position, so the
    DRP's near/far choice is a tie and both its line and continuum
    sky come from the substituted telescope; SKYEW/SKYWW become
    0.5/0.5.

    Rerun lvmdrp.functions.skyMethod.quick_sky_subtraction on
    outfile (e.g. with RunSky.py) to produce a corrected lvmSFrame.

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
    260929 ksl Now produces a fully self-consistent CFrame for use with
        any sky-subtraction method, not just the DRP's: also replaces
        the target's SKY_EAST/SKY_WEST (+_IVAR) extension and
        recomputes SKYEW/SKYWW. Different telescopes on the two sides
        (e.g. SkyE -> SkyW within one exposure) now print a warning and
        fill the rows in order from the source's good fibers, instead
        of failing. New SKY SUBST_MAP / SUBST_SKYEXT provenance
        keywords.

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

# maps a telescope name to its per-fiber sky-model extension in the lvmCFrame
SKY_EXTENSIONS = {'SkyE': 'SKY_EAST', 'SkyW': 'SKY_WEST'}

# SkyE/SkyW combination weights: a joint property of both telescopes
# (set by skyMethod.combine_skies), not a per-telescope one, so never
# block-copied -- recomputed below from the updated pointings instead
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
    Replace the FLUX/IVAR/MASK/LSF fiber rows and SKY_EAST/SKY_WEST
    extension for target_tel in target_cframe with source_tel's from
    source_cframe (rows matched by fiberid for the same telescope, in
    order from the source's good fibers otherwise), copy the
    associated telescope header block, recompute SKYEW/SKYWW, and
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

        if target_tel == source_tel:
            # same telescope: fiber assignment is fixed hardware, identical
            # across exposures, so match fiber-for-fiber by fiberid
            mapping = 'fiberid'
            src_tel_at_rows = np.asarray(slit_s['telescope'])[rows]
            bad = src_tel_at_rows != source_tel
            if bad.any():
                example = fiberids[bad][0]
                example_tel = src_tel_at_rows[bad][0]
                raise ValueError(
                    f"{bad.sum()} of the target's {target_tel} fibers are not tagged '{source_tel}' "
                    f"in the source SLITMAP (e.g. fiberid {example} is '{example_tel}' there)"
                )
            src_rows = rows
        else:
            # different telescopes: there is no fiber-for-fiber correspondence,
            # so fill the target rows in order from the source telescope's good
            # fibers, cycling through them again if there are fewer. The DRP
            # (and the XCframe summaries) average over a telescope's fibers,
            # so which source fiber lands in which target row does not matter.
            mapping = 'order'
            good_src = (slit_s['telescope'] == source_tel) & (slit_s['fibstatus'] == 0)
            src_fiberids = np.asarray(slit_s['fiberid'][good_src])
            if src_fiberids.size == 0:
                raise ValueError(f"no good (fibstatus==0) fibers tagged '{source_tel}' found in source SLITMAP")
            src_rows = np.resize(src_fiberids - 1, rows.size)
            print(f"Warning: target_tel {target_tel} != source_tel {source_tel}; filling the "
                  f"{rows.size} {target_tel} rows in order from the {src_fiberids.size} good "
                  f"{source_tel} fibers (not matched by fiberid)")

        for ext in EXTENSIONS_TO_COPY:
            ht[ext].data[rows] = hs[ext].data[src_rows]

        # SKY_EAST/SKY_WEST (+ _IVAR) are the target telescope's sky model
        # evaluated at every fiber; many non-DRP tools (SummarizeCframe.py,
        # GetTelData.py, lvm_skyfit.py, ...) read them, so they must be
        # replaced along with the fiber rows to keep the CFrame consistent
        src_sky = SKY_EXTENSIONS[source_tel]
        dst_sky = SKY_EXTENSIONS[target_tel]
        for suffix in ('', '_IVAR'):
            if src_sky + suffix in hs and dst_sky + suffix in ht:
                ht[dst_sky + suffix].data = hs[src_sky + suffix].data.copy()
            else:
                print(f"Warning: {src_sky + suffix} or {dst_sky + suffix} missing; not replaced")

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

        # SKYEW/SKYWW: recompute with skyMethod.combine_skies's own formula
        # (inverse angular distance to the science field, normalized), so
        # they describe the sky data actually now in the file
        if 'SKYEW' in hdr_t and 'SKYWW' in hdr_t:
            ad_e = ang_distance(hdr_t['SKYERA'], hdr_t['SKYEDEC'], hdr_t['SCIRA'], hdr_t['SCIDEC'])
            ad_w = ang_distance(hdr_t['SKYWRA'], hdr_t['SKYWDEC'], hdr_t['SCIRA'], hdr_t['SCIDEC'])
            w_e = 1 / (ad_e if ad_e > 0 else 1)
            w_w = 1 / (ad_w if ad_w > 0 else 1)
            w_norm = w_e + w_w
            hdr_t['SKYEW'] = (float(w_e / w_norm), hdr_t.comments['SKYEW'])
            hdr_t['SKYWW'] = (float(w_w / w_norm), hdr_t.comments['SKYWW'])

        hdr_t['HIERARCH SKY SUBST_TEL'] = (target_tel, 'telescope whose fiber data was substituted')
        hdr_t['HIERARCH SKY SUBST_SRC'] = (os.path.basename(source_cframe), 'source CFrame for substituted sky')
        hdr_t['HIERARCH SKY SUBST_SRCTEL'] = (source_tel, 'telescope in source CFrame supplying the data')
        hdr_t['HIERARCH SKY SUBST_SRCEXP'] = (hdr_s.get('EXPOSURE', -1), 'source exposure number')
        hdr_t['HIERARCH SKY SUBST_NFIB'] = (int(rows.size), 'number of fibers substituted')
        hdr_t['HIERARCH SKY SUBST_MAP'] = (mapping, 'fiber mapping: fiberid or order')
        hdr_t['HIERARCH SKY SUBST_SKYEXT'] = (f'{src_sky}->{dst_sky}', 'sky model ext replaced')

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
