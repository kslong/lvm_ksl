#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Perform sky subtraction on an XCframe summary file using the PALACE
    (Ivan Katkov / XSkySepIvan) spectral decomposition to separate the
    sky into emission lines and continuum.  For each row the relevant
    spectra are decomposed with the PALACE model; the resulting line and
    continuum components are then combined and subtracted from the
    science spectrum without any additional scaling.

    Three methods are supported::

        scilines_nearcont lines fitted directly from the science
                          spectrum itself, continuum from the near
                          sky telescope (default)
        nearest           continuum and lines both from the nearest
                          sky telescope
        farlines_nearcont lines from the far sky, continuum from the
                          near sky

    Produces a FITS file with WAVE, FLUX (sky-subtracted), SKY, and
    DRP_ALL extensions.  The DRP_ALL table carries a QA_FLAGS column
    that records per-row quality issues.

Command line usage (if any):

    usage: SkySubDev2.py [-method METHOD] [-delta N] [-lsf FWHM]
                         [-out ROOT] filename

    Arguments::

        filename    XCframe FITS file to process

    Options::

        -method METHOD   sky subtraction method: scilines_nearcont |
                         nearest | farlines_nearcont
                         (default: scilines_nearcont)
        -delta N         process every N-th row; useful for quick tests
                         (default: 1 = all rows)
        -lsf FWHM        constant LSF FWHM in Angstroms (default: 1.3),
                         used only if the input file has no LSF extension
        -out ROOT        output filename root; default is <stem>_dev2_<method>

Description:

    The PALACE decomposer (SkyDecomp) is built from the full wavelength
    grid and an assumed LSF, chosen in priority order:

    1. This file's own per-row, per-wavelength LSF, if it has an LSF
       extension (written by SummarizeCframe.py's make_med_spec from the
       source CFrame's own LSF, one row of FWHM-vs-wavelength per
       exposure) -- the decomposer is rebuilt whenever a row's LSF
       differs from the previous row's (see Notes for the cost).
    2. Otherwise, a representative wavelength-dependent reference curve
       (data/lsf.fits, derived once from real per-row LSF data), if
       found -- the same curve every row, so the decomposer is built
       only once for the whole file (same cost as case 3 below).
    3. Otherwise, the constant -lsf FWHM, built once and reused for
       every row.

    A flat 1.3 A FWHM (case 3's default) is normally too narrow for at
    least part of the real spectrum, making the fitted LINES component
    systematically too narrow -- cases 1 and 2 both capture the true
    wavelength dependence instead.  For each row:

    1. Science and sky spectra are read for the row.
    2. Near/far sky is determined from sci_ra/dec, skye_ra/dec,
       skyw_ra/dec in DRP_ALL.
    3. The near-sky spectrum, and whichever other spectrum the chosen
       method needs (far-sky for farlines_nearcont, the science flux
       itself for scilines_nearcont), are decomposed by the PALACE
       model into::

           LINES = oh + atom + orc + o2
           CONT  = moon + diffuse

    4. The sky model is assembled without scaling::

           scilines_nearcont:  sky = CONT_near + LINES_sci
           farlines_nearcont:  sky = CONT_near + LINES_far
           nearest:            sky = CONT_near + LINES_near

       CONT always comes from the near-sky telescope in all three
       methods -- fitting a continuum from the science spectrum itself
       would be confounded by real astrophysical continuum, so PALACE's
       CONT component is never taken from there, only its LINES
       component (scilines_nearcont).  scilines_nearcont is the default:
       fitting the sky lines directly from the fiber being corrected
       avoids relying on a sky telescope's lines matching the science
       fiber's actual sky-line amplitude, which farlines_nearcont/
       nearest implicitly assume.

    5. sky-subtracted science = flux_sci - sky

    No scale factor is applied anywhere in this method (unlike
    SkySubOrig.py/SkySubDev1.py/SkySepESO.py), so there is no
    DRP_ALL['LINE_SCALE'] column.  Continuum-fit quality against the RAW
    (pre-subtraction) science and sky spectra is still evaluated and
    recorded per row, since PALACE's CONT component is a genuine
    continuum fit even though it is never scaled: if sky_mask.fits is
    found (same convention as SkySub_eval.py), GetSkyCont.arm_continuum_stats
    computes, per spectrograph arm (B/R/Z), the median/NMAD/RMS/skew of
    FLUX-CONT in clean (sky-line-free) pixels -- for the science spectrum
    (SCI_MED_B etc.; for scilines_nearcont this reuses the science-side
    PALACE decomposition already required by the subtraction itself, at
    no extra cost -- for the other two methods it is an extra
    decomposition run purely for this check) and for the near-sky
    spectrum whose CONT is used in the final SKY (SKY_MED_B etc.).  This
    is unrelated to the final sky-subtracted SCI_FLUX -- see
    SkySub_eval.py's Figures 5/6 for that (post-subtraction) residual,
    which measures leftover source continuum rather than fit quality.

    QA flag bits stored in DRP_ALL['QA_FLAGS']::

        0x01  NANDATA   NaN/inf found in input flux or sky data
        0x08  FAILED    row raised an exception; spectrum filled with NaN

Notes:

    Output filename is ``<ROOT>.fits``.  If -out is omitted the name
    is derived as ``<stem>_dev2_<method>.fits`` where ``<stem>`` is the
    input filename without extension.

    Requires the PALACE library, vendored in py_progs/sky_decomp/ (260709;
    previously an external dependency on the lvmsky repository and its
    lvmdrp26-only lvmdrp import) and its reference data in data/palace_ref/;
    paths are taken from XSkySepIvan.py's DEFAULT_BASE_DIR.

    Continuum-quality columns (SCI_*/SKY_* per arm) require sky_mask.fits;
    searched for in the current directory, then in the lvm_ksl data/
    directory.  If not found, do_all() prints a warning and those columns
    are skipped entirely (everything else is unaffected).  For
    farlines_nearcont/nearest, the extra per-row science-side PALACE
    decomposition that feeds them is also skipped in that case; for
    scilines_nearcont that decomposition happens regardless, since the
    subtraction itself needs it.

    Per-row LSF (when the input file has an LSF extension): rebuilding the
    PALACE decomposer costs ~2s (measured) each time a row's LSF differs
    from the previous row's, on top of the ~1s/decomposition already spent
    in one_drp() -- real added cost for a full run, not just a one-time
    setup cost as with the constant-LSF or reference-curve paths.  Any
    non-finite or non-positive FWHM value at a given wavelength in a row's
    own LSF array falls back to the constant -lsf default at that
    wavelength only (not the whole row).  DRP_ALL gains LSF_FWHM_MED
    (this row's median LSF FWHM actually used, whichever source it came
    from); the primary header gains LSFSRC recording which of the three
    sources was used for the whole run.

    Reference LSF curve (data/lsf.fits): a FITS bintable with WAVE, LSF
    (FWHM in Angstroms), LSF_STD columns, searched for in the current
    directory then in the lvm_ksl data/ directory (same convention as
    sky_mask.fits).  Interpolated (linear, edge-extrapolated) onto this
    file's own wavelength grid.  Not derived here -- expected to already
    exist as a representative curve built once from real per-row LSF
    data, analogous to how sky_mask.fits is a pre-built reference rather
    than something each script derives for itself.

History::

    260630 ksl Coding begun; imports PALACE decomposer from XSkySepIvan.py
    260706 ksl DRP_ALL['mjd'] now recomputed precisely from 'obstime' via
               SkySubOrig.obstime_to_mjd(), instead of the truncated
               integer carried through from the input file.
    260708 ksl Added per-arm continuum-fit-quality columns (SCI_*/SKY_* for
               med/nmad/rms/skew x b/r/z), evaluated against the raw
               pre-subtraction science and near-sky spectra using
               GetSkyCont.arm_continuum_stats() (same optional sky_mask.fits
               convention as SkySubOrig.py).  Requires one extra PALACE
               decomposition per row (on the science flux) that the
               subtraction itself does not otherwise need; skipped when no
               mask is available.  Still no LINE_SCALE column -- this
               method genuinely applies no scale factor.
    260709 ksl Added a third method, scilines_nearcont, now the default:
               fits the sky LINES component directly from the science
               spectrum itself (sky = CONT_near + LINES_sci) instead of
               borrowing LINES from a sky telescope fiber.  CONT is still
               always taken from the near-sky telescope in every method.
               one_drp() now decomposes the science flux whenever the
               method needs it (scilines_nearcont) or the continuum-
               quality check needs it (clean_mask given), whichever
               applies -- never both, so there's still only ever one
               science-side PALACE decomposition per row.
    260709 ksl Use the input file's own per-row, per-wavelength LSF (the
               LSF extension written by SummarizeCframe.py's make_med_spec)
               when present, instead of always assuming the constant -lsf
               FWHM (default 1.3 A) -- user found the fitted lines were
               systematically too narrow on real data, traced to this flat
               default not matching the true wavelength-dependent LSF.
               do_all() rebuilds the PALACE decomposer per row when this
               row's LSF differs from the previous row's (_get_decomposer()
               already only caches one instance at a time, so this is
               memory-safe, just slower: ~2s/rebuild, measured). Falls back
               to the constant -lsf default (unchanged behaviour) if no LSF
               extension is found. New DRP_ALL['LSF_FWHM_MED'] and primary
               header LSFSRC record what was actually used per row/run.
    260709 ksl Added a middle tier between the per-row LSF extension and
               the flat -lsf constant: a representative wavelength-
               dependent reference curve (data/lsf.fits), for the common
               case of a file with no LSF extension of its own.  Since the
               same curve is used for every row, _get_decomposer()'s
               single-instance cache means it's still only ever built
               once for the whole file -- no per-row cost, unlike the
               true per-row-LSF path, but still captures the real
               wavelength dependence a flat FWHM misses.  Implemented by
               broadcasting the reference curve to every row and reusing
               the same per-row code path as the true per-row case
               (have_lsf_ext=True), so no special-casing was needed
               downstream (NaN fallback, LSF_FWHM_MED, etc. all just work).
    260709 ksl Even with a real per-row LSF extension, fitted lines were
               still measurably too narrow: tested on a real file with its
               own LSF extension by comparing the fitted model's peak
               height against its integrated area at three well-isolated
               sky lines ([OI] 5577/6300/6364) -- area matched the data to
               within a few percent (as expected, least-squares conserves
               flux), but the model peak came out 4-7% too tall for that
               area, i.e. the effective fitted width is 4-7% narrower than
               the true line, varying somewhat by wavelength. Ruled out
               PALACE's own per-channel LSF refit (n_lsf_refits) as the
               cause -- identical result for 0 through 3 refit iterations.
               Added -lsf_boost (default 1.08, the middle of the measured
               4-7% range), a multiplicative correction applied to
               whatever LSF FWHM is used (any of the three sources) --
               a first-pass empirical fix, not a rigorously
               wavelength-dependent one, pending further investigation of
               whether the correction itself should vary with wavelength.
    260709 ksl Reverted -lsf_boost's default to 1.0 (no correction). The
               260709 entry above was based on a single row's peak/area
               ratio at three isolated lines; a proper test comparing
               boost=1.00/1.08/1.20 over a 280-row median showed the
               opposite of what that single-row test predicted -- the
               noise-corrected HF RMS residual (SkySub_eval.py's own
               diagnostic) got monotonically *worse* with more boost at
               every diagnostic window, and the oversubtracted-center/
               adjacent-bump pattern at [OI] 5577 got monotonically
               deeper and larger, not smaller, as the boost increased.
               So the assumed LSF width is evidently not the (or not the
               only) cause of that residual pattern, and widening it
               further actively hurts. -lsf_boost is left in place for
               further investigation but must not be re-enabled by
               default without re-verifying against a multi-row sample,
               not a single row.

'''

import sys
import os
from pathlib import Path

# ensure py_progs siblings are importable when running directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
from scipy.interpolate import interp1d

from XSkySepIvan import _get_decomposer, estimate_ivar, DEFAULT_BASE_DIR
from SkySubOrig import obstime_to_mjd

try:
    from GetSkyCont import (load_mask, _interp_mask_to_wave,
                            arm_continuum_stats, flatten_arm_stats)
    _HAVE_MASK = True
except ImportError:
    _HAVE_MASK = False

# ──────────────────────────────────────────────────────────────
# QA flag bits
# ──────────────────────────────────────────────────────────────
QA_NANDATA = 1   # NaN or inf values found in input flux or sky data
QA_FAILED  = 8   # row failed entirely; FLUX and SKY are NaN

_QA_FLAG_NAMES = {
    QA_NANDATA: 'NANDATA',
    QA_FAILED:  'FAILED',
}

_USAGE = '''Usage:
  SkySubDev2.py [-method METHOD] [-delta N] [-lsf FWHM] [-lsf_boost FACTOR]
                [-out ROOT] filename

Arguments:
  filename         XCframe FITS file to process

Options:
  -method METHOD   scilines_nearcont | nearest | farlines_nearcont
                   (default: scilines_nearcont)
  -delta N         step size through rows for quick tests (default: 1)
  -lsf FWHM        constant LSF FWHM in Angstroms (default: 1.3), used only
                   as the final fallback, or to fill bad pixels in a real
                   per-row LSF array
  -lsf_boost FACTOR  multiplicative correction applied to whatever LSF is
                   used, any source (default: 1.0, i.e. no correction --
                   a boost was tried as a first-pass fix but a proper
                   multi-row test showed it makes the residual worse,
                   not better; see docstring)
  -out ROOT        output filename root (default: <stem>_dev2_<method>)
'''


# ──────────────────────────────────────────────────────────────
# Reference wavelength-dependent LSF curve (fallback for files with no
# per-row LSF extension of their own)
# ──────────────────────────────────────────────────────────────

def _load_ref_lsf(lsf_file):
    '''
    Load a reference wavelength-dependent LSF FWHM curve: a FITS bintable
    with WAVE, LSF (FWHM in Angstroms), LSF_STD columns (e.g. data/lsf.fits,
    a representative curve derived once from real per-row LSF data, for use
    when an individual input file carries no LSF extension of its own).

    Returns (wave, fwhm) 1-D arrays.
    '''
    with fits.open(lsf_file) as hdul:
        tab = Table(hdul[1].data)
    return np.asarray(tab['WAVE'], dtype=float), np.asarray(tab['LSF'], dtype=float)


def _interp_lsf_to_wave(ref_wave, ref_fwhm, spec_wave):
    '''Linearly interpolate a reference FWHM(wave) curve onto another wavelength grid.'''
    f = interp1d(ref_wave, ref_fwhm, kind='linear', bounds_error=False,
                fill_value=(ref_fwhm[0], ref_fwhm[-1]))
    return f(spec_wave)


# ──────────────────────────────────────────────────────────────
# PALACE decomposition of a raw flux array
# ──────────────────────────────────────────────────────────────

def _decompose(flux, wave, decomposer):
    '''
    Run the PALACE decomposer on a single flux array and return
    (lines, cont) where lines = oh+atom+orc+o2 and cont = moon+diffuse.
    '''
    flux = np.asarray(flux, float)
    ivar = estimate_ivar(flux, wave)

    finite = np.isfinite(flux)
    flux_scale = float(np.sqrt(np.nanmean(flux[finite] ** 2))) if finite.any() else 1.0
    flux_scale = max(flux_scale, 1e-30)

    result = decomposer.fit(flux / flux_scale, ivar * flux_scale ** 2)
    c      = result.components
    lines  = (c['oh'] + c['atom'] + c['orc'] + c['o2']) * flux_scale
    cont   = (c['moon'] + c['diffuse']) * flux_scale
    return lines, cont


# ──────────────────────────────────────────────────────────────
# Per-row sky subtraction
# ──────────────────────────────────────────────────────────────

def one_drp(xfits, drp_all, row, wave, decomposer,
            method='scilines_nearcont', clean_mask=None):
    '''
    Sky-subtract a single row using PALACE decomposition with no scaling.

    method : 'scilines_nearcont' (default) uses the science spectrum's own
        PALACE-fitted LINES component directly (the sky emission actually
        present in the science fiber, rather than borrowed from a sky
        telescope); 'farlines_nearcont' and 'nearest' instead borrow LINES
        from the far/near sky telescope, as before.  CONT always comes
        from the near-sky telescope in all three methods -- a science-side
        continuum fit is confounded by real astrophysical continuum, so it
        is never used for CONT, only (optionally) evaluated below.

    Returns (result_table, qa_flags, cont_stats).  result_table has
    columns WAVE, SCI_FLUX, SKY.

    cont_stats is a flat dict of per-arm continuum-fit-quality stats
    (GetSkyCont.arm_continuum_stats/flatten_arm_stats), evaluated on the
    RAW pre-subtraction science spectrum (lines_sci -- already computed as
    part of the subtraction itself for 'scilines_nearcont'; from an extra
    PALACE decomposition otherwise) and the near-sky spectrum (lines_near;
    cont_near is always the continuum used in the final SKY, for all three
    methods).  Keys are 'sci_<stat>_<arm>' and 'sky_<stat>_<arm>'.  Empty
    dict if clean_mask is None or GetSkyCont is unavailable -- the extra
    science-side decomposition (when method doesn't already need it) is
    skipped entirely in that case, since it exists only to feed this
    evaluation.

    Returns (None, QA_FAILED, {}) on error.
    '''
    qa_flags = 0

    flux      = np.array(xfits['FLUX'].data[row],     dtype=float)
    skye_flux = np.array(xfits['SKY_EAST'].data[row], dtype=float)
    skyw_flux = np.array(xfits['SKY_WEST'].data[row], dtype=float)

    if not (np.all(np.isfinite(flux)) and
            np.all(np.isfinite(skye_flux)) and
            np.all(np.isfinite(skyw_flux))):
        qa_flags |= QA_NANDATA
        # nan_to_num's default replaces +-inf with +-1.8e308 (float64 max),
        # not 0 -- that "poison" value overflows through the PALACE
        # decomposer's internal fit, corrupting the row.  Zero all
        # non-finite pixels explicitly instead.
        flux      = np.nan_to_num(flux,      nan=0.0, posinf=0.0, neginf=0.0)
        skye_flux = np.nan_to_num(skye_flux, nan=0.0, posinf=0.0, neginf=0.0)
        skyw_flux = np.nan_to_num(skyw_flux, nan=0.0, posinf=0.0, neginf=0.0)

    # determine near/far sky from angular separation
    sci_coord  = SkyCoord(ra=drp_all['sci_ra'][row]   * u.degree,
                          dec=drp_all['sci_dec'][row]  * u.degree)
    skye_coord = SkyCoord(ra=drp_all['skye_ra'][row]  * u.degree,
                          dec=drp_all['skye_dec'][row] * u.degree)
    skyw_coord = SkyCoord(ra=drp_all['skyw_ra'][row]  * u.degree,
                          dec=drp_all['skyw_dec'][row] * u.degree)
    if sci_coord.separation(skye_coord) < sci_coord.separation(skyw_coord):
        sky_near, sky_far = skye_flux, skyw_flux
    else:
        sky_near, sky_far = skyw_flux, skye_flux

    # PALACE decomposition — no scaling applied
    lines_near, cont_near = _decompose(sky_near, wave, decomposer)

    # The science flux is decomposed if the method needs its LINES component
    # directly (scilines_nearcont), or purely for the continuum-fit-quality
    # check below (clean_mask given) -- either way, only ever once per row.
    need_sci_decomp = (method == 'scilines_nearcont'
                       or (clean_mask is not None and _HAVE_MASK))
    lines_sci = None
    if need_sci_decomp:
        lines_sci, _ = _decompose(flux, wave, decomposer)

    if method == 'farlines_nearcont':
        lines_far, _ = _decompose(sky_far, wave, decomposer)
        sky = cont_near + lines_far
    elif method == 'nearest':
        sky = cont_near + lines_near
    elif method == 'scilines_nearcont':
        # Fit the sky emission lines directly from the science spectrum
        # itself (rather than borrowing them from a sky telescope fiber),
        # keeping only the near-sky CONT (a science-fiber continuum fit is
        # confounded by real astrophysical continuum, so CONT still has to
        # come from a sky telescope).
        sky = cont_near + lines_sci
    else:
        print('Error: unknown method "%s"' % method)
        return None, QA_FAILED, {}

    cont_stats = {}
    if clean_mask is not None and _HAVE_MASK:
        cont_stats.update(flatten_arm_stats(
            'sci', arm_continuum_stats(wave, lines_sci, clean_mask)))
        cont_stats.update(flatten_arm_stats(
            'sky', arm_continuum_stats(wave, lines_near, clean_mask)))

    result = Table([wave, flux - sky, sky], names=['WAVE', 'SCI_FLUX', 'SKY'])
    return result, qa_flags, cont_stats


# ──────────────────────────────────────────────────────────────
# Batch processing
# ──────────────────────────────────────────────────────────────

def do_all(filename, method='scilines_nearcont', idelta=1,
           fwhm_lsf=1.3, lsf_boost=1.0, outroot=''):
    '''
    Process every row of an XCframe file and write sky-subtracted output.

    Parameters
    ----------
    filename : str
        Path to the input XCframe FITS file.
    method : str
        Sky subtraction method (scilines_nearcont, nearest, or
        farlines_nearcont -- see one_drp()).
    idelta : int
        Row step size (1 = all rows).
    fwhm_lsf : float
        LSF FWHM in Angstroms, used only as the final fallback (see
        priority order below) or to fill bad pixels in a real per-row
        LSF array.
    lsf_boost : float
        Multiplicative correction applied to whatever LSF FWHM is used
        (real per-row, reference curve, or constant), regardless of
        source.  Default 1.0 (no correction).  A single-row peak/area
        test on real data initially suggested a ~1.08 boost would help
        (see History), but a proper multi-row (280-row median) test
        showed the opposite: increasing the boost made the noise-
        corrected HF RMS residual monotonically *worse* at every
        diagnostic window, and made the oversubtracted-center/adjacent-
        bump pattern at [OI] 5577 monotonically deeper/larger, not
        smaller. The single-row test was not representative -- whatever
        is causing that residual pattern is evidently not a pure
        Gaussian-LSF-width deficit, so widening it further just makes
        things worse. Left available for further investigation, but do
        not re-enable as the default without re-testing against a
        proper multi-row sample first.
    outroot : str
        Output filename root; defaults to <stem>_dev2_<method>.
    '''
    x       = fits.open(filename)
    drp_all = Table(x['DRP_ALL'].data)
    wave    = np.array(x['WAVE'].data, dtype=float)

    # LSF source, in priority order:
    #  1. This file's own per-row, per-wavelength LSF (FWHM in Angstroms,
    #     written by SummarizeCframe.py's make_med_spec from the source
    #     CFrame's LSF extension) -- best available, but costs a decomposer
    #     rebuild (~2s, measured) whenever the row-to-row LSF actually changes.
    #  2. A representative wavelength-dependent reference curve (data/lsf.fits,
    #     derived once from real per-row LSF data), for files with no LSF
    #     extension of their own -- captures the true wavelength dependence
    #     (a flat FWHM is normally too narrow for at least part of the
    #     spectrum) at no extra per-row cost, since it's the SAME array every
    #     row: _get_decomposer()'s single-instance cache key only changes
    #     (and only then rebuilds) when the LSF array actually differs from
    #     the previous call, so this reduces to one build for the whole file,
    #     same as the constant-LSF path below.
    #  3. The constant -lsf FWHM, if neither of the above is available.
    have_lsf_ext = 'LSF' in [h.name for h in x]
    if have_lsf_ext:
        lsf_fwhm_ext = np.array(x['LSF'].data, dtype=float) * lsf_boost
        print('Using the per-row, per-wavelength LSF from the LSF extension '
              '(median %.2f A after a %.2fx boost) instead of the constant '
              '-lsf default (%.2f A).  This rebuilds the PALACE decomposer '
              'for each row whose LSF differs from the previous one '
              '(~2s/rebuild) -- slower than the constant-LSF path, but the '
              'fitted lines should no longer be systematically too narrow.'
              % (np.nanmedian(lsf_fwhm_ext), lsf_boost, fwhm_lsf))
        lsf_source = 'LSF extension (per-row)'
    else:
        ref_fwhm = None
        _data_dir = Path(__file__).parent.parent / 'data'
        for _candidate in [Path('lsf.fits'), _data_dir / 'lsf.fits']:
            if _candidate.exists():
                try:
                    _ref_wave, _ref_fwhm_raw = _load_ref_lsf(str(_candidate))
                    ref_fwhm = _interp_lsf_to_wave(_ref_wave, _ref_fwhm_raw, wave) * lsf_boost
                    print('No LSF extension in this file; using the reference '
                          'wavelength-dependent LSF curve %s (median %.2f A '
                          'after a %.2fx boost) instead of the constant -lsf '
                          'default (%.2f A).'
                          % (_candidate, np.nanmedian(ref_fwhm), lsf_boost, fwhm_lsf))
                except Exception as _e:
                    print('Warning: could not load reference LSF %s (%s)' % (_candidate, _e))
                break
        if ref_fwhm is not None:
            # Broadcast the same curve to every row so the per-row loop below
            # (NaN fallback, LSF_FWHM_MED column, etc.) needs no special case
            # -- and, since the array is identical every row, the decomposer
            # is only ever built once (see priority-2 note above).
            lsf_fwhm_ext = np.tile(ref_fwhm, (len(drp_all), 1))
            have_lsf_ext = True
            lsf_source   = 'reference curve (data/lsf.fits)'
        else:
            lsf_fwhm_ext = None
            lsf_sigma  = (fwhm_lsf * lsf_boost) / 2.355
            decomposer = _get_decomposer(wave, lsf_sigma, base_dir=DEFAULT_BASE_DIR)
            lsf_source = 'constant -lsf default'
            print('Warning: no LSF extension or reference LSF curve found; '
                  'using the constant -lsf %.2f A (after a %.2fx boost) for '
                  'every row and wavelength' % (fwhm_lsf * lsf_boost, lsf_boost))

    # Load sky-line mask once (same convention as SkySub_eval.py) purely for
    # the continuum-fit-quality evaluation below; the subtraction itself is
    # unaffected if the mask is unavailable.
    clean_mask = None
    if _HAVE_MASK:
        _data_dir = Path(__file__).parent.parent / 'data'
        for _candidate in [Path('sky_mask.fits'), _data_dir / 'sky_mask.fits']:
            if _candidate.exists():
                try:
                    _mask_wave, _mask_arr = load_mask(str(_candidate))
                    clean_mask = _interp_mask_to_wave(_mask_wave, _mask_arr, wave)
                    print('Loaded sky mask for continuum-quality evaluation: %s' % _candidate)
                except Exception as _e:
                    print('Warning: could not load mask %s (%s)' % (_candidate, _e))
                break
    if clean_mask is None:
        print('Warning: sky_mask.fits not found; skipping continuum-quality columns')

    nan_spectrum = np.full(len(wave), np.nan)

    final_flux      = []
    final_sky       = []
    select          = []
    qa_flags_list   = []
    cont_stats_list = []

    i = 0
    while i < len(drp_all):
        if have_lsf_ext:
            # Per-row FWHM-in-Angstroms -> sigma; fall back to the constant
            # default at any wavelength where this row's LSF is NaN/non-
            # positive (e.g. a masked pixel in the source CFrame's LSF
            # extension), rather than feeding a bad value into the decomposer.
            row_fwhm = lsf_fwhm_ext[i].copy()
            bad = ~np.isfinite(row_fwhm) | (row_fwhm <= 0)
            if bad.any():
                row_fwhm[bad] = fwhm_lsf * lsf_boost
            decomposer = _get_decomposer(wave, row_fwhm / 2.355, base_dir=DEFAULT_BASE_DIR)
        try:
            ftab, row_flags, cont_stats = one_drp(x, drp_all, i, wave, decomposer,
                                                  method=method, clean_mask=clean_mask)
        except Exception as e:
            print('Row %d: exception (%s)' % (i, e))
            ftab       = None
            row_flags  = 0
            cont_stats = {}

        if ftab is None:
            row_flags |= QA_FAILED
            final_flux.append(nan_spectrum.copy())
            final_sky.append(nan_spectrum.copy())
        else:
            final_flux.append(np.array(ftab['SCI_FLUX']))
            final_sky.append(np.array(ftab['SKY']))

        select.append(i)
        qa_flags_list.append(row_flags)
        cont_stats_list.append(cont_stats)
        i += idelta
        if i % 100 == 0:
            print('Completed %6d of %d in steps of %d' % (i, len(drp_all), idelta))

    n_failed = sum(1 for f in qa_flags_list if f & QA_FAILED)
    n_warned = sum(1 for f in qa_flags_list if f != 0 and not (f & QA_FAILED))
    print('\nProcessed %d rows: %d failed (NaN fill), %d with warnings'
          % (len(select), n_failed, n_warned))

    flagged = [(select[j], qa_flags_list[j])
               for j in range(len(select)) if qa_flags_list[j] != 0]
    if flagged:
        print('Rows with QA flags:')
        for orig_row, flags in flagged:
            active = [name for bit, name in _QA_FLAG_NAMES.items() if flags & bit]
            print('  Row %6d  flags=0x%02x  (%s)' % (orig_row, flags, ', '.join(active)))

    hdu1 = fits.PrimaryHDU(data=None)
    hdu1.header['Title']   = 'SkySubDev2'
    hdu1.header['METHOD']  = method
    hdu1.header['FWHMLSF'] = fwhm_lsf
    hdu1.header['LSFBOOST'] = lsf_boost
    hdu1.header['LSFSRC']  = lsf_source
    hdu2 = fits.ImageHDU(data=wave.astype(np.float32),       name='WAVE')
    hdu3 = fits.ImageHDU(data=np.array(final_flux),          name='FLUX')
    hdu4 = fits.ImageHDU(data=np.array(final_sky),           name='SKY')

    xtab = drp_all[select].copy()
    xtab['QA_FLAGS'] = np.array(qa_flags_list, dtype=np.int32)
    if have_lsf_ext:
        xtab['LSF_FWHM_MED'] = np.array(
            [np.nanmedian(lsf_fwhm_ext[j]) for j in select], dtype=np.float32)
    # Continuum-fit-quality columns (raw pre-subtraction sci/sky spectra,
    # not the final sky-subtracted result -- see one_drp() docstring).
    # No LINE_SCALE column here: this method applies no scale factor.
    _cont_keys = sorted({k for d in cont_stats_list for k in d})
    for _key in _cont_keys:
        xtab[_key.upper()] = np.array(
            [d.get(_key, np.nan) for d in cont_stats_list], dtype=np.float32)
    if 'obstime' in xtab.colnames and 'mjd' in xtab.colnames:
        xtab['mjd'] = obstime_to_mjd(xtab['obstime'])
    hdu5 = fits.BinTableHDU(xtab, name='DRP_ALL')

    dwave = float(wave[1] - wave[0]) if len(wave) > 1 else 0.5
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [1, 1]
    wcs.wcs.crval = [float(wave[0]), 0]
    wcs.wcs.cdelt = [dwave, 1]
    wcs.wcs.ctype = ['WAVE', 'LINE']
    hdu3.header.update(wcs.to_header())
    hdu4.header.update(wcs.to_header())

    hdul = fits.HDUList([hdu1, hdu2, hdu3, hdu4, hdu5])

    if outroot == '':
        stem    = os.path.splitext(os.path.basename(filename))[0]
        outroot = '%s_dev2_%s' % (stem, method)
    outfile = '%s.fits' % outroot
    hdul.writeto(outfile, overwrite=True)
    print('Wrote results to %s' % outfile)
    x.close()


# ──────────────────────────────────────────────────────────────
# Command-line entry point
# ──────────────────────────────────────────────────────────────

if __name__ == '__main__':
    argv = sys.argv[1:]
    if not argv or '-h' in argv or '--help' in argv:
        print(_USAGE)
        sys.exit(0)

    method    = 'scilines_nearcont'
    idelta    = 1
    fwhm_lsf  = 1.3
    lsf_boost = 1.0
    outroot   = ''
    filename  = None

    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == '-method':
            i += 1
            method = argv[i]
        elif arg == '-delta':
            i += 1
            idelta = int(argv[i])
        elif arg == '-lsf':
            i += 1
            fwhm_lsf = float(argv[i])
        elif arg == '-lsf_boost':
            i += 1
            lsf_boost = float(argv[i])
        elif arg == '-out':
            i += 1
            outroot = argv[i]
        elif arg.startswith('-'):
            print('Error: unknown option "%s"' % arg)
            print(_USAGE)
            sys.exit(1)
        else:
            if filename is not None:
                print('Error: unexpected argument "%s" (filename already set to "%s")'
                      % (arg, filename))
                print(_USAGE)
                sys.exit(1)
            filename = arg
        i += 1

    if filename is None:
        print('Error: no filename supplied')
        print(_USAGE)
        sys.exit(1)

    if not os.path.exists(filename):
        print('Error: file not found: %s' % filename)
        sys.exit(1)

    valid_methods = {'scilines_nearcont', 'nearest', 'farlines_nearcont'}
    if method not in valid_methods:
        print('Error: -method must be one of: %s' % ', '.join(sorted(valid_methods)))
        sys.exit(1)

    do_all(filename=filename, method=method, idelta=idelta,
           fwhm_lsf=fwhm_lsf, lsf_boost=lsf_boost, outroot=outroot)
