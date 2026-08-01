#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Evaluate how closely real LVM spectra (science and/or sky-telescope)
    match the predictions of the PALACE airglow model, for continuum and
    lines separately -- a direct analog of SkyObsESOCompare.py, with the
    same disjoint clean/line pixel-set fit, the same ``SEP_``/final split, and
    the same output column convention, so SkyObsESO_analysis.py's plotting
    pipeline (make_standard_plots, plot_fit_quality_correlation, etc.) runs
    unchanged against this script's output for a direct ESO-vs-PALACE
    comparison.

    PALACE (see PalaceObs.py) predicts airglow only -- lines and an
    HO2/FeO/O2 continuum -- with no Moon or Zodiacal-light component at
    all.  Rather than invent a bespoke Moon/Zodi treatment, this script
    borrows the ESO Sky Model's own MOON/ZODI prediction (the SAME
    SkySepESO._get_sky_model call SkyObsESOCompare.py itself uses) and
    combines it with PALACE's own DIFFUSE (airglow continuum) and LINES,
    all three already in the same physical flux units -- see Description.

Command line usage (if any):

    usage: SkyObsPalaceCompare.py [-ext EXT [-ext EXT ...]]
                                  [-engine local|remote|auto] [-mask FILE]
                                  [-delta N] [-out ROOT] filename

    Arguments::

        filename    An XCframe/XSFrame summary FITS file (FLUX, SKY_EAST,
                    SKY_WEST, LSF, DRP_ALL extensions), or a Sky_<name>.fits
                    file produced by GetSky_from_CFrame_sum.py's extraction
                    mode (single FLUX extension, DRP_ALL with ra/dec/
                    obstime/tel columns).  Auto-detected via
                    GetSkyCont._is_sky_file.

    Options::

        -ext EXT         FLUX | SKY_EAST | SKY_WEST -- may be given more
                         than once.  XCframe/XSFrame mode only; ignored
                         (with a warning) for Sky_<name>.fits input.
                         Default: all three (FLUX, SKY_EAST, SKY_WEST).
        -engine E        'local' (default) | 'remote' | 'auto' -- which
                         EsoSkyObs.run_sky_obs engine to fetch the ESO
                         MOON/ZODI prediction from (PALACE itself has no
                         engine concept -- see Notes).  Defaults to
                         'local', same as SkyObsESOCompare.py, for the
                         same reason (mixing local/remote-fallback rows
                         would contaminate a judgement of model quality) --
                         but see Notes for why 'remote' is the more
                         independent choice for this particular script.
        -mask FILE       clean-pixel mask FITS file (WAVE, MASK
                         extensions, MASK=1 for line-free/clean pixels);
                         default data/sky_mask.fits, the same mask
                         SkyObsESOCompare.py/GetSkyCont.py/SkySepESO.py use.
        -delta N         process every N-th row; useful for quick tests
                         (default: 1 = all rows)
        -out ROOT        output filename root; default is
                         <stem>_palacecompare

Description:

    For each (extension, row) spectrum requested, exactly mirrors
    SkyObsESOCompare.py's own Description, with one difference: two
    independent model fetches instead of one.

    1. RA/Dec/obstime are read from DRP_ALL, same columns/convention as
       SkyObsESOCompare.py.

    2. MOON, ZODI: SkySepESO._get_sky_model(ra, dec, obstime, engine=...),
       reused unchanged -- only the MOON/ZODI columns of the returned
       table are used; its DIFFUSE/LINES/FLUX are ignored (PALACE
       supplies those instead).

    3. DIFFUSE, LINES: PalaceObs.build_cline/_call_model_linecont/
       rayleigh_per_nm_to_flux, called directly rather than going through
       PalaceObs.do_one() as a whole (see Notes) -- isatm=False gives
       DIFFUSE (matches ESO's above-atmosphere DIFFUSE convention);
       isatm=True gives LINES_LCO (matches ESO's ground-level LINES
       convention -- ESO's own LINES is, as of this writing, NOT divided
       by transmission while its DIFFUSE is, so PALACE's LINES_LCO/
       DIFFUSE, not LINES/DIFFUSE_LCO, are the pair that actually matches
       ESO's own convention like-for-like).

    4. cont0 = MOON + ZODI + DIFFUSE (all three now real physical
       erg/s/cm^2/Angstrom, so combining them is dimensionally sound, same
       as SkyObsESOCompare.py's own cont0).  From here on -- the 3-param
       separation-stage continuum fit (CONT3), the 1-param final continuum
       fit (CONT1/CONT_SCALE), both line fits (SEP against a CONT3
       baseline, final against a CONT1 baseline), NMAD/RMS, NOISE_PROXY,
       FIT_QUALITY_CONT/LINE, and the per-region LINE_FRAC_* metrics --
       are computed by SkyObsESOCompare.py's own _fit_continuum_3param/
       _fit_continuum_1param/_fit_line_scale/_overall_stats/_noise_proxy/
       _find_segments/_line_region_fail_fraction, imported and reused
       completely unchanged (see Notes: none of them are ESO-specific,
       they only operate on whatever arrays they're given).  MOON_FLUX/
       ZODI_FLUX/DIFFUSE_FLUX/LINE_FLUX are recorded the same way too, so
       SkyObsESO_analysis.add_flux_error_columns's FLUX_ERROR_CONT/
       FLUX_ERROR_LINE formulas work unchanged on this script's output as
       well.

    Output DRP_ALL column set, QA flag bits, and the ``SCI_``/``SKYE_``/
    ``SKYW_`` per-extension prefixing convention: identical to
    SkyObsESOCompare.py's own (see its Description for the full list) --
    deliberately, so every downstream SkyObsESO_analysis.py routine runs
    against this script's output with zero changes.

Primary routines::

    one_row   fit a single (wave, flux) spectrum against ESO's MOON/ZODI
              and PALACE's DIFFUSE/LINES
    do_all    batch-process an input file and write the output FITS file

Notes::

    Not calling PalaceObs.do_one() as a whole: it also computes a full
    9-species breakdown and both isatm variants of the combined FLUX,
    none of which this script needs, and writes+reads+deletes a FITS
    file per call besides.  Calling build_cline/_call_model_linecont
    directly (both already module-level, reusable functions in
    PalaceObs.py) skips all of that -- 2 PALACE calls per row here
    (one per isatm) versus do_one()'s 4 internal calls, and no disk
    round-trip at all.

    Cost: two independent model fetches per row now (ESO's
    _get_sky_model plus PALACE's 2x _call_model_linecont) -- roughly on
    par with or somewhat more than SkyObsESOCompare.py's own per-row
    cost alone.  Use -delta for exploratory runs, same guidance as
    SkyObsESOCompare.py/SkyObsESOStack.py.

    On -engine and independence: a prior investigation (see project
    history) found that this codebase's local ESO Sky Model install has
    its airglow-continuum file (acontname) swapped to palace_cont.dat --
    i.e. local-engine ESO's own DIFFUSE is PALACE's continuum table,
    rescaled, not an independent calculation.  This script never uses
    ESO's DIFFUSE at all (PALACE supplies DIFFUSE directly), so that
    substitution doesn't contaminate THIS script's continuum fit -- but
    it does mean ESO's MOON/ZODI (computed by separate, untouched
    physical modules) are the only ESO-sourced components here, and they
    stay genuinely independent of PALACE regardless of -engine.  -engine
    remote (the real, unmodified eso.org service) remains available as a
    fully independent cross-check if wanted, orthogonal to this script's
    own PALACE dependence.

    Reuses SkyObsESOCompare._fit_continuum_3param/_fit_continuum_1param/
    _fit_line_scale/_overall_stats/_noise_proxy/_find_segments/
    _line_region_fail_fraction/_write_qa_report/_EXT_RADEC_COLS/
    _ALL_EXTS/_EXT_COL_PREFIX and the QA_* flag constants directly rather
    than duplicating them -- importing underscore-prefixed helpers across
    scripts has precedent elsewhere in this codebase (SkySepPalace.py
    imports _load_ref_lsf/_interp_lsf_to_wave from SkySubDev2.py;
    SkySubDev2.py imports _get_decomposer from XSkySepIvan.py; PalaceObs
    fetch functions here mirror SkySepESO.py's own _get_sky_model
    convention).  _fit_one_extension/do_all themselves are NOT reused
    from SkyObsESOCompare.py even though structurally similar -- they
    call one_row by (module-level) name, not as a passed-in parameter, so
    reusing them verbatim would silently call ESO's own one_row instead
    of this script's PALACE-aware one.

    No -model option in this version (unlike SkyObsESOCompare.py) --
    there's no single combined "raw model" to write that would mean the
    same thing across two different source models; could be added later
    as two separate optional extensions if wanted.

History::

    260730 ksl  Coding begun: an analog of SkyObsESOCompare.py using
        PALACE for DIFFUSE/LINES and ESO for MOON/ZODI (PALACE predicts
        no Moon/Zodiacal light at all).  Reuses SkyObsESOCompare.py's
        fit-quality/NMAD/RMS/LINE_FRAC machinery unchanged; only the
        model-fetch step and one_row/do_all are new.

'''

import sys
import os
from pathlib import Path

# ensure py_progs siblings are importable when running directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from astropy.io import fits
from astropy.table import Table

from SkySepESO import _get_sky_model
from GetSkyCont import load_mask, _interp_mask_to_wave, _is_sky_file
import PalaceObs

from SkyObsESOCompare import (
    QA_NANDATA, QA_MODELFAIL, QA_FAILED, QA_ZEROCONT, QA_ZEROLINE,
    _fit_continuum_3param, _fit_continuum_1param, _fit_line_scale,
    _overall_stats, _noise_proxy, _find_segments, _line_region_fail_fraction,
    _write_qa_report, _EXT_RADEC_COLS, _ALL_EXTS, _EXT_COL_PREFIX,
    _LINE_FRAC_SIGMA, _LINE_FRAC_ABS_THRESH,
)

DEFAULT_MASK_FILE = Path(__file__).resolve().parent.parent / 'data' / 'sky_mask.fits'

_USAGE = '''Usage:
  SkyObsPalaceCompare.py [-ext EXT [-ext EXT ...]]
                         [-engine local|remote|auto] [-mask FILE]
                         [-delta N] [-out ROOT] filename

Arguments:
  filename    XCframe/XSFrame summary file, or Sky_<name>.fits from
             GetSky_from_CFrame_sum.py (auto-detected)

Options:
  -ext EXT         FLUX | SKY_EAST | SKY_WEST -- may repeat.
                   Default: all three. Ignored for Sky_<name>.fits input.
  -engine E        local (default) | remote | auto -- for ESO's MOON/ZODI
                   fetch only (PALACE has no engine concept)
  -mask FILE       clean-pixel mask FITS file (default: data/sky_mask.fits)
  -delta N         step size through rows for quick tests (default: 1)
  -out ROOT        output filename root (default: <stem>_palacecompare)
'''

# Padding (Angstroms) at each end of the PALACE request, clipped to
# PALACE's own [0.30, 2.50] micron domain -- mirrors
# SkySepPalace._native_wave_range's logic inline (a 2-line computation,
# not worth importing SkySepPalace.py's own heavier dependency chain for).
_WAVE_PAD_ANG = 5.0


def _to_flux(spec, wave):
    '''
    Convert one PALACE line/continuum spectrum table (lam in micron, flux
    in Rayleighs/nm -- the format _call_model_linecont returns) to
    erg/s/cm^2/Angstrom, interpolated onto wave.  Same pattern
    PalaceObs.do_one()'s own to_flux closure uses.  Returns zeros if spec
    is empty (no contribution in this wavelength range).
    '''
    if len(spec) == 0:
        return np.zeros_like(wave)
    spec_wave_ang = np.array(spec['lam'], dtype=float) * 1e4
    spec_flux = PalaceObs.rayleigh_per_nm_to_flux(
        np.array(spec['flux'], dtype=float), spec_wave_ang)
    return np.interp(wave, spec_wave_ang, spec_flux, left=0.0, right=0.0)


def _get_palace_model(ra, dec, obstime, wave):
    '''
    Fetch PALACE's DIFFUSE (airglow continuum, isatm=False -- matches
    ESO's above-atmosphere DIFFUSE convention) and LINES_LCO (isatm=True
    -- matches ESO's ground-level LINES convention; see module Notes on
    the isatm/LCO naming gotcha), both interpolated onto wave.

    Calls PalaceObs.build_cline/_call_model_linecont/
    rayleigh_per_nm_to_flux directly rather than PalaceObs.do_one() as a
    whole (see module Notes) -- 2 PALACE calls total, no disk round-trip.

    Returns (diffuse, lines, err) -- err is '' on success, else a
    description; diffuse/lines are None on failure (e.g. invalid PALACE
    parameter list for this geometry, such as a source below the
    horizon -- same failure mode ESO's own fetch has).
    '''
    lammin = max(0.30, (float(wave.min()) - _WAVE_PAD_ANG) / 1e4)
    lammax = min(2.50, (float(wave.max()) + _WAVE_PAD_ANG) / 1e4)

    cline = PalaceObs.build_cline(ra, dec, obstime, species='all',
                                  lammin=lammin, lammax=lammax,
                                  dlam=PalaceObs.DEFAULT_DLAM,
                                  resol=PalaceObs.DEFAULT_RESOL, isatm=False)
    _linspec, contspec = PalaceObs._call_model_linecont(cline)
    if contspec is None:
        return None, None, ('Invalid PALACE parameter list (isatm=False) '
                            'ra=%s dec=%s obstime=%s' % (ra, dec, obstime))
    diffuse = _to_flux(contspec, wave)

    cline_lco = PalaceObs.build_cline(ra, dec, obstime, species='all',
                                      lammin=lammin, lammax=lammax,
                                      dlam=PalaceObs.DEFAULT_DLAM,
                                      resol=PalaceObs.DEFAULT_RESOL, isatm=True)
    linspec_lco, _contspec_lco = PalaceObs._call_model_linecont(cline_lco)
    if linspec_lco is None:
        return None, None, ('Invalid PALACE parameter list (isatm=True) '
                            'ra=%s dec=%s obstime=%s' % (ra, dec, obstime))
    lines = _to_flux(linspec_lco, wave)

    return diffuse, lines, ''


# ──────────────────────────────────────────────────────────────
# Per-spectrum fit
# ──────────────────────────────────────────────────────────────

def one_row(wave, flux, ra, dec, obstime, clean, engine='local', verbose=False):
    '''
    Fit ESO's MOON/ZODI and PALACE's DIFFUSE/LINES to a single observed
    spectrum, over disjoint clean/line-affected pixel sets -- a direct
    analog of SkyObsESOCompare.one_row.  See module Description.

    Returns (result, qa_flags, coeffs, err_msg) -- same shape as
    SkyObsESOCompare.one_row.  coeffs has the identical key set (see that
    function's own docstring for the full list).
    '''
    qa_flags = 0
    flux = np.asarray(flux, dtype=float)

    if not np.all(np.isfinite(flux)):
        qa_flags |= QA_NANDATA
        flux = np.nan_to_num(flux, nan=0.0, posinf=0.0, neginf=0.0)

    eso_tab, eso_err = _get_sky_model(ra, dec, obstime, engine=engine, verbose=verbose)
    if eso_tab is None:
        return None, QA_MODELFAIL, {}, eso_err

    moon = np.interp(wave, eso_tab['WAVE'], eso_tab['MOON'])
    zodi = np.interp(wave, eso_tab['WAVE'], eso_tab['ZODI'])

    diffuse, lines, palace_err = _get_palace_model(ra, dec, obstime, wave)
    if diffuse is None:
        return None, QA_MODELFAIL, {}, palace_err

    cont0 = moon + zodi + diffuse
    line_mask = ~clean

    if np.dot(cont0[clean], cont0[clean]) == 0:
        qa_flags |= QA_ZEROCONT
    if np.dot(lines[line_mask], lines[line_mask]) == 0:
        qa_flags |= QA_ZEROLINE

    cont3, _coef3, _err3 = _fit_continuum_3param(flux, clean, moon, zodi, diffuse)
    k, k_err, cont1 = _fit_continuum_1param(flux, clean, cont0)

    denom0 = np.sum(cont0[clean])
    sep_cont_scale = float(np.sum(cont3[clean]) / denom0) if denom0 != 0 else np.nan

    sep_line_input = flux - cont3
    sep_d, sep_d_err = _fit_line_scale(sep_line_input, line_mask, lines)
    sep_lines_scaled = sep_d * lines

    final_line_input = flux - cont1
    d, d_err = _fit_line_scale(final_line_input, line_mask, lines)
    lines_scaled = d * lines

    sep_cont_resid = flux - cont3
    cont_resid = flux - cont1
    sep_line_resid = sep_line_input - sep_lines_scaled
    final_line_resid = final_line_input - lines_scaled

    sep_cont_nmad, sep_cont_rms = _overall_stats(wave, sep_cont_resid, clean)
    cont_nmad, cont_rms = _overall_stats(wave, cont_resid, clean)
    sep_line_nmad, sep_line_rms = _overall_stats(wave, sep_line_resid, line_mask)
    line_nmad, line_rms = _overall_stats(wave, final_line_resid, line_mask)

    noise_proxy = _noise_proxy(flux, clean)
    with np.errstate(divide='ignore', invalid='ignore'):
        sep_fit_quality_cont = float(sep_cont_nmad / noise_proxy)
        fit_quality_cont = float(cont_nmad / noise_proxy)
        sep_fit_quality_line = float(sep_line_nmad / noise_proxy)
        fit_quality_line = float(line_nmad / noise_proxy)

    line_segments = _find_segments(line_mask)
    sig_thresh = _LINE_FRAC_SIGMA * noise_proxy
    sep_line_frac_10sig = _line_region_fail_fraction(sep_line_resid, line_segments, sig_thresh)
    line_frac_10sig = _line_region_fail_fraction(final_line_resid, line_segments, sig_thresh)
    abs1, abs2 = _LINE_FRAC_ABS_THRESH
    sep_line_frac_1e15 = _line_region_fail_fraction(sep_line_resid, line_segments, abs1)
    line_frac_1e15 = _line_region_fail_fraction(final_line_resid, line_segments, abs1)
    sep_line_frac_1e14 = _line_region_fail_fraction(sep_line_resid, line_segments, abs2)
    line_frac_1e14 = _line_region_fail_fraction(final_line_resid, line_segments, abs2)

    moon_flux = float(np.mean(moon[clean]))
    zodi_flux = float(np.mean(zodi[clean]))
    diffuse_flux = float(np.mean(diffuse[clean]))
    line_flux = float(np.mean(lines[line_mask]))

    model_unscaled = cont0 + lines
    result = dict(MODEL=model_unscaled)

    coeffs = dict(
        cont_scale=k, cont_scale_err=k_err,
        sep_cont_scale=sep_cont_scale,
        sep_cont_nmad=sep_cont_nmad, sep_cont_rms=sep_cont_rms,
        cont_nmad=cont_nmad, cont_rms=cont_rms,
        noise_proxy=noise_proxy,
        sep_fit_quality_cont=sep_fit_quality_cont,
        fit_quality_cont=fit_quality_cont,
        sep_line_scale=sep_d, sep_line_scale_err=sep_d_err,
        sep_line_nmad=sep_line_nmad, sep_line_rms=sep_line_rms,
        sep_fit_quality_line=sep_fit_quality_line,
        sep_line_frac_10sig=sep_line_frac_10sig,
        sep_line_frac_1e15=sep_line_frac_1e15, sep_line_frac_1e14=sep_line_frac_1e14,
        line_scale=d, line_scale_err=d_err,
        line_nmad=line_nmad, line_rms=line_rms,
        fit_quality_line=fit_quality_line,
        line_frac_10sig=line_frac_10sig,
        line_frac_1e15=line_frac_1e15, line_frac_1e14=line_frac_1e14,
        moon_flux=moon_flux, zodi_flux=zodi_flux, diffuse_flux=diffuse_flux,
        line_flux=line_flux,
    )

    return result, qa_flags, coeffs, ''


# ──────────────────────────────────────────────────────────────
# Batch processing
# ──────────────────────────────────────────────────────────────

def _fit_one_extension(x, drp_all, wave, ext, rows, clean, engine, nan_spectrum):
    '''
    Fit every row in rows for one XCframe/XSFrame extension against ESO's
    MOON/ZODI and PALACE's DIFFUSE/LINES.  Same shape as
    SkyObsESOCompare._fit_one_extension, minus the -model plumbing (not
    supported here -- see module Notes) -- kept as a separate copy rather
    than importing that function, since it calls one_row by module-level
    name and would silently call ESO's own one_row instead of this
    script's PALACE-aware one.
    '''
    ra_col, dec_col = _EXT_RADEC_COLS[ext]
    final_flux = []
    qa_flags_list, coeffs_list, err_msgs = [], [], []

    for n, row in enumerate(rows):
        try:
            flux = np.array(x[ext].data[row], dtype=float)
            ra = float(drp_all[ra_col][row])
            dec = float(drp_all[dec_col][row])
            obstime = drp_all['obstime'][row]
            result, row_flags, coeffs, err_msg = one_row(
                wave, flux, ra, dec, obstime, clean, engine=engine)
            flux_out = flux
        except Exception as e:
            print('%s row %d: exception (%s)' % (ext, row, e))
            result = None
            row_flags = 0
            coeffs = {}
            err_msg = 'exception: %s' % e
            flux_out = nan_spectrum

        final_flux.append(np.asarray(flux_out, dtype=float))
        if result is None:
            row_flags |= QA_FAILED

        qa_flags_list.append(row_flags)
        coeffs_list.append(coeffs)
        err_msgs.append(err_msg)

        if (n + 1) % 100 == 0 or (n + 1) == len(rows):
            print('  Completed %6d of %d' % (n + 1, len(rows)))

    return dict(flux=final_flux, qa=qa_flags_list, coeffs=coeffs_list, err=err_msgs)


def do_all(filename, exts=None, engine='local', mask_file=DEFAULT_MASK_FILE,
          idelta=1, outroot=''):
    '''
    Fit ESO's MOON/ZODI and PALACE's DIFFUSE/LINES to every requested
    spectrum in an input file and write the results.  Same shape as
    SkyObsESOCompare.do_all, minus the -model option (see module Notes).

    Parameters
    ----------
    filename : str
        XCframe/XSFrame summary file, or a Sky_<name>.fits file.
    exts : list of str or None
        Extensions to process (XCframe/XSFrame mode only); None means all
        three (FLUX, SKY_EAST, SKY_WEST).
    engine : str
        'local' (default), 'remote', or 'auto' -- for ESO's MOON/ZODI
        fetch only.
    mask_file : str or Path
        Clean-pixel mask FITS file; default data/sky_mask.fits.
    idelta : int
        Row step size (1 = all rows).
    outroot : str
        Output filename root; defaults to <stem>_palacecompare.
    '''
    x = fits.open(filename)
    wave = np.array(x['WAVE'].data, dtype=float)
    drp_all = Table(x['DRP_ALL'].data)
    rows = list(range(0, len(drp_all), idelta))

    sky_mode = _is_sky_file(filename)
    print('Input mode: %s' % ('Sky_<name>.fits' if sky_mode else 'XCframe/XSFrame'))

    mask_wave, mask_arr = load_mask(mask_file)
    clean = _interp_mask_to_wave(mask_wave, mask_arr, wave)
    print('Loaded mask %s: %d/%d clean pixels' % (mask_file, clean.sum(), len(clean)))

    n_wave = len(wave)
    nan_spectrum = np.full(n_wave, np.nan)

    if outroot == '':
        stem = os.path.splitext(os.path.basename(filename))[0]
        outroot = '%s_palacecompare' % stem

    if sky_mode:
        if exts:
            print("Warning: -ext is ignored for Sky_<name>.fits input "
                 "(spectrum type is already fixed per row via 'tel').")
        tel = (np.array(drp_all['tel'], dtype=str) if 'tel' in drp_all.colnames
              else np.full(len(drp_all), 'FLUX'))
        print('Processing %d spectra ...' % len(rows))

        final_flux = []
        ext_list = []
        qa_flags_list, coeffs_list, err_msgs = [], [], []

        for n, row in enumerate(rows):
            try:
                flux = np.array(x['FLUX'].data[row], dtype=float)
                ra = float(drp_all['ra'][row])
                dec = float(drp_all['dec'][row])
                obstime = drp_all['obstime'][row]
                result, row_flags, coeffs, err_msg = one_row(
                    wave, flux, ra, dec, obstime, clean, engine=engine)
                flux_out = flux
            except Exception as e:
                print('row %d: exception (%s)' % (row, e))
                result = None
                row_flags = 0
                coeffs = {}
                err_msg = 'exception: %s' % e
                flux_out = nan_spectrum

            final_flux.append(np.asarray(flux_out, dtype=float))
            if result is None:
                row_flags |= QA_FAILED

            ext_list.append(str(tel[row]).strip())
            qa_flags_list.append(row_flags)
            coeffs_list.append(coeffs)
            err_msgs.append(err_msg)

            if (n + 1) % 100 == 0 or (n + 1) == len(rows):
                print('Completed %6d of %d' % (n + 1, len(rows)))

        n_failed    = sum(1 for f in qa_flags_list if f & QA_FAILED)
        n_modelfail = sum(1 for f in qa_flags_list if f & QA_MODELFAIL)
        n_warned    = sum(1 for f in qa_flags_list
                          if f != 0 and not (f & (QA_FAILED | QA_MODELFAIL)))
        print('\nProcessed %d spectra: %d failed (NaN fill), %d model-fetch '
             'failures, %d with warnings'
             % (len(rows), n_failed, n_modelfail, n_warned))

        _write_qa_report(
            [(ext_list[j], rows[j], qa_flags_list[j], err_msgs[j])
             for j in range(len(rows)) if qa_flags_list[j] != 0],
            outroot)

        spectral_hdus = [fits.ImageHDU(data=np.array(final_flux), name='FLUX')]
        ext_hdr_val = 'ALL'

        xtab = drp_all[np.array(rows, dtype=int)].copy()
        xtab['EXT'] = np.array(ext_list)
        xtab['ROW'] = np.array(rows, dtype=np.int32)

        all_keys = sorted({k for d in coeffs_list for k in d})
        for key in all_keys:
            xtab[key.upper()] = np.array(
                [d.get(key, np.nan) for d in coeffs_list], dtype=np.float32)
        xtab['QA_FLAGS']  = np.array(qa_flags_list, dtype=np.int32)
        xtab['ERROR_MSG'] = np.array([m[:200] for m in err_msgs])

    else:
        use_exts = exts if exts else list(_ALL_EXTS)
        for e in use_exts:
            if e not in _EXT_RADEC_COLS:
                raise ValueError('unknown extension "%s"; must be one of %s'
                                 % (e, sorted(_EXT_RADEC_COLS)))

        print('Processing %d exposures x %d extension(s) (%s) ...'
             % (len(rows), len(use_exts), ', '.join(use_exts)))

        per_ext = {}
        for e in use_exts:
            print('Extension %s:' % e)
            per_ext[e] = _fit_one_extension(x, drp_all, wave, e, rows, clean, engine, nan_spectrum)

        all_qa = [f for e in use_exts for f in per_ext[e]['qa']]
        n_failed    = sum(1 for f in all_qa if f & QA_FAILED)
        n_modelfail = sum(1 for f in all_qa if f & QA_MODELFAIL)
        n_warned    = sum(1 for f in all_qa
                          if f != 0 and not (f & (QA_FAILED | QA_MODELFAIL)))
        print('\nProcessed %d exposures x %d extension(s): %d failed (NaN fill), '
             '%d model-fetch failures, %d with warnings'
             % (len(rows), len(use_exts), n_failed, n_modelfail, n_warned))

        _write_qa_report(
            [(e, rows[j], per_ext[e]['qa'][j], per_ext[e]['err'][j])
             for e in use_exts for j in range(len(rows))
             if per_ext[e]['qa'][j] != 0],
            outroot)

        spectral_hdus = []
        for e in use_exts:
            d = per_ext[e]
            spectral_hdus.append(fits.ImageHDU(data=np.array(d['flux']), name=e))
        ext_hdr_val = ','.join(use_exts)

        xtab = drp_all[np.array(rows, dtype=int)].copy()
        xtab['ROW'] = np.array(rows, dtype=np.int32)

        for e in use_exts:
            d = per_ext[e]
            p = _EXT_COL_PREFIX[e]
            all_keys = sorted({k for c in d['coeffs'] for k in c})
            for key in all_keys:
                xtab['%s_%s' % (p, key.upper())] = np.array(
                    [c.get(key, np.nan) for c in d['coeffs']], dtype=np.float32)
            xtab['%s_QA_FLAGS' % p]  = np.array(d['qa'], dtype=np.int32)
            xtab['%s_ERROR_MSG' % p] = np.array([m[:200] for m in d['err']])

    hdr = fits.Header()
    hdr['TITLE']  = 'SkyObsPalaceCompare'
    hdr['INPUT']  = str(filename)
    hdr['ENGINE'] = engine
    hdr['MASK']   = str(mask_file)
    hdr['EXT']    = ext_hdr_val

    hdu1 = fits.PrimaryHDU(header=hdr)
    hdu2 = fits.ImageHDU(data=wave, name='WAVE')
    hdu_drp = fits.BinTableHDU(xtab, name='DRP_ALL')

    hdul = fits.HDUList([hdu1, hdu2] + spectral_hdus + [hdu_drp])

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

    exts      = []
    engine    = 'local'
    mask_file = DEFAULT_MASK_FILE
    idelta    = 1
    outroot   = ''
    filename  = None

    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == '-ext':
            i += 1
            exts.append(argv[i])
        elif arg == '-engine':
            i += 1
            engine = argv[i]
        elif arg == '-mask':
            i += 1
            mask_file = argv[i]
        elif arg == '-delta':
            i += 1
            idelta = int(argv[i])
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
                sys.exit(1)
            filename = arg
        i += 1

    if filename is None:
        print('Error: filename is required')
        print(_USAGE)
        sys.exit(1)

    do_all(filename, exts=exts if exts else None, engine=engine,
          mask_file=mask_file, idelta=idelta, outroot=outroot)
