#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Evaluate how closely real LVM spectra (science and/or sky-telescope)
    match the predictions of the ESO Sky Model, by separating each
    observed spectrum's continuum from its lines and fitting independent
    rescalings of the ESO model's continuum and line templates to each
    piece.  Unlike EsoSkyFit.py (which fits all four ESO templates
    jointly, over all pixels), the continuum and line rescalings here are
    fit over disjoint pixel sets (continuum-only vs line-affected, from
    the same sky_mask.fits clean-pixel mask GetSkyCont.py uses), so a
    flux-calibration error shows up as a genuine continuum-scale outlier
    with a real residual, not something a joint fit could quietly absorb
    by shifting weight onto the line term.

    This is a diagnostics tool, not a subtraction tool: it produces
    per-spectrum numbers (scale factors, formal errors, per-arm residual
    statistics) intended to be used later -- once their distributions are
    understood -- to define an actual data-quality flag.  It does not
    define or apply any such flag itself.

Command line usage (if any):

    usage: SkyObsESOCompare.py [-ext EXT [-ext EXT ...]]
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
                         (with a warning) for Sky_<name>.fits input, which
                         has only one spectrum type per row already tagged
                         via its 'tel' column.  Default: all three
                         (FLUX, SKY_EAST, SKY_WEST).
        -engine E        'local' (default) | 'remote' | 'auto' -- which
                         EsoSkyObs.run_sky_obs engine to fetch the ESO
                         model from.  Defaults to 'local' (not SkySepESO.py's
                         'auto') for the same reason as EsoSkyFit.py: mixing
                         local and remote-fallback rows would contaminate a
                         judgement of model quality.
        -mask FILE       palace_mask FITS file (WAVE, MASK extensions,
                         MASK=1 for line-free/clean pixels); default
                         data/sky_mask.fits, the mask already used
                         throughout this codebase (GetSkyCont.py,
                         SkySepESO.py, SkySubDev2.py, ...).
        -delta N         process every N-th row; useful for quick tests
                         (default: 1 = all rows)
        -out ROOT        output filename root; default is <stem>_esocompare
        -model           also write the raw, UNSCALED ESO model prediction
                         (MOON+ZODI+DIFFUSE+LINES, no fitting applied) as a
                         <EXT>_MODEL extension alongside each observed
                         spectrum.  Off by default: the observed spectra
                         actually analyzed (FLUX/SKY_EAST/SKY_WEST) are
                         always written, but the model prediction roughly
                         doubles output size and DRP_ALL alone already
                         carries every scalar this script computes (see
                         Notes).

Description:

    For each (extension, row) spectrum requested:

    1. RA/Dec (sci_ra/sci_dec for FLUX, skye_ra/skye_dec for SKY_EAST,
       skyw_ra/skyw_dec for SKY_WEST -- or, for Sky_<name>.fits input,
       the already-merged ra/dec/tel columns) and obstime are read from
       DRP_ALL, and the ESO sky model is fetched (SkySepESO._get_sky_model,
       reused unchanged) giving MOON/ZODI/DIFFUSE/LINES templates on the
       model's own wavelength grid, interpolated onto the instrument grid.

    2. Continuum-only pixels are identified from the sky_mask.fits clean
       mask (MASK=1), the same mask GetSkyCont.py restricts its own B-spline
       continuum fit to.  Two independent continuum fits are made against
       the observed flux at those pixels only -- one is the continuum
       SEPARATION stage (currently 3-parameter; whatever replaces it
       later, this is where it plugs in), the other is the FINAL,
       1-parameter fit actually reported as CONT_SCALE::

           SEPARATION-stage fit (best shape match):
               CONT3 = a*MOON + b*ZODI + c*DIFFUSE   (a, b, c >= 0, NNLS)

           FINAL fit (model's own predicted mix, one flux-cal number):
               CONT1 = k * (MOON + ZODI + DIFFUSE)   (k >= 0, closed form)

       CONT3 is the best available continuum estimate.  CONT1's
       coefficient k ("CONT_SCALE") is the number intended for judging
       flux calibration: it is uncontaminated by the shape flexibility of
       the separation-stage fit, so its distribution across many
       exposures under similar conditions is what a later quality flag
       would be built from.  A derived SEP_CONT_SCALE (the flux ratio of
       CONT3 to the model's raw, unscaled MOON+ZODI+DIFFUSE sum, both
       summed over clean pixels) is also recorded, so the two can be
       compared directly: agreement means the model's predicted component
       mix is basically right for that geometry; a large difference means
       the model's mix itself (not just its overall amplitude) is off.

    3. On the line-affected pixels (the complement of the clean mask), a
       single non-negative scale factor for the model's own predicted
       LINES template is fit (closed form) -- TWICE, against two
       different baselines, mirroring the SEP/FINAL continuum split
       above so the same "does the separation stage's extra flexibility
       actually help" question can be asked for lines too::

           SEP fit (baseline = SEPARATION-stage continuum):
               LINES = sep_d * MODEL_LINES   against flux - CONT3

           FINAL fit (baseline = FINAL continuum, the one CONT_SCALE
           describes):
               LINES = d * MODEL_LINES       against flux - CONT1

       This (both of them) is what EsoSkyFit.py's real (post-fix) line
       fit does too -- see its History for why fitting the model's own
       predicted line amplitude, rather than treating "whatever's left
       over" as the line spectrum, is the only way this number means
       anything.

    CONT1/CONT3/LINES themselves are not written to the output -- only
    used internally to get the coefficients and NMAD/RMS values below.
    The optional -model extension is a different thing: the raw,
    UNSCALED sum MOON+ZODI+DIFFUSE+LINES exactly as fetched from the
    model, with none of steps 2-3's fitting applied -- there to let you
    plot the model's own prediction directly against the observed
    spectrum (which is always written) without re-fetching anything from
    the model server, not to exactly reconstruct CONT1/CONT3/LINES
    (which would need the continuum/line split redone).

    4. Two summary statistics -- NMAD (robust, median-based) and RMS --
       are computed for each of four residuals: SEP_CONT (flux - CONT3 at
       clean pixels), CONT (flux - CONT1 at clean pixels, unprefixed
       since it's the FINAL answer), SEP_LINE (the SEP line residual at
       line-affected pixels), and LINE (the FINAL line residual,
       unprefixed) -- using GetSkyCont.arm_continuum_stats reused with a
       single band spanning the full wavelength range.

       CONT3 is a strict superset of CONT1 (CONT1 = CONT3 with a=b=c
       forced) -- literally the same predictors, just a linear
       constraint -- so CONT3's fit, searching a strictly larger feasible
       region that includes CONT1's exact solution, is guaranteed
       SEP_CONT_RMS <= CONT_RMS (RMS is exactly what these least-squares
       fits minimize; verified in practice).  The two LINE fits are NOT
       nested this way, despite both fitting the same single amplitude d
       to residual ~ d*LINES: the residual itself is different data
       (flux - CONT3 vs flux - CONT1), evaluated at line-affected pixels
       that neither continuum fit was ever fit TO -- both continua are
       being extrapolated there, and whether CONT3 or CONT1 extrapolates
       better into that region is a real empirical question with no
       mathematical answer, so SEP_LINE_RMS <= LINE_RMS is NOT guaranteed
       (confirmed empirically: it fails for roughly half of rows).  NMAD
       has no guarantee for either pair (it is not what the fits
       minimize): a more flexible fit typically improves RMS mainly by
       reducing a handful of larger/outlier residuals while barely moving
       the median, so ``SEP_*_NMAD`` can come out larger than the unprefixed
       ``*_NMAD`` even when SEP_CONT_RMS <= CONT_RMS holds.  Keep both: RMS
       to compare "is the separation stage actually better, and by how
       much" (a guaranteed answer for continuum, an empirical one for
       lines), NMAD as a "how big is a typical pixel's residual"
       diagnostic that both stages' resulting FIT_QUALITY happens to use.

       Each NMAD is in absolute flux units, so it's naturally larger for
       brighter fields regardless of how good the shape match actually
       is -- not directly comparable across rows on its own.  NOISE_PROXY
       (a per-row, model-independent noise estimate from consecutive-
       pixel differences within clean pixels: ``1.4826*median(|diff|)/
       sqrt(2)``) normalizes each NMAD into a reduced-chi-squared-like
       FIT_QUALITY ratio that IS comparable across rows::

           SEP_FIT_QUALITY_CONT = SEP_CONT_NMAD / NOISE_PROXY
           FIT_QUALITY_CONT     = CONT_NMAD     / NOISE_PROXY
           SEP_FIT_QUALITY_LINE = SEP_LINE_NMAD / NOISE_PROXY
           FIT_QUALITY_LINE     = LINE_NMAD     / NOISE_PROXY

       ~1 means the residual is consistent with noise alone (an excellent
       shape match); values much greater than 1 mean real, unmodeled
       shape mismatch.  Both LINE ratios necessarily share the
       continuum's NOISE_PROXY, since there's no noise-only region within
       the line-affected pixels themselves to estimate it from
       independently.

    5. The line-affected pixels are also split into maximal contiguous
       runs ("regions" -- one per distinct airglow feature, or per
       cluster of blended ones; ~344 of them for the standard
       sky_mask.fits).  FIT_QUALITY_LINE (and SEP_FIT_QUALITY_LINE) is
       one blended number, dominated by whichever regions have the
       most/brightest pixels (see Notes -- for a typical spectrum, ~10 of
       these ~344 regions account for most of LINE_SCALE's own fit
       weight).  As a complementary, non-dominated view, each region's
       median ``|residual|`` is compared against three thresholds, and the
       FRACTION of regions exceeding each is recorded -- SEP_LINE_FRAC_*
       from the SEP line residual, LINE_FRAC_* (unprefixed) from the
       FINAL one::

           LINE_FRAC_10SIG = fraction of regions with median|residual| > 10*NOISE_PROXY
           LINE_FRAC_1E15  = fraction of regions with median|residual| > 1e-15
           LINE_FRAC_1E14  = fraction of regions with median|residual| > 1e-14

       (and the same three for SEP_LINE_FRAC_*, from the SEP line
       residual).  3 sigma was tried first and rejected: it's saturated
       (70-90% of regions exceed it on essentially every row tested, so
       it doesn't discriminate row to row); 10 sigma showed real spread
       in practice.  The two fixed absolute thresholds don't rescale with
       how bright/noisy a given exposure is, unlike the sigma-based one --
       deliberately, so they answer "is the residual big enough to matter
       regardless of this exposure's own noise level" rather than "is it
       statistically inconsistent with this exposure's own noise."
       Unlike the continuum RMS pair, there is no guarantee SEP_LINE_FRAC_*
       <= LINE_FRAC_* (same reason SEP_LINE_RMS <= LINE_RMS isn't
       guaranteed, and this is a discrete per-region count besides, not
       even a minimized quantity) -- confirmed empirically to go either
       way row to row.

    6. The mean brightness of each UNSCALED model component is recorded --
       MOON_FLUX/ZODI_FLUX/DIFFUSE_FLUX (mean of MOON/ZODI/DIFFUSE over
       clean pixels) and LINE_FLUX (mean of LINES over line-affected
       pixels) -- each over the same footprint its corresponding scale
       factor is fit from.  Mean, not median: the scale fits are least-
       squares/dot-product (flux-weighted), so mean is the matching
       summary statistic, and median would read low for LINE_FLUX
       specifically (most line-affected pixels are wings/gaps between
       narrow line cores, not the cores themselves).  Together with
       CONT_SCALE/LINE_SCALE these give an actual flux-error estimate
       without needing to infer one indirectly from moon_phase/moon_alt::

           FLUX_ERROR_CONT = (CONT_SCALE - 1) * (MOON_FLUX + ZODI_FLUX + DIFFUSE_FLUX)
           FLUX_ERROR_LINE = (LINE_SCALE - 1) * LINE_FLUX

    DRP_ALL is built from the ORIGINAL input file's DRP_ALL rows, and in
    XCframe/XSFrame mode has exactly one row per exposure (per row
    processed, regardless of how many extensions -ext requests) so every
    pre-existing column -- moon_alt, moon_phase, sun_alt, tileid, mjd,
    expnum, sci_amass, sci_ra/skye_ra/skyw_ra, object, etc. -- survives
    into the output unchanged, with just a ROW column added.  This
    script's own results are added as per-extension columns, prefixed
    ``SCI_`` for FLUX, ``SKYE_`` for SKY_EAST, ``SKYW_`` for SKY_WEST -- matching the
    input file's own sci_*/skye_*/skyw_* column-prefix convention -- so a
    single row carries every requested extension's numbers for that one
    exposure side by side, e.g.::

        ROW
        SCI_CONT_SCALE, SCI_CONT_SCALE_ERR, SCI_SEP_CONT_SCALE
        SCI_SEP_CONT_NMAD, SCI_SEP_CONT_RMS, SCI_CONT_NMAD, SCI_CONT_RMS
        SCI_NOISE_PROXY, SCI_SEP_FIT_QUALITY_CONT, SCI_FIT_QUALITY_CONT
        SCI_SEP_LINE_SCALE, SCI_SEP_LINE_SCALE_ERR, SCI_SEP_LINE_NMAD,
            SCI_SEP_LINE_RMS, SCI_SEP_FIT_QUALITY_LINE
        SCI_SEP_LINE_FRAC_10SIG, SCI_SEP_LINE_FRAC_1E15, SCI_SEP_LINE_FRAC_1E14
        SCI_LINE_SCALE, SCI_LINE_SCALE_ERR, SCI_LINE_NMAD, SCI_LINE_RMS,
            SCI_FIT_QUALITY_LINE
        SCI_LINE_FRAC_10SIG, SCI_LINE_FRAC_1E15, SCI_LINE_FRAC_1E14
        SCI_MOON_FLUX, SCI_ZODI_FLUX, SCI_DIFFUSE_FLUX, SCI_LINE_FLUX
        SCI_QA_FLAGS, SCI_ERROR_MSG
        (and the same set again with SKYE_/SKYW_ prefixes for whichever
        of SKY_EAST/SKY_WEST were also requested)

    Sky_<name>.fits mode does not have this per-exposure merge available
    (a physical exposure may contribute a SKY_EAST row, a SKY_WEST row,
    both, or neither, depending on which telescope(s) matched the
    requested sky field -- see GetSky_from_CFrame_sum.py), so it keeps one
    row per input row instead, with unprefixed columns (CONT_SCALE, ...,
    QA_FLAGS, ERROR_MSG) and an EXT column (copied from the existing 'tel'
    column) identifying SKY_EAST vs SKY_WEST per row.

    The column set has changed several times already based on what a
    first large-scale run actually needed -- see History: a first,
    richer version (per-extension MOON/ZODI/DIFFUSE coefficients and
    their formal errors, plus a full B/R/Z x med/nmad/rms/skew breakdown
    for three separate residuals) was cut down to just CONT_SCALE/
    LINE_SCALE/NMAD to see whether the comparison was useful at all; once
    it clearly was, MOON_FLUX/ZODI_FLUX/DIFFUSE_FLUX/LINE_FLUX and the
    FIT_QUALITY ratios above were added back in because specific
    follow-up questions (an actual flux-error estimate; whether the
    model's shape, not just its scale, matches the data) needed them; and
    the ``SEP_``/unprefixed naming (dropping "FINAL", extending "SEP" to
    SEP_CONT_SCALE and to a full second SEP line fit) plus the RMS
    columns were added once NMAD by itself turned out not to reliably
    answer "is the separation stage actually better" (see one_row's
    _overall_stats docstring for why).

    QA flag bits (stored in DRP_ALL['QA_FLAGS'] for Sky_<name>.fits mode,
    or DRP_ALL['<PREFIX>_QA_FLAGS'] per requested extension in XCframe/
    XSFrame mode)::

        0x01  NANDATA    NaN/inf found in the input flux
        0x02  MODELFAIL  the ESO sky model could not be fetched for this
                         row with the requested -engine
        0x04  FAILED     row raised an exception; spectra filled with NaN
                         (always set alongside MODELFAIL)
        0x08  ZEROCONT   model's continuum (MOON+ZODI+DIFFUSE) is exactly
                         zero at every clean pixel; CONT_SCALE unreliable
        0x10  ZEROLINE   model's LINES template is exactly zero at every
                         line-affected pixel; LINE_SCALE unreliable

    These are technical-failure flags only -- there is deliberately no
    "this row looks like bad calibration" flag here; see Notes.

Primary routines::

    one_row   fit a single (wave, flux) spectrum against the ESO model
    do_all    batch-process an input file and write the output FITS file

Notes::

    -model and output size: the observed spectra actually analyzed
    (FLUX/SKY_EAST/SKY_WEST, whichever were requested) are always written
    -- that's the point of the tool.  Without -model (the default), that's
    all the per-pixel data in the output; WAVE + those spectra + DRP_ALL.
    -model adds the raw unscaled ESO model prediction alongside each one
    (<EXT>_MODEL), roughly doubling output size: a full, all-three-
    extension run over an ~14000-row XCframe file is roughly 14000 rows x
    3 extensions x 12401 pixels x 4 bytes (float32) per array -- a few GB
    for the observed spectra alone, another few GB if -model is added --
    versus DRP_ALL's roughly hundred scalar columns at a few tens of MB.

    Output spectral extension names: in XCframe/XSFrame mode, each
    requested source extension gets its own set of output image
    extensions -- e.g. the default of all three produces FLUX/FLUX_CONT1/
    FLUX_CONT3/FLUX_LINES/FLUX_RESID, SKY_EAST/SKY_EAST_CONT1/..., and
    SKY_WEST/SKY_WEST_CONT1/..., each holding exactly one spectrum per
    exposure processed (-ext SKY_EAST alone produces just the SKY_EAST
    set).  Row j of every one of these image extensions and row j of
    DRP_ALL all correspond to the same exposure, since DRP_ALL has one row
    per exposure too (see Description).  Sky_<name>.fits input keeps a
    single generic FLUX/CONT1/CONT3/LINES/RESID block instead, since one
    such file can mix SKY_EAST/SKY_WEST rows at arbitrary row order (see
    its own 'tel' column) rather than in ext-grouped blocks.

    No population-level outlier/quality flag is computed or applied here.
    CONT_SCALE, SEP_CONT_SCALE, LINE_SCALE, and the NMAD/RMS residual
    summaries are numbers only; a companion _eval.py script (following
    this codebase's eval_sky.py / SkySub_eval.py / GetSkyCont_eval.py
    convention of separating batch-fit scripts from their plotting/
    threshold-setting companions) would define an actual quality flag once
    the real distributions of these numbers have been examined.

    FLUX rows carry real astrophysical source flux on top of sky, not just
    sky -- a real target's own continuum/line shape at clean/line pixels
    will bias CONT_SCALE/LINE_SCALE regardless of calibration or ESO model
    quality.  This is included by default (rather than defaulting to
    SKY_EAST/SKY_WEST only) specifically so FLUX's numbers can be compared
    against SKY_EAST/SKY_WEST for the same exposure.

    Each row requires exactly one ESO sky model fetch (as EsoSkyFit.py;
    compare SkySepESO.py's two or three per row for real sky subtraction).
    Processing all three extensions of a full XCframe file is therefore
    three times EsoSkyFit.py's already-slow per-row cost -- use -delta or
    -ext for exploratory runs (see SkySepESO.py's Notes for the same
    caveat).

    Reuses SkySepESO._get_sky_model and GetSkyCont.load_mask/
    _interp_mask_to_wave/arm_continuum_stats/_is_sky_file directly rather
    than duplicating them -- importing an underscore-prefixed helper across
    scripts has precedent elsewhere in this codebase (SkySepPalace.py
    imports _load_ref_lsf/_interp_lsf_to_wave from SkySubDev2.py;
    SkySubDev2.py imports _get_decomposer from XSkySepIvan.py; EsoSkyFit.py
    imports _get_sky_model from SkySepESO.py).

History::

    260727 ksl  Coding begun: fits the ESO sky model's continuum (a
        3-parameter NNLS separation-stage fit CONT3, and the model's own
        1-parameter predicted-mix fit CONT1) and lines (a single
        amplitude fit to the model's own LINES template) to each
        requested extension (FLUX/SKY_EAST/SKY_WEST) of an XCframe/
        XSFrame summary file, or to a Sky_<name>.fits file from
        GetSky_from_CFrame_sum.py.  DRP_ALL output has one row per
        exposure (not per extension), with results as SCI_/SKYE_/SKYW_-
        prefixed columns matching the input file's own convention.
        Observed spectra are always written; the raw unscaled ESO model
        prediction (<EXT>_MODEL) is optional (-model).  ESO-model-fetch
        progress (SkySepESO._get_sky_model, given a new verbose=True/
        False switch for this) defaults to off, since batch-processing
        many rows is the whole point; the same diagnostic text is still
        captured into ERROR_MSG/<outroot>_errors.txt on failure.
    260728 ksl  Every scale/fit-quality quantity now comes in two
        versions: SEP_<name> from the 3-parameter separation-stage fit
        (CONT3, or a line fit against a flux-CONT3 baseline), and
        unprefixed <name> from the FINAL fit actually reported (CONT1/
        CONT_SCALE, or a line fit against a flux-CONT1 baseline) -- a
        full second line fit, not just a relabeling, since changing the
        line fit's baseline changes what it measures.  For both
        continuum and lines: NMAD and RMS of the residual (RMS is what
        the least-squares fits actually minimize, so SEP_CONT_RMS <=
        CONT_RMS is a real guaranteed relationship, since CONT1 is a
        literal constrained special case of CONT3; NMAD has no such
        guarantee for either, and SEP_LINE_RMS <= LINE_RMS isn't
        guaranteed at all, since the two line fits use different
        residual baselines rather than being nested the way the two
        continuum fits are); NOISE_PROXY (a per-row, model-independent
        noise estimate from consecutive-pixel differences) turning each
        NMAD into a reduced-chi-squared-like FIT_QUALITY ratio comparable
        across rows regardless of brightness; MOON_FLUX/ZODI_FLUX/
        DIFFUSE_FLUX/LINE_FLUX (mean unscaled-component brightness) so
        CONT_SCALE/LINE_SCALE can be turned into an actual flux-error
        estimate instead of an indirect moon_phase/moon_alt proxy; and
        LINE_FRAC_10SIG/_1E15/_1E14 (fraction of ~344 line-affected
        regions whose median residual exceeds a threshold), added after
        finding LINE_SCALE itself is dominated by a handful of the
        brightest regions ([OI]5577 and a few NIR OH bandheads, not
        [OI]6300) -- FRAC_* gives every region an equal vote regardless
        of brightness or width.

'''

import sys
import os
from pathlib import Path

# ensure py_progs siblings are importable when running directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import nnls

from SkySepESO import _get_sky_model
from GetSkyCont import load_mask, _interp_mask_to_wave, arm_continuum_stats, _is_sky_file

# ──────────────────────────────────────────────────────────────
# QA flag bits
# ──────────────────────────────────────────────────────────────
QA_NANDATA   = 1    # NaN or inf values found in the input flux
QA_MODELFAIL = 2    # ESO sky model could not be fetched for this row
QA_FAILED    = 4    # row raised an exception; spectra filled with NaN
QA_ZEROCONT  = 8    # model continuum is zero at every clean pixel
QA_ZEROLINE  = 16   # model LINES template is zero at every line pixel

_QA_FLAG_NAMES = {
    QA_NANDATA:   'NANDATA',
    QA_MODELFAIL: 'MODELFAIL',
    QA_FAILED:    'FAILED',
    QA_ZEROCONT:  'ZEROCONT',
    QA_ZEROLINE:  'ZEROLINE',
}

DEFAULT_MASK_FILE = Path(__file__).resolve().parent.parent / 'data' / 'sky_mask.fits'

_ALL_EXTS = ['FLUX', 'SKY_EAST', 'SKY_WEST']

# Which DRP_ALL columns carry ra/dec for each XCframe/XSFrame extension.
_EXT_RADEC_COLS = {
    'FLUX':     ('sci_ra',  'sci_dec'),
    'SKY_EAST': ('skye_ra', 'skye_dec'),
    'SKY_WEST': ('skyw_ra', 'skyw_dec'),
}

# Output DRP_ALL column-name prefix for each XCframe/XSFrame extension's
# results, matching the input file's own sci_*/skye_*/skyw_* convention.
_EXT_COL_PREFIX = {
    'FLUX':     'SCI',
    'SKY_EAST': 'SKYE',
    'SKY_WEST': 'SKYW',
}

_USAGE = '''Usage:
  SkyObsESOCompare.py [-ext EXT [-ext EXT ...]] [-engine local|remote|auto]
                      [-mask FILE] [-delta N] [-out ROOT] [-model] filename

Arguments:
  filename         XCframe/XSFrame summary file, or Sky_<name>.fits from
                   GetSky_from_CFrame_sum.py (auto-detected)

Options:
  -ext EXT         FLUX | SKY_EAST | SKY_WEST -- may repeat.
                   Default: all three. Ignored for Sky_<name>.fits input.
  -engine E        local (default) | remote | auto
  -mask FILE       clean-pixel mask FITS file (default: data/sky_mask.fits)
  -delta N         step size through rows for quick tests (default: 1)
  -out ROOT        output filename root (default: <stem>_esocompare)
  -model           also write the raw unscaled ESO model prediction
                   (<EXT>_MODEL) alongside each observed spectrum, which
                   is always written (off by default -- see Notes)
'''


# ──────────────────────────────────────────────────────────────
# Continuum and line fits
# ──────────────────────────────────────────────────────────────

def _fit_continuum_3param(flux, clean, moon, zodi, diffuse):
    '''
    Best-shape non-negative 3-component fit (a*MOON + b*ZODI + c*DIFFUSE)
    to flux, restricted to clean (continuum-only) pixels.

    Returns (cont3, (a, b, c), (a_err, b_err, c_err)).  Errors are the
    ordinary (unconstrained) least-squares formal errors at the clean
    pixels -- see module Notes on why NNLS's non-negativity bound is not
    reflected in them.
    '''
    A = np.column_stack([moon, zodi, diffuse])
    coef, _ = nnls(A[clean, :], flux[clean])
    cont3 = A @ coef

    n_clean = int(clean.sum())
    dof = max(n_clean - 3, 1)
    resid_clean = flux[clean] - cont3[clean]
    sigma2 = np.sum(resid_clean ** 2) / dof
    try:
        cov = sigma2 * np.linalg.inv(A[clean, :].T @ A[clean, :])
        errs = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        errs = np.full(3, np.nan)

    return cont3, tuple(float(v) for v in coef), tuple(float(v) for v in errs)


def _fit_continuum_1param(flux, clean, cont0):
    '''
    Single flux-cal-style non-negative scale k*cont0, restricted to clean
    pixels, where cont0 = MOON + ZODI + DIFFUSE (the model's own predicted
    mix, unscaled).

    Returns (k, k_err, cont1).
    '''
    denom = np.dot(cont0[clean], cont0[clean])
    if denom == 0:
        return 0.0, np.nan, np.zeros_like(flux)

    k = max(np.dot(flux[clean], cont0[clean]) / denom, 0.0)
    cont1 = k * cont0

    n_clean = int(clean.sum())
    dof = max(n_clean - 1, 1)
    resid_clean = flux[clean] - cont1[clean]
    sigma2 = np.sum(resid_clean ** 2) / dof
    k_err = float(np.sqrt(sigma2 / denom))

    return float(k), k_err, cont1


def _fit_line_scale(resid, line_mask, lines):
    '''
    Single non-negative scale d*lines fit to resid, restricted to
    line-affected pixels.  Returns (d, d_err).
    '''
    denom = np.dot(lines[line_mask], lines[line_mask])
    if denom == 0:
        return 0.0, np.nan

    d = max(np.dot(resid[line_mask], lines[line_mask]) / denom, 0.0)

    n_line = int(line_mask.sum())
    dof = max(n_line - 1, 1)
    resid_clean = resid[line_mask] - d * lines[line_mask]
    sigma2 = np.sum(resid_clean ** 2) / dof
    d_err = float(np.sqrt(sigma2 / denom))

    return float(d), d_err


def _overall_stats(wave, resid, mask):
    '''
    Single overall (not per-arm) NMAD and RMS of resid at mask pixels --
    reuses GetSkyCont.arm_continuum_stats (which already computes both in
    one pass) with one band spanning the full wavelength range, rather
    than duplicating its formulas.  Returns (nmad, rms).

    NMAD (robust, median-based) and RMS (what a least-squares fit
    actually minimizes) are NOT guaranteed to move together: a more
    flexible fit is guaranteed to have RMS <= a nested, more constrained
    fit's RMS, but can have EITHER a smaller or larger NMAD, since the
    extra flexibility often improves RMS mainly by reducing a handful of
    larger/outlier residuals while barely moving the median.  Keep both:
    RMS for a reliable "is fit A actually better than fit B" comparison,
    NMAD as a "how big is a typical pixel's residual" diagnostic.
    '''
    stats = arm_continuum_stats(wave, resid, clean=mask,
                                arm_ranges=[('ALL', float(wave.min()), float(wave.max()))])
    return stats['ALL']['nmad'], stats['ALL']['rms']


def _noise_proxy(flux, clean):
    '''
    Robust per-spectrum noise estimate from consecutive-pixel differences
    within clean pixels, independent of any continuum model::

        sigma = 1.4826 * median(|diff|) / sqrt(2)

    Used to normalize NMAD residuals into a reduced-chi-squared-like
    FIT_QUALITY ratio (see module Description) that's directly comparable
    across rows regardless of brightness -- unlike raw NMAD, which is in
    absolute flux units and so is naturally larger for brighter fields
    even when the model's shape match is equally good.
    '''
    vals = flux[clean]
    vals = vals[np.isfinite(vals)]
    if len(vals) < 3:
        return np.nan
    diffs = np.diff(vals)
    return float(1.4826 * np.median(np.abs(diffs)) / np.sqrt(2))


# Thresholds for the per-region line-residual "fail fraction" metrics
# below: one relative to NOISE_PROXY (chosen empirically -- 3 sigma
# saturates near 100% almost everywhere and has little discriminating
# power; 10 sigma has real row-to-row spread), and two fixed absolute
# flux thresholds spanning roughly the noise floor to the typical line
# brightness scale seen in practice.
_LINE_FRAC_SIGMA = 10.0
_LINE_FRAC_ABS_THRESH = (1e-15, 1e-14)


def _find_segments(mask):
    '''
    Return a list of (start, end) [end exclusive] index pairs, one per
    maximal contiguous run of True in mask -- e.g. each run of
    line-affected pixels belonging to one airglow feature (or a cluster
    of blended ones).  Vectorized (diff-based), so cheap enough to call
    fresh every row rather than threading a precomputed list through the
    whole do_all/one_row call chain.
    '''
    padded = np.concatenate(([0], mask.astype(int), [0]))
    idx = np.flatnonzero(np.diff(padded))
    return list(zip(idx[0::2].tolist(), idx[1::2].tolist()))


def _line_region_fail_fraction(resid, segments, threshold):
    '''
    Fraction of segments (as returned by _find_segments) whose median
    |resid| within that segment exceeds threshold -- i.e. the fraction of
    DISTINCT line-affected regions (not pixels) where the model's line
    fit is off by more than threshold, treating a bright, wide feature
    (many pixels) and a faint, narrow one (few pixels) as equally weighted
    votes.  Complements FIT_QUALITY_LINE (a single blended number, and
    itself dominated by whichever regions have the most/brightest pixels)
    with "how many separate features does the model actually get wrong."
    '''
    if not segments:
        return np.nan
    n_bad = sum(1 for s, e in segments if np.median(np.abs(resid[s:e])) > threshold)
    return float(n_bad) / len(segments)


# ──────────────────────────────────────────────────────────────
# Per-spectrum fit
# ──────────────────────────────────────────────────────────────

def one_row(wave, flux, ra, dec, obstime, clean, engine='local', verbose=False):
    '''
    Fit the ESO sky model's continuum and line templates to a single
    observed spectrum, over disjoint clean/line-affected pixel sets.

    Returns (result, qa_flags, coeffs, err_msg).

    result is a dict with one key, MODEL: the raw, unscaled ESO model
    prediction (MOON+ZODI+DIFFUSE+LINES, interpolated onto wave, no
    fitting applied) -- an array the same length as wave.

    coeffs is a flat dict with keys cont_scale, cont_scale_err,
    sep_cont_scale, sep_cont_nmad, sep_cont_rms, cont_nmad, cont_rms,
    noise_proxy, sep_fit_quality_cont, fit_quality_cont, sep_line_scale,
    sep_line_scale_err, sep_line_nmad, sep_line_rms, sep_fit_quality_line,
    sep_line_frac_10sig, sep_line_frac_1e15, sep_line_frac_1e14, line_scale,
    line_scale_err, line_nmad, line_rms, fit_quality_line, line_frac_10sig,
    line_frac_1e15, line_frac_1e14, moon_flux, zodi_flux, diffuse_flux,
    line_flux (see module Description for what each means).

    err_msg is '' on success, else a description of what failed.  On
    error, returns (None, QA_FAILED, {}, err_msg) or
    (None, QA_MODELFAIL, {}, err_msg).

    verbose : bool, default False -- passed straight through to
        SkySepESO._get_sky_model, which defaults to True itself; False
        here (unlike that default) since this script's whole purpose is
        processing many rows in a batch, where echoing each individual
        model fetch's routine progress would flood the terminal without
        adding anything -- the same diagnostic text is still returned in
        err_msg (and stored in DRP_ALL / <outroot>_errors.txt) on failure
        regardless of verbose.
    '''
    qa_flags = 0
    flux = np.asarray(flux, dtype=float)

    if not np.all(np.isfinite(flux)):
        qa_flags |= QA_NANDATA
        # see SkySepESO.one_drp for why this is zeroed explicitly rather
        # than left to nan_to_num's default (+-inf -> +-1.8e308)
        flux = np.nan_to_num(flux, nan=0.0, posinf=0.0, neginf=0.0)

    model_tab, err = _get_sky_model(ra, dec, obstime, engine=engine, verbose=verbose)
    if model_tab is None:
        return None, QA_MODELFAIL, {}, err

    moon    = np.interp(wave, model_tab['WAVE'], model_tab['MOON'])
    zodi    = np.interp(wave, model_tab['WAVE'], model_tab['ZODI'])
    diffuse = np.interp(wave, model_tab['WAVE'], model_tab['DIFFUSE'])
    lines   = np.interp(wave, model_tab['WAVE'], model_tab['LINES'])
    cont0   = moon + zodi + diffuse

    line_mask = ~clean

    if np.dot(cont0[clean], cont0[clean]) == 0:
        qa_flags |= QA_ZEROCONT
    if np.dot(lines[line_mask], lines[line_mask]) == 0:
        qa_flags |= QA_ZEROLINE

    # a/b/c (MOON/ZODI/DIFFUSE coefficients) are only used here to build
    # CONT3, the best available continuum estimate -- not reported
    # individually (see module History on trimming reported columns).
    cont3, _coef3, _err3 = _fit_continuum_3param(flux, clean, moon, zodi, diffuse)
    k, k_err, cont1 = _fit_continuum_1param(flux, clean, cont0)

    # SEP_CONT_SCALE: the 3-parameter (separation-stage) fit's flux ratio
    # to the model's raw, unscaled continuum -- compare directly to
    # CONT_SCALE (the 1-parameter/final fit's own scale) to see whether
    # the model's predicted MOON:ZODI:DIFFUSE mix needed reweighting to
    # get a good shape match.
    denom0 = np.sum(cont0[clean])
    sep_cont_scale = float(np.sum(cont3[clean]) / denom0) if denom0 != 0 else np.nan

    # Two full line fits, over the SAME line-affected pixels but against
    # different continuum baselines -- SEP uses CONT3 (the separation
    # stage's own continuum estimate), FINAL/unprefixed uses CONT1 (the
    # continuum CONT_SCALE actually describes).  This mirrors the
    # continuum comparison below: it directly measures how much the
    # separation-stage's extra flexibility helps the DOWNSTREAM line
    # answer too, not just the continuum answer.
    sep_line_input = flux - cont3
    sep_d, sep_d_err = _fit_line_scale(sep_line_input, line_mask, lines)
    sep_lines_scaled = sep_d * lines

    final_line_input = flux - cont1
    d, d_err = _fit_line_scale(final_line_input, line_mask, lines)
    lines_scaled = d * lines

    # SEP_CONT_NMAD/RMS: residual of the continuum SEPARATION stage
    # (currently the 3-parameter CONT3 fit; whatever replaces it later,
    # this is where its residual goes).  CONT_NMAD/RMS (unprefixed):
    # residual of the FINAL scaled continuum (CONT1, the thing CONT_SCALE
    # actually describes).  CONT3 is a strict superset of CONT1 (CONT1 =
    # CONT3 with a=b=c forced), so SEP_CONT_RMS <= CONT_RMS always (RMS is
    # what the least-squares fits actually minimize) -- but NMAD is a
    # different, robust statistic with no such guarantee (see
    # _overall_stats).  Same SEP-vs-final pairing for the two line fits.
    sep_cont_resid  = flux - cont3
    cont_resid      = flux - cont1
    sep_line_resid   = sep_line_input - sep_lines_scaled
    final_line_resid = final_line_input - lines_scaled

    sep_cont_nmad, sep_cont_rms = _overall_stats(wave, sep_cont_resid, clean)
    cont_nmad, cont_rms         = _overall_stats(wave, cont_resid, clean)
    sep_line_nmad, sep_line_rms = _overall_stats(wave, sep_line_resid, line_mask)
    line_nmad, line_rms         = _overall_stats(wave, final_line_resid, line_mask)

    # Reduced-chi-squared-like ratios: each NMAD above (absolute flux
    # units, naturally larger for brighter fields) divided by a per-row,
    # model-independent noise estimate, so FIT_QUALITY is comparable
    # across rows regardless of brightness.  ~1 = residual consistent
    # with noise (excellent shape match); >>1 = real shape mismatch.
    # Both LINE ratios necessarily share the continuum's noise proxy (no
    # noise-only region exists within the line-affected pixels
    # themselves to estimate it independently).
    noise_proxy = _noise_proxy(flux, clean)
    with np.errstate(divide='ignore', invalid='ignore'):
        sep_fit_quality_cont = float(sep_cont_nmad / noise_proxy)
        fit_quality_cont     = float(cont_nmad / noise_proxy)
        sep_fit_quality_line = float(sep_line_nmad / noise_proxy)
        fit_quality_line     = float(line_nmad / noise_proxy)

    # Per-region ("fail fraction") line-quality metrics: FIT_QUALITY_LINE
    # is one blended number, itself dominated by whichever line-affected
    # regions have the most/brightest pixels (see module Notes); these
    # instead ask what FRACTION of the distinct line-affected regions
    # (each an equally-weighted vote, regardless of brightness or width)
    # have a typical (median) residual worse than a threshold -- one
    # relative to noise (_LINE_FRAC_SIGMA), two fixed absolute flux
    # levels (_LINE_FRAC_ABS_THRESH).
    line_segments = _find_segments(line_mask)
    sig_thresh = _LINE_FRAC_SIGMA * noise_proxy
    sep_line_frac_10sig = _line_region_fail_fraction(sep_line_resid, line_segments, sig_thresh)
    line_frac_10sig     = _line_region_fail_fraction(final_line_resid, line_segments, sig_thresh)
    abs1, abs2 = _LINE_FRAC_ABS_THRESH
    sep_line_frac_1e15 = _line_region_fail_fraction(sep_line_resid, line_segments, abs1)
    line_frac_1e15     = _line_region_fail_fraction(final_line_resid, line_segments, abs1)
    sep_line_frac_1e14 = _line_region_fail_fraction(sep_line_resid, line_segments, abs2)
    line_frac_1e14     = _line_region_fail_fraction(final_line_resid, line_segments, abs2)

    # Mean brightness of each unscaled model component, each over the same
    # pixel footprint its corresponding scale factor is fit from (clean
    # pixels for the continuum components, line-affected pixels for
    # LINES) -- mean rather than median so that (SCALE-1)*<flux> is a
    # physically meaningful flux-error estimate (the scale fits themselves
    # are least-squares/dot-product, i.e. flux-weighted, so mean is the
    # matching summary statistic; median would be biased low for LINES in
    # particular, since most line-affected pixels are wings/gaps between
    # narrow line cores, not the cores themselves).
    moon_flux    = float(np.mean(moon[clean]))
    zodi_flux    = float(np.mean(zodi[clean]))
    diffuse_flux = float(np.mean(diffuse[clean]))
    line_flux    = float(np.mean(lines[line_mask]))

    # The raw, UNSCALED ESO model prediction (no fitting applied) -- the
    # only spectral array this script writes optionally (-model).  CONT1/
    # CONT3/LINES/RESID (the various fitted/scaled quantities) are not
    # returned at all any more: they're fully reconstructable from MODEL's
    # components (MOON+ZODI+DIFFUSE+LINES) and the scalar coefficients in
    # coeffs below, so keeping them as separate output arrays would only
    # duplicate what DRP_ALL + MODEL already contain.
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

def _write_qa_report(flagged, outroot):
    '''Print, and save to <outroot>_errors.txt, a report of every (ext,
    row, flags, msg) tuple in flagged whose QA_FLAGS is non-zero.'''
    if not flagged:
        return
    report_lines = ['Rows with QA flags:',
                    '  %-10s %-8s %-22s %s' % ('Ext', 'Row', 'Flags', 'Reason')]
    for ext_label, row, flags, msg in flagged:
        active = [name for bit, name in _QA_FLAG_NAMES.items() if flags & bit]
        report_lines.append('  %-10s %-8d %-22s %s'
                            % (ext_label, row, ','.join(active), msg or '(no message)'))
    report = '\n'.join(report_lines)
    print(report)

    errfile = '%s_errors.txt' % outroot
    with open(errfile, 'w') as fh:
        fh.write(report + '\n')
    print('Wrote error summary to %s' % errfile)


def _fit_one_extension(x, drp_all, wave, ext, rows, clean, engine, nan_spectrum,
                       save_model=False):
    '''
    Fit every row in rows for one XCframe/XSFrame extension against the
    ESO model.  Returns a dict with qa/coeffs/err (len(rows)-long lists;
    element j corresponds to rows[j]), the observed flux spectra actually
    analyzed (also a len(rows)-long list, always kept), and, only if
    save_model is True, the raw unscaled ESO model prediction per row
    (one_row computes it regardless of save_model since it's cheap; this
    only controls whether it's kept and written out afterward).
    '''
    ra_col, dec_col = _EXT_RADEC_COLS[ext]
    final_flux, final_model = [], []
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
            if save_model:
                final_model.append(nan_spectrum.copy())
        elif save_model:
            final_model.append(result['MODEL'])

        qa_flags_list.append(row_flags)
        coeffs_list.append(coeffs)
        err_msgs.append(err_msg)

        if (n + 1) % 100 == 0 or (n + 1) == len(rows):
            print('  Completed %6d of %d' % (n + 1, len(rows)))

    return dict(flux=final_flux, model=final_model,
               qa=qa_flags_list, coeffs=coeffs_list, err=err_msgs)


def do_all(filename, exts=None, engine='local', mask_file=DEFAULT_MASK_FILE,
          idelta=1, outroot='', save_model=False):
    '''
    Fit the ESO sky model's continuum and line templates to every
    requested spectrum in an input file and write the results.

    Parameters
    ----------
    filename : str
        XCframe/XSFrame summary file, or a Sky_<name>.fits file.
    exts : list of str or None
        Extensions to process (XCframe/XSFrame mode only); None means all
        three (FLUX, SKY_EAST, SKY_WEST).
    engine : str
        'local' (default), 'remote', or 'auto'.
    mask_file : str or Path
        Clean-pixel mask FITS file; default data/sky_mask.fits.
    idelta : int
        Row step size (1 = all rows).
    outroot : str
        Output filename root; defaults to <stem>_esocompare.
    save_model : bool
        The observed spectra actually analyzed (FLUX/SKY_EAST/SKY_WEST,
        whichever were requested) are always written. If save_model is
        also True, the raw unscaled ESO model prediction (MOON+ZODI+
        DIFFUSE+LINES, no fitting applied) is written alongside each one
        as a <EXT>_MODEL extension; off by default since it roughly
        doubles output size and DRP_ALL alone already carries every
        scalar this script computes.
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
        outroot = '%s_esocompare' % stem

    if sky_mode:
        # One row per input row (a physical exposure may contribute a
        # SKY_EAST row, a SKY_WEST row, both, or neither -- see
        # Description -- so there is no further per-exposure merge here).
        if exts:
            print("Warning: -ext is ignored for Sky_<name>.fits input "
                  "(spectrum type is already fixed per row via 'tel').")
        tel = (np.array(drp_all['tel'], dtype=str) if 'tel' in drp_all.colnames
              else np.full(len(drp_all), 'FLUX'))
        print('Processing %d spectra ...' % len(rows))

        final_flux, final_model = [], []
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
                if save_model:
                    final_model.append(nan_spectrum.copy())
            elif save_model:
                final_model.append(result['MODEL'])

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
        if save_model:
            spectral_hdus.append(fits.ImageHDU(data=np.array(final_model), name='MODEL'))
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
        # One row per exposure (per entry of rows), regardless of how many
        # extensions are requested: each extension is fit independently
        # over all of rows, then folded into the SAME xtab as a block of
        # <PREFIX>_<COLNAME> columns (SCI_/SKYE_/SKYW_) -- see Description.
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
            per_ext[e] = _fit_one_extension(
                x, drp_all, wave, e, rows, clean, engine, nan_spectrum,
                save_model=save_model)

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
            if save_model:
                spectral_hdus.append(fits.ImageHDU(data=np.array(d['model']), name='%s_MODEL' % e))
        ext_hdr_val = ','.join(use_exts)

        # One row per exposure: start from the original DRP_ALL rows (all
        # existing sci_/skye_/skyw_ metadata, moon geometry, etc. survive
        # unchanged), then add this script's results as <PREFIX>_<COLNAME>
        # columns per requested extension.
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
    hdr['TITLE']  = 'SkyObsESOCompare'
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

    exts     = []
    engine   = 'local'
    mask_file = DEFAULT_MASK_FILE
    idelta   = 1
    outroot  = ''
    save_model = False
    filename = None

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
        elif arg == '-model':
            save_model = True
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

    valid_engines = {'local', 'remote', 'auto'}
    if engine not in valid_engines:
        print('Error: -engine must be one of: %s' % ', '.join(sorted(valid_engines)))
        sys.exit(1)

    for e in exts:
        if e not in _ALL_EXTS:
            print('Error: -ext must be one of: %s' % ', '.join(_ALL_EXTS))
            sys.exit(1)

    if not os.path.exists(str(mask_file)):
        print('Error: mask file not found: %s' % mask_file)
        sys.exit(1)

    do_all(filename=filename, exts=(exts or None), engine=engine,
          mask_file=mask_file, idelta=idelta, outroot=outroot,
          save_model=save_model)
