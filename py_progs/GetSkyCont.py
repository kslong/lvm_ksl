#!/usr/bin/env python
# coding: utf-8
'''
                    Space Telescope Science Institute

Synopsis:

    Fit a smooth continuum to LVM sky spectra using a two-component B-spline
    design matrix, restricting the fit to wavelength regions unaffected by sky
    emission lines.  Operates independently of the lvmsky and lvmdrp packages.

Command line usage (if any):

    # Sky-file mode (auto-detected from file content):
    usage: GetSkyCont.py sky_file.fits -mask mask.fits [row_no ...] [-delta N]

    # XCframe / XSFrame summary file mode:
    usage: GetSkyCont.py xframe.fits ext -mask mask.fits [row_no ...] [-delta N]

    Arguments:

    sky_file.fits  Sky_<name>.fits from GetSky_from_CFrame_sum.py; all spectra processed if no row numbers are given.
    xframe.fits    an XCframe or XSFrame summary FITS file.
    ext            FITS extension for the spectrum (FLUX, SKY_EAST, SKY_WEST, ...).
    row_no         zero or more 0-based row indices; if omitted all rows are processed (subject to -delta).

    Options:

    -mask file     (required) palace_mask FITS file from palace_make_mask.py (MASK extension: 1=clean 0=line-affected).
    -delta N       process every N-th row (0, N, 2N, ...) instead of all; ignored when explicit row numbers are given.
    -kstep N       B-spline knot spacing in Angstroms (default 100).
    -out outroot   set the output filename root.

Description:

    Builds a two-component cubic B-spline design matrix spanning the full LVM
    wavelength range with knots every kstep Angstroms (default 100 A, giving
    ~66 basis functions per component):

    DIFFUSE component: plain B-spline basis functions.  These absorb smooth
    diffuse airglow emission (HO2, FeO, scattered light) whose spectral shape
    is not described by the solar spectrum.

    MOON component: the same B-spline basis functions each multiplied
    pointwise by a rebinned solar spectrum.  The Fraunhofer absorption
    structure is fixed; the B-spline envelope adjusts the colour and
    amplitude to fit the Moon-reflected and zodiacal-light continuum.

    The design matrix is evaluated only at pixels flagged as clean by the
    palace_make_mask.py mask (MASK=1).  Non-negative least squares (scipy
    nnls) enforces positive component amplitudes.  The fitted components are
    then evaluated at all wavelengths so both interpolate smoothly across
    masked line regions.

    The residual FLUX - CONT = FLUX - MOON - DIFFUSE isolates line emission
    free of continuum bias and is saved for subsequent line fitting.

Primary routines:

    build_design_matrix  construct the two-component B-spline design matrix.
    fit_continuum        fit and decompose continuum for one spectrum.
    process_skyfile      batch process a Sky_<name>.fits file.
    process_many         batch process an XCframe/XSFrame file.

Notes:

    The solar spectrum is read from
    /Users/long/Projects/lvm_sky2606/skysub_ivan/Spectre_HR_LATMOS_Meftah_V1_350_1000nm.txt
    (wavelengths in nm, converted to Angstroms internally).  If the file is
    absent the MOON component is omitted and CONT = DIFFUSE only; a warning
    is printed.

    Output FITS extensions: PRIMARY, WAVE, FLUX, CONT, MOON, DIFFUSE, RESID, MASK, DRP_ALL.
    Default output name: <stem>_cont.fits or <stem>_<ext>_cont.fits.

History:

    260628  ksl  Written; solar component optional.
    260628  ksl  Solar made default; MOON and DIFFUSE always written as separate extensions.
'''

import sys
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import BSpline, interp1d
from scipy.optimize import nnls


DEFAULT_SOLAR_FILE = Path('/Users/long/Projects/lvm_sky2606/skysub_ivan/'
                          'Spectre_HR_LATMOS_Meftah_V1_350_1000nm.txt')

_USAGE = '''Usage:
  GetSkyCont.py sky_file.fits -mask mask.fits [row_no ...] [-delta N]
  GetSkyCont.py xframe.fits ext -mask mask.fits [row_no ...] [-delta N]

Arguments:
  sky_file.fits  Sky_<name>.fits from GetSky_from_CFrame_sum.py
  xframe.fits    XCframe or XSFrame summary FITS file
  ext            FITS extension (FLUX, SKY_EAST, SKY_WEST, ...)
  row_no         0-based row indices (default: all rows)

Options:
  -mask file   (required) palace_mask FITS file from palace_make_mask.py
  -delta N     process every N-th row instead of all
  -kstep N     B-spline knot spacing in Angstroms (default 100)
  -out ROOT    output filename root
'''

_sky_fits_cache = {}
_xcframe_fits   = None


# ──────────────────────────────────────────────────────────────
# File-type detection
# ──────────────────────────────────────────────────────────────

def _is_sky_file(filename):
    fn = str(filename)
    if fn not in _sky_fits_cache:
        _sky_fits_cache[fn] = fits.open(filename)
    try:
        cols = [c.name.lower() for c in _sky_fits_cache[fn]['DRP_ALL'].columns]
        return 'tel' in cols
    except Exception:
        return False


# ──────────────────────────────────────────────────────────────
# Mask and solar loading
# ──────────────────────────────────────────────────────────────

def load_mask(mask_file):
    '''Load a palace_make_mask.py output and return (wave, mask).

    Parameters
    ----------
    mask_file : str or Path
        FITS file with WAVE and MASK extensions.

    Returns
    -------
    wave : ndarray
        Wavelength array (Angstroms).
    mask : ndarray of bool
        True where the pixel is clean (unaffected by sky lines).
    '''
    hdul = fits.open(mask_file)
    wave = hdul['WAVE'].data.astype(float)
    mask = hdul['MASK'].data.astype(bool)
    hdul.close()
    return wave, mask


def _interp_mask_to_wave(mask_wave, mask_bool, spec_wave):
    '''Resample a boolean mask onto a different wavelength grid (nearest-neighbour).'''
    f = interp1d(mask_wave, mask_bool.astype(float), kind='nearest',
                 bounds_error=False, fill_value=0.0)
    return f(spec_wave) > 0.5


def load_solar(wave_out):
    '''Load and resample the solar spectrum onto wave_out (Angstroms).

    Parameters
    ----------
    wave_out : ndarray
        Target wavelength grid in Angstroms.

    Returns
    -------
    solar : ndarray or None
        Solar flux resampled to wave_out, normalised to unit median.
        Returns None if the solar file is not found.
    '''
    if not DEFAULT_SOLAR_FILE.exists():
        print(f'Warning: solar file not found at {DEFAULT_SOLAR_FILE}; MOON component omitted.')
        return None
    data     = np.loadtxt(DEFAULT_SOLAR_FILE, comments=';')
    wave_sol = data[:, 0] * 10.0   # nm → Angstroms
    flux_sol = data[:, 1]
    f        = interp1d(wave_sol, flux_sol, kind='linear',
                        bounds_error=False, fill_value=0.0)
    solar    = f(wave_out)
    med      = np.nanmedian(solar[solar > 0]) if np.any(solar > 0) else 1.0
    return solar / med


# ──────────────────────────────────────────────────────────────
# Design matrix and fitting
# ──────────────────────────────────────────────────────────────

def build_design_matrix(wave, knot_step=100.0, solar=None):
    '''Build a two-component cubic B-spline design matrix.

    Parameters
    ----------
    wave : ndarray
        Wavelength array at which to evaluate the basis functions.
    knot_step : float
        Knot spacing in Angstroms.
    solar : ndarray or None
        Solar spectrum (same length as wave).  If given, a second block of
        columns equal to solar * B-spline_basis is appended for the MOON
        component.  If None only the DIFFUSE block is built.

    Returns
    -------
    A : ndarray, shape (len(wave), n_cols)
        Full design matrix.  n_cols = n_basis (DIFFUSE only) or 2*n_basis (with MOON).
    n_basis : int
        Number of columns in each component block.
    '''
    wmin, wmax = wave[0], wave[-1]
    k     = 3
    t_int = np.arange(wmin + knot_step, wmax, knot_step)
    t     = np.r_[[wmin] * (k + 1), t_int, [wmax] * (k + 1)]
    n_b   = len(t) - k - 1

    A_plain = np.zeros((len(wave), n_b), dtype=float)
    for i in range(n_b):
        c = np.zeros(n_b)
        c[i] = 1.0
        A_plain[:, i] = np.nan_to_num(BSpline(t, c, k, extrapolate=False)(wave))

    if solar is not None:
        A_moon = A_plain * solar[:, np.newaxis]
        return np.hstack([A_plain, A_moon]), n_b
    return A_plain, n_b


def fit_continuum(flux, clean, A, n_b):
    '''Fit and decompose continuum for a single sky spectrum.

    Parameters
    ----------
    flux : ndarray
        Observed sky spectrum (N_pix,).
    clean : ndarray of bool
        True at pixels included in the fit.
    A : ndarray, shape (N_pix, n_cols)
        Pre-built design matrix from build_design_matrix.
    n_b : int
        Number of plain B-spline columns (first block = DIFFUSE).

    Returns
    -------
    cont : ndarray
        Total fitted continuum (DIFFUSE + MOON) at all wavelengths.
    diffuse : ndarray
        Smooth diffuse component (plain B-splines).
    moon : ndarray
        Moon/zodi component (solar x B-splines).  Zero array when A has only one block.
    coef : ndarray
        Non-negative NNLS coefficients.
    n_clean : int
        Number of clean pixels used in the fit.
    '''
    coef, _  = nnls(A[clean, :], flux[clean])
    diffuse  = A[:, :n_b] @ coef[:n_b]
    moon     = (A[:, n_b:] @ coef[n_b:]) if A.shape[1] > n_b else np.zeros_like(diffuse)
    cont     = diffuse + moon
    return cont, diffuse, moon, coef, int(clean.sum())


# ──────────────────────────────────────────────────────────────
# Spectrum I/O
# ──────────────────────────────────────────────────────────────

def _read_sky_row(filename, row):
    fn = str(filename)
    if fn not in _sky_fits_cache:
        _sky_fits_cache[fn] = fits.open(filename)
    hdul = _sky_fits_cache[fn]
    return (hdul['WAVE'].data.astype(float),
            hdul['FLUX'].data[row].astype(float))


def _read_xcframe_row(filename, row, ext):
    global _xcframe_fits
    if _xcframe_fits is None:
        _xcframe_fits = fits.open(filename)
    return (_xcframe_fits['WAVE'].data.astype(float),
            _xcframe_fits[ext].data[row].astype(float))


# ──────────────────────────────────────────────────────────────
# Batch processing
# ──────────────────────────────────────────────────────────────

def _rows_to_process(n_total, rows, delta):
    if rows is not None:
        return list(rows)
    if delta is not None:
        return list(range(0, n_total, int(delta)))
    return list(range(n_total))


def _write_output(outfile, wave, flux_arr, cont_arr, moon_arr, diffuse_arr,
                  resid_arr, clean_pixels, qtab, hdr):
    hdu_list = [
        fits.PrimaryHDU(header=hdr),
        fits.ImageHDU(data=wave.astype(np.float32),          name='WAVE'),
        fits.ImageHDU(data=flux_arr.astype(np.float32),      name='FLUX'),
        fits.ImageHDU(data=cont_arr.astype(np.float32),      name='CONT'),
        fits.ImageHDU(data=moon_arr.astype(np.float32),      name='MOON'),
        fits.ImageHDU(data=diffuse_arr.astype(np.float32),   name='DIFFUSE'),
        fits.ImageHDU(data=resid_arr.astype(np.float32),     name='RESID'),
        fits.ImageHDU(data=clean_pixels.astype(np.uint8),    name='MASK'),
        fits.BinTableHDU(qtab,                                name='DRP_ALL'),
    ]
    fits.HDUList(hdu_list).writeto(outfile, overwrite=True)
    print(f'Wrote {outfile}')


def _get_wave(filename, ext=None):
    '''Return the wavelength array from an open file cache.'''
    if ext is not None:
        global _xcframe_fits
        if _xcframe_fits is None:
            _xcframe_fits = fits.open(filename)
        return _xcframe_fits['WAVE'].data.astype(float)
    fn = str(filename)
    if fn not in _sky_fits_cache:
        _sky_fits_cache[fn] = fits.open(filename)
    return _sky_fits_cache[fn]['WAVE'].data.astype(float)


def process_skyfile(filename, rows=None, delta=None, mask_wave=None, mask_arr=None,
                    knot_step=100.0, outroot=''):
    '''Fit continuum to spectra in a Sky_<name>.fits file.

    Parameters
    ----------
    filename : str or Path
        Sky file from GetSky_from_CFrame_sum.py.
    rows : list of int or None
        Explicit row indices; None means all (or every delta-th).
    delta : int or None
        Row stride; ignored if rows is not None.
    mask_wave : ndarray or None
    mask_arr : ndarray of bool or None
        Boolean mask (True=clean); None means all pixels.
    knot_step : float
    outroot : str
    '''
    fn = str(filename)
    if fn not in _sky_fits_cache:
        _sky_fits_cache[fn] = fits.open(filename)
    hdul    = _sky_fits_cache[fn]
    wave    = hdul['WAVE'].data.astype(float)
    n_total = hdul['FLUX'].data.shape[0]
    xtab    = Table(hdul['DRP_ALL'].data)

    rlist = _rows_to_process(n_total, rows, delta)

    clean = (_interp_mask_to_wave(mask_wave, mask_arr, wave)
             if mask_wave is not None else np.ones(len(wave), dtype=bool))
    print(f'Clean pixels: {clean.sum()} / {len(wave)}')

    solar   = load_solar(wave)
    A, n_b  = build_design_matrix(wave, knot_step=knot_step, solar=solar)
    print(f'Design matrix: {A.shape[1]} columns ({n_b} diffuse + {A.shape[1]-n_b} moon)')

    flux_list, cont_list, moon_list, diff_list, resid_list, n_clean_list = \
        [], [], [], [], [], []

    for i, row in enumerate(rlist):
        _, flux = _read_sky_row(filename, row)
        cont, diffuse, moon, _, n_cl = fit_continuum(flux, clean, A, n_b)
        flux_list.append(flux)
        cont_list.append(cont)
        moon_list.append(moon)
        diff_list.append(diffuse)
        resid_list.append(flux - cont)
        n_clean_list.append(n_cl)
        if (i + 1) % 10 == 0 or (i + 1) == len(rlist):
            print(f'Processed {i + 1}/{len(rlist)} (row {row})')

    qtab            = xtab[rlist].copy()
    qtab['row']     = rlist
    qtab['n_clean'] = n_clean_list

    if outroot == '':
        outroot = f'{Path(filename).stem}_cont'

    hdr = fits.Header()
    hdr['INPUT']  = str(filename)
    hdr['N_PROC'] = len(rlist)
    hdr['KSTEP']  = knot_step
    hdr['SOLAR']  = solar is not None
    if delta is not None:
        hdr['DELTA'] = delta

    _write_output(f'{outroot}.fits', wave,
                  np.array(flux_list), np.array(cont_list),
                  np.array(moon_list), np.array(diff_list),
                  np.array(resid_list), clean, qtab, hdr)


def process_many(filename, ext, rows=None, delta=None, mask_wave=None, mask_arr=None,
                 knot_step=100.0, outroot=''):
    '''Fit continuum to rows of an XCframe/XSFrame file.

    Parameters
    ----------
    filename : str or Path
    ext : str
        FITS extension name (e.g. SKY_EAST).
    rows : list of int or None
    delta : int or None
    mask_wave : ndarray or None
    mask_arr : ndarray of bool or None
    knot_step : float
    outroot : str
    '''
    global _xcframe_fits
    if _xcframe_fits is None:
        _xcframe_fits = fits.open(filename)
    wave    = _xcframe_fits['WAVE'].data.astype(float)
    n_total = _xcframe_fits[ext].data.shape[0]
    xtab    = Table(_xcframe_fits['DRP_ALL'].data)

    rlist = _rows_to_process(n_total, rows, delta)

    clean = (_interp_mask_to_wave(mask_wave, mask_arr, wave)
             if mask_wave is not None else np.ones(len(wave), dtype=bool))
    print(f'Clean pixels: {clean.sum()} / {len(wave)}')

    solar   = load_solar(wave)
    A, n_b  = build_design_matrix(wave, knot_step=knot_step, solar=solar)
    print(f'Design matrix: {A.shape[1]} columns ({n_b} diffuse + {A.shape[1]-n_b} moon)')

    flux_list, cont_list, moon_list, diff_list, resid_list, n_clean_list = \
        [], [], [], [], [], []

    for i, row in enumerate(rlist):
        _, flux = _read_xcframe_row(filename, row, ext)
        cont, diffuse, moon, _, n_cl = fit_continuum(flux, clean, A, n_b)
        flux_list.append(flux)
        cont_list.append(cont)
        moon_list.append(moon)
        diff_list.append(diffuse)
        resid_list.append(flux - cont)
        n_clean_list.append(n_cl)
        if (i + 1) % 10 == 0 or (i + 1) == len(rlist):
            print(f'Processed {i + 1}/{len(rlist)} (row {row})')

    qtab            = xtab[rlist].copy()
    qtab['row']     = rlist
    qtab['n_clean'] = n_clean_list

    if outroot == '':
        outroot = f'{Path(filename).stem}_{ext}_cont'

    hdr = fits.Header()
    hdr['INPUT']  = str(filename)
    hdr['EXT']    = ext
    hdr['N_PROC'] = len(rlist)
    hdr['KSTEP']  = knot_step
    hdr['SOLAR']  = solar is not None
    if delta is not None:
        hdr['DELTA'] = delta

    _write_output(f'{outroot}.fits', wave,
                  np.array(flux_list), np.array(cont_list),
                  np.array(moon_list), np.array(diff_list),
                  np.array(resid_list), clean, qtab, hdr)


# ──────────────────────────────────────────────────────────────
# Command-line interface
# ──────────────────────────────────────────────────────────────

def steer(argv):
    filename        = ''
    ext             = ''
    mask_file       = ''
    outroot         = ''
    delta           = None
    knot_step       = 100.0
    positional_ints = []

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(_USAGE)
            return
        elif argv[i] == '-mask':
            i += 1
            mask_file = argv[i]
        elif argv[i] == '-out':
            i += 1
            outroot = argv[i]
        elif argv[i] == '-delta':
            i += 1
            delta = int(argv[i])
        elif argv[i] == '-kstep':
            i += 1
            knot_step = float(argv[i])
        elif argv[i][0] == '-':
            print('Error: cannot parse command line:', argv)
            return
        elif filename == '':
            filename = argv[i]
        elif ext == '' and not argv[i].lstrip('-').isdigit():
            ext = argv[i]
        elif argv[i].lstrip('-').isdigit():
            positional_ints.append(int(argv[i]))
        else:
            print('Error: cannot parse command line:', argv)
            return
        i += 1

    if not filename:
        print(_USAGE)
        return

    if not mask_file:
        print('Error: -mask mask.fits is required.')
        print(_USAGE)
        return

    mask_wave, mask_arr = load_mask(mask_file)
    print(f'Loaded mask: {mask_arr.sum()} / {len(mask_arr)} clean pixels')

    rows = positional_ints if positional_ints else None

    sky_mode    = (ext == '') and _is_sky_file(filename)
    xcframe_mode = (ext != '')

    if not sky_mode and not xcframe_mode:
        print('Error: specify a FITS extension name (e.g. SKY_EAST) for XCframe files.')
        print(_USAGE)
        return

    if sky_mode:
        process_skyfile(filename, rows=rows, delta=delta,
                        mask_wave=mask_wave, mask_arr=mask_arr,
                        knot_step=knot_step, outroot=outroot)
    else:
        process_many(filename, ext, rows=rows, delta=delta,
                     mask_wave=mask_wave, mask_arr=mask_arr,
                     knot_step=knot_step, outroot=outroot)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
