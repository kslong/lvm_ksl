#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Analyze wavelength offsets in LVM CFrame/SFrame files by fiber,
using FFT cross-correlation of absorption features.

Command line usage (if any)::

    fourier_offset.py [-h] [-wmin 3900] [-wmax 4000] [-file xfile]
                      [-expmin N] [-expmax N] [-every N] [-ver drp_ver] [-plot]
                      filename [filename ...]

    -file xfile      reads the files to process from a Table containing a
                     column named filename or Filename.
    filename         one or more lvmCFrame or lvmSFrame FITS files.
    -wmin, -wmax     wavelength window for the cross-correlation
                     (default 3900-4000 AA).
    -expmin, -expmax inclusive bounds on the exposure number parsed from
                     the filename; exposures outside this range are skipped.
    -every N         process only every Nth exposure after range filtering
                     (default 1 = all exposures).
    -ver drp_ver     DRP version used to locate files via drpall
                     (default 1.2.1).
    -plot            force diagnostic plots even when more than 100 exposures
                     are in the list.  Plots are always made for <= 100 files.

Description:

    Reads the FLUX and WAVE extensions of each file, selects science fibers
    from the SLITMAP, then cross-correlates each spectrum against a median
    reference within [wmin, wmax] to estimate per-fiber wavelength offsets.

    Results for all processed exposures are collected into a FITS table named
    Fourier_<wmin>_<wmax>.<YYMMDD>.fits.  If that file already exists the new
    rows are merged in (duplicate exposures replaced) and the table is
    re-sorted by exposure number.

    Each row of the output table contains::

        Exposure             exposure number
        mjd                  modified Julian date of the observation
        dw_5pct              5th-percentile wavelength offset, all fibers (AA)
        dw_95pct             95th-percentile wavelength offset (AA)
        dw_med               median wavelength offset, all fibers (AA)
        dw_med_sp1/2/3       median wavelength offset per spectrograph (AA)
        dw_mad               scaled MAD of wavelength offsets, all fibers (AA)
        dw_mad_sp1/2/3       scaled MAD per spectrograph (AA)
        med_quality          median cross-correlation quality [0,1]
        flux_max_med         median of per-fiber peak flux in [wmin,wmax]
        flux_max_med_sp1/2/3 median peak flux per spectrograph
        flux_max_mad         scaled MAD of per-fiber peak flux, all fibers
        flux_max_mad_sp1/2/3 scaled MAD of peak flux per spectrograph
        dw_r01_04_sp1/2/3 .. dw_r25_sp1/2/3
                             median offset per ring group per spectrograph

    Diagnostic PNG plots (one per exposure) are written to Fig_Qual/.

Primary routines:

    doit            analyze a single file, return a dict of summary statistics
    do_all          loop over a list of files and write the summary FITS table
    filter_files    filter a file list by exposure number range and stride
    find_top        return the machine-specific top-level data directory
    resolve_filename  resolve a bare location path to a full CFrame path

Notes:

    The cross-correlation uses FFTs with parabolic sub-pixel refinement.
    A quality value in [0,1] is returned alongside each offset; values
    below ~0.5 indicate unreliable measurements (featureless window,
    bad pixels, low S/N).

    The scaled MAD (1.4826 * median(abs(x - median(x)))) equals the standard
    deviation for Gaussian data but is robust against outliers.

History:

260417 ksl Coding begun; wavelength_offset.py incorporated into this file
260501 ksl doit now returns a dict; do_all added; output switched to FITS; flux_max and dw MAD columns added; exposure range/stride filtering added

'''

import sys
import warnings
import datetime
import re
from astropy.io import ascii, fits
from astropy.table import Table, join
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os


# ── wavelength offset estimation (formerly wavelength_offset.py) ──────────────

def estimate_wavelength_offsets(
    w, flux, wmin, wmax,
    ref="median",
    normalize=True,
    continuum_deg=1,
    max_shift=None,
):
    """
    Estimate per-spectrum wavelength offsets from absorption features.

    Cross-correlates each spectrum with a reference template inside the
    window [wmin, wmax].  Sub-pixel precision is obtained by parabolic
    interpolation of the cross-correlation peak.

    Parameters
    ----------
    w : (N_wave,) array_like
        Wavelength array (approximately uniform spacing).
    flux : (N_spectra, N_wave) array_like
        Flux array.
    wmin, wmax : float
        Wavelength range that brackets the absorption features.
    ref : {'median', 'mean'} or (N_wave,) array_like
        Reference template.
        'median' (default) - median of all spectra in the window.
        'mean'             - mean of all spectra.
        array              - external reference with the same length as w.
    normalize : bool
        If True (default), divide each spectrum by a polynomial continuum
        fit before cross-correlating.  Removes the effect of varying
        mean flux levels.
    continuum_deg : int
        Degree of the polynomial continuum (default 1 = linear tilt).
    max_shift : float or None
        Maximum allowed shift in wavelength units (same as w).
        None (default) = no restriction.

    Returns
    -------
    dw : (N_spectra,) ndarray
        Wavelength offset per spectrum, in the same units as w.
        Positive -> features are shifted to longer wavelengths
        (spectrum is red-shifted relative to the reference).
        Negative -> features are shifted to shorter wavelengths.
    quality : (N_spectra,) ndarray
        Normalised cross-correlation peak in [0, 1].
        Values below ~0.5 indicate unreliable measurements
        (featureless window, bad pixels, low S/N, etc.).
    """
    w    = np.asarray(w,    dtype=float)
    flux = np.asarray(flux, dtype=float)

    mask  = (w >= wmin) & (w <= wmax)
    if mask.sum() < 5:
        raise ValueError("Fewer than 5 pixels in [wmin, wmax] – widen the range.")
    w_win  = w[mask]
    f_win  = np.where(np.isnan(flux[:, mask]), 0.0, flux[:, mask])
    N_pix  = w_win.size
    dw_pix = float(np.median(np.diff(w_win)))

    if isinstance(ref, str):
        if ref == "median":
            template = np.nanmedian(f_win, axis=0)
        elif ref == "mean":
            template = np.nanmean(f_win, axis=0)
        else:
            raise ValueError(f"ref must be 'median', 'mean', or an array; got {ref!r}")
    else:
        ref_arr  = np.asarray(ref, dtype=float)
        template = ref_arr[mask] if ref_arr.size == w.size else ref_arr
        if template.size != N_pix:
            raise ValueError("External ref array has wrong length.")

    def _normalise(arr2d):
        P = arr2d.shape[1]
        x = np.linspace(-1.0, 1.0, P)
        V = np.vander(x, continuum_deg + 1, increasing=True)
        coeffs, *_ = np.linalg.lstsq(V, arr2d.T, rcond=None)
        cont = (V @ coeffs).T
        safe = np.where(cont > 0, cont, np.finfo(float).tiny)
        return arr2d / safe

    if normalize:
        f_norm = _normalise(f_win)
        t_norm = _normalise(template[None, :])[0]
    else:
        f_norm = f_win.copy()
        t_norm = template.copy()

    N_fft   = N_pix
    T_fft   = np.fft.rfft(t_norm, n=N_fft)
    F_fft   = np.fft.rfft(f_norm, n=N_fft, axis=1)
    xcorr   = np.fft.irfft(F_fft * np.conj(T_fft), n=N_fft, axis=1)
    xcorr_s = np.fft.fftshift(xcorr, axes=1)

    centre = N_fft // 2
    lags   = np.arange(N_fft, dtype=float) - centre

    if max_shift is not None:
        max_p = int(np.ceil(abs(max_shift) / dw_pix))
        lo = max(0, centre - max_p)
        hi = min(N_fft, centre + max_p + 1)
    else:
        lo, hi = 0, N_fft

    peak_in_win  = np.argmax(xcorr_s[:, lo:hi], axis=1)
    peak_in_full = peak_in_win + lo

    idx    = np.arange(flux.shape[0])
    pk     = peak_in_full
    valid  = (pk > 0) & (pk < N_fft - 1)
    pk_lo  = np.where(valid, pk - 1, pk)
    pk_hi  = np.where(valid, pk + 1, pk)

    y0 = xcorr_s[idx, pk_lo]
    y1 = xcorr_s[idx, pk]
    y2 = xcorr_s[idx, pk_hi]

    denom      = 2.0 * (2.0 * y1 - y0 - y2)
    safe_denom = np.where(np.abs(denom) > 1e-30, denom, 1.0)
    sub_pix    = np.where(valid & (np.abs(denom) > 1e-30), (y2 - y0) / safe_denom, 0.0)
    sub_pix = np.clip(sub_pix, -0.5, 0.5)

    shift_pix = lags[pk] + sub_pix
    dw        = -shift_pix * dw_pix

    # Scale before squaring to avoid float64 overflow; ratio is scale-invariant
    f_scale = np.abs(f_norm).max(axis=1, keepdims=True)
    f_scale = np.where(f_scale > 0, f_scale, 1.0)
    t_scale = float(np.abs(t_norm).max()) or 1.0
    ac_t    = float(np.sum((t_norm / t_scale) ** 2))
    ac_f    = np.sum((f_norm / f_scale) ** 2, axis=1)
    norm    = np.sqrt(ac_f * ac_t)
    peak_scaled = xcorr_s[idx, pk] / (f_scale[:, 0] * t_scale)
    quality = np.where(norm > 0, peak_scaled / norm, 0.0)

    return dw, quality


# ── time conversion (copied from GetSolar.py) ────────────────────────────────

def _detect_time_format(time_input):
    if isinstance(time_input, datetime.datetime):
        return 'datetime'
    time_str = str(time_input) if isinstance(time_input, (int, float)) else time_input
    if re.search(r'\d{4}-\d{2}-\d{2}', time_str) or re.search(r'\d{4}/\d{2}/\d{2}', time_str):
        return 'iso'
    try:
        value = float(time_str)
        return 'jd' if value > 2400000 else 'mjd'
    except ValueError:
        if re.search(r'\d+[/:\-]\d+', time_str):
            return 'iso'
    return 'iso'

def _convert_to_time_object(time_input, input_format):
    try:
        if input_format == 'datetime':
            return Time(time_input, scale='utc')
        elif input_format == 'iso':
            return Time(time_input, format='isot', scale='utc')
        elif input_format == 'mjd':
            return Time(float(time_input), format='mjd', scale='utc')
        elif input_format == 'jd':
            return Time(float(time_input), format='jd', scale='utc')
    except Exception as e:
        raise ValueError('Failed to convert %r as %s format: %s' % (time_input, input_format, e))

def _convert_to_output_format(time_obj, output_format):
    try:
        if output_format == 'datetime':
            return time_obj.to_datetime()
        elif output_format == 'iso':
            return time_obj.iso
        elif output_format == 'mjd':
            return time_obj.mjd
        elif output_format == 'jd':
            return time_obj.jd
        else:
            raise ValueError('Unsupported output format: %s' % output_format)
    except Exception as e:
        raise ValueError('Failed to convert to %s format: %s' % (output_format, e))

def convert_time(time_input, output_format='datetime'):
    '''Convert between time formats (ISO string, MJD, JD, datetime).'''
    input_format = _detect_time_format(time_input)
    t = _convert_to_time_object(time_input, input_format)
    return _convert_to_output_format(t, output_format)


# ── file location ──────────────────────────────────────────────────────────────

def build_drp_lookup(drp_ver='1.2.1'):
    '''
    Build a dict mapping exposure number to full CFrame file path.

    Reads drpall via SummarizeCframe.read_drpall(), prepends the
    machine-specific top directory from SummarizeCframe.find_top(),
    and substitutes SFrame with CFrame in all locations.
    Returns an empty dict if drpall cannot be read.
    '''
    try:
        from SummarizeCframe import find_top, read_drpall
    except ImportError:
        print('Warning: could not import SummarizeCframe; file lookup unavailable.')
        return {}
    topdir = find_top()
    if not topdir:
        return {}
    drp_tab = read_drpall(drp_ver)
    if not drp_tab or len(drp_tab) == 0:
        return {}
    lookup = {}
    for row in drp_tab:
        loc = str(row['location']).replace('SFrame', 'CFrame')
        lookup[int(row['expnum'])] = os.path.join(topdir, loc)
    return lookup


def resolve_filename(filename, lookup=None):
    '''
    Return the full path to a CFrame file, or None if not found.

    First substitutes SFrame with CFrame.  If the file exists at that
    path it is returned immediately.  Otherwise the exposure number is
    extracted from the filename and used to look up the full path in
    *lookup* (a dict built by build_drp_lookup()).
    '''
    filename = filename.replace('SFrame', 'CFrame')
    if os.path.isfile(filename):
        return filename
    if lookup:
        expnum = expnum_from_filename(filename)
        if expnum > 0 and expnum in lookup:
            full = lookup[expnum]
            if os.path.isfile(full):
                return full
    return None


# ── main routines ──────────────────────────────────────────────────────────────

# Rings 1-4 combined (too few fibers individually); rings 5-25 separate
RING_GROUPS = [(1, 4)] + [(r, r) for r in range(5, 26)]


def ring_col_names():
    '''Return the ordered list of per-ring-group per-spectrograph column names.'''
    names = []
    for rmin, rmax in RING_GROUPS:
        rlabel = 'r%02d_%02d' % (rmin, rmax) if rmin != rmax else 'r%02d' % rmin
        for sp in (1, 2, 3):
            names.append('dw_%s_sp%d' % (rlabel, sp))
    return names


def ring_sp_medians(dw, slit_tab):
    '''
    Return median dw for each ring group / spectrograph combination.
    Order matches ring_col_names().  Returns nan where no fibers exist.
    '''
    ring = np.array(slit_tab['ringnum'])
    spec = np.array(slit_tab['spectrographid'])
    vals = []
    for rmin, rmax in RING_GROUPS:
        ring_mask = (ring >= rmin) & (ring <= rmax)
        for sp in (1, 2, 3):
            mask = ring_mask & (spec == sp)
            vals.append(float(np.median(dw[mask])) if mask.sum() > 0 else np.nan)
    return vals

def clean_slitmap(slit_tab, telescope='Sci'):
    '''
    Clean slit_tab and add a column for the row number.
    '''
    nstart = len(slit_tab)
    slit_tab['Row'] = np.arange(nstart)
    colnames = slit_tab.colnames
    for one in colnames:
        if one == 'telescope':
            slit_tab = slit_tab[slit_tab['telescope'] == telescope]
        if one == 'fibstatus':
            slit_tab = slit_tab[slit_tab['fibstatus'] == 0]
    nstop = len(slit_tab)
    return slit_tab


def get_band_flux(xwave, xflux, wmin=4000, wmax=4500):
    '''
    Return the per-fiber median flux in the wavelength band [wmin, wmax].

    Parameters
    ----------
    xwave : (N_wave,) array
        Wavelength array.
    xflux : (N_fibers, N_wave) array
        Flux array.
    wmin, wmax : float
        Wavelength band limits in the same units as xwave.

    Returns
    -------
    band : (N_fibers,) ndarray
        Median flux of each fiber within the band.
    '''
    id_min = np.searchsorted(xwave, wmin, side='right')
    id_max = np.searchsorted(xwave, wmax, side='right')
    xflux = xflux[:, id_min:id_max]
    band = np.nanmedian(xflux, axis=1)
    return band


def plot_offsets(dw, quality, xwave, xflux, xdrp, wmin=3900, wmax=4000, outroot='', expnum=-1):
    '''
    Plot the wavelength-offset diagnostic figure for a single exposure.

    Panels: (1) flux percentiles vs wavelength, (2) offset histogram,
    (3) spatial map of 90th-percentile flux in [wmin, wmax], (4) spatial offset map.
    Figure is saved to Fig_Qual/.
    '''
    x5  = np.percentile(dw, 5)
    x95 = np.percentile(dw, 95)

    plt.figure(1, (16, 4))
    plt.clf()

    plt.subplot(1, 4, 1)
    win_mask = (xwave >= wmin) & (xwave <= wmax)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        p5  = np.nanpercentile(xflux, 5,  axis=0)
        p50 = np.nanmedian(xflux,         axis=0)
        p95 = np.nanpercentile(xflux, 95, axis=0)
        ymin = np.nanmin(p5[win_mask])
        ymax = np.nanmax(p95[win_mask])
    plt.semilogy(xwave, p5,  label='5%')
    plt.semilogy(xwave, p50, label='50%')
    plt.semilogy(xwave, p95, label='95%')
    plt.legend()
    plt.xlim(wmin, wmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(r'$\lambda$ ($\AA$)')
    plt.ylabel('Flux')
    plt.title('Spectra by Brightness')

    plt.subplot(1, 4, 2)
    delta = 0.2
    plt.hist(dw, 100, range=(-delta, delta))
    plt.xlabel(r'$\delta\lambda$ ($\AA$)')
    plt.ylabel('N')
    plt.title('Wavelength Offsets')
    plt.text(0.05, 0.95, '5 - 95 %% offsets: %.3f %.3f' % (x5, x95),
             transform=plt.gca().transAxes)

    # Panel 3 — spatial map of per-fiber maximum flux in the wavelength window
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        fmax_flux = np.nanmax(xflux[:, win_mask], axis=1)
    fmax_lo = np.nanpercentile(fmax_flux, 5)
    fmax_hi = np.nanpercentile(fmax_flux, 95)
    ax3 = plt.subplot(1, 4, 3)
    sc3 = ax3.scatter(xdrp['ra'], xdrp['dec'], c=fmax_flux, cmap='viridis',
                      vmin=fmax_lo, vmax=fmax_hi)
    divider3 = make_axes_locatable(ax3)
    cax3 = divider3.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(sc3, cax=cax3, label=r'$F_{\rm max}$')
    ax3.set_xlabel('RA')
    ax3.set_ylabel('Dec')
    ax3.invert_xaxis()
    ax3.set_title(r'$F_{\rm max}$')

    ax4 = plt.subplot(1, 4, 4)
    sc = ax4.scatter(xdrp['ra'], xdrp['dec'], c=dw, cmap='viridis', vmin=x5, vmax=x95)
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(sc, cax=cax, label=r'$\delta\lambda$ ($\AA$)')
    ax4.set_xlabel('RA')
    ax4.set_ylabel('Dec')
    ax4.invert_xaxis()
    ax4.set_title('Wavelength Offset')
    plt.suptitle(r'Exposure %d  –  Fourier Wavelength Offsets  %d$-$%d $\AA$  (%d spectra)'
                 % (expnum, wmin, wmax, len(xflux)))
    plt.tight_layout()
    if outroot == '':
        outroot = 'fourier_offset'
    fig_dir = 'Fig_Qual'
    os.makedirs(fig_dir, exist_ok=True)
    plt.savefig('%s/%s_%d_%d.png' % (fig_dir, outroot, wmin, wmax))
    plt.close()
    return x5, x95


def doit(filename, wmin=3900, wmax=4000, do_plot=True):
    '''
    Analyze wavelength offsets in a single CFrame/SFrame file.

    Reads the FLUX, WAVE, and SLITMAP extensions, selects science fibers,
    computes per-fiber wavelength offsets via FFT cross-correlation in
    [wmin, wmax], and optionally saves a diagnostic PNG to Fig_Qual/.

    Returns a dict with the following keys (suitable for direct use as
    a row in an astropy Table)::

        Exposure             int   exposure number from the primary header
        mjd                  float modified Julian date
        dw_5pct              float 5th-percentile offset across all fibers (AA)
        dw_95pct             float 95th-percentile offset (AA)
        dw_med               float median offset, all fibers (AA)
        dw_med_sp1/2/3       float median offset per spectrograph (AA)
        dw_mad               float scaled MAD of offsets, all fibers (AA)
        dw_mad_sp1/2/3       float scaled MAD per spectrograph (AA)
        med_quality          float median cross-correlation quality [0,1]
        flux_max_med         float median of per-fiber max flux in [wmin,wmax]
        flux_max_med_sp1/2/3 float median peak flux per spectrograph
        flux_max_mad         float scaled MAD of per-fiber max flux
        flux_max_mad_sp1/2/3 float scaled MAD of peak flux per spectrograph
        dw_r*_sp*            float median offset per ring group/spectrograph
                             (see ring_col_names() for the full column list)
    '''
    x = fits.open(filename)
    expnum  = x[0].header.get('EXPOSURE', -1)
    obstime = x[0].header.get('OBSTIME', '')
    mjd     = convert_time(obstime, 'mjd') if obstime else np.nan
    outroot = os.path.basename(filename).replace('.fits', '')
    slit_tab = Table(x['SLITMAP'].data)
    slit_tab = clean_slitmap(slit_tab)
    wave = x['WAVE'].data
    flux = x['FLUX'].data[slit_tab['Row']]
    slit_tab['B'] = get_band_flux(wave, flux)
    dw, quality = estimate_wavelength_offsets(wave, flux, wmin, wmax)
    x5  = np.percentile(dw, 5)
    x95 = np.percentile(dw, 95)
    if do_plot==True:
        plot_offsets(dw, quality, xwave=wave, xflux=flux, xdrp=slit_tab,
                           wmin=wmin, wmax=wmax, outroot=outroot, expnum=expnum)
    dw_med      = np.median(dw)
    med_quality = np.median(quality)
    spec_id     = np.array(slit_tab['spectrographid'])
    dw_med_sp   = [np.median(dw[spec_id == sp]) for sp in (1, 2, 3)]
    ring_vals   = ring_sp_medians(dw, slit_tab)
    win_mask = (wave >= wmin) & (wave <= wmax)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        fmax_flux = np.nanmax(flux[:, win_mask], axis=1)
    mad = lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))
    flux_max_med    = np.median(fmax_flux)
    flux_max_med_sp = [np.median(fmax_flux[spec_id == sp]) for sp in (1, 2, 3)]
    dw_mad_sp       = [mad(dw[spec_id == sp]) for sp in (1, 2, 3)]
    flux_max_mad_sp = [mad(fmax_flux[spec_id == sp]) for sp in (1, 2, 3)]
    row = {
        'Exposure':           expnum,
        'mjd':                mjd,
        'dw_5pct':            x5,
        'dw_95pct':           x95,
        'dw_med':             dw_med,
        'dw_med_sp1':         dw_med_sp[0],
        'dw_med_sp2':         dw_med_sp[1],
        'dw_med_sp3':         dw_med_sp[2],
        'dw_mad':             mad(dw),
        'dw_mad_sp1':         dw_mad_sp[0],
        'dw_mad_sp2':         dw_mad_sp[1],
        'dw_mad_sp3':         dw_mad_sp[2],
        'med_quality':        med_quality,
        'flux_max_med':       flux_max_med,
        'flux_max_med_sp1':   flux_max_med_sp[0],
        'flux_max_med_sp2':   flux_max_med_sp[1],
        'flux_max_med_sp3':   flux_max_med_sp[2],
        'flux_max_mad':       mad(fmax_flux),
        'flux_max_mad_sp1':   flux_max_mad_sp[0],
        'flux_max_mad_sp2':   flux_max_mad_sp[1],
        'flux_max_mad_sp3':   flux_max_mad_sp[2],
    }
    row.update(zip(ring_col_names(), ring_vals))
    return row


def do_all(files, wmin=3900, wmax=4000, do_plot=False, drp_ver='1.2.1'):
    '''
    Process a list of CFrame/SFrame files, collect per-exposure wavelength-offset
    statistics, and write (or update) the summary FITS table.

    Calls doit() for each file, accumulates results, then writes
    Fourier_<wmin>_<wmax>.<date>.fits.  If that file already exists the new rows
    are merged in (duplicate exposures replaced) and the table is re-sorted by
    exposure number.

    Files that cannot be found directly are resolved by exposure number via
    the drpall table (drp_ver selects which DRP version to search).
    '''
    run_date = datetime.date.today().strftime('%y%m%d')
    lookup   = build_drp_lookup(drp_ver)
    rows = []

    ntot = len(files)
    for i, one in enumerate(files, 1):
        full = resolve_filename(one, lookup)
        if full is None:
            print('Skipping %s (%d/%d): file not found' % (one, i, ntot))
            continue
        print('Processing %s (%d/%d)' % (full, i, ntot))
        row = doit(full, wmin, wmax, do_plot)
        print('  95th pct offset: %.4f AA   median quality: %.3f' % (row['dw_95pct'], row['med_quality']))
        rows.append(row)

    if not rows:
        return

    summary = Table(rows)

    outname = 'Fourier_%d_%d.%s.fits' % (wmin, wmax, run_date)
    if os.path.isfile(outname):
        from astropy.table import vstack
        existing = Table.read(outname)
        existing = existing[~np.isin(existing['Exposure'], summary['Exposure'])]
        summary = vstack([existing, summary])
        summary.sort('Exposure')
    summary.write(outname, overwrite=True)
    print('Wrote %s' % outname)


def expnum_from_filename(filename):
    '''Extract the exposure number embedded in an lvmCFrame/SFrame filename, or -1.'''
    m = re.search(r'(\d+)\.fits', os.path.basename(filename))
    return int(m.group(1)) if m else -1


def filter_files(files, expmin=-1, expmax=-1, every=1):
    '''
    Filter a list of FITS filenames by exposure number range and stride.

    expmin/expmax are inclusive bounds on the exposure number parsed from
    the filename; -1 means no limit.  every=N keeps every Nth file after
    range filtering (1 = keep all).
    '''
    if expmin > 0:
        files = [f for f in files if expnum_from_filename(f) >= expmin]
    if expmax > 0:
        files = [f for f in files if expnum_from_filename(f) <= expmax]
    if every > 1:
        files = files[::every]
    return files


def steer(argv):
    wmin       = 3900
    wmax       = 4000
    files      = []
    input_list = ''
    do_plot    = False
    expmin     = -1
    expmax     = -1
    every      = 1
    drp_ver    = '1.2.1'

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-wmin':
            i += 1
            wmin = eval(argv[i])
        elif argv[i] == '-wmax':
            i += 1
            wmax = eval(argv[i])
        elif argv[i] == '-file':
            i += 1
            input_list = argv[i]
        elif argv[i] == '-plot':
            do_plot = True
        elif argv[i] == '-expmin':
            i += 1
            expmin = int(argv[i])
        elif argv[i] == '-expmax':
            i += 1
            expmax = int(argv[i])
        elif argv[i] == '-every':
            i += 1
            every = int(argv[i])
        elif argv[i] == '-ver':
            i += 1
            drp_ver = argv[i]
        elif argv[i][0] == '-':
            print('Error: Could not interpret command: ', argv[i])
            return
        elif argv[i].count('.fits'):
            files.append(argv[i])
        i += 1

    if input_list != '':
        try:
            xtab = Table.read(input_list)
        except Exception:
            print('Could not read : %s' % input_list)
            return
        if 'filename' in xtab.colnames:
            xfiles = list(xtab['filename'])
        elif 'Filename' in xtab.colnames:
            xfiles = list(xtab['Filename'])
        else:
            print('Read %s but could not find a column named filename or Filename' % input_list)
            return
        files = files + xfiles

    files = filter_files(files, expmin, expmax, every)

    if len(files) < 100:
        do_plot = True

    do_all(files, wmin, wmax, do_plot, drp_ver)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
