#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Analyze wavelength offsets in LVM CFrame/SFrame files by fiber,
using FFT cross-correlation of absorption features.

Command line usage (if any):

    usage: fourier_offset.py [-h] [-wmin 3900] [-wmax 4000] filename [filename ...]

    where filename is one or more lvmCFrame or lvmSFrame FITS files.
    -wmin and -wmax set the wavelength window used for the cross-correlation
    (default 3900-4000 AA).

Description:

    Reads the FLUX and WAVE extensions of each file, selects science fibers
    from the SLITMAP, then cross-correlates each spectrum against a median
    reference within [wmin, wmax] to estimate per-fiber wavelength offsets.
    Results are saved as a PNG showing the offset distribution and a spatial
    map of the offsets.

Primary routines:

    doit

Notes:

    The cross-correlation uses FFTs with parabolic sub-pixel refinement.
    A quality value in [0,1] is returned alongside each offset; values
    below ~0.5 indicate unreliable measurements (featureless window,
    bad pixels, low S/N).

History:

260417 ksl Coding begun; wavelength_offset.py incorporated into this file

'''

import sys
import warnings
from astropy.io import ascii, fits
from astropy.table import Table, join
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

    denom   = 2.0 * (2.0 * y1 - y0 - y2)
    sub_pix = np.where(valid & (np.abs(denom) > 1e-30), (y2 - y0) / denom, 0.0)
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
    id_min = np.searchsorted(xwave, wmin, side='right')
    id_max = np.searchsorted(xwave, wmax, side='right')
    xflux = xflux[:, id_min:id_max]
    band = np.nanmedian(xflux, axis=1)
    return band


def plot_offsets(dw, quality, xwave, xflux, xdrp, wmin=3900, wmax=4000, outroot='', expnum=-1):
    '''
    Plot the wavelength-offset diagnostic figure for a single exposure.

    Panels: (1) flux percentiles vs wavelength, (2) offset histogram,
    (3) spectra at the 5th/95th-percentile offset, (4) spatial offset map.
    Figure is saved to Fig_Qual/.
    '''
    x5  = np.percentile(dw, 5)
    x95 = np.percentile(dw, 95)

    plt.figure(1, (16, 4))
    plt.clf()

    plt.subplot(1, 4, 1)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        p5  = np.nanpercentile(xflux, 5,  axis=0)
        p50 = np.nanmedian(xflux,         axis=0)
        p95 = np.nanpercentile(xflux, 95, axis=0)
    plt.semilogy(xwave, p5,  label='5%')
    plt.semilogy(xwave, p50, label='50%')
    plt.semilogy(xwave, p95, label='95%')
    plt.legend()
    plt.xlim(wmin, wmax)
    ymin, ymax = plt.ylim()
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

    n5  = np.abs(dw - x5).argmin()
    n95 = np.abs(dw - x95).argmin()

    plt.subplot(1, 4, 3)
    plt.semilogy(xwave, xflux[n5],  label='Subtract: %.3f' % -x5)
    plt.semilogy(xwave, xflux[n95], label='       Add: %.3f' % x95)
    plt.xlim(wmin, wmax)
    plt.ylim(ymin, ymax)
    plt.legend()
    plt.xlabel(r'$\lambda$ ($\AA$)')
    plt.ylabel('Flux')
    plt.title('Spectra with offsets')

    ax4 = plt.subplot(1, 4, 4)
    sc = ax4.scatter(xdrp['ra'], xdrp['dec'], c=dw, cmap='viridis', vmin=x5, vmax=x95)
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(sc, cax=cax, label=r'$\delta\lambda$ ($\AA$)')
    ax4.set_xlabel('RA')
    ax4.set_ylabel('Dec')
    plt.suptitle('Exposure %d  –  Fourier Wavelength Offsets  (%d spectra)' % (expnum, len(xflux)))
    plt.tight_layout()
    if outroot == '':
        outroot = 'fourier_offset'
    fig_dir = 'Fig_Qual'
    os.makedirs(fig_dir, exist_ok=True)
    plt.savefig('%s/%s_%d_%d.png' % (fig_dir, outroot, wmin, wmax))
    plt.close()
    return x5, x95


def doit(filename, wmin=3900, wmax=4000):
    '''
    Analyze wavelength offsets in a single CFrame/SFrame file.

    Reads the file, computes per-fiber wavelength offsets via FFT
    cross-correlation, produces the diagnostic plot, and returns the
    exposure number, offset arrays, and summary percentiles.
    '''
    x = fits.open(filename)
    expnum  = x[0].header.get('EXPOSURE', -1)
    outroot = os.path.basename(filename).replace('.fits', '')
    slit_tab = Table(x['SLITMAP'].data)
    slit_tab = clean_slitmap(slit_tab)
    wave = x['WAVE'].data
    flux = x['FLUX'].data[slit_tab['Row']]
    slit_tab['B'] = get_band_flux(wave, flux)
    dw, quality = estimate_wavelength_offsets(wave, flux, wmin, wmax)
    x5, x95 = plot_offsets(dw, quality, xwave=wave, xflux=flux, xdrp=slit_tab,
                           wmin=wmin, wmax=wmax, outroot=outroot, expnum=expnum)
    dw_med  = np.median(dw)
    spec_id = np.array(slit_tab['spectrographid'])
    dw_med_sp = [np.median(dw[spec_id == sp]) for sp in (1, 2, 3)]
    ring_vals = ring_sp_medians(dw, slit_tab)
    return expnum, x5, x95, dw_med, dw_med_sp[0], dw_med_sp[1], dw_med_sp[2], ring_vals


def steer(argv):
    wmin  = 3900
    wmax  = 4000
    files = []

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
        elif argv[i][0] == '-':
            print('Error: Could not interpret command: ', argv[i])
            return
        elif argv[i].count('.fits'):
            files.append(argv[i])
        i += 1

    col_names = ring_col_names()
    exposures, x5_vals, x95_vals = [], [], []
    med_vals, med_sp1_vals, med_sp2_vals, med_sp3_vals = [], [], [], []
    ring_cols = {name: [] for name in col_names}

    ntot = len(files)
    for i, one in enumerate(files, 1):
        print('Processing %s (%d/%d)' % (one, i, ntot))
        expnum, x5, x95, dw_med, dw_sp1, dw_sp2, dw_sp3, ring_vals = doit(one, wmin, wmax)
        exposures.append(expnum)
        x5_vals.append(x5)
        x95_vals.append(x95)
        med_vals.append(dw_med)
        med_sp1_vals.append(dw_sp1)
        med_sp2_vals.append(dw_sp2)
        med_sp3_vals.append(dw_sp3)
        for name, val in zip(col_names, ring_vals):
            ring_cols[name].append(val)

    if exposures:
        summary = Table([exposures, x5_vals, x95_vals,
                         med_vals, med_sp1_vals, med_sp2_vals, med_sp3_vals],
                        names=['Exposure', 'dw_5pct', 'dw_95pct',
                               'dw_med', 'dw_med_sp1', 'dw_med_sp2', 'dw_med_sp3'])
        for name in col_names:
            summary[name] = ring_cols[name]
        for col in summary.colnames[1:]:
            summary[col].format = '%.4f'
        outname = 'Fourier_offsets_%d_%d.txt' % (wmin, wmax)
        summary.write(outname, format='ascii.fixed_width_two_line', overwrite=True)
        print('Wrote %s' % outname)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
