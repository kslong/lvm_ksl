#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Evaluate sky subtraction quality for one or more FITS files produced
    by SkySubOrig.py, SkySubDrp.py, SkySubDev1.py, or SkySubDev2.py.
    All four scripts write files with WAVE, FLUX (sky-subtracted), and
    SKY extensions that this script reads directly.  When multiple files
    are supplied they are overlaid in the same figures, making it easy
    to compare sky subtraction methods on a common footing.

    Writes a single interactive HTML file containing four Plotly figures
    and an inline statistics table.

Command line usage (if any):

    usage: SkySub_eval.py [wmin wmax] [-num N] [-out outroot]
                          filename [filename ...]

    Arguments:

    filename         one or more SkySub output FITS files to evaluate.
    wmin             minimum wavelength in Angstroms (default 3600).
    wmax             maximum wavelength in Angstroms (default 9800).

    Options:

    -num N           overlay N randomly-selected individual spectra (default 20).
    -out outroot     combine all files into one HTML (default: per-file <stem>_eval.html).

Description:

    **Figure 1 — spectral overview (3 panels, shared x-axis, linear scale):**

    Panel 1 (Original Flux): FLUX + SKY, showing the pre-subtraction
    spectrum.  Median and 10th/90th percentile band, with N individual
    spectra overlaid in grey.  One colour per file when multiple files
    are supplied.

    Panel 2 (Sky-subtracted Flux): FLUX extension on the same linear
    scale.

    Panel 3 (Sky Model): SKY extension on the same scale.

    All three panels share a common y-range of -1e-14 to 1e-13
    erg s⁻¹ cm⁻² Å⁻¹, set to make sky-subtracted residuals visible.
    The x-axis is downsampled to ~2000 display pixels for performance.

    **Figure 2 — residual histograms at three diagnostic sky-line windows:**

    Diagnostic windows: [OI] 5577 (5560–5594 Å), [OI] 6300
    (6280–6320 Å), and IR OH (9300–9500 Å).  For each window the
    distribution of all per-pixel sky-subtracted flux values across every
    spectrum is plotted as a bar histogram.  A Gaussian with the same
    median and NMAD is overlaid as a dotted curve.  A statistics box
    (N, median, NMAD, skewness) is placed inside each panel; boxes stack
    vertically when multiple files are given.  The histogram range is
    clipped to median ± 5·NMAD.

    **Figure 3 (Figure 4 in HTML) — diagnostic window median spectra:**

    For each of the three diagnostic windows a wider search region
    (5400–5750, 6100–6500, 9000–9800 Å) is plotted showing the median
    and 10th/90th percentile band of the original (dotted) and
    sky-subtracted (solid) spectra.  Two types of shading show which
    pixels are used by the HF RMS analysis:

    - Red shading — signal window (diagnostic window, used to measure
      sky-line RMS)
    - Green shading — background pixels (mask-selected sky-line-free
      pixels in the search region, used to estimate noise RMS)

    The sky-line mask is read from sky_mask.fits (palace_make_mask.py
    output) if available; otherwise all pixels outside the diagnostic
    window are treated as background.

    **Statistics table:**

    An HTML table between Figure 3 and Figure 4 reports per-file per-window:
    noise RMS (median background RMS), sky-line RMS before/after
    subtraction, and the 50th, 90th, and 95th percentiles of the
    HF RMS ratio (sky_sub / sky_orig) across all spectra, plus the
    fraction of spectra with ratio < 0.5.  The same information is
    printed to the terminal.

    **Figure 4 (Figure 3 in HTML) — HF RMS ratio per spectrum:**

    For each diagnostic window the noise-corrected high-frequency RMS
    ratio (sky-subtracted / original) is plotted against spectrum index.
    The ratio is computed as follows:

    1. A Gaussian smooth (σ = 50 pixels ≈ 25 Å) using only mask-selected
       clean pixels is subtracted from each spectrum to produce a
       high-frequency (HF) residual array that retains sky-line-scale
       (≤ LSF FWHM) structure.
    2. The RMS of the HF residuals in the diagnostic window gives
       ``rms_diag`` (sky lines + noise).
    3. The RMS of the HF residuals in the background region gives
       ``rms_bg`` (noise only).
    4. The sky-line-only amplitude is ``sky_rms = √max(0, rms_diag² −
       rms_bg²)``.
    5. The ratio ``sky_sub / sky_orig`` is plotted.  A value of 0 means
       the sky lines were completely removed; a value of 1 means no
       change.  The panels share a linked x-axis (spectrum index) so
       zooming one panel pans all three.

Notes:

    Requires plotly (pip install plotly) and scipy.
    If GetSkyCont.py is on the Python path, sky_mask.fits is loaded
    automatically from the current directory or the lvm_ksl data/
    directory; without the mask the Gaussian smooth is unweighted.

History:

    260630 ksl  Initial version (matplotlib, mode-based: eval/xplot/bigplot).
    260630 ksl  Rewritten to Plotly HTML with multi-file overlay.
    260630 ksl  Added HF RMS ratio analysis, diagnostic window spectra, stats table.

'''

import sys
import os
import warnings
import numpy as np
from pathlib import Path
from astropy.io import fits
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from GetSkyCont import load_mask, _interp_mask_to_wave
    _HAVE_MASK = True
except ImportError:
    _HAVE_MASK = False


# Diagnostic sky-line windows: (label, win_lo, win_hi, search_lo, search_hi)
# The search range is the broader region scanned for clean background pixels
# (mask=1) outside the diagnostic window itself.
_DIAG_WINDOWS = [
    ('OI 5577',  5560.0, 5594.0,  5400.0, 5750.0),
    ('OI 6300',  6280.0, 6320.0,  6100.0, 6500.0),
    ('IR 9300',  9300.0, 9500.0,  9000.0, 9800.0),
]

# Per-file colour pairs: (solid line, shaded band)
_FILE_COLORS = [
    ('rgb(31,119,180)',   'rgba(31,119,180,0.20)'),
    ('rgb(214,39,40)',    'rgba(214,39,40,0.20)'),
    ('rgb(44,160,44)',    'rgba(44,160,44,0.20)'),
    ('rgb(148,103,189)', 'rgba(148,103,189,0.20)'),
    ('rgb(255,127,14)',  'rgba(255,127,14,0.20)'),
    ('rgb(23,190,207)',  'rgba(23,190,207,0.20)'),
]

_USAGE = '''Usage: SkySub_eval.py [wmin wmax] [-num N] [-out outroot] filename [filename ...]

Arguments:
  filename   one or more SkySub FITS files (SkySubOrig/Drp/Dev1/Dev2 output)
  wmin       spectral overview minimum wavelength in Angstroms (default 3600)
  wmax       spectral overview maximum wavelength in Angstroms (default 9800)

Options:
  -num N        number of individual spectra to overlay (default 20; 0 = band only)
  -out outroot  overlay all files into one HTML with this output root;
                without -out, each file produces its own <stem>_eval.html
'''


# ──────────────────────────────────────────────────────────────
# Analysis helpers
# ──────────────────────────────────────────────────────────────

def _window_stats(flux, wave, wmin, wmax):
    '''Pixel-level residual statistics within a wavelength window.

    Parameters
    ----------
    flux : 2-D array (n_spec, n_wave), sky-subtracted flux
    wave : 1-D array (n_wave,)
    wmin, wmax : float

    Returns
    -------
    dict with keys n, med, nmad, skew, p10, p90, vals — or None if no data.
    '''
    sel  = (wave >= wmin) & (wave <= wmax)
    if not sel.any():
        return None
    vals = flux[:, sel].astype(float).flatten()
    vals = vals[np.isfinite(vals)]
    if len(vals) < 10:
        return None
    med  = float(np.median(vals))
    nmad = float(1.4826 * np.median(np.abs(vals - med)))
    mn   = float(np.mean(vals))
    sd   = float(np.std(vals, ddof=0))
    skew = float(np.mean(((vals - mn) / sd) ** 3)) if sd > 0 else 0.0
    p10, p90 = (float(v) for v in np.percentile(vals, [10, 90]))
    return dict(n=len(vals), med=med, nmad=nmad, skew=skew,
                p10=p10, p90=p90, vals=vals)


def _hf_residual(flux_2d, clean_mask=None, sigma_pix=50):
    '''High-frequency residual: each spectrum minus a Gaussian-smoothed version.

    sigma_pix=50 ≈ 25 Å at LVM sampling, well above the 1.3 Å LSF FWHM,
    so the residual retains only sky-line-scale and finer structure.

    clean_mask : 1-D boolean array (n_wave,), True = sky-line-free pixel.
        Sky-line pixels are zeroed out before smoothing so their flux does
        not leak into the continuum estimate (same masking as GetSkyCont).
        If None, all pixels are treated as clean.

    Uses a weighted Gaussian smooth so masked pixels are reconstructed from
    their neighbours rather than pulled toward zero.
    '''
    filled = np.where(np.isfinite(flux_2d), flux_2d, 0.0)
    if clean_mask is None:
        smooth = gaussian_filter1d(filled, sigma=sigma_pix, axis=1)
    else:
        w = clean_mask.astype(float)               # (n_wave,)
        w2d = np.broadcast_to(w, filled.shape).copy()
        num = gaussian_filter1d(filled * w2d, sigma=sigma_pix, axis=1)
        den = gaussian_filter1d(w2d,           sigma=sigma_pix, axis=1)
        smooth = num / np.where(den > 1e-6, den, 1.0)
    hf = flux_2d - smooth
    return hf


def _per_spec_hf_rms(hf_2d, wave, wmin, wmax):
    '''Per-spectrum RMS of high-frequency residuals within a wavelength window.

    RMS = sqrt(mean(hf^2)) responds to squared amplitude, so even a few
    bright sky-line pixels dominate over the many clean pixels — unlike NMAD
    which is suppressed by the majority of non-line pixels.  The quadrature
    formula rms_sky = sqrt(rms_diag^2 - rms_clean^2) then isolates the
    sky-line contribution.  Equally sensitive to positive and negative
    deviations.

    Returns 1-D array of length n_spec, or None if no pixels fall in window.
    '''
    sel = (wave >= wmin) & (wave <= wmax)
    if not sel.any():
        return None
    r = hf_2d[:, sel].astype(float)
    r[~np.isfinite(r)] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        rms = np.sqrt(np.nanmean(r ** 2, axis=1))
    return rms


def _per_spec_hf_rms_sel(hf_2d, sel):
    '''Per-spectrum RMS of HF residuals for an arbitrary pixel selection.

    sel : 1-D boolean array (n_wave,)
    Returns 1-D array of length n_spec, or None if sel is empty.
    '''
    if not sel.any():
        return None
    r = hf_2d[:, sel].astype(float)
    r[~np.isfinite(r)] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        rms = np.sqrt(np.nanmean(r ** 2, axis=1))
    return rms


# ──────────────────────────────────────────────────────────────
# Plotly trace helpers
# ──────────────────────────────────────────────────────────────

def _band_traces(wave, arr, name, color_med, color_band,
                 legend_ref='legend', show_legend=True):
    '''Median line + shaded 10/90-pct band traces.'''
    p10 = np.nanpercentile(arr, 10, axis=0)
    p90 = np.nanpercentile(arr, 90, axis=0)
    med = np.nanmedian(arr, axis=0)
    return [
        go.Scatter(x=wave, y=p10, mode='lines', line=dict(width=0),
                   showlegend=False, hoverinfo='skip'),
        go.Scatter(x=wave, y=p90, mode='lines', line=dict(width=0),
                   fill='tonexty', fillcolor=color_band,
                   name=f'{name} 10–90 pct', showlegend=False,
                   hoverinfo='skip'),
        go.Scatter(x=wave, y=med, mode='lines',
                   line=dict(color=color_med, width=1.5),
                   name=name, showlegend=show_legend, legend=legend_ref),
    ]


def _sample_traces(wave, arr, idx):
    '''Semi-transparent individual spectrum traces.'''
    if len(idx) == 0:
        return []
    alpha = min(0.7, 3.0 / max(1.0, len(idx) ** 0.5))
    return [
        go.Scatter(x=wave, y=arr[i], mode='lines',
                   line=dict(color=f'rgba(120,120,120,{alpha:.2f})', width=0.6),
                   showlegend=False, hoverinfo='skip')
        for i in idx
    ]


def _fmt(v):
    '''Format a flux value compactly in scientific notation.'''
    if not np.isfinite(v) or abs(v) < 1e-30:
        return '0'
    return f'{v:.2e}'


# ──────────────────────────────────────────────────────────────
# Main evaluation function
# ──────────────────────────────────────────────────────────────

def plot_eval(filenames, wmin=3600.0, wmax=9800.0, n_sample=20, outroot=''):
    '''Build a three-figure Plotly HTML evaluation file.

    Parameters
    ----------
    filenames : str or list of str
        SkySub output FITS file(s).
    wmin, wmax : float
        Wavelength range for the spectral overview panels.
    n_sample : int
        Number of randomly-selected spectra to overlay on band panels.
    outroot : str
        Output filename root.  Default: <first stem>_eval.
    '''
    if isinstance(filenames, str):
        filenames = [filenames]

    # Validate up front: warn and drop missing or unreadable files
    valid = []
    for fn in filenames:
        if not Path(fn).exists():
            print(f'Warning: file not found, skipping: {fn}')
        else:
            valid.append(fn)
    filenames = valid
    if not filenames:
        print('Error: no valid input files.')
        return

    # Load sky-line mask (same as GetSkyCont) to avoid leaking lines into smooth.
    # Raw mask is interpolated to each file's wave grid inside the per-file loop.
    _raw_mask = None
    if _HAVE_MASK:
        _data_dir = Path(__file__).parent.parent / 'data'
        for _candidate in [Path('sky_mask.fits'), _data_dir / 'sky_mask.fits']:
            if _candidate.exists():
                try:
                    _raw_mask = load_mask(str(_candidate))  # (mwave, marr)
                    print(f'Loaded sky mask: {_candidate}')
                except Exception as _e:
                    print(f'Warning: could not load mask {_candidate}: {_e}')
                break
    if _raw_mask is None:
        print('Warning: sky_mask.fits not found; smoothing without sky-line mask')

    rng = np.random.default_rng(42)

    _ax = dict(showline=True, linewidth=1, linecolor='black',
               mirror=True, ticks='outside', ticklen=4, showticklabels=True)
    _leg_style = dict(xanchor='right', yanchor='top', x=0.99,
                      bgcolor='rgba(255,255,255,0.7)')

    # ══════════════════════════════════════════════════════════
    # Figure 1 — spectral overview (orig flux + sky-sub + sky, 3 rows)
    # ══════════════════════════════════════════════════════════
    fig1 = make_subplots(
        rows=3, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.04,
        subplot_titles=[
            f'Original Flux  (sky-sub + sky model)  ({len(filenames)} file(s))',
            'Sky-subtracted Flux',
            'Sky Model',
        ],
    )

    # ══════════════════════════════════════════════════════════
    # Figure 2 — per-window residual histograms (1 row × 3 cols)
    # ══════════════════════════════════════════════════════════
    fig2 = make_subplots(
        rows=1, cols=3,
        horizontal_spacing=0.08,
        subplot_titles=[
            f'{lbl}  ({int(lo)}–{int(hi)} Å)'
            for lbl, lo, hi, *_ in _DIAG_WINDOWS
        ],
    )

    # ══════════════════════════════════════════════════════════
    # Figure 3 — per-spectrum HF RMS ratio at each diagnostic window
    # ══════════════════════════════════════════════════════════
    fig3 = make_subplots(
        rows=1, cols=3,
        horizontal_spacing=0.08,
        subplot_titles=[f'{lbl}  HF RMS ratio (sub/orig)' for lbl, *_ in _DIAG_WINDOWS],
    )

    # ══════════════════════════════════════════════════════════
    # Figure 4 — median spectra per diagnostic window with shading
    # ══════════════════════════════════════════════════════════
    fig4 = make_subplots(
        rows=1, cols=3,
        horizontal_spacing=0.08,
        subplot_titles=[f'{lbl}  ({int(slo)}–{int(shi)} Å)'
                        for lbl, *_, slo, shi in _DIAG_WINDOWS],
    )

    all_n_spec   = []
    all_hf_ratio = {lbl: [] for lbl, *_ in _DIAG_WINDOWS}  # for y-range
    _fig4_shaded = False   # add background/signal shading once (first file)
    _stats_rows  = []      # collected stats for screen + HTML output

    for k, filename in enumerate(filenames):
        c_med, c_band = _FILE_COLORS[k % len(_FILE_COLORS)]
        # Semi-transparent solid colour for histogram bars
        c_hist = c_med.replace('rgb(', 'rgba(').replace(')', ',0.55)')
        label  = Path(filename).stem

        try:
            hdul = fits.open(filename)
            wave = hdul['WAVE'].data.astype(float)
            flux = hdul['FLUX'].data.astype(float)
            sky  = hdul['SKY'].data.astype(float)
            hdul.close()
        except Exception as e:
            print(f'Warning: could not read {filename} ({e}), skipping.')
            continue

        if flux.ndim == 1:
            flux = flux[np.newaxis, :]
            sky  = sky[np.newaxis, :]

        n_spec = len(flux)
        all_n_spec.append(n_spec)

        # Wavelength window + downsample to ~2000 display pixels for Figure 1
        wmask  = (wave >= wmin) & (wave <= wmax)
        ds     = max(1, int(wmask.sum()) // 2000)
        pidx   = np.where(wmask)[0][::ds]
        w      = wave[pidx]
        f_disp = flux[:, pidx]
        s_disp = sky[:,  pidx]

        sample_idx = (rng.choice(n_spec, size=min(n_sample, n_spec), replace=False)
                      if n_sample > 0 else [])

        # ── Figure 1, Panel 1: original flux (FLUX+SKY, linear) ─
        orig_disp = f_disp + s_disp
        for tr in _sample_traces(w, orig_disp, sample_idx):
            fig1.add_trace(tr, row=1, col=1)
        for tr in _band_traces(w, orig_disp, label, c_med, c_band,
                               legend_ref='legend', show_legend=True):
            fig1.add_trace(tr, row=1, col=1)

        # ── Figure 1, Panel 2: sky-subtracted FLUX (linear) ───
        for tr in _sample_traces(w, f_disp, sample_idx):
            fig1.add_trace(tr, row=2, col=1)
        for tr in _band_traces(w, f_disp, label, c_med, c_band,
                               legend_ref='legend2', show_legend=True):
            fig1.add_trace(tr, row=2, col=1)

        # ── Figure 1, Panel 3: SKY (linear) ───────────────────
        for tr in _band_traces(w, s_disp, label, c_med, c_band,
                               legend_ref='legend3', show_legend=True):
            fig1.add_trace(tr, row=3, col=1)

        # ── Figure 2: diagnostic histograms ───────────────────
        for col_idx, (win_lbl, win_lo, win_hi, *_) in enumerate(_DIAG_WINDOWS, start=1):
            st = _window_stats(flux, wave, win_lo, win_hi)
            if st is None:
                continue

            lo = st['med'] - 5.0 * st['nmad']
            hi = st['med'] + 5.0 * st['nmad']
            if lo >= hi:
                lo, hi = st['p10'] - 1e-20, st['p90'] + 1e-20

            counts, edges = np.histogram(st['vals'], bins=60, range=(lo, hi))
            centers = 0.5 * (edges[:-1] + edges[1:])
            bw      = edges[1] - edges[0]

            fig2.add_trace(
                go.Bar(x=centers, y=counts,
                       name=label, legendgroup=label,
                       showlegend=(col_idx == 1),
                       marker_color=c_hist, marker_line_width=0),
                row=1, col=col_idx,
            )

            # Gaussian overlay (dotted, same solid colour)
            xg  = np.linspace(lo, hi, 300)
            sig = max(st['nmad'], 1e-30)
            yg  = (st['n'] * bw / (np.sqrt(2 * np.pi) * sig)
                   * np.exp(-0.5 * ((xg - st['med']) / sig) ** 2))
            fig2.add_trace(
                go.Scatter(x=xg, y=yg, mode='lines',
                           line=dict(color=c_med, width=1.5, dash='dot'),
                           showlegend=False, hoverinfo='skip'),
                row=1, col=col_idx,
            )

            # Stats annotation; stack vertically for multiple files
            ax_sfx = '' if col_idx == 1 else str(col_idx)
            y_top  = 0.97 - k * 0.23
            if y_top > 0.08:
                fig2.add_annotation(
                    x=0.97, y=y_top,
                    xref=f'x{ax_sfx} domain', yref=f'y{ax_sfx} domain',
                    xanchor='right', yanchor='top',
                    showarrow=False, align='right',
                    font=dict(size=9, family='monospace'),
                    bgcolor='rgba(255,255,255,0.80)',
                    bordercolor=c_med, borderwidth=1,
                    text=(f'<b>{label}</b><br>'
                          f'N    = {st["n"]}<br>'
                          f'Med  = {_fmt(st["med"])}<br>'
                          f'NMAD = {_fmt(st["nmad"])}<br>'
                          f'Skew = {st["skew"]:.2f}'),
                )

        # ── Figures 3 & 4: HF RMS ratio and median spectra ────
        # Build clean_mask on the full wave grid for this file
        clean_mask_full = (_interp_mask_to_wave(_raw_mask[0], _raw_mask[1], wave)
                           if _raw_mask is not None and _HAVE_MASK else None)
        # Also on the downsampled grid used for Figure 1
        clean_mask = (_interp_mask_to_wave(_raw_mask[0], _raw_mask[1], w)
                      if _raw_mask is not None and _HAVE_MASK else None)

        hf_orig = _hf_residual(f_disp + s_disp, clean_mask=clean_mask)
        hf_sub  = _hf_residual(f_disp,           clean_mask=clean_mask)
        spec_x  = np.arange(n_spec)

        for col_idx, (win_lbl, win_lo, win_hi, slo, shi) in \
                enumerate(_DIAG_WINDOWS, start=1):

            # ── pixel selections on the downsampled grid ──────
            diag_sel  = (w >= win_lo) & (w <= win_hi)
            # background: search range outside diagnostic window, mask=clean
            bg_sel = (w >= slo) & (w <= shi) & ~diag_sel
            if clean_mask is not None:
                bg_sel &= clean_mask

            rms_orig_diag = _per_spec_hf_rms_sel(hf_orig, diag_sel)
            rms_sub_diag  = _per_spec_hf_rms_sel(hf_sub,  diag_sel)
            rms_orig_bg   = _per_spec_hf_rms_sel(hf_orig, bg_sel)
            rms_sub_bg    = _per_spec_hf_rms_sel(hf_sub,  bg_sel)
            if rms_orig_diag is None or rms_sub_diag is None:
                continue
            noise_orig = rms_orig_bg if rms_orig_bg is not None else np.zeros(n_spec)
            noise_sub  = rms_sub_bg  if rms_sub_bg  is not None else np.zeros(n_spec)
            sky_orig = np.sqrt(np.maximum(0.0, rms_orig_diag**2 - noise_orig**2))
            sky_sub  = np.sqrt(np.maximum(0.0, rms_sub_diag**2  - noise_sub**2))
            with np.errstate(invalid='ignore', divide='ignore'):
                ratio = np.where(sky_orig > 0, sky_sub / sky_orig, np.nan)
            good = ratio[np.isfinite(ratio) & (ratio >= 0)]
            all_hf_ratio[win_lbl].extend(good.tolist())

            # Figure 3 trace
            hover = [f'spec {i}<br>ratio = {ratio[i]:.3f}<br>'
                     f'sky orig = {_fmt(sky_orig[i])}<br>'
                     f'sky sub  = {_fmt(sky_sub[i])}'
                     for i in range(n_spec)]
            fig3.add_trace(
                go.Scatter(
                    x=spec_x, y=ratio,
                    mode='lines',
                    line=dict(color=c_med, width=1.5),
                    name=label, legendgroup=label,
                    showlegend=(col_idx == 1),
                    text=hover,
                    hovertemplate='%{text}<extra></extra>',
                ),
                row=1, col=col_idx,
            )

            # Collect statistics for screen + HTML output
            noise_med    = float(np.nanmedian(noise_sub))
            sky_orig_med = float(np.nanmedian(sky_orig))
            sky_sub_med  = float(np.nanmedian(sky_sub))
            if len(good) > 3:
                p50, p90, p95 = np.nanpercentile(good, [50, 90, 95])
                frac_lt_half  = float(np.mean(good < 0.5)) * 100
            else:
                p50 = p90 = p95 = frac_lt_half = float('nan')
            n_bg = int(bg_sel.sum())
            _stats_rows.append(dict(
                file=label, window=win_lbl,
                n_spec=int(n_spec), n_valid=len(good), n_bg_pix=n_bg,
                noise_med=noise_med,
                sky_orig_med=sky_orig_med, sky_sub_med=sky_sub_med,
                p50=p50, p90=p90, p95=p95, frac_lt_half=frac_lt_half,
            ))

            # ── Figure 4: median spectra in search range ──────
            # Use full-resolution arrays for the search range
            s4_sel = (wave >= slo) & (wave <= shi)
            if not s4_sel.any():
                continue
            w4      = wave[s4_sel]
            orig4   = (flux + sky)[:, s4_sel]
            sub4    = flux[:, s4_sel]

            for arr4, dash in [(orig4, 'dot'), (sub4, 'solid')]:
                p10_4 = np.nanpercentile(arr4, 10, axis=0)
                p90_4 = np.nanpercentile(arr4, 90, axis=0)
                med4  = np.nanmedian(arr4, axis=0)
                fig4.add_trace(go.Scatter(
                    x=w4, y=p10_4, mode='lines', line=dict(width=0),
                    showlegend=False, hoverinfo='skip'), row=1, col=col_idx)
                fig4.add_trace(go.Scatter(
                    x=w4, y=p90_4, mode='lines', line=dict(width=0),
                    fill='tonexty', fillcolor=c_band,
                    showlegend=False, hoverinfo='skip'), row=1, col=col_idx)
                fig4.add_trace(go.Scatter(
                    x=w4, y=med4, mode='lines',
                    line=dict(color=c_med, width=1.2, dash=dash),
                    showlegend=False, hoverinfo='skip'), row=1, col=col_idx)

            # Add shading on first file only (shapes are per-figure)
            if not _fig4_shaded:
                # diagnostic window — red/pink
                fig4.add_vrect(x0=win_lo, x1=win_hi,
                               fillcolor='rgba(220,50,50,0.30)',
                               line_width=1, line_color='rgba(180,0,0,0.5)',
                               row=1, col=col_idx)
                # background: mask-selected clean pixels outside diagnostic window
                if clean_mask_full is not None:
                    bg_full = (wave >= slo) & (wave <= shi) & \
                              ~((wave >= win_lo) & (wave <= win_hi)) & \
                              clean_mask_full
                else:
                    bg_full = (wave >= slo) & (wave <= shi) & \
                              ~((wave >= win_lo) & (wave <= win_hi))
                # shade contiguous clean runs
                idx = np.where(np.diff(np.concatenate([[False], bg_full, [False]])))[0]
                for r0, r1 in zip(idx[::2], idx[1::2]):
                    fig4.add_vrect(x0=wave[r0], x1=wave[r1 - 1],
                                   fillcolor='rgba(30,160,30,0.25)',
                                   line_width=0, row=1, col=col_idx)

        _fig4_shaded = True   # only add shapes once

    # ── Post-loop axis clipping ────────────────────────────────
    # Figure 3: linear y-range 0 to 99th-pct + margin; reference line at 1
    for col_idx, (win_lbl, *_) in enumerate(_DIAG_WINDOWS, start=1):
        vals = all_hf_ratio[win_lbl]
        if vals:
            v = np.array(vals)
            v = v[np.isfinite(v) & (v >= 0)]
            if len(v):
                y_hi = float(np.nanpercentile(v, 99)) + 0.05
                fig3.update_yaxes(range=[-0.05, y_hi], row=1, col=col_idx)
        # reference line at ratio = 1 (no improvement)
        fig3.add_hline(y=1.0, line_dash='dot', line_color='grey',
                       line_width=1, row=1, col=col_idx)

    # ── Figure 1 layout ───────────────────────────────────────
    _lin_range = [-1e-14, 1e-13]   # common linear scale for all three panels

    _title_stem = Path(filenames[0]).stem if len(filenames) == 1 \
                  else f'{len(filenames)} files'
    fig1.update_layout(
        title=f'Sky Subtraction Evaluation — {_title_stem}  {int(wmin)}–{int(wmax)} Å',
        height=1575, template='simple_white',
        legend  = dict(**_leg_style, y=0.99),
        legend2 = dict(**_leg_style, y=0.66),
        legend3 = dict(**_leg_style, y=0.33),
    )
    fig1.update_xaxes(**_ax)
    fig1.update_yaxes(**_ax, exponentformat='e', showexponent='all')
    fig1.update_xaxes(range=[wmin, wmax])
    fig1.update_xaxes(showticklabels=False, row=1, col=1)
    fig1.update_xaxes(showticklabels=False, row=2, col=1)
    fig1.update_xaxes(title_text='Wavelength (Å)', row=3, col=1)
    fig1.update_yaxes(title_text='Flux',         range=_lin_range, row=1, col=1)
    fig1.update_yaxes(title_text='Sky-sub Flux', range=_lin_range, row=2, col=1)
    fig1.update_yaxes(title_text='Sky',          range=_lin_range, row=3, col=1)

    # ── Figure 2 layout ───────────────────────────────────────
    fig2.update_layout(
        height=480, template='simple_white',
        barmode='overlay',
        legend=dict(**_leg_style, y=0.99),
    )
    fig2.update_xaxes(**_ax, exponentformat='e', showexponent='all',
                      title_text='Sky-subtracted flux')
    fig2.update_yaxes(**_ax, exponentformat='none', title_text='N')

    # ── Figure 3 layout ───────────────────────────────────────
    max_n = max(all_n_spec) if all_n_spec else 1
    fig3.update_layout(
        height=420, template='simple_white',
        legend=dict(**_leg_style, y=0.99),
    )
    fig3.update_xaxes(**_ax, title_text='Spectrum index',
                      range=[-0.5, max_n - 0.5], matches='x')
    fig3.update_yaxes(**_ax, title_text='HF RMS ratio (sub/orig)')

    # ── Figure 4 layout ───────────────────────────────────────
    fig4.update_layout(
        height=420, template='simple_white',
        showlegend=False,
    )
    fig4.update_xaxes(**_ax, title_text='Wavelength (Å)')
    fig4.update_yaxes(**_ax, exponentformat='e', showexponent='all',
                      title_text='Flux')

    # ── Print statistics to screen and build HTML stats block ─
    _hdr = (f'{"File":<55} {"Window":<10} {"N":>5} {"Noise":>10} '
            f'{"Sky orig":>10} {"Sky sub":>10} {"Med":>6} '
            f'{"p90":>6} {"p95":>6} {"<50%":>6}')
    print()
    print(_hdr)
    print('-' * len(_hdr))
    _html_rows = []
    for r in _stats_rows:
        line = (f'{r["file"]:<55} {r["window"]:<10} {r["n_valid"]:>5} '
                f'{r["noise_med"]:>10.2e} {r["sky_orig_med"]:>10.2e} '
                f'{r["sky_sub_med"]:>10.2e} {r["p50"]:>6.3f} '
                f'{r["p90"]:>6.3f} {r["p95"]:>6.3f} {r["frac_lt_half"]:>5.0f}%')
        print(line)
        _html_rows.append(
            f'<tr><td>{r["file"]}</td><td>{r["window"]}</td>'
            f'<td>{r["n_valid"]}</td>'
            f'<td>{r["noise_med"]:.2e}</td>'
            f'<td>{r["sky_orig_med"]:.2e}</td>'
            f'<td>{r["sky_sub_med"]:.2e}</td>'
            f'<td>{r["p50"]:.3f}</td><td>{r["p90"]:.3f}</td>'
            f'<td>{r["p95"]:.3f}</td><td>{r["frac_lt_half"]:.0f}%</td></tr>'
        )
    print()

    _stats_html = '''
<div style="font-family:monospace; font-size:13px; margin:10px 20px;">
<h3>Sky subtraction diagnostics — HF RMS ratio (sky-subtracted / original)</h3>
<p>Noise: median per-spectrum background RMS (clean pixels near each window).<br>
Sky orig/sub: median noise-corrected sky-line RMS before/after subtraction.<br>
Ratio columns: 50th (median), 90th, 95th percentile across spectra.<br>
&lt;50%: fraction of spectra where the ratio is below 0.5.</p>
<table border="1" cellpadding="4" cellspacing="0" style="border-collapse:collapse;">
<tr style="background:#ddd;">
  <th>File</th><th>Window</th><th>N valid</th>
  <th>Noise</th><th>Sky orig</th><th>Sky sub</th>
  <th>p50</th><th>p90</th><th>p95</th><th>&lt;50%</th>
</tr>
''' + '\n'.join(_html_rows) + '\n</table></div>\n'

    _fig4_html_header = '''
<div style="margin:10px 20px;">
<b>Diagnostic window spectra</b> — median with 10–90% band.
Dotted = original flux (FLUX+SKY); solid = sky-subtracted (FLUX).
Each file shown with its own colour.<br>
<span style="background:rgba(220,50,50,0.35); padding:2px 8px;">&#9632;</span>
Signal window (used for sky-line RMS)&nbsp;&nbsp;
<span style="background:rgba(30,160,30,0.30); padding:2px 8px;">&#9632;</span>
Background (mask-selected clean pixels, used for noise RMS).
</div>
'''

    # ── Write combined HTML ───────────────────────────────────
    if outroot == '':
        outroot = Path(filenames[0]).stem + '_eval'
    outfile = f'{outroot}.html'

    html1 = fig1.to_html(full_html=False, include_plotlyjs=True)
    html2 = fig2.to_html(full_html=False, include_plotlyjs=False)
    html3 = fig3.to_html(full_html=False, include_plotlyjs=False)
    html4 = fig4.to_html(full_html=False, include_plotlyjs=False)

    with open(outfile, 'w') as fh:
        fh.write('<!DOCTYPE html>\n<html>\n<body>\n')
        fh.write(html1)
        fh.write('\n')
        fh.write(html2)
        fh.write('\n')
        fh.write(_fig4_html_header)
        fh.write(html4)
        fh.write('\n')
        fh.write(_stats_html)
        fh.write(html3)
        fh.write('\n</body>\n</html>\n')

    print(f'Wrote {outfile}')


# ──────────────────────────────────────────────────────────────
# Command-line entry point
# ──────────────────────────────────────────────────────────────

def steer(argv):
    filenames       = []
    positional_nums = []
    n_sample        = 20
    outroot         = ''

    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg in ('-h', '--help'):
            print(_USAGE)
            return
        elif arg == '-num':
            i += 1
            n_sample = int(argv[i])
        elif arg == '-out':
            i += 1
            outroot = argv[i]
        elif arg.startswith('-'):
            print('Error: unknown option "%s"' % arg)
            print(_USAGE)
            return
        else:
            try:
                positional_nums.append(float(arg))
            except ValueError:
                filenames.append(arg)
        i += 1

    if not filenames:
        print(_USAGE)
        return

    wmin = positional_nums[0] if len(positional_nums) > 0 else 3600.0
    wmax = positional_nums[1] if len(positional_nums) > 1 else 9800.0

    if outroot:
        # Explicit output root: overlay all files in one HTML
        plot_eval(filenames, wmin=wmin, wmax=wmax, n_sample=n_sample, outroot=outroot)
    else:
        # Default: one HTML per file
        for filename in filenames:
            plot_eval(filename, wmin=wmin, wmax=wmax, n_sample=n_sample, outroot='')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
