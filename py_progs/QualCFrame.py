#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Create an html file, with various plots, which can be used
as a tool to assess the quality of the lvmdrp flux calibration
and sky measurements in an lvmCFrame exposure.


Command line usage (if any):

    usage: QualCFrame.py [-h] CFrame1 CFrame2 ...

    where
        -h prints this documentation and exits
        CFrame1 CFrame2 ... are files to be analyzed.


Description:

    This routine reads an lvmCFrame file and constructs an html
    file that summarizes:

        - how the SCI/STD/MOD flux calibration methods compare to
          each other (when more than one is available), read
          directly from the FLUXCAL_STD/FLUXCAL_SCI/FLUXCAL_MOD
          extensions already present in the CFrame -- no Gaia
          network access is needed for this comparison

        - how consistent the SkyE and SkyW per-telescope sky
          models (SKY_EAST/SKY_WEST) are with each other

    Unlike QuickLook.py (which works on the sky-subtracted
    lvmSFrame), this routine works on the lvmCFrame, before sky
    subtraction, so it can look directly at the two independent
    sky models and at all three flux calibration methods' own
    sensitivity curves, none of which survive into the lvmSFrame.

Primary routines:

    make_html is the primary driving routine
    steer handles the inputs

Notes:

    The html file is created in the current working directory
    and the various plots are in a subdirectory figs_qual_cf
    (kept separate from QuickLook.py's figs_qual, so both tools
    can be run in the same directory on a matching CFrame/SFrame
    pair without clobbering each other's images).

    Shared helpers (header access, angular distance, moon/sun
    info, percentile-based y-scaling, fiber selection) are
    imported from QuickLook.py rather than duplicated, so the two
    tools don't drift apart on shared logic.

History::

    260905 ksl Coding begun

'''

import os
import re
import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import xhtml
from lvm_ksl import QuickLook
from lvm_ksl import eval_standard
from lvmdrp.core.fluxcal import GaiaXPSpectra


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


SENS_BANDS = ('B', 'R', 'Z')
SENS_METHODS = ('STD', 'SCI', 'MOD')
SENS_COLORS = {'STD': 'tab:blue', 'SCI': 'tab:orange', 'MOD': 'tab:green'}
SENS_DISAGREE_WARN = 0.2  # fractional spread across methods that triggers a WARN note


def get_fluxcal_curve(hdul, ext_name):
    '''
    Read the mean/rms sensitivity curve from a FLUXCAL_STD/FLUXCAL_SCI/
    FLUXCAL_MOD extension of an lvmCFrame.

    Returns (wave, mean, rms, valid) where valid is False if the
    extension is missing or the mean curve is entirely non-finite/zero
    -- i.e. that method produced no usable calibration for this exposure.
    '''
    try:
        wave = hdul['WAVE'].data
        table = hdul[ext_name].data
        mean = np.asarray(table['mean'], dtype=float)
        rms = np.asarray(table['rms'], dtype=float)
    except KeyError:
        return None, None, None, False

    finite = np.isfinite(mean) & (mean != 0)
    valid = np.sum(finite) > 0.5 * len(mean)
    return wave, mean, rms, valid


def sensitivity_summary_table(hdr):
    '''
    Build a small html table (as rows for xhtml.table) comparing the
    band-averaged STD/SCI/MOD sensitivities from the *SENM{band}
    header keywords, flagging the method actually applied (FLUXCAL
    header) and any missing (-999.9 sentinel) or wildly discrepant
    values.
    '''
    method = QuickLook.get_header_string(hdr, 'FLUXCAL', 'Unknown')
    rows = [['Band', 'STD', 'SCI', 'MOD', 'Note']]

    for band in SENS_BANDS:
        vals = {name: QuickLook.get_header_value(hdr, '%sSENM%s' % (name, band))
                for name in SENS_METHODS}
        ok = {name: (vals[name] is not None and vals[name] > -900 and vals[name] > 0)
              for name in SENS_METHODS}
        good_vals = [vals[name] for name in SENS_METHODS if ok[name]]

        note = ''
        if len(good_vals) >= 2:
            spread = (max(good_vals) - min(good_vals)) / np.median(good_vals)
            if spread > SENS_DISAGREE_WARN:
                note = 'WARN: methods disagree by %.0f%%' % (spread * 100)
        note = (note + (', ' if note else '') + 'chosen: %s' % method)

        def fmt(name):
            if not ok[name]:
                return 'FAILED'
            marker = ' *' if name == method else ''
            return '%.3e%s' % (vals[name], marker)

        rows.append([band, fmt('STD'), fmt('SCI'), fmt('MOD'), note])

    return rows


fluxcal_comment = '''
Comparison of the three possible flux-calibration methods (STD=dedicated standard-star fibers,
SCI=Gaia-matched field stars in the science IFU, MOD=stellar atmosphere models fit to the standard
stars), read directly from the FLUXCAL_STD/FLUXCAL_SCI/FLUXCAL_MOD extensions already present in this
lvmCFrame -- these are computed independently for all three methods regardless of which one ends up
applied, so no external network access is needed for this comparison. The top panel overlays whichever
methods produced usable sensitivity curves for this exposure; the thicker line marks the method
actually applied to the delivered FLUX (see the FLUXCAL header, and the table above). The bottom panel
shows the ratio of each available method to MOD (or to whichever pair is available if MOD failed), to
reveal wavelength-dependent disagreement rather than just an overall offset. A '*' in the table above
marks the applied method; FAILED marks a method with no usable data for this exposure/band.
'''


def eval_sensitivity_comparison(filename, outroot=''):
    '''
    Compare the SCI, STD, and MOD flux-calibration sensitivity curves
    stored in the lvmCFrame's FLUXCAL_STD/FLUXCAL_SCI/FLUXCAL_MOD
    extensions.

    Returns (figname, note): figname is None (with an explanatory
    note) only if *no* method has usable data at all. With just one
    valid method, the sensitivity curve is still plotted (no ratio
    panel, since there's nothing to compare it to).
    '''
    try:
        x = fits.open(filename)
    except Exception as e:
        return None, 'Could not open %s (%s)' % (filename, e)

    hdr = x['PRIMARY'].header
    method = QuickLook.get_header_string(hdr, 'FLUXCAL', 'Unknown')

    curves = {}
    for name in SENS_METHODS:
        wave, mean, rms, valid = get_fluxcal_curve(x, 'FLUXCAL_%s' % name)
        if valid:
            curves[name] = (wave, mean, rms)

    if len(curves) == 0:
        return None, 'No flux-cal method has usable data for this exposure'

    have_comparison = len(curves) >= 2
    if have_comparison:
        fig = plt.figure(1, (9, 9))
        plt.clf()
        gs = GridSpec(2, 1, figure=fig, height_ratios=[2, 1])
        ax1 = fig.add_subplot(gs[0])
    else:
        fig = plt.figure(1, (9, 6))
        plt.clf()
        ax1 = fig.add_subplot(1, 1, 1)

    for name, (wave, mean, rms) in curves.items():
        lw = 2.5 if name == method else 1.2
        label = '%s%s' % (name, ' (applied)' if name == method else '')
        ax1.semilogy(wave, mean, label=label, color=SENS_COLORS[name], lw=lw)
    ax1.set_xlim(3600, 9600)
    ax1.set_ylabel('Sensitivity [erg / (ct cm2)]')
    ax1.legend()
    ax1.set_title('Flux calibration comparison, %s' % os.path.basename(filename))
    if not have_comparison:
        ax1.set_xlabel('Wavelength [Angstrom]')

    note = ''
    if have_comparison:
        ax2 = fig.add_subplot(gs[1], sharex=ax1)
        ratios = []
        if 'MOD' in curves:
            _, mean_mod, _ = curves['MOD']
            for name in ('SCI', 'STD'):
                if name in curves:
                    wave, mean, _ = curves[name]
                    ratio = mean / mean_mod
                    ax2.plot(wave, ratio, label='%s / MOD' % name, color=SENS_COLORS[name])
                    ratios.append(ratio)
        if not ratios:
            names = list(curves.keys())
            wave_a, mean_a, _ = curves[names[0]]
            _, mean_b, _ = curves[names[1]]
            ratio = mean_a / mean_b
            ax2.plot(wave_a, ratio, label='%s / %s' % (names[0], names[1]), color='k')
            ratios.append(ratio)
        ax2.axhline(1.0, ls=':', color='0.4')
        # auto-scale around 1.0, never zooming in tighter than +/-0.5 (the normal-
        # agreement case) but widening for a real large disagreement (like SCI vs
        # MOD/STD above) instead of silently clipping it off-screen
        all_ratios = np.concatenate(ratios)
        lo, hi = np.nanpercentile(all_ratios, [1, 99])
        ax2.set_ylim(min(lo, 0.5), max(hi, 1.5))
        ax2.set_xlim(3600, 9600)
        ax2.set_xlabel('Wavelength [Angstrom]')
        ax2.set_ylabel('Ratio')
        ax2.legend()
    else:
        note = ('Only %s has usable data for this exposure -- no comparison possible.'
                % list(curves.keys())[0])

    plt.tight_layout()

    location = './figs_qual_cf/'
    if not os.path.isdir(location):
        os.mkdir(location)
    if outroot == '':
        outroot = os.path.basename(filename).replace('.fits', '')

    figname = '%s/%s_fluxcal.png' % (location, outroot)
    plt.savefig(figname)
    plt.close(fig)

    return figname, note


def get_sci_calibration_fibers(hdr, max_n=15):
    '''
    Retrieve (index, fiberid, gaia_id) for the SCI method's Gaia-
    matched field stars, from the SCI{i}FIB/SCI{i}ID header keywords
    written by science_sensitivity(). SCI{i}FIB is already a raw
    fiberid (int), unlike STD{n}FIB (see get_std_calibration_fibers).
    Not every slot 1..max_n is populated -- a star that failed
    acquisition/matching leaves a gap.
    '''
    out = []
    for i in range(1, max_n + 1):
        fib = hdr.get('SCI%dFIB' % i)
        if fib is None:
            continue
        gaia_id = hdr.get('SCI%dID' % i)
        out.append((i, int(fib), gaia_id))
    return out


def get_std_calibration_fibers(hdr, xtab, max_n=15):
    '''
    Retrieve (index, fiberid, gaia_id) for the STD/MOD methods'
    standard stars. STD{n}FIB is an orig_ifulabel string (e.g.
    "P1-2"), not a raw fiberid, so it has to be matched against the
    SLITMAP to get the numeric fiberid used to index FLUX. Not every
    slot 1..max_n is populated (ACQ=False, or excluded for lacking a
    Gaia XP spectrum, leaves a gap).
    '''
    out = []
    for n in range(1, max_n + 1):
        label = hdr.get('STD%dFIB' % n)
        if label is None or str(label) == 'None':
            continue
        match = xtab[xtab['orig_ifulabel'] == label]
        if len(match) == 0:
            continue
        gaia_id = hdr.get('STD%dID' % n)
        out.append((n, int(match['fiberid'][0]), gaia_id))
    return out


calib_spectra_comment = '''
Approximately sky-subtracted FLUX spectra of the individual stars used by the SCI and STD/MOD flux
calibration methods, read directly from the lvmCFrame's FLUX extension at each star's own fiber, in
color, lightly boxcar-smoothed for display (same convention as eval_standard.py) so residual sky
noise doesn't obscure the continuum shape. The sky estimate subtracted is the median of the ~1800
Sci-telescope science fibers' own FLUX in this same CFrame -- when SKYSRC=SCIMED (the current
default, shown in the overview above), this is the same quantity (and the same broadcast-to-every-
fiber design) apply_fluxcal's own sensitivity derivation actually subtracted internally, so this
should closely match what the DRP itself compared to Gaia; if SKYSRC is not SCIMED, this is only an
approximation and may not match. Top panel: Gaia-matched field stars used by the SCI method
(SCI#FIB headers). Bottom panel: the dedicated standard-star fibers used by both the STD and MOD
methods (STD#FIB headers, which is why one panel serves both). Overplotted in solid black is each
used star's own Gaia XP spectrum, read from the same $LVM_MASTER_DIR/gaia_cache directory the DRP's
own flux calibration uses (so no network access is needed for a star lvmdrp already downloaded
during this exposure's reduction). A star can be acquired (ACQ=True) but still be excluded by the
pipeline itself if its FLUXCAL_STD/FLUXCAL_MOD/FLUXCAL_SCI column is all-NaN (e.g. a low-signal
cut) -- those are drawn dashed and grey, labeled "[excluded]", get no Gaia overlay, and are left out
of the axis auto-scaling so one bad/noisy fiber can't compress the scale for every other,
genuinely-used star. A dotted line labeled "[STD only]"/"[MOD only]" means the star was used by one
of those two methods but dropped by the other.
'''


def _overlay_gaia_spectrum(ax, gaia, gaia_id):
    '''
    Overplot a single star's cached/downloaded Gaia XP spectrum in solid
    black, high-contrast against any observed spectrum's color/noise.
    gaia_id is assumed not None. Returns True on success, False if the
    spectrum could not be retrieved.
    '''
    try:
        gaia.fetch_xp_spectra([gaia_id])
        gwave, gflux = gaia.load_xp_spectra(gaia_id)
        ax.semilogy(gwave, gflux[0], color='k', lw=1.3, alpha=0.85, zorder=10)
        return True
    except Exception as e:
        print('Error: could not retrieve Gaia spectrum for %s: %s' % (gaia_id, e))
        return False


def eval_calibration_spectra(filename, outroot=''):
    '''
    Plot the individual raw FLUX spectra of the stars used by the SCI
    method (field stars) and the STD/MOD methods (dedicated standard-
    star fibers), directly from the lvmCFrame, overlaid with each
    star's own Gaia XP spectrum.

    Returns (figname, note): figname is None (with an explanatory
    note) if no calibration-star header keywords are found at all.
    note reports how many Gaia spectra (of those with a Gaia ID) could
    not be retrieved, if any.
    '''
    try:
        x = fits.open(filename)
    except Exception as e:
        return None, 'Could not open %s (%s)' % (filename, e)

    hdr = x['PRIMARY'].header
    xtab = Table(x['SLITMAP'].data)

    sci_stars = get_sci_calibration_fibers(hdr)
    std_stars = get_std_calibration_fibers(hdr, xtab)

    if len(sci_stars) == 0 and len(std_stars) == 0:
        return None, 'No SCI or STD/MOD calibration-star header keywords were found'

    gaia = GaiaXPSpectra(cache_dir=eval_standard.get_gaia_cache_dir())

    wav = x['WAVE'].data
    flux = x['FLUX'].data
    mask = x['MASK'].data.astype(bool)

    # a star's column in these tables is all-NaN if the pipeline itself
    # excluded it (e.g. a low-signal cut in standard_sensitivity/
    # science_sensitivity), even when it was successfully acquired
    # (ACQ=True) -- flag those distinctly rather than plotting them as if
    # they were normal, on-scale data
    sci_sen = x['FLUXCAL_SCI'].data if 'FLUXCAL_SCI' in x else None
    std_sen = x['FLUXCAL_STD'].data if 'FLUXCAL_STD' in x else None
    mod_sen = x['FLUXCAL_MOD'].data if 'FLUXCAL_MOD' in x else None

    def _col_valid(table, colname):
        if table is None or colname not in table.columns.names:
            return None
        return np.isfinite(np.asarray(table[colname])).any()

    def _spectrum(fiberid):
        return np.ma.masked_array(flux[fiberid - 1], mask[fiberid - 1])

    # Approximate the flux-calibration sky and subtract it, for a cleaner
    # comparison to each star's Gaia spectrum. When SKYSRC=SCIMED (the
    # current default for essentially every exposure), the actual
    # flux-cal sky *is* the median of the Sci-telescope science fibers
    # (see combine_skies) -- computed here directly from the CFrame's own
    # FLUX, at the same fiber selection. Since apply_fluxcal puts every
    # fiber (Sci and STD alike) on the same physical flux-density
    # footing, this median is a valid sky estimate to subtract from any
    # star's spectrum regardless of its own exposure time. This is only
    # an approximation, though: it's computed post-rectification on the
    # final CFrame, not pre-rectification on raw counts the way
    # combine_skies actually did it, and if SKYSRC is not SCIMED (an
    # explicit sky_weights override) this isn't the quantity that was
    # actually used at all.
    sky_src = QuickLook.get_header_string(hdr, 'SKYSRC', 'Unknown')
    sci_fibers_tab = QuickLook.scifib(xtab, select='science', telescope='Sci')
    sci_sky_mask = mask[sci_fibers_tab['fiberid'] - 1]
    sci_sky_flux = np.ma.masked_array(flux[sci_fibers_tab['fiberid'] - 1], sci_sky_mask)
    approx_sky = np.ma.median(sci_sky_flux, axis=0)

    def _sky_subtracted(fiberid):
        return _spectrum(fiberid) - approx_sky

    def _smooth(spec):
        # boxcar-smooth for display only, same convention (and function)
        # eval_standard.py already uses for star spectra -- residual sky
        # noise otherwise obscures the continuum shape needed to compare
        # against Gaia by eye
        filled = np.ma.filled(spec, 0)
        return eval_standard.xsmooth(filled, smooth=11)

    def _set_percentile_ylim(ax, arrays):
        if not arrays:
            return
        combined = np.concatenate(arrays)
        positive = combined[np.isfinite(combined) & (combined > 0)]
        if len(positive) > 0:
            lo, hi = np.nanpercentile(positive, [0.5, 99.9])
            ax.set_ylim(lo, hi * 3)

    fig = plt.figure(3, (9, 9))
    plt.clf()
    gs = GridSpec(2, 1, figure=fig)

    ax1 = fig.add_subplot(gs[0])
    sci_specs_for_scale = []
    sci_tried = sci_failed = 0
    for i, fiberid, gaia_id in sci_stars:
        spec = _smooth(_sky_subtracted(fiberid))
        used = _col_valid(sci_sen, 'SCI%dSEN' % i)
        label = 'SCI%d (fiber %d)' % (i, fiberid)
        if used:
            sci_specs_for_scale.append(spec)
            ax1.semilogy(wav, spec, lw=1.0, label=label)
            if gaia_id is not None:
                sci_tried += 1
                if not _overlay_gaia_spectrum(ax1, gaia, gaia_id):
                    sci_failed += 1
        else:
            label += ' [excluded]'
            ax1.semilogy(wav, spec, lw=0.8, ls='--', alpha=0.5, color='0.5', label=label)
    ax1.set_xlim(3600, 9600)
    _set_percentile_ylim(ax1, sci_specs_for_scale)
    ax1.set_ylabel('FLUX (sky-sub., approx., smoothed)')
    ax1.set_title('SCI method field-star spectra, %s' % os.path.basename(filename))
    if sci_stars:
        ax1.legend(fontsize=7, ncol=3)
    else:
        ax1.text(0.5, 0.5, 'No SCI stars available for this exposure', transform=ax1.transAxes, ha='center')

    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    std_specs_for_scale = []
    std_tried = std_failed = 0
    for n, fiberid, gaia_id in std_stars:
        spec = _smooth(_sky_subtracted(fiberid))
        used_std = _col_valid(std_sen, 'STD%dSEN' % n)
        used_mod = _col_valid(mod_sen, 'STD%dSEN' % n)
        label = 'STD%d (fiber %d)' % (n, fiberid)
        if used_std and used_mod:
            std_specs_for_scale.append(spec)
            ax2.semilogy(wav, spec, lw=1.0, label=label)
            if gaia_id is not None:
                std_tried += 1
                if not _overlay_gaia_spectrum(ax2, gaia, gaia_id):
                    std_failed += 1
        elif not used_std and not used_mod:
            label += ' [excluded]'
            ax2.semilogy(wav, spec, lw=0.8, ls='--', alpha=0.5, color='0.5', label=label)
        else:
            # used by one method but not the other -- still real data, keep
            # it in the axis scaling, just flag which method dropped it
            std_specs_for_scale.append(spec)
            label += (' [MOD only]' if used_mod else ' [STD only]')
            ax2.semilogy(wav, spec, lw=1.0, ls=':', label=label)
            if gaia_id is not None:
                std_tried += 1
                if not _overlay_gaia_spectrum(ax2, gaia, gaia_id):
                    std_failed += 1
    ax2.set_xlim(3600, 9600)
    _set_percentile_ylim(ax2, std_specs_for_scale)
    ax2.set_xlabel('Wavelength [Angstrom]')
    ax2.set_ylabel('FLUX (sky-sub., approx., smoothed)')
    ax2.set_title('STD/MOD method standard-star spectra')
    if std_stars:
        ax2.legend(fontsize=7, ncol=3)
    else:
        ax2.text(0.5, 0.5, 'No STD/MOD stars available for this exposure', transform=ax2.transAxes, ha='center')

    plt.tight_layout()

    location = './figs_qual_cf/'
    if not os.path.isdir(location):
        os.mkdir(location)
    if outroot == '':
        outroot = os.path.basename(filename).replace('.fits', '')

    figname = '%s/%s_calibspec.png' % (location, outroot)
    plt.savefig(figname)
    plt.close(fig)

    notes = []
    if sky_src != 'SCIMED':
        notes.append('SKYSRC=%s (not SCIMED) -- the sky subtracted here (Sci-fiber median) is only '
                     'an approximation and may not match what the DRP actually used internally.' % sky_src)
    total_failed = sci_failed + std_failed
    if total_failed > 0:
        notes.append('%d of %d Gaia spectra could not be retrieved (no network access, and nothing '
                     'cached locally for those stars).' % (total_failed, sci_tried + std_tried))

    return figname, ' '.join(notes)


sky_comment = '''
Comparison of the SkyE and SkyW per-telescope sky models (SKY_EAST/SKY_WEST), before sky subtraction.
Unlike the equivalent panel in QuickLook.py (which compares sky-subtraction residuals in the final
lvmSFrame), this compares the two telescopes' own raw sky estimates directly, evaluated at the
science-telescope fibers. The middle panel shows the difference in total sky flux between the two
telescopes (nearer minus further, matching QuickLook's convention); the bottom row shows the same
difference in three diagnostic wavelength windows, corresponding to Hbeta-[OII], Halpha-[SII], and
[SIII]9071.
'''


def eval_sky_comparison(filename, outroot=''):
    '''
    Compare the SkyE and SkyW per-telescope sky models (SKY_EAST/
    SKY_WEST) stored in an lvmCFrame.

    Returns (figname, note): note is a warning string (possibly empty)
    explaining a caveat about the data, not a failure -- figname is
    None only if the file couldn't be read.
    '''
    try:
        x = fits.open(filename)
    except Exception as e:
        return None, 'Could not open %s (%s)' % (filename, e)

    hdr = x['PRIMARY'].header
    ra = QuickLook.get_header_value(hdr, 'SCIRA')
    dec = QuickLook.get_header_value(hdr, 'SCIDEC')
    ra_sky_e = QuickLook.get_header_value(hdr, 'SKYERA')
    dec_sky_e = QuickLook.get_header_value(hdr, 'SKYEDEC')
    ra_sky_w = QuickLook.get_header_value(hdr, 'SKYWRA')
    dec_sky_w = QuickLook.get_header_value(hdr, 'SKYWDEC')

    distance_sky_e = QuickLook.distance(ra, dec, ra_sky_e, dec_sky_e)
    distance_sky_w = QuickLook.distance(ra, dec, ra_sky_w, dec_sky_w)

    sky_ew = QuickLook.get_header_value(hdr, 'SKYEW')
    sky_ww = QuickLook.get_header_value(hdr, 'SKYWW')

    xtab = Table(x['SLITMAP'].data)
    sci_fibers = QuickLook.scifib(xtab, select='science', telescope='Sci')
    skye_fibers = QuickLook.scifib(xtab, select='SKY', telescope='SkyE')
    skyw_fibers = QuickLook.scifib(xtab, select='SKY', telescope='SkyW')

    wav = x['WAVE'].data
    mask = x['MASK'].data[sci_fibers['fiberid'] - 1].astype(bool)

    sky_e = np.ma.masked_array(x['SKY_EAST'].data[sci_fibers['fiberid'] - 1], mask)
    sky_w = np.ma.masked_array(x['SKY_WEST'].data[sci_fibers['fiberid'] - 1], mask)

    sky_e_med = np.ma.median(sky_e, axis=0)
    sky_w_med = np.ma.median(sky_w, axis=0)

    # Raw FLUX at the actual dedicated SkyE/SkyW telescope fibers -- what
    # those telescopes really observed, independent of whatever combine_skies
    # later did with it. This is the only place a per-telescope problem (e.g.
    # one sky telescope sitting near the Moon) is visible in a file where
    # SKY_EAST/SKY_WEST don't yet reflect genuine independent per-telescope
    # models (e.g. a pre-Option-A SCIMED broadcast, where SKY_EAST/SKY_WEST
    # are identical to each other and disconnected from either telescope's
    # own data) -- the model comparison above would show no disagreement at
    # all in that case even if one telescope's raw sky is badly contaminated.
    skye_mask = x['MASK'].data[skye_fibers['fiberid'] - 1].astype(bool)
    skyw_mask = x['MASK'].data[skyw_fibers['fiberid'] - 1].astype(bool)
    skye_flux = np.ma.masked_array(x['FLUX'].data[skye_fibers['fiberid'] - 1], skye_mask)
    skyw_flux = np.ma.masked_array(x['FLUX'].data[skyw_fibers['fiberid'] - 1], skyw_mask)
    skye_flux_med = np.ma.median(skye_flux, axis=0)
    skyw_flux_med = np.ma.median(skyw_flux, axis=0)

    fig = plt.figure(2, (8, 12))
    plt.clf()
    gs = GridSpec(3, 3, figure=fig)

    ax1 = fig.add_subplot(gs[0, :])
    ax1.semilogy(wav, sky_e_med, label='SKY_EAST (model)', zorder=1)
    ax1.semilogy(wav, sky_w_med, label='SKY_WEST (model)', zorder=0)
    ax1.semilogy(wav, skye_flux_med, label='SkyE fiber FLUX (raw)', zorder=1, ls='--')
    ax1.semilogy(wav, skyw_flux_med, label='SkyW fiber FLUX (raw)', zorder=0, ls='--')
    ax1.set_xlim(3600, 9600)
    # a log-scale autoscale gets distorted by near-zero/negative noise dips;
    # use the continuum-dominated 1st-99.9th percentile of the positive
    # values instead, so genuine emission-line peaks stay visible without a
    # handful of near-zero pixels blowing the range out to 1e-18
    combined = np.ma.filled(np.ma.concatenate([sky_e_med, sky_w_med, skye_flux_med, skyw_flux_med]), np.nan)
    positive = combined[combined > 0]
    if len(positive) > 0:
        ax1.set_ylim(*np.nanpercentile(positive, [1, 99.9]))
    ax1.legend()
    ax1.set_title('SkyE/SkyW sky model vs. raw fiber flux (pre-sky-subtraction)')

    ax2 = fig.add_subplot(gs[1, :])
    delta = sky_w_med - sky_e_med
    if distance_sky_w < distance_sky_e:
        ax2.plot(wav, delta, label='SkyW-SkyE (Nearer-Further)')
    else:
        delta = -delta
        ax2.plot(wav, delta, label='SkyE-SkyW (Nearer-Further)')
    ymin, ymax = QuickLook.get_percentile_yscale(delta, 1, 99.9, min_half_range=2 * QuickLook.MW_5SIGMA)
    ax2.set_ylim(ymin, ymax)
    ax2.set_xlim(3600, 9600)
    ax2.legend()

    windows = [(4800, 5100), (6500, 6800), (9480, 9586)]
    for i, (wmin, wmax) in enumerate(windows):
        ax = fig.add_subplot(gs[2, i])
        xwav, xdelta = QuickLook.limit_spectrum(wav, delta, wmin, wmax)
        xdelta = xdelta - np.ma.median(xdelta)
        ax.plot(xwav, xdelta)
        ax.plot([wmin, wmax], [QuickLook.MW_5SIGMA] * 2, ':r')
        ax.plot([wmin, wmax], [-QuickLook.MW_5SIGMA] * 2, ':r')
        ax.set_xlim(wmin, wmax)
        ymin, ymax = QuickLook.get_percentile_yscale(xdelta, 1, 99, min_half_range=2 * QuickLook.MW_5SIGMA)
        ax.set_ylim(ymin, ymax)

    plt.tight_layout()

    location = './figs_qual_cf/'
    if not os.path.isdir(location):
        os.mkdir(location)
    if outroot == '':
        outroot = os.path.basename(filename).replace('.fits', '')

    figname = '%s/%s_sky.png' % (location, outroot)
    plt.savefig(figname)
    plt.close(fig)

    notes = []
    if sky_ew is not None and sky_ww is not None and abs(sky_ew - 1.0) < 1e-6 and abs(sky_ww) < 1e-6:
        notes.append('SKYEW=1.0, SKYWW=0.0 exactly -- SKY_EAST/SKY_WEST may be the flux-cal SCIMED '
                     'broadcast rather than genuine independent per-telescope sky models. Check that '
                     'this exposure was reduced with the "Option A" combine_skies fix.')

    # flag a badly contaminated sky telescope (e.g. one sitting near the Moon)
    # directly from the raw fiber flux, independent of whatever the derived
    # SKY_EAST/SKY_WEST models show -- this is the check that would have
    # caught issue #250's original motivating case (SkyW 4.5 deg from a
    # near-full Moon) even on a file where the models above show no
    # disagreement at all
    e_level = np.ma.median(skye_flux_med)
    w_level = np.ma.median(skyw_flux_med)
    if e_level is not np.ma.masked and w_level is not np.ma.masked and min(e_level, w_level) > 0:
        raw_ratio = max(e_level, w_level) / min(e_level, w_level)
        if raw_ratio > 3:
            brighter = 'SkyE' if e_level > w_level else 'SkyW'
            notes.append('WARN: raw fiber flux in %s is %.1fx brighter than the other sky telescope -- '
                         'check %s_MOON_SEP/MOON_FLI in the header for possible Moon contamination.'
                         % (brighter, raw_ratio, brighter.upper()))

    return figname, ' '.join(notes)


DIAGNOSTIC_LINES = [
    # (label, center wavelength, y-scale sub-window or None to use the full plotted range)
    ('[OII]3727', 3727.0, None),
    ('Hbeta4861', 4861.0, None),
    ('[OIII]4959,5007', (4959.0 + 5007.0) / 2, None),
    ('Halpha6563', 6563.0, None),
    ('[SII]6717,6731', (6717.0 + 6731.0) / 2, None),
    # this window is dominated by strong OH airglow (see the sharp peaks
    # outside 9525-9540 in the plot) that would otherwise blow out the
    # y-scale -- base it on just the actual [SIII] line region instead,
    # while still showing the full +/-50A window on the x-axis
    ('[SIII]9533', 9533.1, (9525.0, 9540.0)),
]
LINE_WINDOW_HALF_WIDTH = 50.0  # Angstrom, +/- around each line center -- wide enough to show line + local continuum

field_vs_sky_comment = '''
For each diagnostic emission-line window ([OII], Hbeta, the [OIII] doublet, Halpha, and the [SII]
doublet), real spectra as a function of wavelength: the grey band spans the 10th-90th percentile of
the observed field brightness at each wavelength pixel (computed across all ~1800 Sci-telescope
science fibers), with the field median as a black line; SKY_EAST and SKY_WEST are overlaid as their
own spectra (median across Sci fibers at each pixel, same as the SkyE/SkyW comparison above, just
zoomed to these lines), labeled "(near)"/"(far)" for whichever is angularly closer to the science
field this exposure. This asks a different question from that comparison: whether each sky estimate
is *plausible* relative to what the field actually looks like, not just whether the two telescopes
agree with each other. If Vela's nebular emission fills the whole science IFU, a sky spectrum sitting
below even the field's 10th percentile at Halpha/[OIII] isn't necessarily wrong -- it can mean there's
no truly blank sky within the field itself, which is exactly why combine_skies moved to a Sci-fiber-
median sky in the first place. [OII]3727 sits at LVM's blue edge, where S/N is lowest, so expect that
panel to be noisier than the others. The [OIII] window is centered between 4959/5007 and the [SII]
window between 6717/6731 so both doublet components are visible.
'''


def eval_field_vs_sky_lines(filename, outroot=''):
    '''
    For a few diagnostic emission-line windows, plot the field's
    per-wavelength brightness distribution (10th/50th/90th percentile
    across all Sci-telescope science fibers) alongside the SKY_EAST/
    SKY_WEST sky model spectra -- a plausibility check on the sky
    estimate itself, rather than just a SkyE-vs-SkyW agreement check.

    Returns (figname, note).
    '''
    try:
        x = fits.open(filename)
    except Exception as e:
        return None, 'Could not open %s (%s)' % (filename, e)

    hdr = x['PRIMARY'].header
    ra = QuickLook.get_header_value(hdr, 'SCIRA')
    dec = QuickLook.get_header_value(hdr, 'SCIDEC')
    ra_sky_e = QuickLook.get_header_value(hdr, 'SKYERA')
    dec_sky_e = QuickLook.get_header_value(hdr, 'SKYEDEC')
    ra_sky_w = QuickLook.get_header_value(hdr, 'SKYWRA')
    dec_sky_w = QuickLook.get_header_value(hdr, 'SKYWDEC')
    distance_sky_e = QuickLook.distance(ra, dec, ra_sky_e, dec_sky_e)
    distance_sky_w = QuickLook.distance(ra, dec, ra_sky_w, dec_sky_w)
    e_tag = 'near' if distance_sky_e < distance_sky_w else 'far'
    w_tag = 'far' if distance_sky_e < distance_sky_w else 'near'

    xtab = Table(x['SLITMAP'].data)
    sci_fibers = QuickLook.scifib(xtab, select='science', telescope='Sci')

    wav = x['WAVE'].data
    fmask = x['MASK'].data[sci_fibers['fiberid'] - 1].astype(bool)
    flux = np.ma.masked_array(x['FLUX'].data[sci_fibers['fiberid'] - 1], fmask)
    sky_e = np.ma.masked_array(x['SKY_EAST'].data[sci_fibers['fiberid'] - 1], fmask)
    sky_w = np.ma.masked_array(x['SKY_WEST'].data[sci_fibers['fiberid'] - 1], fmask)

    nrows, ncols = 2, 3
    fig, axs = plt.subplots(nrows, ncols, figsize=(16, 9))
    axs = axs.flatten()
    for extra_ax in axs[len(DIAGNOSTIC_LINES):]:
        extra_ax.set_visible(False)

    for i, (ax, (name, wl, yscale_window)) in enumerate(zip(axs, DIAGNOSTIC_LINES)):
        wmin, wmax = wl - LINE_WINDOW_HALF_WIDTH, wl + LINE_WINDOW_HALF_WIDTH
        idx = (wav > wmin) & (wav < wmax)
        if idx.sum() == 0:
            ax.set_title('%s\n(out of range)' % name)
            continue

        xwav = wav[idx]
        filled = np.ma.filled(flux[:, idx], np.nan)
        p10 = np.nanpercentile(filled, 10, axis=0)
        p50 = np.nanpercentile(filled, 50, axis=0)
        p90 = np.nanpercentile(filled, 90, axis=0)

        sky_e_spec = np.ma.filled(np.ma.median(sky_e[:, idx], axis=0), np.nan)
        sky_w_spec = np.ma.filled(np.ma.median(sky_w[:, idx], axis=0), np.nan)

        ax.fill_between(xwav, p10, p90, color='0.85', label='field 10-90%ile' if i == 0 else None)
        ax.plot(xwav, p50, color='k', lw=1.5, label='field median' if i == 0 else None)
        ax.plot(xwav, sky_e_spec, color='tab:blue', lw=1.2,
               label='SKY_EAST (%s)' % e_tag if i == 0 else None)
        ax.plot(xwav, sky_w_spec, color='tab:orange', lw=1.2,
               label='SKY_WEST (%s)' % w_tag if i == 0 else None)

        ax.set_xlim(wmin, wmax)
        if yscale_window is not None:
            y_wmin, y_wmax = yscale_window
            yidx = (xwav >= y_wmin) & (xwav <= y_wmax)
            combined = np.concatenate([p10[yidx], p90[yidx], sky_e_spec[yidx], sky_w_spec[yidx]])
            combined = combined[np.isfinite(combined)]
            if len(combined) > 0:
                ax.set_ylim(min(0, np.nanmin(combined)), np.nanmax(combined) * 1.15)
        ax.set_title('%s (%.0f A)' % (name, wl))
        if i % ncols == 0:
            ax.set_ylabel('FLUX')
        if i >= len(DIAGNOSTIC_LINES) - ncols:
            ax.set_xlabel('Wavelength [Angstrom]')

    axs[0].legend(fontsize=8, loc='best')
    fig.suptitle('Field brightness (10-90%%ile) vs sky estimate, %s' % os.path.basename(filename))
    plt.tight_layout()

    location = './figs_qual_cf/'
    if not os.path.isdir(location):
        os.mkdir(location)
    if outroot == '':
        outroot = os.path.basename(filename).replace('.fits', '')

    figname = '%s/%s_fieldsky.png' % (location, outroot)
    plt.savefig(figname)
    plt.close(fig)

    return figname, ''


def create_overview(filename):
    '''
    Summarize header information about the CFrame exposure
    '''
    try:
        x = fits.open(filename)
    except Exception as e:
        print('Error: Could not open %s: %s' % (filename, e))
        return []

    hdr = x['PRIMARY'].header

    exposure = QuickLook.get_header_value(hdr, 'EXPOSURE')
    mjd = QuickLook.get_header_value(hdr, 'MJD')
    object_name = QuickLook.get_header_string(hdr, 'OBJECT')
    obs_time = QuickLook.get_header_string(hdr, 'OBSTIME')
    drp_version = QuickLook.get_header_string(hdr, 'DRPVER')
    fluxcal_method = QuickLook.get_header_string(hdr, 'FLUXCAL', 'Unknown')
    sky_src = QuickLook.get_header_string(hdr, 'SKYSRC', 'Unknown')

    ra = QuickLook.get_header_value(hdr, 'SCIRA')
    dec = QuickLook.get_header_value(hdr, 'SCIDEC')
    ra_sky_e = QuickLook.get_header_value(hdr, 'SKYERA')
    dec_sky_e = QuickLook.get_header_value(hdr, 'SKYEDEC')
    ra_sky_w = QuickLook.get_header_value(hdr, 'SKYWRA')
    dec_sky_w = QuickLook.get_header_value(hdr, 'SKYWDEC')

    distance_sky_e = QuickLook.distance(ra, dec, ra_sky_e, dec_sky_e)
    distance_sky_w = QuickLook.distance(ra, dec, ra_sky_w, dec_sky_w)

    moon_info = QuickLook.get_moon_info_las_campanas(obs_time)

    xlist = []
    xlist.append('Exposure : %d' % exposure)
    xlist.append('MJD      : %d' % mjd)
    xlist.append('Obs. time: %s' % obs_time)
    xlist.append('Object.  : %s' % object_name)
    xlist.append('DRP Version : %s' % drp_version)
    xlist.append('Flux-cal method applied : %s' % fluxcal_method)
    xlist.append('Sky source (flux-cal)   : %s' % sky_src)
    xlist.append('Science RA  Dec. : %8.2f %8.2f' % (ra, dec))
    xlist.append('SkyE    RA  Dec. (ang distance): %8.2f %8.2f (%8.2f)' %
                 (ra_sky_e, dec_sky_e, distance_sky_e))
    xlist.append('SkyW    RA  Dec. (ang distance): %8.2f %8.2f (%8.2f)' %
                 (ra_sky_w, dec_sky_w, distance_sky_w))
    xlist.append('Moon    RA  Dec. Alt.  Ill:  %8.2f %8.2f %8.2f %8.2f' %
                 (moon_info['MoonRA'], moon_info['MoonDec'], moon_info['MoonAlt'], moon_info['MoonIll']))
    xlist.append('Sun.    RA  Dec. Alt.:  %8.2f %8.2f %8.2f' %
                 (moon_info['SunRA'], moon_info['SunDec'], moon_info['SunAlt']))

    return xlist


def make_html(filename, outroot=''):
    '''
    Create an html file that summarizes the flux-calibration and
    sky-consistency quality of an lvmdrp-processed lvmCFrame
    '''

    if outroot == '':
        outroot = os.path.basename(filename).replace('.fits', '')

    string = xhtml.begin('LVMDRP CFrame Quality Assessment for %s' % filename)
    string += xhtml.hline()

    overview_list = create_overview(filename)
    string += xhtml.add_list(overview_list)

    string += xhtml.hline()
    string += xhtml.h2('Flux Calibration Comparison (STD / SCI / MOD)')

    hdr = fits.getheader(filename, 0)
    string += xhtml.table(sensitivity_summary_table(hdr))
    string += xhtml.paragraph(fluxcal_comment)

    figname, note = eval_sensitivity_comparison(filename, outroot)
    if figname:
        string += xhtml.image(figname, width=900, height=900)
        if note:
            string += xhtml.paragraph('Note: %s' % note)
    else:
        string += xhtml.paragraph('Could not compare flux-cal methods: %s' % note)

    string += xhtml.hline()
    string += xhtml.h2('Flux Calibration Input Spectra (SCI / STD-MOD)')
    string += xhtml.paragraph(calib_spectra_comment)

    calibspec_figname, calibspec_note = eval_calibration_spectra(filename, outroot)
    if calibspec_figname:
        string += xhtml.image(calibspec_figname, width=900, height=900)
        if calibspec_note:
            string += xhtml.paragraph('Note: %s' % calibspec_note)
    else:
        string += xhtml.paragraph('Could not plot calibration-star spectra: %s' % calibspec_note)

    string += xhtml.hline()
    string += xhtml.h2('SkyE / SkyW Consistency')
    string += xhtml.paragraph(sky_comment)

    sky_figname, sky_note = eval_sky_comparison(filename, outroot)
    if sky_figname:
        string += xhtml.image(sky_figname, width=900, height=1200)
        if sky_note:
            string += xhtml.paragraph('Note: %s' % sky_note)
    else:
        string += xhtml.paragraph('Could not compare SkyE/SkyW: %s' % sky_note)

    string += xhtml.hline()
    string += xhtml.h2('Field Brightness vs Sky Estimate (Line Diagnostics)')
    string += xhtml.paragraph(field_vs_sky_comment)

    fieldsky_figname, fieldsky_note = eval_field_vs_sky_lines(filename, outroot)
    if fieldsky_figname:
        string += xhtml.image(fieldsky_figname, width=1000, height=563)
        if fieldsky_note:
            string += xhtml.paragraph('Note: %s' % fieldsky_note)
    else:
        string += xhtml.paragraph('Could not plot field-vs-sky diagnostics: %s' % fieldsky_note)

    string += xhtml.hline()

    g = open(outroot + '.html', 'w')
    g.write(string)
    g.close()


def steer(argv):
    '''
    Just a steering routine
    '''

    i = 1
    files = []
    while i < len(argv):
        if argv[i][0:2] == '-h':
            print(_usage_from_doc(__doc__))
            return
        elif argv[i][0] == '-':
            print('Error: Unknown optional parameter; improperly formatted command line: ', argv)
            return
        else:
            files.append(argv[i])
        i += 1

    for one in files:
        make_html(one)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
