#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Plot how well the standards that are observed in an LVM exposure
are calibrated.


Command line usage (if any):

    usage: eval_standard.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History::

    240318 ksl Coding begun
    260727 ksl Rewrote compare_with_gaia() -- it was silently failing
        on every current SFrame file for two independent reasons: (1)
        it called ancillary_func.retrive_gaia_star(), which no longer
        exists in the installed lvmdrp (an AttributeError swallowed by
        a bare except); (2) it read STD#BIN/STD#ID/STD#FIB header
        keywords matched against an orig_ifulabel slit label in a
        dedicated 'standard'-targettype fiber set, none of which exist
        any more -- the DRP's current flux calibration
        (fluxCalMethod.py's science_sensitivity) instead identifies
        Gaia-matched field stars among the ordinary science-telescope
        fibers and writes SCI#ID/SCI#FIB (Gaia source id / raw
        fiberid) header keywords, up to 15 slots with gaps. get_standard()
        now indexes FLUX/WAVE directly by that fiberid. GAIA XP spectra
        are now fetched via lvmdrp.core.fluxcal.GaiaXPSpectra, cached
        under $LVM_MASTER_DIR/gaia_cache (the same directory the DRP's
        own flux calibration populates during reduction, so spectra it
        already downloaded are reused instead of re-queried). Retrieval
        is now per-star and failure-tolerant rather than all-or-
        nothing. compare_with_gaia()/qual_eval() now return
        (outfile_or_None, message), so a caller (QuickLook.py) can show
        *why* the comparison failed or partially failed instead of a
        generic could-not-do message.
    260907 ksl Updated a docstring reference for QuickLook.py's rename
        to QualSFrame.py -- no functional change.

'''



from astropy.io import fits
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from astropy.table import Table
from lvmdrp.core.fluxcal import GaiaXPSpectra



from scipy.signal.windows import boxcar
from scipy.signal import convolve


import re


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


def xsmooth(flux,smooth=21):
    '''
    boxcar smooth the flux
    '''
    if (smooth)>1:
        q=convolve(flux,boxcar(smooth)/float(smooth),mode='same')
        return(q)
    else:
        return(flux)



NSCI_MAX=15


def get_header_value(header, key, default_value=-999.0, verbose=False):
    '''
    Robust way to get a header value if it exists
    '''

    try:
        value = header[key]
        if value==None:
            value=default_value
        elif isinstance(value, str):
            try:
                value = float(value)  # or int(value) if it's an integer
            except ValueError as e:
                if verbose:
                    print(f"Failed to convert '{value}' to a number for key '{key}': {e}")
                value = default_value
    except KeyError as e:
        if verbose:
            print(f"Key '{key}' not found in header: {e}")
        value = default_value
    return value


def get_header_string(header, key, default_string='Unknown', verbose=False):
    '''
    Robust way to get a header value if it exists
    '''

    try:
        value = header[key]
        if value==None:
            value=default_string
        elif isinstance(value, str):
            if value=='':
                return default_string
            return value
        else:
            if verbose:
                print(f"Key '{key}' found, but not string")
            return default_string
    except KeyError as e:
        if verbose:
            print(f"Key '{key}' not found in header: {e}")
        value = default_string
    return value


SENS_BANDS = ('B', 'R', 'Z')
SENS_METHODS = ('STD', 'SCI', 'MOD')
SENS_COLORS = {'STD': 'tab:blue', 'SCI': 'tab:orange', 'MOD': 'tab:green'}
SENS_DISAGREE_WARN = 0.2  # fractional spread across methods that triggers a WARN note


def get_fluxcal_curve(hdul, ext_name):
    '''
    Read the mean/rms sensitivity curve from a FLUXCAL_STD/FLUXCAL_SCI/
    FLUXCAL_MOD extension of an lvmCFrame or lvmSFrame -- both carry
    the same tables, written once during flux calibration and passed
    through unchanged by quick_sky_subtraction.

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
    method = get_header_string(hdr, 'FLUXCAL', 'Unknown')
    rows = [['Band', 'STD', 'SCI', 'MOD', 'Note']]

    for band in SENS_BANDS:
        vals = {name: get_header_value(hdr, '%sSENM%s' % (name, band))
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
file -- these are computed independently for all three methods regardless of which one ends up
applied, so no external network access is needed for this comparison (unlike the Gaia comparison
below). The top panel overlays whichever methods produced usable sensitivity curves for this
exposure; the thicker line marks the method actually applied to the delivered FLUX (see the FLUXCAL
header, and the table above). The bottom panel shows the ratio of each available method to MOD (or
to whichever pair is available if MOD failed), to reveal wavelength-dependent disagreement rather
than just an overall offset. A '*' in the table above marks the applied method; FAILED marks a
method with no usable data for this exposure/band.
'''


def eval_sensitivity_comparison(filename, outroot='', fignum=1, outdir='./figs_qual/'):
    '''
    Compare the SCI, STD, and MOD flux-calibration sensitivity curves
    stored in the FLUXCAL_STD/FLUXCAL_SCI/FLUXCAL_MOD extensions.
    fignum/outdir let callers avoid a matplotlib figure-number clash
    with their own other plots and keep each tool's PNGs in its own
    directory (e.g. QualSFrame.py's figs_qual/ vs QualCFrame.py's
    figs_qual_cf/).

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
    method = get_header_string(hdr, 'FLUXCAL', 'Unknown')

    curves = {}
    for name in SENS_METHODS:
        wave, mean, rms, valid = get_fluxcal_curve(x, 'FLUXCAL_%s' % name)
        if valid:
            curves[name] = (wave, mean, rms)

    if len(curves) == 0:
        return None, 'No flux-cal method has usable data for this exposure'

    have_comparison = len(curves) >= 2
    if have_comparison:
        fig = plt.figure(fignum, (9, 9))
        plt.clf()
        gs = GridSpec(2, 1, figure=fig, height_ratios=[2, 1])
        ax1 = fig.add_subplot(gs[0])
    else:
        fig = plt.figure(fignum, (9, 6))
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
        # agreement case) but widening for a real large disagreement instead of
        # silently clipping it off-screen
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

    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    if outroot == '':
        outroot = os.path.basename(filename).replace('.fits', '')

    figname = '%s/%s_fluxcal.png' % (outdir, outroot)
    plt.savefig(figname)
    plt.close(fig)

    return figname, note


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


def plot_diagnostic_line_panels(axs, wav, band_flux, overlays=None, refline=None):
    '''
    Fill in a 2x3 (or shorter) grid of axes, one per DIAGNOSTIC_LINES
    window, each showing the 10th/50th/90th percentile band of
    band_flux (fibers x wave, masked or plain) across Sci-telescope
    fibers at each wavelength pixel, plus optional overlay spectra and
    an optional +/-refline pair of reference lines.

    overlays: list of (label, spectrum_1d, color) plotted on top of the
    band -- used by the CFrame's SKY_EAST/SKY_WEST plausibility check;
    pass None/[] where there's nothing to overlay (e.g. the SFrame case,
    where the corresponding SKY_EAST/SKY_WEST extensions don't exist
    post-subtraction).

    refline: if given, draws +/-refline dashed reference lines (e.g.
    MW_5SIGMA) and auto-scales each panel symmetrically around zero
    from the percentile band's own scatter (never so tight the
    reference lines themselves fall off-panel) -- used by the SFrame's
    sky-subtraction-residual check. Left None (no y-limit override,
    matplotlib autoscales) for the CFrame's field-brightness check,
    which isn't residual-shaped and has no natural zero point.

    Only fills axs[:len(DIAGNOSTIC_LINES)]; any extra axes are left
    alone for the caller to hide or reuse.
    '''
    overlays = overlays or []
    filled = np.ma.filled(band_flux, np.nan) if isinstance(band_flux, np.ma.MaskedArray) else np.asarray(band_flux)

    ncols = 3
    for i, (ax, (name, wl, yscale_window)) in enumerate(zip(axs, DIAGNOSTIC_LINES)):
        wmin, wmax = wl - LINE_WINDOW_HALF_WIDTH, wl + LINE_WINDOW_HALF_WIDTH
        idx = (wav > wmin) & (wav < wmax)
        if idx.sum() == 0:
            ax.set_title('%s\n(out of range)' % name)
            continue

        xwav = wav[idx]
        sub = filled[:, idx]
        p10 = np.nanpercentile(sub, 10, axis=0)
        p50 = np.nanpercentile(sub, 50, axis=0)
        p90 = np.nanpercentile(sub, 90, axis=0)

        ax.fill_between(xwav, p10, p90, color='0.85', label='field 10-90%ile' if i == 0 else None)
        ax.plot(xwav, p50, color='k', lw=1.5, label='field median' if i == 0 else None)

        overlay_specs = []
        for label, spec, color in overlays:
            spec_win = np.asarray(spec)[idx]
            ax.plot(xwav, spec_win, color=color, lw=1.2, label=label if i == 0 else None)
            overlay_specs.append(spec_win)

        if refline is not None:
            ax.axhline(refline, ls=':', color='r', lw=1, label=(r'$\pm$ MW 5$\sigma$' if i == 0 else None))
            ax.axhline(-refline, ls=':', color='r', lw=1)

        ax.set_xlim(wmin, wmax)
        if yscale_window is not None:
            y_wmin, y_wmax = yscale_window
            yidx = (xwav >= y_wmin) & (xwav <= y_wmax)
            parts = [p10[yidx], p90[yidx]] + [s[yidx] for s in overlay_specs]
            if refline is not None:
                parts.append(np.array([refline, -refline]))
            combined = np.concatenate(parts)
            combined = combined[np.isfinite(combined)]
            if len(combined) > 0:
                ax.set_ylim(min(0, np.nanmin(combined)), np.nanmax(combined) * 1.15)
        elif refline is not None:
            # residual data has no natural full-range scale like the CFrame's
            # raw field brightness does -- auto-scale symmetrically around
            # zero from the band's own percentile spread, widened if needed
            # so the reference lines never fall off-panel
            combined = np.concatenate([p10, p90])
            combined = combined[np.isfinite(combined)]
            if len(combined) > 0:
                half = max(np.nanmax(np.abs(combined)), abs(refline) * 1.2)
                ax.set_ylim(-half, half)

        ax.set_title('%s (%.0f A)' % (name, wl))
        if i % ncols == 0:
            ax.set_ylabel('FLUX')
        if i >= len(DIAGNOSTIC_LINES) - ncols:
            ax.set_xlabel('Wavelength [Angstrom]')


def get_gaia_cache_dir():
    '''
    The Gaia XP spectra cache directory the lvmdrp flux-calibration
    itself uses (LVM_MASTER_DIR/gaia_cache), so spectra it already
    downloaded during reduction are reused instead of re-querying the
    archive. LVM_MASTER_DIR is normally set as soon as lvmdrp is
    imported (via the sdss_access/tree setup in lvmdrp/__init__.py);
    the local './gaia_cache' fallback only applies if that failed.
    '''
    master_dir=os.getenv('LVM_MASTER_DIR')
    if master_dir:
        return os.path.join(master_dir,'gaia_cache')
    return './gaia_cache'


def get_standard(xx,fiberid):
    '''
    Get the spectrum of a single fiber (by fiberid) from an lvmSFrame/
    lvmCFrame file
    '''
    wave=xx['WAVE'].data
    flux=xx['FLUX'].data[fiberid-1]
    return wave,flux


def get_header_stars(header):
    '''
    Retrieve the (slot, fiberid, gaia_id) triples for the science-
    telescope fibers that land on a Gaia-matched field star, from the
    SCI#ID/SCI#FIB header keywords written by lvmdrp's flux
    calibration (science_sensitivity in fluxCalMethod.py). Not every
    slot 1..15 is populated -- stars that failed acquisition or
    matching leave gaps -- so the original slot number is kept
    alongside each star rather than renumbered sequentially: the
    FLUXCAL_SCI table's columns are named SCI<slot>SEN by slot, not by
    position in this list, and a skipped slot would otherwise silently
    mismatch a star against the wrong column.
    '''
    stars=[]
    for i in range(1,NSCI_MAX+1):
        try:
            gaia_id=header['SCI%dID' % i]
            fiber=header['SCI%dFIB' % i]
        except KeyError:
            continue
        stars.append((i,fiber,gaia_id))
    return stars


def get_std_header_stars(header,xtab):
    '''
    Retrieve the (slot, fiberid, gaia_id) triples for the dedicated
    standard-star fibers used by the STD/MOD flux-calibration methods,
    from the STD#ID/STD#FIB header keywords. Unlike SCI#FIB, STD#FIB
    is an orig_ifulabel string (e.g. "P1-2"), not a raw fiberid, so it
    has to be matched against the SLITMAP to get the numeric fiberid.
    Not every slot 1..15 is populated (ACQ=False, or excluded for
    lacking a Gaia XP spectrum, leaves a gap) -- see get_header_stars
    for why the original slot number must be kept.
    '''
    stars=[]
    for n in range(1,NSCI_MAX+1):
        label=header.get('STD%dFIB' % n)
        if label is None or str(label)=='None':
            continue
        match=xtab[xtab['orig_ifulabel']==label]
        if len(match)==0:
            continue
        gaia_id=header.get('STD%dID' % n)
        if gaia_id is None:
            continue
        stars.append((n,int(match['fiberid'][0]),gaia_id))
    return stars


def _col_valid(table,colname):
    '''
    True if colname exists in table and has at least one finite value
    -- i.e. the pipeline itself didn't exclude this star (e.g. a
    low-signal cut) even if it was successfully acquired (ACQ=True).
    Returns None if the column/table doesn't exist at all.
    '''
    if table is None or colname not in table.columns.names:
        return None
    return np.isfinite(np.asarray(table[colname])).any()


def _plot_star_panel(ax,x,stars,sen_tables,colprefix,gaia):
    '''
    Plot one panel's worth of stars (either the SCI field stars or the
    STD/MOD standard stars): each star's smoothed observed spectrum in
    color, its Gaia XP spectrum overlaid in solid black (skipped for a
    star the pipeline itself excluded, even if the Gaia fetch would
    have succeeded). stars is a list of (slot, fiberid, gaia_id)
    triples (see get_header_stars) -- slot is the original header slot
    number, used (not a renumbered position) to look up each star's
    column, since header slots can have gaps. sen_tables is a list of
    (table, label) pairs used to check exclusion via the column
    colprefix+slot+"SEN" (e.g. "SCI3SEN" or "STD3SEN" -- note
    FLUXCAL_STD and FLUXCAL_MOD share the same STD#SEN column names,
    so they're distinguished by which *table* has a finite value, not
    by column name) -- a star excluded from *every* table in the list
    is drawn dashed/grey and labeled "[excluded]"; used by only some is
    labeled with which.

    Also auto-scales the y-axis from the percentiles of the non-
    excluded stars' spectra only -- an excluded star's near-zero/noisy
    flux would otherwise blow out the log-scale range by many decades.

    Returns (ntried, nfailed) Gaia-retrieval counts, for stars that
    were not excluded.
    '''
    ntried=0
    nfailed=0
    used_flux=[]
    excluded_flux=[]
    for slot,fiber,gaia_id in stars:
        valid=[(name,_col_valid(table,'%s%d%s' % (colprefix,slot,'SEN'))) for table,name in sen_tables]
        used=[v for _,v in valid if v is not None]
        excluded=len(used)>0 and not any(used)
        label='%s%d (fiber %d)' % (colprefix,slot,fiber)

        swave,sflux=get_standard(x,fiber)
        sflux=xsmooth(sflux)

        if excluded:
            ax.semilogy(swave,sflux,lw=0.8,ls='--',alpha=0.5,color='0.5',label=label+' [excluded]')
            excluded_flux.append(sflux)
            continue

        used_names=[name for name,v in valid if v]
        if used_names and len(used_names)<len(valid):
            label+=' [%s only]' % '/'.join(used_names)

        ax.semilogy(swave,sflux,lw=1.0,label=label)
        used_flux.append(sflux)

        ntried+=1
        try:
            gaia.fetch_xp_spectra([gaia_id])
            gwave,gflux=gaia.load_xp_spectra(gaia_id)
            ax.semilogy(gwave,gflux[0],color='k',lw=1.3,alpha=0.85,zorder=10)
            used_flux.append(gflux[0])
        except Exception as e:
            nfailed+=1
            print('Error: Failed on GAIA object %s (fiber %s): %s' % (gaia_id,fiber,e))

    scale_from=used_flux if used_flux else excluded_flux
    if scale_from:
        allvals=np.concatenate(scale_from)
        pos=allvals[np.isfinite(allvals)&(allvals>0)]
        if pos.size:
            zlo=np.nanpercentile(pos,1)
            zhi=np.nanpercentile(pos,99)
            ax.set_ylim(zlo*0.3,zhi*3)

    return ntried,nfailed


def compare_with_gaia(filename='lvmSFrame-00005059.fits',outroot=''):
    '''
    Compare the flux-calibrated spectra of the SCI-method Gaia-matched
    field stars and the STD/MOD-method dedicated standard stars in
    filename to their Gaia XP spectra, in two panels.

    Returns (outfile, message): outfile is the plot filename on
    success and None on failure; message explains why on failure, and
    is a non-empty warning (but still returns an outfile) if only some
    of the stars could be retrieved/plotted.
    '''
    try:
        x=fits.open(filename)
    except Exception as e:
        return None,'Could not open %s (%s)' % (filename,e)

    header=x[0].header
    exposure=header['EXPOSURE']
    mjd=header['MJD']

    sci_stars=get_header_stars(header)

    std_stars=[]
    if 'SLITMAP' in x:
        xtab=Table(x['SLITMAP'].data)
        std_stars=get_std_header_stars(header,xtab)

    if len(sci_stars)==0 and len(std_stars)==0:
        return None,'No SCI#ID/SCI#FIB or STD#ID/STD#FIB header keywords were found in this file'

    sci_sen=x['FLUXCAL_SCI'].data if 'FLUXCAL_SCI' in x else None
    std_sen=x['FLUXCAL_STD'].data if 'FLUXCAL_STD' in x else None
    mod_sen=x['FLUXCAL_MOD'].data if 'FLUXCAL_MOD' in x else None

    gaia=GaiaXPSpectra(cache_dir=get_gaia_cache_dir())

    fig=plt.figure(1,(8,10))
    plt.clf()

    ntried=nfailed=0
    ax1=plt.subplot(2,1,1)
    if sci_stars:
        n_t,n_f=_plot_star_panel(ax1,x,sci_stars,[(sci_sen,'SCI')],'SCI',gaia)
        ntried+=n_t; nfailed+=n_f
        ax1.legend(fontsize=7,ncol=2)
        if n_t==0:
            ax1.text(0.5,0.9,'No SCI star had a valid calibration (all excluded)',
                     transform=ax1.transAxes,ha='center')
    else:
        ax1.text(0.5,0.5,'No SCI stars available',transform=ax1.transAxes,ha='center')
    ax1.set_xlim(3500,9500)
    ax1.set_title('SCI field stars, MJD %d Exposure %d' % (mjd,exposure))

    ax2=plt.subplot(2,1,2)
    if std_stars:
        n_t,n_f=_plot_star_panel(ax2,x,std_stars,[(std_sen,'STD'),(mod_sen,'MOD')],'STD',gaia)
        ntried+=n_t; nfailed+=n_f
        ax2.legend(fontsize=7,ncol=2)
        if n_t==0:
            ax2.text(0.5,0.9,'No STD/MOD star had a valid calibration (all excluded)',
                     transform=ax2.transAxes,ha='center')
    else:
        ax2.text(0.5,0.5,'No STD/MOD stars available',transform=ax2.transAxes,ha='center')
    ax2.set_xlim(3500,9500)
    ax2.set_xlabel(r'Wavelength [\AA]')
    ax2.set_title('STD/MOD standard stars')

    if ntried==0:
        plt.close(1)
        return None,('Failed to retrieve/plot any GAIA-matched star '
                      '(no network access to the GAIA archive, and nothing cached locally)')

    plt.tight_layout()

    if outroot=='':
        word=filename.split('/')
        outroot=word[-1].replace('.fits','')
        outfile='standard_%s.png' % outroot
    else:
        outfile=outroot

    message='' if nfailed==0 else '%d of %d Gaia spectra could not be retrieved/plotted' % (nfailed,ntried)
    return outfile,message


def qual_eval(filename,outname):
    '''
    This is an extra call so this routine can be run from the qual
    evaluation routines.

    Returns (status, message), see compare_with_gaia.
    '''
    outfile,message=compare_with_gaia(filename,outname)
    if outfile==None:
        return False,message
    plt.savefig(outname)
    plt.close()
    return True,message

                
def steer(argv):

    files=[]

    i=1
    while i<len(argv):
        if argv[i].count('-h'):
            print(_usage_from_doc(__doc__))
            return
        elif argv[i][0]=='-':
            print('Error: could not process command line: ',argv)
        else:
            files.append(argv[i])
        i+=1

    for one in files:
        outfile,message=compare_with_gaia(one)
        if outfile!=None:
            plt.savefig(outfile)
            plt.close()
            if message:
                print('Warning: %s' % message)
        else:
            print('Error: %s' % message)
                          
        



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
