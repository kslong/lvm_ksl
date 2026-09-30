#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Create an html file, with various plots, which can be used
as a tool to assess the quality of the reduction of an
exposure, including its sky subtraction.


Command line usage (if any):

    usage: QualSFrame.py [-h] SFrame1 SFrame2 ...

    where 
        -h prints this documentation and exits
        SFrame1 SFrame2 ... are files to be analyzed.


Description::

    This routine reads an lvmSFrame file (or an SFrame-layout
    file from another sky-subtraction method) and constructs an
    html file that contains information from the headers and
    various plots and tables to indicate what the quality of
    the data is: the science and sky spectra (interactive and
    static), the line emission in the subtracted sky, the
    continuum and sky-line subtraction quality, line,
    continuum and [OI] images, and flux-calibration checks.
    The sky-subtraction checks use only FLUX, SKY, IVAR and
    MASK, so any method's file is judged the same way.

Primary routines::

    make_html is the primary driving routine
    make_plotly_spectra   interactive science and sky spectra
    eval_sky_emission     line emission in the subtracted sky
    eval_continuum        continuum subtraction quality
    eval_sky_lines        sky-line subtraction quality
    make_images           line and continuum images
    make_oi_images        [OI]6300 and 5577 maps
    steer handles the inputs

Notes:

    The html file is created in the current working directory
    and the various plots are in a subdirectory figs_qual.

    The links in the html file are relative to the html file
                                       
History::

    240427 ksl Coding begun
    260727 ksl Switched the science/sky telescope pointing header
        keywords read in eval_qual_sframe() and create_overview() from
        the stale commanded/reported pairs (TESCIRA/TESCIDE, POSCIPA,
        POSKYERA/POSKYEDE/POSKYEPA, POSKYWRA/POSKYWDE/POSKYWPA) to
        SCIRA/SCIDEC/SCIPA, SKYERA/SKYEDEC/SKYEPA, SKYWRA/SKYWDEC/
        SKYWPA -- the keywords actually populated by current SFrame
        files, matching the migration already done in rss_combine.py.
        Fixed a 'partition will ignore the mask' UserWarning by
        switching several np.median/np.nanmedian calls on masked
        arrays to np.ma.median. Fixed a 'figure with num: N already
        exists' UserWarning by closing each matplotlib figure right
        after it's saved, in both eval_qual_sframe panels and
        plot_fits_image (eval_standard.py's compare_with_gaia had the
        same issue plus a leak on its GAIA-lookup-failure path, fixed
        there too).
        eval_qual_sframe()'s sky-comparison figure: the top panel
        (SkyE/SkyW noise floor) used a crude fixed 0.2x rescale of the
        full autoscaled range, which let a few cosmic-ray/bad-sky-line
        spikes dominate the visible range. The middle panel (SkyE-SkyW
        total-flux delta) had a real bug -- its y-limits were set from
        a stale `ymax` left over from the top panel's pre-rescale
        get_ylim(), a copy-paste leftover from the science figure's
        semilogy panel (dead giveaway: a commented-out line still
        referencing the science-panel variable) -- which silently
        clipped away every negative excursion of the delta. Both are
        now set via new get_percentile_yscale() (1st/99th percentile
        for the top panel, 1st/99.9th for the middle, since the real
        sky-mismatch spikes there only emerge past the 99th), clamped
        to span at least +/-2x the MW 5 sigma reference line
        (MW_5SIGMA constant) so the axis never zooms in tighter than
        the scale at which the subtraction is already considered good.
        make_html()'s standard-star comparison section now shows the
        specific reason from eval_standard.qual_eval() when the
        comparison fails or partially fails, instead of a generic
        "could not do" message.
    260906 ksl Ported QualCFrame.py's STD/SCI/MOD flux-calibration
        sensitivity comparison (table + overlay/ratio plot) into
        make_html() as a new "Flux Calibration Comparison" section.
        Replaced eval_qual_sframe()'s three hardcoded science-zoom
        panels with the same doublet-aware, 10-90th-percentile-band
        line diagnostics QualCFrame.py uses (six line windows instead
        of three), via new shared eval_standard.
        plot_diagnostic_line_panels(). Moved SENS_BANDS/SENS_METHODS/
        SENS_COLORS/SENS_DISAGREE_WARN, get_fluxcal_curve,
        sensitivity_summary_table, fluxcal_comment,
        eval_sensitivity_comparison, get_header_value/get_header_string,
        and _col_valid into eval_standard.py so this file and
        QualCFrame.py share one implementation instead of two that
        could silently drift apart (see eval_standard.py's History for
        the real header-slot bug this consolidation fixed).
        create_overview() now also reports the FLUXCAL method applied
        and SKYSRC. Fixed the report title's "Asssessment" typo.
    260907 ksl eval_qual_sframe()'s science figure gains a "Continuum
        check" panel: the same sky-subtracted science spectrum as the
        top panel, but y-limits set to the spectrum's own median
        +/-1e-14 instead of a fixed near-zero window -- the near-zero
        window is tuned for sky-line residuals and clips a real
        continuum off the top, hiding broad continuum-level over/
        under-subtraction (a slope or offset spanning the whole band).
        Added a solid orange zero-reference line (a plain black line
        was invisible against the blue spectrum trace) to that panel
        and to every other residual/delta panel in both figures (the
        SkyE/SkyW delta panel and all diagnostic-line sub-panels).
        The SkyE/SkyW figure's bottom row now uses the same 6
        diagnostic-line windows (and matching title format) as the
        science figure's residual check, replacing 3 broader combined
        windows, so the two checks line up panel-for-panel; the figure
        was also resized so its grid cells match the science figure's
        exactly (same column width and row height), making the two
        figures directly comparable side by side.
        create_overview() now returns (xlist, pointing_rows): the
        pointing/Moon/Sun block is a table (Target, RA, Dec., PA, Ang.
        dist., Alt., Illum., Astrometry Src, Shadow Ht) instead of a
        flat text list. Alt. and Shadow Ht come straight from the
        DRP's own SCIALT/SKYEALT/SKYWALT and SKY ..._SH_HGHT header
        keywords rather than being recomputed. Astrometry Src surfaces
        SCIASRC/SKYEASRC/SKYWASRC ('GDR coadd' vs. 'CMD position') so
        it's clear when SkyE/SkyW's PA is just the commanded value
        (SkyE/SkyW aren't actively guided) rather than looking like an
        unexplained inconsistency next to Sci's guider-derived PA.
        Ang. dist. now also covers Moon/Sun separation from the
        science field, not just SkyE/SkyW. make_html() renders the new
        table via xhtml.table() right after the existing bullet list.
        QualCFrame.py's create_overview() and eval_sky_comparison()
        (its analogous pre-sky-subtraction SkyE/SkyW figure) were
        brought up to the same pointing-table/6-window/zero-line/
        sizing conventions -- see its own History entry.
        Removed dead code found along the way: get_yscale()'s debug
        print, eval_qual_sframe()'s unused `xtype` variable (guarded
        by a filename check that could never match real lvmCFrame
        filenames anyway), and an unused plt.ylim() call.
    260907 ksl Renamed QuickLook.py to QualSFrame.py, matching
        QualCFrame.py's naming and making explicit that the two are a
        matched SFrame/CFrame quality-assessment pair. Updated all
        cross-file references (QualCFrame.py, eval_standard.py) and
        the Sphinx docs accordingly; no functional change.
    260930 ksl The line/continuum images now show the SFrame's own,
        sky-subtracted data (quick_map had been reading the matching
        CFrame via GetTelData since 260913).  New "[OI] 6300 Emission
        and Sky-Line Residual Controls" section (make_oi_images): maps
        of [OI]6300 and the 5577 airglow residual on one symmetric
        colour scale, and a table of each map's rms (also as a fraction of that line
        in the subtracted sky) and its correlation with [SII] and 5577.
        New "Sky-Line Subtraction Quality" section (eval_sky_lines): for
        8 bright, isolated sky lines, per-fiber residual rms, noise,
        signed integral and red-blue asymmetry relative to the
        subtracted sky line, summarized in a table, as residual profiles
        (median and 10-90 percentile) and as maps on the sky; uses only
        FLUX/SKY/IVAR/MASK so any sky-subtraction method's SFrame can be
        judged the same way.  Interactive Plotly versions of the full
        science and SkyE/SkyW spectra (make_plotly_spectra) head those
        two sections, so they can be zoomed instead of relying on fixed
        y scales; the static figures follow them.  The Plotly library is
        written once to figs_qual/plotly.min.js, so reports work
        offline.  New "Line Emission in the Subtracted Sky" section
        (eval_sky_emission): [OII], Hbeta, [OIII], Halpha, [NII], [SII],
        [SIII]9531 fitted
        in the subtracted sky, the science total and the raw SkyE/SkyW
        spectra, with the fraction of the field's line flux removed.
        New "Continuum Subtraction Quality" section (eval_continuum):
        per-arm continuum residual in line-free pixels (median, scatter,
        % of sky continuum, vs MW 5 sigma, plane-fit change across the
        field), the step at the b/r and r/z junctions, and per-arm maps.  The new maps use the other images' hot colour map and
        5-95 percentile stretch.
        plot_fits_image() now draws on WCS axes: its RA/Dec tick labels
        had been interpolated linearly between two image corners, which
        put features at the wrong coordinates (0.25 deg off in RA for
        exposure 16998).

'''


import sys
from glob import glob
import os
from astropy.io import ascii, fits
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import xhtml
from astropy.table import Table
from matplotlib.gridspec import GridSpec
from astropy.wcs import WCS
from astropy.coordinates import get_body, solar_system_ephemeris, AltAz, EarthLocation
from astropy.time import Time
import astropy.units as u
from lvm_ksl import quick_map
import plotly.graph_objects as go
from scipy.optimize import curve_fit
from lvm_ksl.GetSkyCont import load_mask, _interp_mask_to_wave
from plotly.subplots import make_subplots
from lvm_ksl import eval_standard
from lvm_ksl.eval_standard import get_header_value, get_header_string


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


def get_moon_info_las_campanas(datetime_utc,verbose=False):
    '''
    Get information about the moon (and sun) as a fuction of UT
    '''
    # Las Campanas Observatory coordinates
    observatory_location = EarthLocation(lat=-29.0089*u.deg, lon=-70.6920*u.deg, height=2281*u.m)

    # Specify the observation time in UT
    obs_time = Time(datetime_utc)
    # print(obs_time.mjd)

    # Set the solar system ephemeris to 'builtin' for faster computation
    with solar_system_ephemeris.set('builtin'):
        # Get the Moon's and Sun's coordinates at the specified time
        moon_coords = get_body('moon', obs_time,location=observatory_location)
        sun_coords = get_body('sun',obs_time,location=observatory_location)

    # Calculate the phase angle (angle between the Sun, Moon, and observer)
    phase_angle = moon_coords.separation(sun_coords).radian

    # Calculate the illuminated fraction of the Moon
    illumination_fraction = (1 - np.cos(phase_angle))/2
    # print('separation',phase_angle,phase_angle*57.29578,illumination_fraction)
    moon_sun_longitude_diff = (moon_coords.ra - sun_coords.ra).wrap_at(360 * u.deg).value
    if moon_sun_longitude_diff>0:
        moon_phase=illumination_fraction/2.
    else:
        moon_phase=1-illumination_fraction/2.

    illumination_fraction*=100.

    # Calculate the Altitude and Azimuth of the Moon from Las Campanas Observatory
    altaz_frame = AltAz(obstime=obs_time, location=observatory_location)
    moon_altaz = moon_coords.transform_to(altaz_frame)
    sun_altaz=sun_coords.transform_to(altaz_frame)



    # Calculate the difference in ecliptic longitudes between the moon and the sun
    delta_longitude = (moon_coords.spherical.lon - sun_coords.spherical.lon).to_value('deg')
    # print('delta_long',delta_longitude)

    # Normalize the difference in ecliptic longitudes to get the moon's phase



    # Print the moon's phase
    # print("Moon's phase:", moon_phase)



    xreturn={
        'SunRA':sun_coords.ra.deg,
        'SunDec':sun_coords.dec.deg,
        'SunAlt': sun_altaz.alt.deg,
        'SunAz': sun_altaz.az.deg,
        'MoonRA': moon_coords.ra.deg,
        'MoonDec': moon_coords.dec.deg,
        'MoonAlt': moon_altaz.alt.deg,
        'MoonAz': moon_altaz.az.deg,
        'MoonPhas': moon_phase,
        'MoonIll': illumination_fraction
    }

    # print(xreturn)

    if verbose:
        for key, value in xreturn.items():
            print(f'{key}: {value}')
    # Return the information
    return xreturn

RADIAN=57.29578
MW_5SIGMA=5.9e-15 # Milky Way sky background 5 sigma sensitivity limit (see the red dotted reference lines)

def distance(r1,d1,r2,d2):
    '''
    distance(r1,d1,r2,d2)
    Return the angular offset between two ra,dec positions
    All variables are expected to be in degrees.
    Output is in degrees

    Note - This routine could easily be made more general
    '''
#    print 'distance',r1,d1,r2,d2
    r1=r1/RADIAN
    d1=d1/RADIAN
    r2=r2/RADIAN
    d2=d2/RADIAN
    xlambda=np.sin(d1)*np.sin(d2)+np.cos(d1)*np.cos(d2)*np.cos(r1-r2)
#    print 'xlambda ',xlambda
    if xlambda>=1.0:
        xlambda=0.0
    else:
        xlambda=np.arccos(xlambda)

    xlambda=xlambda*RADIAN
#    print 'angle ',xlambda
    return xlambda

def scifib(xtab,select='all',telescope=''):
    '''
    Select good fibers from a telescope, of a spefic
    type or all from a telescope from the slitmap table
    of a calbrated file
    '''
    # print(np.unique(xtab['fibstatus']))
    # print(np.unique(xtab['targettype']))
    ztab=xtab[xtab['fibstatus']==0]
    if select=='all':
        ztab=ztab[ztab['targettype']!='standard']
    else:
        ztab=ztab[ztab['targettype']==select]

    if telescope!='' and telescope!='all':
        ztab=ztab[ztab['telescope']==telescope]


    # print('Found %d fibers' % len(ztab))
    return ztab



def limit_spectrum(wave,flux,wmin,wmax):
    '''
    Get a section of the spectrum
    '''

    f=flux[wave>wmin]
    w=wave[wave>wmin]

    f=f[w<wmax]
    w=w[w<wmax]

    return w,f

def get_yscale(f,ymin,ymax):
    new_mask=np.isnan(f.data)
    f.mask = np.logical_or(f.mask, new_mask)
    med=np.ma.median(f)
    zmin=ymin+med
    zmax=ymax+med
    return zmin,zmax


def get_percentile_yscale(arr,low,high,min_half_range=None):
    '''
    Get y-axis limits from the low/high percentiles of arr, ignoring
    NaNs and masked entries. Trims outlier spikes (cosmic rays, bad
    sky lines) directly from the data, rather than via an ad hoc
    rescale factor applied to the full autoscaled range.

    If min_half_range is given, the limits are widened (never
    narrowed) so they span at least +/-min_half_range -- there is no
    diagnostic value in zooming in tighter than the scale at which the
    subtraction is already considered good.
    '''
    if isinstance(arr,np.ma.MaskedArray):
        data=np.ma.filled(arr.astype(float),np.nan)
    else:
        data=np.asarray(arr,dtype=float)
    zmin=np.nanpercentile(data,low)
    zmax=np.nanpercentile(data,high)
    if min_half_range is not None:
        zmin=min(zmin,-min_half_range)
        zmax=max(zmax,min_half_range)
    return zmin,zmax


# SENS_BANDS/SENS_METHODS/SENS_COLORS/SENS_DISAGREE_WARN, get_fluxcal_curve,
# sensitivity_summary_table, fluxcal_comment, and eval_sensitivity_comparison
# now live in eval_standard.py, shared with QualCFrame.py, so the two tools'
# flux-cal comparison logic can't silently drift apart. get_header_value/
# get_header_string are imported from there too (see the import block above)
# so every existing unqualified call site below keeps working unchanged.


def _tel_stats(x,fibers,pct=False):
    '''masked median FLUX and SKY over fibers (and 10/90 percentiles of FLUX if pct)'''
    rows=fibers['fiberid']-1
    bad=x['MASK'].data[rows]!=0
    flux=np.where(bad,np.nan,x['FLUX'].data[rows]).astype(float)
    sky=np.where(bad,np.nan,x['SKY'].data[rows]).astype(float)
    with np.errstate(all='ignore'):
        out=[np.nanmedian(flux,axis=0),np.nanmedian(sky,axis=0)]
        if pct:
            out+=list(np.nanpercentile(flux,[10,90],axis=0))
    return out


def _linear_range(arr):
    '''y range from the 1st/99th percentiles, at least +-2x the MW 5 sigma line'''
    lo,hi=get_percentile_yscale(np.ma.masked_invalid(arr),1,99,min_half_range=2*MW_5SIGMA)
    return [lo,hi]


def make_plotly_spectra(filename='data/lvmSFrame-00011061.fits'):
    '''
    Interactive (Plotly) versions of the full-spectrum plots, which can be
    zoomed rather than relying on a fixed y scale.

    Returns (science html, sky-telescope html): html fragments to embed in
    the report.  The first loads the Plotly library from
    figs_qual/plotly.min.js, written there (once per directory) from the
    installed plotly package, so the report works offline; the second
    relies on it.

    Science: the median sky-subtracted spectrum of the science fibers with
    its 10-90 percentile range across fibers and the +-MW 5 sigma lines
    (top), and the median total (FLUX+SKY) and sky on a log scale
    (bottom).  Sky telescopes: the median sky-subtracted SkyE and SkyW
    spectra (top) and the difference of their total spectra, nearer minus
    further from the science field (bottom).
    '''
    x=fits.open(filename)
    hdr=x['PRIMARY'].header
    xtab=Table(x['SLITMAP'].data)
    wav=np.asarray(x['WAVE'].data,dtype=float)
    f32=lambda a: np.asarray(a,dtype=np.float32)

    sci_med,sci_sky,sci_p10,sci_p90=_tel_stats(x,scifib(xtab,select='science',telescope='Sci'),pct=True)
    zero=dict(color='orange',width=1)
    mw=dict(color='red',width=1,dash='dot')

    fig=make_subplots(rows=2,cols=1,shared_xaxes=True,vertical_spacing=0.08,
                      subplot_titles=('Sky-subtracted science fibers: median and 10-90 percentile range',
                                      'Total (FLUX+SKY) and sky, median over science fibers'))
    fig.add_trace(go.Scatter(x=wav,y=f32(sci_p90),line=dict(width=0),showlegend=False,hoverinfo='skip'),row=1,col=1)
    fig.add_trace(go.Scatter(x=wav,y=f32(sci_p10),line=dict(width=0),fill='tonexty',
                             fillcolor='rgba(31,119,180,0.25)',name='10-90 percentile'),row=1,col=1)
    fig.add_trace(go.Scatter(x=wav,y=f32(sci_med),line=dict(color='rgb(31,119,180)',width=1),name='median'),row=1,col=1)
    fig.add_hline(y=0,line=zero,row=1,col=1)
    fig.add_hline(y=MW_5SIGMA,line=mw,row=1,col=1)
    fig.add_hline(y=-MW_5SIGMA,line=mw,row=1,col=1)
    total=np.clip(sci_med+sci_sky,1e-18,None)
    skyc=np.clip(sci_sky,1e-18,None)
    fig.add_trace(go.Scatter(x=wav,y=f32(total),line=dict(width=1),name='total'),row=2,col=1)
    fig.add_trace(go.Scatter(x=wav,y=f32(skyc),line=dict(width=1),name='sky'),row=2,col=1)
    ceiling=np.nanpercentile(total,99)
    fig.update_yaxes(range=_linear_range(sci_med),exponentformat='e',title_text='FLUX',row=1,col=1)
    fig.update_yaxes(type='log',range=[np.log10(1e-3*ceiling),np.log10(2*ceiling)],exponentformat='e',
                     title_text='FLUX',row=2,col=1)
    fig.update_xaxes(range=[3600,9600],title_text='Wavelength (A)',row=2,col=1)
    fig.update_layout(height=800,margin=dict(l=70,r=20,t=50,b=50),legend=dict(orientation='h',y=1.08))
    location='./figs_qual/'
    if os.path.isdir(location)==False:
        os.mkdir(location)
    jsfile=location+'plotly.min.js'
    if not os.path.isfile(jsfile):
        from plotly.offline import get_plotlyjs
        with open(jsfile,'w') as g:
            g.write(get_plotlyjs())
    sci_html='<script src="%s"></script>\n' % jsfile + fig.to_html(full_html=False,include_plotlyjs=False)

    # sky telescopes
    e_med,e_sky=_tel_stats(x,scifib(xtab,select='SKY',telescope='SkyE'))
    w_med,w_sky=_tel_stats(x,scifib(xtab,select='SKY',telescope='SkyW'))
    try:
        ra,dec=get_header_value(hdr,'SCIRA'),get_header_value(hdr,'SCIDEC')
        near_w=distance(ra,dec,get_header_value(hdr,'SKYWRA'),get_header_value(hdr,'SKYWDEC')) < \
            distance(ra,dec,get_header_value(hdr,'SKYERA'),get_header_value(hdr,'SKYEDEC'))
    except Exception:
        near_w=False
    if near_w:
        delta,dlabel=(w_med+w_sky)-(e_med+e_sky),'SkyW - SkyE (nearer - further)'
    else:
        delta,dlabel=(e_med+e_sky)-(w_med+w_sky),'SkyE - SkyW (nearer - further)'

    fig=make_subplots(rows=2,cols=1,shared_xaxes=True,vertical_spacing=0.08,
                      subplot_titles=('Sky-subtracted SkyE and SkyW fibers (median)',
                                      'Difference of the total spectra: '+dlabel))
    fig.add_trace(go.Scatter(x=wav,y=f32(e_med),line=dict(width=1),name='SkyE'),row=1,col=1)
    fig.add_trace(go.Scatter(x=wav,y=f32(w_med),line=dict(width=1),name='SkyW'),row=1,col=1)
    fig.add_trace(go.Scatter(x=wav,y=f32(delta),line=dict(width=1),name=dlabel),row=2,col=1)
    for r in (1,2):
        fig.add_hline(y=0,line=zero,row=r,col=1)
        fig.add_hline(y=MW_5SIGMA,line=mw,row=r,col=1)
        fig.add_hline(y=-MW_5SIGMA,line=mw,row=r,col=1)
    fig.update_yaxes(range=_linear_range(np.concatenate([e_med,w_med])),exponentformat='e',title_text='FLUX',row=1,col=1)
    fig.update_yaxes(range=_linear_range(delta),exponentformat='e',title_text='FLUX',row=2,col=1)
    fig.update_xaxes(range=[3600,9600],title_text='Wavelength (A)',row=2,col=1)
    fig.update_layout(height=700,margin=dict(l=70,r=20,t=50,b=50),legend=dict(orientation='h',y=1.1))
    sky_html=fig.to_html(full_html=False,include_plotlyjs=False)
    return sci_html,sky_html



# Nebular lines measured in the subtracted sky: name, rest wavelength(s)
# (A).  [OII]3726,3729 is fitted as a doublet (fixed separation, common
# shift and width) and reported as the sum; [SIII]9531 is the brighter of
# the [SIII] pair.  SKY_NEB_MOONLIT are the lines whose fit is unreliable
# when the Moon is up: the solar absorption spectrum in scattered
# moonlight has structure on the scale of the line (Balmer absorption;
# strong absorption either side of [OII]).
SKY_NEB_LINES=[['[OII]3727',[3726.03,3728.82]],['Hbeta',[4861.33]],['[OIII]5007',[5006.84]],
               ['Halpha',[6562.80]],['[NII]6583',[6583.45]],['[SII]6716',[6716.44]],
               ['[SII]6731',[6730.82]],['[SIII]9531',[9530.6]]]
SKY_NEB_MOONLIT=['[OII]3727','Hbeta','Halpha']
SKY_NEB_HALF=6.0          # fit window +-6 A
SKY_NEB_SHIFT=1.5         # line centre allowed within +-1.5 A (~70 km/s)


def _fit_line(wav,spec,lsf,centres):
    '''
    Integrated flux of one emission line, or the summed flux of a doublet:
    Gaussians (common shift within +-SKY_NEB_SHIFT of the rest
    wavelengths, fixed separation, common FWHM between 0.7 and 1.5 times
    the LSF -- nebular lines are barely resolved, and a wider limit lets
    the fit absorb continuum structure) on a linear background, fitted within +-SKY_NEB_HALF of the
    line(s).  The narrow shift range keeps the fit off neighbouring sky
    lines (e.g. the OH lines at 6553.6 and 6568.8 either side of Halpha).
    NaN on failure.
    '''
    ref=np.mean(centres)
    offs=np.array(centres)-ref
    win=(wav>=min(centres)-SKY_NEB_HALF)&(wav<=max(centres)+SKY_NEB_HALF)&np.isfinite(spec)
    if win.sum()<8:
        return np.nan
    xw=wav[win]-ref
    y=spec[win]*1e16
    fw=np.nanmedian(lsf[win]) if lsf is not None else 1.5
    nl=len(centres)

    def g(x,*p):
        sig=p[nl+1]/2.3548
        out=p[nl+2]+p[nl+3]*x
        for k in range(nl):
            out=out+p[k]/(sig*np.sqrt(2*np.pi))*np.exp(-0.5*((x-offs[k]-p[nl])/sig)**2)
        return out

    edge=np.min(np.abs(xw[:,None]-offs[None,:]),axis=1)>3
    b0=np.median(y[edge]) if edge.any() else np.median(y)
    f0=max(np.sum(y-b0)*np.median(np.diff(xw))/nl,1e-3)
    p0=[f0]*nl+[0.,fw,b0,0.]
    lo=[-np.inf]*nl+[-SKY_NEB_SHIFT,0.7*fw,-np.inf,-np.inf]
    hi=[np.inf]*nl+[SKY_NEB_SHIFT,1.5*fw,np.inf,np.inf]
    try:
        p,_=curve_fit(g,xw,y,p0=p0,bounds=(lo,hi),maxfev=5000)
    except Exception:
        return np.nan
    return np.sum(p[:nl])/1e16


def eval_sky_emission(filename='data/lvmSFrame-00011061.fits',outroot='test'):
    '''
    How much nebular line emission is in the sky that was subtracted.

    Each line in SKY_NEB_LINES is fitted (_fit_line) in four median
    spectra: the subtracted sky (SKY, median over science fibers), the
    science fibers' total (FLUX+SKY), and the raw SkyE and SkyW
    spectra (FLUX+SKY of their fibers).  "sky / total" is the fraction of
    the field's median line flux that the sky subtraction removed from
    every fiber.  For a sky from the sky telescopes this measures nebular
    emission (or geocoronal Halpha) in their fields; for a sky taken from
    the science field itself it is the emission "floor" subtracted.

    Returns (table rows, figure name).  The figure shows each line region
    in the four spectra, local continuum removed.
    '''
    x=fits.open(filename)
    xtab=Table(x['SLITMAP'].data)
    wav=np.asarray(x['WAVE'].data,dtype=float)
    lsf=None
    if 'LSF' in x:
        lsf=np.nanmedian(x['LSF'].data[scifib(xtab,select='science',telescope='Sci')['fiberid']-1],axis=0)

    sci_med,sci_sky=_tel_stats(x,scifib(xtab,select='science',telescope='Sci'))
    e_med,e_sky=_tel_stats(x,scifib(xtab,select='SKY',telescope='SkyE'))
    w_med,w_sky=_tel_stats(x,scifib(xtab,select='SKY',telescope='SkyW'))
    spectra=[('subtracted sky',sci_sky),('science total',sci_med+sci_sky),
             ('SkyE raw',e_med+e_sky),('SkyW raw',w_med+w_sky)]

    try:
        moonlit=float(x['PRIMARY'].header['SKY MOON_ALT'])>0
    except (KeyError,ValueError):
        moonlit=False
    table=[['Line','Subtracted sky','Science total (median)','Sky / total (%)','SkyE raw','SkyW raw']]
    for name,wl in SKY_NEB_LINES:
        fl={lab:_fit_line(wav,spec,lsf,wl) for lab,spec in spectra}
        frac=100*fl['subtracted sky']/fl['science total']
        label=name+(' (Moon up: uncertain)' if moonlit and name in SKY_NEB_MOONLIT else '')
        table.append([label,'%.2e' % fl['subtracted sky'],'%.2e' % fl['science total'],
                      '%.0f' % frac if np.isfinite(frac) else '--','%.2e' % fl['SkyE raw'],'%.2e' % fl['SkyW raw']])

    location='./figs_qual/'
    if os.path.isdir(location)==False:
        os.mkdir(location)
    if outroot=='':
        outroot='test'
    figname=location+outroot+'.skyneb.png'
    regions=[('[OII]3726,3729',3715,3740),('Hbeta',4850,4872),('[OIII]5007',4995,5018),
             ('[NII]+Halpha',6540,6595),('[SII]',6705,6742),('[SIII]9531',9518,9543)]
    styles=[dict(color='#1f3b99',ls='-'),dict(color='black',ls='-'),
            dict(color='#1a7d1a',ls='--'),dict(color='#b22222',ls='--')]
    fig,axes=plt.subplots(2,3,figsize=(18,9))
    for ax,(title,lo,hi) in zip(axes.flat,regions):
        r=(wav>=lo)&(wav<=hi)
        for (lab,spec),st in zip(spectra,styles):
            y=spec[r]
            ax.plot(wav[r],y-np.nanpercentile(y,20),lw=2,label=lab,**st)
        ax.axhline(0,color='orange',lw=1)
        ax.set_title(title,fontsize=13)
        ax.set_xlabel('Wavelength (A)')
    for ax in axes[:,0]:
        ax.set_ylabel('FLUX - local continuum')
    axes[0,0].legend(fontsize=10)
    fig.tight_layout()
    fig.savefig(figname)
    plt.close(fig)
    return table,figname


# Continuum check: arms (line-free part used), boundary windows either side
# of the b/r and r/z junctions, and nebular lines masked (+-10 A).
CONT_ARMS=[['b',3700.,5750.],['r',5810.,7450.],['z',7650.,9500.]]
CONT_JUNCTIONS=[['b/r',(5700.,5750.),(5810.,5860.)],['r/z',(7395.,7450.),(7680.,7735.)]]
CONT_NEB_MASK=[3727.,3869.,4102.,4340.,4861.,4959.,5007.,5876.,6300.,6364.,6548.,6563.,6583.,6716.,6731.,
               7136.,7320.,7330.,9069.,9531.]


def eval_continuum(filename='data/lvmSFrame-00011061.fits',outroot='test'):
    '''
    Quantify the continuum left after sky subtraction, per arm, using
    pixels that are free of sky lines (data/sky_mask.fits) and of the
    nebular lines in CONT_NEB_MASK.  Each science fiber's continuum is the
    median FLUX over those pixels in each arm.

    Returns (table rows, figure name).  Per arm the table gives the median
    over fibers and its robust scatter, the same as a percentage of the
    subtracted sky's continuum and in units of the MW 5 sigma level, and
    the change across the field of a plane fitted to the fibers'
    continua (a gradient the sky subtraction did not remove, or a real
    one in the source).  A second block gives the median step at the b/r
    and r/z junctions (red side minus blue side), where the flux
    calibration is weakest.  The figure maps each arm's continuum.
    Real source continuum (stars, nebular continuum) is included: the
    median over fibers and the plane fit are robust to a few stars, not
    to a genuinely bright extended continuum.
    '''
    x=fits.open(filename)
    xtab=Table(x['SLITMAP'].data)
    sci=scifib(xtab,select='science',telescope='Sci')
    rows=sci['fiberid']-1
    wav=np.asarray(x['WAVE'].data,dtype=float)
    bad=x['MASK'].data[rows]!=0
    flux=np.where(bad,np.nan,x['FLUX'].data[rows]).astype(float)
    sky=np.nanmedian(np.where(bad,np.nan,x['SKY'].data[rows]),axis=0)

    mask_file=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'data','sky_mask.fits')
    mw,mb=load_mask(mask_file)
    clean=_interp_mask_to_wave(mw,mb,wav)
    clean&=np.min(np.abs(wav[:,None]-np.array(CONT_NEB_MASK)[None,:]),axis=1)>10

    ra=np.array(sci['ra'],float); dec=np.array(sci['dec'],float)
    xx=(ra-np.nanmean(ra))*np.cos(np.radians(np.nanmean(dec)))*60
    yy=(dec-np.nanmean(dec))*60

    def _rstd(a):
        a=a[np.isfinite(a)]
        return 1.4826*np.median(np.abs(a-np.median(a))) if a.size else np.nan

    def _plane(v):
        '''robust plane fit; returns the max-min of the plane over the fibers'''
        ok=np.isfinite(v)
        if ok.sum()<20:
            return np.nan
        med,r=np.median(v[ok]),_rstd(v)
        ok&=np.abs(v-med)<4*r
        A=np.column_stack([np.ones(ok.sum()),xx[ok],yy[ok]])
        c,*_=np.linalg.lstsq(A,v[ok],rcond=None)
        p=c[0]+c[1]*xx+c[2]*yy
        return np.nanmax(p)-np.nanmin(p)

    table=[['Arm','Median residual','Scatter (fibers)','% of sky continuum','Median / MW 5 sigma',
            'Change across field (plane fit)']]
    arm_cont={}
    with np.errstate(all='ignore'):
        for arm,lo,hi in CONT_ARMS:
            pix=clean&(wav>=lo)&(wav<=hi)
            c=np.nanmedian(flux[:,pix],axis=1)
            arm_cont[arm]=c
            skyc=np.nanmedian(sky[pix])
            med=np.nanmedian(c)
            table.append([arm,'%.2e' % med,'%.2e' % _rstd(c),'%+.2f' % (100*med/skyc),
                          '%+.2f' % (med/MW_5SIGMA),'%.2e' % _plane(c)])
        table.append(['Junction','Median step (red - blue)','Scatter (fibers)','% of sky continuum','Step / MW 5 sigma',''])
        for name,(b0,b1),(r0,r1) in CONT_JUNCTIONS:
            bp=clean&(wav>=b0)&(wav<=b1)
            rp=clean&(wav>=r0)&(wav<=r1)
            step=np.nanmedian(flux[:,rp],axis=1)-np.nanmedian(flux[:,bp],axis=1)
            skyc=np.nanmedian(sky[bp|rp])
            med=np.nanmedian(step)
            table.append([name,'%.2e' % med,'%.2e' % _rstd(step),'%+.2f' % (100*med/skyc),
                          '%+.2f' % (med/MW_5SIGMA),''])

    location='./figs_qual/'
    if os.path.isdir(location)==False:
        os.mkdir(location)
    if outroot=='':
        outroot='test'
    figname=location+outroot+'.continuum.png'
    fig,axes=plt.subplots(1,3,figsize=(20,6.2))
    for ax,(arm,lo,hi) in zip(axes,CONT_ARMS):
        v=arm_cont[arm]
        vmin,vmax=np.nanpercentile(v,[5,95])
        ax.set_facecolor((0.5,0.5,0.5,0.2))
        sc=ax.scatter(ra,dec,c=v,s=14,marker='h',cmap='hot',vmin=vmin,vmax=vmax)
        ax.set_aspect(1/np.cos(np.radians(np.nanmean(dec))))
        ax.invert_xaxis()
        ax.set_xlabel('RA (deg)')
        ax.set_ylabel('Dec (deg)')
        ax.set_title('%s arm continuum residual (%.0f-%.0f A)' % (arm,lo,hi),fontsize=10)
        plt.colorbar(sc,ax=ax,shrink=0.8,label='median FLUX')
    fig.tight_layout()
    fig.savefig(figname)
    plt.close(fig)
    return table,figname



def eval_qual_sframe(filename='data/lvmSFrame-00011061.fits',ymin=-0.2e-13,ymax=1e-13,xmin=3600,xmax=9500,outroot=''):
    '''
    Provide a standard plot for looking at how well the sky subtraction has worked overall
    '''

    try:
        x=fits.open(filename)
    except:
        print('Error: eval_qual: Could not open %s' % filename)
        return

    hdr=x['PRIMARY'].header
    mjd=hdr['MJD']
    exposure=hdr['EXPOSURE']

    ra=get_header_value(hdr,'SCIRA')
    dec=get_header_value(hdr,'SCIDEC')

    ra_sky_e=get_header_value(hdr,'SKYERA')
    dec_sky_e=get_header_value(hdr,'SKYEDEC')

    ra_sky_w=get_header_value(hdr,'SKYWRA')
    dec_sky_w=get_header_value(hdr,'SKYWDEC')


    distance_sky_w=distance(ra,dec,ra_sky_w,dec_sky_w)
    distance_sky_e=distance(ra,dec,ra_sky_e,dec_sky_e)


    xtab=Table(x['SLITMAP'].data)

    science_fibers=scifib(xtab,select='science',telescope='Sci')
    skye_fibers=scifib(xtab,select='SKY',telescope='SkyE')
    skyw_fibers=scifib(xtab,select='SKY',telescope='SkyW')

    wav=x['WAVE'].data
    sci_flux=x['FLUX'].data[science_fibers['fiberid']-1]
    sci_sky=x['SKY'].data[science_fibers['fiberid']-1]
    sci_mask=x['MASK'].data[science_fibers['fiberid']-1]
    sci_flux=np.ma.masked_array(sci_flux,sci_mask)
    sci_sky=np.ma.masked_array(sci_sky,sci_mask)
    sci_flux_med=np.ma.median(sci_flux,axis=0)
    sci_sky_med=np.ma.median(sci_sky,axis=0)


    skye_flux= x['FLUX'].data[skye_fibers['fiberid']-1]
    skye_sky=x['SKY'].data[skye_fibers['fiberid']-1]
    skye_mask=x['MASK'].data[skye_fibers['fiberid']-1]
    skye_flux=np.ma.masked_array(skye_flux,skye_mask)
    skye_sky=np.ma.masked_array(skye_sky,skye_mask)
    skye_flux_med=np.ma.median(skye_flux,axis=0)
    skye_sky_med=np.ma.median(skye_sky,axis=0)

    skyw_flux= x['FLUX'].data[skyw_fibers['fiberid']-1]
    skyw_sky=x['SKY'].data[skyw_fibers['fiberid']-1]
    skyw_mask=x['MASK'].data[skyw_fibers['fiberid']-1]
    skyw_flux=np.ma.masked_array(skyw_flux,skyw_mask)
    skyw_sky=np.ma.masked_array(skyw_sky,skyw_mask)
    skyw_flux_med=np.ma.median(skyw_flux,axis=0)
    skyw_sky_med=np.ma.median(skyw_sky,axis=0)

    fig=plt.figure(1,(12,16))
    plt.clf()
    gs= GridSpec(5, 3, figure=fig)

    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(wav,sci_flux_med,label='Sky-Subtracted Science',zorder=2)
    # ax1.plot(wav,skye_flux_med,label='SkyE-Subtracted SkyE',zorder=1)
    # ax1.plot(wav,skyw_flux_med,label='SkyW-Subtracted SkyW',zorder=0)
    ax1.plot([3600,9600],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    ax1.plot([3600,9600],[-5.9e-15,-5.9e-15],':r')
    ax1.set_xlim(3600,9600)
    ymin,ymax=ax1.get_ylim()
    ax1.set_ylim(-1e-14,ymax)
    ax1.legend()

    # Continuum-subtraction check: same sky-subtracted science spectrum as
    # ax1, but y-limits are set to the spectrum's own median +/-1e-14
    # instead of a fixed near-zero window. ax1's near-zero window is tuned
    # to show sky-line residuals and clips a real (non-zero) continuum off
    # the top, so a broad continuum-level over/under-subtraction -- a slope
    # or offset spanning the whole band -- is invisible there. Centering on
    # the spectrum's own median keeps the panel on-scale regardless of the
    # field's overall brightness while still using a fixed window width, so
    # panels are visually comparable across exposures.
    axc = fig.add_subplot(gs[1, :])
    axc.plot(wav,sci_flux_med,label='Sky-Subtracted Science',zorder=2)
    axc.axhline(0,color='orange',lw=1.5,ls='-',zorder=3)
    axc.plot([3600,9600],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    axc.plot([3600,9600],[-5.9e-15,-5.9e-15],':r')
    axc.set_xlim(3600,9600)
    ymin,ymax=get_yscale(sci_flux_med,-1e-14,1e-14)
    axc.set_ylim(ymin,ymax)
    axc.set_title('Continuum check (Median +/- 1e-14)')
    axc.legend()

    ax2 = fig.add_subplot(gs[2, :])
    ax2.semilogy(wav,sci_flux_med+sci_sky_med,label='Science Total',zorder=2)
    ax2.semilogy(wav,sci_sky_med,label='Science Sky',zorder=1)
    ymax=np.nanmax(sci_flux_med+sci_sky_med)
    ax2.set_ylim(1e-3*ymax,1.1*ymax)
    ax2.set_xlim(3600,9600)
    ax2.legend()


    # doublet-aware sky-subtraction residual check: for each diagnostic line
    # window, the 10-90%ile band (not just the population median) of the
    # sky-subtracted Sci-fiber FLUX across all Sci fibers, which should
    # hover near zero within the +/-MW_5SIGMA reference lines if the sky
    # subtraction is clean. Shares its line-window definitions and per-panel
    # plotting with QualCFrame.py's eval_field_vs_sky_lines (eval_standard.
    # plot_diagnostic_line_panels) -- that CFrame check instead compares raw
    # (pre-subtraction) field brightness to the SKY_EAST/SKY_WEST models,
    # which don't exist as SFrame extensions once the sky is subtracted.
    line_axs = [fig.add_subplot(gs[3 + i // 3, i % 3]) for i in range(6)]
    eval_standard.plot_diagnostic_line_panels(line_axs, wav, sci_flux, refline=MW_5SIGMA)
    line_axs[0].legend(fontsize=8, loc='best')

    plt.tight_layout()

    location='./figs_qual/'

    if os.path.isdir(location)==False:
        os.mkdir(location)

    words=filename.split('/')
    root=words[-1].replace('.fits','')
    figname='%s/%s.png' % (location,root)
    plt.savefig(figname)
    plt.close(fig)

    # Now make another plot for the sky fibers

    # Matches the science figure's per-row/per-column size (12x16 over 5
    # rows/3 cols -> 3.2in tall, 4in wide per cell) so the two figures'
    # panels are directly comparable side by side rather than differing in
    # aspect just because this figure has one fewer row.
    fig=plt.figure(2,(12,12.8))
    plt.clf()
    gs= GridSpec(4, 3, figure=fig)

    ax1 = fig.add_subplot(gs[0, :])
    # ax1.plot(wav,sci_flux_med,label='Sky-Subtracted Science',zorder=2)
    ax1.plot(wav,skye_flux_med,label='SkyE-Subtracted SkyE',zorder=1)
    ax1.plot(wav,skyw_flux_med,label='SkyW-Subtracted SkyW',zorder=0)
    ax1.plot([3600,9600],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    ax1.plot([3600,9600],[-5.9e-15,-5.9e-15],':r')
    ax1.set_xlim(3600,9600)
    ymin,ymax=get_percentile_yscale(np.ma.concatenate([skye_flux_med,skyw_flux_med]),1,99,min_half_range=2*MW_5SIGMA)
    ax1.set_ylim(ymin,ymax)
    ax1.legend()

    ax2 = fig.add_subplot(gs[1, :])
    delta=skyw_flux_med+skyw_sky_med-(skye_flux_med+skye_sky_med)
    if distance_sky_w<distance_sky_e:
        ax2.plot(wav,delta,label='SkyW-SkyE (Nearer-Further)',zorder=1)
    else:
        delta=-delta
        ax2.plot(wav,delta,label='SkyE-SkyW (Nearer-Further)',zorder=1)

    ax2.axhline(0,color='orange',lw=1.5,ls='-',zorder=3)
    ymin,ymax=get_percentile_yscale(delta,1,99.9,min_half_range=2*MW_5SIGMA)
    ax2.set_ylim(ymin,ymax)
    ax2.set_xlim(3600,9600)
    ax2.legend()


    # Same 6 diagnostic line windows (and half-width) as the science
    # figure's per-line residual check (eval_standard.DIAGNOSTIC_LINES /
    # plot_diagnostic_line_panels), so the two figures' bottom panels line
    # up one-to-one instead of the 3 broader, differently-chosen windows
    # this used to show. delta here is already a single difference
    # spectrum (not a fibers x wave array), so each panel is filled
    # directly rather than via plot_diagnostic_line_panels (which expects
    # to compute a per-fiber percentile band).
    line_axs = [fig.add_subplot(gs[2 + i // 3, i % 3]) for i in range(6)]
    for i,(ax,(name,wl,_yscale_window)) in enumerate(zip(line_axs,eval_standard.DIAGNOSTIC_LINES)):
        wmin=wl-eval_standard.LINE_WINDOW_HALF_WIDTH
        wmax=wl+eval_standard.LINE_WINDOW_HALF_WIDTH

        xwav,delta_limit=limit_spectrum(wav,delta,wmin,wmax)
        delta_median=np.ma.median(delta_limit)
        delta_limit=delta_limit-delta_median

        ax.plot(xwav,delta_limit,zorder=1)
        ax.axhline(0,color='orange',lw=1.5,ls='-',zorder=3)
        ax.plot([wmin,wmax],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' if i==0 else None)
        ax.plot([wmin,wmax],[-5.9e-15,-5.9e-15],':r')
        ax.set_xlim(wmin,wmax)
        ymin,ymax=get_yscale(delta_limit,-2e-14,2e-14)
        ax.set_ylim(ymin,ymax)
        ax.set_title('%s (%.0f A)' % (name,wl))
    line_axs[0].legend(fontsize=8,loc='best')

    plt.tight_layout()

    words=filename.split('/')
    root=words[-1].replace('.fits','')
    sky_figname='%s/%s_sky.png' % (location,root)
    plt.savefig(sky_figname)
    plt.close(fig)


    return figname,sky_figname
                 


# get_header_value/get_header_string now live in eval_standard.py (imported
# above), shared with QualCFrame.py.


def create_overview(filename='data/lvmSFrame-00011061.fits'):
    ''' Sumarize information about the  processed data file
    '''
    try:
        x=fits.open(filename)
    except:
        print('Error: Could not open %s' % filename)
        return [],[]

    hdr=x['PRIMARY'].header

    exposure=get_header_value(hdr,'EXPOSURE')
    mjd=get_header_value(hdr,'MJD')
    object_name=get_header_string(hdr,'OBJECT')
    obs_time=get_header_string(hdr,'OBSTIME')
    drp_version=get_header_string(hdr,'DRPVER')
    drp_commit=get_header_string(hdr,'COMMIT')
    fluxcal_method=get_header_string(hdr,'FLUXCAL','Unknown')
    sky_src=get_header_string(hdr,'SKYSRC','Unknown')
    ra=get_header_value(hdr,'SCIRA')
    dec=get_header_value(hdr,'SCIDEC')
    pa=get_header_value(hdr,'SCIPA',default_value=0)
    alt=get_header_value(hdr,'SCIALT')
    sh_hght=get_header_value(hdr,'SKY SCI_SH_HGHT')
    # ASRC records which of the two ways set_telescope_astrometry() (lvmdrp
    # core/astrometry.py) can fill RA/Dec/PA actually applied: 'GDR coadd'
    # means a real guider astrometric solution, 'CMD position' means it
    # fell back to the commanded pointing (guider coadd missing/unsolved).
    # SkyE/SkyW are not actively guided, so their PA is normally the
    # commanded value (often exactly 0) rather than a measured one -- this
    # is usually why a SkyE/SkyW PA looks surprising next to Sci's.
    asrc=get_header_string(hdr,'SCIASRC','Unknown')

    ra_sky_e=get_header_value(hdr,'SKYERA')
    dec_sky_e=get_header_value(hdr,'SKYEDEC')
    pa_sky_e=get_header_value(hdr,'SKYEPA')
    alt_sky_e=get_header_value(hdr,'SKYEALT')
    sh_hght_sky_e=get_header_value(hdr,'SKY SKYE_SH_HGHT')
    asrc_sky_e=get_header_string(hdr,'SKYEASRC','Unknown')

    ra_sky_w=get_header_value(hdr,'SKYWRA')
    dec_sky_w=get_header_value(hdr,'SKYWDEC')
    pa_sky_w=get_header_value(hdr,'SKYWPA')
    alt_sky_w=get_header_value(hdr,'SKYWALT')
    sh_hght_sky_w=get_header_value(hdr,'SKY SKYW_SH_HGHT')
    asrc_sky_w=get_header_string(hdr,'SKYWASRC','Unknown')


    distance_sky_w=distance(ra,dec,ra_sky_w,dec_sky_w)
    distance_sky_e=distance(ra,dec,ra_sky_e,dec_sky_e)

  

    moon_info=get_moon_info_las_campanas(obs_time)
    #for key, value in moon_info.items():
    #    print(f'{key}: {value}')

    # Moon/Sun ang. distance from the science field. lvmdrp's own sky-model
    # header block carries the Moon one (SKY SCI_MOON_SEP) but no Sun
    # equivalent, so both are computed the same way as the SkyE/SkyW
    # distances above for consistency (verified to match SKY SCI_MOON_SEP
    # to the precision reported there).
    distance_moon=distance(ra,dec,moon_info['MoonRA'],moon_info['MoonDec'])
    distance_sun=distance(ra,dec,moon_info['SunRA'],moon_info['SunDec'])

   

    xlist=[]
    xlist.append('Exposure : %d' % exposure)
    xlist.append('MJD      : %d' % mjd)
    xlist.append('Obs. time: %s' % obs_time)
    xlist.append('Object.  : %s' % object_name)
    xlist.append('DRP Version : %s' % drp_version)
    xlist.append('DRP Commit  : %s' % drp_commit)
    xlist.append('Flux-cal method applied : %s' % fluxcal_method)
    xlist.append('Sky source (flux-cal)   : %s' % sky_src)

    # Pointing/moon/sun geometry, as a table (one row per target) rather
    # than one fixed-column-format line per target -- RA/Dec apply to every
    # row, but PA/astrometry source only make sense for the science/sky
    # telescopes, and Illum. only for the Moon, so a shared table with
    # blank cells where a quantity doesn't apply reads more clearly than
    # five differently-shaped printed lines. Ang. dist. is the separation
    # from the science field throughout (for Sci itself, blank). Alt. for
    # Sci/SkyE/SkyW comes straight from the SCIALT/SKYEALT/SKYWALT header
    # keywords (already computed by the DRP), not recomputed here.
    # Column order groups Ang. dist./Alt./Illum. together since together
    # they indicate how bright the general sky background should be;
    # Astrometry source ('GDR coadd' vs. 'CMD position', from SCIASRC/
    # SKYEASRC/SKYWASRC) is included because SkyE/SkyW are not actively
    # guided -- their PA is normally just the commanded value (frequently
    # exactly 0) rather than a measured one, which otherwise looks like an
    # inconsistency next to Sci's guider-derived PA. Shadow height is last
    # since, unlike the rest of the table, it speaks to geocoronal emission
    # rather than general sky brightness.
    pointing_rows=[['Target','RA','Dec.','PA','Ang. dist.','Alt.','Illum. (%)','Astrometry Src','Shadow Ht (km)']]
    pointing_rows.append(['Science','%.2f' % ra,'%.2f' % dec,'%.2f' % pa,'','%.2f' % alt,'',asrc,'%.1f' % sh_hght])
    pointing_rows.append(['SkyE','%.2f' % ra_sky_e,'%.2f' % dec_sky_e,'%.2f' % pa_sky_e,'%.2f' % distance_sky_e,'%.2f' % alt_sky_e,'',asrc_sky_e,'%.1f' % sh_hght_sky_e])
    pointing_rows.append(['SkyW','%.2f' % ra_sky_w,'%.2f' % dec_sky_w,'%.2f' % pa_sky_w,'%.2f' % distance_sky_w,'%.2f' % alt_sky_w,'',asrc_sky_w,'%.1f' % sh_hght_sky_w])
    pointing_rows.append(['Moon','%.2f' % moon_info['MoonRA'],'%.2f' % moon_info['MoonDec'],'','%.2f' % distance_moon,'%.2f' % moon_info['MoonAlt'],'%.2f' % moon_info['MoonIll'],'',''])
    pointing_rows.append(['Sun','%.2f' % moon_info['SunRA'],'%.2f' % moon_info['SunDec'],'','%.2f' % distance_sun,'%.2f' % moon_info['SunAlt'],'','',''])

    return xlist,pointing_rows



def calculate_percentiles(arr, percentiles):
    '''
    Calculate specified percentiles, handling NaN values
    '''
    result = np.nanpercentile(arr, percentiles)
    
    return result



def plot_fits_image(filename,title='Cont.(5000-8000)',outname='test.png'):
    '''
    Create an image of an LVM RSS image in a specific wavelentth range
    '''
    # Read the FITS file
    hdul = fits.open(filename)
    data = hdul[0].data
    header = hdul[0].header
    
    # Get the WCS information
    wcs = WCS(header)
    
    # Get the 5th and 95th percentiles of the image data
    min_val, max_val = np.nanpercentile(data, [5, 95])

    # Plot the image on WCS axes, so RA/Dec labels follow the image's
    # own projection and rotation
    fig=plt.figure(figsize=(8, 8))
    ax=fig.add_subplot(projection=wcs)
    cmap=plt.get_cmap('hot').copy()
    cmap.set_bad(color='gray', alpha=0.2)
    im=ax.imshow(data, cmap=cmap, vmin=min_val, vmax=max_val, origin='lower')
    plt.colorbar(im, ax=ax, label='Intensity', shrink=0.8)
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    ax.set_title(title)
    ax.grid(color='white', ls='dotted')
    

    plt.savefig(outname)
    plt.close(fig)




def make_images(filename='data/llvmSFrame-00011061.fits',outroot='test'):
    '''
    Make multiple images of the CFrame of SFrame data
    '''
    xha=['Ha',[6560.,6566.],[6590.,6630.]]
    xs2=['SII',[6710.,6735.],[6740.,6760.]]
    cont=['Cont',[5299.,6200.],None]

    ha_file=quick_map.doit(filename,xha[0],xha[1],xha[2])
    s2_file=quick_map.doit(filename,xs2[0],xs2[1],xs2[2])
    c_file=quick_map.doit(filename,cont[0],cont[1],cont[2])

    location='./figs_qual/'

    if os.path.isdir(location)==False:
        os.mkdir(location)

    if outroot=='':
        outroot='test'

    

    cont_plot=location+outroot+'.cont.png'
    ha_plot=location+outroot+'.ha.png'
    s2_plot=location+outroot+'.s2.png'

    plot_fits_image(filename=c_file,title='Cont.(5200-6200)',outname=cont_plot)
    plot_fits_image(filename=ha_file,title=r'H$\alpha$',outname=ha_plot)
    plot_fits_image(filename=s2_file,title='SII',outname=s2_plot)
    return ha_plot,s2_plot,cont_plot


# [OI]6300 map and its airglow control: name, line window, continuum
# window (A).  Continuum windows are sky-line-free in data/sky_mask.fits;
# 6315.5-6320 also stays clear of [SIII]6312.  [OI]6364 is not mapped:
# atomic physics fixes it at 1/3 of 6300, so it adds no information.
OI_BANDS=[['OI6300',[6297.,6304.],[6315.5,6320.]],
          ['Sky5577',[5574.,5581.],[5587.5,5593.]]]
SII_BAND=['SII',[6710.,6735.],[6740.,6760.]]


def _band_level(wav,spec,band,cont):
    '''band-mean of spec minus the continuum-window mean (same measure as quick_map)'''
    b=(wav>=band[0])&(wav<=band[1])
    c=(wav>=cont[0])&(wav<=cont[1])
    return np.nanmean(spec[b])-np.nanmean(spec[c])


def _show_map(ax,data,vmin,vmax,cmap,title):
    '''one map panel on WCS axes (ax created with projection=wcs)'''
    cm=plt.get_cmap(cmap).copy()
    cm.set_bad(color='gray', alpha=0.2)
    im=ax.imshow(data, cmap=cm, vmin=vmin, vmax=vmax, origin='lower')
    ax.coords[0].set_ticks(number=4)
    ax.coords[1].set_ticks(number=4)
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    ax.set_title(title)
    ax.grid(color='black', ls='dotted', alpha=0.3)
    plt.colorbar(im, ax=ax, shrink=0.8, label='band mean - continuum (FLUX units)')


def make_oi_images(filename='data/lvmSFrame-00011061.fits',outroot='test'):
    '''
    Map of [OI]6300 -- to look for [OI] emission from the source -- next to
    the 5577 airglow line as a control: 5577 has no nebular contribution,
    so its map is the sky-subtraction residual pattern alone.  Both are
    shown relative to their own median on one colour scale.

    Returns (figure name, table rows).  The table gives, for each map,
    the median and robust rms over the IFU, the same rms as a fraction
    of that line's own level in the subtracted sky (SKY, median over
    science fibers), and the pixel correlation with the [SII] map (made
    but not shown: it duplicates the [SII] image above) and with 5577.
    '''
    x=fits.open(filename)
    xtab=Table(x['SLITMAP'].data)
    sci=scifib(xtab,select='science',telescope='Sci')
    wav=x['WAVE'].data
    rows=sci['fiberid']-1
    sky=np.ma.median(np.ma.masked_array(x['SKY'].data[rows],x['MASK'].data[rows]!=0),axis=0).filled(np.nan)

    maps={}
    wcs=None
    for name,band,cont in OI_BANDS+[SII_BAND]:
        mapfile=quick_map.doit(filename,name,list(band),list(cont))
        if mapfile is None:
            print('Error: make_oi_images: could not make the %s map of %s' % (name,filename))
            return None,[]
        with fits.open(mapfile) as m:
            maps[name]=np.array(m[0].data,dtype=float)
            if wcs is None:
                wcs=WCS(m[0].header)

    location='./figs_qual/'
    if os.path.isdir(location)==False:
        os.mkdir(location)
    if outroot=='':
        outroot='test'
    figname=location+outroot+'.oi.png'

    # The sky-line maps can carry a large uniform offset (e.g. a DRP sky
    # telescope with brighter [OI] than the science field), which would
    # hide any structure: each is shown relative to its own median, both on
    # one scale (the [OI]6300 map's 5-95 percentiles, the colour map and
    # stretch of the other images).
    meds={k:np.nanmedian(maps[k]) for k in maps}
    lo,hi=np.nanpercentile(maps['OI6300']-meds['OI6300'],[5,95])
    fig=plt.figure(figsize=(14,6.5))
    for i,(name,title) in enumerate((('OI6300','[OI] 6300'),('Sky5577','5577 airglow residual'))):
        ax=fig.add_subplot(1,2,i+1,projection=wcs)
        _show_map(ax,maps[name]-meds[name],lo,hi,'hot',
                  '%s\nminus its median, %.2e' % (title,meds[name]))
    fig.tight_layout()
    fig.savefig(figname)
    plt.close(fig)

    def _corr(a,b):
        ok=np.isfinite(a)&np.isfinite(b)
        return np.corrcoef(a[ok],b[ok])[0,1] if ok.sum()>10 else np.nan

    table=[['Map','Median','Robust rms','rms / sky line','Corr. with [SII]','Corr. with 5577']]
    for name,band,cont in OI_BANDS:
        d=maps[name]
        med=np.nanmedian(d)
        rms=1.4826*np.nanmedian(np.abs(d-med))
        table.append([name,'%.2e' % med,'%.2e' % rms,'%.4f' % (rms/_band_level(wav,sky,band,cont)),
                      '%.2f' % _corr(d,maps['SII']),'%.2f' % _corr(d,maps['Sky5577'])])
    return figname,table


# Bright, isolated sky lines for the sky-line subtraction check (air, A),
# chosen from the SKY spectrum to avoid nebular lines.  [OI]6300 is also
# emitted by shocked gas, so it is reported but left out of the summaries.
SKYLINE_CHECK=[5577.34,6300.30,6863.96,7340.89,7993.33,8399.18,8885.85,9375.98]
SKYLINE_SOURCE=[6300.30]
SKYLINE_HALF=4.0                 # line window: +-4 A
SKYLINE_SIDE=(6.0,15.0)          # local continuum: 6-15 A either side


def eval_sky_lines(filename='data/lvmSFrame-00011061.fits',outroot='test'):
    '''
    Quantify how well bright sky lines were subtracted, from FLUX, SKY,
    IVAR and MASK alone -- so any sky-subtraction method writing an
    SFrame-layout file can be judged the same way.

    For each line in SKYLINE_CHECK and each good science fiber, the
    residual (FLUX minus its local continuum) within +-SKYLINE_HALF of the
    line is compared with the subtracted sky line (median SKY over the
    science fibers, minus its local continuum)::

        rms       rms of the residual / sky-line peak
        noise     rms expected from IVAR / sky-line peak
        integral  summed residual / summed sky line (signed: > 0 means
                  under-subtracted, e.g. a throughput too low)
        asym      (red half - blue half) of the residual / summed sky line
                  (a wavelength offset gives an antisymmetric residual)

    Returns (table rows, profile figure, map figure).  The table gives the
    median over fibers of each quantity per line, the systematic part
    sqrt(rms^2 - noise^2), and the robust scatter of the integral; the
    profile figure shows the median and 10-90 percentile residual profile
    at each line; the map figure shows each fiber's median rms and median
    integral over the lines that are not also nebular.
    '''
    x=fits.open(filename)
    xtab=Table(x['SLITMAP'].data)
    sci=scifib(xtab,select='science',telescope='Sci')
    rows=sci['fiberid']-1
    wav=x['WAVE'].data
    dw=np.median(np.diff(wav))
    bad=x['MASK'].data[rows]!=0
    flux=np.where(bad,np.nan,x['FLUX'].data[rows]).astype(float)
    ivar=np.where(bad,np.nan,x['IVAR'].data[rows]).astype(float)
    sky=np.nanmedian(np.where(bad,np.nan,x['SKY'].data[rows]),axis=0)

    nfib=len(rows)
    res={k:np.full((len(SKYLINE_CHECK),nfib),np.nan) for k in ('rms','noise','integral','asym')}
    profiles=[]
    with np.errstate(all='ignore'):
        for i,line in enumerate(SKYLINE_CHECK):
            win=np.abs(wav-line)<=SKYLINE_HALF
            side=(np.abs(wav-line)>SKYLINE_SIDE[0])&(np.abs(wav-line)<SKYLINE_SIDE[1])
            sky_line=sky[win]-np.nanmedian(sky[side])
            peak=np.nanmax(sky_line)
            total=np.nansum(sky_line)*dw
            r=flux[:,win]-np.nanmedian(flux[:,side],axis=1)[:,None]
            res['rms'][i]=np.sqrt(np.nanmean(r**2,axis=1))/peak
            res['noise'][i]=np.sqrt(np.nanmean(np.where(ivar[:,win]>0,1/ivar[:,win],np.nan),axis=1))/peak
            res['integral'][i]=np.nansum(r,axis=1)*dw/total
            red=wav[win]>line
            res['asym'][i]=(np.nansum(r[:,red],axis=1)-np.nansum(r[:,~red],axis=1))*dw/total
            profiles.append((wav[win]-line,np.nanpercentile(r/peak,[10,50,90],axis=0),sky_line/peak))

    def _rstd(a):
        a=a[np.isfinite(a)]
        return 1.4826*np.median(np.abs(a-np.median(a))) if a.size else np.nan

    table=[['Line (A)','rms / peak','noise / peak','systematic / peak','integral (%)',
            'integral scatter (%)','asymmetry (%)']]
    clean=np.array([l not in SKYLINE_SOURCE for l in SKYLINE_CHECK])
    for i,line in enumerate(SKYLINE_CHECK):
        rms,noise=np.nanmedian(res['rms'][i]),np.nanmedian(res['noise'][i])
        label='%.2f' % line+('' if clean[i] else ' (also nebular [OI])')
        table.append([label,'%.4f' % rms,'%.4f' % noise,'%.4f' % np.sqrt(max(rms**2-noise**2,0)),
                      '%+.2f' % (100*np.nanmedian(res['integral'][i])),'%.2f' % (100*_rstd(res['integral'][i])),
                      '%+.2f' % (100*np.nanmedian(res['asym'][i]))])
    rms=np.nanmedian(res['rms'][clean]); noise=np.nanmedian(res['noise'][clean])
    table.append(['median, sky-only lines','%.4f' % rms,'%.4f' % noise,'%.4f' % np.sqrt(max(rms**2-noise**2,0)),
                  '%+.2f' % (100*np.nanmedian(res['integral'][clean])),
                  '%.2f' % (100*np.nanmedian([_rstd(v) for v in res['integral'][clean]])),
                  '%+.2f' % (100*np.nanmedian(res['asym'][clean]))])

    location='./figs_qual/'
    if os.path.isdir(location)==False:
        os.mkdir(location)
    if outroot=='':
        outroot='test'

    # residual profiles, one panel per line, all on the same y scale
    prof_name=location+outroot+'.skylines.png'
    ncol=4
    nrow=int(np.ceil(len(SKYLINE_CHECK)/ncol))
    fig,axes=plt.subplots(nrow,ncol,figsize=(16,4*nrow),sharey=True)
    for ax,line,(dx,pct,shape) in zip(axes.flat,SKYLINE_CHECK,profiles):
        ax.fill_between(dx,pct[0],pct[2],color='C0',alpha=0.3,label='10-90%')
        ax.plot(dx,pct[1],color='C0',label='median')
        ax.plot(dx,0.05*shape,'k:',label='5% of sky line')
        ax.axhline(0,color='orange',lw=1)
        ax.set_title('%.2f' % line+('' if line not in SKYLINE_SOURCE else '  (also nebular)'))
        ax.set_xlabel(r'$\Delta\lambda$ (A)')
    for ax in axes[:,0]:
        ax.set_ylabel('residual / sky-line peak')
    for ax in list(axes.flat)[len(SKYLINE_CHECK):]:
        ax.set_visible(False)
    axes.flat[0].set_ylim(-0.08,0.08)
    axes.flat[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(prof_name)
    plt.close(fig)

    # per-fiber maps on the sky
    map_name=location+outroot+'.skyline_map.png'
    ra=np.array(sci['ra'],float); dec=np.array(sci['dec'],float)
    fib_rms=np.nanmedian(res['rms'][clean],axis=0)
    fib_int=np.nanmedian(res['integral'][clean],axis=0)
    fig,axes=plt.subplots(1,2,figsize=(15,6.5))
    for ax,val,title in (
            (axes[0],fib_rms,'rms / sky-line peak (median over sky-only lines)'),
            (axes[1],100*fib_int,'integrated residual, % of sky line (median over sky-only lines)')):
        vmin,vmax=np.nanpercentile(val,[5,95])
        cmap='hot'
        ax.set_facecolor((0.5,0.5,0.5,0.2))
        sc=ax.scatter(ra,dec,c=val,s=18,marker='h',cmap=cmap,vmin=vmin,vmax=vmax)
        ax.set_aspect(1/np.cos(np.radians(np.nanmean(dec))))
        ax.invert_xaxis()
        ax.set_xlabel('RA (deg)')
        ax.set_ylabel('Dec (deg)')
        ax.set_title(title,fontsize=10)
        plt.colorbar(sc,ax=ax,shrink=0.8)
    fig.tight_layout()
    fig.savefig(map_name)
    plt.close(fig)
    return table,prof_name,map_name

plotly_comment='''
Interactive plot (drag to zoom, double-click to reset, click a legend entry to hide it): the median
sky-subtracted spectrum of the science fibers with its 10-90 percentile range across fibers, and,
below, the median total spectrum and sky on a log scale.  The x axes are linked.  The static figure
below shows the same spectra at fixed scales plus close-ups of the diagnostic line regions.
'''

science_plot_comment='''
The median sky subtracted spectrum from the science fibers.  The top panel shows the median spectrum,
scaled to highlight sky-line residuals near zero.  The second panel shows the same spectrum but scaled
to the spectrum's own median +/- 1e-14, to make broad continuum-level over/under-subtraction (a slope or
offset spanning the whole band) visible without being clipped by the top panel's near-zero window.
The third panel shows the sum of the flux and sky, and just the sky.  The three panels at the
bottom show the median scence spectra (after a crude contiumm subtraction) in three wavelength ranges, corresponding to Hbeta-[0II], Halpha-[SII], and [SIII]9071.
Several of the plots also have dashed lines which outline the 5 sigma
sensitivity limit in the Milky Way.
'''

sky_plot_comment='''
Comparisons of the median spectra in the SkyE and SkyW telescopes.  The top panel shows the median spectrum 
in each of the two sky telescopes (after sky subtraction.) The middle panel shows the difference in the the 
total fluxes (with sky included) in the two sky telescopes. The difference is computed by subtracting the 
spectrum of the sky telescope that is furthers from the science target to from the spectrum of the nearer sky telescope.
The bottom panels show the differences in sky subtracted spectra in the same six diagnostic line regions used for the
science-spectrum check above ([OII]3727, Hbeta4861, [OIII]4959,5007, Halpha6563, [SII]6717,6731, [SIII]9533), so the
two checks line up panel-for-panel. (Note that at present,
the lvmdrp uses the sky calculated for the science telescope 
for subtracting sky from the sky telescopes. This implies that what is presented in this figure tells one 
mostly about the differences in the sky in the two telescopes.)
'''

oi_image_comment='''
A map of [OI]6300, to look for [OI] emission from the source.  [OI]6300 is also a bright airglow
line, so its map in sky-subtracted data contains a residual from the sky subtraction as well as any
real emission.  The 5577 airglow line, which has no nebular contribution, is shown alongside as a
control: its map is the sky-subtraction residual pattern by itself.  Real [OI] emission should
follow the [SII] image above (in shocks) rather than 5577.  ([OI]6364 is not shown: it is always
1/3 of 6300.)  Both maps are shown relative to their own median (given in the panel title; a large
median means the subtracted sky's line was brighter or fainter than the science field's), on one
shared colour scale, displayed like the images above (5th to 95th percentile of the [OI]6300 map).  Each map is the mean FLUX in a narrow
line window minus the mean in a nearby sky-line-free continuum window.  In the table, "rms / sky
line" is the map's robust rms as a fraction of the same line measured in the subtracted sky (SKY,
median over science fibers); "Corr." are pixel correlations with the [SII] image above and with the
5577 map.
'''

sky_emission_comment='''
Whether the sky that was subtracted contains nebular line emission.  Each line is fitted in four
median spectra: the subtracted sky (SKY, median over science fibers), the science fibers' total
(FLUX+SKY, i.e. before subtraction), and the raw spectra of the SkyE and SkyW telescopes.  "Sky /
total" is the percentage of the field's median line flux that the sky subtraction removed from every
fiber.  If the sky came from the sky telescopes, a non-zero value means nebular emission (or
geocoronal Halpha) in their fields; if it came from the science field itself, it is the emission
"floor" removed along with the sky, so fluxes in the sky-subtracted data are relative to it.  The
panels show the same four spectra around each line, with a local continuum removed.  When the Moon
is up, scattered moonlight puts the solar absorption spectrum into every spectrum (Balmer absorption
at Hbeta and Halpha, strong absorption either side of [OII], and a Ca I line at 6717.6 A next to
[SII]6716), so a fitted line "flux" can be negative.  With the Moon up, [OII], Hbeta and Halpha are
marked uncertain: their values change by tens of percent with reasonable changes to the fit, while
the other lines change by much less; check the panels.
'''

continuum_comment='''
How much continuum is left after sky subtraction.  For each science fiber the continuum is the
median FLUX over pixels free of sky lines and nebular lines, separately in each arm (b 3700-5750,
r 5810-7450, z 7650-9500 A; the arm junctions are excluded).  The table gives the median over fibers
and its fiber-to-fiber scatter, the median as a percentage of the subtracted sky's continuum and in
units of the Milky Way 5 sigma level, and the change across the field of a plane fitted to the
fibers' continua -- a gradient left by the sky subtraction, or a real one in the source.  The second
block gives the median step across the b/r and r/z junctions (red side minus blue side), where the
flux calibration is weakest.  Real continuum from the source (stars, nebular continuum) is included
in all of these; the medians are robust to a few stars.  The maps show each fiber's continuum per arm.
'''

skyline_comment='''
How well bright sky lines were subtracted, measured from the file's FLUX, SKY and IVAR only, so
files from any sky-subtraction method can be compared.  For each line and each science fiber the
residual within +-4 A of the line (after removing the fiber's local continuum) is compared with the
sky line that was subtracted (median SKY over science fibers).  Table columns (medians over
fibers): "rms / peak", the residual rms as a fraction of the sky-line peak; "noise / peak", the rms
expected from the IVAR alone; "systematic / peak", the part of the rms not explained by noise;
"integral", the summed residual as a percentage of the sky line (positive = under-subtracted, as
from too low a throughput; negative = over-subtracted); "integral scatter", its fiber-to-fiber
scatter; "asymmetry", red half minus blue half of the residual (a wavelength offset gives an
antisymmetric residual).  [OI]6300 is also emitted by shocked gas, so it is listed but left out of
the summary row and the maps.  The profile plot shows the median residual and its 10-90 percentile
range across fibers at each line, on a common scale, with 5% of the sky line (dotted) for
reference: a symmetric bump or dip indicates a throughput mismatch, an S-shape a wavelength offset,
and a W or M shape a difference in line width.  The maps show each fiber's median rms and median
integrated residual over the sky-only lines.
'''

image_comment='''
Line and continuum images of fluxes the science telescope. The emission line emission images use a nearby band pass for
continuum subtraction.
Note that the emisison lime bandpasses for the LMC and SMC images are adjusted 
for the red shifts of these galaxies.  In all other cases, no shifts are assumed.  The 
The images are displayed linearly between the 5th and 95th percentile
'''

def make_html(filename='data/lvmSFrame-00011061.fits', outroot=''):
    '''
    Create an html file that contains summary information about the
    data quality of an lvmdrp analyzed exposrue
    '''

    if outroot=='':
        words=filename.split('/')
        outroot=words[-1].replace('.fits','')

    string=xhtml.begin('LVMDRP SFrame Quality Assessment for %s' % filename)
    string+=xhtml.hline()

    overview_list,pointing_rows=create_overview(filename)
    string+=xhtml.add_list(overview_list)
    string+=xhtml.table(pointing_rows)

    string+=xhtml.hline()
    string+=xhtml.h2('Science Spectrum')
    sci_html,sky_html=make_plotly_spectra(filename)
    string+=xhtml.paragraph(plotly_comment)
    string+=sci_html
    string+=xhtml.paragraph(science_plot_comment)
    

    figname,sky_figname= eval_qual_sframe(filename,ymin=-0.2e-13,ymax=1e-13,xmin=3600,xmax=9500)

    string+=xhtml.image('%s' % (figname),width=900,height=1500)
    string+=xhtml.hline()
    string+=xhtml.h2('SkyE and SkyW  Spectra')
    string+=sky_html
    string+=xhtml.paragraph(sky_plot_comment)

    string+=xhtml.image('%s' % (sky_figname),width=900,height=960)
    string+=xhtml.hline()
    string+=xhtml.h2('Line and Continuum images')

    ha_plot,s2_plot,cont_plot=make_images(filename,outroot)

    string+=xhtml.paragraph(image_comment)

    string+=xhtml.image('%s' % (ha_plot),width=900,height=900)
    string+=xhtml.image('%s' % (s2_plot),width=900,height=900)
    string+=xhtml.image('%s' % (cont_plot),width=900,height=900)

    string+=xhtml.hline()
    string+=xhtml.h2('Line Emission in the Subtracted Sky')
    string+=xhtml.paragraph(sky_emission_comment)
    se_table,se_fig=eval_sky_emission(filename,outroot)
    string+=xhtml.table(se_table)
    string+=xhtml.image('%s' % (se_fig),width=1100,height=550)

    string+=xhtml.hline()
    string+=xhtml.h2('Continuum Subtraction Quality')
    string+=xhtml.paragraph(continuum_comment)
    ct_table,ct_fig=eval_continuum(filename,outroot)
    string+=xhtml.table(ct_table)
    string+=xhtml.image('%s' % (ct_fig),width=1100,height=341)

    string+=xhtml.hline()
    string+=xhtml.h2('Sky-Line Subtraction Quality')
    string+=xhtml.paragraph(skyline_comment)
    sl_table,sl_prof,sl_map=eval_sky_lines(filename,outroot)
    string+=xhtml.table(sl_table)
    string+=xhtml.image('%s' % (sl_prof),width=1000,height=500)
    string+=xhtml.image('%s' % (sl_map),width=1000,height=433)

    string+=xhtml.hline()
    string+=xhtml.h2('[OI] 6300 Emission and Sky-Line Residual Control')
    string+=xhtml.paragraph(oi_image_comment)
    oi_plot,oi_table=make_oi_images(filename,outroot)
    if oi_plot:
        string+=xhtml.table(oi_table)
        string+=xhtml.image('%s' % (oi_plot),width=1000,height=464)
    else:
        string+=xhtml.paragraph('Could not make the [OI] maps')

    string+=xhtml.hline()
    string+=xhtml.h2('Flux Calibration Comparison (STD / SCI / MOD)')

    hdr = fits.getheader(filename, 0)
    string+=xhtml.table(eval_standard.sensitivity_summary_table(hdr))
    string+=xhtml.paragraph(eval_standard.fluxcal_comment)

    fluxcal_figname,fluxcal_note = eval_standard.eval_sensitivity_comparison(filename, outroot, fignum=3, outdir='./figs_qual/')
    if fluxcal_figname:
        string+=xhtml.image(fluxcal_figname,width=900,height=900)
        if fluxcal_note:
            string+=xhtml.paragraph('Note: %s' % fluxcal_note)
    else:
        string+=xhtml.paragraph('Could not compare flux-cal methods: %s' % fluxcal_note)

    string+=xhtml.hline()

    string+=xhtml.paragraph('Comparision between the flux calibrated star fibers to the GAIA spectra of the stars')

    outname='figs_qual/standard_%s.png' % outroot
    status,message=eval_standard.qual_eval(filename,outname)

    if status==True:
        string+=xhtml.image(outname,width=900,height=900)
        if message:
            string+=xhtml.paragraph('Warning: %s' % message)
    else:
        string+=xhtml.paragraph('Could not compare standard stars to GAIA: %s' % message)
    
    string+=xhtml.hline()

    # Finally write out the html file
    # print(string)
    g=open(outroot+'.html','w')
    g.write(string)
    g.close()


def steer(argv):
    '''
    Just a steering routine
    '''

    i=1
    files=[]
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(_usage_from_doc(__doc__))
            return
        elif argv[i][0]=='-':
            print('Error: Unknown optional parameter; improperly formatted command line: ',argv)
            return
        else:
            files.append(argv[i])
        i+=1

    for one in files:
        make_html(one)




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
