#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Create an html file, with various plo5s, which can be used
as a tool to assess the quality of of the lvmdrp reduction
of an exposure.


Command line usage (if any):

    usage: QuickLook.py [-h] SFrame1 SFrame2 ...

    where 
        -h prints this documentation and exits
        SFrame1 SFrame2 ... are files to be analyzed.


Description:  

    This routine reads an lvmSFrame file and constructs
    an html file, that contains information from the
    headers and various plots to indicate what the 
    quality of the data is

Primary routines:

    make_html is the primary driving routine
    steer handles the inputs

Notes:

    The html file is created in the current working directory
    and the various plots are in a subdirectory fig_quals.

    The links in the html file are relative to the html file
                                       
History:

240427 ksl Coding begun

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
from astropy.coordinates import get_body, solar_system_ephemeris, get_sun, AltAz, EarthLocation, get_moon
from astropy.time import Time
import astropy.units as u
from lvm_ksl import quick_map
from lvm_ksl import eval_standard


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
    med=np.median(f)
    zmin=ymin+med
    zmax=ymax+med
    print(zmin,zmax,med,ymin,ymax)
    return zmin,zmax



def eval_qual_sframe(filename='data/lvmSFrame-00011061.fits',ymin=-0.2e-13,ymax=1e-13,xmin=3600,xmax=9500,outroot=''):
    '''
    Provide a standard plot for looking at how well the sky subtraction has worked overall
    '''

    xtype='SkySubtracted'
    if filename.count('CFRAME'):
        print('This is an LCFrame')
        xtype='NotSkySubtracted'


    try:
        x=fits.open(filename)
    except:
        print('Error: eval_qual: Could not open %s' % filename)
        return

    hdr=x['PRIMARY'].header
    mjd=hdr['MJD']
    exposure=hdr['EXPOSURE']

    ra=get_header_value(hdr,'TESCIRA')
    dec=get_header_value(hdr,'TESCIDE')

    ra_sky_e=get_header_value(hdr,'POSKYERA')
    dec_sky_e=get_header_value(hdr,'POSKYEDE')

    ra_sky_w=get_header_value(hdr,'POSKYWRA')
    dec_sky_w=get_header_value(hdr,'POSKYWDE')


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

    fig=plt.figure(1,(8,12))
    plt.clf()
    gs= GridSpec(3, 3, figure=fig)

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

    ax2 = fig.add_subplot(gs[1, :])
    ax2.semilogy(wav,sci_flux_med+sci_sky_med,label='Science Total',zorder=2)
    ax2.semilogy(wav,sci_sky_med,label='Science Sky',zorder=1)
    ymin,ymax=plt.ylim()
    ymax=np.max(sci_flux_med+sci_sky_med)
    ax2.set_ylim(1e-3*ymax,1.1*ymax)
    ax2.set_xlim(3600,9600)
    ax2.legend()


    ax3 = fig.add_subplot(gs[2, 0])
    wmin=4650
    wmax=5100

    xwav,xsci_flux_med=limit_spectrum(wav,sci_flux_med,wmin,wmax)
    xmedian=np.median(xsci_flux_med)
    xsci_flux_med-=xmedian
    ax3.plot(xwav,xsci_flux_med,label='Sky-Subtracted Science',zorder=2)

    ax3.plot([wmin,wmax],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    ax3.plot([wmin,wmax],[-5.9e-15,-5.9e-15],':r')
    ax3.set_xlim(wmin,wmax)
    ymin,ymax=get_yscale(xsci_flux_med,-2e-14,2e-14)
    ax3.set_ylim(ymin,ymax)


    ax4 = fig.add_subplot(gs[2, 1])
    wmin=6250
    wmax=6800
    wmin=6500
    wmax=6800

    xwav,xsci_flux_med=limit_spectrum(wav,sci_flux_med,wmin,wmax)
    xmedian=np.median(xsci_flux_med)
    xsci_flux_med-=xmedian
    ax4.plot(xwav,xsci_flux_med,label=r'H$\alpha$/[NII]/[SII]',zorder=2)

    ax4.plot([wmin,wmax],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    ax4.plot([wmin,wmax],[-5.9e-15,-5.9e-15],':r')
    ax4.set_xlim(wmin,wmax)
    ymin,ymax=get_yscale(xsci_flux_med,-2e-14,2e-14)
    ax4.set_ylim(ymin,ymax)
    # ax4.legend()
    

    ax5 = fig.add_subplot(gs[2, 2])
    wmin=9450
    wmax=9600
    wmin=9480
    wmax=9586
    xwav,xsci_flux_med=limit_spectrum(wav,sci_flux_med,wmin,wmax)
    xmedian=np.median(xsci_flux_med)
    xsci_flux_med-=xmedian
    

    ax5.plot(xwav,xsci_flux_med,label='Sky-Subtracted Science',zorder=2)

    ax5.plot([wmin,wmax],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    ax5.plot([wmin,wmax],[-5.9e-15,-5.9e-15],':r')
    ax5.set_xlim(wmin,wmax)
    ymin,ymax=get_yscale(xsci_flux_med,-2e-14,2e-14)
    ax5.set_ylim(ymin,ymax)

    plt.tight_layout()

    location='./figs_qual/'

    if os.path.isdir(location)==False:
        os.mkdir(location)

    words=filename.split('/')
    root=words[-1].replace('.fits','')
    print(root)
    figname='%s/%s.png' % (location,root)
    plt.savefig(figname)

    # Now make another plot for the sky fibers

    fig=plt.figure(2,(8,12))
    plt.clf()
    gs= GridSpec(3, 3, figure=fig)

    ax1 = fig.add_subplot(gs[0, :])
    # ax1.plot(wav,sci_flux_med,label='Sky-Subtracted Science',zorder=2)
    ax1.plot(wav,skye_flux_med,label='SkyE-Subtracted SkyE',zorder=1)
    ax1.plot(wav,skyw_flux_med,label='SkyW-Subtracted SkyW',zorder=0)
    ax1.plot([3600,9600],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    ax1.plot([3600,9600],[-5.9e-15,-5.9e-15],':r')
    ax1.set_xlim(3600,9600)
    ymin,ymax=ax1.get_ylim()
    ax1.set_ylim(0.2*ymin,0.2*ymax)
    ax1.legend()

    ax2 = fig.add_subplot(gs[1, :])
    # ax2.semilogy(wav,sci_flux_med+sci_sky_med,label='Science Total',zorder=2)
    # ax2.semilogy(wav,sci_sky_med,label='Science Sky',zorder=1)
    delta=skyw_flux_med+skyw_sky_med-(skye_flux_med+skye_sky_med)
    if distance_sky_w<distance_sky_e:
        ax2.plot(wav,delta,label='SkyW-SkyE (Nearer-Further)',zorder=1)
    else:
        delta=-delta
        ax2.plot(wav,delta,label='SkyE-SkyW (Nearer-Further)',zorder=1)

    # ymin,ymax=plt.ylim()
    # ymax=np.max(sci_flux_med+sci_sky_med)
    ax2.set_ylim(1e-3*ymax,1.1*ymax)
    ax2.set_xlim(3600,9600)
    ax2.legend()


    ax3 = fig.add_subplot(gs[2, 0])
    wmin=4650
    wmax=5100
    wmin=4800
    wmax=5100



    xwav,delta_limit=limit_spectrum(wav,delta,wmin,wmax)
    delta_median=np.median(delta_limit)
    print('a',delta_median)
    delta_limit-=delta_median
    delta_median=np.median(delta_limit)
    print('b',delta_median)
    ax3.plot(xwav,delta_limit,label='SkyE-Subtracted SkyE',zorder=1)

    # ax3.plot(xwav,xskyw_flux_med,label='SkyW-Subtracted SkyW',zorder=0)

    ax3.plot([wmin,wmax],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    ax3.plot([wmin,wmax],[-5.9e-15,-5.9e-15],':r')

    ax3.set_xlim(wmin,wmax)
    ymin,ymax=get_yscale(delta_limit,-2e-14,2e-14)
    ax3.set_ylim(ymin,ymax)


    ax4 = fig.add_subplot(gs[2, 1])
    wmin=6250
    wmax=6800
    wmin=6500
    wmax=6800

    xwav,delta_limit=limit_spectrum(wav,delta,wmin,wmax)

    delta_median=np.median(delta_limit)
    print('a',delta_median)
    delta_limit-=delta_median
    delta_median=np.median(delta_limit)
    print('b',delta_median)

    ax4.plot(xwav,delta_limit,zorder=1)

    ax4.plot([wmin,wmax],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    ax4.plot([wmin,wmax],[-5.9e-15,-5.9e-15],':r')
    ax4.set_xlim(wmin,wmax)
    # ymin,ymax=get_yscale(xsci_flux_med,-1e-14,1e-14)
    ymin,ymax=get_yscale(delta_limit,-2e-14,2e-14)
    ax4.set_ylim(ymin,ymax)

    ax5 = fig.add_subplot(gs[2, 2])
    wmin=9450
    wmax=9600
    wmin=9480
    wmax=9586
    xwav,xsci_flux_med=limit_spectrum(wav,sci_flux_med,wmin,wmax)
    xwav,delta_limit=limit_spectrum(wav,delta,wmin,wmax)
    delta_median=np.median(delta_limit)
    delta_limit-=delta_median
    delta_median=np.median(delta_limit)


    ax5.plot(xwav,delta_limit,label='SkyE-Subtracted SkyE',zorder=1)
    ax5.plot([wmin,wmax],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    ax5.plot([wmin,wmax],[-5.9e-15,-5.9e-15],':r')
    ax5.set_xlim(wmin,wmax)
    ymin,ymax=get_yscale(delta_limit,-2e-14,2e-14)
    ax5.set_ylim(ymin,ymax)

    plt.tight_layout()

    words=filename.split('/')
    root=words[-1].replace('.fits','')
    print(root)
    sky_figname='%s/%s_sky.png' % (location,root)
    plt.savefig(sky_figname)




    return figname,sky_figname
                 


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
            value=default_value
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
        value = default_value
    return value



def create_overview(filename='data/lvmSFrame-00011061.fits'):
    ''' Sumarize information about the  processed data file
    '''
    try:
        x=fits.open(filename)
    except:
        print('Error: Could not open %s' % filename)
        return

    hdr=x['PRIMARY'].header

    exposure=get_header_value(hdr,'EXPOSURE')
    mjd=get_header_value(hdr,'MJD')
    object_name=get_header_string(hdr,'OBJECT')
    obs_time=get_header_string(hdr,'OBSTIME')
    ra=get_header_value(hdr,'TESCIRA')
    dec=get_header_value(hdr,'TESCIDE')
    try:
        pa=get_header_value(hdr,'POSCIPA')
    except:
        print('Error: PSCIIPA keyword is missing, assuming PA is zero')
        pa=0

    ra_sky_e=get_header_value(hdr,'POSKYERA')
    dec_sky_e=get_header_value(hdr,'POSKYEDE')
    pa_sky_e=get_header_value(hdr,'POSKYEPA')

    ra_sky_w=get_header_value(hdr,'POSKYWRA')
    dec_sky_w=get_header_value(hdr,'POSKYWDE')
    pa_sky_w=get_header_value(hdr,'POSKYWPA')


    distance_sky_w=distance(ra,dec,ra_sky_w,dec_sky_w)
    distance_sky_e=distance(ra,dec,ra_sky_e,dec_sky_e)

  
    # print(obstime)

    moon_info=get_moon_info_las_campanas(obs_time)
    #for key, value in moon_info.items():
    #    print(f'{key}: {value}')

   

    xlist=[]
    xlist.append('Exposure : %d' % exposure)
    xlist.append('MJD      : %d' % mjd)
    xlist.append('Obs. time: %s' % obs_time)
    xlist.append('Object.  : %s' % object_name)
    xlist.append('Science RA  Dec. PA : %8.2f %8.2f %8.2f' % (ra,dec,pa))
    xlist.append('SkyE    RA  Dec. PA (ang distance): %8.2f %8.2f %8.2f (%8.2f)' % (ra_sky_e,dec_sky_e,pa_sky_e,distance_sky_e))
    xlist.append('SkyW    RA  Dec. PA (ang distance): %8.2f %8.2f %8.2f (%8.2f)' % (ra_sky_w,dec_sky_w,pa_sky_w,distance_sky_w))
    xlist.append('Moon    RA  Dec. Alt.  Ill:  %8.2f %8.2f %8.2f %8.2f' % 
                 (moon_info['MoonRA'],moon_info['MoonDec'],moon_info['MoonAlt'],moon_info['MoonIll']))
    xlist.append('Sun.    RA  Dec. Alt.:  %8.2f %8.2f %8.2f' %
                 (moon_info['SunRA'],moon_info['SunDec'],moon_info['SunAlt']))

 
    return xlist



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
    
    # Create a grid of pixel coordinates
    y, x = np.indices(data.shape)
    
    # Convert corner pixel coordinates to celestial coordinates
    ra_min, dec_min = wcs.all_pix2world(0, 0, 0)
    ra_max, dec_max = wcs.all_pix2world(data.shape[1], data.shape[0], 0)
    
    # Plot the image
    plt.figure(figsize=(8, 8))
    cmap=plt.get_cmap('hot')
    cmap.set_bad(color='gray', alpha=0.2)
    plt.imshow(data, cmap=cmap, vmin=min_val, vmax=max_val, origin='lower', extent=(0, data.shape[1], 0, data.shape[0]))
    plt.colorbar(label='Intensity', shrink=0.8)  # Adjust the shrink parameter as needed
    
    # Set up RA and DEC axis labels and ticks
    plt.xlabel('RA')
    plt.ylabel('DEC')
    
    # Set RA tick positions and labels
    ra_ticks = np.linspace(ra_min, ra_max, 5)
    ra_tick_labels = [f'{ra:.2f}' for ra in ra_ticks]
    plt.xticks(np.linspace(0, data.shape[1], 5), ra_tick_labels)
    
    # Set DEC tick positions and labels
    dec_ticks = np.linspace(dec_min, dec_max, 5)
    dec_tick_labels = [f'{dec:.2f}' for dec in dec_ticks]
    plt.yticks(np.linspace(0, data.shape[0], 5), dec_tick_labels)
    
    plt.title(title)
    plt.grid(color='white', ls='dotted')
    

    plt.savefig(outname)




def make_images(filename='data/llvmSFrame-00011061.fits',outroot='test'):
    '''
    Make multiple images of the CFrame of SFrame data
    '''
    xha=['Ha',[6560.,6566.],[6590.,6630.]]
    xs2=['[SII]',[6710.,6735.],[6740.,6760.]]
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
    plot_fits_image(filename=s2_file,title='[SII]',outname=s2_plot)
    return ha_plot,s2_plot,cont_plot

science_plot_comment='''
The median sky subtracted spectrum from the science fibers.  The top panel shows the median spectrum.
The middle panel shows the sum of the flux and sky, and just the sky.  The three panels at the 
bottom show the median scence spectra (after a crude contiumm subtraction) in three wavelength ranges, corresponding to Hbeta-[0II], Halpha-[SII], and [SIII]9071.
Several of the plots also have dashed lines which outline the 5 sigma 
sensitivity limit in the Milky Way.
'''

sky_plot_comment='''
Comparisons of the median spectra in the SkyE and SkyW telescopes.  The top panel shows the median spectrum 
in each of the two sky telescopes (after sky subtraction.) The middle panel shows the difference in the the 
total fluxes (with sky included) in the two sky telescopes. The difference is computed by subtracting the 
spectrum of the sky telescope that is furthers from the science target to from the spectrum of the nearer sky telescope.
The bottom panel shows the differences in sky subtracted spectra in three spectral regions. (Note that at present, 
the lvmdrp uses the sky calculated for the science telescope 
for subtracting sky from the sky telescopes. This implies that what is presented in this figure tells one 
mostly about the differences in the sky in the two telescopes.)
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
    print(outroot)

    string=xhtml.begin('LVMDRP Quality Asssessment for %s' % filename)
    string+=xhtml.hline()
    
    overview_list=create_overview(filename)
    string+=xhtml.add_list(overview_list)



    string+=xhtml.hline()
    string+=xhtml.h2('Science Spectrum')
    string+=xhtml.paragraph(science_plot_comment)
    

    figname,sky_figname= eval_qual_sframe(filename,ymin=-0.2e-13,ymax=1e-13,xmin=3600,xmax=9500)

    string+=xhtml.image('file:%s' % (figname),width=900,height=1200)
    string+=xhtml.hline()
    string+=xhtml.h2('SkyE and SkyW  Spectra')
    string+=xhtml.paragraph(sky_plot_comment)

    string+=xhtml.image('file:%s' % (sky_figname),width=900,height=1200)
    string+=xhtml.hline()
    string+=xhtml.h2('Line and Continuum images')

    ha_plot,s2_plot,cont_plot=make_images(filename,outroot)

    string+=xhtml.paragraph(image_comment)

    string+=xhtml.image('file:%s' % (ha_plot),width=900,height=900)
    string+=xhtml.image('file:%s' % (s2_plot),width=900,height=900)
    string+=xhtml.image('file:%s' % (cont_plot),width=900,height=900)

    string+=xhtml.hline()

    string+=xhtml.paragraph('Comparision between the flux calibrated star fibers to the GAIA spectra of the stars')

    outname='figs_qual/standard_%s.png' % outroot
    status=eval_standard.qual_eval(filename,outname)

    if status==True:
        string+=xhtml.image(outname,width=900,height=900)
    else:
        string+=xhtml.paragraph('Could not do standard standard star comparision')
    
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
            print(__doc__)
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
