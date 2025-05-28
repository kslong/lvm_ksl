#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Generate a prediction of the sky using the ESO Sky Model
given an RA and Dec and a time


Command line usage (if any):

    usage: SkyModelObs.py [-h] [-config] [-data data_dir] [-out whatever] ra dec time

    where
        -h prints the _doc__ string and exits
        -config forces a new set up of the directories needed to 
            run the routines. (Regardless of this switch the directories will be set up the
            first time the routine is run.)
        -data data_dir  Sets the data directory for sky model to a location
            not pointed to by ESO_SKY_MODEL
        -out whatever changes the root name of the output file
        ra, dec  a position in the sky in degrees
        time time in one of serveral coordiantes, either a date_time string, mjd, or jd

Description:  

    The routine uses the routine calcskymodel which is part of the ESO Sky model to
    estimate the sky from LCO at a certain time and position on the sky.  If reconstructs
    the outputs to resemble the data files created by skycalc, the web application.



Primary routines:

    doit

Notes:

    There must be conistency between the various inputs, and in particular one must have
    the agreement between the calculated sky models in the data/lib directory and the 
    parameters described in confit/sky_model.par
                                       
History:

250510 ksl Coding begun

'''

# # Generate inputs for running the sky model
# 
# Also attempt to actually run model, and then reconstuct a fits file that resembles SkyCalc


import os
import sys
from datetime import date

import subprocess
from astropy.io import ascii,fits
import matplotlib.pyplot as plt


from astropy.table import Table, join
from astropy.wcs import WCS
from astropy.coordinates import get_body, solar_system_ephemeris, get_sun, AltAz, EarthLocation
from astropy.coordinates import SkyCoord, Distance
from astropy.coordinates import GeocentricTrueEcliptic
from astropy.time import Time
import astropy.units as u
import numpy as np
from astropy.coordinates import get_body

import shutil


import requests
import inspect



# This section is a fairly flexiable way to convert times
# The main routine is convert_time

from astropy.time import Time
import datetime
import re

def convert_time(time_input, output_format='datetime'):
    """
    Convert between different time formats: date strings, MJD, JD, and datetime objects.

    Parameters
    ----------
    time_input : str, float, or int
        Time input in various formats:
        - ISO date string (e.g., '2025-04-21T03:00:00')
        - MJD (e.g., 60394.125 or '60394.125')
        - JD (e.g., 2460394.625 or '2460394.625')
        - datetime object
    output_format : str
        Desired output format:
        - 'datetime': Python datetime object
        - 'iso': ISO 8601 format string (e.g., '2025-04-21T03:00:00')
        - 'iso_ms': ISO 8601 format with milliseconds (e.g., '2023-08-29T03:20:43.668')
        - 'mjd': Modified Julian Date (float)
        - 'jd': Julian Date (float)

    Returns
    -------
    Time in the requested format
    """
    # Detect input format
    input_format = _detect_time_format(time_input)

    # Convert to Time object (intermediate representation)
    t = _convert_to_time_object(time_input, input_format)

    # Convert to desired output format
    return _convert_to_output_format(t, output_format)

def _detect_time_format(time_input):
    """
    Detect the format of the input time.

    Parameters
    ----------
    time_input : str, float, int, or datetime.datetime
        Time input in various formats

    Returns
    -------
    str
        Detected format: 'iso', 'mjd', 'jd', or 'datetime'
    """
    # If it's already a datetime object
    if isinstance(time_input, datetime.datetime):
        return 'datetime'

    # Convert to string for pattern matching
    if isinstance(time_input, (int, float)):
        time_str = str(time_input)
    else:
        time_str = time_input

    # Check for ISO format (contains date separators and possibly time separators)
    # Enhanced to detect milliseconds/microseconds patterns
    if (re.search(r'\d{4}-\d{2}-\d{2}', time_str) or 
        re.search(r'\d{4}/\d{2}/\d{2}', time_str) or
        re.search(r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.?\d*', time_str)):
        return 'iso'

    # Try to convert to float
    try:
        value = float(time_str)
        # JD is typically > 2400000
        if value > 2400000:
            return 'jd'
        # MJD is typically < 100000
        else:
            return 'mjd'
    except ValueError:
        # If it's not convertible to a number but has date-like patterns, treat as ISO
        if re.search(r'\d+[/:\-]\d+', time_str):
            return 'iso'

    # Default case - try as ISO
    return 'iso'

def _convert_to_time_object(time_input, input_format):
    """
    Convert input to an astropy Time object.

    Parameters
    ----------
    time_input : str, float, int, or datetime.datetime
        Time input
    input_format : str
        Format of the input ('iso', 'mjd', 'jd', or 'datetime')

    Returns
    -------
    astropy.time.Time
        Time object representing the input
    """
    try:
        if input_format == 'datetime':
            return Time(time_input, scale='utc')
        elif input_format == 'iso':
            # Handle various ISO formats including those with milliseconds/microseconds
            return Time(time_input, format='isot', scale='utc')
        elif input_format == 'mjd':
            return Time(float(time_input), format='mjd', scale='utc')
        elif input_format == 'jd':
            return Time(float(time_input), format='jd', scale='utc')
    except Exception as e:
        raise ValueError(f"Failed to convert {time_input!r} as {input_format} format: {str(e)}")

def _convert_to_output_format(time_obj, output_format):
    """
    Convert Time object to the desired output format.

    Parameters
    ----------
    time_obj : astropy.time.Time
        Time object
    output_format : str
        Desired output format ('datetime', 'iso', 'iso_ms', 'mjd', or 'jd')

    Returns
    -------
    The time in the requested format
    """
    try:
        if output_format == 'datetime':
            return time_obj.to_datetime()
        elif output_format == 'iso':
            return time_obj.iso
        elif output_format == 'iso_ms':
            # Format with millisecond precision: YYYY-MM-DDTHH:MM:SS.sss
            dt = time_obj.to_datetime()
            # Format with 3 decimal places for milliseconds
            return dt.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]
        elif output_format == 'mjd':
            return time_obj.mjd
        elif output_format == 'jd':
            return time_obj.jd
        else:
            raise ValueError(f"Unsupported output format: {output_format}")
    except Exception as e:
        raise ValueError(f"Failed to convert to {output_format} format: {str(e)}")

# End of routines to convert time


def safe_remove(path):
    if os.path.islink(path):
        os.unlink(path)  # ✅ Just remove the symlink
    elif os.path.isdir(path):
        shutil.rmtree(path)  # ✅ Remove actual directory (even if not empty)
    else:
        return

# safe_remove("your_path_here")

def setup(eso_sky_dir='',config=True):
    '''
    Set up the directories that are needed for the routine to run
    '''

    if config==False:
        icheck=True
        # We just check to see if the directories exist
        if os.path.isdir('config')==False:
            icheck=False
        if os.path.isdir('output')==False:
            icheck=False
        if os.path.isdir('data')==False and os.path.islink('data')==False :
            icheck=False
        if icheck==True:
            return



    xdir=os.getenv('ESO_SKY_MODEL')
    if eso_sky_dir=='':
        eso_sky_dir=xdir

    data_dir='%s/sm-01_mod2/data' % eso_sky_dir
    if os.path.isdir(data_dir)==False:
        print('Error: %s does not appear to exist')
        return 
    safe_remove('data')
    os.symlink(data_dir,'data')
    os.makedirs('output', exist_ok=True)
    os.makedirs('config', exist_ok=True)


#   output_dir='%s/sm-01_mod1/output' % eso_sky_dir
#   config_dir='%s/sm-01_mod2/config' % eso_sky_dir
#   if os.path.exists(data_dir) == False and os.path.islink(data_dir)==False:
#       os.symlink(data_dir, 'data')
#   if os.path.exists(output_dir) ==False and os.path.islink(output_dir)==False :
#       os.symlink(output_dir, 'output')
#   if os.path.exists(config_dir)==False  and os.path.islink(config_dir)==False :
#       os.symlink(config_dir, 'config')
    return



def get_info_las_campanas(datetime_utc, ra, dec, verbose=False):
    '''
    Get information about the sun, moon, and a source at given RA and Dec as a function of UT
    
    Parameters:
    -----------
    datetime_utc : str or datetime
        UTC time for the observation
    ra : float
        Right ascension of the source in degrees
    dec : float
        Declination of the source in degrees
    verbose : bool, optional
        If True, print the information
        
    Returns:
    --------
    dict
        Dictionary containing information about the sun, moon, and source
    '''

    if verbose:
        print('get_info_las_campanas,Start: ',datetime_utc, ra, dec)
    # Las Campanas Observatory coordinates
    observatory_location = EarthLocation(lat=-29.0089*u.deg, lon=-70.6920*u.deg, height=2281*u.m)
    
    # Specify the observation time in UT
    obs_time = Time(datetime_utc)
    
    # Create source coordinates
    source_coords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    
    # Set the solar system ephemeris to 'builtin' for faster computation
    with solar_system_ephemeris.set('builtin'):
        # Get the Moon's and Sun's coordinates at the specified time
        moon_coords = get_body('moon', obs_time, location=observatory_location)
        sun_coords = get_body('sun', obs_time, location=observatory_location)
    
    # Calculate the phase angle (angle between the Sun, Moon, and observer)
    phase_angle = moon_coords.separation(sun_coords).radian
    
    # Calculate the illuminated fraction of the Moon
    illumination_fraction = (1 - np.cos(phase_angle))/2
    
    moon_sun_longitude_diff = (moon_coords.ra - sun_coords.ra).wrap_at(360 * u.deg).value
    if moon_sun_longitude_diff > 0:
        moon_phase = illumination_fraction/2.
    else:
        moon_phase = 1-illumination_fraction/2.
    illumination_fraction *= 100.
    
    # Calculate the Altitude and Azimuth from Las Campanas Observatory
    altaz_frame = AltAz(obstime=obs_time, location=observatory_location)
    moon_altaz = moon_coords.transform_to(altaz_frame)
    sun_altaz = sun_coords.transform_to(altaz_frame)
    source_altaz = source_coords.transform_to(altaz_frame)
    if source_altaz.alt.deg<0:
        print('Error: Source altitude is negative :', source_altaz.alt.deg,ra,dec,datetime_utc)
    
    # Calculate ecliptic coordinates
    moon_ecliptic = moon_coords.transform_to(GeocentricTrueEcliptic(equinox=obs_time))
    sun_ecliptic = sun_coords.transform_to(GeocentricTrueEcliptic(equinox=obs_time))
    source_ecliptic = source_coords.transform_to(GeocentricTrueEcliptic(equinox=obs_time))
    
    # Convert ecliptic longitudes to -180 to 180 range
    moon_eclip_lon = moon_ecliptic.lon.deg
    if moon_eclip_lon > 180:
        moon_eclip_lon -= 360

    if verbose:
        print('XXX %.1f %.1f -> %.1f ' % (source_ecliptic.lon.deg,sun_ecliptic.lon.deg,source_ecliptic.lon.deg-sun_ecliptic.lon.deg))
        
    sun_eclip_lon = sun_ecliptic.lon.deg
    # if sun_eclip_lon > 180:
    #     sun_eclip_lon -= 360
        
    source_eclip_lon = source_ecliptic.lon.deg
    # if source_eclip_lon > 180:
    #   source_eclip_lon -= 360

        
    # Calculate Moon-Earth distance in units of mean distance
    # Mean Earth-Moon distance is 384,400 km
    mean_moon_distance = 384400 * u.km
    moon_distance = moon_coords.distance.to(u.km)
    moon_distance_in_mean = moon_distance / mean_moon_distance
    
    # Calculate separations between objects
    moon_sun_separation = moon_coords.separation(sun_coords).deg
    moon_source_separation = moon_coords.separation(source_coords).deg
    sun_source_separation = sun_coords.separation(source_coords).deg
    
    # Create result dictionary
    xreturn = {
        'SunRA': sun_coords.ra.deg,
        'SunDec': sun_coords.dec.deg,
        'SunAlt': sun_altaz.alt.deg,
        'SunAz': sun_altaz.az.deg,
        'SunEclipLon': sun_eclip_lon,
        'SunEclipLat': sun_ecliptic.lat.deg,
        'MoonRA': moon_coords.ra.deg,
        'MoonDec': moon_coords.dec.deg,
        'MoonAlt': moon_altaz.alt.deg,
        'MoonAz': moon_altaz.az.deg,
        'MoonEclipLon': moon_eclip_lon,
        'MoonEclipLat': moon_ecliptic.lat.deg,
        'MoonPhas': moon_phase,
        'MoonIll': illumination_fraction,
        'MoonDistance': moon_distance.value,
        'MoonDistanceInMeanUnits': moon_distance_in_mean.value,
        'SourceRA': source_coords.ra.deg,
        'SourceDec': source_coords.dec.deg,
        'SourceAlt': source_altaz.alt.deg,
        'SourceAz': source_altaz.az.deg,
        'SourceEclipLon': source_eclip_lon,
        'SourceEclipLat': source_ecliptic.lat.deg,
        'Moon-Sun_Separation': moon_sun_separation,
        'Moon-Source_Separation': moon_source_separation,
        'Sun-Source_Separation': sun_source_separation
    }
    
    if verbose:
        for key, value in xreturn.items():
            print(f'{key}: {value}')
    
    # Return the information
    return xreturn


inst_base='''
# Wavelength grid:

# minimum and maximum wavelength [mum]
limlam     = 0.36 0.98

# step size [mum]
dlam       = 0.00005


# Line-spread function:

# radius of convolution kernel [pixels] (N_pixel = 2 x kernrad + 1)
kernrad    = 3

# FWHM of boxcar kernel [pixels]
wbox       = 0.8

# FWHM of Gaussian kernel [pixels]
wgauss     = 0.8

# FWHM of Lorentzian kernel [pixels]
wlorentz   = 0.8

# variable kernel (width proportional to wavelength)? -> 1 = yes; 0 = no
# if varkern = 1: kernel radius and FWHM for central wavelength
varkern    = 1

# output file for kernel ("stdout": screen; "null": no output)
kernelfile = output/kernel.dat
'''


obs_base='''
# observatory height in km [2.4, 3.06] (default: 2.64)
# sm_h = 2.64
sm_h = 2.5

# lower height limit in km (default: 2.0)
sm_hmin = 2.0

# altitude of object above horizon [0,90]
alt      = %.1f

# separation of Sun and Moon as seen from Earth [0,360]
# (> 180 for waning Moon)
alpha    = %.1f

# separation of Moon and object [0,180]
rho      = %.1f

# altitude of Moon above horizon [-90,90]
altmoon  = %.1f

# distance to Moon (mean distance = 1; [0.91,1.08])
moondist = %.2f

# pressure at observer altitude in hPa (default: 744)
pres     = 744.

# single scattering albedo for aerosols [0,1] (default: 0.97)
ssa      = 0.97

# calculation of double scattering of moonlight ('Y' or 'N')
calcds   = N

# relative UV/optical ozone column density (1 -> 258 DU)
o3column = 1.

# scaling factor for scattered moonlight (default: 1.0)
moonscal = 1.0

# heliocentric ecliptic longitude of object [-180,180]
lon_ecl  = %.1f

# ecliptic latitude of object [-90,90]
lat_ecl  = %.1f

# grey-body emissivity (comma-separated list)
emis_str = 0.2

# grey-body temperature in K (comma-separated list)
temp_str = 290.

# monthly-averaged solar radio flux [sfu]
msolflux = 101.

# bimonthly period (1: Dec/Jan, ..., 6: Oct/Nov; 0: entire year)
season   = 0

# period of the night (x/3 of night, x = 1,2,3; 0: entire night)
time     = 0 

# vac[uum] or air wavelengths
vac_air  = air

# precipitable water vapour in mm (-1: bimonthly mean)
pwv      = 3.5

# radiative transfer code L(BLRTM) or R(FM) for molecular spectra
rtcode   = L

# resolution of molecular spectra in library (crucial for run time)
# resol    = 1e6
resol    = 6e4

# path to file sm_filenames.dat for data paths and file names
filepath = data

# inclusion of sky model components
# format: "xxxxxxx" where x = "Y" (yes) or x = "N" (no)
# pos. 1: scattered moonlight
#      2: scattered starlight
#      3: zodiacal light
#      4: thermal emission by telescope/instrument
#      5: molecular emission of lower atmosphere
#      6: sky emission lines of upper atmosphere
#      7: airglow con
incl     = YYYYYYY
'''




def create_inputs(ra=296.242608,dec=-14.811007,obstime='2023-08-29T03:20:43.668',verbose=False):
    '''
    calcskymodel reads inputs for the actual source location etc from a fixed file
    called sky_model_etc.par
    '''
    info=get_info_las_campanas(obstime, ra=ra,dec=dec,verbose=verbose)
    # print(info)

    longitude=info['SourceEclipLon']-info['SunEclipLon']
    # longitude must be between -180 and 180
    if longitude>180:
        longitude=longitude-360.
    if longitude< -180.:
        longitude=longitude+360.


    # print('ZZZ ',info['SourceEclipLon']-info['SunEclipLon'],longitude)


    xout=open('config/skymodel_etc.par','w')
    out_string=obs_base % (info['SourceAlt'],info['Moon-Sun_Separation'],info['Moon-Source_Separation'],info['MoonAlt'],info['MoonDistanceInMeanUnits'],
                           longitude,info['SourceEclipLat'])
    xout.write(out_string)
    xout.close()

    keys=['RA','Dec','ObsTime']
    values=[ra,dec,obstime]
    
    lines=out_string.split('\n')
    for one_line in lines:
        word=one_line.split()
        if len(word)>2 and word[1]=='=':
            keys.append(word[0])
            try:
                v=eval(word[2])
            except:
                v=word[2]
            values.append(v)

    xinst=open('config/instrument_etc.par','w')
    ibase=inst_base
    xinst.write(ibase)
    xinst.close()
    
    lines=ibase.split('\n')
    for one_line in lines:
        word=one_line.split()
        if len(word)>2 and word[1]=='=':
            keys.append(word[0])
            try:
                v=eval(word[2])
            except:
                v=word[2]
            values.append(v)


    return keys,values
        


def reformat_model(rfile='output/radspec.fits',tfile='output/transspec.fits',xkey=[],xval=[],outfile='foo.fits'):
    '''
    Reformat the outputs of the sky model into something resembling the SkyCalc fits file
    '''
    rad=fits.open(rfile)
    trans=fits.open(tfile)
    rtab=Table(rad[1].data)
    ttab=Table(trans[1].data)
    ztab=join(rtab,ttab,join_type='left')
    ztab['lam']*=1000.

    ztab.rename_column('lam','WAVE')
    ztab.rename_column('flux','FLUX')
    ztab.rename_column('flux_sml','MOON')
    ztab.rename_column('flux_zl','ZODI')
    ztab.rename_column('flux_ael','LINES')
    ztab.rename_column('flux_arc','DIFFUSE')
    area=np.pi*(37/2)**2

    q=1.98644586e-17*area/ztab['WAVE']
    ztab['FLUX']*=q
    ztab['MOON']*=q
    ztab['LINES']*=q
    ztab['ZODI']*=q
    ztab['DIFFUSE']*=q
    ztab['WAVE']*=10.
    ztab['CONT']=ztab['MOON']+ztab['ZODI']+ztab['DIFFUSE']
    ztab['CONT']/=ztab['trans']
    ztab['MOON']/=ztab['trans']
    ztab['ZODI']/=ztab['trans']
    ztab['DIFFUSE']/=ztab['trans']
    ztab['DIFFUSE']/=ztab['trans']
    new_hdu = fits.BinTableHDU(data=ztab)
    
    primary_hdu = fits.PrimaryHDU(header=rad[0].header)
    table_hdu=table_hdu = fits.BinTableHDU(data=ztab, header=rad[1].header)
    new_hdul = fits.HDUList([primary_hdu, table_hdu])
    i=0
    while i<len(xkey):
        new_hdul['PRIMARY'].header[xkey[i]]=xval[i]
        i+=1
    new_hdul.writeto(outfile, overwrite=True)


def do_one(ra=296.242608,dec=-14.811007,obstime='2023-08-29T03:20:43.668',xdata='',config=False,outroot=''):

    obstime=convert_time(obstime,'iso_ms')
    # print(obstime)

    setup(xdata,config)
    # print('Got Here')

    key,value=create_inputs(ra=ra,dec=dec,obstime=obstime)

    result=subprocess.run(['/Users/long/SDSS/skymodel/sm-01_mod2/bin/calcskymodel'],capture_output=True,text=True)

    # print("stdout:", result.stdout)
    if len(result.stderr):
        print("stderr:", result.stderr)
        print('Could not create model due to errors')
        return

    if outroot=='':
        mjd=convert_time(obstime,'mjd')
        print(mjd)
        if dec>0:
            outroot='SkyM_%8.2f_%05.1f_+%04.1f' % (mjd,ra,dec)
        else:
            outroot='SkyM_%8.2f_%05.1f_%.1f' % (mjd,ra,dec)
    if outroot.count('.fits')==0:
        outname='%s.fits' % outroot
    else:
        outname=outroot
    reformat_model(rfile='output/radspec.fits',tfile='output/transspec.fits',xkey=key,xval=value,outfile=outname)
    return outroot


def is_number(x):
    try:
        eval(x)
        return True
    except:
        return False


def steer(argv):
    '''
    SkyModelObs.py [-h]  RA Dec time
    '''

    ra=-99
    dec=-99.
    xtime=''
    data_dir=''
    config=False
    outroot=''

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-config':
            config=True    
        elif argv[i]=='-data':
            i+=1
            data_dir=argv[i]
        elif argv[i][:4]=='-out':
            i+=1
            outroot=argv[i]
        elif argv[i][0]=='-' and is_number(argv[i])==False:
            print('Error: unknown option: ',argv)
            return
        elif ra < 0:
            ra=eval(argv[i])
        elif dec <-90:
            dec=eval(argv[i])
        elif xtime=='':
            xtime=argv[i]
        else:
            print('Error: Too many argments: ',argv)
            return


        i+=1

    print(config)
    do_one(ra=ra,dec=dec,obstime=xtime,xdata=data_dir,config=config,outroot=outroot)



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)





