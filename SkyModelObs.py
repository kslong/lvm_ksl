#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Generate a prediction of the sky using the ESO Sky Model
given an RA and Dec and a time


Command line usage (if any):

    usage: SkyModelObs.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
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

def setup(eso_sky_dir=''):
    '''
    Make sure the directies we need exist
    '''

    xdir=os.getenv('ESO_SKY_MODEL')
    if eso_sky_dir=='':
        eso_sky_dir=xdir

    data_dir='%s/sm-01_mod2/data' % eso_sky_dir
    output_dir='%s/sm-01_mod1/output' % eso_sky_dir
    config_dir='%s/sm-01_mod2/config' % eso_sky_dir
    if os.path.exists(data_dir) == False and os.path.islink(data_dir)==False:
        os.symlink(data_dir, 'data')
    if os.path.exists(output_dir) ==False and os.path.islink(output_dir)==False :
        os.symlink(output_dir, 'output')
    if os.path.exists(config_dir)==False  and os.path.islink(config_dir)==False :
        os.symlink(config_dir, 'config')




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
    
    # Calculate ecliptic coordinates
    moon_ecliptic = moon_coords.transform_to(GeocentricTrueEcliptic(equinox=obs_time))
    sun_ecliptic = sun_coords.transform_to(GeocentricTrueEcliptic(equinox=obs_time))
    source_ecliptic = source_coords.transform_to(GeocentricTrueEcliptic(equinox=obs_time))
    
    # Convert ecliptic longitudes to -180 to 180 range
    moon_eclip_lon = moon_ecliptic.lon.deg
    if moon_eclip_lon > 180:
        moon_eclip_lon -= 360
        
    sun_eclip_lon = sun_ecliptic.lon.deg
    if sun_eclip_lon > 180:
        sun_eclip_lon -= 360
        
    source_eclip_lon = source_ecliptic.lon.deg
    if source_eclip_lon > 180:
        source_eclip_lon -= 360
        
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
season   = 4

# period of the night (x/3 of night, x = 1,2,3; 0: entire night)
time     = 0

# vac[uum] or air wavelengths
vac_air  = air

# precipitable water vapour in mm (-1: bimonthly mean)
pwv      = -1

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




def create_inputs(ra=296.242608,dec=-14.811007,obstime='2023-08-29T03:20:43.668'):
    '''
    calcskymodel reads inputs for the actual source location etc from a fixed file
    called sky_model_etc.par
    '''
    info=get_info_las_campanas(obstime, ra=ra,dec=dec,verbose=True)
    print(info)

    xout=open('config/skymodel_etc.par','w')
    out_string=obs_base % (info['SourceAlt'],info['Moon-Sun_Separation'],info['Moon-Source_Separation'],info['MoonAlt'],info['MoonDistanceInMeanUnits'],
                           info['SourceEclipLon'],info['SourceEclipLat'])
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
    ftab=join(rtab,ttab,join_type='left')
    ftab['lam']*=1000.
    
    primary_hdu = fits.PrimaryHDU(header=rad[0].header)
    table_hdu=table_hdu = fits.BinTableHDU(data=ftab, header=rad[1].header)
    new_hdul = fits.HDUList([primary_hdu, table_hdu])
    i=0
    while i<len(xkey):
        new_hdul['PRIMARY'].header[xkey[i]]=xval[i]
        i+=1
    new_hdul.writeto(outfile, overwrite=True)


def do_one(ra=296.242608,dec=-14.811007,obstime='2023-08-29T03:20:43.668',outroot=''):

    setup()

    key,value=create_inputs(ra=ra,dec=dec,obstime=obstime)

    result=subprocess.run(['/Users/long/SDSS/skymodel/sm-01_mod2/bin/calcskymodel'],capture_output=True,text=True)

    print("stdout:", result.stdout)
    print("stderr:", result.stderr)

    if outroot=='':
        outroot='SkyM_%5.1f_%5.1f' % (ra,dec)
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

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i][0]=='-' and is_number(argv[i])==False:
            print('Error unknown option: ',argv)
            return
        elif ra < 0:
            ra=eval(argv[i])
        elif dec <-90:
            dec=eval(argv[i])
        elif xtime=='':
            xtime=argv[i]
        else:
            print('Erorr Too many argments: ',argv)
            return


        i+=1

    do_one(ra=ra,dec=dec,obstime=xtime,outname='')



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)





