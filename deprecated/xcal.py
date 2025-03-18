#!/usr/bin/env python 
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Add RA and Dec  and moon/sun location 
information to calibrated data from LVM


Command line usage (if any):

    usage: xcal.py filename

    where one or more calibrated filenames can be provided



Description:  

    The rouine looks in a number of locations for the 
    caliibrated files.  In general if the file is
    in the directory on in standard data structure
    for LVM the file will be found.  

    Normally one will not provide the full path
    name to the file, unless one wants the output
    file to be written to the directory. Without
    a directory name, the output file will
    be written to the current working directory.

    The routine looks in the standard place for the
    guider information and adds that infomation to
    the header, to allow kslmap to be indepencedn
    of still having the guider data.

    The routine calculatates the postions of the fibers
    in RA and Dec and adds the positions to the SLITMAP.
    The routine also uses the UTC to calculate information
    about the sun and moon (including the illumination)
    which are added to the header).




Primary routines:

    doit

Notes:

    The routine makes use of routines in fib2radec
                                       
History:

231226 ksl Coding begun

'''

import sys

from astropy.coordinates import Angle
import numpy as np
from astropy.io import fits
from astropy.table import Table
import fib2radec
from astropy.coordinates import get_body, solar_system_ephemeris, get_sun, AltAz, EarthLocation, get_moon
from astropy.time import Time
import astropy.units as u
import numpy as np




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



def doit(filename='lvmCFrame-00006162.fits',outfile=''):
    '''
    Add fiber postions, best astometry, and fiber positions
    to calibrated science data
    '''
    xcal=fib2radec.locate_file(filename)
    if xcal==None:
        print('Could not locate calibrated file: %s' % filename)
        return
    else:
        print('Found %s' % xcal)

    x=fits.open(xcal, mode='update')
    exposure=x[0].header['Exposure']
    obstime=x['PRIMARY'].header['Obstime']
    # print(obstime)

    moon_info=get_moon_info_las_campanas(obstime)

    for key, value in moon_info.items():
        # print(f'{key}: {value}')
        x['PRIMARY'].header[f'{key}']=value
    
    guider_file='lvm.sci.coadd_s%08d.fits' % exposure

    xguide=fib2radec.locate_file(guider_file)
    if xguide==None:
        print('Could not locate guider file: %s' % guider_file)
        xguide=''    
    else:
        print('Found guider file %s' % xguide)

    xguide,racen,deccen,pa,xtab=fib2radec.get_ra_dec(xcal,xguide,outname='Fib2radec.txt')
    # print(xguide,racen,deccen,pa)

    if xguide:
        print('Using guider information for astometry')
        x['PRIMARY'].header['Best']='FromGuiTel'
    else:
        print('Using science  image for astometry')
        x['PRIMARY'].header['Best']='FromSciTel'

    x['PRIMARY'].header['BestRA']=racen
    x['PRIMARY'].header['BestDec']=deccen
    x['PRIMARY'].header['BestPA']=pa

    
    # print(xtab)

    ztab=Table(x['SLITMAP'].data)
    ztab['RA']=-99.
    ztab['Dec']=-99.

    for one in xtab:
        ztab['RA'][one['fiberid']-1]=one['RA']
        ztab['Dec'][one['fiberid']-1]=one['Dec']   

    # ztab.info()


    # Convert the modified table back to a binary table
    binary_table_hdu = fits.table_to_hdu(ztab)

    # Update the original FITS file with the modified binary table
    x['SLITMAP'].data = binary_table_hdu.data

    if outfile=='':
        outfile=filename.replace('lvmCFrame','XCFrame')
    
    x.writeto(outfile,overwrite=True)

    
def steer(argv):
    '''
    This is just a simple steering routine
    '''

    xfiles=[]


    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print (__doc__)
            return
        elif argv[i]=='-':
            print('Error: Unknown option:',argv)
            return
        else:
            xfiles.append(argv[i])
        i+=1

    for one in xfiles:

        try:
            xfoo=eval(one)
            one='lvmCFrame-%08d.fits' % xfoo
        except:
            pass

        print('\nBeginning %s' % one)
        doit(one)






# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
