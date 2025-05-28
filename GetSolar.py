#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Get the latest version of the Solar flux file and store them 
in the data directory for lvm_ksl

See

https://www.spaceweather.gc.ca/solar_flux_data/daily_flux_values/fluxtable.txt

as descibed by

https://lasp.colorado.edu/lisird/data/penticton_radio_flux


And store this in the data directory for under lvm_ksl

Command line usage (if any):

    usage:  GetSolar.py [-h] [-retrieve]  time1 ...


    where 
        -h displays the __doc__ file and quites
        -retrieve causes a new fluxtable.txt to be retrieved

    and then the solar fluxes are retrived for one or more times, given various formats

Description:  




Primary routines:

    doit  gets a new solar data file and stores it in the lvm_ksl/data directory
    get_flux obtains the solar flux for a specific time

Notes:
                                       
History:

250520 ksl Coding begun

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt


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
    if re.search(r'\d{4}-\d{2}-\d{2}', time_str) or re.search(r'\d{4}/\d{2}/\d{2}', time_str):
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
        Desired output format ('datetime', 'iso', 'mjd', or 'jd')

    Returns
    -------
    The time in the requested format
    """
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
            raise ValueError(f"Unsupported output format: {output_format}")
    except Exception as e:
        raise ValueError(f"Failed to convert to {output_format} format: {str(e)}")

# End of section that translates times.

# Replace with the actual URL of the text file


def get_text(url,location='.'):

    

    # Send a GET request to the URL
    response = requests.get(url)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
      # Open a file in write binary mode ('wb')
      outname="%s/solar.txt" % location
      with open(outname, "wb") as file:
        # Write the content of the response to the file
        file.write(response.content)
      print("File downloaded successfully! to %s" % outname)
    else:
      print(f"Error downloading file. Status code: {response.status_code}")


def get_source_location():
    file_path = inspect.getfile(get_text)
    return file_path


def doit(xurl='https://www.spaceweather.gc.ca/solar_flux_data/daily_flux_values/fluxtable.txt'):

    print('Getting solar data from:\n %s' % xurl)
    location=get_source_location()
    location=location.replace('GetSolar.py','data')
    get_text(xurl,location)


def get_flux(xtime='2025-04-21T03:00:00',filename='solar.txt'):

    location=get_source_location()
    print(location)
    location=location.replace('GetSolar.py','data')
    print(location)
    
    xtab=ascii.read('%s/%s' % (location,filename))
    jd=convert_time(xtime,'jd')
    # print(jd)
    ztab=xtab[xtab['fluxjulian']>jd-5]
    if len(ztab)==0:
        print('Error: Time jd %f string %s later than data available, using last set of data' % (jd,convert_time(xtime,'iso')))
        ztab=xtab[-10:]
    else:
        ztab=ztab[ztab['fluxjulian']<jd+5]

    solar_flux=np.median(ztab['fluxadjflux'])
    return solar_flux


def convert_to_number_if_possible(one_number):
    try:
        x=eval(one_number)
    except:
        x=on_number
    return x


def steer(argv):
    '''
    usage:  GetSolar.py -h -retrieve xxx
    '''

    times=[]
    xret=False
    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-retrieve':
            xret=True
        elif argv[i][0]=='-':
            print(__doc__)
            print ('Error: Unkown switch: ',argv)
            return 
        else:
            one_time=argv[i]
            one_time=convert_to_number_if_possible(one_time)
            times.append(one_time)

        i+=1

        if xret:
            doit()

        for one_time in times:
            flux=get_flux(one_time)
            print(one_time,' :  ',  flux)

        
# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
