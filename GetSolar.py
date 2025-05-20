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

    usage: GetSolar.py 

Description:  

    There are no switches for this

Primary routines:

    doit

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
    location=location.replace('./getSolar.py','data')
    get_text(xurl,location)

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    doit()
