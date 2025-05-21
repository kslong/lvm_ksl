#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Use SkyCalc to find the sky for a given RA and DEC and time


Command line usage (if any):

    Usage:  SkyCalObs.py  RA  DEC time

Description:  

Primary routines:

    just_run_SkyCalc_from_observation(xinput='test.json',almanac='',outroot='',msol=0,print_output=False):
    run_SkyCalc_from_observation(xdefault=default,ra=121.75,dec=-29.7,xtime="2012-07-17T21:12:14", outroot='test',msol=137,print_output=False):

Notes:
Run SkyCalc based on RA, Dec, and time
# 
Skycalc can be found [here](https://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC).  There is a [command line interface](https://www.eso.org/observing/etc/doc/skycalc/helpskycalccli.html) for it as well, which should allow one to create a series of models for it.  There is a text and json file interface that one can use to work with the cli.  
 
As indicated in the documentation above, one must use pip to install the skycalc_cli.  If you are using astroconda, then the command is 
 
pip install skycalc_cli
If you use the command described in the documenation, then you must make sure that wherever it is installed is in one's path.
In this script SkyCalc is run based on RA, a Dec, and an observation time.  (The observation time can be given in terms of MJD as a float, or as a time string).  If one one provides a time which is non-physical, then the final routine will produce and error, but will not really tell you what is wrong.
# 

                                       
History:

250421 ksl Coding begun

'''

import sys
import os
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess
from astropy.time import Time
from datetime import datetime


default='''
{
    "ra": 121.75,
    "dec": -29.7,
    "date": "2012-07-17T21:12:14",
    "observatory": "lasilla"
}
'''


def update(xdict,param='wmax',value=1000.):
    '''
    Update a paremeter in a dictionary.  Note that this
    occurs in place so that one does not have to return anything
    '''
    xdict.update({param : value} )
    return



def write_obs_inputs(xdefault=default,ra=121.75,dec=-29.7,xtime="2012-07-17T21:12:14", msol=-1,outroot='test'):
    '''
    ra  - in decimal degress
    dec - in decimal degrees
    time - eitherr string (timestamp, YYYY-MM-DDThh:mm:ss) | float (Modified Julian Day)
    
    '''
    
    xdict=json.loads(xdefault)
    update(xdict,'ra',ra)
    
    update(xdict,'dec',dec)
    
    print('A ',msol)

    xtime=Time(parse_to_datetime(xtime))
    isot_whole_seconds = xtime.datetime.strftime('%Y-%m-%dT%H:%M:%S')
    update(xdict,'date',isot_whole_seconds)
    # print('XXX',xtime)
    year=eval(xtime.datetime.strftime('%Y'))
    fudge_sol=False
    if year>2024:
        fudge_sol=True
        msol=0

    print('B ',msol,year)

    if msol>0:
        xdict['msolflux']=msol
    
    
    # print(xdict)
    
    jsonString=json.dumps(xdict, indent = 4, sort_keys=False)
    
    # print('Writing sky calc inputs to %s.json'  % outroot)
    jsonFile = open("%s.json" % outroot, "w")
    jsonFile.write(jsonString)
    jsonFile.close()
    

import time

xdefaults='''
{
    "vacair": "air",
    "wmin": 360.0,
    "wmax": 980.0,
    "wgrid_mode": "fixed_wavelength_step",
    "wdelta": 0.05,
    "lsf_type": "Gaussian",
    "lsf_gauss_fwhm": 2.0,
    "lsf_boxcar_fwhm": 2.0,
    "observatory": "lasilla"
}
'''


ydefaults='''
{
    "vacair": "air",
    "wmin": 360.0,
    "wmax": 980.0,
    "wgrid_mode": "fixed_wavelength_step",
    "wdelta": 0.05,
    "lsf_type": "Gaussian",
    "lsf_gauss_fwhm": 2.0,
    "lsf_boxcar_fwhm": 2.0,
    "observatory": "lasilla"
}
'''



def just_run_SkyCalc_from_observation(xinput='test.json',almanac='',outroot='',msol=0,print_output=False):
    
    '''
    Run skycalc on a valid set of inputs
    '''
    
    print('xsol', msol)
    if xinput.count('.json')==0:
        xroot=xinput
        xinput='%s.json' % xinput
    else:
        xroot=xinput.replace('.json','')
    
    if outroot=='':
        outroot=xroot

    xstring= 'skycalc_cli -i xsky.json -a %s.json -o %s.fits' % (xroot,outroot)
    g=open('xsky.json','w')
    g.write(xdefaults)
    g.flush()
    g.close()
    time.sleep(1)
    
    # print(xstring)
    command_line=xstring.split()
    print(command_line)
    Xerror=False
        
    try:
        result=subprocess.run(command_line,capture_output=True,text=True)
    except subprocess.CalledProcessError as e:
        print('Error: subprocess error')
        Xerror==True
    except FileNotFoundError:
        print("Error: 'skycalc_cli' not found in your PATH.")
        return ''
    except Exception:
        print("Unexpected Error")
        Xerror=True

    # Detect hidden errors
    if result.returncode != 0 or "Traceback" in result.stderr:
        print("ERROR: skycalc_cli failed to execute successfully.")
        Xerror=True

    if print_output or Xerror:
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)


    if Xerror:
        return ''


    if is_recent('%s.fits' % xroot)==False:
        print('Cannot verify that file %s.fits was modified or possibly even created' % xroot)
        return ''

 
    return outroot

def is_recent(filepath, minutes=10):
    """Return True if file was modified in the last `minutes`."""
    if not os.path.isfile(filepath):
        return False
    mtime = os.path.getmtime(filepath)  # modification time in seconds
    age_seconds = time.time() - mtime
    return age_seconds < (minutes * 60)

def parse_to_datetime(time_input):
    """
    Convert a time input (ISO string or MJD number/string) to a Python datetime object.

    Parameters
    ----------
    time_input : str or float
        ISO time string (e.g. '2025-04-21T03:00:00') or MJD (e.g. 60394.125 or '60394.125')

    Returns
    -------
    dt : datetime.datetime
        Corresponding UTC datetime object
    """
    try:
        # If input looks like a number or numeric string, treat as MJD
        if isinstance(time_input, (int, float)) or str(time_input).replace('.', '', 1).isdigit():
            t = Time(float(time_input), format='mjd', scale='utc')
        else:
            t = Time(time_input, format='isot', scale='utc')
        return t.to_datetime()
    except Exception as e:
        raise ValueError(f"Unrecognized time format: {time_input!r}") from e

def run_SkyCalc_from_observation(xdefault=default,ra=121.75,dec=-29.7,xtime="2012-07-17T21:12:14", outroot='',msol=137,print_output=False):
    if outroot=='':
        outroot='SkyC_%5.1f_%5.1f' % (ra,dec)
    # print('XXX',outroot,xtime)
    # isot_whole_seconds = xtime.datetime.strftime('%Y-%m-%dT%H:%M:%S')
    # print('XXX',outroot,xtime.isot,isot_whole_seconds)

    # write_obs_inputs(xdefault,ra,dec,isot_whole_seconds, msol, outroot)
    write_obs_inputs(xdefault,ra,dec,xtime, msol, outroot)
    xroot=just_run_SkyCalc_from_observation(outroot,almanac='',outroot=outroot,msol=msol,print_output=print_output)
    return xroot


def steer(argv):
    '''
    Usage:  SkyCalObs RA  DEC time
    '''

    ra=-1.0
    dec=-100
    xtime=-1
    i=1
    outroot=''
    msol=137

    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][:4]=='-out':
            i+=1
            outroot=argv[i]
        # elif argv[i][0]=='-':
        #     print('Error: Unknown option: ',argv)
        #     return
        elif ra<0.0:
            ra=eval(argv[i])
        elif dec<-90:
            dec=eval(argv[i])
        elif xtime<0:
            xtime=argv[i]
        i+=1


    
    run_SkyCalc_from_observation(xdefault=default,ra=ra,dec=dec,xtime=xtime, outroot=outroot)





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
