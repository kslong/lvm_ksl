#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Run the ESO sky model remotely given alt-az and the postion of the moon


Command line usage (if any):

    usage: RunSkyCalcAltAz.py filename

Description:  

Primary routines:

    doit

Notes:

Generate sky spectra with Sky_Calc over a range of azimuths to see how much things vary
 
Skycalc can be found [here](https://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC).  There is a [command line interface](https://www.eso.org/observing/etc/doc/skycalc/helpskycalccli.html) for it as well, which should allow one to create a series of models for it.  There is a text and json file interface that one can use to work with the cli.  
 
As indicated in the documentation above, one must use pip to install the skycalc_cli.  If you are using astroconda, then the command is 

pip install skycalc_cli
 
If you use the command described in the documenation, then you must make sure that wherever it is installed is in ones path

History:

250421 ksl Coding begun

'''

import sys
from astropy.io import ascii,fits
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess


# The values below can be set to the defaults that one wishes to have


default='''
{
        "airmass": 1.0,
        "pwv_mode": "season",
        "season": 0,
        "time": 0,
        "pwv": 3.5,
        "msolflux": 130.0,
        "incl_moon": "Y",
        "moon_sun_sep": 90.0,
        "moon_target_sep": 45.0,
        "moon_alt": 45.0,
        "moon_earth_dist": 1.0,
        "incl_starlight": "Y",
        "incl_zodiacal": "Y",
        "ecl_lon": 135.0,
        "ecl_lat": 90.0,
        "incl_loweratm": "Y",
        "incl_upperatm": "Y",
        "incl_airglow": "Y",
        "incl_therm": "N",
        "therm_t1": 0.0,
        "therm_e1": 0.0,
        "therm_t2": 0.0,
        "therm_e2": 0.0,
        "therm_t3": 0.0,
        "therm_e3": 0.0,
        "vacair": "air",
        "wmin": 360.0,
        "wmax": 980.0,
        "wgrid_mode": "fixed_wavelength_step",
        "wdelta": 0.015,
        "wres": 20000,
        "lsf_type": "Gaussian",
        "lsf_gauss_fwhm": 5.0,
        "lsf_boxcar_fwhm": 5.0,
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


def write_inputs(xdefault=default,alt=90,moon_alt=-20,date=2021.5, time_of_night=2, solflux=200, outroot='test'):
    '''
    alt is the altitude above the horizon
    date is a year date in decimal degees it is used to determine the season, if date is 0, the sky glow for the
    average of the year will be produced
    
    time of night indicates whether to caculate the airflow for the average (0) or first, second, third (1,2,3)
    part of the night
    
    
    '''
    
    xdict=json.loads(xdefault)
    
    airmass=1./(np.sin(alt/57.29578))
    update(xdict,'airmass',airmass)
    
    update(xdict,'moon_alt',moon_alt)
    
    update(xdict,'time',time_of_night)
    

    if date!=0:
        date-=int(date)
        season=6*date+1
        season=round(season)
        update(xdict,'season',season)
    else:
        update(xdict,'pwv_mode','pwv')
        update(xdict,'season',0)
        
    
        
    
    moon_targ=xdict.get('moon_target_sep')
    print(np.abs(moon_alt-alt),moon_targ)
    if np.abs(moon_alt-alt)> moon_targ:
        update(xdict,'moon_target_sep',float(np.abs(moon_alt-alt)))
        
    update(xdict,'msolflux',solflux)
    
    # print(xdict)
    
    jsonString=json.dumps(xdict, indent = 4, sort_keys=False)
    
    print('Writing sky calc inputs to %s.json'  % outroot)
    jsonFile = open("%s.json" % outroot, "w")
    jsonFile.write(jsonString)
    jsonFile.close()
    


def just_run_SkyCalc(xinput='test.json',almanac='',outroot='',print_output=True):
    
    '''
    Run skycalc on a valid set of inputs
    '''
    
    if xinput.count('.json')==0:
        xroot=xinput
        xinput='%s.json' % xinput
    else:
        xroot=xinput.replace('.json','')
    
    if outroot=='':
        outroot=xroot
    
    xstring= 'skycalc_cli -i %s.json -o %s.fits' % (xroot,outroot)
    command_line=xstring.split()
    print(command_line)
        
    result=subprocess.run(command_line,capture_output=True,text=True)
    if print_output:
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)



def run_SkyCalc(xdefault=default,alt=90,moon_alt=-20,date=2021.5, time_of_night=2, solflux=200, outroot='zenith'):
    write_inputs(xdefault,alt,moon_alt,date,time_of_night,solflux, outroot)
    just_run_SkyCalc(outroot,almanac='',outroot='',print_output=True)
    return





def steer(argv):
    '''
    '''

    alt=33
    moon_alt=30
    solflus=200
    outroot=''

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-alt':
            i+=1
            alt=eval(argv[i])
        elif argv[i]=='-moon_alt':
            i+=1
            moon_alt=eval(argv[i])
        elif argv[i][0]=='-':
            print('Error: Unknown option: ',argv)
            return
        i+=1

    outroot='Mod_%03d_%03d' % (alt,moon_alt)
    run_SkyCalc(alt=alt,moon_alt=moon_alt,outroot=outroot)



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)

