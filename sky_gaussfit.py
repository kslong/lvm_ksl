#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Extract fluxes etc from a standard set of emission lines from
the scince fibers of a sky-subtracted LVM exposure


Command line usage (if any):

    usage: lvm_bootstrap.py [-lmc] [-smc] [-out root] filename

    where 

    -lmc or -smc applies a velocity offset for fitting
    -out root sets the rootname for the output file
    filename is the aname of an SFrame compatiable file


Description:  

Primary routines:

    doit

Notes:
                                       
History:

240604 ksl Coding begun

'''

import sys


import os 
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.table import Table,vstack, hstack
from astropy.io import fits
import numpy as np
from astropy.modeling import models, fitting
from glob import glob


from scipy.optimize import curve_fit
import astropy.units as u
from datetime import datetime
from lvm_ksl.lvm_gaussfit import fit_gaussian_to_spectrum,fit_double_gaussian_to_spectrum,save_fit,scifib,check_for_nan,plot_one,plot_all,clean


def do_one(spectrum_table,vel=0.,xplot=False):
    '''
    Completely process a single spectrum
    '''

    clean()
    
    zz=1.+ (vel/3e5)

    records=[]

    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='oi_a',init_wavelength=6300, init_fwhm=1., wavelength_min=6290, wavelength_max=6305)
        records.append(results)
        if xplot:
            save_fit('oi_a',xspec)
    except Exception as e:
        print(f"Fitting oi_a; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    

    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='oi_b',init_wavelength=6300, init_fwhm=1., wavelength_min=6290, wavelength_max=6305)
        records.append(results)
        if xplot:
            save_fit('oi_b',xspec)
    except Exception as e:
        print(f"Fitting oi_b; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    

    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='ha', init_wavelength=zz*6563, init_fwhm=1., wavelength_min=zz*6555, wavelength_max=zz*6570)
        records.append(results)
        if xplot:
            save_fit('ha',xspec)
    except Exception as e:
        print(f"Fitting Ha An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")

    
            
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='sii_a',init_wavelength=zz*6716, init_fwhm=1., wavelength_min=zz*6706, wavelength_max=zz*6726)
        records.append(results)
        if xplot:
            save_fit('sii_a',xspec)
    except Exception as e:
        print(f"Fitting [SII]6716:  An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='sii_b',init_wavelength=zz*6731, init_fwhm=1., wavelength_min=zz*6721, wavelength_max=zz*6741)
        records.append(results)
        if xplot:
            save_fit('sii_b',xspec)
    except Exception as e:
        print(f"Fitting [SII]6731:  An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")


    # These are airglow lines and so should not be redshifted
    
    try:
        xname='sky6553'
        xwcen=6553.
        xwmin=6549.
        xwmax=6556.
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line=xname,init_wavelength=xwcen, init_fwhm=1., wavelength_min=xwmin, wavelength_max=xwmax)
        records.append(results)
        if xplot:
            save_fit(xname,xspec)
    except Exception as e:
        print(f"Fitting %s:  An exception occurred: {e}" % xname)
        print(f"Exception type: {type(e).__name__}")

    
    try:
        xname='sky6533'
        xwcen=6533.04
        xwmin=6528.
        xwmax=6538.
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line=xname,init_wavelength=xwcen, init_fwhm=1., wavelength_min=xwmin, wavelength_max=xwmax)
        records.append(results)
        if xplot:
            save_fit(xname,xspec)
    except Exception as e:
        print(f"Fitting %s:  An exception occurred: {e}" % xname)
        print(f"Exception type: {type(e).__name__}")

    
    try:
        xname='sky6577'
        xwcen=6577.2
        xwmin=6572.
        xwmax=6582.
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line=xname,init_wavelength=xwcen, init_fwhm=1., wavelength_min=xwmin, wavelength_max=xwmax)
        records.append(results)
        if xplot:
            save_fit(xname,xspec)
    except Exception as e:
        print(f"Fitting %s:  An exception occurred: {e}" % xname)
        print(f"Exception type: {type(e).__name__}")


    
    xname='sky6912'
    xwcen=6912.623
    xwmin=6907.
    xwmax=6917.
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line=xname,init_wavelength=xwcen, init_fwhm=1., wavelength_min=xwmin, wavelength_max=xwmax)
        records.append(results)
        if xplot:
            save_fit(xname,xspec)
    except Exception as e:
        print(f"Fitting %s:  An exception occurred: {e}" % xname)
        print(f"Exception type: {type(e).__name__}")

    
    xname='sky6923'
    xwcen=6923.220
    xwmin=6918.
    xwmax=6928.
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line=xname,init_wavelength=xwcen, init_fwhm=1., wavelength_min=xwmin, wavelength_max=xwmax)
        records.append(results)
        if xplot:
            save_fit(xname,xspec)
    except Exception as e:
        print(f"Fitting %s:  An exception occurred: {e}" % xname)
        print(f"Exception type: {type(e).__name__}")

    
    xname='sky6939'
    xwcen=6939.521
    xwmin=6934.
    xwmax=6944.
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line=xname,init_wavelength=xwcen, init_fwhm=1., wavelength_min=xwmin, wavelength_max=xwmax)
        records.append(results)
        if xplot:
            save_fit(xname,xspec)
    except Exception as e:
        print(f"Fitting %s:  An exception occurred: {e}" % xname)
        print(f"Exception type: {type(e).__name__}")


    try:
        ztab=hstack(records)
        return ztab
    except:
        print('Nothing fit for this spectrum')
        return []

    return ztab


def do_individual(filenames,vel,outname,xplot=False):
    '''
    This is to process individual spectra from an astropy table
    containing a WAVE and FLUX column
    '''
    xresults=[]
    xbad=[]
    xgood=[]

    plot_dir='Gauss_plot'

    os.makedirs(plot_dir,exist_ok=True)

    for one_file in filenames:
        try:
            xtab=ascii.read(one_file)
            foo=xtab['WAVE']
            foo=xtab['FLUX']
        except:
            xbad.append(one_file)
            continue
        try:
            xtab=ascii.read(one_file)
            results=do_one(xtab,vel,xplot)
            word=one_file.split('/')
            root=word[-1]
            root=root.replace('.txt','')
            print(results)
            results['Object']=root
            xresults.append(results)
            xgood.append(one_file)
            plot_all()
            word=one_file.split('/')
            proot=word[-1].replace('.txt','')
            plot_name='%s/%s.png' % (plot_dir,proot)
            print(plot_name)
            plt.savefig(plot_name)
        except:
            print('Could not analyze %s' % one_file)
            xbad.append(one_file)

    if len(xgood)>0:
        for one_file in xgood:
            print('Successfully fit %s' % one_file)
    if len(xbad)>0:
        print('!! Failed to fit %s' % one_file)
    if len(xresults)==0:
        print('Duh')
        return
    ftab=vstack(xresults)

    current_date = datetime.now()
    formatted_date = current_date.strftime("%d%m%y")
    if outname=='':
        outname=formatted_date

    if len(filenames)==1:
        outname='Gauss_%s.txt' % root
        ftab.write(outname,format='ascii.fixed_width_two_line',overwrite=True)
    else:
        outname='Gauss_%s.txt' % outname
        ftab.write(outname,format='ascii.fixed_width_two_line',overwrite=True)

def do_all(filename='data/lvmSFrame-00009088.fits',vel=0.0,outname=''):
    '''
    Do all of the spectra in a rss fits file
    '''
    try:
        x=fits.open(filename)
    except:
        print('Error: Could not open ',filename)
        return

    wave=x['WAVE'].data
    flux=x['FLUX'].data

        
    # Now get the good fibers from the science telecsope
    slittab=Table(x['SLITMAP'].data)
    good=scifib(slittab,'science')
    # print(len(good))
    # print(good)
    results=[]
    records=[]
    for i in range(len(good)):
        j=good['fiberid'][i]-1
        one_spec=Table([wave,flux[j]],names=['WAVE','FLUX'])
        if check_for_nan(flux[j])==False:
            rtab=do_one(spectrum_table=one_spec,vel=vel,xplot=False)
            if len(rtab)>0:
                rtab['fiberid']=good['fiberid'][i]
                rtab['ra']=good['ra'][i]
                rtab['dec']=good['dec'][i]
                records.append(rtab)
            else:
                print('Nothing fit for fiber %d at %.2f %.2f'  % (good['fiberid'][i],good['ra'][i],good['dec'][i]))
        else:
            print('Too many nans for  fiber %d at %.2f %.2f'  % (good['fiberid'][i],good['ra'][i],good['dec'][i]))

    results=vstack(records)

    if outname=='':
        outname=filename.split('/')[-1]
        outname=outname.replace('.fits','.txt')
    else:
        outname=outname+'.txt'

    columns=results.colnames
    for one in columns:
        if one.count('flux'):
            results[one].format='.2e'
        elif one.count('wave'):
            results[one].format='.2f'
        elif one.count('fwhm'):
            results[one].format='.2f'
        elif one.count('rmse'):
            results[one].format='.2e'
        elif one.count('back'):
            results[one].format='.2e'
        elif one.count('ra'):
            results[one].format='.5f'
        elif one.count('dec'):
            results[one].format='.5f'
        elif one.count('fiberid'):
            results[one].format='d'
        else:
            print('Did not reformat ',one)



    results.write(outname,format='ascii.fixed_width_two_line',overwrite=True)
    return results
        
def analyze(xtab):
    print('There are %d spectra' % len(xtab))
    xx=xtab[xtab['eflux_ha']!=999.]
    print('There are %d spectra with ha measured' % len(xx))
    xx=xx[xx['eflux_sii_a']!=999.]
    print('Of these, %d spectra also have s2 measured' % len(xx))
    xx=xx[xx['eflux_oi_a']!=999.]
    print('Of these, %d spectra also have o1 measured' % len(xx))
    print('Measured means having an error calculated')

def steer(argv):
    filename=''
    lmc=262.
    smc=146.
    outname=''
    fitsfiles=[]
    specfiles=[]

    vel=0
    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-lmc':
            vel=lmc
        elif argv[i]=='-smc':
            vel=smc
        elif argv[i][0:4]=='-out':
            i+=1
            outname=argv[i]
        elif argv[i]=='-v':
            i+=1
            vel=eval(argv[i])
        elif argv[i][0]=='-':
            print('Unknown options :',argv)
            return
        elif argv[i].count('.fits'):
            fitsfiles.append(argv[i])
        elif argv[i].count('.txt'):
            specfiles.append(argv[i])
        else:
            print('Unknown options :',argv)
            return
        i+=1

    for one_file in fitsfiles:
        results=do_all(one_file,vel,outname)
        analyze(results)

    if len(specfiles)>0:
        do_individual(specfiles,vel,outname)

    return




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
