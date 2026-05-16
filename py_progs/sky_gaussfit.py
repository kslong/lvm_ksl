#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Fit Gaussians to a fixed set of nebular emission lines and airglow lines
in the science fibers of a sky-subtracted LVM SFrame exposure, producing
one output row per fiber.


Command line usage::

    sky_gaussfit.py [-lmc] [-smc] [-v vel] [-out root] filename [filename ...]

    where

    -lmc         Apply the LMC radial velocity (~262 km/s) when fitting nebular lines
    -smc         Apply the SMC radial velocity (~146 km/s) when fitting nebular lines
    -v vel       Apply an arbitrary radial velocity (km/s) when fitting nebular lines
    -out root    Set the root name for the output file (default: derived from input filename)
    filename     One or more SFrame FITS files or ASCII spectrum files (WAVE, FLUX columns)


Description:

For each good science fiber (telescope == 'Sci', fibstatus == 0) the script fits
single Gaussians to the following lines:

    Nebular lines (wavelengths shifted by the supplied velocity):
        [OI]  6300.304  (oi_a, oi_b)
        [NII] 6548.04   (nii_a)
        Ha    6562.80   (ha)
        [NII] 6583.46   (nii_b)
        [SII] 6716.440  (sii_a)
        [SII] 6730.815  (sii_b)

    Airglow lines (fitted at fixed, unshifted wavelengths):
        sky5577   5577.34 A
        sky6300   6300.31 A  (sky OI; nebular OI also fit separately above)
        sky6363   6363.78 A
        sky6533   6533.04 A
        sky6553   6553.0  A
        sky6577   6577.2  A
        sky6912   6912.62 A
        sky6923   6923.22 A
        sky6939   6939.52 A
        sky7358   7358.68 A
        sky7392   7392.21 A
        sky7914   7913.72 A
        sky8344   8344.61 A
        sky8399   8399.18 A
        sky8827   8827.11 A
        sky8988   8988.38 A
        sky9552   9552.55 A
        sky9719   9719.84 A

For each line the fit returns the flux, centroid wavelength, FWHM, background
level, and RMSE.  Fibers with too many NaN pixels are skipped.

Primary routines:

    do_one    - fit all lines for a single spectrum table
    do_all    - iterate over all science fibers in an SFrame FITS file
    do_individual - iterate over a list of ASCII spectrum files


Output:

When the input is an SFrame FITS file, the output is an ASCII fixed-width table
with one row per successfully fit fiber.  Columns cover the fit parameters
(flux, wave, fwhm, back, rmse) for each line together with fiberid, ra, and dec.
The output filename defaults to the input filename with .fits replaced by .txt,
or <root>.txt if -out is supplied.

When the input is one or more ASCII spectrum files, the output is written to
Gauss_<root>.txt (multiple files) or Gauss_<stem>.txt (single file), where
<stem> is the input filename without path or extension.


Notes:

The airglow lines are present in sky-subtracted SFrame data regardless of the
target velocity.  Their measured fluxes reflect whatever signal remains at those
wavelengths after sky subtraction.

History:

240604 ksl Coding begun
260516 ksl Add 12 airglow lines from gauss_offset (ESO UVES wavelengths);
           sky6300 fit at fixed wavelength alongside existing nebular OI
260516 ksl Add nii_a (6548.04) and nii_b (6583.46); correct SII wavelengths
           to 6716.440 and 6730.815

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
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='nii_a', init_wavelength=zz*6548.04, init_fwhm=1., wavelength_min=zz*6543, wavelength_max=zz*6553)
        records.append(results)
        if xplot:
            save_fit('nii_a',xspec)
    except Exception as e:
        print(f"Fitting [NII]6548: An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")

    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='ha', init_wavelength=zz*6562.80, init_fwhm=1., wavelength_min=zz*6555, wavelength_max=zz*6570)
        records.append(results)
        if xplot:
            save_fit('ha',xspec)
    except Exception as e:
        print(f"Fitting Ha An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")

    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='nii_b', init_wavelength=zz*6583.46, init_fwhm=1., wavelength_min=zz*6578, wavelength_max=zz*6593)
        records.append(results)
        if xplot:
            save_fit('nii_b',xspec)
    except Exception as e:
        print(f"Fitting [NII]6583: An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")

    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='sii_a',init_wavelength=zz*6716.440, init_fwhm=1., wavelength_min=zz*6706, wavelength_max=zz*6726)
        records.append(results)
        if xplot:
            save_fit('sii_a',xspec)
    except Exception as e:
        print(f"Fitting [SII]6716:  An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")

    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='sii_b',init_wavelength=zz*6730.815, init_fwhm=1., wavelength_min=zz*6721, wavelength_max=zz*6741)
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

    # Additional airglow lines (ESO UVES wavelengths).
    # sky6300 is the sky OI line; the nebular OI doublet is fit separately above.
    for xname, xwcen, xwmin, xwmax in [
        ('sky5577', 5577.34668,  5572., 5582.),
        ('sky6300', 6300.308594, 6295., 6305.),
        ('sky6363', 6363.782715, 6358., 6368.),
        ('sky7358', 7358.680176, 7353., 7363.),
        ('sky7392', 7392.209961, 7387., 7397.),
        ('sky7914', 7913.717773, 7908., 7918.),
        ('sky8344', 8344.613281, 8339., 8349.),
        ('sky8399', 8399.175781, 8394., 8404.),
        ('sky8827', 8827.112305, 8822., 8832.),
        ('sky8988', 8988.383789, 8983., 8993.),
        ('sky9552', 9552.546875, 9547., 9557.),
        ('sky9719', 9719.838867, 9714., 9724.),
    ]:
        try:
            results, xspec = fit_gaussian_to_spectrum(
                spectrum_table, line=xname,
                init_wavelength=xwcen, init_fwhm=1.,
                wavelength_min=xwmin, wavelength_max=xwmax)
            records.append(results)
            if xplot:
                save_fit(xname, xspec)
        except Exception as e:
            print('Fitting %s: An exception occurred: %s' % (xname, e))

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
