#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Extract fluxes etc from a standard set of emission lines from
the scince fibers of a sky-subtracted LVM exposure


Command line usage (if any):

    usage: lvm_gaussfit.py --h] [-lmc] [-smc] [-v vel] [-stype SOURCE][-out root] filename ...

    where 

    -h prient this tdocumeantiaon and fits
    -lmc or -smc applies a velocity offset for fitting
    -vel whatever applies a velocity offset that the user specifies  
    -stype SOURCE or BACK - used only for spectra that has been created as text files
    -out root sets the rootname for the output file
    filename is the a name of an SFrame compatiable file or one or more txt files


Description:  

    The code can be used to perform gaussian to prominente lines from either a 
    rss spectrum file (including but not limited to a standard SFrame file) or
    alternative to an spectrum that has been extracted in the a table containing
    at least a column for WAVE, FLUX, ERROR.

    For an extracted spectrum of this type, which may have additional columns, at
    present one containing the spectrum before local subtraction and another containing
    the back ground spectrum which was subtracted, one can redirect the fitting with
    the stype option.  

Primary routines:

    doit

Notes:
                                       
History:

240604 ksl Coding begun

'''

import sys


import os 
from astropy.io import ascii, fits
import matplotlib.pyplot as plt
from astropy.table import Table,vstack, hstack
import numpy as np
from glob import glob


from scipy.optimize import curve_fit
import astropy.units as u
from datetime import datetime
from astropy.table import Table


import numpy as np
from astropy.table import Table
from lmfit import Model, Parameters

import warnings
import traceback

def custom_warning_handler(message, category, filename, lineno, file=None, line=None):
    if "std_dev==0" in str(message):
        print(f"\n{category.__name__}: {message}")
        traceback.print_stack()

warnings.showwarning = custom_warning_handler

def patch_stderr_fraction(result, min_frac=0.01, verbose=False):
    for name, par in result.params.items():
        # Calculate fallback stderr
        abs_val = abs(par.value)
        min_stderr = min_frac * abs_val if abs_val != 0 else min_frac

        if par.stderr is None:
            if verbose:
                print(f"⚠️  Parameter '{name}' has stderr=None → setting to {min_stderr:.3e}")
            par.stderr = min_stderr
        elif par.stderr == 0.0:
            if verbose:
                print(f"⚠️  Parameter '{name}' has stderr=0 → setting to {min_stderr:.3e}")
            par.stderr = min_stderr
        # elif par.stderr < min_stderr:
        #     if verbose
        #       print(f"⚠️  Parameter '{name}' has very small stderr={par.stderr:.3e} "
        #           f"(< {min_frac:.0%} of value) → setting to {min_stderr:.3e}")
        #     par.stderr = min_stderr

    return result


def fit_gaussian_to_spectrum(spectrum_table, line, init_wavelength, init_fwhm, wavelength_min, wavelength_max, verbose=False):
    """
    Fit a Gaussian with a constant background to a spectral line using lmfit.
    
    Parameters:
    spectrum_table (astropy.table.Table): Table with 'WAVE', 'FLUX', and 'ERROR' columns
    line (str): Name of the line being fit (for output column names)
    init_wavelength (float): Initial guess for the center wavelength
    init_fwhm (float): Initial guess for the FWHM
    wavelength_min, wavelength_max (float): Wavelength range for fitting
    
    Returns:
    qtab (astropy.table.Table): Table with fit results and uncertainties
    fit_table (astropy.table.Table): Table with original and fitted values
    """
    # Define model functions
    def gaussian(x, flux, center, fwhm, background):
        """Gaussian function parameterized by flux, center, FWHM, and background"""
        sigma = fwhm / 2.355
        amp = flux / (sigma * np.sqrt(2 * np.pi))
        return amp * np.exp(-0.5 * ((x - center) / sigma)**2) + background
    
    # Extract and clean data
    mask = (spectrum_table['WAVE'] >= wavelength_min) & (spectrum_table['WAVE'] <= wavelength_max)
    x = spectrum_table['WAVE'][mask]
    y = spectrum_table['FLUX'][mask]
    
    # Get errors if available
    if 'ERROR' in spectrum_table.colnames:
        errors = spectrum_table['ERROR'][mask]
    else:
        errors = None
    
    # Filter for finite values
    finite_mask = np.isfinite(y)
    if errors is not None:
        finite_mask = finite_mask & np.isfinite(errors)
    
    x = x[finite_mask]
    y = y[finite_mask]
    if errors is not None:
        errors = errors[finite_mask]
        # Ensure errors are not too small to prevent numerical issues
        min_error_value = np.max(np.abs(y)) * 1e-6
        errors = np.maximum(errors, min_error_value)
    
    # Initial parameter estimates
    dx = x[1] - x[0]
    init_flux = np.sum(y * dx)
    init_background = np.median(y)
    
    # Set up lmfit model
    gmodel = Model(gaussian)
    params = gmodel.make_params(
        flux=init_flux,
        center=init_wavelength,
        fwhm=init_fwhm,
        background=init_background
    )
    
    # Set parameter bounds
    params['fwhm'].min = 0.1 * init_fwhm  # Prevent unrealistically narrow lines
    params['fwhm'].max = 5.0 * init_fwhm   # Prevent unrealistically broad lines
    params['center'].min = wavelength_min
    params['center'].max = wavelength_max
    
    # Perform the fit 
    if verbose:
        print('Starting ...')
    if errors is not None:
        if verbose:
            print('With Errors')
        if np.any(errors <= 0):
            print("Some Errors are zero or negative:", errors[errors <= 0])
        with warnings.catch_warnings():
            warnings.filterwarnings( "ignore",
            message="Using UFloat objects with std_dev==0 may give unexpected results.",
            category=UserWarning,
            module="uncertainties.core"
            )
            result = gmodel.fit(y, params, x=x, weights=1.0 / errors)
            patch_stderr_fraction(result, min_frac=0.01)
    else:
        if verbose:
            print('No Errors!!!')
        result = gmodel.fit(y, params, x=x)
    if verbose:
        print('Finishing ...')
    
    # Extract fitted parameters and uncertainties
    flux = result.params['flux'].value
    flux_err = result.params['flux'].stderr if result.params['flux'].stderr is not None else np.nan
    
    wavelength = result.params['center'].value
    wave_err = result.params['center'].stderr if result.params['center'].stderr is not None else np.nan
    
    fwhm = result.params['fwhm'].value
    fwhm_err = result.params['fwhm'].stderr if result.params['fwhm'].stderr is not None else np.nan
    
    background = result.params['background'].value
    back_err = result.params['background'].stderr if result.params['background'].stderr is not None else np.nan
    
    # Calculate statistics
    y_fit = result.best_fit
    rmse = np.sqrt(np.mean((y - y_fit) ** 2))
    
    # Use lmfit's built-in reduced chi-squared calculation
    reduced_chi2 = result.redchi
    
    # Create output tables
    fit_table = Table([x, y, y_fit], names=('WAVE', 'FLUX', 'Fit'))
    if errors is not None:
        fit_table['ERROR'] = errors
    
    col_names = [f'{param}_{line}' for param in 
                 ['flux', 'eflux', 'wave', 'ewave', 'fwhm', 'efwhm', 'back', 'eback', 'rmse', 'chi2']]
    
    qtab = Table(names=col_names)
    qtab.add_row([flux, flux_err, wavelength, wave_err, fwhm, fwhm_err, background, back_err, rmse, reduced_chi2])
    
    return qtab, fit_table

def fit_double_gaussian_to_spectrum(spectrum_table, line, init_wavelength1, init_wavelength2, init_fwhm, wavelength_min, wavelength_max):
    """
    Fit two Gaussians with a shared FWHM and a constant background using lmfit.
    
    Parameters:
    spectrum_table (astropy.table.Table): Table with 'WAVE', 'FLUX', and 'ERROR' columns
    line (str): Name of the line being fit (for output column names)
    init_wavelength1, init_wavelength2 (float): Initial guesses for the center wavelengths
    init_fwhm (float): Initial guess for the shared FWHM
    wavelength_min, wavelength_max (float): Wavelength range for fitting
    
    Returns:
    qtab (astropy.table.Table): Table with fit results and uncertainties
    fit_table (astropy.table.Table): Table with original and fitted values
    """
    # Define model function
    def double_gaussian(x, flux1, flux2, center1, center2, fwhm, background):
        """Two Gaussians with shared FWHM and a background"""
        sigma = fwhm / 2.355
        amp1 = flux1 / (sigma * np.sqrt(2 * np.pi))
        amp2 = flux2 / (sigma * np.sqrt(2 * np.pi))
        g1 = amp1 * np.exp(-0.5 * ((x - center1) / sigma)**2)
        g2 = amp2 * np.exp(-0.5 * ((x - center2) / sigma)**2)
        return g1 + g2 + background
    
    # Extract and clean data
    mask = (spectrum_table['WAVE'] >= wavelength_min) & (spectrum_table['WAVE'] <= wavelength_max)
    x = spectrum_table['WAVE'][mask]
    y = spectrum_table['FLUX'][mask]
    
    # Get errors if available
    if 'ERROR' in spectrum_table.colnames:
        errors = spectrum_table['ERROR'][mask]
    else:
        errors = None
    
    # Filter for finite values
    finite_mask = np.isfinite(y)
    if errors is not None:
        finite_mask = finite_mask & np.isfinite(errors)
    
    x = x[finite_mask]
    y = y[finite_mask]
    if errors is not None:
        errors = errors[finite_mask]
        # Ensure errors are not too small to prevent numerical issues
        min_error_value = np.max(np.abs(y)) * 1e-6
        errors = np.maximum(errors, min_error_value)
    
    # Initial parameter estimates
    dx = x[1] - x[0]
    total_flux = np.sum(y * dx)
    init_flux1 = total_flux / 2
    init_flux2 = total_flux / 2
    init_background = np.median(y)
    
    # Set up lmfit model
    gmodel = Model(double_gaussian)
    params = gmodel.make_params(
        flux1=init_flux1,
        flux2=init_flux2,
        center1=init_wavelength1,
        center2=init_wavelength2,
        fwhm=init_fwhm,
        background=init_background
    )
    
    # Set parameter bounds
    params['fwhm'].min = 0.1 * init_fwhm
    params['fwhm'].max = 5.0 * init_fwhm
    params['center1'].min = wavelength_min
    params['center1'].max = wavelength_max
    params['center2'].min = wavelength_min
    params['center2'].max = wavelength_max
    
    # Ensure center1 < center2 to avoid swapping during fitting
    if init_wavelength1 < init_wavelength2:
        params['center1'].max = init_wavelength2
        params['center2'].min = init_wavelength1
    else:
        params['center1'].max = init_wavelength1
        params['center2'].min = init_wavelength2
    
    # Perform the fit
    if errors is not None:
        result = gmodel.fit(y, params, x=x, weights=1.0/errors)  # Using error as weights, not squared
    else:
        result = gmodel.fit(y, params, x=x)
    
    # Extract fitted parameters and uncertainties
    flux1 = result.params['flux1'].value
    flux1_err = result.params['flux1'].stderr if result.params['flux1'].stderr is not None else np.nan
    
    flux2 = result.params['flux2'].value
    flux2_err = result.params['flux2'].stderr if result.params['flux2'].stderr is not None else np.nan
    
    wave1 = result.params['center1'].value
    wave1_err = result.params['center1'].stderr if result.params['center1'].stderr is not None else np.nan
    
    wave2 = result.params['center2'].value
    wave2_err = result.params['center2'].stderr if result.params['center2'].stderr is not None else np.nan
    
    fwhm = result.params['fwhm'].value
    fwhm_err = result.params['fwhm'].stderr if result.params['fwhm'].stderr is not None else np.nan
    
    background = result.params['background'].value
    back_err = result.params['background'].stderr if result.params['background'].stderr is not None else np.nan
    
    # Calculate statistics
    y_fit = result.best_fit
    rmse = np.sqrt(np.mean((y - y_fit) ** 2))
    
    # Use lmfit's built-in reduced chi-squared calculation
    reduced_chi2 = result.redchi
    
    # Create output tables
    fit_table = Table([x, y, y_fit], names=('WAVE', 'FLUX', 'Fit'))
    if errors is not None:
        fit_table['ERROR'] = errors
    
    # Define column names using the same naming convention as the original code
    xflux = f'flux_{line}_a'
    xwave = f'wave_{line}_a'
    xfwhm = f'fwhm_{line}_a'
    
    exflux = f'eflux_{line}_a'
    exwave = f'ewave_{line}_a'
    exfwhm = f'efwhm_{line}_a'
    
    yflux = f'flux_{line}_b'
    ywave = f'wave_{line}_b'
    yfwhm = f'fwhm_{line}_b'
    
    eyflux = f'eflux_{line}_b'
    eywave = f'ewave_{line}_b'
    eyfwhm = f'efwhm_{line}_b'
    
    xback = f'back_{line}_ab'
    eback = f'eback_{line}_ab'
    xrmse = f'rmse_{line}_ab'
    xchi2 = f'chi2_{line}_ab'
    
    col_names = [xflux, exflux, yflux, eyflux, xwave, exwave, ywave, eywave, 
                xfwhm, exfwhm, xback, eback, xrmse, xchi2]
    
    qtab = Table(names=col_names)
    qtab.add_row([flux1, flux1_err, flux2, flux2_err, wave1, wave1_err, wave2, wave2_err, 
                fwhm, fwhm_err, background, back_err, rmse, reduced_chi2])
    
    return qtab, fit_table

def save_fit(line='oi',ftab=None,xdir='Gauss_dir'):
    '''
    Save the fit to a directory
    '''
    # print('saving to  %s/%s.txt' % (xdir,line))
    if os.path.isdir(xdir)==False:
        os.makedirs(xdir,exist_ok=True)
    ftab.write('%s/%s.txt' % (xdir,line), format='ascii.fixed_width_two_line',overwrite=True)

def clean(xdir='Gauss_dir',wild='*.txt'):
    '''
    Clean a certain group of files from a directory
    '''
    xfiles=glob('%s/%s' % (xdir,wild))
    if len(xfiles)>0:
        for one in xfiles:
            os.remove(one)
    return

def do_one(spectrum_table,vel=0.,xplot=False,outroot=''):
    '''
    Completely process a single spectrum
    '''

    # First make sure all of the files in directory containing the fits are removed

    clean()



    zz=1.+ (vel/3e5)

    records=[]
    try:
        results,xspec=fit_double_gaussian_to_spectrum(spectrum_table, line='oii',init_wavelength1=zz*3726.092, init_wavelength2=zz*3729.875, init_fwhm=1, wavelength_min=zz*3717, wavelength_max=zz*3737)
        records.append(results)
        if xplot:
            save_fit('oii',xspec)
    except Exception as e:
        print(f"Fitting OII; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")

    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='hb',init_wavelength=zz*4861.325, init_fwhm=1., wavelength_min=zz*4855, wavelength_max=zz*4870)
        records.append(results)
        if xplot:
            save_fit('hb',xspec)
    except Exception as e:
        print(f"Fitting Hb; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='oiii_a',init_wavelength=zz*4958.91, init_fwhm=1., wavelength_min=zz*4949, wavelength_max=zz*4969)
        records.append(results)
        if xplot:
            save_fit('oiii_a',xspec)
    except Exception as e:
        print(f"Fitting [OIII]4959; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='oiii_b', init_wavelength=zz*5006.843, init_fwhm=1., wavelength_min=zz*4997, wavelength_max=zz*5017)
        if xplot:
            save_fit('oiii_b',xspec)
        records.append(results)
    except Exception as e:
        print(f"Fitting [OIII]5007; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")

    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='ha', init_wavelength=zz*6562.8, init_fwhm=1., wavelength_min=zz*6555, wavelength_max=zz*6570)
        if xplot:
            save_fit('ha',xspec)
        records.append(results)
    except Exception as e:
        print(f"Fitting Ha An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")

    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='nii_a',init_wavelength=zz*6548.04, init_fwhm=1., wavelength_min=zz*6538, wavelength_max=zz*6558)
        records.append(results)
        if xplot:
            save_fit('nii_a',xspec)
    except Exception as e:
        print(f"Fitting [NII]6548:  An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='nii_b',init_wavelength=zz*6583.46, init_fwhm=1., wavelength_min=zz*6574, wavelength_max=zz*6594)
        records.append(results)
        if xplot:
            save_fit('nii_b',xspec)
    except Exception as e:
        print(f"Fitting [NII]6584:  An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
            
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='sii_a',init_wavelength=zz*6716.44, init_fwhm=1., wavelength_min=zz*6706, wavelength_max=zz*6726)
        records.append(results)
        if xplot:
            save_fit('sii_a',xspec)
    except Exception as e:
        print(f"Fitting [SII]6716:  An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='sii_b',init_wavelength=zz*6730.81, init_fwhm=1., wavelength_min=zz*6721, wavelength_max=zz*6741)
        records.append(results)
        if xplot:
            save_fit('sii_b',xspec)
    except Exception as e:
        print(f"Fitting [SII]6731:  An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    

    #250110 - Add more lines
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='heii',init_wavelength=zz*4685.71, init_fwhm=1., wavelength_min=zz*4666, wavelength_max=zz*4706)
        records.append(results)
        if xplot:
            save_fit('heii',xspec)
    except Exception as e:
        print(f"Fitting heii; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='hei',init_wavelength=zz*5875.6, init_fwhm=1., wavelength_min=zz*5856, wavelength_max=zz*5886)
        records.append(results)
        if xplot:
            save_fit('hei',xspec)
    except Exception as e:
        print(f"Fitting hei; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='oiii_4363',init_wavelength=zz*4363.21, init_fwhm=1., wavelength_min=zz*4343, wavelength_max=zz*4383)
        records.append(results)
        if xplot:
            save_fit('oiii_4363',xspec)
    except Exception as e:
        print(f"Fitting oiii_4363; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    

    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='siii_a',init_wavelength=zz*9068.6, init_fwhm=1., wavelength_min=zz*9055, wavelength_max=zz*9090)
        records.append(results)
        if xplot:
            save_fit('siii_a',xspec)
    except Exception as e:
        print(f"Fitting siii_a; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='siii_b',init_wavelength=zz*9530.6, init_fwhm=1., wavelength_min=zz*9525, wavelength_max=zz*9545)
        records.append(results)
        if xplot:
            save_fit('siii_b',xspec)
    except Exception as e:
        print(f"Fitting siii_b; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    

    try:
        ztab=hstack(records)
        # return ztab
    except:
        print('Nothing fit for this spectrum')
        return []


    if xplot and outroot!='':
        plot_all(title=outroot)
        plt.savefig('Gauss_dir/%s.png' % outroot)

    return ztab
    
def do_individual(filenames,vel,stype,outname,xplot=True):
    '''
    This is to process individual spectra from an astropy table
    containing a WAVE and FLUX column

    This routine does read the spectra, and so one
    can modify the what is read in at this point.
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
            if stype=='SOURCE':
                xtab['FLUX']=xtab['SOURCE_FLUX']
                xtab['ERROR']=xtab['SOURCE_ERROR']
            elif stype=='BACK':
                xtab['FLUX']=xtab['BACK_FLUX']
                xtab['ERROR']=xtab['BACK_ERROR']
            elif stype!='':
                print('Error: Unkown spectrum type: %s' % stype)
                raise ValueError



            if np.nanmedian(xtab['FLUX'])<1:
                xtab['FLUX']*=1e16
                xtab['ERROR']*=1e16
            efactor=2.25
            xtab['ERROR']/=efactor

            results=do_one(xtab,vel,xplot)
        
            word=one_file.split('/')
            root=word[-1]
            root=root.replace('.txt','')
            # print(results)
            results['Object']=root
            xresults.append(results)
            xgood.append(one_file)
            word=one_file.split('/')
            proot=word[-1].replace('.txt','')
            if stype!='':
                proot='%s.%s' % (proot,stype)
            plot_name='%s/%s.png' % (plot_dir,proot)
            print('Plot created: ',plot_name)
            plot_all(title=proot)
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
    formatted_date = current_date.strftime("%y%m%d")
    if outname=='':
        outname=formatted_date


    if len(filenames)==1:
        if stype!='':
            outname='Gauss_%s.%s.txt' % (root,stype)
        else:
            outname='Gauss_%s.txt' % root
        ftab.write(outname,format='ascii.fixed_width_two_line',overwrite=True)
    else:
        if stype!='':
            outname='Gauss_%s.%s.txt' % (outname,stype)
        else:
            outname='Gauss_%s.txt' % outname
        ftab.write(outname,format='ascii.fixed_width_two_line',overwrite=True)

    return ftab



    
   
    


def scifib(xtab,select='all',telescope=''):
    '''
    Select good fibers from a telescope, of a spefic
    type or all from a telescope from the slitmap table
    of a calbrated file
    '''
    try:
        ztab=xtab[xtab['fibstatus']==0]
    except:
        print('This table does not contain fibstatus, returning all rows')
        return xtab
    if select=='all':
        ztab=ztab[ztab['targettype']!='standard']
    else:
        ztab=ztab[ztab['targettype']==select]

    if telescope!='' and telescope!='all':
        ztab=ztab[ztab['telescope']==telescope]


    print('Found %d fibers' % len(ztab))
    return ztab

def check_for_nan(flux,max_frac=0.5):
    '''
    Check an array for nans. Some number of which are allowed
    '''
    unique_elements, counts = np.unique(np.isnan(flux), return_counts=True)
    if True in unique_elements:
        nan_count = counts[unique_elements == True][0]
    else:
        return False

    if nan_count>max_frac*len(flux):
        return  True
    else:
        return False


def do_all(filename='data/lvmSFrame-00009088.fits',vel=0.0,outname='',xplot=False):
    '''
    Do all of the spectra in a rss fits file.  This multiplies everything
    by 1e16 and corrects the errors as well
    '''
    try:
        x=fits.open(filename)
    except:
        print('Error: Could not open ',filename)
        return

    wave=x['WAVE'].data
    flux=x['FLUX'].data*1e16
    # error=1./np.sqrt(x['IVAR'].data)*1e16
    with np.errstate(divide='ignore', invalid='ignore'):
        error = 1. / np.sqrt(x['IVAR'].data) * 1e16
        error = np.where(np.isfinite(error), error, 1e30)

    efactor=2.25
    error/=efactor

        
    # Now get the good fibers from the science telecsope
    slittab=Table(x['SLITMAP'].data)
    good=scifib(slittab,'science')
    # print(len(good))
    # print(good)
    results=[]
    records=[]
    for i in range(len(good)):
        j=good['fiberid'][i]-1
        one_spec=Table([wave,flux[j],error[j]],names=['WAVE','FLUX','ERROR'])
        if check_for_nan(flux[j])==False:
            rtab=do_one(spectrum_table=one_spec,vel=vel,xplot=xplot,outroot='Fib%04d' % good['fiberid'][i])
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
        outname=outname.replace('.fits','.gauss.txt')
    else:
        outname=outname+'.gauss.txt'

    columns=results.colnames
    for one in columns:
        if one.count('flux'):
            results[one].format='.3e'
        elif one.count('wave'):
            results[one].format='.3f'
        elif one.count('fwhm'):
            results[one].format='.3f'
        elif one.count('rmse'):
            results[one].format='.3e'
        elif one.count('back'):
            results[one].format='.3e'
        elif one.count('chi'):
            results[one].format='.3f'
        elif one.count('ra'):
            results[one].format='.3f'
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
    xx=xx[xx['eflux_oiii_b']!=999.]
    print('Of these, %d spectra also have o3 measured' % len(xx))
    xx=xx[xx['eflux_oii_b']!=999.]
    print('Of these, %d spectra also have o2 measured' % len(xx))
    print('Measured means having an error calculated')

def plot_one(qtab,name='ha'):
    '''
    Make a very simple comparison of the quality of a fit
    '''
    plt.plot(qtab['WAVE'],qtab['FLUX'],label=name)
    plt.plot(qtab['WAVE'],qtab['Fit'])
    plt.plot(qtab['WAVE'],qtab['ERROR'],'.',alpha=0.5)
    plt.legend()
        
def plot_all(qdir='Gauss_dir',title=None):
    '''
    Read in a bunch of fits and plot them all in one big
    plot, assuming all of the fits have been saved to a directoryt

    This produces a single plot of all of the files in the directory
    so one should make sure the directory is emptied at the right times

    This does not save the plot


    '''

    files=glob('%s/*.txt' % qdir)
    files=np.sort(files)
    if len(files)<5:
        nq=2
    elif len(files)<10:
        nq=3
    elif len(files)<17:
        nq=4
    elif len(files)<26:
        nq=5
    elif len(files)<37:
        nq=7
    else:
        print('Too many files to plot at once')
        return
    plt.figure(1,(12,12))
    plt.clf()
    i=0
    while i<len(files):
        plt.subplot(nq,nq,i+1)
        one_file=files[i]
        xtab=ascii.read(one_file)
        word=one_file.split('/')
        name=word[-1].replace('.txt','')
        plot_one(xtab,name)
        i+=1
    plt.tight_layout()

    # Add an overall title if provided
    if title:
        plt.suptitle(title, fontsize=16, y=0.98)
        # Adjust the layout to make room for the title
        plt.subplots_adjust(top=0.92)

    return 

def steer(argv):
    '''
    usage: lvm_gaussfit.py --h] [-lmc] [-smc] [-v vel] [-stype SOURCE][-out root] filename
    '''
    filename=''
    lmc=262.
    smc=146.
    outname=''
    fitsfiles=[]
    specfiles=[]
    xplot=False
    stype=''

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
        elif argv[i]=='-stype':
            i+=1
            stype=argv[i]
        elif argv[i]=='-plot':
            xplot=True
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
        results=do_all(one_file,vel,outname,xplot)
        analyze(results)

    if len(specfiles)>0:
        do_individual(specfiles,vel,stype, outname)

    print('Errors have been decreased by a factor of 2.25; this should be removed after a new processing')

    return




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
