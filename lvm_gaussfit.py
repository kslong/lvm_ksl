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
from astropy.io import ascii, fits
import matplotlib.pyplot as plt
from astropy.table import Table,vstack, hstack
import numpy as np
from glob import glob


from scipy.optimize import curve_fit
import astropy.units as u
from datetime import datetime
from astropy.modeling.parameters import Parameter
from astropy.modeling import models, fitting, Fittable1DModel
from astropy.table import Table


def fit_gaussian_to_spectrum(spectrum_table, line, init_wavelength, init_fwhm, wavelength_min, wavelength_max, n_bootstrap=100):
    """
    Fit a Gaussian with a constant background to a spectral line, parameterized by total flux, wavelength, and FWHM.
    Falls back to bootstrapping if covariance matrix is unavailable.
    
    Parameters:
    spectrum_table (astropy.table.Table): Table with 'WAVE' and 'FLUX' columns
    line (str): Name of the line being fit (for output column names)
    init_wavelength (float): Initial guess for the center wavelength
    init_fwhm (float): Initial guess for the FWHM
    wavelength_min, wavelength_max (float): Wavelength range for fitting
    n_bootstrap (int): Number of bootstrap iterations if needed
    
    Returns:
    qtab (astropy.table.Table): Table with fit results and uncertainties
    fit_table (astropy.table.Table): Table with original and fitted values
    """
    
    class FluxGaussian(Fittable1DModel):
        total_flux = Parameter()
        mean = Parameter()
        stddev = Parameter()
        
        @staticmethod
        def evaluate(x, total_flux, mean, stddev):
            amplitude = total_flux / (stddev * np.sqrt(2 * np.pi))
            return amplitude * np.exp(-0.5 * ((x - mean) / stddev) ** 2)
    
    def bootstrap_fit(x, y, init_model):
        """Perform a single bootstrap fit"""


        indices = np.random.randint(0, len(x), len(x))
        x_resampled = x[indices]
        y_resampled = y[indices]
        
        fitter = fitting.LevMarLSQFitter()
        fitted_model = fitter(init_model.copy(), x_resampled, y_resampled)

        return (
            fitted_model.total_flux_0.value,
            fitted_model.mean_0.value,
            fitted_model.stddev_0.value * 2.355,  # Convert to FWHM
            fitted_model.amplitude_1.value  # Background amplitude
        )
    
    # Extract and clean data
    mask = (spectrum_table['WAVE'] >= wavelength_min) & (spectrum_table['WAVE'] <= wavelength_max)
    x = spectrum_table['WAVE'][mask]
    y = spectrum_table['FLUX'][mask]
    
    finite_mask = np.isfinite(y)
    x = x[finite_mask]
    y = y[finite_mask]
    
    # Initial parameter estimates
    init_stddev = init_fwhm / 2.355
    dx = x[1] - x[0]
    init_flux = np.sum(y * dx)
    
    # Initialize models
    g_init = FluxGaussian(total_flux=init_flux, mean=init_wavelength, stddev=init_stddev)
    const_init = models.Const1D(amplitude=np.median(y))
    combined_init = g_init + const_init
    
    # Perform initial fit
    fitter = fitting.LevMarLSQFitter(calc_uncertainties=True)
    combined_fit = fitter(combined_init, x, y)
    
    # Get best-fit parameters
    flux = combined_fit.total_flux_0.value
    wavelength = combined_fit.mean_0.value
    fwhm = combined_fit.stddev_0.value * 2.355
    background = combined_fit.amplitude_1.value
    
    # print('xtest %s' %  line)
    # bootstrap_params = np.array([bootstrap_fit(x, y, combined_init) for _ in range(n_bootstrap)])
    # print('ztest %s' % line)
    
    # Handle parameter errors
    if fitter.fit_info['param_cov'] is not None:
        try:
            param_errors = np.sqrt(np.diag(fitter.fit_info['param_cov']))
            flux_err, wave_err, fwhm_err, back_err = param_errors
            # print(f"Using covariance matrix for {line} errors")
        except Exception:
            print(f"\nSingle1: Covariance matrix failed for {line}, using bootstrap")
            bootstrap_params = np.array([bootstrap_fit(x, y, combined_init) for _ in range(n_bootstrap)])
            flux_err, wave_err, fwhm_err, back_err = np.std(bootstrap_params, axis=0)
    else:
        print(f"\nSingle2: No covariance matrix for {line}, using bootstrap")
        print('Flux Wave FWHM, Back:      %.3e %6.1f %4.1f %.3e' % (flux,wavelength,fwhm,background))
        bootstrap_params = np.array([bootstrap_fit(x, y, combined_init) for _ in range(n_bootstrap)])
        # print("Bootstrap results:", bootstrap_params.shape)
        # print("First bootstrap fit:", bootstrap_params[0])
        flux_err, wave_err, fwhm_err, back_err = np.std(bootstrap_params, axis=0)
        print('eFlux eWave eFWHM, eBack:  %.3e %6.1f %4.1f %.3e' % (flux_err,wave_err,fwhm_err,back_err))
    
    # Calculate RMSE
    y_fit = combined_fit(x)
    rmse = np.sqrt(np.mean((y - y_fit) ** 2))
    
    # Create output tables
    fit_table = Table([x, y, y_fit], names=('WAVE', 'FLUX', 'Fit'))
    
    col_names = [f'{param}_{line}' for param in 
                 ['flux', 'eflux', 'wave', 'ewave', 'fwhm', 'efwhm', 'back', 'rmse']]
    
    qtab = Table(names=col_names)
    qtab.add_row([flux, flux_err, wavelength, wave_err, fwhm, fwhm_err, background, rmse])
    
    return qtab, fit_table

def fit_double_gaussian_to_spectrum(spectrum_table, line, init_wavelength1, init_wavelength2, init_fwhm, wavelength_min, wavelength_max, n_bootstrap=100):
    """
    Fit two Gaussians with a shared FWHM and a constant background to a spectral region.
    Falls back to bootstrapping if covariance matrix is unavailable.
    
    Parameters:
    spectrum_table (astropy.table.Table): Table with 'WAVE' and 'FLUX' columns
    line (str): Name of the line being fit (for output column names)
    init_wavelength1, init_wavelength2 (float): Initial guesses for the center wavelengths of the two Gaussians
    init_fwhm (float): Initial guess for the shared FWHM
    wavelength_min, wavelength_max (float): Wavelength range for fitting
    n_bootstrap (int): Number of bootstrap iterations if needed
    
    Returns:
    qtab (astropy.table.Table): Table with fit results and uncertainties
    fit_table (astropy.table.Table): Table with original and fitted values
    """
    class SharedFWHMGaussian(Fittable1DModel):
        flux1 = Parameter()
        flux2 = Parameter()
        mean1 = Parameter()
        mean2 = Parameter()
        stddev = Parameter()
        
        @staticmethod
        def evaluate(x, flux1, flux2, mean1, mean2, stddev):
            amplitude1 = flux1 / (stddev * np.sqrt(2 * np.pi))
            amplitude2 = flux2 / (stddev * np.sqrt(2 * np.pi))
            g1 = amplitude1 * np.exp(-0.5 * ((x - mean1) / stddev) ** 2)
            g2 = amplitude2 * np.exp(-0.5 * ((x - mean2) / stddev) ** 2)
            return g1 + g2
    
    def bootstrap_fit(x, y, init_model):
        indices = np.random.randint(0, len(x), len(x))
        x_resampled = x[indices]
        y_resampled = y[indices]
        
        fitter = fitting.LevMarLSQFitter()
        fitted_model = fitter(init_model.copy(), x_resampled, y_resampled)
        
        return (
            fitted_model.flux1_0.value,
            fitted_model.flux2_0.value,
            fitted_model.mean1_0.value,
            fitted_model.mean2_0.value,
            fitted_model.stddev_0.value * 2.355,
            fitted_model.amplitude_1.value
        )
    
    mask = (spectrum_table['WAVE'] >= wavelength_min) & (spectrum_table['WAVE'] <= wavelength_max)
    x = spectrum_table['WAVE'][mask]
    y = spectrum_table['FLUX'][mask]
    
    finite_mask = np.isfinite(y)
    x = x[finite_mask]
    y = y[finite_mask]
    
    init_stddev = init_fwhm / 2.355
    dx = x[1] - x[0]
    init_flux1 = np.sum(y * dx) / 2
    init_flux2 = np.sum(y * dx) / 2
    
    g_init = SharedFWHMGaussian(flux1=init_flux1, flux2=init_flux2,
                                mean1=init_wavelength1, mean2=init_wavelength2,
                                stddev=init_stddev)
    const_init = models.Const1D(amplitude=np.median(y))
    combined_init = g_init + const_init
    
    fitter = fitting.LevMarLSQFitter(calc_uncertainties=True)
    combined_fit = fitter(combined_init, x, y)
    
    flux1 = combined_fit.flux1_0.value
    flux2 = combined_fit.flux2_0.value
    wave1 = combined_fit.mean1_0.value
    wave2 = combined_fit.mean2_0.value
    fwhm = combined_fit.stddev_0.value * 2.355
    background = combined_fit.amplitude_1.value
    
    if fitter.fit_info['param_cov'] is not None:
        try:
            param_errors = np.sqrt(np.diag(fitter.fit_info['param_cov']))
            flux1_err, flux2_err, wave1_err, wave2_err, fwhm_err, back_err = param_errors
            # print(f"Using covariance matrix for {line} errors")
        except Exception:
            print(f"\ndouble1: Covariance matrix failed for {line}, using bootstrap")
            bootstrap_params = np.array([bootstrap_fit(x, y, combined_init) for _ in range(n_bootstrap)])
            print('flux1 flux2 wave1 wave2 fwhm background:       %8.3e %8.3e %8.3f %8.3f %8.3f %8.3e' % (flux1,flux2,wave1,wave2,fwhm,background))
            flux1_err, flux2_err, wave1_err, wave2_err, fwhm_err, back_err = np.std(bootstrap_params, axis=0)
            print('eflux1 eflux2 ewave1 ewave2 efwhm ebackground: %8.3e %8.3e %8.3f %8.3f %8.3f %8.3e' % (flux1_err,flux2_err,wave1_err,wave2_err,fwhm_err,back_err))
    else:
        print(f"\ndouble2: No covariance matrix for {line}, using bootstrap")
        print('flux1 flux2 wave1 wave2 fwhm background:       %8.3e %8.3e %8.3f %8.3f %8.3f %8.3e' % (flux1,flux2,wave1,wave2,fwhm,background))
        bootstrap_params = np.array([bootstrap_fit(x, y, combined_init) for _ in range(n_bootstrap)])
        flux1_err, flux2_err, wave1_err, wave2_err, fwhm_err, back_err = np.std(bootstrap_params, axis=0)
        print('eflux1 eflux2 ewave1 ewave2 efwhm ebackground: %8.3e %8.3e %8.3f %8.3f %8.3f %8.3e' % (flux1_err,flux2_err,wave1_err,wave2_err,fwhm_err,back_err))
    
    y_fit = combined_fit(x)
    rmse = np.sqrt(np.mean((y - y_fit) ** 2))
    
    fit_table = Table([x, y, y_fit], names=('WAVE', 'FLUX', 'Fit'))

    xflux='flux_%s_a' % line
    xwave='wave_%s_a' % line
    xfwhm='fwhm_%s_a' % line

    exflux='eflux_%s_a' % line
    exwave='ewave_%s_a' % line
    exfwhm='efwhm_%s_a' % line

    xback='back_%s_ab' % line
    xrmse='rmse_%s_ab' % line

    yflux='flux_%s_b' % line
    ywave='wave_%s_b' % line
    yfwhm='fwhm_%s_b' % line

    eyflux='eflux_%s_b' % line
    eywave='ewave_%s_b' % line
    eyfwhm='efwhm_%s_b' % line


    yback='back_%s_b' % line
    yrmse='rmse_%s_b' % line
    

    col_names= [xflux, exflux, yflux, eyflux, xwave, exwave, ywave, eywave, xfwhm, exfwhm, xback, xrmse]
    
    qtab = Table(names=col_names)
    qtab.add_row([flux1, flux1_err, flux2, flux2_err, wave1, wave1_err, wave2, wave2_err, fwhm, fwhm_err, background, rmse])
    
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

def do_one(spectrum_table,vel=0.,xplot=False):
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
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='hb',init_wavelength=zz*4861, init_fwhm=1., wavelength_min=zz*4855, wavelength_max=zz*4870)
        records.append(results)
        if xplot:
            save_fit('hb',xspec)
    except Exception as e:
        print(f"Fitting Hb; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='oiii_a',init_wavelength=zz*4959, init_fwhm=1., wavelength_min=zz*4949, wavelength_max=zz*4969)
        records.append(results)
        if xplot:
            save_fit('oiii_a',xspec)
    except Exception as e:
        print(f"Fitting [OIII]4959; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='oiii_b', init_wavelength=zz*5007, init_fwhm=1., wavelength_min=zz*4997, wavelength_max=zz*5017)
        if xplot:
            save_fit('oiii_b',xspec)
        records.append(results)
    except Exception as e:
        print(f"Fitting [OIII]5007; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")

    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='ha', init_wavelength=zz*6563, init_fwhm=1., wavelength_min=zz*6555, wavelength_max=zz*6570)
        if xplot:
            save_fit('ha',xspec)
        records.append(results)
    except Exception as e:
        print(f"Fitting Ha An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")

    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='nii_a',init_wavelength=zz*6548, init_fwhm=1., wavelength_min=zz*6538, wavelength_max=zz*6558)
        records.append(results)
        if xplot:
            save_fit('nii_a',xspec)
    except Exception as e:
        print(f"Fitting [NII]6548:  An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='nii_b',init_wavelength=zz*6584, init_fwhm=1., wavelength_min=zz*6574, wavelength_max=zz*6594)
        records.append(results)
        if xplot:
            save_fit('nii_b',xspec)
    except Exception as e:
        print(f"Fitting [NII]6584:  An exception occurred: {e}")
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
    

    #250110 - Add more lines
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='heii',init_wavelength=zz*4686, init_fwhm=1., wavelength_min=zz*4666, wavelength_max=zz*4706)
        records.append(results)
        if xplot:
            save_fit('heii',xspec)
    except Exception as e:
        print(f"Fitting heii; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='hei',init_wavelength=zz*5876, init_fwhm=1., wavelength_min=zz*5856, wavelength_max=zz*5896)
        records.append(results)
        if xplot:
            save_fit('hei',xspec)
    except Exception as e:
        print(f"Fitting hei; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='oiii_4363',init_wavelength=zz*4363, init_fwhm=1., wavelength_min=zz*4343, wavelength_max=zz*4383)
        records.append(results)
        if xplot:
            save_fit('oiii_4363',xspec)
    except Exception as e:
        print(f"Fitting oiii_4363; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    

    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='siii_a',init_wavelength=zz*9068, init_fwhm=1., wavelength_min=zz*9048, wavelength_max=zz*9088)
        records.append(results)
        if xplot:
            save_fit('siii_a',xspec)
    except Exception as e:
        print(f"Fitting siii_a; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='siii_b',init_wavelength=zz*9531, init_fwhm=1., wavelength_min=zz*9511, wavelength_max=zz*9551)
        records.append(results)
        if xplot:
            save_fit('siii_b',xspec)
    except Exception as e:
        print(f"Fitting siii_b; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    

    try:
        ztab=hstack(records)
        return ztab
    except:
        print('Nothing fit for this spectrum')
        return []

    return ztab
    
def do_individual(filenames,vel,outname,xplot=True):
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
            # print(results)
            results['Object']=root
            xresults.append(results)
            xgood.append(one_file)
            plot_all()
            word=one_file.split('/')
            proot=word[-1].replace('.txt','')
            plot_name='%s/%s.png' % (plot_dir,proot)
            print('Plot created: ',plot_name)
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
            results[one].format='.3e'
        elif one.count('wave'):
            results[one].format='.3f'
        elif one.count('fwhm'):
            results[one].format='.3f'
        elif one.count('rmse'):
            results[one].format='.3e'
        elif one.count('back'):
            results[one].format='.3e'
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
    plt.legend()
        
def plot_all(qdir='Gauss_dir'):
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
    return 

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
