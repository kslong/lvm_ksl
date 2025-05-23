#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Extract fluxes etc from a standard set of emission lines from
the scince fibers of a sky-subtracted LVM exposure


Command line usage (if any):

    usage: lvm_gaussfit.py filename

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




from astropy.io import fits, ascii
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt



import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.modeling import models, fitting

def fit_braod_narrow_to_spectrum(spectrum_table, line, init_wavelength, init_fwhm1, init_fwhm2, wavelength_min, wavelength_max, plot=False):
    # Extract WAVE and FLUX columns from the table
    x = spectrum_table['WAVE']
    y = spectrum_table['FLUX']


    # Limit the wavelengths over which the fit is carried out
    mask = (x >= wavelength_min) & (x <= wavelength_max)
    x_limited = x[mask]
    y_limited = y[mask]

    finite_mask=np.isfinite(y_limited)
    x_limited=x_limited[finite_mask]
    y_limited=y_limited[finite_mask]


    # Initial guesses for Gaussian parameters
    init_amplitude1 = (np.max(y_limited) - np.min(y_limited)) / 2
    init_amplitude2 = init_amplitude1 / 2
    init_mean1 = init_wavelength
    init_mean2 = init_wavelength + init_fwhm1  # Initial guess for the second peak slightly offset
    init_stddev1 = init_fwhm1 / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to standard deviation
    init_stddev2 = init_fwhm2 / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to standard deviation
    init_constant = np.min(y_limited)

    # Define the model: Two Gaussians + constant
    gauss1 = models.Gaussian1D(amplitude=init_amplitude1, mean=init_mean1, stddev=init_stddev1)
    gauss2 = models.Gaussian1D(amplitude=init_amplitude2, mean=init_mean2, stddev=init_stddev2)
    const = models.Const1D(amplitude=init_constant)
    model = gauss1 + gauss2 + const

    # Perform the first fit
    fitter = fitting.LevMarLSQFitter()
    fitted_model = fitter(model, x_limited, y_limited)

    param_errors = None
    if fitter.fit_info['param_cov'] is not None:
        try:
            # Calculate errors as square root of diagonal elements of covariance matrix
            param_errors = np.sqrt(np.diag(fitter.fit_info['param_cov']))
            # print('gotcha ',param_errors)
            amp_err1=param_errors[0]
            wave_err1=param_errors[1]
            fwhm_err1=param_errors[2]
            amp_err2=param_errors[4]
            wave_err2=param_errors[3]
            fwhm_err2=param_errors[4]
        except Exception as e:
            print(f"Error calculating parameter errors: {e}")
            param_errors = None
            amp_err1=wave_err1=fwhm_err1=999.
            amp_err2=wave_err2=fwhm_err2=999.
    else:
        print("Covariance matrix is None, cannot calculate parameter errors.")
        amp_err1=wave_err1=fwhm_err1=999.
        amp_err2=wave_err2=fwhm_err2=999.       
    

    
    # Extract fit parameters
    amplitude1 = fitted_model[0].amplitude.value
    mean1 = fitted_model[0].mean.value
    stddev1 = fitted_model[0].stddev.value 
    amplitude2 = fitted_model[1].amplitude.value
    mean2 = fitted_model[1].mean.value
    stddev2= fitted_model[1].stddev.value  # Shared standard deviation
    background = fitted_model[2].amplitude.value
    
    # Convert stddev back to FWHM
    fwhm1 = stddev1 * 2.355
    fwhm2 = stddev2 * 2.355
    

    gauss_one = models.Gaussian1D(amplitude=fitted_model.amplitude_0.value, mean=fitted_model.mean_0.value, stddev=fitted_model.stddev_0.value)
    gauss_two = models.Gaussian1D(amplitude=fitted_model.amplitude_1.value, mean=fitted_model.mean_1.value, stddev=fitted_model.stddev_1.value)
    const_fin = models.Const1D(amplitude=fitted_model.amplitude_2.value)
    zgauss_1=gauss_one(x_limited)
    zgauss_2=gauss_two(x_limited)
    zconst=const_fin(x_limited)

    # Plot the data and the first fit if requested
    if plot:
        # Print the fitted parameters for the first fit
        print(f"Fitted amplitude 1: {fitted_model.amplitude_0.value}")
        print(f"Fitted mean 1: {fitted_model.mean_0.value}")
        print(f"Fitted stddev 1: {fitted_model.stddev_0.value}")
        print(f"Fitted amplitude 2: {fitted_model.amplitude_1.value}")
        print(f"Fitted mean 2: {fitted_model.mean_1.value}")
        print(f"Fitted stddev 2: {fitted_model.stddev_1.value}")
        print(f"Fitted constant: {fitted_model.amplitude_2.value}")

        plt.figure()
        plt.plot(x_limited, y_limited, 'ko', label='Data')
        plt.plot(x_limited, fitted_model(x_limited), label='Fit')
        plt.plot(x_limited,  y_limited-fitted_model(x_limited), label='delta')
        plt.plot(x_limited,zgauss_1+zconst,label='Narrow')
        plt.plot(x_limited,zgauss_2+zconst,label='Broad')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.legend()
        plt.title(f'Two Gaussian Fit to Line {line}')
        plt.show()


    # Create a new table with WAVE, FLUX, FIT1, and FIT2
    rmse = np.sqrt(np.mean(y_limited - fitted_model(x_limited))**2)
    fit1 = fitted_model(x_limited)
 
    result_table = Table([x_limited, y_limited, fit1], names=('WAVE', 'FLUX', 'FIT'))



    xflux='flux_%s_a' % line
    xwave='wave_%s_a' % line
    xfwhm='fwhm_%s_a' % line
    
    exflux='eflux_%s_a' % line
    exwave='ewave_%s_a' % line
    exfwhm='efwhm_%s_a' % line
        
    xback='back_%s_a' % line
    xrmse='rmse_%s_a' % line    

    yflux='flux_%s_b' % line
    ywave='wave_%s_b' % line
    yfwhm='fwhm_%s_b' % line

    eyflux='eflux_%s_b' % line
    eywave='ewave_%s_b' % line
    eyfwhm='efwhm_%s_b' % line
   
    
    yback='back_%s_b' % line
    yrmse='rmse_%s_b' % line  

    qtab=Table(names=[xflux,exflux,xwave,exwave,xfwhm,exfwhm,xback,xrmse,yflux,eyflux,ywave,eywave,yfwhm,eyfwhm, yback,yrmse])

    for one_col in qtab.colnames:
        qtab[one_col].format='.5e'


    record=[amplitude1,amp_err1,mean1,wave_err1,fwhm1,fwhm_err1,background,rmse,amplitude2,amp_err2,mean2,wave_err2,fwhm2,fwhm_err2,background,rmse]
    qtab.add_row(record)
    return qtab,result_table






def fit_gaussian_to_spectrum(spectrum_table, line, init_wavelength, init_fwhm, wavelength_min, wavelength_max, plot=False):
    """
    Fit a Gaussian with a constant background to a spectral line in a given spectrum table with specified initial parameters and wavelength limits.
    
    Parameters:
    spectrum_table (astropy.table.Table): Table with 'Wavelength' and 'Flux' columns.
    init_wavelength (float): Initial guess for the mean of the Gaussian.
    init_fwhm (float): Initial guess for the FWHM of the Gaussian.
    wavelength_min (float): Minimum wavelength to include in the fitting.
    wavelength_max (float): Maximum wavelength to include in the fitting.
    plot (bool): If True, plot the data and the fit.
    
    Returns:
    amplitude (float): Amplitude of the fitted Gaussian.
    mean (float): Mean of the fitted Gaussian.
    fwhm (float): FWHM of the fitted Gaussian.
    background (float): Amplitude of the constant background.
    rmse (float): Root Mean Square Error between the fit and the data.
    fit_table (astropy.table.Table): Table containing the wavelengths of interest and the fitted values.
    """
    # Extract data from the table within the specified wavelength limits
    mask = (spectrum_table['WAVE'] >= wavelength_min) & (spectrum_table['WAVE'] <= wavelength_max)
    x = spectrum_table['WAVE'][mask]
    y = spectrum_table['FLUX'][mask]

    finite_mask=np.isfinite(y)
    x=x[finite_mask]
    y=y[finite_mask]

    # Convert FWHM to standard deviation: FWHM = 2.355 * sigma
    init_stddev = init_fwhm / 2.355
    
    # Initialize the Gaussian model with the provided initial guesses
    g_init = models.Gaussian1D(amplitude=max(y), mean=init_wavelength, stddev=init_stddev)
    
    # Initialize the constant background model
    const_init = models.Const1D(amplitude=np.median(y))
    
    # Combine the Gaussian and constant models
    combined_init = g_init + const_init
    
    # Initialize the fitter
    fit_g = fitting.LevMarLSQFitter(calc_uncertainties=True)
    
    # Perform the fit
    combined_fit = fit_g(combined_init, x, y)
    param_errors = None
    if fit_g.fit_info['param_cov'] is not None:
        try:
            # Calculate errors as square root of diagonal elements of covariance matrix
            param_errors = np.sqrt(np.diag(fit_g.fit_info['param_cov']))
            # print('gotcha ',param_errors)
            amp_err=param_errors[0]
            wave_err=param_errors[1]
            fwhm_err=param_errors[2]
        except Exception as e:
            print(f"Error calculating parameter errors: {e}")
            param_errors = None
            amp_err=wave_err=fwhm_err=999.
    else:
        print("Covariance matrix is None, cannot calculate parameter errors.")
        amp_err=wave_err=fwhm_err=999.
        
    # Extract fit parameters
    amplitude = combined_fit[0].amplitude.value
    mean = combined_fit[0].mean.value
    stddev = combined_fit[0].stddev.value
    background = combined_fit[1].amplitude.value
    
    # Convert stddev back to FWHM
    fwhm = stddev * 2.355
    
    # Calculate RMSE
    y_fit = combined_fit(x)
    rmse = np.sqrt(np.mean((y - y_fit)**2))
    
    # Create a new table with the wavelengths of interest and the fitted values
    fit_table = Table([x, y,y_fit], names=('WAVE', 'FLUX','Fit'))
    
    # Plot the results if requested
    if plot:
        plt.figure()
        plt.plot(x, y, 'ko', label='Data')
        plt.plot(x, y_fit, 'r-', label='Gaussian + Background Fit')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.legend()
        plt.show()

    xflux='flux_%s' % line
    xwave='wave_%s' % line
    xfwhm='fwhm_%s' % line
    
    eflux='eflux_%s' % line
    ewave='ewave_%s' % line
    efwhm='efwhm_%s' % line
    
    xback='back_%s' % line
    xrmse='rmse_%s' % line

    qtab=Table(names=[xflux,eflux,xwave,ewave,xfwhm,efwhm,xback,xrmse])
    qtab.add_row([amplitude, amp_err, mean, wave_err, fwhm, fwhm_err, background, rmse])
    # qtab.info()
    return qtab, fit_table


def fit_double_gaussian_to_spectrum(spectrum_table, line, init_wavelength1, init_wavelength2, init_fwhm, wavelength_min, wavelength_max, plot=False):
    """
    Fit two Gaussians with a constant background to a spectral line in a given spectrum table with specified initial parameters and wavelength limits.
    
    Parameters:
    spectrum_table (astropy.table.Table): Table with 'WAVE' and 'FLUX' columns.
    init_wavelength1 (float): Initial guess for the mean of the first Gaussian.
    init_wavelength2 (float): Initial guess for the mean of the second Gaussian.
    init_fwhm (float): Initial guess for the FWHM of the Gaussians.
    wavelength_min (float): Minimum wavelength to include in the fitting.
    wavelength_max (float): Maximum wavelength to include in the fitting.
    plot (bool): If True, plot the data and the fit.
    
    Returns:
    amplitude1 (float): Amplitude of the first fitted Gaussian.
    mean1 (float): Mean of the first fitted Gaussian.
    amplitude2 (float): Amplitude of the second fitted Gaussian.
    mean2 (float): Mean of the second fitted Gaussian.
    fwhm (float): FWHM of the fitted Gaussians.
    background (float): Amplitude of the constant background.
    rmse (float): Root Mean Square Error between the fit and the data.
    fit_table (astropy.table.Table): Table containing the wavelengths of interest and the fitted values.
    """
    # Extract data from the table within the specified wavelength limits
    mask = (spectrum_table['WAVE'] >= wavelength_min) & (spectrum_table['WAVE'] <= wavelength_max)
    x = spectrum_table['WAVE'][mask]
    y = spectrum_table['FLUX'][mask]

    finite_mask=np.isfinite(y)
    xx=x[finite_mask]
    yy=y[finite_mask]

    # Extract side regions for background estimation
    side_region_mask = ((spectrum_table['WAVE'] >= wavelength_min - (wavelength_max - wavelength_min)) & 
                        (spectrum_table['WAVE'] < wavelength_min)) | \
                       ((spectrum_table['WAVE'] > wavelength_max) & 
                        (spectrum_table['WAVE'] <= wavelength_max + (wavelength_max - wavelength_min)))
    y_side = spectrum_table['FLUX'][side_region_mask]
    
    # Estimate background as the median flux of the side regions
    background = np.nanmedian(y_side)
    
    # Subtract the estimated background from the flux
    yy = yy - background
    
    # Convert FWHM to standard deviation: FWHM = 2.355 * sigma
    init_stddev = init_fwhm / 2.355
    
    # Initialize the Gaussian models with the provided initial guesses
    g1_init = models.Gaussian1D(amplitude=max(yy), mean=init_wavelength1, stddev=init_stddev)
    g2_init = models.Gaussian1D(amplitude=max(yy), mean=init_wavelength2, stddev=init_stddev)
    
    # Tie the standard deviations of both Gaussians to the same value
    g2_init.stddev.tied = lambda model: model.stddev_0
    
    # Combine the two Gaussians (background is already subtracted)
    combined_init = g1_init + g2_init
    
    # Initialize the fitter
    fit_g = fitting.LevMarLSQFitter(calc_uncertainties=True)
    
    # Perform the fit
    combined_fit = fit_g(combined_init, xx, yy)

    param_errors = None
    if fit_g.fit_info['param_cov'] is not None:
        try:
            # Calculate errors as square root of diagonal elements of covariance matrix
            param_errors = np.sqrt(np.diag(fit_g.fit_info['param_cov']))
            # print('gotcha ',param_errors)
            amp_err1=param_errors[0]
            wave_err1=param_errors[1]
            fwhm_err1=param_errors[2]
            amp_err2=param_errors[4]
            wave_err2=param_errors[3]
            fwhm_err2=param_errors[4]
        except Exception as e:
            print(f"Error calculating parameter errors: {e}")
            param_errors = None
            amp_err1=wave_err1=fwhm_err1=999.
            amp_err2=wave_err2=fwhm_err2=999.
    else:
        print("Covariance matrix is None, cannot calculate parameter errors.")
        amp_err1=wave_err1=fwhm_err1=999.
        amp_err2=wave_err2=fwhm_err2=999.       
    
    # Extract fit parameters
    amplitude1 = combined_fit[0].amplitude.value
    mean1 = combined_fit[0].mean.value
    amplitude2 = combined_fit[1].amplitude.value
    mean2 = combined_fit[1].mean.value
    stddev = combined_fit[0].stddev.value  # Shared standard deviation
    
    # Convert stddev back to FWHM
    fwhm = stddev * 2.355
    
    # Calculate RMSE
    y_fit = combined_fit(x) + background
    rmse = np.sqrt(np.mean((spectrum_table['FLUX'][mask] - y_fit)**2))
    
    # Create a new table with the wavelengths of interest and the fitted values
    fit_table = Table([x, y_fit], names=('WAVE', 'Fit'))
    
    # Plot the results if requested
    if plot:
        plt.figure()
        plt.plot(x, spectrum_table['FLUX'][mask], 'ko', label='Data')
        plt.plot(x, y_fit, 'r-', label='Double Gaussian Fit')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.legend()
        plt.show()

    xflux='flux_%s_a' % line
    xwave='wave_%s_a' % line
    xfwhm='fwhm_%s_a' % line
    
    exflux='eflux_%s_a' % line
    exwave='ewave_%s_a' % line
    exfwhm='efwhm_%s_a' % line
        
    xback='back_%s_a' % line
    xrmse='rmse_%s_a' % line    

    yflux='flux_%s_b' % line
    ywave='wave_%s_b' % line
    yfwhm='fwhm_%s_b' % line

    eyflux='eflux_%s_b' % line
    eywave='ewave_%s_b' % line
    eyfwhm='efwhm_%s_b' % line
   
    
    yback='back_%s_b' % line
    yrmse='rmse_%s_b' % line  

    qtab=Table(names=[xflux,exflux,xwave,exwave,xfwhm,exfwhm,xback,xrmse,yflux,eyflux,ywave,eywave,yfwhm,eyfwhm, yback,yrmse])

    for one_col in qtab.colnames:
        qtab[one_col].format='.5e'

    record=[amplitude1,amp_err1,mean1,wave_err1,fwhm,fwhm_err1,background,rmse,amplitude2,amp_err2,mean2,wave_err2,fwhm,fwhm_err2,background,rmse]
    qtab.add_row(record)
    
    return qtab, fit_table



def do_one(spectrum_table,xplot=False):

    records=[]
    results,xspec=fit_double_gaussian_to_spectrum(spectrum_table, line='oii',init_wavelength1=3726.092, init_wavelength2=3729.875, init_fwhm=2, wavelength_min=3717, wavelength_max=3737, plot=xplot)
    
    records.append(results)
    
    results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='hb',init_wavelength=4863, init_fwhm=1., wavelength_min=4855, wavelength_max=4870, plot=xplot)
    records.append(results)
    
    results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='oiii_a',init_wavelength=4959, init_fwhm=1., wavelength_min=4949, wavelength_max=4969, plot=xplot)
    records.append(results)
    
    results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='oiii_b', init_wavelength=5007, init_fwhm=1., wavelength_min=4997, wavelength_max=5017, plot=xplot)
    records.append(results)

    
    results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='ha', init_wavelength=6563, init_fwhm=1., wavelength_min=6555, wavelength_max=6570, plot=xplot)
    records.append(results)

    
    results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='nii_a',init_wavelength=6548, init_fwhm=1., wavelength_min=6538, wavelength_max=6558, plot=xplot)
    records.append(results)
    
    results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='nii_b',init_wavelength=6584, init_fwhm=1., wavelength_min=6574, wavelength_max=6594, plot=xplot)
    records.append(results)
            
    results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='sii_a',init_wavelength=6716, init_fwhm=1., wavelength_min=6706, wavelength_max=6726, plot=xplot)
    records.append(results)
    
    results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='sii_b',init_wavelength=6731, init_fwhm=1., wavelength_min=6721, wavelength_max=6741, plot=xplot)
    records.append(results)


    results,xspec=fit_braod_narrow_to_spectrum(spectrum_table, line='ha', 
        init_wavelength=6563, init_fwhm1=1.4, init_fwhm2=7, wavelength_min=6552, wavelength_max=6575, plot=xplot)
    records.append(results)


    results,xspec=fit_braod_narrow_to_spectrum(spectrum_table, line='s6716', 
        init_wavelength=6716, init_fwhm1=1.4, init_fwhm2=7, wavelength_min=6706, wavelength_max=6726, plot=xplot)
    records.append(results)


    results,xspec=fit_braod_narrow_to_spectrum(spectrum_table, line='s6731', 
        init_wavelength=6731, init_fwhm1=1.4, init_fwhm2=7, wavelength_min=6721, wavelength_max=6741, plot=xplot)
    records.append(results)


    ztab=hstack(records)
    # print(ztab)
    return ztab
    


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




def do_all(filename='data/lvmSFrame-00009088.fits'):
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
        rtab=do_one(spectrum_table=one_spec,xplot=False)
        rtab['fiberid']=good['fiberid'][i]
        rtab['ra']=good['ra'][i]
        rtab['dec']=good['dec'][i]

        records.append(rtab)

    results=vstack(records)
    outname=filename.split('/')[-1]
    outname=outname.replace('.fits','.double.txt')

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
        
def doit(filename):
    results=do_all(filename)
    analyze(results)
    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        doit(sys.argv[1])
    else:
        print (__doc__)
