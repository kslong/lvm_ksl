#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This is an experimental routine to attempt to fit a single 
line in lvm spectra to a single and double gaussian
and provide information on whethe a complest fit is required.


Command line usage (if any):

    usage: xgauus.py -h -wmin whatever -wmax whatever file.tx or file.fits  

Description:  

Primary routines:

    doit

Notes:
                                       
History:

250323 ksl Coding begun

'''






import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from lmfit import Model
from lmfit.models import GaussianModel
from scipy import integrate
from astropy.io import ascii, fits
from astropy.table import Table,vstack, hstack
from glob import glob
import astropy.units as u
from datetime import datetime


def sanitize_table(input_table, default_value=-999.99):
    """
    Returns a sanitized copy of the input Astropy Table where:
      - All None, NaN, and inf values are replaced with default_value
      - Object-type columns are converted to float
    """
    sanitized = input_table.copy()

    for name in sanitized.colnames:
        col = sanitized[name]

        # Convert object columns to float if possible
        if col.dtype.kind == 'O':
            # Try converting to float, handling None
            new_data = [float(val) if val not in (None, '') else default_value for val in col]
            sanitized[name] = np.array(new_data, dtype=float)
            continue

        # Handle masked or numeric columns
        if np.issubdtype(col.dtype, np.number):
            # Replace non-finite values (nan, inf)
            bad = ~np.isfinite(col)
            col[bad] = default_value

        # If masked, also handle masked values
        if hasattr(col, 'mask') and np.any(col.mask):
            col[col.mask] = default_value
            col.mask = np.zeros(len(col), dtype=bool)

    return sanitized

def print_and_visualize_results(model_row, x_data, y_data, y_errors, single_result, double_result,label=''):
    """
    Print detailed results and create visualization for model comparison.

    Parameters:
        model_row (astropy.table.Row): Table row containing model comparison results.
        x_data (array-like): X-axis data (e.g., wavelength).
        y_data (array-like): Observed y-axis data.
        y_errors (array-like or None): Errors in y_data.
        single_result (lmfit.ModelResult): Result of single Gaussian model fit.
        double_result (lmfit.ModelResult): Result of double Gaussian model fit.
        label (str): A plot label. Default is empty string.

    Returns:
        None: Creates a matplotlib figure showing the fit comparison.
    """

    model_row=sanitize_table(model_row)
    # Unpack statistical values from the table row
    aic1 = model_row['aic_single'][0]
    bic1 = model_row['bic_single'][0]
    aic2 = model_row['aic_double'][0]
    bic2 = model_row['bic_double'][0]
    f_stat = model_row['f_statistic'][0]
    p_value = model_row['p_value'][0]
    preferred_model = model_row['preferred_model'][0]

    # Prepare parameter dictionaries for printing
    single_params_dict = {
        'Amplitude': f"{model_row['s1_amplitude'][0]:.4f} ± {model_row['s1_amplitude_err'][0]:.4f}",
        'Wave     ': f"{model_row['s1_wave'][0]:.4f} ± {model_row['s1_wave_err'][0]:.4f}",
        'Sigma   ': f"{model_row['s1_sigma'][0]:.4f} ± {model_row['s1_sigma_err'][0]:.4f}"
    }

    double_params_dict = {
        'Component 1 Amplitude': f"{model_row['d1_amplitude'][0]:.4f} ± {model_row['d1_amplitude_err'][0]:.4f}",
        'Component 1      Wave': f"{model_row['d1_wave'][0]:.4f} ± {model_row['d1_wave_err'][0]:.4f}",
        'Component 1     Sigma': f"{model_row['d1_sigma'][0]:.4f} ± {model_row['d1_sigma_err'][0]:.4f}",
        'Component 2 Amplitude': f"{model_row['d2_amplitude'][0]:.4f} ± {model_row['d2_amplitude_err'][0]:.4f}",
        'Component 2      Wave': f"{model_row['d2_wave'][0]:.4f} ± {model_row['d2_wave_err'][0]:.4f}",
        'Component 2     Sigma': f"{model_row['d2_sigma'][0]:.4f} ± {model_row['d2_sigma_err'][0]:.4f}"
    }

    # Print parameter values
    print("=== Single Gaussian Parameters ===")
    for name, value in single_params_dict.items():
        print(f"{name}: {value}")

    print("\n=== Double Gaussian Parameters ===")
    for name, value in double_params_dict.items():
        print(f"{name}: {value}")

    # Print statistical comparison
    print("\n=== Model Comparison ===")
    print(f"Single Gaussian: AIC = {aic1:.4f}, BIC = {bic1:.4f}, Reduced χ² = {model_row['redchi_single'][0]:.4f}")
    print(f"Double Gaussian: AIC = {aic2:.4f}, BIC = {bic2:.4f}, Reduced χ² = {model_row['redchi_double'][0]:.4f}")
    print(f"F-test: F = {f_stat:.4f}, p-value = {p_value:.6f}")
    print(f"Preferred model: {preferred_model}")

    # Print total flux information for both models
    print("\n=== Total Flux Information ===")
    print(f"Preferred model ({preferred_model}):")
    if preferred_model == 'single_gaussian':
        print(f"  Total flux: {model_row['single_total_flux'][0]:.4e} ± {model_row['single_total_flux_err'][0]:.4e} ergs/cm²/s")
    else:
        print(f"  Total flux: {model_row['double_total_flux'][0]:.4e} ± {model_row['double_total_flux_err'][0]:.4e} ergs/cm²/s")

    print("\nSingle Gaussian model:")
    print(f"  Total flux: {model_row['single_total_flux'][0]:.4e} ± {model_row['single_total_flux_err'][0]:.4e} ergs/cm²/s")

    print("\nDouble Gaussian model:")
    print(f"   Total flux: {model_row['double_total_flux'][0]:.4e} ± {model_row['double_total_flux_err'][0]:.4e} ergs/cm²/s")
    print(f"  Component 1: {model_row['d1_flux'][0]:.4e} ± {model_row['d1_flux_err'][0]:.4e} ergs/cm²/s")
    print(f"  Component 2: {model_row['d2_flux'][0]:.4e} ± {model_row['d2_flux_err'][0]:.4e} ergs/cm²/s")

    # Plotting (most of the previous plotting code remains the same)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Plot the data and fits
    ax1.errorbar(x_data, y_data, yerr=y_errors, fmt='o', alpha=0.5, label='Data')
    ax1.plot(x_data, single_result.best_fit, 'r-', label='Single Gaussian')
    ax1.plot(x_data, double_result.best_fit, 'g-', label='Double Gaussian')

    # If double Gaussian is preferred, plot individual components
    if preferred_model == 'double_gaussian':
        # Create component models
        g1 = GaussianModel(prefix='g1_')
        g2 = GaussianModel(prefix='g2_')

        # Set parameters from the fit
        g1_params = g1.make_params()
        g2_params = g2.make_params()

        for param_name, param in double_result.params.items():
            if param_name.startswith('g1_'):
                g1_params[param_name].set(value=param.value)
            elif param_name.startswith('g2_'):
                g2_params[param_name].set(value=param.value)

        # Plot components
        ax1.plot(x_data, g1.eval(g1_params, x=x_data), 'g--', alpha=0.5, label='Component 1')
        ax1.plot(x_data, g2.eval(g2_params, x=x_data), 'g:', alpha=0.5, label='Component 2')

    # Set title
    if label=='':
        label='Model Comparison'

    if preferred_model == 'single_gaussian':
        ax1.set_title('%s : Single Gaussian Preferred' % label)

    else:
        ax1.set_title('%s : Double Gaussian Preferred' % label)

    ax1.legend()
    ax1.set_ylabel('Flux (ergs/cm²/s/Å)')

    # Plot residuals
    ax2.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    ax2.plot(x_data, single_result.residual, 'r-', label='Single Gaussian Residuals')
    ax2.plot(x_data, double_result.residual, 'g-', label='Double Gaussian Residuals')
    ax2.set_xlabel('Wavelength (Å)')
    ax2.set_ylabel('Residuals')
    ax2.legend()

    # Add text with statistical information and total flux
    preferred_flux = (model_row['double_total_flux'][0] if preferred_model == 'double_gaussian' 
                      else model_row['single_total_flux'][0])
    info_text = (
        f"Single Gaussian: AIC={aic1:.2f}, BIC={bic1:.2f}, χ²/dof={model_row['redchi_single'][0]:.2f}\n"
        f"Double Gaussian: AIC={aic2:.2f}, BIC={bic2:.2f}, χ²/dof={model_row['redchi_double'][0]:.2f}\n"
        f"F-test: F={f_stat:.2f}, p={p_value:.4f} ({'' if p_value < 0.05 else 'not '}significant)\n"
        f"Total flux ({preferred_model}): {preferred_flux:.4e} ergs/cm²/s"
    )
    plt.figtext(0.5, 0.01, info_text, ha='center', fontsize=10, bbox={"facecolor":"white", "alpha":0.8, "pad":5})

    plt.tight_layout(rect=[0, 0.05, 1, 0.98])
    return



def compare_models(x_data, y_data, y_errors=None, do_plot=False,label=''):
    """
    Compare single and double Gaussian models and return results as an Astropy Table row
    with distinct prefixes for single and double gaussian parameters
    """
    # Set up weights
    weights = 1/y_errors if y_errors is not None else None
    
    # Single Gaussian model (previous implementation remains the same)
    single_gauss = GaussianModel(prefix='g1_')
    single_params = single_gauss.guess(y_data, x=x_data)
    
    # Add amplitude constraint: must be positive
    single_params['g1_amplitude'].set(min=0.0)
    single_params['g1_sigma'].set(min=0.4,max=2.0)
    single_params['g1_center'].set(min=np.min(x_data),max=np.max(x_data))
    
    # Fit single Gaussian
    single_result = single_gauss.fit(y_data, single_params, x=x_data, weights=weights)
    
    # Double Gaussian model (previous implementation remains the same)
    double_gauss = GaussianModel(prefix='g1_') + GaussianModel(prefix='g2_')
    
    # Guess parameters for double Gaussian
    double_params = double_gauss.make_params()
    
    # First Gaussian similar to single Gaussian result
    for param_name in single_result.params:
        double_params[param_name].set(value=single_result.params[param_name].value)
    
    # Residual-based guess for second Gaussian
    residuals = y_data - single_result.best_fit

    # Find positive residual peaks that are sufficiently separated from the first Gaussian
    g1_center = single_result.params['g1_center'].value
    g1_sigma = single_result.params['g1_sigma'].value
    min_separation = 1.5 * g1_sigma  # Require peaks to be at least 1.5 sigma away from first peak

    # Create a mask for points that are far enough from the first Gaussian
    distance_mask = np.abs(x_data - g1_center) > min_separation

    # Find the highest positive residual peak that satisfies the distance constraint
    masked_residuals = residuals.copy()
    masked_residuals[~distance_mask] = 0  # Zero out residuals too close to first peak

    # Only consider positive residuals - we're looking for a second emission peak
    positive_mask = masked_residuals > 0
    if np.any(positive_mask):
        masked_residuals[~positive_mask] = 0  # Zero out negative residuals
        
        # Find the index of the maximum positive residual outside the exclusion zone
        if np.max(masked_residuals) > 0:
            max_resid_idx = np.argmax(masked_residuals)
        else:
            # Fall back to a default position if no clear second peak
            max_resid_idx = np.argmax(np.abs(residuals))
    else:
        # Fall back if no positive residuals found outside exclusion zone
        max_resid_idx = np.argmax(np.abs(residuals))
    
    # Store initial parameter values for diagnosis
    g2_init_center = x_data[max_resid_idx]
    g2_init_amplitude = max(residuals[max_resid_idx] * 2, 0.1)
    g2_init_sigma = single_result.params['g1_sigma'].value * 0.8
    
    # Set the parameters
    double_params['g2_center'].set(g2_init_center)
    double_params['g2_amplitude'].set(g2_init_amplitude)  # Ensure positive
    double_params['g2_sigma'].set(g2_init_sigma)  # Slightly narrower
    
    # Add positive constraints
    double_params['g1_amplitude'].set(min=0.0)
    double_params['g2_amplitude'].set(min=0.0)
    double_params['g1_sigma'].set(min=0.4,max=2.0)
    double_params['g2_sigma'].set(min=0.4,max=2.0)
    double_params['g1_center'].set(min=np.min(x_data),max=np.max(x_data))
    double_params['g2_center'].set(min=np.min(x_data),max=np.max(x_data))
    
    # Fit double Gaussian
    double_result = double_gauss.fit(y_data, double_params, x=x_data, weights=weights)
    
    
    # Model comparison statistics
    n = len(x_data)
    k1 = len(single_result.params)  # Number of parameters in model 1
    k2 = len(double_result.params)  # Number of parameters in model 2
    
    # Calculate AIC 
    aic1 = n * np.log(single_result.chisqr/n) + 2 * k1
    aic2 = n * np.log(double_result.chisqr/n) + 2 * k2
    
    # Calculate BIC
    bic1 = n * np.log(single_result.chisqr/n) + k1 * np.log(n)
    bic2 = n * np.log(double_result.chisqr/n) + k2 * np.log(n)
    
    # F-test
    dof1 = n - k1
    dof2 = n - k2
    
    if single_result.chisqr < double_result.chisqr:
        f_stat = 0
        p_value = 1.0
    else:
        f_stat = ((single_result.chisqr - double_result.chisqr) / (k2 - k1)) / (double_result.chisqr / dof2)
        p_value = 1 - stats.f.cdf(f_stat, k2 - k1, dof2)
    
    # Calculate flux for single Gaussian model
    amp = single_result.params['g1_amplitude'].value
    amp_err = single_result.params['g1_amplitude'].stderr
    sigma = single_result.params['g1_sigma'].value
    sigma_err = single_result.params['g1_sigma'].stderr
    
    single_total_flux = amp * sigma * np.sqrt(2 * np.pi)
    
    # Error propagation for single Gaussian
    if amp_err is not None and sigma_err is not None:
        single_total_flux_err = np.sqrt(
            (sigma * np.sqrt(2 * np.pi) * amp_err) ** 2 + 
            (amp * np.sqrt(2 * np.pi) * sigma_err) ** 2
        )
    else:
        single_total_flux_err = np.nan
    
    # Calculate flux for double Gaussian model
    amp1 = double_result.params['g1_amplitude'].value
    amp1_err = double_result.params['g1_amplitude'].stderr
    sigma1 = double_result.params['g1_sigma'].value
    sigma1_err = double_result.params['g1_sigma'].stderr
    
    amp2 = double_result.params['g2_amplitude'].value
    amp2_err = double_result.params['g2_amplitude'].stderr
    sigma2 = double_result.params['g2_sigma'].value
    sigma2_err = double_result.params['g2_sigma'].stderr
    
    # Calculate flux for each component of double Gaussian
    def calculate_component_flux(amp, amp_err, sigma, sigma_err):
        """
        Calculate flux and error for a single Gaussian component
        """
        flux = amp * sigma * np.sqrt(2 * np.pi)
        
        if amp_err is not None and sigma_err is not None:
            flux_err = np.sqrt(
                (sigma * np.sqrt(2 * np.pi) * amp_err) ** 2 + 
                (amp * np.sqrt(2 * np.pi) * sigma_err) ** 2
            )
        else:
            flux_err = np.nan
        
        return flux, flux_err
    
    flux1, flux1_err = calculate_component_flux(amp1, amp1_err, sigma1, sigma1_err)
    flux2, flux2_err = calculate_component_flux(amp2, amp2_err, sigma2, sigma2_err)
    
    double_total_flux = flux1 + flux2
    
    # Error propagation for double Gaussian total flux
    if not np.isnan(flux1_err) and not np.isnan(flux2_err):
        double_total_flux_err = np.sqrt(flux1_err**2 + flux2_err**2)
    else:
        double_total_flux_err = np.nan
    
    # Determine preferred model
    preferred_model = 'double_gaussian' if (aic2 < aic1 and bic2 < bic1) or p_value < 0.05 else 'single_gaussian'
    
    # Check if second Gaussian is negligible
    if preferred_model == 'double_gaussian':
        g1_amp = double_result.params['g1_amplitude'].value
        g2_amp = double_result.params['g2_amplitude'].value
        g2_amp_stderr = double_result.params['g2_amplitude'].stderr
        
        if g2_amp < 0.05 * g1_amp or (g2_amp_stderr is not None and g2_amp_stderr > g2_amp):
            preferred_model = 'single_gaussian'
    # Create Astropy Table row
    model_row = Table()
    
    # Model selection metrics
    model_row['preferred_model'] = [preferred_model]
    model_row['aic_single'] = [aic1]
    model_row['aic_double'] = [aic2]
    model_row['bic_single'] = [bic1]
    model_row['bic_double'] = [bic2]
    model_row['f_statistic'] = [f_stat]
    model_row['p_value'] = [p_value]
    
    # Single Gaussian parameters - using 's1' prefix
    model_row['s1_amplitude'] = [single_result.params['g1_amplitude'].value]
    model_row['s1_amplitude_err'] = [single_result.params['g1_amplitude'].stderr]
    model_row['s1_wave'] = [single_result.params['g1_center'].value]
    model_row['s1_wave_err'] = [single_result.params['g1_center'].stderr]
    model_row['s1_sigma'] = [single_result.params['g1_sigma'].value]
    model_row['s1_sigma_err'] = [single_result.params['g1_sigma'].stderr]
    model_row['single_total_flux'] = [single_total_flux]
    model_row['single_total_flux_err'] = [single_total_flux_err]
    model_row['redchi_single'] = [single_result.redchi]
    
    # Double Gaussian parameters - using 'd1' and 'd2' prefixes
    model_row['d1_amplitude'] = [double_result.params['g1_amplitude'].value]
    model_row['d1_amplitude_err'] = [double_result.params['g1_amplitude'].stderr]
    model_row['d1_wave'] = [double_result.params['g1_center'].value]
    model_row['d1_wave_err'] = [double_result.params['g1_center'].stderr]
    model_row['d1_sigma'] = [double_result.params['g1_sigma'].value]
    model_row['d1_sigma_err'] = [double_result.params['g1_sigma'].stderr]
    
    model_row['d2_amplitude'] = [double_result.params['g2_amplitude'].value]
    model_row['d2_amplitude_err'] = [double_result.params['g2_amplitude'].stderr]
    model_row['d2_wave'] = [double_result.params['g2_center'].value]
    model_row['d2_wave_err'] = [double_result.params['g2_center'].stderr]
    model_row['d2_sigma'] = [double_result.params['g2_sigma'].value]
    model_row['d2_sigma_err'] = [double_result.params['g2_sigma'].stderr]
    
    # Add individual component fluxes for double Gaussian
    model_row['d1_flux'] = [flux1]
    model_row['d1_flux_err'] = [flux1_err]
    model_row['d2_flux'] = [flux2]
    model_row['d2_flux_err'] = [flux2_err]
    
    model_row['double_total_flux'] = [double_total_flux]
    model_row['double_total_flux_err'] = [double_total_flux_err]
    model_row['redchi_double'] = [double_result.redchi]
    
    # Add initial parameter values for the second Gaussian fit
    model_row['d2_init_wave'] = [g2_init_center]
    model_row['d2_init_amplitude'] = [g2_init_amplitude]
    model_row['d2_init_sigma'] = [g2_init_sigma]
    
    # Plotting if requested
    if do_plot:
        print_and_visualize_results(model_row, x_data, y_data, y_errors, single_result, double_result,label)
    
    # Return just the table row
    return model_row

# Example usage
def example(do_plot=True):
    # Create synthetic data with two Gaussians
    np.random.seed(42)
    x = np.linspace(-10, 10, 100)
    
    # True parameters
    amp1, cen1, sigma1 = 5.0, -1.5, 1.0
    amp2, cen2, sigma2 = 3.0, 2.0, 1.3
    
    # Create a single or double Gaussian (to test both cases)
    double_peak = True
    
    if double_peak:
        y_true = amp1 * np.exp(-(x-cen1)**2/(2*sigma1**2)) + amp2 * np.exp(-(x-cen2)**2/(2*sigma2**2))
        print("True model: Double Gaussian")
        print(f"  Peak 1: amplitude={amp1}, center={cen1}, sigma={sigma1}")
        print(f"  Peak 2: amplitude={amp2}, center={cen2}, sigma={sigma2}")
    else:
        y_true = amp1 * np.exp(-(x-cen1)**2/(2*sigma1**2))
        print("True model: Single Gaussian")
        print(f"  Peak: amplitude={amp1}, center={cen1}, sigma={sigma1}")
    
    # Add noise
    noise_level = 0.2
    y_noise = y_true + noise_level * np.random.normal(size=len(x))
    y_errors = np.ones_like(x) * noise_level
    
    # Compare models
    results = compare_models(x, y_noise, y_errors,do_plot)
    
    return results

from astropy.stats import sigma_clipped_stats

def do_one(ztab,wmin=6715,wmax=6750,do_plot=True,label=''):
    '''
    Prooces a single spectrum which has already been converted to an
    astropy table, if that was not already the case q

    '''
    if len(ztab)<10:
        print('Houston we have a problem')
        return []

    scale_factor=1e16

    stab=ztab.copy()
    stab=stab[stab['WAVE']>wmin]
    stab=stab[stab['WAVE']<wmax]
    xdata=np.array(stab['WAVE'])
    ydata=np.array(stab['FLUX'])*scale_factor
    # ydata=ydata-np.min(ydata)

    mean_flux, median_flux, std_flux = sigma_clipped_stats(ydata, sigma_upper=1,sigma_lower=3)
    ydata-=median_flux

    error=np.array(stab['ERROR'])*scale_factor

    results=compare_models(xdata,ydata,error,do_plot=do_plot,label=label)
    return results

def do_individual(filenames,wmin=6710,wmax=6750,outname='',xplot=False):
    '''
    This is to process individual spectra from an astropy table
    containing a WAVE and FLUX column
    '''
    xresults=[]
    xbad=[]
    xgood=[]

    plot_dir='Double_plot'
    os.makedirs(plot_dir,exist_ok=True)

    for one_file in filenames:
        try:
            xtab=ascii.read(one_file)
            foo=xtab['WAVE']
            foo=xtab['FLUX']
            foo=xtab['ERROR']
        except:
            xbad.append(one_file)
            continue
        try:
            xtab=ascii.read(one_file)
            results=do_one(xtab,wmin,wmax,xplot)
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


def do_all(filename='data/lvmSFrame-00009088.fits',wmin=6710,wmax=6750,outname='',do_plot=False):
    '''
    Do all of the spectra in a rss fits file
    
    '''
    try:
        x=fits.open(filename)
    except:
        print('Error: Could not open ',filename)
        return

    if do_plot:
        plot_dir='Double_plot'
        if outname!='':
            plot_dir='Double_plot_%s' % outname
        os.makedirs(plot_dir,exist_ok=True)


    wave=x['WAVE'].data
    flux=x['FLUX'].data
    b=np.sqrt(x['IVAR'].data)
    with np.errstate(divide='ignore', invalid='ignore'):
        error = np.where(b != 0, 1. / b, 1e10) # Replace with a large number
    mask=x['MASK'].data

        
    # Now get the good fibers from the science telecsope
    slittab=Table(x['SLITMAP'].data)
    good=scifib(slittab,'science')
    # print(len(good))
    # print(good)
    results=[]
    records=[]
    for i in range(len(good)):
        j=good['fiberid'][i]-1
        one_spec=Table([wave,flux[j],error[j],mask[j]],names=['WAVE','FLUX','ERROR','MASK'])
        one_spec=one_spec[one_spec['MASK']==0]
        if check_for_nan(flux[j])==False:
            rtab=do_one(one_spec,wmin=wmin,wmax=wmax,do_plot=do_plot,label='Fib %04d' % good['fiberid'][i])
            if len(rtab)>0:
                rtab['fiberid']=good['fiberid'][i]
                rtab['ra']=good['ra'][i]
                rtab['dec']=good['dec'][i]
                records.append(rtab)

                if do_plot:
                    plt.savefig('%s/p%05d_%05d_%05d.png' % (plot_dir,good['fiberid'][i],wmin,wmax))
            else:
                nan_count=np.isnan(flux[j]).sum()
                print('Nothing (%d/%d) fit for fiber %d at %.2f %.2f'  % (nan_count,len(flux[j]),good['fiberid'][i],good['ra'][i],good['dec'][i]))
        else:
            nan_count=np.isnan(flux[j]).sum()
            print('Too many nans (%d/%d) for  fiber %d at %.2f %.2f'  % (nan_count,len(flux[j]),good['fiberid'][i],good['ra'][i],good['dec'][i]))

    # Now clean the rows
    for table in records:
        table.meta.clear() 
    for table in records:
        colnames=table.colnames
        for col_name in colnames:
        
            if table[col_name].dtype.kind == 'O':  # Check if the column is of object type
                try:
                    table[col_name] = table[col_name].astype(float)  # Convert if possible
                except ValueError:
                    table[col_name] = np.nan  # Assign NaN if conversion fails

    results=vstack(records)

    if outname=='':
        outname=filename.split('/')[-1]
        outname=outname.replace('.fits','.double.txt')
    else:
        outname=outname+'.double.txt'

    outname=outname.replace('.txt','.%d_%d.txt' % (wmin,wmax))
    print('Output File:  %s ' % outname)

    columns=results.colnames
    for one in columns:
        if one.count('d1_'):
            results[one].format='.3f'
        elif one.count('d2_'):
            results[one].format='.3f'
        elif one.count('s1_'):
            results[one].format='.3f'
        elif one.count('aic_'):
            results[one].format='.3f'
        elif one.count('bic_'):
            results[one].format='.3f'
        elif one.count('redchi_'):
            results[one].format='.3f'
        elif one.count('f_stat'):
            results[one].format='.3e'
        elif one.count('p_value'):
            results[one].format='.3e'
        elif one.count('flux'):
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
    wmin=6710
    wmax=6730
    do_plot=False

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
        elif argv[i]=='-wmin':
            i+=1
            wmin=eval(argv[i])
        elif argv[i]=='-wmax':
            i+=1
            wmax=eval(argv[i])
        elif argv[i]=='-plot':
            do_plot=True
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
        results=do_all(one_file,wmin,wmax,outname,do_plot)
        # analyze(results)

    if len(specfiles)>0:
        do_individual(specfiles,wmin,wmax,outname,do_plot)


    return





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)

