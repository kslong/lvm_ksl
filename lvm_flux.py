#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Calculate fluxes from a line complex given the parameters 
from another line complex.


Command line usage (if any):

    usage: lvm_flux.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

250412 ksl Coding begun

'''

import sys
import os


from astropy.io import ascii,fits
from astropy.table import Table, join
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, Column, vstack
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

import matplotlib
from scipy import stats
from astropy.stats import sigma_clipped_stats


def gaussian(x, amp, x0, sigma):
    """
    Gaussian function: amp * exp(-(x - x0)^2 / (2 * sigma^2))
    
    Parameters:
    -----------
    x : array-like
        Wavelength values
    amp : float
        Amplitude of the Gaussian
    x0 : float
        Center of the Gaussian
    sigma : float
        Width of the Gaussian (sigma)
    
    Returns:
    --------
    array-like
        Gaussian values at points x
    """
    return amp * np.exp(-(x - x0)**2 / (2 * sigma**2))



def fit_gaussians(spectrum_table, lines_table):
    """
    Fit multiple Gaussian profiles to a spectrum with positive flux constraint.

    Parameters:
    -----------
    spectrum_table : astropy.table.Table
        Table containing the spectrum with columns 'WAVE', 'FLUX', 'ERROR'
    lines_table : astropy.table.Table
        Table containing line parameters with columns 'WAVE' (center), 'FWHM'

    Returns:
    --------
    tuple
        (updated_spectrum_table, updated_lines_table)
    """
    # Extract data from tables
    wave = spectrum_table['WAVE']
    flux = spectrum_table['FLUX']
    error = spectrum_table['ERROR']

    # Extract line parameters
    line_centers = lines_table['WAVE']
    line_fwhms = lines_table['FWHM']

    # Convert FWHM to sigma (FWHM = 2.355 * sigma)
    sigmas = line_fwhms / 2.355

    # Number of lines
    n_lines = len(line_centers)

    # Initial amplitude guesses - use max flux near each line center
    initial_amps = []
    for center in line_centers:
        # Find the wavelength point closest to line center
        idx = np.argmin(np.abs(wave - center))
        # Use this flux as initial amplitude guess, ensure positive
        initial_amps.append(max(0.1, flux[idx]))

    # Create the model function for multiple Gaussians
    def multi_gaussian(x, *params):
        """Model function: sum of multiple Gaussians"""
        y = np.zeros_like(x)
        for i in range(n_lines):
            # Each line has 1 parameter (amplitude)
            # Center and sigma are fixed
            amp = params[i]
            center = line_centers[i]
            sigma = sigmas[i]
            y += gaussian(x, amp, center, sigma)
        return y

    # Initial parameter guesses (just the amplitudes)
    initial_params = initial_amps

    # Set bounds for amplitudes (all must be positive)
    # Lower bound: 0 for all parameters
    # Upper bound: infinity for all parameters
    lower_bounds = [0.0] * n_lines
    upper_bounds = [np.inf] * n_lines
    bounds = (lower_bounds, upper_bounds)

    try:
        # Perform the fit with bounds to ensure positive amplitudes
        params, pcov = curve_fit(
            multi_gaussian,
            wave,
            flux,
            p0=initial_params,
            sigma=error,
            absolute_sigma=True,
            bounds=bounds,
            method='trf'  # Trust Region Reflective algorithm works well with bounds
        )

        # Calculate parameter errors from covariance matrix
        # Note: when using bounds, errors might be less reliable near the bounds
        param_errors = np.sqrt(np.diag(pcov))

        # Calculate individual line fluxes and update lines_table
        fitted_amps = params
        fitted_amp_errors = param_errors

        # Calculate the flux for each line (integral of Gaussian = amp * sigma * sqrt(2*pi))
        line_fluxes = []
        line_flux_errors = []

        for i in range(n_lines):
            # Calculate integrated flux for this line
            amp = fitted_amps[i]
            sigma = sigmas[i]
            flux_value = amp * sigma * np.sqrt(2 * np.pi)

            # Calculate error propagation
            amp_err = fitted_amp_errors[i]
            # Error propagation for amp * sigma * sqrt(2*pi)
            flux_error = flux_value * (amp_err / amp) if amp > 0 else 0.0

            line_fluxes.append(flux_value)
            line_flux_errors.append(flux_error)

        # Add flux and error columns to the lines_table
        lines_table['FLUX'] = line_fluxes
        lines_table['FLUX_ERROR'] = line_flux_errors

        # Calculate individual line contributions to the spectrum
        line_contributions = []
        for i in range(n_lines):
            amp = fitted_amps[i]
            center = line_centers[i]
            sigma = sigmas[i]
            contribution = gaussian(wave, amp, center, sigma)
            line_contributions.append(contribution)

            # Add individual line flux to spectrum_table
            spectrum_table[f'FLUX_LINE_{i+1}'] = contribution

        # Add total model flux to spectrum_table
        spectrum_table['FLUX_MODEL'] = multi_gaussian(wave, *params)

        # Compute residuals
        residuals=flux - spectrum_table['FLUX_MODEL']
        spectrum_table['RESIDUALS'] = residuals



        # Calculate reduced chi-squared statistic
        # Number of data points
        n_points = len(wave)
        # Number of free parameters: just n_lines amplitudes
        n_params = n_lines
        # Degrees of freedom

        dof = n_points - n_params

        if dof > 0:  # Ensure we don't divide by zero
            # Chi-squared: sum of ((observed - expected)/error)^2
            chi_squared = np.sum((residuals / error)**2)
            # Reduced chi-squared: chi-squared / degrees of freedom
            reduced_chi_squared = chi_squared / dof
        else:
            reduced_chi_squared = np.nan

        lines_table['REDUCED_CHI2']=reduced_chi_squared

        return spectrum_table, lines_table

    except Exception as e:
        print(f"Fitting failed: {e}")
        return spectrum_table, lines_table

def plot_results(spectrum_table, lines_table):
    """
    Plot the spectrum, fitted lines, and residuals

    Parameters:
    -----------
    spectrum_table : astropy.table.Table
        Updated spectrum table with model fluxes
    lines_table : astropy.table.Table
        Updated lines table with fitted fluxes
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]})

    # Plot the data and total model
    wave = spectrum_table['WAVE']
    flux = spectrum_table['FLUX']
    error = spectrum_table['ERROR']
    model = spectrum_table['FLUX_MODEL']
    residuals = spectrum_table['RESIDUALS']

    ax1.errorbar(wave, flux, yerr=error, fmt='o', alpha=0.5, label='Data')
    ax1.plot(wave, model, 'r-', lw=2, label='Total Model')

    # Plot individual lines
    n_lines = len(lines_table)
    for i in range(n_lines):
        line_flux = spectrum_table[f'FLUX_LINE_{i+1}']
        # Convert wavelength to a scalar value before string formatting
        wave_value = float(lines_table["WAVE"][i].item())
        ax1.plot(wave, line_flux, '--', lw=1,
                 label=f'Line {i+1} ({wave_value:.2f} Å)')

    # Plot residuals
    ax2.errorbar(wave, residuals, yerr=error, fmt='o', alpha=0.5)
    ax2.axhline(y=0, color='r', linestyle='-', lw=1)

    # Labels and legend
    ax1.set_xlabel('Wavelength')
    ax1.set_ylabel('Flux')
    ax2.set_xlabel('Wavelength')
    ax2.set_ylabel('Residuals')
    ax1.legend()

    plt.tight_layout()
    return

def fit_gaussians_with_shift(spectrum_table, lines_table, allow_shift=True, initial_shift=0.0):
    """
    Fit multiple Gaussian profiles to a spectrum with positive flux constraint
    and optional uniform wavelength shift.

    Parameters:
    -----------
    spectrum_table : astropy.table.Table
        Table containing the spectrum with columns 'WAVE', 'FLUX', 'ERROR'
    lines_table : astropy.table.Table
        Table containing line parameters with columns 'WAVE' (center), 'FWHM'
    allow_shift : bool, optional
        Whether to allow a uniform wavelength shift in the fit (default: True)
    initial_shift : float, optional
        Initial guess for the wavelength shift in the same units as WAVE (default: 0.0)

    Returns:
    --------
    tuple
        (updated_spectrum_table, updated_lines_table, wavelength_shift)
    """
    # Extract data from tables
    wave = spectrum_table['WAVE']
    flux = spectrum_table['FLUX']
    error = spectrum_table['ERROR']

    # Extract line parameters
    line_centers = lines_table['WAVE']
    line_fwhms = lines_table['FWHM']

    # Convert FWHM to sigma (FWHM = 2.355 * sigma)
    sigmas = line_fwhms / 2.355

    # Number of lines
    n_lines = len(line_centers)

    # Initial amplitude guesses - use max flux near each line center
    initial_amps = []
    for center in line_centers:
        # Find the wavelength point closest to line center
        idx = np.argmin(np.abs(wave - center))
        # Use this flux as initial amplitude guess, ensure positive
        initial_amps.append(max(0.1, flux[idx]))

    if allow_shift:
        # Create the model function for multiple Gaussians with a wavelength shift
        def multi_gaussian_with_shift(x, *params):
            """
            Model function: sum of multiple Gaussians with a wavelength shift
            First n_lines parameters are amplitudes
            Last parameter is the wavelength shift
            """
            y = np.zeros_like(x)
            # Extract the wavelength shift (last parameter)
            wave_shift = params[-1]

            for i in range(n_lines):
                # Each line has 1 parameter (amplitude)
                # Center is shifted by wave_shift
                # Sigma is fixed
                amp = params[i]
                center = line_centers[i] + wave_shift
                sigma = sigmas[i]
                y += gaussian(x, amp, center, sigma)
            return y

        # Initial parameter guesses (amplitudes + wavelength shift)
        initial_params = initial_amps + [initial_shift]

        # Set bounds for parameters
        # Lower bound: 0 for all amplitudes, reasonable negative shift
        # Upper bound: infinity for all amplitudes, reasonable positive shift
        lower_bounds = [0.0] * n_lines + [-20.0]  # Adjust shift bounds as needed
        upper_bounds = [np.inf] * n_lines + [20.0]  # Adjust shift bounds as needed
        bounds = (lower_bounds, upper_bounds)

        try:
            # Perform the fit with bounds
            params, pcov = curve_fit(
                multi_gaussian_with_shift,
                wave,
                flux,
                p0=initial_params,
                sigma=error,
                absolute_sigma=True,
                bounds=bounds,
                method='trf'  # Trust Region Reflective algorithm works well with bounds
            )

            # Extract fitted parameters
            fitted_amps = params[:-1]  # All but the last parameter
            wave_shift = params[-1]    # Last parameter is the wavelength shift

            # Calculate parameter errors from covariance matrix
            param_errors = np.sqrt(np.diag(pcov))
            fitted_amp_errors = param_errors[:-1]
            wave_shift_error = param_errors[-1]

            # Store the wavelength shift in the lines_table
            # Extract scalar values if the values are array-like
            scalar_shift = wave_shift
            scalar_shift_error = wave_shift_error

            lines_table['WAVE_SHIFT'] = [scalar_shift] * n_lines
            lines_table['WAVE_SHIFT_ERROR'] = [scalar_shift_error] * n_lines

            # Calculate the fitted wavelengths
            fitted_waves = []
            for i in range(n_lines):
                if hasattr(line_centers[i], 'item'):
                    fitted_waves.append(line_centers[i].item() + scalar_shift)
                else:
                    fitted_waves.append(line_centers[i] + scalar_shift)

            # Store the shifted line centers
            lines_table['WAVE_FITTED'] = fitted_waves

            # Calculate the flux for each line
            line_fluxes = []
            line_flux_errors = []

            for i in range(n_lines):
                # Calculate integrated flux for this line
                amp = fitted_amps[i]
                sigma = sigmas[i]
                flux_value = amp * sigma * np.sqrt(2 * np.pi)

                # Calculate error propagation
                amp_err = fitted_amp_errors[i]
                # Error propagation for amp * sigma * sqrt(2*pi)
                flux_error = flux_value * (amp_err / amp) if amp > 0 else 0.0

                line_fluxes.append(flux_value)
                line_flux_errors.append(flux_error)

            # Add flux and error columns to the lines_table
            lines_table['FLUX'] = line_fluxes
            lines_table['FLUX_ERROR'] = line_flux_errors

            # Calculate individual line contributions to the spectrum
            for i in range(n_lines):
                amp = fitted_amps[i]
                center = line_centers[i] + wave_shift  # Apply the shift
                sigma = sigmas[i]
                contribution = gaussian(wave, amp, center, sigma)

                # Add individual line flux to spectrum_table
                spectrum_table[f'FLUX_LINE_{i+1}'] = contribution

            # Add total model flux to spectrum_table
            spectrum_table['FLUX_MODEL'] = multi_gaussian_with_shift(wave, *params)

            # Compute residuals
            residuals=flux - spectrum_table['FLUX_MODEL']
            spectrum_table['RESIDUALS'] = residuals



            # Calculate reduced chi-squared statistic
            # Number of data points
            n_points = len(wave)
            # Number of free parameters: just n_lines amplitudes
            n_params = n_lines
            # Degrees of freedom

            dof = n_points - n_params

            if dof > 0:  # Ensure we don't divide by zero
                # Chi-squared: sum of ((observed - expected)/error)^2
                chi_squared = np.sum((residuals / error)**2)
                # Reduced chi-squared: chi-squared / degrees of freedom
                reduced_chi_squared = chi_squared / dof
            else:
                 reduced_chi_squared = np.nan


            lines_table['REDUCED_CHI2']=reduced_chi_squared

            return spectrum_table, lines_table, wave_shift

        except Exception as e:
            print(f"Fitting failed: {e}")
            return spectrum_table, lines_table, 0.0

    else:
        # Use the original fit_gaussians function without shift
        spectrum_table, lines_table = fit_gaussians(spectrum_table, lines_table)
        return spectrum_table, lines_table, 0.0

def plot_results_with_shift(spectrum_table, lines_table):
    """
    Plot the spectrum, fitted lines, and residuals with shift information
    
    Parameters:
    -----------
    spectrum_table : astropy.table.Table
        Updated spectrum table with model fluxes
    lines_table : astropy.table.Table
        Updated lines table with fitted fluxes and wavelength shift
    
    Returns:
    --------
    matplotlib.figure.Figure
        The figure with the plots
    """
    # Clear any existing plots to avoid duplication
    plt.close('all')
    
    # Create a new figure with two panels
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]})
    
    # Extract data
    wave = spectrum_table['WAVE']
    flux = spectrum_table['FLUX']
    error = spectrum_table['ERROR']
    model = spectrum_table['FLUX_MODEL']
    residuals = spectrum_table['RESIDUALS']
    
    # TOP PANEL: Data and fit
    ax1.errorbar(wave, flux, yerr=error, fmt='o', alpha=0.5, label='Data')
    ax1.plot(wave, model, 'r-', lw=2, label='Total Model')
    
    # Plot individual lines
    n_lines = len(lines_table)
    for i in range(n_lines):
        if f'FLUX_LINE_{i+1}' in spectrum_table.colnames:
            line_flux = spectrum_table[f'FLUX_LINE_{i+1}']
            
            # Get the correct wavelength (shifted if available)
            if 'WAVE_FITTED' in lines_table.colnames:
                wave_value = lines_table["WAVE_FITTED"][i]
                # Convert to float if it's an array or numpy type
                if hasattr(wave_value, 'item'):
                    wave_value = float(wave_value.item())
                else:
                    wave_value = float(wave_value)
            else:
                wave_value = lines_table["WAVE"][i]
                if hasattr(wave_value, 'item'):
                    wave_value = float(wave_value.item())
                else:
                    wave_value = float(wave_value)
                
            ax1.plot(wave, line_flux, '--', lw=1, 
                     label=f'Line {i+1} ({wave_value:.2f} Å)')
    
    # BOTTOM PANEL: Residuals
    ax2.errorbar(wave, residuals, yerr=error, fmt='o', alpha=0.5)
    ax2.axhline(y=0, color='r', linestyle='-', lw=1)
    
    # Add shift information to the title if available
    if 'WAVE_SHIFT' in lines_table.colnames and len(lines_table['WAVE_SHIFT']) > 0:
        shift_value = lines_table['WAVE_SHIFT'][0]
        if hasattr(shift_value, 'item'):
            shift_value = float(shift_value.item())
        else:
            shift_value = float(shift_value)
            
        if 'WAVE_SHIFT_ERROR' in lines_table.colnames:
            shift_error = lines_table['WAVE_SHIFT_ERROR'][0]
            if hasattr(shift_error, 'item'):
                shift_error = float(shift_error.item())
            else:
                shift_error = float(shift_error)
        else:
            shift_error = 0.0
            
        fig.suptitle(f'Fitted Spectrum (Wavelength Shift: {shift_value:.3f} ± {shift_error:.3f} Å)')
    
    # Labels and legend
    ax1.set_xlabel('Wavelength (Å)')
    ax1.set_ylabel('Flux')
    ax2.set_xlabel('Wavelength (Å)')
    ax2.set_ylabel('Residuals')
    ax1.legend()
    
    plt.tight_layout()
    if 'WAVE_SHIFT' in lines_table.colnames:
        plt.subplots_adjust(top=0.93)  # Make room for suptitle
    
    return

def fits2table(fib=300,wmin=6721,wmax=6741,filename='W28_all_ave.fits'):
    # print(fib,wmin,wmax,filename)

    e_factor=2.25

    x=fits.open(filename)
    wave=x['WAVE'].data
    flux=x['FLUX'].data[fib-1]*1e16
    error=x['IVAR'].data[fib-1]

    mask=x['MASK'].data[fib-1]
    spectrum=Table([wave,flux,error,mask],names=['WAVE','FLUX','ERROR','MASK'])
    spectrum=spectrum[spectrum['WAVE']>wmin]
    spectrum=spectrum[spectrum['WAVE']<wmax]
    spectrum=spectrum[spectrum['MASK']==0]
    spectrum['ERROR']=1/np.sqrt(spectrum['ERROR'])*1e16/e_factor
    mean_flux, median_flux, std_flux = sigma_clipped_stats(spectrum['FLUX'], sigma_upper=2,sigma_lower=3)
    spectrum['FLUX']-=median_flux
    return spectrum


def construct_inputs(fib=500, w_in=6730.8, w_out=6716.4, fit_tab='lvmSFrame-00003470.triple.6707_6727.txt'):
    '''
    Read a table which caondatins the outputs of lvm_triple and use it to fit the line complex of another line

    '''
    try:
        xresult = ascii.read(fit_tab)
    except:
        print('Error: Could not read %s' % fib_tab)
        return []

    xresult = xresult[xresult['fiberid'] == fib]


    # If no rows match, return empty table
    if len(xresult) == 0:
        print('No match for fib %5d' % fib)
        return []

    # Create simple lists to store scalar values
    w = []
    fwhm = []
    # print(xresult['preferred_model'])
    if xresult['preferred_model']=='single_gaussian':
        w.append(xresult['s1_wave'][0])
        fwhm.append(xresult['s1_sigma'][0])
    elif xresult['preferred_model']=='double_gaussian':
        w.append(xresult['d1_wave'][0])
        w.append(xresult['d2_wave'][0])
        fwhm.append(xresult['d1_sigma'][0])
        fwhm.append(xresult['d2_sigma'][0])
    elif xresult['preferred_model']=='triple_gaussian':
        w.append(xresult['t1_wave'][0])
        w.append(xresult['t2_wave'][0])
        w.append(xresult['t3_wave'][0])
        fwhm.append(xresult['t1_sigma'][0])
        fwhm.append(xresult['t2_sigma'][0])
        fwhm.append(xresult['t3_sigma'][0])
    else:
        print('Preferreed model type is unknown :', xresult['preferred_model'])



    # Create table with scalar values
    xtab = Table()
    xtab['WAVE'] = w
    xtab['FWHM'] = fwhm

    # Apply FWHM adjustment and wavelength shift
    xtab['FWHM'] *= 2.35
    xtab['WAVE'] -= (w_in - w_out)

    xtab['FWHM'].format='.3f'
    xtab['WAVE'].format='.3f'

    xtab['preferred_model']=xresult['preferred_model']
    xtab['fiberid']=xresult['fiberid']

    # xtab.info()

    return xtab

def do_one(fib, w_in, w_out, wmin, wmax, fit_tab, rss_file='W28_all_ave.fits', do_plot=True):
    '''
    Process a single spectrum and ensure all table values are scalar
    '''
    spec = fits2table(fib, wmin, wmax, filename=rss_file)
    if len(spec) == 0:
        print('Spectrum had all zeros for %f' % fib)
        return []
    
    linetab = construct_inputs(fib, w_in, w_out, fit_tab)
    if len(linetab) == 0:
        # print('linetab empty for %4d' % fib)
        return []
    
    try:
        spectrum_table, lines_table, shift = fit_gaussians_with_shift(
            spectrum_table=spec, 
            lines_table=linetab, 
            allow_shift=True, 
            initial_shift=0.0
        )
        
        # Create a new table with guaranteed scalar values
        clean_lines_table = Table()
        
        # Copy each column, ensuring scalar values
        for col_name in lines_table.colnames:
            if col_name in ['fiberid']:
                # For columns that should be the same for all rows
                clean_lines_table[col_name] = fib
                continue
                
            # Extract values as scalars
            clean_values = []
            for i in range(len(lines_table)):
                val = lines_table[col_name][i]
                # Convert to scalar if needed
                if hasattr(val, 'item'):
                    val = val.item()
                clean_values.append(val)
            
            # Add column with scalar values
            clean_lines_table[col_name] = clean_values
        
        plot_dir = 'Flux_plot'
        if os.path.isdir(plot_dir) == False:
            os.makedirs(plot_dir, exist_ok=True)
        
        outname = 'P%04d_%4d_%04d' % (fib, w_in, w_out)
        if do_plot:
            plot_results_with_shift(spectrum_table, clean_lines_table)
            plt.savefig('%s/%s.png' % (plot_dir, outname))
        
        clean_lines_table['WAVE'].format='.3f'
        clean_lines_table['FWHM'].format='.3f'
        clean_lines_table['WAVE_SHIFT'].format='.3f'
        clean_lines_table['WAVE_SHIFT_ERROR'].format='.3f'
        clean_lines_table['WAVE_FITTED'].format='.3f'
        clean_lines_table['FLUX'].format='.3f'
        clean_lines_table['FLUX_ERROR'].format='.3f'
        clean_lines_table['REDUCED_CHI2'].format='.3f'
        return clean_lines_table
    
    except Exception as e:
        print(f'failure to fit: %04d - {str(e)}' % fib)
        return []

def do_all(w_in,w_out,wmin,wmax,fit_tab,rss_file,do_plot=True):
    x=fits.open(rss_file)
    fib_tab=Table(x['SLITMAP'].data)
    fib_tab=fib_tab[fib_tab['fibstatus']==0]
    good_fibers=np.array(fib_tab['fiberid'])
    x.close()
    results=[]
    icount=0
    failed = []
    for one_fib in good_fibers:
        one_result=do_one(one_fib,w_in,w_out,wmin,wmax,fit_tab,rss_file,do_plot)
        if len(one_result)>0:
            results.append(one_result)
        else:
            failed.append(one_fib)
        icount+=1
        if icount % 50 == 0:
            print('Finished %d of %d fits' % (icount,len(good_fibers)))
    xresults=vstack(results)
    outname=rss_file.replace('.fits','')
    outname='%s.flux.%4d_%4d.txt' % (outname,w_in,w_out)
    print('writing to ',outname)
    print(xresults.info())
    xresults.write(outname,format='ascii.fixed_width_two_line',overwrite=True)

    print('Of %d attempts, %d produced no results' % (len(good_fibers),len(failed)))
    for one in failed:
        print(one)

    return xresults


    
def steer(argv):

    w_in=-1.    
    w_out=-1.    
    wmin=-1  
    wmax=-1
    rssfile=''
    specfile=''
    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-w_in':
            i+=1
            w_in=eval(argv[i])
        elif argv[i]=='-w_out':
            i+=1
            w_out=eval(argv[i])
        elif argv[i]=='-wmin':
            i+=1
            wmin=eval(argv[i])
        elif argv[i]=='-wmax':
            i+=1
            wmax=eval(argv[i])
        elif argv[i][0]=='1':
            print('Error: Unknown switch: ', argv)
            return
        elif argv[i].count('.txt'):
            specfile=argv[i]
        elif argv[i].count('.fits'):
            rssfile=argv[i]
        elif w_in<0:
            w_in=eval(argv[i])
        elif w_out<9:
            w_out=eval(argv[i])
        i+=1

    if specfile=='' or rssfile=='':
        print('Not enough inputs - specfile or rss file not specified: ',argv)
        return

    if w_out==-1 or w_out==-1:
        print('Not enough inputs - wavelengths w_in and w_out not  specified: ',argv)
        return

    if wmin<0:
        wmin=w_out-10.
    if wmax<0:
        wmax=w_out+10


    do_all(w_in,w_out,wmin,wmax,fit_tab=specfile,rss_file=rssfile,do_plot=True)

    print('Warning:  Errors have been rescaled by 2.25, which should change with a new reduction')
    return
    
# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
