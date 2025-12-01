#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This is a routine to find peaks in spectra (lvm spectra).


Command line usage (if any):

    usage: peaks.py [-h] -wmin 6720 -wmax 6740  [-all]   filename.fits   one or more fiber numbers

    where 

    -h prints this documentation and exits
    -wmin   the minimum wavelength use for searching for peaks
    -wmax   the maxium wavelength used for seach 
    -all    an option to process all fibers, in which case the fibers are all used



Description:  

    The routine looks for peaks in a wavelength interval and creates an output table, and 
    figures to show what peaks were found.

Primary routines:

    do_one

Notes:

    Currently there are a lot of hardwired portions of the routine, including the
    criteria for finding peaks.  It thise routine is used extensively it likely 
    needs tuning.
                                       
History:

250408 ksl Coding begun

'''

import os

from astropy.io import ascii,fits
from astropy.table import Table, join, vstack
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
from astropy.stats import sigma_clipped_stats
from scipy.signal import find_peaks, peak_widths
from scipy import stats
from scipy.optimize import curve_fit
from datetime import date




def gaussian(x, amp, cen, sigma, offset):
    """Gaussian function with offset"""
    return amp * np.exp(-(x - cen)**2 / (2 * sigma**2)) + offset

def multi_gaussian(x, *params):
    """
    Multiple Gaussian model for fitting closely spaced lines.
    
    Parameters:
    x: wavelength array
    params: list of parameters [amp1, cen1, sigma1, amp2, cen2, sigma2, ..., offset]
    """
    n_gaussians = (len(params) - 1) // 3
    offset = params[-1]
    result = np.zeros_like(x) + offset
    
    for i in range(n_gaussians):
        amp = params[i*3]
        cen = params[i*3 + 1]
        sigma = params[i*3 + 2]
        result += amp * np.exp(-(x - cen)**2 / (2 * sigma**2))
        
    return result




def plot_emission_lines(wave, flux, continuum, emission_flux, peaks_table, valid_peaks, fitted_regions=None,name=''):
    """
    Plot the detected emission lines in the spectrum.
    
    Parameters
    ----------
    wave : array-like
        Wavelength array
    flux : array-like
        Original flux array
    continuum : float
        Estimated continuum level
    emission_flux : array-like
        Continuum-subtracted flux array
    peaks_table : astropy.table.Table
        Table of detected emission peaks
    valid_peaks : array-like
        Indices of valid peaks in the original arrays
    fitted_regions : list of tuples, optional
        List of (region_slice, fitted_flux, param_values) for fitted regions
        
    Returns
    -------
    None
    """
    plt.figure(figsize=(12, 8))
    
    # Extract data from peaks table
    peak_waves = peaks_table['center_wavelength']
    peak_fluxes = peaks_table['flux']
    valid_snr = peaks_table['snr']
    fwhm_wavelengths = peaks_table['fwhm']
    
    # Plot original spectrum
    plt.subplot(2, 1, 1)
    plt.plot(wave, flux, 'k-', alpha=0.7, label='Spectrum')
    plt.plot(wave, np.ones_like(wave) * continuum, 'r--', label='Continuum')
    
    # Highlight detected peaks
    plt.plot(peak_waves, peak_fluxes, 'ro', label='Detected Lines')
    
    # Plot any fitted regions
    print(fitted_regions)
    if fitted_regions:
        for region_slice, fitted_flux, params, n_peaks in fitted_regions:
            region_wave = wave[region_slice]
            plt.plot(region_wave, fitted_flux + continuum, 'g-', linewidth=2, alpha=0.7, 
                    label=f'{n_peaks}-Component Fit' if n_peaks > 1 else 'Gaussian Fit')
    
    # Add peak labels
    for i, (w, f, s) in enumerate(zip(peak_waves, peak_fluxes, valid_snr)):
        plt.text(w, f*1.05, f"{w:.2f}\nSNR={s:.1f}", 
                 ha='center', va='bottom', fontsize=8)
    
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    # if name=='':
    #    plt.title('Spectrum with Detected Emission Lines')
    # else:
    #     plt.title(name)

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.grid(alpha=0.3)
    
    # Plot continuum-subtracted spectrum
    plt.subplot(2, 1, 2)
    plt.plot(wave, emission_flux, 'b-', alpha=0.7, label='Continuum Subtracted')
    plt.plot(peak_waves, emission_flux[valid_peaks], 'ro', label='Detected Lines')
    
    # Plot fitted regions in subtracted spectrum
    if fitted_regions:
        for region_slice, fitted_flux, params, n_peaks in fitted_regions:
            region_wave = wave[region_slice]
            plt.plot(region_wave, fitted_flux, 'g-', linewidth=2, alpha=0.7, 
                    label=f'{n_peaks}-Component Fit' if n_peaks > 1 else 'Gaussian Fit')
    
    # Plot zero line
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    # Show FWHM for each line
    wave_spacing = np.median(np.diff(wave))
    for i, (peak_idx, fwhm) in enumerate(zip(valid_peaks, fwhm_wavelengths)):
        half_max = emission_flux[peak_idx] / 2
        fwhm_points = fwhm / wave_spacing
        left_idx = int(peak_idx - fwhm_points / 2)
        right_idx = int(peak_idx + fwhm_points / 2)
        if left_idx < 0:
            left_idx = 0
        if right_idx >= len(wave):
            right_idx = len(wave) - 1
            
        plt.plot([wave[left_idx], wave[right_idx]], 
                 [half_max, half_max], 'g-', alpha=0.6)
        
    plt.xlabel('Wavelength')
    plt.ylabel('Emission Flux')
    if name=='':
        plt.title('Continuum-Subtracted Spectrum')
    else:
        plt.title(name)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.grid(alpha=0.3)
    
    plt.tight_layout()



def find_emission_lines(spec_table, min_snr=3.0, typical_fwhm=None, 
                        plot_results=True, continuum_percentile=50,
                        resolve_close_lines=True, prominence=0.1,name=''):
    """
    Find emission lines in a spectrum given as an Astropy Table,
    with dramatically improved detection of closely spaced lines.
    
    Parameters
    ----------
    spec_table : astropy.table.Table
        Table containing columns 'WAVE', 'FLUX', and 'ERROR'
    min_snr : float, optional
        Minimum signal-to-noise ratio for peak detection. Default is 3.0.
    typical_fwhm : float, optional
        Typical FWHM of emission lines in same units as wavelength.
        If None, will estimate from data.
    plot_results : bool, optional
        Whether to plot the results. Default is True.
    continuum_percentile : float, optional
        Percentile used to estimate the continuum. Default is 50.
    resolve_close_lines : bool, optional
        Whether to attempt to resolve closely spaced lines. Default is True.
    prominence : float, optional
        Minimum prominence for peak detection, as a fraction of the max emission.
        
    Returns
    -------
    peaks_table : astropy.table.Table
        Table containing detected peaks ordered by significance
    """
    # Extract data from the input table
    wave = spec_table['WAVE']
    flux = spec_table['FLUX']
    error = spec_table['ERROR']
    
    # Estimate the continuum level using a percentile
    continuum = np.percentile(flux, continuum_percentile)
    
    # Subtract continuum to get emission-only spectrum
    emission_flux = flux - continuum
    
    # Calculate the typical width in data points if not provided
    wave_spacing = np.median(np.diff(wave))
    if typical_fwhm is None:
        # Assume typical FWHM is ~10 data points if not specified
        typical_fwhm_points = 10
    else:
        typical_fwhm_points = int(typical_fwhm / wave_spacing)
    
    # Calculate sigma (std dev) from FWHM
    typical_sigma_points = typical_fwhm_points / 2.355
    
    # DRAMATICALLY IMPROVED PEAK DETECTION:
    # ====================================
    
    # 1. First, try to detect all potential peaks with very sensitive settings
    # Use much more relaxed parameters to catch all potential peaks
    initial_width = max(1, int(typical_fwhm_points / 8))  # Much smaller width requirement
    initial_distance = max(1, int(initial_width / 2))     # Minimum distance between peaks
    initial_prominence = prominence * np.max(emission_flux) / 3  # Lower prominence threshold

    
    # Use very sensitive settings for initial peak detection
    potential_peaks, _ = find_peaks(emission_flux, 
                                   height=0,  # Just require positive peaks
                                   width=initial_width,
                                   prominence=initial_prominence,
                                   distance=initial_distance)
    
    print(f"Initial sensitive peak detection found {len(potential_peaks)} potential peaks")
    
    # 2. Refine these potential peaks with more restrictive criteria
    if len(potential_peaks) == 0:
        print("No potential emission lines detected with even the most sensitive parameters.")
        return Table()
    
    # Calculate SNR for each potential peak
    potential_snr = emission_flux[potential_peaks] / error[potential_peaks]
    valid_peaks = potential_peaks[potential_snr >= min_snr]
    valid_snr = potential_snr[potential_snr >= min_snr]
    
    print(f"After SNR filtering, {len(valid_peaks)} peaks remain")
    
    if len(valid_peaks) == 0:
        print(f"No emission lines with SNR >= {min_snr} detected.")
        return Table()
    
    # 3. If initial detection didn't find enough peaks, try derivative-based detection
    if len(valid_peaks) <= 1 and resolve_close_lines:
        print("Only one peak detected. Trying derivative-based detection to find close lines...")
        
        # Calculate first and second derivatives
        first_deriv = np.gradient(emission_flux)
        second_deriv = np.gradient(first_deriv)
        
        # Find zero-crossings in first derivative (potential peak locations)
        # A peak has first_deriv = 0 and second_deriv < 0
        zero_crossings = []
        for i in range(1, len(first_deriv)):
            if (first_deriv[i-1] > 0 and first_deriv[i] <= 0) and second_deriv[i] < 0:
                zero_crossings.append(i)
        
        print(f"Derivative analysis found {len(zero_crossings)} potential peak locations")
        
        # Add these to our valid peaks if they exceed SNR threshold
        for idx in zero_crossings:
            if idx not in valid_peaks:  # Avoid duplicates
                peak_snr = emission_flux[idx] / error[idx]
                if peak_snr >= min_snr:
                    valid_peaks = np.append(valid_peaks, idx)
                    valid_snr = np.append(valid_snr, peak_snr)
        
        # Sort valid_peaks by position
        sorted_indices = np.argsort(valid_peaks)
        valid_peaks = valid_peaks[sorted_indices]
        valid_snr = valid_snr[sorted_indices]
        
        print(f"After derivative analysis, {len(valid_peaks)} peaks found")
        
    # Calculate peak properties
    peak_waves = wave[valid_peaks]
    peak_fluxes = flux[valid_peaks]
    peak_emissions = emission_flux[valid_peaks]
    
    # 4. Try multi-component deconvolution for the entire spectrum using expected peak positions
    fwhm_wavelengths = []
    all_fitted_regions = []
    
    # Process peaks either individually or in groups
    if resolve_close_lines:
        # Group nearby peaks with a more careful approach
        grouped_peaks = []
        min_separation = typical_fwhm_points * 0.7  # Use 0.7x FWHM as threshold
        
        if len(valid_peaks) > 0:
            current_group = [valid_peaks[0]]
            
            for i in range(1, len(valid_peaks)):
                if valid_peaks[i] - current_group[-1] < min_separation:
                    current_group.append(valid_peaks[i])
                else:
                    if current_group:
                        grouped_peaks.append(current_group)
                    current_group = [valid_peaks[i]]
            
            if current_group:
                grouped_peaks.append(current_group)
                
        print(f"Grouped peaks into {len(grouped_peaks)} groups")
        for i, group in enumerate(grouped_peaks):
            print(f"  Group {i+1}: {len(group)} peaks at wavelengths:", 
                  [f"{wave[idx]:.2f}" for idx in group])
        
        # Process each group or individual peak
        new_peaks = []
        new_peak_waves = []
        new_peak_fluxes = []
        new_peak_emissions = []
        new_snr = []
        
        for group in grouped_peaks:
            if len(group) == 1:
                # Single isolated peak - use standard width measurement
                peak_idx = group[0]
                
                # Measure FWHM directly from data
                half_max = emission_flux[peak_idx] / 2
                left_idx = peak_idx
                while left_idx > 0 and emission_flux[left_idx] > half_max:
                    left_idx -= 1
                    
                right_idx = peak_idx
                while right_idx < len(emission_flux)-1 and emission_flux[right_idx] > half_max:
                    right_idx += 1
                
                fwhm_pts = right_idx - left_idx
                fwhm_wavelengths.append(fwhm_pts * wave_spacing)
                
                # Store the peak
                new_peaks.append(peak_idx)
                new_peak_waves.append(wave[peak_idx])
                new_peak_fluxes.append(flux[peak_idx])
                new_peak_emissions.append(emission_flux[peak_idx])
                new_snr.append(emission_flux[peak_idx] / error[peak_idx])
            
            else:
                # Multiple close peaks - fit with multiple Gaussians
                # Define a wider region around the group for better fitting
                min_idx = max(0, min(group) - int(typical_fwhm_points * 2))
                max_idx = min(len(wave) - 1, max(group) + int(typical_fwhm_points * 2))
                region_slice = slice(min_idx, max_idx + 1)
                region_wave = wave[region_slice]
                region_flux = emission_flux[region_slice]
                region_error = error[region_slice]
                
                # Ensure errors are positive
                if np.any(region_error <= 0):
                    region_error = np.ones_like(region_error) * np.mean(region_error[region_error > 0])
                    if np.mean(region_error) <= 0:
                        region_error = np.ones_like(region_error) * np.mean(region_flux) * 0.01
                
                # Setup initial parameters and bounds
                params = []
                lower_bounds = []
                upper_bounds = []
                
                for idx in group:
                    # Get precise peak parameters
                    amp = emission_flux[idx]
                    cen = wave[idx]
                    sigma = typical_sigma_points * wave_spacing * 0.8  # Start with narrower lines
                    
                    params.extend([amp, cen, sigma])
                    
                    # Set tight constraints on line position
                    lower_bounds.append(0.2 * amp)  # Min amplitude
                    upper_bounds.append(5.0 * amp)  # Max amplitude
                    
                    # Allow limited movement around detected positions
                    lower_bounds.append(cen - 0.3 * typical_fwhm_points * wave_spacing)
                    upper_bounds.append(cen + 0.3 * typical_fwhm_points * wave_spacing)
                    
                    # Critical: TIGHT constraints on line width
                    # This prevents fitting a single broad component instead of multiple narrow ones
                    lower_bounds.append(0.3 * sigma)  # Minimum width (narrower than typical)
                    upper_bounds.append(1.2 * sigma)  # Maximum width (only slightly broader than typical)
                
                # Add offset parameter
                params.append(0)
                lower_bounds.append(-0.1 * np.max(region_flux))
                upper_bounds.append(0.1 * np.max(region_flux))
                
                try:
                    # Fit multiple Gaussians with strict constraints
                    popt, pcov = curve_fit(multi_gaussian, region_wave, region_flux, 
                                          p0=params, sigma=region_error,
                                          maxfev=10000, bounds=(lower_bounds, upper_bounds))
                    
                    # Calculate fitted flux
                    fitted_flux = multi_gaussian(region_wave, *popt)
                    all_fitted_regions.append((region_slice, fitted_flux, popt))
                    
                    # Evaluate fit quality
                    residuals = region_flux - fitted_flux
                    rel_residuals = np.abs(residuals) / (np.abs(region_flux) + np.median(np.abs(region_flux))/10)
                    fit_quality = 1.0 - np.mean(rel_residuals)
                    
                    print(f"Group at {wave[group[0]]:.2f}: fit quality = {fit_quality:.2f}")
                    
                    # Process the fitted parameters for each component
                    n_gaussians = (len(popt) - 1) // 3
                    for i in range(n_gaussians):
                        amp = popt[i*3]
                        cen = popt[i*3 + 1]
                        sigma = popt[i*3 + 2]
                        
                        print(f"  Component {i+1}: center={cen:.2f}, amp={amp:.2f}, sigma={sigma:.2f}")
                        
                        # Calculate FWHM from sigma
                        fwhm = 2.355 * sigma
                        fwhm_wavelengths.append(fwhm)
                        
                        # Find the nearest wavelength point
                        nearest_idx = np.argmin(np.abs(wave - cen))
                        
                        # Calculate SNR at fitted peak
                        fitted_peak_emission = amp  # Pure amplitude, offset already included
                        fit_error = np.interp(cen, wave[region_slice], region_error)
                        fitted_snr = fitted_peak_emission / fit_error
                        
                        if fitted_snr >= min_snr:
                            new_peaks.append(nearest_idx)
                            new_peak_waves.append(cen)  # Use fitted center
                            new_peak_fluxes.append(amp + continuum)  # Add continuum back
                            new_peak_emissions.append(fitted_peak_emission)
                            new_snr.append(fitted_snr)
                
                except (RuntimeError, ValueError) as e:
                    print(f"Fitting failed for group at {wave[group[0]]:.2f}: {str(e)}")
                    print("Falling back to individual peak processing")
                    
                    # Process each peak individually
                    for peak_idx in group:
                        # Try to fit a single Gaussian to this peak
                        try:
                            # Define a small region around just this peak
                            mini_min_idx = max(0, peak_idx - int(typical_fwhm_points))
                            mini_max_idx = min(len(wave) - 1, peak_idx + int(typical_fwhm_points))
                            mini_slice = slice(mini_min_idx, mini_max_idx + 1)
                            mini_wave = wave[mini_slice]
                            mini_flux = emission_flux[mini_slice]
                            mini_error = error[mini_slice]
                            
                            # Ensure errors are positive
                            if np.any(mini_error <= 0):
                                mini_error = np.ones_like(mini_error) * np.mean(mini_error[mini_error > 0])
                                if np.mean(mini_error) <= 0:
                                    mini_error = np.ones_like(mini_error) * np.mean(mini_flux) * 0.01
                            
                            # Initial parameters for single Gaussian
                            p0 = [emission_flux[peak_idx], wave[peak_idx], 
                                  typical_sigma_points * wave_spacing, 0]
                            
                            # Bounds for single Gaussian
                            l_bounds = [0.2 * emission_flux[peak_idx],  # amp min
                                       wave[peak_idx] - 0.3 * typical_fwhm_points * wave_spacing,  # center min
                                       0.3 * typical_sigma_points * wave_spacing,  # sigma min
                                       -0.1 * emission_flux[peak_idx]]  # offset min
                            
                            u_bounds = [5.0 * emission_flux[peak_idx],  # amp max
                                       wave[peak_idx] + 0.3 * typical_fwhm_points * wave_spacing,  # center max
                                       1.2 * typical_sigma_points * wave_spacing,  # sigma max
                                       0.1 * emission_flux[peak_idx]]  # offset max
                            
                            # Fit single Gaussian
                            mini_popt, _ = curve_fit(gaussian, mini_wave, mini_flux,
                                                    p0=p0, sigma=mini_error,
                                                    bounds=(l_bounds, u_bounds))
                            
                            # Extract fitted parameters
                            amp, cen, sigma, offset = mini_popt
                            fwhm = 2.355 * sigma
                            fwhm_wavelengths.append(fwhm)
                            
                            # Find nearest index
                            nearest_idx = np.argmin(np.abs(wave - cen))
                            
                            # Calculate SNR
                            fitted_snr = amp / np.mean(mini_error)
                            
                            if fitted_snr >= min_snr:
                                new_peaks.append(nearest_idx)
                                new_peak_waves.append(cen)
                                new_peak_fluxes.append(amp + continuum)
                                new_peak_emissions.append(amp)
                                new_snr.append(fitted_snr)
                                
                        except (RuntimeError, ValueError) as e2:
                            print(f"  Single peak fitting also failed for peak at {wave[peak_idx]:.2f}: {str(e2)}")
                            # Last resort: just use the raw peak
                            fwhm_pts = typical_fwhm_points  # Use typical value
                            fwhm_wavelengths.append(fwhm_pts * wave_spacing)
                            
                            new_peaks.append(peak_idx)
                            new_peak_waves.append(wave[peak_idx])
                            new_peak_fluxes.append(flux[peak_idx])
                            new_peak_emissions.append(emission_flux[peak_idx])
                            new_snr.append(emission_flux[peak_idx] / error[peak_idx])
        
        # Update peaks with our refined measurements
        if len(new_peaks) > 0:
            valid_peaks = np.array(new_peaks)
            peak_waves = np.array(new_peak_waves)
            peak_fluxes = np.array(new_peak_fluxes)
            peak_emissions = np.array(new_peak_emissions)
            valid_snr = np.array(new_snr)
            fwhm_wavelengths = np.array(fwhm_wavelengths)
        else:
            print("Warning: No valid peaks found after processing")
            return Table()
    
    else:
        # Just use standard FWHM calculation for all peaks without multi-component fitting
        widths_result = peak_widths(emission_flux, valid_peaks, rel_height=0.5)
        fwhm_points = widths_result[0]
        fwhm_wavelengths = fwhm_points * wave_spacing
    
    # Calculate significance as detection SNR * peak flux / continuum ratio
    peak_significance = valid_snr * (peak_fluxes / continuum)
    
    # Create a result table
    peaks_table = Table()
    peaks_table['center_wavelength'] = peak_waves
    peaks_table['flux'] = peak_fluxes
    peaks_table['emission_flux'] = peak_emissions
    peaks_table['snr'] = valid_snr
    peaks_table['fwhm'] = fwhm_wavelengths
    peaks_table['significance'] = peak_significance

    
    # Sort by significance
    peaks_table.sort('significance', reverse=True)
    
    # Plot results if requested
    if plot_results:
        plot_emission_lines(wave, flux, continuum, emission_flux, peaks_table, valid_peaks, all_fitted_regions,name=name)
    
    return peaks_table


def fits2table(fib,wmin,wmax,filename,plot=False):
    try:
        x=fits.open(filename)
    except:
        print('Error: Could not open %s' % filename)
        raise IOError
    wave=x['WAVE'].data
    flux=x['FLUX'].data[fib-1]*1e16
    error=x['IVAR'].data[fib-1]


    mask=x['MASK'].data[fib-1]
    spectrum=Table([wave,flux,error,mask],names=['WAVE','FLUX','ERROR','MASK'])
    spectrum=spectrum[spectrum['WAVE']>wmin]
    spectrum=spectrum[spectrum['WAVE']<wmax]
    spectrum=spectrum[spectrum['MASK']==0]
    spectrum['ERROR']=1/np.sqrt(spectrum['ERROR'])*1e16
    # mean_flux, median_flux, std_flux = sigma_clipped_stats(spectrum['FLUX'], sigma_upper=2,sigma_lower=3)
    # spectrum['FLUX']-=median_flux
    if plot==True:
        plt.clf()
        plt.plot(spectrum['WAVE'],spectrum['FLUX'])
        plt.plot(spectrum['WAVE'],spectrum['ERROR'])
    return spectrum


def generate_sample_spectrum(delta_wave=1.0, wave_start=4000, wave_end=7000, 
                             n_lines=5, noise_level=0.1, n_doublets=2,
                             fwhm_range=(5, 20), amplitude_range=(0.5, 2.0),
                             doublet_separation_range=(2.0, 10.0)):
    """
    Generate a sample spectrum with emission lines and return both the spectrum
    and a table of line parameters.
    
    Parameters
    ----------
    delta_wave : float
        Wavelength spacing in Angstroms.
    wave_start, wave_end : float
        Wavelength range of the spectrum in Angstroms.
    n_lines : int
        Number of isolated emission lines.
    noise_level : float
        Standard deviation of Gaussian noise added to the continuum.
    n_doublets : int
        Number of emission line doublets to include.
    fwhm_range : tuple (min, max)
        Range of FWHM values in Angstroms for the emission lines.
    amplitude_range : tuple (min, max)
        Range of peak amplitudes for the emission lines.
    doublet_separation_range : tuple (min, max)
        Range of separations (in Angstroms) between doublet components.
    
    Returns
    -------
    spec_table : astropy.table.Table
        Simulated spectrum with columns: WAVE, FLUX, ERROR.
    line_table : astropy.table.Table
        Table of line properties with columns: TYPE, WAVELENGTH, FLUX, FWHM.
    """
    
    wave = np.arange(wave_start, wave_end + delta_wave, delta_wave)
    n_points = len(wave)
    
    continuum = 1.0
    noise = np.random.normal(0, noise_level, n_points)
    flux = continuum + noise
    
    line_table = Table(names=['TYPE', 'WAVELENGTH', 'FLUX', 'FWHM'], dtype=['str', 'f8', 'f8', 'f8'])
    
    # Generate single emission lines
    line_centers = np.random.uniform(wave_start + 200, wave_end - 200, n_lines)
    line_strengths = np.random.uniform(*amplitude_range, n_lines)
    line_widths = np.random.uniform(*fwhm_range, n_lines)
    
    for i in range(n_lines):
        center = line_centers[i]
        strength = line_strengths[i]
        fwhm = line_widths[i]
        sigma = fwhm / 2.355
        
        emission = strength * np.exp(-0.5 * ((wave - center) / sigma) ** 2)
        flux += emission
        line_table.add_row(['single', center, strength, fwhm])
    
    # Generate doublets
    if n_doublets > 0:
        doublet_centers = np.random.uniform(wave_start + 300, wave_end - 300, n_doublets)
        
        for i in range(n_doublets):
            center = doublet_centers[i]
            fwhm = np.random.uniform(*fwhm_range)
            sigma = fwhm / 2.355
            separation = np.random.uniform(*doublet_separation_range)
            
            strength1 = np.random.uniform(*amplitude_range)
            strength2 = np.random.uniform(*amplitude_range)
            
            emission1 = strength1 * np.exp(-0.5 * ((wave - center) / sigma) ** 2)
            emission2 = strength2 * np.exp(-0.5 * ((wave - (center + separation)) / sigma) ** 2)
            
            flux += emission1 + emission2
            
            line_table.add_row(['doublet', center, strength1, fwhm])
            line_table.add_row(['doublet', center + separation, strength2, fwhm])
    
    error = np.sqrt(np.abs(flux)) * noise_level
    
    spec_table = Table()
    spec_table['WAVE'] = wave
    spec_table['FLUX'] = flux
    spec_table['ERROR'] = error
    
    return spec_table, line_table




def do_one(rss_file,fib,wmin,wmax,plot_results=False):

    spec=fits2table(fib=fib,wmin=wmin,wmax=wmax,filename=rss_file)

    if len(spec)==0:
        print('Nothing to do')
        return []

    print(spec)

    results=find_emission_lines(spec_table=spec, min_snr=0.5, typical_fwhm=0.6, 
                        plot_results=plot_results, continuum_percentile=20,
                        resolve_close_lines=True, prominence=0.02,name='Fiber %5d' % fib)
    results['fib']=fib

    if plot_results:
        dir_name='Peaks.%s' % (date.today().strftime("%y%m%d"))
        if os.path.isdir(dir_name)==False:
            os.mkdir(dir_name)
        name='%s/p%05d_%d_%d.png' % (dir_name,fib,wmin,wmax)
        plt.savefig(name)
    return results

def steer(argv):
    '''
    peaks.py -h -wmin 6720 -wmax 6740  -all   filename.fits   one or more fiber numbers
    '''
    filename=''
    wmin=6720
    wmax=6740
    fib_no=[]
    xplot=True
    all_results=[]
    do_all=False

    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-wmin':
            i+=1
            wmin=eval(argv[i])
        elif argv[i]=='-wmax':
            wmax=eval(argv[i])
        elif argv[i]=='-all':
            do_all=True       
        elif argv[i][0]=='-':
            print('Error:Unknown option: ',argv)
            return
        elif argv[i].count('fits'):
            filename=argv[i]
        else:
            fib_no.append(int(argv[i]))
        i+=1

    if do_all==True:
        x=fits.open(filename)
        f=x['FLUX'].data
        fib_no=range(len(f))
        x.close()

    i=0
    while i<len(fib_no):
        results=do_one(filename,fib_no[i],wmin,wmax,xplot)
        if len(results):
            all_results.append(results)
            print('Got %d' % fib_no[i])
        else:
            print('failed %d' % fib_no[i])

        i+=1

    if len(all_results):
        all_results=vstack(all_results)
        out_name='Peaks_%04d_%04d.%s.txt' % (wmin,wmax,date.today().strftime("%y%m%d"))
        print('Writing %d to %s ' % (len(all_results),out_name ))
        all_results.write(out_name,format='ascii.fixed_width_two_line',overwrite=True)




    

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
#!/usr/bin/env python
# coding: utf-8

