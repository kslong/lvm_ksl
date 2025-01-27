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
    fit_table = Table([x, y,y_fit], names=('WAVE', 'FLUX','Fit'))
    
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
        qtab[one_col].format='.3e'

    record=[amplitude1,amp_err1,mean1,wave_err1,fwhm,fwhm_err1,background,rmse,amplitude2,amp_err2,mean2,wave_err2,fwhm,fwhm_err2,background,rmse]
    qtab.add_row(record)
    
    return qtab, fit_table

def save_fit(line='oi',ftab=None,xdir='Gauss_dir'):
    '''
    Save the fit to a directory
    '''
    if os.path.isdir(xdir)==False:
        os.makedirs(xdir,exist_ok=True)
    ftab.write('%s/%s.txt' % (xdir,line), format='ascii.fixed_width_two_line',overwrite=True)

def do_one(spectrum_table,vel=0.,xplot=False,xfit=False):
    '''
    Completely process a single spectrum
    '''


    
    zz=1.+ (vel/3e5)

    records=[]
    try:
        results,xspec=fit_double_gaussian_to_spectrum(spectrum_table, line='oii',init_wavelength1=zz*3726.092, init_wavelength2=zz*3729.875, init_fwhm=1, wavelength_min=zz*3717, wavelength_max=zz*3737, plot=xplot)
        records.append(results)
        save_fit('oii',xspec)
    except Exception as e:
        print(f"Fitting OII; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")

    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='hb',init_wavelength=zz*4861, init_fwhm=1., wavelength_min=zz*4855, wavelength_max=zz*4870, plot=xplot)
        records.append(results)
        save_fit('hb',xspec)
    except Exception as e:
        print(f"Fitting Hb; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='oiii_a',init_wavelength=zz*4959, init_fwhm=1., wavelength_min=zz*4949, wavelength_max=zz*4969, plot=xplot)
        records.append(results)
        save_fit('oiii_a',xspec)
    except Exception as e:
        print(f"Fitting [OIII]4959; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='oiii_b', init_wavelength=zz*5007, init_fwhm=1., wavelength_min=zz*4997, wavelength_max=zz*5017, plot=xplot)
        save_fit('oiii_b',xspec)
        records.append(results)
    except Exception as e:
        print(f"Fitting [OIII]5007; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")

    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='ha', init_wavelength=zz*6563, init_fwhm=1., wavelength_min=zz*6555, wavelength_max=zz*6570, plot=xplot)
        save_fit('ha',xspec)
        records.append(results)
    except Exception as e:
        print(f"Fitting Ha An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")

    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='nii_a',init_wavelength=zz*6548, init_fwhm=1., wavelength_min=zz*6538, wavelength_max=zz*6558, plot=xplot)
        records.append(results)
        save_fit('nii_a',xspec)
    except Exception as e:
        print(f"Fitting [NII]6548:  An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='nii_b',init_wavelength=zz*6584, init_fwhm=1., wavelength_min=zz*6574, wavelength_max=zz*6594, plot=xplot)
        records.append(results)
        save_fit('nii_b',xspec)
    except Exception as e:
        print(f"Fitting [NII]6584:  An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
            
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='sii_a',init_wavelength=zz*6716, init_fwhm=1., wavelength_min=zz*6706, wavelength_max=zz*6726, plot=xplot)
        records.append(results)
        save_fit('sii_a',xspec)
    except Exception as e:
        print(f"Fitting [SII]6716:  An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='sii_b',init_wavelength=zz*6731, init_fwhm=1., wavelength_min=zz*6721, wavelength_max=zz*6741, plot=xplot)
        records.append(results)
        save_fit('sii_b',xspec)
    except Exception as e:
        print(f"Fitting [SII]6731:  An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    

    #250110 - Add more lines
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='heii',init_wavelength=zz*4686, init_fwhm=1., wavelength_min=zz*4666, wavelength_max=zz*4706, plot=xplot)
        records.append(results)
        save_fit('heii',xspec)
    except Exception as e:
        print(f"Fitting heii; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='hei',init_wavelength=zz*5876, init_fwhm=1., wavelength_min=zz*5856, wavelength_max=zz*5896, plot=xplot)
        records.append(results)
        save_fit('hei',xspec)
    except Exception as e:
        print(f"Fitting hei; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    

    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='oiii_4363',init_wavelength=zz*4363, init_fwhm=1., wavelength_min=zz*4343, wavelength_max=zz*4383, plot=xplot)
        records.append(results)
        save_fit('oiii_4363',xspec)
    except Exception as e:
        print(f"Fitting hei; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    

    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='siii_a',init_wavelength=zz*9068, init_fwhm=1., wavelength_min=zz*9048, wavelength_max=zz*9088, plot=xplot)
        records.append(results)
        save_fit('siii_a',xspec)
    except Exception as e:
        print(f"Fitting siii_a; An exception occurred: {e}")
        print(f"Exception type: {type(e).__name__}")
    
    
    try:
        results,xspec=fit_gaussian_to_spectrum(spectrum_table, line='siii_b',init_wavelength=zz*9531, init_fwhm=1., wavelength_min=zz*9511, wavelength_max=zz*9551, plot=xplot)
        records.append(results)
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
    
def do_individual(filenames,vel,outname,xplot=False,xfit=False):
    '''
    This is to process individual spectra from an astropy table
    containing a WAVE and FLUX column
    '''
    xresults=[]
    xbad=[]
    xgood=[]
    if len(filenames)==1:
        xfit=True

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
            results=do_one(xtab,vel,xplot,xfit)
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
    '''
    files=glob('%s/*.txt' % qdir)
    files=np.sort(files)
    print(len(files))
    plt.figure(1,(12,12))
    plt.clf()
    i=0
    while i<len(files):
        plt.subplot(4,4,i+1)
        one_file=files[i]
        xtab=ascii.read(one_file)
        word=one_file.split('/')
        name=word[-1].replace('.txt','')
        print(name)
        plot_one(xtab,name)
        i+=1
    plt.tight_layout()

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
