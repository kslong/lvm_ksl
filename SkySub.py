#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Carry out simple sky subtraction of a science exposure
and a appropriate sky exposure.


Command line usage (if any):

    usage: SkySub.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

240330 ksl Coding begun

'''

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table,vstack
from scipy.optimize import bisect
from astropy.stats import sigma_clip



def fit_func(r, sci, sky):
    '''
    The function that is minimized, basically
    the absolute value of the difference in science and
    sky squared, which basically weights regions with
    the brightest sky lines the most.
    '''

    delta=sci-r*sky

    xx=sci*delta
    xx=np.abs(xx)

    xresult=np.sum(xx)
    xnorm=np.dot(sky,sky)

    
    return xresult/xnorm



def ksl_bisection(func, a, b, tol=1e-2, args=(), maxiter=30):
    """
    Custom implementation of bisection method to find the minimum of a function within an interval.
    
    Parameters:
        func (callable): The objective function.
        a (float): The lower bound of the interval.
        b (float): The upper bound of the interval.
        tol (float, optional): The tolerance for the minimum value. Default is 1e-6.
        maxiter (int, optional): The maximum number of iterations. Default is 100.
    
    Returns:
        float: The x-coordinate of the estimated minimum.
    """


    j=0
    while j<maxiter:
        interval = [a, (3./4.*a+1/4.*b), (a + b) / 2, 1/4* a + 3/4*b, b]
        print(interval)
        values = [func(x,*args) for x in interval]
        print(values)
        min_index = values.index(min(values))
        print(min_index)


        if min_index == 0:
            a, b = interval[0], interval[2]
        elif min_index == 1:
            a, b = interval[0], interval[2]
        elif min_index == 2:
            a, b = interval[1], interval[3]
        elif min_index == 3:
            a, b = interval[2], interval[4]
        elif min_index == 4:
            a, b = interval[2], interval[4]
        print(a,b)
        if abs(a-b)<tol:
            print('Accuracy achieved on interation %d' % j)
            break
        j+=1


    # Return the midpoint after the maximum number of iterations
    return (a + b) / 2



def polynomial_fit_with_outliers(spectrum_table, degree=3, sigma_lower=3, sigma_upper=3, grow=0, max_iter=10):
    """
    Perform a polynomial fit to the total spectrum while iteratively removing outliers.

    Parameters:
        spectrum_table (astropy.table.Table): Table containing columns for wavelength, total flux, line flux, and continuum flux.
        degree (int): Degree of the polynomial fit.
        sigma_lower (float): Lower sigma value for sigma-clipping to identify outliers.
        sigma_upper (float): Upper sigma value for sigma-clipping to identify outliers.
        grow (int): Number of pixels to grow the mask of clipped values.
        max_iter (int): Maximum number of iterations for outlier removal.

    Returns:
        astropy.table.Table: Table containing columns for wavelength, total flux, line flux, continuum flux, fitted spectrum, and mask.
    """
    # Copy input spectrum table
    output_table = spectrum_table.copy()

    # Initial fit
    coefficients = np.polyfit(spectrum_table['WAVE'], spectrum_table['FLUX'], degree)
    fitted_flux = np.polyval(coefficients, spectrum_table['WAVE'])

    # Iterate until convergence or max_iter
    for _ in range(max_iter):
        # Compute residuals
        residuals = spectrum_table['FLUX'] - fitted_flux

        # Perform sigma clipping
        mask = sigma_clip(residuals, sigma_lower=sigma_lower, sigma_upper=sigma_upper, grow=grow).mask

        # Update masked spectrum table
        masked_spectrum = spectrum_table[~mask]

        # Perform polynomial fit
        coefficients = np.polyfit(masked_spectrum['WAVE'], masked_spectrum['FLUX'], degree)

        # Evaluate polynomial fit on the wavelength grid
        fitted_flux = np.polyval(coefficients, spectrum_table['WAVE'])

    # Add fitted spectrum and mask columns to output table
    output_table['CONT'] = fitted_flux
    output_table['MASK'] = mask

    return output_table




def sky_sub_lines(sci_file='Sci_5554.fits',sky_file='SkyE_5554.fits',wmin=3600,wmax=9000):
    sci=fits.open(sci_file)
    sky=fits.open(sky_file)
    sci_tab=Table(sci[1].data)
    sky_tab=Table(sky[1].data)
    # sky_tab.info()
    plt.figure(1,(6,6))    

    xsci=sci_tab[sci_tab['WAVE']>wmin]
    xsci=xsci[xsci['WAVE']<wmax]

    xsky=sky_tab[sky_tab['WAVE']>wmin]
    xsky=xsky[xsky['WAVE']<wmax]

    xsci= polynomial_fit_with_outliers(xsci, degree=4, sigma_lower=3, sigma_upper=1, grow=5)
    xsky= polynomial_fit_with_outliers(xsky, degree=4, sigma_lower=3, sigma_upper=1, grow=5)
    
    xnorm=np.dot(xsky['FLUX'],xsky['FLUX'])

    xsci['LINES']=xsci['FLUX']-xsci['CONT']
    xsky['LINES']=xsky['FLUX']-xsky['CONT']
    
    rmin=0.5
    rmax=1.5

    minimum =ksl_bisection(fit_func, rmin, rmax, tol=0.001, maxiter=8,args=(xsci['LINES'],xsky['LINES']))


    print('Results %f' % (minimum))

    plt.figure(1,(8,6))
    plt.subplot(4,1,1)
    plt.plot(xsci['WAVE'],xsci['FLUX'],label='Science')
    plt.legend()
    plt.xlim(wmin,wmax)
    plt.tight_layout()
    
    plt.subplot(4,1,2)
    plt.plot(xsky['WAVE'],xsky['FLUX'],label='Sky')

    plt.xlim(wmin,wmax)
    plt.legend()
    plt.tight_layout()

    plt.subplot(4,1,3)
    plt.plot(xsci['WAVE'],xsci['FLUX']-minimum*xsky['LINES'],label='Sky-subtracted')
    plt.legend()
    plt.xlim(wmin,wmax)
    plt.tight_layout()

    plt.subplot(4,1,4)
    plt.plot(xsci['WAVE'],xsci['FLUX']-xsci['CONT']-minimum*xsky['LINES'],label='Sky and Cont subtracted')
    plt.legend()
    ymin,ymax=plt.ylim()
    if ymin<-2e-14:
        ymin=-2e-14
    if ymax>1e-13:
        ymax=1e-13
    plt.ylim(ymin,ymax)
    plt.tight_layout()

    outname='%s_%s_lines.png' % (sci_file.replace('.fits',''),sky_file.replace('.fits',''))
    plt.savefig(outname)




def sky_sub_all(sci_file='Sci_4852.fits',sky_file='SkyE_4852.fits',wmin=6200,wmax=6450):
    sci=fits.open(sci_file)
    sky=fits.open(sky_file)
    sci_tab=Table(sci[1].data)
    sky_tab=Table(sky[1].data)
    plt.figure(1,(6,6))    

    xsci=sci_tab[sci_tab['WAVE']>wmin]
    xsci=xsci[xsci['WAVE']<wmax]

    xsky=sky_tab[sky_tab['WAVE']>wmin]
    xsky=xsky[xsky['WAVE']<wmax]
    
    xnorm=np.dot(xsky['FLUX'],xsky['FLUX'])
    
    rmin=0.5
    rmax=1.5

    minimum =ksl_bisection(fit_func, rmin, rmax, tol=0.001, maxiter=8,args=(xsci['FLUX'],xsky['FLUX']))


    print('Results %f' % (minimum))

    plt.figure(1,(8,6))
    plt.subplot(4,1,1)
    plt.plot(xsci['WAVE'],xsci['FLUX'],label='Science')
    plt.legend()
    plt.xlim(wmin,wmax)
    plt.tight_layout()
    
    plt.subplot(4,1,2)
    plt.plot(xsky['WAVE'],xsky['FLUX'],label='Sky')
    plt.legend()
    plt.xlim(wmin,wmax)
    plt.tight_layout()

    plt.subplot(4,1,3)
    plt.plot(xsci['WAVE'],xsci['FLUX']-minimum*xsky['FLUX'],label='Sky-subtracted')
    plt.legend()
    plt.xlim(wmin,wmax)
    ymin,ymax=plt.ylim()
    if ymin<-2e-14:
        ymin=-2e-14
    if ymax>1e-13:
        ymax=1e-13
    plt.ylim(ymin,ymax)
    plt.tight_layout()

    outname='%s_%s_all.png' % (sci_file.replace('.fits',''),sky_file.replace('.fits',''))
    plt.savefig(outname)





def steer(argv):
    '''
    The simple directs the flow
    '''

    science=''
    sky=''

    lines_only=True

    i=1
    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-all':
            lines_only=False
        elif argv[i][0]=='-':
            print('Error: Cannot parse command line: ',argv)
            return
        elif science=='':
            science=argv[i]
        elif sky=='':
            sky=argv[i]
        else:
            print('Error: Too many arguments in command line: ',argv)
            return
        i+=1

    if sky=='':
        print('Error: not enoubh argmeents: ',argv)

    if lines_only:
        sky_sub_lines(sci_file=science,sky_file=sky,wmin=3600,wmax=9000)
    else:
        sky_sub_all(sci_file=science,sky_file=sky,wmin=3600,wmax=9000)



    return


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
