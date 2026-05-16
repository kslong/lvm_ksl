#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Fit Gaussians to a fixed set of nebular emission lines and airglow lines
in the science fibers of a sky-subtracted LVM SFrame exposure, producing
one output row per fiber.


Command line usage::

    sky_gaussfit.py [-lmc] [-smc] [-v vel] [-out root] [-np nproc] filename [filename ...]

    where

    -lmc         Apply the LMC radial velocity (~262 km/s) when fitting nebular lines
    -smc         Apply the SMC radial velocity (~146 km/s) when fitting nebular lines
    -v vel       Apply an arbitrary radial velocity (km/s) when fitting nebular lines
    -out root    Set the root name for the output file (default: derived from input filename)
    -np nproc    Number of parallel processes for fiber fitting (default: 8)
    filename     One or more SFrame FITS files or ASCII spectrum files (WAVE, FLUX columns)


Description:

For each good science fiber (telescope == 'Sci', fibstatus == 0) the script fits
single Gaussians to the following lines:

    Nebular lines (wavelengths shifted by the supplied velocity):
        [OIII] 4958.911  (oiii_a)
        [OIII] 5006.843  (oiii_b)
        [OI]   6300.309  (oi_a)
        [OI]   6363.783  (oi_b)
        [NII]  6548.04   (nii_a)
        Ha     6562.80   (ha)
        [NII]  6583.46   (nii_b)
        [SII]  6716.440  (sii_a)
        [SII]  6730.815  (sii_b)

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
260516 ksl Add 12 airglow lines (ESO UVES wavelengths); sky6300 fit at fixed wavelength alongside nebular OI
260516 ksl Add nii_a (6548.04) and nii_b (6583.46); correct SII to 6716.440 and 6730.815
260516 ksl Correct oi_a/oi_b: use ESO wavelengths (6300.309, 6363.783) with velocity shift
260516 ksl Add oiii_a (4958.911) and oiii_b (5006.843); refactor to NEBULAR_LINES/SKY_LINES constants; pre-trim and pre-compute per-line indices in do_all

'''

import sys
import os
import multiprocessing
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


# Nebular emission lines: (name, rest_wavelength, window_min, window_max) in Angstroms.
# All values are multiplied by zz = 1 + vel/c at runtime to apply the radial velocity.
NEBULAR_LINES = [
    ('oiii_a', 4958.911,     4953.,  4964.),
    ('oiii_b', 5006.843,     5001.,  5012.),
    ('oi_a',   6300.308594,  6295.,  6305.),
    ('oi_b',   6363.782715,  6358.,  6368.),
    ('nii_a',  6548.04,      6543.,  6553.),
    ('ha',     6562.80,      6555.,  6570.),
    ('nii_b',  6583.46,      6578.,  6593.),
    ('sii_a',  6716.440,     6706.,  6726.),
    ('sii_b',  6730.815,     6721.,  6741.),
]

# Airglow lines: (name, center_wavelength, window_min, window_max) in Angstroms.
# Fitted at fixed, unshifted wavelengths (ESO UVES sky spectrum atlas).
SKY_LINES = [
    ('sky5577', 5577.34668,  5572.,  5582.),
    ('sky6300', 6300.308594, 6295.,  6305.),
    ('sky6363', 6363.782715, 6358.,  6368.),
    ('sky6533', 6533.04,     6528.,  6538.),
    ('sky6553', 6553.0,      6549.,  6556.),
    ('sky6577', 6577.2,      6572.,  6582.),
    ('sky6912', 6912.623,    6907.,  6917.),
    ('sky6923', 6923.220,    6918.,  6928.),
    ('sky6939', 6939.521,    6934.,  6944.),
    ('sky7358', 7358.680176, 7353.,  7363.),
    ('sky7392', 7392.209961, 7387.,  7397.),
    ('sky7914', 7913.717773, 7908.,  7918.),
    ('sky8344', 8344.613281, 8339.,  8349.),
    ('sky8399', 8399.175781, 8394.,  8404.),
    ('sky8827', 8827.112305, 8822.,  8832.),
    ('sky8988', 8988.383789, 8983.,  8993.),
    ('sky9552', 9552.546875, 9547.,  9557.),
    ('sky9719', 9719.838867, 9714.,  9724.),
]


def do_one(spectrum_table, vel=0., xplot=False):
    '''
    Completely process a single spectrum
    '''
    clean()
    zz = 1.0 + vel / 3e5
    records = []

    for name, center, wmin, wmax in NEBULAR_LINES:
        try:
            results, xspec = fit_gaussian_to_spectrum(
                spectrum_table, line=name,
                init_wavelength=zz*center, init_fwhm=1.,
                wavelength_min=zz*wmin, wavelength_max=zz*wmax)
            records.append(results)
            if xplot:
                save_fit(name, xspec)
        except Exception as e:
            print('Fitting %s: An exception occurred: %s' % (name, e))

    for name, center, wmin, wmax in SKY_LINES:
        try:
            results, xspec = fit_gaussian_to_spectrum(
                spectrum_table, line=name,
                init_wavelength=center, init_fwhm=1.,
                wavelength_min=wmin, wavelength_max=wmax)
            records.append(results)
            if xplot:
                save_fit(name, xspec)
        except Exception as e:
            print('Fitting %s: An exception occurred: %s' % (name, e))

    try:
        ztab = hstack(records)
        return ztab
    except:
        print('Nothing fit for this spectrum')
        return []


def _fit_one_fiber(args):
    '''
    Worker function for parallel fiber fitting.  Must be at module level
    so multiprocessing can pickle it.

    args: (line_specs, line_waves, line_idxs, trim_flux_row, fiberid, ra, dec)
    where line_specs is a list of (name, center, wmin, wmax) with velocity
    already applied, line_waves[name] is the pre-sliced wavelength sub-array,
    and line_idxs[name] is the integer index array into trim_flux_row.

    Returns: (rtab_or_None, fiberid, ra, dec, status)
    where status is 'ok', 'nothing', 'nan', or 'error'.
    '''
    line_specs, line_waves, line_idxs, trim_flux_row, fiberid, ra, dec = args
    if check_for_nan(trim_flux_row):
        return None, int(fiberid), float(ra), float(dec), 'nan'
    records = []
    for name, center, wmin, wmax in line_specs:
        idx = line_idxs[name]
        if len(idx) == 0:
            continue
        sub_table = Table([line_waves[name], trim_flux_row[idx]],
                          names=['WAVE', 'FLUX'])
        try:
            results, xspec = fit_gaussian_to_spectrum(
                sub_table, line=name,
                init_wavelength=center, init_fwhm=1.,
                wavelength_min=wmin, wavelength_max=wmax)
            records.append(results)
        except Exception as e:
            pass
    if not records:
        return None, int(fiberid), float(ra), float(dec), 'nothing'
    rtab = hstack(records)
    rtab['fiberid'] = int(fiberid)
    rtab['ra']      = float(ra)
    rtab['dec']     = float(dec)
    return rtab, int(fiberid), float(ra), float(dec), 'ok'


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

def do_all(filename='data/lvmSFrame-00009088.fits', vel=0.0, outname='', nproc=8):
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

    slittab=Table(x['SLITMAP'].data)
    good=scifib(slittab,'science')

    # Build the full list of line specs with velocity applied to nebular lines.
    zz = 1.0 + vel / 3e5
    line_specs = (
        [(name, zz*center, zz*wmin, zz*wmax) for name, center, wmin, wmax in NEBULAR_LINES] +
        [(name, center, wmin, wmax) for name, center, wmin, wmax in SKY_LINES]
    )

    # Pre-trim the wavelength array to only the range covered by any line window.
    trim_wmin = min(s[2] for s in line_specs) - 2.
    trim_wmax = max(s[3] for s in line_specs) + 2.
    trim_mask = (wave >= trim_wmin) & (wave <= trim_wmax)
    trim_wave = wave[trim_mask]

    # Pre-compute per-line index arrays and wave sub-arrays once for the whole file.
    line_idxs = {}
    line_waves = {}
    for name, center, wmin, wmax in line_specs:
        idx = np.where((trim_wave >= wmin) & (trim_wave <= wmax))[0]
        line_idxs[name] = idx
        line_waves[name] = trim_wave[idx]

    args_list = [
        (line_specs, line_waves, line_idxs,
         flux[good['fiberid'][i]-1][trim_mask],
         good['fiberid'][i], good['ra'][i], good['dec'][i])
        for i in range(len(good))
    ]

    print('Fitting %d fibers using %d processes' % (len(good), nproc))
    if nproc > 1:
        with multiprocessing.Pool(nproc) as pool:
            raw = pool.map(_fit_one_fiber, args_list)
    else:
        raw = [_fit_one_fiber(a) for a in args_list]

    records = []
    for rtab, fiberid, ra, dec, status in raw:
        if status == 'nan':
            print('Too many nans for fiber %d at %.2f %.2f' % (fiberid, ra, dec))
        elif status == 'nothing':
            print('Nothing fit for fiber %d at %.2f %.2f' % (fiberid, ra, dec))
        elif status == 'error':
            print('Error fitting fiber %d at %.2f %.2f' % (fiberid, ra, dec))
        else:
            records.append(rtab)

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
    nproc=8

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
        elif argv[i]=='-np':
            i+=1
            nproc=int(argv[i])
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
        results=do_all(one_file,vel,outname,nproc)
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
