#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Fit Gaussians to a fixed set of airglow lines in each LVM exposure,
producing one output row per exposure.  Two input modes are supported:
individual lvmCFrame files (supplied directly, via a file list, or
located via drpall), or a pre-built SummarizeCframe FITS file.  Both
modes fit the same quantity — the median spectrum across all good science
fibers — and produce equivalent results for the same exposures.


Command line usage:

    gauss_offset.py [-ext FLUX] [-out root] [-ver drp_ver] [-file xfile]
                    [exp_start [exp_stop [delta]]]
                    filename [filename ...]

    where

    -ext extname  Extension to fit: FLUX (default), SKY_EAST, or SKY_WEST
    -out root     Root name for the output file.  If omitted the output
                  filename is derived from the inputs (see Output section).
    -ver drp_ver  DRP version for drpall lookup (default 1.2.1).
    -file xfile   Read lvmCFrame filenames from a table with a column named
                  filename or Filename.  SFrame paths are converted to CFrame
                  automatically.  These are added to any CFrame files from
                  exp_start/exp_stop.
    exp_start     First exposure number.  Combined with exp_stop, queries
                  drpall to build a list of CFrame files to process.
    exp_stop      Last exposure number (inclusive).
    delta         Use every Nth exposure in [exp_start, exp_stop] (default 1).
    filename      One or more SummarizeCframe FITS files (XCframe_*.fits).
                  Direct filename arguments select SummarizeCframe mode
                  (see Description).


Description:

In both modes the script fits single Gaussians to the following airglow
lines and returns one output row per exposure:

    sky5577   5577.34 A
    sky6300   6300.31 A
    sky6363   6363.78 A
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
level, and RMSE.

Two input modes are supported.  Both fit the same quantity — the median
spectrum across all good science fibers — and give equivalent results for
the same set of exposures.

  SummarizeCframe mode (triggered by direct .fits filename arguments):
    A SummarizeCframe file (XCframe_*.fits) stores, for each exposure, a
    pre-computed median science-fiber spectrum as one row of its FLUX
    extension (or SKY_EAST / SKY_WEST for sky-telescope spectra).  This
    script reads each row directly and fits it.  Use this mode when you
    have already run SummarizeCframe or one of its variants.

  CFrame mode (triggered by -file or exp_start/exp_stop):
    Individual lvmCFrame files are processed one at a time.  For each
    file, the median spectrum is computed on the fly across all good
    science fibers (telescope == 'Sci', fibstatus == 0), with bad pixels
    masked via the MASK extension.  Gaussians are then fit to that median.
    SFrame paths in the input table are converted to CFrame paths
    automatically; drpall is used as a fallback if a path cannot be
    resolved directly.  Use this mode when you have a drpall-derived file
    list or want to skip the SummarizeCframe step.

Primary routines:

    do_one        - fit all airglow lines for a single spectrum table
    doit          - fit one CFrame file (median science-fiber spectrum)
    do_all_files  - iterate over a list of CFrame files
    do_all        - iterate over rows in a SummarizeCframe FITS file


Output:

A FITS binary table with one row per exposure.  Columns cover the fit
parameters (flux, wave, fwhm, back, rmse) for each line together with
expnum and mjd.

CFrame mode: output is named Gauss_<ext>.<YYMMDD>.fits by default.

SummarizeCframe mode: output is named <input_stem>_gauss.<ext>.fits by
default.

If -out root is supplied the output is written to root.fits (adding .fits
if not already present).


Notes:

Wavelengths are from the ESO UVES sky spectrum atlas.

Because CFrame spectra are before sky subtraction, the airglow lines are
bright and well-suited to measuring centroid wavelengths and fluxes as a
function of exposure, revealing wavelength-calibration drifts or changes
in sky brightness over a night or survey.

History:

240604 ksl Coding begun (as sky_gaussfit.py)
260504 ksl Adapted for SummarizeCframe input; restricted to airglow lines only
260504 ksl FITS output; selectable extension (-ext)
260515 ksl Add -file, -ver, exp_start/exp_stop, and CFrame direct-file mode

'''

import sys
import os
import re
import datetime
import numpy as np
from astropy.table import Table, vstack, hstack
from astropy.io import fits
from astropy.time import Time

from lvm_ksl.lvm_gaussfit import fit_gaussian_to_spectrum, check_for_nan


# Wavelengths taken from
# https://www.eso.org/observing/dfo/quality/UVES/pipeline/sky_spectrum.html
AIRGLOW_LINES = [
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
]


def build_drp_lookup(drp_ver='1.2.1'):
    '''
    Build a dict mapping exposure number to full CFrame file path.

    Reads drpall via SummarizeCframe.read_drpall(), prepends the
    machine-specific top directory, and converts SFrame to CFrame paths.
    Returns an empty dict if drpall cannot be read.
    '''
    try:
        from SummarizeCframe import find_top, read_drpall
    except ImportError:
        print('Warning: could not import SummarizeCframe; drpall lookup unavailable.')
        return {}
    topdir = find_top()
    if not topdir:
        return {}
    drp_tab = read_drpall(drp_ver)
    if drp_tab is None or len(drp_tab) == 0:
        return {}
    lookup = {}
    for row in drp_tab:
        loc = str(row['location']).replace('SFrame', 'CFrame')
        lookup[int(row['expnum'])] = os.path.join(topdir, loc)
    return lookup


def resolve_filename(filename, lookup=None):
    '''
    Return the full path to a CFrame file, or None if not found.

    Converts any SFrame path to CFrame first.  If the file exists at that
    path it is returned immediately.  Otherwise the exposure number is
    extracted from the filename and used to look up the full path in lookup
    (a dict built by build_drp_lookup()).
    '''
    filename = filename.replace('SFrame', 'CFrame')
    if os.path.isfile(filename):
        return filename
    if lookup:
        m = re.search(r'(\d+)\.fits', os.path.basename(filename))
        if m:
            expnum = int(m.group(1))
            if expnum in lookup and os.path.isfile(lookup[expnum]):
                return lookup[expnum]
    return None


def do_one(spectrum_table):
    '''
    Fit all airglow lines for a single spectrum supplied as an astropy Table
    with WAVE and FLUX columns.  Returns an hstacked single-row Table of fit
    results, or an empty list if nothing could be fit.
    '''
    records = []

    for xname, xwcen, xwmin, xwmax in AIRGLOW_LINES:
        try:
            results, xspec = fit_gaussian_to_spectrum(
                spectrum_table, line=xname,
                init_wavelength=xwcen, init_fwhm=1.,
                wavelength_min=xwmin, wavelength_max=xwmax)
            records.append(results)
        except Exception as e:
            print('Fitting %s: An exception occurred: %s' % (xname, e))

    if len(records) == 0:
        print('Nothing fit for this spectrum')
        return []

    return hstack(records)


def doit(filename, ext='FLUX'):
    '''
    Fit airglow lines for a single CFrame file.

    Computes the median spectrum across all good science fibers
    (telescope == 'Sci', fibstatus == 0) using the MASK extension,
    then fits Gaussians to each airglow line.

    Returns a single-row astropy Table, or None if the file cannot be
    processed.
    '''
    try:
        x = fits.open(filename)
    except Exception:
        print('Error: Could not open %s' % filename)
        return None

    expnum  = x[0].header.get('EXPOSURE', -1)
    obstime = x[0].header.get('OBSTIME', '')
    mjd     = Time(obstime, format='isot', scale='utc').mjd if obstime else np.nan

    wave    = x['WAVE'].data
    slitmap = Table(x['SLITMAP'].data)
    sci     = slitmap[(slitmap['telescope'] == 'Sci') & (slitmap['fibstatus'] == 0)]
    rows    = sci['fiberid'] - 1   # fiberid is 1-based

    try:
        flux_2d = x[ext].data[rows]
    except KeyError:
        print('Error: extension %s not found in %s' % (ext, filename))
        return None

    mask_2d  = x['MASK'].data[rows]
    flux_ma  = np.ma.masked_array(flux_2d, mask_2d.astype(bool))
    med_flux = np.ma.median(flux_ma, axis=0).data

    if check_for_nan(med_flux):
        print('Too many NaNs in %s, skipping' % filename)
        return None

    one_spec = Table([wave, med_flux], names=['WAVE', 'FLUX'])
    rtab = do_one(one_spec)
    if not len(rtab):
        return None

    rtab['expnum'] = expnum
    rtab['mjd']    = mjd
    return rtab


def do_all_files(files, ext='FLUX', outname='', drp_ver='1.2.1'):
    '''
    Fit airglow lines for a list of CFrame files.

    Parameters
    ----------
    files : list of str
        Paths to lvmCFrame FITS files.  SFrame paths are converted
        automatically; drpall is used as a fallback if a path does not
        exist at its literal location.
    ext : str
        Extension to fit (FLUX, SKY_EAST, or SKY_WEST).
    outname : str
        Output filename root.  If empty, Gauss_<ext>.<YYMMDD>.fits is used.
    drp_ver : str
        DRP version used for the fallback drpall path lookup.
    '''
    lookup  = build_drp_lookup(drp_ver)
    records = []

    for one in files:
        full = resolve_filename(one, lookup)
        if full is None:
            print('Skipping %s: file not found' % one)
            continue
        print('Processing %s' % full)
        rtab = doit(full, ext)
        if rtab is not None and len(rtab) > 0:
            records.append(rtab)

    if not records:
        print('No spectra successfully fit.')
        return

    results = vstack(records)

    if outname == '':
        today   = datetime.date.today().strftime('%y%m%d')
        outname = 'Gauss_%s.%s.fits' % (ext, today)
    elif not outname.endswith('.fits'):
        outname = outname + '.fits'

    results.write(outname, format='fits', overwrite=True)
    print('Wrote results to %s' % outname)
    return results


def do_all(filename, ext='FLUX', outname=''):
    '''
    Fit airglow lines for every exposure row in a SummarizeCframe FITS file.

    Parameters
    ----------
    filename : str
        Path to a SummarizeCframe output FITS file.
    ext : str
        Name of the flux extension to fit (FLUX, SKY_EAST, or SKY_WEST).
    outname : str
        Output filename.  If empty, auto-generated from the input stem.

    Returns
    -------
    results : astropy Table or None
    '''
    try:
        x = fits.open(filename)
    except Exception:
        print('Error: Could not open %s' % filename)
        return

    wave = x['WAVE'].data
    try:
        flux = x[ext].data
    except KeyError:
        print('Error: extension %s not found in %s' % (ext, filename))
        return

    drp   = Table(x['drp_all'].data)
    n_exp = flux.shape[0]

    records = []
    for i in range(n_exp):
        one_spec = Table([wave, flux[i]], names=['WAVE', 'FLUX'])
        if check_for_nan(flux[i]):
            print('Too many NaNs for row %d (expnum %d), skipping'
                  % (i, drp['expnum'][i]))
            continue
        rtab = do_one(one_spec)
        if len(rtab) > 0:
            rtab['expnum'] = int(drp['expnum'][i])
            rtab['mjd']    = Time(drp['obstime'][i], format='isot', scale='utc').mjd
            records.append(rtab)
        else:
            print('Nothing fit for row %d (expnum %d)' % (i, drp['expnum'][i]))

    if len(records) == 0:
        print('No spectra successfully fit in %s' % filename)
        return

    results = vstack(records)

    if outname == '':
        stem    = os.path.basename(filename).replace('.fits', '')
        outname = '%s_gauss.%s.fits' % (stem, ext)
    elif not outname.endswith('.fits'):
        outname = outname + '.fits'

    results.write(outname, format='fits', overwrite=True)
    print('Wrote results to %s' % outname)
    return results


def steer(argv):
    outname    = ''
    ext        = 'FLUX'
    drp_ver    = '1.2.1'
    input_list = ''
    cfiles     = []   # individual CFrame files (from -file or exp range)
    xfiles     = []   # SummarizeCframe files (direct .fits args)
    exp_start  = -1
    exp_stop   = -1
    delta      = 1

    i = 1
    while i < len(argv):
        if argv[i][:2] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-ext':
            i += 1
            ext = argv[i]
        elif argv[i][0:4] == '-out':
            i += 1
            outname = argv[i]
        elif argv[i] == '-ver':
            i += 1
            drp_ver = argv[i]
        elif argv[i] == '-file':
            i += 1
            input_list = argv[i]
        elif argv[i][0] == '-':
            print('Unknown option: %s' % argv[i])
            return
        elif argv[i].endswith('.fits'):
            xfiles.append(argv[i])
        else:
            try:
                val = int(argv[i])
            except ValueError:
                print('Unrecognised argument: %s' % argv[i])
                return
            if exp_start < 0:
                exp_start = val
            elif exp_stop < 0:
                exp_stop = val
            else:
                delta = val
        i += 1

    # Build CFrame file list from drpall when exp_start/exp_stop are given
    if exp_start >= 0 and exp_stop >= 0:
        try:
            from SummarizeCframe import find_top, read_drpall, select
        except ImportError:
            print('Error: could not import SummarizeCframe')
            return
        topdir  = find_top()
        drp_tab = read_drpall(drp_ver)
        if drp_tab is None or len(drp_tab) == 0:
            print('Error: could not read drpall for version %s' % drp_ver)
            return
        ztab = select(drp_tab, exp_start, exp_stop, delta)
        print('Selected %d exposures from drpall' % len(ztab))
        for row in ztab:
            loc = str(row['location']).replace('SFrame', 'CFrame')
            cfiles.append(os.path.join(topdir, loc))

    # Add files from -file table
    if input_list != '':
        try:
            xtab = Table.read(input_list)
        except Exception:
            print('Could not read: %s' % input_list)
            return
        if 'filename' in xtab.colnames:
            cfiles += list(xtab['filename'])
        elif 'Filename' in xtab.colnames:
            cfiles += list(xtab['Filename'])
        else:
            print('Read %s but could not find a column named filename or Filename'
                  % input_list)
            return

    if not cfiles and not xfiles:
        print('Nothing to do: provide filename(s), -file xfile, or exp_start exp_stop.')
        return

    if cfiles:
        do_all_files(cfiles, ext, outname, drp_ver)

    for one_file in xfiles:
        do_all(one_file, ext, outname)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
