#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Drpall-driven, remote-friendly version of SkySubSci.py: estimate a
    sky spectrum for a range of LVM exposures directly from the science
    IFU, selecting exposures the same way SummarizeCframe.py does
    (exposure-number range against a drpall table) instead of taking an
    explicit list of filenames.

Command line usage (if any):

    usage: SummarizeSciSky.py [-h] [-emin 900] [-ver 1.2.1] [-drp_all FILE]
                              [-low 10] [-high 90] [-navg 10] [-sigma 3.0]
                              [-maxiters 5] [-mask FILE] [-stat median|mean]
                              [-out ROOT] exp_start exp_stop [delta]

    Arguments::

        exp_start   starting exposure number
        exp_stop    stopping exposure number
        delta       process every delta-th exposure in range (default 1)

    Options::

        -emin N       minimum exposure time to include (default 900)
        -ver VER      DRP version, used to locate drpall-VER.fits (default 1.2.1)
        -drp_all FILE explicit drpall table to read instead of drpall-VER.fits
                      (FITS, or ascii if the name contains "txt"/".tab")
        -low PCT      percentile rank (0-100) of the faint/sky-like fiber
                      window (default 10)
        -high PCT     percentile rank (0-100) of the bright/science-like
                      fiber window (default 90)
        -navg N       number of nearest-rank fibers to combine per percentile,
                      via a sigma-clipped mean (default 10; use 1 for a
                      single nearest-rank fiber)
        -sigma S      sigma-clipping threshold for the robust mean (default 3.0)
        -maxiters K   sigma-clipping iteration limit (default 5)
        -mask FILE    palace_mask FITS file (default: sky_mask.fits searched
                      in cwd then the lvm_ksl data/ directory)
        -stat STAT    median (default) or mean, ranking statistic
        -out ROOT     output filename root; default is
                      SummarizeSciSky_<ver>_<exp_start>_<exp_stop>_<delta>

Description:

    This script exists so the science-fiber-based sky estimate in
    SkySubSci.py can be run unattended over many exposures (e.g. at
    Utah) without having to hand it a file list. Exposure selection and
    file-path resolution follow the same pattern as SummarizeCframe.py's
    read_drpall/select/find_top, re-implemented locally here (with an
    added -drp_all override) rather than imported from SumCframe.py or
    SummarizeCframe.py, so this script has no dependency on the optional
    "dask" package that SumCframe.py requires for an unrelated function.
    read_drpall/select_exps/find_top pick out rows in [exp_start,
    exp_stop] with exptime >= -emin, every delta-th one, and resolve
    each row's "location" column to an actual file under find_top()'s
    data directory (Utah / Rainbow / Muskie), renaming SFrame -> CFrame.

    For each resolved exposure, the same algorithm as SkySubSci.py is
    applied (re-implemented locally here -- this script is intentionally
    standalone and does not import from SkySubSci.py):

    1. Select science fibers from SLITMAP (scifib, imported from
       SummarizeCframe.py).
    2. Measure each fiber's continuum flux as the median (or mean) FLUX
       over sky-line-free pixels, using the mask from sky_mask.fits
       (resampled onto this file's wavelength grid via
       GetSkyCont._interp_mask_to_wave).
    3. Sort fibers by that continuum level and, around each of the -low
       and -high percentile ranks, take a window of the -navg fibers
       whose rank is closest to that target and combine them
       pixel-by-pixel with a sigma-clipped mean.
    4. The -low window is the sky estimate; the -high window is a
       bright/science-like reference; their difference is the
       sky-subtracted result.

    Rather than building a fresh metadata table, the calculated values
    (continuum flux, fiber IDs used, positions, spectrograph IDs) are
    added as new columns directly onto the selected drpall rows, which
    become the DRP_ALL extension of the output.

    Output FITS structure, extension names chosen to be readable
    directly by SkySub_eval.py (FLUX = sky-subtracted spectrum, SKY =
    sky model, so FLUX + SKY reconstructs the original)::

        PRIMARY   header with run parameters (LOWPCT, HIGHPCT, NAVG,
                  SIGCLIP, MAXITER, MASKFILE, STAT, DRPVER, EMIN,
                  EXPSTART, EXPSTOP, DELTA, N_PROC)
        WAVE      float32 (Npix,)          wavelength grid (from the
                                            first successfully processed
                                            exposure)
        SCI       float32 (Nobs, Npix)     robust mean of the high-
                                            percentile fiber window
                                            (reference only; not read by
                                            SkySub_eval)
        SKY       float32 (Nobs, Npix)     robust mean of the low-
                                            percentile fiber window (sky
                                            model)
        FLUX      float32 (Nobs, Npix)     SCI - SKY (sky-subtracted)
        DRP_ALL   BinTable (Nobs rows)     the selected drpall rows, with
                      added columns: n_sci_fibers, n_avg_sci, n_avg_sky,
                      fiberid_sci_list, fiberid_sky_list, ra_sci,
                      dec_sci, ra_sky, dec_sky, spectrographid_sci,
                      spectrographid_sky, contflux_sci, contflux_sky

Primary routines:

    pick_sky_sci     process one CFrame file
    process_drpall    select exposures from a drpall table and write output

Notes:

    Exposures that cannot be opened, or that have too few usable science
    fibers, are skipped with a warning; the run continues with the
    remaining exposures.

    This script deliberately duplicates a small amount of logic from
    SkySubSci.py (rank-window selection and the sigma-clipped robust
    mean) and from SummarizeCframe.py/SumCframe.py (drpall selection and
    file-path resolution) instead of importing it, so it can run
    standalone with no dependency beyond astropy/numpy/scipy -- in
    particular it avoids SumCframe.py's "dask" import, which is only
    needed there for a different (unused) function.

History::

    260703  ksl  Coding begun
    260706  ksl  DRP_ALL['mjd'] now recomputed precisely from 'obstime' via
                 SkySubOrig.obstime_to_mjd(), instead of the truncated
                 integer carried through from the master drpall table.

'''

import sys
import os

# ensure py_progs siblings are importable when running directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from astropy.io import fits
from astropy.table import Table

from astropy.io import ascii as apy_ascii
from SummarizeCframe import scifib, _rank_window, _robust_mean
from GetSkyCont import load_mask, _interp_mask_to_wave
from SkySubOrig import obstime_to_mjd

_USAGE = '''Usage:
  SummarizeSciSky.py [-emin N] [-ver VER] [-drp_all FILE] [-low PCT]
                     [-high PCT] [-navg N] [-sigma S] [-maxiters K]
                     [-mask FILE] [-stat median|mean] [-out ROOT]
                     exp_start exp_stop [delta]

Arguments:
  exp_start   starting exposure number
  exp_stop    stopping exposure number
  delta       process every delta-th exposure in range (default 1)

Options:
  -emin N       minimum exposure time to include (default 900)
  -ver VER      DRP version, used to locate drpall-VER.fits (default 1.2.1)
  -drp_all FILE explicit drpall table (FITS or ascii) instead of
                drpall-VER.fits
  -low PCT      percentile rank of the sky-like fiber window (default 10)
  -high PCT     percentile rank of the science-like fiber window (default 90)
  -navg N       number of nearest-rank fibers to combine per percentile,
                via a sigma-clipped mean (default 10)
  -sigma S      sigma-clipping threshold for the robust mean (default 3.0)
  -maxiters K   sigma-clipping iteration limit (default 5)
  -mask FILE    palace_mask FITS file (default: sky_mask.fits searched in
                cwd then data/)
  -stat STAT    median (default) or mean, ranking statistic
  -out ROOT     output filename root (default:
                SummarizeSciSky_<ver>_<exp_start>_<exp_stop>_<delta>)
'''


# ──────────────────────────────────────────────────────────────
# Standalone copies of the SummarizeCframe.py / SumCframe.py drpall
# selection logic (avoids SumCframe.py's "dask" dependency, which is
# only needed there for a different, unused function)
# ──────────────────────────────────────────────────────────────

def read_drpall(filename='', drp_ver='1.2.1'):
    '''Read a drpall FITS file, or an ascii table, and return the table.

    Parameters
    ----------
    filename : str
        Explicit drpall table to read (FITS, or ascii if the name
        contains "txt" or ".tab"). If empty, drpall-<drp_ver>.fits is
        looked for locally, then under the Utah redux tree.
    drp_ver : str
        DRP version, used to build the default filename.

    Returns
    -------
    astropy.table.Table, or an empty list if the file cannot be found
    or read.
    '''
    if filename.count('txt') or filename.count('.tab'):
        try:
            return apy_ascii.read(filename)
        except Exception:
            print('Error: Could not locate : ', filename)
            return []

    if filename == '':
        DRPFILE = 'drpall-%s.fits' % (drp_ver)
    else:
        DRPFILE = filename

    if os.path.isfile(DRPFILE):
        xfile = DRPFILE
    else:
        BASEDIR = '/uufs/chpc.utah.edu/common/home/sdss51/sdsswork/lvm/spectro/redux/%s/' % (drp_ver)
        xfile = '%s/%s' % (BASEDIR, DRPFILE)
        if not os.path.isfile(xfile):
            print('Error: Could not locate : ', xfile)
            return []

    try:
        drpall = fits.open(xfile)
        print('Successfully opened ', xfile)
    except Exception:
        print('Error: Located but could not read : ', xfile)
        return []

    return Table(drpall[1].data)


def select_exps(ztab, exp_start=4000, exp_stop=8000, delta=5, exp_min=900.):
    '''Select every delta-th exposure between exp_start and exp_stop.'''
    xtab = ztab[ztab['expnum'] >= exp_start]
    xtab = xtab[xtab['expnum'] <= exp_stop]
    if exp_min > 0:
        xtab = xtab[xtab['exptime'] >= exp_min]
    if delta > 1:
        xtab = xtab[::delta]
    return xtab


_XTOP     = '/uufs/chpc.utah.edu/common/home/sdss51/'
_XRAINBOW = '/Users/long/Projects/lvm_data/sas'
_XMUSKIE  = '/home/long/Projects/lvm_data/sas'


def find_top():
    '''Locate the top of the local redux data tree (Utah / Rainbow / Muskie).'''
    if os.path.isdir(_XTOP):
        loc, topdir = 'Utah', _XTOP
    elif os.path.isdir(_XRAINBOW):
        loc, topdir = 'Rainbow', _XRAINBOW
    elif os.path.isdir(_XMUSKIE):
        loc, topdir = 'Muskie', _XMUSKIE
    else:
        print('Error: I do not know where I am:', os.getcwd())
        return ''
    print('We are on : ', loc)
    return topdir


def pick_sky_sci(filename, low=10, high=90, navg=10, sigma=3.0, maxiters=5,
                 mask_wave=None, mask_bool=None, stat='median'):
    '''Select sky-like and science-like fiber spectra from one CFrame file.

    Parameters
    ----------
    filename : str
        lvmCFrame FITS file.
    low, high : float
        Percentile ranks (0-100) of the sky-like and science-like fiber
        windows, after sorting science fibers by continuum flux.
    navg : int
        Number of nearest-rank fibers averaged (via a sigma-clipped mean)
        around each of the -low/-high targets.
    sigma, maxiters : float, int
        Sigma-clipping parameters for the robust mean.
    mask_wave, mask_bool : ndarray
        Sky-line mask (from GetSkyCont.load_mask), resampled onto this
        file's wave grid internally.
    stat : str
        'median' or 'mean' -- statistic used to rank fibers.

    Returns
    -------
    dict with keys wave, sci_flux, sky_flux, meta (dict of calculated
    scalar metadata), or None if the file could not be processed.
    '''
    try:
        x = fits.open(filename)
    except Exception as e:
        print('Warning: could not open %s (%s), skipping.' % (filename, e))
        return None

    try:
        slitmap = Table(x['SLITMAP'].data)
        sci = scifib(slitmap, 'all', 'Sci')
        if len(sci) < 10:
            print('Warning: only %d science fibers found in %s, skipping.'
                  % (len(sci), filename))
            return None

        wave = x['WAVE'].data.astype(np.float64)
        flux = x['FLUX'].data[sci['fiberid'] - 1].astype(np.float64)
        bad  = x['MASK'].data[sci['fiberid'] - 1] != 0
        flux[bad] = np.nan

        clean = _interp_mask_to_wave(mask_wave, mask_bool, wave)
        if stat == 'mean':
            cont = np.nanmean(flux[:, clean], axis=1)
        else:
            cont = np.nanmedian(flux[:, clean], axis=1)

        good = np.isfinite(cont)
        if good.sum() < 10:
            print('Warning: too few fibers with valid continuum flux in %s, skipping.'
                  % filename)
            return None
        sci_tab = sci[good]
        cont    = cont[good]
        flux    = flux[good]

        order = np.argsort(cont)
        n     = len(order)
        i_lo  = int(round(low  / 100.0 * (n - 1)))
        i_hi  = int(round(high / 100.0 * (n - 1)))
        win_sky = _rank_window(order, i_lo, navg)
        win_sci = _rank_window(order, i_hi, navg)

        sky_flux = _robust_mean(flux[win_sky], sigma=sigma, maxiters=maxiters)
        sci_flux = _robust_mean(flux[win_sci], sigma=sigma, maxiters=maxiters)

        def _mode_int(arr):
            arr = np.asarray(arr, int)
            return int(np.bincount(arr).argmax())

        meta = dict(
            n_sci_fibers       = n,
            n_avg_sci          = len(win_sci),
            n_avg_sky          = len(win_sky),
            fiberid_sci_list   = ','.join(str(v) for v in sci_tab['fiberid'][win_sci]),
            fiberid_sky_list   = ','.join(str(v) for v in sci_tab['fiberid'][win_sky]),
            ra_sci             = float(np.mean(sci_tab['ra'][win_sci])),
            dec_sci            = float(np.mean(sci_tab['dec'][win_sci])),
            ra_sky             = float(np.mean(sci_tab['ra'][win_sky])),
            dec_sky            = float(np.mean(sci_tab['dec'][win_sky])),
            spectrographid_sci = _mode_int(sci_tab['spectrographid'][win_sci]),
            spectrographid_sky = _mode_int(sci_tab['spectrographid'][win_sky]),
            contflux_sci       = float(np.mean(cont[win_sci])),
            contflux_sky       = float(np.mean(cont[win_sky])),
        )

        return dict(wave=wave, sci_flux=sci_flux, sky_flux=sky_flux, meta=meta)
    finally:
        x.close()


# ──────────────────────────────────────────────────────────────
# Drpall-driven batch processing
# ──────────────────────────────────────────────────────────────

def process_drpall(exp_start, exp_stop, delta=1, exp_min=900., drp_ver='1.2.1',
                   drp_all_file='', low=10, high=90, navg=10, sigma=3.0,
                   maxiters=5, mask_file='', stat='median', outroot=''):
    '''Select exposures from a drpall table and write a combined
    SummarizeSciSky FITS file.

    Parameters
    ----------
    exp_start, exp_stop, delta : int
        Exposure-number selection, same semantics as SummarizeCframe.py.
    exp_min : float
        Minimum exposure time to include.
    drp_ver : str
        DRP version, used to locate drpall-<drp_ver>.fits if drp_all_file
        is not given.
    drp_all_file : str
        Explicit drpall table (FITS or ascii) to use instead.
    low, high : float
        Percentile ranks of the sky-like and science-like fiber windows.
    navg : int
        Number of nearest-rank fibers combined per percentile (robust mean).
    sigma, maxiters : float, int
        Sigma-clipping parameters for the robust mean.
    mask_file : str
        sky_mask.fits path (palace_make_mask.py output).
    stat : str
        'median' or 'mean'.
    outroot : str
        Output filename root; default derived from drp_ver/exposure range.
    '''
    xtop = find_top()
    xtab = read_drpall(drp_all_file, drp_ver)
    if len(xtab) == 0:
        print('Error: could not read a drpall table; nothing to do.')
        return
    ztab = select_exps(xtab, exp_start, exp_stop, delta, exp_min)
    print('Selected %d exposures from the drpall table' % len(ztab))

    mask_wave, mask_bool = load_mask(mask_file)

    filenames = []
    keep = []
    for i in range(len(ztab)):
        xfile = '%s/%s' % (xtop, ztab['location'][i])
        if xfile.count('SFrame'):
            xfile = xfile.replace('SFrame', 'CFrame')
        if os.path.isfile(xfile):
            keep.append(i)
            filenames.append(xfile)
    print('There are %d files to process' % len(filenames))
    ztab = ztab[keep]

    ref_wave = None
    sci_list, sky_list, meta_list, good_rows = [], [], [], []

    for i, filename in enumerate(filenames):
        result = pick_sky_sci(filename, low=low, high=high, navg=navg,
                              sigma=sigma, maxiters=maxiters,
                              mask_wave=mask_wave, mask_bool=mask_bool, stat=stat)
        if result is None:
            continue

        wave     = result['wave']
        sci_flux = result['sci_flux']
        sky_flux = result['sky_flux']

        if ref_wave is None:
            ref_wave = wave
        elif len(wave) != len(ref_wave) or not np.allclose(wave, ref_wave, rtol=0, atol=1e-3):
            print('Warning: %s has a different wavelength grid; interpolating onto '
                  'the reference grid.' % filename)
            sci_flux = np.interp(ref_wave, wave, sci_flux)
            sky_flux = np.interp(ref_wave, wave, sky_flux)

        sci_list.append(sci_flux)
        sky_list.append(sky_flux)
        meta_list.append(result['meta'])
        good_rows.append(i)

        if (i + 1) % 10 == 0 or (i + 1) == len(filenames):
            print('Finished %d of %d' % (i + 1, len(filenames)))

    if not sci_list:
        print('Error: no exposures were successfully processed; nothing written.')
        return

    sci_arr  = np.array(sci_list, dtype=np.float32)
    sky_arr  = np.array(sky_list, dtype=np.float32)
    flux_arr = sci_arr - sky_arr

    ztab_final = ztab[good_rows]
    for col in meta_list[0]:
        ztab_final[col] = [m[col] for m in meta_list]
    if 'obstime' in ztab_final.colnames and 'mjd' in ztab_final.colnames:
        ztab_final['mjd'] = obstime_to_mjd(ztab_final['obstime'])

    if outroot == '':
        outroot = 'SummarizeSciSky_%s_%d_%d_%d' % (drp_ver, exp_start, exp_stop, delta)
    outfile = outroot if outroot.endswith('.fits') else outroot + '.fits'

    hdr = fits.Header()
    hdr['ROUTINE']  = ('SummarizeSciSky', 'Script that produced this file')
    hdr['LOWPCT']   = low
    hdr['HIGHPCT']  = high
    hdr['NAVG']     = navg
    hdr['SIGCLIP']  = sigma
    hdr['MAXITER']  = maxiters
    hdr['MASKFILE'] = os.path.basename(mask_file)
    hdr['STAT']     = stat
    hdr['DRPVER']   = drp_ver
    hdr['EMIN']     = exp_min
    hdr['EXPSTART'] = exp_start
    hdr['EXPSTOP']  = exp_stop
    hdr['DELTA']    = delta
    hdr['N_PROC']   = len(meta_list)

    hdul = fits.HDUList([
        fits.PrimaryHDU(header=hdr),
        fits.ImageHDU(data=ref_wave.astype(np.float32), name='WAVE'),
        fits.ImageHDU(data=sci_arr,  name='SCI'),
        fits.ImageHDU(data=sky_arr,  name='SKY'),
        fits.ImageHDU(data=flux_arr, name='FLUX'),
        fits.BinTableHDU(ztab_final, name='DRP_ALL'),
    ])
    hdul.writeto(outfile, overwrite=True)
    print('Wrote results to %s' % outfile)


def steer(argv):
    exp_start = -1
    exp_stop  = -1
    delta     = -1
    exp_min   = 900.
    drp_ver   = '1.2.1'
    drp_all_file = ''
    low       = 10.0
    high      = 90.0
    navg      = 10
    sigma     = 3.0
    maxiters  = 5
    mask_file = ''
    stat      = 'median'
    outroot   = ''

    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg in ('-h', '--help'):
            print(_USAGE)
            return
        elif arg == '-emin':
            i += 1
            exp_min = float(argv[i])
        elif arg == '-ver':
            i += 1
            drp_ver = argv[i]
        elif arg == '-drp_all':
            i += 1
            drp_all_file = argv[i]
        elif arg == '-low':
            i += 1
            low = float(argv[i])
        elif arg == '-high':
            i += 1
            high = float(argv[i])
        elif arg == '-navg':
            i += 1
            navg = int(argv[i])
        elif arg == '-sigma':
            i += 1
            sigma = float(argv[i])
        elif arg == '-maxiters':
            i += 1
            maxiters = int(argv[i])
        elif arg == '-mask':
            i += 1
            mask_file = argv[i]
        elif arg == '-stat':
            i += 1
            stat = argv[i]
        elif arg == '-out':
            i += 1
            outroot = argv[i]
        elif arg.startswith('-'):
            print('Error: unknown option "%s"' % arg)
            print(_USAGE)
            return
        elif exp_start < 0:
            exp_start = int(arg)
        elif exp_stop < 0:
            exp_stop = int(arg)
        elif delta < 0:
            delta = int(arg)
        else:
            print('Error: unexpected argument "%s"' % arg)
            print(_USAGE)
            return
        i += 1

    if exp_start < 0 or exp_stop < 0:
        print(_USAGE)
        return
    if delta < 0:
        delta = 1

    if stat not in ('median', 'mean'):
        print('Error: -stat must be median or mean')
        return

    if not mask_file:
        _data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
        for _candidate in [os.path.join(os.getcwd(), 'sky_mask.fits'),
                           os.path.join(_data_dir,   'sky_mask.fits')]:
            if os.path.exists(_candidate):
                mask_file = _candidate
                print('Using default mask: %s' % mask_file)
                break
        if not mask_file:
            print('Error: no mask file supplied and sky_mask.fits not found '
                 'in the current directory or data/')
            print(_USAGE)
            return
    elif not os.path.exists(mask_file):
        print('Error: mask file not found: %s' % mask_file)
        return

    process_drpall(exp_start, exp_stop, delta=delta, exp_min=exp_min, drp_ver=drp_ver,
                   drp_all_file=drp_all_file, low=low, high=high, navg=navg,
                   sigma=sigma, maxiters=maxiters, mask_file=mask_file,
                   stat=stat, outroot=outroot)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
