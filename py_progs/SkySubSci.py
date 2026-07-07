#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Estimate a sky spectrum for one or more LVM CFrame exposures directly
    from the science IFU, without using the dedicated sky telescopes.

Command line usage (if any):

    usage: SkySubSci.py [-h] [-low PCT] [-high PCT] [-navg N] [-sigma S]
                        [-maxiters K] [-mask FILE] [-stat median|mean]
                        [-out ROOT] filename [filename ...]

    Arguments::

        filename    one or more lvmCFrame FITS files (one row is written per
                    file, in the order given)

    Options::

        -low PCT    percentile rank (0-100) of the faint/sky-like fiber
                    (default 10)
        -high PCT   percentile rank (0-100) of the bright/science-like fiber
                    (default 90)
        -navg N     number of fibers, ranked closest to -low/-high, to combine
                    with a sigma-clipped (robust) mean (default 10; use 1 to
                    reproduce the original single-fiber behaviour)
        -sigma S    sigma-clipping threshold for the robust mean (default 3.0)
        -maxiters K sigma-clipping iteration limit (default 5)
        -mask FILE  palace_mask FITS file from palace_make_mask.py (WAVE/MASK
                    extensions, MASK=1 means clean/sky-line-free).  If omitted
                    the script searches for sky_mask.fits in the current
                    directory, then in the lvm_ksl data/ directory.
        -stat STAT  statistic used to rank fibers by continuum flux: median
                    (default) or mean
        -out ROOT   output filename root; default is
                    SkySubSci_<first_expnum>_<last_expnum> (or
                    SkySubSci_<expnum> for a single file)

Description:

    At a typical LVM science pointing, most of the 1801 science fibers on
    the Sci telescope see mostly sky rather than an astronomical source.
    For each exposure this script:

    1. Selects science fibers from SLITMAP (scifib, imported from
       SummarizeCframe.py).
    2. Measures each fiber's continuum flux as the median (or mean) FLUX
       over sky-line-free pixels, using the mask from sky_mask.fits
       (palace_make_mask.py output, resampled onto this file's own
       wavelength grid via GetSkyCont._interp_mask_to_wave).
    3. Sorts fibers by that continuum level and, around each of the -low
       and -high percentile ranks (default 10th / 90th), takes a window of
       the -navg fibers whose rank is closest to that target and combines
       them pixel-by-pixel with a sigma-clipped mean (astropy.stats.
       sigma_clipped_stats).  With -navg 1 this reduces to picking the
       single nearest-rank fiber's raw spectrum, the original behaviour.
    4. The -low window is treated as the sky estimate; the -high window is
       treated as a bright/science-like reference; their difference is the
       sky-subtracted result.

    Averaging several fibers per percentile trades spatial locality in
    rank-space for lower per-pixel noise (~1/sqrt(navg) for clean pixels);
    the sigma clipping keeps a fiber with a faint source or a defect from
    dominating the average.

    Multiple exposures are stacked into one output file, one row per
    exposure, in the order given on the command line.

    This prototype originates from
    /Users/long/Projects/lvm_sky2607/MakeSky/Today260702.ipynb.

    Output FITS structure, extension names chosen to be readable directly
    by SkySub_eval.py (which expects FLUX = sky-subtracted spectrum and
    SKY = sky model, so that FLUX + SKY reconstructs the original)::

        PRIMARY   header with run parameters (LOWPCT, HIGHPCT, NAVG,
                  SIGCLIP, MAXITER, MASKFILE, STAT, N_PROC)
        WAVE      float32 (Npix,)          wavelength grid (from the first
                                            input file)
        SCI       float32 (Nobs, Npix)     robust mean of the high-
                                            percentile fiber window
                                            (reference only; not read by
                                            SkySub_eval)
        SKY       float32 (Nobs, Npix)     robust mean of the low-
                                            percentile fiber window (sky
                                            model)
        FLUX      float32 (Nobs, Npix)     SCI - SKY (sky-subtracted)
        DRP_ALL   BinTable (Nobs rows)     one row per exposure:
                      filename, expnum, exptime, obstime, mjd,
                      n_sci_fibers, n_avg_sci, n_avg_sky,
                      fiberid_sci_list, fiberid_sky_list,
                      ra_sci, dec_sci, ra_sky, dec_sky,
                      spectrographid_sci, spectrographid_sky,
                      contflux_sci, contflux_sky

Primary routines:

    pick_sky_sci     process one CFrame file
    process_files    batch-process a list of CFrame files and write output

Notes:

    Files that cannot be opened, or that have too few usable science
    fibers, are skipped with a warning; the run continues with the
    remaining files.

History::

    260702  ksl  Coding begun
    260702  ksl  Added -navg sigma-clipped robust mean over a window of
                 fibers around each percentile rank, instead of a single
                 nearest-rank fiber
    260706  ksl  DRP_ALL['mjd'] now computed precisely from 'obstime' via
                 SkySubOrig.obstime_to_mjd(), instead of the truncated
                 integer header keyword MJD.

'''

import sys
import os
import re
import warnings

# ensure py_progs siblings are importable when running directly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from astropy.utils.exceptions import AstropyWarning

from SummarizeCframe import scifib
from GetSkyCont import load_mask, _interp_mask_to_wave
from SkySubOrig import obstime_to_mjd

_USAGE = '''Usage:
  SkySubSci.py [-low PCT] [-high PCT] [-navg N] [-sigma S] [-maxiters K]
              [-mask FILE] [-stat median|mean] [-out ROOT]
              filename [filename ...]

Arguments:
  filename    one or more lvmCFrame FITS files (one output row per file)

Options:
  -low PCT    percentile rank of the sky-like fiber (default 10)
  -high PCT   percentile rank of the science-like fiber (default 90)
  -navg N     number of nearest-rank fibers to combine per percentile,
              via a sigma-clipped mean (default 10; use 1 for the
              original single-fiber behaviour)
  -sigma S    sigma-clipping threshold for the robust mean (default 3.0)
  -maxiters K sigma-clipping iteration limit (default 5)
  -mask FILE  palace_mask FITS file (default: sky_mask.fits searched in
              cwd then data/)
  -stat STAT  median (default) or mean, ranking statistic
  -out ROOT   output filename root (default: SkySubSci_<first>_<last>)
'''

_EXPNUM_RE = re.compile(r'(\d+)')


def _expnum_from_filename(filename):
    '''Best-effort exposure number parsed from an lvmCFrame filename.'''
    stem = os.path.basename(filename)
    m = _EXPNUM_RE.findall(stem)
    return int(m[-1]) if m else -1


def _rank_window(order, i_target, navg):
    '''Return up to navg entries of order centred on rank i_target.

    The window is shifted inward at the ends of the array so it still has
    navg entries where possible, rather than being truncated.
    '''
    n = len(order)
    navg = max(1, min(navg, n))
    lo = i_target - navg // 2
    lo = max(0, min(lo, n - navg))
    hi = lo + navg
    return order[lo:hi]


def _robust_mean(flux_window, sigma=3.0, maxiters=5):
    '''Per-pixel sigma-clipped mean across a window of fiber spectra.

    flux_window : ndarray, shape (n_fib_in_window, n_pix)

    Returns a 1-D array of length n_pix.  Pixels where every fiber in the
    window is NaN (e.g. a detector column masked bad for all fibers) come
    back as NaN; numpy's "empty slice"/"all-NaN slice" RuntimeWarnings for
    those columns are expected and suppressed here, along with the
    AstropyWarning sigma_clipped_stats raises for the same reason.
    '''
    if flux_window.shape[0] == 1:
        return flux_window[0].copy()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', AstropyWarning)
        warnings.simplefilter('ignore', RuntimeWarning)
        mean, _, _ = sigma_clipped_stats(flux_window, sigma=sigma,
                                         maxiters=maxiters, axis=0)
    return np.asarray(mean)


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
        around each of the -low/-high targets.  navg=1 picks the single
        nearest-rank fiber's raw spectrum.
    sigma, maxiters : float, int
        Sigma-clipping parameters for the robust mean (see
        astropy.stats.sigma_clipped_stats).
    mask_wave, mask_bool : ndarray
        Sky-line mask (from GetSkyCont.load_mask), resampled onto this
        file's wave grid internally.
    stat : str
        'median' or 'mean' -- statistic used to rank fibers.

    Returns
    -------
    dict with keys wave, sci_flux, sky_flux, meta (dict of scalar metadata),
    or None if the file could not be processed.
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

        hdr = x[0].header
        obstime_str = str(hdr.get('OBSTIME', ''))
        meta = dict(
            filename           = os.path.basename(filename),
            expnum             = _expnum_from_filename(filename),
            exptime            = float(hdr.get('EXPTIME', np.nan)),
            obstime            = obstime_str,
            mjd                = float(obstime_to_mjd(obstime_str)) if obstime_str else -1.0,
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


def process_files(filenames, low=10, high=90, navg=10, sigma=3.0, maxiters=5,
                  mask_file='', stat='median', outroot=''):
    '''Process a list of CFrame files and write a combined SkySubSci FITS file.

    Parameters
    ----------
    filenames : list of str
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
        Output filename root; default derived from exposure numbers.
    '''
    mask_wave, mask_bool = load_mask(mask_file)

    ref_wave = None
    sci_list, sky_list, meta_list = [], [], []

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

        if (i + 1) % 10 == 0 or (i + 1) == len(filenames):
            print('Finished %d of %d' % (i + 1, len(filenames)))

    if not sci_list:
        print('Error: no files were successfully processed; nothing written.')
        return

    sci_arr  = np.array(sci_list, dtype=np.float32)
    sky_arr  = np.array(sky_list, dtype=np.float32)
    flux_arr = sci_arr - sky_arr

    drp_all = Table(rows=[tuple(m.values()) for m in meta_list],
                    names=list(meta_list[0].keys()))

    if outroot == '':
        expnums = [m['expnum'] for m in meta_list]
        if len(expnums) == 1:
            outroot = 'SkySubSci_%d' % expnums[0]
        else:
            outroot = 'SkySubSci_%d_%d' % (expnums[0], expnums[-1])
    if not outroot.endswith('.fits'):
        outfile = outroot + '.fits'
    else:
        outfile = outroot

    hdr = fits.Header()
    hdr['LOWPCT']   = low
    hdr['HIGHPCT']  = high
    hdr['NAVG']     = navg
    hdr['SIGCLIP']  = sigma
    hdr['MAXITER']  = maxiters
    hdr['MASKFILE'] = os.path.basename(mask_file)
    hdr['STAT']     = stat
    hdr['N_PROC']   = len(meta_list)

    hdul = fits.HDUList([
        fits.PrimaryHDU(header=hdr),
        fits.ImageHDU(data=ref_wave.astype(np.float32), name='WAVE'),
        fits.ImageHDU(data=sci_arr,  name='SCI'),
        fits.ImageHDU(data=sky_arr,  name='SKY'),
        fits.ImageHDU(data=flux_arr, name='FLUX'),
        fits.BinTableHDU(drp_all, name='DRP_ALL'),
    ])
    hdul.writeto(outfile, overwrite=True)
    print('Wrote results to %s' % outfile)


def steer(argv):
    filenames = []
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
        else:
            filenames.append(arg)
        i += 1

    if not filenames:
        print(_USAGE)
        return

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

    process_files(filenames, low=low, high=high, navg=navg, sigma=sigma,
                 maxiters=maxiters, mask_file=mask_file, stat=stat, outroot=outroot)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(_USAGE)
