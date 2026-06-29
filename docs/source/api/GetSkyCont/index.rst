GetSkyCont
==========

.. py:module:: GetSkyCont

.. autoapi-nested-parse::

                       Space Telescope Science Institute

   Synopsis:

       Fit a smooth continuum to LVM sky spectra using a two-component B-spline
       design matrix, restricting the fit to wavelength regions unaffected by sky
       emission lines.  Operates independently of the lvmsky and lvmdrp packages.

   Command line usage (if any):

       # Sky-file mode (auto-detected from file content):
       usage: GetSkyCont.py sky_file.fits -mask mask.fits [row_no ...] [-delta N]

       # XCframe / XSFrame summary file mode:
       usage: GetSkyCont.py xframe.fits ext -mask mask.fits [row_no ...] [-delta N]

       Arguments:

       sky_file.fits  Sky_<name>.fits from GetSky_from_CFrame_sum.py; all spectra processed if no row numbers are given.
       xframe.fits    an XCframe or XSFrame summary FITS file.
       ext            FITS extension for the spectrum (FLUX, SKY_EAST, SKY_WEST, ...).
       row_no         zero or more 0-based row indices; if omitted all rows are processed (subject to -delta).

       Options:

       -mask file     (required) palace_mask FITS file from palace_make_mask.py (MASK extension: 1=clean 0=line-affected).
       -delta N       process every N-th row (0, N, 2N, ...) instead of all; ignored when explicit row numbers are given.
       -kstep N       B-spline knot spacing in Angstroms (default 100).
       -out outroot   set the output filename root.

   Description:

       Builds a two-component cubic B-spline design matrix spanning the full LVM
       wavelength range with knots every kstep Angstroms (default 100 A, giving
       ~66 basis functions per component):

       DIFFUSE component: plain B-spline basis functions.  These absorb smooth
       diffuse airglow emission (HO2, FeO, scattered light) whose spectral shape
       is not described by the solar spectrum.

       MOON component: the same B-spline basis functions each multiplied
       pointwise by a rebinned solar spectrum.  The Fraunhofer absorption
       structure is fixed; the B-spline envelope adjusts the colour and
       amplitude to fit the Moon-reflected and zodiacal-light continuum.

       The design matrix is evaluated only at pixels flagged as clean by the
       palace_make_mask.py mask (MASK=1).  Non-negative least squares (scipy
       nnls) enforces positive component amplitudes.  The fitted components are
       then evaluated at all wavelengths so both interpolate smoothly across
       masked line regions.

       The residual FLUX - CONT = FLUX - MOON - DIFFUSE isolates line emission
       free of continuum bias and is saved for subsequent line fitting.

   Primary routines:

       build_design_matrix  construct the two-component B-spline design matrix.
       fit_continuum        fit and decompose continuum for one spectrum.
       process_skyfile      batch process a Sky_<name>.fits file.
       process_many         batch process an XCframe/XSFrame file.

   .. rubric:: Notes

   The solar spectrum is read from
   /Users/long/Projects/lvm_sky2606/skysub_ivan/Spectre_HR_LATMOS_Meftah_V1_350_1000nm.txt
   (wavelengths in nm, converted to Angstroms internally).  If the file is
   absent the MOON component is omitted and CONT = DIFFUSE only; a warning
   is printed.

   Output FITS extensions: PRIMARY, WAVE, FLUX, CONT, MOON, DIFFUSE, RESID, MASK, DRP_ALL.
   Default output name: skycont_<stem>.fits or skycont_<ext>_<stem>.fits.

   History:

       260628  ksl  Written; solar component optional.
       260628  ksl  Solar made default; MOON and DIFFUSE always written as separate extensions.



Attributes
----------

.. autoapisummary::

   GetSkyCont.DEFAULT_SOLAR_FILE


Functions
---------

.. autoapisummary::

   GetSkyCont.build_design_matrix
   GetSkyCont.fit_continuum
   GetSkyCont.load_mask
   GetSkyCont.load_solar
   GetSkyCont.process_many
   GetSkyCont.process_skyfile
   GetSkyCont.steer


Module Contents
---------------

.. py:data:: DEFAULT_SOLAR_FILE

.. py:function:: build_design_matrix(wave, knot_step=100.0, solar=None)

   Build a two-component cubic B-spline design matrix.

   :param wave: Wavelength array at which to evaluate the basis functions.
   :type wave: ndarray
   :param knot_step: Knot spacing in Angstroms.
   :type knot_step: float
   :param solar: Solar spectrum (same length as wave).  If given, a second block of
                 columns equal to solar * B-spline_basis is appended for the MOON
                 component.  If None only the DIFFUSE block is built.
   :type solar: ndarray or None

   :returns: * **A** (*ndarray, shape (len(wave), n_cols)*) -- Full design matrix.  n_cols = n_basis (DIFFUSE only) or 2*n_basis (with MOON).
             * **n_basis** (*int*) -- Number of columns in each component block.


.. py:function:: fit_continuum(flux, clean, A, n_b)

   Fit and decompose continuum for a single sky spectrum.

   :param flux: Observed sky spectrum (N_pix,).
   :type flux: ndarray
   :param clean: True at pixels included in the fit.
   :type clean: ndarray of bool
   :param A: Pre-built design matrix from build_design_matrix.
   :type A: ndarray, shape (N_pix, n_cols)
   :param n_b: Number of plain B-spline columns (first block = DIFFUSE).
   :type n_b: int

   :returns: * **cont** (*ndarray*) -- Total fitted continuum (DIFFUSE + MOON) at all wavelengths.
             * **diffuse** (*ndarray*) -- Smooth diffuse component (plain B-splines).
             * **moon** (*ndarray*) -- Moon/zodi component (solar x B-splines).  Zero array when A has only one block.
             * **coef** (*ndarray*) -- Non-negative NNLS coefficients.
             * **n_clean** (*int*) -- Number of clean pixels used in the fit.


.. py:function:: load_mask(mask_file)

   Load a palace_make_mask.py output and return (wave, mask).

   :param mask_file: FITS file with WAVE and MASK extensions.
   :type mask_file: str or Path

   :returns: * **wave** (*ndarray*) -- Wavelength array (Angstroms).
             * **mask** (*ndarray of bool*) -- True where the pixel is clean (unaffected by sky lines).


.. py:function:: load_solar(wave_out)

   Load and resample the solar spectrum onto wave_out (Angstroms).

   :param wave_out: Target wavelength grid in Angstroms.
   :type wave_out: ndarray

   :returns: **solar** -- Solar flux resampled to wave_out, normalised to unit median.
             Returns None if the solar file is not found.
   :rtype: ndarray or None


.. py:function:: process_many(filename, ext, rows=None, delta=None, mask_wave=None, mask_arr=None, knot_step=100.0, outroot='')

   Fit continuum to rows of an XCframe/XSFrame file.

   :param filename:
   :type filename: str or Path
   :param ext: FITS extension name (e.g. SKY_EAST).
   :type ext: str
   :param rows:
   :type rows: list of int or None
   :param delta:
   :type delta: int or None
   :param mask_wave:
   :type mask_wave: ndarray or None
   :param mask_arr:
   :type mask_arr: ndarray of bool or None
   :param knot_step:
   :type knot_step: float
   :param outroot:
   :type outroot: str


.. py:function:: process_skyfile(filename, rows=None, delta=None, mask_wave=None, mask_arr=None, knot_step=100.0, outroot='')

   Fit continuum to spectra in a Sky_<name>.fits file.

   :param filename: Sky file from GetSky_from_CFrame_sum.py.
   :type filename: str or Path
   :param rows: Explicit row indices; None means all (or every delta-th).
   :type rows: list of int or None
   :param delta: Row stride; ignored if rows is not None.
   :type delta: int or None
   :param mask_wave:
   :type mask_wave: ndarray or None
   :param mask_arr: Boolean mask (True=clean); None means all pixels.
   :type mask_arr: ndarray of bool or None
   :param knot_step:
   :type knot_step: float
   :param outroot:
   :type outroot: str


.. py:function:: steer(argv)

