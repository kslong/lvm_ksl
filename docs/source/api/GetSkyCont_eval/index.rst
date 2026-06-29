GetSkyCont_eval
===============

.. py:module:: GetSkyCont_eval

.. autoapi-nested-parse::

                  Space Telescope Science Institute

   Synopsis:

       Evaluate continuum fits produced by GetSkyCont.py by plotting flux,
       continuum, components, and residuals over a chosen wavelength window.
       Reads a FITS file produced by GetSkyCont.py, writes an interactive
       Plotly HTML file, and updates the DRP_ALL table in the input FITS
       file with per-spectrum fit-quality statistics.

   Command line usage (if any):

       usage: GetSkyCont_eval.py skycont_file.fits [wmin wmax] [-num N] [-out outroot]

       Arguments:

       skycont_file.fits  FITS file from GetSkyCont.py (process_skyfile or process_many).
       wmin               minimum wavelength in Angstroms (default 3600).
       wmax               maximum wavelength in Angstroms (default 9800).

       Options:

       -num N             overlay N randomly-selected individual spectra on the band (default 20; 0 = band only).
       -out outroot       output filename root; default is <stem>_<wmin>_<wmax>.

   Description:

       Writes a single HTML file containing three interactive Plotly figures
       and updates the DRP_ALL table in the input FITS file with per-spectrum
       fit-quality statistics.

       Figure 1 — four-panel spectral overview:

       Panel 1 (Flux + Continuum, log): observed FLUX with the fitted CONT median
       overlaid in red.  Blue vertical bands mark the clean (fit) regions;
       grey bands mark the masked (line-affected, excluded) regions.

       Panel 2 (Total Continuum, log): the total CONT band, showing the overall
       smoothness and level of the B-spline model.

       Panel 3 (Components, log): MOON (red; B-spline x solar) and DIFFUSE
       (orange; plain B-splines) plotted together with the total CONT median
       (purple) for reference.  The relative amplitude of the components shows
       how much of the continuum is Moon/zodiacal vs. diffuse airglow.

       Panel 4 (Residual, linear): RESID = FLUX - CONT; should be near zero in
       clean regions and show sky line emission in the masked regions.

       Both log panels share a floor at 1e-18 to suppress near-zero noise.

       Figure 2 — per-arm residual histograms:

       Three panels (Blue 3600-5900 A, Red 5900-7600 A, NIR 7600-9800 A) each
       showing the distribution of all residual flux values (all spectra, all
       clean pixels) in that wavelength range.  The x-axis is residual flux;
       the y-axis is count N.  A Gaussian with center = median and sigma = NMAD
       is overlaid in black.  An annotation box reports N, the median, NMAD,
       skewness, and the 10th/90th percentiles.  The histogram range is clipped
       to median +/- 5*NMAD to suppress extreme outliers.

       Figure 3 — per-spectrum fit quality:

       Row 1: three scatter plots (Blue, Red, NIR) of per-spectrum median
       residual (x) vs NMAD (y).  Each point is one spectrum; hover shows the
       spectrum row index and its median/NMAD.  Axis limits are clipped to the
       2nd-98th percentile range so that extreme outliers do not compress the
       scale; the annotation reports how many spectra are off-scale.  A dashed
       vertical line marks x=0 (ideal median); a dotted horizontal line marks
       the ensemble NMAD for reference.

       Row 2: NMAD for all three arms plotted against original spectrum number
       (log y-scale) so that clusters of temporally adjacent bad observations
       are immediately visible.  Hover shows spectrum index and NMAD value.
       A dotted horizontal line marks the ensemble NMAD per arm.

       Per-spectrum statistics written to DRP_ALL:

       For each arm (blue, red, nir) and for all arms combined (all), four
       float32 columns are added (or updated) in the DRP_ALL table:
           resid_med_<arm>   — median residual in clean pixels
           resid_nmad_<arm>  — NMAD of residuals in clean pixels
           resid_rms_<arm>   — RMS of residuals in clean pixels
           resid_skew_<arm>  — skewness of residuals in clean pixels

   Primary routines:

       plot_eval   build and write the Plotly HTML evaluation plot.

   .. rubric:: Notes

   Requires plotly (pip install plotly).
   Wavelength window defaults to the full LVM range (3600-9800 A) when
   wmin and wmax are not given.
   The input FITS file is updated in place; a backup is not created.

   History:

       260628  ksl  Written.
       260628  ksl  Redesigned to four panels: Flux, Total Cont, Components, Residual.
       260629  ksl  Added per-arm residual histogram figure (Figure 2).
       260629  ksl  Added per-spectrum fit-quality figure (Figure 3) and DRP_ALL update.



Functions
---------

.. autoapisummary::

   GetSkyCont_eval.plot_eval
   GetSkyCont_eval.steer


Module Contents
---------------

.. py:function:: plot_eval(filename, wmin=3600.0, wmax=9800.0, n_sample=20, outroot='')

   Build spectral evaluation plots and per-spectrum quality diagnostics
   from a GetSkyCont output FITS file.

   Writes a three-figure HTML file and updates the DRP_ALL table in the
   input FITS file with per-spectrum fit-quality statistics (med, nmad, rms,
   skew) for each spectrograph arm and for all arms combined.

   :param filename: skycont_*.fits file produced by GetSkyCont.py.
   :type filename: str or Path
   :param wmin: Minimum wavelength in Angstroms.
   :type wmin: float
   :param wmax: Maximum wavelength in Angstroms.
   :type wmax: float
   :param n_sample: Number of randomly-selected individual spectra to overlay.  0 = band only.
   :type n_sample: int
   :param outroot: Output filename root.  Default: <stem>_<wmin>_<wmax>.
   :type outroot: str


.. py:function:: steer(argv)

