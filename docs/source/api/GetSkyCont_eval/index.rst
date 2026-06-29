GetSkyCont_eval
===============

.. py:module:: GetSkyCont_eval

.. autoapi-nested-parse::

                       Space Telescope Science Institute

   Synopsis:

       Evaluate continuum fits produced by GetSkyCont.py by plotting flux,
       continuum, components, and residuals over a chosen wavelength window.
       Reads a FITS file produced by GetSkyCont.py and writes an interactive
       Plotly HTML file.

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

       Creates a four-panel interactive HTML plot from a skycont_*.fits file:

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

   Primary routines:

       plot_eval   build and write the Plotly HTML evaluation plot.

   .. rubric:: Notes

   Requires plotly (pip install plotly).
   Wavelength window defaults to the full LVM range (3600-9800 A) when
   wmin and wmax are not given.

   History:

       260628  ksl  Written.
       260628  ksl  Redesigned to four panels: Flux, Total Cont, Components, Residual.



Functions
---------

.. autoapisummary::

   GetSkyCont_eval.plot_eval
   GetSkyCont_eval.steer


Module Contents
---------------

.. py:function:: plot_eval(filename, wmin=3600.0, wmax=9800.0, n_sample=20, outroot='')

   Build a four-panel Plotly evaluation plot from a GetSkyCont output FITS file.

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

