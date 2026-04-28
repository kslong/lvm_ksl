fourier_offset_check
====================

.. py:module:: fourier_offset_check

.. autoapi-nested-parse::

                       Space Telescope Science Institute

   Synopsis:

   Add a synthetic Gaussian emission line to the FLUX extension of an
   lvmCFrame file, for use in testing fourier_offset.py.

   Command line usage (if any):

       usage: check_fourier.py [-h] [-wave 4000] [-flux 1e-13]
                               [-off1 0.0] [-off2 0.1] [-off3 -0.1] filename

       where filename is an lvmCFrame FITS file.
       -wave sets the central wavelength of the line in Angstroms (default 4200).
       -flux sets the integrated line flux in erg/s/cm^2 (default 1e-11).
       -off1, -off2, -off3 set wavelength offsets in Angstroms applied to all
       fibers in spectrographs 1, 2, and 3 respectively (defaults 0, +0.1, -0.1).

   Description:

       Reads the WAVE, FLUX, LSF, and SLITMAP extensions of the input file.
       For each fiber the FWHM at the requested wavelength is taken from the
       LSF extension and converted to a Gaussian sigma.  The line centre is
       shifted by the per-spectrograph offset before the Gaussian is computed.
       The peak amplitude is derived from the integrated flux as
       A = flux / (sigma * sqrt(2*pi)), and the Gaussian is added to the FLUX
       array.  The modified file is written to test_<basename>.fits in the
       current directory.

   Primary routines:

       do_one

   Notes:

       The LSF extension is assumed to contain FWHM values in Angstroms on
       the same wavelength grid as the WAVE extension.

   History:

   260428 ksl Coding begun



Functions
---------

.. autoapisummary::

   fourier_offset_check.do_one
   fourier_offset_check.steer


Module Contents
---------------

.. py:function:: do_one(filename, wave=4200.0, flux=1e-11, off1=0.0, off2=0.1, off3=-0.1)

   Add a synthetic Gaussian emission line to the FLUX extension of a
   CFrame file and write the result to test_<basename>.fits.

   Parameters
   ----------
   filename : str
       Path to the input lvmCFrame FITS file.
   wave : float
       Central wavelength of the line in Angstroms.
   flux : float
       Integrated line flux in erg/s/cm^2.
   off1, off2, off3 : float
       Wavelength offset in Angstroms applied to all fibers in
       spectrographs 1, 2, and 3 respectively.


.. py:function:: steer(argv)

