Paranal Airglow Line And Continuum Emission (PALACE) model: data and code of
v1.0

Authors: Stefan Noll (st@noll-x.de), Carsten Schmidt, Patrick Hannawald,
         Wolfgang Kausch, and Stefan Kimeswenger
Year: 2024

This is a brief description of the PALACE model data (PMD) release. It is
connected to the paper "PALACE v1.0: Paranal Airglow Line And Continuum
Emission model" of the same authors that should be published in Geoscientific
Model Development (GMD). The purpose of this release is to provide data that
were used to build and evaluate the model. In this way, it should also be
possible to reproduce the figures of the paper. All data are provided as ASCII
tables.

PALACE produces spectra in the wavelength range from 0.3 to 2.5 µm that show
the nightglow line and continuum emission for different chemical species in
the middle and upper atmosphere above Cerro Paranal in Chile for different
observing conditions. The latter include zenith angle, local time, month, and
solar activity. The main data source of the model are 10 years of spectra of
the astronomical echelle spectrograph X-shooter of the Very Large Telescope
at Cerro Paranal. These data were partly supplemented by data from the UVES
echelle spectrograph that is located at the same site. Apart from
observational data, the development of the model also relied on theoretical
molecular and atomic properties such as wavelengths, energy levels, Einstein-A
coefficients, and recombination coefficients.

pmd_refspec_air.dat:
PALACE reference spectra with atmospheric extinction for the total emission
and for different chemical species

pmd_refspec_vac.dat:
PALACE reference spectra without atmospheric extinction for the total emission
and for different chemical species

pmd_reflineint_mol.dat:
PALACE reference intensities for molecular lines

pmd_corrA_OHbands.dat:
Correction of OH Einstein-A coefficients from HITRAN2020 for each branch and
band

pmd_popdata_OH.dat:
OH level populations from measured reference intensities of Lambda doublets
(based on X-shooter/UVES data) that were used for the OH population model

pmd_popmodel_OH.dat:
Reference populations and intensities from the OH population model for lines
with HITRAN2020 data

pmd_popdata_O2.dat:
O2 level populations from measured reference intensities (based on
X-shooter/UVES data) that were used for the O2 population model

pmd_popmodel_O2.dat:
Reference populations and intensities from the O2 population model for lines
with HITRAN2020 data

pmd_reflineint_atom.dat:
PALACE reference intensities for atomic lines

pmd_intdata_atom.dat:
List of reference intensities of atomic lines/multiplets with measured
X-shooter or UVES time series that were used for the derivation of the PALACE
line intensity model

pmd_intmodel_Orc.dat:
Illustration of the calculation of the PALACE reference intensities for O
recombination lines

pmd_refcont.dat:
Reference continuum components of PALACE based on X-shooter-related
measurements

pmd_refclim:
Folder with 2D reference climatologies of relative intensity, solar cycle
effect, and residual variability for 23 PALACE variability classes (23 files)

pmd_climfeat.dat:
List of lines/multiplets/features with measured X-shooter/UVES time series
that were used for the derivation of the PALACE reference climatologies

pmd_climdata:
Folder with 23 subdirectories for each variability class that include the
intensity time series from X-shooter data (UVES in the case of K) that were
used for the derivation of the PALACE reference climatologies (396 files)

pmd_reseval_UVB.dat, pmd_reseval_VIS.dat, pmd_reseval_NIR.dat:
Results of the evaluation of PALACE and the airglow component of the ESO Sky
Model in comparison to X-shooter UVB/VIS/NIR-arm spectra

Note: Some file headers refer to the files palace_lines.fits,
palace_cont.fits, and palace_var.fits. These are the input files for the
PALACE code, which is also part of the release.
