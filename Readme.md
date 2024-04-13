This repository contains tooks tools ksl has or is developing for anaylyzing LVM data or for develping S/W for better sky subtraction.

Some are intended for scince and others are intended to aide in the development of better sky subtraction procecures

To use the various routines that are in this repository, one should add
this directory to your PYTHONPATH and PATH

Most of the programs in this directory accept -h options to inindicate
what they do, and the various options associated with them

To have a good coneption of the overall concept, one should understand that a data reduction pipeline exists for LVM, which one can run on the raw science exposures after downloading the raw scince data.  The data reduction pipepline expects a well-defined, but compliated directory structure. The data reduction pipeline produces flux-calibrated row-stacked spectra of each exposure with names like lvmCFrame-00009240.fits which contain the row stacked flux-calibrated spectra, where 9240 is the exposure number.  

For doing, science, the current versions of the row-stacked spectra are lacking information about where the various fibers are positioned on the sky.  To address this problem, ksl has written S/W to obtain this information, the quality of which depends somewhat on how much ancilary information is available about each exposure, and in particular whether the data from the guide star camera exists.  If it does, it can be used to update information abut the pointing positions from those that were desired; if not the desired postions are the only information that is available.

To obtain the best information available currently, ksl has written a subsidiary routine xcal.py, which (a) if one has a local version of the lvmdrp instaled, and (b) which one has used to process the data with the lvmdrp, then it adds information about the pointing positions (and for purpose of sky subraction information about the moon) to a modified version of the calibrated RSS spectra.  This produces a file with names like  XCFrame-00009240.fits, which are intended for further analysis.  

Here is a brief description of each program:

* xcal.py adds the ra and decs for individual fibers in a calibrated exposure. It also adds information needed to create a WCS that is used by kslmap to create an image of an RSS file with good astrometry. The LVM DRP (data analysis pipelines) crates files with names like lvmCFrame-00009240.fits which contain the row stacked flux-calibrated spectra, where 9240 is the exposure number. xcal.py files will have names like XCFrame-00009240.fits

* kslmap.py creates 'images" of an exposure from either the lvm drp created reduced data or the xcal data. 

* PlotClosest.py creates an ansci table fo the fiber closest to a given ra and dec.  It can run on either the direct output of xcal.py or on the direct output of the lvmdrp.  

* fib2radec.py creates an astropy table giving the RA and Dec of each fibeer in a calbirated LVM imaage; if the guide star data are avalible the program will use information from that ,but if not it will use data in the calbrated file to determin the ra's and dec's of each fiber.

The remaining programs are mainly intended for use when one has a local version of the LVM DRP installed, or for work on sky subtraction

* Reduce.py downlaods data from Utah and processes it with the quick-reduction version of the LVM drp.

* LocateData.py is intended to find a lvmCFrame file locally or in the file structure that is part of LVMDRP.  If multiple versions of a file exist, it returns the name of
the last one created.  

* PrepSkyCorr.py creates data in a format that can be used by SkyCorr, specificaly, it creates 3 different fits files for the science, SkyE and SkyW telescopes

* RunSkyCorr.py allows one to run SkyCoor on calibrated data.  The routine creates the .par files that SkyCorr needs and then runs SkyCorr on the Prepped files.

* SkySub.py allows one to run a simple sky subtraction (separating the continuum from the sky lines, and rescaling the sky lines to minimize the differences between the science and sky spectrum.  This uses the same input files as SkyCorr.

* SummarizeData.py - Summarizes the raw data that has been dowwnloaded.  Used on the Utah cluster one can make a summary of all of the data that has been obtained

* SkyID.py is intended to identify what SkyFields are associated with what pointings.  It uses the output of SummarizeData.py.

* CleanAncillary.py - removes the ancillary data files that have been created by processing (which can take space if one is not intereted in them)


The repository will likely contain other routines as time goes forward.  Hopefully, they will be docuemnted, but usually there is at least some indication of what they do in the headers information at the top of the routines.
