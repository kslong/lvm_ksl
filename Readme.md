This repository contains tools ksl has or is developing for analyzing LVM data or for developing S/W for better sky subtraction.

Some are intended for science and others are intended to aide in the development of better sky subtraction procedures.

To use the various routines that are in this repository, one should add
this directory to your PYTHONPATH and PATH.

Most of the programs in this directory accept -h options to indicate
what they do, and the various options associated with them.

To have a good conception of the overall concept, one should understand that a data reduction pipeline exists for LVM, which one can run on the raw science exposures after downloading the raw science data.  The data reduction pipeline expects a well-defined, but complicated directory structure. The data reduction pipeline produces flux-calibrated row-stacked spectra of each exposure with names like lvmCFrame-00009240.fits which contain the row stacked flux-calibrated spectra, where 9240 is the exposure number.

For doing science, the current versions of the row-stacked spectra are lacking information about where the various fibers are positioned on the sky.  To address this problem, ksl has written S/W to obtain this information, the quality of which depends somewhat on how much ancillary information is available about each exposure, and in particular whether the data from the guide star camera exists.  If it does, it can be used to update information about the pointing positions from those that were desired; if not the desired positions are the only information that is available.

To obtain the best information available currently, ksl has written a subsidiary routine xcal.py, which (a) if one has a local version of the lvmdrp installed, and (b) which one has used to process the data with the lvmdrp, then it adds information about the pointing positions (and for the purpose of sky subtraction information about the moon) to a modified version of the calibrated RSS spectra.  This produces a file with names like XCFrame-00009240.fits, which are intended for further analysis.

## Quick Start

### Requirements

- Python 3.8 or later
- [numpy](https://numpy.org/), [scipy](https://scipy.org/), [astropy](https://www.astropy.org/), [matplotlib](https://matplotlib.org/)

Optional, for full functionality:
- **lvmdrp** — LVM Data Reduction Pipeline (needed for `Reduce.py` and local reprocessing)
- **sdss_access** — for retrieving data directly from the Utah archive
- **skycorr** — ESO sky-subtraction tool (needed for `RunSkyCorr.py`)

### Installation

Clone the repository:

```bash
git clone https://github.com/kslong/lvm_ksl.git
cd lvm_ksl
```

Add `py_progs` to your shell environment so Python can import the modules
and the scripts can be run directly from the command line:

**bash / zsh** — add to `~/.bashrc` or `~/.zshrc`:
```bash
export PYTHONPATH="/path/to/lvm_ksl/py_progs:$PYTHONPATH"
export PATH="/path/to/lvm_ksl/py_progs:$PATH"
```

**csh / tcsh** — add to `~/.cshrc`:
```csh
setenv PYTHONPATH /path/to/lvm_ksl/py_progs:$PYTHONPATH
setenv PATH /path/to/lvm_ksl/py_progs:$PATH
```

Replace `/path/to/lvm_ksl` with the directory where you cloned the repository.

### Verifying the installation

Most scripts accept `-h` to print usage information:

```bash
GetFromUtah.py -h
rss_snap.py -h
```

### Documentation

Full documentation is built with Sphinx:

```bash
pip install sphinx sphinx-rtd-theme sphinx-autoapi
cd /path/to/lvm_ksl/docs
make html
```

Open `docs/html/index.html` in a browser.  The documentation covers all
major workflows including data retrieval, sky subtraction, spectral
analysis, and the Magellanic Cloud SNR analysis pipeline.

---

## Overview

The repository contains around 60 independent Python scripts in `py_progs/`,
organised around the following workflows.  Full documentation for each is
in the Sphinx docs (see above).

**Locating observations** — cross-match a source catalog against the LVM
drpall file to find which exposures cover each target and how much
integration time has accumulated (`find_obs.py`).

**Data retrieval and reduction** — retrieve reduced data from the Utah
archive or run the LVM DRP locally, and locate processed files on disk
(`GetFromUtah.py`, `Reduce.py`, `LocateData.py`, `LocateReduced.py`).

**Summarizing data** — create compact spectral summaries across many
exposures for quality evaluation, broken down by ring position or by
spectrograph (`SummarizeData.py`, `SummarizeCframe.py`,
`SummarizeSframe.py`, `SummarizeRings.py`, `SummarizeSpec.py`).

**Sky subtraction** — prepare inputs, run SkyCorr, and evaluate sky
subtraction residuals (`Prep4SkyCorr.py`, `RunSkyCorr.py`, `SkySub.py`,
`eval_sky.py`).

**Spectral analysis** — extract spectra, fit Gaussian emission-line
profiles, and measure fluxes (`GetSpec.py`, `GetRegSpec.py`,
`lvm_gaussfit.py`, `lvm_flux.py`).

**Visualization** — produce line maps, broadband images, and spectral
overview plots (`kslmap.py`, `line_map.py`, `quick_map.py`, `PlotSpec.py`,
`PlotSpec3.py`, `SNRPlot.py`).

**RSS combining** — combine multiple dithered SFrame exposures into a
single row-stacked spectrum on a regular or original fiber grid
(`rss_combine.py`, `rss_snap.py`).

**Region files and background subtraction** — define source and
background regions, generate annular backgrounds, and extract
background-subtracted spectra (`MakeLVMReg.py`,
`GenAnnularBackground.py`).

**End-to-end science workflow** — a documented pipeline for analysing
Magellanic Cloud SNRs from exposure identification through emission-line
fitting is described in full in the Sphinx documentation.

Scripts not yet covered by detailed documentation still include a usage
summary in their module docstring; run any script with `-h` to see it.
Obsolete scripts have been moved to the `deprecated/` subdirectory.
