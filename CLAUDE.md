# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**lvm_ksl** is an astronomical data analysis toolkit for the SDSS-V Local Volume Mapper (LVM) survey. It provides tools for data retrieval, reduction pipeline integration, sky subtraction, spectral analysis, and visualization. Developed at STScI by Knox Long and Sean Points.

## Setup

No formal package installation. Add to Python path:
```bash
export PYTHONPATH="/path/to/lvm_ksl/py_progs:$PYTHONPATH"
```

## Documentation

```bash
cd docs && make html
```
Builds Sphinx documentation with AutoAPI. Output in `docs/html/`.

## Architecture

**Script-based design**: 58 independent executable Python programs in `py_progs/`, each handling a specific task with command-line interface. Most scripts accept `-h` for usage help.

### Key Workflows

1. **Data Pipeline**: `GetFromUtah.py` → `Reduce.py` → analysis scripts
2. **Sky Subtraction**: `Prep4SkyCorr.py` → `RunSkyCorr.py` → `eval_sky.py`
3. **Spectral Analysis**: `GetSpec.py` → `lvm_gaussfit.py` → `lvm_flux.py`
4. **Visualization**: `kslmap.py`, `quick_map.py`, `line_map.py`

### Core Scripts by Function

| Category | Scripts |
|----------|---------|
| Data retrieval | `GetFromUtah.py`, `Reduce.py`, `GetDAP.py` |
| File location | `LocateData.py`, `SummarizeData.py` |
| Sky subtraction | `SkySub.py`, `RunSkyCorr.py`, `Prep4SkyCorr.py` |
| Spectral fitting | `lvm_gaussfit.py`, `lvm_flux.py`, `lvm_double.py`, `lvm_triple.py` |
| Astrometry | `fib2radec.py` (fiber → RA/Dec conversion) |
| Imaging | `kslmap.py`, `quick_map.py`, `line_map.py` |
| Data quality | `CheckReduced.py`, `CheckData.py` |

### Key Dependencies

- **astropy**: Tables, FITS I/O, WCS, coordinates, units
- **scipy**: `curve_fit` for Gaussian fitting, `bisect` for sky subtraction, `griddata` for interpolation
- **External**: lvmdrp (LVM DRP), sdss_access (Utah data), skycorr (sky subtraction tool)

### Data Formats

- Works with lvmCFrame and SFrame files from LVM DRP
- Uses astropy Tables extensively for data manipulation
- Configuration via `lvm_base.par` (SkyCorr parameters)

### Directory Structure

- `py_progs/` - All executable scripts
- `data/` - Reference data (sky tiles, solar flux, SNR observations)
- `deprecated/` - Legacy code kept for reference
- `docs/` - Sphinx documentation

### Important Constants

- Plate scale: 112.36748321030637 arcsec/mm
- Environment variables: `LVM_MASTER_DIR`, `LVMAGCAM_DIR`

## Code Style

- Docstrings with Synopsis, Description, Notes, History sections
- Version dates in YYMMDD format
- argparse for command-line parsing
- Heavy use of astropy Table operations
