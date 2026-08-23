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

**Rule: rebuild the docs and finalize History entries at commit time, not after every iterative edit.** During a work session, docstrings/History may go through several drafts as the approach changes — don't rebuild docs or polish History after each one; that's wasted effort if the design isn't settled yet. Once the change is actually ready to commit: rebuild the docs and fix any new warnings/errors, and consolidate History into a clean final entry (or entries) that reflects what was actually built, not a blow-by-blow of intermediate attempts.

```bash
cd docs && rm -rf source/api html doctrees && make html
```
(`docs/source/api/` is AutoAPI-generated and gitignored — force a clean rebuild rather than trusting a stale copy on disk, since `autoapi_keep_files` leaves old generated `.rst` around otherwise.) Watch for `ERROR`/`WARNING` lines from `docutils` in the output.

The dominant, recurring cause is multi-line docstring sections (`History:`, `Parameters:`, numbered/columnar lists) where a wrapped continuation line is indented differently from its first line — RST reads that as a malformed nested block ("Unexpected indentation" / "Block quote ends without a blank line"). **The fix is always the same, apply it directly on sight — do not hand-tune indentation or re-derive a fix by trial and error:** change the section header to end in `::` (e.g. `History::` not `History:`) and indent the whole body at least one level deeper than the header. This makes it an RST literal block, which is never reparsed, so internal wrapping/indentation can't break it. RST's tolerance for wrapped/deeper-indented continuation lines under a plain `:` header is inconsistent and context-dependent — near-identical patterns sometimes parse and sometimes don't — so don't trust a passing test on one such section as a reason to leave it as `:`; any multi-line section gets `::`. Applied repo-wide 260707–260716.

Other docstring RST gotchas:
- A literal `|word|` (e.g. `sqrt(|flux|)`) is parsed as an undefined substitution reference — reword without pipes or wrap in double backticks.
- A bare word ending in `_` in free prose (e.g. `SkyE_`) is parsed as an implicit hyperlink target — wrap in double backticks.
- A single unmatched `*` (e.g. `drpall-*.fits`) is parsed as the start of emphasis markup — wrap in double backticks, or put it inside a section already marked `::`.
- A `Parameters:`/`Returns:` docstring without a colon after each name (relying on column alignment) confuses `sphinx.ext.napoleon`'s Google/NumPy parser — either use `name: description` per line, or mark the whole section `::` like `History`.

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
| Spectral fitting | `lvm_gaussfit.py`, `lvm_snrfit.py`, `lvm_flux.py`, `lvm_double.py`, `lvm_triple.py` |
| Astrometry | `fib2radec.py` (fiber → RA/Dec conversion) |
| Imaging | `kslmap.py`, `quick_map.py`, `line_map.py` |
| Spectrum overview plots | `PlotSpec.py`, `PlotSpec3.py`, `PlotSpecI.py` (interactive/Plotly) |
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
