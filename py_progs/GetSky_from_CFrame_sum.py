#!/usr/bin/env python
"""
GetSky_from_CFrame_sum.py
=========================

Synopsis
--------
Summarise sky field observation counts from an LVM XCframe summary file, or
extract all sky spectra for a named sky field into a new FITS file.

Description
-----------
LVM exposures record both an east and a west sky field alongside each science
pointing.  The DRP_ALL table inside an XCframe summary file carries the field
names in the columns skye_name and skyw_name.  This program has two operating
modes depending on whether a source name is supplied:

Summary mode (no source name)
    Counts how many times each sky field appears across both the east and west
    sky telescopes, prints the top --n-print entries sorted by observation
    count (most observed first), and writes the full count table to a FITS
    file named sky_summary_<input_stem>.fits.

Extraction mode (source name given)
    Collects every row in the DRP_ALL table where skye_name or skyw_name
    matches the requested source, extracts the corresponding sky spectra from
    the SKY_EAST and SKY_WEST extensions, and writes them to a FITS file
    named Sky_<source_name>.fits (or the path given by --output).

Output FITS structure (extraction mode)
----------------------------------------
  PRIMARY    No spectral data; header records Source and input filename.
  WAVE       float32  (Npix,)            Wavelength array (Å).
  FLUX       float32  (Nobs, Npix)       Sky spectra; one row per observation.
  DRP_ALL    BinTable (Nobs rows)        Subset of the input DRP_ALL table for
                                         the selected source, with three extra
                                         columns added:
                                           line_no      — original row index in
                                                          the input DRP_ALL table
                                           Source_name  — the requested name
                                           tel          — 'SKY_EAST' or 'SKY_WEST'

Usage
-----
    # Summary mode — count observations per sky field
    python GetSky_from_CFrame_sum.py <fits_file>

    # Extraction mode — extract spectra for one sky field
    python GetSky_from_CFrame_sum.py <fits_file> <source_name>

    # Custom output filename
    python GetSky_from_CFrame_sum.py <fits_file> <source_name> --output my_sky.fits

Arguments
---------
    fits_file     Path to the LVM XCframe summary FITS file.  Must contain
                  a DRP_ALL BinTableHDU with skye_name and skyw_name columns,
                  a WAVE ImageHDU, and SKY_EAST / SKY_WEST ImageHDUs.
    source_name   (Optional) Name of the sky field to extract.  Must match
                  the values in the skye_name or skyw_name columns exactly.

Options
-------
    --output PATH   Output FITS path for extraction mode.
                    (default: Sky_<source_name>.fits)
    --n-print N     Number of top sky fields to print in summary mode.
                    (default: 10)

Notes
-----
History
    260628  Written based on LocateSky.ipynb (Knox Long / Claude)
"""

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack


# ---------------------------------------------------------------------------
# Summary mode
# ---------------------------------------------------------------------------

def summarise_sky_fields(fits_file, n_print=10):
    """Count observations per sky field and return a sorted astropy Table.

    Both east (skye_name) and west (skyw_name) sky telescope pointings are
    counted together.  Empty or blank field names are excluded.

    Parameters
    ----------
    fits_file : str or Path
        Path to the XCframe summary FITS file.
    n_print : int
        Number of top entries to print to stdout.

    Returns
    -------
    Table
        Columns Source_name (str) and Count (int), sorted descending by Count.
    """
    with fits.open(fits_file) as hdul:
        drp = Table(hdul["DRP_ALL"].data)

    skye = np.array(drp["skye_name"], dtype=str)
    skyw = np.array(drp["skyw_name"], dtype=str)
    sky  = np.concatenate([skye, skyw])

    # remove blank / whitespace-only entries
    sky = sky[np.char.strip(sky) != ""]

    names, counts = np.unique(sky, return_counts=True)
    tab = Table([names, counts], names=["Source_name", "Count"])
    tab.sort("Count", reverse=True)

    print(f"\nSky field observation counts from {Path(fits_file).name}")
    print(f"Total unique sky fields: {len(tab)}")
    print(f"\nTop {min(n_print, len(tab))} sky fields by observation count:\n")
    print(f"  {'Source_name':<20}  {'Count':>6}")
    print("  " + "-" * 28)
    for row in tab[:n_print]:
        print(f"  {row['Source_name']:<20}  {row['Count']:>6d}")

    return tab


def save_summary(tab, fits_file):
    """Write the sky field count table to a FITS BinTable file."""
    outpath = f"sky_summary_{Path(fits_file).stem}.fits"
    hdr = fits.Header()
    hdr["INPUT"] = (str(fits_file), "source XCframe summary file")
    fits.HDUList([
        fits.PrimaryHDU(header=hdr),
        fits.BinTableHDU(tab, name="SKY_COUNTS"),
    ]).writeto(outpath, overwrite=True)
    print(f"\nFull summary table written to {outpath}")


# ---------------------------------------------------------------------------
# Extraction mode
# ---------------------------------------------------------------------------

def get_sky_spectra(fits_file, source_name, outpath=None):
    """Extract sky spectra for a named sky field from an XCframe summary file.

    Rows where skye_name == source_name are drawn from the SKY_EAST extension;
    rows where skyw_name == source_name are drawn from SKY_WEST.  The combined
    spectra and metadata are written to a FITS file.

    Parameters
    ----------
    fits_file : str or Path
        Path to the XCframe summary FITS file.
    source_name : str
        Exact sky field name to extract (matched against skye_name / skyw_name).
    outpath : str or Path or None
        Output FITS path.  Defaults to Sky_<source_name>.fits in the current
        directory.
    """
    if outpath is None:
        outpath = f"Sky_{source_name}.fits"

    print(f"Reading {fits_file} ...")
    with fits.open(fits_file) as hdul:
        drp      = Table(hdul["DRP_ALL"].data)
        wave     = hdul["WAVE"].data
        sky_east = hdul["SKY_EAST"].data   # shape (n_fibers, n_wave) in numpy
        sky_west = hdul["SKY_WEST"].data

    # add a row-index column so we can select spectra by position
    drp["line_no"] = np.arange(len(drp))
    drp["Source_name"] = source_name

    skye_mask = np.char.strip(np.array(drp["skye_name"], dtype=str)) == source_name
    skyw_mask = np.char.strip(np.array(drp["skyw_name"], dtype=str)) == source_name

    skye_rows = drp[skye_mask].copy()
    skyw_rows = drp[skyw_mask].copy()
    skye_rows["tel"] = "SKY_EAST"
    skyw_rows["tel"] = "SKY_WEST"

    n_east = int(skye_mask.sum())
    n_west = int(skyw_mask.sum())
    n_total = n_east + n_west

    if n_total == 0:
        print(f"ERROR: source '{source_name}' not found in skye_name or skyw_name.")
        return

    print(f"Source: {source_name}")
    print(f"  SKY_EAST observations: {n_east}")
    print(f"  SKY_WEST observations: {n_west}")
    print(f"  Total:                 {n_total}")

    skye_spec = sky_east[skye_rows["line_no"]] if n_east > 0 else np.empty((0, wave.size))
    skyw_spec = sky_west[skyw_rows["line_no"]] if n_west > 0 else np.empty((0, wave.size))

    tab  = vstack([skye_rows, skyw_rows])
    spec = np.vstack([skye_spec, skyw_spec])   # shape (n_total, n_wave)

    hdr = fits.Header()
    hdr["TITLE"]  = "Sky_Summary"
    hdr["SOURCE"] = (source_name, "sky field name")
    hdr["INPUT"]  = (str(fits_file), "source XCframe summary file")
    hdr["N_EAST"] = (n_east,  "number of SKY_EAST observations")
    hdr["N_WEST"] = (n_west,  "number of SKY_WEST observations")
    hdr["N_TOTAL"] = (n_total, "total sky observations")

    fits.HDUList([
        fits.PrimaryHDU(header=hdr),
        fits.ImageHDU(data=wave,               name="WAVE"),
        fits.ImageHDU(data=spec.astype(np.float32), name="FLUX"),
        fits.BinTableHDU(tab,                  name="DRP_ALL"),
    ]).writeto(outpath, overwrite=True)

    print(f"\nOutput written to {outpath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=("Summarise sky field counts or extract sky spectra from "
                     "an LVM XCframe summary FITS file."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("fits_file",
                   help="LVM XCframe summary FITS file")
    p.add_argument("source_name", nargs="?", default=None,
                   help="Sky field name to extract (omit for summary mode)")
    p.add_argument("--output", default=None,
                   help="Output FITS path (extraction mode only; "
                        "default: Sky_<source_name>.fits)")
    p.add_argument("--n-print", type=int, default=10,
                   help="Number of top sky fields to print in summary mode")
    args = p.parse_args()

    if args.source_name is None:
        tab = summarise_sky_fields(args.fits_file, n_print=args.n_print)
        save_summary(tab, args.fits_file)
    else:
        get_sky_spectra(args.fits_file, args.source_name, outpath=args.output)


if __name__ == "__main__":
    main()
