#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

    Summarise sky field observation counts from an LVM XCframe summary file,
    or extract all sky spectra for a named sky field into a new FITS file.

Command line usage (if any):

    usage: GetSky_from_CFrame_sum.py [-h] [--output PATH] [--n-print N]
                                     fits_file [source_name]

    where

    fits_file     is the path to an LVM XCframe summary FITS file containing
                  a DRP_ALL table (with skye_name and skyw_name columns),
                  a WAVE extension, and SKY_EAST / SKY_WEST extensions.

    source_name   (optional) exact sky field name to extract.  If omitted
                  the routine runs in summary mode.

    --output PATH sets the output FITS filename (extraction mode only;
                  default: Sky_<source_name>.fits).

    --n-print N   number of top sky fields to print in summary mode
                  (default: 10).

Description:

    LVM exposures observe both an east and a west sky field alongside each
    science pointing.  The DRP_ALL table inside an XCframe summary file
    records which sky fields were observed in the skye_name and skyw_name
    columns.  This program operates in two modes:

    Summary mode (no source_name given):
        Counts how many times each sky field appears across both the east and
        west sky telescopes, prints the top N entries sorted by count (most
        observed first), and writes the full count table to a FITS file named
        sky_summary_<input_stem>.fits.  This is useful for identifying which
        sky fields have been observed enough times to build a high-quality
        stacked sky spectrum.

    Extraction mode (source_name given):
        Collects every row in DRP_ALL where skye_name or skyw_name matches
        the requested source, extracts the corresponding spectra from the
        SKY_EAST and SKY_WEST extensions, and writes them to a FITS file.
        Paired skye_*/skyw_* metadata columns (ra, dec, amass, alt, etc.)
        are merged into single columns reflecting the correct telescope for
        each spectrum row.  Science-pointing columns (sci_ra, sci_dec, etc.)
        are removed as they do not pertain to the sky field.

    Output FITS structure (extraction mode):
        PRIMARY   header: SOURCE, INPUT, N_EAST, N_WEST, N_TOTAL
        WAVE      float32 (Npix,)       wavelength array (Å)
        FLUX      float32 (Nobs, Npix)  sky spectra; one row per observation
        DRP_ALL   BinTable (Nobs rows)  metadata; added cols: line_no, Source_name, tel

Primary routines:

    summarise_sky_fields
    get_sky_spectra

Notes:

    If source_name is not found in the file, the program exits with a
    non-zero status and prints the top 10 available sky fields.

History:

    260628  ksl  Coding begun based on LocateSky.ipynb
'''

import argparse
import sys
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
# Sky column merging
# ---------------------------------------------------------------------------

# Paired skye_*/skyw_* columns that are collapsed into a single column in
# the output table.  For each row, the value is taken from the skye_* column
# when tel == 'SKY_EAST' and from the skyw_* column when tel == 'SKY_WEST'.
_PAIRED_COLS = [
    ("skye_ra",       "skyw_ra",       "ra"),
    ("skye_dec",      "skyw_dec",      "dec"),
    ("skye_pa",       "skyw_pa",       "pa"),
    ("skye_amass",    "skyw_amass",    "amass"),
    ("skye_astsrc",   "skyw_astsrc",   "astsrc"),
    ("skye_kmpos",    "skyw_kmpos",    "kmpos"),
    ("skye_focpos",   "skyw_focpos",   "focpos"),
    ("skye_name",     "skyw_name",     "name"),
    ("skye_alt",      "skyw_alt",      "alt"),
    ("skye_sh_hght",  "skyw_sh_hght",  "sh_hght"),
    ("skye_moon_sep", "skyw_moon_sep", "moon_sep"),
]

# Columns removed from the output table because they describe the science
# pointing rather than the sky field (sci_mean_sen* are kept).
_DROP_COLS = [
    "sci_ra", "sci_dec", "sci_pa", "sci_amass", "sci_astsrc",
    "sci_kmpos", "sci_focpos", "sci_alt", "sci_sh_hght", "sci_moon_sep",
    "sci_skye_sep", "sci_skyw_sep",
]


def _merge_sky_columns(tab):
    """Replace paired skye_*/skyw_* columns with a single merged column.

    For each (skye_col, skyw_col, merged_name) triplet in _PAIRED_COLS, a new
    column named merged_name is added whose value is taken from skye_col for
    rows where tel == 'SKY_EAST' and from skyw_col for rows where tel ==
    'SKY_WEST'.  Both original columns are then removed.

    Parameters
    ----------
    tab : astropy Table
        Combined table from vstack([skye_rows, skyw_rows]).  Must contain a
        'tel' column with values 'SKY_EAST' or 'SKY_WEST'.

    Returns
    -------
    astropy Table
        Table with merged columns substituted for the paired originals.
    """
    is_east = np.array(tab["tel"]) == "SKY_EAST"

    for col_e, col_w, col_out in _PAIRED_COLS:
        if col_e not in tab.colnames or col_w not in tab.colnames:
            continue
        vals_e = np.array(tab[col_e])
        vals_w = np.array(tab[col_w])
        tab[col_out] = np.where(is_east, vals_e, vals_w)
        tab.remove_columns([col_e, col_w])

    to_drop = [c for c in _DROP_COLS if c in tab.colnames]
    if to_drop:
        tab.remove_columns(to_drop)

    return tab


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

    n_east  = int(skye_mask.sum())
    n_west  = int(skyw_mask.sum())
    n_total = n_east + n_west

    if n_total == 0:
        # build the available source list so the user knows what to ask for
        sky = np.concatenate([
            np.array(drp["skye_name"], dtype=str),
            np.array(drp["skyw_name"], dtype=str),
        ])
        sky = sky[np.char.strip(sky) != ""]
        available, counts = np.unique(sky, return_counts=True)
        order = np.argsort(counts)[::-1]

        print(f"\nERROR: '{source_name}' not found in skye_name or skyw_name.")
        print(f"       {len(available)} sky fields are present in {Path(fits_file).name}.")
        print(f"       Top 10 by observation count:")
        print(f"\n         {'Source_name':<20}  {'Count':>6}")
        print(f"         " + "-" * 28)
        for i in order[:10]:
            print(f"         {available[i]:<20}  {counts[i]:>6d}")
        print(f"\n       Run without a source name to see the full list.")
        sys.exit(1)

    skye_rows = drp[skye_mask].copy()
    skyw_rows = drp[skyw_mask].copy()
    skye_rows["tel"] = "SKY_EAST"
    skyw_rows["tel"] = "SKY_WEST"

    print(f"Source: {source_name}")
    print(f"  SKY_EAST observations: {n_east}")
    print(f"  SKY_WEST observations: {n_west}")
    print(f"  Total:                 {n_total}")

    skye_spec = sky_east[skye_rows["line_no"]] if n_east > 0 else np.empty((0, wave.size))
    skyw_spec = sky_west[skyw_rows["line_no"]] if n_west > 0 else np.empty((0, wave.size))

    tab  = _merge_sky_columns(vstack([skye_rows, skyw_rows]))
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
