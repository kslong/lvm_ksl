#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

List the emission lines a DAP output file was configured to fit, by
reading the ``e_flux_*`` column names out of its NP_ELINES_B/_R/_I
extensions, and write them to a table in the same format as
data/dap_lines.txt.

Command line usage (if any):

    usage: ListDapLines.py [-h] [-out dap_lines.txt] dapfile

Description:

    The DAP's input yaml configuration controls which emission lines it
    fits, so the line set (and even individual wavelengths) can differ
    between DAP runs.  data/dap_lines.txt is not derived from atomic
    physics -- it is a transcription of the DAP's own column names -- so
    it can go stale if the yaml config changes.  This script writes its
    output locally (never touching data/dap_lines.txt) so the result can
    be diffed against the checked-in table before deciding to replace it.

    For each of NP_ELINES_B/NP_ELINES_R/NP_ELINES_I present in the file,
    every ``e_flux_<Ion>_<Wave_air>`` column names one fitted line.  Ion and
    Wave_air are split off the tail of that string (mirroring DAP2tab.py's
    get_lines()); Wave_vac is added via the standard Morton 2000 air-to-
    vacuum conversion (the same formula used in IDL astrolib's
    airtovac.pro).  LineID is written as a copy of Ion, giving a separate
    column that can be hand-edited afterwards (e.g. to shorten or
    disambiguate display labels) without touching the physical Ion id.

    If any line comes back with Ion == '[SIII]', a comment note is
    prepended flagging that the DAP's [SIII] wavelengths are its own
    legacy values and do not match those currently used in Mappings or
    Cloudy (see dap.rst) -- written with whatever wavelengths this run
    actually produced, in case a yaml change has moved them.

    -out file.txt
        Output path (default: ./dap_lines.txt in the current directory --
        deliberately not data/dap_lines.txt, so a run never overwrites the
        checked-in table; compare the two and copy over by hand if wanted).

Primary routines:

    make_dap_lines

Notes:

    Output format matches data/dap_lines.txt: fixed-width columns (widths
    5/13/13/11/11/20), left-justified, Arm/Ion/LineID/Wave_air/Wave_vac/
    DAP_name, sorted by arm (B, R, I) then by Wave_air.  Readable back
    with astropy.io.ascii.read() same as before.

History::

    260822 ksl Coding begun

'''

import os
import sys

from astropy.io import fits

from PlotSpec import _usage_from_doc


DEFAULT_OUTPUT = 'dap_lines.txt'

ARMS = (('B', 'NP_ELINES_B'), ('R', 'NP_ELINES_R'), ('I', 'NP_ELINES_I'))

WIDTHS = (5, 13, 13, 11, 11, 20)


def airtovac(wave_air):
    '''
    Convert an air wavelength (Angstroms) to vacuum, via the Morton 2000
    formula also used by IDL astrolib's airtovac.pro.
    '''
    sigma2 = (1e4 / wave_air) ** 2
    n = 1 + 6.4328e-5 + 2.94981e-2 / (146 - sigma2) + 2.5540e-4 / (41 - sigma2)
    return wave_air * n


def get_dap_names(colnames):
    '''
    Strip the ``e_flux_`` prefix off column names, returning the bare
    ``<Ion>_<Wave_air>`` strings.  Mirrors DAP2tab.py's get_lines().
    '''
    return [c[len('e_flux_'):] for c in colnames if c.startswith('e_flux_')]


def make_dap_lines(dapfile):
    '''
    Read dapfile and return a list of row dicts (Arm, Ion, LineID,
    Wave_air, Wave_vac, DAP_name), one per fitted line, sorted by arm
    (B, R, I) then by Wave_air.
    '''
    rows = []
    with fits.open(dapfile) as hdul:
        extnames = [hdu.name for hdu in hdul]
        for arm, extname in ARMS:
            if extname not in extnames:
                continue
            colnames = hdul[extname].columns.names
            for dap_name in get_dap_names(colnames):
                ion, wave_str = dap_name.rsplit('_', 1)
                wave_air = float(wave_str)
                wave_vac = airtovac(wave_air)
                rows.append(dict(Arm=arm, Ion=ion, LineID=ion,
                                  Wave_air=wave_air, Wave_vac=wave_vac,
                                  DAP_name=dap_name))

    arm_order = {arm: i for i, (arm, _) in enumerate(ARMS)}
    rows.sort(key=lambda r: (arm_order[r['Arm']], r['Wave_air']))
    return rows


def _format_row(fields):
    return ''.join(str(f).ljust(w) for f, w in zip(fields, WIDTHS))


def _siii_note(rows):
    siii = [r for r in rows if r['Ion'] == '[SIII]']
    if not siii:
        return []
    waves = ', '.join('%.2f' % r['Wave_air'] for r in siii)
    return [
        '# NOTE: [SIII] wavelengths here (%s air) are as returned' % waves,
        '# by the DAP and are NOT the currently accepted values -- they do not match',
        '# those used in Mappings or Cloudy. Do not treat them as authoritative without',
        '# checking against current atomic data.',
    ]


def write_dap_lines(rows, outfile):
    lines = _siii_note(rows)
    lines.append(_format_row(['Arm', 'Ion', 'LineID', 'Wave_air', 'Wave_vac', 'DAP_name']))
    for r in rows:
        lines.append(_format_row([r['Arm'], r['Ion'], r['LineID'],
                                   '%.2f' % r['Wave_air'], '%.2f' % r['Wave_vac'],
                                   r['DAP_name']]))
    with open(outfile, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def steer(argv):
    '''
    Steering routine: parse the command line and write the line table.
    '''
    outfile = DEFAULT_OUTPUT
    dapfile = None

    i = 1
    while i < len(argv):
        if argv[i][0:2] == '-h':
            print(_usage_from_doc(__doc__))
            return
        elif argv[i] == '-out':
            i += 1
            outfile = argv[i]
        elif argv[i][0] == '-':
            print('Error: Unknown switch: ', argv[i])
            return
        else:
            dapfile = argv[i]
        i += 1

    if dapfile is None:
        print('Error: No DAP file specified')
        return

    if not os.path.isfile(dapfile):
        print('Error: Could not find %s' % dapfile)
        return

    rows = make_dap_lines(dapfile)
    if not rows:
        print('Error: No e_flux_* columns found in NP_ELINES_B/R/I extensions of %s' % dapfile)
        return

    write_dap_lines(rows, outfile)
    print('Wrote %d lines to %s' % (len(rows), outfile))


if __name__ == "__main__":
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
