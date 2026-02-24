#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

This is intended to create a standardized plot for looking at spectra
of SNRs and SNR candidates extracted from LVM data.


Command line usage (if any):

    usage: SNRPlot.py [-v_offset km/s] [-smc] [-lmc] file [files ...]

    Plots are written to SNRPlot_fig/<basename>.png.
    A README.txt is written to SNRPlot_fig/ describing the run.

    -v_offset  systemic velocity offset in km/s (default 280)
    -smc       set v_offset to 146 km/s (SMC)
    -lmc       set v_offset to 280 km/s (LMC, default)

Description:

Each plot has a top panel covering the full wavelength range with
FLUX and ERROR overlaid, and six sub-panels zoomed around key
emission lines: [OII] 3728, [OIII] 5007, [OI] 6300, Ha 6563,
[SII] 6720, and [SIII] 9531.

Primary routines:

    xplot

Notes:

History:

241111 ksl Coding begun
260224 ksl Add v_offset flag, SNRPlot_fig output directory, README

'''


import os
import sys
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from astropy.io import ascii


OUT_DIR = 'SNRPlot_fig'


def limit_spectrum(xtab, w_cen=3728, vlim=1000, v_offset=280):
    ztab = xtab.copy()
    wmin = w_cen * (1 - vlim / 3e5)
    wmax = w_cen * (1 + vlim / 3e5)
    delta = w_cen * v_offset / 3e5
    ztab['WAVE'] -= delta
    ztab = ztab[ztab['WAVE'] > wmin]
    ztab = ztab[ztab['WAVE'] < wmax]
    return ztab


def xplot(filename='Spec_09444_test_ave.txt', outroot='', vlim=1200, v_offset=280):
    xtab = ascii.read(filename)
    plt.close(1)
    fig = plt.figure(1, (8, 8))
    gs = GridSpec(3, 3, figure=fig)

    ax1 = fig.add_subplot(gs[0, :])
    ax1.set_title(os.path.basename(filename))
    ax1.plot(xtab['WAVE'], xtab['FLUX'], label='Spectrum', zorder=2)
    ax1.plot(xtab['WAVE'], xtab['ERROR'], label='Error', zorder=2)
    ax1.set_xlim(3600, 9700)
    ymax = np.nanmax(xtab['FLUX'])
    ax1.set_ylim(-0.1 * ymax, 1.1 * ymax)
    ax1.legend()

    line_panels = [
        (gs[1, 0], 3728,   '[OII]'),
        (gs[1, 1], 5007,   '[OIII]'),
        (gs[1, 2], 6300,   '[OI]'),
        (gs[2, 0], 6563,   r'H$\alpha$'),
        (gs[2, 1], 6720,   '[SII]'),
        (gs[2, 2], 9530.6, '[SIII]'),
    ]
    for gs_cell, w_cen, label in line_panels:
        ax2 = fig.add_subplot(gs_cell)
        ztab = limit_spectrum(xtab, w_cen, vlim, v_offset)
        ax2.plot(ztab['WAVE'], ztab['FLUX'], label=label)
        ax2.plot(ztab['WAVE'], ztab['ERROR'])
        if len(ztab) > 0:
            ax2.set_xlim(ztab['WAVE'][0], ztab['WAVE'][-1])
            ymax = np.max(ztab['FLUX'])
            ymin = np.min(ztab['FLUX'])
            if ymin < 0 and np.median(ztab['FLUX']) > 0:
                ax2.set_ylim(-0.1 * ymax, 1.1 * ymax)
        ax2.legend()

    plt.tight_layout()

    if outroot == '':
        outroot = filename.replace('.txt', '').replace('.tab', '')
    plt.savefig(outroot + '.png')


def steer(argv):
    '''
    Make a plot of spectra extracted from a region of LVM data, with
    an emphasis on lines useful for looking at SNR data.
    '''
    files = []
    v_offset = 280
    vlim = 1200

    i = 1
    while i < len(argv):
        if argv[i][0:2] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-v_offset':
            i += 1
            v_offset = eval(argv[i])
        elif argv[i] == '-smc':
            v_offset = 146
        elif argv[i] == '-lmc':
            v_offset = 280
        elif argv[i][0] == '-':
            print('Unknown switch: ', argv[i])
            return
        else:
            files.append(argv[i])
        i += 1

    if not files:
        print('Error: No input files specified')
        return

    os.makedirs(OUT_DIR, exist_ok=True)

    readme = os.path.join(OUT_DIR, 'README.txt')
    with open(readme, 'w') as f:
        f.write('SNRPlot_fig\n')
        f.write('===========\n')
        f.write('Overview plots produced by SNRPlot.py\n\n')
        f.write('Generated : %s\n' % datetime.date.today().isoformat())
        f.write('Command   : %s\n' % ' '.join(argv))
        f.write('v_offset  : %g km/s\n' % v_offset)
        f.write('vlim      : %g km/s\n' % vlim)
        f.write('\nEach PNG shows the full spectrum in the top panel and\n')
        f.write('six key emission lines ([OII] 3728, [OIII] 5007, [OI] 6300,\n')
        f.write('Ha 6563, [SII] 6720, [SIII] 9531) in sub-panels, with\n')
        f.write('FLUX and ERROR overlaid in each panel.\n')

    for one in files:
        basename = os.path.basename(one).replace('.txt', '').replace('.tab', '')
        outroot = os.path.join(OUT_DIR, basename)
        xplot(filename=one, outroot=outroot, vlim=vlim, v_offset=v_offset)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        steer(sys.argv)
    else:
        print(__doc__)
