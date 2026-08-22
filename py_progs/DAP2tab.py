#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Parse the DAP output file or files, and produce an astropy 
table with the information needed to make simple plots


Command line usage (if any):

    usage: DAP2tab.py [-h] filename(s)

Description:  

Primary routines:

    doit

Notes:

History::

    241224 ksl Coding begun
    260822 ksl get_radec_fluxes(): replaced 18 near-duplicate try/except
        blocks with a single (arm_table, DAP_name, gauss_name) list and a
        loop over it, and expanded the line set to 37 (SNR-relevant lines
        identified by comparing against Mapping_v_DAP.txt).  Fixed a
        naming bug
        (HeII_4685.68 was 'hii', now 'heii') and renamed two existing
        lines for a consistent doublet convention (hei->hei_b,
        caii_7291->caii_a) alongside their new partners (hei_a, caii_b).
        See docs/source/dap.rst for the full name table and the naming
        rules (species+ion lowercased; _a/_b for doublets, shorter
        wavelength first; a wavelength suffix otherwise).

'''

# # Develop tools for visualising the DAP results (using Vela as an example)



from astropy.io import ascii,fits
from astropy.table import Table, join
import matplotlib.pyplot as plt
import numpy as np


import re


def _usage_from_doc(doc):
    '''
    __doc__ truncated just before a line consisting of "History:" (or
    "History::"/"Version History" -- whitespace/colon-insensitive), so
    -h stays short even as that section grows -- without hand-
    duplicating the Synopsis/Options text in a second string.  Anchored
    to a whole line (not a bare substring search) so it can't misfire on
    "History:" appearing mid-sentence, and returns doc unchanged if no
    such line is present.
    '''
    m = re.search(r'^\s*(?:Version\s+)?History:{0,2}\s*$', doc, re.MULTILINE)
    return doc[:m.start()].rstrip() + '\n' if m else doc


def get_lines(colnames):
    names=[]
    for one in colnames:
        if one.count('e_flux'):
            names.append(one.replace('e_flux_',''))
    return names


def eval_sig2noise(xtab,sig=3):
    colnames=xtab.colnames
    lines=get_lines(colnames)
    for one_line in lines:
        flux=xtab['flux_%s' % one_line] 
        eflux=xtab['e_flux_%s' % one_line]
        sn=flux/eflux
        num=np.sum(sn>sig)
        print('%15s  %3d' % (one_line,num))

def get_one_line(xtab,name='Halpha_6562.85',xname='ha'):
    xline=xtab['id','flux_%s' % name,'e_flux_%s' % name,'vel_%s' % name,'e_vel_%s' % name,'disp_%s' % name,'e_disp_%s' % name]
    xline.rename_column('flux_%s' % name,'flux_%s' % xname)
    xline.rename_column('e_flux_%s' % name,'eflux_%s' % xname)
    xline['flux_%s' % xname]*=1e-16
    xline['eflux_%s' % xname]*=1e-16
    word=name.split('_')
    wave=eval(word[-1])
    xline.rename_column('vel_%s' % name,'vel_%s' % xname)
    xline.rename_column('e_vel_%s' % name,'e_vel_%s' % xname)

    xline['wave_%s' % xname] = wave*(1.+ xline['vel_%s' % xname]/2.997e5)
    xline['ewave_%s' % xname] = wave*(1+ xline['e_vel_%s' % xname]/2.997e5)

    xline.rename_column('disp_%s' % name,'fwhm_%s' % xname)
    xline['fwhm_%s' % xname]*=2.355
    xline.rename_column('e_disp_%s' % name,'efwhm_%s' % xname)
    xline['efwhm_%s' % xname]*=2.355
    return xline



def get_radec_fluxes(filename='DAP/dap-rsp108-sn20-00009083.dap.fits.gz'):
    '''
    Return information about he most prominent lines
    '''
    x=fits.open(filename)
    pt=Table(x['PT'].data)
    B=Table(x['NP_ELINES_B'].data)
    R=Table(x['NP_ELINES_R'].data)
    I=Table(x['NP_ELINES_I'].data)

    # (arm table, DAP_name from data/dap_lines.txt, output gauss_name),
    # sorted by wavelength within each arm.  See docs/source/dap.rst for
    # the naming convention (species+ion lowercased; _a/_b for doublets,
    # shorter wavelength first; a wavelength suffix instead when a line
    # isn't part of a fixed-ratio pair).
    lines=[
        (B,'[OII]_3726.03','oii_a'),
        (B,'[OII]_3728.82','oii_b'),
        (B,'[NeIII]_3868.75','neiii_a'),
        (B,'HeI_3888.65','hei_a'),
        (B,'HI_3889.05','hi_3889'),
        (B,'CaII_3933.66','caii_3933'),
        (B,'[NeIII]_3967.46','neiii_b'),
        (B,'Hepsilon_3970.07','hepsilon'),
        (B,'Hdelta_4101.77','hdelta'),
        (B,'[FeV]_4227.2','fev_4227'),
        (B,'Hgamma_4340.49','hgamma'),
        (B,'[OIII]_4363.21','oiii_4363'),
        (B,'[FeIII]_4658.1','feiii_a'),
        (B,'HeII_4685.68','heii'),
        (B,'Hbeta_4861.36','hb'),
        (B,'[OIII]_4958.91','oiii_a'),
        (B,'[FeIII]_4985.9','feiii_b'),
        (B,'[OIII]_5006.84','oiii_b'),
        (B,'[NI]_5197.9','ni_a'),
        (B,'[NI]_5200.26','ni_b'),
        (R,'HeI_5876.0','hei_b'),
        (R,'[OI]_6300.3','oi_a'),
        (R,'[OI]_6363.78','oi_b'),
        (R,'[FeX]_6374.51','fex_6374'),
        (R,'[NII]_6548.05','nii_a'),
        (R,'Halpha_6562.85','ha'),
        (R,'[NII]_6583.45','nii_b'),
        (R,'[SII]_6716.44','sii_a'),
        (R,'[SII]_6730.82','sii_b'),
        (R,'[CaII]_7291.46','caii_a'),
        (R,'[OII]_7318.92','oii_7320'),
        (R,'[CaII]_7323.88','caii_b'),
        (R,'[NiII]_7377.83','ni_ii_7377'),
        (I,'[FeII]_8616.96','feii_a'),
        (I,'[FeII]_8891.88','feii_b'),
        (I,'[SIII]_9069.0','siii_a'),
        (I,'[SIII]_9531.1','siii_b'),
    ]

    for arm_tab,dap_name,gauss_name in lines:
        try:
            xline=get_one_line(arm_tab,dap_name,gauss_name)
            pt=join(pt,xline)
        except:
            print('Could not get %s -> %s' % (dap_name,gauss_name))


    # OK that this end now wrap up
    for one_name in pt.colnames:
        if one_name.count('flux'):
            pt[one_name].format='.3e'
        elif one_name.count('fwhm') or one_name.count('vel') or one_name.count('wave'):
            pt[one_name].format='.3f'

    pt['ra'].format='.5f'
    pt['dec'].format='.5f'


    outname='DAPsum_test.txt'
    words=filename.split('/')
    root=words[-1]
    root=root.replace('.dap.fits.gz','')
    root=root.replace('dap-rsp108-sn20-','')
    outname='DAPsum_%s.txt' % root
    pt.write(outname,format='ascii.fixed_width_two_line',overwrite=True)

    return pt


def steer(argv):
    files=[]
    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(_usage_from_doc(__doc__))
            return
        elif argv[i][0]=='-':
            print('Error: Could not interet command line',argv)
        elif argv[i].count('fits'):
            files.append(argv[i])
        i+=1
    if len(files)==0:
        print('Apperently nothiing to do')

    for one_file in files:
        spec_sum=get_radec_fluxes(filename=one_file)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
