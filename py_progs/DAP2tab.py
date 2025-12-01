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
                                       
History:

241224 ksl Coding begun

'''

# # Develop tools for visualising the DAP results (using Vela as an example)



from astropy.io import ascii,fits
from astropy.table import Table, join
import matplotlib.pyplot as plt
import numpy as np

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

    dap_name='[OII]_3726.03'
    gauss_name='oii_a'
    # xline=get_one_line(B,dap_name,gauss_name)
    try:
        xline=get_one_line(B,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))

    dap_name='[OII]_3728.82'
    gauss_name='oii_b'
    
    try:
        xline=get_one_line(B,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))

    dap_name='Hgamma_4340.49'
    gauss_name='hgamma'
    
    try:
        xline=get_one_line(B,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))


    dap_name='HeII_4685.68'
    gauss_name='hii'
    
    try:
        xline=get_one_line(B,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))

    dap_name='Hbeta_4861.36'
    gauss_name='hb'
    
    try:
        xline=get_one_line(B,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))

    dap_name='[OIII]_4958.91'
    gauss_name='oiii_a'
    
    try:
        xline=get_one_line(B,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))


    dap_name='[OIII]_5006.84'
    gauss_name='oiii_b'
    
    try:
        xline=get_one_line(B,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))

    # Switch to R

    R=Table(x['NP_ELINES_R'].data)

    dap_name='HeI_5876.0'
    gauss_name='hei'
    
    try:
        xline=get_one_line(R,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))


    dap_name='[OI]_6300.3'
    gauss_name='oi_a'
    
    try:
        xline=get_one_line(R,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))

    dap_name='[NII]_6548.05'
    gauss_name='nii_a'
        
    
    try:
        xline=get_one_line(R,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))

    dap_name='Halpha_6562.85'
    gauss_name='ha'
        
    
    try:
        xline=get_one_line(R,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))

    dap_name='[NII]_6583.45' 
    gauss_name='nii_b'

    
    try:
        xline=get_one_line(R,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))

    dap_name='[SII]_6716.44'
    gauss_name='sii_a'
    
    
    try:
        xline=get_one_line(R,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))

    dap_name='[SII]_6730.82'
    gauss_name='sii_b'
    
    try:
        xline=get_one_line(R,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))

    dap_name='[CaII]_7291.46'
    gauss_name='caii_7291'
    
    try:
        xline=get_one_line(R,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))



    dap_name='[OII]_7318.92'
    gauss_name='oii_7320'

    
    try:
        xline=get_one_line(R,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))



    # Switch to I 

    I=Table(x['NP_ELINES_I'].data)

    dap_name='[SIII]_9069.0' 
    gauss_name='siii_a'
    
    try:
        xline=get_one_line(I,dap_name,gauss_name)
        pt=join(pt,xline)
    except:
        print('Could not get %s -> %s' % (dap_name,gauss_name))

    
    dap_name='[SIII]_9531.1'
    gauss_name='siii_b'

    
    try:
        xline=get_one_line(I,dap_name,gauss_name)
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
    
    # Now do the lines that are in the R changell
    return pt


def steer(argv):
    files=[]
    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
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
