#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Parse the DAP output file, and produce an astropy 
table with the information needed to make simple plots


Command line usage (if any):

    usage: DAP2tab.py filename

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



def get_radec_fluxes(filename='DAP/dap-rsp108-sn20-00009083.dap.fits.gz'):
    '''
    Return information about he most prominent lines
    '''
    x=fits.open(filename)
    pt=Table(x['PT'].data)
    B=Table(x['NP_ELINES_B'].data)
    
    name='[OII]_3726.03'
    xline=B['id','flux_%s' % name,'e_flux_%s' % name]
    xname='o2a'
    xline.rename_column('flux_%s' % name,'flux_%s' % xname)
    xline.rename_column('e_flux_%s' % name,'e_flux_%s' % xname)    
    
    pt=join(pt,xline)
    
    name='[OII]_3728.82'

    xline=B['id','flux_%s' % name,'e_flux_%s' % name]
    xname='o2b'
    xline.rename_column('flux_%s' % name,'flux_%s' % xname)
    xline.rename_column('e_flux_%s' % name,'e_flux_%s' % xname)    
    
    pt=join(pt,xline)

    name='HeII_4685.68'
    xline=B['id','flux_%s' % name,'e_flux_%s' % name]
    xname='he2'
    xline.rename_column('flux_%s' % name,'flux_%s' % xname)
    xline.rename_column('e_flux_%s' % name,'e_flux_%s' % xname)    
    
    pt=join(pt,xline)

    name='Hbeta_4861.36'
    xline=B['id','flux_%s' % name,'e_flux_%s' % name]
    xname='hb'
    xline.rename_column('flux_%s' % name,'flux_%s' % xname)
    xline.rename_column('e_flux_%s' % name,'e_flux_%s' % xname)    
    
    pt=join(pt,xline)    

    name='[OIII]_5006.84'
    xline=B['id','flux_%s' % name,'e_flux_%s' % name]
    xname='o3'
    xline.rename_column('flux_%s' % name,'flux_%s' % xname)
    xline.rename_column('e_flux_%s' % name,'e_flux_%s' % xname)    
    pt=join(pt,xline)

    R=Table(x['NP_ELINES_R'].data)

    name='[OI]_6300.3'
    xline=R['id','flux_%s' % name,'e_flux_%s' % name]
    xname='o1'
    xline.rename_column('flux_%s' % name,'flux_%s' % xname)
    xline.rename_column('e_flux_%s' % name,'e_flux_%s' % xname)    
    pt=join(pt,xline)
    
    name='[NII]_6548.05'
    xline=R['id','flux_%s' % name,'e_flux_%s' % name]
    xname='n2a'
    xline.rename_column('flux_%s' % name,'flux_%s' % xname)
    xline.rename_column('e_flux_%s' % name,'e_flux_%s' % xname)    
    pt=join(pt,xline)
        
    name='Halpha_6562.85'
    xline=R['id','flux_%s' % name,'e_flux_%s' % name,'vel_%s' % name,'e_vel_%s' % name,'disp_%s' % name,'e_disp_%s' % name]
    xname='ha'
    xline.rename_column('flux_%s' % name,'flux_%s' % xname)
    xline.rename_column('e_flux_%s' % name,'e_flux_%s' % xname)
    xline.rename_column('vel_%s' % name,'vel_%s' % xname)
    xline.rename_column('e_vel_%s' % name,'e_vel_%s' % xname)
    xline.rename_column('disp_%s' % name,'disp_%s' % xname)
    xline.rename_column('e_disp_%s' % name,'e_disp_%s' % xname)
    
    pt=join(pt,xline)
        
    name='[NII]_6583.45' 
    xline=R['id','flux_%s' % name,'e_flux_%s' % name]
    xname='n2b'
    xline.rename_column('flux_%s' % name,'flux_%s' % xname)
    xline.rename_column('e_flux_%s' % name,'e_flux_%s' % xname)    
    pt=join(pt,xline)

    name='[SII]_6716.44'
    xline=R['id','flux_%s' % name,'e_flux_%s' % name]
    xname='s2a'
    xline.rename_column('flux_%s' % name,'flux_%s' % xname)
    xline.rename_column('e_flux_%s' % name,'e_flux_%s' % xname)    
    pt=join(pt,xline)
    
    name='[SII]_6730.82'
    xline=R['id','flux_%s' % name,'e_flux_%s' % name]
    xname='s2b'
    xline.rename_column('flux_%s' % name,'flux_%s' % xname)
    xline.rename_column('e_flux_%s' % name,'e_flux_%s' % xname)    
    pt=join(pt,xline)

    I=Table(x['NP_ELINES_I'].data)

    name='[SIII]_9069.0' 
    xline=I['id','flux_%s' % name,'e_flux_%s' % name]
    xname='s3a'
    xline.rename_column('flux_%s' % name,'flux_%s' % xname)
    xline.rename_column('e_flux_%s' % name,'e_flux_%s' % xname)    
    pt=join(pt,xline)
    
    name='[SIII]_9531.1'
    xline=I['id','flux_%s' % name,'e_flux_%s' % name]
    xname='s3b'
    xline.rename_column('flux_%s' % name,'flux_%s' % xname)
    xline.rename_column('e_flux_%s' % name,'e_flux_%s' % xname)    
    pt=join(pt,xline)

    for one_name in pt.colnames:
        if one_name.count('flux') or one_name.count('vel') or one_name.count('disp'):
            pt[one_name].format='.2f'

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
