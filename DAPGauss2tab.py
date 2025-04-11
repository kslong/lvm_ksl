#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:  

Read one or more DAP output files and produce a 
table in the format of various other routines
that contains the fluxes etc of the lines
that were fit with Gaussians in the DAP


Command line usage (if any):

    usage: DAPGauss2tab.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

250301 ksl Coding begun

'''

import sys
from astropy.io import fits
from astropy.table import Table, join
import numpy as np
import matplotlib.pyplot as plt



def get_one_line(ptab,wavelength=3726.03,name='oii_a'):
    '''
    Return a table containging the data from one line, with column names that are those
    that I commonly use
    '''
    ztab=ptab[ptab['wl']==wavelength]
    print (len(ztab))
    mtab=ztab['id','id_fib']
    mtab['flux_%s' % name]=ztab['flux']
    mtab['eflux_%s' % name]=ztab['e_flux']
    mtab['wave_%s' % name]=ztab['wl'] * (1.+ztab['vel']/2.997e5)
    mtab['ewave_%s' % name]=ztab['wl'] * (1.+ztab['e_vel']/2.997e5)    
    mtab['fwhm_%s' % name] = ztab['disp']* 2.355
    mtab['efwhm_%s' % name] = ztab['e_disp']* 2.355
    
    return mtab




def make_primary_emission_line_table(xfile='dap-rsp108-sn20-00003470.dap.fits.gz'):
    try:
        x=fits.open(xfile)
    except:
        print('Error could not open :',xfile)
        return

    pt=Table(x['PT'].data)
    xtab=Table(x['PM_ELINES'].data)
    print(np.unique(xtab['wl']))
    pt.info()

    o2a=get_one_line(ptab=xtab,wavelength=3726.03,name='oii_a')
    o2a.info()
    prime=join(pt,o2a,join_type='left')

    o2b=get_one_line(ptab=xtab,wavelength=3728.82,name='oii_b')
    prime=join(prime,o2b,join_type='left')

    hb=get_one_line(ptab=xtab,wavelength=4861.36,name='hb')
    prime=join(prime,hb,join_type='left')
    
    o3a=get_one_line(ptab=xtab,wavelength=4958.91,name='oiii_a')
    prime=join(prime,o3a,join_type='left')

    o3b=get_one_line(ptab=xtab,wavelength=5006.84,name='oiii_b')
    prime=join(prime,o3b,join_type='left')

    o1a=get_one_line(ptab=xtab,wavelength=6300.3,name='oi_a')
    prime=join(prime,o1a,join_type='left')

    n2a=get_one_line(ptab=xtab,wavelength=6548.05,name='nii_a')
    prime=join(prime,n2a,join_type='left')

    ha=get_one_line(ptab=xtab,wavelength=6562.85,name='ha')
    prime=join(prime,ha,join_type='left')
    
    n2b=get_one_line(ptab=xtab,wavelength=6583.45,name='nii_b')
    prime=join(prime,n2b,join_type='left')    


    s2a=get_one_line(ptab=xtab,wavelength=6716.44,name='sii_a')
    prime=join(prime,s2a,join_type='left')

    s2b=get_one_line(ptab=xtab,wavelength=6730.82,name='sii_b')
    prime=join(prime,s2b,join_type='left')    

    s3a=get_one_line(ptab=xtab,wavelength=9069.0,name='siii_a')
    prime=join(prime,s3a,join_type='left')

    s3b=get_one_line(ptab=xtab,wavelength=9531.1,name='siii_b')
    prime=join(prime,s3b,join_type='left')    

    # OK that this end now wrap up
    for one_name in prime.colnames:
        if one_name.count('flux'):
            prime[one_name]*=1e-16
            prime[one_name].format='.3e'
        elif one_name.count('fwhm') or one_name.count('vel') or one_name.count('wave'):
            prime[one_name].format='.3f'

    prime['ra'].format='.5f'
    prime['dec'].format='.5f'

    
    prime.info()
    words=xfile.split('/')
    root=words[-1]
    root=root.replace('.dap.fits.gz','')
    root=root.replace('dap-rsp108-sn20-','')
    outname='DAP_PME_%s.txt' % root
    prime.write(outname,format='ascii.fixed_width_two_line',overwrite=True)
    return prime



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
        spec_sum=make_primary_emission_line_table(one_file)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)



