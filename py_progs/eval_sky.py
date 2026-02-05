#!/usr/bin/env python
# coding: utf-8



'''
                    Space Telescope Science Institute

Synopsis:  

Create a simple plot to evaluate how well
sky subtraction is working on individual
SFrame files


Command line usage (if any):

    usage: eval_sky.py filename1 filename2 

Description:  

Primary routines:

    doit

Notes:

    This does not work on the various 
    median type spectra of multiple
    observations that ksl has written
                                       
History:

240318 ksl Coding begun

'''

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

def scifib(xtab,select='all',telescope=''):
    # print(np.unique(xtab['fibstatus']))
    # print(np.unique(xtab['targettype']))
    ztab=xtab[xtab['fibstatus']==0]
    if select=='all':
        ztab=ztab[ztab['targettype']!='standard']
    else:
        ztab=ztab[ztab['targettype']==select]

    if telescope!='' and telescope!='all':
        ztab=ztab[ztab['telescope']==telescope]

    print('Found %d fibers' % len(ztab))
    return ztab
    
    


def get_masked(filename='/Users/long/Projects/lvm_science/N103b/XCFrame-00007755.fits',ext='ERROR',select='science',telescope=''):
    x=fits.open(filename)
    xtab=Table(x['SLITMAP'].data)
    xgood=scifib(xtab,select,telescope)
    wav=x['WAVE'].data
    fwhm=x[ext].data[xgood['fiberid']-1]
    print('ok',x['MASK'].data.shape)
    # wav=np.ma.masked_array(wav,x['MASK'].data[xgood['fiberid']-1])
    foo=x['MASK'].data[xgood['fiberid']-1]
    print('mask',np.unique(foo,return_counts=True))
    fwhm=np.ma.masked_array(fwhm,x['MASK'].data[xgood['fiberid']-1])
    ave_fwhm=np.ma.average(fwhm,axis=0)
    min_fwhm=np.ma.min(fwhm,axis=0)
    max_fwhm=np.ma.max(fwhm,axis=0)
    med_fwhm=np.ma.median(fwhm,axis=0)
    plt.figure(2,(8,6))
    plt.clf()
    plt.plot(wav,ave_fwhm,label='Ave')
    plt.plot(wav,med_fwhm,label='Med')
    plt.plot(wav,min_fwhm,label='Min')
    plt.plot(wav,max_fwhm,label='Max')
    plt.ylim(-1e-13,1e-13)
    plt.legend()
    plt.xlabel('Wavelength',size=16)
    plt.ylabel(ext,size=16)





def eval_qual(filename='/Users/long/Projects/lvm_science/N103b/XCFrame-00007755.fits',ymin=-0.2e-13,ymax=1e-13,xmin=3600,xmax=9500,outroot=''):
    '''
    Provide a standard plot for looking at how well the sky subtraction has worked overall
    '''
    
    
    
    try:
        x=fits.open(filename)
    except:
        print('Error: eval_qual: Could not open %s' % filename)
        return

    header=x[0].header
    mjd=header['MJD']
    exposure=header['EXPOSURE']
    
    xtab=Table(x['SLITMAP'].data)

    plt.figure(1,(12,16))
    plt.clf()

    # science
    xgood=scifib(xtab,select='science',telescope='all')
    wav=x['WAVE'].data
    flux=x['FLUX'].data[xgood['fiberid']-1]
    sky=x['SKY'].data[xgood['fiberid']-1]
    tot=flux+sky    
    
    tot=np.ma.masked_array(tot,x['MASK'].data[xgood['fiberid']-1])
    flux=np.ma.masked_array(flux,x['MASK'].data[xgood['fiberid']-1])
    sky=np.ma.masked_array(sky,x['MASK'].data[xgood['fiberid']-1])

    med_tot=np.ma.median(tot,axis=0)
    med_flux=np.ma.median(flux,axis=0)
    med_sky=np.ma.median(sky,axis=0)

    plt.subplot(4,1,1)

    # plt.plot(wav,med_tot,label='Flux+Sky')

    xmed=np.median(med_flux)
    ymin=ymin+xmed
    ymax=ymax+xmed

    #plt.plot(wav,med_sky,label='Sky')
    plt.plot(wav,med_flux,'k',label='Flux')
    plt.plot([xmin,xmax],[5.9e-15+xmed,5.9e-15+xmed],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    plt.plot([xmin,xmax],[-5.9e-15+xmed,-5.9e-15+xmed],':r')
    
    plt.ylim(ymin,ymax)
    plt.xlim(xmin,xmax)
    plt.legend()
    plt.xlabel('Wavelength',size=16)
    plt.ylabel('Median',size=16)  
    plt.title('Science for MJD %d exposure %d' %(mjd,exposure))
    plt.tight_layout()
    
    plt.subplot(4,1,2)

    plt.plot(wav,med_tot,label='Flux+Sky')

    plt.plot(wav,med_sky,label='Sky')
    plt.plot(wav,med_flux,'k',label='Flux')
    
    plt.ylim(ymin,ymax)
    plt.xlim(xmin,xmax)
    plt.legend()
    plt.xlabel('Wavelength',size=16)
    plt.ylabel('Median',size=16)  
    plt.title('Science')
    plt.tight_layout()

    plt.subplot(4,1,3)
    
    xgood=scifib(xtab,select='SKY',telescope='SkyE')
    wav=x['WAVE'].data
    flux=x['FLUX'].data[xgood['fiberid']-1]
    sky=x['SKY'].data[xgood['fiberid']-1]
    tot=flux+sky    
    
    tot=np.ma.masked_array(tot,x['MASK'].data[xgood['fiberid']-1])
    flux=np.ma.masked_array(flux,x['MASK'].data[xgood['fiberid']-1])
    sky=np.ma.masked_array(sky,x['MASK'].data[xgood['fiberid']-1])

    med_tot=np.ma.median(tot,axis=0)
    med_flux=np.ma.median(flux,axis=0)
    med_sky=np.ma.median(sky,axis=0)
    
 

    plt.plot(wav,med_tot,label='Flux+Sky')

    plt.plot(wav,med_sky,label='Sky')
    plt.plot(wav,med_flux,'k',label='Flux')
    
    plt.ylim(ymin,ymax)
    plt.xlim(xmin,xmax)
    plt.legend()
    plt.xlabel('Wavelength',size=16)
    plt.ylabel('Median',size=16)  
    plt.title('SkyE')
    plt.tight_layout()

    plt.subplot(4,1,4)
    
    xgood=scifib(xtab,select='SKY',telescope='SkyW')
    wav=x['WAVE'].data
    flux=x['FLUX'].data[xgood['fiberid']-1]
    sky=x['SKY'].data[xgood['fiberid']-1]
    tot=flux+sky    
    
    tot=np.ma.masked_array(tot,x['MASK'].data[xgood['fiberid']-1])
    flux=np.ma.masked_array(flux,x['MASK'].data[xgood['fiberid']-1])
    sky=np.ma.masked_array(sky,x['MASK'].data[xgood['fiberid']-1])

    med_tot=np.ma.median(tot,axis=0)
    med_flux=np.ma.median(flux,axis=0)
    med_sky=np.ma.median(sky,axis=0)
    
 

    plt.plot(wav,med_tot,label='Flux+Sky')

    plt.plot(wav,med_sky,label='Sky')
    plt.plot(wav,med_flux,'k',label='Flux')
    
    plt.ylim(ymin,ymax)
    plt.xlim(xmin,xmax)
    plt.legend()
    plt.xlabel('Wavelength',size=16)
    plt.ylabel('Median',size=16)  
    plt.title('SkyW')
    plt.tight_layout()

    if outroot=='':
        word=filename.split('/')
        outroot=word[-1].replace('.fits','')
        
    plt.savefig('sky_%s.png' % outroot)



def steer(argv):

    files=[]

    i=1
    while i<len(argv):
        if argv[i].count('-h'):
            print(__doc__)
            return
        elif argv[i][0]=='-':
            print('Error: could not process command line: ',argv)
        else:
            files.append(argv[i])
        i+=1

    for one in files:
        eval_qual(one)





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)


