#!/usr/bin/env python
# coding: utf-8



'''
                    Space Telescope Science Institute

Synopsis:  

Create a simple plot to evaluate how well
sky subraction is working


Command line usage (if any):

    usage: eval_sky.py filename1 filename2 

Description:  

Primary routines:

    doit

Notes:
                                       
History:

240318 ksl Coding begun

'''

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.gridspec import GridSpec

def scifib(xtab,select='all',telescope=''):
    '''
    Select good fibers from a telescope, of a spefic
    type or all from a telescope from the slitmap table
    of a calbrated file
    '''
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



def limit_spectrum(wave,flux,wmin,wmax):
    '''
    Get a section of the spectrum
    '''

    f=flux[wave>wmin]
    w=wave[wave>wmin]

    f=f[w<wmax]
    w=w[w<wmax]

    return w,f

def get_yscale(f,ymin,ymax):
    med=np.median(f)
    zmin=ymin+med
    zmax=ymax+med
    print(zmin,zmax,med,ymin,ymax)
    return zmin,zmax

def eval_qual(filename='/Users/long/Projects/lvm_science/N103b/XCFrame-00007755.fits',ymin=-0.2e-13,ymax=1e-13,xmin=3600,xmax=9500,outroot=''):
    '''
    Provide a standard plot for looking at how well the sky subtraction has worked overall
    '''
    
    xtype='SkySubtracted'
    if filename.count('CFRAME'):
        print('This is an LCFrame')
        xtype='NotSkySubtracted'
    
    
    try:
        x=fits.open(filename)
    except:
        print('Error: eval_qual: Could not open %s' % filename)
        return

    header=x[0].header
    mjd=header['MJD']
    exposure=header['EXPOSURE']
    
    xtab=Table(x['SLITMAP'].data)

    science_fibers=scifib(xtab,select='science',telescope='Sci')
    skye_fibers=scifib(xtab,select='SKY',telescope='SkyE')
    skyw_fibers=scifib(xtab,select='SKY',telescope='SkyW')

    wav=x['WAVE'].data
    sci_flux=x['FLUX'].data[science_fibers['fiberid']-1]
    sci_sky=x['SKY'].data[science_fibers['fiberid']-1]
    sci_mask=x['MASK'].data[science_fibers['fiberid']-1]
    sci_flux=np.ma.masked_array(sci_flux,sci_mask)
    sci_sky=np.ma.masked_array(sci_sky,sci_mask)
    sci_flux_med=np.ma.median(sci_flux,axis=0)
    sci_sky_med=np.ma.median(sci_sky,axis=0)



    skye_flux= x['FLUX'].data[skye_fibers['fiberid']-1]
    skye_sky=x['SKY'].data[skye_fibers['fiberid']-1]
    skye_mask=x['MASK'].data[skye_fibers['fiberid']-1]
    skye_flux=np.ma.masked_array(skye_flux,skye_mask)
    skye_sky=np.ma.masked_array(skye_sky,skye_mask)
    skye_flux_med=np.ma.median(skye_flux,axis=0)
    skye_sky_med=np.ma.median(skye_sky,axis=0)

    skyw_flux= x['FLUX'].data[skyw_fibers['fiberid']-1]
    skyw_sky=x['SKY'].data[skyw_fibers['fiberid']-1]
    skyw_mask=x['MASK'].data[skyw_fibers['fiberid']-1]
    skyw_flux=np.ma.masked_array(skyw_flux,skyw_mask)
    skyw_sky=np.ma.masked_array(skyw_sky,skyw_mask)
    skyw_flux_med=np.ma.median(skyw_flux,axis=0)
    skyw_sky_med=np.ma.median(skyw_sky,axis=0)

    fig=plt.figure(1,(8,12))
    plt.clf()
    gs= GridSpec(3, 3, figure=fig)

    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(wav,sci_flux_med,label='Sky-Subtracted Science',zorder=2)
    ax1.plot(wav,skye_flux_med,label='SkyE-Subtracted SkyE',zorder=1)
    ax1.plot(wav,skyw_flux_med,label='SkyW-Subtracted SkyW',zorder=0)
    ax1.plot([3600,9600],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    ax1.plot([3600,9600],[-5.9e-15,-5.9e-15],':r')
    ax1.set_xlim(3600,9600)
    ymin,ymax=ax1.get_ylim()
    ax1.set_ylim(-1e-14,ymax)
    ax1.legend()

    ax2 = fig.add_subplot(gs[1, :])
    ax2.semilogy(wav,sci_flux_med+sci_sky_med,label='Science Total',zorder=2)
    ax2.semilogy(wav,sci_sky_med,label='Science Sky',zorder=1)
    ymin,ymax=plt.ylim()
    ymax=np.max(sci_flux_med+sci_sky_med)
    ax2.set_ylim(1e-3*ymax,1.1*ymax)
    ax2.set_xlim(3600,9600)
    ax2.legend()


    ax3 = fig.add_subplot(gs[2, 0])
    wmin=4650
    wmax=5100

    xwav,xsci_flux_med=limit_spectrum(wav,sci_flux_med,wmin,wmax)
    ax3.plot(xwav,xsci_flux_med,label='Sky-Subtracted Science',zorder=2)
    # ax5.plot(wav,sci_flux_med,label='Sky-Subtracted Science',zorder=2)

    xwav,xskye_flux_med=limit_spectrum(wav,skye_flux_med,wmin,wmax)
    ax3.plot(xwav,xskye_flux_med,label='SkyE-Subtracted SkyE',zorder=1)

    xwav,xskyw_flux_med=limit_spectrum(wav,skyw_flux_med,wmin,wmax)
    ax3.plot(xwav,xskyw_flux_med,label='SkyW-Subtracted SkyW',zorder=0)
    ax3.plot([wmin,wmax],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    ax3.plot([wmin,wmax],[-5.9e-15,-5.9e-15],':r')
    ax3.set_xlim(wmin,wmax)
    ymin,ymax=get_yscale(xsci_flux_med,-1e-14,1e-14)
    ax3.set_ylim(ymin,ymax)
    
    #ax3.plot(wav,sci_flux_med,label='Sky-Subtracted Science',zorder=2)
    #ax3.plot(wav,skye_flux_med,label='SkyE-Subtracted SkyE',zorder=1)
    #ax3.plot(wav,skyw_flux_med,label='SkyW-Subtracted SkyW',zorder=0)
    #ax3.plot([3600,9600],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    #ax3.plot([3600,9600],[-5.9e-15,-5.9e-15],':r')
    #ax3.set_xlim(4600,5200)
    # ax3.legend()



    ax4 = fig.add_subplot(gs[2, 1])
    wmin=6250
    wmax=6800

    xwav,xsci_flux_med=limit_spectrum(wav,sci_flux_med,wmin,wmax)
    ax4.plot(xwav,xsci_flux_med,label='Sky-Subtracted Science',zorder=2)
    # ax5.plot(wav,sci_flux_med,label='Sky-Subtracted Science',zorder=2)

    xwav,xskye_flux_med=limit_spectrum(wav,skye_flux_med,wmin,wmax)
    ax4.plot(xwav,xskye_flux_med,label='SkyE-Subtracted SkyE',zorder=1)

    xwav,xskyw_flux_med=limit_spectrum(wav,skyw_flux_med,wmin,wmax)
    ax4.plot(xwav,xskyw_flux_med,label='SkyW-Subtracted SkyW',zorder=0)
    ax4.plot([wmin,wmax],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    ax4.plot([wmin,wmax],[-5.9e-15,-5.9e-15],':r')
    ax4.set_xlim(wmin,wmax)
    ymin,ymax=get_yscale(xsci_flux_med,-1e-14,1e-14)
    ax4.set_ylim(ymin,ymax)
    
    


    ax5 = fig.add_subplot(gs[2, 2])
    wmin=9450
    wmax=9600
    xwav,xsci_flux_med=limit_spectrum(wav,sci_flux_med,wmin,wmax)
    ax5.plot(xwav,xsci_flux_med,label='Sky-Subtracted Science',zorder=2)
    # ax5.plot(wav,sci_flux_med,label='Sky-Subtracted Science',zorder=2)

    xwav,xskye_flux_med=limit_spectrum(wav,skye_flux_med,wmin,wmax)
    ax5.plot(xwav,xskye_flux_med,label='SkyE-Subtracted SkyE',zorder=1)

    xwav,xskyw_flux_med=limit_spectrum(wav,skyw_flux_med,wmin,wmax)
    ax5.plot(xwav,xskyw_flux_med,label='SkyW-Subtracted SkyW',zorder=0)
    ax5.plot([wmin,wmax],[5.9e-15,5.9e-15],':r',label=r'$Med \pm$ MW 5 $\sigma$' )
    ax5.plot([wmin,wmax],[-5.9e-15,-5.9e-15],':r')
    ax5.set_xlim(wmin,wmax)
    ymin,ymax=get_yscale(xsci_flux_med,-1e-14,1e-14)
    ax5.set_ylim(ymin,ymax)
    # ax5.legend()

    #plt.plot(wav,sci_flux_med+sci_sky_med,label='Science Total')
    #plt.plot(wav,skye_flux_med+skye_sky_med,label='SkyE Total')
    #plt.plot(wav,skyw_flux_med+skyw_sky_med,label='SkyWe Total')
    #plt.legend()

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


