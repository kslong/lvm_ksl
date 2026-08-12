#!/usr/bin/env python 

'''

Synopsis:  
    This is a simple program to that makes a master file from the 
    a region file          

Description:  

    Usage:

    reg2master.py regionfile  [masterfile]

    The program requires 1 parameters:

        The name of the regionfile
    
    The program has 1 optional parameter which is the
    name of the master file.  If the optional parameter
    is missing the name will be the same as the regionfile
    with .txt attached


Notes:

    General-purpose region/table conversion utility used across several
    of the author's projects (not LVM-specific); vendored into lvm_ksl so
    the region-based spectral extraction workflow (see
    docs/source/region_spectra.rst) is reproducible without an external
    dependency. See also master2reg.py, the inverse conversion.

History::

    090109 ksl Coded
    111111 ksl Fixed so would work with pydocs. Routine might be
        rewritten to make the main routine smaller
    190507 ksl Modified to read the color in, because Bill often uses
        this to indicate something about the source. This is not
        currently incorporated into master2reg, because I am not sure
        exactly how I want to do this, and because I do not need it
        currently.
    241230 ksl Add box as a type of region. Note that for a box Theta
        corresponds to the angle E of N for the 'Minor' axis. This is
        also true for ellipses.
    260812 ksl Vendored into lvm_ksl/py_progs from the author's personal
        py_xray utilities (unchanged apart from this History entry and
        docstring RST formatting).

'''

import sys
import os
from astropy.table import Table
import numpy


def radec2deg(ra,dec):
    ''' Convert an ra dec string to degrees.  The string can already
    be in degrees in which case all that happens is a conversion to
    a float'''

    r=ra.split(':')
    d=dec.split(':')

    rr=float(r[0])
    if len(r)>1:
        rr=rr+float(r[1])/60.
    if len(r)>2:
        rr=rr+float(r[2])/3600.
    if len(r)>1:
        rr=15.*rr  # Since we assume ra was in hms
    
    dd=float(d[0])
    x=0
    if len(d)>1:
        x=x+float(d[1])/60.
    if len(d)>2:
        x=x+float(d[2])/3600.
    
    if dd > 0:
        dd=dd+x
    else:
        dd=dd-x
    
    return rr,dd  


def size2arcsec(word):
    '''
    Convert a string to arcsec if possible

    Otherwise return the value
    '''
    if word.count('"')==1: # We have arcsec
        value=float(word.rstrip('"'))
    elif word.count("'")==1: # We have arcsec
        value=60.*float(word.rstrip("'"))
    else:
        print('Error : Could not parse ',word,' to arcsec')
        value=float(word)
    return value

    


def read_regions(filename):
    '''
    read_regions(filename) reads and parses a region file

    170511  ksl Give sources names if they do not have one
    '''
    f=open(filename,'r')
    type='unknown'
    xcolor='unknown'

    records=[]
    source_no=1

    lines=f.readlines()
    for line in lines:
        line=line.replace('(',' ')
        line=line.replace(')',' ')
        line=line.replace('{',' ')
        line=line.replace('}',' ')
        line=line.replace('=',' ')
        line=line.replace(',',' ')
        line=line.split()
        if len(line)==0:
            continue

        # process a generic line
        # print(line)
 

        if line[0]=='global':
            # print('Global:',line)
            j=0
            while j<len(line)-1:
                if line[j]=='color':
                    xcolor=line[j+1]
                j+=1

        # find the name
        j=0
        name='unknown'
        color='unknown'
        while j<len(line)-1:
            if line[j]=='text':
                name=line[j+1]
            if line[j]=='color':
                color=line[j+1]
            j=j+1
        if name=='unknown':
            name='zzz%03d' % source_no
        if color=='unknown':
            color=xcolor



        # Now process circles
        if len(line)==0:
            pass
        elif line[0]=='fk5':
            type='fk5'
        elif line[0]=='image':
            type='image'
        elif line[0]=='physical':
            type='physical'
        elif line[0]=='circle':
            if type=='unknown' or type=='fk5':
                rr,dd=radec2deg(line[1],line[2])
            else:
                rr=float(line[1])
                dd=float(line[2])
            x1=size2arcsec(line[3])
            record=[name,rr,dd,'circle',x1,0.0,0.0,color]
            records=records+[record]
            source_no+=1
        elif line[0]=='ellipse':
            if type=='unknown' or type=='fk5':
                rr,dd=radec2deg(line[1],line[2])
            else:
                rr=float(line[1])
                dd=float(line[2])
            x1=size2arcsec(line[3])
            x2=size2arcsec(line[4])
            theta=float(line[5])
            record=[name,rr,dd,'ellipse',x1,x2,theta,color]
            records=records+[record]
            source_no+=1
        elif line[0]=='box':
            if type=='unknown' or type=='fk5':
                rr,dd=radec2deg(line[1],line[2])
            else:
                rr=float(line[1])
                dd=float(line[2])
            x1=size2arcsec(line[3])
            x2=size2arcsec(line[4])
            theta=float(line[5])
            record=[name,rr,dd,'box',x1,x2,theta,color]
            records=records+[record]
            source_no+=1
        elif line[0]=='annulus':
            if type=='unknown' or type=='fk5':
                rr,dd=radec2deg(line[1],line[2])
            else:
                rr=float(line[1])
                dd=float(line[2])
            x1=size2arcsec(line[3])
            x2=size2arcsec(line[4])
            theta=0.0
            record=[name,rr,dd,'annulus',x1,x2,theta,color]
            records=records+[record]
            source_no+=1
        else:
            print('Did not use :',line)
        # print(line)
    f.close()
    return type,records




def write_masterfile(masterfile,records,source='unknown',type='unknown'):
    '''
    write a master file

    Notes::

        141102 ksl Set so that first line was written in a better
            format for astropy
    '''

    names=['Source_name','RA','Dec','RegType','Major','Minor','Theta','Color']

    x=Table()

    i=0
    while i<len(names):
        value=[]
        for record in records:
            # print(record)
            value.append(record[i])
        x[names[i]]=value
        i+=1

    x['RA'].format='10.6f'
    x['Dec'].format='10.6f'
    x['Major'].format='8.2f'
    x['Minor'].format='8.2f'
    x['Theta'].format='8.2f'
    x['Color'].format='10s'

    x.write(masterfile,format='ascii.fixed_width_two_line',overwrite=True)

    return

def doit(regionfile,masterfile=''):
    '''
    Create a master table from a region file
    '''

    if masterfile=='':
        masterfile=regionfile+'.txt'
    elif masterfile[-4:].count('.txt')==0:
        masterfile=masterfile+'.txt'
    

    print('Regionfile          :', regionfile)
    print('Masterfile          :', masterfile)

    # Finished getting input information

    # read the regionfile

    type,records=read_regions(regionfile)

    print('No of records parsed ', len(records))

    write_masterfile(masterfile,records,regionfile,type)


def steer(argv):
    regionfile=''
    root=''
    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][0]=='':
            print('Could not parse command line: ',argv)
            return
        elif regionfile=='':
            regionfile=argv[i]
        elif root=='':
            root=argv[i]
        else:
            print('Too many arguments in command line: ',argv[i])
            return

        i+=1

    doit(regionfile,masterfile=root)






# Subroutines and function calls are above.  The driving routine is below
if __name__ == "__main__":
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print(__doc__)


