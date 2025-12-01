#!/usr/bin/env python
# coding: utf-8


'''
                    Space Telescope Science Institute

Synopsis:  

Match RAs and Decs of observations made with the various talescopes to
positions in a catalog, usually the locations of sky fields


Command line usage (if any):

    usage: SkyID.py [-h] [-sep 60] [-out whatever ] Observations.txt  skyfields.txt

    where 

        Observation.txt is a file produced by the routine SummarizeData
        skyfields.txt is file containing the skyfield positions (or any file with the
            columns Source_name, RA and Dec

        -h prints this documentation an quits
        -sep whatever sets the maximum separation for a match
        -out whatever changes a portion of the names of the output file

Description:  

    This routine reads a file produced by the routine SummarizeData.py that contains
    the positions of all of the science and sky telescopes, and first produces
    a table Xpos that has the RAs and DECs of each of the telescopes (when they exist)
    for each observation.  A column in this table indicates what telescope the
    position applies to.

    Then the routine finds the closet object in the second table, asumming that
    object is within the agreed separation (in arcsec).

    The results of this are wirtten to a table (the XMatch) table.

    Finally, the routine looks at the Xmatch table and determines how many times
    the locations in the skyfield table thabe been observed.  


Primary routines:

    doit

Notes:

    Although designed to work with a table of the skyfields, the second file can containing
    any file with an RA and Dec column indicating the positions, e.g for a set of SNRs.
                                       
History:

240402 ksl Coding begun

'''


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.table import Table,vstack,hstack
import astropy.units as u



def get_sky(filename='All_data.txt'):
    '''
    Clean and reformat the observation file
    '''
    xtab=ascii.read(filename)
    xtab=xtab[xtab['Type']=='object']
    sci=xtab['RA','Dec','Exposure','MJD','Object']
    sci['Telescope']='Science'
    sky_e=xtab['RA_E','Dec_E','Exposure','MJD','Object']
    sky_e.rename_column('RA_E','RA')
    sky_e.rename_column('Dec_E','Dec')
    sky_e['Telescope']='SkyE'
    sky_w=xtab['RA_W','Dec_W','Exposure','MJD','Object']
    sky_w['Telescope']='SkyW'
    sky_w.rename_column('RA_W','RA')
    sky_w.rename_column('Dec_W','Dec')
    xall=vstack([sci,sky_e,sky_w])
    print(len(xall))
    xall=xall[xall['RA']>0.0]
    print(len(xall))
    return xall 


def find_closest_objects(table1_path, table2_path, max_sep=0.5):
    '''
    Find the objects with a given distance given two astropy tables.  The
    routine returns only the closest object that satisvies this criterion.
    '''
    # Read the two Astropy tables
    try:
        table1 = Table.read(table1_path,format='ascii.fixed_width_two_line')
    except:
        print('Error: find_closest_objects: could not read %s' % table1_path)
        return []

    try:
        table2 = Table.read(table2_path,format='ascii.fixed_width_two_line')
    except:
        print('Error: find_closest_objects: could not read %s' % table2_path)
        return []


    print('get_closest_objects: Beginning x-match of %s and %s' % (table1_path,table2_path))

    # Convert RA and Dec columns to SkyCoord objects
    coords_table1 = SkyCoord(ra=table1['RA'] * u.degree, dec=table1['Dec'] * u.degree)
    coords_table2 = SkyCoord(ra=table2['RA'] * u.degree, dec=table2['Dec'] * u.degree)

    # Find the closest objects in table2 to each object in table1
    closest_indices = []
    closest_separations = []

    for coord1 in coords_table1:
        separations = coord1.separation(coords_table2)
        min_index = np.argmin(separations)
        closest_indices.append(min_index)
        closest_separations.append(separations[min_index].to(u.arcsec).value)

    # Create a new table to store the closest objects and their separations
    table1['Sep']=closest_separations*u.arcsec
    xtab=hstack([table1,table2[closest_indices]])

    xtab=xtab[xtab['Sep']<max_sep]
    xtab['Sep'].format='.2f'

    print('Of %d objects in %s and %d objects in %s, found %d matches' % (len(table1),table1_path,len(table2),table2_path,len(xtab)))

    return xtab


def doit( obs_file, obj_file, max_sep=60,outroot=''):

    if outroot=='':
        obs_words=obs_file.split('/')
        obj_words=obj_file.split('/')
        posroot='%s' % (obs_words[-1].replace('.txt',''))
        outroot='%s_%s' % (obs_words[-1].replace('.txt',''),obj_words[-1].replace('.txt',''))
    else:
        posroot=outroot

    xsky=get_sky(filename=obs_file)
    good_name='XPos.%s.txt' % posroot
    xsky.write(good_name,format='ascii.fixed_width_two_line',overwrite=True)
    print('Wrote good telescope positions to %s' % good_name)

    xmatch=find_closest_objects(good_name,obj_file,max_sep)


    xmatch_name='Xmatch.%s.txt' %outroot
    xmatch.write(xmatch_name ,format='ascii.fixed_width_two_line',overwrite=True)
    print('Wrote xmatches  to %s' % xmatch_name)
    source,counts=np.unique(xmatch['Source_name'],return_counts=True)

    qtab=Table([source,counts],names=['Source_name','NumObs'])
    qtab.sort(['NumObs'],reverse=True)

    num_name='Xnum.%s.txt'% outroot
    qtab.write(num_name,format='ascii.fixed_width_two_line',overwrite=True)
    print('Wrote number counts to %s' % num_name)

    return


def steer(argv):
    '''
    This is primarily intended to gather inputs
    '''

    obs_file=''
    obj_file=''
    outroot=''
    sep=60.

    i=1
    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-sep':
            i+=1
            try:
                sep=eval(argv[i])
            except:
                print('Error: Could not parse separation of : ',argv[i])
                return
        elif argv[i][0:4]=='-out':
            i+=1
            outroot=argv[i]
        elif obs_file=='':
            obs_file=argv[i]
        elif obj_file=='':
            obj_file=argv[i]
        else:
            print('Error: Cannot parse command line: ',argv)
            return

        i+=1

    if obj_file=='':
        print('Error: Not enough arguments: ', argv)
        return

    doit(obs_file,obj_file,sep,outroot)
    return



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
