#!/usr/bin/env python

"""

Synopsis:  
    This is a simple program to that makes a region file from the 
    so-called master file.  At present, it can accepts properly
    only circles or ellipses, or boxes

Description:  

    Usage:

    master2reg.py [-h] [-r 4.0] [-color red] masterfile  [regionfile]

    The program requires 1 parameters:

        The name of the masterfile
    
    The program has several optional parameters  including the
    name of the region file.  If the optional parameter
    is missing the name will be the same as the masterfile
    with .reg attached

    The masterfile is an ascii file which should consist of a

    source_name ra dec regiontype one or more numbers 

    or

    source_man ra dec

    or and astropy table in which case one needs columns named

    Source_name
    RA
    Dec
    RegType
    Major
    Minor
    Theta

    If only the first three of these are supplied as for example in
    a coin file then, the RegType is assumed to be a circle and the 
    defaults are applied


    
    The numbers depend on the region type

    for a circle there is only a radius (which is assumed to be in arcsec)

    for an ellipse, one has a major and minor axis and an angle (in degrees)

    -color    allows one to set the color of the region files
    -r    allows one to set the radius of the circular regions, when one 
        provides only a name and position
    -h    prints this information and quits


Notes:

    The routine should also process the output file of ae.py, e.g.
    all.collated.txt.  130830 - I don't beleive the flag for all.collated
    is needed anymore.

    The case where a circular apertue is used but a minor axis is not
    supplied is probably not handled correctly

    General-purpose region/table conversion utility used across several
    of the author's projects (not LVM-specific); vendored into lvm_ksl so
    the region-based spectral extraction workflow (see
    docs/source/region_spectra.rst) is reproducible without an external
    dependency. See also reg2master.py, the inverse conversion.

History::

    090109 ksl Coded
    111111 ksl Modified so pydocs would work. Routine should be
        rewritten to move most of this material out of main
    111216 ksl Modified so that if it cannot interpret the last few
        columns, that it still produces a region file based on the
        first 3 columns
    121212 ksl Modified to increase the input options, notably color
    130830 ksl Modified to allow one to specify a radius for the region
        file in instances where it is not defined.
    150220 ksl Modified so uses astropy to read the masterfile
    251210 ksl Modified for situation where the input file is in fits
        format, and also to allow for the possibility that we have a
        Color, but not necessarily a size
    260812 ksl Vendored into lvm_ksl/py_progs from the author's personal
        py_xray utilities (unchanged apart from this History entry and
        docstring RST formatting).

"""

import sys
import os
from astropy.io import ascii
import numpy
from astropy.table import Table


def radec2deg(ra, dec):
    """ Convert an ra dec string to degrees.  The string can already
    be in degrees in which case all that happens is a conversion to
    a float

    If what is transferred is a float, the routine assumes it has been
    given ra and dec in degrees and just returns ra,dec
    
    """

    try:
        r = ra.split(":")
        d = dec.split(":")
    except AttributeError:
        return ra, dec

    rr = float(r[0])
    if len(r) > 1:
        rr = rr + float(r[1]) / 60.0
    if len(r) > 2:
        rr = rr + float(r[2]) / 3600.0
    if len(r) > 1:
        rr = 15.0 * rr  # Since we assume ra was in hms

    dd = float(d[0])
    x = 0
    if len(d) > 1:
        x = x + float(d[1]) / 60.0
    if len(d) > 2:
        x = x + float(d[2]) / 3600.0

    if dd > 0:
        dd = dd + x
    else:
        dd = dd - x

    return rr, dd


def read_masterfile(filename, xtype="circle", xmajor=3, xminor=3, xtheta=0.0):
    """
    Read the masterfile

    History::

        100705 ksl Modified so would accept a file that just contained
            three columns, a name, and ra and a dec
        111113 ksl Comment - While this works it looks like I only
            partially completed an upgrade of this to make it more
            general
        121212 ksl Modified so will also read a all.collated.txt file
            and produce a region file form this
        141102 ksl Modified so reads master file that can be read by
            astropy.io.ascii. Note that I did not do this by switching
            to astropy.io which would have been a better option. This
            is very much a kluge now.
        150220 ksl Modified to use astropy
    """

    try:
        data = ascii.read(filename)
    except:
        try:
            data=Table.read(filename)
        except:
            print("Failed to read %s", filename)
            return 0

    colnames = data.colnames

    print(colnames)

    # Check whether the column names are defined or whether this a generic file

    if colnames[0] == "col1":
        if len(colnames) == 2:
            # Assume we have only a list of RA and Decs, in which
            # case we need to generate the source names
            data.rename_column("col1", "RA")
            data.rename_column("col2", "Dec")
            x = []
            i = 0
            while i < len(data["RA"]):
                x.append("x%03d" % i)
                i += 1
            data["Source_name"] = x
            # Add a check for nan's which cropped up in tweakreg files
            i = 0
            nan_row = []
            while i < len(data):
                if numpy.isnan(data["RA"][i]) or numpy.isnan(data["Dec"][i]):
                    print("gotcha %d" % i)
                    nan_row.append(i)
                i += 1
            if len(nan_row):
                data.remove_rows(nan_row)
        else:
            data.rename_column("col1", "Source_name")
            data.rename_column("col2", "RA")
            data.rename_column("col3", "Dec")
        if len(colnames) > 3:
            try:
                data.rename_column("col4", "RegType")
                data.rename_column("col5", "Major")
                data.rename_column("col6", "Minor")
                data.rename_column("col7", "Theta")
                complete = True
            except ValueError:
                complete = False
        else:
            complete = False
    else:  # This is an astropy table
        # Check that the essential files exist
        print("This is an astropy table with headers")
        ok = True
        if ("Source_name" in colnames) == False:
            print("There is no column named Source_name, so manufacturing one")
            names = []
            i = 0
            while i < len(data):
                names.append("x%03d" % (i + 1))
                i += 1
            data["Source_name"] = names
            # ok=False
        if ("RA" in colnames) == False:
            ok = False
        if ("Dec" in colnames) == False:
            ok = False
        if ok == False:
            print("Error: File read but column names are not correct:", colnames)
            print("       Minimally need Source_name,RA,Dec")
            return
        ok = True
        # print('test',colnames)
        # print('test','RegType' in colnames)

        if ("RegType" in colnames) == False:
            print("Did not get RegType")
            ok = False
        if ("Major" in colnames) == False:
            print("Did not get Major")
            ok = False
        if ("Minor" in colnames) == False:
            print("Did not get Minor")
            ok = False
        if ("Theta" in colnames) == False:
            print("Did not get Theta")
            ok = False

        if ("Color" in colnames) == False:
            print("Did not get  Color")
            data["Color"] = "Unknown"

        print("OK", ok)
        if ok == False:
            print("Colnames in table:\n", colnames)
            print(
                "Warning: Creating region file but with generic values for sizes:",
                colnames,
            )
            complete = False
        else:
            complete = True

    # So at this point we should have everything we need to continue

    print("Were we complete?", complete)

    if complete == False:
        # Create the region type
        data["RegType"] = xtype
        data["Major"] = xmajor
        data["Minor"] = 0.0
        data["Theta"] = 0.0
        # data["Color"] = "Unknown"

    # Make sure it is in the correct order
    xdata = data[
        "Source_name", "RA", "Dec", "RegType", "Major", "Minor", "Theta", "Color"
    ]

    print('XXX')
    print(xdata[:10])
    print('XXX')

    ttype = "unknown"
    return xdata, ttype


def write_regionfile(
    regionfile, records, source="unknown", type="unknown", color="red"
):
    """
    write_regionfile(regionfile,records):

    History::

        101231 ksl Increased the precision of the radii or major/minor
            axes to be hundredths of an arcsec for HST images where
            regions can be very small
        150220 ksl Modified to use astrpy tables.
    """

    # print(records[len(records)/2])

    # print('ok: got here')

    # First we need to understand how we are going to treat colors, namely whether
    # we need to write this out individually or collectively.  Our assumption is
    # that we will not change colors if the colors are different for different
    # lines of the master table

    records.info()

    xcolor = False
    i = 0
    while i < len(records):
        if records["Color"][i] != records["Color"][i - 1]:
            xcolor = True
            break
        i += 1

    if xcolor == True:
        print(
            "There are different colors for different region files, so we will preserve colors"
        )

    if color == "Unknown":
        color = "red"

    f = open(regionfile, "w")

    # write the header for the region file

    f.write("# Region file format: DS9 version 4.1\n")
    f.write("# Filename: %s\n" % regionfile)
    f.write(
        'global color=%s width=3 font="helvetica 14 bold" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
        % color
    )
    if type != "unknown":
        f.write("%s\n" % type)
    else:
        f.write("fk5\n")

    i = 0
    while i < len(records):
        z = records[i]
        name = z["Source_name"]
        ra = z["RA"]
        dec = z["Dec"]
        regtype = z["RegType"]
        x1 = z["Major"]
        zcolor = z["Color"]
        if type == "physical" or type == "image":
            if z["RegType"] == "ellipse" and z["Minor"] > 0.0:
                x2 = z["Minor"]
                theta = z["Theta"]
                if xcolor:
                    f.write(
                        "ellipse(%f,%f,%.2f,%.2f,%.1f) # color=%s text={%s}\n"
                        % (ra, dec, x1, x2, theta, zcolor, name)
                    )
                else:
                    f.write(
                        "ellipse(%f,%f,%.2f,%.2f,%.1f) # text={%s}\n"
                        % (ra, dec, x1, x2, theta, name)
                    )
            else:
                if xcolor:
                    f.write(
                        "circle(%f,%f,%.2f)  # color=%s text={%s}\n"
                        % (ra, dec, x1, zcolor, name)
                    )
                else:
                    f.write("circle(%f,%f,%.2f)  # text={%s}\n" % (ra, dec, x1, name))
        else:  # We assume this is a normal masterfile with RA, DECs and sizes in arcsec
            if (z["RegType"] == "ellipse" or z["RegType"] == "box") and z["Minor"] > 0.0:
                zztype=z["RegType"]
                x2 = z["Minor"]
                theta = z["Theta"]
                try:
                    if xcolor:
                        f.write(
                            '%s(%f,%f,%.2f",%.2f",%.1f) # color=%s text={%s}\n'
                            % (zztype,ra, dec, x1, x2, theta, zcolor, name)
                        )
                    else:
                        f.write(
                            '%s(%f,%f,%.2f",%.2f",%.1f) # text={%s}\n'
                            % (zztype,ra, dec, x1, x2, theta, name)
                        )
                except:
                    if xcolor:
                        f.write(
                            '%s(%s,%s,%.2f",%.2f",%.1f) # color=%s,text={%s}\n'
                            % (zztype,ra, dec, x1, x2, theta, zcolor, name)
                        )
                    else:
                        f.write(
                            '%s(%s,%s,%.2f",%.2f",%.1f) # text={%s}\n'
                            % (zztype,ra, dec, x1, x2, theta, name)
                        )
            elif z["RegType"] == "annulus":
                zztype=z["RegType"]
                x2 = z["Minor"]
                if xcolor:
                    f.write( '%s(%f,%f,%.2f",%.2f") # color=%s text={%s}\n'
                            % (zztype,ra, dec, x1, x2, zcolor, name))
                else:
                    f.write( '%s(%f,%f,%.2f",%.2f") #  text={%s}\n'
                            % (zztype,ra, dec, x1, x2, name))
            else:
                try:
                    if xcolor:
                        f.write(
                            'circle(%f,%f,%.2f")  # color=%s text={%s}\n'
                            % (ra, dec, x1, zcolor, name)
                        )
                    else:
                        f.write(
                            'circle(%f,%f,%.2f")  # text={%s}\n' % (ra, dec, x1, name)
                        )
                except TypeError:
                    if xcolor:
                        f.write(
                            'circle(%s,%s,%.2f")  # color=%s text={%s}\n'
                            % (ra, dec, x1, zcolor, name)
                        )
                    else:
                        f.write(
                            'circle(%s,%s,%.2f")  # text={%s}\n' % (ra, dec, x1, name)
                        )

        i = i + 1
    f.close()


def doit(argv):
    """
    Process the command line and call the main routines
    which are read_masterfile and write_masterfile  

    130830 Added extra parameter to allow changing the default radius
    """

    set_color = False
    color = "red"
    collated = "no"
    masterfile = ""
    regionfile = ""
    rad = 3.0

    i = 1
    while i < len(argv):
        if argv[i] == "-h":
            print(__doc__)
            return
        elif argv[i] == "-color":
            i = i + 1
            color = argv[i]
            set_color = True
        elif argv[i] == "-r":
            i = i + 1
            rad = eval(argv[i])
        elif argv[i] == "-collated":
            collated = "yes"
        else:
            if masterfile == "":
                masterfile = argv[i]
            elif regionfile == "":
                regionfile = argv[i]
        i = i + 1

    # Use the default region file name if the name is not yet assigned.
    if regionfile == "":
        regionfile = masterfile + ".reg"

    print("Masterfile          :", masterfile)
    print("Regionfile          :", regionfile)

    # Finished getting input information

    # read the masterfile

    data, ttype = read_masterfile(masterfile, xmajor=rad, xminor=rad)

    if data["Color"][0] == "Unknown":
        xcolor = color
    elif set_color:
        xcolor = color
    else:
        xcolor = data["Color"][0]

    # print('After return', data.colnames)
    # print('Got type',ttype)

    write_regionfile(regionfile, data, masterfile, ttype, xcolor)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":

    # Read the command line

    argc = len(sys.argv)

    if argc < 2:
        print(__doc__)
        sys.exit()

    doit(sys.argv)
