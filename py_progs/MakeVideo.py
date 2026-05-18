#!/usr/bin/env python3

'''
                    Space Telescope Science Institute

Synopsis:  

Create an Mpeg from a lot of pngs


Command line usage (if any):

    usage: MakeVideo.py [-h] [-out outroot] [-fps N] boat_load of pngs

    where:
        -h        print this help and exit
        -out root root name for output mp4 file (default: test.mp4)
        -fps N    frames per second (default: 4 = 0.25 s/frame); lower = slower, higher = faster

Description:  

    This is a simple routine that convets a bung of pngs to 
    an mpeg file for displaying.  It is very basic, so
    the pngs need to be presented in a logical order


Primary routines:

    doit

Notes:
                                       
History:

250409 ksl Coding begun

'''


import sys
import subprocess
import os
import tempfile

def create_mpeg_from_png_list(png_list, output_file='output.mp4', fps=4):
    """
    Create an MPEG video from a list of PNG files using ffmpeg.

    Parameters:
    - png_list: List of PNG file paths
    - output_file: Output video file name
    - fps: Frames per second
    """
    # Sort the list to ensure frame order
    png_list = sorted(png_list)

    # Duration per frame in seconds
    duration = 1.0 / fps

    # Concat file with explicit per-frame duration (avoids DTS timestamp errors)
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as list_file:
        for png in png_list:
            list_file.write(f"file '{os.path.abspath(png)}'\n")
            list_file.write(f"duration {duration}\n")
        list_path = list_file.name

    try:
        command = [
            'ffmpeg',
            '-y',
            '-f', 'concat',
            '-safe', '0',
            '-i', list_path,
            # Pad to even dimensions required by libx264
            '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2',
            '-c:v', 'libx264',
            '-pix_fmt', 'yuv420p',
            output_file
        ]
        subprocess.run(command, check=True)
    finally:
        os.remove(list_path)

def steer(argv):
    '''
    make_video.py -h -out outroot boat_load of pngs
    '''

    pngs=[]
    outname='test'
    fps = 4
    
    i=1
    while i<len(argv):
        if argv[i][:2]=='-h':
            print(__doc__)
            return
        elif argv[i][:4]=='-out':
            i+=1
            outname=argv[i]
        elif argv[i][:4]=='-fps':
            i+=1
            fps=int(argv[i])
        elif argv[i][0]=='-':
            print('Error: Unknown switch: ',argv[:10])
            return
        elif argv[i].count('png'):
            pngs.append(argv[i])

        i+=1

    if outname.count('.mp4')==0:
        outname='%s.mp4' % outname

    create_mpeg_from_png_list(pngs, output_file=outname, fps=fps)
    return

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
