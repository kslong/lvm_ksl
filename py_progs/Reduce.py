#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Reduce one or more  LVM datasets from one or more MJDs


Command line usage (if any):

    usage: Reduce.py [-h] [-keep] [-cp] [-np N] [-get] [-force] exposures_to_process

    where -h prints this help message, -keep retains the ancillary files which are otherwise
    deleted (the default is to delete these files to save disk space), -cp causes the routine
    to copy the reduced frames to a directory ./data relative to where the program is being run,
    -np N is the number of threads to use to process the data, -get downloads the raw data
    from Utah for the requested exposures without running the pipeline on them (useful for
    pre-staging data ahead of time; -keep/-cp/-np are ignored in this mode), and -force
    attempts the reduction even when raw spectrograph frames or the science-telescope
    astrometry (agcam coadd) are missing -- by default such exposures are not reduced, and
    the reason is written to xlog/log_<exposure>.txt instead.

    The exposures_to_process is a string of words interpreted as follows: a word greater than
    50000 is treated as an MJD, a word less than 500 means process this exposure, 500-510 means
    process exposures 500 to 510 inclusive, and 500,505,513 means process exposures 500, 505,
    and 513.

    An example set of words might be: 60188 4155-4157 60189 4321 60190 5011,5012,5013 which would
    mean to process 4155-4157 from MJD 60188, 4321 from 60189, and 5011, 5012, and 5013 from 60190.




Description:  

    This routine downloades raw data to one's local verion of the lvm reduction
    data structure and process the data with the pipleline.
    One can optionally optionally copy the
    reduced data to a local ./data directory.  The normal log fils
    produced by the pipeline can be found in a local ./xlog directory.


Primary routines:


Notes:
                                       
History::

    240404 ksl Coding begun
    240526 ksl Adapt to new version of the DRP
    240527 ksl Add multiprocessing and allow for a more complicated
    input spectrum.
    260714 ksl get_data() rewritten to use sdss_access (Access.add for the
    9 raw camspec frames + agcam coadd) instead of hand-rolled rsync
    subprocess calls against ~/.sdss_rsync_password -- dtn.sdss.org now
    requires 2FA for that rsync auth path, which sdss_access sidesteps
    via .netrc. Also fixed a latent bug where the agcam coadd was
    written under the original mjd instead of the resolved qmjd when
    the exposure rolled over to mjd+1. Dropped unused ascii/fits/np/plt
    imports. The download half now works in any sdss_access env (e.g.
    ksl); process_one() still shells out to `drp run`, so the actual
    reduction step still requires lvmdrp26.
    260716 ksl Added ensure_metadata_store(), called once per distinct MJD
    at the top of do_many() before any parallel `drp run` processes for
    that MJD start. Fixes a FileExistsError race in lvmdrp's
    raw_metadata.hdf5 store: h5py opens that file with mode='a', which
    is a check-then-act "open, and if not found, create exclusively" --
    if two `drp run` processes for exposures on the same new MJD hit
    this at once, both can see "not found" and race to create it, and
    the loser crashes. Downloading one exposure and running `drp
    metadata regenerate -m` serially up front creates the store before
    the parallel processes start, so they only ever open an existing
    file.
    260717 ksl get_data() was only ever fetching the Sci telescope's agcam
    coadd (tel='sci'), never SkyE/SkyW. That meant load_guider_header() in
    lvmdrp always found "file not found" for the sky telescopes on locally
    reduced exposures, forcing the CMD-position fallback (ASRC='CMD
    position') every time regardless of what guider solution Utah actually
    had. Now fetches tel in ('sci', 'skye', 'skyw').
    260717 ksl Fixed do_many() only joining its last-created
    multiprocessing.Process instead of all of them: the loop that builds
    `jobs` leaves `p` referring to the final job, and the old `p.join()`
    at the end of do_many() waited on just that one. If an earlier job
    outlived the last one (e.g. it was still downloading/reducing while a
    later job started and quickly failed), do_many() -- and so
    Reduce.py's own process -- could return while that earlier job was
    still running in the background, letting a `source`d sequence of
    Reduce.py calls race ahead into the next MJD while the previous one
    was still active. This was observed directly: MJD 60194's exposure
    4415 was still writing its log at 12:11 while the next DoW28 line's
    (MJD 60203) logs had already started at 12:18. Now every job in
    `jobs` is explicitly joined before do_many() returns.
    260717 ksl process_one() now times the `drp run` subprocess call and
    reports the elapsed wall-clock seconds alongside the
    completed/FAILED message, so per-exposure reduction time is visible
    in the (unredirected) Reduce.py output, not just inferred from xlog
    file timestamps.
    260719 ksl Added -get, a download-only mode threaded through
    do_one()/do_many() as get_only: it calls get_data() as usual but
    returns before process_one() (the `drp run` call), so exposures can
    be pre-staged from Utah without reducing them. In do_many(), -get
    also skips ensure_metadata_store()'s serial priming step, since that
    exists only to protect concurrent `drp run` calls and isn't needed
    when nothing is being reduced.
    260719 ksl Fixed get_data() silently downloading nothing when an
    exposure has no agcam coadd (e.g. a calibration exposure taken
    without active guiding): it used to queue all 9 raw camspecs and all
    3 coadds in one sdss_access batch, and a.set_stream() raises before
    a.commit() ever runs if any single queued file is missing remotely,
    so one absent coadd aborted the whole batch and blocked every raw
    camspec too, even though they existed remotely. Raw camspecs and
    coadds are now fetched as two independent batches. Also added
    check_local_files(), which verifies -- via sdss_access's own local
    path templates -- which of the 9 raw camspecs and 3 coadds actually
    landed on disk, and prints a clear status summary before
    process_one() runs (or, under -get, before returning). A missing sci
    coadd is called out explicitly, since it makes lvmdrp's
    load_guider_header() fall back to the commanded (CMD) position
    instead of a real astrometric solution for the science telescope.
    get_data() now returns (qmjd, raw_ok, missing_coadd) instead of
    either qmjd or the string 'Failed', so callers get an honest
    raw_ok/missing_coadd status instead of an unconditional success.
    260719 ksl Added -force. By default, do_one() now refuses to run
    process_one() (the `drp run` call) if get_data() reports any raw
    camspec missing or the science-telescope agcam coadd missing --
    since a missing sci coadd means load_guider_header() would fall back
    to the commanded (CMD) position instead of a real astrometric
    solution. Rather than silently skipping, do_one() writes the reason
    to xlog/log_<exp>.txt via the new write_skip_log(), using "ERROR:"
    lines so CheckReduced.py's existing scan for ERROR lines across
    xlog/log*.txt flags it exactly like a failed reduction. -force
    overrides this block and attempts the reduction anyway. Threaded
    through do_many() and steer() alongside get_only.
    260719 ksl Trimmed parse()/steer()'s start-of-run diagnostic prints
    (the raw words, the parsed MJD/ExpNo lists, and the echoed xtab) down
    to a single "Beginning retrieval for ..." line, since they cluttered
    the top of every run without adding information beyond what the
    user just typed.

'''

from astropy.table import Table
import subprocess
import timeit
import time
import multiprocessing
multiprocessing.set_start_method("spawn",force=True)
import os
import sky_plot
import re
from collections import Counter
from sdss_access import Access


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


def clean_and_count_lines_with_keywords(filename):
    # Regular expression to remove non-ASCII characters
    non_ascii_regex = re.compile(r'[^\x00-\x7F]+')

    # Dictionary to store cleaned lines and their counts
    line_counts = Counter()

    # Read the file line by line
    with open(filename, 'r', encoding='utf-8') as file:
        for line in file:
            # Check if the line contains 'WARNING' or 'ERROR'
            if 'WARNING' in line or 'ERROR' in line:
                # Remove non-ASCII characters
                cleaned_line = non_ascii_regex.sub('', line)

                # Count the cleaned line
                line_counts[cleaned_line.strip()] += 1

    # Display results sorted by frequency (most common first)
    print('### WARNING/ERROR Summary ###')
    for line, count in line_counts.most_common():
        print(f"{count}: {line}")


def format_number(number):
    return "{:08d}".format(number)

def process_one(mjd,i,clean):

        # Regenerate metadata
        # metadata_process = subprocess.run(["drp", "metadata", "regenerate", "-m", mjd])
        # if metadata_process.returncode == 0:
        #     print("Metadata regenerated successfully.")
        # else:
        #     print("Failed to regenerate metadata.")

        print("Starting reduction with log in xlog/log_{}.txt".format(i))
        with open("xlog/log_{}.txt".format(i), "w") as logfile:
            # xcommand=["drp", "run", "-m",str(mjd),"-e", str(i)]
            xcommand=["drp", "run", "-e", str(i)]
            if clean==True:
                # xcommand=["drp", "run", "-c","-m",str(mjd),"-e", str(i)]
                xcommand=["drp", "run", "-c","-e", str(i)]

            print("Begin processing of ",i,"with command: ",xcommand)
            start_time = timeit.default_timer()
            reduction_process = subprocess.run(xcommand, stdout=logfile, stderr=subprocess.STDOUT)
            elapsed = timeit.default_timer() - start_time
        if reduction_process.returncode == 0:
            print(f"Reduction for {i} completed successfully in {elapsed:.1f} seconds.")
        else:
            print(f"FAILED to complete reduction for {i} on {mjd} after {elapsed:.1f} seconds, check log for errors.")

        print('Type for logfile', type(logfile))
        print('Path ',os.path.isfile(logfile.name))
        print('name', logfile.name)

        clean_and_count_lines_with_keywords(logfile.name)

        print("Finished reduction of", i)

        return  reduction_process.returncode


# Raw frames come from 3 cameras (b, r, z) x 3 spectrographs (1, 2, 3)
CAMSPECS = [f'{cam}{spec}' for cam in ('b', 'r', 'z') for spec in (1, 2, 3)]

# The telescopes with an agcam coadd; 'sci' is the one whose astrometry
# actually matters for fiber-to-RA/Dec, see check_local_files().
AGCAM_TELS = ('sci', 'skye', 'skyw')


def check_local_files(mjd,i):
    '''
    Verify, using sdss_access's own local path templates, that the raw
    camspec frames and agcam coadds for exposure i on mjd actually landed
    on disk, and print a clear summary of what is/isn't there.

    Returns raw_ok, missing_raw, missing_coadd:: the 9 raw camspecs are
    considered essential (raw_ok is False if any are missing), while a
    missing coadd is reported but does not by itself make raw_ok False.
    A missing sci coadd is called out explicitly, since
    load_guider_header() in lvmdrp falls back to the commanded (CMD)
    position instead of a real guider astrometric solution for the
    science telescope when that file is absent.
    '''

    a = Access(release='sdsswork')

    missing_raw = [camspec for camspec in CAMSPECS
                   if not os.path.isfile(a.full('lvm_raw', mjd=mjd, hemi='s', camspec=camspec, expnum=i))]

    missing_coadd = [tel for tel in AGCAM_TELS
                      if not os.path.isfile(a.full('lvm_agcam_coadd', mjd=mjd, tel=tel, specframe=i))]

    print('--- Local file check for exposure %s (MJD %s) ---' % (i,mjd))
    if not missing_raw:
        print('  Raw camspecs: all %d present' % len(CAMSPECS))
    else:
        print('  Raw camspecs: MISSING %d of %d: %s' % (len(missing_raw),len(CAMSPECS),', '.join(missing_raw)))

    if not missing_coadd:
        print('  Agcam coadds: all present (%s)' % ', '.join(AGCAM_TELS))
    else:
        print('  Agcam coadds: MISSING: %s' % ', '.join(missing_coadd))
        if 'sci' in missing_coadd:
            print('  WARNING: science-telescope agcam coadd is missing -- fiber astrometry')
            print('  for this exposure will fall back to the commanded (CMD) position')
            print('  instead of a real guider solution when it is reduced.')

    return len(missing_raw)==0, missing_raw, missing_coadd


def get_data(mjd,i):
    '''
    qmjd is a string

    Uses sdss_access (HTTPS + .netrc) rather than rsync, since dtn.sdss.org
    now requires 2FA for rsync access.

    Raw camspecs and agcam coadds are fetched as two independent
    sdss_access batches, then checked against what actually landed on
    disk (see check_local_files()). They used to be queued together in a
    single a.set_stream()/a.commit() pair -- but a.set_stream() raises
    before a.commit() ever runs if any one queued file is missing
    remotely, so one absent coadd (e.g. a calibration exposure with no
    agcam data) silently aborted the whole batch and prevented every raw
    camspec from being fetched, even when they existed remotely and would
    otherwise have downloaded fine.

    Always returns (qmjd, raw_ok, missing_coadd) -- it never hard-fails
    itself. Even when the exposure can't be found remotely under either
    mjd or mjd+1, it still falls through to check_local_files() (which
    will simply report everything missing) so callers have one single,
    consistent source of truth for whether there is enough data to
    process this exposure. Deciding whether to actually proceed anyway
    (-force) is left to the caller, not to get_data().
    '''

    os.environ["LVMAGCAM_DIR"] = os.path.join(os.environ["SAS_BASE_DIR"], "sdsswork/data/agcam/lco/")
    mjd='%s' % mjd
    xmjd='%d' % (int(mjd)+1)
    xnumb = format_number(i)

    a = Access(release='sdsswork')

    # A single representative camspec is enough to tell whether this exposure
    # was recorded under mjd or rolled over into mjd+1.
    if a.exists('lvm_raw', remote=True, mjd=mjd, hemi='s', camspec='b1', expnum=i):
        print('All is OK with %s so proceeding' % mjd)
        qmjd=mjd
        found_remotely=True
    elif a.exists('lvm_raw', remote=True, mjd=xmjd, hemi='s', camspec='b1', expnum=i):
        print('Failed on orginal %s, but succeeded with  %s' % (mjd,xmjd))
        qmjd=xmjd
        found_remotely=True
    else:
        print('Failed with both %s and %s -- no raw data found remotely for %s' % (mjd,xmjd,xnumb))
        qmjd=mjd
        found_remotely=False

    if found_remotely:
        try:
            a_raw = Access(release='sdsswork')
            a_raw.remote()
            for camspec in CAMSPECS:
                a_raw.add('lvm_raw', mjd=qmjd, hemi='s', camspec=camspec, expnum=i)
            a_raw.set_stream()
            a_raw.commit()
            print(f"Raw frames for {xnumb} successfully downloaded.")
        except Exception as e:
            print(f"Failed to download raw frames for {xnumb}: {e}")

        try:
            a_coadd = Access(release='sdsswork')
            a_coadd.remote()
            for tel in AGCAM_TELS:
                a_coadd.add('lvm_agcam_coadd', mjd=qmjd, tel=tel, specframe=i)
            a_coadd.set_stream()
            a_coadd.commit()
            print(f"Agcam coadds for {xnumb} successfully downloaded.")
        except Exception as e:
            print(f"Failed to download agcam coadds for {xnumb}: {e}")

    raw_ok, missing_raw, missing_coadd = check_local_files(qmjd,i)

    return qmjd, raw_ok, missing_coadd


def ensure_metadata_store(mjd,exp):
    '''
    Force-create the raw_metadata.hdf5 store for mjd before any parallel
    `drp run` processes start.

    h5py opens that store with mode='a', which is implemented as "try to
    open read/write, and if that raises FileNotFoundError, create the file
    exclusively" -- a check-then-act race. If several `drp run` processes
    for exposures on the same new MJD hit this at once, all of them can see
    "not found" and race to create it exclusively, and the losers crash
    with FileExistsError. Downloading one exposure and running `drp
    metadata regenerate` here, serially, creates the store up front so the
    parallel `drp run` calls only ever open an existing file.
    '''
    print('Pre-creating metadata store for MJD %s using exposure %s' % (mjd,exp))
    qmjd,raw_ok,missing_coadd=get_data(mjd,exp)
    if not raw_ok:
        print('WARNING: could not download exposure %s to prime metadata store for MJD %s' % (exp,mjd))
        return
    regen_process=subprocess.run(["drp","metadata","regenerate","-m",str(qmjd)])
    if regen_process.returncode != 0:
        print('WARNING: drp metadata regenerate failed for MJD %s' % qmjd)


def write_skip_log(exp,reason_lines):
    '''
    Write xlog/log_<exp>.txt recording why this exposure was not
    processed, in the same spot process_one() would otherwise have used,
    so CheckReduced.py's scan of xlog/log*.txt for lines containing
    "ERROR" flags it just like any other failed reduction.
    '''
    if os.path.isdir('./xlog')==False:
        os.mkdir('./xlog')
    logname="xlog/log_{}.txt".format(exp)
    with open(logname,"w") as logfile:
        for line in reason_lines:
            logfile.write(line+'\n')
    print('Reason for not processing %s written to %s' % (exp,logname))


def do_one(mjd,exp,clean=True,xcopy=True,get_only=False,force=False):
    '''
    mjd is a string, as is qmjd

    By default, if any raw camspec frame or the science-telescope agcam
    coadd is missing, the exposure is not reduced -- pass force=True to
    attempt the reduction anyway.
    '''

    qmjd,raw_ok,missing_coadd=get_data(mjd,exp)

    if mjd!=qmjd:
        print('Although these data were taken on %s, the are in SMJD %s' % (mjd,qmjd))
        mjd=qmjd

    if get_only:
        print("Downloaded %s only, as requested (-get); skipping reduction" % str(exp))
        return 0

    sci_missing = 'sci' in missing_coadd
    blocked = (not raw_ok) or sci_missing

    if blocked and not force:
        reason=['ERROR: exposure %s was NOT processed' % str(exp)]
        if not raw_ok:
            reason.append('ERROR: one or more raw spectrograph camspec frames are missing for this exposure')
        if sci_missing:
            reason.append('ERROR: the science-telescope agcam coadd (astrometry) is missing for this exposure')
        reason.append('Rerun Reduce.py with -force to attempt reduction anyway.')
        for line in reason:
            print(line)
        write_skip_log(exp,reason)
        return 1

    if blocked and force:
        print('WARNING: -force set -- attempting reduction of %s despite the missing data noted above' % str(exp))

    process_one(mjd,exp,clean)



    if xcopy:
        print("Running LocateData and copying files on %s" % str(exp))
        locate_data_process = subprocess.run(["LocateData.py", "-cp", str(exp), str(exp)])

        if locate_data_process.returncode == 0:
            print("LocateData.py executed successfully.")
        else:
            print("Failed to execute LocateData.py.")

    return 0
    
def get_no_jobs(jobs):
    '''
    Check how many jobs are running
    '''
    njobs=0
    for one in jobs:
        if one.is_alive():
            njobs+=1
    return njobs

def do_many(xtab,clean=True,xcopy=True,nproc=8,get_only=False,force=False):
    '''
    A routine to run the lvmdrp in parallel

    If get_only is True, this only downloads the raw data from Utah for
    each exposure and skips both the metadata-store priming (which exists
    solely to protect concurrent `drp run` calls) and the reduction itself.

    If force is False (the default), do_one() skips reduction (and writes
    a reason to xlog/log_<exp>.txt) for any exposure missing raw camspec
    frames or the science-telescope agcam coadd; force=True attempts the
    reduction regardless.
    '''



    if os.path.isdir('./xlog')==False:
        os.mkdir('./xlog')


    start_time = timeit.default_timer()

    if not get_only:
        # Prime the metadata store for each distinct MJD serially, before any
        # parallel `drp run` processes for that MJD start (see
        # ensure_metadata_store for why).
        seen_mjds=set()
        for one in xtab:
            mjd=one['MJD']
            if mjd not in seen_mjds:
                ensure_metadata_store(mjd,one['ExpNo'])
                seen_mjds.add(mjd)

    jobs=[]
    for one in xtab:
        if not get_only and int(one['MJD']) < 60177:
            print('WARNING: THESE DATA ARE UNLIKELY BE CALIBRATABLE WITHOUT SPECIAL EFFORT, AS THEY WERE OBTAINED BEFORE MJD 60177')

        p=multiprocessing.Process(target=do_one,args=[one['MJD'],one['ExpNo'],clean,xcopy,get_only,force])
        jobs.append(p)

    i=0
    while i<nproc and i<len(jobs):
        t = time.localtime()
        print('STARTING  %d of %d:\n' % (i+1,len(xtab)))
        one=jobs[i]
        one.start()
        time.sleep(60) # Space out starts to reduce concurrency issue with hdf metadate file
        i+=1

    njobs=get_no_jobs(jobs)

    while i<len(jobs):
        time.sleep(2)
        njobs=get_no_jobs(jobs)

        while njobs<nproc and i<len(jobs):
            t = time.localtime()
            print('STARTING  %d of %d:\n' % (i+1,len(xtab)))
            one=jobs[i]
            one.start()
            time.sleep(60)  # Space out starts to reduceconcurrency issues with hdf metadata file 
            njobs+=1
            i+=1

    for one in jobs:
        one.join()
        one.close()

    elapsed = timeit.default_timer() - start_time
    print('Completed multiprocessing of  %d exposures  ' % (len(xtab)))

    return



def doit(mjd,first_exp,last_exp,clean=True,xcopy=False):
    '''
    Process a consecutive sequence of exposures from the same MJD
    with our without sky subtraction.
    '''

    if int(mjd) < 60177:
        print('WARNING: THESE DATA ARE UNLIKELY BE CALIBRATABLE WITHOUT SPECIAL EFFORT, AS THEY WERE OBTAINED BEFORE MJD 60177')

    mjd='%s' % mjd
    scani=first_exp
    scanf=last_exp

    if os.path.isdir('./xlog')==False:
        os.mkdir('./xlog')



    # Get the necessary FITS files
    for i in range(scani, scanf + 1):
        do_one(mjd,i,clean)
    
    if xcopy:
        print("Running LocateData and copying files")
        locate_data_process = subprocess.run(["LocateData.py", "-cp", str(scani), str(scanf)])
    else:
        print("Running LocateData w/o copying files")
        locate_data_process = subprocess.run(["LocateData.py", str(scani), str(scanf)])

    if locate_data_process.returncode == 0:
        print("LocateData.py executed successfully.")
    else:
        print("Failed to execute LocateData.py.")


def parse(words):

    mjd=[]
    xexp=[]
    mjd_now=0
    for word in words:
        
        if word.count('-'):
            xword=word.split('-')
            if len(xword)==2:
                imin=int(xword[0])
                imax=int(xword[1])
                i=imin
                while i<=imax:
                    mjd.append(mjd_now)
                    xexp.append(i)
                    i+=1
            else:
                print('Badly formated inputs',words)
                return []
        elif word.count(','):
            xword=word.split(',')
            for one in xword:
                mjd.append(mjd_now)
                xexp.append(int(one))
        else:
            try:
                value=int(word)
                if value>50000:
                    mjd_now=value
                else:
                    mjd.append(mjd_now)
                    xexp.append(value)
            
 
            except:
                print('Badly formated inputs',words)
                return []
    xtab=Table([mjd,xexp],names=['MJD','ExpNo'])
    return xtab
    


def steer(argv):
    '''
    Control the flow of the program
    '''

    i=1
    mjd=-1
    first_exp=-1
    last_exp=-1
    clean=True
    xcopy=False
    clean=True
    nproc=1
    get_only=False
    force=False

    words=[]

    while i<len(argv):
        if argv[i][0:2]=='-h':
            print(_usage_from_doc(__doc__))
            return
        if argv[i]=='-np':
            i+=1
            nproc=int(argv[i])
        elif argv[i]=='-cp':
            xcopy=True
        elif argv[i]=='-keep':
            clean=False
        elif argv[i]=='-get':
            get_only=True
        elif argv[i]=='-force':
            force=True
        elif argv[i][0]=='-':
            print('Error: Could not parse command line: ', argv)
            return
        else:
            words.append(argv[i])
        i+=1

    print('Beginning retrieval for %s' % ' '.join(words))

    xtab=parse(words)

    if len(xtab)==0:
        print('Error: Could not parse command line: ', argv)
        return
        
    if nproc==1:
        if os.path.isdir('./xlog')==False:
            os.mkdir('./xlog')

        for one in xtab:
            do_one(one['MJD'],one['ExpNo'],clean,get_only=get_only,force=force)
    else:
        do_many(xtab,clean,xcopy,nproc,get_only,force)

if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        steer(sys.argv)
    else:
        print (__doc__)
