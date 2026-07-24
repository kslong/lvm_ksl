#!/usr/bin/env python
# coding: utf-8

'''
                    Space Telescope Science Institute

Synopsis:

Wire up and run lvm-dap-conf (LVM-DAP) on a FITS spectrum from the current
directory, without having to hand-edit output_path or recreate the _legacy
symlink for every new run directory.


Command line usage (if any):

    usage: RunDap.py [-h] [-force] [-template file] [-config name] [-out label] fits_file

    where fits_file is the input FITS spectrum, with or without a .fits
    suffix, and may be a relative path (e.g. xdata/foo/whatever.fits).
    -out label sets the label that becomes the root of the filenames
    lvm-dap-conf creates (the same label DoDap/lvm-dap-conf normally
    takes); if omitted, label defaults to fits_file's basename with any
    directory and .fits/.fits.gz suffix stripped, e.g.
    xdata/foo/whatever.fits becomes whatever. -force regenerates the
    yaml config from the template even if one already exists in the
    current directory (by default an existing config is left untouched,
    so hand-tuned parameters survive a second run). -template lets you
    point at a different yaml template than the bundled one; -config
    lets you name the generated config file something other than the
    default.

Description:

    lvm-dap-conf takes a yaml config file whose output_path must be an
    absolute path -- a relative output_path silently corrupts (writes end
    up under a bogus path rooted at "/" instead of the run directory).
    auto_ssp_LVM.config (referenced by the yaml's config-file key) also
    contains bare relative paths like "_legacy/Ha_LVM.config" that are
    opened relative to whatever directory the fit is run from, which is
    why every run directory needs its own _legacy symlink pointing at
    $LVM_DAP/_legacy.

    This routine does both of those automatically in the current
    directory -- creates _legacy if missing, creates output_dap/, writes
    ksl-dap_v110.yaml from the template in ../data/dap_ref with
    output_path/lvmdap_dir filled in -- and then runs lvm-dap-conf.

Primary routines:

    default_label
    wire_run_dir
    run_dap
    steer

Notes:

    See ../data/dap_ref/ksl-dap_v110.template.yaml for the template that
    gets copied into each run directory. Edit that once to change fitting
    parameters (RSP file, w-range, etc.) for all future runs.

History::

    260723 ksl Coding begun
    260723 ksl label is now optional, via -out, instead of a required
    positional argument; when omitted it defaults to fits_file's
    basename with any directory and .fits/.fits.gz suffix stripped.

'''

import os
import re
import subprocess
import sys

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE_DEFAULT = os.path.join(REPO_DIR, 'data', 'dap_ref', 'ksl-dap_v110.template.yaml')
CONFIG_DEFAULT = 'ksl-dap_v110.yaml'


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


def resolve_fits(fits_file):
    '''
    Accept fits_file with or without the .fits suffix (matching the old
    DoDap convention of `lvm-dap-conf $1.fits ...`), and return an
    absolute path. Returns None (after printing why) if neither form
    exists.
    '''
    if os.path.isfile(fits_file):
        return os.path.abspath(fits_file)
    with_suffix = fits_file + '.fits'
    if os.path.isfile(with_suffix):
        return os.path.abspath(with_suffix)
    print('Error: could not find %s or %s' % (fits_file, with_suffix))
    return None


def default_label(fits_file):
    '''
    Derive a label from fits_file when -out is not given: strip any
    directory and a trailing .fits or .fits.gz suffix, e.g.
    xdata/foo/whatever.fits -> whatever.
    '''
    base = os.path.basename(fits_file)
    if base.endswith('.fits.gz'):
        return base[:-len('.fits.gz')]
    if base.endswith('.fits'):
        return base[:-len('.fits')]
    return base


def ensure_legacy_symlink():
    '''
    Create ./_legacy -> $LVM_DAP/_legacy if it does not already exist.
    Required because auto_ssp_LVM.config's own contents reference
    "_legacy/..." paths that lvm-dap-conf opens relative to the current
    working directory, not relative to lvmdap_dir. Refuses to touch
    _legacy if it already exists and points somewhere else, rather than
    silently replacing it.
    '''
    target = os.path.join(os.environ['LVM_DAP'].rstrip('/'), '_legacy')
    if os.path.islink('_legacy') or os.path.exists('_legacy'):
        if os.path.islink('_legacy') and os.path.realpath('_legacy') == os.path.realpath(target):
            return True
        print('Error: ./_legacy already exists and does not point at %s -- fix or remove it by hand' % target)
        return False
    os.symlink(target, '_legacy')
    print('Created _legacy -> %s' % target)
    return True


def write_config(template, config_name, force):
    '''
    Generate config_name from template, filling in output_path (must be
    absolute -- see module docstring) and lvmdap_dir with the current
    directory and $LVM_DAP respectively. Leaves an existing config_name
    untouched unless force is True, so hand-tuned parameters from a
    previous run survive.
    '''
    if os.path.isfile(config_name) and not force:
        print('Using existing %s (pass -force to regenerate from the template)' % config_name)
        return True

    if not os.path.isfile(template):
        print('Error: template not found: %s' % template)
        return False

    output_path = os.path.join(os.getcwd(), 'output_dap') + '/'
    lvmdap_dir = os.environ['LVM_DAP'].rstrip('/') + '/'

    with open(template) as f:
        text = f.read()
    text = text.replace('__OUTPUT_PATH__', output_path)
    text = text.replace('__LVMDAP_DIR__', lvmdap_dir)

    with open(config_name, 'w') as f:
        f.write(text)
    print('Wrote %s (output_path=%s)' % (config_name, output_path))
    return True


def wire_run_dir(template, config_name, force):
    '''
    Set up everything lvm-dap-conf needs in the current directory:
    output_dap/, _legacy, and the yaml config. Returns True if the
    directory is ready to run.
    '''
    if not os.path.isdir('output_dap'):
        os.makedirs('output_dap')
        print('Created ./output_dap')

    if not ensure_legacy_symlink():
        return False

    return write_config(template, config_name, force)


def run_dap(fits_path, label, config_name):
    '''
    Invoke lvm-dap-conf on fits_path/label/config_name, streaming its
    output directly, and return its exit code.
    '''
    cmd = ['lvm-dap-conf', fits_path, label, config_name]
    print('Running: %s' % ' '.join(cmd))
    result = subprocess.run(cmd)
    return result.returncode


def steer(argv):

    force = False
    template = TEMPLATE_DEFAULT
    config_name = CONFIG_DEFAULT
    label = None
    words = []

    i = 1
    while i < len(argv):
        if argv[i][0:2] == '-h':
            print(_usage_from_doc(__doc__))
            return
        elif argv[i] == '-force':
            force = True
        elif argv[i] == '-out':
            i += 1
            label = argv[i]
        elif argv[i] == '-template':
            i += 1
            template = argv[i]
        elif argv[i] == '-config':
            i += 1
            config_name = argv[i]
        elif argv[i][0] == '-':
            print('Unknown option :', argv)
            return
        else:
            words.append(argv[i])
        i += 1

    if len(words) != 1:
        print('Error: expected a single fits_file, got:', words)
        print(_usage_from_doc(__doc__))
        return

    fits_file = words[0]
    if label is None:
        label = default_label(fits_file)
        print('No -out given, using label: %s' % label)

    fits_path = resolve_fits(fits_file)
    if fits_path is None:
        return

    if not wire_run_dir(template, config_name, force):
        return

    return run_dap(fits_path, label, config_name)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        retcode = steer(sys.argv)
        sys.exit(retcode if isinstance(retcode, int) else 0)
    else:
        print(__doc__)
