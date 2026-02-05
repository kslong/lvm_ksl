Installation
============

This page describes how to install the lvm_ksl package from GitHub.


Requirements
------------

lvm_ksl requires Python 3.8 or later and the following packages:

- numpy
- scipy
- astropy
- matplotlib
- psutil

Optional dependencies for full functionality:

- lvmdrp (LVM Data Reduction Pipeline)
- sdss_access (for retrieving data from Utah)
- skycorr (ESO sky correction tool)


Installation from GitHub
------------------------

Clone the repository from GitHub using HTTPS::

    git clone https://github.com/kslong/lvm_ksl.git

Or, if you have a GitHub account with SSH keys configured::

    git clone git@github.com:kslong/lvm_ksl.git

The package is script-based and does not require a formal installation.
Add the ``py_progs`` directory to your Python path:

**Bash/Zsh** (add to your ``.bashrc`` or ``.zshrc``)::

    export PYTHONPATH="/path/to/lvm_ksl/py_progs:$PYTHONPATH"

**Csh/Tcsh** (add to your ``.cshrc``)::

    setenv PYTHONPATH /path/to/lvm_ksl/py_progs:$PYTHONPATH

Replace ``/path/to/lvm_ksl`` with the actual path where you cloned the
repository.


Running Scripts
---------------

Scripts can be run directly from the command line::

    python /path/to/lvm_ksl/py_progs/rss_combine.py -h

Or, if ``py_progs`` is in your PATH::

    rss_combine.py -h

Most scripts accept ``-h`` to display usage information.


Importing as a Module
---------------------

With the PYTHONPATH set, you can import modules in Python::

    from lvm_ksl import rss_combine
    from lvm_ksl import lvm_gaussfit


Updating
--------

To update to the latest version::

    cd /path/to/lvm_ksl
    git pull


Environment Variables
---------------------

Some scripts use environment variables to locate data:

LVM_MASTER_DIR
    Directory containing LVM master calibration files.

LVMAGCAM_DIR
    Directory containing LVM AG camera data.

These are optional and only needed for specific workflows.


Building the Documentation
--------------------------

The documentation is built using Sphinx with the AutoAPI extension,
which automatically generates API documentation from the Python source
code.

**Required packages**

Install the documentation dependencies::

    pip install sphinx sphinx-rtd-theme sphinx-autoapi

**Building**

To build the HTML documentation::

    cd /path/to/lvm_ksl/docs
    make html

The output will be in ``docs/html/``. Open ``docs/html/index.html`` in
a browser to view the documentation.

**Rebuilding**

To rebuild from scratch (useful if you've added new modules)::

    cd /path/to/lvm_ksl/docs
    make clean
    make html
