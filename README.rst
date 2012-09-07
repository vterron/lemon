LEMON
=====

LEMON is a CCD differential-photometry pipeline, written in Python, developed at the `Institute of Astrophysics of Andalusia (CSIC) <http://www.iaa.es/>`_ and originally designed for its use at the `1.23m CAHA telescope <http://www.caha.es/telescopes-overview-and-instruments-manuals.html/>`_ for automated variable stars detection and analysis. The aim of this tool is to make it possible to *completely reduce thousands of images of time series* in a matter of hours and with minimal user interaction, if not none at all, automatically detecting variable stars and presenting the results to the astronomer.

A first overview of LEMON, now slightly outdated, was presented some time ago at `<http://adsabs.harvard.edu/abs/2011hsa6.conf..755T>`_.


Modules
=======

The pipeline consists of **ten command-line scripts**, which implement the data reduction and analysis steps and are usually run sequentially, although depending on your needs only a specific subset of them may be used. In this sense, and following the Unix
tools philosophy ("*write programs that do one thing and do it well*"), LEMON can be viewed as a set of tasks that *may* be used as a pipeline.

* import.py — group all the images of a campaign
* seeing.py — find image with best seeing
* offsets.py — determine translation offsets
* mosaic.py — create master frame
* astrometry.py — astrometric calibration
* annuli.py — find best parameters for photometry
* photometry.py — simply do photometry
* diffphot.py — generate light curves
* periods.py — string-length method
* mining.py — identify *interesting* stars

Installation
============

LEMON stands on the shoulders of many giants, using excellent, robust programs developed by people much more skilled than us to detect sources, do aperture photometry and compute astrometric solutions on the FITS images. The disadvantage, however, is that for many of them there are not (yet?) GNU/Debian packages available, so they have to be installed manually — the configuration of IRAF and PyRAF, although heavily simplified in recent versions, is particularly tedious and painful.

These are the steps to install LEMON on a clean GNU/Debian machine:

1. ``apt-get install python-numpy python-scipy python-pyfits python-lxml python-uncertainties plplot-bin sextractor``
#. Install `IRAF <http://iraf.noao.edu/>`_
#. Install `STSDAS/TABLES <http://www.stsci.edu/institute/software_hardware/stsdas/download-stsdas/>`_, followed by `PyRAF <http://www.stsci.edu/institute/software_hardware/pyraf/current/download/>`_ version 2.0 or newer.
#. Install the `CDSClient <http://cdsarc.u-strasbg.fr/doc/cdsclient.html>`_ package
#. Install `SCAMP <http://www.astromatic.net/software/scamp>`_ and `SWarp <http://www.astromatic.net/software/swarp>`_
#. ``git clone git://github.com/vterron/lemon.git ~/lemon``
#. ``python ~/lemon/setup.py``

Note that, starting from version 2.16, IRAF is now released `under a free software license <ftp://iraf.noao.edu/iraf/v216/v216revs.txt>`_. There is, thus, reasonable hope that it may be packaged for drop-in installation in GNU/Linux systems in the near future, which would enormously simplify the process of installing LEMON. Until then, please bear with us.

