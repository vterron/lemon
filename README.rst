|logo| LEMON
============

.. image:: https://travis-ci.org/vterron/lemon.png?branch=master
  :target: https://travis-ci.org/vterron/lemon

LEMON is a CCD differential-photometry pipeline, written in Python, developed at the `Institute of Astrophysics of Andalusia (CSIC) <http://www.iaa.es/>`_ and originally designed for its use at the `1.23m CAHA telescope <http://www.caha.es/telescopes-overview-and-instruments-manuals.html/>`_ for automated variable stars detection and analysis. The aim of this tool is to make it possible to *completely reduce thousands of images of time series* in a matter of hours and with minimal user interaction, if not none at all, automatically detecting variable stars and presenting the results to the astronomer.

A first overview of LEMON, now slightly outdated, was presented some time ago at `<http://adsabs.harvard.edu/abs/2011hsa6.conf..755T>`_.

Commands
========

The pipeline consists of **nine commands**, six of which are considered to be *essential* because they implement the data reduction and analysis steps and are usually run sequentially. However, depending on your needs only a specific subset of them may be used. In this sense, LEMON can be viewed as a set of tasks that *may* be used as a pipeline. The other three commands are *auxiliary* in the sense that they provide features that are convenient in some, but not all, scenarios.

.. code::

  usage: lemon [--help] [--version] [--update] COMMAND [ARGS]

  The essential commands are:
     astrometry   Calibrate the images astrometrically
     mosaic       Assemble the images into a mosaic
     photometry   Perform aperture photometry
     diffphot     Generate light curves
     periods      Dworetsky's string-length method
     juicer       LEMONdB browser and variability analyzer

  The auxiliary, not-always-necessary commands are:
     import       Group the images of an observing campaign
     seeing       Discard images with bad seeing or elongated
     annuli       Find optimal parameters for photometry

  See 'lemon COMMAND' for more information on a specific command.

Installation
============

- Current version: **0.2**
- View `CHANGELOG <./Misc/CHANGES>`_

LEMON stands on the shoulders of many giants, using excellent, robust programs developed by people much more skilled than us to detect sources, do aperture photometry and compute astrometric solutions of the FITS images. The disadvantage, however, is that for several of them there are not (yet?) Debian packages available, so they have to be installed manually. But fret not — none of them should be particularly difficult to get up and running on your system.

These are the steps to install LEMON on a fresh Debian 7 (`Wheezy <https://www.debian.org/releases/wheezy/>`_) machine:

1. ``apt-get install git python-pip csh``
#. ``apt-get build-dep python-matplotlib python-scipy``
#. ``apt-get install openmpi-dev`` # you may need this to compile Montage

#. ``git clone git://github.com/vterron/lemon.git ~/lemon``
#. ``cd ~/lemon``
#. ``pip install "numpy>=1.7.1"``
#. ``pip install -r pre-requirements.txt``
#. ``pip install -r requirements.txt``

#. Install `IRAF <http://iraf.noao.edu/>`_
#. Install `SExtractor <http://www.astromatic.net/software/sextractor>`_ (version 2.19.5 or newer) [#]_
#. Install `Astrometry.net <http://astrometry.net/use.html>`_
#. Install the MPI-enabled `Montage <http://montage.ipac.caltech.edu/docs/download2.html>`_ binaries [#]_
#. ``python ./setup.py``
#. ``echo 'PATH=$PATH:~/lemon' >> ~/.bashrc``
#. ``echo "source ~/lemon/lemon-completion.sh" >> ~/.bashrc``
#. ``./run_tests.py`` — optional, although recommended!

Note that, starting from version 2.16, IRAF is now released `under a free software license <ftp://iraf.noao.edu/iraf/v216/v216revs.txt>`_. There is, thus, reasonable hope that it may be packaged for drop-in installation in GNU/Linux systems in the near future, which would enormously simplify the process of installing LEMON. Until then, please bear with us.

.. |logo| image:: ./Misc/lemon-icon_200px.png
          :width: 200 px
          :alt: LEMON icon

.. [#] The important thing to `keep in mind <http://www.astromatic.net/forum/showthread.php?tid=587>`_ is that SExtractor does not rely on the `CLAPACK <http://www.netlib.org/clapack/>`_ implementation of `LAPACK <http://www.netlib.org/lapack/>`_ — instead, it only uses the subset of the LAPACK functions available in `ATLAS <http://math-atlas.sourceforge.net/>`_. That is the reason why, in case the ``liblapack-dev`` package is installed, you may encounter an error such as :code:`configure: error: CBLAS/LAPack library files not found at usual locations! Exiting`. If that is your case, you may need to do something like this:

.. code:: bash

  cd ./sextractor-2.19.5
  apt-get install fftw3-dev libatlas-base-dev
  update-alternatives --set liblapack.so /usr/lib/atlas-base/atlas/liblapack.so
  ./configure --with-atlas-incdir=/usr/include/atlas
  make
  make install

.. [#] Edit these two lines in ``Montage/Makefile.LINUX`` before doing ``make``

.. code:: bash

  # uncomment the next two lines to build MPI modules
  # MPICC  =	mpicc
  # BINS = 	$(SBINS) $(MBINS)
