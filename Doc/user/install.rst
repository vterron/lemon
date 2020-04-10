.. _install:

Installation
============

Release v\ |version|.

The installation of LEMON is a somewhat tedious —although not particularly difficult— process, involving several dependencies for which there not exist Debian packages.

These are the steps to install LEMON on a fresh `Debian 7`_ machine:

1. ``apt-get install git python-pip csh``
#. ``apt-get build-dep python-matplotlib python-scipy``
#. ``apt-get install openmpi-dev``
#. ``easy_install -U distribute``
#. ``git clone --branch v0.3 git://github.com/vterron/lemon.git ~/lemon``
#. ``cd ~/lemon``
#. ``git submodule add git@github.com:vterron/barycorr.git``
#. ``pip install "numpy>=1.7.1"``
#. ``pip install -r pre-requirements.txt`` # :download:`[View] <../../pre-requirements.txt>`
#. ``pip install -r requirements.txt`` # :download:`[View] <../../requirements.txt>`

#. Install IRAF_.
#. Install SExtractor_ (version 2.19.5 or newer) [#]_
#. Install `Astrometry.net`_.
#. Install the MPI-enabled Montage_ binaries [#]_
#. ``python ./setup.py``
#. ``echo 'PATH=$PATH:~/lemon' >> ~/.bashrc``
#. ``echo "source ~/lemon/lemon-completion.sh" >> ~/.bashrc``
#. ``./run_tests.py`` — optional, although recommended!

Note that in October 2018 NOAO began the transition of IRAF to an end-of-support state. IRAF is now maintained by the `iraf-community project <https://iraf-community.github.io/>`_ on GitHub, which integrates any available patches into the source code. You can download the IRAF 2.16.1 snapshot `here <https://github.com/iraf-community/iraf/releases/tag/v2.16.1+2018.11.01>`_.

.. [#] The important thing to `keep in mind <http://www.astromatic.net/forum/showthread.php?tid=587>`_ is that SExtractor does not rely on the CLAPACK_ implementation of LAPACK_ — instead, it only uses the subset of the LAPACK functions available in ATLAS_. That is the reason why, in case the ``liblapack-dev`` package is installed, you may encounter an error such as :code:`configure: error: CBLAS/LAPack library files not found at usual locations! Exiting`. If that is your case, you may need to do something like this:

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


.. note::

   LEMON is **not** yet available on PyPI_, but we intend to package it soon. This will enormously simplify the installation process, which should consist of a single ``pip install lemon`` command — provided that IRAF_, SExtractor_, `Astrometry.net`_ and Montage_ are already installed on your system.

.. _Debian 7: https://www.debian.org/releases/wheezy/
.. _IRAF: http://iraf.noao.edu/
.. _SExtractor: http://www.astromatic.net/software/sextractor
.. _Astrometry.net: http://astrometry.net/use.html
.. _Montage: http://montage.ipac.caltech.edu/docs/download2.html
.. _CLAPACK: http://www.netlib.org/clapack/
.. _LAPACK: http://www.netlib.org/lapack/
.. _ATLAS: http://math-atlas.sourceforge.net/
.. _PyPI: https://pypi.python.org/pypi
