.. _quickstart:

Quickstart
==========

This page gives an introduction in how to get started with LEMON,
showing the steps that would normally be necessary to reduce a data
set. In particular, this example assumes that we have a series of FITS
images from a transit of exoplanet `WASP-10b`_.

.. important::

  Your **data is assumed to be calibrated**. Bias and dark
  subtraction, flat-fielding correction and any other necessary steps
  should have been performed before any data is fed to the pipeline.

Let's take a look at each of the steps of the data reduction.

Do astrometry
-------------

Your FITS images need to be astrometrically solved before we can do
photometry. The ``astrometry`` command takes multiple FITS files as
input and writes a copy of them, containing the WCS header, to the
output directory specified in the last argument.

::

    $ lemon astrometry data/*.fits WASP10/
    >> The output directory 'WASP10' did not exist, so it had to be created.
    >> Using a local build of Astrometry.net.
    >> Doing astrometry on the 193 paths given as input.
    >> 100%[======================================================================>]
    >> You're done ^_^

As all LEMON commands, ``astrometry`` spawns multiple processes in
order to fully leverage all the processors your machine has. This
usually results in significant speed-ups on multi-core systems. If
that is not your case, this can always be modified via the ``--cores``
option.

The ``astrometry`` command is built on top of `Astrometry.net`_.

.. _WASP-10b: http://exoplanet.eu/catalog/wasp-10_b/
.. _Astrometry.net: http://astrometry.net/
