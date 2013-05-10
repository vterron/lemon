.. _commands-import:

######
import
######

The purpose of :command:`import` is to automatically detect all the FITS files
belonging to an observation campaign and copy them to our working directory.
This is particularly needed when reducing data taken at observatories that
enforce :ref:`specific naming conventions <import-historical-note>`, where the
astronomer may easily end up with hundreds of images scattered over dozens of
directories, many times with duplicate filenames.


Overview
========

This command walks down a series of directory trees, detects the FITS files,
gets rid of those that are saturated or irrelevant for the data reduction and
saves a copy of the remaining to an output directory. But this is just a very
succinct explanation. Let's take a closer look at what :command:`import` does:

#. The input directories are recursively walked down and all the FITS files
   contained in them are automatically detected. LEMON considers a FITS file
   anything that conforms to `the FITS standard`_, regardless of its
   extension. Note that, among many other things, this means that the `SIMPLE
   keyword`_, containing a logical constant with the value ``T``, is required
   to be the first one in the primary header of your files — do not worry, this
   is almost certainly the case.

   .. _SIMPLE keyword: http://archive.stsci.edu/fits/fits_standard/node39.html#SECTION00941110000000000000

#. Then, saturated files are discarded. The saturation is defined in terms of
   the entire FITS file: if the median number of counts (:abbr:`ADUs
   (Analog-to-digital unit)`) of all the pixels is above a certain saturation
   threshold (see the ``--counts`` option), the file is marked as saturated and
   discarded. This allows for a reasonable number of saturated pixels while at
   the same time excluding those FITS files that are essentially useless for
   any scientific purpose, such as skyflats that went out of hand.

#. Although not used by default, is it possible to define a Unix-style pattern
   that the object name of the FITS files must match in order to be
   imported. This is useful if you are interested in working only with specific
   sets of FITS files contained in a directory tree. For example, the pattern
   ``NGC2264*`` only imports those whose object name starts with *NGC2264*. For
   more information, see the ``--pattern`` option.

#. The remaining FITS files are sorted by their date of observation, which is
   defined as the date at the start of observation plus half of the total
   exposure time — in other words, the exact date at the middle of the
   observation.

#. The FITS files are copied to the output directory and renamed
   sequentially. The first FITS file in chronological order is `assigned the
   number zero`_ and the following numbers go from there. The most common
   filename among the imported files is used to determine the basename to which
   these sequence numbers are appended, but this behavior can be changed with
   the ``--filename`` option. Please refer to its documentation for further
   details.

   .. _assigned the number zero: `would have disagreed`_

.. note::

   Due to the algorithms used for image alignment and source detection, LEMON
   requires that all the FITS files of an observation campaign have the same
   size (number of pixels along the x- and y-axes). In case there are multiple
   sizes among the input FITS files, only those with the most common size will
   be imported.


Usage
=====

The :command:`import` command accepts a variable number of arguments: FITS
files are detected by recursively walking through all the paths except the
last one, which is used as the output directory: ::

  $ lemon import [OPTION]... INPUT_DIRS... OUTPUT_DIR

The extension of the FITS files is irrelevant, as LEMON does not pay attention
to it. The only thing that matters is whether they conform to `the FITS
standard`_. If the output directory does not exist, it is created for you.

Input paths may be directories or individual files. Therefore, it is possible
to import whole directories, specific sets of FITS files or a combination of
both. Non-existent paths and non-standard FITS files are silently ignored. For
example: ::

  $ lemon import ~/2013-01-22/ ~/2012-01-23/M101_*fits ~/M101_raw

This will import all the FITS files in the directory tree that starts at
``~/2013-01-22/`` and, in addition, those in ``~/2012-01-23`` whose name begins
with *M101_* and have the *.fits* extension, provided that they conform to the
FITS standard. The files are copied to ``~/M101_raw``, which is created if it
does not exist.

.. _the FITS standard: http://fits.gsfc.nasa.gov/fits_standard.html


.. _import-historical-note:

A historical note
=================

The origins of this command trace back to our observation campaigns with the
optical CCD at the `1.23m CAHA telescope`_, whose `manual`_ instructs observers
to follow a strict naming convention: a separate directory, *yymmdd*, must be
used for each night's data, and FITS filenames must have the form *nnnF_*, with
*nnn* being the first three letters of the surname of the :abbr:`PI (Principal
Investigator)` and *F* the first letter, in upper case, of the first name.

.. _1.23m CAHA telescope: http://www.caha.es/telescopes-overview-and-instruments-manuals.html/
.. _manual: http://www.caha.es/CAHA/Instruments/IA123/ObsManual.pdf

In our case, where the campaigns lasted a full month, this meant that the files
were spread over thirty different directories, each one of them containing
images that always had the same filename. These names were not only as cryptic
as this misguided convention doomed them to be [#]_ but, even worse, they
always started from one (something with which Dijkstra, by the way, `would have
disagreed`_). Therefore, when the time to reduce the data came, we had thirty
different FITS files named *ferM_0001.fits* (as per our PI, `Matilde
Fernández`_), another thirty named *ferM_0002.fits*, and so forth::

     $ ls
     111115  111119  111123  111127  111201  111205  111209  111213
     111116  111120  111124  111128  111202  111206  111210  111214
     111117  111121  111125  111129  111203  111207  111211
     111118  111122  111126  111130  111204  111208  111212
     $ ls 111115/
     ferM_0001.fits  ferM_0010.fits  ferM_0019.fits  ferM_0028.fits
     ferM_0002.fits  ferM_0011.fits  ferM_0020.fits  ferM_0029.fits
     ferM_0003.fits  ferM_0012.fits  ferM_0021.fits  ferM_0030.fits
     ferM_0004.fits  ferM_0013.fits  ferM_0022.fits  ferM_0031.fits
     ferM_0005.fits  ferM_0014.fits  ferM_0023.fits  ferM_0032.fits
     ferM_0006.fits  ferM_0015.fits  ferM_0024.fits  ferM_0033.fits
     ferM_0007.fits  ferM_0016.fits  ferM_0025.fits  ferM_0034.fits
     ferM_0008.fits  ferM_0017.fits  ferM_0026.fits
     ferM_0009.fits  ferM_0018.fits  ferM_0027.fits

We could not simply move the FITS files to our working directory, as the names
would collide, so in 2009 we wrote a :download:`Bash script <./rename.sh>` to
rename them sequentially. It was eventually rewritten in Python, gradually
incorporated additional functionality and was finally merged into the code of
LEMON.

.. _would have disagreed: http://www.cs.utexas.edu/~EWD/transcriptions/EWD08xx/EWD831.html
.. _Matilde Fernández: http://www.iaa.es/~matilde/

.. [#] This convention may seem reasonable at first glance, but it comes at the
   expense of stripping the filename of rather useful information that would
   greatly simplify the job of astronomers.  Imagine, for example, how
   different it would be to come across a file named
   *ferM_0056_OrionF1_20minV.fits* instead of one that just says
   *ferM_0056.fits*. As the Python mantra goes, explicit is better than
   implicit.

