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
   threshold (see the :option:`--counts` option), the file is marked as
   saturated and discarded. This allows for a reasonable number of saturated
   pixels while at the same time excluding those FITS files that are
   essentially useless for any scientific purpose, such as skyflats that went
   out of hand.

#. Although not used by default, is it possible to define a Unix-style pattern
   that the object name of the FITS files must match in order to be
   imported. This is useful if you are interested in working only with specific
   sets of FITS files contained in a directory tree. For example, the pattern
   ``NGC2264*`` only imports those whose object name starts with *NGC2264*. For
   more information, see the :option:`--pattern` option.

#. The remaining FITS files are sorted by their date of observation, which is
   defined as the date at the start of observation plus half of the total
   exposure time — in other words, the exact date at the middle of the
   observation.

#. The FITS files are copied to the output directory and renamed
   sequentially. The first FITS file in chronological order is `assigned the
   number zero`_ and the following numbers go from there. The most common
   filename among the imported files is used to determine the basename to which
   these sequence numbers are appended, but this behavior can be changed with
   the :option:`--filename` option.

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


Options
=======

.. cmdoption:: --object <patterns>

  List of case-insensitive patterns, according to `the rules used by the Unix
  shell`_, and separated by commas, of the object names of the FITS files to
  import. Those FITS files whose object name (see the :option:`--objectk`
  option) matches one or more of these patterns are imported, while the rest
  are ignored. There cannot be spaces around the commas that separate the
  patterns, as if that is case what follows a whitespace is considered a
  different argument to the program. By default, all object names are matched
  (pattern ``*``).

  Examples: ``--object Andromeda`` imports only those FITS files whose object
  name is exactly that. A pattern such as ``Trumpler 37`` must be either
  quoted, ``'Trumpler 37'``, or have its whitespace escaped, ``Trumpler\
  37``. Finally, ``--object 'skyflat*,lampflat*'`` imports those FITS files
  whose object name starts with ``skyflat`` or ``lamplat``.

.. cmdoption:: --pattern <pattern>

  The case-insensitive, `Unix-style pattern`_ that the filename of a FITS file
  must match to be detected when the input directory trees are walked down.
  Files with a non-matching filename are ignored. The pattern must be quoted or
  escaped to prevent wildcard expansion, if any. We use the term *filename*
  because that is what it means for an end user, but the technical term would
  be *basename*: the name of the file along with its extension, such as
  ``GJ436-006V.fits``. By default, all filenames are imported.

  Examples: ``--object '*.fit'`` imports only those FITS files with the
  ``.fit`` extension, while ``'GJ436*.fits'`` imports those whose name starts
  with ``GJ436`` and have the ``.fits`` extension.

  .. _the rules used by the Unix shell:
  .. _Unix-style pattern:
     https://en.wikipedia.org/wiki/Glob_(programming)#Syntax

.. cmdoption:: --counts <ADUs>

   Number of :abbr:`ADUs (Analog-to-digital unit)` at which saturation occurs.
   The median of the pixel distribution is computed for each FITS file, and
   those with a value above this threshold are discarded. If this option is not
   used, no file is discarded no matter what its median number of ADUs is.

   Example: ``--counts 50000`` imports only those FITS files whose median
   number of ADUs is equal to or less than 50,000.

.. cmdoption:: --filename <prefix>

   The base name common to the copies, made in the output directory, of all the
   imported files. The sequence number, once the FITS files are sorted by their
   date of observation, is appended to this prefix before the extension.
   Leading zeros are used so that the filenames of the copies of the imported
   FITS files are all of equal length. If we import 100 FITS files, for
   example, sequence numbers can be written with no more than two digits, so
   the first file will be assigned the sequence number ``00`` and the last
   ``99``.

   Example: ``--filename WASP-44b_``, assuming that we are importing a total of
   437 files with the ``.fit`` extension, makes the first file copied to the
   output directory have the name ``WASP-44b_000.fit``, while the last one is
   named ``WASP-44b_436.fit``.

.. cmdoption:: --follow

   By default, when detecting FITS files we do not walk down into symbolic
   links that resolve to directories. Use this option to visit directories
   pointed to by symlinks, on systems that support them. This can lead to
   infinite recursion if a link points to a parent directory of itself.

.. cmdoption:: --exact

   For each imported FITS file, the `HISTORY`_ keyword is used to store both
   the path to the original file and the date at which it was imported. In
   addition, the copy of each imported file has its own path stored in the
   keyword specified with the :option:`--uik` option.

   Use this option in case you do not want to modify the FITS files, but
   instead prefer to work with an exact copy. The FITS headers will be left
   untouched and, so that even the most paranoid among us can rest assured that
   the copy of each file is identical, the `SHA-1 hash`_ is used to verify
   their integrity.

   .. _SHA-1 hash: https://en.wikipedia.org/wiki/SHA-1


.. _import-keywords:

FITS keywords
-------------

In order to correctly process the FITS files, access is needed to some of the
information stored in their headers. The default keywords where the necessary
values are looked for are those used by `PANIC`_, so you do not need to tinker
with these options if your data was taken at that instrument. If that is not
the case, you **must** make sure that the keywords here defined exactly match
those present in your FITS files.  Failing to do so will result in apocalyptic
consequences — the least severe of them being LEMON aborting its execution.

.. _PANIC: https://w3.iaa.csic.es/PANIC/

.. cmdoption:: --datek <keyword>

   The date of the observation, in the Y2K compliant date format specified in
   `the FITS standard`_: ``yyyy-mm-dd`` or ``yyyy-mm-ddTHH:MM:SS[.sss]``
   (default: ``DATE-OBS``)

   .. _the FITS standard: http://fits.gsfc.nasa.gov/fits_standard.html

.. cmdoption:: --expk <keyword>

   The exposure time in seconds (default: ``EXPTIME``)

.. cmdoption:: --objectk <keyword>

   The name of the object observed (default: ``OBJECT``)

.. cmdoption:: --uik <keyword>

   Along with some book-keeping information using the `HISTORY`_ keyword, the
   copies of the imported FITS files also have their own path stored in their
   headers, using the keyword defined by this option. This provides, since
   keywords propagate when FITS file are manipulated, a means of getting the
   path to the original file in case they are calibrated or modified in any
   other way. In case you do not want the path to be saved to the header, set
   this option to an empty string (``''``) to disable it (default: ``UNCIMG``)

   .. note::

      The path to the original FITS file is needed when we do aperture
      photometry, in order to check which pixels are saturated: calibration
      steps such as bias subtraction or, particularly, flat-fielding may cause
      a pixel to go below the saturation level when, in actuality, before the
      calibration took place it was above, or vice versa. This keyword, thus,
      allows LEMON to do this check in the original file, which is the one that
      matters.

   .. _HISTORY: http://archive.stsci.edu/fits/fits_standard/node40.html#SECTION00942420000000000000


Examples
========

Let's see three examples very similar to what you may need to use::

  lemon import /root/data/2013-03-28/ ~/Trumpler37/ --counts 45000 --filename "Trumpler_37_" --exact

This first example detects all the FITS files in ``/root/data/2013-03-28/`` or
any of its subdirectories, and saves a copy of them to ``~/Trumpler37/``, which
is created if it does not exist. Because of ``--counts 45000``, files whose
median number of ADUs is greater than this number are discarded. The name of
the files copied to the output directory starts with the value of
:option:`--filename`, to which the sequence number is appended: the first one,
therefore, could be named, for example, ``Trumpler_37_0000.fits`` — note that
how many zeros are used depends on the number of files imported. Finally,
:option:`--exact` guarantees that the FITS files copied to the output directory
are an identical copy of the original, with no book-keeping information added
to their headers.

::

  lemon import ~/ ~/exoplanets --object "WASP*b,HD*b,Gliese*b" --follow

In this second example, LEMON scans your entire home directory, copying to
``~/exoplanets/`` those FITS files whose object name matches any of these
patterns: ``WASP*b``, ``HD*b`` and ``Gliese*b``. Examples of object names that
would be matched are ``WASP-44b``, ``WASP-1 B`` (patterns are
case-insensitive), ``HD 100655 b`` and ``Gliese 876 d``, while others such as
``HAT-P-30-WASP-51 b``, ``HD 10180 g`` or ``Gliese 876 e`` would not be so. Due
to the presence of :option:`--follow`, LEMON will walk down into symbolic links
that point to directories.

::

  lemon import /disk-b/obs12_images/ /data/ --pattern "*.fit[s]"

Here LEMON walks down the directory ``/disk-b/obs12_images``, detecting all the
FITS files contained there or in any of its subdirectories and making a copy of
them to ``/data/``. Thanks to ``--pattern "*.fit[s]"``, the search for these
FITS files is restricted to those with the extensions ``.fit`` and ``.fits``.
This illustrates how :option:`--pattern` may be used to considerably speed up
the execution time of this command, as by default it checks whether all the
regular files it comes across are standard-conforming FITS files.


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

