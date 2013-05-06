.. _commands-import:

######
import
######

The purpose of :command:`import` is to automatically detect all the FITS files
belonging to an observation campaign, get rid of those that are saturated or
irrelevant for the data reduction (e.g., test images that were taken when
centering the field) and copy them sequentially to the working directory. To
put it more simply, this command walks down a directory tree and makes a copy
of all the FITS files that are of our interest.

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

