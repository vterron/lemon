#! /usr/bin/env python

# Copyright (c) 2012 Victor Terron. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of LEMON.
#
# LEMON is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import division

description = """
This module does aperture photometry on a series of images. In its most basic
usage, it only receives two parameters: (a) an image of the field, usually
taken under excellent atmospheric conditions, from which the reference
instrumental magnitude of each object is extracted and (b) a LEMON XML file
which lists the horizontal and vertical translation offsets of all the images
with respect to the master frame. In this manner, after each object has been
detected on the master frame, photometry can be done on the same exact position
in all the images, as the translation objects are known.

The output of the module is a LEMON database, already indexed in order to
considerably speed up the queries that are most likely to be done.

Most, if not all, of the aspects of the module are configurable. For example,
which aperture and sky annuli (instead of the default values) are used or on
which exact pixels photometry is done (instead of on all the sources that are
automatically detected on the master frame) can be specified. Please refer to
the manual for further information.

"""

import collections
import hashlib
import logging
import multiprocessing
import numpy
import optparse
import os
import os.path
import pwd
import socket
import sys
import time
import warnings

# LEMON modules
import astromatic
import customparser
import database
import defaults
import fitsimage
import keywords
import methods
import passband
import qphot
import seeing
import style
import xmlparse

# The Queue is global -- this works, but note that we could have
# passed its reference to the function managed by pool.map_async.
# See http://stackoverflow.com/a/3217427/184363
queue = multiprocessing.Queue()

class InputFITSFiles(collections.defaultdict):
    """ Map each photometric filter to a list of FITS files.

    A convenience class to simplify how we work with the images on which
    photometry has to be done: it is a defaultdict, mapping each photometric
    filter to a list of the corresponding FITS files, but that iterates over
    the values (the FITS files, independently of their filter). This way, it
    may be viewed as a sequence of FITS files that also allows access by the
    photometric filter.

    """

    def __init__(self):
        super(InputFITSFiles, self).__init__(list)

    def __iter__(self):
        """ Iterate over the FITS files, regardless or their filter """
        for pfilter in self.itervalues():
            for img in pfilter:
                yield img

    def __len__(self):
        """ Return the number of FITS files, in any filter """
        return sum(len(pfilter) for pfilter in self.itervalues())

    def remove(self, img):
        """ Remove all occurrences of the FITS file, in any filter.

        Loop over the different photometric filters and remove, from each
        associated list of FITS files, all occurrences of 'img'. Returns the
        number of elements that were deleted. Emits a warning for each file
        that is removed, stating that is has been 'excluded'.

        """

        discarded = 0
        for pfilter_imgs in self.itervalues():
            # Iterate backwards, modify original list in situ
            for index in reversed(xrange(len(pfilter_imgs))):
                if pfilter_imgs[index] == img:
                    del pfilter_imgs[index]
                    discarded += 1
                    msg = "%s %s excluded."
                    warnings.warn(msg % (style.prefix, img))
        return discarded


def parallel_photometry(args):
    """ Method argument of map_async to do photometry in parallel.

    Functions defined in classes don't pickle, so we have moved this code here
    in order to be able to use it with multiprocessing's map_async. As it
    receives a single argument, values are passed in a tuple which is then
    unpacked.

    """

    img_offset, reference_pixels, qphot_params, options = args
    aperture, annulus, dannulus = qphot_params

    photometry_image = \
        seeing.FITSeeingImage(img_offset.shifted, options.maximum,
                                 options.margin, coaddk = options.coaddk)
    img_path = photometry_image.path
    logging.debug("Doing photometry on image %s" % img_path)

    saturation_level = photometry_image.saturation
    qphot_result = qphot.run(img_offset, aperture, annulus, dannulus,
                             saturation_level, options.exptimek,
                             options.uncimgk, reference_pixels)
    logging.info("Done photometry on image %s " % img_path)
    logging.debug("%s: photometry returned %d records (expected = %d)" % \
                 (img_path, len(qphot_result), len(reference_pixels)))
    assert len(qphot_result) == len(reference_pixels)

    # The filter name and qphot parameters are stored in their own
    # table and referred to using their ID (primary key) in order not
    # to have duplicated information in the database. The aperture,
    # annulus and dannulus, on the other hand, were received as input
    # by the program.

    img_pfilter = img_offset.filter
    img_unix_time = img_offset.date
    img_object = img_offset.object
    img_airmass = img_offset.airmass
    img_gain = options.gain or photometry_image.read_keyword(options.gaink)
    img_xoffset = img_offset.x
    img_xoverlap = img_offset.x_overlap
    img_yoffset = img_offset.y
    img_yoverlap = img_offset.y_overlap

    logging.debug("%s: filter = %s" % (img_path, img_pfilter))
    logging.debug("%s: observation date: %.2f (%s)" % \
                  (img_path, img_unix_time, methods.utctime(img_unix_time)))
    logging.debug("%s: object = %s" % (img_path, img_object))
    logging.debug("%s: airmass = %.4f" % (img_path, img_airmass))
    gain_comes_from_msg = "given by user" if options.gain else "read from header"
    logging.debug("%s: gain (%s) = %.4f" % \
                  (img_path, gain_comes_from_msg, img_gain))
    logging.debug("%s: x-offset: %.4f" % (img_path, img_xoffset))
    logging.debug("%s: x-overlap: %.d" % (img_path, img_xoverlap))
    logging.debug("%s: y-offset: %.4f" % (img_path, img_yoffset))
    logging.debug("%s: y-overlap: %.d" % (img_path, img_yoverlap))

    logging.debug("%s: photometry aperture: %.3f" % (img_path, aperture))
    logging.debug("%s: photometry annulus: %.3f" % (img_path, annulus))
    logging.debug("%s: photometry dannulus: %.3f" % (img_path, dannulus))
    pparams = database.PhotometricParameters(aperture, annulus, dannulus)

    pimage = \
        database.Image(img_path, img_pfilter, pparams, img_unix_time,
                       img_object, img_airmass, img_gain, img_xoffset,
                       img_yoffset, img_xoverlap, img_yoverlap)

    queue.put((pimage, img_offset, qphot_result))


parser = customparser.get_parser(description)
parser.usage = "%prog [OPTION]... OFFSETS_XML_FILE"

parser.add_option('--overwrite', action = 'store_true', dest = 'overwrite',
                  help = "overwrite output database if it already exists")

parser.add_option('--filter', action = 'store', type = 'passband',
                  dest = 'filter', default = None,
                  help = "do not do photometry on all the FITS files given as "
                  "input, but only on those taken in this photometric filter. "
                  "The supported systems are Johnson, Cousins, Gunn, SDSS, "
                  "2MASS, Stromgren and H-alpha, but letters, designating a "
                  "particular section of the electromagnetic spectrum, may "
                  "also be used without a system (e.g., 'V')")

parser.add_option('--maximum', action = 'store', type = 'int',
                  dest = 'maximum', default = defaults.maximum,
                  help = defaults.desc['maximum'])

parser.add_option('--margin', action = 'store', type = 'int',
                  dest = 'margin', default = defaults.margin,
                  help = defaults.desc['margin'])

parser.add_option('--gain', action = 'store', type = 'float',
                  dest = 'gain', default = None,
                  help = "the gain of the CCD, in e-/ADU. Needed in order to "
                  "accurately calculate the SNR of each measurement. In case "
                  "this option is given, the value will not be read from the "
                  "FITS header (--gaink option)")

parser.add_option('--annuli', action = 'store', type = str,
                  dest = 'xml_annuli', default = None,
                  help = "ignore the 'Aperture Photometry' (FWHM and pixels) "
                  "sections below an instead read the aperture and sky annuli "
                  "to be used from a XML file outputted by the LEMON "
                  "annuli.py module. This file should, of course, have been "
                  "generated for the same set of images on which photometry "
                  "is now being done. For this same reason, the execution "
                  "will be aborted if the XML file does not have information "
                  "for all the filters listed in OFFSETS_XML_FILE. Note that "
                  "the parameters in the 'Aperture Photometry' sections are "
                  "still used when photometry is done on the reference image "
                  "to extract the instrumental magnitude of each image.")

parser.add_option('--cores', action = 'store', type = 'int',
                  dest = 'ncores', default = defaults.ncores,
                  help = defaults.desc['ncores'])

parser.add_option('-v', '--verbose', action = 'count',
                  dest = 'verbose', default = defaults.verbosity,
                  help = defaults.desc['verbosity'])

pixels_group = optparse.OptionGroup(parser, "List of Pixels",
               "By default, the module runs SExtractor on the reference image "
               "in order to determine which stars must have photometry done. "
               "Alternatively, photometry may be done on a fixed number of "
               "pixels by means of this option. It should be noted that this "
               "method is exclusive with the automatic detection of sources "
               "in the image; that is, if this option is used, photometry "
               "will be done only on these pixels.")

pixels_group.add_option('--pixels', action = 'store', type = str,
                        dest = 'list', default = None,
                        help = "path to the file which lists the pixels on "
                        "which photometry must be done. These coordinates "
                        "must be listed one per line: x_coordinate, "
                        "whitespace, y-coordinate. Improperly-formatted "
                        "lines are ignored.")
parser.add_option_group(pixels_group)

qphot_group = optparse.OptionGroup(parser, "Aperture Photometry (FWHM)",
              "Photometry is done on the images by IRAF's qphot, the quick "
              "aperture photometer, which computes accurate sky values and "
              "magnitudes for a series of objects. Instead of using "
              "absolute values (such as, for example, 8, 11 or 13 pixels), "
              "the value of the following parameters are defined in terms of "
              "the *median* FWHM of the images in each band.\n\n"

              "As the full width at half-maximum of an image is defined here "
              "as the median of FWHM (as reported by SExtractor) for all the "
              "detected stars, you may think of the median 'FWHM of the "
              "images' as the 'median of the median of the FWHM of all the "
              "stars, if you find it easier. Note also that a different "
              "'median FWHM' is computed for each passband on which "
              "photometry is done.")

qphot_group.add_option('--aperture', action = 'store', type = 'float',
                       dest = 'aperture', default = 3.0,
                       help = "the aperture radius, in number of times the "
                       "median FWHM [default: %default]")

qphot_group.add_option('--annulus', action = 'store', type = 'float',
                       dest = 'annulus', default = 4.5,
                       help = "the inner radius of the sky annulus, in "
                       "number of times the median FWHM [default: %default]")

qphot_group.add_option('--dannulus', action = 'store', type = 'float',
                       dest = 'dannulus', default = 1.0,
                       help = "the width of the sky annulus, in number "
                       "of times the median FWHM [default: %default]")

qphot_group.add_option('--min-sky', action = 'store', type = 'float',
                       dest = 'min', default = 3.0,
                       help = "the minimum width of the sky annulus, in "
                       "pixels, regardless of the value specified in the "
                       "above parameter. This option is intended to prevent "
                       "small FWHMs from resulting in too thin an sky "
                       "annulus. [default = %default]")

qphot_group.add_option('--individual-fwhm', action = 'store_true',
                       dest = 'individual_fwhm',
                       help = "consider FITS images individually when the "
                       "FWHM (and, therefore, the derived aperture and sky "
                       "annuli) is computed. That is, instead of using the "
                       "median FWHM of *all* the images in the same filter, "
                       "when this option is present the aperture and sky "
                       "annuli used to do photometry on an image depend "
                       "exclusively on the FWHM of *that* image.")
parser.add_option_group(qphot_group)

qphot_fixed = optparse.OptionGroup(parser, "Aperture Photometry (pixels)",
              "In case the exact size of the aperture and sky annuli is "
              "known, their dimensions can be specified in pixels. The three "
              "parameters are required, so if one or more are not specified "
              "the module will ignore them and revert back to the calculation "
              "of the sizes based on the FWHM of the images. In other and "
              "more concise, emphatic words: use them all or none at all!\n\n"

              "Note that setting the aperture photometry parameters to a "
              "fixed size, as these options do, means that the same values "
              "are used regardless of the filter with which the images were "
              "observed, and even for the reference image, so you should "
              "probably always use these options in conjunction with "
              "--filter, in order to to photometry one filter at a time.")

qphot_fixed.add_option('--aperture-pix', action = 'store', type = 'float',
                       dest = 'aperture_pix', default = None,
                       help = "the aperture radius, in pixels")

qphot_fixed.add_option('--annulus-pix', action = 'store', type = 'float',
                       dest = 'annulus_pix', default = None,
                       help = "the inner radius of the sky annulus, in pixels")

qphot_fixed.add_option('--dannulus-pix', action = 'store', type = 'float',
                       dest = 'dannulus_pix', default = None,
                       help = "the width of the sky annulus, in pixels")
parser.add_option_group(qphot_fixed)

key_group = optparse.OptionGroup(parser, "FITS Keywords",
                                 keywords.group_description)

key_group.add_option('--filterk', action = 'store', type = 'str',
                     dest = 'filterk', default = keywords.filterk,
                     help = keywords.desc['filterk'])

key_group.add_option('--datek', action = 'store', type = 'str',
                     dest = 'datek', default = keywords.datek,
                     help = keywords.desc['datek'])

key_group.add_option('--timek', action = 'store', type = 'str',
                     dest = 'timek', default = keywords.timek,
                     help = keywords.desc['timek'])

key_group.add_option('--expk', action = 'store', type = 'str',
                     dest = 'exptimek', default = keywords.exptimek,
                     help = keywords.desc['exptimek'])

key_group.add_option('--coaddk', action = 'store', type = 'str',
                     dest = 'coaddk', default = keywords.coaddk,
                     help = keywords.desc['coaddk'])

key_group.add_option('--gaink', action = 'store', type = 'str',
                     dest = 'gaink', default = keywords.gaink,
                     help = keywords.desc['gaink'])

key_group.add_option('--uik', action = 'store', type = 'str',
                     dest = 'uncimgk', default = keywords.uncimgk,
                     help = keywords.desc['uncimgk'])

parser.add_option_group(key_group)
customparser.clear_metavars(parser)

def main(arguments = None):
    """ main() function, encapsulated in a method to allow for easy invokation.

    This method follows Guido van Rossum's suggestions on how to write Python
    main() functions in order to make them more flexible. By encapsulating the
    main code of the script in a function and making it take an optional
    argument the script can be called not only from other modules, but also
    from the interactive Python prompt.

    Guido van van Rossum - Python main() functions:
    http://www.artima.com/weblogs/viewpost.jsp?thread=4829

    Keyword arguments:
    arguments - the list of command line arguments passed to the script.

    """

    if arguments is None:
        arguments = sys.argv[1:] # ignore argv[0], the script name

    # All the parameters must be cast to string before beging passed to the
    # option parser, as apparently it only expects str data types. Otherwise,
    # optparse may raise a TypeError exception and complains about how int or
    # float objects are unsuscriptable.
    arguments = [str(param) for param in arguments]
    (options, args) = parser.parse_args(args = arguments)

    # Adjust the logger level to WARNING, INFO or DEBUG, depending on the
    # given number of -v options (none, one or two or more, respectively)
    logging_level = logging.WARNING
    if options.verbose == 1:
        logging_level = logging.INFO
    elif options.verbose >= 2:
        logging_level = logging.DEBUG
    logging.basicConfig(format = style.LOG_FORMAT, level = logging_level)

    # Print the help and abort the execution if there are not two positional
    # arguments left after parsing the options, as the user must specify at
    # least one (only one?) input FITS file and the output LEMON database.
    if len(args) < 2:
        parser.print_help()
        return 2     # 2 is generally used for command line syntax errors
    else:
        sources_img_path = args[0]
        input_paths = set(args[1:-1])
        output_db_path = args[-1]

    # If the user gives an empty string as the FITS keyword which stores the
    # path to the original image, it is understood as meaning that we want
    # saturation to be checked for in the same images on which photometry is
    # done. If that is the case, we need to set the option to None (although an
    # empty string would also work), as that is what qphot.run expects to
    # receive in these cases.

    if not options.uncimgk:
        options.uncimgk = None

    # The annuli XMl file (that generated by the annuli command, and specified
    # with the --annuli option) must exist. The use of this file automatically
    # discards whathever was specified with the Aperture Photometry (FWHM and
    # pixels) options.
    xml_annuli = None
    fixed_annuli = False

    if options.xml_annuli:
        if not os.path.exists(options.xml_annuli):
            print "%sError. The file '%s' does not exist." % \
                  (style.prefix, options.xml_annuli)
            print style.error_exit_message
            return 1
        else:
            xml_annuli = \
                xmlparse.CandidateAnnuli.xml_load(options.xml_annuli)
            print "%sPhotometric paramaters read from the '%s' file." % \
                  (style.prefix, os.path.basename(options.xml_annuli))

    # Even if the annuli XML is used by the module, the aperture and sky annuli
    # parameters (whether FWHM-based or given in pixels) are still needed in
    # order to do photometry on the reference image and extract the
    # instrumental magnitude that for each star is stored in the database.

    # This warning is printed if one or more, but not the three, options are given
    if 0 < bool(options.aperture_pix) + bool(options.annulus_pix) + \
           bool(options.dannulus_pix) < 3:

        msg = ("%(p)sWarning: the --aperture-pix, --annulus-pix and "
               "--dannulus-pix options must\n"
               "%(p)salways be used in conjunction. The sizes of the apertures "
               "and sky annuli\n"
               "%(p)swill still be determined by the FWHMs.\n"
               "%(p)s") % {'p': style.prefix}
        warnings.warn(msg)

    if options.aperture_pix and options.annulus_pix and options.dannulus_pix:
        fixed_annuli = True

    # Abort the execution if the user gives, at the same time, the
    # --aperture-pix, --annulus-pix, --dannulus-pix and --annuli options. It
    # does not make sense to set the aperture and sky annuli to a fixed value
    # and simultaneously specify that these very values must be read from the
    # XML file. Does not compute.
    if fixed_annuli and xml_annuli:
        print ("%sError. The --aperture-pix, --annulus-pix and --dannulus-pix "
              "options are incompatible with --annuli." % style.prefix)
        print style.error_exit_message
        return 1

    if options.individual_fwhm:

        # If the photometric parameters are set to a fixed value, they cannot
        # be also derived from the FWHM of each image.
        if fixed_annuli:
            print "%sError. The --aperture-pix, --annulus-pix and " \
                  "--dannulus-pix options are incompatible with " \
                  "--individual-fwhm." % style.prefix
            print style.error_exit_message
            return 1

        # The same applies to --annuli: if the photometric parameters are read
        # from the XML file, they cannot also depend on the FWHM of the images.
        if xml_annuli:
            print "%sError. The --annuli option is incompatible with " \
                  "--individual-fwhm." % style.prefix
            print style.error_exit_message
            return 1

    # The aperture, annulus and dannulus values, whether expressed in number of
    # times the median FWHM or by a fixed number of pixels, must be positive
    # numbers. By definition, also, the inner radius of the sky annulus must be
    # greater than or equal to the aperture radius. Obviously!

    fwhm_options = (options.aperture, options.annulus, options.dannulus)
    pixel_options = (options.aperture_pix, options.annulus_pix, options.dannulus_pix)

    if (not fixed_annuli and min(fwhm_options)  <= 0) or \
           (fixed_annuli and min(pixel_options) <= 0):
        print "%sError. The aperture, annulus and dannulus values must be " \
              "positive numbers." % style.prefix
        print style.error_exit_message
        return 1

    if (not fixed_annuli and options.aperture > options.annulus) or \
           (fixed_annuli and options.aperture_pix > options.annulus_pix):
        print "%sError. The aperture radius (%.2f) must be smaller than or equal\n" \
              "%sto the inner radius of the sky annulus (%.2f)" % \
              (style.prefix, options.aperture_pix if fixed_annuli else options.aperture,
               style.prefix, options.annulus_pix if fixed_annuli else options.annulus)

        print style.error_exit_message
        return 1

    # If a list of pixels has been specified, read the file and extract them.
    # We also need to make sure that the file exists and complain otherwise.
    # The execution must be also aborted if the list of pixels is empty.
    if options.list:

        if not os.path.exists(options.list):
            print "%sError. The pixels list '%s' does not exist." % \
                  (style.prefix, options.list)
            print style.error_exit_message
            return 1

        else:
            # methods.load_file_list(path) returns a list of two-element tuples
            # which we cast to astromatic.Pixel instances for later convenience
            list_of_pixels = [astromatic.Pixel(*coords) for coords \
                              in methods.load_file_list(options.list, warn = True)]
            if not list_of_pixels:
                print "%sError. The pixels list '%s' is empty." % \
                      (style.prefix, options.list)
                print style.error_exit_message
                return 1

    # Each campaign must be saved to its own LEMON database, as it would not
    # make much sense to merge data (as the same tables would be used) of stars
    # that belong to different fields. Thus, by default the module will refuse
    # to work with an existing database (which is what the LEMONdB class would
    # do otherwise, as other modules shall need to update its contents) unless
    # the overwrite option is set, in which case it will be deleted and thus
    # created again from scratch.

    if os.path.exists(output_db_path):
        if not options.overwrite:
            print "%sError. The output database '%s' already exists." % \
                  (style.prefix, output_db_path)
            print style.error_exit_message
            return 1
        else:
            os.unlink(output_db_path)

    # Loop over all the input FITS files, mapping (a) each photometric filter
    # to a list of the FITS images that were observed in it, and (b) each FITS
    # image to its date of observation (UTC), in Unix time.

    msg = "%sExamining the headers of the %s FITS files given as input..."
    print msg % (style.prefix, len(input_paths))

    def get_date(img):
        """ Return the date() of a FITSImage object """
        return img.date(date_keyword = options.datek,
                        time_keyword = options.timek,
                        exp_keyword = options.exptimek)

    files = InputFITSFiles()
    img_dates = {}

    methods.show_progress(0.0)
    for index, img_path in enumerate(input_paths):
        img = fitsimage.FITSImage(img_path)
        pfilter = img.pfilter(options.filterk)
        files[pfilter].append(img_path)

        date = get_date(img)
        img_dates[img_path] = date

        percentage = (index + 1) / len(input_paths) * 100
        methods.show_progress(percentage)

    print # progress bar doesn't include newline

    msg = "%s%d different photometric were detected:"
    print msg % (style.prefix, len(files))

    for pfilter, images in files.iteritems():
        msg = "%s %s: %d files (%.2f %%)"
        percentage = len(images) / len(files) * 100
        print msg % (style.prefix, pfilter, len(images), percentage)

    # Light curves, which are our ultimate goal, can only have one magnitude
    # for each point in time. Therefore, we cannot do photometry on two or more
    # images with the same observation date. This may seem (and is) unlikely to
    # happen, but astronomical instruments also have software errors -- we have
    # already come across this while reducing images taken with Omega 2000, a
    # camera for the 3.5m CAHA.
    #
    # We might be tempted to keep one of them, such as, for example, that with
    # the highest number of sources, or the best FWHM. However, the safest bet
    # is to discard them all, because we cannot know which image has the right
    # observation date. The only certain thing is that an error occurred. Thus,
    # we better forget about these images.

    msg = "%sMaking sure there are no images with the same observation date..."
    print msg % style.prefix ,
    sys.stdout.flush()

    # Map each date of observation (UTC), in Unix time, to a list of the
    # corresponding FITS files. If there are no duplicate dates, all the values
    # of the dictionary will have a length of one. We do not use the Counter
    # class, from the collections module, for Python 2.6 compatibility.

    dates_counter = collections.defaultdict(list)
    for img_path, date in img_dates.iteritems():
        dates_counter[date].append(img_path)

    if max(len(imgs) for imgs in dates_counter.itervalues()) == 1:
        print 'done.'

    else:

        # "Making sure..." message printed above does not include newline.
        print

        # Find the dates for which there is more than one FITS file. Then,
        # remove from the InputFITSFiles object each of these files with a
        # duplicate date.
        discarded = 0
        for date, images in dates_counter.iteritems():
            nimgs = len(images)
            if nimgs > 1:
                msg = "%sWarning! Multiple images have date %s"
                warnings.warn(msg % (style.prefix, methods.utctime(date)))
                for img in images:
                    discarded += files.remove(img)

        # There should be no duplicate observation dates anymore, and at least
        # two FITS files should have been discarded. Otherwise, how did we get
        # to the 'else' clause in the first place?
        if __debug__:
            dates = []
            for image in files:
                dates.append(img_dates[image])
            assert len(set(dates)) == len(dates)
            assert discarded >= 2

        msg = "%s%d images had duplicate dates and were discarded, %d remain."
        print msg % (style.prefix, discarded, len(files))

    # The --filter option allows the user to specify on which FITS files, among
    # all those received as input, photometry must be done: only those files in
    # the options.filter photometric filter. The Passband class, which supports
    # comparison operations, makes it possible to compare filters for what they
    # really are, not how they were written: "Johnson V" and "johnson_v", for
    # example, are the same filter after all, but if we just compared the two
    # strings we would consider them to be different.
    if options.filter:

        msg = "%sIgnoring images not taken in the '%s' photometric filter..."
        print msg % (style.prefix, options.filter) ,
        sys.stdout.flush()

        discarded = 0
        for pfilter, images in files.items():
            if pfilter != options.filter:
                discarded += len(images)
                del files[pfilter]

        if not files:
            print
            msg = "%sError. No image was taken in the '%s' filter."
            print msg % (style.prefix, options.filter)
            print style.error_exit_message
            return 1

        else:
            print 'done.'
            msg = "%s%d images taken in the '%s' filter, %d were discarded."
            print msg % (style.prefix, len(files), options.filter, discarded)

    # If an XML file is specified with --annuli, it must list the photometric
    # parameters for all the filters on which photometry is to be done.
    if xml_annuli:
        for pfilter in files.iterkeys():
            if pfilter not in xml_annuli.iterkeys():
                msg = ("%sError. Photometric parameters for the '%s' "
                       "filter not listed in '%s. Wrong file, maybe?")
                xml_basename = os.path.basename(options.xml_annuli)
                print msg % (style.prefix, pfilter, xml_basename)
                print style.error_exit_message
                return 1

    # But, hey, shouldn't we make sure that the reference image about which we
    # are talking that much does really _exist_ and is a FITS, standard-
    # conforming file? Just to make sure, you know. In order to do that, simply
    # pass it to FITSImage, which will throw the appropiate exceptions if an
    # error is encountered.
    reference_img = fitsimage.FITSImage(xml_offsets.reference['path'])

    # Just a reminder of which the reference image (i.e., the image on which
    # sources are detected and to which the offsets in the XML file refer) and
    # the output database are.
    print "%sReference image: %s" % (style.prefix, reference_img.path)
    print "%sOutput database: %s" % (style.prefix, output_db_path)

    # Let the user know on how many coordinates, that were listed in the pixels
    # file, photometry will be done. If no list of pixels was given to the
    # module, run SExtractor on the reference image, discard those sources too
    # close to the edges of the image and create a list of pixels with the
    # coordinates of the remaining objects. Note that, therefore, we are also
    # working with a list of pixels, whether given by the user or generated by
    # us using SExtractor.

    if options.list:
        print "%sPhotometry will be done on the %d pixels listed in '%s'." % \
              (style.prefix, len(list_of_pixels), options.list)
    else:
        print "%sDetecting sources on the reference image..." % style.prefix ,
        sys.stdout.flush()
        reference_img = \
             seeing.FITSeeingImage(str(reference_img.path), options.maximum,
                                   options.margin, coaddk = options.coaddk)
        print 'done.'

        # Get the list of coordinates (Pixels) for the valid detections
        list_of_pixels = reference_img.coordinates

        ipercentage = reference_img.ignored / reference_img.total * 100
        rpercentage = len(reference_img) / reference_img.total * 100
        print "%s%d detections (%.2f %%) within %d pixels of the edge " \
              "were removed." % (style.prefix, reference_img.ignored,
                                 ipercentage, options.margin)
        print "%sThere remain %d sources (%.2f %%) on which to do " \
              "photometry." % (style.prefix, len(reference_img), rpercentage)

    print style.prefix
    print "%sPhotometry has to be done on the reference image in order to " \
          "extract" % style.prefix
    print "%sthe instrumental magnitude of each star, using the following " \
          "parameters:" % style.prefix

    # The photometry method always receives an XMLOffset instance, so we need
    # to create one to stupidly indicate that the offset between the reference
    # image and itself is... zero. Who would have guessed that?

    keys = ('path', 'path', 'object', 'filter', 'date', 'fwhm', 'airmass')
    args = [xml_offsets.reference[k] for k in keys]
    args += [0.0, 0.0, float('inf'), float('inf')]
    null_offset = xmlparse.XMLOffset(*args)

    if not fixed_annuli:

        reference_fwhm = xml_offsets.reference['fwhm']
        reference_aperture = options.aperture * reference_fwhm
        reference_annulus = options.annulus  * reference_fwhm
        reference_dannulus = options.dannulus * reference_fwhm

        print "%sFWHM (reference image) = %.3f pixels, therefore:" % \
              (style.prefix, reference_fwhm)
        print "%sAperture radius = %.3f x %.2f = %.3f pixels" % \
              (style.prefix, reference_fwhm, options.aperture, reference_aperture)
        print "%sSky annulus, inner radius = %.3f x %.2f = %.3f pixels" % \
              (style.prefix, reference_fwhm, options.annulus, reference_annulus)
        print "%sSky annulus, width = %.3f x %.2f = %.3f pixels" % \
              (style.prefix, reference_fwhm, options.dannulus, reference_dannulus)

        if reference_dannulus < options.min:
            reference_dannulus = options.min
            msg = ("%sWhoops! Sky annulus too thin, setting it to the minimum "
                   "of %.2f pixels") %  (style.prefix, reference_dannulus)
            warnings.warn(msg)

    else:
        reference_aperture = options.aperture_pix
        reference_annulus = options.annulus_pix
        reference_dannulus = options.dannulus_pix

        print "%sAperture radius = %.3f pixels" % \
              (style.prefix, reference_aperture)
        print "%sSky annulus, inner radius = %.3f pixels" % \
              (style.prefix, reference_annulus)
        print "%sSky annulus, width = %.3f pixels" % \
              (style.prefix, reference_dannulus)

    print "%sExtracting the instrumental magnitude of each star..." % style.prefix ,
    sys.stdout.flush()

    # The path to the uncalibrated, original image can not be used here, as for
    # the reference image, most probably being the result of combining multiple
    # FITS files, there is no such thing as a 'original' file.
    #
    # What is the saturation level of the reference image? We don't know, and
    # in any case it would depend on the saturation level of the images that
    # were combined. Anyway, it does not matter at all, as here we are only
    # detecting sources, saturated or not. Hence the usage of sys.maxint, which
    # returns the largest positive integer supported by the regular integer
    # type. It is not the largest integer supported by Python, but being at
    # least 2**31-1, as a saturation level it is quite close to infinity.

    qphot_result = qphot.run(null_offset, reference_aperture,
                             reference_annulus, reference_dannulus,
                             sys.maxint, options.exptimek, None,
                             list_of_pixels)

    # From the reference catalog, ignore those stars which are INDEF (None)
    #
    # Also: unless the --pixels option is given, the photometry function does
    # photometry on all the sources detected on the reference image. However,
    # this is not what we want if one or more stars have been discarded because
    # they were INDEF in the reference image. Solution: we use list_of_pixels
    # even if it was not given by the user, listing the coordinates of only the
    # non-INDEF stars, and pass it to 'photometry'. This not only solves the
    # problem; as a side effect, the execution time is sped up as sources do no
    # longer have to be detected (i.e., the SExtractor loaded from disk, as it
    # is cached on a file by LEMON) on the reference image over and over.

    reference_stars = [] # a tuple for each non-INDEF star
    ignored_counter = 0  # number of stars that are INDEF
    assert len(qphot_result) == len(list_of_pixels)

    for star_index, star_photometry in enumerate(qphot_result):
        imag = star_photometry.magnitude
        if imag is not None:  # NoneType returned for INDEF mags
            x = star_photometry.x
            y = star_photometry.y

            if options.list:
                ra, dec = None, None
            else:
                # RA / DEC extracted from sextractor catalog (that is why
                # we cannot currently determine them if --pixels is used)
                db_star = reference_img[star_index]
                ra, dec = db_star.alpha, db_star.delta

                if __debug__:
                    epsilon = 0.001  # a thousandth of a pixel!
                    assert abs(db_star.x - x) < epsilon
                    assert abs(db_star.y - y) < epsilon

            reference_stars.append((x, y, ra, dec, imag))
        else:
            ignored_counter += 1

    # Make sure no INDEF star made it to the list of stars with which we'll
    # work. The instrumental magnitude is the last element of each tuple
    assert all([x[-1] for x in reference_stars])
    assert len(reference_stars) + ignored_counter == \
           len(list_of_pixels or reference_img)
    print 'done.'

    if ignored_counter:
        print "%s%d stars were INDEF in the reference image and therefore " \
              "ignored." % (style.prefix, ignored_counter)

    if not reference_stars:
        print "%sError. There are no stars on which to do photometry." % \
              style.prefix
        print style.error_exit_message
        return 1
    if ignored_counter:
        print "%sThere are %d detections left on which to do photometry." % \
              (style.prefix, len(reference_stars))

    output_db = database.LEMONdB(output_db_path)

    # In the first place, we extract the image and sky coordinates of each
    # star, immediately storing this information in the database. But this is
    # only possible if sources were automatically detected on the intrumental
    # magnitudes image. Otherwise, we (currently, at least) have no way of
    # knowing neither the celestial coordinates that correspond to a pixel on
    # the image nor the magnitude reported by SExtractor for that pixel.
    #
    # This means that if sources were automatically detected, all of this
    # information is taken from the SExtractor catalog, while if a list of
    # pixels was given we only store the x- and y-image coordinates, leaving
    # alpha, delta and the reference magnitude set to minus one. This may seem
    # useless, but think that the pipeline will mainly (and even maybe only)
    # use the --pixels option when the light curve of a specific set of stars
    # has to be evaluated, and in that case we do not care (or, at least, we do
    # not currently need to know) about the declination, right ascension or
    # instrumental magnitude of these pixels in the reference image.

    print "%sStoring the information of each star in the database..." % style.prefix,
    sys.stdout.flush()
    for star_id, star_info in enumerate(reference_stars):
        x, y, ra, dec, imag = star_info
        output_db.add_star(star_id, x, y, ra, dec, imag)
    print 'done.'

    # Although not used for photometry, we also need to store the path to the
    # reference image so that the 'x' and 'y' coordinates in the STARS table
    # can be correlated to something.
    #
    # The date of the image may not make much sense in many cases (for example,
    # when the reference image is the output of the mosaic.py module, i.e., the
    # result of combining multiple images in order to maximize the SNR), but it
    # is still valid when a single image is used for sources detection.

    pfilter = xml_offsets.reference['filter']
    date = xml_offsets.reference['date']
    object_ = xml_offsets.reference['object']
    airmass = xml_offsets.reference['airmass']
    gain = options.gain or reference_img.read_keyword(options.gaink)
    args = reference_img.path, pfilter, date, object_, airmass, gain
    output_db.rimage = database.ReferenceImage(*args)

    # Insert FITS file (as blob) into the LEMONdB
    output_db.mosaic = reference_img.path

    # Determine how many different filters there are among the images contained
    # in the list of offsets. Then, sort the filters by their wavelength, so
    # that, e.g., photometry for the images observed in the Johnson V band (551
    # nanometers) is done before that for the images observed in Johnson I (806
    # nms), even although lexicographically they go the other way around.

    photometry_filters = sorted(set([x.filter for x in xml_offsets]))

    for pfilter in photometry_filters:
        band_offsets = [x for x in xml_offsets if x.filter == pfilter]

        print style.prefix
        print "%sPhotometry will now be done on the %d images taken in the " \
              "%s filter." % (style.prefix, len(band_offsets), pfilter)

        # The procedure if the dimensions of the aperture and sky annuli are to
        # be extracted from the --annuli file is simple: just take the first
        # CandidateAnnuli instance, as they are sorted in increasing order by
        # the standard deviation (which means that the best one is the first
        # element of the list) and use it.
        #
        # Alternatively, if the dimensions of the annuli are to be determined
        # by the median FWHM of the images, this has to be done for each
        # different filter in which images were taken. This contrasts with when
        # specific sizes (in pixels) are given for the annuli, which are used
        # for all the filters.

        if xml_annuli:

            # Store all the CandidateAnnuli objects in the LEMONdB
            assert len(xml_annuli[pfilter])
            for cand in xml_annuli[pfilter]:
                output_db.add_candidate_pparams(cand, pfilter)

            filter_annuli = xml_annuli[pfilter][0]
            aperture = filter_annuli.aperture
            annulus  = filter_annuli.annulus
            dannulus = filter_annuli.dannulus

            print "%sUsing the photometric parameters listed in the XML " \
                  "file, which are:" % style.prefix
            print "%sAperture radius = %.3f pixels" % \
                  (style.prefix, aperture)
            print "%sSky annulus, inner radius = %.3f pixels" % \
                  (style.prefix, annulus)
            print "%sSky annulus, width = %.3f pixels" % \
                  (style.prefix, dannulus)

        elif options.individual_fwhm:
            print "%sUsing photometric parameters derived from the FWHM of " \
                  "each image:" % style.prefix
            print "%sAperture radius = %.2f x FWHM pixels" % \
                  (style.prefix, options.aperture)
            print "%sSky annulus, inner radius = %.2f x FWHM pixels" % \
                  (style.prefix, options.annulus)
            print "%sSky annulus, width = %.2f x FWHM pixels" % \
                  (style.prefix, options.dannulus)

        elif not fixed_annuli:
            print "%sCalculating the median FWHM for this filter..." %  \
                  style.prefix ,
            sys.stdout.flush()

            img_fwhms = []
            for img in band_offsets:
                logging.debug("%s: FWHM = %.3f" % (img.shifted, img.fwhm))
                img_fwhms.append(img.fwhm)

            fwhm = numpy.median(img_fwhms)
            print 'done.'

            aperture = fwhm * options.aperture
            annulus  = fwhm * options.annulus
            dannulus = fwhm * options.dannulus

            print "%sFWHM (%s passband) = %.3f pixels, therefore:" % \
                  (style.prefix, pfilter, fwhm)
            print "%sAperture radius = %.3f x %.2f = %.3f pixels" % \
                  (style.prefix, fwhm, options.aperture, aperture)
            print "%sSky annulus, inner radius = %.3f x %.2f = %.3f pixels" % \
                  (style.prefix, fwhm, options.annulus, annulus)
            print "%sSky annulus, width = %.3f x %.2f = %.3f pixels" % \
                  (style.prefix, fwhm, options.dannulus, dannulus)

            if dannulus < options.min:
                dannulus = options.min
                msg = ("%sWhoops! Sky annulus too thin, setting it to the "
                      "minimum of %.2f pixels") % (style.prefix, dannulus)
                warnings.warn(msg)

        else: # fixed aperture and sky annuli directly specified in pixels
            aperture = options.aperture_pix
            annulus  = options.annulus_pix
            dannulus = options.dannulus_pix

            print "%sAperture radius = %.3f pixels" % \
                  (style.prefix, aperture)
            print "%sSky annulus, inner radius = %.3f pixels" % \
                  (style.prefix, annulus)
            print "%sSky annulus, width = %.3f pixels" % \
                  (style.prefix, dannulus)

        # The task of doing photometry on a series of images is inherently
        # parallelizable; use a pool of workers to which to assign the images.
        pool = multiprocessing.Pool(options.ncores)

        # Stars are stored in 'reference_stars' as database.DBStar instances.
        # However, the photometry function receives the pixels as a sequence
        # of astromatic.Pixel instances; hence the cast. Note that in
        # 'reference_stars' we have a five-element tuple, where the x and y
        # coordinates are the first two values.
        reference_pixels = [astromatic.Pixel(*x[:2]) for x in reference_stars]

        def fwhm_derived_params(offset):
            """ Return the FWHM-derived aperture and sky annuli parameters.

            Return a three-element tuple with (1) the aperture radius, (2) the
            sky annulus inner radius and (3) its width. These are equal to the
            FWHM of the image (offset.fwhm) times the --aperture, --annulus and
            --dannulus options, respectively.

            """

            fwhm = offset.fwhm

            aperture = fwhm * options.aperture
            annulus  = fwhm * options.annulus
            dannulus = fwhm * options.dannulus

            path = offset.shifted
            logging.debug("%s: FWHM = %.3f" % (path, fwhm))
            msg = "%s: FWHM-derived aperture: %.3f x %.2f = %.3f pixels"
            logging.debug(msg % (path, fwhm, options.aperture, aperture))
            msg = "%s: FWHM-derived annulus: %.3f x %.2f = %.3f pixels"
            logging.debug(msg % (path, fwhm, options.annulus, annulus))
            msg = "%s: FWHM-derived dannulus: %.3f x %.2f = %.3f pixels"
            logging.debug(msg % (path, fwhm, options.dannulus, dannulus))

            return aperture, annulus, dannulus

        # Define qphot_params either as a function that always returns the same
        # aperture, annulus and dannulus (since the same photometric parameters
        # are to be used for all the images in the same filter) *or*, if the
        # --individual-fwhm option was given, derives them from the FWHM of
        # each of the images. This allows us to, in both cases, populate
        # map_async_args by looping over the images on which photometry is to
        # be done and, for each one of them, calling qphot_params to get the
        # parameters.

        if not options.individual_fwhm:
            qphot_params = lambda x: (aperture, annulus, dannulus)
        else:
            qphot_params = fwhm_derived_params

        map_async_args = \
            ((offset, reference_pixels, qphot_params(offset), options)
             for offset in band_offsets)
        result = pool.map_async(parallel_photometry, map_async_args)

        methods.show_progress(0.0)
        while not result.ready():
            time.sleep(1)
            methods.show_progress(queue.qsize() / len(band_offsets) * 100)
            # Do not update the progress bar when debugging; instead, print it
            # on a new line each time. This prevents the next logging message,
            # if any, from being printed on the same line that the bar.
            if logging_level < logging.WARNING:
                print

        result.get() # reraise exceptions of the remote call, if any
        methods.show_progress(100) # in case the queue was ready too soon
        print

        # And now store the photometry for all the stars in the database
        print "%sStoring photometric information in the database..." % style.prefix
        methods.show_progress(0)
        photometric_results = (queue.get() for x in xrange(queue.qsize()))
        for index, (pimage, img_offset, qphot_result) in enumerate(photometric_results):

            logging.debug("Storing image %s in database" % pimage.path)
            output_db.add_image(pimage)
            logging.debug("Image %s successfully stored" % pimage.path)

            # Now store the photometry for each star
            for star_id, star_photometry in enumerate(qphot_result):

                if __debug__:
                    # Make sure that photometry was done on the correct pixels:
                    # apply the image offset to the original coordinates of the
                    # star (those in the reference image) and verify that this
                    # value equals the coordinates were photometry was done.
                    star_x, star_y = reference_stars[star_id][:2]
                    expected_x = star_x + img_offset.x
                    expected_y = star_y + img_offset.y

                    logging.debug("%s: star %d x-coord = %.3f (expected = %.3f)" % \
                                  (pimage.path, star_id, star_photometry.x, expected_x))
                    logging.debug("%s: star %d y-coord = %.3f (expected = %.3f)" % \
                                  (pimage.path, star_id, star_photometry.y, expected_y))

                    epsilon = 0.001  # a thousandth of a pixel!
                    assert abs(star_photometry.x - expected_x) < epsilon
                    assert abs(star_photometry.y - expected_y) < epsilon

                # INDEF stars will have a magnitude of None, and those with at
                # least one saturated pixel in the aperture will have infinity.
                # In both cases the measure is useless for our photometric
                # purposes and can be ignored -- it won't even be saved.
                magnitude = star_photometry.magnitude
                if magnitude is None:
                    logging.debug("%s: star %d is INDEF (None)" % \
                                  (pimage.path, star_id))
                    continue

                elif magnitude == float('infinity'):
                    logging.debug("%s: star %d is saturated (infinity)" % \
                                  (pimage.path, star_id))
                    continue
                else:
                    logging.debug("%s: star %d magnitude = %f" % \
                                  (pimage.path, star_id, magnitude))

                # Stars with a signal-to-noise ratio less than or equal to one
                # are ignored -- not only because these measures are anything
                # but reliable, but also because such values are outside of the
                # domain function that converts SNRs to errors in magnitudes.
                snr = star_photometry.snr(pimage.gain)
                if snr <= 1:
                    logging.debug("%s: star %d ignored (SNR = %f <= 1.0)" % \
                                  (pimage.path, star_id, snr))
                    continue
                else:
                    logging.debug("%s: star %d SNR = %f" % \
                                  (pimage.path, star_id, snr))
                    logging.debug("%s: storing photometry for star %d in database" % \
                                  (pimage.path, star_id))
                    args = star_id, pimage.unix_time, pimage.pfilter, magnitude, snr
                    output_db.add_photometry(*args)
                    logging.debug("%s: photometry for star %d successfully stored" % \
                                  (pimage.path, star_id))

            methods.show_progress(100 * (index + 1) / len(band_offsets))
            if logging_level < logging.WARNING:
                print

        else:
            logging.info("Photometry for %s completed" % pfilter)
            logging.debug("Committing database transaction")
            output_db.commit()
            logging.info("Database transaction commited")

            methods.show_progress(100.0)
            print

    # Collect information that can be used by the query optimizer to help make
    # better query planning choices. In the absence of ANALYZE information,
    # SQLite assumes that each table contains one million records when deciding
    # between doing a full table scan and constructing an automatic index.

    print "%sGathering statistics about tables and indexes..." % style.prefix ,
    sys.stdout.flush()
    output_db.analyze()
    print 'done.'

    # Store into the METADATA table of the LEMONdB the current time (in seconds
    # since the Unix epoch), the login name of the currently effective user id
    # and the hostname of the machine where Python is currently executing.

    output_db.date = time.time()
    output_db.author = pwd.getpwuid(os.getuid())[0]
    output_db.hostname = socket.gethostname()
    output_db.commit()

    # Use as unique identifier of the LEMONdB a 32-digit hexadecimal number:
    # the MD5 hash of the concatenation, in this order, of the LEMONdB.date,
    # author and hostname properties, which we have just set above.

    md5 = hashlib.md5()
    md5.update(str(output_db.date))
    md5.update(str(output_db.author))
    md5.update(str(output_db.hostname))
    output_db.id = md5.hexdigest()
    output_db.commit()

    methods.owner_writable(output_db_path, False) # chmod u-w
    print "%sYou're done ^_^" % style.prefix
    return 0

if __name__ == "__main__":
    sys.exit(main())

