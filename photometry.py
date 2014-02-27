#! /usr/bin/env python
# -*- coding: utf-8 -*-

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

import atexit
import collections
import hashlib
import itertools
import logging
import multiprocessing
import numpy
import optparse
import os
import os.path
import pwd
import socket
import sys
import tempfile
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

# Message of the warning that is issued when the width of the sky annulus is
# smaller than the number of pixels specified with the --min-sky option.
DANNULUS_TOO_THIN_MSG = \
"Whoops! Sky annulus too thin, setting it to the minimum of %.2f pixels"

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

def clean_tmp_coords_file(path):
    """ Try to unlink the coordinates file 'path'. """

    msg = "Cleaning up temporary file '%s'"
    logging.debug(msg % path)

    try:
        os.unlink(path)
    except OSError, e:
        msg = "Cannot delete '%s' (%s)"
        logging.debug(msg % (path, e))
    else:
        msg = "Temporary coordinates file '%s' removed"
        logging.debug(msg % path)

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

coords_group = optparse.OptionGroup(parser, "List of Coordinates",
               "By default, we run SExtractor on the sources image in order "
               "to detect the astronomical objects on which to do photometry. "
               "Alternatively, with this option it is possible to skip the "
               "detection stage and directly do photometry on the objects "
               "whose celestial coordinates are specified in a text file. "
               "Note that this option is exclusive with the automatic "
               "detection of sources, so photometry will be done *only* on "
               "these coordinates.")

coords_group.add_option('--coordinates', action = 'store', type = str,
                        dest = 'coordinates', default = None,
                        help = "path to the file containing the celestial "
                        "coordinates of the objects on which photometry must "
                        "be done. These coordinates must be given in degrees "
                        "and listed one per line, in two columns, right "
                        "ascension and declination, respectively.")
parser.add_option_group(coords_group)

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

key_group.add_option('--objectk', action = 'store', type = 'str',
                     dest = 'objectk', default = keywords.objectk,
                     help = keywords.desc['objectk'])

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

key_group.add_option('--saturk', action = 'store', type = 'str',
                     dest = 'saturk', default = keywords.saturk,
                     help = keywords.desc['saturk'])

key_group.add_option('--gaink', action = 'store', type = 'str',
                     dest = 'gaink', default = keywords.gaink,
                     help = keywords.desc['gaink'])

key_group.add_option('--fwhmk', action = 'store', type = 'str',
                     dest = 'fwhmk', default = keywords.fwhmk,
                     help = keywords.desc['fwhmk'])

key_group.add_option('--airmk', action = 'store', type = 'str',
                     dest = 'airmassk', default = keywords.airmassk,
                     help = keywords.desc['airmassk'])

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

    # Abort the execution if one or more of the {aperture,{,d}annulus}-pix:
    # options are given, but not all three. If we are going to use fixed sizes
    # for the aperture and sky annuli we need to know all of them.
    fixed_pix_count = bool(options.aperture_pix) + \
                      bool(options.annulus_pix)  + \
                      bool(options.dannulus_pix)

    if fixed_pix_count:
        if fixed_pix_count < 3:
            assert 1 <= fixed_pix_count <= 2
            msg = "%sError. The --aperture-pix, --annulus-pix and " + \
                  "--dannulus-pix options must be used together."
            print msg % style.prefix
            print style.error_exit_message
            return 1
        else:
            assert fixed_pix_count == 3
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

    # If the --coordinates option has been given, read the text file and
    # extract the right ascensions and declinations as two-element tuples,
    # storing them in a list of astromatic.Coordinates objects. Abort the
    # execution if the coordinates file is empty.

    if options.coordinates:

        sources_coordinates = []
        for ra, dec in methods.load_coordinates(options.coordinates):
            coords = astromatic.Coordinates(ra, dec)
            sources_coordinates.append(coords)

        if not sources_coordinates:
            msg = "%sError. Coordinates file '%s' is empty."
            print msg % (style.prefix, options.coordinates)
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

    msg = "%s%d different photometric filters were detected:"
    print msg % (style.prefix, len(files.keys()))

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

    print "%sSources image: %s" % (style.prefix, sources_img_path)
    print "%sRunning SExtractor on the sources image..." % style.prefix ,
    sys.stdout.flush()
    args = (sources_img_path, options.maximum, options.margin)
    kwargs = dict(coaddk = options.coaddk,
                  saturk = options.saturk)
    sources_img = seeing.FITSeeingImage(*args, **kwargs)
    print 'done.'

    msg = "%sCalculating coordinates of field center..."
    print msg % style.prefix ,
    sys.stdout.flush()

    ra, dec = sources_img.center_wcs()
    sources_img_ra = ra
    sources_img_dec = dec
    print 'done.'

    # Print coordinates, in degrees and sexagesimal
    print "%sα = %11.7f" % (style.prefix, sources_img_ra) ,
    msg = " (%.02d %.02d %05.2f)"
    args = methods.DD_to_HMS(sources_img_ra)
    print msg % args

    print "%sδ = %11.7f" % (style.prefix, sources_img_dec) ,
    msg = "(%+.02d %.02d %05.2f)"
    args = methods.DD_to_DMS(sources_img_dec)
    print msg % args

    # If --coordinates was given, let the user know on how many celestial
    # coordinates we are going to do photometry. If not, run SExtractor on the
    # sources image, discard those detections too close to the edges and create
    # a list of Coordinates objects with the right ascension and declination of
    # the remaining detections. Note that internally we always work with a list
    # of coordinates, whether given by the user or generated by us.

    if options.coordinates:
        msg = "%sPhotometry will be done on the %d coordinates listed in '%s'."
        args = (style.prefix, len(sources_coordinates), options.coordinates)
        print msg % args

    else:
        sources_coordinates = sources_img.coordinates
        assert len(sources_coordinates) == len(sources_img)
        ipercentage = sources_img.ignored / sources_img.total * 100
        rpercentage = len(sources_img) / sources_img.total * 100

        msg = "%s%d detections (%.2f %%) within %d pixels of the edge were removed."
        print msg % (style.prefix, sources_img.ignored, ipercentage, options.margin)
        msg = "%sThere remain %d sources (%.2f %%) on which to do photometry."
        print msg % (style.prefix, len(sources_img), rpercentage)

        # Save the coordinates of the astronomical objects detected by
        # SExtractor to a text file, needed to do photometry later on. Note
        # that we store the absolute pathname of the temporary file in the
        # options.coordinates attribute, so that the coordinates file can be
        # always accessed there, independently of whether the user gave the
        # --coordinates option or not.

        kwargs = dict(prefix = '%s_' % sources_img.basename_woe,
                      suffix = '_sextractor.coords',
                      text = True)

        coords_fd, options.coordinates = tempfile.mkstemp(**kwargs)
        atexit.register(clean_tmp_coords_file, options.coordinates)
        for ra, dec in sources_coordinates:
            os.write(coords_fd, "%.10f\t%.10f\n" % (ra, dec))
        os.close(coords_fd)

    print style.prefix
    msg = "%sNeed to determine the instrumental magnitude of each source."
    print msg % style.prefix
    msg = "%sDoing photometry on the sources image, using the parameters:"
    print msg % style.prefix

    # Unless the photometric parameters are given in pixels, the sizes of the
    # aperture and sky annulus are determined by the FWHM of the sources image.
    if not fixed_annuli:

        sources_img_fwhm = sources_img.read_keyword(options.fwhmk)
        sources_aperture = options.aperture * sources_img_fwhm
        sources_annulus  = options.annulus  * sources_img_fwhm
        sources_dannulus = options.dannulus * sources_img_fwhm

        t = (style.prefix, sources_img_fwhm)
        msg = "%sFWHM (sources image) = %.3f pixels, therefore:"
        print msg % t
        msg = "%sAperture radius = %.3f x %.2f = %.3f pixels"
        print msg % (t + (options.aperture, sources_aperture))
        msg = "%sSky annulus, inner radius = %.3f x %.2f = %.3f pixels"
        print msg % (t + (options.annulus, sources_annulus))
        msg = "%sSky annulus, width = %.3f x %.2f = %.3f pixels"
        print msg % (t + (options.dannulus, sources_dannulus))

        if sources_dannulus < options.min:
            sources_dannulus = options.min
            msg = style.prefix + DANNULUS_TOO_THIN_MSG
            warnings.warn(msg % sources_dannulus)

    else:
        sources_aperture = options.aperture_pix
        sources_annulus  = options.annulus_pix
        sources_dannulus = options.dannulus_pix

        msg = "%sAperture radius = %.3f pixels"
        print msg % (style.prefix, sources_aperture)
        msg = "%sSky annulus, inner radius = %.3f pixels"
        print msg % (style.prefix, sources_annulus)
        msg = "%sSky annulus, width = %.3f pixels"
        print msg % (style.prefix, sources_dannulus)

    print style.prefix
    msg = "%sRunning IRAF's qphot..."
    print msg % style.prefix ,
    sys.stdout.flush()

    # Some (or even many) astronomical objects may be saturated in the sources
    # image, but (a) there is nothing we can really do about it and, anyway,
    # (b) this fact is irrelevant for our purposes. The instrumental magnitude
    # computed by IRAF's qphot in the sources image is exclusively intended to
    # serve as a very rough estimate of how bright each object is, allowing us
    # to compare its intensity to that of other objects, but nothing more.
    # Because of their saturation, there is no guarantee that the instrumental
    # magnitudes of the brightest objects will be the right ones: they may
    # appear less bright than they actually are, we hypothesize that following
    # a non-linear distribution.
    #
    # The number of ADUs at which saturation arises must be sufficiently large
    # so that qphot.run() does not mark any object as saturated. An approach
    # could be using float('infinity'), but the function expects an integer.
    # That is why we instead use sys.maxint, which returns the largest positive
    # integer supported by the regular integer type. Being at least 2 ** 31 -
    # 1, as a saturation level this value is sufficiently close to infinity.

    qphot_args = [sources_img, options.coordinates,
                  sources_aperture, sources_annulus, sources_dannulus,
                  sys.maxint, options.exptimek, None]
    sources_phot = qphot.run(*qphot_args)
    print 'done.'

    # Remove those astronomical objects so faint that they are INDEF in the
    # sources image. After all, if they are not even visible in this image,
    # which ideally should be as deep as possible, they will not be visible in
    # the individual images either. This may happen, for example, with false
    # positive detections by SExtractor, or if incorrect coordinates, that do
    # not correspond to any object, are given with the --coordinates option.
    #
    # Make a copy of the options.coordinates file and delete the lines of the
    # objects that are INDEF (i.e., whose magnitude is None). This is possible
    # because the order of the QPhotResult objects in the QPhot object returned
    # by qphot.run() is guaranteed to preserve that of the astronomical objects
    # listed in the coordinates file.

    msg = "%sDetecting INDEF objects..."
    print msg % style.prefix ,
    sys.stdout.flush()

    with open(options.coordinates, 'rt') as fd:
        lines = fd.readlines()

    kwargs = dict(prefix = 'PID_%d' % os.getpid(),
                  suffix = '_NO_INDEF.coords',
                  text = True)

    # A second coordinates file, listing only non-INDEF objects
    coords_fd, options.coordinates = tempfile.mkstemp(**kwargs)
    atexit.register(clean_tmp_coords_file, options.coordinates)

    ignored_counter = 0
    non_ignored_counter = 0
    original_size = len(sources_phot)

    for line, object_phot in itertools.izip(lines, sources_phot):
        if object_phot.mag is not None:
            os.write(coords_fd, line)
            non_ignored_counter += 1
        else:
            ignored_counter += 1
    os.close(coords_fd)

    # Delete INDEF photometric measurements, in-place
    for index in xrange(len(sources_phot) - 1, -1, -1):
        if sources_phot[index].mag is None:
            sources_phot.pop(index)

    assert non_ignored_counter == len(sources_phot)
    assert ignored_counter + non_ignored_counter == original_size
    print 'done.'

    if ignored_counter:
        msg = "%s%s objects" % (style.prefix, ignored_counter)
    else:
        msg = "%sNo objects" % style.prefix
    print msg + " are INDEF in the sources image."

    if not non_ignored_counter:
        msg = "%sError. There are no objects left on which to do photometry."
        print msg % style.prefix
        print style.error_exit_message
        return 1

    elif ignored_counter:
        msg = "%sThere are %d objects left on which to do photometry."
        print msg % (style.prefix, len(sources_phot))

    if __debug__:

        msg = "%sMaking sure INDEF objects were removed..."
        print msg % style.prefix ,
        sys.stdout.flush()

        # Do photometry again, use the non-INDEF coordinates file
        qphot_args[1] = options.coordinates
        non_INDEF_phot = qphot.run(*qphot_args)
        assert sources_phot == non_INDEF_phot
        print 'done.'

    print style.prefix
    msg = "%sInitializing output LEMONdB..."
    print msg % style.prefix ,
    sys.stdout.flush()

    output_db = database.LEMONdB(output_db_path)

    # The fact that the QPhot object returned by qphot.run() preserves the
    # order of the astronomical objects listed in options.coordinates proves to
    # be useful again: it allows us to match the celestial coordinates of each
    # object (a row in the coordinates file) to the corresponding QPhotResult
    # object. Note that qphot.run() accepts celestial coordinates but returns
    # the x- and y-coordinates of their centers, as IRAF's qphot does.

    sources_coords = list(methods.load_coordinates(options.coordinates))
    assert len(sources_coords) == len(sources_phot)
    it = itertools.izip(sources_coords, sources_phot)
    for id_, (object_coords, object_phot) in enumerate(it):
        x, y = object_phot.x, object_phot.y
        ra, dec = object_coords
        imag = object_phot.mag

        args = (id_, x, y, ra, dec, imag)
        output_db.add_star(*args)

    output_db.commit()
    print 'done.'

    # Store some relevant information about the sources image in the LEMONdB.
    # Do this by creating a database.Image object, which encapsulates a FITS
    # file, and assign it to the LEMONdB.simage attribute. The image is also
    # stored as a blob and is available through the LEMONdB.mosaic attribute.

    path = sources_img.path
    pfilter = sources_img.pfilter(options.filterk)

    kwargs = dict(date_keyword = options.datek,
                  time_keyword = options.timek,
                  exp_keyword = options.exptimek)
    unix_time = sources_img.date(**kwargs)

    object_ = sources_img.read_keyword(options.objectk)
    airmass = sources_img.read_keyword(options.airmassk)
    # If not given with --gaink, read it from the FITS header
    gain = options.gain or sources_img.read_keyword(options.gaink)
    ra, dec = sources_img_ra, sources_img_dec

    args = (path, pfilter, unix_time, object_, airmass, gain, ra, dec)
    simage = database.Image(*args)
    output_db.simage = simage
    output_db.commit()

    for pfilter, images in files.iteritems():
        print style.prefix
        msg = "%sLet's do photometry on the %d images taken in the %s filter."
        args = (style.prefix, len(images), pfilter)
        print msg % args

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

            msg = "%sUsing the parameters listed in the XML file, which are:"
            print msg % style.prefix
            msg = "%sAperture radius = %.3f pixels"
            print msg % (style.prefix, aperture)
            msg = "%sSky annulus, inner radius = %.3f pixels"
            print msg % (style.prefix, annulus)
            msg = "%sSky annulus, width = %.3f pixels"
            print msg % (style.prefix, dannulus)

        elif options.individual_fwhm:
            msg = "%sUsing parameters derived from the FWHM of each image:"
            print msg % style.prefix
            msg = "%sAperture radius = %.2f x FWHM pixels"
            print msg % (style.prefix, options.aperture)
            msg = "%sSky annulus, inner radius = %.2f x FWHM pixels"
            print msg % (style.prefix, options.annulus)
            msg = "%sSky annulus, width = %.2f x FWHM pixels"
            print msg % (style.prefix, options.dannulus)

        elif not fixed_annuli:
            msg = "%sCalculating the median FWHM for this filter..."
            print msg % style.prefix ,
            sys.stdout.flush()

            pfilter_fwhms = []
            for path in images:
                img = fitsimage.FITSImage(path)
                img_fwhm = img.read_keyword(options.fwhmk)
                logging.debug("%s: FWHM = %.3f" % (img.path, img_fwhm))
                pfilter_fwhms.append(img_fwhm)

            fwhm = numpy.median(pfilter_fwhms)
            print 'done.'

            aperture = fwhm * options.aperture
            annulus  = fwhm * options.annulus
            dannulus = fwhm * options.dannulus

            msg = "%sFWHM (%s) = %.3f pixels, therefore:"
            print msg % (style.prefix, pfilter, fwhm)
            msg = "%sAperture radius = %.3f x %.2f = %.3f pixels"
            print msg % (style.prefix, fwhm, options.aperture, aperture)
            msg = "%sSky annulus, inner radius = %.3f x %.2f = %.3f pixels"
            print msg % (style.prefix, fwhm, options.annulus, annulus)
            msg = "%sSky annulus, width = %.3f x %.2f = %.3f pixels"
            print msg % (style.prefix, fwhm, options.dannulus, dannulus)

            if dannulus < options.min:
                dannulus = options.min
                msg = style.prefix + DANNULUS_TOO_THIN_MSG
                warnings.warn(msg % dannulus)

        else: # fixed aperture and sky annuli directly specified in pixels
            aperture = options.aperture_pix
            annulus  = options.annulus_pix
            dannulus = options.dannulus_pix

            msg = "%sAperture radius = %.3f pixels"
            print msg % (style.prefix, aperture)
            msg = "%sSky annulus, inner radius = %.3f pixels"
            print msg % (style.prefix, annulus)
            msg = "%sSky annulus, width = %.3f pixels"
            print msg % (style.prefix, dannulus)

        # The task of doing photometry on a series of images is inherently
        # parallelizable; use a pool of workers to which to assign the images.
        pool = multiprocessing.Pool(options.ncores)

        def fwhm_derived_params(img):
            """ Return the FWHM-derived aperture and sky annuli parameters.

            Return a database.PhotometricParameters object (a three-element
            named tuple) containing (1) the aperture radius, (2) sky annulus
            inner radius and (3) its width, in pixels, which with to do
            photometry. These are equal to the FWHM of the FITS file (a
            fitsimage.FITSImage object) times the --aperture, --annulus
            and --dannulus options, respectively.

            """

            fwhm = img.read_keyword(options.fwhmk)
            aperture = fwhm * options.aperture
            annulus  = fwhm * options.annulus
            dannulus = fwhm * options.dannulus

            path = img.path
            logging.debug("%s: FWHM = %.3f" % (path, fwhm))
            msg = "%s: FWHM-derived aperture: %.3f x %.2f = %.3f pixels"
            logging.debug(msg % (path, fwhm, options.aperture, aperture))
            msg = "%s: FWHM-derived annulus: %.3f x %.2f = %.3f pixels"
            logging.debug(msg % (path, fwhm, options.annulus, annulus))
            msg = "%s: FWHM-derived dannulus: %.3f x %.2f = %.3f pixels"
            logging.debug(msg % (path, fwhm, options.dannulus, dannulus))

            args = aperture, annulus, dannulus
            return database.PhotometricParameters(*args)

        # Define qphot_params either as a function that always returns the same
        # PhotometricParameters object (since identical photometric parameters
        # are to be used for all the images in this photometric filter) or, if
        # the --individual-fwhm option was used, derives them from the FWHM of
        # each of the FITS images. This allows us to, in both cases, populate
        # map_async_args by looping over the images on which photometry is to
        # be done and, for each one of them, calling qphot_params to get the
        # parameters that have to be used.

        if not options.individual_fwhm:
            args = aperture, annulus, dannulus
            pparams = database.PhotometricParameters(*args)
            qphot_params = lambda x: pparams
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

