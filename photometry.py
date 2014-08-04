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
This module does aperture photometry on all the FITS images that it receives as
arguments. Astronomical objects are automatically detected, using SExtractor,
on the first image (which will be referred to from now on as the 'sources
image'), and then photometry is done for their celestial coordinates on the
rest of the images. The output is a LEMON database. The sources image is also
used to extract the instrumental magnitude that for each astronomical object
is stored in the LEMONdB, in order to allow us to approximately estimate how
bright each object is.

By default, the sizes of the aperture and sky annulus are determined by the
median FWHM of the images in each photometric filter, but different options
make it possible, among others, to directly specify those sizes in pixels, or
to use apertures and sky annuli that depend on the FWHM of each individual
image. The celestial coordinates of the objects on which to do photometry,
if known, can be listed in a text file, in this manner skipping the sources
detection step and working exclusively with the specified objects.

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
import json_parse
import keywords
import methods
import qphot
import seeing
import style

# Message of the warning that is issued when the width of the sky annulus is
# smaller than the number of pixels specified with the --min-sky option.
DANNULUS_TOO_THIN_MSG = \
"Whoops! Sky annulus too thin, setting it to the minimum of %.2f pixels"

# The Queue is global -- this works, but note that we could have
# passed its reference to the function managed by pool.map_async.
# See http://stackoverflow.com/a/3217427/184363
queue = methods.Queue()

def get_fwhm(img, options):
    """ Return the FWHM of the FITS image.

    Attempt to read the full width at half maximum from the header of the FITS
    image (keyword options.fwhmk). If the keyword cannot be found, then compute
    the FWHM by calling FITSeeingImage.fwhm(). In this manner, we can always
    call this method to get the FWHM of each image, without having to worry
    about whether it is in the header already. The 'img' argument must be a
    fitsimage.FITSImage object, while 'options' must be the optparse.Values
    object returned by optparse.OptionParser.parse_args().

    """

    try:
        msg = "%s: reading FWHM from keyword '%s'"
        args = img.path, options.fwhmk
        logging.debug(msg % args)

        fwhm = img.read_keyword(options.fwhmk)

        msg = "%s: FWHM = %.3f (keyword '%s')"
        args = img.path, fwhm, options.fwhmk
        logging.debug(msg % args)
        return fwhm

    except KeyError:

        msg = "%s: keyword '%s' not found in header"
        args = img.path, options.fwhmk
        logging.debug(msg % args)

        if not isinstance(img, seeing.FITSeeingImage):

            msg = "%s: type of argument 'img' is not FITSeeingImage ('%s')"
            args = img.path, type(img)
            logging.debug(msg % args)

            msg = "%s: calling FITSeeingImage.__init__() with 'img'"
            logging.debug(msg % img.path)

            args = (img.path, options.maximum, options.margin)
            kwargs = dict(coaddk = options.coaddk)
            img = seeing.FITSeeingImage(*args, **kwargs)

        msg = "%s: calling FITSeeingImage.fwhm() to compute FWHM"
        logging.debug(msg % img.path)

        mode = 'mean' if options.mean else 'median'
        kwargs = dict(per = options.per, mode = mode)
        fwhm = img.fwhm(**kwargs)

        msg = "%s: FITSeeingImage.fwhm() returned %.3f"
        args = img.path, fwhm
        logging.debug(msg % args)
        return fwhm

@methods.print_exception_traceback
def parallel_photometry(args):
    """ Function argument of map_async() to do photometry in parallel.

    This will be the first argument passed to multiprocessing.Pool.map_async(),
    which chops the iterable into a number of chunks that are submitted to the
    process pool as separate tasks. 'args' must be a three-element tuple with
    (1) a fitsimage.FITSImage object, (2) a database.PhotometricParameters
    object and (3) 'options', the optparse.Values object returned by
    optparse.OptionParser.parse_args().

    This function does photometry (qphot.run()) on the astronomical objects of
    the FITS image listed in options.coordinates, using the aperture, annulus
    and dannulus defined by the PhotometricParameters object. The result is
    another three-element tuple, which is put into the module-level 'queue'
    object, a process shared queue. This tuple contains (1) a database.Image
    object, (2) a database.PhotometricParameters object and (3) a qphot.QPhot
    object -- therefore mapping each FITS file and the parameters used for
    photometry to the measurements returned by qphot.

    """

    image, pparams, options = args

    logging.debug("Doing photometry on %s" % image.path)
    msg = "%s: qphot aperture: %.3f"
    logging.debug(msg % (image.path, pparams.aperture))
    msg = "%s: qphot annulus: %.3f"
    logging.debug(msg % (image.path, pparams.annulus))
    msg = "%s: qphot dannulus: %.3f"
    logging.debug(msg % (image.path, pparams.dannulus))

    maximum = image.saturation(options.maximum, coaddk = options.coaddk)
    msg = "%s: saturation level = %d ADUs'"
    args = (image.path, maximum)
    logging.debug(msg % args)

    logging.info("Running qphot on %s" % image.path)
    args = (image, options.coordinates,
            pparams.aperture, pparams.annulus, pparams.dannulus,
            maximum, options.exptimek, options.uncimgk)
    img_qphot = qphot.run(*args)
    logging.info("Finished running qphot on %s" % image.path)

    msg = "%s: qphot.run() returned %d records"
    args = (image.path, len(img_qphot))
    logging.debug(msg % args)

    pfilter = image.pfilter(options.filterk)
    logging.debug("%s: filter = %s" % (image.path, pfilter))

    kwargs = dict(date_keyword = options.datek,
                  time_keyword = options.timek,
                  exp_keyword = options.exptimek)
    unix_time = image.date(**kwargs)
    msg = "%s: observation date: %.2f (%s)"
    args = (image.path, unix_time, methods.utctime(unix_time))
    logging.debug(msg % args)

    object_ = image.read_keyword(options.objectk)
    logging.debug("%s: object = %s" % (image.path, object_))

    airmass = image.read_keyword(options.airmassk)
    logging.debug("%s: airmass = %.4f" % (image.path, airmass))

    # If not given with --gaink, read it from the FITS header
    gain = options.gain or image.read_keyword(options.gaink)
    msg = "%s: gain (%s) = %.4f"
    gain_msg = "given by user" if options.gain else "read from header"
    args = image.path, gain_msg, gain
    logging.debug(msg % args)

    msg = "%s: calculating coordinates of field center"
    logging.debug(msg % image.path)
    ra, dec = image.center_wcs()
    logging.debug("%s: RA (field center) = %.8f" % (image.path, ra))
    logging.debug("%s: DEC (field center) = %.8f" % (image.path, dec))

    args = (image.path, pfilter, unix_time, object_, airmass, gain, ra, dec)
    db_image = database.Image(*args)
    queue.put((db_image, pparams, img_qphot))
    msg = "%s: photometry result put into global queue"
    logging.debug(msg % image.path)


parser = customparser.get_parser(description)
parser.usage = "%prog [OPTION]... SOURCES_IMG INPUT_IMGS... OUTPUT_DB"

parser.add_option('--overwrite', action = 'store_true', dest = 'overwrite',
                  help = "overwrite output database if it already exists")

parser.add_option('--filter', action = 'store', type = 'passband',
                  dest = 'filter', default = None,
                  help = "do not do photometry on all the FITS files given "
                  "as input, but only on those taken in this photometric "
                  "filter. " + defaults.desc['filter'])

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
                  dest = 'json_annuli', default = None,
                  help = "ignore the Aperture Photometry (FWHM and Pixels) "
                  "sections below an instead read the apertures and sky annuli "
                  "to use from a JSON file output by the 'annuli' command. "
                  "This file should, of course, have been generated for the "
                  "same set of images on which photometry is now being done. "
                  "For this same reason, the execution will be aborted if the "
                  "JSON file does not have information for the photometric "
                  "filters of all the input FITS files. Even when this option "
                  "is used, the aperture and sky annulus used for the sources "
                  "image are determined by the 'Aperture Photometry' sections.")

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
              "the values of the following parameters are defined in terms "
              "of the *median* FWHM of the images in each filter.\n\n"

              "The full width at half-maximum of an image is defined here as "
              "the median of the FWHM of all the astronomical objects detected "
              "by SExtractor. Therefore, you may think of the median 'FWHM of "
              "the images' as the 'median of the median of the FWHM of all the "
              "astronomical objects, if you find it easier. Note also that a "
              "different median FWHM is computed for each photometric filter "
              "on which photometry is done.")

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
                       "small FWHMs from resulting in too thin a sky "
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
              "In case the exact sizes of the aperture and sky annulus are "
              "known, their dimensions can be specified in pixels. If used, "
              "these three options must be used together.\n\n"

              "Note that setting the aperture photometry parameters to a "
              "fixed size, as these options do, means that the same values "
              "are used regardless of the photometric filter of the images, "
              "including the sources image. You probably want to use these "
              "options in conjunction with --filter, in order to do "
              "photometry one photometric filter at a time.")

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

fwhm_group = optparse.OptionGroup(parser, "FWHM",
              "The full width at half maximum of each FITS image is written "
              "to its header by the 'seeing' command, which is generally run "
              "before photometry is done. However, it may be the case that "
              "not all the images with which we have to work now are the "
              "output of 'seeing' (for example, if the sources image is that "
              "output by the 'mosaic' command). In these cases, we need to "
              "compute the FWHM of these images, in the same way the "
              "'seeing' command does.")

fwhm_group.add_option('--snr-percentile', action = 'store', type = 'float',
                      dest = 'per', default = defaults.snr_percentile,
                      help = defaults.desc['snr_percentile'])

fwhm_group.add_option('--mean', action = 'store_true', dest = 'mean',
                      help = defaults.desc['mean'])
parser.add_option_group(fwhm_group)

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

    # Print the help and abort the execution if there are not three positional
    # arguments left after parsing the options, as the user must specify the
    # sources image, at least one (only one?) image on which to do photometry
    # and the output LEMON database.
    if len(args) < 3:
        parser.print_help()
        return 2     # 2 is generally used for command line syntax errors
    else:
        sources_img_path = args[0]
        input_paths = set(args[1:-1])
        output_db_path = args[-1]

    assert input_paths

    # If the user gives an empty string as the FITS keyword which stores the
    # path to the original image, it is understood as meaning that we want
    # saturation to be checked for in the same images on which photometry is
    # done. If that is the case, we need to set the option to None (although an
    # empty string would also work), as that is what qphot.run expects to
    # receive in these cases.

    if not options.uncimgk:
        options.uncimgk = None

    # The annuli JSON file (that generated by the annuli command, and specified
    # with the --annuli option) must exist. The use of this file automatically
    # discards whathever was specified with the Aperture Photometry (FWHM and
    # pixels) options.
    json_annuli = None
    fixed_annuli = False

    if options.json_annuli:
        if not os.path.exists(options.json_annuli):
            print "%sError. The file '%s' does not exist." % \
                  (style.prefix, options.json_annuli)
            print style.error_exit_message
            return 1
        else:
            json_annuli = json_parse.CandidateAnnuli.load(options.json_annuli)
            msg = "%sPhotometric parameters read from the '%s' file."
            args = (style.prefix, os.path.basename(options.json_annuli))
            print msg % args

    # Even if the annuli JSON file is used, the aperture and sky annuli
    # parameters (whether FWHM-based or given in pixels) are still needed in
    # order to do photometry on the sources image and extract the instrumental
    # magnitude that for each astronomical object is stored in the database.

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
    # JSON file. Does not compute.
    if fixed_annuli and json_annuli:
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
        # from the JSON file, they cannot also depend on the FWHM of the images.
        if json_annuli:
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
    # make much sense to merge data (since the same tables would be used) of
    # astronomical objects that belong to different fields. Thus, we refuse to
    # work with an existing database (which is what the LEMONdB class would do
    # otherwise) unless the --overwrite option is given, in which case it is
    # deleted and created again from scratch.

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

    files = fitsimage.InputFITSFiles()
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

    for pfilter, images in sorted(files.iteritems()):
        msg = "%s %s: %d files (%.2f %%)"
        percentage = len(images) / len(files) * 100
        print msg % (style.prefix, pfilter, len(images), percentage)

    # Light curves, which are our ultimate goal, can only have one magnitude
    # for each point in time. Therefore, we cannot do photometry on two or more
    # images with the same observation date and photometric filter. This may
    # seem (and, indeed, is) unlikely to happen, but astronomical instruments
    # also have software errors - we have already come across this while
    # reducing images taken with Omega 2000, a camera for the 3.5m CAHA.
    #
    # We might be tempted to keep one of them, such as, for example, that with
    # the highest number of sources, or the best FWHM. However, the safest bet
    # is to discard them all, because we cannot know which image has the right
    # observation date. The only certain thing is that an error occurred. Thus,
    # we better forget about these images.

    msg = "%sMaking sure there are no images with the same date and filter..."
    print msg % style.prefix ,
    sys.stdout.flush()

    # Two-level dictionary: map each date of observation (UTC), in Unix time,
    # to a photometric filter to a list of the corresponding FITS files. If
    # there are no two or more images with the same date and filter, all the
    # second-level values of the dictionary will have a length of one. We do
    # not use the Counter class, from the collections module, for Python 2.6
    # compatibility.

    get_dict = lambda: collections.defaultdict(list)
    dates_counter = collections.defaultdict(get_dict)

    # 'files' maps each filter to a list of images, while 'img_dates' maps each
    # image to its Unix date. There is no need, therefore, to read the date or
    # photometric filter of the images from their FITS headers again.

    for pfilter, images in files.iteritems():
        for img_path in images:
            date = img_dates[img_path]
            dates_counter[date][pfilter].append(img_path)

    # Find the dates and filters for which there is more than one FITS file.
    # Then, remove from the InputFITSFiles object each of these images with
    # a duplicate Unix time and photometric filter.

    discarded = 0
    for date, date_images in dates_counter.items():
        for pfilter, images in date_images.items():

            if len(images) > 1:

                # "Making sure..." message above does not include newline.
                # We need to print it, but only for the first issued warning.
                if not discarded:
                    print

                msg = "%sWarning! Multiple images have date %s and filter %s"
                args = style.prefix, methods.utctime(date), pfilter
                warnings.warn(msg % args)
                del dates_counter[date][pfilter]
                for img in images:
                    discarded += files.remove(img)

    if not discarded:
        print 'done.'

    else:

        # There should be no FITS files with the same observation date and
        # filter anymore, and at least two of them should have been discarded.
        # Otherwise, how did we get to the 'else' clause in the first place?
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

    # If a JSON file is specified with --annuli, it must list the photometric
    # parameters for all the filters on which photometry is to be done.
    if json_annuli:
        for pfilter in files.iterkeys():
            if pfilter not in json_annuli.iterkeys():
                msg = ("%sError. Photometric parameters for the '%s' "
                       "filter not listed in '%s. Wrong file, maybe?")
                json_basename = os.path.basename(options.json_annuli)
                print msg % (style.prefix, pfilter, json_basename)
                print style.error_exit_message
                return 1

    print "%sSources image: %s" % (style.prefix, sources_img_path)
    print "%sRunning SExtractor on the sources image..." % style.prefix ,
    sys.stdout.flush()

    # Remove from the FITS header the path to the on-disk catalog, if present,
    # thus forcing SExtractor to detect sources on the image. This is necessary
    # because, if SExtractor (via the seeing.FITSeeingImage class) were run on
    # the image before it was calibrated astrometrically, the on-disk catalog
    # would only contain the X and Y image coordinates of the astronomical
    # objects, using zero for both their right ascensions and declinations.
    img = fitsimage.FITSImage(sources_img_path)
    img.delete_keyword(keywords.sex_catalog)

    args = (sources_img_path, options.maximum, options.margin)
    kwargs = dict(coaddk = options.coaddk)
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
        atexit.register(methods.clean_tmp_files, options.coordinates)
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

        sources_img_fwhm = get_fwhm(sources_img, options)
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

    # The options.exptimek FITS keyword is allowed to be missing from the
    # header of the sources image (for example, a legitimate scenario: we
    # detect sources on a mosaic created with IPAC's Montage, combining several
    # images). In those cases, qphot() uses the default value, an empty string.
    # We can ignore the MissingFITSKeyword warning for (and only for) the
    # sources image: it is not critical if magnitudes cannot be normalized to
    # an exposure time of one time unit, as these values are only expected to
    # serve as an estimate of how bright each astronomical object is.

    with warnings.catch_warnings():
        kwargs = dict(category = qphot.MissingFITSKeyword)
        warnings.filterwarnings('ignore', **kwargs)
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
    atexit.register(methods.clean_tmp_files, options.coordinates)

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

        with warnings.catch_warnings():
            kwargs = dict(category = qphot.MissingFITSKeyword)
            warnings.filterwarnings('ignore', **kwargs)
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
    #
    # In the case of the sources image, unlike for the images on which we do
    # photometry, there are several fields that are allowed to be None. This
    # is because we may detect sources on an image resulting from assembling
    # several ones into a custom mosaic: the resulting image does not have a
    # proper (a) photometric filter, (b) observation date, (c) airmass or (d)
    # gain. Therefore, we use None, which SQLite interprets as NULL.

    path = sources_img.path
    pfilter = methods.func_catchall(sources_img.pfilter, options.filterk)

    kwargs = dict(date_keyword = options.datek,
                  time_keyword = options.timek,
                  exp_keyword = options.exptimek)
    unix_time = methods.func_catchall(sources_img.date, **kwargs)

    # In theory, sources should be detected on the result on mosaicking several
    # FITS images, in order to improve the signal-to-noise ratio and allow for
    # a more accurate detection of faint astronomical objects. However, and as
    # Javier Blasco pointed out in issue #19, not all users need to do this: it
    # may be enough for them to use to detect sources one of the FITS images on
    # which they also want to do photometry.
    #
    # Allow to do photometry on the sources FITS image
    # [URL] https://github.com/vterron/lemon/issues/19
    #
    # In order to make this possible, ignore the Unix time and photometric
    # filter of the sources image (using None instead, regardless of what we
    # read from the FITS header) if there is an image with the same date and
    # filter among those on which we are going to do photometry. This prevents
    # the database.DuplicateImageError exception, with a message such as "Image
    # with Unix time 1325631812.2045 (Tue Jan 3 23:03:32 2012 UTC) and filter J
    # already in database"), from being raised. The idea is to store in the
    # output database as much information as possible about the sources image,
    # but if needed we can get by without these two values. After all, the data
    # about the sources image is mostly stored for book-keeping purposes, in
    # order to simplify future analysis and debugging.

    # Nested defaultdict, always returns a list
    if dates_counter[unix_time][pfilter]:

        # There can only be one FITS file with the same observation date and
        # photometric filter, as duplicate images were previously discarded.
        assert len(dates_counter[unix_time][pfilter]) == 1
        img = fitsimage.FITSImage(dates_counter[unix_time][pfilter][0])
        if pfilter == img.pfilter(options.filterk):

            msg1 = ("%s has the same date (%.4f, %s) and filter (%s) as the "
                    "sources image (%s)")
            date_str = methods.utctime(unix_time)
            args = (img.path, unix_time, date_str, pfilter, path)
            logging.debug(msg1 % args)

            msg2 = ("This must mean you are doing photometry on the FITS image "
                    "that you are also using to detect astronomical sources")
            logging.debug(msg2)

            msg3 = ("Avoid collision: ignore date and filter of the sources "
                    "image (store in the LEMONdB a None instead)")
            logging.debug(msg3)

            unix_time = None
            pfilter   = None

    object_ = methods.func_catchall(sources_img.read_keyword, options.objectk)
    airmass = methods.func_catchall(sources_img.read_keyword, options.airmassk)
    # If not given with --gaink, read it from the FITS header
    if options.gain:
        gain = options.gain
    else:
        gain = methods.func_catchall(sources_img.read_keyword, options.gaink)

    ra, dec = sources_img_ra, sources_img_dec

    args = (path, pfilter, unix_time, object_, airmass, gain, ra, dec)
    simage = database.Image(*args)
    output_db.simage = simage
    output_db.commit()

    for pfilter, images in sorted(files.iteritems()):
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

        if json_annuli:
            # Store all the CandidateAnnuli objects in the LEMONdB
            assert len(json_annuli[pfilter])
            for cand in json_annuli[pfilter]:
                output_db.add_candidate_pparams(cand, pfilter)

            filter_annuli = json_annuli[pfilter][0]
            aperture = filter_annuli.aperture
            annulus  = filter_annuli.annulus
            dannulus = filter_annuli.dannulus

            msg = "%sUsing the parameters listed in the JSON file, which are:"
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
                img_fwhm = get_fwhm(img, options)
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

            fwhm = get_fwhm(img, options)
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
        # each of the FITS images. This allows us to, in both cases, make the
        # map_async_args() generator loop over the images on which photometry
        # is to be done and, for each one of them, call qphot_params() to get
        # the parameters that have to be used.

        if not options.individual_fwhm:
            args = aperture, annulus, dannulus
            pparams = database.PhotometricParameters(*args)
            qphot_params = lambda x: pparams
        else:
            qphot_params = fwhm_derived_params

        def map_async_args():
            for path in images:
                img = fitsimage.FITSImage(path)
                yield (img, qphot_params(img), options)

        # Unlike the sources image, the options.exptimek FITS keyword is *not*
        # optional for the images on which we do photometry: qphot() needs it
        # to normalize the computed magnitudes to an exposure time of one time
        # unit. However, this point cannot be reached if one of the images does
        # not contain this keyword, as it was needed in order to make sure that
        # there are no duplicate observation dates. There is no need to turn
        # the MissingFITSKeyword warning into an exception.

        result = pool.map_async(parallel_photometry, map_async_args())
        methods.show_progress(0.0)
        while not result.ready():
            time.sleep(1)
            methods.show_progress(queue.qsize() / len(images) * 100)
            # Do not update the progress bar when debugging; instead, print it
            # on a new line each time. This prevents the next logging message,
            # if any, from being printed on the same line that the bar.
            if logging_level < logging.WARNING:
                print

        result.get() # reraise exceptions of the remote call, if any
        methods.show_progress(100) # in case the queue was ready too soon
        print

        msg = "%sStoring photometric measurements in the database..."
        print msg % style.prefix
        sys.stdout.flush()

        methods.show_progress(0)
        qphot_results = (queue.get() for x in xrange(queue.qsize()))
        for index, args in enumerate(qphot_results):
            db_image, pparams, img_qphot = args
            logging.debug("Storing image %s in database" % db_image.path)
            output_db.add_image(db_image)
            logging.debug("Image %s successfully stored" % db_image.path)

            # Now store each photometric measurement
            for object_id, object_phot in enumerate(img_qphot):
                # INDEF photometric measurements have a magnitude of None, and
                # those with at least one saturated pixel in the aperture have
                # a magnitude of infinity. In both cases the measurement is
                # useless for our photometric purposes and can be ignored.
                if object_phot.mag is None:
                    msg = "%s: object %d is INDEF (None)"
                    args = db_image.path, object_id
                    logging.debug(msg % args)
                    continue

                elif object_phot.mag == float('infinity'):
                    msg = "%s: object %d is saturated (infinity)"
                    args = db_image.path, object_id
                    logging.debug(msg % args)
                    continue

                else:
                    msg = "%s: object %d magnitude = %f"
                    args = db_image.path, object_id, object_phot.mag
                    logging.debug(msg % args)

                # Photometric measurements with a signal-to-noise ratio less
                # than or equal to one are ignored -- not only because these
                # measurements are anything but reliable, but also because such
                # values are outside of the domain of the function that
                # converts SNRs to errors in magnitudes.
                object_snr = object_phot.snr(db_image.gain)
                if object_snr <= 1:
                    msg = "%s: object %d ignored (SNR = %f <= 1)"
                    args = db_image.path, object_id, object_snr
                    logging.debug(msg % args)
                    continue

                else:
                    msg = "%s: object %d SNR = %f"
                    args = db_image.path, object_id, object_snr
                    logging.debug(msg % args)

                    msg = "%s: storing measurement for object %d in database"
                    args = db_image.path, object_id
                    logging.debug(msg % args)

                    args = (object_id,
                            db_image.unix_time,
                            db_image.pfilter,
                            object_phot.mag,
                            object_snr)

                    output_db.add_photometry(*args)

                    msg = "%s: measurement for object %d successfully stored"
                    args = db_image.path, object_id
                    logging.debug(msg % args)

            methods.show_progress(100 * (index + 1) / len(images))
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

