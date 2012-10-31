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
import logging
import multiprocessing
import numpy
import optparse
import os
import os.path
import sys
import time
import warnings

# LEMON modules
import astromatic
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
    img_unix_time = photometry_image.date(options.datek, options.exptimek)
    img_object = img_offset.object
    img_airmass = photometry_image.read_keyword(options.airmassk)
    img_gain = options.gain or photometry_image.read_keyword(options.gaink)
    img_xoffset  = img_offset.x
    img_xoverlap = img_offset.x_overlap
    img_yoffset = img_offset.y
    img_yoverlap = img_offset.y_overlap

    logging.debug("%s: filter = %s" % (img_path, img_pfilter))
    logging.debug("%s: observation date: %.2f (%s)" % \
                  (img_path, img_unix_time, time.ctime(img_unix_time)))
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


parser = optparse.OptionParser(description = description,
                               formatter = style.NewlinesFormatter())

parser.usage = "%prog [OPTION]... OFFSETS_XML_FILE"
parser.add_option('-o', action = 'store', type = 'str',
                  dest = 'output_db', default = 'photometry.LEMONdB',
                  help = "path to output database [default: %default]")

parser.add_option('-w', action = 'store_true', dest = 'overwrite',
                  help = "overwrite output database if it already exists")

parser.add_option('--passband', action = 'store', type = str,
                  dest = 'passband', default = None,
                  help = "do not do photometry on all the images whose "
                  "offset is listed in OFFSETS_XML_FILE, but instead only on "
                  "those whose passband (i.e., the filter with which they were "
                  "observed) matches the value specified by this option")

parser.add_option('-m', action = 'store', type = 'int',
                  dest = 'maximum', default = defaults.maximum,
                  help = defaults.desc['maximum'])

parser.add_option('--margin', action = 'store', type = 'int',
                  dest = 'margin', default = defaults.margin,
                  help = defaults.desc['margin'])

parser.add_option('-g', action = 'store', type = 'float',
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

qphot_group.add_option('-p', action = 'store', type = 'float',
                       dest = 'aperture', default = 3.0,
                       help = "the aperture radius, in number of times the "
                       "median FWHM [default: %default]")

qphot_group.add_option('-i', action = 'store', type = 'float',
                       dest = 'annulus', default = 4.5,
                       help = "the inner radius of the sky annulus, in "
                       "number of times the median FWHM [default: %default]")

qphot_group.add_option('-k', action = 'store', type = 'float',
                       dest = 'dannulus', default = 1.0,
                       help = "the width of the sky annulus, in number "
                       "of times the median FWHM [default: %default]")

qphot_group.add_option('--min', action = 'store', type = 'float',
                       dest = 'min', default = 3.0,
                       help ="the minimum width of the sky annulus, in "
                       "pixels, regardless of the value specified in the "
                       "above parameter. This option is intended to prevent "
                       "small FWHMs from resulting in too thin an sky "
                       "annulus. [default = %default]")
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
              "--passband, in order to to photometry one filter at a time.")

qphot_fixed.add_option('--pp', action = 'store', type = 'float',
                       dest = 'aperture_pix', default = None,
                       help = "the aperture radius, in pixels")

qphot_fixed.add_option('--ip', action = 'store', type = 'float',
                       dest = 'annulus_pix', default = None,
                       help = "the inner radius of the sky annulus, in pixels")

qphot_fixed.add_option('--kp', action = 'store', type = 'float',
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

key_group.add_option('--rak', action = 'store', type = 'str',
                     dest = 'rak', default = keywords.rak,
                     help = keywords.desc['rak'])

key_group.add_option('--deck', action = 'store', type = 'str',
                     dest = 'deck', default = keywords.deck,
                     help = keywords.desc['deck'])

key_group.add_option('--datek', action = 'store', type = 'str',
                     dest = 'datek', default = keywords.datek,
                     help = keywords.desc['datek'])

key_group.add_option('--expk', action = 'store', type = 'str',
                     dest = 'exptimek', default = keywords.exptimek,
                     help = keywords.desc['exptimek'])

key_group.add_option('--airmk', action = 'store', type = 'str',
                     dest = 'airmassk', default = keywords.airmassk,
                     help = keywords.desc['airmassk'])

key_group.add_option('--coaddk', action = 'store', type = 'str',
                     dest = 'coaddk', default = keywords.coaddk,
                     help = keywords.desc['coaddk'])

key_group.add_option('--gaink', action = 'store', type = 'str',
                     dest = 'gaink', default = keywords.gaink,
                     help = keywords.desc['gaink'])

key_group.add_option('--uik', action = 'store', type = 'str',
                     dest = 'uncimgk', default = keywords.uncimgk,
                     help = keywords.desc['uncimgk'])

key_group.add_option('--fwhmk', action = 'store', type = 'str',
                     dest = 'fwhmk', default = keywords.fwhmk,
                     help = "keyword for the FWHM of the image, stored by "
                     "LEMON at the seeing.py stage [default: %default]")
parser.add_option_group(key_group)

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

    # All the parameters must be casted to string before beging passed to the
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

    if len(args) != 1:
        parser.print_help()
        return 2 # used for command line syntax errors
    else:
        assert len(args) == 1
        xml_path      = args[0]

    # If the user gives an empty string as the FITS keyword which stores the
    # path to the original image, it is understood as meaning that we want
    # saturation to be checked for in the same images on which photometry is
    # done. If that is the case, we need to set the option to None (although an
    # empty string would also work), as that is what qphot.run expects to
    # receive in these cases.

    if not options.uncimgk:
        options.uncimgk = None

    # It goes without saying that the input offsets XML file must exist...
    if not os.path.exists(xml_path):
        print "%sError. The file '%s' does not exist." % (style.prefix, xml_path)
        print style.error_exit_message
        return 1

    # ... and the same goes for the annuli XMl file (that generated by the
    # annuli.py module, whose path may have been passed with the --annuli
    # option). The use of this file automatically discards whathever was
    # specified with the Aperture Photometry (FWHM and pixels) options.
    xml_annuli   = None
    fixed_annuli = False

    if options.xml_annuli:
        if not os.path.exists(options.xml_annuli):
            print "%sError. The file '%s' does not exist." % \
                  (style.prefix, xml_path)
            print style.error_exit_message
            return 1
        else:
            xml_annuli = \
                xmlparse.CandidateAnnuli.xml_load(options.xml_annuli, best_only = True)
            print "%sPhotometric paramaters read from the '%s' file." % \
                  (style.prefix, os.path.basename(options.xml_annuli))

    # Even if the annuli XML is used by the module, the aperture and sky annuli
    # parameters (whether FWHM-based or given in pixels) are still needed in
    # order to do photometry on the reference image and extract the
    # instrumental magnitude that for each star is stored in the database.

    # This warning is printed if one or more, but not the three, options are given
    if 0 < bool(options.aperture_pix) + bool(options.annulus_pix) + \
           bool(options.dannulus_pix) < 3:

        msg = ("%(p)sWarning: the --pp, --ip and --kp options must always be "
               "used in conjunction;\n%(p)sthe sizes of the apertures and sky "
               "annuli will still be given by the FWHM.") % \
               {'p': style.prefix}
        warnings.warn(msg)

    if options.aperture_pix and options.annulus_pix and options.dannulus_pix:
        fixed_annuli = True

    # The aperture, annulus and dannulus values, whether expressed in number of
    # times the median FWHM or by a fixed number of pixels, must be positive
    # numbers. By definition, also, the inner radius of the sky annulus must be
    # greater than or equal to the aperture radius. Obviously!

    fwhm_options  = (options.aperture, options.annulus, options.dannulus)
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

    if os.path.exists(options.output_db):
        if not options.overwrite:
            print "%sError. The output database '%s' already exists." % \
                  (style.prefix, options.output_db)
            print style.error_exit_message
            return 1
        else:
            os.unlink(options.output_db)

    print "%sValidating input XML file '%s'..." % \
          (style.prefix, os.path.basename(xml_path)) ,
    xmlparse.validate_dtd(xml_path)
    print 'done.'

    # Load all the offsets, saving them as a list of XMLOffset instances.
    # Note that the XML file is initially read as a XMLOffsetsFile instance,
    # from which we extract all the XMLOffsets and save them in a list.
    print "%sLoading input XML file into memory..." % style.prefix,
    sys.stdout.flush()
    xml_offsets = list(xmlparse.XMLOffsetFile(xml_path))

    # Also, make sure that all of them pertain to the same reference frame.
    # Otherwise, whiningly complain and abort the execution. Note the cast to
    # map and then to list, in order to remove duplicated names.
    list_of_reference_imgs = list(set((x.reference for x in xml_offsets)))
    if len(list_of_reference_imgs) != 1:
        print style.prefix
        print "%sError. The offsets must all refer to the same reference " \
              "image." % style.prefix
        print "%sHowever, multiple reference images were found in the " \
              "input XML file:" % style.prefix
        print "%s[%s]" % (style.prefix, ", ".join(list_of_reference_imgs))
        print style.error_exit_message
        return 1
    else:
        print 'done.'
        print "%s%d offsets were read." % (style.prefix, len(xml_offsets))
        reference_img_path = list_of_reference_imgs[0]


    # Although all the offsets listed in the XML file are loaded into memory,
    # the --passband option allows the user to specify which must be taken into
    # account: only those offsets for which the shifted image was taken with
    # the specified filter. The comparison is made by looking at the wavelength
    # of the filter -- after all, 'Johnson V' and 'V' are the same filter, but
    # if we just compared both strings they would not be so. The rest will be
    # discarded and lost like tears in rain.

    if options.passband:

        print "%sIgnoring offsets for shifted images with a passband other " \
              "than '%s'..." % (style.prefix, options.passband) ,
        sys.stdout.flush()

        pfilter = passband.Passband(options.passband)
        xml_offsets = [x for x in xml_offsets if x.filter.wavelength == pfilter.wavelength]

        if not xml_offsets:
            print style.prefix
            print "%sError. No shifted image was taken with the '%s' " \
                  "passband." % (style.prefix, options.passband)
            print style.error_exit_message
            return 1

        else:
            print 'done.'
            print "%s%d offsets matched the '%s' passband; the rest were " \
                  "discarded." % (style.prefix, len(xml_offsets),
                                  options.passband)


    # If an annuli XML file is being used, it must list the photometric
    # parameters for all the filters on which photometry is to be done, whether
    # they are all the input images or only those taken in the filter specified
    # by the --passband option.

    if xml_annuli:
        xml_offsets_bands = set(x.filter for x in xml_offsets)
        for pfilter, annuli in xml_annuli.iteritems():
            # Only the best CandidateAnnuli must have been loaded
            assert len(annuli) == 1
            if pfilter not in xml_offsets_bands:
                print "%sError. The photometric parameters for the %s filter " \
                      "are not listed in" % (style.prefix, pfilter)
                print "%sthe '%s' file. Wrong file, maybe?" % (style.prefix,
                       os.path.basename(options.xml_annuli))
                print style.error_exit_message
                return 1

    # But, hey, shouldn't we make sure that the reference image about which we
    # are talking that much does really _exist_ and is a FITS, standard-
    # conforming file? Just to make sure, you know. In order to do that, simply
    # pass it to FITSImage, which will throw the appropiate exceptions if an
    # error is encountered.
    reference_img = fitsimage.FITSImage(reference_img_path)

    # Light curves, which are our ultimate goal, can only have one differential
    # magnitude for each point in time. Therefore, we cannot have two or more
    # images with the same observation date. When this happens, all the images
    # with the duplicate date must be discarded. This may seem unlikely to
    # happen, but astronomical instruments also have software errors, and we
    # have already encountered this while reducing images taken with Omega2000,
    # a wide-field camera for the 3.5m telescope at Calar Alto.
    #
    # We may be tempted to keep one of them, such as, for example, that which
    # the highest number of sources, or the one with the best full width at
    # half maximum. However, as Matilde Fernandez thoughtfully pointed out, the
    # safest bet is to discard them all: how do we know which of the images has
    # the right observation date, and which have it wrong? We have no way of
    # finding it out; the only certain thing is that an error occurred. Thus,
    # we better forget about all these images.

    msg = "%sMaking sure there are no images with the same observation date..."
    print msg % style.prefix ,
    sys.stdout.flush()

    # Use a defaultdict, defaulting to an empty list, to map each Unix time to
    # the list of XMLOffsets (and the corresponding FITS images) that have it
    # as the observation date. If there are no duplicate dates, all the values
    # of the dictionary will have a length of one.
    dates_counter = collections.defaultdict(list)
    for offset in xml_offsets:
        dates_counter[offset.date].append(offset)

    if max([len(offsets) for offsets in dates_counter.itervalues()]) == 1:
        print 'done.'

    else:
        # Find the Unix times for which there is more than one XMLOffset
        # (and therefore a FITS image) and remove them all from the list.
        discarded = 0
        for unix_time, utime_offsets in dates_counter.iteritems():
            ntimes = len(utime_offsets)
            if ntimes > 1:

                # Newline needed only the first time, since the "Making sure
                # there are no images..." printed above did not include it.
                if not discarded:
                    print

                msg = "%sWarning! Multiple images have date %s"
                warnings.warn(msg % (style.prefix, time.ctime(unix_time)))
                for offset in utime_offsets:
                    xml_offsets.remove(offset)
                    discarded += 1
                    msg = "%sImage %s excluded"
                    warnings.warn(msg % (style.prefix, offset.shifted))

        # There should be no duplicate observation dates anymore, and at least
        # two offsets should have been discarded -- otherwise, how did we get
        # into the 'else' clause in the first place?
        assert len(set(x.date for x in xml_offsets)) == len(xml_offsets)
        assert discarded >= 2

        msg = "%s%d images had duplicate dates and were discarded, %d remain"
        print msg % (style.prefix, discarded, len(xml_offsets))

    # Just a reminder of which the reference image (i.e., the image on which
    # sources are detected and to which the offsets in the XML file refer) and
    # the output database are.
    print "%sReference image: %s" % (style.prefix, reference_img.path)
    print "%sOutput database: %s" % (style.prefix, options.output_db)

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
    # image and itself is... zero. Who would have guessed that? Do not bother
    # reading the object name from the FITS header as it is irrelevant here.
    reference_date = reference_img.date(date_keyword = options.datek,
                                        exp_keyword = options.exptimek)
    null_offset = xmlparse.XMLOffset(reference_img.path, reference_img.path, 'N/A',
                                     'N/A', reference_date, 0.0, 0.0, 'inf', 'inf')

    if not fixed_annuli:

        try:
            reference_fwhm = reference_img.read_keyword(options.fwhmk)
        except KeyError:
            # The FWHM of the image was not found at the FITS header. This may
            # happen, even if all the images went through the seeing.py module,
            # if the FWHM keyword has whitespaces (such as "LEMON FWHM"), since
            # SWarp, or at least its version 2.19.1, seems not to propagate
            # them. Better than aborting the execution, and until this bug (?)
            # is fixed, we compute it now, using the default parameters of
            # seeing.py.

            seeing_default_opts = seeing.parser.defaults
            fwhm_per  = seeing_default_opts['per']
            fwhm_mode = seeing_default_opts['mean'] and 'mean' or 'median'
            reference_fwhm = reference_img.fwhm(per  = fwhm_per,
                                                mode = fwhm_mode)

        assert isinstance(reference_fwhm, float)
        reference_aperture = options.aperture * reference_fwhm
        reference_annulus  = options.annulus  * reference_fwhm
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
        reference_annulus  = options.annulus_pix
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
            x    = star_photometry.x
            y    = star_photometry.y

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

    output_db = database.LEMONdB(options.output_db)

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

    rimage_filter = reference_img.read_keyword(options.filterk)
    pfilter = passband.Passband(rimage_filter)
    unix_time = reference_img.date(options.datek, options.exptimek)
    object_ = reference_img.read_keyword(options.objectk)
    airmass  = reference_img.read_keyword(options.airmassk)
    gain = options.gain or reference_img.read_keyword(options.gaink)
    args = reference_img.path, pfilter, unix_time, object_, airmass, gain
    output_db.rimage = database.ReferenceImage(*args)

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
        # be extracted from the --annuli file is simple: just take the
        # CandidateAnnuli instance with the best photometric paramaters for
        # this filter (as we previously made sure that there is an instance for
        # each passband on which photometry is to be done) and use them.
        #
        # Alternatively, if the dimensions of the annuli are to be determined
        # by the median FWHM of the images, this has to be done for each
        # different filter in which images were taken. This contrasts with when
        # specific sizes (in pixels) are given for the annuli, which are used
        # for all the filters.

        if xml_annuli:

            assert len(xml_annuli[pfilter]) == 1
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

        elif not fixed_annuli:
            print "%sCalculating the median FWHM for this filter..." %  \
                  style.prefix ,
            sys.stdout.flush()

            img_fwhms = []
            for img in band_offsets:
                logging.debug("%s: reading value of '%s' keyword..." % \
                              (img.shifted, options.fwhmk))
                image_fwhm = fitsimage.FITSImage(img.shifted).read_keyword(options.fwhmk)
                logging.debug("%s: FWHM = %.3f (read from FITS header)" % \
                              (img.shifted, image_fwhm))
                img_fwhms.append(image_fwhm)


            assert all([isinstance(x, float) for x in img_fwhms])
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

        # Stars are stored in 'reference_stars' as database.DBStar instances.
        # However, the photometry function receives the pixels as a sequence
        # of astromatic.Pixel instances; hence the cast. Note that in
        # 'reference_stars' we have a five-element tuple, where the x and y
        # coordinates are the first two values.
        reference_pixels = [astromatic.Pixel(*x[:2]) for x in reference_stars]
        qphot_params = float(aperture), float(annulus), float(dannulus)

        # The task of doing photometry on a series of images is inherently
        # parallelizable; use a pool of workers to which assign the images.
        pool = multiprocessing.Pool(options.ncores)
        map_async_args = \
            ((offset, reference_pixels, qphot_params, options)
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
                magnitude  = star_photometry.magnitude
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
                    output_db.add_photometry(star_id, pimage.unix_time, magnitude, snr)
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

    methods.owner_writable(options.output_db, False) # chmod u-w
    print "%sYou're done ^_^" % style.prefix
    return 0

if __name__ == "__main__":
    sys.exit(main())

