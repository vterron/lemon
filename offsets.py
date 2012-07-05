#! /usr/bin/env python

# Copyright (c) 2012 Victor Terron. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of LEMON.
#
# LEMON is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
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
This module receives a series of FITS images and calculates the translation
offsets between the first one, which is considered to be the reference image,
and the rest of them. The offsets are calculated by means of cross-correlating
the 1-D profiles of the images in both axes (horizontal and vertical) and
finding the exact point at which most stars overlap. You may think of this as
if we were 'sliding' one of the images over the other and finding the point at
which most stars were aligned. The offsets are saved to a XML file which can be
parsed by LEMON at the next stages.

"""

import multiprocessing
import os
import os.path
import sys
import time

# LEMON modules
import defaults
import keywords
import fitsimage
import optparse
import methods
import seeing
import style
import xmlparse

def offset(reference_path, shifted_path, maximum, margin, per,
           filter_keyword, date_keyword, exp_keyword, coadd_keyword):
    """ Calculate the offset between two FITS images.

    The method returns an instance of xmlparse.XMLOffset, which encapsulates
    the translation offset between the reference image and the shifted image.
    This is just a high-level wrapper around FITSeeingImage.offset method, as
    XMLOffset instances, which give the offset of the 'shifted' image with
    respect to 'reference' image, also contain information on the photometric
    filter and date of observation of the former.

    FITSeeingImage.offset estimates the offset by generating a matrix
    representation (or 'single-point' mask) of both images and
    cross-correlating the 1-D profiles of the images in both axes (horizontal
    and vertical). The correct translation offset is that which maximizes the
    number of images that overlap. Although the matrix representation of the
    images rounds each star center to the nearest integer, sub-pixel precision
    is still achievable by fitting a parabola to the peak of each 1-D profile
    and finding its maximum. See FITSeeingImage.offset for further information.

    ValueError is raised in case no star makes it to the array representation
    of one of the two FITS images involved. This usually happens because
    'margin' is set to a too large a value, although it may also occur if its
    value is reasonable but there are few stars, located too close to the image
    borders. The same exception is also raised if the two images do not have
    the same dimensions. KeyError is raised if one of the keywords cannot be
    found in the FITS header.

    Arguments:
    reference_path - the path to the reference FITS image.
    shifted_path - the path to the shifted ('moved') FITS image.
    maximum - level at which arises saturation, in ADUs. For coadded
              observations, the effective saturation level is obtained by
              multiplying this value by the number of coadded images.
    margin - the width, in pixels, of the areas adjacent to the borders
             that will not be considered for the matrix representation of
             the images. Stars whose whose center is fewer than 'margin'
             pixels from any border (either horizontal or vertical) of the
             FITS image are not included in the matrix representation.
    per - the score at the given percentile of the signal-to-noise ratio of the
          stars that are to be included in the matrix representation of the
          images. In other words: stars whose SNR is below the 'per' percentile
          are excluded from the matrix representation. Note that the percentile
          is calculated taking into account only the stars within the image
          margins (see 'margin' argument).
    filter_keyword - keyword that stores the photometric filter in which the
                     FITS images was taken.
    date_keyword - FITS keyword that stores the date of the observation, in the
                   format specified in the FITS Standard. The old date format
                   was 'yy/mm/dd' and may be used only for dates from 1900
                   through 1999.  The new Y2K compliant date format is
                   'yyyy-mm-dd' or 'yyyy-mm-ddTHH:MM:SS[.sss]'.
    exp_keyword - the FITS keyword in which the duration of the exposure is
                  stored. It is expected to be a floating-point number which
                  gives the duration in seconds. The exact definition of
                  'exposure time' is mission dependent and may, for example,
                  include corrections for shutter open and close duration,
                  detector dead time, vignetting, or other effects.
    coaddk_keyword - FITS keyword that stores the number of effective coadds.
                     If the keyword is missing, we assume a value of one (that
                     is, that the observation consisted in a single exposure).

    """

    reference = seeing.FITSeeingImage(reference_path, maximum,
                                      margin, coaddk = coadd_keyword)
    shifted = seeing.FITSeeingImage(shifted_path, maximum,
                                    margin, coaddk = coadd_keyword)

    # Get the filter in which the shifted image was taken and the date of the
    # observation. Raise KeyError if it could not be found in the FITS header.
    shifted_filter = shifted.read_keyword(filter_keyword)
    shifted_date   = shifted.date(date_keyword = date_keyword,
                                  exp_keyword = exp_keyword)

    x_offset, y_offset, x_overlap, y_overlap = \
        reference.offset(shifted, per = per)

    return xmlparse.XMLOffset(reference.path, shifted.path,
                              shifted_filter, shifted_date,
                              x_offset, y_offset, x_overlap, y_overlap)

# The Queue is global -- this works, but note that we could have
# passed its reference to the function managed by pool.map_async.
# See http://stackoverflow.com/a/3217427/184363
queue = multiprocessing.Queue()

def parallel_offset(args):
    """ Method argument of map_async to compute offsets in parallel.

    Functions defined in classes don't pickle, so we have moved this code here
    in order to be able to use it with multiprocessing's map_async. As it
    receives a single argument, values are passed in a tuple which is then
    unpacked.

    """

    reference, shifted_path, options = args

    # We did not check that the given paths do actually refer to existing, FITS
    # files before starting the work in parallel. Thus, IOError may be raised
    # if 'shifted_path' does not exist, as well as the self-explanatory
    # NonFITSFile and NonStandardFITS exceptions. In these cases, the path is
    # silently ignored and the worker (CPU core) can move to the next path.

    try:
        fitsimage.FITSImage(shifted_path)
    except (IOError, fitsimage.NonFITSFile, fitsimage.NonStandardFITS):
        return

    img_offset = offset(reference.path, shifted_path, options.maximum,
                        options.margin, options.percentile, options.filterk,
                        options.datek, options.exptimek, options.coaddk)
    queue.put(img_offset)


parser = optparse.OptionParser(description = description,
                               formatter = style.NewlinesFormatter())

parser.usage = "%prog [OPTION]... REFERENCE_IMAGE  REST_OF_IMAGES..."
parser.add_option('-o', action = 'store', type = 'str',
                  dest = 'output_xml', default = 'offsets.xml',
                  help = "path to the output XML file [default: %default]")

parser.add_option('-w', action = 'store_true', dest = 'overwrite',
                  help = "overwrite output XML file if it already exists")

parser.add_option('--cores', action = 'store', type = 'int',
                  dest = 'ncores', default = multiprocessing.cpu_count(),
                  help = "the maximum number of cores available to the "
                  "module. This option defaults to the number of CPUs in "
                  "the system, which are automatically detected "
                  "[default: %default]")

mask_group = optparse.OptionGroup(parser, "Single-point-masks",
             "The offsets between the FITS images are not computed using them "
             "directly, but instead with a two-dimensional NumPy array in "
             "which each object detected by SExtractor is represented by a "
             "single, non-zero value, located at the coordinates specified by "
             "its X_IMAGE and Y_IMAGE parameters in the SExtractor catalog. "
             "Unfortunately, the center of each star is practically never an "
             "exact value. Instead, centers will most likely be something "
             "like (311.046, 887.279), so they have to be rounded to the "
             "nearest integer.")

# Note for developers: this value is not per se needed to compute the offsets,
# but it is required by the instantiation method of the FITSeeingImage class.
# Although we could pass any other value, using the same saturation level as in
# other stages of the pipeline allows us to reuse the on-disk cached SExtractor
# catalog.
parser.add_option('-m', action = 'store', type = 'int',
                  dest = 'maximum', default = defaults.maximum,
                  help = defaults.desc['maximum'])

parser.add_option('--margin', action = 'store', type = 'int',
                  dest = 'margin', default = defaults.margin,
                  help = defaults.desc['margin'])

mask_group.add_option('--per', action = 'store', type = 'int',
                      dest = 'percentile', default = 50,
                      help = "the score at the given percentile of the "
                      "signal-to-noise ratio of the stars that are to be "
                      "included in the matrix representation of the images. "
                      "In other words: stars whose SNR is below this "
                      "percentile are excluded from the matrix "
                      "representation. Note that the percentile is "
                      "calculated taking into account only the stars within "
                      "the image margins (see the --margin option) "
                      "[default: %default]")
parser.add_option_group(mask_group)

key_group = optparse.OptionGroup(parser, "FITS Keywords",
                                 keywords.group_description)

key_group.add_option('--filterk', action = 'store', type = 'str',
                     dest = 'filterk', default = keywords.filterk,
                     help = keywords.desc['filterk'])

key_group.add_option('--datek', action = 'store', type = 'str',
                     dest = 'datek', default = keywords.datek,
                     help = keywords.desc['datek'])

key_group.add_option('--expk', action = 'store', type = 'str',
                     dest = 'exptimek', default = keywords.exptimek,
                     help = keywords.desc['exptimek'])

key_group.add_option('--coaddk', action = 'store', type = 'str',
                     dest = 'coaddk', default = keywords.coaddk,
                     help = keywords.desc['coaddk'])

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
    (options, args) = parser.parse_args(args = arguments)

    if not args:
        parser.print_help()
        return 2 # used for command line syntax errors
    elif len(args) < 2:
        print "%sError. At least one image other than REFERENCE_IMAGE must " \
              "be given." % style.prefix
        print style.error_exit_message
        return 2
    else:
        reference_path       = args[0]
        rest_of_images_paths = args[1:]

    # Do not overwrite an existing XML offsets file unless the user actually
    # (-w option) intended to do so. We could not bear the loss of hours of
    # careful computation because of a negligent typo by the astronomer.
    if os.path.exists(options.output_xml):
        if not options.overwrite:
            print "%sError. The output XML file '%s' already exists." % \
                  (style.prefix, options.output_xml)
            print style.error_exit_message
            return 1
        else:
            os.unlink(options.output_xml)

    print "%sIndexing FITS files in REST_OF_IMAGES..." % style.prefix ,
    indexed_paths = fitsimage.find_files(rest_of_images_paths)
    print 'done.'

    if not indexed_paths:
        raise TypeError("REST_OF_IMAGES is empty")

    print "%sDetecting sources on %s..." % (style.prefix, reference_path) ,
    sys.stdout.flush()
    reference_img = seeing.FITSeeingImage(reference_path, options.maximum,
                                          options.margin, options.coaddk)
    print 'done.'

    # The original code used the FITSeeingSet.fitsimage_cast() to convert at
    # once all the images to FITSeeingImage instances, thus running SExtractor
    # consecutively on all the images that needed it. Although this methodical
    # approach simplified the code, it has the unsurmountable problem (given
    # the current, barbaric technology of our times, readers of the future)
    # that causes the module to eat up *lots* of memory: hundreds of
    # FITSeeingImage instances, each one of them with thousands of Star
    # objects, result in the module unnecessary needing gigabytes of RAM.
    #
    # What we do now, instead, is to cast the images to FITSeeingImage one by
    # one, deleting each instance as we calculate the offsets as they are no
    # longer needed. An easy to implement improvement, yes, but which lowers
    # the memory requirements to what you could expect: almost nothing.

    offsets_list = []
    print "%sCalculating offsets between REFERENCE_IMAGE and the %d " \
          "indexed FITS images... " % (style.prefix, len(indexed_paths))

    # Use a pool of workers to which to assign the calculation of the offsets
    pool = multiprocessing.Pool(options.ncores)

    # map_async supports only one iterable argument, so we need to pass all the
    # arguments as a tuple which is then unpacked by the parallel_offset method
    map_async_args = ((reference_img, path, options) for path in indexed_paths)
    result = pool.map_async(parallel_offset, map_async_args)

    methods.show_progress(0.0)
    while not result.ready():
        time.sleep(1)
        methods.show_progress(queue.qsize() / len(indexed_paths) * 100)

    result.get() # reraise exceptions of the remote call, if any
    methods.show_progress(100) # in case the queue was ready too soon
    print

    # The offsets are sorted by their observation date
    offsets_list = sorted(queue.get() for x in xrange(queue.qsize()))

    print  "%sGenerating the output XML tree..." % style.prefix ,
    sys.stdout.flush()
    xml_tree = xmlparse.XMLOffsetFile()

    for xml_offset in offsets_list:
        xml_tree.add(xml_offset)
    print 'done.'

    print "%sSaving the translation offsets to '%s'..." % \
          (style.prefix, options.output_xml) ,
    sys.stdout.flush()
    xml_tree.dump(options.output_xml)
    print 'done.'

    print "%sYou're done ^_^" % style.prefix
    return 0

if __name__ == "__main__":
    sys.exit(main())

