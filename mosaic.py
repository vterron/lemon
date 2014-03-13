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

description = \
""" This module receives a series of offsets, aligns their corresponding FITS
    images and finally combines them all into a single image by averaging them
    after rejecting a certain fraction of the lowest and highest pixels. The
    offsets are read from the XML file received as input, which is expected to
    have been outputted by the 'offset.py' module. Note that, therefore, all
    the offsets must refer to the same reference image, as only one reference
    FITS image is created.

    A very important note: as you surely already know, shifting images, as is
    needed in order to have them aligned, distorts slightly both signal and
    noise. This means that you should not (even worse, cannot!) do photometry
    on the resulting image. It should be only used in order to detect sources,
    as it maximizes the signal-to-noise ratio. Please stand on the shoulders of
    the photometrist giants and adhere to the standards in our guild.

    Note that it is assumed that, for offsets, North is up and East is left.
    Images for which North is down or East is right are not yet supported, and
    using them shall have catastrophic results. Yes, this is the unexpected
    consequence of having most of your astronomers and astronomical software
    developers based in countries in the Northern hemisphere.

"""

import atexit
import logging
import math
import montage_wrapper as montage
import numpy
import os
import os.path
import optparse
import pyfits
import shutil
import sys
import tempfile
import traceback

# LEMON modules
import customparser
import keywords
import fitsimage
import methods
import style
import xmlparse

def clean_tmp_dir(dir_path):
    """ Try to delete an entire directory tree. """

    msg = "Cleaning up temporary directory '%s'"
    logging.debug(msg % dir_path)

    error_count = []

    def log_error(function, path, excinfo):
        """ Error handler for shutil.tree() """

        # nonlocal is not available in Python 2.x so, being it outside of the
        # local scope, we cannot rebind 'error_count' each time we come across
        # an error. Instead of initializing it to zero and incrementing it by
        # one every time this function is called, use it as a list, appending
        # an element for each error.
        error_count.append(1)
        msg = "%s: error deleting '%s' (%s)"
        args = function, path, excinfo[1]
        logging.debug(msg % args)

    try:
        kwargs = dict(ignore_errors = False, onerror = log_error)
        shutil.rmtree(dir_path,**kwargs)

    finally:
        msg = "Temporary directory '%s' deleted"
        if error_count:
            msg += " (but there were failed removals)"
        logging.debug(msg % dir_path)


parser = customparser.get_parser(description)
parser.usage = "%prog [OPTION]... OFFSETS_XML_FILE"

parser.add_option('--overwrite', action = 'store_true', dest = 'overwrite',
                  help = "overwrite output image if it already exists")

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
    (options, args) = parser.parse_args(args = arguments)

    # Print the help and abort the execution if there are fewer than three
    # positional arguments left, as the user must specify at least two FITS
    # images and the output mosaic into which they are assembled.
    if len(args) < 3:
        parser.print_help()
        return 2 # used for command line syntax errors
    else:
        assert len(args) >= 3
        input_paths = set(args[1:-1])
        output_path = args[-1]

    # Refuse to overwrite the output FITS file unless explicitly instructed to
    # do so. Note that, if the --overwritten option is given, we do not need to
    # delete the existing file: it will be silently overwritten when the output
    # of montage.mosaic() is shutil.move()'d to the output path.

    if os.path.exists(output_path):
        if not options.overwrite:
            msg = "%sError. The output file '%s' already exists."
            print msg % (style.prefix, output_path)
            print style.error_exit_message
            return 1

    # montage.mosaic() requires as first argument the directory containing the
    # input FITS images but, in order to maintain the same syntax across all
    # LEMON commands, we receive them as command-line arguments. Thus, create a
    # temporary directory and symlink from it the input images. Hard links are
    # not an option because os.link() will raise "OSError: [Errno 18] Invalid
    # cross-device link" if the temporary directory is created in a different
    # partition.

    pid = os.getpid()
    suffix = "_LEMON_%d_mosaic" % pid
    kwargs = dict(suffix = suffix + '_input')
    input_dir = tempfile.mkdtemp(**kwargs)
    atexit.register(clean_tmp_dir, input_dir)

    for path in input_paths:
        source = os.path.abspath(path)
        basename = os.path.basename(path)
        link_name = os.path.join(input_dir, basename)
        os.symlink(source, link_name)

    # The output of montage.mosaic() is another directory, to which several
    # files are written, so we need the path to a second temporary directory.
    # Delete it before calling mosaic(), as otherwise it will raise IOError
    # ("Output directory already exists").

    kwargs = dict(suffix = suffix + '_output')
    output_dir = tempfile.mkdtemp(**kwargs)
    atexit.register(clean_tmp_dir, output_dir)
    os.rmdir(output_dir)

    montage.mosaic(input_dir, output_dir)

    # montage.mosaic() writes several files to the output directory, but we are
    # only interested in one of them: 'mosaic.fits', the mosaic FITS image. We
    # need to move it to the output path specified by the user.

    MOSAIC_OUTPUT = 'mosaic.fits'
    src = os.path.join(output_dir, MOSAIC_OUTPUT)
    shutil.move(src, output_path)

    print "%sYou're done ^_^" % style.prefix
    return 0

if __name__ == "__main__":
    sys.exit(main())

