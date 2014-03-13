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
Use the Montage (Montage Astronomical Image Mosaic Engine) toolkit [1] to
assemble the input FITS images into a composite mosaic that preserves their
flux calibration and positional fidelity. This is a high-level interface to
mosaic(), a convenience function of the montage-wrapper module which runs the
mosaicking process from start to end. The input FITS images, all of which must
have been astrometrically calibrated, are reprojected onto a common coordinate
system and combined into a mosaic.

Montage is an extremely powerful toolkit, whose algorithms preserve the
astrometric and photometric accuracy of the input images and perform background
rectification in such a fashion that its impact on the photometric quality of
the data is almost negligible [2]. For example, according to the results of an
accuracy testing, 99.7% of the sources in re-projected synthetic images were
within 0.1% of the original flux [3]. There are, however, certain assumptions
of which you should be aware. For example, Montage assumes that the input
images are all calibrated to an absolute energy scale and that any
discrepancies between the images are due to variations in their background
levels that are terrestrial or instrumental in origin [4]

Note that montage_wrapper is not a replacement for the IPAC Montage mosaicking
software, whose commands (such as mAdd or mProject) must be present in PATH.

[1]_http://montage.ipac.caltech.edu/
[2]_http://adsabs.harvard.edu/abs/2003ASPC..295..343B
[3]_http://montage.ipac.caltech.edu/docs/accuracy.html
[4]_http://montage.ipac.caltech.edu/docs/algorithms.html

"""

import atexit
import logging
import montage_wrapper as montage
import os
import os.path
import shutil
import sys
import tempfile

# LEMON modules
import customparser
import fitsimage
import methods
import style

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
        shutil.rmtree(dir_path, **kwargs)

    finally:
        msg = "Temporary directory '%s' deleted"
        if error_count:
            msg += " (but there were failed removals)"
        logging.debug(msg % dir_path)


parser = customparser.get_parser(description)
parser.usage = "%prog [OPTION]... INPUT_IMGS... OUTPUT_IMG"

parser.add_option('--overwrite', action = 'store_true', dest = 'overwrite',
                  help = "overwrite output image if it already exists")

parser.add_option('--background-match', action = 'store_true',
                  dest = 'background_match',
                  help = "include a background-matching step, thus removing "
                  "any discrepancies in brightness or background. Note that, "
                  "although an amazing feature of Montage, this makes the "
                  "assembling of the images take remarkably longer.")

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

    msg = "%sMaking sure the %d input paths are FITS images..."
    print msg % (style.prefix, len(input_paths))

    methods.show_progress(0.0)
    for index, path in enumerate(input_paths):
        # fitsimage.FITSImage.__init__() raises fitsimage.NonStandardFITS if
        # one of the paths is not a standard-conforming FITS file. We do not
        # need the FITSImage object that is created.
        try:
            fitsimage.FITSImage(path)
        except fitsimage.NonStandardFITS:
            print
            msg = "'%s' is not a standard FITS file"
            raise fitsimage.NonStandardFITS(msg % path)

        percentage = (index + 1) / len(input_paths) * 100
        methods.show_progress(percentage)
    print # progress bar doesn't include newline

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

    kwargs = dict(background_match = options.background_match)
    montage.mosaic(input_dir, output_dir, **kwargs)

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

