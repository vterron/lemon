#! /usr/bin/env python
#encoding:UTF-8

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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division

import atexit
import logging
import optparse
import os
import shutil
import subprocess
import sys
import tempfile

# LEMON modules
import customparser
import defaults
import fitsimage
import methods
import style

description = """
This module uses a local build of the Astrometry.net software in order to
compute the astrometric solution of a FITS image. This is, essentially, a mere
simple interface to solve-field, Astrometry.net's command-line high-level user
interface, which must be present in PATH. Keep in mind that for Astrometry.net
to work it is also necessary to download the index files.

"""

ASTROMETRY_COMMAND = 'solve-field'

class AstrometryNetNotInstalled(StandardError):
    """ Raised if Astrometry.net is not installed on the system """
    pass

class AstrometryNetError(subprocess.CalledProcessError):
    """ Raised if the execution of Astrometry.net fails """
    pass

def astrometry_net(path, ra = None, dec = None, radius = 1, verbosity = 0):
    """ Do astrometry on a FITS image using Astrometry.net.

    Use a local build of the amazing Astrometry.net software [1] in order to
    compute the astrometric solution of a FITS image. This software has many,
    many advantages over the well-respected SCAMP, but the most important one
    is that it is a blind astrometric calibration service. We do not need to
    know literally anything about the image, including approximate coordinates,
    scale and equinox. It just works, giving us a new FITS file containing the
    WCS header.

    In order for this function to work, you must have built and installed the
    Astrometry.net code in your machine [2]. The main high-level command-line
    user interface, 'solve-field', is expected to be available in your PATH;
    otherwise, the AstrometryNetNotInstalled exception is raised. Note that you
    also need to download the appropriate index files, which are considerably
    heavy. At the time of this writing, the entire set of indexes built from
    the 2MASS catalog [4] has a total size of ~32 gigabytes.

    [1] http://astrometry.net/
    [2] http://astrometry.net/doc/build.html
    [3] http://astrometry.net/doc/readme.html#getting-index-files
    [4] http://data.astrometry.net/4200/

    Keyword arguments:

    ra,
    dec,
    radius - restrict the Astrometry.net search to those indexes within
             'radius' degrees of the field center given by ('ra', 'dec').
             Both the right ascension and declination must be given in order
             for this feature to work. The three arguments must be expressed
             in degrees.

    verbosity - the verbosity level. The higher this value, the 'chattier'
                Astrometry.net will be. Most of the time, a verbosity other
                than zero, the default value, is only needed for debugging.

    """

    emsg = "'%s' not found in the current environment"
    if not methods.which(ASTROMETRY_COMMAND):
        raise AstrometryNetNotInstalled(emsg % ASTROMETRY_COMMAND)

    img = fitsimage.FITSImage(path)
    tempfile_prefix = '%s_' % img.basename_woe
    # Place all output files in this directory
    kwargs = dict(prefix = tempfile_prefix, suffix = '_astrometry.net')
    output_dir = tempfile.mkdtemp(**kwargs)

    # Path to the temporary FITS file containing the WCS header
    root, ext = os.path.splitext(img.basename)
    kwargs = dict(prefix = '%s_astrometry_' % root, suffix = ext)
    with tempfile.NamedTemporaryFile(**kwargs) as fd:
        output_path = fd.name

    # --dir: place all output files in the specified directory.
    # --no-plots: don't create any plots of the results.
    # --new-fits: the new FITS file containing the WCS header.
    # --no-fits2fits: don't sanitize FITS files; assume they're already valid.
    # --overwrite: overwrite output files if they already exist.

    args = [ASTROMETRY_COMMAND, path,
            '--dir', output_dir,
            '--no-plots',
            '--new-fits', output_path,
            '--no-fits2fits',
            '--overwrite']

    # -3 / --ra <degrees or hh:mm:ss>: only search in indexes within 'radius'
    # of the field center given by 'ra' and 'dec'
    # -4 / --dec <degrees or [+-]dd:mm:ss>: only search in indexes within
    # 'radius' of the field center given by 'ra' and 'dec'
    # -5 / --radius <degrees>: only search in indexes within 'radius' of the
    # field center given by ('ra', 'dec')

    if ra is not None:
        args += ['--ra', '%f' % ra]

    if dec is not None:
        args += ['--dec', '%f' % dec]

    if radius is not None:
        args += ['--radius', '%f' % radius]

    # -v / --verbose: be more chatty -- repeat for even more verboseness
    if verbosity:
        args.append('-%s' % ('v' * verbosity))

    try:
        subprocess.check_call(args)
        return output_path
    except subprocess.CalledProcessError, e:
        raise AstrometryNetError(e.returncode, e.cmd)
    finally:
        shutil.rmtree(output_dir, ignore_errors = True)


parser = customparser.get_parser(description)
parser.usage = "%prog [OPTION]... FITS_IMAGE"

parser.add_option('--output', action = 'store', type = 'str',
                  dest = 'output_path', default = 'astrometry.fits',
                  help = "path to the output image [default: %default]")

parser.add_option('--overwrite', action = 'store_true', dest = 'overwrite',
                  help = "overwrite output image if it already exists")

parser.add_option('--update', action = 'store_true', dest = 'update',
                  help = "do not output a FITS image; instead, update the "
                  "input image with the astrometric solution. This causes "
                  "the --output and --overwrite options to be ignored")

parser.add_option('-v', '--verbose', action = 'count',
                  dest = 'verbose', default = defaults.verbosity,
                  help = defaults.desc['verbosity'] + " The verbosity "
                  "level is also passed down to Astrometry.net, causing "
                  "it to be increasingly chattier as more -v flags are "
                  "given")

coords_group = optparse.OptionGroup(parser, "Approximate coordinates",
               "Although one of its main advantages is that it is a blind "
               "calibration tool, the execution of Astrometry.net can be "
               "enormously sped up (in more than an order of magnitude, in "
               "fact) if we know the approximate coordinates of the image. "
               "If that is your case, you may use these options to restrict "
               "the search to those indexes close to the field center")

coords_group.add_option('--ra', action = 'store', type = 'float',
                        dest = 'ra', help = "Right ascension, in degrees")

coords_group.add_option('--dec', action = 'store', type = 'float',
                        dest = 'dec', help = "Declination, in degrees")

parser.add_option_group(coords_group)
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
        img_path = args[0]

    # Images cannot be directly updated with the astrometric solution. Instead,
    # what we do is to save it to a temporary file and then overwrite the input
    # image. This is what the user views as an "update" of the original file.
    if options.update:
        tmp_fd, tmp_path = tempfile.mkstemp(suffix = '.fits')
        os.close(tmp_fd)
        options.output_path = tmp_path
        options.overwrite = True

        # Use a simple cleanup function to guarante that the temporary file is
        # removed -- we want it to be deleted even if the program crashes or
        # the execution is aborted before the input image is overwritten.
        @atexit.register
        def clean_tempfile():
            if os.path.exists(options.output_path):
                try:
                    os.unlink(options.output_path)
                except (IOError, OSError):
                    pass

    elif os.path.exists(options.output_path) and not options.overwrite:
        print "%sError. The output image '%s' already exists." % \
              (style.prefix, options.output_path)
        print style.error_exit_message
        return 1


    print "%sInput FITS image: %s" % (style.prefix, img_path)
    msg = "%sUsing a local build of Astrometry.net to solve the FITS image."
    print msg % style.prefix

    msg = "%sLines not starting with '%s' come from Astrometry.net."
    print msg % (style.prefix, style.prefix.strip())
    print

    kwargs = dict(ra = options.ra,
                  dec = options.dec,
                  verbosity = options.verbose)

    output_path = astrometry_net(img_path, **kwargs)

    try:
        shutil.move(output_path, options.output_path)
    except (IOError, OSError):
        try: os.unlink(output_path)
        except (IOError, OSError): pass

    output_img = fitsimage.FITSImage(options.output_path)

    msg1 = "Astrometry done via LEMON on %s" % methods.utctime()
    msg2 = "[Astrometry] WCS solution found by Astrometry.net"
    msg3 = "[Astrometry] Original image: %s" % img_path

    output_img.add_history(msg1)
    output_img.add_history(msg2)
    output_img.add_history(msg3)

    if not options.update:
        print "%sImage with astrometry saved to '%s'." % \
              (style.prefix, options.output_path)
    else:
        shutil.move(options.output_path, img_path)
        print "%sImage '%s' was updated with the astrometric solution." % \
              (style.prefix, img_path)

    print "%sYou're done ^_^" % style.prefix
    return 0

if __name__ == "__main__":
    sys.exit(main())

