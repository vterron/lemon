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
import astromatic
import customparser
import defaults
import keywords
import fitsimage
import methods
import seeing
import style

description = """
This module does astrometry on an image entirely using Emmanuel Bertin's tools,
namely SExtractor, SCAMP and SWarp. First, SExtractor is run on the image and
the output catalog saved in the FITS_LDAC binary format. This is then read by
SCAMP, which computes the astrometic solution and saves it to a FITS-like image
header. Finally, this header file is merged with the original FITS image using
SWarp, thus updating it with the astrometric information.

"""

ASTROMETRY_COMMAND = 'solve-field'

class AstrometryNetNotInstalled(StandardError):
    """ Raised if Astrometry.net is not installed on the system """
    pass

class AstrometryNetError(subprocess.CalledProcessError):
    """ Raised if the execution of Astrometry.net fails """
    pass

def astrometry_net(path):
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
    try:
        subprocess.check_call(args)
        return output_path
    except subprocess.CalledProcessError, e:
        raise AstrometryNetError(e.returncode, e.cmd)
    finally:
        shutil.rmtree(output_dir, ignore_errors = True)

def astrometry(img_path, scale, equinox, radecsys, saturation,
               copy_keywords = None, ra_keyword = 'RA',
               dec_keyword = 'DEC', stdout = None, stderr = None):
    """ Do astrometry on a FITS image.

    This method chains the execution of SExtractor, SCAMP and SWarp, as
    explained in the description of this module, computing the astrometric
    solution of the input FITS image. Returns the path to the output image,
    which is saved to a temporary file and for whose deletion when it is no
    longer needed the user is responsible.

    scale - scale of the image, in degrees per pixel
    equinox - equinox in years (e.g., 2000)
    radecsys - reference system (e.g., ICRS)
    saturation - number of ADUs at which arises saturation. Note that for
                 coadded images this value is the result of multiplying the CCD
                 saturation level by the number of images that were coadded.

    Keyword arguments:
    copy_keywords - FITS keywords, and their values, to propagate from the
                    input image header to the resampled and coadded image
                    header produced by SWarp.  Needed since not all FITS
                    keywords are automatically copied to the output image
                    header, as many of them become irrelevant.
    ra_keyword - FITS keyword for the right ascension, in decimal degrees.
    dec_keyword - FITS keyword for the declination, in decimal degrees.
    stdout - the SExtractor, SCAMP and SWarp standard output file handle.
             If set to None, no redirection will occur.
    stderr - the SExtractor, SCAMP and SWarp standard error file handle.
             If set to None, no redirection will occur.

    """

    # This does astrometry on the image and returns the path to the
    # temporary file to which the '.head' file has been saved
    head_path = astromatic.scamp(img_path, scale, equinox,
                                 radecsys, saturation,
                                 ra_keyword = ra_keyword,
                                 dec_keyword = dec_keyword,
                                 stdout = stdout, stderr = stderr)
    try:
        # Now merge the '.head' file with the original image
        return astromatic.swarp(img_path, head_path,
                                copy_keywords = copy_keywords,
                                stdout = stdout, stderr = stderr)
    finally:
        try: os.unlink(head_path)
        except (IOError, OSError): pass


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

# CCD SITE#2b focus scale = 0.50209205 arcsec/pixel
parser.add_option('--scale', action = 'store', type = 'float',
                  dest = 'scale', default = 0.50,
                  help = "the scale of the images, in arcsec/pixel "
                  "[default: %default]")

parser.add_option('--equinox', action = 'store', type = 'int',
                  dest = 'equinox', default = '2000',
                  help = "mean equinox, in years [default: %default]")

parser.add_option('--radecsys', action = 'store', type = 'str',
                  dest = 'radecsys', default = 'ICRS',
                  help = "WCS astrometric system [default: %default]")

parser.add_option('--maximum', action = 'store', type = 'int',
                  dest = 'maximum', default = defaults.maximum,
                  help = defaults.desc['maximum'])

# Note for developers: we are not doing anything with the stars of the image,
# so the width of the margin is irrelevant here. However, it is required by the
# __init__ method of the FITSeeingImage class. We could use any value, but we
# prefer to use the same as in other stages of the pipeline.
parser.add_option('--margin', action = 'store', type = 'int',
                  dest = 'margin', default = defaults.margin,
                  help = defaults.desc['margin'])

parser.add_option('-v', '--verbose', action = 'count',
                  dest = 'verbose', default = defaults.verbosity,
                  help = defaults.desc['verbosity'])

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

key_group.add_option('--timek', action = 'store', type = 'str',
                     dest = 'timek', default = keywords.timek,
                     help = keywords.desc['timek'])

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
                     help = keywords.desc['fwhmk'])
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

    print "%sReading WCS information from the FITS header..." % style.prefix,
    img = seeing.FITSeeingImage(img_path, options.maximum, options.margin)
    ra  = img.read_keyword(options.rak)
    dec = img.read_keyword(options.deck)
    print 'done.'

    print "%sAstrometic system: %s" % (style.prefix, options.radecsys)
    print "%sRight ascension = %f (%.2dh %.2dm %.4fs)" % \
          ((style.prefix, ra) + methods.DD_to_HMS(ra))
    print "%sDeclination = %f (%.2d° %.2d′ %.4f″)" % \
          ((style.prefix, dec) + methods.DD_to_DMS(dec))

    print "%sMean equinox: %d" % (style.prefix, options.equinox)
    print "%sScale: %.3f arcsec/pixel" % (style.prefix, options.scale)
    msg = "%sImage saturation level: %d x %d = %d ADUs"
    print msg % (style.prefix, options.maximum, img.ncoadds, img.saturation)

    # All these keywords have to be propagated to the resulting FITS
    # image, as subsequent modules of the pipeline need access to them.
    propagated = \
        [options.objectk, options.filterk, options.rak, options.deck,
         options.datek, options.timek, options.exptimek, options.airmassk,
         keywords.coaddk, options.gaink, options.uncimgk, options.fwhmk]

    print "%sRunning SExtractor, SCAMP and SWarp on the image..." % \
          style.prefix ,
    sys.stdout.flush()

    with open(os.devnull, 'wt') as fd:
        output_path = astrometry_net(img_path)
        try:
            shutil.move(output_path, options.output_path)
        except (IOError, OSError):
            try: os.unlink(output_path)
            except (IOError, OSError): pass

    print 'done.'
    output_img = fitsimage.FITSImage(options.output_path)

    msg1 = "Astrometry done by LEMON on %s" % methods.utctime()
    msg2 = "[Astrometry] Implemented using Emmanuel Bertin's SCAMP and SWarp"
    msg3 = "[Astrometry] Astrometric system: %s" % options.radecsys
    msg4 = "[Astrometry] Mean equinox: %d" % options.equinox
    msg5 = "[Astrometry] Original image: %s" % img_path

    output_img.add_history(msg1)
    output_img.add_history(msg2)
    output_img.add_history(msg3)
    output_img.add_history(msg4)
    output_img.add_history(msg5)

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

