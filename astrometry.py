#! /usr/bin/env python2
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

import logging
import multiprocessing
import optparse
import os
import os.path
import shutil
import sys
import tempfile
import time
import warnings

# The 'timeout' argument of subprocess.call() was added in version 3.3.
# In previous versions we need to use 'subprocess32', a backport of the
# subprocess module from Python 3.2/3.3 for use on 2.x.

if sys.version_info < (3, 3):
    import subprocess32 as subprocess
else:
    import subprocess

# LEMON modules
import customparser
import defaults
import fitsimage
import keywords
import methods
import style

description = """
This module uses a local build of the Astrometry.net software in order to
compute the astrometric solution of the input FITS files, saving the new files,
containing the WCS header, to the output directory. This is, in essence, a mere
simple interface to solve-field, Astrometry.net's command-line high-level user
interface, which must be present in PATH. Keep in mind that for Astrometry.net
to work it is also necessary to download the index files.

"""

ASTROMETRY_COMMAND = 'solve-field'

# The Queue is global -- this works, but note that we could have
# passed its reference to the function managed by pool.map_async.
# See http://stackoverflow.com/a/3217427/184363
queue = methods.Queue()

class AstrometryNetNotInstalled(StandardError):
    """ Raised if Astrometry.net is not installed on the system """
    pass

class AstrometryNetError(subprocess.CalledProcessError):
    """ Raised if the execution of Astrometry.net fails """
    pass

class AstrometryNetUnsolvedField(subprocess.CalledProcessError):
    """ Raised if Astrometry.net could not solve the field """

    def __init__(self, path):
        self.path = path

    def __str__(self):
        return "%s: could not solve field" % self.path

class AstrometryNetTimeoutExpired(AstrometryNetUnsolvedField):
    """ Raised if the Astrometry.net timeout was reached """

    def __init__(self, path, timeout):
        self.path = path
        self.timeout = timeout

    def __str__(self):
        msg = "%s: could not solve field in less than %d seconds"
        return msg % (self.path, self.timeout)

def astrometry_net(path, ra = None, dec = None, radius = 1,
                   verbosity = 0, timeout = None, options = None):
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

    Raises AstrometryNetError if Astrometry.net exits with a non-zero status
    code, AstrometryNetTimeoutExpired if the 'timeout' limit is exceeded and
    AstrometryNetUnsolvedField if the CPU time limit, set in the backend.cfg
    file (by default located in /usr/local/astrometry/etc/) is hit.

    [1] http://astrometry.net/
    [2] http://astrometry.net/doc/build.html
    [3] http://astrometry.net/doc/readme.html#getting-index-files
    [4] http://data.astrometry.net/4200/
    [5] https://groups.google.com/d/msg/astrometry/ORVkOk0jSZg/PeCMeAJodyAJ

    Keyword arguments:

    ra,
    dec,
    radius - restrict the Astrometry.net search to those indexes within
             'radius' degrees of the field center given by ('ra', 'dec').
             Both the right ascension and declination must be given in order
             for this feature to work. The three arguments must be expressed
             in degrees.
    verbosity - the verbosity level. The default value is zero, meaning that
                the function executes silently. A value of one makes both the
                standard output and standard error of Astrometry.net visible.
                Above that, the number of -v flags send to it equals the value
                of the argument minus one. For example: verbosity = 3 allows us
                to see stdout and stderr, and calls Astrometry.net with two -v
                flags. Most of the time, verbosities greater than one are only
                needed for debugging.
    timeout - the maximum number of seconds that Astrometry.net spends on the
              image before giving up and raising AstrometryNetTimeoutExpired.
              Note that the backend configuration file (astrometry.cfg) puts a
              limit on the CPU time that is spent on an image: this can reduce
              that value but not increase it.
    options - a dictionary, containing additional options to be passed to
              solve-field. Each option must map to the corresponding argument
              (for example, {'--downsample' : '2'}), except in case they do not
              take any, when they must map to None (e.g., {'--invert' : None}).
              Both options and values should be given as strings, but they will
              be automatically cast to string just to be safe.

    """

    emsg = "'%s' not found in the current environment"
    if not methods.which(ASTROMETRY_COMMAND):
        raise AstrometryNetNotInstalled(emsg % ASTROMETRY_COMMAND)

    basename = os.path.basename(path)
    root, ext = os.path.splitext(basename)
    # Place all output files in this directory
    kwargs = dict(prefix = root + '_', suffix = '_astrometry.net')
    output_dir = tempfile.mkdtemp(**kwargs)

    # Path to the temporary FITS file containing the WCS header
    kwargs = dict(prefix = '%s_astrometry_' % root, suffix = ext)
    with tempfile.NamedTemporaryFile(**kwargs) as fd:
        output_path = fd.name

    # If the field solved, Astrometry.net creates a <base>.solved output file
    # that contains (binary) 1. That is: if this file does not exist, we know
    # that an astrometric solution could not be found.
    solved_file = os.path.join(output_dir, root + '.solved')

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

    # -v / --verbose: be more chatty -- repeat for even more verboseness. A
    # value of 'verbosity' equal to zero means that both the standard output
    # and error of Astrometry.net and redirected to the null device. Above
    # that, we send 'verbosity' minus one -v flags to Astrometry.net.

    if verbosity > 1:
        args.append('-%s' % ('v' * (verbosity - 1)))

    # If additional options for solve-field have been specified, append them to
    # the argument list. All options are assumed to take an argument, except if
    # they are mapped to None. In this manner, {'--downsample' : 2} is appended
    # to the argument list as ['--downsample', '2'] (note the automatic cast to
    # string), while {'--invert' : None} appends only '--invert'.

    if options:
        for opt, value in options.iteritems():
            opt = str(opt)
            if value is None:
                args.append(opt)
            else:
                args += [opt, str(value)]

    # Needed when 'verbosity' is 0
    null_fd = open(os.devnull, 'w')

    try:
        kwargs = dict(timeout = timeout)
        if not verbosity:
            kwargs['stdout'] = kwargs['stderr'] = null_fd

        subprocess.check_call(args, **kwargs)

        # .solved file must exist and contain a binary one
        with open(solved_file, 'rb') as fd:
            if ord(fd.read()) != 1:
                raise AstrometryNetUnsolvedField(path)

        return output_path

    except subprocess.CalledProcessError, e:
        raise AstrometryNetError(e.returncode, e.cmd)
    # If .solved file doesn't exist or contain one
    except (IOError, AstrometryNetUnsolvedField):
        raise AstrometryNetUnsolvedField(path)
    except subprocess.TimeoutExpired:
        raise AstrometryNetTimeoutExpired(path, timeout)
    finally:
        null_fd.close()
        methods.clean_tmp_files(output_dir)

@methods.print_exception_traceback
def parallel_astrometry(args):
    """ Function argument of map_async() to do astrometry in parallel.

    This will be the first argument passed to multiprocessing.Pool.map_async(),
    which chops the iterable into a number of chunks that are submitted to the
    process pool as separate tasks. 'args' must be a three-element tuple with
    (1) a string with the path to the FITS image, (2) a string with the path to
    the output directory and (3) 'options', the optparse.Values object returned
    by optparse.OptionParser.parse_args().

    This function does astrometry on each FITS image with the astrometry_net()
    function. The output FITS files, containing the WCS headers calculated by
    Astrometry.net, are written to the output directory with the same basename
    as the original files but with the string options.suffix appended before
    the file extension.

    The path to each solved image is put, as a string, into the module-level
    'queue' object, a process shared queue. If the image cannot be solved, None
    is put instead. Note that the contents of the shared queue are necessary so
    that the progress bar can be updated to reflect the number of input images
    that have been processed so far. Apart from that, you most probably do not
    need to do anything with these paths, as the output files are written to
    the output directory by astrometry_net().

    """

    path, output_dir, options = args

    img = fitsimage.FITSImage(path)
    # Add the suffix to the basename of the FITS image
    root, ext = os.path.splitext(os.path.basename(path))
    output_filename = root + options.suffix + ext
    dest_path = os.path.join(output_dir, output_filename)

    if options.blind:
        msg = "%s: solving the image blindly (--blind option)"
        logging.debug(msg % img.path)
        ra = dec = None
        msg = "%s: using α = δ = None"
        logging.debug(msg % img.path)

    else:

        try:
            ra  = img.ra (options.rak)
            dec = img.dec(options.deck)
        except (ValueError, KeyError), e:
            msg = "%s: %s" % (img.path, str(e))
            logging.debug(msg)
            ra = dec = None
            msg = "%s: could not read coordinates from FITS header"
            logging.debug(msg % img.path)
            msg = "%s: using α = δ = None"
            logging.debug(msg % img.path)

    kwargs = dict(ra = ra,
                  dec = dec,
                  radius = options.radius,
                  verbosity = options.verbose,
                  timeout = options.timeout,
                  options = options.solve_field_options)

    try:
        output_path = astrometry_net(img.path, **kwargs)

    except AstrometryNetUnsolvedField, e:

        # A subclass of AstrometryNetUnsolvedField
        if isinstance(e, AstrometryNetTimeoutExpired):
            msg = "%s exceeded the timeout limit. Ignored."
        else:
            msg = "%s did not solve. Ignored."

        msg %= img.path
        warnings.warn(msg, RuntimeWarning)
        queue.put(None)
        logging.debug("%s: None put into global queue" % path)
        return

    try:
        shutil.move(output_path, dest_path)
        logging.debug("%s: solved image saved to %s" % (path, dest_path))
    except (IOError, OSError), e:
        logging.debug("%s: can't solve image (%s)" % (path, str(e)))
        methods.clean_tmp_files(output_path)

    output_img = fitsimage.FITSImage(dest_path)

    debug_args = path, output_img.path
    logging.debug("%s: updating header of output image (%s)" % debug_args)
    msg1 = "Astrometry done via LEMON on %s" % methods.utctime()
    msg2 = "[Astrometry] WCS solution found by Astrometry.net"
    msg3 = "[Astrometry] Original image: %s" % img.path

    output_img.add_history(msg1)
    output_img.add_history(msg2)
    output_img.add_history(msg3)
    logging.debug("%s: header of output image (%s) updated" % debug_args)

    queue.put(output_img.path)
    msg = "{0}: astrometry result ({1!r}) put into global queue"
    logging.debug(msg.format(*debug_args))


parser = customparser.get_parser(description)
parser.usage = "%prog [OPTION]... INPUT_IMGS... OUTPUT_DIR"

parser.add_option('--radius', action = 'store', type = 'float',
                  dest = 'radius', default = 1,
                  help = "only search in indexes within this number of "
                  "degrees of the field center, whose coordinates are read "
                  "from the FITS header of each image (see --rak and --deck). "
                  "In case these keywords cannot be read or contain invalid "
                  "values, the image is solved blindly, as if the --blind "
                  "option had been used. [default: %default]")

parser.add_option('--blind', action = 'store_true', dest = 'blind',
                  help = "ignore --radius, --rak and --deck and solve the "
                  "images blindly. A necessity in case the FITS headers of "
                  "your data have no information about the telescope "
                  "pointing, or when they do but it is deemed to be "
                  "entirely unreliable")

parser.add_option('--timeout', action = 'store', type = 'int',
                  dest = 'timeout', default = 600,
                  help = "the maximum number of seconds that may be spent "
                  "attempting to find the astrometric solution of a FITS "
                  "image. If this time limit is exceeded, we give up and no "
                  "solution for the image is saved to the output directory. "
                  "Note, however, that Astrometry.net's backend configuration "
                  "file (astrometry.cfg) puts a limit on the CPU time that is "
                  "spent on an image: this option can reduce this value but "
                  "not increase it. [default: %default]")

parser.add_option('--suffix', action = 'store', type = 'str',
                  dest = 'suffix', default = 'a',
                  help = "string to be appended to output images, before "
                  "the file extension, of course [default: %default]")

parser.add_option('--cores', action = 'store', type = 'int',
                  dest = 'ncores', default = defaults.ncores,
                  help = defaults.desc['ncores'])

parser.add_option('-o', action = 'callback', type = 'str',
                  dest = 'solve_field_options', default = {},
                  callback = customparser.additional_options_callback,
                  help = "additional options to pass to Astrometry.net's "
                  "solve-field. If the option and the corresponding value, "
                  "if any, contain any whitespace they must be enclosed in "
                  "quotes. For example: -o=--invert, -o '--downsample 2', "
                  "-o --sigma=3. This option may be used multiple times.")

parser.add_option('-v', '--verbose', action = 'count',
                  dest = 'verbose', default = defaults.verbosity,
                  help = defaults.desc['verbosity'] + " By default, the "
                  "standard output and error of Astrometry.net are ignored. "
                  "The first -v flag is necessary to be able to see them; the "
                  "rest are passed down to Astrometry.net, causing it to be "
                  "increasingly chattier as more -v flags are given.")

key_group = optparse.OptionGroup(parser, "FITS Keywords",
                                 keywords.group_description)

key_group.add_option('--rak', action = 'store', type = 'str',
                     dest = 'rak', default = keywords.rak,
                     help = keywords.desc['rak'])

key_group.add_option('--deck', action = 'store', type = 'str',
                     dest = 'deck', default = keywords.deck,
                     help = keywords.desc['deck'])

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

    # Print the help and abort the execution if there are not two positional
    # arguments left after parsing the options, as the user must specify at
    # least one (only one?) input FITS file and the output directory
    if len(args) < 2:
        parser.print_help()
        return 2     # 2 is generally used for command line syntax errors
    else:
        input_paths = args[:-1]
        output_dir = args[-1]

    # No index can be within the search area if the radius is not > 0
    if options.radius <= 0:
        msg = "%sError: --radius must a positive number of degrees"
        print msg % style.prefix
        sys.exit(style.error_exit_message)

    # Make sure that the output directory exists; create it if it doesn't.
    methods.determine_output_dir(output_dir)

    print "%sUsing a local build of Astrometry.net." % style.prefix
    msg = "%sDoing astrometry on the %d paths given as input."
    print msg % (style.prefix, len(input_paths))

    pool = multiprocessing.Pool(options.ncores)
    map_async_args = ((path, output_dir, options) for path in input_paths)
    result = pool.map_async(parallel_astrometry, map_async_args)

    while not result.ready():
        time.sleep(1)
        methods.show_progress(queue.qsize() / len(input_paths) * 100)
        # Do not update the progress bar when debugging; instead, print it
        # on a new line each time. This prevents the next logging message,
        # if any, from being printed on the same line that the bar.
        if logging_level < logging.WARNING:
            print

    result.get() # reraise exceptions of the remote call, if any
    methods.show_progress(100) # in case the queue was ready too soon
    print

    # Results in the process shared queue were only necessary to accurately
    # update the progress bar. They are no longer needed, so empty it now.
    queue.clear()

    print "%sYou're done ^_^" % style.prefix
    return 0

if __name__ == "__main__":
    sys.exit(main())

