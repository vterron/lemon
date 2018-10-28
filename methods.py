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

import contextlib
import functools
import logging
import math
import multiprocessing
import multiprocessing.queues
import numpy
import os
import os.path
import platform
import re
import shutil
import stat
import sys
import tempfile
import traceback
import time
import warnings

# LEMON modules
import style

def percentage_change(old, new):
    """ Return the relative change between the old value and the new one.

    Note we need to use an absolute value for V1 in the denominator regarding
    values with V1 being a negative and V2 being positive, as well as V1 being
    negative, and V2 being greater than V1 but still negative. """

    if old < 0 and (new > 0 or old < new < 0):
        return (new - old) / float(abs(old))
    else:
        return (new - old) / float(old)

def show_progress(percentage):
    """ Print a progress bar strikingly similar to that of the wget command.

    Displays a progress bar with the format used as the Unix wget command:
    51%[===============================>                                 ]

    The whole bar, including the percentage and the surrounding square
    brackets, has a length of 79 characters, as recommended by the Style Guide
    for Python Code (PEP 8). It should also be noted that the progress bar is
    printed in the current line, so a new one should be started before calling
    the method if the current line is not to be overwritten. It also does not
    include a newline character either. The percentage must be in the range
    [0, 100].

    """

    if not 0 <= percentage <= 100:
        raise ValueError("The value of 'percentage' must be in the range [0,100]")

    length = 79     # the 79 character line recommendation
    occupied = 3 + len(style.prefix)  # "%" + "[" + "]" + style.prefix


    percent = str(int(percentage))
    available = length-occupied-len(percent)
    bar = int((available)*int(percent)/100.0)
    spaces = available-bar

    sys.stdout.write("\r" + style.prefix + percent + "%[" + bar * "=" + \
                     ">" + spaces*" " + "]")
    sys.stdout.flush()

def determine_output_dir(output_directory, dir_suffix = None, quiet = False):
    """ Ensure that the specified directory does exist.

    The method receives the intented output directory and guarantees that it
    can be used. If None was specified, it creates a temporary directory
    (readable, writable, and searchable only by the creating user ID), while if
    the specified directory does not exists, it creates it. The method always
    return the path to the output directory, whether it is the path to the
    just-created temporary directory or that originally specified.

    Keyword arguments:
    dir_suffix - the string to be appended to the output directory if it is to
                 be a temporary directory (that is, if the specified output
                 directory is None.
    quiet - if True, no messages will be printed to standard output.

    """

    # If no output directory was specified, create a temporary directory
    if not output_directory:
        if not quiet:
            print "%sNo output directory was specified." % style.prefix
        if dir_suffix is not None:
            temporary_directory = tempfile.mkdtemp(suffix = dir_suffix)
        else:
            temporary_directory = tempfile.mkdtemp()
        if not quiet:
            print "%sImages will be saved to temporary directory %s" % \
                  (style.prefix, temporary_directory)
        return temporary_directory

    # If the path to the output directory was specified, check that it exists
    # and that it is writable. If it does not exist, try to create it.
    else:
        if not os.path.exists(output_directory):
            try:
                os.makedirs(output_directory)
            except OSError:
                msg = "The output directory '%s' could not be created. " \
                      "Is the directory writable?" % output_directory
                raise IOError(msg)

            if not quiet:
                print "%sThe output directory '%s' did not exist, so it " \
                      "had to be created." % (style.prefix, output_directory)
        else:
            if not os.path.isdir(output_directory):
                raise IOError("%s is not a directory" % output_directory)
            if not os.access(output_directory, os.W_OK):
                raise IOError("The output directory '" + output_directory + "' is not writable.")
            if not quiet:
                print "%sImages will be saved to directory '%s'" % \
                      (style.prefix, output_directory)

        return output_directory

def DMS_to_DD(degrees, arcminutes, arcseconds):
    """ Degrees, arcminutes, arcseconds to decimal degrees conversion. """

    decimal = abs(degrees) + float(arcminutes)/60 + float(arcseconds)/3600
    if degrees < 0:
        decimal = -decimal
    return decimal

def DD_to_DMS(decimal_degrees):
    """ Decimal degrees to degrees, arcminutes, arcseconds conversion. """

    # math.modf(x) returns the fractional and integer parts of x
    arcminutes, degrees = math.modf(decimal_degrees)
    arcminutes = abs(arcminutes) * 60  # do not propagate the minus sign, if any
    arcseconds, arcminutes = math.modf(arcminutes)
    arcseconds *= 60
    return int(degrees), int(arcminutes), arcseconds

def HMS_to_DD(hours, minutes, seconds):
    """ Hours, minutes, seconds to decimal degrees conversion. """

    decimal = abs(hours)*15 + float(minutes)/60*15 + float(seconds)/3600*15
    if hours < 0:
        decimal = -decimal
    return decimal

def DD_to_HMS(decimal_degrees):
    """ Decimal degrees to hours, minutes, seconds conversion. """

    # math.modf(x) returns the fractional and integer parts of x
    minutes, hours = math.modf(decimal_degrees / 15.0)
    hours = abs(hours)    # do not propagate the minus sign, if any
    minutes = abs(minutes) * 60
    seconds, minutes = math.modf(minutes)
    seconds *= 60
    return int(hours), int(minutes), seconds

def ra_str(decimal_degrees):
    """ Return the string representation of a right ascension.
    Example: 14h 03m 12.64s """

    ra_hours, ra_min, ra_sec = DD_to_HMS(decimal_degrees)
    ra_hours = '%02d' % ra_hours
    ra_min = '%02d' % ra_min
    ra_sec = '%05.2f' % ra_sec
    return "%sh %sm %ss" % (ra_hours, ra_min, ra_sec)

def dec_str(decimal_degrees):
    """ Return the string representation of a declination.
    Example: +57d 33m 10.75s """

    dec_deg, dec_arcmin, dec_arcsec = DD_to_DMS(decimal_degrees)
    dec_deg = '%+03d' % dec_deg
    dec_arcmin = '%02d' % dec_arcmin
    dec_arcsec = '%05.2f' % dec_arcsec
    return "%sd %sm %ss" % (dec_deg, dec_arcmin, dec_arcsec)

def str_split_callback(option, opt, value, parser):
    """ opt-parse callback function to parse a list of values.

    This method is intended to be used in order to parse opt-parse options that
    contain a list of values. In other words, the received comma-separated
    values are converted to a list, so that when the user specifies, i.e.,
    '--groups one,two' the value returned by opt-parse is ['one', 'two'].
    [URL] http://stackoverflow.com/questions/392041/python-optparse-list

    option - the Option instance that is calling the callback.
    opt_str - the option string seen on the command-line.
    value - the argument to this option seen on the command-line.
    parser - the OptionParser instance driving the whole thing.

    """

    setattr(parser.values, option.dest, value.split(','))

def owner_writable(path, add):
    """ Make the file owner writeable or unwriteable.

    The method either adds or removes write permission for the file's owner,
    depending on the value of the 'add' parameter. In other words, if 'add' is
    True, the method is equivalent to chmod u+w on the file, while if the
    parameter evaluates to False it is equivalent to a chmod u-w on 'path'.

    Note that doing this is not as straight-forward as you may think, as
    os.chmod(path, mode) _sets_ the mode of path to a numeric mode, which means
    that before doing that we need to determine which the current permissions
    of the file add.

    """

    # stat.ST_MODE returns the protection bits of the file.
    # stat.S_IMODE returns the portion of the file's mode that can be set by
    # os.chmod(), i.e., the file's permission bits, plus the sticky bit,
    # set-group-id, and set-user-id bits (on systems that support them).
    # stat.S_IWUSR = owner has write permission.

    mode = stat.S_IMODE(os.stat(path)[stat.ST_MODE])

    if add:
        mode |= stat.S_IWUSR    # bitwise (inclusive) OR
    else:
        mode ^= stat.S_IWUSR    # bitwise XOR (exclusive OR)

    os.chmod(path, mode)

def load_coordinates(path):
    """ Load a list of celestial coordinates from a text file.

    Parse a text file containing the celestial coordinates of a series of
    astronomical objects, one per line, and return a generator that yields (ra,
    dec, pm_ra, pm_dec) tuples. The file must have two columns, with the right
    ascension and declination (in decimal degrees) and, optionally, two other
    columns with the proper motion in right ascension and declination (in
    seconds of arc per year) surrounded by brackets. For example:

      269.456271 4.665281
      269.452075 4.693391 [-0.79858] [10.32812] # Barnard's Star
      269.466450 4.705625 [0.0036] [-.0064]     # TYC 425-262-1

    The four-element tuples contain the right ascension, declination, proper
    motion in right ascension and proper motion in declination, respectively.
    Nones are used in case the proper motion of an astronomical object is not
    specified. Empty lines, as well as comments (which start with the hash
    character, #, and extend to the end of the physical line) are ignored.

    ValueError is raised (a) if in any line there is a number of coordinates
    other than two (right ascension and declination) or the proper motions are
    not surrounded by brackets, (b) if any right ascension or declination is
    out of range or (c) if the proper motion in right ascension is specified
    but not that in declination, or vice versa.

    """

    with open(path, 'rt') as fd:
        for line in fd:

            # Ignore comments
            line = line.split('#')[0]

            # Ignore empty lines
            regexp = "^\s*$"
            if re.match(regexp, line):
                continue

            kwargs = dict(float = "([+-]?\d+(\.\d+)*)")
            regexp = ("^\s*"
                      "(?P<ra>{float})"
                      "\s+"
                      "(?P<dec>{float})"
                      "("
                        "\s+"
                        "\[\s*(?P<pm_ra>{float})\s*\]"
                        "\s+"
                        "\[\s*(?P<pm_dec>{float})\s*\]"
                      ")?"
                      "\s*$")

            match = re.match(regexp.format(**kwargs), line)

            if not match:
                msg = ("Unable to parse line %r. Astronomical objects must be "
                       "listed one per line with coordinate values in columns one "
                       "(right ascension) and two (declination). Proper motions "
                       "may be optionally specified in columns three (ra) and "
                       "four (dec), surrounded by brackets -- but, in that case, "
                       "both of them are required." % line)
                raise ValueError(msg)

            ra = float(match.group('ra'))
            if not 0 <= ra < 360:
                msg = "Right ascension '%s' not in range [0, 360[ degrees"
                raise ValueError(msg % ra)

            dec = float(match.group('dec'))
            if not -90 <= dec <= 90:
                msg = "Declination '%s' not in range [-90, 90] degrees"
                raise ValueError(msg % dec)

            pm_ra = match.group('pm_ra')
            if pm_ra is not None:
                pm_ra = float(pm_ra)

            pm_dec = match.group('pm_dec')
            if pm_dec is not None:
                pm_dec = float(pm_dec)

            yield ra, dec, pm_ra, pm_dec

def which(*names):
    """ Search PATH for executable files with the given names.

    Replicate the functionality of Unix 'which', returning a list of the full
    paths to the executables that would be executed in the current environment
    if the arguments were given as commands in a POSIX-conformant shell. This
    is done by searching, in the directories listed in the environment variable
    PATH, for executable files matching the names of the arguments. If all the
    command are nonexistent or not executable, an empty list is returned.

    The code in this function is largely ripped from Twister's repository:
    https://twistedmatrix.com/trac/browser/trunk/twisted/python/procutils.py

    """

    result = []
    for directory in os.environ.get('PATH', '').split(os.pathsep):
        for name in names:
            path = os.path.join(directory, name)
            if os.access(path, os.X_OK) and os.path.isfile(path):
                result.append(path)
    return result

def split_by_diff(iterable, delta = 3):
    """ Split a sequence by the difference between consecutive elements.

    The method returns an interator over the result of splitting the input
    sequence, stopping each sub-sequence at the element at which its difference
    with the next one is greater than 'delta'. In other words: the difference
    between consecutive elements of each of the returned sublists will be
    smaller than or equal to delta.

    For example, split_by_diff([1, 2, 3, 8, 9, 15], delta = 3) returns an
    iterator over three lists: [1, 2, 3], [8, 9] and [15]

    """

    differences = numpy.diff(iterable)
    sublist_indexes = numpy.where(differences > delta)[0]

    sublists = []
    iterable = list(iterable)  # work on a copy
    for index in reversed(sublist_indexes):
        sublists.append(iterable[index + 1:])
        del iterable[index + 1:]
    sublists.append(iterable)
    return reversed(sublists)

def memoize(f):
    """ Minimalistic memoization decorator (*args / **kwargs)
    Based on: http://code.activestate.com/recipes/577219/ """

    cache = {}
    @functools.wraps(f)
    def memf(*args, **kwargs):
        fkwargs = frozenset(kwargs.iteritems())
        if (args, fkwargs) not in cache:
            cache[args, fkwargs] = f(*args, **kwargs)
        return cache[args, fkwargs]
    return memf

def utctime(seconds = None, suffix = True):
    """ UTC version of time.ctime.

    Convert a time expressed in seconds since the Unix epoch to a 28-character
    string, representing Coordinated Universal Time, of the following form:
    'Sun Jun 20 23:21:05 1993 UTC'. Fractions of a second are ignored. The last
    part of the string, ' UTC', is omitted (and therefore a 24-character string
    returned) if the 'suffix' argument evaluates to False. If 'seconds' is not
    provided or None, the current time as returned by time.time() is used.

    """

    utc_ctime = time.asctime(time.gmtime(seconds))
    if suffix:
        utc_ctime += ' UTC'
    return utc_ctime

def log_uncaught_exceptions(func):
    """ Decorator to log uncaught exceptions with level DEBUG.

    This decorator catches any exception raised by the decorated function and
    logs it with level DEBUG on the root logger. Only subclasses of Exception
    are caught, as we do not want to log SystemExit or KeyboardInterrupt. The
    usage of this decorator makes probably only sense when the function raising
    the uncaught exception cannot be fixed, for example when working with a
    third-party library.

    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception:
            type, value, traceback = sys.exc_info()
            msg = "%s raised %s('%s')" % (func.__name__, type.__name__, value)
            logging.debug(msg)
    return wrapper

@contextlib.contextmanager
def tmp_chdir(path):
    """ A context manager to temporarily change the working directory.

    This is a rather simple context manager to change the current working
    directory within a with statement, restoring the original one upon exit.

    """

    cwd = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(cwd)

def clean_tmp_files(*paths):
    """ Try to remove multiple temporary files and directories.

    Loop over the provided positional arguments, calling os.unlink() on files
    and shutil.rmtree() on directories. Errors never raise an exception, but
    are logged at DEBUG level. These files are considered to be 'temporary' in
    the sense that, being no longer necessary, they must be cleaned up, but
    they are not important enough as to require special handling if they cannot
    be deleted. After all, if they are located in /tmp/, as they are expected,
    they will eventually get cleared.

    """

    for path in paths:

        if os.path.isdir(path):

            msg = "Cleaning up temporary directory '%s'"
            logging.debug(msg % path)

            error_count = [0]

            def log_error(function, path, excinfo):
                """ Error handler for shutil.rmtree() """

                # nonlocal is not available in Python 2.x so, being it outside
                # of the local scope, we cannot use 'error_count' as a counter
                # and rebind it each time we come across an error. But we can
                # make it a list, which being mutable allows us to modify its
                # elements inside the function.

                error_count[0] += 1
                msg = "%s: error deleting '%s' (%s)"
                args = function, path, excinfo[1]
                logging.debug(msg % args)

            try:
                kwargs = dict(ignore_errors = False, onerror = log_error)
                shutil.rmtree(path, **kwargs)

            finally:
                msg = "Temporary directory '%s' deleted"
                if max(error_count) > 0:
                    msg += " (but there were failed removals)"
                logging.debug(msg % path)

        else:

            msg = "Cleaning up temporary file '%s'"
            logging.debug(msg % path)

            try:
                os.unlink(path)
            except OSError, e:
                msg = "Cannot delete '%s' (%s)"
                logging.debug(msg % (path, e))
            else:
                msg = "Temporary file '%s' removed"
                logging.debug(msg % path)

class SharedCounter(object):
    """ A synchronized shared counter.

    The locking done by multiprocessing.Value ensures that only a single
    process or thread may read or write the in-memory ctypes object. However,
    in order to do n += 1, Python performs a read followed by a write, so a
    second process may read the old value before the new one is written by the
    first process. The solution is to use a multiprocessing.Lock to guarantee
    the atomicity of the modifications to Value.

    This class comes almost entirely from Eli Bendersky's blog:
    http://eli.thegreenplace.net/2012/01/04/shared-counter-with-pythons-multiprocessing/

    """

    def __init__(self, n = 0):
        self.count = multiprocessing.Value('i', n)

    def increment(self, n = 1):
        """ Increment the counter by n (default = 1) """
        with self.count.get_lock():
            self.count.value += n

    @property
    def value(self):
        """ Return the value of the counter """
        return self.count.value


class Queue(multiprocessing.queues.Queue):
    """ A portable implementation of multiprocessing.Queue.

    Because of multithreading / multiprocessing semantics, Queue.qsize() may
    raise the NotImplementedError exception on Unix platforms like Mac OS X
    where sem_getvalue() is not implemented. This subclass addresses this
    problem by using a synchronized shared counter (initialized to zero) and
    increasing / decreasing its value every time the put() and get() methods
    are called, respectively. This not only prevents NotImplementedError from
    being raised, but also allows us to implement a reliable version of both
    qsize() and empty().

    """

    def __init__(self, *args, **kwargs):
        super(Queue, self).__init__(*args, **kwargs)
        self.size = SharedCounter(0)

    def put(self, *args, **kwargs):
        super(Queue, self).put(*args, **kwargs)
        self.size.increment(1)

    def get(self, *args, **kwargs):
        item = super(Queue, self).get(*args, **kwargs)
        self.size.increment(-1)
        return item

    def qsize(self):
        """ Reliable implementation of multiprocessing.Queue.qsize() """
        return self.size.value

    def empty(self):
        """ Reliable implementation of multiprocessing.Queue.empty() """
        return not self.qsize()

    def clear(self):
        """ Remove all elements from the Queue. """
        while not self.empty():
            self.get()

def print_exception_traceback(func):
    """ Decorator to print the stack trace of an exception.

    This decorator catches any exception raised by the decorated function,
    prints its information to the standard output (traceback.print_exc()) and
    re-raises it. In particular, this may be used to decorate functions called
    through the multiprocessing module, since exceptions in spawned child
    processes do not print stack traces.

    The idea for this decorator comes from Chuan Ji's blog:
    http://seasonofcode.com/posts/python-multiprocessing-and-exceptions.html

    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception:
            traceback.print_exc()
            print
            raise
    return wrapper

@contextlib.contextmanager
def tempinput(data):
    """ A context manager to work with StringIO-like temporary files.

    The open() built-in command only takes filenames, so we cannot feed to it a
    StringIO. This context manager writes 'data' to a temporary file, ensuring
    that it is cleaned up afterwards. In this manner, we can work in a similar
    manner to what we would have done with StringIO, as the data is written to
    disk only temporarily and transparently to us.

    with tempinput('some data\more data') as path:
        open(path)

    Taken from Martijn Pieters's answer on Stack Overflow:
    [URL] https://stackoverflow.com/a/11892712/184363

    """

    fd, path = tempfile.mkstemp()
    os.write(fd, data)
    os.close(fd)
    yield path
    os.unlink(path)

def func_catchall(func, *args, **kwargs):
    """ Return func(*args, **kwargs), or None if an exception is raised.

    Return whatever func() returns when called with the positional arguments
    'args' and keyword arguments 'keywords'. If any exception is raised by
    func(), catch it, log it and return None. This is a convenience function
    useful when we need to call a function that may raise an exception,
    scenario in which we want to use None instead of what would have been
    normally returned.

    """

    try:
        return func(*args, **kwargs)
    except Exception as e:
        exc_type = sys.exc_info()[0]
        args = (func.__name__, exc_type.__name__, str(e))
        msg = "%s() raised %s (%s), None returned instead" % args
        logging.debug(msg)
        return None

class StreamToWarningFilter(object):
    """ A file-like class that matches strings and issues them as warnings.

    This class creates a file-like object that normally writes a string to
    'fd', a file type object. However, those lines that match 'regexp' are
    issued as warnings of the 'category' class and logged at INFO level. The
    message used for the warning is that matched by the 'msg' named group that
    must be present in the regular expression; otherwise, IndexError is raised.
    This class may be used (and, in fact, was coded with this purpose in mind),
    for example, to capture a message written by a third-party library to the
    standard error and issue it as a warning instead.

    For example, consider the scenario where the following object is created:
    fd = StreamToWarningFilter(sys.stdout, 'v(?P<msg>(\d\.?)+)', UserWarning)
    After this, fd.write('v2.19.5') will raise UserWarning (with the message
    '2.19.5'), while fd.write('Nobody expects the Spanish inquisition') will
    not match the regexp and therefore print the string to sys.stdout.

    """

    def __init__(self, fd, regexp, category):
        self.fd = fd
        self.regexp = regexp
        self.category = category

    def write(self, str_):
        """ Write str_ to the file; issue as a warning if regexp is matched"""

        match = re.match(self.regexp, str_)
        if match:
            msg = match.group('msg')
            logging.info(msg)
            warnings.warn(msg, self.category)
        else:
            self.fd.write(str_)

    def flush(self):
        self.fd.flush()

    def close(self):
        self.fd.close()

class LoggerWriter(object):
    """ Wrap a logger with a file-like API.

    Sometimes, we need to interface to a third-party API which expects a
    file-like object to write to, but we want to direct the API's output to a
    logger. This can be done using this class, based on that written by Vinay
    Sajip [https://stackoverflow.com/a/9422332/184363].

    """


    def __init__(self, level):
        """ Initialize the LoggerWritter object.

        The argument 'level' must be a valid logging level: 'debug', 'info',
        'warning', 'error' or 'critical', and will correspond to the function
        of the logging module that will be used every time the method write()
        is called. For example, LoggerWriter('info').write(msg) is equivalent
        to logging.info(msg).

        """

        self.log_func = getattr(logging, level)

    def write(self, msg):
        self.log_func(msg)

    def flush(self):
        """ Do nothing -- required by PyRAF.

        This method does nothing, but it is required in order to be able to
        redirect the output of the PyRAF tasks to the logger. Otherwise, the
        AttributeError ('LoggerWriter' object has no attribute 'flush')
        exception is raised.

        """

        pass

def get_nbits():
    """ Return the bit architecture of the Python interpreter binary. """

    bits = platform.architecture()[0]
    for n in (32, 64):
        if str(n) in bits:
            return n

    # The Python documentation warns us that "On Mac OS X (and perhaps other
    # platforms), executable files may be universal files containing multiple
    # architectures. To get at the "64-bitness" of the current interpreter, it
    # is more reliable to query the sys.maxsize attribute".

    if sys.maxsize > 2**32:
        return 64
    else:
        return 32 # safe guess

