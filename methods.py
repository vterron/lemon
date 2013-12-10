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
import numpy
import os
import os.path
import stat
import sys
import tempfile
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
    included a newline character either. The percentage must be in the range
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

def deprecated(func):
    """ Mark a function as deprecated.

    This is a decorator which can be used to mark functions as deprecated.
    It results in a DeprecationWarning being issued when the function is used.
    [URL] http://code.activestate.com/recipes/391367-deprecated/

    """

    @functools.wraps(func)
    def depre_func(*args, **kwargs):
        msg = "call to deprecated function '%s'." % func.__name__
        warnings.warn(msg, category = DeprecationWarning, stacklevel = 2)
        return func(*args, **kwargs)

    return depre_func

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


class WrongIntervalError(ValueError):
    pass

def intervals_split_callback(option, opt, value, parser):
    """ opt-parse callback function to parse a list of intervals.

    Intended to be used to parse opt-parse options that contain list of
    intervals, the method receives, comma-separated, dash-separated integers
    and returns the values as a list of two-element tuples. '1-4,5-6', for
    example, returns [(1,4),(5,6)].

    As a special case, the one-character '*' string would translate into (-inf,
    inf), which means that all numbers belong to the interval. Also, missing
    values are understood as positive or negative infinity, depending on
    whether it is the lower or upper bound of the range. Thus, '-5' translates
    into (-inf, 5), while '6-' is equivalent to (6, inf). This means that '-'
    is a synonym for '*', as both are equal to (-inf, inf).

    The WrongIntervalError exception is raised for improperly-formatted
    intervals and also for ranges which overlap with the preceding one.

    """

    intervals = []
    numbers = value.split(',')

    if len(numbers) == 1 and numbers[0] == '*':
        intervals.append((float('-inf'), float('inf'))) # all the elements
    else:
        for pair in numbers:
            try:
                # May raise "need more than 1 value to unpack" (ValueError)
                lower, upper = pair.split('-')
            except ValueError:
                raise WrongIntervalError(pair)

            lower = float('-inf') if not lower else int(lower)
            upper = float('inf') if not upper else int(upper)

            if lower > upper:
                msg = "error in '%s': upper bound must be greater " \
                      "than or equal to the lower bound" % pair
                raise WrongIntervalError(msg)

            if intervals and lower <= intervals[-1][-1]:
                msg = "error in '%s': cannot overlap with preceding " \
                      "interval" % pair
                raise WrongIntervalError(msg)

            intervals.append((lower, upper))

    setattr(parser.values, option.dest, intervals)

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
    """ Load a list of Pixels from a file.

    The method parses a text file which shall contain two values (the x and y
    coordinates of a star) per line. These pixels are then returned in a list
    of tuples. Note that improperly-formatted lines, such as those with tree
    values or a non-real value, are ignored.

    """

    coordinates = []
    with open(path, 'rt') as fd:
        for line in fd:

            words = line.split()

            if not words:
                continue

            if len(words) != 2:
                msg = ("Unable to parse line '%r'. Objects must be listed one "
                       "per line with coordinate values in columns one (right "
                       "ascension) and two (declination)" % line)
                raise ValueError(msg)

            try:
                ra = float(words[0])
                if not 0 <= ra < 360:
                    raise ValueError
            except ValueError:
                msg = "Right ascension '%r' not in range [0, 360[ degrees"
                raise ValueError(msg % ra)

            try:
                dec = float(words[1])
                if not -90 <= dec <= 90:
                    raise ValueError
            except ValueError:
                msg = "Declination '%r' not in range [-90, 90] degrees"
                raise ValueError(msg % dec)

            coordinates.append((ra, dec))

    return coordinates

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
    """ Minimalistic memoization decorator.
    http://code.activestate.com/recipes/577219-minimalistic-memoization/ """

    cache = {}
    @functools.wraps(f)
    def memf(*x):
        if x not in cache:
            cache[x] = f(*x)
        return cache[x]
    return memf

def int_to_roman(i):
    """ Convert an integer to Roman numerals.

    This is a slightly modified version of the recipe published by Tim Valenta
    at ActiveState Code Recipes. As it was posted in a comment, it did not list
    its exact license in the sidebar, but, according to the FAQ, all of the
    code licenses allowed for recipes are OSI approved open source licenses.
    http://code.activestate.com/recipes/81611-roman-numerals/#c4

    """

    numeral_map = (
      (1000, 'M'), (900, 'CM'), (500, 'D'), (400, 'CD'), (100, 'C'),
      (90, 'XC'), (50, 'L'), (40, 'XL'),  (10, 'X'), (9, 'IX'),
      (5, 'V'), (4, 'IV'), (1, 'I')
    )

    result = []
    for integer, numeral in numeral_map:
        count = int(i / integer)
        result.append(numeral * count)
        i -= integer * count
    return ''.join(result)

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

