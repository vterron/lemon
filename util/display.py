#! /usr/bin/env python

# Copyright (c) 2019 Victor Terron. All rights reserved.
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

import functools
import sys
import time

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
