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
import platform
import re
import shutil
import stat
import sys
import tempfile
import traceback
import time
import warnings

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
