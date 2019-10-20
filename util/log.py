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

import logging
import os.path
import re
import sys
import warnings

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
