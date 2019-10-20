#! /usr/bin/env python

# Copyright (c) 2014 Victor Terron. All rights reserved.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import division

import StringIO
import operator
import os
import random
import time
import warnings


# LEMON modules
from test import unittest
from astromatic import Coordinates
import methods

class MethodsFunctionsTests(unittest.TestCase):

    def test_func_catchall(self):

        # Returns func(*args, **kwargs) ...
        self.assertEqual(3,  methods.func_catchall(operator.div, 9, 3))
        self.assertEqual(4,  methods.func_catchall(int, '4'))
        self.assertEqual(-5, methods.func_catchall(max, -1, -5, 4, key=abs))

        # ... unless the function raises an exception. In that case, it is
        # catched and None is returned instead.

        def foo_except():
            raise ValueError

        self.assertEqual(None, methods.func_catchall(foo_except))
        self.assertEqual(None, methods.func_catchall(operator.div, 1, 0))


class StreamToWarningFilterTest(unittest.TestCase):

    def test_filter_stream(self):

        # A filter that, if the 'foa', 'fooa' or 'foooa' strings are matched,
        # issues UserWarning using the string as the warning message. If not,
        # writes the string to the output stream (here, a StringIO)

        output = StringIO.StringIO()
        expected_output = 'fo{1,3}a'
        regexp = '(?P<msg>{0})'.format(expected_output)
        category = UserWarning
        args = output, regexp, category
        stdout_filter = methods.StreamToWarningFilter(*args)

        with warnings.catch_warnings():
            warnings.filterwarnings('error')

            # 'fooa' matches the regexp, so issue the warning
            with self.assertRaisesRegexp(category, expected_output):
                stdout_filter.write("fooa")

            # 'spam' does not match, so write to the output stream
            str_ = "spam"
            stdout_filter.write(str_)
            self.assertEqual(output.getvalue(), str_)


        # A filter that, if "Warning: Keyword: EXPTIME not found" is matched,
        # issues RuntimeWarning with "EXPTIME not matched" as its message. If
        # not, writes the string to the output stream (again, a StringIO).

        output = StringIO.StringIO()
        expected_output = "EXPTIME not found"
        regexp = "Warning: Keyword: (?P<msg>{0})".format(expected_output)
        category = RuntimeWarning
        args = output, regexp, category
        stdout_filter = methods.StreamToWarningFilter(*args)

        with warnings.catch_warnings():
            warnings.filterwarnings('error')

            with self.assertRaisesRegexp(category, expected_output):
                str_ = "Warning: Keyword: EXPTIME not found"
                stdout_filter.write(str_)

            str_ = str_.replace('EXPTIME', 'OBJECT')
            stdout_filter.write(str_)
            self.assertEqual(output.getvalue(), str_)


        # The example mentioned in the class docstring

        output = StringIO.StringIO()
        expected_output = '2.19.5'
        regexp = 'v(?P<msg>(\d\.?)+)'
        category = UserWarning
        args = output, regexp, category
        stdout_filter = methods.StreamToWarningFilter(*args)

        with warnings.catch_warnings():
            warnings.filterwarnings('error')

            with self.assertRaisesRegexp(category, expected_output):
                stdout_filter.write('v2.19.5')

            str_ = "Nobody expects the Spanish inquisition"
            stdout_filter.write(str_)
            self.assertEqual(output.getvalue(), str_)


        # The regular expression does not include the mandatory 'msg' named
        # group, so IndexError is raised when the string is matched, as the
        # StreamToWarningFilter class tries to refer to a non-existent group.

        output = StringIO.StringIO()
        regexp = 'spam'
        args = output, regexp, UserWarning
        stdout_filter = methods.StreamToWarningFilter(*args)

        with self.assertRaisesRegexp(IndexError, "no such group"):
            stdout_filter.write("spam")

    def test_filter_close(self):

        fd = open(os.devnull, 'wt')
        regexp = "Keyword (?P<msg>EXPTIME) not found"
        args = fd, regexp, RuntimeWarning
        devnull_filter = methods.StreamToWarningFilter(*args)
        # Must close the underlying file object
        devnull_filter.close()
        self.assertTrue(fd.closed)
