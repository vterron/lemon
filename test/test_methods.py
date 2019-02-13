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
import Queue
import random
import time
import warnings


# LEMON modules
from test import unittest
from astromatic import Coordinates
import methods

NITERS = 100  # How many times each test case is run with random data

class LoadCoordinatesTest(unittest.TestCase):

    # A series of two-element tuples, one per line. The first element is a
    # string containing the name of an astronomical object. The second is a
    # four-element tuple with the right ascension, declination and proper
    # motions. Nones are used if the proper motions are not known.
    TEST_DATA_DIR = './test/test_data/SIMBAD_objects'

    # Parse the SIMBAD file and map each astronomical object (a string) to
    # its right ascension and declination (an astromatic.Coordinates object)
    COORDINATES = {}
    with open(TEST_DATA_DIR, 'rt') as fd:
        for line in fd:
            line = line.strip()
            if line and not line.startswith("#"):
                object_, coords = eval(line)
                COORDINATES[object_] = Coordinates(*coords)

    NCOORDS = (1, len(COORDINATES))  # number of objects in each file
    NEMPTY = (1, 50)    # minimum and maximum number of empty lines
    NCOMMENTS = (1, 50) # minimum and maximum number of comment lines
    COMMENT_PROB = 0.35 # probability of inline comments
    SEPS = [' ', '\t'] # separators randomly added to the coords file
    MAX_SEPS = 5       # maximum number of consecutive separators

    @classmethod
    def get_seps(cls, minimum):
        """ Return a string containing a random number of separators.

        The separators are randomly chosen from cls.SEPS. The returned string
        contains N of them, where N is a random integer such that: minimum <=
        N <= cls.MAX_SEPS.

        """

        n = random.randint(minimum, cls.MAX_SEPS)
        return ''.join(random.choice(cls.SEPS) for _ in range(n))

    @classmethod
    def get_comment(cls):
        """ Return a random string that starts with '#'. """

        # Use the name of one of the SIMBAD objects
        object_ = random.choice(cls.COORDINATES.keys())
        sep1 = cls.get_seps(0)
        sep2 = cls.get_seps(0)
        return "#" + sep1 + object_ + sep2

    @classmethod
    def get_coords_data(cls, coords):
        """ Format 'coords' as the contents of a coordinates file.

        Return a string that contains a line for each astronomical object in
        'coords', an iterable argument, listing their right ascensions and
        declinations in two columns and, if available, their proper motions
        in two additional columns, surrounded by brackets. For example:

          269.466450 4.705625 [0.0036] [-.0064]

        These columns and brackets are surrounded by a random number (up to the
        value of the MAX_SEPS class attribute) or random separators (SEPS class
        attribute). The returned string, after written to disk, is expected to
        be successfully parsed by methods.load_coordinates().

        """

        lines = []
        for ra, dec, pm_ra, pm_dec in coords:

            sep0 = cls.get_seps(0)
            sep1 = cls.get_seps(1)
            sep2 = cls.get_seps(0)
            line = "%s%.8f%s%.8f%s" % (sep0, ra, sep1, dec, sep2)

            if None not in (pm_ra, pm_dec):

                sep3 = cls.get_seps(0)
                sep4 = cls.get_seps(0)
                pm_ra_column = "[%s%.6f%s]" % (sep3, pm_ra, sep4)

                sep5 = cls.get_seps(0)
                sep6 = cls.get_seps(0)
                pm_dec_column = "[%s%.6f%s]" % (sep5, pm_dec, sep6)

                sep7 = cls.get_seps(1)
                sep8 = cls.get_seps(1)
                sep9 = cls.get_seps(0)
                line += sep7 + pm_ra_column + sep8 + pm_dec_column + sep9

            lines.append(line)

        return '\n'.join(lines)

    def test_load_coordinates(self):

        for _ in xrange(NITERS):

            # Randomly choose some of the SIMBAD astronomical objects and write
            # them to a temporary file, formatting their coordinates and proper
            # motions in four columns (or just two, if the proper motions are
            # not known) and inserting a random number of separators before,
            # between and after the columns and brackets. Then make sure that
            # methods.load_coordinates() returns the same astronomical objects,
            # in the same order and with the same coordinates that we wrote.

            n = random.randint(*self.NCOORDS)
            objects = random.sample(self.COORDINATES.values(), n)
            data = self.get_coords_data(objects)
            with methods.tempinput(data) as path:
                coordinates = methods.load_coordinates(path)
                for coords, expected in zip(coordinates, objects):
                    self.assertEqual(coords, expected)

    def test_load_coordinates_empty_lines_and_comments(self):

        # The same as test_load_coordinates(), but randomly inserting a few
        # empty and comment lines, as well as inline comments, all of which
        # must be ignored by load_coordinates().

        for _ in xrange(NITERS):
            n = random.randint(*self.NCOORDS)
            objects = random.sample(self.COORDINATES.values(), n)
            data = self.get_coords_data(objects)

            lines = data.split('\n')

            # Randomly insert inline comments
            for index in range(len(lines)):
                if random.random() < self.COMMENT_PROB:
                    sep = self.get_seps(0)
                    comment = self.get_comment()
                    lines[index] += sep + comment

            # Randomly insert empty lines
            for _ in range(random.randint(*self.NEMPTY)):
                index = random.randint(0, len(lines))
                empty = self.get_seps(0)
                lines.insert(index, empty)

            # Randomly insert comment lines
            for _ in range(random.randint(*self.NCOMMENTS)):
                index = random.randint(0, len(lines))
                sep = self.get_seps(0)
                comment = self.get_comment()
                lines.insert(index, sep + comment)

            data = '\n'.join(lines)

            with methods.tempinput(data) as path:
                coordinates = methods.load_coordinates(path)
                for coords, expected in zip(coordinates, objects):
                    self.assertEqual(coords, expected)

    def test_load_coordinates_empty_file(self):
        # If the file is empty, nothing is returned
        with methods.tempinput('') as path:
            self.assertEqual([], list(methods.load_coordinates(path)))

    def test_load_coordinates_invalid_data(self):

        def check_raise(data, exception, regexp):
            """ Make sure that methods.load_coordinates() raises 'exception'
            when a file containing 'data' is parsed. 'regexp' is the regular
            expression that must be matched by the string representation of
            the raised exception"""

            with methods.tempinput(data) as path:
                with self.assertRaisesRegexp(exception, regexp):
                    list(methods.load_coordinates(path))

        def get_coords():
            """ Return an element from COORDINATES with known proper motions """

            coords = []
            for c in self.COORDINATES.itervalues():
                if None not in (c.pm_ra, c.pm_dec):
                    coords.append(c)
            return random.choice(coords)

        # (1) Lines with other than (a) two floating-point numbers (right
        # ascension and declination) or (b) four floating-point numbers (alpha,
        # delta and proper motions, the last two surrounded by brackets).

        c = get_coords()
        unparseable_data = [

            # The names of three objects, no coordinates
            '\n'.join(["NGC 4494", "11 Com b", "TrES-1"]),
            # String + float
            "foo %.8f" % c.dec,
            # Three floating-point numbers
            ("%.8f " * 3) % c[:3],
            # Proper motions not in brackets
            ("%.8f " * 4) % c,
            # Missing declination proper motion
            (("%.8f " * 2) + "[%.6f]") % c[:-1],
            # Three proper motions
            (("%.8f " * 2) + ("[%.6f] " * 3)) % (c + (c.pm_dec,))]

        regexp = "Unable to parse line"
        for data in unparseable_data:
            check_raise(data, ValueError, regexp)

        # (2) An object with right ascension out of range
        c = get_coords()
        regexp = "Right ascension .* not in range"
        fmt = "%.8f %.8f [%.6f] [%.6f]"
        c = c._replace(ra = -24.19933)
        data1 = fmt % c  # RA < 0
        check_raise(data1, ValueError, regexp)
        data2 = "%.8f %.8f" % (360, c.dec)     # RA >= 360
        check_raise(data2, ValueError, regexp)
        data3 = "%.8f %.8f" % (417.993, c.dec) # RA >= 360
        check_raise(data3, ValueError, regexp)

        # (3) An object with declination out of range
        c = get_coords()
        regexp = "Declination .* not in range"
        c = c._replace(dec = -90.21)
        data1 = fmt % c # DEC < -90
        check_raise(data1, ValueError, regexp)
        data2 = "%.8f %.8f" % (c.ra,  113.93) # DEC > +90
        check_raise(data2, ValueError, regexp)

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


class MethodsQueueTest(unittest.TestCase):

    def shortDescription(self):
        """Don't use the first line of the test method's docstring."""
        pass

    def test_put_and_get(self):
        """Test that this is a FIFO queue."""

        q = methods.Queue()
        q.put(5)
        q.put(7)
        q.put(9)

        self.assertEqual(5, q.get())
        self.assertEqual(7, q.get())
        self.assertEqual(9, q.get())

    def test_put_when_queue_is_full(self):
        """Test that if Queue.Full is raised, qsize() doesn't break."""

        q = methods.Queue(maxsize=1)
        q.put(5)
        self.assertEqual(1, q.qsize())
        with self.assertRaises(Queue.Full):
            q.put(7, block=False)

        # Internal counter doesn't change, see:
        # https://github.com/vterron/lemon/pull/103
        self.assertEqual(1, q.qsize())

        # Pause to allow queue to complete operation before terminating
        # Otherwise sometimes results in 'IOError: [Errno 32] Broken pipe'
        time.sleep(0.1)

    def test_get_when_queue_is_empty(self):
        """Test that if Queue.Empty is raised, qsize() doesn't break."""

        q = methods.Queue()
        self.assertEqual(0, q.qsize())
        with self.assertRaises(Queue.Empty):
            q.get(block=False)

        # Internal counter doesn't change, see:
        # https://github.com/vterron/lemon/pull/103
        self.assertEqual(0, q.qsize())
        time.sleep(0.1)

    def test_qsize(self):
        """Test that qsize() gets updated correctly as we put() and get()."""

        q = methods.Queue()
        self.assertEqual(0, q.qsize())
        q.put(5)
        self.assertEqual(1, q.qsize())
        q.put(7)
        self.assertEqual(2, q.qsize())
        q.put(9)
        self.assertEqual(3, q.qsize())

        q.get()
        self.assertEqual(2, q.qsize())
        q.get()
        self.assertEqual(1, q.qsize())
        q.get()
        self.assertEqual(0, q.qsize())

    def test_empty(self):
        """Test empty()."""

        q = methods.Queue()
        self.assertTrue(q.empty())
        q.put(1)
        self.assertFalse(q.empty())
        time.sleep(0.1)

    def test_clear(self):
        """Test that clear() empties the queue."""

        q = methods.Queue()
        q.put(5)
        q.put(7)
        q.put(9)

        q.clear()
        self.assertEqual(0, q.qsize())
