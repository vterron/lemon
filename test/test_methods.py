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

import random

# LEMON modules
from test import unittest
from astromatic import Coordinates
import methods

NITERS = 100  # How many times each test case is run with random data

class LoadCoordinatesTest(unittest.TestCase):

    # A series of two-element tuples, one per line. The first element is a
    # string containing the name of an astronomical object. The second is
    # another two-element tuple with its sky coordinates (right ascension
    # and declination), in decimal degrees.
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
    NEMPTY = (1, 50) # number of empty lines in each coordinates file
    SEPS = [' ', '\t'] # separators randomly added to the coords file
    MAX_SEPS = 5       # maximum number of consecutive separators

    @classmethod
    def get_seps(cls, n):
        """ Return a string containing 'n' random separators """
        return ''.join(random.choice(cls.SEPS) for _ in range(n))

    @classmethod
    def get_coords_data(cls, coords):
        """ Format 'coords' as the contents of a coordinates file.

        Return a string that contains a line for each astronomical object in
        'coords', an iterable argument, listing their right ascensions and
        declinations in two columns and decimal degrees. These columns are
        surrounded by a random number (up to the value of the MAX_SEPS class
        attribute) or random separators (SEPS class attribute). The returned
        string should, after written to disk, is expected to be successfully
        parsed by methods.load_coordinates().

        """

        lines = []
        for ra, dec in coords:
            sep1 = cls.get_seps(random.randint(0, cls.MAX_SEPS))
            sep2 = cls.get_seps(random.randint(1, cls.MAX_SEPS))
            sep3 = cls.get_seps(random.randint(0, cls.MAX_SEPS))
            line = "%s%.8f%s%.8f%s" % (sep1, ra, sep2, dec, sep3)
            lines.append(line)
        return '\n'.join(lines)

    def test_load_coordinates(self):

        for _ in xrange(NITERS):

            # Randomly choose some of the SIMBAD astronomical objects and write
            # them to a temporary file, formatting their right ascensions and
            # coordinates in two columns and inserting a random number of
            # separators before, between and after the columns. Make sure that
            # methods.load_coordinates() returns the same astronomical objects,
            # in the same order and with the same coordinates that we wrote.

            n = random.randint(*self.NCOORDS)
            objects = random.sample(self.COORDINATES.values(), n)
            data = self.get_coords_data(objects)
            with methods.tempinput(data) as path:
                coordinates = methods.load_coordinates(path)
                for index, (ra, dec) in enumerate(coordinates):
                    self.assertAlmostEqual(ra,  objects[index].ra)
                    self.assertAlmostEqual(dec, objects[index].dec)

    def test_load_coordinates_empty_lines(self):

        # The same as test_load_coordinates(), but randomly inserting a
        # few empty lines. These must be ignored by load_coordinates()

        for _ in xrange(NITERS):
            n = random.randint(*self.NCOORDS)
            objects = random.sample(self.COORDINATES.values(), n)
            data = self.get_coords_data(objects)

            # Randomly insert empty lines
            for _ in range(*self.NEMPTY):
                lines = data.split('\n')
                index = random.randint(0, len(lines))
                empty = self.get_seps(random.randint(0, self.MAX_SEPS))
                lines.insert(index, empty)

            data = '\n'.join(lines)
            # Empty lines must be ignored
            with methods.tempinput(data) as path:
                coordinates = methods.load_coordinates(path)
                for index, (ra, dec) in enumerate(coordinates):
                    self.assertAlmostEqual(ra,  objects[index].ra)
                    self.assertAlmostEqual(dec, objects[index].dec)

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

        # Lines with other than two floating-point numbers
        data1 = '\n'.join(["NGC 4494", "11 Com b", "TrES-1"])
        regexp = "Unable to parse line"
        check_raise(data1, ValueError, regexp)
        ra, dec = random.choice(self.COORDINATES.values())
        data2 = "foo %.8f" % dec
        check_raise(data1, ValueError, regexp)

        # An object with right ascension out of range
        ra, dec = random.choice(self.COORDINATES.values())
        regexp = "Right ascension .* not in range"
        data1 = "%.8f %.8f" % (-24.19933, dec) # RA < 0
        check_raise(data1, ValueError, regexp)
        data2 = "%.8f %.8f" % (360, dec)     # RA >= 360
        check_raise(data2, ValueError, regexp)
        data3 = "%.8f %.8f" % (417.993, dec) # RA >= 360
        check_raise(data3, ValueError, regexp)

        # An object with declination out of range
        ra, dec = random.choice(self.COORDINATES.values())
        regexp = "Declination .* not in range"
        data1 = "%.8f %.8f" % (ra,  -90.21) # DEC < -90
        check_raise(data1, ValueError, regexp)
        data2 = "%.8f %.8f" % (ra,  113.93) # DEC > +90
        check_raise(data2, ValueError, regexp)
