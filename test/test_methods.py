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

    SEPS = [' ', '\t'] # separators randomly added to the coords file

    @classmethod
    def get_seps(cls, n):
        """ Return a string containing 'n' random separators """
        return ''.join(random.choice(cls.SEPS) for _ in range(n))

