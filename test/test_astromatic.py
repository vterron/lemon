#! /usr/bin/env python

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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import copy
import random
import unittest

from astromatic import Pixel

NITERS = 100

class PixelTest(unittest.TestCase):

    X_COORD_RANGE = (1, 2048)
    Y_COORD_RANGE = (1, 2048)

    @classmethod
    def random(cls):
        """ Return a random Pixel object """
        x = random.uniform(*cls.X_COORD_RANGE)
        y = random.uniform(*cls.Y_COORD_RANGE)
        return Pixel(x, y)

    @classmethod
    def different(cls, pixel):
        """ Return a Pixel with a different x- or y-coordinate value, or both"""
        while True:
            another = cls.random()
            assert id(pixel) != id(another)
            if pixel != another:
                return another

    def test_init(self):
        for _ in xrange(NITERS):
            x = random.uniform(*self.X_COORD_RANGE)
            y = random.uniform(*self.Y_COORD_RANGE)
            pixel = Pixel(x, y)
            self.assertEqual(pixel.x, x)
            self.assertEqual(pixel.y, y)

    def test_repr(self):
        for _ in xrange(NITERS):
            pixel = self.random()
            repr_pixel = eval(`pixel`)

            # We need to use TestCase.assertAlmostEqual, instead of a simple
            # equality comparison, because the precision of the coordinates in
            # 'repr_pixel' is limited by the number of decimal places printed
            # by __repr__. The Pixel returned by eval(`pixel`), therefore, may
            # not be exactly equal to 'pixel' when the coordinates are real
            # numbers, but we just want to verify that it __repr__ computes a
            # valid, approximate-enough string representation of the object.

            kwargs = dict(places = 5)
            self.assertAlmostEqual(pixel.x, repr_pixel.x, **kwargs)
            self.assertAlmostEqual(pixel.y, repr_pixel.y, **kwargs)

    def test_eq(self):

        pix1 = Pixel(13.2, 15.4)
        pix2 = Pixel(67.3, 12.5)
        self.assertEqual(pix1, pix1)
        self.assertEqual(pix2, pix2)
        self.assertNotEqual(pix1, pix2)

        for _ in xrange(NITERS):

            # Three objects: 'pixel' == 'identical' != 'different'
            pixel = self.random()
            identical = copy.deepcopy(pixel)

            self.assertEqual(pixel, pixel)
            self.assertEqual(identical, identical)
            self.assertEqual(pixel, identical)

            different = self.different(pixel)
            self.assertEqual(different, different)
            self.assertNotEqual(pixel, different)
            self.assertNotEqual(identical, different)

    def test_hash(self):

        pix1 = Pixel(83.1, 18.2)
        pix2 = Pixel(33.1, 13.4)
        self.assertEqual(hash(pix1), hash(pix1))
        self.assertEqual(hash(pix2), hash(pix2))
        self.assertNotEqual(hash(pix1), hash(pix2))

        for _ in xrange(NITERS):

            pixel = self.random()
            self.assertEqual(hash(pixel), hash(pixel))

            # Half of these random cases use 'pixel' and another Pixel to which
            # it compares equal (so they must have the same hash value), where
            # the other half uses a different Pixel (the hashes must differ).

            if random.choice([True, False]):
                identical = copy.deepcopy(pixel)
                self.assertEqual(hash(identical), hash(identical))
                self.assertEqual(hash(pixel), hash(identical))
            else:
                different = self.different(pixel)
                self.assertEqual(hash(different), hash(different))
                self.assertNotEqual(hash(pixel), hash(different))

