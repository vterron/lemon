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
import numpy
import random
import unittest

from astromatic import Pixel, Star

NITERS = 100

class PixelTest(unittest.TestCase):

    X_COORD_RANGE = (1, 2048)
    Y_COORD_RANGE = (1, 2048)

    @classmethod
    def rargs(cls):
        """ Return the arguments needed to instantiate a random Pixel """
        x = random.uniform(*cls.X_COORD_RANGE)
        y = random.uniform(*cls.Y_COORD_RANGE)
        return x, y

    @classmethod
    def random(cls):
        """ Return a random Pixel object """
        x, y = cls.rargs()
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
            x, y = self.rargs()
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

    def test_distance(self):

        pixel1 = Pixel(4, 6)
        pixel2 = Pixel(8, 5)
        distance = pixel1.distance(pixel2)
        # ((4 - 8) ^ 2 + (6 - 5) ^ 2) ^ 0.5
        self.assertAlmostEqual(distance, 4.1231056256176606)

        pixel1 = Pixel(2.14, 1.89)
        pixel2 = Pixel(3.34, 8.01)
        distance = pixel1.distance(pixel2)
        # ((2.14 - 3.34) ^ 2 + (1.89 - 8.01) ^ 2) ^ 0.5
        self.assertAlmostEqual(distance, 6.2365375008894155)

        for _ in xrange(NITERS):

            pixel1 = self.random()
            pixel2 = self.random()
            distance = pixel1.distance(pixel2)

            # Compute the expected Euclidean distance with NumPy
            array1 = numpy.array([pixel1.x, pixel1.y])
            array2 = numpy.array([pixel2.x, pixel2.y])
            expected = numpy.linalg.norm(array1 - array2)
            self.assertAlmostEqual(distance, expected)


class StarTest(unittest.TestCase):

    X_COORD_RANGE = (1, 2048)
    Y_COORD_RANGE = (1, 2048)
    RIGHT_ASCENSION_RANGE = (0, 360)
    DECLINATION_RANGE = (-90, 90)
    ISOAREAF_RANGE = (1, 100)
    MAGNITUDE_RANGE = (1.47, 25)
    SNR_RANGE = (2, 10000)
    FWHM_RANGE = (0.5, 4.5)
    ELONGATION_RANGE = (1, 3.5)

    @classmethod
    def rargs(cls):
        """ Return the arguments needed to instantiate a random Star """

        x = random.uniform(*cls.X_COORD_RANGE)
        y = random.uniform(*cls.Y_COORD_RANGE)
        alpha = random.uniform(*cls.RIGHT_ASCENSION_RANGE)
        delta = random.uniform(*cls.DECLINATION_RANGE)
        isoareaf = random.uniform(*cls.ISOAREAF_RANGE)
        mag = random.uniform(*cls.MAGNITUDE_RANGE)
        saturated = random.choice([True, False])
        snr = random.uniform(*cls.SNR_RANGE)
        fwhm = random.uniform(*cls.FWHM_RANGE)
        elong = random.uniform(*cls.ELONGATION_RANGE)
        return x, y, alpha, delta, isoareaf, mag, saturated, snr, fwhm, elong

    @classmethod
    def random(cls):
        """ Return a random Star object """
        return Star(*cls.rargs())

    def test_init(self):
        for _ in xrange(NITERS):
            args = self.rargs()
            star = Star(*args)

            x, y, ra, dec, area, mag, satur, snr, fwhm, elong = args
            self.assertEqual(star.x, x)
            self.assertEqual(star.y, y)
            self.assertEqual(star.alpha, ra)
            self.assertEqual(star.delta, dec)
            self.assertEqual(star.area, area)
            self.assertEqual(star.mag, mag)
            self.assertEqual(star.saturated, satur)
            self.assertEqual(star.snr, snr)
            self.assertEqual(star.fwhm, fwhm)
            self.assertEqual(star.elongation, elong)

    def test_angular_distance(self):

        # The first case is taken from Sten Odenwald's Ask the Astronomer
        # [http://www.astronomycafe.net/qadir/q1890.html]. Except for the right
        # ascension and declination, which we set manually, the two Stars are
        # instantiated with random data.

        args1 = list(self.rargs())
        args1[2:4] = 100.2, -16.58
        star1 = Star(*args1)

        args2 = list(self.rargs())
        args2[2:4] = 87.5, 7.38
        star2 = Star(*args2)

        sky_distance = star1.angular_distance(star2)
        self.assertAlmostEqual(sky_distance, 27.054384870767787)

        # The second case is an example by David Oesper, taken from The Sky
        # This Week [http://www.skythisweek.info/angsep.pdf]. It computes the
        # angular distance between Merak and Dubhe, in the Big Dipper.

        args3 = list(self.rargs())
        args3[2:4] = 165.458, 56.3825
        star3 = Star(*args3)

        args4 = list(self.rargs())
        args4[2:4] = 165.933, 61.7511
        star4 = Star(*args4)

        sky_distance = star3.angular_distance(star4)
        self.assertAlmostEqual(sky_distance, 5.374111607543190)

