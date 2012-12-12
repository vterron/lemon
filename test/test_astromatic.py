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
import functools
import numpy
import os
import os.path
import random
import tempfile
import unittest

from astromatic import Pixel, Star, Catalog

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

    def test_immutability(self):
        """ Make sure that the Pixel class is immutable """

        pixel = self.random()
        for name in pixel._asdict().iterkeys():
            value = random.random()
            self.assertRaises(AttributeError, setattr, pixel, name, value)
            self.assertRaises(AttributeError, delattr, pixel, name)

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

            # Make sure the 'x', 'y', 'alpha' and 'delta' properties work
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

    def test_distance(self):

        args1 = list(self.rargs())
        args1[0:2] = 134, 345
        star1 = Star(*args1)

        args2 = list(self.rargs())
        args2[0:2] = 178, 471
        star2 = Star(*args2)

        # ((134 - 178) ^ 2 + (345 - 471) ^ 2) ^ 0.5
        distance = star1.distance(star2)
        self.assertAlmostEqual(distance, 133.46160496562297)

        for _ in xrange(NITERS):

            star1 = self.random()
            star2 = self.random()
            distance = star1.distance(star2)

            # Compute the expected Euclidean distance with NumPy
            array1 = numpy.array([star1.x, star1.y])
            array2 = numpy.array([star2.x, star2.y])
            expected = numpy.linalg.norm(array1 - array2)
            self.assertAlmostEqual(distance, expected)


class CatalogTest(unittest.TestCase):

    args = os.path.dirname(__file__), './test_data/'
    get_path = functools.partial(os.path.join, *args)
    SAMPLE_CATALOG_PATH = get_path('sextractor.cat')
    SAMPLE_INCOMPLETE_PATH = get_path('sextractor_incomplete.cat')
    SAMPLE_NOASCIIHEAD_PATH = get_path('sextractor_noasciihead.cat')

    def test_flag_saturated(self):

        # Numbers in the [0, 255] range whose 3rd least significant bit is one
        saturated_flags = set([4, 5, 6, 7, 12, 13, 14, 15, 20, 21, 22, 23, 28,
        29, 30, 31, 36, 37, 38, 39, 44, 45, 46, 47, 52, 53, 54, 55, 60, 61, 62,
        63, 68, 69, 70, 71, 76, 77, 78, 79, 84, 85, 86, 87, 92, 93, 94, 95,
        100, 101, 102, 103, 108, 109, 110, 111, 116, 117, 118, 119, 124, 125,
        126, 127, 132, 133, 134, 135, 140, 141, 142, 143, 148, 149, 150, 151,
        156, 157, 158, 159, 164, 165, 166, 167, 172, 173, 174, 175, 180, 181,
        182, 183, 188, 189, 190, 191, 196, 197, 198, 199, 204, 205, 206, 207,
        212, 213, 214, 215, 220, 221, 222, 223, 228, 229, 230, 231, 236, 237,
        238, 239, 244, 245, 246, 247, 252, 253, 254, 255])

        for flag in xrange(0, 256):
            is_saturated = Catalog.flag_saturated(flag)
            if flag in saturated_flags:
                self.assertTrue(is_saturated)
            else:
                self.assertFalse(is_saturated)

        # Flags outside of the range raise ValueError
        self.assertRaises(ValueError, Catalog.flag_saturated, -2)
        self.assertRaises(ValueError, Catalog.flag_saturated, -1)
        self.assertRaises(ValueError, Catalog.flag_saturated, 256)
        self.assertRaises(ValueError, Catalog.flag_saturated, 257)

    def test_init(self):

        path = self.SAMPLE_CATALOG_PATH
        catalog = Catalog(path)
        self.assertEqual(catalog.path, path)
        self.assertEqual(len(catalog), 127)

        star = catalog[0]
        self.assertEqual(star.x, 844.359)
        self.assertEqual(star.y, 434.109)
        self.assertEqual(star.alpha, 100.2910553)
        self.assertEqual(star.delta, +9.4697032)
        self.assertEqual(star.area, 8085)
        self.assertEqual(star.mag, 8.2460)
        self.assertFalse(star.saturated)
        # S/N = FLUX_ISO / FLUXERR_ISO = 8545557 / 12811.09
        self.assertAlmostEqual(star.snr, 667.043709786)
        self.assertAlmostEqual(star.fwhm, 12.668) # FLUX_RADIUS x 2
        self.assertEqual(star.elongation, 1.562)

        star = catalog[68]
        self.assertEqual(star.x, 4.189)
        self.assertEqual(star.y, 1995.196)
        self.assertEqual(star.alpha, 100.1844554)
        self.assertEqual(star.delta, +9.2744660)
        self.assertEqual(star.area, 16)
        self.assertEqual(star.mag, 14.8506)
        self.assertFalse(star.saturated)
        # S/N = FLUX_ISO / FLUXERR_ISO = 10539.52 / 569.9098
        self.assertAlmostEqual(star.snr, 18.49331245)
        self.assertAlmostEqual(star.fwhm, 2.126) # FLUX_RADIUS x 2
        self.assertEqual(star.elongation, 1.339)

        star = catalog[126]
        self.assertEqual(star.x, 107.033)
        self.assertEqual(star.y, 1953.500)
        self.assertEqual(star.alpha, 100.1971009)
        self.assertEqual(star.delta, +9.2796855)
        self.assertEqual(star.area, 10)
        self.assertEqual(star.mag, 15.9567)
        self.assertTrue(star.saturated)
        # S/N = FLUX_ISO / FLUXERR_ISO = 3833.773 / 450.5532
        self.assertAlmostEqual(star.snr, 8.509035115)
        self.assertAlmostEqual(star.fwhm, 2.412) # FLUX_RADIUS x 2
        self.assertEqual(star.elongation, 1.613)

        # IOError is raised if the SExtractor catalog does not exist; in order
        # to get the path to a non-existent file, we use NamedTemporaryFile and
        # close the file (deleting it) immediately after.
        with tempfile.NamedTemporaryFile(suffix = '.cat') as fd: pass
        self.assertFalse(os.path.exists(fd.name))
        self.assertRaises(IOError, Catalog, fd.name)

        # ValueError is raised if one of the required SExtractor parameters is
        # not in the catalog, or if it was not saved in the ASCII_HEAD format
        # and therefore there are no comment lines listing column labels, which
        # are needed in order to determine in which column each parameter is.
        self.assertRaises(ValueError, Catalog, self.SAMPLE_INCOMPLETE_PATH)
        self.assertRaises(ValueError, Catalog, self.SAMPLE_NOASCIIHEAD_PATH)

    def test_get_image_coordinates(self):

        catalog = Catalog(self.SAMPLE_CATALOG_PATH)
        self.assertEqual(len(catalog), 127)
        pixels = catalog.get_image_coordinates()

        # The returned list must contain as many Pixel objects, and in the same
        # order, as detected sources are in the SExtractor catalog. That is to
        # say that the x- and y-coordinates of the first Pixel must match those
        # of the first object in the Catalog, and so on.

        self.assertEqual(len(pixels), len(catalog))
        for pixel, star in zip(pixels, catalog):
            self.assertEqual(pixel.x, star.x)
            self.assertEqual(pixel.y, star.y)

