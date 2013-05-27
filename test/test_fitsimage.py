#! /usr/bin/env python

# Copyright (c) 2013 Victor Terron. All rights reserved.
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

import numpy.random
import os
import pyfits
import random
import tempfile
import unittest

# LEMON modules
import fitsimage

class FITSTestImage(fitsimage.FITSImage):
    """ Delete FITSImage on exit from the body of the with statement.

    This class does not alter the functionality of the fitsimage.FITSImage
    class. The only difference is that it adds support for the with statement
    (the __enter__() and __exit__() magic methods), deleting the FITS file on
    exit from the body of the with statement. In this manner, we can focus on
    writing our unit tests without having to worry about deleting the temporary
    FITS files that we use.

    """

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.unlink()


class FITSImageTest(unittest.TestCase):

    MIN_PIXEL_VALUE = 0
    MAX_PIXEL_VALUE = 65535
    MIN_SIZE = 100
    MAX_SIZE = 2048

    @classmethod
    def mkfits(cls, x_size, y_size, **keywords):
        """ Return a two-dimensional FITS image with random data.

        Create a temporary FITS image of size 'x_size' times 'y_size' pixels
        and populate it with random values in the range [MIN_PIXEL_VALUE,
        MAX_PIXEL_VALUE]. Keyword/value pairs are stored in the header of the
        FITS image: e.g., passing object='Salvor Hardin' sets the value of the
        OBJECT keyword in the FITS header to this very value, 'Salvor Hardin'.
        The method returns the absolute pathname of the temporary FITS image.

        """

        kwargs = dict(low = cls.MIN_PIXEL_VALUE,
                      high = cls.MAX_PIXEL_VALUE,
                      size = (x_size, y_size))
        pixels = numpy.random.random_integers(**kwargs)
        hdu = pyfits.PrimaryHDU(pixels)

        for keyword, value in keywords.iteritems():
            if len(keyword) > 8:
                keyword = 'HIERARCH' + keyword
            hdu.header.update(keyword, value)

        fd, path = tempfile.mkstemp(suffix = '.fits')
        os.close(fd)
        hdu.writeto(path)
        return path

    @classmethod
    def random(cls, **keywords):
        """ Return a random, temporary FITS image.

        This method returns a temporary FITS image (see FITSImageTest.mkfits),
        encapsulated as a FITSTestImage object -- which means that the file is
        automatically deleted on exit from the body of a with statement. The
        dimensions of the FITS image along the x- and y- axes are random
        numbers in the range [MIN_SIZE, MAX_SIZE]. Keyword/value pairs are
        passed down to mkfits() and stored in the header of the FITS image.

        """

        x_size = random.randint(cls.MIN_SIZE, cls.MAX_SIZE)
        y_size = random.randint(cls.MIN_SIZE, cls.MAX_SIZE)
        path = cls.mkfits(x_size, y_size, **keywords)
        return FITSTestImage(path)

    def test_read_keyword(self):

        # Read a character string
        keyword = 'OBSERVER'
        observer = 'Winnie-the-Pooh'
        kwargs = {keyword : observer}
        with self.random(**kwargs) as img:
            self.assertEqual(img.read_keyword(keyword), observer)

        # Read a floating-point number
        keyword = 'AIRMASS'
        airmass = 1.231
        kwargs = {keyword : airmass}
        with self.random(**kwargs) as img:
            self.assertAlmostEqual(img.read_keyword(keyword), airmass)

        # Read a keyword longer than eight characters. These have to use the
        # special keyword HIERARCH internally, with the actual long keyword
        # following, but can be accessed by using the keyword name directly,
        # with or without the 'hierarch' prepended.

        keyword = "CAHA GEN AMBI WIND SPEED"
        wind_speed = 6.1
        kwargs = {keyword : wind_speed}
        with self.random(**kwargs) as img:
            # Try both ways, with and without 'HIERARCH'
            self.assertAlmostEqual(img.read_keyword(keyword), wind_speed)
            keyword = 'HIERARCH' + keyword
            self.assertAlmostEqual(img.read_keyword(keyword), wind_speed)

        # Keywords are case-insensitive
        keyword = 'FILTER'
        img_filter = 'Johnson R'
        kwargs = {keyword.lower() : img_filter}
        with self.random(**kwargs) as img:
            self.assertEqual(img.read_keyword(keyword), img_filter)

        # Read several keywords from the FITS header
        instrument = "3.5m CAHA"
        image_type = 'science'
        ccd_temperature = -115.4
        keywords = {'INSTRUMENT' : instrument,
                    'IMAGE_TYPE' : image_type,
                    'CCDTEMP' : ccd_temperature}
        with self.random(**keywords) as img:
            self.assertEqual(img.read_keyword('INSTRUMENT'), instrument)
            self.assertEqual(img.read_keyword('IMAGE_TYPE'), image_type)
            self.assertEqual(img.read_keyword('CCDTEMP'), ccd_temperature)

        # TypeError is raised if the value of the 'keyword' argument is None,
        # ValueError if it is an empty string and, finally, KeyError if the
        # keyword cannot be found in the FITS header.

        with self.random() as img:
            self.assertRaises(TypeError, img.read_keyword, None)
            self.assertRaises(ValueError, img.read_keyword, '')
            self.assertRaises(KeyError, img.read_keyword, 'EXPTIME')

