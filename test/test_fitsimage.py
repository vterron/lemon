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
import tempfile
import unittest

class FITSImageTest(unittest.TestCase):

    MIN_PIXEL_VALUE = 0
    MAX_PIXEL_VALUE = 65535

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

