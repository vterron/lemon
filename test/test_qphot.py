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

import astropy.wcs
import itertools
import operator
import os.path
import pyfits
import random
import tempfile

# LEMON modules
from test import unittest
import astromatic
import fitsimage
import methods
import qphot
import test.test_fitsimage

NITERS = 100  # How many times random-data tests cases are run

def fix_DSS_image(path):
    """ Fix images downloaded from STScI DSS before we can do photometry.

    The images that we download from the STScI Digitized Sky Survey (DSS) to
    use in our unit tests are multi-extension FITS files. This causes IRAF's
    qphot to complain that 'FXF: must specify which FITS extension' and exit.
    Also, the exposure time is stored in the header in minutes, instead of in
    seconds as our code expects.

    This function accepts as input a DSS image and returns the path to a
    temporary copy of it with only the primary HDU, which is the one we are
    interested in; and the exposure time in seconds. You are responsible for
    deleting this temporary FITS file when done with it.

    """

    kwargs = dict(suffix = os.path.splitext(path)[-1])
    with tempfile.NamedTemporaryFile(**kwargs) as fd:
        output_path = fd.name

    with pyfits.open(path) as hdu:
        # Keep only the primary header
        del hdu[1:]
        header = hdu[0].header
        # From minutes to seconds
        header['EXPOSURE'] *= 60
        hdu.writeto(output_path)

    return output_path


class QPhotTest(unittest.TestCase):

    QPHOT_KWARGS = dict(
        epoch = 2000,
        aperture = 11,
        annulus  = 13,
        dannulus = 8,
        maximum = 30000,
        datek = 'DATE-OBS',
        timek = None,
        exptimek = 'EXPOSURE',
        uncimgk = None)

    def test_get_coords_file(self):

        def c(*args):
            """ Return an astromatic.Coordinates object """
            return astromatic.Coordinates(*args)

        coordinates = (
            c(316.724802,  38.74944,   4.16831,   3.2692),    # 61 Cygni A
            c(316.730266,  38.742056,  4.1069,    3.14468),   # 61 Cygni B
            c(165.834142,  35.96988,  -0.58027,  -4.76585),   # Lalande 21185
            c(152.11717,   12.3065,   -0.000114, -0.000126),  # Leo I
            c(346.466817, -35.853069,  6.7682,    1.32752),   # HD 217987
            c(348.992913,  31.462856,  None,      None),      # WASP-10 b
            c( 97.19046,   38.962963)                         # HD 45350 b
            )

        # There is no need to write code to test that proper-motion correction
        # is correctly calculated: we do already have a unit test for exactly
        # that: Coordinates.get_exact_coordinates(). Parse the output file of
        # get_coords_file() and make sure that the coordinates written to it
        # match those returned by get_exact_coordinates() for the same year
        # and epoch.

        for _ in xrange(NITERS):

            year  = random.uniform(1900, 2050)
            epoch = random.choice([1950, 2000])

            # The proper-motion corrected astromatic.Coordinates objects. We
            # only correct those objects with proper motions (that is, those
            # for which both 'pm_ra' and 'pm_dec' are not None neither zero).

            func = operator.methodcaller('get_exact_coordinates', year, epoch)
            has_pm = lambda c: c.pm_ra and c.pm_dec
            expected = (func(x) if has_pm(x) else x for x in coordinates)

            output_path = qphot.get_coords_file(coordinates, year, epoch)

            try:
                # The coordinates written by get_coords_file(), returned as
                # four-element tuples (ra, dec, pm_ra, pm_dec)
                output = methods.load_coordinates(output_path)
                for c1, c2 in zip(output, expected):
                    self.assertAlmostEqual(c1[0], c2.ra)
                    self.assertAlmostEqual(c1[1], c2.dec)
            finally:
                os.unlink(output_path)

    def test_qphot_run(self):

        # A simple test: do photometry on the DSS image of NGC 2264, measuring
        # ten astronomical objects whose celestial coordinates (right ascension
        # and declination) we know. Then, compare the photometric measurements,
        # returned as qphot.QPhotResult objects, to the expected values. This
        # is a simple case in which the objects do not have proper motions: we
        # are mostly testing here that IRAF's qphot runs without any problems.

        ngc2264_path = './test/test_data/fits/NGC_2264.fits'
        ngc2264_input_coords = (
            astromatic.Coordinates(100.1543316, 9.7909363),
            astromatic.Coordinates(100.1597762, 9.7878795),
            astromatic.Coordinates(100.2147546, 9.8636567),
            astromatic.Coordinates(100.2502955, 9.8714701),
            astromatic.Coordinates(100.2933265, 9.8838196),
            astromatic.Coordinates(100.1191901, 9.8177770),
            astromatic.Coordinates(100.1598790, 9.9627296),
            astromatic.Coordinates(100.2446191, 9.8962391),
            astromatic.Coordinates(100.2579343, 9.8802548),
            astromatic.Coordinates(100.3635468, 9.8540181))

        ngc2264_expected_output = [
            #                 x        y        mag     sum      flux     stdev
            qphot.QPhotResult(755.241, 75.308,  17.821, 6015410, 3124231, 579.0784),
            qphot.QPhotResult(736.138, 64.345,  17.777, 6035459, 3254578, 572.4526),
            qphot.QPhotResult(542.399, 334.835, 18.01,  5197306, 2624735, 374.3998),
            qphot.QPhotResult(417.419, 362.544, 18.305, 4541155, 2001931, 362.4446),
            qphot.QPhotResult(266.116, 406.437, 18.12,  4804557, 2372792, 378.5256)]

        # IRAF returns different values depending on the architecture (32 vs 64
        # bits). This is a very strange and confusing behavior, since most of
        # the time only the *flux* is different, while the magnitude, total
        # number of counts and even the standard deviation are the same. The
        # difference looks too big and specific to be only caused by the
        # limitations of floating-point arithmetic. Until the transition to
        # Astropy's photutils, take into account both cases (32- vs 64-bit
        # architecture). See issue #56 for futher information:
        # https://github.com/vterron/lemon/issues/57

        nbits = methods.get_nbits()
        if nbits == 64:
            ngc2264_expected_output += [                             # <<>>
                qphot.QPhotResult(878.62,  171.496, 17.623, 6266626, 3749740, 484.9941),
                qphot.QPhotResult(734.601, 689.319, 17.631, 6102620, 3723123, 555.338) ,
                qphot.QPhotResult(437.227, 451.111, 18.54,  9239881, 1611601, 3583.002),
                qphot.QPhotResult(390.527, 393.897, 17.768, 5925938, 3280150, 450.3236),
                qphot.QPhotResult(19.449,  299.555, 17.409, 7143799, 4567058, 540.6945)]
        else:
            assert nbits == 32
            ngc2264_expected_output += [                             # <<>>
                qphot.QPhotResult(878.62,  171.496, 17.623, 6266626, 3750371, 484.9941),
                qphot.QPhotResult(734.601, 689.319, 17.631, 6102620, 3723439, 555.338) ,
                qphot.QPhotResult(437.227, 451.111, 18.54,  9239881, 1611932, 3583.002),
                #                                             <<>>
                qphot.QPhotResult(390.527, 393.897, 17.768, 5925937, 3280362, 450.3236),
                qphot.QPhotResult(19.449,  299.555, 17.409, 7143799, 4567332, 540.6945)]

        path = fix_DSS_image(ngc2264_path)
        with test.test_fitsimage.FITSImage(path) as img:
            result = qphot.run(img, ngc2264_input_coords, **self.QPHOT_KWARGS)
            for phot, expected_phot in zip(result, ngc2264_expected_output):
                self.assertEqual(phot, expected_phot)

        # When none of the coordinates have a known proper motion, KeyError is
        # *not* raised if 'datek' / 'timek' cannot be found in the FITS header:
        # since there are no proper motions to correct, qphot.run() does not
        # even bother reading the keywords from the header.

        kwargs = self.QPHOT_KWARGS.copy()
        kwargs['datek'] = 'MISSING-KWD'
        kwargs['timek'] = 'MISSING-KWD'

        path = fix_DSS_image(ngc2264_path)
        with test.test_fitsimage.FITSImage(path) as img:

            # Make sure that both keywords are missing from the FITS header
            with self.assertRaises(KeyError):
                img.read_keyword(kwargs['datek'])
            with self.assertRaises(KeyError):
                img.read_keyword(kwargs['timek'])

            # No exception raised, although both keywords are missing
            qphot.run(img, ngc2264_input_coords, **kwargs)

    def test_qphot_run_nonzero_cbox(self):

        # A call to run() with a 'cbox' other than zero. Although the four
        # input coordinates are the same as in the previous test case, the
        # output image coordinates are slightly different, as qphot has
        # computed the accurate center for each object using the centroid
        # centering algorithm.

        ngc2264_path = './test/test_data/fits/NGC_2264.fits'
        ngc2264_input_coords = (
            astromatic.Coordinates(100.1191901, 9.8177770),
            astromatic.Coordinates(100.1543316, 9.7909363),
            astromatic.Coordinates(100.1597762, 9.7878795),
            astromatic.Coordinates(100.1598790, 9.9627296))

        ngc2264_expected_output = (
            #                 x        y        mag     sum      flux     stdev
            qphot.QPhotResult(877.992, 171.373, 17.627, 6248653, 3737111, 484.8067),
            qphot.QPhotResult(752.713, 76.071,  17.904, 5820262, 2894284, 548.6724),
            qphot.QPhotResult(735.552, 65.18,   17.785, 6019633, 3231141, 575.0005),
            qphot.QPhotResult(734.0,   689.025, 17.639, 6075931, 3694561, 567.6364))

        kwargs = self.QPHOT_KWARGS.copy()
        kwargs['cbox'] = 5

        path = fix_DSS_image(ngc2264_path)
        with test.test_fitsimage.FITSImage(path) as img:
            result = qphot.run(img, ngc2264_input_coords, **kwargs)
            for phot, expected_phot in zip(result, ngc2264_expected_output):
                self.assertEqual(phot, expected_phot)

    def test_qphot_run_proper_motions(self):

        # Do photometry on Barnard's Star, the star with the largest-known
        # proper motion relative to the Solar System. Its coordinates are
        # J2000, but the DSS image dates back from 1993: qphot.run() must
        # correct the right ascension and declination before doing photometry,
        # adjusting for the change in its position over those seven years.

        barnard_path = "./test/test_data/fits/Barnard's_Star.fits"
        barnard = astromatic.Coordinates(269.452075, 4.693391, -0.79858, 10.32812)
        expected_output = (   #  x        y       mag     sum     flux      stdev
            qphot.QPhotResult(440.947, 382.595, 17.245, 8563646, 5311413, 1039.844))

        path = fix_DSS_image(barnard_path)
        with test.test_fitsimage.FITSImage(path) as img:

            # Fix incorrect date in FITS header: '1993-07-26T04:87:00'.
            # Subtract sixty minutes and add one hour: 4h87m == 5h27m
            keyword = 'DATE-OBS'
            assert img.read_keyword(keyword) == '1993-07-26T04:87:00'
            img.update_keyword(keyword, '1993-07-26T05:27:00')

            # The proper-motion corrected coordinates
            year = img.year(exp_keyword = 'EXPOSURE')
            expected_coordinates = barnard.get_exact_coordinates(year)

            result = qphot.run(img, [barnard], **self.QPHOT_KWARGS)[0]
            self.assertEqual(result, expected_output)

            # Transform the pixel coordinates that IRAF's qphot outputs for
            # each measured object to celestial coordinates. This allows us to
            # make sure that photometry has been effectively done on the right,
            # proper-motion corrected coordinates.

            wcs = astropy.wcs.WCS(img._header)
            ra, dec = wcs.all_pix2world(result.x, result.y, 1)
            f = self.assertAlmostEqual
            f(ra,  expected_coordinates.ra,  delta = 1e-4) # delta = 0.024 arcsec
            f(dec, expected_coordinates.dec, delta = 1e-7) # delta = 0.00036 arsec
