#! /usr/bin/env python2

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

import absl.testing.parameterized as parameterized
import datetime
import calendar
import numpy.random
import os
import pyfits
import random
import stat
import tempfile
import warnings

# LEMON modules
from test import unittest
import fitsimage

NITERS = 10  # How many times random-data tests case are run

class FITSImage(fitsimage.FITSImage):
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
        os.unlink(self.path)


class FITSImageTest(parameterized.TestCase):

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
                      size = (y_size, x_size))
        pixels = numpy.random.random_integers(**kwargs)
        hdu = pyfits.PrimaryHDU(pixels)

        for keyword, value in keywords.iteritems():
            hdu.header[keyword] = value

        fd, path = tempfile.mkstemp(suffix = '.fits')
        os.close(fd)
        hdu.writeto(path)
        return path

    @classmethod
    def random_data(cls, **keywords):
        """ Return a random, temporary FITS image.

        Return a three-element tuple containing a random temporary FITS image
        (see the FITSImageTest.mkfits() class method) and its dimensions along
        the x- and y- axes -- which are random numbers in the range [MIN_SIZE,
        MAX_SIZE]. Keyword/value pairs are passed down to mkfits() and stored
        in the header of the FITS image.

        """

        # Ignore keywords mapping to None.
        for k, v in keywords.items():
            if v is None:
                keywords.pop(k)

        x_size = random.randint(cls.MIN_SIZE, cls.MAX_SIZE)
        y_size = random.randint(cls.MIN_SIZE, cls.MAX_SIZE)
        path = cls.mkfits(x_size, y_size, **keywords)
        return path, x_size, y_size

    @classmethod
    def random(cls, **keywords):
        """ Return a random FITSImage object.

        Return a temporary FITS image encapsulated as a FITSImage object,
        which means that the file is automatically deleted on exit from the
        body of a with statement. See the FITSImageTest.random_data() class
        method for further information on how the random image is created.

        """

        path = cls.random_data(**keywords)[0]
        return FITSImage(path)

    def test_init(self):
        for _ in xrange(NITERS):
            path, x_size, y_size = self.random_data()
            with FITSImage(path) as img:
                self.assertEqual(img.path, path)
                self.assertEqual(img.size, (x_size, y_size))

        # IOError raised if we do not have permission to open the file...
        with self.random() as img:
            nonreadable_path = img.path
            mode = stat.S_IMODE(os.stat(nonreadable_path)[stat.ST_MODE])
            mode ^= stat.S_IRUSR
            os.chmod(nonreadable_path, mode)
            with self.assertRaises(IOError):
                FITSImage(nonreadable_path)

        # ... or if it simply does not exist
        nonexistent_path = path
        self.assertFalse(os.path.exists(nonexistent_path))
        with self.assertRaises(IOError):
            FITSImage(nonexistent_path)

        # NonStandardFITS raised if the FITS file does not conform to the FITS
        # standard. When this is the case, the 'SIMPLE' keyword is expected to
        # be False (as, according to the Standard, "A [logical constant] of F
        # signifies that the file does not conform to this standard").
        nonstandard_path = self.random_data(SIMPLE = False)[0]
        with self.assertRaises(fitsimage.NonStandardFITS):
            FITSImage(nonstandard_path)
        os.unlink(nonstandard_path)

        # It may also happen that the 'SIMPLE' keyword does not exist
        nonstandard_path = self.random_data()[0]
        handler = pyfits.open(nonstandard_path, mode = 'update')
        del handler[0].header['SIMPLE']
        handler.close(output_verify = 'ignore')

        # Ignore PyFITS <= 3.2 warning: "Error validating header for HDU #0
        # (note: PyFITS uses zero-based indexing). Block does not begin with
        # SIMPLE or XTENSION. There may be extra bytes after the last HDU or
        # the file is corrupted".

        with warnings.catch_warnings():
            msg = "(?s)Error validating header for .+ file is corrupted."
            warnings.filterwarnings('ignore', message=msg)

            with self.assertRaises(fitsimage.NonStandardFITS):
                FITSImage(nonstandard_path)
            os.unlink(nonstandard_path)

        # NonStandardFITS exception must also be raised when we try to open
        # anything that is not a FITS file. Among the countless kinds of file
        # types that could be used for this, try to open (a) an empty file...
        with tempfile.NamedTemporaryFile(suffix = '.fits') as fd:
            empty_path = fd.name
            with self.assertRaises(fitsimage.NonStandardFITS):
                FITSImage(empty_path)

        # ... and (b) a non-empty text file.
        with tempfile.NamedTemporaryFile(suffix = '.fits') as fd:
            fd.write("Lorem ipsum dolor sit amet,\n")
            fd.write("consectetur adipisicing elit\n")
            fd.flush()
            text_path = fd.name
            with self.assertRaises(fitsimage.NonStandardFITS):
                FITSImage(text_path)

    def test_repr(self):
        with self.random() as img1:
            self.assertEqual(img1.path, eval(repr(img1)).path)

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
            self.assertAlmostEqual(img.read_keyword('HIERARCH ' + keyword), wind_speed)
            # KeyError raised if 'HIERARCH' does not include a whitespace
            with self.assertRaises(KeyError):
                img.read_keyword('HIERARCH' + keyword)

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
            with self.assertRaises(TypeError):
                img.read_keyword(None)
            with self.assertRaises(ValueError):
                img.read_keyword('')
            with self.assertRaises(KeyError):
                img.read_keyword('EXPTIME')

    def test_update_keyword(self):

        def get_comment(image, keyword):
            """ Return the comment associated with a FITS keyword """
            return image._header.comments[keyword]

        keywords = {'INSTRUMENT' : "2.2m CAHA", 'CCD_TEMP' : -129.12}
        with self.random(**keywords) as img:

            # The keyword is not in the FITS header, so add it.
            keyword = 'OBSERVER'
            observer = "Subrahmanyan Chandrasekhar"
            observer_comment = "1983 Nobel Prize for Physics"
            kwargs = dict(comment = observer_comment)
            img.update_keyword(keyword, observer, **kwargs)
            self.assertEqual(img.read_keyword(keyword), observer)
            self.assertEqual(get_comment(img, keyword), observer_comment)

            # Update a floating-point number
            keyword = 'CCD_TEMP'
            new_temperature = -89.133
            temperature_comment = "degrees on the Celsius scale"
            kwargs = dict(comment = temperature_comment)
            img.update_keyword(keyword, new_temperature, **kwargs)
            self.assertEqual(img.read_keyword(keyword), new_temperature)
            self.assertEqual(get_comment(img, keyword), temperature_comment)

            # Update a keyword longer than eight characters.
            keyword = 'INSTRUMENT'
            new_instrument = "PANIC (2.2m CAHA)"
            img.update_keyword(keyword, new_instrument)
            # Try both ways, with and without 'HIERARCH'.
            self.assertEqual(img.read_keyword(keyword), new_instrument)
            self.assertEqual(img.read_keyword('HIERARCH ' + keyword), new_instrument)
            self.assertEqual(get_comment(img, keyword), '') # no comment

            # Keywords are case-insensitive. Also, note that only the name of
            # the observer is updated here: we are not specifying the comment
            # associated with the keyword. Thus, it retains its value ("1983
            # Nobel Prize for Physics", the one we used for Chandrasekha)

            keyword = 'oBSeRVer'
            new_observer = "William Alfred Fowler"
            img.update_keyword(keyword.upper(), new_observer)
            self.assertEqual(img.read_keyword(keyword), new_observer)
            self.assertEqual(get_comment(img, keyword), observer_comment)

            # Create a second FITSImage object referring to the same FITS file.
            # This reads the image from disk again and thus we can verify that
            # not only the in-memory copy of the FITS header was updated.

            img2 = FITSImage(img.path)
            self.assertEqual(img2.read_keyword(keyword), new_observer)
            self.assertEqual(get_comment(img2, keyword), observer_comment)
            self.assertEqual(img._header.items(), img2._header.items())

    def test_delete_keyword(self):

        keywords = {'OBJECT': 'M101', 'INSTRUMENT': 'CanariCam'}
        with self.random(**keywords) as img:

            # Deleting a non-existing keyword raises no error
            with warnings.catch_warnings():
                # Warning emitted until PyFITS 3.2 or 3.3, apparently
                msg = "(?s).+ Please update your code so that KeyErrors are " \
                      "caught and handled when deleting non-existent keywords."
                warnings.filterwarnings('ignore', message=msg)
                img.delete_keyword('TELESCOPE')

            # Keywords are case-insensitive
            keyword = 'oBJeCT'
            img.delete_keyword(keyword)
            with self.assertRaises(KeyError):
                img.read_keyword(keyword.upper())

            # Create a second FITSImage object referring to the same file. This
            # re-reads the image from disk and we can verify that the keyword
            # was not only removed from the in-memory copy of the FITS header.
            img2 = FITSImage(img.path)
            with self.assertRaises(KeyError):
                img2.read_keyword(keyword)

            # Remove a keyword longer than eight characters.
            keyword = 'INSTRUMENT'
            img.delete_keyword(keyword)
            # Try both ways, with and without 'HIERARCH'.
            with self.assertRaises(KeyError):
                img.read_keyword(keyword)
            with self.assertRaises(KeyError):
                img.read_keyword('HIERARCH ' + keyword)

    def test_read_barycentric_date(self):
        keywords = {
            'BJD_TDB': 2458902.321777873,
            'RA': '03:47:24.00',
            'DEC': '+24:15:19.0',
        }
        # From http://astroutils.astronomy.ohio-state.edu/time/bjd2utc.html,
        # we see that for these coordinates, the input BJD_TDB corresponds to
        # JD UTC = 2458902.321287484, which we convert this to a UTC timestamp.
        with self.random(**keywords) as img:
            self.assertAlmostEqual(img.date(barycentric=True), 1582400559.23860979)

    @parameterized.named_parameters(
        {
            'testcase_name': 'yyyy-mm-ddTHH:MM:SS.sss',
            'date': '2009-06-12T21:12:47.238',
            'time': None,
            'exptime': 1505,
            'want_date': 1244841919.738,  # 2009-06-12 21:25:19.738
            'want_year': 2009.44628104,
        },{
            'testcase_name': 'yyyy-mm-ddTHH:MM:SS',
            'date': '1998-11-29T02:53:21',
            'time': None,
            'exptime': 1834,
            'want_date': 912308918.0,  # 1998-11-29 03:08:38.0
            'want_year': 1998.90994793,
        },{
            'testcase_name': 'yyyy-mm-ddTHH:MM:SS (leading and trailing whitespaces)',
            'date': ' 1998-11-29T02:53:21  ',
            'time': None,
            'exptime': 1834,
            'want_date': 912308918.0,  # 1998-11-29 03:08:38.0
            'want_year': 1998.90994793,
        },{
            'testcase_name': 'yyyy-mm-ddTHH:MM:SS.sss (midpoint falls in the following year)',
            'date': '1999-12-31T23:53:12.565',
            'time': None,
            'exptime': 4312.56,
            'want_date': 946686548.845,  # 2000-01-01 00:29:08.845
            'want_year': 2000.00005528,
        },{
            'testcase_name': 'yyyy-mm-dd and HH:MM:SS.sss',
            'date': '2003-03-28',
            'time': '04:23:17.87',
            'exptime': 1357.25,
            'want_date': 1048826076.495,  # 2003-03-28 04:34:36.495
            'want_year': 2003.23613889,
        },{
            'testcase_name': 'yyyy-mm-dd and HH:MM:SS.sss (leading and trailing whitespaces)',
            'date': ' 2003-03-28  ',
            'time': ' 04:23:17.87 ',
            'exptime': 1357.25,
            'want_date': 1048826076.495,  # 2003-03-28 04:34:36.495
            'want_year': 2003.23613889,
        },{
            'testcase_name': 'yyyy-mm-dd_and_HH:MM:SS',
            'date': '2011-06-13',
            'time': '19:04:28',
            'exptime': 2341,
            'want_date': 1307993038.5,  # 2011-06-13 19:23:58.5
            'want_year': 2011.44878989,
        },{
            'testcase_name': 'yy/mm/dd',
            'date': '89/08/30',  # interpreted as 1989
            'time': '01:32:22.14',
            'exptime': 83.55,
            'want_date': 620443983.915,  # 1989-08-30 01:33:03.915
            'want_year': 1989.66045101,
        },{
            'testcase_name': 'yy/mm/dd (another example)',
            'date': '13/04/20',  # interpreted as 1913
            'time': '22:58:12.1',
            'exptime': 1328,
            'want_date': -1789260643.9,  # 1913-04-20 23:09:16.1
            'want_year': 1913.30127334,
        },
    )
    def test_date_and_year(self, date, time, exptime, want_date, want_year):
        keywords = {
            'DATE-OBS': date,
            'TIME-OBS': time,
            'EXPTIME': exptime,
        }
        with self.random(**keywords) as img:
            self.assertAlmostEqual(img.date(), want_date)
            self.assertAlmostEqual(img.year(), want_year)

    def test_date_and_year_non_standard_keywords(self):
        # The FITS image does not use standard keywords, but date() and year()
        # can be pointed to read the date, time and exposure duration from the
        # right keywords.
        keywords = {
            'DATE': '2011-06-13',
            'TIME': '19:04:28',
            'EXPOSURE': 2341,
        }
        with self.random(**keywords) as img:
            got_date = img.date(
                date_keyword='DATE',
                time_keyword='TIME',
                exp_keyword='EXPOSURE',
            )
            got_year = img.year(
                date_keyword='DATE',
                time_keyword='TIME',
                exp_keyword='EXPOSURE',
            )
            self.assertAlmostEqual(got_date, 1307993038.5)  # 2011-06-13 19:23:58.5
            self.assertAlmostEqual(got_year, 2011.44878989)

    @parameterized.named_parameters(
        {
            'testcase_name': 'Missing DATE-OBS',
            'date': None,
            'exptime': 1500,
            'want': KeyError,
        },{
            'testcase_name': 'Missing EXPTIME',
            'date': '1988-09-13T21:45:19',
            'exptime': None,
            'want': KeyError,
        },{
            # DATE-OBS in the 'yyyy-mm-dd' format, so the time at the start of
            # the observation must be read from TIME-OBS (which is not found).
            'testcase_name': 'Missing TIME-OBS',
            'date': '1993-04-29',
            'exptime': 1200,
            'want': KeyError,
        },
    )
    def test_date_and_year_missing_keyword(self, date, exptime, want):
        keywords = {
            'DATE-OBS': date,
            'EXPTIME': exptime,
        }
        with self.random(**keywords) as img:
            with self.assertRaises(want):
                img.date()
            with self.assertRaises(want):
                img.year()

    @parameterized.named_parameters(
        {
            'testcase_name': 'Slashes instead of dashes',
            'date': '1993/04/29T12:42:23',
        },{
            'testcase_name': "'T' surrounded by whitespaces",
            'date': '1993-04-29 T 12:42:23',
        },{
            'testcase_name': "No 'T' between date and time",
            'date': '1993-04-29 12:42:23',
        },{
            'testcase_name': 'Whitespaces instead of colons',
            'date': '1993-04-29T12 42 23',
        },{
            'testcase_name': "Nothing after the 'T'",
            'date': '1993-04-29T',
        },{
            'testcase_name': "Year missing",
            'date': '04-29T12:42:23',
        },{
            'testcase_name': 'Time present but date missing',
            'date': '12:42:23',
        },{
            'testcase_name': 'ctime() format',
            'date': 'Thu Apr 29 12:42:23 1993',
        },
    )
    def test_date_and_year_nonstandard_DATE_OBS(self, date):
        keywords = {
            'DATE-OBS': date,
            'EXPTIME': 1800,
        }
        with self.random(**keywords) as img:
            with self.assertRaisesRegexp(fitsimage.NonStandardFITS, 'DATE-OBS'):
                img.date()
            with self.assertRaisesRegexp(fitsimage.NonStandardFITS, 'DATE-OBS'):
                img.year()

    @parameterized.named_parameters(
        {
            'testcase_name': 'Slashes instead of colons',
            'time': '12/42/23',
        },{
            'testcase_name': 'Whitespaces instead of colons',
            'time': '12 42 23',
        },{
            'testcase_name': 'Seconds missing',
            'time': '12:42',
        },{
            'testcase_name': "'h', 'm' and 's' instead of colons",
            'time': '12h42m23s',
        },{
            'testcase_name': 'Empty',
            'time': '',
        },
    )
    def test_date_and_year_nonstandard_TIME_OBS(self, time):
        keywords = {
            'DATE-OBS': '1998-11-29',
            'TIME-OBS': time,
            'EXPTIME': 1800,
        }
        with self.random(**keywords) as img:
            with self.assertRaisesRegexp(fitsimage.NonStandardFITS, 'TIME-OBS'):
                img.date()
            with self.assertRaisesRegexp(fitsimage.NonStandardFITS, 'TIME-OBS'):
                img.year()

    @parameterized.named_parameters(
        {
            'testcase_name': 'Not a number',
            'exptime': 'unknown',
        },{
            'testcase_name': 'Empty',
            'exptime': '',
        },
    )
    def test_date_and_year_nonstandard_EXPTIME(self, exptime):
        keywords = {
            'DATE-OBS': '2003-03-28',
            'TIME-OBS': '04:23:17.87',
            'EXPTIME': exptime,
        }
        with self.random(**keywords) as img:
            with self.assertRaisesRegexp(fitsimage.NonStandardFITS, 'EXPTIME'):
                img.date()
            with self.assertRaisesRegexp(fitsimage.NonStandardFITS, 'EXPTIME'):
                img.year()

    def test_ra(self):

        ra_kwd = 'RA'

        keywords = {ra_kwd : 344.412916667}
        with self.random(**keywords) as img:
            self.assertAlmostEqual(keywords[ra_kwd], img.ra(ra_kwd))

        keywords = {ra_kwd.lower() : '19:43:56.01'}
        with self.random(**keywords) as img:
            self.assertAlmostEqual(295.983375, img.ra(ra_kwd))

        keywords = {ra_kwd : '17:58:17          '}
        with self.random(**keywords) as img:
            self.assertAlmostEqual(269.570833333, img.ra(ra_kwd))

        keywords = {ra_kwd.lower() : '    00:42:30.997 '}
        with self.random(**keywords) as img:
            self.assertAlmostEqual(10.6291541667, img.ra(ra_kwd))

        keywords = {ra_kwd : '13:00:01.'}
        with self.random(**keywords) as img:
            self.assertAlmostEqual(195.004166667, img.ra(ra_kwd))

        def assertValueErrorRaised(ra_value):
            """ Assert that a random FITS image whose right ascension is
            'ra_value' raises ValueError when its ra() method is called. """

            keywords = {ra_kwd : ra_value}
            with self.random(**keywords) as img:
                regexp = "not in decimal degrees or 'HH:MM:SS.sss' format"
                with self.assertRaisesRegexp(ValueError, regexp):
                    img.ra(ra_kwd)

        # Not in decimal or dd:mm:ss[.sss] format
        assertValueErrorRaised('00h42m30s')
        assertValueErrorRaised('00:42:30.9997')
        assertValueErrorRaised('3:24:57.12')
        assertValueErrorRaised('N/A')

        # 'RA' keyword not in FITS header
        with self.assertRaises(KeyError):
            with self.random() as img:
                img.ra(ra_kwd)

    def test_dec(self):

        dec_kwd = 'DEC'

        keywords = {dec_kwd : -52.6956605556}
        with self.random(**keywords) as img:
            self.assertAlmostEqual(keywords[dec_kwd], img.dec(dec_kwd))

        keywords = {dec_kwd.lower() : '45:59:52.768'}
        with self.random(**keywords) as img:
            self.assertAlmostEqual(45.9979911111, img.dec(dec_kwd))

        keywords = {dec_kwd : '+19:10:56          '}
        with self.random(**keywords) as img:
            self.assertAlmostEqual(19.1822222222, img.dec(dec_kwd))

        keywords = {dec_kwd.lower() : '    -08:12:05.8  '}
        with self.random(**keywords) as img:
            self.assertAlmostEqual(-8.20161111111, img.dec(dec_kwd))

        keywords = {dec_kwd : '26:25:08.'}
        with self.random(**keywords) as img:
            self.assertAlmostEqual(26.4188888889, img.dec(dec_kwd))

        def assertValueErrorRaised(dec_value):
            """ Assert that a random FITS image whose declination is
            'dec_value' raises ValueError when its dec() method is called."""

            keywords = {dec_kwd : dec_value}
            with self.random(**keywords) as img:
                regexp = "not in decimal degrees or 'DD:MM:SS.sss' format"
                with self.assertRaisesRegexp(ValueError, regexp):
                    img.dec(dec_kwd)

        # Not in decimal or dd:mm:ss[.sss] format
        assertValueErrorRaised('+41d16m9s')
        assertValueErrorRaised('45:59:52.7685')
        assertValueErrorRaised('-8:12:05.8')
        assertValueErrorRaised('N/A')

        # 'DEC' keyword not in FITS header
        with self.assertRaises(KeyError):
            with self.random() as img:
                img.dec(dec_kwd)
