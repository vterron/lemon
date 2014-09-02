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

import datetime
import calendar
import numpy.random
import os
import pyfits
import random
import shutil
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

        # Ignore PyFITS warning: "Error validating header for HDU #0 (note:
        # PyFITS uses zero-based indexing). Block does not begin with SIMPLE or
        # XTENSION. There may be extra bytes after the last HDU or the file is
        # corrupted".

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

    def test_unlink(self):
        for _ in xrange(NITERS):
            img = self.random()
            path = img.path
            self.assertTrue(os.path.exists(path))
            img.unlink()
            self.assertFalse(os.path.exists(path))
            self.assertEqual(img.path, None)

    def test_repr(self):
        with self.random() as img1:
            self.assertEqual(img1, eval(repr(img1)))

    def test_eq_and_ne(self):

        with self.random() as img1:
            self.assertEqual(img1, img1)
            self.assertFalse(img1 != img1)

            # A different FITS file
            with self.random() as img2:
                self.assertNotEqual(img1, img2)

            # The same FITS file
            img3 = FITSImage(img1.path)
            self.assertEqual(img1, img3)

            # An exact copy of the FITS file
            copy_path = self.random_data()[0]
            shutil.copy2(img1.path, copy_path)
            with FITSImage(copy_path) as img4:
                self.assertNotEqual(img1.path, img4.path)
                self.assertEqual(img1.sha1sum, img4.sha1sum)
                self.assertNotEqual(img1, img4)

            # A symbolic link to the FITS file
            currdir = os.path.dirname(os.path.realpath(__file__))
            symlink_path = os.path.join(currdir, img1.basename)
            os.symlink(img1.path, symlink_path)
            with FITSImage(symlink_path) as img5:
                self.assertNotEqual(img1.path, img5.path)
                self.assertEqual(img1, img5)

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

    def test_date_and_year(self):

        def strptime_utc(date_string):
            """ Parse a string representing a UTC time according to strptime's
            format '%Y-%m-%d %H:%M:%S.%f' and return it in seconds since the
            Unix epoch. Note that '%s' is not a documented option to strftime()
            as it is OS dependent, but it should work on any decent system."""

            format_ = "%Y-%m-%d %H:%M:%S.%f"
            dt = datetime.datetime.strptime(date_string, format_)
            struct_time = dt.utctimetuple()
            # strftime('.%f') because struct_time does not store fractions
            return calendar.timegm(struct_time) + float(dt.strftime('.%f'))

        def unstrip(str_, max_whitespaces = 5):
            """ Return a copy of 'str' with leading and trailing whitespaces.
            The number of whitespaces on each side is a random number in the
            range [1, max_whitespaces] """

            leading  = ' ' * random.randint(1, max_whitespaces)
            trailing = ' ' * random.randint(1, max_whitespaces)
            return '%s%s%s' % (leading, str_, trailing)

        # (1) The ideal case: date format 'yyyy-mm-ddTHH:MM:SS.sss'
        # Start of observation = 2009-06-12 21:12:47.238, EXPTIME = 1505
        # 2009-06-12 21:12:47.238 + (00:25:05 / 2) =
        # 2009-06-12 21:12:47.238 + 00:12:32.5 =
        # 2009-06-12 21:25:19.738 UTC ~=
        # 2009.44628104

        keywords = {'DATE-OBS' : '2009-06-12T21:12:47.238',
                    'EXPTIME' : 1505}
        with self.random(**keywords) as img:
            expected_date = strptime_utc("2009-06-12 21:25:19.738")
            expected_year = 2009.44628104
            self.assertAlmostEqual(expected_date, img.date())
            self.assertAlmostEqual(expected_year, img.year())

        # (2) Format 'yyyy-mm-ddTHH:MM:SS', without fractions of a second
        # Start of observation = 1998-11-29 02:53:21, EXPTIME = 1834
        # 1998-11-29 02:53:21 + (0:30:34 / 2) =
        # 1998-11-29 02:53:21 + 0:15:17 =
        # 1998-11-29 03:08:38 UTC ~=
        # 1998.90994793

        keywords = {'DATE-OBS' : '1998-11-29T02:53:21',
                    'EXPTIME' : 1834}
        with self.random(**keywords) as img:
            expected_date = strptime_utc("1998-11-29 03:08:38.0")
            expected_year = 1998.90994793
            self.assertAlmostEqual(expected_date, img.date())
            self.assertAlmostEqual(expected_year, img.year())

        # Repeat (2) with leading and trailing whitespaces
        keywords['DATE-OBS'] = unstrip(keywords['DATE-OBS'])
        with self.random(**keywords) as img:
            self.assertAlmostEqual(expected_date, img.date())
            self.assertAlmostEqual(expected_year, img.year())

        # (3) Another one with fractions of a second, turn to a new year
        # Start of observation = 1999-12-31 23:53:12.565, EXPTIME = 4312.56
        # 1999-12-31 23:53:12.565 + (1:11:52.56 / 2) =
        # 1999-12-31 23:53:12.565 + 0:35:56.28 =
        # 2000-01-01 00:29:08.845 UTC =~
        # 2000.00005528

        keywords = {'DATE-OBS' : '1999-12-31T23:53:12.565',
                    'EXPTIME' : 4312.56}
        with self.random(**keywords) as img:
            expected_date = strptime_utc("2000-01-01 00:29:08.845")
            expected_year = 2000.00005528
            self.assertAlmostEqual(expected_date, img.date())
            self.assertAlmostEqual(expected_year, img.year())

        # (4) Format 'yyyy-mm-dd', HH:MM:SS.sss read from TIME-OBS
        # Start of observation = 2003-03-28 04:23:17.87, EXPTIME = 1357.25
        # 2003-03-28 04:23:17.87 + (0:22:37.25 / 2) =
        # 2003-03-28 04:23:17.87 + 0:11:18.625 =
        # 2003-03-28 04:34:36.495 UTC ~=
        # 2003.23613889

        keywords = {'DATE-OBS' : '2003-03-28',
                    'TIME-OBS' : '04:23:17.87',
                    'EXPTIME' : 1357.25}
        with self.random(**keywords) as img:
            expected_date = strptime_utc("2003-03-28 04:34:36.495")
            expected_year = 2003.23613889
            self.assertAlmostEqual(expected_date, img.date())
            self.assertAlmostEqual(expected_year, img.year())

        # Repeat (4) with leading and trailing whitespaces
        keywords['DATE-OBS'] = unstrip(keywords['DATE-OBS'])
        keywords['TIME-OBS'] = unstrip(keywords['TIME-OBS'])
        with self.random(**keywords) as img:
            self.assertAlmostEqual(expected_date, img.date())
            self.assertAlmostEqual(expected_year, img.year())

        # (5) Another 'yyyy-mm-dd', HH:MM:SS (no fractions) read from TIME-OBS
        # Start of observation = 2011-06-13 19:04:28, EXPTIME = 2341
        # 2011-06-13 19:04:28 + (0:39:01 / 2) =
        # 2011-06-13 19:04:28 + 0:19:30.5 =
        # 2011-06-13 19:23:58.5 UTC ~=
        # 2011.44878989

        keywords = {'DATE-OBS' : '2011-06-13',
                    'TIME-OBS' : '19:04:28',
                    'EXPTIME' : 2341}
        with self.random(**keywords) as img:
            expected_date = strptime_utc("2011-06-13 19:23:58.5")
            expected_year = 2011.44878989
            self.assertAlmostEqual(expected_date, img.date())
            self.assertAlmostEqual(expected_year, img.year())

        # (6) Repeat (5), but using different, non-standard keywords
        keywords['DATE'] = keywords.pop('DATE-OBS')
        keywords['TIME'] = keywords.pop('TIME-OBS')
        keywords['EXPOSURE'] = keywords.pop('EXPTIME')
        with self.random(**keywords) as img:
            kwargs = dict(date_keyword = 'DATE',
                          time_keyword = 'TIME',
                          exp_keyword = 'EXPOSURE')
            self.assertAlmostEqual(expected_date, img.date(**kwargs))
            self.assertAlmostEqual(expected_year, img.year(**kwargs))

        # (7) Format 'yy/mm/dd' (only for dates from 1900 through 1999)
        # The year 13, when this format is used, is interpreted as 1913
        # Start of observation = 1913-04-20 22:58:12.1, EXPTIME = 1328
        # 1913-04-20 22:58:12.1 + (0:22:08 / 2) =
        # 1913-04-20 22:58:12.1 + 0:11:04 =
        # 1913-04-20 23:09:16.1 UTC ~=
        # 1913.30127334

        keywords = {
            'DATE-OBS' : '13/04/20',
            'TIME-OBS' : '22:58:12.1',
            'EXPTIME' : 1328}
        with self.random(**keywords) as img:
            expected_date = strptime_utc("1913-04-20 23:09:16.1")
            expected_year = 1913.30127334
            self.assertAlmostEqual(expected_date, img.date())
            self.assertAlmostEqual(expected_year, img.year())

        # (8) Another 'yy/mm/dd' format: year 89 means here 1989
        # Start of observation = 1989-08-30 01:32:22.14, EXPTIME = 83.55
        # 1989-08-30 01:32:22.14 + (0:01:23.55 / 2) =
        # 1989-08-30 01:32:22.14 + 0:00:41.775 =
        # 1989-08-30 01:33:03.915 ~=
        # 1989.66045101

        keywords =  {
            'DATE-OBS' : '89/08/30',
            'TIME-OBS' : '01:32:22.14',
            'EXPTIME' : 83.55}
        with self.random(**keywords) as img:
            expected_date = strptime_utc("1989-08-30 01:33:03.915")
            expected_year = 1989.66045101
            self.assertAlmostEqual(expected_date, img.date())
            self.assertAlmostEqual(expected_year, img.year())

        def assertRaised(exception, keywords):
            """ Assert that a random FITS image with the 'keywords' FITS
            keywords raises 'exception' when both its date() and year()
            methods are called."""
            with self.random(**keywords) as img:
                with self.assertRaises(exception):
                    img.date()
                with self.assertRaises(exception):
                    img.year()

        # Keyword DATE-OBS not found
        keywords = {'EXPTIME' : 1500}
        assertRaised(KeyError, keywords)

        # Keyword EXPTIME not found
        keywords = {'DATE-OBS' : '1988-09-13T21:45:19'}
        assertRaised(KeyError, keywords)

        # DATE-OBS in the 'yyyy-mm-dd' format, so the time at the start of
        # the observation must be read from TIME-OBS (which is not found)
        keywords = {'DATE-OBS' : '1993-04-29', 'EXPTIME' : 1200}
        assertRaised(KeyError, keywords)

        # Both DATE-OBS and EXPTIME exist, but the latter does not contain a
        # floating-point number giving the exposure time in units of seconds.
        keywords = {'DATE-OBS' : '1988-09-13T21:45:19', 'EXPTIME' : 'unknown'}
        assertRaised(fitsimage.NonStandardFITS, keywords)

        # Several non-standard DATE-OBS formats
        nonstandard_dates = [
            '1993/04/29T12:42:23',      # slashes, not dashes
            '1993-04-29 T 12:42:23',    # 'T' surrounded by whitespaces
            '1993-04-29 12:42:23',      # no 'T' between date and time
            '1993-04-29T12 42 23',      # whitespaces, not colons
            '1993-04-29T',              # nothing after the 'T'
            '04-29T12:42:23',           # year missing
            '12:42:23',                 # time, but not date
            'Thu Apr 29 12:42:23 1993'] # ctime() format

        for wrong_date in nonstandard_dates:
            keywords['DATE-OBS'] = wrong_date
            assertRaised(fitsimage.NonStandardFITS, keywords)

        # Several non-standard TIME-OBS formats
        nonstandard_times = [
            '12/42/23',  # slashes, not colons
            '12 42 23',  # whitespaces, not colons
            '12:42',     # seconds missing
            '12h42m23s', # 'h', 'm' and 's' instead of colons
            '']          # empty

        keywords = {'DATE-OBS' : '1993-04-29', 'EXPTIME' : 1200}
        for wrong_time in nonstandard_times:
            keywords['TIME-OBS'] = wrong_time
            assertRaised(fitsimage.NonStandardFITS, keywords)

