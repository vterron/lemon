#! /usr/bin/env python2
# encoding: UTF-8

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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division

# Tell PyRAF to skip all graphics initialization and run in terminal-only mode
# via. Otherwise we will get annoying warning messages (such as "could not open
# XWindow display" or "No graphics display available for this session") when
# wprking at a remote terminal or a terminal without any X Windows support. Any
# tasks which attempt to display graphics will fail, of course, but we are not
# going to make use of any of them, anyway.


import astropy.time
import astropy.wcs
import barycorr
import calendar
import collections
import datetime
import fnmatch
import hashlib
import itertools
import logging
import numpy
import numbers
import os
import os.path
import pyfits
import re
import warnings

# LEMON modules
import util
import keywords
import passband
import style

# Ignore the "adding a HIERARCH keyword" PyFITS warning that is emitted when
# the keyword name does not start with 'HIERARCH' and is greater than eight
# characters or contains spaces. This solution is taken from:
# https://github.com/geminiutil/geminiutil/commit/9aa46fd9cd3
warnings.filterwarnings('ignore', message=".+ a HIERARCH card will be created.")

class NonStandardFITS(IOError):
    """ Raised when a non-standard file is attempted to be opened."""
    pass

class NoWCSInformationError(ValueError):
    """ Raised if WCS information is not found in a FITS header. """
    pass

class FITSImage(object):
    """ Encapsulates a FITS image located in the filesystem. """

    def __init__(self, path):
        """ Instantiation method for the FITSImage class.

        A copy of the header of the FITS file is kept in memory for fast access
        to its keywords. IOError is raised if 'path' does not exist or is not
        readable, and NonStardardFITS in case that it does not conform to the
        FITS standard or it simply is not a FITS file at all.

        FITS Standard Document:
        http://fits.gsfc.nasa.gov/fits_standard.html

        An image is considered to follow the FITS standard if it has 'SIMPLE'
        as its first keyword of the primary header. According to the standard,
        it contains a logical constant with the value 'T' if the file conforms
        to it. This keyword is mandatory for the primary header and is not
        permitted in extension headers. A value of 'F', on the other hand,
        means that the file does not conform to the standard.

        We trust the 'SIMPLE' keyword blindly: if is says that the FITS image
        follows the standard, we believe it. Period. We do not consider the
        possibility (although this may change in the future, if we begin to
        work with much less reliable data) that the keyword has a value of 'T'
        while at the same time there are violations of the standard. That is
        why we instruct PyFITS to ignore any FITS standard violations we come
        across (output_verify = 'ignore').

        """

        if not os.path.exists(path):
            raise IOError("file '%s' does not exist" % path)

        self.path = path

        try:
            # The file must be opened to make sure it is a standard FITS.
            # We would rather use the with statement, but in that case we
            # would not be able to set the output verification of close() to
            # 'ignore'. Thus, the default option, 'exception', which raises
            # an exception if any FITS standard is violated, would be used.
            handler = pyfits.open(self.path, mode = 'readonly')
            try:

                # This is kind of an ugly hack to make sure that all PyFITS
                # versions are supported. Up to version 3.2, PyFITS returned
                # the HDUs exactly as it saw them. This allowed us to check the
                # existence and value of the 'SIMPLE' keyword in order to make
                # sure the file conforms to the FITS standard. However, PyFITS
                # 3.3 changed this, and now adds the keywords required for a
                # minimal viable primary HDU: this means that the header will
                # always contain the 'SIMPLE' keyword. Refer to this link for
                # more info: https://github.com/spacetelescope/PyFITS/issues/94

                try:
                    type_ = pyfits.info(self.path, output=False)[0][2]
                    if type_ == 'NonstandardHDU':
                        # 'SIMPLE' exists but does not equal 'T'
                        msg = "%s: value of 'SIMPLE' keyword is not 'T'"
                        raise NonStandardFITS(msg % self.path)

                except AttributeError as e:
                    # 'SIMPLE' keyword does not exist
                    error_msg = "'_ValidHDU' object has no attribute '_summary'"
                    assert error_msg in str(e)
                    msg = "%s: 'SIMPLE' keyword missing from header"
                    raise NonStandardFITS(msg % self.path)

                # A copy of the FITS header is kept in memory and the file is
                # closed; otherwise we may run into trouble when working with
                # thousands of images ("too many open files" and such). This
                # approach gives us fast read-only access to the image header;
                # if modified, we will have to take care of 'reloading' (call
                # it synchronize, if you wish) the header.

                self.size = handler[0].data.shape[::-1]
                self._header = handler[0].header
            finally:
                handler.close(output_verify = 'ignore')

        # PyFITS raises IOError if we do not have permission to open the file,
        # if we attempt to open a non-FITS file, and also if we open one whose
        # first keyword is not either SIMPLE or XTENSION. Nothing is raised if
        # the value of SIMPLE is 'F'; that is why we had to specifically make
        # sure it was 'T' a few lines above.
        except IOError, e:
            pyfits_msg = "Block does not begin with SIMPLE or XTENSION"
            if str(e) == pyfits_msg:
                msg = "%s: 'SIMPLE' keyword missing from primary header"
                raise NonStandardFITS(msg % self.path)
            elif "Permission denied" in str(e):
                raise
            else:
                msg = "%s (%s)" % (self.path, str(e))
                raise NonStandardFITS(msg)

    def __repr__(self):
        """ The unambiguous string representation of a FITSImage object """
        return "%s(%r)" % (self.__class__.__name__, self.path)

    def read_keyword(self, keyword):
        """ Read a keyword from the header of the FITS image.

        Note that, although always upper-case in the FITS header, keywords are
        here case-insensitive, for user's convenience. TypeError is raised if
        the keyword is None, ValueError if it is left empty and KeyError if the
        keyword cannot be found in the header.

        There is no need to prepend 'HIERARCH' to the keyword name when it is
        longer than eight characters, since PyFITS handles this transparently.
        You may nevertheless use it, but note in that case there *must* be a
        whitespace between 'HIERARCH' and the keyword name: e.g., you must
        write 'HIERARCH AMBI WIND SPEED', never 'HIERARCHAMBI WIND SPEED'.

        """

        if keyword is None:
            raise TypeError("keyword cannot be None")
        if not keyword:
            raise ValueError("keyword cannot be empty")
        try:
            return self._header[keyword.upper()]
        except KeyError:
            msg = "%s: keyword '%s' not found" % (self.path, keyword)
            raise KeyError(msg)

    def update_keyword(self, keyword, value, comment = None):
        """ Updates the value of a FITS keyword, adding it if it does not exist.

        The method updates the value of a keyword in the FITS header, replacing
        it with the specified value or simply adding it in case if does not yet
        exist. Note that, although always upper-case inside the FITS file,
        keywords are here case-insensitive, for user's convenience.

        Raises ValueError if a HIERARCH keyword (that is, a keyword longer than
        eight characters or that contains spaces) and its value exceed eighty
        characters. The reason for this limitation is that PyFITS does not
        support CONTINUE for HIERARCH. If the value is too long, therefore,
        make sure that the keyword does not need to be HIERARCH-ed.

        Keyword arguments:
        comment - the comment to be added to the keyword.

        """

        if len(keyword) > 8:
            msg = "%s: keyword '%s' is longer than eight characters or " \
                  "contains spaces; a HIERARCH card will be created"
            logging.debug(msg % (self.path, keyword))

        handler = pyfits.open(self.path, mode = 'update')
        msg = "%s: file opened to update '%s' keyword" % (self.path, keyword)
        logging.debug(msg)

        try:
            header = handler[0].header

            # Ignore the 'card is too long, comment is truncated' warning
            # printed by PyRAF in case, well, the comment is too long.
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                header[keyword] = (value, comment)
                args = self.path, keyword, value
                msg = "%s: keyword '%s' updated to '%s'" % args
                if comment:
                    msg += " with comment '%s'" % comment
                logging.debug(msg)

            # Update in-memory copy of the FITS header
            self._header = header

        except ValueError, e:

            # ValueError is raised if a HIERARCH keyword is used and the total
            # length (keyword, equal sign string and value) is greater than 80
            # characters. The default exception message is a bit cryptic ("The
            # keyword {...} with its value is too long"), so add some more
            # information to help the user understand what went wrong.

            pattern = "The keyword .*? with its value is too long"
            if re.match(pattern, str(e)):
                assert len(keyword) > 8
                msg = ("%s: keyword '%s' could not be updated (\"%s\"). Note "
                       "that PyFITS does not support CONTINUE for HIERARCH. "
                       "In other words: if your keyword has more than eight "
                       "characters or contains spaces, the total length of "
                       "the keyword with its value cannot be longer than %d "
                       "characters.")
                args = self.path, keyword, str(e), pyfits.Card.length
                logging.warning(msg % args)
                raise ValueError(msg % args)
            else:
                # Different ValueError, re-raise it
                msg = "%s: keyword '%s' could not be updated (%s)"
                args = self.path, keyword, e
                logging.warning(msg % args)
                raise

        except Exception, e:
            msg = "%s: keyword '%s' could not be updated (%s)"
            args = self.path, keyword, e
            logging.warning(msg % args)
            raise

        finally:
            handler.close(output_verify = 'ignore')
            msg = "%s: file closed" % self.path
            logging.debug(msg)

    def delete_keyword(self, keyword):
        """ Delete a keyword from the header of the FITS image.

        The method removes from the header of the image the specified keyword.
        Unlike read_keyword(), no exception is raised if the keyword is not
        present in the header. Keywords are case-insensitive.

        """

        handler = pyfits.open(self.path, mode = 'update')
        try:
            header = handler[0].header
            try:
                # Ignore DeprecationWarning: "Deletion of non-existent
                # keyword [...] In a future PyFITS version Header.__delitem__
                # may be changed so that this raises a KeyError just like a
                # dict would. Please update your code so that KeyErrors are
                # caught and handled when deleting non-existent keywords.
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    del header[keyword]
                # Update in-memory copy of the FITS header
                self._header = header

            # Future versions of PyFITS (by 3.2 or 3.3, most probably) will
            # raise KeyError when a non-existent keyword is deleted, just
            # like a dictionary would, so we better get ready for this.
            except KeyError:
                pass
        finally:
            handler.close(output_verify = 'ignore')

    def add_history(self, history):
        """ Add another record to the history of the FITS image.

        The 'HISTORY' keyword shall have, according to the FITS Standard, no
        associated value; columns 9-80 may contain any ASCII text. The text
        should contain a history of steps and procedures associated with the
        processing of the associated data. Any number of HISTORY card images
        may appear in a header. The in-memory copy of the header is not
        updated.

        """

        handler = pyfits.open(self.path, mode = 'update')
        try:
            header = handler[0].header
            header.add_history(history)
        finally:
            handler.close(output_verify = 'ignore')

    def date(self,
             date_keyword=keywords.datek,
             time_keyword=keywords.timek,
             exp_keyword=keywords.exptimek,
             barycentric=False,
             bjd_keyword=keywords.bjdk,
        ):
        """ Return the date at mid-exposure (UTC), in Unix time.

        This function always returns a UTC Unix timestamp. By default, the date
        at mid-exposure is derived using the traditional, standard FITS keywords
        (e.g. DATE-OBS, EXPTIME). However, if 'barycentric' is set to True, the
        value used is the Barycentric Julian Date in Barycentric Dynamical Time
        (BJD_TDB)). This is a non-standard FITS keyword, but commonly used in
        the exoplanet community for high-accuracy observations. For more
        information, see http://astroutils.astronomy.ohio-state.edu/time/.
        """

        if barycentric:
            return self._read_barycentric_date(bjd_keyword)
        return self._read_uncorrected_date(
            date_keyword=date_keyword,
            time_keyword=time_keyword,
            exp_keyword=exp_keyword)

    def _read_barycentric_date(self, bjd_keyword):
        """Returns the Barycentric Julian Date as a UTC timestamp."""

        bjd_tdb = self.read_keyword(bjd_keyword)
        jd_utc = barycorr.bjd2utc(bjd_tdb, self.ra(), self.dec())
        t = astropy.time.Time(jd_utc, format='jd', scale='utc')
        return t.unix

    def _read_uncorrected_date(self, date_keyword, time_keyword, exp_keyword):
        """ Return the date of observation (UTC), in Unix time.

        This method returns, in seconds after the Unix epoch (the time 00:00:00
        UTC on 1 January 1970), the date at the 'midpoint' of the observation.
        This is defined as the time of the start of the exposure plus one-half
        the exposure duration.

        This function makes no assumptions or corrections for the time reference
        frame (e.g. geocentric, heliocentric or barycentric; i.e. the different
        geometric locations from which one could measure time, differing by the
        light-travel time between them). It simply just returns the date found
        in the header as a UTC timestamp.

        The KeyError exception is raised if any of the specified keywords
        cannot be found in the FITS header. NonStandardFITS is thrown if the
        keywords exist but they do not follow the FITS standard. See the URLs
        below for more information on the standard and other popular keywords:

        Definition of the Flexible Image Transport System (FITS):
        http://archive.stsci.edu/fits/fits_standard/

        Dictionary of Commonly Used FITS Keywords:
        http://heasarc.gsfc.nasa.gov/docs/fcg/common_dict.html

        Keyword arguments:
        date_kewyord - the FITS keyword in which the date of the observation is
                       stored, in the format specified in the FITS Standard. The
                       old date format was 'yy/mm/dd' and may be used only for
                       dates from 1900 through 1999. The new Y2K compliant date
                       format is 'yyyy-mm-dd' or 'yyyy-mm-ddTHH:MM:SS[.sss]'.
        time_keyword - FITS keyword storing the time at which the observation
                       started, in the format HH:MM:SS[.sss]. This keyword is
                       ignored (and, thus, should not be used) if the time is
                       included directly as part of the 'date_keyword' keyword
                       value with the format 'yyyy-mm-ddTHH:MM:SS[.sss]'.
        exp_keyword - the FITS keyword in which the duration of the exposure is
                      stored. It is expected to be a floating-point number which
                      gives the duration in seconds. The exact definition of
                      'exposure time' is mission dependent and may, for example,
                      include corrections for shutter open and close duration,
                      detector dead time, vignetting, or other effects.

        """

        # Throws KeyError is the specified keyword is not found in the header
        start_date_str = self.read_keyword(date_keyword).strip()

        # Time format string, needed by strptime()
        format_str = '%Y-%m-%dT%H:%M:%S'

        # There are two formats which only store the date, without the time:
        # 'yy/dd/dd' (deprecated, may be used only for dates 1900-1999) and
        # 'yyyy-mm-dd' (new, Y2K-compliant format).
        old_date_regexp = '\d{2}/\d{2}/\d{2}'
        date_regexp = '\d{4}-\d{2}-\d{2}' # 'yyyy-mm-dd'

        # 'HH:MM:SS[.sss]': the time format. Note that we allow up to four
        # decimals in the seconds [.ssss], instead of three. This is not
        # standard, but used anyway in the header of O2K CAHA FITS images.
        time_regexp = '\d{2}:\d{2}:\d{2}(?P<secs_fraction>\.\d{0,4})?'

        # 'yyyy-mm-ddTHH:MM:SS[.sss]: the format that ideally we would always
        # come across. It contains both the date and time at the start of the
        # observation, so there is no need to read a second keyword (TIME-OBS,
        # for example) to extract the time.
        complete_regexp = '%sT%s' % (date_regexp, time_regexp)
        match = re.match(complete_regexp, start_date_str)
        if match:
            if match.group('secs_fraction'):
                format_str += '.%f'

        # Non-ideal scenario: 'date_keyword' (e.g., DATE-OBS) does not include
        # the time, so we need to read it from a second keyword, 'time_keyword'
        # (e.g., TIME-OBS).

        else:

            # Must be either 'yy/mm/dd' or 'yyyy-mm-dd'
            args = old_date_regexp, date_regexp
            regexp = '^((?P<old>%s)|(?P<new>%s))$' % args
            match = re.match(regexp, start_date_str)

            if match:

                # Translate 'yy/mm/dd' dates to the 'yyyy-mm-dd' format. Note
                # that, according to the standard, dates using the old format
                # may be used only for years [1900, 1999]. Therefore, a date
                # such as '02/04/21' *will* be interpreted as '1902/04/21'.
                if match.group('old'):
                    start_date_str = '19' + start_date_str.replace('/', '-')

                # At this point the format of the date must be 'yyyy-mm-dd'
                assert re.match('^%s$' % date_regexp, start_date_str)

                # Read the time from its keyword. Does it follow the standard?
                start_time_str = self.read_keyword(time_keyword).strip()
                regexp = '^%s$' % time_regexp
                time_match = re.match(regexp, start_time_str)

                if not time_match:
                    args = time_keyword, start_time_str
                    msg = ("'%s' keyword (%s) does not follow the FITS "
                           "standard: 'HH:MM:SS[.sss]'")
                    raise NonStandardFITS(msg % args)

                start_date_str += 'T%s' % start_time_str
                if time_match.group('secs_fraction'):
                    format_str += '.%f'

            else:
                args = date_keyword, start_date_str
                msg = ("'%s' keyword (%s) does not follow the FITS standard: "
                       "yyyy-mm-dd[THH:MM:SS[.sss]] or yy/dd/mm (deprecated)")
                raise NonStandardFITS(msg % args)

        start_date = datetime.datetime.strptime(start_date_str, format_str)

        try:
            # Divide the duration of the exposure, in seconds, by two
            exposure_time = self.read_keyword(exp_keyword)
            half_exp_time = float(exposure_time) / 2
        except KeyError:    # if the keyword is not found
            raise
        except ValueError:  # if the casting to float fails
            args = exp_keyword, exposure_time
            msg = "'%s' keyword (%s) is not a floating-point number"
            raise NonStandardFITS(msg % args)

        # strftime('.%f') is needed because struct_time does not store
        # fractions of second: we need to add it after the conversion from
        # datetime to seconds [http://stackoverflow.com/a/698279/184363]
        seconds_fraction = float(start_date.strftime('.%f'))
        start_struct_time = start_date.utctimetuple()
        start_date = calendar.timegm(start_struct_time) + seconds_fraction
        return start_date + half_exp_time

    def year(self, **kwargs):
        """ Return the date of observation (UTC) as a fractional year.

        Run FITSImage.date() and convert the returned timestamp to a fractional
        year. The start of the year is at time t0, while the end of the year is
        at time t1: if we take any time between those two, we have a fraction.
        This method is probably accurate to within the second and it also works
        correctly during leap years. The function signature is the same as that
        of FITSImage.date()

        Taken from ninjagecko's answer on Stack Overflow:
        [URL] https://stackoverflow.com/a/6451892/184363

        """

        def since_epoch(date):
            """ From datetime to seconds since epoch """
            return calendar.timegm(date.timetuple())

        timestamp = self.date(**kwargs)
        date = datetime.datetime.utcfromtimestamp(timestamp)

        year = date.year
        t0 = datetime.datetime(year = year,     month = 1, day = 1)
        t1 = datetime.datetime(year = year + 1, month = 1, day = 1)

        year_elapsed  = since_epoch(date) - since_epoch(t0)
        year_duration = since_epoch(t1)   - since_epoch(t0)
        fraction = year_elapsed / year_duration
        return year + fraction

    def pfilter(self, keyword):
        """ Return the photometric filter of the image as a Passband instance.

        Read the 'keyword' keyword from the header of the FITS image and
        encapsulate it as a passband.Passband object. For user's convenience,
        the keyword is case-insensitive. This method can raise four different
        exceptions: TypeError if 'keyword' is None, ValueError if it is left
        empty, KeyError if it cannot be found in the FITS header and, lastly,
        passband.NonRecognizedPassband if the name of the photometric filter
        cannot be parsed.

        """

        try:
            pfilter_str = self.read_keyword(keyword)
            return passband.Passband(pfilter_str)
        except passband.NonRecognizedPassband:
            kwargs = dict(path = self.path, keyword = keyword)
            raise passband.NonRecognizedPassband(pfilter_str, **kwargs)

    def ra(self, ra_keyword = keywords.rak):
        """ Return the right ascension, in decimal degrees.

        Return the value in decimal degrees of the 'ra_keyword' FITS keyword.
        The right ascension must be expressed either as a floating point number
        in units of decimal degrees, or as a string in the 'hh:mm:ss[.sss]'
        format; otherwise, ValueError is raised. If the keyword cannot be
        found in the FITS header, raise KeyError.

        """

        msg = "{0}: reading α from FITS header (keyword '{1}')"
        logging.debug(msg.format(self.path, ra_keyword))
        ra_str = self.read_keyword(ra_keyword)

        try:
            ra = float(ra_str)
            logging.debug("{0}: α = {1:.5f}".format(self.path, ra))
            return ra

        except ValueError, e:

            msg1 = "{0}: {1}".format(self.path, str(e))
            msg2 = "{0}: parsing α as sexagesimal".format(self.path)
            msg3 = "{0}: α = '{1}'".format(self.path, ra_str)

            logging.debug(msg1)
            logging.debug(msg2)
            logging.debug(msg3)

            # HH:MM:SS[.sss]
            regexp = '^(?P<hh>\d{2}):(?P<mm>\d{2}):(?P<ss>\d{2}(\.\d{0,3})?)$'
            match = re.match(regexp, ra_str.strip())
            if match:
                hh =   int(match.group('hh'))
                mm =   int(match.group('mm'))
                ss = float(match.group('ss'))

                logging.debug(  "{0}: hours = {1} (α)".format(self.path, hh))
                logging.debug("{0}: minutes = {1} (α)".format(self.path, mm))
                logging.debug("{0}: seconds = {1} (α)".format(self.path, ss))

                ra = util.HMS_to_DD(hh, mm, ss)
                logging.debug("{0}: α = {1:.5f}".format(self.path, ra))
                return ra
            else:
                msg = "{0}: '{1}' not in decimal degrees or 'HH:MM:SS.sss' format (got = {2!r})"
                raise ValueError(msg.format(self.path, ra_keyword, ra_str))

    def dec(self, dec_keyword = keywords.deck):
        """ Return the declination, in decimal degrees.

        Return the value in decimal degrees of the 'dec_keyword' FITS keyword.
        The declination must be expressed either as a floating point number
        in units of decimal degrees, or as a string in the 'dd:mm:ss[.sss]'
        format; otherwise, ValueError is raised. If the keyword cannot be
        found in the FITS header, raise KeyError.

        """

        msg = "{0}: reading δ from FITS header (keyword '{1}')"
        logging.debug(msg.format(self.path, dec_keyword))
        dec_str = self.read_keyword(dec_keyword)

        try:
            dec = float(dec_str)
            logging.debug("{0}: δ = {1:.5f}".format(self.path, dec))
            return dec

        except ValueError, e:

            msg1 = "{0}: {1}".format(self.path, str(e))
            msg2 = "{0}: parsing α as sexagesimal".format(self.path)
            msg3 = "{0}: δ = '{1}'".format(self.path, dec_str)

            logging.debug(msg1)
            logging.debug(msg2)
            logging.debug(msg3)

            # DD:MM:SS[.sss]
            regexp = '^(?P<dd>([-+])?\d{2}):(?P<mm>\d{2}):(?P<ss>\d{2}(\.\d{0,3})?)$'
            match = re.match(regexp, dec_str.strip())
            if match:
                dd =   int(match.group('dd'))
                mm =   int(match.group('mm'))
                ss = float(match.group('ss'))

                logging.debug("{0}: degrees = {1} (δ)".format(self.path, dd))
                logging.debug("{0}: minutes = {1} (δ)".format(self.path, mm))
                logging.debug("{0}: seconds = {1} (δ)".format(self.path, ss))

                dec = util.DMS_to_DD(dd, mm, ss)
                logging.debug("{0}: δ = {1:.5f}".format(self.path, dec))
                return dec
            else:
                msg = "{0}: '{1}' not in decimal degrees or 'DD:MM:SS.sss' format (got = {2!r}"
                raise ValueError(msg.format(self.path, dec_keyword, dec_str))

    @property
    def prefix(self):
        """ Extract the leftmost non-numeric substring of the image base name.

        The method returns the leftmost substring, entirely composed of
        non-numeric characters, that can be found in the extension-stripped
        basename of the FITS image. For example, for the image
        /caha/ferM_0013.fits, the returned string would be 'ferM_'.
        Note that, for an image whose filename contained no numbers, such as
        'ferM_no_number.fit', only 'ferM_no_number' would be returned.

        """

        str_char = ''
        basename = os.path.basename(self.path)
        root = os.path.splitext(basename)[0]
        for character in root:
            try:
                int(character)
                return str_char # stop as soon as an integer is found
            except ValueError:
                # Raised because the character is not an integer
                str_char += character
        return str_char

    @property
    def x_size(self):
        """ Return the number of pixels of the image in the x-axis."""
        return self.size[0]

    @property
    def y_size(self):
        """ Return the number of pixels of the image in the y-axis."""
        return self.size[1]

    @property
    def center(self):
        """ Returns the x, y coordinates of the central pixel of the image. """
        return list(int(round(x / 2)) for x in self.size)

    @util.memoize
    def _get_wcs(self):
        """ Return the astropy.wcs.WCS object for the header of this image. """

        # astropy.wcs.WCS() is extremely slow (in the order of minutes) if we
        # work with the in-memory FITS header (self._header). I cannot fathom
        # the reason, but the problem goes away if we use astropy.io.fits to
        # load the FITS header, as illustrated in the Astropy documentation:
        # http://docs.astropy.org/en/stable/wcs/index.html
        with astropy.io.fits.open(self.path) as hdulist:
            header = hdulist[0].header

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return astropy.wcs.WCS(header)

    def pix2world(self, x, y):
        """ Transform pixel coordinates to world coordinates.

        Return a two-element tuple with the right ascension and declination to
        which the specified x- and y-coordinates correspond in the FITS image.
        Raises NoWCSInformationError if the header of the FITS image does not
        contain an astrometric solution -- i.e., if the astropy.wcs.WCS class
        is unable to recognize it as such. This is something that should very
        rarely happen, and almost positively caused by non-standard systems or
        FITS keywords.

        """

        coords = (x, y)
        wcs = self._get_wcs()
        pixcrd = numpy.array([coords])
        ra, dec = wcs.all_pix2world(pixcrd, 1)[0]

        # We could use astropy.wcs.WCS.has_celestial for this, but as of today
        # [Tue Jan 20 2015] it is only available in the development version of
        # Astropy. Therefore, do a simple (but in theory enough) check: if the
        # header does not contain an astrometric solution, WCS.all_pix2world()
        # will not be able to transform the pixel coordinates, and therefore
        # will return the same coordinates as the FITSImage.center attribute.

        if (ra, dec) == coords:
            msg = ("{0}: the header of the FITS image does not seem to "
                   "contain WCS information. You may want to make sure that "
                   "the image has been solved astrometrically, for example "
                   "with the 'astrometry' LEMON command.".format(self.path))
            raise NoWCSInformationError(msg)

        return ra, dec

    def center_wcs(self):
        """ Return the world coordinates of the central pixel of the image.

        Transform the pixel coordinates of the center of the image to world
        coordinates. Return a two-element tuple with the right ascension and
        declination. Most of the time the celestial coordinates of the field
        center can be found in the FITS header, but (a) these keywords are
        non-standard and (b) it would not be the first time that we come
        across incorrect (or, at least, not as accurate as we would expect)
        coordinates. Instead of blindly trusting the FITS header, compute these
        values ourselves -- provided, of course, that our images are calibrated
        astrometrically.

        Raises NoWCSInformationError if the header of the FITS image does not
        contain an astrometric solution -- i.e., if the astropy.wcs.WCS class
        is unable to recognize it as such. This is something that should very
        rarely happen, and almost positively caused by non-standard systems or
        FITS keywords.

        """

        center = tuple(self.center)
        return self.pix2world(*center)

    def has_wcs(self):
        """ Check whether the header of the image contains WCS information. """

        try:
            self.center_wcs()
            return True
        except NoWCSInformationError:
            return False

    def saturation(self, maximum, coaddk = keywords.coaddk):
        """ Return the effective saturation level, in ADUs.

        Return the number of ADUs at which saturation arises. This is the
        result of multiplying 'maximum', the saturation level of a single
        frame, by the number of coadded frames, which is read from the keyword
        'coaddk'. If the keyword cannot be found in the FITS header, a value of
        one is assumed. TypeError must be an integer or real number, while the
        number of coadded frames must be an integer; otherwise, TypeError is
        raised.

        Keyword arguments:
        coaddk - FITS keyword for the number of coadded frames.

        """

        if not isinstance(maximum, numbers.Real):
            msg = "'maximum' must be an integer or real number"
            raise TypeError(msg)

        msg = "%s: saturation level of a single image: %d ADUs"
        logging.debug(msg % (self.path, maximum))

        try:
            ncoadds = self.read_keyword(coaddk)
        except KeyError:
            ncoadds = 1
            msg = "%s: keyword '%s' not in header, assuming value of one"
            logging.debug(msg % (self.path, coaddk))

        msg = "%s: value of keyword '%s' must be a positive integer"
        args = self.path, coaddk
        error_msg = msg % args

        if not isinstance(ncoadds, numbers.Integral):
            raise TypeError(error_msg)

        if ncoadds < 1:
            raise ValueError(error_msg)

        msg = "%s: number of effective coadds (%s keyword) = %d"
        args = self.path, coaddk, ncoadds
        logging.debug(msg % args)

        saturation = maximum * ncoadds
        msg = "%s: effective saturation level = %d x %d = %d"
        args = self.path, maximum, ncoadds, saturation
        logging.debug(msg % args)
        return saturation

    @property
    def sha1sum(self):
        """ Return the hexadecimal SHA-1 checksum of the FITS image """

        sha1 = hashlib.sha1()
        with open(self.path, 'rb') as fd:
            for line in fd:
                sha1.update(line)
            return sha1.hexdigest()


class InputFITSFiles(collections.defaultdict):
    """ Map each photometric filter to a list of FITS files.

    A convenience class to simplify the manipulation of a series of FITS
    images: it is a defaultdict, mapping each photometric filter to a list of
    the corresponding FITS files, but that iterates over the values (the FITS
    files, independently of their filter). This way, it may be viewed as a
    sequence of FITS files that also allows access by the photometric filter.

    """

    def __init__(self):
        super(InputFITSFiles, self).__init__(list)

    def __iter__(self):
        """ Iterate over the FITS files, regardless or their filter """
        return itertools.chain.from_iterable(self.itervalues())

    def __len__(self):
        """ Return the number of FITS files, in any filter """
        return sum(len(pfilter) for pfilter in self.itervalues())

    def remove(self, img):
        """ Remove all occurrences of the FITS file, in any filter.

        Loop over the different photometric filters and remove, from each
        associated list of FITS files, all occurrences of 'img'. Returns the
        number of elements that were deleted. Emits a warning for each file
        that is removed, stating that is has been 'excluded'.

        """

        discarded = 0
        for pfilter_imgs in self.itervalues():
            # Iterate backwards, modify original list in situ
            for index in reversed(xrange(len(pfilter_imgs))):
                if pfilter_imgs[index] == img:
                    del pfilter_imgs[index]
                    discarded += 1
                    msg = "%s %s excluded."
                    warnings.warn(msg % (style.prefix, img))
        return discarded


def find_files(paths, followlinks = True, pattern = None):
    """ Find all the regular files that can be found in the given paths.

    The method receives a variable number of paths and returns a list with all
    the existing regular files that were found at these locations. If a path
    corresponds to a regular file, it is simply added to the list, while if it
    points to a directory it is recursively walked top-down in search of
    regular files. In other words: if the path to a directory is given, all
    the regular files in the directory tree are included in the returned list.

    Keyword arguments:
    followlinks - by default, the method will walk down into symbolic links
                  that resolve to directories. You may set this to False to
                  disable visiting directories pointed to by symlinks. Note
                  that setting followlinks to True can lead to infinite
                  recursion if a link points to a parent directory of itself.
    pattern - the pattern, according to the rules used by the Unix shell (which
              are not the same as regular expressions) that the base name of a
              regular file must match to be considered when scanning the
              paths. Non-matching files are ignored.

    """

    files_paths = []
    for path in sorted(paths):
        if os.path.isfile(path):
            basename = os.path.basename(path)
            if not pattern or fnmatch.fnmatch(basename, pattern):
                files_paths.append(path)

        elif os.path.isdir(path):
            tree = os.walk(path, followlinks = followlinks)
            for dirpath, dirnames, filenames in tree:
                dirnames.sort()
                for basename in sorted(filenames):
                    abs_path = os.path.join(dirpath, basename)
                    files_paths += find_files([abs_path],
                                              followlinks = followlinks,
                                              pattern = pattern)
    return files_paths
