#! /usr/bin/env python
# encoding:UTF-8

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

import os
os.environ['PYRAF_NO_DISPLAY'] = '1'
import pyraf.iraf
from pyraf.iraf import noao, artdata  # 'noao.artdata' package

import calendar
import fnmatch
import hashlib
import logging
import os.path
import pyfits
import re
import shutil
import stat
import tempfile
import time
import warnings

# LEMON modules
import astromatic
import methods
import style

class NonFITSFile(TypeError):
    """ Raised when a non-FITS file is opened by the FITSImage class """
    pass

class NonStandardFITS(ValueError):
    """ Raised when a non-standard file is attempted to be opened."""
    pass

class FITSImage(object):
    """ Encapsulates a FITS image located in the filesystem. """

    def __init__(self, path):
        """ Instantiation method for the FITSImage class.

        The method creates an instance of a FITS image, whose header is cached
        into memory for fast access to its keywords. The NonFITSFile exception
        is raised in case 'path' is not a FITS file, while NonStardarFITS is
        raised if it is but does not conform to the FITS standard.

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

                # The SIMPLE keyword is required to be the first keyword in
                # the primary header of all FITS files. It contains a logical
                # constant with the value 'T' if the file conforms to this
                # standard. This keyword is mandatory for the primary header
                # and is not permitted in extension headers. A value of 'F'
                # signifies that the file does not conform to this standard.
                # [URL] http://archive.stsci.edu/fits/fits_standard/node39.html

                if not handler[0].header['SIMPLE']:
                    msg = "file '%s' does not conform to the FITS standard" % self.path
                    raise NonStandardFITS(msg)

                # A copy of the FITS header is kept in memory and the file is
                # closed; otherwise we may run into trouble when working with
                # thousands of images ("too many open files" and such). This
                # approach gives us fast read-only access to the image header;
                # if modified, we will have to take care of 'reloading' (call
                # it synchronize, if you wish) the header.

                self.size = handler[0].data.shape[::-1]
                self.header = handler[0].header
            finally:
                handler.close(output_verify = 'ignore')

        # PyFITS raises IOError if we attempt to open a non-FITS file. However,
        # catching the IndexError exception is also necessary, as it seems that
        # PyFITS does not throw the IOError exception when opening an empty
        # file, but instead waits until the header or data are accessed to
        # raise IndexError.

        except (IOError, IndexError, AttributeError):
            raise NonFITSFile("file '%s' is not a FITS file" % self.path)
        except NonStandardFITS:
            raise

    def __del__(self):
        pass

    def unlink(self):
        """ Remove the FITS Image to which the instance is linked.

        The method deletes the file to which the instance is currently linked.
        Be careful! It is not the instance that is deleted, but the _file_!
        This method should _only_ by called right before the instance deletion,
        as otherwise further operations may, understandably, fail.

        """

        os.unlink(self.path)
        self.path = None

    def __repr__(self):
        """ The 'official' string representation of the FITS image """
        return "%s('%s')" % (self.__class__.__name__, self.path)

    def __str__(self):
        """ The 'informal' string representation of the FITS image. """
        return repr(self)

    def __eq__(self, other):
        """ Test if both FITSImage instances encapsulate the same FITS image.

        The method returns True if both instances point to the same FITS image,
        and False otherwise. In order to take into account symbolic links, this
        is not checked by simply comparing the path of both instances, but by
        actually testing if both paths point to the same file on disk. If this
        happens, the inode numbers will be equal, although we must also check
        that the device the inodes reside on are the same, as two files on
        different devices might have the same inode number.

        """

        sinode = os.stat(self.path)[stat.ST_INO]
        oinode = os.stat(other.path)[stat.ST_INO]
        sdevice = os.stat(self.path)[stat.ST_DEV]
        odevice = os.stat(other.path)[stat.ST_DEV]
        return sinode == oinode and sdevice == odevice

    def read_keyword(self, keyword):
        """ Read a keyword from the header of the FITS image.

        Note that, although always upper-case in the FITS header, keywords are
        here case-insensitive, for user's convenience. TypeError is raised if
        the keyword is None, ValueError if it is left empty and KeyError if the
        keyword cannot be found in the header.

        """

        if keyword is None:
            raise TypeError("keyword cannot be None")
        if not keyword:
            raise ValueError("keyword cannot be empty")
        return self.header[keyword.upper()]

    def update_keyword(self, keyword, value, comment = None):
        """ Updates the value of a FITS keyword, adding it if it does not exist.

        The method updates the value of a keyword in the FITS header, replacing
        it with the specified value or simply adding it in case if does not yet
        exist. Note that, although always upper-case inside the FITS file,
        keywords are here case-insensitive, for user's convenience.

        Keyword arguments:
        comment - the comment to be added to the keyword.

        """

        # [From the PyFits Handbook, page 20] For keywords longer than 8
        # characters, there is a convention originated at ESO to facilitate
        # such use. It uses a special keyword HIERARCH with the actual long
        # keyword following. [...] When creating or updating using the
        # header.update() method, it is necessary to prepend 'hierarch'.

        if len(keyword) > 8:
            msg = "%s: keyword '%s' is longer than eight characters"
            logging.debug(msg % (self.path, keyword))
            hierarch_keyword = 'HIERARCH ' + keyword
            msg = "%s: keyword '%s' became '%s'"
            logging.debug(msg % (self.path, keyword, hierarch_keyword))
            keyword = hierarch_keyword

        handler = pyfits.open(self.path, mode = 'update')
        msg = "%s: file opened to update '%s' keyword" % (self.path, keyword)
        logging.debug(msg)

        try:
            header = handler[0].header

            # Ignore the 'card is too long, comment is truncated' warning
            # printed by PyRAF in case, well, the comment is too long.
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                header.update(keyword, value, comment = comment)
                args = self.path, keyword, value
                msg = "%s: keyword '%s' updated to '%s'" % args
                if comment:
                    msg += " with comment '%s'" % comment
                logging.debug(msg)

            # Update in-memory copy of the FITS header
            self.header = header

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
            del header[keyword]
            # Update in-memory copy of the FITS header
            self.header = header
        finally:
            handler.close(output_verify = 'ignore')

    def add_comment(self, comment):
        """ Add a descriptive comment to the header of the image.

        According to the FITS standard, the 'COMMENTS' keyword has no
        associated value and columns 9-80 may contain any ASCII text. Any
        number of COMMENT card images may appear in a header. Ironically,
        there is no comment in a commentary card, only a string value.
        This method does not update the in-memory copy of the header.

        """

        handler = pyfits.open(self.path, mode = 'update')
        try:
            header = handler[0].header
            header.add_comment(comment)
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

    def date(self, date_keyword = 'DATE-OBS', exp_keyword = 'EXPTIME'):
        """ Return the date of observation, in Unix time.

        This method returns, in seconds after the Unix epoch (the time 00:00:00
        UTC on 1 January 1970), the date of the 'midpoint' of the exposure;
        i.e., the time of the start of the observation of plus one-half the
        duration of the exposure.

        The KeyError exception is raised if any of the specified keywords
        cannot be found in the FITS header. NonStandardFITS is thrown if the
        keywords exist but the date they represent does not follow the FITS
        standard. See the URL below for more information on the standard.

        Definition of the Flexible Image Transport System (FITS):
        http://archive.stsci.edu/fits/fits_standard/

        Keyword arguments:
        date_kewyord - the FITS keyword in which the date of the observation is
                       stored, in the format specified in the FITS Standard. The
                       old date format was 'yy/mm/dd' and may be used only for
                       dates from 1900 through 1999. The new Y2K compliant date
                       format is 'yyyy-mm-dd' or 'yyyy-mm-ddTHH:MM:SS[.sss]'.
        exp_keyword - the FITS keyword in which the duration of the exposure is
                      stored. It is expected to be a floating-point number which
                      gives the duration in seconds. The exact definition of
                      'exposure time' is mission dependent and may, for example,
                      include corrections for shutter open and close duration,
                      detector dead time, vignetting, or other effects.

        """

        # Throws KeyError is the specified keyword is not found in the header
        obs_date_string = self.read_keyword(date_keyword)

        # Parse the string representation of the observation date according to
        # the format specified in the FITS Standard. The old date format was
        # 'yy/mm/dd' and may be used only for dates from 1900 through 1999.
        # The new Y2K compliant date format is 'yyyy-mm-dd' or
        # 'yyyy-mm-ddTHH:MM:SS[.sss]'.
        #
        # time.strptime raises ValueError if the string cannot be parsed
        # according to format, or if it has excess data after parsing.

        try:
            # yy/mm/dd; deprecated, old date format
            if re.match('^\d{2}/\d{2}/\d{2}$', obs_date_string):
                obs_date = time.strptime(obs_date_string, '%y-%m-%d')
            # yyyy-mm-dd
            elif re.match('^\d{4}-\d{2}-\d{2}$', obs_date_string):
                obs_date = time.strptime(obs_date_string, '%Y-%m-%d')

            # yyyy-mm-ddTHH:MM:SS[.sss] Note that we discard the decimals of
            # the seconds (by parsing the first parenthesized group instead of
            # the entire matching string) as time.strptime raises ValueError if
            # the string has excess data after the parsing
            #
            # [Fix Wed Feb 8 2012: allow up to four decimals in the seconds
            # [.ssss]. Although this is not standard, it is the information
            # used in the header of O2K (Calar Alto Observatory) images]
            matched = re.match("^(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})"
                               "(\.\d{0,4})?$", obs_date_string)
            if matched:
                obs_date = time.strptime(matched.group(1), '%Y-%m-%dT%H:%M:%S')
            else:
                raise NonStandardFITS

        except (TypeError, ValueError, NonStandardFITS):
            msg = "'%s' keyword does not follow the FITS standard" % date_keyword
            raise NonStandardFITS(msg)

        try:
            # Divide the duration of the exposure, in seconds, by two
            half_exp_time = self.read_keyword(exp_keyword) / 2
        except KeyError:    # if the keyword is not found
            raise
        except ValueError:  # if the casting to float fails
            msg = "'%s' keyword must be a real number" % exp_keyword
            raise NonStandardFITS(msg)

        # conversion from struct_time in UTC to seconds since the Unix epoch
        return calendar.timegm(obs_date) + half_exp_time

    @property
    def basename(self):
        """ Return the base name of the FITS image.

        The method returns the string that results from deleting from the path
        of the image any prefix up to the last slash. For example, for an image
        located in /caha/ferM_0013.fits, 'ferM_0013.fits' is returned.

        """

        return os.path.basename(self.path)

    @property
    def basename_woe(self):
        """ Return the base name of the FITS image, without extension.

        The method returns the basename of the image with its extension
        stripped. That is, all the characters of the base name starting from
        the rightmost dot are discarded. For example, for an image located in
        /caha/ferM_0013.fits, the string 'ferM_0013' is returned.

        """

        return os.path.splitext(self.basename)[0]

    @property
    def extension(self):
        """ Return the extension of the FITS image.

        The method returns the extension of the filename, that is, the
        characters that follow the rightmost dot ('.'), which is also
        included, in the filename. For example, for an image located in
        /caha/ferM_0013.fits, the string '.fits' is returned.
        An empty string is returned if the filename has no extension.

        """

        return os.path.splitext(self.path)[1]

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
        for character in self.basename_woe:
            try:
                int(character)
                return str_char # stop as soon as an integer is found
            except ValueError:
                # Raised because the character is not an integer
                str_char += character
        return str_char

    @property
    def number(self):
        """ Extract the number that follows the filename prefix.

        The method strips the filename of its prefix and, after that, loops
        over the resulting substring until a non-numeric character is found,
        moment at which the seen (and contiguous) numbers seen so far are
        returned as a string. For example, for the image /caha/ferM_0013.fits,
        the string '0013' would be returned.

        The reason why the number is returned as a string is to not remove
        the leading zeros, if any. The value returned by the method may, thus,
        need to be cast to integer.

        """

        # Strip the prefix, if any, from the filename
        f_prefix = self.prefix
        if f_prefix:
            basename_woprefix = self.basename[len(f_prefix):]
        else:
            basename_woprefix = self.basename

        characters = basename_woprefix
        for index in xrange(len(characters) - 1):
            if not characters[:index+1].isdigit():
                break

        # The variable 'index' will not be defined when the basename of the
        # image, without its prefix, is less than two characters long. If that
        # is the case, the for loop does not iterate over any value, and thus
        # the return statement in the try clause raises UnboundLocalError.

        try:
            return characters[:index]
        except UnboundLocalError:
            assert len(characters) < 2
            return ''

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

    def _subscript(self, x1, x2, y1, y2):
        """ Return the string representation of the image subscript.

        The method receives four coordinates and returns the path wich would
        have to be passed to IRAF in order only to work with a section of an
        image. For the image '/home/123m/ferM_0001.fits' and x1 = 1, x2 = 100,
        y1 = 1 and y2 = 200, the string '/home/123m/ferM_0001.fiys[1:100,1:200]'
        would be returned.

        """
        return "%s[%d:%d,%d:%d]" % (self.path, x1, x2, y1, y2)

    def _subscript_instance(self, x1, x2, y1, y2):
        """ Return a FITSImage instance of the subscripted version of the image.

        The method returns an instance of FITSImage that encapsulates the
        'subscripted' path of the image. Although still pointing to the same
        FITS image, these instances can be used to perform IRAF tasks on a
        section of the image, such as ferM_0001.fits[0:100,400:450].

        """

        simage = FITSImage(self.path)
        simage.path = simage._subscript(x1, x2, y1, y2)
        return simage

    def imstat(self, statistic, x1 = None, x2 = None, y1 = None, y2 = None):
        """ Compute image pixel statistics.

        The method invokes IRAF's 'imstat' on the image and returns the value
        computed by the task. The available statistics, according to the IRAF
        help page, are the following:

        image - the image name
        npix - the number of pixels used to do the statistics
        mean - the mean of the pixel distribution
        midpt - estimate of the median of the pixel distribution
        mode - the mode of the pixel distribution
        stddev - the standard deviation of the pixel distribution
        skew - the skew of the pixel distribution
        kurtosis - the kurtosis of the pixel distribution
        min - the minimum pixel value
        max - the maximum pixel value

        ValueError is raised if an invalid statistic is specified.

        By default, the statistics of the entire image are calculated. However,
        the results can be restricted to a section of the image if the four x1,
        x2, y1 and y2 coordinates are specified. The four coordinates must be
        given; otherwise, the entire image is used.

        Keyword arguments:
        x1 - the lower bound in the x-axis, if only a section is to be considered.
        x2 - the upper bound in the x-axis, if only a section is to be considered.
        y1 - the lower bound in the y-axis, if only a section is to be considered.
        y2 - the upper bound in the y-axis, if only a section is to be considered.

        """

        stats = ['image', 'npix', 'mean', 'midpt', 'mode',
                 'stddev', 'skew', 'kurtosis', 'min', 'max']

        if statistic.lower() not in stats:
            msg = "invalid IRAF's imstat statistic '%s'" % statistic
            raise ValueError(msg)

        # If a region has been specified, create a subscripted instance (which
        # does not duplicate the FITS image, but simply 'points' to the region
        # of the image in which we are interested) and call the method on it.
        if x1 is not None and \
           x2 is not None and \
           y1 is not None and \
           y2 is not None:
            sinstance = self._subscript_instance(x1, x2, y1, y2)
            return sinstance.imstat(statistic)

        with tempfile.TemporaryFile() as temp_file:

            pyraf.iraf.imstat(self.path, fields = statistic, Stdout = temp_file)
            temp_file.flush()
            temp_file.seek(0)

            # The contents of the file, as of IRAF v2.14, look like this:
            #  MEAN
            #  14217.
            #
            # Therefore, we need to read the value in the second line, casting
            # it to float. Note that IRAF's imstat returns nothing if the
            # specified field does not exist, a case in which the second line
            # of the file would be empty. This, however, can never happen, as
            # we previously made sure that the specified statistic is valid.

            stripped_line = temp_file.readlines()
            assert stripped_line
            field_value = stripped_line[1].strip()

            # All the fields except the image name must be parsed to float
            if statistic != 'image':
                field_value = float(field_value)

            return field_value

    def stable_overscan(self, x1, x2, y1, y2, threshold, offset, verbose = True):
        """ Determine the coordinates of the stable section of the overscan area.

        The method automatically finds the portion of the overscan section of
        the frame that is considered to be 'stable'. The initial area, that
        delimited by the coordinates [x1:x2, y1:y2], is iteratively shrunk by
        'offset' pixels until the percentage change between (1) the arithmetic
        mean and (2) the standard deviation of the pixels in the current area
        and those in the previous area are smaller or equal than the specified
        threshold. In other words, the stopping criterion of the overscan area
        reduction sequence is given as the maximum allowed percentage change
        between the arithmetic mean and the standard deviation of the pixels in
        the two last considered overscan areas.

        Upon convergence, the coordinates of the determined 'stable' overscan
        area of the frame are returned as a four-element tuple,
        (x1, x2, y1, y2).

        If the overscan area is exahusted (that is, completely shrunk in one or
        both axes) without convergence, the method is run again with twice a
        threshold value. That is, the detection process starts over with a
        'relaxed' stopping condition. This is repeated until the method
        converges and a stable region of the overscan area is found.

        theshold - the maximum allowed percentage change between the arithmetic
                   mean and standard deviation of the pixels in the current
                   overscan area and those in the previous one if they are
                   considered to have converged. Note that, e.g., 0.1 means 10%
        offset - the number of pixels by which each coordinate is increased (x1
                 and y1) or decreased (x2 and y2) every time the overscan area
                 is shrunk.
        verbose - if True, information will be printed to the standard output
                  as the 'stable' overscan area of the frame is determined.
                  Otherwise, the method will run silently.

        """

        # Calculate the arithmetic mean and standard deviation
        # of the initial overscan portion of the frame
        mean   = self.imstat('mean', x1, x2, y1, y2)
        stddev = self.imstat('stddev', x1, x2, y1, y2)

        if verbose:
            print "%s%s --> mean = %.3f, stddev = %.3f" % \
                  (style.prefix, self._subscript(x1, x2, y1, y2), mean, stddev)

        while True:

            # Shrink overscan area
            x1 += offset
            x2 -= offset
            y1 += offset
            y2 -= offset

            # Test that x1 and y1 are smaller or equal than y1 and y2, respectively.
            # Otherwise we have to stop, as the overscan area has been completely
            # shrunk, in one or both axes. In that case, we try to determine the
            # stable overscan area again, but this time with threshold * 2 (that
            # is, we 'relax' the stopping condition).

            if x1 > x2 or y1 > y2:

                if verbose:
                    if x1 > x2 and y1 > y2:
                        which_axis = "x- and y-axes"
                    elif x1 > x2:
                        which_axis = "x-axis"
                    else:
                        which_axis = "y-axis"

                    msg1 = "The overscan area was completely explored in " \
                           "the %s without convergence." % which_axis
                    msg2 = "It was not possible to determine the overscan " \
                           "area with threshold = %.3f." % threshold
                    msg3 = "Trying to determine the stable overscan area " \
                           "with threshold = %.3f..." % (threshold * 0.2)

                    warnings.warn(msg1, RuntimeWarning)
                    warnings.warn(msg2, RuntimeWarning)
                    warnings.warn(msg3, RuntimeWarning)

                else:
                    msg = "WARNING: the overscan area did not converge, " \
                          "trying with threshold = %.3f..." % (threshold * 0.2)
                    warnings.warn(msg, RuntimeWarning)

                return self.stable_overscan(x1, x2, y1, y2, threshold * 2,
                                            offset, verbose = verbose)

            # Now calculate the arithmetic mean and standard deviation of the
            # new section of the overscan area
            new_mean = self.imstat('mean', x1, x2, y1, y2)
            new_stddev = self.imstat('stddev', x1, x2, y1, y2)

            # We stop when the percentage change between the current and former
            # arithmetic means and standard deviations is small enough
            mean_pc = abs(methods.percentage_change(mean, new_mean))
            stddev_pc = abs(methods.percentage_change(stddev, new_stddev))

            if verbose:
                print "%s%s --> mean = %.3f, stddev = %.3f" % (style.prefix,
                      self._subscript(x1, x2, y1, y2), new_mean, new_stddev) ,
                print "(changes: %.6f and %.6f)" % (mean_pc, stddev_pc)

            if stddev_pc <= threshold and mean_pc <= threshold:
                if verbose:
                    print "%sThe overscan area of %s converged at coordinates " \
                          "[%d:%d,%d:%d]." % \
                           (style.prefix, self.path, x1, x2, y1, y2)

                # Return the coordinates of the stable area as a four-element tuple
                return x1, x2, y1, y2

            # Update the value of the current standard deviation and arithmetic
            # mean before moving to the next iteration
            stddev = new_stddev
            mean = new_mean

    def imcopy(self, output_path = None, suffix = '.fits',
               x1 = None, x2 = None, y1 = None, y2 = None):
        """ Return the FITSImage instance of a copy of the FITS image.

        The method makes a copy of the FITS image and returns the FITS instance
        that encapsulates it. If no output name is specified, the FITS image
        will be copied to a temporary file (most likely in the /tmp/ directory)
        for whose deletion when it is no longer needed the user is responsible.

        By default, a copy of the full image is done. However, the copy can be
        restricted to a section of the image if the four x1, x2, y1 and y2
        coordinates are specified.

        Keyword arguments:
        output_path - the path to which the copy of the FITS image will be
                      saved. If not specified, the image will be copied to a
                      temporary file.
        suffix - if no output path is specified, use this suffix when saving
                 the resulting image to a temporary file.
        x1 - the lower bound in the x-axis, if only a section is to be copied.
        x2 - the upper bound in the x-axis, if only a section is to be copied.
        y1 - the lower bound in the y-axis, if only a section is to be copied.
        y2 - the upper bound in the y-axis, if only a section is to be copied.

        """

        if output_path:
            copied_path = output_path
        else:
            # If no output path was given, copy the image to a temporary file.
            copy_fd, copied_path = tempfile.mkstemp(suffix = suffix)
            os.close(copy_fd)

        # If the method is called so that the image overwrites itself (i.e.,
        # if the output path matches the path to the FITS image to which the
        # instance is linked), make a temporary copy and work with it.

        if self.path == output_path:
            source_image = self.imcopy()
            image_copied = True
        else:
            source_image = self
            image_copied = False

        # The four coordinates must have been specified if only a portion of
        # the image is to be copied.
        if x1 is not None and \
           x2 is not None and \
           y1 is not None and \
           y2 is not None:
            subscript_path = source_image._subscript(x1, x2, y1, y2)
            if os.path.exists(copied_path):
                os.unlink(copied_path) # remove silently

            # 'Stdout' prevents IRAF from printing a message like this:
            # ferM_0001.fits[500:1000,300:1300] -> output_image.fits
            with open(os.devnull, 'wt') as fd:
                pyraf.iraf.imcopy(subscript_path, copied_path, Stdout = fd)

        else:
            shutil.copy2(source_image.path, copied_path)

        # Delete the temporary of the FITS image, if any
        if image_copied:
            os.unlink(source_image.path)
            del source_image

        return FITSImage(copied_path)

    def imarith(self, operator, other, output_path = None):
        """ Binary image arithmetic.

        The method invokes the IRAF task 'imarith' on the two operators, which
        can be FITSImage instances or numeric values (if, for example, the
        operation ferM_0001.fits * 1.25 were to be calculated), performing the
        specified operation and returning a new instance of FITSImage which
        encapsulates the output image. Note that, while IRAF does not overwrite
        files and instead aborts the execution, this method prefers to delete
        the existing file before proceeding, thus effectively overwriting it.

        'operator' must be a string with the operator to be applied to the
        operands. Must be one of +, -, *, /, min or max; otherwise, ValueError
        is raised.

        If no output path is specified, the resulting image will be saved to a
        temporary file, most likely in the directory /tmp.

        """

        # Test that the specified operand is a valid one
        if operator not in ['+', '-', '*', '/', 'min', 'max']:
            raise ValueError("illegal operand when invoking 'imarith'")

        # Determine the path to the output file, if none was specified
        if not output_path:
            output_fd, output_path = tempfile.mkstemp(suffix = '.fits')
            os.close(output_fd)

        if not isinstance(other, FITSImage) and \
           not isinstance(other, int) and not isinstance(other, float):
            msg = "second operand must be either a FITSImage or a real value"
            raise ValueError(msg)

        # If the specified output will overwrite one (or both, e.g.,
        # ferM_0001.fits = ferM_0001.fits * ferM0001.fits) of the input images,
        # make a temporary copy of those images to be overwritten

        self_copied, other_copied = False, False
        if self.path == output_path:
            first_operand = self.imcopy().path
            self_copied = True
        else:
            first_operand = self.path

        if isinstance(other, FITSImage):
            if other.path == output_path:
                second_operand = other.imcopy().path
                other_copied = True
            else:
                second_operand = other.path

        else:
            # Either a float or an integer
            second_operand = other

       # If the output file already exists, delete it silently
        if os.path.exists(output_path):
            os.unlink(output_path)

        pyraf.iraf.imarith(first_operand, operator, second_operand, output_path)

        # If temporary copies of the FITS images were made, delete them
        if self_copied:
            os.unlink(first_operand)
        if other_copied:
            os.unlink(second_operand)

        return FITSImage(output_path)

    def normalize(self, path = None, fraction = 0.5, method = 'mean'):
        """ Normalize the image, by dividing it by its central value.

        This method, originally written for normalizing flat-field images,
        calculates the central value of the image (the median by default,
        although other parameters accepted by imstat, such as 'mean', may be
        used) and divides all the pixels by the result. The dimensions of the
        central area are defined in terms of the fraction of the smallest axis:
        a square of side fraction * min(x,y) will be centered in the image and
        those pixels that fall within it used to compute the central value used
        for the normalization.

        In case the value of 'fraction' causes the square area to be bigger than
        the image itself, the method defaults back to the full frame dimensions,
        that is, the entire image will be used for the normalization. Note that
        this may happen for values smaller than (but very close to) one, as for
        images with sides of an odd number of pixels the coordinates of its
        center are rounded to the nearest integer.

        If no output path is given, the normalized image is saved to a temporary
        file created in the most secure manner possible -- readable and writable
        only by the creating user ID. The method returns a FITSImage instance
        which encapsulates the resulting image.

        """

        # Determine the coordinates of the central area
        x_size, y_size = self.size
        x_center, y_center = self.center
        square_width = min(x_size, y_size)
        half_width = int((square_width * fraction) / 2)
        x1 = x_center - half_width
        x2 = x_center + half_width
        y1 = y_center - half_width
        y2 = y_center + half_width

        if min(x1, y1) < 1 or x2 > self.x_size or y2 > self.y_size:
            msg = "exceeded frame dimensions, using full frame instead"
            warnings.warn(msg, RuntimeWarning)
            x1 = y1 = 1
            x2 = self.x_size
            y2 = self.y_size

        subscripted_img = self._subscript_instance(x1, x2, y1, y2)
        central_value = subscripted_img.imstat(method)

        if not path:
            prefix = '%s_normalized_%s_' % (os.path.basename(self.path), fraction)
            fd, path = tempfile.mkstemp(prefix = prefix, suffix = '.fits')
            os.close(fd)

        return self.imarith('/', central_value, output_path = path)

    def ahead_file(self, path_to_ahead_file, scale, equinox, radecsys,
                   ra_keyword = 'RA', dec_keyword = 'DEC'):
        """ Generate the .ahead file needed by SCAMP in order to do astrometry.

        The method creates the .ahead file that ensures that the required FITS
        keywords, which define an initial guess of the astrometic solution,
        needed by Emmanuel Bertin's SCAMP are present. These keywords could
        also be direclty added or updated in the FITS image, but using the
        .ahead file allows us not to modify the original file. On successful
        completion, the method returns the path to the generated .ahead file.

        [From the SCAMP user guide, page 4] 'The binary catalogues in 'FITS
        LDAC' format read by SCAMP contain a copy of the original FITS image
        headers. These headers provide fundamental information such as frame
        dimensions, World Coordinate System (WCS) data and many other FITS
        keywords which SCAMP uses to derive a full astrometric and photometric
        calibration. It is often needed to change or add keywords in some
        headers. Editing FITS files is not convenient, so SCAMP provides read
        (and write) support for 'external' header files. External headers may
        either be real FITS header cards (no carriage-return), or ASCII files
        containing lines in FITS-like format, with the final line starting with
        'END⊔⊔⊔⊔⊔'. Multiple extensions must be separated by an 'END⊔⊔⊔⊔⊔'
        line. External “headers” need not contain all the FITS keywords
        normally required. The keywords present in external headers are only
        there to override their counterparts in the original image headers or
        to add new ones. Hence for every input (say, xxxx.cat) FITS catalogue,
        SCAMP looks for a xxxx.ahead header file, loads it if present, and
        overrides or adds to image header keywords those found there.'

        For further information on the World Coordinate System (WCS), refer to
        http://tdc-www.harvard.edu/software/wcstools/wcstools.wcs.html

        path_to_ahead_file - path to wich the ahead file will be saved.
        scale - scale of the image, in degrees per pixel
        equinox - equinox in years (e.g., 2000)
        radecsys - reference system (e.g., ICRS)

        Keyword arguments:
        ra_keyword - FITS keyword for the right ascension, in decimal degrees.
        dec_keyword - FITS keyword for the declination, in decimal degrees.

        """

        with open(path_to_ahead_file, 'wt') as fd:

            # Pixel coordinates of the reference point
            xcenter, ycenter = self.center
            fd.write("CRPIX1  = %d\n" % xcenter)
            fd.write("CRPIX2  = %d\n" % ycenter)

            ra  = self.read_keyword(ra_keyword)
            dec = self.read_keyword(dec_keyword)

            fd.write("CRVAL1  = %s\n" % ra)    # RA at reference point
            fd.write("CRVAL2  = %s\n" % dec)   # DEC at reference point
            fd.write("CTYPE1  = 'RA---TAN'\n") # Gnomonic projection
            fd.write("CTYPE2  = 'DEC--TAN'\n")

            # The FITS WCS standard uses a rotation matrix, CD1_1, CD1_2,
            # CD2_1, and CD2_2 to indicate both rotation and scale, allowing a
            # more intuitive computation if the axes are skewed. This model has
            # been used by HST and IRAF for several years.

            fd.write("CD1_1   = %f\n" % scale)
            fd.write("CD1_2   = 0.0\n")
            fd.write("CD2_1   = 0.0\n")
            fd.write("CD2_2   = %f\n" % -scale)
            fd.write("EQUINOX = %d\n" % equinox)  # equinox in years
            fd.write("RADECSYS= '%s'\n" % radecsys)  # reference system

            fd.write("END⊔⊔⊔⊔⊔")
            return path_to_ahead_file

    def add_margin(self, left_margin, right_margin, bottom_margin, top_margin):
        """ Add a blank margin to the FITS image.

        The method adds zero-valued margins of the given width, given in
        pixels, to the image. Internally, this is achieved by creating a
        larger, blank canvas with IRAF's mkpattern and then replacing the
        section of the cavas within those margins with the original image.
        Returns a FITSImage instance that encapsulates the resulting image,
        saved to a temporary file and for whose deletion when done with it
        the user is responsible .

        """

        # No duck typing here
        if not isinstance(left_margin, int) or \
           not isinstance(right_margin, int) or \
           not isinstance(bottom_margin, int) or \
           not isinstance(top_margin, int):
            raise ValueError("width of margings must be an integer")

        # Get path to the output temporary file
        output_fd, output_path = tempfile.mkstemp(suffix = '.fits')
        os.close(output_fd)

        # Determine the dimensions of the canvas, and therefore those of the
        # resulting image, by adding the width of the margings in each axis to
        # the dimensions of the original image.
        x_size, y_size = self.size
        dest_x_size = x_size + left_margin + right_margin
        dest_y_size = y_size + bottom_margin + top_margin

        os.unlink(output_path) # IRAF will refuse to overwrite it
        artdata.mkpattern(output_path, pattern = 'constant', v1 = 0,
                          ncols = dest_x_size, nlines = dest_y_size,
                          header = self.path) # also copy FITS header

        output_image = FITSImage(output_path)
        assert dest_x_size == output_image.x_size
        assert dest_y_size == output_image.y_size

        # Finally, copy the original image at the corresponding position in
        # the blank canvas. IRAF uses one indexing, so we need to transform
        # our zero-indexed coordinates.
        x1 = left_margin + 1
        x2 = x1 + x_size - 1
        y1 = bottom_margin + 1
        y2 = y1 + y_size - 1

        assert (x2 - x1 + 1) == x_size
        assert (y2 - y1 + 1) == y_size

        with open(os.devnull, 'wt') as fd:
            dest_path = output_image._subscript(x1, x2, y1, y2)
            pyraf.iraf.imcopy(self.path, dest_path, Stdout = fd)

        return output_image

    def imshift(self, xshift, yshift, interp_type = 'linear'):
        """ Shift an image in the x- and y-axes.

        This method, a high-level wrapper around IRAF's imshift, shifts an
        image in x and y such that xout = xin + xshift and yout = yin + yshift.
        The output image gray levels are determined by interpolating in the
        input image at the positions of the shifted output pixels.

        Returns a FITSImage instance that encapsulates the resulting image,
        saved to a temporary file and for whose deletion when done with it the
        user is responsible.

        Keyword arguments:
        interp_type - the interpolant type used to compute the output shifted
                      image. Defaults to 'linear' (bilinear interpolation in
                      x and y, although any of the other values accepted by
                      IRAF's imshift may be used.

        """

        output_fd, output_path = tempfile.mkstemp(suffix = '.fits')
        os.close(output_fd)

        os.unlink(output_path) # IRAF will refuse to overwrite it
        pyraf.iraf.imshift(self.path, output_path, xshift, yshift,
                           interp_type = interp_type,
                           boundary_type = 'constant')

        return FITSImage(output_path)

    def check_image(self, check_type = 'OBJECTS'):
        """ Return one of the multiple SExtractor's check-image.

        The method returns a FITSImage instance with the SExtractor check-image
        generated for this image. The check-image is saved to a temporary file,
        for whose deletion when done with it the user is responsible .

        The types of check-image currently accepted by SExtractor are:
        1. IDENTICAL: identical to input image (useful for converting formats)
        2. BACKGROUND: full-resolution interpolated background map
        3. BACKGROUND_RMS: full-resolution interpolated background noise map
        4. MINIBACKGROUND: low-resolution background map
        5. MINIBACK_RMS: low-resolution background noise map
        6. -BACKGROUND: background-subtracted image
        7. FILTERED: background-subtracted filtered image (requires FILTER = Y)
        8. OBJECTS: detected objects
        9. -OBJECTS: background-subtracted image with detected objects blanked
        10. APERTURES: MAG_APER and MAG_AUTO integration limits.
        11. SEGMENTATION: display patches corresponding to pixels attributed to
                          each object.

        """

        prefix = '%s_checkimage_%s_' % (self.basename, check_type.lower())
        check_fd, check_path = tempfile.mkstemp(prefix = prefix, suffix = '.fits')
        os.close(check_fd)

        try:

            sextractor_options = ['-CHECKIMAGE_NAME', check_path,
                                  '-CHECKIMAGE_TYPE', check_type]

            # Redirect SExtractor output to the null device
            with open(os.devnull, 'wt') as fd:
                retcode, catalog_path = \
                astromatic.sextractor(self.path, options = sextractor_options,
                                      stdout = fd, stderr = fd)
            return FITSImage(check_path)

        finally:
            try: os.unlink(catalog_path)
            except NameError: pass

    @property
    def sha1sum(self):
        """ Return the hexadecimal SHA-1 checksum of the FITS image """

        sha1 = hashlib.sha1()
        with open(self.path, 'rb') as fd:
            for line in fd:
                sha1.update(line)
            return sha1.hexdigest()


class FITSet(list):
    """ Encapsulates a set of FITS images (although it subclasses 'list') """

    def __init__(self, iterable = None):
        """ Instantiation method for the FITSet class.

        Return a new instance of FITSet whose elements are taken from iterable.
        Those elements in iterable which are not FITSImages are silently
        ignored. If iterable is not specified, a new empty set is returned.

        Keyword arguments:
        iterable - from which the initial FITSImages of the set are taken.

        """

        if iterable is None:
            iterable = []

        super(list, self).__init__()
        for img in iterable:
            if isinstance(img, FITSImage):
                self.append(img)

    def imcombine(self, output_path, flow, fhigh):
        """ Mim-max average combine all the images in the set.

        The method invokes the IRAF task 'imcombine' on all the FITS images
        contained in the set. The type of combining operation performed on the
        pixels (the 'combine' parameter when 'imcombine' is invoked) is
        'average', while the rejection operation (the 'reject' parameter) is
        'minmax'. That is the reason why 'flow' and 'fhigh', which determine
        the fraction of low and high pixels, respectively, that are to be
        rejected, are required as arguments.

        If existing, the output image is silently overwritten. The method
        returns a FITSImage instance which encapsulated the output image. The
        ValueError is raised if the set is empty or if not all the images in
        the set have the same dimensions.

        """

        if not len(self):
            raise ValueError("FITSet is empty")

        if not 0 <= flow <= 1 or not 0 <= fhigh <= 1:
            raise ValueError("both 'nlow' and 'fhigh' must be in the range [0,1]")

        # Make sure that all the images are of the same size
        if len(set(img.size for img in self)) != 1:
            raise ValueError("images in the FITSet differ in size")

        # Number of low and high pixels are to reject; round to nearest integer
        nlow  = int(round(len(self) * flow))
        nhigh = int(round(len(self) * fhigh))

        # If only two images are combined, we always discard the highest pixel
        # (since some stars could be discernible in the case of flat-fields,
        # for example). When there are three or more stars we also discard the
        # lowest and highest pixels, unless the specified fraction is zero: in
        # that case, it is understood as to mean "do not discard anything".
        if len(self) == 2 and not nhigh and fhigh:
            nhigh = 1
        elif len(self) >= 3:
            if not nlow and flow:
                nlow = 1
            if not nhigh and fhigh:
                nhigh = 1

        # There is what seems to be a bug in PyRAF that makes it impossible to
        # imcombine more than approximately forty images, so we save the list
        # of input images to a temporary file, use it as the input to IRAF's
        # imcombine and delete it before leaving the method.

        input_fd, input_path = tempfile.mkstemp(suffix = '.lst', text = True)

        try:
            for img in self:
                os.write(input_fd, '%s\n' % img.path)
            os.close(input_fd)

            if os.path.exists(output_path):
                os.unlink(output_path)

            with open(os.devnull, 'wt') as fd:
                pyraf.iraf.imcombine('@' + input_path, output_path,
                                     combine = 'average', reject = 'minmax',
                                     nlow = nlow, nhigh = nhigh,
                                     Stdout = fd)
        finally:
            os.unlink(input_path)

        # Delete the keywords that store the cached SExtractor catalog and the
        # hash of the configuration files, if present. Otherwise, it would mean
        # that the catalog of one of the input images (most likely the first in
        # the FITSet, and thus the first passed to IRAF's imcombine) would be
        # that of the resulting image, which would be obviously wrong! See
        # seeing.FITSeeingImage.__init__ for further information on these two
        # SExtractor-related keywords.

        output_img = FITSImage(output_path)
        output_img.delete_keyword('SEX_CATALOG')
        output_img.delete_keyword('SEX_MD5SUM')
        return output_img

    def sort(self):
        """ Sort the images in the set by their number.

        The method returns a new instance of FITSet that only contains those
        images that have a number, as returned by FITSImage.number, criterion
        by which they are sorted in increasing order. In other words the method
        returns a set whose images are sorted by their number, and from which
        those without one have been excluded.

        """

        imgs = [img for img in self if len(img.number)]
        imgs.sort(key = lambda x: int(x.number))
        return FITSet(imgs)

    def number_interval_filter(self, interval_lower, interval_upper):
        """ Return those images whose number is within the specified interval.

        The method returns a new instance of FITSet that contains only those
        images in the set whose number, as returned by FITSImage.number,
        belongs to the interval [interval_lower, interval_upper]. Images
        are sorted in increasing order in the returned FITSet.

        For example: if the set contains the images './ultra2/ferM_0098.fits',
        './ultra2/ferM_0099.fits', './ultra2/ferM_0100.fits',
        './ultra2/ferM_0101.fits' and './ultra2/ferM_0102.fits', interval_lower
        = 70 and interval_upper = 100, the returned FITSet would contain the
        images './ultra2/ferM_0098.fits', './ultra2/ferM_0099.fits' and
        './ultra2/ferM_0098.fits'], that is, those whose number x is
        interval_lower <= x <= interval_upper.

        """

        selected = FITSet()
        for image in self:

            # Ignore image if it has no number
            if image.number:
                if interval_lower <= int(image.number) <= interval_upper:
                    selected.append(image)
            else:
                msg = "%s has no number, ignored" % image.path
                warnings.warn(msg, RuntimeWarning)

        return selected.sort()

    def date_sort(self, date_keyword = 'DATE-OBS', exp_keyword = 'EXPTIME'):
        """ Sort the images in the set by their observation date.

        The method returns a new instance of FITSet, which images are sorted
        into ascending order by the date in which they were taken. In other
        words, the order of the FITS images in the returned set is determined
        by their observation date, as returned by FITSImage.date. Images for
        which one of the keywords cannot be found in the FITS header of that
        do not follow the FITS standard are excluded from the returned set.

        Keyword arguments:
        date_kewyord - the FITS keyword in which the date of the observation is
                       stored, in the format specified in the FITS Standard. The
                       old date format was 'yy/mm/dd' and may be used only for
                       dates from 1900 through 1999. The new Y2K compliant date
                       format is 'yyyy-mm-dd' or 'yyyy-mm-ddTHH:MM:SS[.sss]'.
        exp_keyword - the FITS keyword in which the duration of the exposure is
                      stored. It is expected to be a floating-point number which
                      gives the duration in seconds. The exact definition of
                      'exposure time' is mission dependent and may, for example,
                      include corrections for shutter open and close duration,
                      detector dead time, vignetting, or other effects.

        """

        # Store the date (as Unix time) of observation of each image in a
        # dictionary, which is used as a cache of dates. It is also here when
        # we discard those images for which the date cannot be calculated:
        # KeyError is raised if any of the keywords is not in the header, while
        # NonStandardFITS is raised if they do not conform to the FITS
        # specification.

        dates_cache = {}
        for img in self:
            try:
                obs_date = img.date(date_keyword = date_keyword,
                                    exp_keyword = exp_keyword)
                dates_cache[img] = obs_date
            except (KeyError, NonStandardFITS):
                pass

        imgs = sorted(dates_cache.iterkeys(), key = lambda x: dates_cache[x])
        return FITSet(imgs)

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

