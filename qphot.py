#! /usr/bin/env python

# Copyright (c) 2012 Victor Terron. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of LEMON.
#
# LEMON is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
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

from __future__ import division

"""
This module implements a wrapper for the IRAF tasks 'qphot' (quick aperture
photometer) and 'txdump' (which print fields from selected records in an
APPHOT/DAOPHOT text database). The routines implemented in this module provide
a way to automatically do photometry on an image returning the records in a
qphot.QPhot instance. Temporary files are used for both qphot and txdump, but
they are automatically removed, so the entire process takes place in memory,
from the user's perspective.

"""

import collections
import itertools
import logging
import math
import os
import os.path
import sys
import tempfile
import warnings

# LEMON modules
import fitsimage
import methods

# Tell PyRAF to skip all graphics initialization and run in terminal-only mode.
# Otherwise we will get annoying warning messages (such as "could not open
# XWindow display" or "No graphics display available for this session") when
# working at a remote terminal or at a terminal without any X Windows support.
# Any tasks which attempt to display graphics will fail, of course, but we are
# not going to make use of any of them, anyway.

os.environ['PYRAF_NO_DISPLAY'] = '1'

# Avoid PyRAF "Warning: no login.cl found" messages and do not clutter the
# filesystem with pyraf/ cache directories. The current working directory is
# also temporarily changed in fitsimage.py, for the same reason, so refer to
# that module for a more detailed explanation.

with methods.tmp_chdir(os.path.dirname(os.path.abspath(__file__))):
    import pyraf.iraf
    from pyraf.iraf import digiphot, apphot  # 'digiphot.apphot' package

# Turn PyRAF process caching off; otherwise, if we spawn multiple processes
# and run them in parallel, each one of them would use the same IRAF running
# executable, which could sometimes result in the most arcane of errors.
pyraf.iraf.prcacheOff()

# Decorate pyraf.subproc.Subprocess.__del__() to catch the SubprocessError
# exception that it occasionally raises (when the process is not gone after
# sending the TERM and KILL signals to it) and log it with level DEBUG on the
# root logger. As explained in the Python data model, uncaught exceptions in
# __del__() are ignored and a warning message, which we want to get rid of,
# such as the following, printed to the standard error instead:
#
# Exception pyraf.subproc.SubprocessError: SubprocessError("Failed
# kill of subproc 24600, '/iraf/iraf/bin.linux/x_images.e -c', with
# signals ['TERM', 'KILL']",) in <bound method Subprocess.__del__
# of <Subprocess '/iraf/iraf/bin.linux/x_images.e -c', at
# 7f9f3f408710>> ignored

func = methods.log_uncaught_exceptions(pyraf.subproc.Subprocess.__del__)
pyraf.subproc.Subprocess.__del__ = func

class MissingFITSKeyword(RuntimeWarning):
    """ Warning about keywords that cannot be read from a header (non-fatal) """
    pass

typename = 'QPhotResult'
field_names = "x, y, mag, sum, flux, stdev"
class QPhotResult(collections.namedtuple(typename, field_names)):
    """ Encapsulate the photometry of an astronomical object. In other words,
    each one of the lines that for each object are output by IRAF's qphot are
    parsed and saved as an object of this class.

    x, y - the x- and y-coordinates of the center of the object.
    mag - the instrumental magnitude of the object in the aperture.
    sum - the total number of counts in the aperture *including* the sky.
    flux - the total number of counts in the aperture *excluding* the sky.
    stdev - the standard deviation of the best estimate of the sky value,
            per pixel.

    """

    def snr(self, gain):
        """ Return the signal-to-noise ratio of the photometric measurement.

        The method returns the S/N, a quantitative measurement of the accuracy
        with which the object was observed. The signal-to-noise ratio tells us
        the relative size of the desired signal to the underlying noise or
        background light. The noise is defined as the standard deviation of a
        single measurement from the mean of the measurements made on an object.

        For photon arrivals, the statistical noise fluctuation is represented
        by the Poisson distribution. For bright sources where the sky
        background is negligible, S/N = total counts / sqrt(total counts).
        When the sky is not insignificant, the formula becomes S/N = (total
        counts - sky counts) / sqrt(total counts).

        Astronomers typically consider a good photoelectric measurement as one
        that has a signal-to-noise ratio of 100, or in other words, the noise
        is one percent of the signal. This implies an observational error of
        0.01 magnitude. [Henden, Astronomical Photometry, 1982]

        The 'gain' parameter gives the gain of the CCD in e-/ADU, as the above
        formulas apply to the number of electrons. Using the ADUs instead would
        introduce an error equal to the square root of the actual gain. In case
        the gain is unknown you may use a value of one (so as many electrons as
        ADUs will be used in the formula), as long as you understand that,
        although one is a common gain among many instruments, results may be
        only approximate.

        As counterintuitive as it seems, IRAF's qphot may return negative
        values of 'sum': e.g., if one of the steps of the data calibration
        (such as the bias subtraction) resulted in pixels with negative values.
        When that is the case, this method returns a negative SNR, which should
        be viewed as a red flag that something went wrong with this photometric
        measurement.

        """

        if gain <= 0:
            raise ValueError("CCD gain must be a positive value")

        # Division by zero must be avoided, as sum and flux may both be zero if
        # photometry is done on celestial coordinates that do not correspond to
        # any astronomical object, or if it is so faint that it is not visible.
        if not self.sum:
            return 0.0

        elif self.sum < 0.0:
            return -(abs(self.flux * gain) / math.sqrt(abs(self.sum * gain)))

        else:
            return (self.flux * gain) / math.sqrt(self.sum * gain)


class QPhot(list):
    """ The photometry of an image, as returned by IRAF's qphot.

    This class stores the result of the photometry done by IRAF's qphot (quick
    aperture photometer) on an image. A QPhotResult object is created for each
    object listed in the text file: after calling QPhot.run(), this subclass of
    the built-in list contains the photometric measurement of each astronomical
    object. The order of these QPhotResult objects is guaranteed to respect
    that in which coordinates are listed in the text file. In other words: the
    i-th QPhotResult object corresponds to the i-th astronomical object.

    """

    def __init__(self, img_path, coords_path):
        """ Instantiation method for the QPhot class.

        img_path - path to the FITS image on which to do photometry.
        coords_path - path to the text file with the celestial coordinates
                      (right ascension and declination) of the astronomical
                      objects to be measured. These objects must be listed
                      one per line, in two columns.

        """

        super(list, self).__init__()
        self.image = fitsimage.FITSImage(img_path)
        self.coords_path = coords_path

        for ra, dec in methods.load_coordinates(self.coords_path):
            if ra == 0 and dec == 0:
                msg = (
                  "the right ascension and declination of one or more "
                  "astronomical objects in '%s' is zero. This is a very bad "
                  "sign: these are the celestial coordinates that SExtractor "
                  "uses for sources detected on a FITS image that has not been "
                  "calibrated astrometrically (may that be your case?), and "
                  "without that it is impossible to do photometry on the "
                  "desired coordinates" % self.coords_path)
                # Emit warning only once
                warnings.warn(msg)
                break

    @property
    def path(self):
        """ Return the path to the FITS image. """
        return self.image.path

    def clear(self):
        """ Remove all the photometric measurements. """
        del self[:]

    def run(self, annulus, dannulus, aperture, exptimek):
        """ Run IRAF's qphot on the FITS image.

        This method is a wrapper, equivalent to (1) running 'qphot' on a FITS
        image and (2) using 'txdump' in order to extract some fields from the
        resulting text database. This subroutine, therefore, allows to easily
        do photometry on a FITS image, storing as a sequence of QPhotResult
        objects the photometric measurements. The method returns the number
        of astronomical objects on which photometry has been done.

        No matter how convenient, QPhot.run() should still be seen as a
        low-level method: INDEF objects will have a magnitude and standard
        deviation of None, but it does not pay attention to the saturation
        levels -- the sole reason for this being that IRAF's qphot, per se,
        provides no way of knowing if one or more pixels in the aperture are
        above some saturation level. For this very reason, the recommended
        approach to doing photometry is to use photometry(), a convenience
        function defined below.

        In the first step, photometry is done, using qphot, on the astronomical
        objects whose coordinates have been listed, one per line, in the text
        file passed as an argument to QPhot.__init__(). The output of this IRAF
        task is saved to temporary file (a), an APPHOT text database from which
        'txdump' extracts the fields to another temporary file, (b). Then this
        file is parsed, the information of each of its lines, one per object,
        used in order to create a QPhotResult object. All previous photometric
        measurements are lost every time this method is run. All the temporary
        files (a, b) are guaranteed to be deleted on exit, even if an error is
        encountered.

        An important note: you may find extremely confusing that, although the
        input that this method accepts are celestial coordinates (the right
        ascension and declination of the astronomical objects to be measured),
        the resulting QPhotResult objects store the x- and y-coordinates of
        their centers. The reason for this is that IRAF's qphot supports
        'world' coordinates as the system of the input coordinates, but the
        output coordinate system options are 'logical', 'tv', and 'physical':
        http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?qphot

        The x- and y-coordinates of the centers of the astronomical objects may
        be reported as -1 if their celestial coordinates fall considerably off
        those of the FITS image on which photometry is done. This is caused by
        a bug in IRAF, which, as of v.2.16.1, outputs invalid floating-point
        numbers (such as '-299866.375-58') when the astronomical object is
        approximately >= 100 degrees off the image boundaries. In those cases,
        the fallback value of -1 does not give us any other information than
        that the object falls off the image. This must probably be one of the
        least important bugs ever, considering how wrong the input celestial
        coordinates must be.

        Bug report that we have submitted to the IRAF development team:
        [URL] http://iraf.net/forum/viewtopic.php?showtopic=1468373

        Arguments:
        annulus - the inner radius of the sky annulus, in pixels.
        dannulus - the width of the sky annulus, in pixels.
        aperture - the aperture radius, in pixels.
        exptimek - the image header keyword containing the exposure time. Needed
                   by qphot in order to normalize the computed magnitudes to an
                   exposure time of one time unit.

        """

        self.clear() # empty object

        try:
            # Temporary file to which the APPHOT text database produced by
            # qphot will be saved. Even if empty, it must be deleted before
            # calling qphot. Otherwise, an error message, stating that the
            # operation "would overwrite existing file", will be thrown.
            output_fd, qphot_output = \
                tempfile.mkstemp(prefix = os.path.basename(self.path),
                                 suffix = '.qphot_output', text = True)
            os.close(output_fd)
            os.unlink(qphot_output)

            # In case the FITS image does not contain the keyword for the
            # exposure time, qphot() shows a message such as: "Warning: Image
            # NGC2264-mosaic.fits Keyword:  EXPTIME not found". This, however,
            # is not a warning, but a simple message written to standard error.
            # Filter this stream so that, if this message is written, it is
            # captured and issued as a MissingFITSkeyword warning instead.

            # Note the two whitespaces before 'Keyword'
            regexp = ("Warning: Image (?P<msg>{0}  Keyword: {1} "
                      "not found)".format(self.path, exptimek))

            args = sys.stderr, regexp, MissingFITSKeyword
            stderr = methods.StreamToWarningFilter(*args)

            # Run qphot on the image and save the output to our temporary
            # file. Note that cbox *must* be set to zero, as otherwise qphot
            # will compute accurate centers for each object using the centroid
            # centering algorithm. This is generally a good thing, but in our
            # case we want the photometry to be done exactly on the specified
            # coordinates.

            kwargs = dict(cbox = 0, annulus = annulus, dannulus = dannulus,
                          aperture = aperture, coords = self.coords_path,
                          output = qphot_output, exposure = exptimek,
                          wcsin = 'world', interactive = 'no',
                          Stderr = stderr)

            apphot.qphot(self.path, **kwargs)

            # Make sure the output was written to where we said
            assert os.path.exists(qphot_output)

            # Now extract the records from the APPHOT text database. We need
            # the path of another temporary file to which to save them. Even
            # if empty, we need to delete the temporary file created by
            # mkstemp(), as IRAF will not overwrite it.
            txdump_fd, txdump_output = \
                tempfile.mkstemp(prefix = os.path.basename(self.path),
                                 suffix ='.qphot_txdump', text = True)
            os.close(txdump_fd)
            os.unlink(txdump_output)

            # The type casting of Stdout to string is needed as txdump will not
            # work with unicode, if we happen to come across it: PyRAF requires
            # that redirection be to a file handle or string.

            txdump_fields = ['xcenter', 'ycenter', 'mag', 'sum', 'flux', 'stdev']
            pyraf.iraf.txdump(qphot_output, fields = ','.join(txdump_fields),
                              Stdout = str(txdump_output), expr = 'yes')

            # Now open the output file again and parse the output of txdump,
            # creating a QPhotResult object for each record.
            with open(txdump_output, 'rt') as txdump_fd:
                for line in txdump_fd:

                    fields = line.split()

                    # As of IRAF v.2.16.1, the qphot task may output an invalid
                    # floating-point number (such as "-299866.375-58") when the
                    # coordinates of the object to be measured fall considerably
                    # off (>= ~100 degrees) the image. That raises ValueError
                    # ("invalid literal for float()") when we attempt to convert
                    # the output xcenter or ycenter to float. In those cases, we
                    # use -1 as the fallback value.

                    try:
                        xcenter_str = fields[0]
                        xcenter     = float(xcenter_str)
                        msg = "%s: xcenter = %.8f" % (self.path, xcenter)
                        logging.debug(msg)
                    except ValueError, e:
                        msg = "%s: can't convert xcenter = '%s' to float (%s)"
                        logging.debug(msg % (self.path, xcenter_str, str(e)))
                        msg = "%s: xcenter set to -1 (fallback value)"
                        logging.debug(msg % self.path)
                        xcenter = -1

                    try:
                        ycenter_str = fields[1]
                        ycenter     = float(ycenter_str)
                        msg = "%s: ycenter = %.8f" % (self.path, ycenter)
                        logging.debug(msg)
                    except ValueError, e:
                        msg = "%s: can't convert ycenter = '%s' to float (%s)"
                        logging.debug(msg % (self.path, ycenter_str, str(e)))
                        msg = "%s: ycenter set to -1 (fallback value)"
                        logging.debug(msg % self.path)
                        ycenter = -1

                    try:
                        mag_str = fields[2]
                        mag     = float(mag_str)
                    except ValueError:  # float("INDEF")
                        assert mag_str == 'INDEF'
                        mag = None

                    sum_ = float(fields[3])
                    flux = float(fields[4])

                    try:
                        stdev_str = fields[5]
                        stdev = float(stdev_str)
                    except ValueError:  # float("INDEF")
                        assert stdev_str == 'INDEF'
                        stdev = None

                    args = xcenter, ycenter, mag, sum_, flux, stdev
                    self.append(QPhotResult(*args))

        finally:

            # Remove temporary files. The try-except is necessary because the
            # deletion may fail (OSError), or something could go wrong and an
            # exception be raised before the variables are defined (NameError).

            try:
                os.unlink(qphot_output)
            except (NameError, OSError):
                pass

            try:
                os.unlink(txdump_output)
            except (NameError, OSError):
                pass

        return len(self)


def run(img, coords_path,
        aperture, annulus, dannulus,
        maximum, exptimek, uncimgk):
    """ Do photometry on a FITS image.

    This convenience function does photometry on a FITSImage object, measuring
    the astronomical objects listed in the 'coords_path' text file. Note that
    the FITS image *must* have been previously calibrated astrometrically, so
    that the right ascensions and declinations listed in the text file are
    meaningful. Returns a QPhot object, using None as the magnitude of those
    astronomical objects that are INDEF (i.e., so faint that qphot could not
    measure anything) and positive infinity if they are saturated (i.e., if
    one or more pixels in the aperture are above the saturation level).

    Arguments:
    img - the fitsimage.FITSImage object on which to do photometry.
    coords_path - the text file containing the celestial coordinates for the
                  astronomical objects to be measured. These objects must be
                  listed one per line, in two columns: right ascension and
                  declination, in this order.
    aperture - the aperture radius, in pixels.
    annulus - the inner radius of the sky annulus, in pixels.
    dannulus - the width of the sky annulus, in pixels.
    maximum - number of ADUs at which saturation arises. If one or more pixels
              in the aperture are above this value, the magnitude of the
              astronomical object is set to positive infinity. For coadded
              observations, the effective saturation level is obtained by
              multiplying this value by the number of coadded images.
    exptimek - the image header keyword containing the exposure time. Needed
               by qphot in order to normalize the computed magnitudes to an
               exposure time of one time unit.
    uncimgk - the image header keyword containing the path to the image used to
              check for saturation. It is expected to be the original FITS file
              (that is, before any calibration step, since corrections such as
              flat-fielding may move a saturated pixel below the saturation
              level) of the very image on which photometry is done. If this
              argument is set to an empty string or None, saturation is checked
              for on the same FITS image used for photometry, 'img'.

    """

    img_qphot = QPhot(img.path, coords_path)
    img_qphot.run(annulus, dannulus, aperture, exptimek)

    # How do we know whether one or more pixels in the aperture are above a
    # saturation threshold? As suggested by Frank Valdes at the IRAF.net
    # forums, we can make a mask of the saturated values, on which we can do
    # photometry using the same aperture. If we get a non-zero flux, we know it
    # has saturation: http://iraf.net/forum/viewtopic.php?showtopic=1466068

    if not uncimgk:
        orig_img_path = img.path

    else:
        orig_img_path = img.read_keyword(uncimgk)
        if not os.path.exists(orig_img_path):
            msg = "image %s (keyword '%s' of image %s) does not exist"
            args = orig_img_path, uncimgk, img.path
            raise IOError(msg % args)

    try:
        # Temporary file to which the saturation mask is saved
        basename = os.path.basename(orig_img_path)
        mkstemp_prefix = "%s_satur_mask_%d_ADUS_" % (basename, maximum)
        kwargs = dict(prefix = mkstemp_prefix,
                      suffix = '.fits', text = True)
        mask_fd, satur_mask_path = tempfile.mkstemp(**kwargs)
        os.close(mask_fd)

        # IRAF's imexpr won't overwrite the file. Instead, it will raise an
        # IrafError exception stating that "IRAF task terminated abnormally
        # ERROR (1121, "FXF: EOF encountered while reading FITS file".
        os.unlink(satur_mask_path)

        # The expression that will be given to 'imexpr'. The space after the
        # colon is needed to avoid sexigesimal interpretation. 'a' is the first
        # and only operand, linked to our image at the invokation of the task.
        expr = "a>%d ? 1 : 0" % maximum
        logging.debug("%s: imexpr = '%s'" % (img.path, expr))
        logging.debug("%s: a = %s" % (img.path, orig_img_path))
        pyraf.iraf.images.imexpr(expr, a = orig_img_path,
                                 output = satur_mask_path, verbose = 'no')

        # Now we just do photometry again, on the same pixels, but this time on
        # the saturation mask. Those objects for which we get a non-zero flux
        # will be known to be saturated and their magnitude set to infinity.
        mask_qphot = QPhot(satur_mask_path, coords_path)
        mask_qphot.run(annulus, dannulus, aperture, exptimek)

        assert len(img_qphot) == len(mask_qphot)
        for object_phot, object_mask in itertools.izip(img_qphot, mask_qphot):
            assert object_phot.x == object_mask.x
            assert object_phot.y == object_mask.y

            if object_mask.flux > 0:
                object_phot = object_phot._replace(mag = float('infinity'))
    finally:

        # Remove saturation mask. The try-except is necessary because the
        # deletion may fail (OSError), or something could go wrong and an
        # exception be raised before the variable is defined (NameError).

        try:
            os.unlink(satur_mask_path)
        except (NameError, OSError):
            pass

    return img_qphot

