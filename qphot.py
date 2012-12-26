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

# Tell PyRAF to skip all graphics initialization and run in terminal-only mode.
# Otherwise we will get annoying warning messages (such as "could not open
# XWindow display" or "No graphics display available for this session") when
# working at a remote terminal or at a terminal without any X Windows support.
# Any tasks which attempt to display graphics will fail, of course, but we are
# not going to make use of any of them, anyway.

import os
os.environ['PYRAF_NO_DISPLAY'] = '1'
import pyraf.iraf
from pyraf.iraf import digiphot, apphot  # 'digiphot.apphot' package

# PyRAF versions older than 2.0 caused IRAF's imexpr to raise an obscure
# exception ("Error in parameter 'a' for task imexpr\n'Key a not found") from
# time to time. This only happened, at least on our machines, when the work was
# divided up between multiple processes using the multiprocessing module, and
# it did not seem to occur with pools of only one worker.

# This line converts a string like '2.0-r1785' to the tuple (2, 0)
pyraf_version = tuple(int(x) for x in pyraf.__version__[:3].split('.'))
if pyraf_version < (2, 0):
    raise ImportError("PyRAF version 2.0 or newer is required")

# Turn PyRAF process caching off; otherwise, if we spawn multiple processes
# and run them in parallel, each one of them would use the same IRAF running
# executable, which could sometimes result in the most arcane of errors.
pyraf.iraf.prcacheOff()

import copy
import logging
import math
import os.path
import tempfile

# LEMON modules
import astromatic
import fitsimage
import methods
import seeing

class QPhotResult(object):
    """ This class simply encapsulates the photometry of a star. In other
        words, each one of the lines that for each star will be outputted by
        IRAF's qphot will be parsed and saved as an instance of this class."""

    def __init__(self, xcenter, ycenter, mag, sum_, flux, stdev):
        """ Instantiation method for the QPhotResult class.

        xcenter - x coordinate of the center of the star.
        ycenter - y coordinate of the center of the star.
        mag - the instrumental magnitude of the star in the aperture.
        sum_ - the total number of counts in the aperture _including_ the sky.
        flux - the total number of counts in the aperture _excluding_ the sky.
        stdev - the standard deviation of the best estimate of the sky value,
                per pixel.

        """

        # Historical note: here there used to be an assertion that made sure
        # that sum >= flux, as the total number of counts in the aperture
        # including the sky should never, ever be smaller than when the sky is
        # excluded. However, this may not be true if we attempt to do
        # photometry on something that is not an astronomical object.
        #
        # How could this possibly happen, you may ask. Well, in our case, at
        # least, this was due to SExtractor, when sources were detected in the
        # reference image, which identified some points in the extremely noisy
        # (because of the flat-field reduction step) margins as stars. Anyway,
        # I suspect that the assertion may also fail for large (bigger than the
        # aperture) objects.
        #
        # The moral of the story, my developer fellow, is that the S/N of the
        # star also provides information about, well, whether it is _actually_
        # a star: negative signal-to-noise ratios clearly indicate that what
        # had photometry done on was definitely not a star.

        self.x         = xcenter
        self.y         = ycenter
        self.magnitude = mag
        self.sum       = sum_
        self.flux      = flux
        self.stdev     = stdev

    def snr(self, gain):
        """ Return the signal-to-noise ratio of the star.

        The method returns the S/N, a quantitative measure of the accuracy with
        which the star was observed. The signal-to-noise ratio tells us the
        relative size of the desired signal to the underlying noise or
        background light. The noise is defined as the standard deviation of a
        single measurement from the mean of the measurements made on a star.

        For photon arrivals, the statistical noise fluctuation is represented
        by the Poisson distribution, and for bright sources where the sky
        background is negligible, S/N = total counts / sqrt(total counts).
        Note that when the sky is not insignificant, the formula becomes
        S/N = (total counts - sky counts) / sqrt(total counts).

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

        """

        if gain <= 0:
            raise ValueError("CCD gain must be a positive value")

        # Division by zero must be avoided, as sum and flux can both be zero if
        # photometry is done on a pair of coordinates that do not correspond to
        # any star (or, in other words, if the star on which we are doing
        # photometry is so faint that it is not visible in the image).

        if not self.sum:
            return 0.0

        # Also, sum and flux could be negative if photometry, as already
        # explained, is done on the margins (overscan area and such) of the
        # image, where the calibration of the images may have resulted in
        # pixels with negative values. In that case we do not return zero, but
        # instead a _negative_ SNR, so that it is clear that almost surely what
        # had photometry done on was not a star.

        elif self.sum < 0.0:
            return -(abs(self.flux * gain) / math.sqrt(abs(self.sum * gain)))

        else:
            return (self.flux * gain) / math.sqrt(self.sum * gain)


class QPhot(list):
    """ The photometry of an image, as returned by IRAF's qphot.

    This class encapsulates the photometry done by the qphot (quick aperture
    photometer) task implemented in IRAF. An instance of QPhotResult will be
    created for each star detected in the image, except for those stars that
    were INDEF or saturated -- for these, None is used.

    """

    def __init__(self, path):
        """ Instantiation method for the QPhot class.

        path - path to the image whose photometry is to be done.

        """

        super(list, self).__init__()
        self.image = fitsimage.FITSImage(path)

    @property
    def path(self):
        """ Return the path to the FITS image. """
        return self.image.path

    def clear(self):
        """ Remove from the instance the photometry of all the stars. """
        del self[:]

    def run(self, annulus, dannulus, aperture, exptimek, pixels):
        """ Run IRAF's qphot on the image.

        This method is a wrapper, equivalent to (1) running 'qphot' on an image
        and (2) using 'txdump' in order to extract some fields from the text
        database that the former outputs. The goal of this routine is to make
        it possible to easily do photometry on an image, saving in the current
        instance the information for all the objects (that is, when running
        txdump the parameter 'expr' is set to 'yes').

        This method may be seen as a low-level routine: INDEF stars will have
        *a magnitude of None*, but it does not pay attention to the saturation
        levels. For that, you will have to use the photometry function, defined
        in this very module. That, and not this one, is the routine you are
        expected to use to do photometry on your images.

        In the first step, photometry will be done for the objects whose
        coordinates are listed in 'pixels', a sequence of two-element tuples
        (x- and y- values). These coordinates are saved to a temporary file,
        which is given as input to 'qphot' to do photometry on these pixels.

        The output of 'qphot' is saved to a second temporary file (b), on which
        'txdump' is run to extract the fields from the APPHOT text database
        that was produced by 'qphot' to a third temporary file (c). This last
        file is then parsed, the information of each of its lines used in order
        to create an instance of QPhotResult, one per star, which is added to
        the current instance. Note that previous results will be lost if this
        method is run twice on the same instance.

        All the temporary files (a, b, c) are guaranteed to be deleted before
        the method exits, even if an error is encountered. The method returns
        the number of stars whose photometry has been done.

        Arguments:
        annulus - the inner radius of the sky annulus, in pixels.
        dannulus - the width of the sky annulus, in pixels.
        aperture - the single aperture radius, in pixels.
        exptimek - the image header keyword containing the exposure time. Needed
                   by qphot in order to normalize the computed magnitudes to an
                   exposure time of one time unit.
        pixels - list of instances of Pixel which encapsulate the image
                 coordinates that will be passed to IRAF's qphot. It set to None
                 (the default value) or left empty, these coordinates will be
                 obtained by running SExtractor on the image.

        """

        self.clear() # empty the current instance

        try:
            # The task qphot will receive a text file with the coordinates for
            # the objects to be measured. We need, then, to create a temporary
            # file, with the coordinates of the stars, one per line.
            coords_fd, stars_coords = \
                tempfile.mkstemp(prefix = os.path.basename(self.path),
                                 suffix = '.qphot_coords', text = True)
            for pixel in pixels:
                os.write(coords_fd, "%f\t%f\n" % (pixel.x, pixel.y))
            os.close(coords_fd)

            # Temporary file to which the text database produced by qphot will
            # be saved. Even if empty, it must be deleted before calling
            # qphot. Otherwise, an error message, stating that the operation
            # "would overwrite existing file", will be thrown.
            output_fd, qphot_output = \
                tempfile.mkstemp(prefix = os.path.basename(self.path),
                                 suffix = '.qphot_output', text = True)
            os.close(output_fd)
            os.unlink(qphot_output)

            # The list of the fields passed to qphot
            qphot_fields = ['xcenter', 'ycenter', 'mag', 'sum', 'flux', 'stdev']

            # Run IRAF's qphot on the image and save the output to our
            # temporary file. Note that cbox *must* be set to zero, as
            # otherwise IRAF's qphot will compute accurate centers for each
            # object using the centroid centering algorithm. This is generally
            # a good thing, but in this case we need the photometry to be done
            # exactly on the specified coordinates.

            kwargs = dict(cbox = 0, annulus = annulus, dannulus = dannulus,
                          aperture = aperture, coords = stars_coords,
                          output = qphot_output, exposure = exptimek,
                          interactive = 'no')

            # Ignore the annoying warning message that is, from time to time,
            # printed to the standard error, and that look like what follows:
            #
            # Exception pyraf.subproc.SubprocessError: SubprocessError("Failed
            # kill of subproc 24600, '/iraf/iraf/bin.linux/x_images.e -c', with
            # signals ['TERM', 'KILL']",) in <bound method Subprocess.__del__
            # of <Subprocess '/iraf/iraf/bin.linux/x_images.e -c', at
            # 7f9f3f408710>> ignored
            #
            # This message is printed because of an uncaught exception in the
            # __del__() method of the PyRAF's Subprocess class. As explained in
            # the Python data model, those exceptions are ignored and a warning
            # message printed to stderr instead. Using this context manager we
            # can get rid of these messages.

            with methods.ignore_del_exceptions():
                apphot.qphot(self.path, **kwargs)

            # Make sure the outpout was written to where we said
            assert os.path.exists(qphot_output)

            # Now extract the specified records from the qphot output. We need
            # the path (the file, even empty, must be deleted, as IRAF won't
            # overwrite it) another temporary file to which save the output;
            # IRAF's txdump 'Stdout' (why the upper case?) redirection must
            # be to a file handle or string.

            txdump_fd, txdump_output = \
                tempfile.mkstemp(prefix = os.path.basename(self.path),
                                 suffix ='.qphot_txdump', text = True)
            os.close(txdump_fd)
            os.unlink(txdump_output)

            # The cast of Stdout to string is needed as txdump won't work with
            # unicode, if we happen to come across it -- it would insist on
            # that redirection must be to a file handle or string.
            pyraf.iraf.txdump(qphot_output, fields = ','.join(qphot_fields),
                              Stdout = str(txdump_output), expr = 'yes')

            # Now open the outut file again and parse the output of 'txdump',
            # creating an instance of QPhotResult and saving it in the current
            # instance for each line of the file.
            txdump_fd = open(txdump_output, 'rt')
            txdump_lines = txdump_fd.readlines()

            assert len(txdump_lines) == len(pixels)
            for line, pixel in zip(txdump_lines, pixels):
                try:
                    # The i-th line in the file corresponds to the i-th pixel.
                    splitted_line = line.split()
                    xcenter = float(splitted_line[0]) # same value as pixel.x
                    ycenter = float(splitted_line[1]) # same value as pixel.y

                    if __debug__:
                        epsilon = 0.001  # a thousandth of a pixel!
                        assert abs(xcenter - pixel.x) < epsilon
                        assert abs(ycenter - pixel.y) < epsilon

                    try:
                        mag_str = splitted_line[2]
                        mag     = float(mag_str)
                    except ValueError:  # raised by float("INDEF")
                        assert mag_str == 'INDEF'
                        mag = None

                    sum_ = float(splitted_line[3])
                    flux = float(splitted_line[4])

                    try:
                        stdev_str = splitted_line[5]
                        stdev = float(stdev_str)
                    except ValueError:  # raised by float("INDEF")
                        assert stdev_str == 'INDEF'
                        stdev = None

                    star_phot = QPhotResult(xcenter, ycenter, mag, sum_, flux, stdev)
                    self.append(star_phot)

                except IndexError:
                    raise  # we aren't expecting any improperly formatted line

        finally:

            # Remove the temporary files, as we do not longer need them. The
            # try-except is needed to prevents raised exceptions in case a file
            # does not get to be created -- that is, if an IRAF task fails or
            # the execution of the module is aborted by the user. NameError is
            # needed because an exception may be raised before the variables
            # are declared.

            try:
                os.unlink(stars_coords)
            except (NameError, OSError):
                pass

            try:
                os.unlink(qphot_output)
            except (NameError, OSError):
                pass

            try:
                os.unlink(txdump_output)
            except (NameError, OSError):
                pass

        return len(self)

def run(xml_offset, aperture, annulus, dannulus,
        maximum, exptimek, uncimgk, pixels):
    """ Do photometry on an image.

    The method receives a xmlparse.XMLOffset instance and runs IRAF's qphot in
    the shifted image on those those stars whose coordinates are listed in the
    'pixels' keyword argument. Where the stars are in the shifted image is
    determined by applying to their cooredinates the translation offset given
    by the XMLOffset instance, thus aligning both images and 'correcting' the
    star coordinates.

    Returns a qphot.QPhot instance, which contains a magnitude of None for
    INDEF stars and positive infinity for saturated stars (those with one or
    more pixels in the aperture above the saturation level).

    Arguments:
    xml_offset - a XMLOffset instance, encapsulating the translation offset
                 between the reference and the shifted FITS images.
    aperture - the single aperture radius, in pixels.
    annulus - the inner radius of the sky annulus, in pixels.
    dannulus - the width of the sky annulus, in pixels.
    maximum - level at which arises saturation in the shifted image, in
              ADUs. If one or more pixels in the aperture of the star are above
              this value, it will be considered to be saturated. For coadded
              observations, the effective saturation level is obtained by
              multiplying this value by the number of coadded images.
    exptimek - the image header keyword containing the exposure time. Needed
               by qphot in order to normalize the computed magnitudes to an
               exposure time of one time unit.
    uncimgk - keyword for the relative path to the uncalibrated shifted image;
              which will be one used to check whether pixels are saturated.
              This value may be set to an empty string or None if saturation is
              to be checked for on the same image in which photometry is done.
    pixels - sequence of instances of Pixel which encapsulate the image
             coordinates that will be passed to IRAF's qphot.

    """

    # The instantiation method of FITSImage makes sure that the image on which
    # photometry will be done does exist and is standard-conforming.
    shifted_image = fitsimage.FITSImage(xml_offset.shifted)

    # Apply the offset to the coordinates of each star.
    shifted_pixels = []
    for pixel in pixels:
        x = pixel.x + xml_offset.x
        y = pixel.y + xml_offset.y
        shifted = astromatic.Pixel(x, y)
        shifted_pixels.append(shifted)

    if __debug__:
        assert len(pixels) == len(shifted_pixels)
        epsilon = 0.001  # a thousandth of a pixel!
        for pixel, shifted in zip(pixels, shifted_pixels):
            assert abs(shifted.x - xml_offset.x - pixel.x) < epsilon
            assert abs(shifted.y - xml_offset.y - pixel.y) < epsilon

    # And, finally, do photometry on the shifted image
    img_qphot = QPhot(shifted_image.path)
    img_qphot.run(annulus, dannulus, aperture, exptimek, shifted_pixels)

    # How do we know whether one or more pixels in the aperture are above a
    # saturation threshold? As suggested by Frank Valdes at the IRAF.net
    # forums, we can make a mask of the saturated values, on which we can do
    # photometry using the same aperture. If we get a non-zero flux, we know
    # it has saturation. [URL http://iraf.net/phpBB2/viewtopic.php?t=90621]
    #
    # In order to check for the saturation of a star, we can use the same image
    # on which we are doing photometry or, as Matilde Fernandez suggested, the
    # _original_ image -- the raw one, before any calibration step, as
    # processes such as flat-fielding may move a saturated pixel below the
    # saturation level. Far-sighted as we are, the path to this primeval image
    # was stored in the FITS header by the import.py module, the first of the
    # pipeline. However, in case we want to use the same image, or if the
    # keyword is not available for these images, 'uncimgk' may be set to None.

    # Note that, if 'uncimgk' is given, the path to the original image will be
    # loaded and the execution of the script aborted. The path stored in the
    # keyword is trusted blindly, so you better do not modify it to a different
    # value.

    if not uncimgk:
        orig_img_path = shifted_image.path

    else:
        orig_img_path = shifted_image.read_keyword(uncimgk)
        if not os.path.exists(orig_img_path):
            msg = "image %s (keyword '%s' of image %s) does not exist" % \
                  (orig_img_path, uncimgk, shifted_image.path)
            raise IOError(msg)

    try:
        # The temporary file to which the saturation mask is saved
        basename = os.path.basename(orig_img_path)
        mkstemp_prefix = "%s_satur_mask_%d_ADUS_" % (basename, maximum)
        mask_fd, satur_mask_path = \
            tempfile.mkstemp(prefix = mkstemp_prefix,
                             suffix = '.fits', text = True)

        # IRAF's imexpr won't overwrite the file. Instead, it will raise an
        # IrafError exception stating that "IRAF task terminated abnormally
        # ERROR (1121, "FXF: EOF encountered while reading FITS file", with the
        # path to the output file. It seems it is attempting to read it (and
        # failing, since it's empty) even although it is the output path.
        # Anyway, we can fix this by simply deleting it.
        os.close(mask_fd)
        os.unlink(satur_mask_path)

        # The expression that will be given to 'imexpr'. The space after the
        # colon is needed to avoid sexigesimal intepretation. 'a' is the first
        # and only operand, linked to our image at the invokation of the task.
        expr = "a>%d ? 1 : 0" % maximum
        logging.debug("%s: imexpr = '%s'" % (shifted_image.path, expr))
        logging.debug("%s: a = %s" % (shifted_image.path, orig_img_path))
        pyraf.iraf.images.imexpr(expr, a = orig_img_path,
                                 output = satur_mask_path, verbose = 'no')

        # Now we just do photometry again, on the same pixels, but this time on
        # the saturation mask. Those stars for which we get a non-zero flux
        # will be known to be saturated and their magnitude set to infinity.
        mask_qphot = QPhot(satur_mask_path)
        mask_qphot.run(annulus, dannulus, aperture, exptimek, shifted_pixels)

        assert len(img_qphot) == len(mask_qphot)
        for star_phot, star_mask in zip(img_qphot, mask_qphot):
            assert star_phot.x == star_mask.x
            assert star_phot.y == star_mask.y

            if star_mask.flux > 0:
                star_phot.mag = float('infinity')
    finally:
        try:
            os.unlink(satur_mask_path)
        # NameError needed as something may go wrong before it's declared
        except (NameError, OSError):
            pass

    return img_qphot

