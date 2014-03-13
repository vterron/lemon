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

description = \
""" This module receives a series of offsets, aligns their corresponding FITS
    images and finally combines them all into a single image by averaging them
    after rejecting a certain fraction of the lowest and highest pixels. The
    offsets are read from the XML file received as input, which is expected to
    have been outputted by the 'offset.py' module. Note that, therefore, all
    the offsets must refer to the same reference image, as only one reference
    FITS image is created.

    A very important note: as you surely already know, shifting images, as is
    needed in order to have them aligned, distorts slightly both signal and
    noise. This means that you should not (even worse, cannot!) do photometry
    on the resulting image. It should be only used in order to detect sources,
    as it maximizes the signal-to-noise ratio. Please stand on the shoulders of
    the photometrist giants and adhere to the standards in our guild.

    Note that it is assumed that, for offsets, North is up and East is left.
    Images for which North is down or East is right are not yet supported, and
    using them shall have catastrophic results. Yes, this is the unexpected
    consequence of having most of your astronomers and astronomical software
    developers based in countries in the Northern hemisphere.

"""

import atexit
import logging
import math
import montage_wrapper as montage
import numpy
import os
import os.path
import optparse
import pyfits
import shutil
import sys
import tempfile
import traceback

# LEMON modules
import customparser
import keywords
import fitsimage
import methods
import style
import xmlparse

def clean_tmp_dir(dir_path):
    """ Try to delete an entire directory tree. """

    msg = "Cleaning up temporary directory '%s'"
    logging.debug(msg % dir_path)

    error_count = []

    def log_error(function, path, excinfo):
        """ Error handler for shutil.tree() """

        # nonlocal is not available in Python 2.x so, being it outside of the
        # local scope, we cannot rebind 'error_count' each time we come across
        # an error. Instead of initializing it to zero and incrementing it by
        # one every time this function is called, use it as a list, appending
        # an element for each error.
        error_count.append(1)
        msg = "%s: error deleting '%s' (%s)"
        args = function, path, excinfo[1]
        logging.debug(msg % args)

    try:
        kwargs = dict(ignore_errors = False, onerror = log_error)
        shutil.rmtree(dir_path,**kwargs)

    finally:
        msg = "Temporary directory '%s' deleted"
        if error_count:
            msg += " (but there were failed removals)"
        logging.debug(msg % dir_path)

class Mosaic(object):
    """ Encapsulates the canvas needed to mosaic a set of images.

    The purpose of this class is to be used in order to determine the size of
    the canvas that is needed so that all the images of the same field can be
    combined without any information being lost. In simple terms, what we are
    trying to do is to shift all the images, so that they become perfectly
    aligned, and combine them, as if we were composing a mosaic.

    Ideally, we would simply do what is described above: align and combine.
    However, the underlying software that does the shifting, IRAF's imshift,
    does not alter the size of the FITS image: it simply moves the pixels of
    data. This means that some information will be unavoidably lost, as the
    shift causes a region of the image to fall outside of it. This, clearly,
    is unacceptable for out purposes here.

    The solution, thus, is to resize each image by adding to it a blank margin
    big enough so that no information is lost. That's easy, you may argue: if
    the image is going to be shifted, say, fifty pixels to the left, we just
    have to add a margin of that same size, fifty pixels, to the left margin of
    the image, don't we? Well, no.

    Unfortunately, it is not that easy, because all the images that are to be
    combined must be of the _same_size_. The images on which the pipeline works
    must meet this requirement, and so must do the images after being aligned.
    Therefore, we need to resize all the images to the same dimensions. We do
    this by means of calculating what the size of the final mosaic will be, and
    then resizing all the images to that size before the shifting is done.

    The size of the mosaic is expressed in terms of the margin that has to be
    added to the reference (the non-shifted, to which all the offsets refer)
    image. This may see counterintuitive, but actually makes perfect sense, as
    the exact margin to be added to each image can be derived from each offset.

    As an example, assume that the margin that has to be added to the reference
    image is 50 pixels on the left and 100 on the right (meaning that the size
    of the mosaic will be 50 + 100 = 150 pixels larger on the x-axis than the
    original images). Say that we want to align one of the images, whose offset
    in the x-axis (the one that concerns this example) is 15 pixels (that is,
    it is shifted to the right). Which margin must be added before shiting it?
    It's simple: 50 + 15 = 65 on the left and 100 - 15 = 85 on the right.

    In summary, this class is used in order to determine the margins that will
    have to be added to each image being combined. Note that the size of the
    canvas is determined by the maximum offset (on each axis) among the images
    that are to be combined. Too shifted images, thus, may result in rather
    larger than expected mosaics.

    """

    @staticmethod
    def _maximum_offsets(offsets):
        """ Determine what are the maximum offsets among the shifted images.

        This method receives a sequence of offsets (XMLOffset instances) and
        calculates the maximum offset that will have to be applied on both
        directions, for both axes when the images are aligned. In other words,
        the method determines the maximum offsets that, among the shifted
        images, there are to the left (x-axis, negative), right (x-axis,
        positive), downwards (y-axis, negative) and upwards (y-axis, negative).
        These values are returned, in that order, as a four-element tuple, and
        establish the size of the margins to be added to each image, and by
        extension also the dimensions of the mosaic.

        Implementation detail: the XMLOffset class used by the pipeline stores
        'how much' each shifted image is, well, shifted from the reference
        image, but in order to align the images we need to negate these values.
        In other words: the offset (-15.6, 56.9), which means that the shifted
        image is moved 15.6 pixels to the left and 56.9 upwwards, becomes for
        our purposes (15.6, -56.9), as in order to align the image with the
        reference one we need to shift it 15.6 pixels to the right and 56.9
        downwards. This is absolutely transparent to the user, but it is worth
        explaining here in case you wonder why these values are negated below.

        """

        # The zeros appended to both lists are the 'default' values. If we did
        # not use them, the 'max' or 'min' built-in functions would fail below
        # if all the offsets in one of the axes happened to be positive or
        # negative, as they would receive an empty sequence.

        x_offsets = [-offset.x for offset in offsets] + [0]
        y_offsets = [-offset.y for offset in offsets] + [0]

        right  = max(x for x in x_offsets if x >= 0)
        left   = min(x for x in x_offsets if x <= 0)
        top    = max(y for y in y_offsets if y >= 0)
        bottom = min(y for y in y_offsets if y <= 0)

        return left, right, bottom, top

    def __init__(self, offsets):
        """ Instantiation method for the Mosaic class.

        The Mosaic constructor receives a sequence of offsets (encapsulated as
        XMLOffset instances) and determines what the final dimensions will be.
        Due to the computation cost, the shifted images are neither aligned not
        combined until the methods Mosaic.align(index) and Mosaic.combine() are
        invoked. Right after instantiation, however, it can be known what the
        size of the mosaic, when the images are aligned and combined, will be.

        Although the previous methods of the pipeline should have taken care of
        this, the method checks that all the shifted images are of the same
        size and that their offsets are given with respect to the same image.

        """

        if not offsets:
            raise TypeError("__init__() arg is an empty sequence")

        self.reference = fitsimage.FITSImage(offsets[0].reference)

        for offset in offsets:
            if offset.reference != self.reference.path:
                raise ValueError("not all the offsets refer to the same image")
            if fitsimage.FITSImage(offset.shifted).size != self.reference.size:
                raise ValueError("not all the images have the same size")

        # The offsets are stored in a dictionary, which maps each offset to the
        # path to the aligned image. The initial value, None, indicates that
        # the shifted image for that offset has not yet been aligned. The class
        # shall make sure that all the images have been aligned (i.e., that
        # there is no key in the dictionary with a value of None) before
        # combining all the images and generating the resulting mosaic.

        self.offsets = dict((x, None) for x in offsets)

        # Compute and store the maximum offset that will have to be applied on
        # both directions, for both axes when the images are aligned. These
        # determine the width of the margins that will have to be added to the
        # images, and therefore the size of the mosaic.

        self.min_x_offset, self.max_x_offset, \
        self.min_y_offset, self.max_y_offset = Mosaic._maximum_offsets(self)

        # This private attribute is used in order to keep the track of how many
        # images have been overlapped for each pixel of the resulting image.
        # Here, at instantiation time, a two-dimensional matrix with the same
        # size (that is, with an equal number of rows and columns) and filled
        # with zeros is created. Each position of the matrix, thus, corresponds
        # exactly to a pixel of the mosaic.
        #
        # Later on, every time an image is aligned, this attribute is updated
        # and the positions corresponding to the non-blank pixels on the
        # aligned image will be incremented by one. This will allow us, when
        # the images are combined, to remove the 'irrelevant' parts, namely
        # those pixels which result from having combined margins for the most
        # part. In this manner, it becomes possible to remove from the output
        # mosaic those parts that are the product of combining less than ten
        # images, for example.

        self._pixels_counter = numpy.zeros(self.size[::-1], dtype = float)

    def __enter__(self):
        """ Enter the runtime context related to this Mosaic object.

        The context manager (that is, the with statement) simplifies the task
        of managing resources that need to be cleaned up. In our case, the
        images that may have been aligned are no longer needed once the mosaic
        is generated. In case we want to have them automatically removed once
        they have been combined, regardless of any error that we may run into,
        the context manager can be used.

        """
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """ Exit the runtime context related to this Mosaic object.

        The method automatically removes the temporary images that have been
        aligned during the lifetime of the object, and that are no longer
        needed once the mosaic has been created.

        """

        for aligned_img in self.offsets.values():
            try:
                os.unlink(aligned_img.path)
            # Nones cannot be deleted; other errors may occur
            except (TypeError, IOError):
                pass

    @property
    def right(self):
        """ Return the right margin needed to mosaic the images. """
        return int(math.ceil(abs(self.max_x_offset)))

    @property
    def left(self):
        """ Return the left margin needed to mosaic the images. """
        return int(math.ceil(abs(self.min_x_offset)))

    @property
    def top(self):
        """ Return the top margin needed to mosaic the images. """
        return int(math.ceil(abs(self.max_y_offset)))

    @property
    def bottom(self):
        """ Return the bottom margin needed to mosaic the images. """
        return int(math.ceil(abs(self.min_y_offset)))

    @property
    def size(self):
        """ Return the size of the mosaic.

        The method returns the dimensions of the resulting mosaic, once all the
        shifted images are aligned and combined. The size is returned in a
        two-element tuple, whose first element is the number of columns (that
        is, width of the image) and the second the number of rows (height).

        """

        x_rsize, y_rsize = self.reference.size
        mosaic_x_size = x_rsize + self.left + self.right
        mosaic_y_size = y_rsize + self.top  + self.bottom
        return mosaic_x_size, mosaic_y_size

    def __str__(self):
        """ The 'informal' string representation of a Mosaic. """
        return '%s(%s)' % (self.__class__.__name__, self.reference.path)

    def __len__(self):
        """ Return the number of shifted images that are to be mosaicked.

        The method returns the total number of images that will be aligned with
        the reference image a combined into what we call the 'mosaic'. What the
        method is internally doing is returning the number of offsets that were
        passed to the instance at construction time, and from which the shifted
        images are aligned, one by one.

        """

        return len(self.offsets)

    def __getitem__(self, index):
        """ Return the offset of the index-th shifted image.

        The method returns the index-th offset that was passed to the instance
        at construction time. In other words, what is returned is the object
        that encapsulates the translation offset of the index-th shifted image
        that is to be aligned and combined. For the aligned image, see the
        Mosaic.aligned(index) image.

        """

        return self.offsets.keys()[index]

    def aligned(self, index):
        """ Return the index-th shifted image, already aligned.

        The method returns the FITSImage instance of the already-aligned
        version of the index-th shifted image. That is, for the index-th
        shifted image that is to be combined, the instance of the FITSImage
        class which encapsulates it after its alignment with the reference
        image is returned. In case it has not yet been aligned -- see the
        Mosaic.align(index) method --, None is returned.

        """

        return self.offsets[self.offsets.keys()[index]]

    def _set_aligned(self, index, aligned):
        """ Set the aligned-version of the index-th shifted image. """
        self.offsets[self.offsets.keys()[index]] = aligned

    def align(self, index, checkimage_type = None):
        """ Align the index-ith shifted image.

        The method aligns the index-th shifted image with the reference image,
        saving the resulting image to a temporary file whose path is internally
        stored. In order to retrieve it, use the Mosaic.aligned(index) method.
        Note that in order to have these aligned files automatically deleted,
        as they are no longer need once they have been combined and the mosaic
        obtained, the context manager of the class must be used.

        The 'checkimage_type' determines the SExtractor check-image that will
        be used in the mosaic for each image. In other words: the images are not
        directly combined, but their check-image computed and taken instead. The
        default value, '-BACKGROUND', means that the background-subtracted
        version of each image will be used. Please refer to the SExtractor
        documentation or the fitsimage.check_image() method for further
        information on the different types of check-image available. In case no
        check-image is needed, the parameter should be set to None.

        Note for developers: the check-image generation was originally done at
        Mosaic.combine, but it was moved here as the module shows a progress
        bar as the images are aligned. This is better than the initial approach,
        in which the user had no way of knowing how much progress the task of
        combining the images had made -- as they were not only aligned, but also
        their check-images generated, a very CPU-intensive task. Although now a
        parameter in this method, you should, of course, use the same kind of
        check-image across all the images in the mosaic.

        """

        # The offsets in both axes have to be negated as here we are not
        # interested in how much the shifted image is moved with respect to the
        # reference once: what we need to know is how much we have to shift the
        # shifted image so that it becomes aligned with the reference. In other
        # words, we are correcting the misalignment.

        shifted_path = self[index].shifted
        x_offset     = -self[index].x
        y_offset     = -self[index].y

        # Before the image can be shifted, and as it has already been explained
        # several times, a blank margin has to be added to the image, so that
        # its new size matches that of the mosaic. This file with the margins
        # is only useful as the input to FITSImage.imshift(xshift, yshift);
        # hence the try-finally, to make sure it is always deleted, under all
        # circumstances, no matter which errors we may encounter.

        shifted_img = fitsimage.FITSImage(shifted_path)

        # If needed, compute the SExtractor check-image before aligning
        if checkimage_type:
            shifted_img = shifted_img.check_image(check_type = checkimage_type)

        try:
            margin_img = shifted_img.add_margin(self.left, self.right,
                                                self.bottom, self.top)

            prefix = '%s_aligned_' % shifted_img.basename_woe
            kwargs = dict(prefix = prefix)
            aligned_img = margin_img.imshift(x_offset, y_offset, **kwargs)
            self._set_aligned(index, aligned_img)

            comment_msgs = (
                "Image shifted by LEMON on %s" % methods.utctime(),
                "[imshift] xout = xin + xshift",
                "[imshift] yout = yin + yshift",
                "[imshift] Blank frames were added before the shift")

            history_msgs = (
                "[imshift] Original image: %s" % shifted_img.path,
                "[imshift] Pixel shift in x-axis: %.3f" % x_offset,
                "[imshift] Pixel shift in y-axis: %.3f" % y_offset,
                "[imshift] Left frame, width (pixels): %d" % self.left,
                "[imshift] Right frame, width (pixels): %d" % self.right,
                "[imshift] Bottom frame, width (pixels): %d" % self.bottom,
                "[imshift] Top frame, width (pixels): %d" % self.top)

            for msg in comment_msgs:
                aligned_img.add_comment(msg)
            for msg in history_msgs:
                aligned_img.add_history(msg)

        finally:
            try:
                margin_img.unlink()
            except:  # such as NameError if not yet instantiated
                pass

        # The shifted image has been aligned and its path internally stored.
        # Now we need to update _pixels_counter, the attribute used to keep the
        # track of how many images have been overlapped for each pixel of the
        # resulting image. In order to do that, we create another matrix of the
        # same dimensions filled with zeros, set to one those positions that do
        # not correspond to a blank margin (or, put another way, we set to one
        # those pixels of the now aligned image for which there is actual data)
        # and add it to _pixels_counter, updating its value.
        #
        # Note that, although pixels are the smallest unit of a picture, images
        # may have to be shifted by a fractional number of pixels in order to
        # become aligned. This poses no problem for us, since, according to its
        # help page, IRAF's imshift determines the value of the output pixels
        # "by interpolating in the input image at the positions of the shifted
        # output pixels". However, this means that the pixel counter cannot be
        # exactly precise, but just an approximation: in order to establish
        # which parts of the aligned image are just blank margins, we round
        # each offset (x- and y-axis) to the nearest integer.

        mosaic_x_size, mosaic_y_size = self.size

        # The 'r' stands for 'rounded'
        x_roffset = int(round(x_offset))
        y_roffset = int(round(y_offset))

        # Here we determine, for each axis, where the beginning and end of the
        # data (that is, non-blank, non-margin) portions of the image are. The
        # initial starting point of each margin is at the n-th pixel, where n
        # means the width of the margin added to the left (x-axis) or bottom
        # (y-axis), and its end is located at the n-th to last pixel, where n
        # is now the width of the margin added to the right (x-axis) or top
        # (y-axis) of the image. These are the 'initial' values, which have to
        # be corrected (using rounded-to-integer values) for the shifting that
        # has been applied to the image.
        #
        # For example, assume that the original image is 500 pixels width, and
        # that the size of the margins that are applied to it, on the x-axis,
        # are 20 pixels to the left and 35 to the right. Now the image is
        # shifted 30 pixels to the right. What is the width of the margins on
        # the aligned image on the x-axis? In other words, for which pixels are
        # there no shifted pixels on the x-axis? The margin on the left was 20
        # pixels width, but it had to be incremented by the 30 pixels by which
        # the image was shifted in the other direction. Analogously, the right
        # margin, initially 35 pixels width, has decreesed to 35-30 = 5 pixels.

        x_init = self.left   + x_roffset
        y_init = self.bottom + y_roffset
        x_end  = mosaic_x_size - self.right + x_roffset
        y_end  = mosaic_y_size - self.top   + y_roffset

        assert x_init >= 0
        assert y_init >= 0
        assert x_end <= mosaic_x_size
        assert y_end <= mosaic_y_size

        image_pixels = numpy.zeros([mosaic_y_size, mosaic_x_size])
        assert aligned_img.size == image_pixels.shape[::-1]
        image_pixels[x_init:x_end, y_init:y_end] = 1
        self._pixels_counter += image_pixels

    def combine(self, output_path, minimum = 0, flow = 0, fhigh = 0):
        """ Create the mosaic by combining the aligned images.

        This method combines the already-aligned images and returns the
        resulting mosaic as a FITSImage instance. Shifted images which have not
        yet been aligned (see the Mosaic.align(index) method) when this method
        is invoked are simply ignored. Furthermore, if no image at all has been
        aligned before the method is run, None is returned.

        The type of combining operation performed on the aligned pixels is
        IRAF's imcombine 'average' ('combine' parameter), which means that the
        arithmetic mean of the pixels is calculated. Before the combination,
        however, a rejection operation is performed on the pixels: IRAF's
        'mixmax' ('reject' parameter), meaning that the specified fraction of
        the highest and lowest pixels are rejected and therefore excluded from
        the combination. Note that IRAF's imcombine receives the number of
        pixels to be rejected as input, while this method uses the fraction of
        the total number of pixels (equal to the number of images) being
        combined ('flow' and 'fhigh' parameters).

        The 'minimum' argument may be used to remove (by means of 'blanking'
        them, in other words, setting them to zero) the pixels that result from
        having combined more blank margins than actual data, so to speak. More
        strictly speaking, this parameter determines the minimum fraction of
        data pixels (not corresponding to the margins that were added before
        the alignment) that must have been combined for each pixel in the
        output image if it is not to be zeroed. For example, a value of 0.75
        means that those pixels of the mosaic for which more than 0.25 of the
        input pixels corresponded to a margin will be set to zero after their
        combination.

        """

        align_imgs = fitsimage.FITSet()
        for index in range(len(self)):
            aligned_img = self.aligned(index)
            if aligned_img:  # ignore Nones
                align_imgs.append(aligned_img)

        combined_img = align_imgs.imcombine(output_path, flow, fhigh)

        try:

            # The _pixels_counter attribute, updated each time an image was
            # aligned, keeps the track of how many data pixels (versus what we
            # call here 'margin' pixels) have been overlapped for each pixel of
            # the mosaic. We need to convert these absolute numbers to a
            # fraction of the total of images being combined, and then set to
            # zero those pixels in the output mosaic for which not enough data
            # pixels were combined.

            pixel_usage_fraction = self._pixels_counter / len(align_imgs)

            img_handler = pyfits.open(combined_img.path, mode = 'update')
            combined_data = img_handler[0].data
            assert combined_data.shape == pixel_usage_fraction.shape
            combined_data[pixel_usage_fraction < minimum] = 0
            img_handler.close(output_verify = 'ignore')

            return combined_img

        # In case of an error during the update of the mosaic coordinates, an
        # exception will be raised, but we want the output image to be
        # deleted. The inner try block only catches subclasses of exception, as
        # we do not want to catch SystemExit or KeyboardInterrupt here. Also,
        # any exception that may occur during the file deletion is reported,
        # instead of silently dropped. Pythonic code courtesy of Sven Marnach
        # at Stack Overflow [http://stackoverflow.com/q/7687012/184363]

        except:
            e = sys.exc_info()
            try:
                combined_img.unlink()
            except Exception:
                traceback.print_exc()
            raise e[0], e[1], e[2]

    def update_coordinates(self, output_mosaic, scale,
                           ra_keyword = 'RA', dec_keyword = 'DEC'):
        """ Update the coordinate keywords of the image.

        After the aligned images have been combined, the resulting mosaic
        (passed to the method as a FITSImage instance) must have its right
        ascension and declination updated. These values indicate the pointing
        position and are expected to refer to the center (both in the x- and
        y-axes) of the image. Although updating these values is not essential,
        but it will facilitate the astrometry that may be done on the mosaic.

        The coordinates of the center of the mosaic are derived from those of
        the center of the reference image (i.e., the image from which all the
        images that have been combined where shifted), and is actually a
        straight-forward calculation: since we know the width of the margins
        that were added to the reference image and the pixel size ('scale'
        parameter, must be given in arcsec/pixel), we can easily determine the
        coordinates of the mosaic.

        For each axis, if the same number of pixels was added on each side
        (i.e., the width of both margins on the axis is exactly the same), the
        center of the mosaic remains at the same coordinates. In all the other
        cases, the value of the right ascension and declination must be updated
        in the direction of the wider added margin, its value increased or
        decreased by scale x difference arcsecs, where 'difference' represents
        the difference between the size of the widest added margin and that of
        the narrowest one on that axis.

        For example, if the right margin (i.e., on the x-axis, right ascension)
        is 12.31 and the left -3.21, we know that the center of the image has
        moved to the right, as 'more margin' has been added on the right than
        on the left. The value of the right ascension, therefore, becomes
        RA = RA + (12.31 - 3.21) * scale --> RA += (9.1 * scale) arcsecs.

        As a second example, let us consider that the top margin (i.e., on the
        y-axis, declinaton) is 7.89 and the bottom -25.92. In this case, the
        center of the image has shifted to the bottom, since 'more margin' has
        been added at the bottom than at the top. The declination of the
        mosaic, then, would become DEC = DEC - (25.92 - 7.89) * scale -->
        DEC -= (18.03 * scale) arcsecs.

        Note that it is assumed that North is up and East is left. Images for
        which North is down or East is right are not yet supported. Although
        the coordinates of the mosaic are updated in-place, the method returns
        a two-element tuple with the right ascension and declination of the
        mosaic, just in case you are curious about which the new values are.

        Keyword arguments:
        ra_keyword - FITS keyword that stores the right ascension at the center
                     of the mosaic (pointing coordinates), in decimal degrees.
        dec_keyword - FITS keyword that stores the declination at the center of
                     the mosaic (pointing coordinates), in decimal degrees.

        """

        # If the margin on both sides has exactly the same width, the
        # coordinate on that axis (right ascension for the x-axis and
        # declination for the y-axis) remains the same.

        if abs(self.left) != self.right:
            x_diff = self.left + self.right
            right_ascension = output_mosaic.read_keyword(ra_keyword)
            offset_in_dec_degrees = methods.DMS_to_DD(0, 0, x_diff * scale)
            updated_ra = right_ascension + offset_in_dec_degrees
            output_mosaic.update_keyword(ra_keyword, updated_ra)
        else:
            updated_ra = output_mosaic.read_keyword(ra_keyword)

        if abs(self.bottom) != self.top:
            y_diff = self.bottom + self.top
            declination = output_mosaic.read_keyword(dec_keyword)
            offset_in_dec_degrees = methods.DMS_to_DD(0, 0, y_diff * scale)

            # Note that here the offset is subtracted instead of added. The
            # reason for this it that, in our Eurocentric view of the world,
            # East is left. Therefore, as right ascension increases towards the
            # East, 'more margin' on the left means that the right ascension
            # has increased -- subtractig a negative value equals addition.

            updated_dec = declination - offset_in_dec_degrees
            output_mosaic.update_keyword(dec_keyword, updated_dec)
        else:
            updated_dec = output_mosaic.read_keyword(dec_keyword)

        return updated_ra, updated_dec

    def dump(self, directory):
        """ Copy the temporary images to the specified directory.

        This method copies to a directory the temporary images that have been
        aligned during the lifetime of the object. No silent errors here: it is
        your responsibility to provide the path to a directory that exists and
        is writable by the user; otherwise, an exception is raised. Note that
        in case of a name collision files are overwritten, so you should make
        sure that the specified directory is empty or, at least, does not
        contain any data you cannot afford to lose.

        """

        for aligned_img in self.offsets.values():
            dst_path = os.path.join(directory, aligned_img.basename)
            shutil.copy2(aligned_img.path, dst_path)


parser = customparser.get_parser(description)
parser.usage = "%prog [OPTION]... OFFSETS_XML_FILE"

parser.add_option('--overwrite', action = 'store_true', dest = 'overwrite',
                  help = "overwrite output image if it already exists")

parser.add_option('--fraction', action = 'store', type = 'float',
                  dest = 'fraction', default = 0.85,
                  help = "the minimum fraction of data pixels (not "
                  "corresponding to the margins that were added to each image "
                  "before the alignment) that must have been combined for "
                  "each pixel if it is not to be zeroed in the output mosaic "
                  "[default: %default]")

parser.add_option('--name', action = 'store', type = 'str',
                  dest = 'name', default = 'Reference Mosaic',
                  help = "name of the observed object. The OBJECT keyword of "
                  "the output mosaic will be set to this value, as according "
                  "to the FITS Standard this keyword must give 'a name for "
                  "the object observed' [default: %default]")

parser.add_option('--tempdir', action = 'store', type = 'str',
                  dest = 'tempdir', default = None,
                  help = "directory to which to save a copy of the temporary, "
                  "aligned images, which are otherwise removed before the "
                  "program exits. This is primarily useful for debugging "
                  "purposes, or if you do not want to have to blindly trust "
                  "the final mosaic but prefer to be able to examine the "
                  "images from which it was generated and verify that they "
                  "all were correctly aligned")

check_group = optparse.OptionGroup(parser, "Check-image",
              "Instead of the original images, after their alignment, the "
              "SExtractor check-images can be used when composing the mosaic. "
              "This provides a built-in, robust image-reduction feature, "
              "already present in SExtractor and of which we can take "
              "advantage in order to interpolate the background map "
              "('BACKGROUND' check-image type) or subtract the sky "
              "('-BACKGROUND' type), for example. As the generation of the "
              "check-images is a CPU-intensive task, this option is None "
              "by default, meaning that the original images are used.")

check_group.add_option('--check-type', action = 'store', type = 'str',
                       dest = 'check_type', default = None,
                       help = "SExtractor check-image to be used in the "
                       "mosaic for each image.")
parser.add_option_group(check_group)

reject_group = optparse.OptionGroup(parser, "Minmax algorithm parameters",
               "After the alignment is done, all the images are combined "
               "into the mosaic that results from min-max averaging them. "
               "Instead of directly combining all the pixels that are in "
               "the same position, some of them may be rejected: these two "
               "parameters define the fraction of the highest and lowest "
               "that are excluded from the averaging (using IRAF's 'minmax' "
               "rejection algorithm).")

reject_group.add_option('--min', action = 'store', type = 'float',
                        dest = 'flow', default = 0.10,
                        help = "fraction of the lowest pixels to be "
                        "rejected [default: %default]")
reject_group.add_option('--max', action = 'store', type = 'float',
                        dest = 'fhigh', default = 0.10,
                        help = "fraction of the highest pixels to be "
                        "rejected [default: %default]")
parser.add_option_group(reject_group)


key_group = optparse.OptionGroup(parser, "FITS Keywords",
                                 keywords.group_description)

key_group.add_option('--rak', action = 'store', type = 'str',
                     dest = 'rak', default = keywords.rak,
                     help = keywords.desc['rak'])

key_group.add_option('--deck', action = 'store', type = 'str',
                     dest = 'deck', default = keywords.deck,
                     help = keywords.desc['deck'])

key_group.add_option('--objectk', action = 'store', type = 'str',
                     dest = 'objectk', default = keywords.objectk,
                     help = keywords.desc['objectk'])
parser.add_option_group(key_group)
customparser.clear_metavars(parser)

def main(arguments = None):
    """ main() function, encapsulated in a method to allow for easy invokation.

    This method follows Guido van Rossum's suggestions on how to write Python
    main() functions in order to make them more flexible. By encapsulating the
    main code of the script in a function and making it take an optional
    argument the script can be called not only from other modules, but also
    from the interactive Python prompt.

    Guido van van Rossum - Python main() functions:
    http://www.artima.com/weblogs/viewpost.jsp?thread=4829

    Keyword arguments:
    arguments - the list of command line arguments passed to the script.

    """

    if arguments is None:
        arguments = sys.argv[1:] # ignore argv[0], the script name
    (options, args) = parser.parse_args(args = arguments)

    # Print the help and abort the execution if there are fewer than three
    # positional arguments left, as the user must specify at least two FITS
    # images and the output mosaic into which they are assembled.
    if len(args) < 3:
        parser.print_help()
        return 2 # used for command line syntax errors
    else:
        assert len(args) >= 3
        input_paths = set(args[1:-1])
        output_path = args[-1]

    # Refuse to overwrite the output FITS file unless explicitly instructed to
    # do so. Note that, if the --overwritten option is given, we do not need to
    # delete the existing file: it will be silently overwritten when the output
    # of montage.mosaic() is shutil.move()'d to the output path.

    if os.path.exists(output_path):
        if not options.overwrite:
            msg = "%sError. The output file '%s' already exists."
            print msg % (style.prefix, output_path)
            print style.error_exit_message
            return 1

    # montage.mosaic() requires as first argument the directory containing the
    # input FITS images but, in order to maintain the same syntax across all
    # LEMON commands, we receive them as command-line arguments. Thus, create a
    # temporary directory and symlink from it the input images. Hard links are
    # not an option because os.link() will raise "OSError: [Errno 18] Invalid
    # cross-device link" if the temporary directory is created in a different
    # partition.

    pid = os.getpid()
    suffix = "_LEMON_%d_mosaic" % pid
    kwargs = dict(suffix = suffix + '_input')
    input_dir = tempfile.mkdtemp(**kwargs)
    atexit.register(clean_tmp_dir, input_dir)

    for path in input_paths:
        source = os.path.abspath(path)
        basename = os.path.basename(path)
        link_name = os.path.join(input_dir, basename)
        os.symlink(source, link_name)

    # The output of montage.mosaic() is another directory, to which several
    # files are written, so we need the path to a second temporary directory.
    # Delete it before calling mosaic(), as otherwise it will raise IOError
    # ("Output directory already exists").

    kwargs = dict(suffix = suffix + '_output')
    output_dir = tempfile.mkdtemp(**kwargs)
    atexit.register(clean_tmp_dir, output_dir)
    os.rmdir(output_dir)

    montage.mosaic(input_dir, output_dir)

    # montage.mosaic() writes several files to the output directory, but we are
    # only interested in one of them: 'mosaic.fits', the mosaic FITS image. We
    # need to move it to the output path specified by the user.

    MOSAIC_OUTPUT = 'mosaic.fits'
    src = os.path.join(output_dir, MOSAIC_OUTPUT)
    shutil.move(src, output_path)

    print "%sYou're done ^_^" % style.prefix
    return 0

if __name__ == "__main__":
    sys.exit(main())

