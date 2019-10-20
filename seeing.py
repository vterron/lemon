#! /usr/bin/env python

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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from __future__ import division

description = """
This module runs on parallel SExtractor on all the input FITS files, fully
leveraging multiple processors available on the machine. Both images with a too
large full width at half maximum) or elongation (ratio between its semi-major
and semi-minor axis lengths) are automatically discarded by fitting a Gaussian
distribution to the values and excluding those which are above a certain
sigma. Finally, the 'best-seeing' image, the optimal one on which to identify
sources, is found by identifying among the non-discarded images that whose
number of sources is at a given percentile and which has the best FWHM.

"""

import atexit
import logging
import multiprocessing
import numpy
import optparse
import os
import os.path
import scipy.stats
import scipy.signal
import shutil
import sys
import tempfile
import time

# LEMON modules
import util
import astromatic
import customparser
import defaults
import fitsimage
import keywords
import methods
import style

class FITSeeingImage(fitsimage.FITSImage):
    """ High-level interface to the SExtractor catalog of each FITS image.

    As running SExtractor on an image is a CPU-expensive operation, in the
    order of seconds, this class uses the local filesystem as a cache, so that
    SExtractor does not have to be unnecessarily run twice. To achieve this,
    the SExtractor catalog is saved to a temporary file, whose path is stored
    in the FITS header of the image. In order to prevent an outdated catalog
    (that is, if the SExtractor configuration files are modified in between),
    the MD5 hash of the configuration files is also stored.

    """

    def __init__(self, path, maximum, margin, coaddk = keywords.coaddk):
        """ Instantiation method for the FITSeeingImage class.

        The path to the SExtractor catalog is read from the FITS header: if the
        keyword is not found, or if it is present but refers to a non-existent
        file, SExtractor has to be executed again. Even if the catalog exists,
        however, the MD5 of the SExtractor configuration files that were used
        when the catalog was created must be the same.

        The 'maximum' parameter determines the pixel value above which it is
        considered saturated. This value depends not only on the CCD, but also
        on the number of coadded images. If three images were coadded, for
        example, the effective saturation level equals three times the value of
        'saturation'. The number of effective coadds is read from the 'coaddk'
        parameter. If the keyword is missing, a value of one (that is, that the
        image consisted of a single exposure) is assumed.

        The 'margin' argument gives the width, in pixels, of the areas adjacent
        to the edges of the image that are to be ignored when detecting sources
        on the reference image. Stars whose center is fewer than this number of
        pixels from any border of the FITS image are not considered.

        """

        super(FITSeeingImage, self).__init__(path)
        self.margin = margin
        msg = "%s: width of margin: %d pixels" % (self.path, self.margin)
        logging.debug(msg)

        # Compute the MD5 hash of the SExtractor configuration files and this
        # saturation level, which overrides the definition of SATUR_LEVEL.
        satur_level = self.saturation(maximum, coaddk = coaddk)
        options = dict(SATUR_LEVEL = str(satur_level))
        sex_md5sum = astromatic.sextractor_md5sum(options = options)
        msg = "%s: SExtractor MD5 hash: %s" % (self.path, sex_md5sum)
        logging.debug(msg)

        try:

            try:
                self.catalog_path = self.read_keyword(keywords.sex_catalog)
            except KeyError:
                msg = "%s: keyword '%s' not found"
                logging.debug(msg % (self.path, keywords.sex_catalog))
                raise

            # The path not only must be stored: the catalog must also exist
            msg = "%s: on-disk catalog %s" % (self.path, self.catalog_path)
            if not os.path.exists(self.catalog_path):
                logging.debug(msg + " does not exist")
                raise IOError
            else:
                logging.debug(msg)

            try:
                img_sex_md5sum = self.read_keyword(keywords.sex_md5sum)
                msg = "%s: on-disk catalog SExtractor MD5 hash (%s): %s"
                args = self.path, keywords.sex_md5sum, img_sex_md5sum
                logging.debug(msg % args)
            except KeyError:
                msg = "%s: keyword '%s' not found"
                logging.debug(msg % (self.path, keywords.sex_md5sum))
                raise

            if img_sex_md5sum != sex_md5sum:
                msg = "%s: outdated catalog (SExtractor hashes do not match)"
                logging.debug(msg % self.path)
                try:
                    os.unlink(self.catalog_path)
                    msg = "%s: outdated catalog removed from disk" % self.path
                    logging.debug(msg)

                except (IOError, OSError), e:
                    msg = "%s: could not remove outdated catalog (%s)"
                    logging.warning(msg % (self.path, e))

                finally:
                    raise ValueError

            # This point only reached if on-disk catalog can be reused
            msg = "%s: on-disk catalog exists and MD5 hashes match. Yay!"
            logging.debug(msg % self.path)

        except (KeyError, IOError, ValueError):

            msg = ("%s: could not reuse an existing, on-disk cached catalog; "
                   "SExtractor must be run") % self.path
            logging.debug(msg)

             # Redirect standard and error outputs to null device
            with open(os.devnull, 'wt') as fd:
                logging.info("%s: running SExtractor" % self.path)

                self.catalog_path = \
                    astromatic.sextractor(self.path, options = options,
                                          stdout = fd, stderr = fd)

                logging.debug("%s: SExtractor OK" % self.path)

                try:
                    # Update the FITS header with the path and MD5; give up
                    # silently in case we do not have permissions to do it
                    # (IOError) or if the length of the HIERARCH keyword, equal
                    # sign and value is longer than 80 characters (ValueError,
                    # see FITSImage.update_keyword() for details). The cast to
                    # str is needed because PyFITS has complained sometimes
                    # about "illegal values" if it receives a Unicode string.
                    self.update_keyword(keywords.sex_catalog, str(self.catalog_path))
                    self.update_keyword(keywords.sex_md5sum, sex_md5sum)
                except (IOError, ValueError):
                    pass


    @property
    @methods.memoize
    def catalog(self):
        """ Return the SExtraxtor catalog of the FITS image.

        The method returns the SExtractor catalog of the FITS image (a
        astromatic.Catalog instance), excluding from it those stars that are
        too close to the image edges, according to the value 'margin' value
        given when the object was instantiated.

        As parsing the catalog file and constructing the astromatic.Catalog
        instance requires a disk access, the catalog is memoized in order to
        speed up our code. Note that this means that on-disk modifications of
        the catalog (although this is something you should not be doing,
        anyway) will not be reflected after the first call to this method.

        """

        # Catalogs are a subclass of tuple and cannot therefore be modified, so
        # we need to use a list comprehension first, to identify the stars that
        # are not too close to the edges, and use from_sequence() to create a
        # new Catalog instance with the selected stars.

        self._ignored_sources = 0
        catalog = astromatic.Catalog(self.catalog_path)
        logging.info("Removing from the catalog objects too close to the "
                     "edges of %s" % self.path)
        logging.debug("Margin width: %d pixels" % self.margin)
        logging.debug("Image size: (%d, %d)" % self.size)

        non_discarded = []
        for star in catalog:
            logger_msg = "Star at %.3f, %.3f " % (star.x, star.y)
            if star.x < self.margin or \
               star.x > (self.x_size - self.margin) or \
               star.y < self.margin or \
               star.y > (self.y_size - self.margin):
                   logging.debug(logger_msg + "ignored -- too close to edges")
                   self._ignored_sources += 1
            else:
                non_discarded.append(star)
                logging.debug(logger_msg + "OK -- inside of margins")

        if __debug__:
            assert len(catalog) == len(non_discarded) + self._ignored_sources
            for star in non_discarded:
                assert star.x >= self.margin
                assert star.x <= (self.x_size - self.margin)
                assert star.y >= self.margin
                assert star.y <= (self.y_size - self.margin)

        return astromatic.Catalog.from_sequence(*non_discarded)

    def __len__(self):
        """ Return the number of stars detected by SExtractor in the image """
        return len(self.catalog)

    @property
    def ignored(self):
        """ Return the number of SExtractor detections ignored because they
        were too close to the image edges"""
        try:
            return self._ignored_sources
        # The _ignored_sources attribute does not exist until the SExtractor
        # catalog is accessed the first time; so if this method is called
        # before it happens, just access (and memoize) it and try again
        except AttributeError:
            self.catalog
            return self.ignored

    @property
    def total(self):
        """ Return the total number of detections in the SExtractor catalog,
        ignored or not"""
        return len(self) + self.ignored

    def __getitem__(self, key):
        """ Return the key-th star detected by SExtractor in the image --
        obviously, ignored stars (too close to the edges) are not considered"""
        return self.catalog[key]

    @property
    def coordinates(self):
        """ Return a list with the celestial coordinates of the stars.

        This method returns a list containing a Coordinates object for each
        astronomical source that was detected by SExtractor and not ignored
        because of its proximity to the image edge.

        """
        return self.catalog.get_sky_coordinates()

    def snr_percentile(self, per):
        """ Return the score at the given percentile of the SNR of the stars.

        This methods answers the question of what is the lowest SNR among a
        given fraction of the stars with the best signal-to-noise ratio. Those
        stars whose ADU counts are above the effective saturation level are
        excluded from the computation. If the desired quantile lies between two
        data points, we interpolate between them. Note that the score at
        per = 50 is, by definition, the median.

        The ValueError exception is raised in case no star can be used to
        compute the percentile of the signal-to-noise ratio. This is usually
        caused by a too low saturation level, although it may also happen with
        small images of sparsely populated fields, for which the value of the
        margin may cause all the stars to be ignored.

        Arguments:
        per - percentile at which to extract the score

        """

        snrs = [star.snr for star in self if not star.saturated]
        if not snrs:
            raise ValueError("no stars available to compute the percentile SNR of '%s'" % self.path)

        # If empty, scoteatpercentile would raise IndexError
        return scipy.stats.scoreatpercentile(snrs, per)

    def fwhm(self, per = 50, mode = 'median'):
        """ Return the median (or mean) FWHM of the stars in the image.

        The method returns the median (or mean) of the full width at half
        maximum of the stars detected by SExtractor in the image. Those stars
        whose ADU counts are above the effective saturation level are excluded
        from the computation. Also are excluded stars with too low a SNR (see
        the 'per' keyword argument), in order to prevent false detections and
        faint objets from distorting the reliability of the result.

        The ValueError exception is raised in case no star can be used to
        compute the FWHM. This is usually caused by a too low saturation level,
        although it may also happen with small images of sparsely populated
        fields, for which the value of the margin may cause all the stars to
        be ignored.

        The FWHM is computed thanks to the FLUX_RADIUS parameter, which
        estimates the radius of the circle centered on the barycenter of the
        star that encloses about half of the total flux. For a Gaussian profile,
        this is equal to 1/2 FWHM, although, according to the author of
        SExtractor, for most astronomical images it will be slightly higher.
        [URL] http://www.astromatic.net/forum/showthread.php?tid=318

        Keyword arguments:
        per - the score at the given percentile of the signal-to-noise ratio of
              the stars that are to be included in the calculation of the
              FWHM. That is: stars whose SNR is below the 'per' percentile are
              excluded when the FWHM of the stars in the image is calculated.
              Note that the percentile is calculated taking into account only
              the stars within the image margins.

        mode - how the FWHM of the image is computed out of the FWHM of the
               stars on it detected. Must by 'median' or 'mean'; otherwise,
               ValueError is raised.

        """

        if mode not in ('median', 'mean'):
            raise ValueError("'mode' must be 'median' or 'mean'")

        # Find the signal-to-noise ratio of the image at the 'per' percentile.
        # This will be the minimum SNR that a star must have for its FWHM to be
        # taken into account when calculating the FWHM of the image as a whole.
        snr = self.snr_percentile(per)

        fwhms = []
        for star in self:
            # Ignore saturated and noisy stars
            if star.saturated or star.snr < snr:
                continue
            else:
                fwhms.append(star.fwhm)

        if not fwhms:
            # Exception needed, NumPy would return NaN for an empty list
            raise ValueError("no stars available to compute the FWHM")

        if mode == 'median':
            return numpy.median(fwhms)
        else:
            assert mode == 'mean'
            return numpy.mean(fwhms)

    def elongation(self, per = 50, mode = 'median'):
        """ Return the median (or mean) elongation of the stars in the image.

        The method returns the median (or mean) of the elongation of the stars
        detected by SEXtractor in the image. Those stars whose ADU counts are
        above the effective saturation level are excluded from the computation.
        Also are excluded stars with too low a SNR (see the 'per' keyword
        argument), in order to prevent false detections and faint objets from
        distorting the reliability of the result.

        The ValueError exception is raised in case no star can be used to
        compute the elongation. This is usually caused by a too low saturation
        level, although it may also happen with small images of sparsely
        populated fields, for which the value of the margin may cause all
        the stars to be ignored.

        The elongation of an object is given in the SExtractor catalogs by the
        ELONGATION parameter, and is defined as A/B, where A and B are its
        semi-major and semi-minor axis lengths, respectively. More precisely,
        A and B represent the maximum and minimum spatial rms of the object
        profile along any direction. Therefore, the average elongation is
        always >= 1.0, as A is known to be greater or equal than B.

        Keyword arguments:
        per - the score at the given percentile of the signal-to-noise ratio of
              the stars that are to be included in the calculation of the
              elongation. That is: stars whose SNR is below the 'per'
              percentile are excluded when the elongation of the stars in the
              image is calculated.  Note that the percentile is calculated
              taking into account only the stars within the image margins.

        mode - how the elongation of the image is computed out of the
               elongation of the stars on it detected. Must by 'median' or
               'mean'; otherwise, ValueError is raised.

        """

        if mode not in ('median', 'mean'):
            raise ValueError("'mode' must be 'median' or 'mean'")

        # Find the signal-to-noise ratio of the image at the 'per' percentile.
        snr = self.snr_percentile(per)

        elongations = []
        for star in self:
            # Ignore saturated and noisy stars
            if star.saturated or star.snr < snr:
                continue
            else:
                elongations.append(star.elongation)

        if not elongations:
            # Exception needed, NumPy would return NaN for an empty list
            raise ValueError("no stars available to compute the elongation")

        if mode == 'median':
            return numpy.median(elongations)
        else:
            assert mode == 'mean'
            return numpy.mean(elongations)

# This Queue is global -- this works, but note that we could have
# passed its reference to the function managed by pool.map_async.
# See http://stackoverflow.com/a/3217427/184363
queue = util.Queue()

parser = customparser.get_parser(description)
parser.usage = "%prog [OPTION]... INPUT_IMGS... OUTPUT_DIR"
parser.add_option('--filename', action = 'store', type = 'str',
                  dest = 'bseeingfn', default = 'best_seeing.fits',
                  help = "filename with which the FITS image with the best "
                  "full width at half-maximum (FWHM) will be saved to "
                  "OUTPUT_DIR. If set to an empty string (''), the image "
                  "will retain its original filename, although --suffix will "
                  "still be appended (in other words, it is treated as all "
                  "the other images) [default: %default]")

parser.add_option('--maximum', action = 'store', type = 'int',
                  dest = 'maximum', default = defaults.maximum,
                  help = defaults.desc['maximum'])

parser.add_option('--margin', action = 'store', type = 'int',
                  dest = 'margin', default = defaults.margin,
                  help = defaults.desc['margin'])

parser.add_option('--snr-percentile', action = 'store', type = 'float',
                  dest = 'per', default = defaults.snr_percentile,
                  help = defaults.desc['snr_percentile'])

parser.add_option('--mean', action = 'store_true', dest = 'mean',
                  help = defaults.desc['mean'])

parser.add_option('--sources-percentile', action = 'store', type = 'float',
                  dest = 'stars_per', default = 75.0,
                  help = "the best-seeing image is identified by finding the "
                  "one with the best (lowest) FWHM among the non-discarded "
                  "images whose number of detected sources is at this "
                  "percentile. Taking directly the image with the best FWHM "
                  "could not work as we need the image to also be one of the "
                  "most populated [default: %default]")

parser.add_option('--suffix', action = 'store', type = 'str',
                  dest = 'suffix', default = 's',
                  help = "string to be appended to output images, before "
                  "the file extension, of course [default: %default]")

parser.add_option('--overwrite', action = 'store_true', dest = 'overwrite',
                  help = "overwrite any output file if it already exists")

parser.add_option('--cores', action = 'store', type = 'int',
                  dest = 'ncores', default = defaults.ncores,
                  help = defaults.desc['ncores'])

parser.add_option('-v', '--verbose', action = 'count',
                  dest = 'verbose', default = defaults.verbosity,
                  help = defaults.desc['verbosity'])

fwhm_group = optparse.OptionGroup(parser, "Full width at half maximum",
             "After SExtractor is run on the FITS files, the FWHM of each "
             "image is computed by taking the median or mean of the stars "
             "which center is not too close to the margins see --mean and "
             "--margin options) and fitting a Gaussian distribution to the "
             "values. Those images whose FWHM exceed the specified number "
             "of standard deviations of the mean (mu + n x sigma) are "
             "discarded")

fwhm_group.add_option('--fsigma', action = 'store', type = 'float',
                      dest = 'fwhm_sigma', default = 3.0,
                      help = "the maximum number of standard deviations of "
                      "the mean; gives the maximum FWHM allowed for images "
                      "[default: %default]")

fwhm_group.add_option('--fwhm_dir', action = 'store', type = 'str',
                      dest = 'fwhm_dir', default = 'fwhm_discarded',
                      help = "subdirectory where automatically discarded "
                      "images because of their FWHM will be saved in "
                      "OUTPUT_DIR [default: %default]")
parser.add_option_group(fwhm_group)

elong_group = optparse.OptionGroup(parser, "Elongation",
              "After images with a bad (excessively large) FWHM are "
              "discarded, those elongated are identified and excluded by "
              "using the same approach. The elongation of a star, given by "
              "SExtractor, is defined as A/B, where A and B are its "
              "semi-major and semi-minor axis lengths, respectively. More "
              "precisely, A and B represent the maximum and minimum spatial "
              "rms of the object profile along any direction. Elongated "
              "images are usually caused by intermittent problems with the "
              "telescope tracking.")

elong_group.add_option('--esigma', action = 'store', type = 'float',
                       dest = 'elong_sigma', default = 3.0,
                       help = "the maximum number of standard deviations of "
                       "the mean; gives the maximum elongation allowed for "
                       "images [default: %default]")

elong_group.add_option('--elong_dir', action = 'store', type = 'str',
                       dest = 'elong_dir', default = 'elong_discarded',
                       help = "subdirectory where automatically discarded "
                       "images because of their elongation will be saved in "
                       "OUTPUT_DIR [default: %default]")
parser.add_option_group(elong_group)

key_group = optparse.OptionGroup(parser, "FITS Keywords",
                                 keywords.group_description)

key_group.add_option('--coaddk', action = 'store', type = 'str',
                     dest = 'coaddk', default = keywords.coaddk,
                     help = keywords.desc['coaddk'])

key_group.add_option('--fwhmk', action = 'store', type = 'str',
                     dest = 'fwhmk', default = keywords.fwhmk,
                     help = "keyword to which to write the estimated FWHM "
                     "of of the image. This value is needed by subsequent "
                     "modules of the pipeline [default: %default]")
parser.add_option_group(key_group)
customparser.clear_metavars(parser)

@methods.print_exception_traceback
def parallel_sextractor(args):
    """ Run SExtractor and compute the FWHM and elongation of a FITS image.

    This method is intended to be used with a multiprocessing' pool of workers.
    It receives a two-element tuple, the path of an image and the 'instance'
    object returned by optparse's parse_args(), containing the values for all
    the options of the program, and runs SExtractor on the image. The median or
    mean (depending whether the --mean option was given) FWHM and elongation of
    the stars in the image is also computed. Nothing is returned; instead, the
    result is saved to the global variable 'queue' as a five-element tuple: (1)
    path of the input image, (2) path of the temporary output image (a copy of
    the input image, with the path to the SExtractor catalog and the MD5 hash
    of the configuration files written to its header), (3) FWHM, (4) elongation
    and (5) number of objects that were detected by SExtractor. Nothing is
    added to 'queue' in case an error is encountered.

    """

    path, options = args
    mode = options.mean and 'mean' or 'median'

    try:

        # FITSeeingImage.__init__() writes to the header of the FITS image the
        # path to the SExtractor catalog and the MD5 hash of the configuration
        # files. Work on a copy of the input image so as not to modify it.
        kwargs = dict(prefix = '%s_' % os.path.basename(path),
                      suffix = '.fits')
        fd, output_path = tempfile.mkstemp(**kwargs)
        os.close(fd)
        shutil.copy2(path, output_path)

        # Allow FITSeeingImage.__init__() to write to the FITS header
        util.owner_writable(output_path, True) # chmod u+w

        args = output_path, options.maximum, options.margin
        kwargs = dict(coaddk = options.coaddk)
        image = FITSeeingImage(*args, **kwargs)
        fwhm = image.fwhm(per = options.per, mode = mode)
        logging.debug("%s: FWHM = %.3f" % (path, fwhm))
        elong = image.elongation(per = options.per, mode = mode)
        logging.debug("%s: Elongation = %.3f" % (path, elong))
        nstars = len(image)
        logging.debug("%s: %d sources detected" % (path, nstars))
        queue.put((path, output_path, fwhm, elong, nstars))

    except fitsimage.NonStandardFITS:
        logging.info("%s ignored (non-standard FITS)" % path)
        return

    # Raised if no star can be used to compute the FWHM or elongation
    except ValueError, e:
        logging.info("%s ignored (%s)" % (path, str(e)))

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

    # Adjust the logger level to WARNING, INFO or DEBUG, depending on the
    # given number of -v options (none, one or two or more, respectively)
    logging_level = logging.WARNING
    if options.verbose == 1:
        logging_level = logging.INFO
    elif options.verbose >= 2:
        logging_level = logging.DEBUG
    logging.basicConfig(format = style.LOG_FORMAT, level = logging_level)

    # Print the help and abort the execution if there are not two positional
    # arguments left after parsing the options, as the user must specify at
    # least one (only one?) input FITS file and the output directory
    if len(args) < 2:
        parser.print_help()
        return 2     # 2 is generally used for command line syntax errors
    else:
        input_paths = args[:-1]
        output_dir = args[-1]

    # Make sure that the output directory exists, and create it if it doesn't.
    # The subdirectories for discarded images are not yet created; we put this
    # off until we know that at least one image is indeed going to be excluded.
    util.determine_output_dir(output_dir)
    fwhm_dir = os.path.join(output_dir, options.fwhm_dir)
    elong_dir = os.path.join(output_dir, options.elong_dir)

    print "%s%d paths given as input, on which sources will be detected." % \
          (style.prefix, len(input_paths))
    print "%sRunning SExtractor on all the FITS images..." % style.prefix

    # Use a pool of workers and run SExtractor on the images in parallel!
    pool = multiprocessing.Pool(options.ncores)
    map_async_args = ((path, options) for path in input_paths if os.path.isfile(path))
    result = pool.map_async(parallel_sextractor, map_async_args)

    util.show_progress(0.0)
    while not result.ready():
        time.sleep(1)
        util.show_progress(queue.qsize() / len(input_paths) * 100)
        # Do not update the progress bar when debugging; instead, print it
        # on a new line each time. This prevents the next logging message,
        # if any, from being printed on the same line that the bar.
        if logging_level < logging.WARNING:
            print

    result.get()      # reraise exceptions of the remote call, if any
    util.show_progress(100) # in case the queue was ready too soon
    print

    # Three sets, to keep the track of all the images on which SExtractor
    # has been run and also of which have been discarded because of their
    # unnaceptable FWHM or elongation ratio.
    all_images = set()
    fwhm_discarded = set()
    elong_discarded = set()

    # Dictionary mapping each input image to the temporary output file: a copy
    # of the input image but whose FITS header has been updated with the path
    # to the SExtractor catalog and the MD5 hash of the configuration files.
    seeing_tmp_paths = dict()

    # Extract the four-element tuples (path to the image, FWHM, elongation and
    # number of sources detected by SExtractor) from the multiprocessing' queue
    # and store the values in three independent dictionaries; these provide
    # fast access, with O(1) lookup, to the data.
    fwhms  = {}
    elongs = {}
    nstars = {}

    for _ in xrange(queue.qsize()):
        path, output_tmp_path, fwhm, elong, stars = queue.get()
        all_images.add(path)
        seeing_tmp_paths[path] = output_tmp_path

        # The clean-up function cannot be registered in parallel_sextractor()
        # because it would remove the temporary FITS file when the process
        # terminates (instead of when our program exits, which is what we
        # need). Do it here, to make sure that whatever happens next these
        # temporary files are always deleted.
        atexit.register(util.clean_tmp_files, output_tmp_path)

        fwhms[path]  = fwhm
        elongs[path] = elong
        nstars[path] = stars

    if not all_images:
        print "%sError. No FITS images were detected." % style.prefix
        print style.error_exit_message
        return 1

    # Let's first discard those images with a bad full width at half maximum.
    # In order to to this, we fit a normal distribution (assuming the FWHMs to
    # be Gaussian distributed) and define the maximum allowed value as that
    # which exceeds the specified number of standard deviations of the mean.

    print "%sFitting a Gaussian distribution to the FWHMs..." % style.prefix ,
    sys.stdout.flush()
    logging.debug("Fitting a Gaussian distribution to the %d FWHMs" % len(fwhms))
    mu, sigma = scipy.stats.norm.fit(fwhms.values())
    logging.debug("FWHMs mean = %.3f" % mu)
    logging.debug("FWHMs sigma = %.3f" % sigma)
    print 'done.'
    sys.stdout.flush()

    print "%sFWHMs mean = %.3f, sigma = %.3f pixels" % (style.prefix, mu, sigma)
    maximum_fwhm = mu + (options.fwhm_sigma * sigma)
    logging.debug("Maximum allowed FWHM = %.3f + %.1f x %.3f = %.3f pixels" % \
                 (mu, options.fwhm_sigma, sigma, maximum_fwhm))
    print "%sDiscarding images with a FWHM > %.3f + %.1f x %.3f = %.3f pixels..." % \
          (style.prefix, mu, options.fwhm_sigma, sigma, maximum_fwhm)

    # Exclude images by adding them to the FWHM-discarded set
    for path, fwhm in sorted(fwhms.iteritems()):
        if fwhm > maximum_fwhm:
            fwhm_discarded.add(path)
            logging.debug("%s discarded (FWHM = %.3f > %.3f" % \
                         (path, fwhm, maximum_fwhm))
            print "%s%s discarded (FWHM = %.3f)" % (style.prefix, path, fwhm)

    logging.info("Images discarded by FWHM: %d" % len(fwhm_discarded))
    if not fwhm_discarded:
        print "%sNo images were discarded because of their FWHM. Hooray!" % style.prefix
    else:
        discarded_fraction = len(fwhm_discarded) / len(all_images) * 100
        nleft = len(all_images) - len(fwhm_discarded)  # non-discarded images
        print "%s%d FITS images (%.2f %%) discarded, %d remain" % \
              (style.prefix,  len(fwhm_discarded), discarded_fraction, nleft)


    # Repeat the same approach, now with the elongation ratios. Images already
    # discarded because of their FWHM are not even considered -- why discard
    # them twice? They can simply be ignored.

    print "%sFitting a Gaussian distribution to the elongations..." % style.prefix ,
    sys.stdout.flush()
    mu, sigma = scipy.stats.norm.fit(elongs.values())
    logging.debug("Elongations mean = %.3f" % mu)
    logging.debug("Elongations sigma = %.3f" % sigma)
    print 'done.'
    sys.stdout.flush()

    print "%sElongation mean = %.3f, sigma = %.3f pixels" % (style.prefix, mu, sigma)
    maximum_elong = mu + (options.elong_sigma * sigma)
    logging.debug("Maximum allowed elongation = %.3f + %.1f x %.3f = %.3f pixels" % \
                 (mu, options.elong_sigma, sigma, maximum_elong))
    print "%sDiscarding images with an elongation > %.3f + %.1f x %.3f = %.3f ..." % \
          (style.prefix, mu, options.elong_sigma, sigma, maximum_elong)

    for path, elong in sorted(elongs.iteritems()):
        # Ignore FWHM-discarded images
        if path in fwhm_discarded:
            logging.debug("%s ignored (already discarded by FWHM)" % path)
            continue
        elif elong > maximum_elong:
            elong_discarded.add(path)
            logging.debug("%s discarded (elongation = %.3f > %.3f" % \
                         (path, fwhm, maximum_elong))
            print "%s%s discarded (elongation = %.3f)" % (style.prefix, path, elong)

    logging.info("Images discarded by elongation: %d" % len(elong_discarded))
    if not elong_discarded:
        print "%sNo images were discarded because of their elongation. Yay!" % style.prefix
    else:
        initial_size = len(all_images) - len(fwhm_discarded)
        discarded_fraction = len(elong_discarded) / initial_size * 100
        nleft = initial_size - len(elong_discarded)
        print "%s%d FITS images (%.2f %%) discarded, %d remain" % \
              (style.prefix,  len(elong_discarded), discarded_fraction, nleft)


    # Finally, take the images whose number of stars is at the 'stars_per'
    # percentile and select the one with the best FWHM. This will be our
    # 'best-seeing' image, in which sources may be detected. Taking directly
    # the image with the best FWHM may not work as we need the best-seeomg
    # image to also be one of the most populated.

    print "%sIdentifying the images whose number of detected sources it at " \
          "the %.2f percentile..." % (style.prefix, options.stars_per) ,
    sys.stdout.flush()
    # Ignore discarded images, for whatever reason
    logging.debug("Finding the %.2f percentile of the number of stars " \
                 "detected by SExtractor"  % options.stars_per)
    for path in fwhm_discarded.union(elong_discarded):
        del nstars[path]
        reason = 'FWHM' if path in fwhm_discarded else 'elongation'
        logging.debug("%s ignored (was discarded by %s)" % (path, reason))
    min_nstars = scipy.stats.scoreatpercentile(nstars.values(), options.stars_per)
    print 'done.'

    print "%sNumber of stars at percentile = %d, taking the images with at " \
          "least this number of sources..." % (style.prefix, min_nstars) ,
    sys.stdout.flush()
    most_populated_images = [path
                             for path, stars in nstars.iteritems()
                             if stars >= min_nstars]

    logging.debug("There are %s images with a number of stars at the %.2f " \
                 "percentile" % (len(most_populated_images), options.stars_per))
    logging.debug("Identifying the image with the lowest FWHM")
    print 'done.'

    print "%sFinally, finding the image with the lowest FWHM among these " \
          "%d images..." % (style.prefix, len(most_populated_images)),
    sys.stdout.flush()

    # Find the image with the best seeing (lowest FWHM)
    best_seeing = min(most_populated_images, key = lambda path: fwhms[path])
    logging.debug("Best-seeing image: %s" % path)
    logging.debug("Best-seeing image FWHM = %.3f" % fwhms[best_seeing])
    logging.debug("Best-seeing image elongation = %.3f" % elongs[best_seeing])
    logging.debug("Best-seeing image sources = %d" % nstars[best_seeing])
    assert best_seeing not in fwhm_discarded
    assert best_seeing not in elong_discarded
    print 'done.'

    print "%sBest-seeing image = %s, with %d sources and a FWHM of %.3f pixels" % \
          (style.prefix, best_seeing, nstars[best_seeing], fwhms[best_seeing])

    # The subdirectories are created only if at least one image is going to be
    # discarded. We do not want empty directories in case no image is discarded
    # because of its full-width at half maximum (FWHM) or elongation.

    if fwhm_discarded:
        util.determine_output_dir(fwhm_dir, quiet = True)

    if elong_discarded:
        util.determine_output_dir(elong_dir, quiet = True)

    # Finally, copy all the FITS images to the output directory
    processed = 0
    for path in sorted(all_images):
        # Add the suffix to the basename of the FITS image
        root, ext = os.path.splitext(os.path.basename(path))
        output_filename = root + options.suffix + ext
        logging.debug("Basename '%s' + '%s' becomes '%s'" % \
                     (path, options.suffix, output_filename))

        if path in fwhm_discarded:
            output_path = os.path.join(fwhm_dir, output_filename)
            logging.debug("%s was discarded because of its FWHM" % path)
            logging.debug("%s to be copied to subdirectory %s" % (path, fwhm_dir))
            history_msg1 = "Image discarded by LEMON on %s" % methods.utctime()
            history_msg2 = "[Discarded] FWHM = %.3f pixels, maximum allowed value = %.3f" % \
                           (fwhms[path], maximum_fwhm)

        elif path in elong_discarded:
            output_path = os.path.join(elong_dir, output_filename)
            logging.debug("%s was discarded because of its elongation ratio" % path)
            logging.debug("%s to be copied to subdirectory %s" % (path, elong_dir))
            history_msg1 = "Image discarded by LEMON on %s" % methods.utctime()
            history_msg2 = "[Discarded] Elongation = %.3f, maximum allowed value = %.3f" % \
                           (elongs[path], maximum_elong)

        elif path == best_seeing:

            # Retain original name if --filename is an empty string
            if not options.bseeingfn:
                filename = output_filename
            else:
                filename = options.bseeingfn

            output_path = os.path.join(output_dir, filename)
            logging.debug("%s is the best-seeing image" % path)
            logging.debug("%s to be copied to directory %s with name %s" % \
                         (path, output_dir, options.bseeingfn))
            history_msg1 = "Image identified by LEMON as the 'best-seeing' one"
            history_msg2 = "FWHM = %.3f | Elongation = %.3f | Sources: %d (at %.2f percentile)" % \
                           (fwhms[path], elongs[path], nstars[path], options.stars_per)

        else:
            output_path = os.path.join(output_dir, output_filename)
            logging.debug("%s to be copied to %s" % (path, output_dir))
            history_msg1 = "Image FWHM = %.3f" % fwhms[path]
            history_msg2 = "Image elongation = %.3f" % elongs[path]

        if os.path.exists(output_path) and not options.overwrite:
            msg = ("%sError. Output FITS file '%s' already exists. "
                   "You need to use --overwrite.")
            args = style.prefix, output_path
            print msg % args
            print style.error_exit_message
            return 1

        else:
            src = seeing_tmp_paths[path]
            shutil.move(src, output_path)

        util.owner_writable(output_path, True) # chmod u+w
        logging.debug("%s copied to %s" % (path, output_path))
        output_img = fitsimage.FITSImage(output_path)
        output_img.add_history(history_msg1)
        output_img.add_history(history_msg2)
        logging.debug("%s: FITS header updated (HISTORY keywords)" % path)

        # Copy the FWHM to the FITS header, for future reference
        comment = "Margin = %d, SNR percentile = %.3f" % (options.margin, options.per)
        output_img.update_keyword(options.fwhmk, fwhms[path], comment = comment)
        logging.debug("%s: FITS header updated (%s keyword)" % (path, options.fwhmk))

        print "%sFITS image %s saved to %s" % (style.prefix, path, output_path)
        processed += 1

    print "%sA total of %d images was saved to directory '%s'." % (style.prefix, processed, output_dir)
    print "%sWe're done ^_^" % style.prefix
    return 0

if __name__ == "__main__":
    sys.exit(main())
