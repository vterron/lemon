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

import collections
import functools
import hashlib
import os
import os.path
import re
import math
import shutil
import tempfile
import subprocess

# LEMON modules
import fitsimage
import methods

ASTROMATIC_FILES = os.path.join(os.path.dirname(__file__), 'astromatic/')
SEXTRACTOR_CONFIG = os.path.join(ASTROMATIC_FILES, 'sextractor.sex')
SEXTRACTOR_PARAMS = os.path.join(ASTROMATIC_FILES, 'sextractor.param')
SEXTRACTOR_FILTER = os.path.join(ASTROMATIC_FILES, 'sextractor.conv')
SEXTRACTOR_STARNNW = os.path.join(ASTROMATIC_FILES, 'sextractor.nnw')
SEXTRACTOR_COMMANDS = 'sextractor', 'sex' # may be any of these
SEXTRACTOR_REQUIRED_VERSION = (2, 8, 6)

SCAMP_CONFIG = os.path.join(ASTROMATIC_FILES, 'scamp.conf')
SCAMP_AHEADER_SUFFIX = '.ahead'
SCAMP_HEADER_SUFFIX = '.head'
SCAMP_COMMAND = 'scamp'
CDSCLIENT_COMMAND = 'aclient'

SWARP_CONFIG = os.path.join(ASTROMATIC_FILES, 'swarp.conf')
SWARP_COMMAND = 'swarp'

class SExtractorNotInstalled(StandardError):
    pass

class SExtractorUpgradeRequired(StandardError):
    """ Raised if a too-old version of SExtractor is installed """
    pass

class SExtractorError(subprocess.CalledProcessError):
    pass

class Pixel(collections.namedtuple('Pixel', "x y")):
    """ A pair of immutable x- and y-coordinates. """

    def distance(self, another):
        """ Return the Euclidean distance between two Pixels """

        x_axis = pow(self.x - another.x, 2)
        y_axis = pow(self.y - another.y, 2)
        return math.sqrt(x_axis + y_axis)

class Star(collections.namedtuple('Pixel', "img_coords, sky_coords, area, "
           "mag, saturated, snr, fwhm, elongation")):
    """ An immutable class with a source detected by SExtractor. """

    def __new__(cls, x, y, alpha, delta, area, mag, satur, snr, fwhm, elong):
        """ Customize the creation of a Star instance: encapsulate the (x, y)
        and (alpha, delta) as Pixel objects and pass them as the first two
        arguments of the named tuple. The other arguments are not modified.

        x - star position along x.
        y - star position along y.
        alpha - right ascension of the star.
        delta - declination of the star.
        area - isophotal area (filtered) above detection threshold (pix^2).
        mag - measure of the brightness as seen by an observer on Earth.
        satur - at least one pixel of the Star is saturated, or very close to.
        snr - the signal-to-noise ratio of the star. This is the only value
              not directly read from the SExtractor catalog, but which has
              to be derived by us from other values.
        fwhm - the full width at half maximum (FWHM) of the star.
        elong - the value of A/B, where A and B are its semi-major and
                semi-minor axis lengths, as reported by SExtractor. More
                precisely, A and B represent the maximum and minimum spatial
                rms of the object profile along any direction.

        """

        img_coords = Pixel(x, y)
        sky_coords = Pixel(alpha, delta)
        args = img_coords, sky_coords, area, mag, satur, snr, fwhm, elong
        return super(Star, cls).__new__(cls, *args)

    @property
    def x(self):
        return self.img_coords.x

    @property
    def y(self):
        return self.img_coords.y

    @property
    def alpha(self):
        return self.sky_coords.x

    @property
    def delta(self):
        return self.sky_coords.y

    def angular_distance(self, another):
        """ Return the angular distance, in degrees, between two Stars. """

        # Formula: cos(A) = sin(d1)sin(d2) + cos(d1)cos(d2)cos(ra1-ra2)
        # http://www.astronomycafe.net/qadir/q1890.html

        ra1  = math.radians(self.alpha)
        dec1 = math.radians(self.delta)
        ra2  = math.radians(another.alpha)
        dec2 = math.radians(another.delta)
        return math.degrees(math.acos(math.sin(dec1) * math.sin(dec2) +
                                      math.cos(dec1) * math.cos(dec2) *
                                      math.cos(ra1-ra2)))

    def distance(self, another):
        """ The Euclidean distance between the image coordinates of two Stars """
        return self.img_coords.distance(another.img_coords)

class Catalog(tuple):
    """ High-level interface to a SExtractor catalog """

    @staticmethod
    def _find_column(contents, parameter):
        """ Return the index of a SExtractor paramater in the catalog.

        The method takes as inputs the contents of a SExtractor catalog and the
        name of a parameter (such as 'X_IMAGE' or 'FLUX_MAX') and returns the
        zero-based index of the column in which it is located. As happens with
        list.index(), ValueError is raised if the parameter cannot be found in
        the catalog. For user's convenience the search is case-insensitive, so
        'x_world' would match 'X_WORLD'.

        The 'contents' parameter must be a list (one element per line in the
        SExtractor catalog) of lists (one element per word in each line, as
        returned by str.split(). This is how 'content' could look like, e.g.:
        [['#', '1', 'X_IMAGE', 'Object', 'position', 'along', 'x', '[pixel]'],
         ['#', '2', 'Y_IMAGE', 'Object', 'position', 'along', 'y', '[pixel]'],
         ...]

        Is is important no note that, for the parameter being able to be found,
        the catalog must have been saved in the SExtractor ASCII_HEAD format,
        so that the file stars with comment lines listing column labels.

        """

        # We need to examine the comment lines (those that start with a '#')
        # and extract the value of the column associated to the parameter.
        # SExtractor comments follow the following format:
        #
        #       # 1 X_IMAGE     Object position along x     [pixel]
        #       # 2 Y_IMAGE     Object position along y     [pixel]
        #
        # The first integer in each line, right after the '#', indicates the
        # column of the parameter. Note that we must subtract one from these
        # indexes, as they are one-based.

        for line in contents:
            if line[0][0] == '#':
                param_name = line[2]
                if param_name.upper() == parameter.upper():
                    param_index = int(line[1]) - 1
                    return param_index
        else:
            msg = "parameter '%s' not found" % parameter
            raise ValueError(msg)

    @staticmethod
    def flag_saturated(flag_value):
        """ Test the value of FLAGS and determine if the object has saturated.

        The method receives the value of the internal flag computed for a star
        and returns True if the decimal value of the flag indicates that at
        least one pixel of the object is saturated or very close to. Otherwise,
        False is returned.

        [From the SExtractor user guide, page 25] The internal flags are always
        computed. They are accessible through the FLAGS catalog parameter, which
        is a short integer. FLAGS contains, coded in decimal, all the extraction
        flags as a sum of powers of 2:

        1 - The object has neighbours, bright and close enough to significantly
        bias the MAG AUTO photometry , or bad pixels (more than 10% of the
        integrated area affected),
        2 - The object was originally blended with another one,
        4 - At least one pixel of the object is saturated (or very close to),
        8 - The object is truncated (too close to an image boundary),
        16 - Object's aperture data are incomplete or corrupted,
        36 - Object's isophotal data are incomplete or corrupted,
        64 - A memory overflow occurred during deblending,
        128 - A memory overflow occurred during extraction.

        For example, an object close to an image border may have FLAGS = 16,
        and perhaps FLAGS = 8+16+32 = 56. [End of quote]

        A flag is saturated, therefore, if 4 was one of the values that were
        added to calculate it. In order to test this, the flag is represented
        as a binary string and we check whether the third bit (2**(3-1) == 4)
        from the left equals one.

        Since the value of the flag is determined by the first eight powers of
        two, its minimum valid value is zero and the maximum (2**8)-1 = 255.
        The ValueError exception is raised if the decimal value of the flag
        is outside of this range.

        """

        if not 0 <= flag_value <= 255:
            msg = "flag value out of range [0, 255]"
            raise ValueError(msg)

        # Convert from, for example, '0b1' to '1', and then fill in with as
        # many zeros are needed to represent the flag as an eight-bit binary
        # number ('00000001'), so that we can always check the value of the
        # third least significant bit.
        return int(bin(flag_value)[2:].zfill(8)[-3]) == 1

    @classmethod
    def _load_stars(cls, path):
        """ Load a SExtractor catalog into memory.

        The method parses a SExtractor catalog and returns a generator of Star
        objects, once for each detected object. It is mandatory, or ValueError
        will be raised otherwise, that the following parameters are present in
        the catalog: X_IMAGE, Y_IMAGE, ALPHA_SKY, DELTA_SKY, ISOAREAF_IMAGE,
        MAG_AUTO, FLUX_ISO, FLUXERR_ISO, FLUX_RADIUS and ELONGATION. Also, the
        catalog must have been saved in the SExtractor ASCII_HEAD format, as
        the comment lines listing column labels are needed in order to detect
        in which column each parameter is.

        The FWHM is derived from the FLUX_RADIUS parameter, which estimates the
        radius of the circle centered on the barycenter that encloses about
        half of the total flux. For a Gaussian profile, this is equal to 1/2
        FWHM, although with most astronomical images it will be slightly
        higher [http://www.astromatic.net/forum/showthread.php?tid=318]

        The signal-to-noise ratio is calculated as FLUX_ISO / FLUXERR_ISO,
        that is, the isophotal flux (photometry derived from the counts above
        the threshold minus the background) divided by the RMS error for the
        isophotal flux; in other words, the signal divided by the noise.

        """

        with open(path, 'rt') as fd:
            contents = [line.split() for line in fd]

        get_index = functools.partial(cls._find_column, contents)
        x_index = get_index('X_IMAGE')
        y_index = get_index('Y_IMAGE')
        alpha_index = get_index('ALPHA_SKY')
        delta_index = get_index('DELTA_SKY')
        area_index = get_index('ISOAREAF_IMAGE')
        mag_index = get_index('MAG_AUTO')
        flux_index = get_index('FLUX_ISO')
        fluxerr_index = get_index('FLUXERR_ISO')
        flux_radius_index = get_index('FLUX_RADIUS')
        flags_index = get_index('FLAGS')
        elong_index = get_index('ELONGATION')

        for line in contents:
            if line[0][0] != '#': # ignore comments

                def get_param(index, type_ = float):
                    """ Get the index-th element of 'line', cast to 'type_'"""
                    return type_(line[index])

                x = get_param(x_index)
                y = get_param(y_index)
                alpha = get_param(alpha_index)
                delta = get_param(delta_index)
                area = get_param(area_index, type_ = int)
                mag = get_param(mag_index)
                flux = get_param(flux_index)
                fluxerr = get_param(fluxerr_index)
                flux_radius = get_param(flux_radius_index)
                flags = get_param(flags_index, type_ = int)
                elongation = get_param(elong_index)

                saturated = Catalog.flag_saturated(flags)
                snr = flux / fluxerr
                fwhm = flux_radius * 2

                args = (x, y, alpha, delta, area, mag, saturated, snr,
                        fwhm, elongation)

                yield Star(*args)

    def __new__(cls, path):
        stars = cls._load_stars(path)
        catalog = super(Catalog, cls).__new__(cls, stars)
        catalog._path = path
        return catalog

    @property
    def path(self):
        """ Read-only 'path' attribute """
        return self._path

    @classmethod
    def from_sequence(cls, *stars):
        """ Create a Catalog from a sequence of Stars.

        Return a Catalog that is not the result of loading a SExtractor catalog
        into memory, but that encapsulates a series of Star objects. Note that,
        being an 'in-memory' catalog, so to speak, the returned instance does
        not have the 'path' attribute, so any attempt to access it will raise
        the AttributeError exception.

        """

        return super(Catalog, cls).__new__(cls, stars)

    def get_image_coordinates(self):
        """ Return as a list the image coordinates of the stars in the catalog.

        The method loops over the stars present in the SExtractor catalog,
        encapsulating the image coordinates of each star in a Pixel instance,
        all of which are returned in a list.

        """
        return [star.img_coords for star in self]

def sextractor_md5sum(options = None):
    """ Return the MD5 hash of the SExtractor configuration.

    This method returns the MD5 hash of the concatenation of the four
    configuration files (.sex, .param, .conv and .nnw) used by SExtractor, as
    well as the command-line configuration parameters (given in 'options', a
    dictionary that maps each parameter to its value, both strings) that
    override the corresponding definition in the configuration files or any
    default value. The hash is returned expressed as a 32-digit hexadecimal
    number.

    Note that the returned MD5 hash is just that of the input SExtractor
    configuration files and the overriding command-line options, not those
    values that were used in the end by SExtractor. This means that, for
    example, a SATUR_LEVEL of 50000 in the configuration file overridden by a
    SATUR_LEVEL option with a value of 45000 returns a different hash than a
    SATUR_LEVEL of 45000 in the configuration file and no command-line option,
    although in practical terms they are the same configuration.

    Although we could use an even more secure hash function, that would be
    overkill. The possibility of a colision is already rather small: a MD5 hash
    is 128 bits long, so assuming all hashes have an equal chance of occuring,
    the odds of any two random strings hashing to the same value are 1 in 2^128
    [http://ask.metafilter.com/50343/MD5-and-the-probability-of-collisions]

    The IOError exception is raised if any of the four SExtractor configuration
    files does not exist or is not readable. TypeError is raised if 'options'
    is not a dictionary or any of its keys or values is not a string. The
    latter means that, to compute the hash overriding the saturation level
    specified in the configuration file, something like {'SATUR_LEVEL' :
    '45000'}, for example, must be used.

    """

    sex_files = (SEXTRACTOR_CONFIG, SEXTRACTOR_PARAMS,
                 SEXTRACTOR_FILTER, SEXTRACTOR_STARNNW)

    md5 = hashlib.md5()
    for path in sex_files:
        with open(path, 'rt') as fd:
            for line in fd:
                md5.update(line)

    if options:
        # CPython returns the elements of a dictionary in an arbitrary order,
        # so it is necessary to sort the items to guarantee that two different
        # dictionaries with the same (key, value) pairs return the same hash.
        try:
            for key, value in sorted(options.items()):
                md5.update(key)
                md5.update(value)
        except AttributeError:
            msg = "'options' must be a dictionary"
            raise TypeError(msg)

    return md5.hexdigest()

def sextractor_version():
    """ Return the SExtractor version as a tuple.

    Run SExtractor with the --version option as its sole argument, capture the
    standard output and parse it. The version number of SExtractor is returned
    as a tuple (major, minor, micro), such as (2, 8, 6). SExtractorNotInstalled
    is raised if its executable cannot be found in the current environment.

    """

    # For example: "SExtractor version 2.8.6 (2009-04-09)"
    PATTERN = "^SExtractor version (\d\.\d\.\d) \(\d{4}-\d{2}-\d{2}\)$"

    for executable in SEXTRACTOR_COMMANDS:
        if methods.which(executable):
            break
    else:
        msg = "SExtractor not found in the current environment"
        raise SExtractorNotInstalled(msg)

    try:
        with tempfile.TemporaryFile() as fd:
            args = [executable, '--version']
            subprocess.check_call(args, stdout = fd)
            fd.seek(0)
            output = fd.readline()
        version = re.match(PATTERN, output).group(1)
        # From, for example, '2.8.6' to (2, 8, 6)
        return tuple(int(x) for x in version.split('.'))

    except subprocess.CalledProcessError, e:
        raise SExtractorError(e.returncode, e.cmd)

def sextractor(path, ext = 0, options = None, stdout = None, stderr = None):
    """ Run SExtractor on the image and return the path to the output catalog.

    This function runs SExtractor on 'path', using the configuration files
    defined in the module-level variables SEXTRACTOR_CONFIG, SEXTRACTOR_PARAMS,
    SEXTRACTOR_FILTER and SEXTRACTOR_STARNNW. It returns the path to the output
    catalog, which is saved to a temporary location and for whose deletion when
    it is no longer needed the user is responsible.

    The SExtractorNotInstalled exception is raised if a SExtractor executable
    cannot be found, and IOError if any of the four SExtractor configuration
    files does not exist or is not readable. If a SExtractor version earlier
    than SEXTRACTOR_REQUIRED_VERSION is installed, SExtractorUpgradeRequired
    is raised; this is necessary because the syntax that allows us to select on
    which extension sources are detected was not added until version 2.8.6. Any
    errors thrown by SExtractor are propagated as SExtractorError exceptions.
    Lastly, TypeEror is raised if (a) 'ext' is not an integer or (b) 'options'
    is not a dictionary or any of its keys or values is not a string.

    Keyword arguments:
    ext - for multi-extension FITS images, the index of the extension on which
          SExtractor will be run. It defaults to zero, meaning that sources are
          detected on the first extension of the FITS image. If a nonexistent
          extension is specified, the execution of SExtractor fails and the
          SExtractorError exception is raised.
    options - a dictionary mapping each SExtractor parameter to its value, and
              that will override their definition in the configuration files or
              any default value. In this manner, it is possible to execute
              SExtractor with different parameters without having to modify the
              configuration files. For example, {'CLEAN' : 'N', 'CLEAN_PARAM' :
              '1.1'}, would make SExtractor run with the parameters 'CLEAN' set
              to 'N' and 'CLEAN_PARAM' set to 1.1, regardless of what the
              configuration files say. All the keys and values in this
              dictionary must be strings.
    stdout - standard output file handle. If None, no redirection will occur.
    stderr - standard error file handle. If None, no redirection will occur.

    """

    # It is easier to ask forgiveness than permission, yes, but checking the
    # type here helps avoid some subtle errors. If, say, 'ext' is assigned a
    # value of 3.8, we do not want it to be silently casted (and truncated)
    # to three; it is much better (and safer) to have TypeError raised and
    # let the user know that an invalid, non-integer index was given.

    if not isinstance(ext, (int, long)):
        raise TypeError("'ext' must be an integer")

    for executable in SEXTRACTOR_COMMANDS:
        if methods.which(executable):
            break
    else:
        msg = "SExtractor not found in the current environment"
        raise SExtractorNotInstalled(msg)

    if sextractor_version() < SEXTRACTOR_REQUIRED_VERSION:
        # From, for example, (2, 8, 6) to '2.8.6'
        version_str = '.'.join(str(x) for x in SEXTRACTOR_REQUIRED_VERSION)
        msg = "SExtractor version %s or newer is needed" % version_str
        raise SExtractorUpgradeRequired(msg)

    # If the loop did not break (and thus SExtractorNotInstalled was not
    # raised), 'executable' contains the first command that was found

    root, _ = os.path.splitext(os.path.basename(path))
    catalog_fd, catalog_path = \
        tempfile.mkstemp(prefix = '%s_' % root, suffix = '.cat')
    os.close(catalog_fd)

    # Raise IOError if any of the configuration files is nonexistent or not
    # readable. We cannot trust that SExtractor will fail when this happens as
    # it may not abort the execution, but instead just issue a warning and use
    # the internal defaults. As of version 2.8.6, only -PARAMETERS_NAME and
    # -FILTER_NAME, if unreadable, cause the execution of SExtractor to fail.

    for config_file in (SEXTRACTOR_CONFIG, SEXTRACTOR_PARAMS,
                        SEXTRACTOR_FILTER, SEXTRACTOR_STARNNW):

        if not os.path.exists(config_file):
            msg = "configuration file %s not found"
            raise IOError(msg % config_file)
        if not os.access(config_file, os.R_OK):
            msg = "configuration file %s cannot be read"
            raise IOError(msg % config_file)

    args = [executable, path + '[%d]' % ext,
            '-c', SEXTRACTOR_CONFIG,
            '-PARAMETERS_NAME', SEXTRACTOR_PARAMS,
            '-FILTER_NAME', SEXTRACTOR_FILTER,
            '-STARNNW_NAME', SEXTRACTOR_STARNNW,
            '-CATALOG_NAME', catalog_path]

    if options:
        try:
            for key, value in options.iteritems():
                args += ['-%s' % key, value]
        except AttributeError:
            msg = "'options' must be a dictionary"
            raise TypeError(msg)

    try:
        subprocess.check_call(args, stdout = stdout, stderr = stderr)
        return catalog_path
    except subprocess.CalledProcessError, e:
        try: os.unlink(catalog_path)
        except (IOError, OSError): pass
        raise SExtractorError(e.returncode, e.cmd)

