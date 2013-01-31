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
import math
import shutil
import tempfile
import subprocess

# LEMON modules
import fitsimage
import methods

ASTROMATIC_FILES = os.path.abspath('./astromatic')
SEXTRACTOR_CONFIG = os.path.join(ASTROMATIC_FILES, 'sextractor.sex')
SEXTRACTOR_PARAMS = os.path.join(ASTROMATIC_FILES, 'sextractor.param')
SEXTRACTOR_FILTER = os.path.join(ASTROMATIC_FILES, 'sextractor.conv')
SEXTRACTOR_STARNNW = os.path.join(ASTROMATIC_FILES, 'sextractor.nnw')
SEXTRACTOR_COMMANDS = 'sextractor', 'sex' # may be any of these

SCAMP_CONFIG = os.path.join(ASTROMATIC_FILES, 'scamp.conf')
SCAMP_HEADER_SUFFIX = '.head' # extension of header files
SCAMP_COMMAND = 'scamp'
ACLIENT_COMMAND = 'aclient'

SWARP_CONFIG = os.path.join(ASTROMATIC_FILES, 'swarp.conf')
SWARP_COMMAND = 'swarp'


class CDSclientNotInstalled(StandardError):
    """ Raised if F.Ochsenbein's CDSclient package is not installed """
    pass

class SExtractorNotInstalled(StandardError):
    pass

class SExtractorError(subprocess.CalledProcessError):
    pass

class SCAMPNotInstalled(StandardError):
    pass

class SCAMPError(subprocess.CalledProcessError):
    pass

class SWarpNotInstalled(StandardError):
    pass

class SWarpError(subprocess.CalledProcessError):
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
    """ An immutable class with a sourced detected by SExtractor. """

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

def sextractor(path, ext = 0, options = None, stdout = None, stderr = None):
    """ Run SExtractor on the image and return the path to the output catalog.

    This function runs SExtractor on 'path', using the configuration files
    defined in the module-level variables SEXTRACTOR_CONFIG, SEXTRACTOR_PARAMS,
    SEXTRACTOR_FILTER and SEXTRACTOR_STARNNW. It returns the path to the output
    catalog, which is saved to a temporary location and for whose deletion when
    it is no longer needed the user is responsible.

    The SExtractorNotInstalled exception is raised if a SExtractor executable
    cannot be found, and IOError if any of the four SExtractor configuration
    files does not exist or is not readable. Any errors thrown by SExtractor
    are propagated as SExtractorError exceptions. Lastly, TypeEror is raised if
    (a) 'ext' is not an integer or (b) 'options' is not a dictionary or any of
    its keys or values is not a string.

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

def ahead_file(img, output_path, scale, equinox, radecsys,
               ra_keyword = 'RA', dec_keyword = 'DEC'):
    """ Generate the .ahead file needed by SCAMP in order to do astrometry.

    This function receives a fitsimage.FITSImage object, which encapsulates
    an astronomical image, and creates the .ahead file that guarantees that
    the required FITS keywords needed by Emmanuel Bertin's SCAMP (defining
    an initial guess of the astrometic solution) are present. Although the
    keywords could also be directly added or updated in the FITS image,
    using the .ahead file allows us not to modify the file.

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
    line. External 'headers' need not contain all the FITS keywords
    normally required. The keywords present in external headers are only
    there to override their counterparts in the original image headers or
    to add new ones. Hence for every input (say, xxxx.cat) FITS catalogue,
    SCAMP looks for a xxxx.ahead header file, loads it if present, and
    overrides or adds to image header keywords those found there.'

    For further information on the World Coordinate System (WCS), refer to
    http://tdc-www.harvard.edu/software/wcstools/wcstools.wcs.html

    img - FITSImage object for which to generate the .ahead file.
    output_path - path to which the .ahead file will be saved.
    scale - scale of the image, in arcseconds per pixel
    equinox - equinox in years (e.g., 2000)
    radecsys - reference system (e.g., ICRS)

    Keyword arguments:
    ra_keyword - FITS keyword for the right ascension, in decimal degrees.
    dec_keyword - FITS keyword for the declination, in decimal degrees.

    """

    with open(output_path, 'wt') as fd:

        # Pixel coordinates of the reference point
        xcenter, ycenter = img.center
        fd.write("CRPIX1  = %d\n" % xcenter)
        fd.write("CRPIX2  = %d\n" % ycenter)

        ra  = img.read_keyword(ra_keyword)
        dec = img.read_keyword(dec_keyword)

        fd.write("CRVAL1  = %s\n" % ra)    # RA at reference point
        fd.write("CRVAL2  = %s\n" % dec)   # DEC at reference point
        fd.write("CTYPE1  = 'RA---TAN'\n") # Gnomonic projection
        fd.write("CTYPE2  = 'DEC--TAN'\n")

        # The FITS WCS standard uses a rotation matrix, CD1_1, CD1_2, CD2_1 and
        # CD2_2 to indicate both rotation and scale, allowing a more intuitive
        # computation if the axes are skewed. This model has been used by HST
        # and IRAF for several years.

        scale /= 3600 # arcsec/pixel to deg/pixel
        fd.write("CD1_1   = %f\n" % scale)
        fd.write("CD1_2   = 0.0\n")
        fd.write("CD2_1   = 0.0\n")
        fd.write("CD2_2   = %f\n" % -scale)
        fd.write("EQUINOX = %d\n" % equinox)    # equinox in years
        fd.write("RADECSYS= '%s'\n" % radecsys) # reference system

        fd.write("END⊔⊔⊔⊔⊔")

def scamp(path, scale, equinox, radecsys, saturation, ext = 0, options = None,
          ra_keyword = 'RA', dec_keyword = 'DEC', stdout = None, stderr = None):
    """ Return a FITS-like image header with updated astrometric information.

    This function runs SExtractor on the input FITS image, creating a catalog
    in the FITS_LDAC binary format which is then read by SCAMP to compute the
    astrometric projection parameters. This is achieved by automatically
    downloading a catalog of astrometric standards from the VizieR database and
    cross-matching them with the sources detected by SExtractor in the image.
    The result of this calibration is a FITS-like image header which contains
    updated astrometric information and that SWarp can read and use for image
    stacking. The function returns the path to the .head file, which is saved
    to a temporary file.

    The process of computing the astrometric solution is very dependent on the
    preliminary, approximate astrometric data of the FITS image. Therefore, it
    is extremely important that the FITS keywords with the right ascension and
    declination have the right values, so that the astrometric standards
    downloaded from VizieR correspond to stars in the image.

    Arguments:
    scale - scale of the image, in degrees per pixel.
    equinox - equinox in years (e.g., 2000)
    radecsys - celestial reference system (e.g., ICRS)
    saturation - the number of ADUs at which arises saturation.

    Keyword arguments:
    ext - for multi-extension FITS images, the index of the extension on which
          SExtractor will be run and that therefore will determine the updated
          astrometric information saved to the FITS-like image header. It
          defaults to zero, which means that the first extension of the FITS
          image is used. If a nonexistent extension is specified, the execution
          of SExtractor fails and the SExtractorError exception is raised.
    options - a dictionary mapping each SCAMP parameter to its value, and that
              will override their definition in the configuration file or any
              default value. In this manner, it is possible to execute SCAMP
              with different parameters without having to modify the
              configuration file. For example, {'POSANGLE_MAXERR' : '2.5'}
              would make SCAMP run with the parameter 'POSANGLE_MAXERR' set to
              2.5, regardless of what the configuration file says. All the keys
              and values in this dictionary must be strings.
    ra_keyword - FITS keyword for the right ascension, in decimal degrees.
    dec_keyword - FITS keyword for the declination, in decimal degrees.
    stdout - standard output file handle. If None, no redirection will occur.
    stderr - standard error file handle. If None, no redirection will occur.

    """

    emsg = "'%s' not found in the current environment"
    if not methods.which(ACLIENT_COMMAND):
        raise CDSclientNotInstalled(emsg % ACLIENT_COMMAND)

    if not methods.which(SCAMP_COMMAND):
        raise SCAMPNotInstalled(emsg % SCAMP_COMMAND)

    if not os.path.exists(SCAMP_CONFIG):
        msg = "configuration file %s not found"
        raise IOError(msg % SCAMP_CONFIG)
    if not os.access(SCAMP_CONFIG, os.R_OK):
        msg = "configuration file %s cannot be read"
        raise IOError(msg % SCAMP_CONFIG)

    # Do all the work in a temporary directory and remove it
    # (and thus all the involved files) upon the method exit.
    tmp_dir = tempfile.mkdtemp(suffix = '_scamp')

    try:
        img = fitsimage.FITSImage(path)

        # Determine the paths to the FITS LDAC SExtractor catalog, the
        # '.ahead' header file and the merged catalogue produced by SCAMP.
        tmp_path = os.path.join(tmp_dir, img.basename_woe)
        ldac_path = '%s.ldac' % tmp_path
        ahead_path = '%s.ahead' % tmp_path
        merged_path = '%s.cat' % tmp_path

        # Use the FITS LDAC' format and image saturation level
        sextractor_options = dict(CATALOG_TYPE = 'FITS_LDAC',
                                  SATUR_LEVEL = str(saturation))

        # The SExtractor catalog is saved to a temporary file, which
        # we have to move to our temporary, working directory.
        kwargs = dict(ext = ext, options = sextractor_options,
                      stdout = stdout, stderr = stderr)
        output_path = sextractor(img.path, **kwargs)
        shutil.move(output_path, ldac_path)

        # Then, create the SCAMP .ahead file
        kwargs = dict(ra_keyword = ra_keyword, dec_keyword = dec_keyword)
        ahead_file(img, ahead_path, scale, equinox, radecsys, **kwargs)

        # Finally, run SCAMP on the image. Those keywords defined in the .ahead
        # file will override their counterparts in the original FITS header.
        # The output catalog will be saved in ASCII format, preceded with
        # comment lines listing column labels.
        args = [SCAMP_COMMAND, ldac_path,
                '-c', SCAMP_CONFIG,
                '-MERGEDOUTCAT_NAME', merged_path,
                '-MERGEDOUTCAT_TYPE', 'ASCII_HEAD',
                '-HEADER_SUFFIX', SCAMP_HEADER_SUFFIX]

        if options:
            try:
                for key, value in options.iteritems():
                    args += ['-%s' % key, value]
            except AttributeError:
                msg = "'options' must be a dictionary"
                raise TypeError(msg)

        try:
            subprocess.check_call(args, stdout = stdout, stderr = stderr)
        except subprocess.CalledProcessError, e:
            raise SCAMPError(e.returncode, e.cmd)

        # Move the output 'head' file to a different, temporary location
        # (as this directory is going to be deleted when we are done)
        src_head_path = tmp_path + SCAMP_HEADER_SUFFIX
        dst_head_fd, dst_head_path = tempfile.mkstemp(suffix = SCAMP_HEADER_SUFFIX)
        os.close(dst_head_fd)
        shutil.move(src_head_path, dst_head_path)
        return dst_head_path

    finally:
        # Needed if it could not be moved to the temporary directory
        try: os.unlink(output_path)
        except (NameError, IOError, OSError): pass

        try: shutil.rmtree(tmp_dir)
        except (IOError, OSError): pass

def swarp(img_path, head_path, copy_keywords = None,
          stdout = None, stderr = None):
    """ Use SWarp to merge a FITS-like image header with an image.

    This method uses the SWarp stacking tool to merge a FITS image with a .head
    header, thus updating it with the astrometric information that was computed
    and saved by SCAMP to this external header file. Returns the path to the
    output FITS image, which is saved to a temporary file and for whose
    deletion when it is no longer needed the user is responsible.

    Keyword arguments:
    copy_keywords - FITS keywords, and their values, to propagate from the
                    input image header to the resampled and coadded image
                    header produced by SWarp.  Needed since not all FITS
                    keywords are automatically copied to the output image
                    header, as many of them become irrelevant.
    stdout - the SWarp standard output file handle. If set to None, no
             redirection will occur.
    stderr - the SWarp standard error file handle. If set to None, no
             redirection will occur.

    """

    # The work will be done in a temporary directory, so we need to save the
    # absolute path to both the current directory (to be able to move back to
    # it once we are finished) and to the input FITS image (as if a relative
    # path were given its meaning would not be the same when we changed our
    # directory).

    cwd = os.getcwdu() # save current directory
    img_path = os.path.abspath(img_path)
    tmp_dir = tempfile.mkdtemp(suffix = '_swarp')

    try:
        # Symlink the image and the header file from the temporary
        # directory. We cannot use hard links as they would not work
        # if the temporary directory is created in a different device
        # ('Invalid cross-device link' error)
        img_basename = os.path.basename(img_path)
        image_dest = os.path.join(tmp_dir, img_basename)
        root, ext = os.path.splitext(img_basename)
        header_dst = os.path.join(tmp_dir, '%s%s' % (root, SCAMP_HEADER_SUFFIX))
        os.symlink(img_path, image_dest)
        os.symlink(head_path, header_dst)

        output_fd, output_path = \
            tempfile.mkstemp(prefix = '%s_' % root, suffix = ext)
        os.close(output_fd)

        os.chdir(tmp_dir)
        args = [SWARP_COMMAND, img_basename,
                '-c', SWARP_CONFIG,
                '-IMAGEOUT_NAME', output_path,
                '-HEADER_SUFFIX', SCAMP_HEADER_SUFFIX]

        if copy_keywords:
            # Coma-separated list of FITS keywords that will be propagated from
            # the input FITS header to the coadded and resampled image header.
            args.append('-COPY_KEYWORDS')
            args.append(','.join(copy_keywords))

        swarp_command = args[0]
        if not methods.which(swarp_command):
            msg = "'%s' not found in the current environment" % swarp_command
            raise SWarpNotInstalled(msg)

        try:
            subprocess.check_call(args, stdout = stdout, stderr = stderr)
        except subprocess.CalledProcessError, e:
            raise SWarpError(e.returncode, e.cmd)

        # The latest stable version (v2.19.1) of SWarp does not propagate FITS
        # keywords with whitespaces (such as "LEMON FWHM"). This is not much of
        # a problem as we can check that all the keywords have been actually
        # propagated, manually copying them to the header of the output image
        # if that is not the case. It is not that we do not trust SWarp; just
        # that we prefer to make sure everything went well.

        input_img = fitsimage.FITSImage(img_path)
        output_img = fitsimage.FITSImage(output_path)

        for keyword in copy_keywords:

            # Read the keyword from the output image; we do not check its
            # value, but just make sure that it is there, as expected. If
            # it is missing, the KeyError exception is raised and the
            # keyword is copied from the header of the input image.

            try:
                output_img.read_keyword(keyword)

            except KeyError:

                try:
                    value = input_img.read_keyword(keyword)
                    kwargs = dict(comment = "Propagated keyword")
                    output_img.update_keyword(keyword, value, **kwargs)

                # Ignore the keyword (it's not in the input image either)
                except KeyError:
                    pass

        # Do not return a FITSImage object, only the path
        return output_img.path

    finally:
        os.chdir(cwd) # move back to original directory
        try: shutil.rmtree(tmp_dir)
        except (IOError, OSError): pass

