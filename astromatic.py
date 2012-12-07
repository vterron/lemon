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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division

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
SCAMP_CONFIG = os.path.join(ASTROMATIC_FILES, 'scamp.conf')
SWARP_CONFIG = os.path.join(ASTROMATIC_FILES, 'swarp.conf')

# Filename extension of header files returned by SCAMP
HEADER_SUFFIX = '.head'

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

class Pixel(object):
    """ A pair of x- and y-coordinates. """

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return '%s(%f, %f)' % (self.__class__.__name__, self.x, self.y)

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __hash__(self):
        return hash((self.x, self.y))

    def distance(self, another):
        """ Return the Euclidean distance between two Pixels """

        x_axis = pow(self.x - another.x, 2)
        y_axis = pow(self.y - another.y, 2)
        return math.sqrt(x_axis + y_axis)

class Star(object):
    """ A Star is primarily made up of two Pixels: one that encapsulates the
    image coordinates of the Star and another that encapsulates its celestial
    coordinates (right ascension and declination).

    """

    def __init__(self, x, y, alpha, delta, isoareaf, magnitude, saturated,
                 snr, fwhm, elongation):
        """ Instantiate a Star with values read from a SExtractor catalog.

        x - star position along x.
        y - star position along y.
        alpha - right ascension of the star.
        delta - declination of the star.
        isoareaf - isophotal area (filtered) above detection threshold (pix^2).
        magnitude - the measure of the brightness of the Star as seen by an
                    observer on Earth.
        saturated - indicated whether at least one pixel of the Star is
                    saturated, or very close to.
        snr - the signal-to-noise ratio of the star. This is the only value not
              directly read from the SExtractor catalog, but which has to be
              derived by us from other values.
        fwhm - the full width at half maximum (FWHM) of the star
        elongation - the value of A/B, where A and B are its semi-major and
                     semi-minor axis lengths, as reported by SExtractor. More
                     precisely, A and B represent the maximum and minimum
                     spatial rms of the object profile along any direction.

        """

        self.img_coords = Pixel(x, y)
        self.sky_coords = Pixel(alpha, delta)
        self.area = isoareaf
        self.mag = magnitude
        self.saturated = saturated
        self.snr = snr
        self.fwhm = fwhm
        self.elongation = elongation

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

    def __str__(self):
        return '%s<(%f, %f)-->(%f, %f)' % (self.__class__.__name__,
               self.x, self.y, self.alpha, self.delta)

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

class Catalog(list):
    """ High-level interface to a SExtractor catalog """

    def _read_catalog(self):
        """ A generator over the contents of the catalog.

        The method returns, one-by-one, the lines in the SExtractor catalog,
        splitting the words in a list where runs of consecutive whitespace are
        regarded as a single separator. ['#', '1', 'X_IMAGE', 'Object',
        'position', 'along', 'x', '[pixel]'] is an example of what may
        be returned.

        """
        with open(self.path, 'rt') as fd:
            for line in fd:
                yield line.split()

    def _find_column(self, parameter):
        """ Return the index of a SExtractor paramater in the catalog.

        The method returns the zero-based index of the column in which a
        SExtractor parameter (such as X_IMAGE or FLUX_MAX) was written to a
        catalog. As happens with list.index, ValueError is raised if the
        parameter cannot be found in the catalog. For user's convenience the
        search is case-insensitive, so 'x_world' would match 'X_WORLD'.

        Is is important no note that, for the parameter being able to be found,
        the catalog must have been saved in the SExtractor ASCII_HEAD format,
        so that the file stars with comment lines listing column labels.

        In order to avoid unnecessary disk accesses, already-found parameter
        indexes internally cached, so that they can be directly returned the
        next time they are needed. I would have rather used our memoization
        decorator, but it is not that simple as, subclassing from list, this
        class is unhashable.

        """

        attr = '__%s_index' % parameter

        try:
            return getattr(self, attr)

        except AttributeError:

            # We need to examine the comment lines (those that start with a
            # '#') and extract the value of the column associated to the
            # parameter.  SExtractpr comments follow the following format:
            #
            #       # 1 X_IMAGE     Object position along x     [pixel]
            #       # 2 Y_IMAGE     Object position along y     [pixel]
            #
            # Where the first integer in each line, right after the '#',
            # indicates the column of the parameter. Note that we must
            # subtract one from these indexes, as they are one-based.

            comments = (line
                        for line in self._read_catalog()
                        if line and line[0] and line[0][0] == '#')

            for line in comments:
                param_name = line[2]
                if param_name.upper() == parameter.upper():
                    param_index = int(line[1]) - 1
                    setattr(self, attr, param_index)
                    return self._find_column(parameter)
            else:
                raise ValueError("'%s' not in %s" % (parameter, self.path))

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

        """

        binary_repr = bin(flag_value)
        return len(binary_repr) >= 3 and binary_repr[-3] == '1'

    def _load_stars(self):
        """ Load a SExtractor catalog into memory.

        The method parses a SExtractor catalog and adds to 'self', after
        emptying it, an instance of the Star class for each detected object.
        It is mandatory, or the ValueError exception will be raised otherwise,
        that the following parameters are present in the catalog: X_IMAGE,
        Y_IMAGE, ALPHA_SKY, DELTA_SKY, ISOAREAF_IMAGE, MAG_AUTO, FLUX_ISO,
        FLUXERR_ISO, FLUX_RADIUS, ELONGATION and CLASS_STAR. Also, the catalog
        must have been saved in the SExtractor ASCII_HEAD format, as the
        comment lines listing column labels are needed in order to detect in
        which column each parameter is.

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

        del self[:]

        for line in self._read_catalog():
            if line[0][0] != '#': # ignore comments
                x           = float(line[self._find_column('X_IMAGE')])
                y           = float(line[self._find_column('Y_IMAGE')])
                alpha       = float(line[self._find_column('ALPHA_SKY')])
                delta       = float(line[self._find_column('DELTA_SKY')])
                isoareaf    =   int(line[self._find_column('ISOAREAF_IMAGE')])
                mag_auto    = float(line[self._find_column('MAG_AUTO')])
                flux_iso    = float(line[self._find_column('FLUX_ISO')])
                fluxerr_iso = float(line[self._find_column('FLUXERR_ISO')])
                flux_radius = float(line[self._find_column('FLUX_RADIUS')])

                flags = int(line[self._find_column('FLAGS')])
                saturated = Catalog.flag_saturated(flags)
                snr = flux_iso / fluxerr_iso
                fwhm = flux_radius * 2

                elongation     = float(line[self._find_column('ELONGATION')])
                self.append(Star(x, y, alpha, delta, isoareaf, mag_auto,
                                 saturated, snr, fwhm, elongation))

        return self

    def __init__(self, path):
        super(list, self).__init__()
        self.path = path
        self._load_stars()

    def get_image_coordinates(self):
        """ Return as a list the image coordinates of the stars in the catalog.

        The method loops over the stars present in the SExtractor catalog,
        encapsulating the image coordinates of each star in a Pixel instance,
        all of which are returned in a list.

        """
        return [star.img_coords for star in self]

def sextractor_md5sum(options):
    """ Return the MD5 hash of the SExtractor configuration.

    This method returns the MD5 hash of the concatenation of the four
    configuration files (.sex, .param, .conv and .nnw) used by SExtractor, as
    well as the command-line configuration parameters (given in 'options' as
    a sequence of strings) that override the corresponding definition in the
    configuration files or any default value.  The hash is returned expressed
    as a 32-digit hexadecimal number.

    Note that the returned MD5 hash is just that of the input SExtractor
    configuration files and the overriding command-line options, not those
    values that were used in the end by SExtractor. This means that, for
    example, a SATUR_LEVEL of 50000 in the configuration file overridden by a
    -SATUR_LEVEL option with a value of 45000 returns a different hash than a
    SATUR_LEVEL of 45000 in the configuration file and no command-line option,
    although in practical terms they are the same configuration.

    Although we could use an even more secure hash function, that would be
    overkill. The possibility of a colision is already rather small: a MD5 hash
    is 128 bits long, so assuming all hashes have an equal chance of occuring,
    the odds of any two random strings hashing to the same value are 1 in 2^128
    [http://ask.metafilter.com/50343/MD5-and-the-probability-of-collisions]

    """

    sex_files = (SEXTRACTOR_CONFIG, SEXTRACTOR_PARAMS,
                 SEXTRACTOR_FILTER, SEXTRACTOR_STARNNW)

    md5 = hashlib.md5()
    for path in sex_files:
        with open(path, 'rt') as fd:
            for line in fd:
                md5.update(line)

    for opt in options:
        md5.update(opt)

    return md5.hexdigest()

def sextractor(path, options = None, stdout = None, stderr = None):
    """ Run SExtractor on the image and return the path to the output catalog.

    The method rus SExtractor (which is assumed to be installed on the system)
    on the input image, using the .sex, .param, .conv and .nnw configuration
    files defined by the global variables SEXTRACTOR_CONFIG, SEXTRACTOR_PARAMS,
    SEXTRACTOR_FILTER and SEXTRACTOR_STARNNW. It returns the path to its output
    catalog, saved to a temporary file, for whose deletion when it is no longer
    needed the user is responsible.

    Keyword arguments:
    options - sequence of parameters and their corresponding values, which
              override their definition in the configuration files or any
              default value. In this manner, it is possible to execute
              SExtractor with different values without having to modify the
              configuration files. For example, ['-CLEAN', 'N', '-CLEAN_PARAM',
              1.1], would make SExtractor run with the paramters 'CLEAN' set to
              'N' and 'CLEAN_PARAM' set to 1.1, regardless of what the
              configuration files say.
    stdout - the SExtractor standard output file handle. If set to None, no
             redirection will occur.
    stderr - the SExtractor standard error file handle. If set to None, no
             redirection will occur.

    """

    for executable in ['sextractor', 'sex']:
        if methods.check_command(executable):
            break
    else:
        msg = "SExtractor not found in the current environment"
        raise SExtractorNotInstalled(msg)

    # If the loop did not break (and thus SExtractorNotInstalled was not
    # raised), 'executable' contains the first command that was found

    root, ext = os.path.splitext(os.path.basename(path))
    catalog_fd, catalog_path = \
        tempfile.mkstemp(prefix = '%s_' % root, suffix = '.cat')
    os.close(catalog_fd)

    args = [executable, path,
            '-c', SEXTRACTOR_CONFIG,
            '-PARAMETERS_NAME', SEXTRACTOR_PARAMS,
            '-FILTER_NAME', SEXTRACTOR_FILTER,
            '-STARNNW_NAME', SEXTRACTOR_STARNNW,
            '-CATALOG_NAME', catalog_path]

    if options:
        args += options

    try:
        subprocess.check_call(args, stdout = stdout, stderr = stderr)
        return catalog_path
    except subprocess.CalledProcessError, e:
        try: os.unlink(catalog_path)
        except (IOError, OSError): pass
        raise SExtractorError(e.returncode, e.cmd)

def scamp(path, scale, equinox, radecsys, saturation, ra_keyword = 'RA',
          dec_keyword = 'DEC', stdout = None, stderr = None):
    """ Run SCAMP to create a FITS-like image header that SWarp can read.

    The method runs SExtractor on the image, outputting a catalog in the
    'FITS LDAC' binary format, and then uses SCAMP to read it and produce
    the FITS-like image header (the '.head' header file), containing updated
    astrometric and photometric information, than SWarp can read and use for
    image stacking. Returns the path to the .head file, which is saved to a
    temporary file and for whose deletion when it is no longer needed the
    user is responsible.

    scale - scale of the image, in degrees per pixel
    equinox - equinox in years (e.g., 2000)
    radecsys - reference system (e.g., ICRS)
    saturation - the number of ADUs at which arises saturation.

    Keyword arguments:
    ra_keyword - FITS keyword for the right ascension, in decimal degrees.
    dec_keyword - FITS keyword for the declination, in decimal degrees.
    stdout - the SExtractor and SCAMP standard output file handle. If set
             to None, no redirection will occur.
    stderr - the SExtractor and SCAMP standard error file handle. If set
             to None, no redirection will occur.

    """

    aclient_command = 'aclient'
    scamp_command = 'scamp'

    emsg = "'%s' not found in the current environment"
    if not methods.check_command(aclient_command):
        raise CDSclientNotInstalled(emsg % aclient_command)

    if not methods.check_command(scamp_command):
        raise CDSclientNotInstalled(emsg % scamp_command)

    # Do all the work in a temporary directory and remove it
    # (and thus all the involved files) upon the method exit.
    tmp_dir = tempfile.mkdtemp(suffix = '_scamp')

    try:
        img = fitsimage.FITSImage(path)

        # Determine the paths to the FITS LDAC SExtractor catalog, the
        # '.ahead' header file and the marged catalogue produced by SCAMP.
        tmp_path = os.path.join(tmp_dir, img.basename_woe)
        ldac_path = '%s.ldac' % tmp_path
        ahead_path = '%s.ahead' % tmp_path
        marged_path = '%s.cat' % tmp_path

        # Use the FITS LDAC' format and image saturation level
        options = ['-CATALOG_TYPE', 'FITS_LDAC',
                   '-SATUR_LEVEL', '%d' % saturation]

        # The SExtractor catalog is saved to a temporary file, which
        # we have to move to our temporary, working directory.
        output_path = sextractor(img.path, options = options,
                                 stdout = stdout, stderr = stderr)
        shutil.move(output_path, ldac_path)

        # Then, create the SCAMP .ahead file. Note the conversion
        # from arcsec/pixel to degrees/pix, as expected by SCAMP.
        img.ahead_file(ahead_path, scale / 3600, equinox, radecsys,
                       ra_keyword = ra_keyword, dec_keyword = dec_keyword)

        # Finally, run SCAMP on the image. Those keywords defined in the .ahead
        # file will override their counterparts in the original FITS header.
        # The output catalog will be saved in ASCII format, preceded with
        # comment lines listing column labels.
        args = ['scamp', ldac_path,
                '-c', SCAMP_CONFIG,
                '-MERGEDOUTCAT_NAME', marged_path,
                '-MERGEDOUTCAT_TYPE', 'ASCII_HEAD',
                '-HEADER_SUFFIX', HEADER_SUFFIX]

        scamp_command = args[0]
        if not methods.check_command(scamp_command):
            raise SCAMPNotInstalled(emsg % scamp_command)

        try:
            subprocess.check_call(args, stdout = stdout, stderr = stderr)
        except subprocess.CalledProcessError, e:
            raise SCAMPError(e.returncode, e.cmd)

        # Move the output 'head' file to a different, temporary location
        # (as this directory is going to be deleted when we are done)
        src_head_path = tmp_path + HEADER_SUFFIX
        dst_head_fd, dst_head_path = tempfile.mkstemp(suffix = HEADER_SUFFIX)
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
        header_dst = os.path.join(tmp_dir, '%s%s' % (root, HEADER_SUFFIX))
        os.symlink(img_path, image_dest)
        os.symlink(head_path, header_dst)

        output_fd, output_path = \
            tempfile.mkstemp(prefix = '%s_' % root, suffix = ext)
        os.close(output_fd)

        os.chdir(tmp_dir)
        args = ['swarp', img_basename,
                '-c', SWARP_CONFIG,
                '-IMAGEOUT_NAME', output_path,
                '-HEADER_SUFFIX', HEADER_SUFFIX]

        if copy_keywords:
            # Coma-separated list of FITS keywords that will be propagated from
            # the input FITS header to the coadded and resampled image header.
            args.append('-COPY_KEYWORDS')
            args.append(','.join(copy_keywords))

        swarp_command = args[0]
        if not methods.check_command(swarp_command):
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

