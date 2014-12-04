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

import copy
import functools
import mock
import numpy
import os
import os.path
import pyfits
import random
import shutil
import stat
import subprocess
import tempfile

from test import unittest
import astromatic
from astromatic import Pixel, Coordinates, Star, Catalog
import dss_images
import fitsimage
import methods

NITERS = 100

class PixelTest(unittest.TestCase):

    X_COORD_RANGE = (1, 2048)
    Y_COORD_RANGE = (1, 2048)

    @classmethod
    def random(cls):
        """ Return a random Pixel object """
        x = random.uniform(*cls.X_COORD_RANGE)
        y = random.uniform(*cls.Y_COORD_RANGE)
        return Pixel(x, y)

    def test_immutability(self):
        """ Make sure that the Pixel class is immutable """

        pixel = self.random()
        for name in pixel._asdict().iterkeys():
            value = random.random()
            self.assertRaises(AttributeError, setattr, pixel, name, value)
            self.assertRaises(AttributeError, delattr, pixel, name)

class CoordinatesTest(unittest.TestCase):

    RIGHT_ASCENSION_RANGE = (0, 360)
    DECLINATION_RANGE = (-90, 90)
    PM_RA_RANGE =  (0, 0.79858)
    PM_DEC_RANGE = (0, 10.32812)

    @staticmethod
    def get_random_pm(pm_range):
        """ Return a positive or negative random number in 'pm_range' """

        value = random.uniform(*pm_range)
        sign = random.choice([-1, 1])
        return sign * value

    @classmethod
    def random(cls):
        """ Return a random Coordinates object """

        ra = random.uniform(*cls.RIGHT_ASCENSION_RANGE)
        dec = random.uniform(*cls.DECLINATION_RANGE)

        # Use None for both 'pm_ra' and 'pm_dec' or none of them. If does not
        # make much sense to know the proper motion in declination but not in
        # right ascension, or vice versa.
        if random.choice([True, False]):
            pm_ra  = cls.get_random_pm(cls.PM_RA_RANGE)
            pm_dec = cls.get_random_pm(cls.PM_DEC_RANGE)
        else:
            pm_ra = pm_dec = None

        return Coordinates(ra, dec, pm_ra, pm_dec)

    def test_immutability(self):
        """ Make sure that the Coordinates class is immutable """

        coords = self.random()
        for name in coords._asdict().iterkeys():
            value = random.random()
            self.assertRaises(AttributeError, setattr, coords, name, value)
            self.assertRaises(AttributeError, delattr, coords, name)

    def test_distance(self):

        # http://www.astronomycafe.net/qadir/q1890.html (Sten Odenwald)
        coords1 = Coordinates(100.2, -16.58)
        coords2 = Coordinates(87.5, 7.38)
        distance = coords1.distance(coords2)
        self.assertAlmostEqual(distance, 27.054384870767787)

        # http://www.skythisweek.info/angsep.pdf (David Oesper)
        coords3 = Coordinates(165.458, 56.3825)
        coords4 = Coordinates(165.933, 61.7511)
        distance = coords3.distance(coords4)
        self.assertAlmostEqual(distance, 5.374111607543190)

    def test_get_exact_coordinates(self):

        # Barnard's Star (J2000): -798.58 10328.12 (mas/yr)
        barnard = Coordinates(269.452075, 4.693391, -0.79858, 10.32812)

        # year = epoch, so coordinates do not change
        coords = barnard.get_exact_coordinates(2000)
        self.assertEqual(coords.ra,  barnard.ra)
        self.assertEqual(coords.dec, barnard.dec)
        self.assertIs(coords.pm_ra,  None)
        self.assertIs(coords.pm_dec, None)

        # 2005 - 2000 = 5 years
        # 269.452075 + (-0.79858 * 5) / 3600 = 269.450965861
        #   4.693391 + (10.32812 * 5) / 3600 =   4.707735611
        coords = barnard.get_exact_coordinates(2005)
        self.assertAlmostEqual(coords.ra, 269.450965861)
        self.assertAlmostEqual(coords.dec,  4.707735611)
        self.assertIs(coords.pm_ra,  None)
        self.assertIs(coords.pm_dec, None)

        # 2014.5 = July 3, 2014
        # 2014.5 - 2000 = 14.5 years
        # 269.452075 + (-0.79858 * 14.5) / 3600 = 269.448858497
        #   4.693391 + (10.32812 * 14.5) / 3600 =   4.734990372
        coords = barnard.get_exact_coordinates(2014.5)
        self.assertAlmostEqual(coords.ra, 269.448858497)
        self.assertAlmostEqual(coords.dec,  4.734990372)
        self.assertIs(coords.pm_ra,  None)
        self.assertIs(coords.pm_dec, None)

        # 1975.35 = May 9, 1975
        # 1975.35 - 2000 = -24.65 years
        # 269.452075 + (-0.79858 * -24.65) / 3600 = 269.457543055
        #   4.693391 + (10.32812 * -24.65) / 3600 =   4.622672067
        coords = barnard.get_exact_coordinates(1975.35)
        self.assertAlmostEqual(coords.ra, 269.457543055)
        self.assertAlmostEqual(coords.dec,  4.622672067)
        self.assertIs(coords.pm_ra,  None)
        self.assertIs(coords.pm_dec, None)

        # Kapteyn's Star (J1950): 6505.08 -5730.84 (mas/yr)
        kapteyn = Coordinates(77.791453, -44.938748, 6.50508, -5.73084)

        # 2008 - 1950 = 58 years
        # 77.791453  + ( 6.50508 * 58) / 3600 =  77.896257067
        # -44.938748 + (-5.73084 * 58) / 3600 = -45.0310782
        coords = kapteyn.get_exact_coordinates(2008, epoch = 1950)
        self.assertAlmostEqual(coords.ra,   77.896257067)
        self.assertAlmostEqual(coords.dec, -45.0310782)
        self.assertIs(coords.pm_ra,  None)
        self.assertIs(coords.pm_dec, None)

        # 1905.49180328 = June 30, 1905
        # 1905.49180328 - 1950 = -44.50819672
        # 77.791453  + ( 6.50508 * -44.50819672) / 3600 =  77.711028172
        # -44.938748 + (-5.73084 * -44.50819672) / 3600 = -44.867895402
        coords = kapteyn.get_exact_coordinates(1905.49180328, epoch = 1950)
        self.assertAlmostEqual(coords.ra,   77.711028172)
        self.assertAlmostEqual(coords.dec, -44.867895402)
        self.assertIs(coords.pm_ra,  None)
        self.assertIs(coords.pm_dec, None)

        # IOK 1 (z = 6.96, 12.88 Gly)
        # No proper motion; object does not move
        iok1 = Coordinates(200.999170, 27.415500)
        coords = iok1.get_exact_coordinates(2061)
        self.assertEqual(coords.ra,  iok1.ra)
        self.assertEqual(coords.dec, iok1.dec)
        self.assertIs(coords.pm_ra,  None)
        self.assertIs(coords.pm_dec, None)


class StarTest(unittest.TestCase):

    X_COORD_RANGE = (1, 2048)
    Y_COORD_RANGE = (1, 2048)
    RIGHT_ASCENSION_RANGE = (0, 360)
    DECLINATION_RANGE = (-90, 90)
    ISOAREAF_RANGE = (1, 100)
    MAGNITUDE_RANGE = (1.47, 25)
    SNR_RANGE = (2, 10000)
    FWHM_RANGE = (0.5, 4.5)
    ELONGATION_RANGE = (1, 3.5)

    @classmethod
    def rargs(cls):
        """ Return the arguments needed to instantiate a random Star """

        x = random.uniform(*cls.X_COORD_RANGE)
        y = random.uniform(*cls.Y_COORD_RANGE)
        alpha = random.uniform(*cls.RIGHT_ASCENSION_RANGE)
        delta = random.uniform(*cls.DECLINATION_RANGE)
        isoareaf = random.uniform(*cls.ISOAREAF_RANGE)
        mag = random.uniform(*cls.MAGNITUDE_RANGE)
        saturated = random.choice([True, False])
        snr = random.uniform(*cls.SNR_RANGE)
        fwhm = random.uniform(*cls.FWHM_RANGE)
        elong = random.uniform(*cls.ELONGATION_RANGE)
        return x, y, alpha, delta, isoareaf, mag, saturated, snr, fwhm, elong

    @classmethod
    def random(cls):
        """ Return a random Star object """
        return Star(*cls.rargs())

    def test_init(self):
        for _ in xrange(NITERS):
            args = self.rargs()
            star = Star(*args)

            # Make sure the 'x', 'y', 'alpha' and 'delta' properties work
            x, y, ra, dec, area, mag, satur, snr, fwhm, elong = args
            self.assertEqual(star.x, x)
            self.assertEqual(star.y, y)
            self.assertEqual(star.alpha, ra)
            self.assertEqual(star.delta, dec)
            self.assertEqual(star.area, area)
            self.assertEqual(star.mag, mag)
            self.assertEqual(star.saturated, satur)
            self.assertEqual(star.snr, snr)
            self.assertEqual(star.fwhm, fwhm)
            self.assertEqual(star.elongation, elong)

    def test_immutability(self):
        """ Make sure that the Star class is immutable """

        star = self.random()
        for name in star._asdict().iterkeys():
            value = random.random()
            self.assertRaises(AttributeError, setattr, star, name, value)
            self.assertRaises(AttributeError, delattr, star, name)

    def test_distance(self):

        args1 = list(self.rargs())
        args1[0:2] = 134, 345
        star1 = Star(*args1)

        args2 = list(self.rargs())
        args2[0:2] = 178, 471
        star2 = Star(*args2)

        # ((134 - 178) ^ 2 + (345 - 471) ^ 2) ^ 0.5
        distance = star1.distance(star2)
        self.assertAlmostEqual(distance, 133.46160496562297)

        for _ in xrange(NITERS):

            star1 = self.random()
            star2 = self.random()
            distance = star1.distance(star2)

            # Compute the expected Euclidean distance with NumPy
            array1 = numpy.array([star1.x, star1.y])
            array2 = numpy.array([star2.x, star2.y])
            expected = numpy.linalg.norm(array1 - array2)
            self.assertAlmostEqual(distance, expected)


class CatalogTest(unittest.TestCase):

    args = os.path.dirname(__file__), './test_data/'
    get_path = functools.partial(os.path.join, *args)
    SAMPLE_CATALOG_PATH = get_path('sextractor.cat')
    SAMPLE_INCOMPLETE_PATH = get_path('sextractor_incomplete.cat')
    SAMPLE_NOASCIIHEAD_PATH = get_path('sextractor_noasciihead.cat')

    def test_flag_saturated(self):

        # Numbers in the [0, 255] range whose 3rd least significant bit is one
        saturated_flags = set([4, 5, 6, 7, 12, 13, 14, 15, 20, 21, 22, 23, 28,
        29, 30, 31, 36, 37, 38, 39, 44, 45, 46, 47, 52, 53, 54, 55, 60, 61, 62,
        63, 68, 69, 70, 71, 76, 77, 78, 79, 84, 85, 86, 87, 92, 93, 94, 95,
        100, 101, 102, 103, 108, 109, 110, 111, 116, 117, 118, 119, 124, 125,
        126, 127, 132, 133, 134, 135, 140, 141, 142, 143, 148, 149, 150, 151,
        156, 157, 158, 159, 164, 165, 166, 167, 172, 173, 174, 175, 180, 181,
        182, 183, 188, 189, 190, 191, 196, 197, 198, 199, 204, 205, 206, 207,
        212, 213, 214, 215, 220, 221, 222, 223, 228, 229, 230, 231, 236, 237,
        238, 239, 244, 245, 246, 247, 252, 253, 254, 255])

        for flag in xrange(0, 256):
            is_saturated = Catalog.flag_saturated(flag)
            if flag in saturated_flags:
                self.assertTrue(is_saturated)
            else:
                self.assertFalse(is_saturated)

        # Flags outside of the range raise ValueError
        self.assertRaises(ValueError, Catalog.flag_saturated, -2)
        self.assertRaises(ValueError, Catalog.flag_saturated, -1)
        self.assertRaises(ValueError, Catalog.flag_saturated, 256)
        self.assertRaises(ValueError, Catalog.flag_saturated, 257)

    def test_init(self):

        path = self.SAMPLE_CATALOG_PATH
        catalog = Catalog(path)
        self.assertEqual(catalog.path, path)
        self.assertEqual(len(catalog), 127)

        star = catalog[0]
        self.assertEqual(star.x, 844.359)
        self.assertEqual(star.y, 434.109)
        self.assertEqual(star.alpha, 100.2910553)
        self.assertEqual(star.delta, +9.4697032)
        self.assertEqual(star.area, 8085)
        self.assertEqual(star.mag, 8.2460)
        self.assertFalse(star.saturated)
        # S/N = FLUX_ISO / FLUXERR_ISO = 8545557 / 12811.09
        self.assertAlmostEqual(star.snr, 667.043709786)
        self.assertAlmostEqual(star.fwhm, 12.668) # FLUX_RADIUS x 2
        self.assertEqual(star.elongation, 1.562)

        star = catalog[68]
        self.assertEqual(star.x, 4.189)
        self.assertEqual(star.y, 1995.196)
        self.assertEqual(star.alpha, 100.1844554)
        self.assertEqual(star.delta, +9.2744660)
        self.assertEqual(star.area, 16)
        self.assertEqual(star.mag, 14.8506)
        self.assertFalse(star.saturated)
        # S/N = FLUX_ISO / FLUXERR_ISO = 10539.52 / 569.9098
        self.assertAlmostEqual(star.snr, 18.49331245)
        self.assertAlmostEqual(star.fwhm, 2.126) # FLUX_RADIUS x 2
        self.assertEqual(star.elongation, 1.339)

        star = catalog[126]
        self.assertEqual(star.x, 107.033)
        self.assertEqual(star.y, 1953.500)
        self.assertEqual(star.alpha, 100.1971009)
        self.assertEqual(star.delta, +9.2796855)
        self.assertEqual(star.area, 10)
        self.assertEqual(star.mag, 15.9567)
        self.assertTrue(star.saturated)
        # S/N = FLUX_ISO / FLUXERR_ISO = 3833.773 / 450.5532
        self.assertAlmostEqual(star.snr, 8.509035115)
        self.assertAlmostEqual(star.fwhm, 2.412) # FLUX_RADIUS x 2
        self.assertEqual(star.elongation, 1.613)

        # IOError is raised if the SExtractor catalog does not exist; in order
        # to get the path to a non-existent file, we use NamedTemporaryFile and
        # close the file (deleting it) immediately after.
        with tempfile.NamedTemporaryFile(suffix = '.cat') as fd: pass
        self.assertFalse(os.path.exists(fd.name))
        self.assertRaises(IOError, Catalog, fd.name)

        # ValueError is raised if one of the required SExtractor parameters is
        # not in the catalog, or if it was not saved in the ASCII_HEAD format
        # and therefore there are no comment lines listing column labels, which
        # are needed in order to determine in which column each parameter is.
        self.assertRaises(ValueError, Catalog, self.SAMPLE_INCOMPLETE_PATH)
        self.assertRaises(ValueError, Catalog, self.SAMPLE_NOASCIIHEAD_PATH)

    def test_immutability(self):
        """ Make sure that the Catalog class is immutable """

        catalog = Catalog(self.SAMPLE_CATALOG_PATH)
        with self.assertRaises(AttributeError):
            catalog.path = self.SAMPLE_INCOMPLETE_PATH
        with self.assertRaises(AttributeError):
            del catalog.path

        for index in xrange(len(catalog)):
            star = StarTest.random()
            with self.assertRaises(TypeError):
                catalog[index] = star
            with self.assertRaises(TypeError):
                del catalog[index]

    def test_from_sequence(self):

        # Load a SExtractor catalog into memory and pass all of its stars to
        # Catalog.from_sequence(): both Catalogs must be equal. This is not the
        # case if the order of the stars passed to from_sequence() is altered,
        # of course, as tuples are compared position by position.

        catalog = Catalog(self.SAMPLE_CATALOG_PATH)
        stars = catalog[:]

        identical = Catalog.from_sequence(*stars)
        self.assertEqual(type(identical), Catalog)
        self.assertEqual(identical, catalog)

        different = Catalog.from_sequence(*stars[::-1])
        self.assertEqual(type(different), Catalog)
        self.assertNotEqual(different, catalog)

        # Now create a second Catalog containing those stars whose magnitude is
        # greater than or equal to sixteen. There are nineteen of these stars,
        # but this number is anyway irrelevant for our purposes. What matters
        # is that all the stars passed to from_sequence() must also be in the
        # returned Catalog, in the same order.

        stars = [star for star in catalog if star.mag > 16]
        self.assertEqual(len(stars), 19)
        faint_catalog = Catalog.from_sequence(*stars)
        self.assertEqual(type(faint_catalog), Catalog)
        self.assertEqual(list(faint_catalog), list(stars))

    def test_get_sky_coordinates(self):

        catalog = Catalog(self.SAMPLE_CATALOG_PATH)
        self.assertEqual(len(catalog), 127)
        coordinates = catalog.get_sky_coordinates()

        # There should be as many Coordinates objects, and in the same order,
        # as sources there are in the SExtractor catalog. Therefore, the right
        # ascension and declination of the i-th Coordinates object must match
        # those of the i-th source in the Catalog.

        self.assertEqual(len(coordinates), len(catalog))
        for coord, star in zip(coordinates, catalog):
            self.assertEqual(coord.ra, star.alpha)
            self.assertEqual(coord.dec, star.delta)


def get_nonexistent_path(ext = None):
    """ Return the path to a nonexistent file """
    with tempfile.NamedTemporaryFile(suffix = ext) as fd:
        return fd.name


class SExtractorFunctionsTest(unittest.TestCase):

    SEXTRACTOR_MODULE_VARS = \
      ('SEXTRACTOR_CONFIG',
       'SEXTRACTOR_PARAMS',
       'SEXTRACTOR_FILTER',
       'SEXTRACTOR_STARNNW')

    def test_sextractor_md5sum(self):

        # Testing that astromatic.sextractor_md5sum() works as expected means,
        # by extension, checking that the four SExtractor configuration files
        # (.sex, .param, .conv and .nnw) defined in the 'astromatic' module
        # exist and are readable.

        # If the SExtractor configuration files are not modified in between,
        # two different executions of the function must yield the same hash.
        checksum  = astromatic.sextractor_md5sum()
        identical = astromatic.sextractor_md5sum()
        self.assertEqual(checksum, identical)

        # The MD5 hash is that of the concatenation of the lines of the four
        # configuration files and the overriding SExtractor options (i.e., the
        # sequence of strings that are optionally passed to the method). If any
        # of the files is modified, or if an overriding option is given, the
        # MD5 hash returned by the function must be different.
        #
        # In order to test this we are not going to temporarily modify the
        # SExtractor configuration files: although we could make a copy of
        # them, something beyond our control (e.g., the SIGKILL signal or a
        # power outage) could prevent the original file from being restored.
        # Instead, what we will temporarily modify are the module-level
        # variables, so that they refer to a modified copy of the files.

        for variable in self.SEXTRACTOR_MODULE_VARS:

            # The variable must exist and refer to an existing file
            path = eval('astromatic.%s' % variable)
            self.assertTrue(os.path.exists(path))

            # Make a temporary copy of the configuration file and modify it,
            # appending a comment to it. Then, mock the module-level variable
            # so that it refers to the modified configuration file. Although
            # such an irrelevant change does not alter how SExtractor works (as
            # the settings are still the same), the MD5 hash must be different.

            ext = os.path.splitext(path)[1]
            copy_path = get_nonexistent_path(ext = ext)

            try:
                shutil.copy2(path, copy_path)
                with open(copy_path, 'at') as fd:
                    fd.write("# useless comment\n")

                with mock.patch.object(astromatic, variable, copy_path):
                    different = astromatic.sextractor_md5sum()
                    self.assertNotEqual(checksum, different)

            finally:
                os.unlink(copy_path)

        # If overriding options are given, they are also used to compute the
        # MD5 hash. Thus, the hash will be different even if the options have
        # the same value as those defined in the configuration file (although
        # that means that we are not actually overriding anything). Test this
        # for all the options defined in the configuration file: in all cases
        # the returned checksum must be different.

        with open(astromatic.SEXTRACTOR_CONFIG) as fd:

            for line in fd:
                stripped = line.strip()
                if stripped and stripped[0] != '#': # ignore comment lines
                    key, value = stripped.split()[:2]
                    options = {key : value}
                    different = astromatic.sextractor_md5sum(options)
                    self.assertNotEqual(checksum, different)

                    # Try also overriding the option with a different value
                    # (the result of incrementing it by one). Again, this must
                    # result in a different MD5 hash.

                    try:
                        options[key] = str(int(options[key]) + 1)
                        also_different = astromatic.sextractor_md5sum(options)
                        self.assertNotEqual(checksum,  also_different)
                        self.assertNotEqual(different, also_different)

                    except ValueError:
                        pass # non-numeric option

        # IOError is raised if any of the SExtractor configuration files is not
        # readable or does not exist. To test that, mock again the module-level
        # variables so that they refer to temporary copies of the configuration
        # files, but which are unreadable by the user, first, and then removed.

        for variable in self.SEXTRACTOR_MODULE_VARS:

            path = eval('astromatic.%s' % variable)
            ext = os.path.splitext(path)[1]
            copy_path = get_nonexistent_path(ext = ext)
            shutil.copy2(path, copy_path)

            with mock.patch.object(astromatic, variable, copy_path):

                # chmod u-r
                mode = stat.S_IMODE(os.stat(copy_path)[stat.ST_MODE])
                mode ^= stat.S_IRUSR
                os.chmod(copy_path, mode)

                args = IOError, astromatic.sextractor_md5sum
                self.assertFalse(os.access(copy_path, os.R_OK))
                self.assertRaises(*args)

                os.unlink(copy_path)
                self.assertFalse(os.path.exists(copy_path))
                self.assertRaises(*args)

        # TypeError raised if 'options' is not a dictionary
        kwargs = dict(options = ['DETECT_MINAREA', '5'])
        with self.assertRaises(TypeError):
            astromatic.sextractor_md5sum(**kwargs)

        # ... or if any of its elements is not a string
        kwargs['options'] = {'DETECT_MINAREA' : 125}
        with self.assertRaises(TypeError):
            astromatic.sextractor_md5sum(**kwargs)

    def test_sextractor_version(self):

        # We have no other way of knowing the SExtractor version that is
        # installed on the system than running it ourselves with the --version
        # option and examine its output: if sextractor_version() returns the
        # tuple (2, 8, 6), for example, that means that the SExtractor version
        # number must contain the string '2.8.6'. What we are doing, basically,
        # is to test the functionality the other way around: we take the tuple
        # that sextractror_version() outputs, transform it back to a string and
        # verify that it corresponds to what the --version option prints.

        try:
            executable = methods.which(*astromatic.SEXTRACTOR_COMMANDS)[0]
        except IndexError:
            msg = "SExtractor not found in the current environment"
            raise astromatic.SExtractorNotInstalled(msg)

        with tempfile.TemporaryFile() as fd:
            args = [executable, '--version']
            subprocess.check_call(args, stdout = fd)
            fd.seek(0)
            # For example: "SExtractor version 2.5.0 (2006-07-14)"
            stdout = '\n'.join(fd.readlines())

        version = astromatic.sextractor_version()
        # From, for example, (2, 5, 0) to '2.5.0'
        version_str = '.'.join(str(x) for x in version)
        self.assertTrue(version_str in stdout)

        # SExtractorNotInstalled must be raised if no SExtractor executable is
        # detected. In order to simulate this, mock os.environ and clear the
        # 'PATH' environment variable. In this manner, as the list of paths to
        # directories where executables may be found is empty, SExtractor (and
        # any other command) will appear as not installed on the system.

        environment_copy = copy.deepcopy(os.environ)
        with mock.patch.object(os, 'environ', environment_copy) as mocked:
            mocked['PATH'] = ''
            with self.assertRaises(astromatic.SExtractorNotInstalled):
                astromatic.sextractor_version()

    def test_sextractor(self):

        # Note that, being a third-party program, we cannot write a unit test
        # to check that SExtractor works as it should (i.e., that it correctly
        # detects astronomical sources on FITS images): that is something that
        # Emmanuel Bertin, its developer, certainly takes care of. What we must
        # test for is whether astromatic.sextractor(), our Python wrapper to
        # call SExtractor, works as we expect and raises the appropriate errors
        # when something goes wrong. Nothing less, nothing more.

        # First, the most probable scenario: our wrapper should run SExtractor
        # on any FITS image that it receives (here we use some downloaded from
        # the STScI Digitized Sky Survey) and return the path to the output
        # catalog, which astromatic.Catalog should be always able to parse.

        kwargs = dict(stdout = open(os.devnull), stderr = open(os.devnull))
        for img_path in dss_images.TEST_IMAGES:
            catalog_path = astromatic.sextractor(img_path, **kwargs)
            try:
                Catalog(catalog_path)
            finally:
                os.unlink(catalog_path)

        # SExtractorNotInstalled must be raised if no SExtractor executable is
        # detected. In order to simulate this, mock os.environ and clear the
        # 'PATH' environment variable. In this manner, as the list of paths to
        # directories where executables may be found is empty, SExtractor (and
        # any other command) will appear as not installed on the system.

        environment_copy = copy.deepcopy(os.environ)
        with mock.patch.object(os, 'environ', environment_copy) as mocked:
            mocked['PATH'] = ''
            with self.assertRaises(astromatic.SExtractorNotInstalled):
                astromatic.sextractor(img_path)

        # IOError is raised if any of the four SExtractor configuration files
        # is not readable or does no exist. To test that, mock the module-level
        # variables so that they first refer to an unreadable file and later to
        # a nonexistent one.

        for variable in self.SEXTRACTOR_MODULE_VARS:

            path = eval('astromatic.%s' % variable)
            ext = os.path.splitext(path)[1]
            copy_path = get_nonexistent_path(ext = ext)
            shutil.copy2(path, copy_path)

            with mock.patch.object(astromatic, variable, copy_path):

                # chmod u-r
                mode = stat.S_IMODE(os.stat(copy_path)[stat.ST_MODE])
                mode ^= stat.S_IRUSR
                os.chmod(copy_path, mode)

                args = IOError, astromatic.sextractor, img_path
                self.assertFalse(os.access(copy_path, os.R_OK))
                self.assertRaises(*args)

                os.unlink(copy_path)
                self.assertFalse(os.path.exists(copy_path))
                self.assertRaises(*args)

        # SExtractorUpgradeRequired must be raised if the version of SExtractor
        # that is installed on the system is older than that defined in the
        # SEXTRACTOR_REQUIRED_VERSION module-level variable. We cannot (easily,
        # at least) change the version that is installed, but we can test for
        # this by changing the value of the required version, setting it to a
        # tuple that is greater than the installed version of SExtractor.

        # From, for example, (2, 8, 6) to (2, 8, 7)
        version = list(astromatic.sextractor_version())
        version[-1] += 1
        version = tuple(version)

        with mock.patch.object(astromatic, 'SEXTRACTOR_REQUIRED_VERSION', version):
            with self.assertRaises(astromatic.SExtractorUpgradeRequired):
                astromatic.sextractor(img_path)

        # The SExtractorError exception is raised if anything goes wrong during
        # the execution of SExtractor (in more technical terms, if its return
        # code is other than zero). For example, the FITS image may not exist,
        # or a non-numerical value may be assigned to a parameter that expects
        # one, or we could try to run SExtractor on a FITS extension that does
        # not exist. These three are just examples: SExtractor may fail for
        # many different reasons that we cannot even foresee.

        # (1) Try to run SExtractor on a non-existent image
        kwargs = dict(stdout = open(os.devnull), stderr = open(os.devnull))
        nonexistent_path = get_nonexistent_path(ext = '.fits')
        with self.assertRaises(astromatic.SExtractorError):
            astromatic.sextractor(nonexistent_path, **kwargs)

        # (2) The DETECT_MINAREA parameter (minimum number of pixels above the
        # threshold needed to trigger detection) expects an integer. Assigning
        # to it a string will cause SExtractor to complain ("keyword out of
        # range") and abort its execution.

        kwargs['options'] = dict(DETECT_MINAREA = 'Y')
        with self.assertRaises(astromatic.SExtractorError):
            astromatic.sextractor(img_path, **kwargs)
        del kwargs['options']

        # (3) Try to run SExtractor on a nonexistent FITS extension.
        hdulist = pyfits.open(img_path, mode = 'readonly')
        nextensions = len(hdulist)
        hdulist.close()

        kwargs['ext'] = nextensions + 1
        with self.assertRaises(astromatic.SExtractorError):
            astromatic.sextractor(img_path, **kwargs)

        # TypeError raised if 'options' is not a dictionary
        kwargs = dict(options = ['DETECT_MINAREA', '5'])
        with self.assertRaises(TypeError):
            astromatic.sextractor(img_path, **kwargs)

        # ... or if any of its elements is not a string.
        kwargs['options'] = {'DETECT_MINAREA' : 125}
        with self.assertRaises(TypeError):
            astromatic.sextractor(img_path, **kwargs)

        # TypeError also raised if 'ext' is not an integer...
        kwargs = dict(ext = 0.56)
        with self.assertRaises(TypeError):
            astromatic.sextractor(img_path, **kwargs)

        # ... even if it is a float but has nothing after the decimal point
        # (for which the built-in is_integer() function would return True).
        kwargs = dict(ext = 0.0)
        with self.assertRaises(TypeError):
            astromatic.sextractor(img_path, **kwargs)

