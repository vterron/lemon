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

import collections
import copy
import itertools
import numpy
import operator
import os
import random
import re
import sqlite3
import string
import tempfile
import time

from test import unittest
import passband
from database import \
  (DBStar,
   DuplicateLightCurvePointError,
   DuplicateImageError,
   DuplicatePeriodError,
   DuplicatePhotometryError,
   DuplicateStarError,
   Image,
   LEMONdB,
   LightCurve,
   PhotometricParameters,
   UnknownImageError,
   UnknownStarError)

from diffphot import Weights
from json_parse import CandidateAnnuli
import test.test_fitsimage
# https://stackoverflow.com/q/12603541/184363
from test.test_astromatic import CoordinatesTest
get_random_coords = CoordinatesTest.random

NITERS = 100      # How many times each test case is run with random data
MIN_NSTARS = 25   # Minimum number of items for random collections of DBStars
MAX_NSTARS = 100  # Maximum number of items for random collections of DBStars
MIN_UNIX_TIME = 0            # Thu Jan  1 01:00:00 1970
MAX_UNIX_TIME = time.time()  # Minimum and maximum random Unix times

def assertSequenceOfTuplesAlmostEqual(cls, first, second):
    """ Assert that the tuples of the two sequences are almost equal.

    The first parameter, 'cls', must be the subclass of unittest.TestCase
    instance in which this test is done. Values are considered to be almost
    equal according to the rules of the TestCase.assertAlmostEqual method.

    """

    if len(first) != len(second):
        cls.fail("sequences must be of the same size")
    for tuple1, tuple2 in zip(first, second):
        if len(tuple1) != len(tuple2):
            cls.fail("tuples must be of the same size")
        for values in zip(tuple1, tuple2):
            cls.assertAlmostEqual(*values)

def runix_times(size):
    """ Return a list of 'size' random, different Unix times """
    rtimes = []
    while len(rtimes) < size:
        utime = random.uniform(MIN_UNIX_TIME, MAX_UNIX_TIME)
        if utime not in rtimes:
            rtimes.append(utime)
    return rtimes

def different_runix_time(utimes):
    """ Return a random Unix time not in 'utimes' """
    while True:
        t = ImageTest.random().unix_time
        if t not in utimes:
            return t


class DBStarTest(unittest.TestCase):

    MIN_ID = 1     # Minimum ID for random DBStars
    MAX_ID = 9999  # Maximum ID for random DBStars
    MIN_SIZE = 1   # Minimum number of photometric records
    MAX_SIZE = 250 # Maximum number of photometric records
    MIN_MAG = 1.47   # Minimum value for random magnitudes (Saturn)
    MAX_MAG = 25     # Maximum value for random magnitudes (Fenrir)
    MIN_SNR = 2      # Minimum value for random signal-to-noise ratios
    MAX_SNR = 10000  # Maximum value for random signal-to-noise ratios

    @classmethod
    def random_data(cls, pfilter = None):
        """ Return the information needed to instantiate a random DBStar.

        If 'pfilter' is given, it is used as the photometric filter in the
        returned tuple (the second element), instead of a random one.

        """

        id_ = random.randint(cls.MIN_ID, cls.MAX_ID)
        if not pfilter:
            pfilter = passband.Passband.random()
        size = random.randint(cls.MIN_SIZE, cls.MAX_SIZE)
        phot_info = numpy.empty((3, size))
        times_indexes = {}
        for index, unix_time in enumerate(runix_times(size)):
            magnitude = random.uniform(cls.MIN_MAG, cls.MAX_MAG)
            snr = random.uniform(cls.MIN_SNR, cls.MAX_SNR)
            phot_info[:,index] = unix_time, magnitude, snr
            times_indexes[unix_time] = index
        return id_, pfilter, phot_info, times_indexes

    @classmethod
    def random(cls, pfilter = None):
        """ Return a random DBStar.

        If 'pfilter' is given, it is used as the photometric filter of the
        DBStar, instead of a random one.

        """
        args = cls.random_data(pfilter = pfilter)
        return DBStar(*args)

    def test_init(self):
        for _ in xrange(NITERS):
            id_, pfilter, phot_info, times_indexes = DBStarTest.random_data()
            star = DBStar(id_, pfilter, phot_info, times_indexes)
            # Test that the attributes are correctly set at instantiation time
            self.assertEqual(star.id, id_)
            self.assertEqual(star.pfilter, pfilter)
            self.assertTrue(numpy.all(numpy.equal(star._phot_info, phot_info)))
            self.assertEqual(star._time_indexes, times_indexes)

        # ValueError must be raised if 'phot_info' does not have three rows
        id_, pfilter, phot_info, times_indexes = DBStarTest.random_data()
        row_index = random.randint(0, 2)
        phot_info = numpy.delete(phot_info, row_index, axis = 0)
        with self.assertRaises(ValueError):
            DBStar(id_, pfilter, phot_info, times_indexes)

    def test_len(self):
        for _ in xrange(NITERS):
            id_, pfilter, phot_info, times_indexes = DBStarTest.random_data()
            size = phot_info.shape[1] # as many photometric records as columns
            star = DBStar(id_, pfilter, phot_info, times_indexes)
            self.assertEqual(len(star), size)

    def test_time_mag_and_snr(self):
        for _ in xrange(NITERS):
            id_, pfilter, phot_info, times_indexes = DBStarTest.random_data()
            star = DBStar(id_, pfilter, phot_info, times_indexes)
            for index in xrange(len(star)):
                unix_time, magnitude, snr = phot_info[:, index]
                self.assertEqual(star.time(index), unix_time)
                self.assertEqual(star.mag(index), magnitude)
                self.assertEqual(star.snr(index), snr)

    def test_time_index(self):
        for _ in xrange(NITERS):
            id_, pfilter, phot_info, times_indexes = DBStarTest.random_data()
            star = DBStar(id_, pfilter, phot_info, times_indexes)
            for unix_time in times_indexes.iterkeys():
                unix_time_index = star._time_index(unix_time)
                self.assertEqual(unix_time_index, times_indexes[unix_time])

    def test_unix_times(self):
        for _ in xrange(NITERS):
            id_, pfilter, phot_info, times_indexes = DBStarTest.random_data()
            star = DBStar(id_, pfilter, phot_info, times_indexes)
            self.assertTrue(numpy.all(numpy.equal(star._unix_times, phot_info[0])))

    @classmethod
    def make_subset(cls, star):
        """ Return a DBstar which is a subset of 'star'. In other words: the
        set of Unix times of the returned DBStar star is a subset of the Unix
        times of 'star'."""

        # For the subset star we make a copy of 'phot_info' and remove from it
        # a random number (which may be zero) of columns, but always leaving at
        # least one column.
        sid = random.randint(cls.MIN_ID, cls.MAX_ID)
        to_delete = random.randint(0, len(star) - 1)
        delete_indexes = random.sample(range(len(star)), to_delete)
        sphot_info = numpy.delete(star._phot_info, delete_indexes, axis = 1)

        # Populate 'sphot_info' with other magnitudes and SNRs
        for index in xrange(sphot_info.shape[1]):
            magnitude = random.uniform(cls.MIN_MAG, cls.MAX_MAG)
            snr = random.uniform(cls.MIN_SNR, cls.MAX_SNR)
            sphot_info[1:, index] = magnitude, snr

        stimes_indexes = {}
        for index, unix_time in enumerate(sphot_info[0, :]):
            stimes_indexes[unix_time] = index

        subset = DBStar(sid, star.pfilter, sphot_info, stimes_indexes)
        assert 0 < len(subset) <= len(star)
        assert subset.issubset(star)
        return subset

    @classmethod
    def make_superset(cls, star):
        """ Return a DBstar which is a superset of 'star'. This means that the
        set of Unix times of the returned DBStar star is a subset of the Unix
        times of the 'star'."""

        # For the superset star we make a copy of 'phot_info' and add to it
        # a random number (which may be zero) of columns. At most we will
        # duplicate the size of the star.
        sid = random.randint(cls.MIN_ID, cls.MAX_ID)
        to_add = random.randint(0, len(star))
        sphot_info = numpy.copy(star._phot_info)
        stimes_indexes = star._time_indexes.copy()

        for unix_time in runix_times(to_add):
            magnitude = random.uniform(cls.MIN_MAG, cls.MAX_MAG)
            snr = random.uniform(cls.MIN_SNR, cls.MAX_SNR)
            row = [unix_time, magnitude, snr]

            size = sphot_info.shape[1]
            sphot_info = numpy.insert(sphot_info, size, row, axis = 1)
            stimes_indexes[unix_time] = size

        # There must be no gaps in the indexes mapped by the dictionary
        assert sorted(stimes_indexes.values()) == range(sphot_info.shape[1])

        superset = DBStar(sid, star.pfilter, sphot_info, stimes_indexes)
        assert star.issubset(superset)
        return superset

    def test_issubset(self):
        for _ in xrange(NITERS):
            # Any star is a subset of itself
            rstar = DBStarTest.random()
            self.assertTrue(rstar.issubset(rstar))

            original = DBStarTest.random()
            subset = DBStarTest.make_subset(original)
            superset = DBStarTest.make_superset(original)

            self.assertTrue(subset.issubset(original))
            self.assertTrue(original.issubset(superset))
            self.assertTrue(subset.issubset(superset)) # transitive relation

            self.assertTrue(subset.issubset(subset))
            self.assertTrue(original.issubset(original))
            self.assertTrue(superset.issubset(superset))

            # If a star has the same size than its superset (i.e., if no
            # columns were removed by DBSTarTest.make_subset/superset), they
            # are effectively equal -- and thus the superset is also a subset
            # of its subset, as any star is a subset of itself
            stars = (superset, original, subset)
            for first, second in itertools.combinations(stars, 2):
                if len(first) == len(second):
                    func = self.assertTrue
                else:
                    func = self.assertFalse
                func(first.issubset(second))

            # Modify one of the Unix times of 'subset', so that it is no longer
            # a subset of 'original', and verify that the method returns False
            index = random.randint(0, len(subset) - 1)
            while True: # +/- one hour, until we get a different Unix time
                unix_time = subset.time(index) + random.uniform(-3600, 3600)
                if unix_time != subset.time(index):
                    break

            # Update the Unix time (a key) in the 'time_indexes' dictionary
            time_indexes = subset._time_indexes # shorter variable name
            time_indexes[unix_time] = time_indexes.pop(subset.time(index))
            subset._phot_info[0][index] = unix_time
            assert subset.time(index) == unix_time
            self.assertFalse(subset.issubset(original))

    @classmethod
    def make_star(cls, records):
        """ Return a DBStar whose photometric information (Unix time, magnitude
        and SNR of each record) is taken from the tree-element tuples contained
        in 'records'). The rest of the information is randomly generated """

        if not records:
            raise ValueError("'records' cannot be empty")
        id_ = random.randint(cls.MIN_ID, cls.MAX_ID)
        pfilter = passband.Passband.random()
        return DBStar.make_star(id_, pfilter, records)

    def test_trim_to(self):

        # A specific, non-random test case...
        original = DBStarTest.make_star(
        [(13000, 15.6, 100),
         (13100, 14.5, 230),
         (14000, 13.4, 200),
         (16700, 15.3, 250)])

        subset = DBStarTest.make_star(
        [(13100, 19.3, 85),
         (14000, 17.4, 115)])

        self.assertTrue(subset.issubset(original))
        trimmed = original._trim_to(subset)
        self.assertEqual(len(trimmed), 2)
        self.assertEqual(trimmed.time(0), 13100)
        self.assertEqual(trimmed.time(1), 14000)
        self.assertTrue(trimmed.issubset(subset))

        # ... and the random test cases
        for _ in xrange(NITERS):
            original = DBStarTest.random()
            subset = DBStarTest.make_subset(original)
            trimmed = original._trim_to(subset)
            self.assertTrue(trimmed.issubset(subset))

            # The Unix times in both DBStars must now be equal...
            self.assertEqual(len(trimmed), len(subset))
            for index in xrange(len(trimmed)):
                unix_time = trimmed.time(index)
                self.assertEqual(unix_time, subset.time(index))
                # ... and each Unix time in 'trimmed' must be associated to the
                # right magnitude and SNR -- those which it had in 'original'.
                oindex = original._time_index(unix_time)
                self.assertEqual(trimmed.mag(index), original.mag(oindex))
                self.assertEqual(trimmed.snr(index), original.snr(oindex))

            # The internal dictionaries mapping each Unix time to its position
            # in '_phot_info' must be identical, as both stars contain now
            # information for exactly the same Unix times, in the same order.
            self.assertTrue(numpy.all(numpy.equal(trimmed._unix_times,
                                                  subset._unix_times)))

        # KeyError is raised if we attempt to trim a star which is not a subset
        star = DBStarTest.random()
        while True:
            nonsubset = DBStarTest.random()
            if not nonsubset.issubset(star):
                break
        with self.assertRaises(KeyError):
            star._trim_to(nonsubset)

    def test_complete_for(self):
        for _ in xrange(NITERS):
            original = DBStarTest.random()
            stars = [original]
            super_ids = [] # the IDs of the superset stars

            while len(stars) < random.randint(MIN_NSTARS, MAX_NSTARS):
                # Half of the stars in the collection will be supersets
                if random.choice((True, False)):
                    superset = DBStarTest.make_superset(original)
                    super_ids.append(superset.id)
                    stars.append(superset)
                else:
                    while True:
                        nonsuperset = DBStarTest.random()
                        if not original.issubset(nonsuperset):
                            break
                    stars.append(nonsuperset)

            complete = original.complete_for(stars)
            complete_ids = [star.id for star in complete]
            self.assertEqual(sorted(complete_ids), sorted(super_ids))
            for cstar in complete:
                self.assertTrue(original.issubset(cstar))

    def test_make_star(self):

        id_ = 1
        pfilter = passband.Passband('V')

        row0 = (13000, 15.6, 100)
        row1 = (14000, 14.5, 230)
        row2 = (13100, 13.4, 200)
        records = [row0, row1, row2]

        # Note how rows are not given sorted by their Unix time -- make_star
        # does not pay attention to it; it just saves them in the internal
        # NumPy array in the order in which they are in 'records'
        star = DBStar.make_star(id_, pfilter, records)
        self.assertEqual(star.id, id_)
        self.assertEqual(star.pfilter, pfilter)

        self.assertEqual(star.time(0), row0[0])
        self.assertEqual(star.mag(0),  row0[1])
        self.assertEqual(star.snr(0),  row0[2])

        self.assertEqual(star.time(1), row1[0])
        self.assertEqual(star.mag(1),  row1[1])
        self.assertEqual(star.snr(1),  row1[2])

        self.assertEqual(star.time(2), row2[0])
        self.assertEqual(star.mag(2),  row2[1])
        self.assertEqual(star.snr(2),  row2[2])

        # Now random test cases
        for _ in xrange(NITERS):
            id_, pfilter, phot_info, times_indexes = self.random_data()
            # Construct the row of three-element tuples, out of this random
            # NumPy array, and then check that the array in the DBStar returned
            # by DBStar.make_star is equal than this input, original array.
            rows = [phot_info[:, index] for index in xrange(phot_info.shape[1])]
            star = DBStar.make_star(id_, pfilter, rows)
            self.assertEqual(star.id, id_)
            self.assertEqual(star.pfilter, pfilter)
            self.assertTrue(numpy.all(numpy.equal(star._phot_info, phot_info)))
            self.assertEqual(star._time_indexes, times_indexes)

    @staticmethod
    def equal(first, second):
        """ Check whether two DBStars are equal """

        if (first.id != second.id or
            first.pfilter != second.pfilter or
            len(first) != len(second)):
              return False

        for index in xrange(len(first)):
            if (first.time(index) != second.time(index) or
                first.mag(index) != second.mag(index) or
                first.snr(index) != second.snr(index)):
                  return False

        return True


class PhotometricParametersTest(unittest.TestCase):

    MIN_APERTURE = 0.1
    MAX_APERTURE = 20
    MIN_ANNULUS = 0.1
    MAX_ANNULUS = 20
    MIN_DANNULUS = 0.1
    MAX_DANNULUS = 20

    @classmethod
    def random_data(cls):
        """ Return the information needed to instantiate a random instance"""

        # These values may be unrealistic, as the inner radius of the sky
        # annulus ('annulus') may be smaller than the aperture radius, but
        # it does not matter here. We only need to verify that the values
        # are correctly set at instantiation time, and the class must work
        # for all numbers anyway, whether realistic or not.
        aperture = random.uniform(cls.MIN_APERTURE, cls.MAX_APERTURE)
        annulus  = random.uniform(cls.MIN_ANNULUS, cls.MAX_ANNULUS)
        dannulus = random.uniform(cls.MIN_DANNULUS, cls.MAX_DANNULUS)
        return aperture, annulus, dannulus

    @classmethod
    def random(cls):
        """ Return a random PhotometricParameters """
        args = cls.random_data()
        return PhotometricParameters(*args)

    def test_init_(self):
        for _ in xrange(NITERS):
            aperture, annulus, dannulus = self.random_data()
            pparams = PhotometricParameters(aperture, annulus, dannulus)
            self.assertEqual(pparams.aperture, aperture)
            self.assertEqual(pparams.annulus, annulus)
            self.assertEqual(pparams.dannulus, dannulus)


class ImageTest(unittest.TestCase):

    MIN_AIRMASS = 1     # Minimum value for random airmasses (zenith)
    MAX_AIRMASS = 2.92  # Maximum value for random airmasses (~70 dec)
    MIN_GAIN = 1   # Minimum value for random CCD gains
    MAX_GAIN = 10  # Maximum value for random CCD gains
    RA_RANGE = (0, 359.99999)  # Minimum and maximum right ascensions (degrees)
    DEC_RANGE = (-90, 90) # Minimum and maximum declinations (degrees)
    OBJECT_LENGTH = (1, 30) # Minimum and maximum length of the random string
                            # used as the name of the observed object

    @classmethod
    def random_data(cls):
        """ Return the information needed to instantiate a random Image """

        fd, path = tempfile.mkstemp(suffix = '.fits')
        os.close(fd)
        os.unlink(path)

        pfilter = passband.Passband.random()
        unix_time = runix_times(1)[0]
        size = random.randint(*cls.OBJECT_LENGTH)
        object_ = ''.join(random.choice(string.letters) for _ in xrange(size))
        airmass = random.uniform(cls.MIN_AIRMASS, cls.MAX_AIRMASS)
        gain = random.uniform(cls.MIN_GAIN, cls.MAX_GAIN)
        ra = random.uniform(*cls.RA_RANGE)
        dec = random.uniform(*cls.DEC_RANGE)

        return (path, pfilter, unix_time, object_, airmass, gain, ra, dec)

    @classmethod
    def nrandom(cls, size, pfilter = None, simult_prob = 0.5):
        """ Return a generator which steps through 'size' random Images.

        The random Image instances may have the same Unix time, but in these
        cases they are guaranteed to have a different photometric filter, to
        simulate a campaign with simultaneous observations. In other words: in
        no case two returned Images may have the same Unix time and photometric
        filter.

        If 'pfilter' is given, it is set as the photometric filter of all the
        Images, instead of using a random one for each image. When this is the
        case it is not possible to have simultaneous observations (as all the
        Images have the same filter), so all the Unix times will be different.

        The 'simult_prob' keyword parameter defines the probability that each
        random generated image has the same Unix time, but a different filter,
        that one of the Images previously returned by the method. A value of
        0.5, for example, means that on average half of the random Images will
        share their Unix time with another Image.

        """

        # Map each filter to the Unix times that have been used for it
        used_utimes = collections.defaultdict(set)

        for index, unix_time in enumerate(runix_times(size)):
            args = list(cls.random_data())
            args[2] = unix_time

            img = Image(*args)

            # No simultaneus observations if all the images will have the same
            # photometric filter, or if there will be simultaneous ones but
            # this image is (a) either the first one, a case in which there is
            # no previous image whose Unix time we can reuse or (b) not chosen
            # randomly to have a duplicate Unix time.

            if pfilter:
                img = img._replace(pfilter = pfilter)

            elif bool(used_utimes) and random.random() < simult_prob:

                # Reuse the Unix time of one the Images already returned, but
                # we need to find a filter for which it has not been used yet.

                all_utimes = set().union(*used_utimes.itervalues())
                duplicate_utime = random.choice(list(all_utimes))

                all_pfilters = passband.Passband.all()
                random.shuffle(all_pfilters)
                for duplicate_pfilter in all_pfilters:
                    if duplicate_utime not in used_utimes[duplicate_pfilter]:
                        img = img._replace(unix_time = duplicate_utime)
                        img = img._replace(pfilter = duplicate_pfilter)
                        break

                # This point is reached only if all the photometric filters
                # have been already used for this Unix time. Therefore, we
                # cannot reuse the Unix time, so give up -- use the original,
                # unique Unix time that obtained with the runix_times function
                # and with which we have instantiated the Image object.

            pfilter_utimes = used_utimes[img.pfilter]
            assert img.unix_time not in pfilter_utimes
            pfilter_utimes.add(img.unix_time)
            yield img

        if __debug__:
            all_utimes = set().union(*used_utimes.itervalues())
            assert len(all_utimes) <= size

    @classmethod
    def random(cls, pfilter = None):
        """ Return a random Image.

        If 'pfilter' is given, it is used as the photometric filter, instead
        of a random one -- though it could be argued then that the Image is
        not completely 'random'.

        """

        return cls.nrandom(1, pfilter = pfilter).next()

    def test_init_(self):
        for _ in xrange(NITERS):
            args = self.random_data()
            img = Image(*args)
            self.assertEqual(img.path, args[0])
            self.assertEqual(img.pfilter, args[1])
            self.assertEqual(img.unix_time, args[2])
            self.assertEqual(img.object, args[3])
            self.assertEqual(img.airmass, args[4])
            self.assertEqual(img.gain, args[5])
            self.assertEqual(img.ra, args[6])
            self.assertEqual(img.dec, args[7])


class LightCurveTest(unittest.TestCase):

    MIN_ID = 1     # Minimum value for random DBStar IDs
    MAX_ID = 9999  # Maximum value for random DBStar IDs
    MIN_MAG = 1.47   # Minimum value for random magnitudes (Saturn)
    MAX_MAG = 25     # Maximum value for random magnitudes (Fenrir)
    MIN_SNR = 2      # Minimum value for random signal-to-noise ratios
    MAX_SNR = 10000  # Maximum value for random signal-to-noise ratios

    @classmethod
    def random_points(cls, size):
        """ Return a generator which steps through 'size' random (Unix
        time, magnitude, snr) tuples, guaranteed to have different times.

        """

        already_used = set()
        while len(already_used) != size:
            unix_time = runix_times(1)[0]
            if unix_time not in already_used:
                already_used.add(unix_time)
                magnitude = random.uniform(cls.MIN_MAG, cls.MAX_MAG)
                snr = random.uniform(cls.MIN_SNR, cls.MAX_SNR)
                yield unix_time, magnitude, snr

    @classmethod
    def random_point(cls):
        """ Return a random (Unix time, magnitude, snr) tuple """
        return cls.random_points(1).next()

    @classmethod
    def random_data(cls, pfilter = None, cstars = None):
        """ Return the arguments needed to instantiate a random LightCurve.

        Return a four-element tuple: (a) a photometric filter and three lists
        with the (b) comparison stars, (c) weights and (b) standard deviations.
        If specified, the photometric filter and comparison stars IDs are not
        random. Instead, the passed values (a Passband instance and a sequence
        of integers, respectively) are used.

        """

        if not pfilter:
            pfilter = passband.Passband.random()
        if cstars:
            size = len(cstars)
        else:
            size = random.randint(MIN_NSTARS, MAX_NSTARS)
            cstars = [random.randint(cls.MIN_ID, cls.MAX_ID) for x in xrange(size)]
        cweights = Weights.random(size)
        cstdevs = Weights.inversely_proportional(cweights)
        return pfilter, cstars, cweights, cstdevs

    @classmethod
    def random(cls, pfilter = None, cstars = None):
        """ Return a random, empty LightCurve.

        If given, the photometric filter and comparison star IDs are not
        random. Instead, the passed values (a Passband instance and a
        sequence of integers, respectively) are used.

        """

        args = cls.random_data(pfilter = pfilter, cstars = cstars)
        return LightCurve(*args)

    @classmethod
    def populate(cls, curve, images):
        """ Populate a curve with a random point for each of the Images """
        for img in images:
            magnitude, snr = cls.random_point()[1:]
            curve.add(img.unix_time, magnitude, snr)
        return curve

    def test_init(self):
        for _ in xrange(NITERS):
            pfilter, cstars, cweights, cstdevs = args = self.random_data()
            curve = LightCurve(*args)
            self.assertEqual(curve.pfilter, pfilter)
            self.assertEqual(curve.cstars, cstars)
            self.assertTrue(numpy.all(numpy.equal(curve.cweights, cweights)))
            self.assertTrue(numpy.all(numpy.equal(curve.cstdevs, cstdevs)))

        # ValueError raised if number of weights != stdevs != comparison stars
        pfilter, cstars, cweights, cstdevs = self.random_data()
        rindex = random.choice(range(len(cweights)))
        cweights = cweights.rescale(rindex) # remove rindex-th coefficient
        assert len(cstars) != len(cweights) != len(cstdevs)
        with self.assertRaises(ValueError):
            LightCurve(pfilter, cstars, cweights, cstdevs)

        # Same exception also raised if there are no comparison stars
        cstars = cweights = cstdevs = []
        with self.assertRaises(ValueError):
            LightCurve(pfilter, cstars, cweights, cstdevs)

    def test_add_len_and_getitem(self):
        for _ in xrange(NITERS):
            curve = self.random()
            size = random.randint(MIN_NSTARS, MAX_NSTARS)
            for index, point in enumerate(self.random_points(size)):
                curve.add(*point)
                self.assertEqual(len(curve), index + 1)
                self.assertEqual(curve[index], point)

    def test_iter(self):

        # A specific, non-random test case...
        curve = self.random()
        point1 = (25000, 21.1, 115)
        point2 = (23000, 19.2, 125)
        point3 = (24000, 19.9, 110)
        point4 = (19000, 20.3, 105)
        for index in xrange(1, 5):
            curve.add(*eval('point%d' % index))
        points = iter(curve)
        self.assertEqual(points.next(), point4)
        self.assertEqual(points.next(), point2)
        self.assertEqual(points.next(), point3)
        self.assertEqual(points.next(), point1)
        with self.assertRaises(StopIteration):
            points.next()

        # ... and the random test cases
        for _ in xrange(NITERS):
            curve = self.random()
            size = random.randint(MIN_NSTARS, MAX_NSTARS)
            curve_points = list(self.random_points(size))
            for point in curve_points:
                curve.add(*point)

            # __iter__ returns the tuples *sorted chronologically*
            curve_points.sort(key = operator.itemgetter(0))
            self.assertEqual(list(curve), curve_points)

    def test_stdev(self):

        curve = self.random()
        assert not len(curve)
        curve.add(15000, 14.5, 100)
        curve.add(16000, 15.6, 125)
        curve.add(21000, 13.1, 200)
        self.assertAlmostEqual(curve.stdev, 1.0230672835481871)

        for _ in xrange(NITERS):
            curve = self.random()
            size = random.randint(MIN_NSTARS, MAX_NSTARS)
            magnitudes = []
            for point in self.random_points(size):
                unix_time, mag, snr = point
                curve.add(*point)
                magnitudes.append(mag)
            self.assertAlmostEqual(curve.stdev, numpy.std(magnitudes))

        # ValueError is raised if the LightCuve is empty
        curve = self.random()
        assert not len(curve)
        with self.assertRaises(ValueError):
            curve.stdev

    def test_weights(self):
        for _ in xrange(NITERS):
            pfilter, cstars, cweights, cstdevs = self.random_data()
            curve = LightCurve(pfilter, cstars, cweights, cstdevs)
            it = enumerate(curve.weights())
            for index, (cstar_id, cweight, cstdev) in it:
                self.assertEqual(cstar_id, cstars[index])
                self.assertEqual(cweight, cweights[index])
                self.assertEqual(cstdev, cstdevs[index])

    def test_amplitude(self):

        curve = self.random()
        assert not len(curve)
        curve.add(10000, 14.5, 100)
        curve.add(11000, 13.1, 110)
        curve.add(16000, 13.4, 125)
        curve.add(19000, 14.3, 150)
        curve.add(21000, 14.45, 125)

        # 14.5 - 13.1 = 1.4 (only one point is used to determine the peak
        # and trough, so whether we use the median or mean is irrelevant)
        amplitude1 = curve.amplitude(npoints = 1, median = True)
        amplitude2 = curve.amplitude(npoints = 1, median = False)
        self.assertEqual(amplitude1, amplitude2)
        self.assertAlmostEqual(amplitude1, 1.4)

        # median([14.5, 14.45]) - median([13.1, 13.4]) = 1.225
        amplitude = curve.amplitude(npoints = 2, median = True)
        self.assertAlmostEqual(amplitude, 1.225)

        # mean([14.5, 14.45, 14.3]) - mean([13.1, 13.4, 14.3]) = ~0.81666
        amplitude = curve.amplitude(npoints = 3, median = False)
        self.assertAlmostEqual(amplitude, 0.81666666666666643)

        for _ in xrange(NITERS):

            curve = self.random()
            size = random.randint(MIN_NSTARS, MAX_NSTARS)
            magnitudes = []
            for point in self.random_points(size):
                unix_time, mag, snr = point
                curve.add(*point)
                magnitudes.append(mag)

            magnitudes.sort()
            npoints = random.choice(xrange(1, len(curve) + 1))

            # First using the median to combine the points...
            maximum = numpy.median(magnitudes[-npoints:])
            minimum = numpy.median(magnitudes[:npoints])
            amplitude = maximum - minimum

            kwargs = dict(npoints = npoints, median = True)
            self.assertAlmostEqual(amplitude, curve.amplitude(**kwargs))

            # ... and now with the arithmetic mean
            maximum = numpy.mean(magnitudes[-npoints:])
            minimum = numpy.mean(magnitudes[:npoints])
            amplitude = maximum - minimum
            kwargs = dict(npoints = npoints, median = False)
            self.assertAlmostEqual(amplitude, curve.amplitude(**kwargs))

        # ValueError is raised if the LightCuve is empty
        curve = self.random()
        assert not len(curve)
        self.assertRaises(ValueError, curve.amplitude)

    def test_ignore_noisy(self):

        curve = self.random()
        assert not len(curve)
        curve.add(10000, 12.1, 100)
        curve.add(11000, 12.7, 250)
        curve.add(12000, 12.9, 400)
        curve.add(13000, 12.6, 375)
        curve.add(14000, 12.2, 95)

        # Two points (mags 12.1 and 12.2) have a SNR < 250
        curve = curve.ignore_noisy(250)
        self.assertEqual(len(curve), 3)
        utimes, mags, snrs = zip(*curve)
        self.assertEqual(utimes, (11000, 12000, 13000))
        self.assertEqual(mags, (12.7, 12.9, 12.6))
        self.assertEqual(snrs, (250, 400, 375))

        for _ in xrange(NITERS):

            curve = self.random()
            size = random.randint(MIN_NSTARS, MAX_NSTARS)
            [curve.add(*point) for point in self.random_points(size)]
            utimes, mags, snrs = zip(*curve)

            # The random minimum SNR could discard none to all of the points
            threshold = random.uniform(0, max(snrs) + 1)
            non_noisy_curve = curve.ignore_noisy(threshold)

            # The values that we expect the new light curve to contain. If the
            # SNR of all the points is below the threshold, ignore_noisy must
            # return an empty LightCurve object.
            npoints = [p for p in curve if p[-1] >= threshold]
            self.assertEqual(len(non_noisy_curve), len(npoints))

            if npoints:
                nutimes, nmags, nsnrs = zip(*npoints)
                self.assertTrue(min(nsnrs) >= threshold)
                self.assertEqual(nutimes, tuple(p[0] for p in non_noisy_curve))
                self.assertEqual(nmags,   tuple(p[1] for p in non_noisy_curve))
                self.assertEqual(nsnrs,   tuple(p[2] for p in non_noisy_curve))

    @staticmethod
    def assertThatAreEqual(cls, first, second):
        """ Assert that two LightCurves are equal.

        Two LightCurve are instances if they have the same photometric filter
        and number of points, with the same values. Furthermore, the comparison
        stars that they use and the corresponding weights must also be equal.
        Except for the photometric filter and the star IDs, values are
        considered to be equal according to the rules of the
        TestCase.assertAlmostEqual method.

        The first parameter, 'cls', must be the subclass of unittest.TestCase
        instance in which this test is done. And we need to make the method
        static so that it does not receives an implicit first argument.
        This is kind of an ugly hack, yes.

        """

        cls.assertEqual(first.pfilter, second.pfilter)
        cls.assertEqual(len(first), len(second))
        for point1, point2 in zip(first, second):
            # Three-element tuples: time, mag and SNR (floats)
            cls.assertAlmostEqual(point1[0], point2[0])
            cls.assertAlmostEqual(point1[1], point2[1])
            cls.assertAlmostEqual(point1[2], point2[2])

        for weight1, weight2 in zip(sorted(first.weights()),
                                    sorted(first.weights())):
            # Two-element tuples: ID and the weight
            cls.assertEqual(weight1[0], weight2[0])
            cls.assertAlmostEqual(weight1[1], weight2[1])


class LEMONdBTest(unittest.TestCase):

    MIN_NIMAGES = 100  # Minimum number of images for random databases
    MAX_NIMAGES = 250  # Maximum number of images for random databases
    MIN_ID = 1     # Minimum value for random star IDs
    MAX_ID = 9999  # Maximum value for random star IDs
    XSIZE = 2000 # Dimensions in pixels of the imaginary reference image on
    YSIZE = 2000 # which random stars were detected. Thus, the values of their
                 # x- and y-coordinates will be in the range [0, X/YSIZE]
    MIN_MAG = 1.47   # Minimum value for random magnitudes (Saturn)
    MAX_MAG = 25     # Maximum value for random magnitudes (Fenrir)
    MIN_SNR = 2      # Minimum value for random signal-to-noise ratios
    MAX_SNR = 10000  # Maximum value for random signal-to-noise ratios

    OBSERVED_PROB = 0.50  # Probability of a star being observed in each image,
                          # needed some test cases which use random sets of data

    # Minimum and maximum number of random json_parse.CandidateAnnuli objects
    # to be added and retrieved from random databases in the test cases of the
    # LEMONdB.add_candidate_pparams and get_candidate_pparams methods.
    MIN_CANDIDATE_PPARAMS = 1
    MAX_CANDIDATE_PPARAMS = 128

    # Minimum and maximum number of photometrics filters in which a
    # random star may have been observed in some of the test cases
    MIN_NFILTERS = 1
    MAX_NFILTERS = len(passband.Passband.all())

    # Minimum and maximum values (in seconds) for random periods
    MIN_PERIOD = 1
    MAX_PERIOD = 3600 * 24 * 365 # a year

    # Minimum and maximum random steps (used with the string-length method)
    MIN_STEP = 1
    MAX_STEP = 3600 # one hour

    PERIOD_PROB = 0.50  # Probability of a star having a period in a filter;
                        # using in test test_and_get_period_and_get_periods

    @staticmethod
    def random_path():
        """ Return to path to a temporary file, without creating it """
        fd, path = tempfile.mkstemp(suffix = '.LEMONdB')
        os.close(fd)
        os.unlink(path)
        return path

    def test_init(self):
        for _ in xrange(NITERS):
            path = self.random_path()
            try:
                db = LEMONdB(path)
                self.assertEqual(db.path, path)
                self.assertTrue(os.path.exists(path))
            finally:
                os.unlink(path)

    def test_add_and_get_candidate_pparams(self):

        for _ in xrange(NITERS):

            db = LEMONdB(':memory:')
            args = self.MIN_CANDIDATE_PPARAMS, self.MAX_CANDIDATE_PPARAMS
            how_many = random.randint(*args) # number of CandidateAnnuli

            # Map each photometric filter to a list of CandidateAnnuli objects;
            # to keep the track of which have been added for each filter. This
            # is what we we expect LEMONdB.get_candidate_pparams to return.
            added = collections.defaultdict(list)

            for index in xrange(how_many):
                args = list(PhotometricParametersTest.random_data())
                args.append(random.random()) # the standard deviation
                candidate = CandidateAnnuli(*args)
                pfilter = passband.Passband.random()

                added[pfilter].append(candidate)
                db.add_candidate_pparams(candidate, pfilter)

            # Add the last CandidateAnnuli object for the same photometric
            # filter again, but using a different standard deviation. This
            # should replace the record we added previously, and for that
            # reason we need to update 'added' accordingly.
            candidate = candidate._replace(stdev = random.random())

            kwargs = dict(stdev = candidate.stdev)
            added[pfilter][-1] = added[pfilter][-1]._replace(**kwargs)
            db.add_candidate_pparams(candidate, pfilter)

            for pfilter, expected in added.iteritems():
                # The returned CandidateAnnuli are sorted by their stdev
                expected.sort(key = operator.attrgetter('stdev'))
                retrieved = db.get_candidate_pparams(pfilter)
                self.assertEqual(retrieved, expected)

    @staticmethod
    def images_filters_tables_status(db):
        """ Return two sorted tuples with all the information stored in the
        IMAGES and PHOTOMETRIC_FILTERS tables of the 'db' LEMONdB"""

        query_images  = "SELECT * FROM images ORDER BY unix_time"
        query_filters = "SELECT * FROM photometric_filters ORDER BY id"
        db._execute(query_images)
        images = tuple(db._rows)
        db._execute(query_filters)
        filters = tuple(db._rows)
        return images, filters

    def test_simage_and_mosaic(self):
        db = LEMONdB(':memory:')

        self.assertEqual(db.simage, None)
        self.assertEqual(db.mosaic, None)

        for index in xrange(NITERS):
            # We need more than the path to a fictitious FITS file: it must
            # actually exist, since LEMONdB.simage stores it in the database.
            # Use test.test_fitsimage.FITSImage.random() to get a subclass of
            # fitsimage.FITSImage which is automatically deleted on exit from
            # the body of the with statement.
            with test.test_fitsimage.FITSImageTest.random() as input:
                img = ImageTest.random()
                img = img._replace(path = input.path)
                db.simage = img
                self.assertEqual(db.simage, img)
                with test.test_fitsimage.FITSImage(db.mosaic) as output:
                    self.assertEqual(input.sha1sum, output.sha1sum)

        # LEMONdB.simage replaces the existing source image, if any, but it
        # raises DuplicateImageError if an image with the same Unix time and
        # photometric filter is already stored in the database. To test this,
        # add a new image with the same filter and observation date than the
        # last Image added in the above loop.

        kwargs = dict(pfilter = img.pfilter, unix_time = img.unix_time)
        duplicate_img = ImageTest.random()._replace(**kwargs)
        before_tables = self.images_filters_tables_status(db)
        with self.assertRaises(DuplicateImageError):
            db.simage = duplicate_img
        after_tables = self.images_filters_tables_status(db)
        self.assertEqual(before_tables, after_tables)

    def test_add_and_get_image(self):
        db = LEMONdB(':memory:')
        size = random.randint(self.MIN_NIMAGES, self.MAX_NIMAGES)

        # Add the random Images to the database. Use a two-level dictionary to
        # map each photometric filter to a Unix time to an Image, so that we
        # know what we have added and can later be able to fetch them from the
        # LEMONdB and check that they are the same.
        added_images = collections.defaultdict(dict)

        for index, img in enumerate(ImageTest.nrandom(size)):
            added_images[img.pfilter][img.unix_time] = img
            db.add_image(img)

        # Add the lengths of the second-level dictionaries (i.e., find the
        # number of images that have been added to the LEMONdB, regardless
        # of their photometric filter and any duplicate Unix time)
        assert sum(len(p) for p in added_images.itervalues()) == size

        # Extract the Images from the LEMONdB and check that the returned
        # objects are exactly equal to what we added previously
        for pfilter in added_images.iterkeys():
            for unix_time, input_img in added_images[pfilter].iteritems():
                output_img = db.get_image(unix_time, pfilter)
                self.assertEqual(input_img, output_img)

        # LEMONdB.get_image raises KeyError if there is no Image with the
        # specified Unix time *and* photometric filter. Try the three possible
        # combinations: missing Unix time, missing photometric filter, or both.

        # (1) Missing Unix time. To do this we need a set with the Unix times
        # of all the Images that have been added to the LEMONdB, regardless of
        # their photometric filter; in other words, we need to extract all the
        # keys of the second-level dictionary.
        all_utimes = itertools.chain(*[x.keys() for x in added_images.values()])
        nonexistent_unix_time = different_runix_time(all_utimes)
        pfilter = random.choice(added_images.keys())
        with self.assertRaises(KeyError):
            db.get_image(nonexistent_unix_time, pfilter)

        # (2) Missing photometric filter. We need a fresh LEMONdB here as all
        # the photometric filters that the Passband class can encapsulate may
        # have been used when we populated the database in the above code.
        db = LEMONdB(':memory:')
        image = ImageTest.random()
        db.add_image(image)

        nonexistent_filter = image.pfilter.different()
        with self.assertRaises(KeyError):
            db.get_image(image.unix_time, nonexistent_filter)

        # (3) Missing Unix time *and* photometric filter
        nonexistent_unix_time = different_runix_time([image.unix_time])
        with self.assertRaises(KeyError):
            db.get_image(nonexistent_unix_time, nonexistent_filter)

        # LEMONdB.add_image raises DuplicateImageError if we add an image with
        # the same date of observation *and* photometric filter that another
        # image that was previously added. In these cases we need to verify
        # that any possible insertions into the tables are rolled back, so
        # that the database is not modified in case of an error.
        db = LEMONdB(':memory:')
        img1 = ImageTest.random()
        db.add_image(img1)

        # Another random Image with the same Unix time filter.
        img2 = ImageTest.random(pfilter = img1.pfilter)
        img2 = img2._replace(unix_time = img1.unix_time)

        assert img1.unix_time == img2.unix_time
        assert img1.pfilter == img2.pfilter

        before_tables = self.images_filters_tables_status(db)
        with self.assertRaises(DuplicateImageError):
            db.add_image(img2)
        after_tables = self.images_filters_tables_status(db)
        self.assertEqual(before_tables, after_tables)

    def test_add_and_get_image_None_fields(self):

        def img_None_attr(attr):
            """ Return a random Image where only this attribute is None """
            kwargs = {attr : None}
            return ImageTest.random()._replace(**kwargs)

        def assert_None_raises(db, attr, is_sources_img):
            """ Test that, when an Image whose 'attr' attribute is None is
            added to the LEMONdB with add_image(), sqlite3.IntegrityError is
            raised. Make also sure that none of the involved tables of the
            database are modified (i.e., that the transacion is rolled back)"""

            before_tables = self.images_filters_tables_status(db)
            regexp = re.compile("%s may not be NULL" % attr, re.IGNORECASE)
            img = img_None_attr(attr)
            with self.assertRaisesRegexp(sqlite3.IntegrityError, regexp):
                db.add_image(img, _is_sources_img = is_sources_img)
            after_tables = self.images_filters_tables_status(db)
            self.assertEqual(before_tables, after_tables)

        # Except for 'object', LEMONdB.add_image() does not allow any of the
        # attributes of the Image object to be None (and therefore NULL in the
        # IMAGES table). Raises sqlite3.IntegrityError if that is the case (for
        # example: "images.ra may not be NULL")

        db = LEMONdB(':memory:')
        attributes = set(Image._fields)
        for attr in attributes.difference(set(['object'])):
            img = img_None_attr(attr)
            assert_None_raises(db, attr, False)

        # The object name can be None
        input_img = img_None_attr('object')
        db.add_image(input_img)
        output_img = db.get_image(input_img.unix_time, input_img.pfilter)
        self.assertEqual(input_img, output_img)

        # Things are different for the sources image (_is_sources_img = True):
        # in those cases, only the path to the image, its right ascension and
        # declination are mandatory. The other fields (filter, date, airmass
        # and gain) may be meaningless if we assemble several images into a
        # custom mosaic and use the resulting image to detect the sources.

        mandatory_fields = set(['path', 'ra', 'dec'])
        for attr in mandatory_fields:
            img = img_None_attr(attr)
            assert_None_raises(db, attr, True)

        for attr in attributes.difference(mandatory_fields):
            img = img_None_attr(attr)
            db.add_image(img, _is_sources_img = True)
            # We do not test LEMONdB.get_image() here: the sources image should
            # be extracted from the database using the LEMONdB.simage property.

    def test_add_image_ra_dec_out_of_range(self):

        def test_img_addition(ra = None, dec = None):
            """ Insert a random Image object, returned by ImageTest.random(),
            in a fresh LEMONdB. In case the 'ra' or 'dec' (or both) keyword
            arguments are specified, these attributes of the random Image
            object are replaced with their values before the insertion.

            """
            db = LEMONdB(':memory:')
            img = ImageTest.random()

            kwargs = {}
            if ra:
                kwargs['ra'] = ra
            if dec:
                kwargs['dec'] = dec
            if kwargs:
                img = img._replace(**kwargs)
            db.add_image(img)

        def test_img_addition_raises(regexp, ra = None, dec = None):
            """ Make sure that, when test_img_addition(ra = ra, dec = dec)
            is called, the sqlite3.IntegrityError exception is raised and
            its string representation matches 'regexp'

            """
            with self.assertRaisesRegexp(sqlite3.IntegrityError, regexp):
                test_img_addition(ra = ra, dec = dec)

        # RA [0, 360]: boundary cases
        regexp = "RA out of range \[0, 360\["
        test_img_addition_raises(regexp, ra = -5)
        test_img_addition_raises(regexp, ra = -0.001)
        test_img_addition(ra = 0) # boundary (inclusive)
        test_img_addition(ra = 0.0025)

        test_img_addition(ra = 120)
        test_img_addition(ra = 240)

        test_img_addition(ra = 359.9999999)
        test_img_addition_raises(regexp, ra = 360) # boundary (exclusive)
        test_img_addition_raises(regexp, ra = 360.13)
        test_img_addition_raises(regexp, ra = 400)

        # DEC [-90, 90]: boundary cases
        regexp = "DEC out of range \[-90, 90\]"
        test_img_addition_raises(regexp, dec = -120)
        test_img_addition_raises(regexp, dec = -90.001)
        test_img_addition(dec = -90) # boundary (inclusive)
        test_img_addition(dec = -45.25)

        test_img_addition(dec = 0)

        test_img_addition(dec = 75.81)
        test_img_addition(dec = 90) # boundary (inclusive)
        test_img_addition_raises(regexp, dec = 90.0001)
        test_img_addition_raises(regexp, dec = 91)

    @classmethod
    def random_stars_info(cls, size):
        """ Return a generator which steps through 'size' random nine-element
        tuples (the ID, x- and y- coordinates of the star in the sources image,
        the right ascension and declination, the astronomical epoch, the two
        proper motions and the instrumental magnitude, in this very order),
        guaranteed to have different IDs.

        """

        stars_ids = range(cls.MIN_ID, cls.MAX_ID + 1)
        if size > len(stars_ids):
            msg = "cannot ask for more stars than available IDs"
            raise ValueError(msg)

        for star_id in random.sample(stars_ids, size):
            x = random.uniform(0, cls.XSIZE)
            y = random.uniform(0, cls.YSIZE)
            ra, dec, pm_ra, pm_dec = get_random_coords()
            epoch = random.choice([1950, 2000])
            imag = random.uniform(cls.MIN_MAG, cls.MAX_MAG)
            yield star_id, x, y, ra, dec, epoch, pm_ra, pm_dec, imag

    @classmethod
    def random_star_info(cls, id_ = None):
        """ Return a nine-element tuple with the information of a random star.

        The method returns a tuple with the ID, x- and y-coordinates, right
        ascension, declination, epoch, both proper motions and instrumental
        magnitude of a random star. The ID is random unless a value for it
        is specified with the 'id_' keyword argument.

        """

        star_info = list(cls.random_stars_info(1).next())
        if id_ is not None:
            star_info[0] = id_
        return star_info

    def test_add_and_get_star(self):
        db = LEMONdB(':memory:')
        size = random.randint(MIN_NSTARS, MAX_NSTARS)
        stars = {} # map each ID to the rest of the star info
        for star_info in self.random_stars_info(size):
            id_ = star_info[0]
            stars[id_] = star_info[1:]
            db.add_star(*star_info)
        self.assertEqual(len(stars), size)

        items = stars.items()
        random.shuffle(items)
        for id_, input_info in items:
            output_info = db.get_star(id_)
            self.assertEqual(input_info, output_info)

        # DuplicateStarError must be raised if the specified ID was
        # already used for another star in the database.
        star_info = list(self.random_star_info())
        star_info[0] = random.choice(stars.keys()) # reuse an ID
        with self.assertRaises(DuplicateStarError):
            db.add_star(*star_info)

        # KeyError must be raised if there is no star for the given
        # ID. We need to generate one not already in the database
        while True:
            id_ = random.randint(self.MIN_ID, self.MAX_ID)
            if id_ not in stars.iterkeys():
                with self.assertRaises(KeyError):
                    db.get_star(id_)
                break

    def test_len(self):
        db = LEMONdB(':memory:')
        size = random.randint(MIN_NSTARS, MAX_NSTARS)
        for index, star_info in enumerate(self.random_stars_info(size)):
            self.assertEqual(len(db), index)
            db.add_star(*star_info)
            self.assertEqual(len(db), index + 1)

    def test_star_ids(self):
        db = LEMONdB(':memory:')
        size = random.randint(MIN_NSTARS, MAX_NSTARS)
        stars_ids = [] # keep track of inserted IDs
        self.assertEqual([], db.star_ids)
        for star_info in self.random_stars_info(size):
            id_ = star_info[0]
            stars_ids.append(id_)
            db.add_star(*star_info)
            self.assertEqual(sorted(stars_ids), db.star_ids)

    @classmethod
    def random_stars(cls, size, unix_times):
        """ Return a generator which steps through 'size' random DBstars.

        The returned stars are observed in the dates given in the sequence
        'unix_times': a random number of unique elements will be chosen among
        these Unix times for each star. In other words: each of the random
        stars may be observed, potentially, in all the Unix times. The
        returned stars are guaranteed to have different IDs.

        """

        stars_ids = range(cls.MIN_ID, cls.MAX_ID + 1)
        if size > len(stars_ids):
            msg = "cannot ask for more stars than available IDs"
            raise ValueError(msg)

        for star_id in random.sample(stars_ids, size):
            pfilter = passband.Passband.random()
            rows = []
            for unix_time in unix_times:
                if random.random() < cls.OBSERVED_PROB:
                    magnitude = random.uniform(cls.MIN_MAG, cls.MAX_MAG)
                    snr = random.uniform(cls.MIN_SNR, cls.MAX_SNR)
                    rows.append((unix_time, magnitude, snr))
            yield DBStar.make_star(star_id, pfilter, rows)

    def test_add_and_get_pm_correction(self):

        # A specific, non-random test case. Three astronomical objects,
        # the first two with known proper-motions, the third without it.

        db = LEMONdB(':memory:')

        star1 = self.random_star_info(id_ = 1)
        star1[6:8] = [0.57095399531758917, -9.0025061305781175]
        db.add_star(*star1)

        star2 = self.random_star_info(id_ = 2)
        star2[6:8] = [-0.016290290457260315, -4.9776521868125601]
        db.add_star(*star2)

        star3 = self.random_star_info(id_ = 3)
        star3[6:8] = [None, None]
        db.add_star(*star3)

        pfilter1 = passband.Passband("Johnson V")
        pfilter2 = passband.Passband("Johnson R")

        img1, img2 = ImageTest.nrandom(2)
        img1 = img1._replace(pfilter = pfilter1)
        img2 = img2._replace(pfilter = pfilter2)

        db.add_image(img1)
        db.add_image(img2)

        utime1 = img1.unix_time
        utime2 = img2.unix_time

        f = db.add_pm_correction
        f(1, utime1, pfilter1,  28.87617936271731, -36.84344057247144)
        f(1, utime2, pfilter2, -34.18950418859132,  10.93687152558971)
        f(2, utime1, pfilter1,  -4.65494684748566, -31.26563816482958)
        # Do not store a proper-motion correction for the second image

        x, y = db.get_pm_correction(1, utime1, pfilter1)
        self.assertAlmostEqual(x,  28.87617936271731)
        self.assertAlmostEqual(y, -36.84344057247144)

        x, y = db.get_pm_correction(1, utime2, pfilter2)
        self.assertAlmostEqual(x, -34.18950418859132)
        self.assertAlmostEqual(y,  10.93687152558971)

        x, y = db.get_pm_correction(2, utime1, pfilter1)
        self.assertAlmostEqual(x,  -4.65494684748566)
        self.assertAlmostEqual(y, -31.26563816482958)

        x, y = db.get_pm_correction(2, utime2, pfilter2)
        self.assertIs(x, None)
        self.assertIs(y, None)

        # If no proper motion was stored for the astronomical object with
        # LEMONdB.add_star(), ValueError is raised if we now try to store
        # a proper-motion correction for it.
        regexp = re.compile("has no proper motion", re.IGNORECASE)
        with self.assertRaises(ValueError):
            db.add_pm_correction(3, utime1, pfilter1, -42.878132, -32.923991)
        with self.assertRaises(ValueError):
            db.add_pm_correction(3, utime2, pfilter2, 9.9565731, -21.9961378)

        # UnknownStarError if 'star_id' does not match the ID of any of the
        # stars in the database
        star_id = 4
        star4 = self.random_star_info(id_ = star_id)
        star4[6:8] = [-8.6902812929289581, -41.295196884794784]
        self.assertTrue(star_id not in db.star_ids)

        args = star_id, utime1, pfilter1, 4.88648990615375, -3.50462606269579
        with self.assertRaises(UnknownStarError):
            db.add_pm_correction(*args)

        # But if the star is added before, it works, of course
        db.add_star(*star4)
        db.add_pm_correction(*args)
        x, y = db.get_pm_correction(star_id, utime1, pfilter1)
        self.assertAlmostEqual(x, args[-2])
        self.assertAlmostEqual(y, args[-1])

        # As with LEMONdB.add_photometry(), the UnknownImageError exception is
        # raised if the image with the specified Unix time *and* photometric
        # filter does not exist. Try the three possible combinations: missing
        # Unix time, missing photometric filter, or both.

        # (1) Missing Unix time
        nonexistent_unix_time = different_runix_time([utime1, utime2])
        args = star_id, nonexistent_unix_time, pfilter1, 29.94581, -74.20631
        with self.assertRaises(UnknownImageError):
            db.add_pm_correction(*args)

        # (2) Missing photometric filter
        nonexistent_filter = passband.Passband("Johnson I")
        self.assertTrue(nonexistent_filter not in db.pfilters)
        args = star_id, utime1, nonexistent_filter, -42.14321, 58.89305
        with self.assertRaises(UnknownImageError):
            db.add_pm_correction(*args)

        # (3) Missing Unix time *and* photometric filter
        args = star_id, nonexistent_unix_time, nonexistent_filter, -57.188, -2.267
        with self.assertRaises(UnknownImageError):
            db.add_pm_correction(*args)

        # sqlite3.IntegrityError is raised if we attempt to add a second record
        # for the same astronomical object, Unix time and photometric filter
        # (that is, if more than one proper-motion correction is added for the
        # same object and image)
        args = [star_id, utime1, pfilter1, -25.592231, -21.32372]
        regexp = re.compile("columns star_id, image_id are not unique")
        with self.assertRaisesRegexp(sqlite3.IntegrityError, regexp):
            db.add_pm_correction(*args)

        # LEMONdB.get_pm_correction() raises KeyError if the ID does not match
        # that of any of the astronomical objects already in the database.
        nonexistent_id = 5
        self.assertTrue(nonexistent_id not in db.star_ids)
        with self.assertRaises(KeyError):
            db.get_pm_correction(nonexistent_id, utime1, pfilter1)

    def test_add_and_get_photometry(self):

        # A specific, non-random test case...
        db = LEMONdB(':memory:')
        johnson_B = passband.Passband('B')
        johnson_V = passband.Passband('V')

        star_ids = [star1_id, star2_id] = range(2)
        for id_ in star_ids:
            star_info =  self.random_star_info(id_ = id_)
            db.add_star(*star_info)

        # Note that the second image precedes the first one in time, and that
        # the third and fourth have the same Unix time. This is fine as long
        # as they have different photometric filters, as the LEMONdB class
        # supports simultaneous observations.

        img1 = ImageTest.random(johnson_B)
        img1 = img1._replace(unix_time = 100000)
        img2 = ImageTest.random(johnson_B)
        img2 = img2._replace(unix_time = 90000)
        img3 = ImageTest.random(johnson_V)
        img3 = img3._replace(unix_time = 150400)
        img4 = ImageTest.random(johnson_B)
        img4 = img4._replace(unix_time = img3.unix_time)

        images = [img1, img2, img3, img4]
        for img in images:
            db.add_image(img)

        db.add_photometry(star1_id, img1.unix_time, img1.pfilter, 10.1, 150)
        db.add_photometry(star1_id, img2.unix_time, img2.pfilter, 10.8, 125)
        db.add_photometry(star1_id, img3.unix_time, img3.pfilter, 9.5, 175)
        db.add_photometry(star1_id, img4.unix_time, img4.pfilter, 7.3, 200)

        db.add_photometry(star2_id, img1.unix_time, img1.pfilter, 7.1, 350)
        db.add_photometry(star2_id, img2.unix_time, img2.pfilter, 7.9, 315)
        db.add_photometry(star2_id, img3.unix_time, img3.pfilter, 5.8, 550)
        db.add_photometry(star2_id, img4.unix_time, img4.pfilter, 5.9, 450)

        # The records in DBStars are sorted by their Unix time; therefore, the
        # magnitudes of the stars in the second Image should come before those
        # in the first one, as Unix time 90,000 < 100,000.
        star1_B = db.get_photometry(star1_id, johnson_B)
        self.assertEqual(len(star1_B), 3)

        self.assertEqual(star1_B.time(0), img2.unix_time)
        self.assertEqual(star1_B.mag(0), 10.8)
        self.assertEqual(star1_B.snr(0), 125)

        self.assertEqual(star1_B.time(1), img1.unix_time)
        self.assertEqual(star1_B.mag(1), 10.1)
        self.assertEqual(star1_B.snr(1), 150)

        self.assertEqual(star1_B.time(2), img4.unix_time)
        self.assertEqual(star1_B.mag(2), 7.3)
        self.assertEqual(star1_B.snr(2), 200)

        # The star ID can have a data type other than Python's built-in int
        star1_V = db.get_photometry(numpy.int16(star1_id), johnson_V)
        self.assertEqual(len(star1_V), 1)

        self.assertEqual(star1_V.time(0), img3.unix_time)
        self.assertEqual(star1_V.mag(0), 9.5)
        self.assertEqual(star1_V.snr(0), 175)

        star2_B = db.get_photometry(numpy.int32(star2_id), johnson_B)
        self.assertEqual(len(star2_B), 3)

        self.assertEqual(star2_B.time(0), img2.unix_time)
        self.assertEqual(star2_B.mag(0), 7.9)
        self.assertEqual(star2_B.snr(0), 315)

        self.assertEqual(star2_B.time(1), img1.unix_time)
        self.assertEqual(star2_B.mag(1), 7.1)
        self.assertEqual(star2_B.snr(1), 350)

        self.assertEqual(star2_B.time(2), img4.unix_time)
        self.assertEqual(star2_B.mag(2), 5.9)
        self.assertEqual(star2_B.snr(2), 450)

        star2_V = db.get_photometry(numpy.int64(star2_id), johnson_V)
        self.assertEqual(star2_V.time(0), img3.unix_time)
        self.assertEqual(star2_V.mag(0), 5.8)
        self.assertEqual(star2_V.snr(0), 550)

        # LEMONdB.add_photometry raises UnknownStarError if 'star_id' does not
        # match the ID of any of the stars in the database
        star_info = self.random_star_info()
        star_id = star_info[0] = 2
        self.assertTrue(star_id not in db.star_ids)

        img = ImageTest.random(pfilter = johnson_B)
        db.add_image(img)

        args = [star_id, img.unix_time, img.pfilter, 14.5, 100]
        with self.assertRaises(UnknownStarError):
            db.add_photometry(*args)
        # But if the star is added before, it works, of course
        db.add_star(*star_info)
        db.add_photometry(*args)

        # The UnknownImageError exception is raised by LEMONdB.add_photometry
        # if there is no Image with the specified Unix time *and* photometric
        # filter. Try the three possible combinations: missing Unix time,
        # missing photometric filter, or both.

        # (1) Missing Unix time
        nonexistent_unix_time = different_runix_time([img.unix_time])
        args = star_id, nonexistent_unix_time, img.pfilter, 15.4, 125
        with self.assertRaises(UnknownImageError):
            db.add_photometry(*args)

        # (2) Missing photometric filter
        nonexistent_filter = passband.Passband('Z')
        self.assertTrue(nonexistent_filter not in db.pfilters)
        args = star_id, img.unix_time, nonexistent_filter, 16.5, 175
        with self.assertRaises(UnknownImageError):
            db.add_photometry(*args)

        # (3) Missing Unix time *and* photometric filter
        args = star_id, nonexistent_unix_time, nonexistent_filter, 17.8, 235
        with self.assertRaises(UnknownImageError):
            db.add_photometry(*args)

        # DuplicatePhotometryError is raised if we attempt to add a second
        # record for the same star, Unix time and photometric filter (that is,
        # if more than one photometric record is added for the same Image).
        args = [star_id, img.unix_time, img.pfilter, 13.2, 125]
        with self.assertRaises(DuplicatePhotometryError):
            db.add_photometry(*args)

        # LEMONdB.get_photometry raises KeyError if the star ID does not
        # match that of any of the stars already stored in the database
        nonexistent_id = 3
        self.assertTrue(nonexistent_id not in db.star_ids)
        with self.assertRaises(KeyError):
            db.get_photometry(nonexistent_id, johnson_V)

        # If the star has no records in a filter, an empty DBStar is returned
        empty_star = db.get_photometry(star_id, johnson_V)
        self.assertEqual(len(empty_star), 0)

    def test_pfilters_and_star_pfilters(self):

        db = LEMONdB(':memory:')

        # FITS file is deleted on exit from the with statement.
        with test.test_fitsimage.FITSImageTest.random() as fits:
            img = ImageTest.random()
            img = img._replace(path = fits.path)
            db.simage = img
            self.assertEqual([], db.pfilters)

        nstars = random.randint(MIN_NSTARS, MAX_NSTARS)
        # Add the information of the stars in the database.
        # Use a list to store the IDs of the stars that are added.
        stars_ids = []
        for star_info in self.random_stars_info(nstars):
            star_id = star_info[0]
            stars_ids.append(star_id)
            db.add_star(*star_info)

        # Add some random images to the database. For each one of them, the
        # probability of each star being observed is OBSERVED_PROB. We want to
        # see how the list of filters in which each star (and also the database
        # as a whole) has been observed is correctly updated as images with
        # different filters are added.

        nimages = random.randint(self.MIN_NIMAGES, self.MAX_NIMAGES)
        images = list(ImageTest.nrandom(nimages))

        for img in images:
            db.add_image(img)

        # No photometric record has yet been added
        self.assertEqual([], db.pfilters)

        # Map each star ID to the filters in which it was observed
        stars_pfilters = collections.defaultdict(set)

        seen_pfilters = set()
        for star_id in stars_ids:
            for img in images:
                if random.random() < self.OBSERVED_PROB:
                    mag = random.uniform(self.MIN_MAG, self.MAX_MAG)
                    snr = random.uniform(self.MIN_SNR, self.MAX_SNR)
                    db.add_photometry(star_id, img.unix_time, img.pfilter, mag, snr)
                    stars_pfilters[star_id].add(img.pfilter)
                    seen_pfilters.add(img.pfilter)

                    star_pfilters = db._star_pfilters(star_id)
                    # The returned Passband objects must be already sorted
                    self.assertEqual(star_pfilters, sorted(stars_pfilters[star_id]))
                    self.assertEqual(db.pfilters, sorted(seen_pfilters))

        # LEMONdB._star_pfilters raises KeyError if the star is not in the
        # database. To test this we need to generate a random, new star ID
        while True:
            star_id = random.randint(self.MIN_ID, self.MAX_ID)
            if star_id not in stars_ids:
                with self.assertRaises(KeyError):
                    db._star_pfilters(star_id)
                break

    def test_add_and_get_light_curve(self):

        db = LEMONdB(':memory:')
        curve = LightCurveTest.random()
        nstars = random.randint(MIN_NSTARS, MAX_NSTARS)
        for star_info in LEMONdBTest.random_stars_info(nstars):
            db.add_star(*star_info)

        # Add random images to the database, and use a dictionary to map
        # each photometric filter to the Images that in it were observed.
        images = collections.defaultdict(list)
        size = random.randint(self.MIN_NIMAGES, self.MAX_NIMAGES)
        for img in ImageTest.nrandom(size):
            images[img.pfilter].append(img)
            db.add_image(img)

        # Two level dictionary: map each photometric filter to a star ID to the
        # light curve which was added to the database. Needed to check that the
        # curves read from the database are exactly equal to the input ones.
        light_curves = collections.defaultdict(dict)
        for pfilter in images.iterkeys():
            for star_id in db.star_ids:

                # Select a random number of the other stars as comparison
                candidate_cstars = set(db.star_ids) - set([star_id])
                assert star_id not in candidate_cstars
                ncstars = random.randint(1, len(candidate_cstars))
                cstars = random.sample(candidate_cstars, ncstars)
                curve = LightCurveTest.random(pfilter = pfilter, cstars = cstars)
                curve = LightCurveTest.populate(curve, images[pfilter])
                assert len(curve) == len(images[pfilter])

                light_curves[pfilter][star_id] = curve
                db.add_light_curve(star_id, curve)

        # Now read all the curves, and check that they are equal
        images_keys = images.keys()
        random.shuffle(images_keys)
        for pfilter in images_keys:

            star_ids = db.star_ids
            random.shuffle(star_ids)
            for star_id in star_ids:
                icurve = light_curves[pfilter][star_id]
                ocurve = db.get_light_curve(star_id, pfilter)
                LightCurveTest.assertThatAreEqual(self, icurve, ocurve)

        # Finally, check that the proper exceptions are raised when
        # needed. Again, we need to ensure that when an error occurs
        # none of the tables of the database is modified.

        def tables_status(db):
            """ Return three sorted tuples with all the information stored in
            the LIGHT_CURVES, CMP_STARS and PHOTOMETRIC_FILTERS tables of the
            'db' LEMONdB """

            query_curves = "SELECT * FROM light_curves ORDER BY id"
            query_cmp_stars = "SELECT * FROM cmp_stars ORDER BY id"
            query_filters = "SELECT * FROM photometric_filters ORDER BY id"

            db._execute(query_curves)
            curves = tuple(db._rows)
            db._execute(query_cmp_stars)
            cmp_stars = tuple(db._rows)
            db._execute(query_filters)
            filters = tuple(db._rows)
            return curves, cmp_stars, filters

        def add_test_star(db, *ids):
            """ Add to 'db' LEMONdB random information for the stars whose IDs
            are given as arguments """
            for id_ in ids:
                db.add_star(*LEMONdBTest.random_star_info(id_ = id_))

        def add_test_img(db, *images):
            """ Add to 'db' LEMONdB the Image instances given as arguments """
            for img in images:
                db.add_image(img)


        db = LEMONdB(':memory:')
        star1_id, star2_id, star3_id = range(3)
        pfilter = passband.Passband.random()
        # Three random Images in the same filter, for use only two
        images = list(ImageTest.nrandom(3, pfilter = pfilter))
        img1, img2, img3 = images
        add_test_img(db, img1, img2)
        kwargs = dict(pfilter = pfilter, cstars = [star2_id, star3_id])
        curve = LightCurveTest.random(**kwargs)
        curve = LightCurveTest.populate(curve, [img1, img2])

        # Raises UnknownStarError ('star1_id' not in database)
        args = star1_id, curve
        before_tables = tables_status(db)
        self.assertTrue(star1_id not in db.star_ids)
        with self.assertRaises(UnknownStarError):
            db.add_light_curve(*args)
        self.assertEqual(tables_status(db), before_tables)

        # Let's add the star and one of the comparison stars. The method must
        # fail, as one of the comparison stars, 'star2_id', is still missing.
        # Once we add it, the light curve can finally be added to the LEMONdB.
        add_test_star(db, star1_id, star2_id)
        before_tables = tables_status(db)
        with self.assertRaises(UnknownStarError):
            db.add_light_curve(*args)
        self.assertEqual(tables_status(db), before_tables)

        add_test_star(db, star3_id)
        db.add_light_curve(*args) # works!

        # UnknownImageError raised if the stars are in the LEMONdB but not one
        # (or more) of the images to which the light curve refers; we use the
        # light curve of the second star, which has the third as comparison but
        # records for the three images (although only the first two have been
        # added to the LEMONdB).

        kwargs = dict(pfilter = pfilter, cstars = [star3_id])
        curve = LightCurveTest.random(**kwargs)
        curve = LightCurveTest.populate(curve, [img1, img2, img3])

        args = star2_id, curve
        before_tables = tables_status(db)
        with self.assertRaises(UnknownImageError):
            db.add_light_curve(*args)
        self.assertEqual(tables_status(db), before_tables)

        add_test_img(db, img3)
        db.add_light_curve(*args) # works!

        # DuplicateLightCurvePointError raised if the light curve has more than
        # one point for the same Unix time, or if more than one light curve per
        # star and filter is attempted to be added.

        # Attempt to insert it again... (and fail)
        before_tables = tables_status(db)
        with self.assertRaises(DuplicateLightCurvePointError):
            db.add_light_curve(*args)
        self.assertEqual(tables_status(db), before_tables)

        # Generate the light curve of the third star (uses both the first and
        # second as comparison), but insert two points with the same Unix time
        kwargs = dict(pfilter = pfilter, cstars = [star1_id, star2_id])
        curve = LightCurveTest.random(**kwargs)
        curve = LightCurveTest.populate(curve, [img1, img2, img3])

        # Reuse the Unix time of the last point
        point = list(LightCurveTest.random_point())
        point[0] = curve[-1][0]
        duplicate_curve = copy.deepcopy(curve)
        duplicate_curve.add(*point)

        before_tables = tables_status(db)
        with self.assertRaises(DuplicateLightCurvePointError):
            db.add_light_curve(star3_id, duplicate_curve)
        self.assertEqual(tables_status(db), before_tables)

        db.add_light_curve(star3_id, curve) # the good one works, of course!

        # ValueError if the star uses itself as one of its comparison stars; to
        # test this, we add to the LEMONdB a new star (the fourth one), and use
        # a light curve where it has itself as one of the comparison stars.
        star4_id = 3
        add_test_star(db, star_id, star4_id)
        kwargs = dict(pfilter = pfilter, cstars = [star2_id, star4_id])
        curve = LightCurveTest.random(**kwargs)
        curve = LightCurveTest.populate(curve, [img1, img2, img3])

        before_tables = tables_status(db)
        with self.assertRaises(ValueError):
            db.add_light_curve(star4_id, curve)
        self.assertEqual(tables_status(db), before_tables)

        # The LEMONdB.get_light_curve method, on the other hand, raises
        # KeyError if the star with the given ID is not in the database...
        nstar_id = 17 # any other new value would work
        self.assertTrue(nstar_id not in db.star_ids)
        with self.assertRaises(KeyError):
            db.get_light_curve(nstar_id, pfilter)

        # ... and returns None if the star exists but there is no curve
        add_test_star(db, nstar_id)
        self.assertEqual(None, db.get_light_curve(nstar_id, pfilter))

        # Finally, let's corrupt the database manually: insert a point for the
        # light curve of the star (using the private method _add_curve_point)
        # and leave the database without information on the comparison stars.
        # After this, LEMONdB.get_light_curve must raise sqlite3.IntegrityError
        args = nstar_id, img1.unix_time, img1.pfilter, 7.4, 100
        db._add_curve_point(*args)
        with self.assertRaises(sqlite3.IntegrityError):
            db.get_light_curve(nstar_id, pfilter)

    def test_add_and_get_period_and_get_periods(self):
        db = LEMONdB(':memory:')
        nstars = random.randint(MIN_NSTARS, MAX_NSTARS)
        for star_info in LEMONdBTest.random_stars_info(nstars):
            db.add_star(*star_info)

        # Add random periods to the database, with the probability of each star
        # having had its period calculated being PERIOD_PROB. Use a two-level
        # dictionary to map each star ID to a photometric filter to a
        # two-element tuple with the period and string-length step that was
        # used. None is used, instead of a tuple, when the period is missing.

        periods = collections.defaultdict(dict)
        nfilters = random.randint(self.MIN_NFILTERS, self.MAX_NFILTERS)
        spfilters = random.sample(passband.Passband.all(), nfilters)

        for pfilter in spfilters:
            for star_id in db.star_ids:
                if random.random() < self.PERIOD_PROB:
                    period = random.uniform(self.MIN_PERIOD, self.MAX_PERIOD)
                    step = random.uniform(self.MIN_STEP, self.MAX_STEP)
                    args = (period, step)
                    db.add_period(star_id, pfilter, *args)
                else:
                    args = None
                periods[star_id][pfilter] = args

        # Now read all the periods, and check that they are equal
        star_ids = db.star_ids
        random.shuffle(star_ids)
        random.shuffle(spfilters)
        for star_id in star_ids:
            for pfilter in spfilters:
                istar_period = periods[star_id][pfilter]
                ostar_period = db.get_period(star_id, pfilter)
                if ostar_period is None:
                    self.assertEqual(istar_period, ostar_period)
                else:
                    iperiod, istep = istar_period
                    operiod, ostep = ostar_period
                    self.assertAlmostEqual(iperiod, operiod)
                    self.assertAlmostEqual(istep, ostep)

            # The LEMONdB.get_periods method returns a NumPy
            # array with all the periods of the star.
            expected_periods = []
            for pfilter in periods[star_id].iterkeys():
                period = periods[star_id][pfilter]
                if period is not None:
                    expected_periods.append(period[0])

            returned_periods = list(db.get_periods(star_id))

            expected_periods.sort()
            returned_periods.sort()

            self.assertEqual(len(expected_periods), len(returned_periods))
            if len(expected_periods):
                for iperiod, operiod in zip(expected_periods, returned_periods):
                    self.assertAlmostEqual(iperiod, operiod)

        # LEMONdB.get_periods returns an empty array if there are no periods
        db = LEMONdB(':memory:')
        star_id = random.randint(self.MIN_ID, self.MAX_ID)
        db.add_star(*LEMONdBTest.random_star_info(id_ = star_id))
        self.assertEqual([], list(db.get_periods(star_id)))

        # Finally, check that the proper exceptions are raised when
        # needed. Again, we need to ensure that when an error occurs
        # none of the tables of the database is modified.

        # LEMONdB.add_period raises UnknownStarError if the star ID does
        # not match that of any of the stars already stored in the database
        db = LEMONdB(':memory:')
        star_id = 1
        pfilter = passband.Passband.random()
        period, step = 1400, 15

        def tables_status(db):
            """ Return three sorted tuples with all the information stored in
            the PHOTOMETRIC_FILTERS and PERIODS tables of the 'db' LEMONdB """

            query_filters = "SELECT * FROM photometric_filters ORDER BY id"
            query_periods = "SELECT * FROM periods ORDER BY id"
            db._execute(query_filters)
            filters = tuple(db._rows)
            db._execute(query_periods)
            periods = tuple(db._rows)
            return filters, periods

        args = star_id, pfilter, period, step
        before_tables = tables_status(db)
        with self.assertRaises(UnknownStarError):
            db.add_period(*args)
        self.assertEqual(tables_status(db), before_tables)
        # If the star is added, the ID is no longer unknown
        db.add_star(*LEMONdBTest.random_star_info(id_ = star_id))
        db.add_period(*args) # it works!

        # DuplicatePeriodError is raised if the same is added twice
        before_tables = tables_status(db)
        with self.assertRaises(DuplicatePeriodError):
            db.add_period(*args)
        self.assertEqual(tables_status(db), before_tables)

        # A period is considered duplicate even if we use a different value;
        # which must be unique is the (star ID, photometric filter) pair
        before_tables = tables_status(db)
        with self.assertRaises(DuplicatePeriodError):
            db.add_period(star_id, pfilter, 3755, 13)
        self.assertEqual(tables_status(db), before_tables)

        # Both LEMONdB.get_period and LEMONdB.get_periods must raise KeyError
        # id the star with the specified ID is not stored in the database
        star2_id = 89
        assert star2_id not in db.star_ids
        with self.assertRaises(KeyError):
            db.get_period(star2_id, pfilter)
        with self.assertRaises(KeyError):
            db.get_periods(star2_id)

    def test_airmasses(self):
        db = LEMONdB(':memory:')

        # Two-level dictionary: map each photometric filter to a dictionary
        # which in turns maps each Unix time to its airmass (we're building
        # here the dictionary LEMONdB.airmasses is expected to return)
        airmasses = collections.defaultdict(dict)

        # Make sure that the sources image is not ignored. This FITS file is
        # automatically deleted on exit from the with statement, once it has
        # been stored in the database.
        with test.test_fitsimage.FITSImageTest.random() as fits:
            simg = ImageTest.random()
            simg = simg._replace(path = fits.path)
            db.simage = simg
            airmasses[simg.pfilter][simg.unix_time] = simg.airmass

        size = random.randint(self.MIN_NIMAGES, self.MAX_NIMAGES)
        for img in ImageTest.nrandom(size):
            airmasses[img.pfilter][img.unix_time] = img.airmass
            db.add_image(img)

        for pfilter, airmasses in airmasses.iteritems():
            self.assertEqual(db.airmasses(pfilter), airmasses)

    def test_get_phase_diagram(self):

        db = LEMONdB(':memory:')
        star_id = 1
        pfilter = passband.Passband.random()
        cstars = [2, 3]
        cweights = Weights.random(2)
        curve = LightCurve(pfilter, cstars, cweights)

        all_ids = [star_id] + cstars
        [db.add_star(*LEMONdBTest.random_star_info(id_ = x)) for x in all_ids]
        images = list(ImageTest.nrandom(3, pfilter = pfilter))
        # Note that img2 precedes img1 in time
        img1, img2, img3 = images
        img1 = img1._replace(unix_time = 500)
        img2 = img2._replace(unix_time = 100)
        img3 = img3._replace(unix_time = 785)

        db.add_image(img1)
        db.add_image(img2)
        db.add_image(img3)

        point1 = (img1.unix_time, 14.5, 100)
        point2 = (img2.unix_time, 10.1, 200)
        point3 = (img3.unix_time, 12.2, 140)

        curve.add(*point1)
        curve.add(*point2)
        curve.add(*point3)

        db.add_light_curve(star_id, curve)

        phase = sorted(db.get_phase_diagram(star_id, pfilter, 55))

        epoint1 = (0.272727272727273, 14.5, 100)
        epoint2 = (0, 10.1, 200)
        epoint3 = (0.454545454545455, 12.2, 140)
        expected_phase = sorted((epoint1, epoint2, epoint3))

        assertSequenceOfTuplesAlmostEqual(self, phase, expected_phase)

        # Try now with another value and a different repeat value
        phase = sorted(db.get_phase_diagram(star_id, pfilter, 175, repeat = 2))

        epoint1 = (0.2857142857142856, 14.5, 100)
        epoint2 = (0, 10.1, 200)
        epoint3 = (0.9142857142857141, 12.2, 140)
        epoint4 = (epoint1[0] + 1, 14.5, 100)
        epoint5 = (epoint2[0] + 1, 10.1, 200)
        epoint6 = (epoint3[0] + 1, 12.2, 140)
        expected_phase = [locals()['epoint%d' % x] for x in xrange(1, 7)]
        expected_phase.sort()

        assertSequenceOfTuplesAlmostEqual(self, phase, expected_phase)

        # None returned if the star has no light curve in the filter
        phase = db.get_phase_diagram(star_id, pfilter.different(), 125)
        self.assertEqual(None, phase)

        # KeyError must be raised if the star ID does not match
        # that of any of the stars already stored in the database
        nstar_id = 13
        assert nstar_id not in db.star_ids
        with self.assertRaises(KeyError):
            db.get_phase_diagram(nstar_id, pfilter, 65)

    def test_most_similar_magnitude(self):

        # The ID, x- and y-image coordinates, right ascension, declination,
        # epoch, proper motions (where None means that they are not known) and
        # instrumental magnitudes of the five stars that will be used in this
        # test case. These values are all required by LEMONdB.add_star(), but
        # we are only interested in the first (star ID) and last (magnitude).

        stars_info = \
        [[1, 115, 235, 145.3, 190.4, 2000,   None,   None, 11.5], # star 1
         [2, 931, 833, 146.2, 191.3, 2000,   None,   None, 12.9], # star 2
         [3, 342, 781, 144.1, 190.9, 2000, -0.021, 0.0013, 10.1], # star 3
         [4, 782, 101, 143.9, 189.2, 2000,   None,   None, 13.4], # star 4
         [5, 689, 288, 144.2, 189.4, 2000,   None,   None, 10.3]] # star 5

        db = LEMONdB(':memory:')
        for info in stars_info:
            db.add_star(*info)

        # These stars, which we have just added to the database, will have
        # light curves in the Johnson V and I photometric filters. We will use
        # ten images for the first filter (V) and five for the second (I) but,
        # of course, any number greater than zero would work too. These images
        # need to be generated all at once by ImageTest.nrandom using the same
        # filter (Johnson V), in order to guarantee that their Unix times are
        # unique, and then we change the filter of the last five to Johnson I.

        johnson_V = passband.Passband('V')
        johnson_I = passband.Passband('I')

        # Map each photometric filter to the list of Images observed in it
        images = collections.defaultdict(list)
        for index, img in enumerate(ImageTest.nrandom(15, pfilter = johnson_V)):
            if index >= 10:
                img = img._replace(pfilter = johnson_I)
            images[img.pfilter].append(img)

        assert len(images[johnson_V]) == 10
        assert len(images[johnson_I]) == 5

        # Add the images, for both filters, to the LEMONdB
        for img in itertools.chain(*images.values()):
            db.add_image(img)

        def add_rcurve(star_id, all_ids, pfilter):
            """ Add a random light curve to the LEMONdB for the star with ID
            'star_id' and the 'pfilter' photometric filter. Use as comparison
            stars those whose IDs are in 'all_ids', excluding 'star_id' in case
            it is contained there -- as a star cannot use itself as comparison.
            """

            cstars = list(set(all_ids) - set([star_id]))
            kwargs = dict(pfilter = pfilter, cstars = cstars)
            curve = LightCurveTest.random(**kwargs)
            curve = LightCurveTest.populate(curve, images[pfilter])
            db.add_light_curve(star_id, curve)

        # The first four stars (IDs 1-4) have a light curve in Johnson V, while
        # only the last two (IDs 4 and 5) have it in Johnson I. Note that only
        # the fourth star (ID = 4) has a curve in both photometric filters.

        ids = [1, 2, 3, 4]
        for star_id in ids:
            add_rcurve(star_id, ids, johnson_V)

        ids = [4, 5]
        for star_id in ids:
            add_rcurve(star_id, ids, johnson_I)

        # ********************* Johnson V *********************

        # Star 1: magnitude = 11.5
        # Star 2: abs(11.5 - 12.9) = 1.4 --> 1st
        # Star 3: abs(11.5 - 10.1) = 1.4 --> 2nd
        # Star 4: abs(11.5 - 13.4) = 2.1 --> 3rd
        # Star 5: no light curve in V, ignored
        #
        # The difference between stars 1-2 and 1-3 is exactly the same: 1.4.
        # In these cases, due to its implementation (and the fact the sorts are
        # guaranteed to be stable in Python), the method should return first
        # the star with the lowest ID, as it comes across it sooner.

        star1_similar_V = list(db.most_similar_magnitude(1, johnson_V))
        self.assertEqual(len(star1_similar_V), 3)
        expected = [(2, 12.9), (3, 10.1), (4, 13.4)]
        self.assertEqual(star1_similar_V, expected)

        # Star 2: magnitude = 12.9
        # Star 1: abs(12.9 - 11.5) = 1.4 --> 2nd
        # Star 3: abs(12.9 - 10.1) = 2.8 --> 3rd
        # Star 4: abs(12.9 - 13.4) = 0.5 --> 1st
        # Star 5: no light curve in V, ignored
        star2_similar_V = list(db.most_similar_magnitude(2, johnson_V))
        self.assertEqual(len(star2_similar_V), 3)
        expected = [(4, 13.4), (1, 11.5), (3, 10.1)]
        self.assertEqual(star2_similar_V, expected)

        # Star 3: magnitude = 10.1
        # Star 1: abs(10.1 - 11.5) = 1.4 --> 1st
        # Star 2: abs(10.1 - 12.9) = 2.8 --> 2nd
        # Star 4: abs(10.1 - 13.4) = 3.3 --> 3rd
        # Star 5: no light curve in V, ignored
        star3_similar_V = list(db.most_similar_magnitude(3, johnson_V))
        self.assertEqual(len(star3_similar_V), 3)
        expected = [(1, 11.5), (2, 12.9), (4, 13.4)]
        self.assertEqual(star3_similar_V, expected)

        # Star 4: magnitude 13.4
        # Star 1: abs(13.4 - 11.5) = 1.9 --> 2nd
        # Star 2: abs(13.4 - 12.9) = 0.5 --> 1st
        # Star 3: abs(13.4 - 10.1) = 3.3 --> 3rd
        # Star 5: no light curve in V, ignored
        star4_similar_V = list(db.most_similar_magnitude(4, johnson_V))
        self.assertEqual(len(star4_similar_V), 3)
        expected = [(2, 12.9), (1, 11.5), (3, 10.1)]
        self.assertEqual(star4_similar_V, expected)

        # Star 5: magnitude = 10.3
        # Star 1: abs(10.3 - 11.5) = 1.2 --> 2nd
        # Star 2: abs(10.3 - 12.9) = 2.6 --> 3rd
        # Star 3: abs(10.3 - 10.1) = 0.2 --> 1st
        # Star 4: abs(10.3 - 13.4) = 3.1 --> 4th
        #
        # Note that the method can be called even it the reference star (i.e.,
        # the star to which the instrumental magnitudes of all the stars are
        # compared) does not have a light curve in this filter.

        star5_similar_V = list(db.most_similar_magnitude(5, johnson_V))
        self.assertEqual(len(star5_similar_V), 4)
        expected = [(3, 10.1), (1, 11.5), (2, 12.9), (4, 13.4)]
        self.assertEqual(star5_similar_V, expected)

        # ********************* Johnson I *********************

        # Star 1: magnitude = 11.5
        # Star 2: no light curve in I, ignored
        # Star 3: no light curve in I, ignored
        # Star 4: abs(11.5 - 13.4) = 1.9 --> 2nd
        # Star 5: abs(11.5 - 10.3) = 1.2 --> 1st
        star1_similar_I = list(db.most_similar_magnitude(1, johnson_I))
        self.assertEqual(len(star1_similar_I), 2)
        expected = [(5, 10.3), (4, 13.4)]
        self.assertEqual(star1_similar_I, expected)

        # Star 2: magnitude = 12.9
        # Star 1: no light curve in I, ignored
        # Star 3: no light curve in I, ignored
        # Star 4: abs(12.9 - 13.4) = 0.5 --> 1st
        # Star 5: abs(12.9 - 10.3) = 2.6 --> 2nd
        star2_similar_I = list(db.most_similar_magnitude(2, johnson_I))
        self.assertEqual(len(star2_similar_I), 2)
        expected = [(4, 13.4), (5, 10.3)]
        self.assertEqual(star2_similar_I, expected)

        # Star 3: magnitude = 10.1
        # Star 1: no light curve in I, ignored
        # Star 2: no light curve in I, ignored
        # Star 4: abs(10.1 - 13.4) = 3.3 --> 2nd
        # Star 5: abs(10.1 - 10.3) = 0.2 --> 1st
        star3_similar_I = list(db.most_similar_magnitude(3, johnson_I))
        self.assertEqual(len(star3_similar_I), 2)
        expected = [(5, 10.3), (4, 13.4)]
        self.assertEqual(star3_similar_I, expected)

        # Star 4: magnitude = 13.4
        # Star 1: no light curve in I, ignored
        # Star 2: no light curve in I, ignored
        # Star 3: no light curve in I, ignored
        # Star 5: abs(13.4 - 10.3) = 3.1 --> 1st
        star4_similar_I = list(db.most_similar_magnitude(4, johnson_I))
        self.assertEqual(len(star4_similar_I), 1)
        expected = [(5, 10.3)]
        self.assertEqual(star4_similar_I, expected)

        # Star 5: magnitude = 10.3
        # Star 1: no light curve in I, ignored
        # Star 2: no light curve in I, ignored
        # Star 3: no light curve in I, ignored
        # Star 4: abs(10.3 - 13.4) = 3.1 --> 1st
        star5_similar_I = list(db.most_similar_magnitude(5, johnson_I))
        self.assertEqual(len(star5_similar_I), 1)
        expected = [(4, 13.4)]
        self.assertEqual(star5_similar_I, expected)

    def test_field_name(self):

        def populate_with_images(names):
            """ Return a LEMONdB that contains Images whose object names are
            taken from 'names', which must be an iterable of strings. """

            db = LEMONdB(':memory:')
            rimages = ImageTest.nrandom(len(names))
            for image, name in zip(rimages, names):
                image = image._replace(object = name)
                db.add_image(image)
            return db

        # (1) An easy case: 'FT_Tau_' is the prefix common to all the object
        # names. However, the trailing '_' is stripped, resulting in 'FT_Tau'.
        names = ['FT_Tau_15secI', 'FT_Tau_15secI', 'FT_Tau_2minV',
                 'FT_Tau_5minV',  'FT_Tau_15secI', 'FT_Tau_30secI',
                 'FT_Tau_20secI', 'FT_Tau_25secI', 'FT_Tau_25secI',
                 'FT_Tau_25secI', 'FT_Tau_5minV',  'FT_Tau_4minV',
                 'FT_Tau_3minV',  'FT_Tau_10minB', 'FT_Tau_25secI',
                 'FT_Tau_15secI', 'FT_Tau_15secI', 'FT_Tau_15secI']
        db = populate_with_images(names)
        self.assertEqual(db.field_name, 'FT_Tau')

        # (2) The method must return 'IC5146', which is common to all of the
        # images except for the last one, where lower-case letters were used.
        names = ['IC5146_30minV', 'IC5146_30minR', 'IC5146_1minI',
                 'IC5146_25minI', 'IC5146_20minV', 'ic5146_1minI']
        db = populate_with_images(names)
        self.assertEqual(db.field_name, 'IC5146')

        # (3) Note the typo ('mgc') in the second object name
        names = ['ngc2264_20minB', 'mgc2264_25minV', 'ngc2264_5minV',
                 'ngc2264_3minI',  'ngc2264_25minI', 'ngc2264_15minV']
        db = populate_with_images(names)
        self.assertEqual(db.field_name, 'ngc2264')

        # (4) Two images of FT Tau amid a series of observations of BD+78_779
        names = ['BD+78_779_30secI', 'BD+78_779_750secI', 'BD+78_779_75secV',
                 'BD+78_779_20minV', 'BD+78_779_2minB',   'BD+78_779_20minB',
                 'BD+78_779_30secI', 'BD+78_779_750secI', 'BD+78_779_3minB',
                 'FT_Tau_5minV',  'FT_Tau_15secI']
        db = populate_with_images(names)
        self.assertEqual(db.field_name, 'BD+78_779')

        # (5) The object name of a single image is also the common prefix
        name = 'BD+78_779_3minB'
        db = populate_with_images([name])
        self.assertEqual(db.field_name, name)

        # (6) None is returned if there is no common prefix
        names = ['IC5146_1minI', 'FT_Tau_5minV', 'BD+78_779_30secI']
        db = populate_with_images(names)
        self.assertEqual(db.field_name, None)

        # (7) Exception raised if there are no images in the LEMONdB
        db = LEMONdB(':memory:')
        with self.assertRaises(ValueError):
            db.field_name

    def test_set_get_and_del_metadata(self):

        db = LEMONdB(':memory:')

        key, value = 'AUTHOR', 'John Doe'
        # AttributeError is raised if the key is not in the table
        regexp = "METADATA table has no record '%s'" % key
        with self.assertRaisesRegexp(AttributeError, regexp):
            db._get_metadata(key)

        db._set_metadata(key, value)
        self.assertEqual(db._get_metadata(key), value)

        # If assigned again, the value of the record is replaced
        value = 'Jane Doe'
        db._set_metadata(key, value)
        self.assertEqual(db._get_metadata(key), value)

        # After deleting the record from the METADATA table, attempting to
        # read it again must raise AttributeError, as it no longer exists.
        db._del_metadata(key)
        with self.assertRaisesRegexp(AttributeError, regexp):
            db._get_metadata(key)

        # Numbers and None are also allowed as values
        key, value = 'ATTEMPTS', 17
        db._set_metadata(key, value)
        self.assertEqual(db._get_metadata(key), value)

        key, value = 'Pi', numpy.pi
        db._set_metadata(key, value)
        self.assertEqual(db._get_metadata(key), value)

        key, value = 'OBSERVER', None
        db._set_metadata(key, value)
        self.assertEqual(db._get_metadata(key), value)

        # These keys and values must not only be kept in memory in the LEMONdB
        # object: they are stored in the METADATA table of the LEMON database.
        # Therefore, a second LEMONdB object must be able to access them.

        key, value = 'VMIN', 12143.281
        fd, path = tempfile.mkstemp(suffix = '.LEMONdB')
        os.close(fd)
        db1 = LEMONdB(path)
        db1._set_metadata(key, value)
        db1.commit() # otherwise changes are lost
        del db1
        db2 = LEMONdB(path)
        self.assertEqual(db2._get_metadata(key), value)
        os.unlink(path)

        def test_raises(regexp, key, value):
            """ Make sure that, when db._set_metadata(key, value) is called,
            the ValueErrror exception is raised and its string representation
            matches 'regexp'.
            """
            with self.assertRaisesRegexp(ValueError, regexp):
                db._set_metadata(key, value)

        regexp = "key must be a string"
        test_raises(regexp, None, value)
        test_raises(regexp, True, value)
        test_raises(regexp, 0, value)
        test_raises(regexp, numpy.pi, value)

        regexp = "value must be a string, number or None"
        test_raises(regexp, key, [])
        test_raises(regexp, key, ())
        test_raises(regexp, key, {})
        test_raises(regexp, key, dict(one = 1, two = 2))

    def test_metadata_properties(self):

        # Expected error message of AttributeError
        regexp = "'LEMONdB' object has no attribute '%s'"

        # Make sure that the METADATA properties are effectively stored in,
        # well, the METADATA table. In order to do that, create a temporary
        # LEMONdB that resides on disk (instead of in memory, as we usually do
        # in the unit tests with ":memory:") and set the value of a property.
        # Then, open the database using a second LEMONdB object and test that
        # the value that it returns is the same that was previously assigned.

        def make_db(name, value):
            """ Return the path to a LEMON database with the 'name' property
            set to 'value'. We are responsible for deleting the database file
            when done with it.
            """
            fd, path = tempfile.mkstemp(suffix = '.LEMONdB')
            os.close(fd)
            db = LEMONdB(path)
            setattr(db, name, value)
            db.commit()
            return path

        # The different values that each property will take
        metaproperties = \
          dict(date = [time.time(), time.time() + 100],
               author = ["Jane Doe", "John Doe"],
               hostname = ['example.com', 'github.com'],
               vmin = [11345.641, None],
               vmax = [20000, 1298820.91])

        for name, values in metaproperties.iteritems():

            # AttributeError is raised if the property is not set
            db1 = LEMONdB(":memory:")
            with self.assertRaisesRegexp(AttributeError, regexp % name):
                getattr(db1, name)

            assert len(values) == 2
            first_value, second_value = values
            path = make_db(name, first_value)

            db2 = LEMONdB(path)
            self.assertEqual(getattr(db2, name), first_value)
            setattr(db2, name, second_value)
            self.assertEqual(getattr(db2, name), second_value)
            db2.commit() # otherwise changes are lost
            del db2

            # Test not only that the new value of the METADATA property has
            # been successfully updated in the database, but also that, when
            # we delete it, it is indeed removed from the table.

            db3 = LEMONdB(path)
            self.assertEqual(getattr(db3, name), second_value)
            delattr(db3, name)
            db3.commit()
            del db3

            db4 = LEMONdB(path)
            with self.assertRaisesRegexp(AttributeError, regexp % name):
                getattr(db4, name)
            del db4

            os.unlink(path)

        # Attempting to access properties that we have not explicitly declared
        # in LEMONdB raises AttributeError, as it should. We include this test
        # in case LEMONdB is refactored using __getattr__() / __setattr__(),
        # which may have dire consequences if not done with care.

        regexp = "'LEMONdB' object has no attribute 'foo'"
        with self.assertRaisesRegexp(AttributeError, regexp):
            db = LEMONdB(':memory:')
            db.foo

    def test_star_closest_to_world_coords(self):

        db = LEMONdB(':memory:')

        # LEMONdB is empty, raises ValueError
        point = (24.19933, 41.40547) # or any other value
        with self.assertRaises(ValueError):
            db.star_closest_to_world_coords(*point)

        star_1 = list(self.random_star_info(id_ = 1))
        star_1[3:5] = (186.612871, 31.223544) # NGC 4414
        db.add_star(*star_1)

        star_2 = list(self.random_star_info(id_ = 2))
        star_2[3:5] = (5.166987, 31.989942) # WASP-1 b
        db.add_star(*star_2)

        star_3 = list(self.random_star_info(id_ = 3))
        star_3[3:5] = (278.064554, 6.945753) # HD 171028 b
        db.add_star(*star_3)

        star_4 = list(self.random_star_info(id_ = 4))
        star_4[3:5] = (97.69629, 58.16264) # 6 Lyn b
        db.add_star(*star_4)

        # First point: ra = 73.68192, dec = 12.35219 (HD 31253 b)
        # Angular distances from (73.68192, 12.35219) to:
        # Star 1 (186.612871, 31.223544) = 102.39093634720719
        # Star 2 (  5.166987, 31.989942) =  65.368742708468304
        # Star 3 (278.064554,  6.945753) = 149.01762956453169
        # Star 4 ( 97.69629,  58.16264)  =  49.274795088346579 <- (shortest)

        point_1 = (73.68192, 12.35219)
        star_id, distance = db.star_closest_to_world_coords(*point_1)
        self.assertEqual(star_id, 4)
        self.assertAlmostEqual(distance, 49.274795088346579)

        # Second point: ra = 346.869646, dec = 21.134251 (HR 8799 d)
        # Angular distances from (346.869646, 21.134251) to:
        # Star 1 (186.612871, 31.223544) = 124.32181703665258
        # Star 2 (  5.166987, 31.989942) =  19.591535371385636 <- (shortest)
        # Star 3 (278.064554,  6.945753) =  67.768448323818134
        # Star 4 ( 97.69629,  58.16264)  =  82.451118082774897

        point_2 = (346.869646, 21.134251)
        star_id, distance = db.star_closest_to_world_coords(*point_2)
        self.assertEqual(star_id, 2)
        self.assertAlmostEqual(distance, 19.591535371385636)

        # Third point: ra = 227.21558, dec = 2.34333 (WASP-24 b)
        # Angular distances from (227.21558, 2.34333) to:
        # Star 1 (186.612871, 31.223544) =  47.939281840122732 <- (shortest)
        # Star 2 (  5.166987, 31.989942) = 127.41779312476839
        # Star 3 (278.064554,  6.945753) =  50.864720140121634
        # Star 4 ( 97.69629,  58.16264)  = 107.49712742084881

        point_3 = (227.21558, 2.34333)
        star_id, distance = db.star_closest_to_world_coords(*point_3)
        self.assertAlmostEqual(star_id, 1)
        self.assertEqual(distance, 47.939281840122732)

