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

import itertools
import numpy
import operator
import os
import random
import tempfile
import time
import unittest

import passband
from database import DBStar
from database import Image
from database import LightCurve
from database import PhotometricParameters
from database import ReferenceImage

from diffphot import Weights

NITERS = 100      # How many times each test case is run with random data
MIN_NSTARS = 10   # Minimum number of items for random collections of DBStars
MAX_NSTARS = 100  # Maximum number of items for random collections of DBStars

class DBStarTest(unittest.TestCase):

    MIN_ID = 1     # Minimum ID for random DBStars
    MAX_ID = 9999  # Maximum ID for random DBStars
    MIN_SIZE = 1   # Minimum number of photometric records
    MAX_SIZE = 250 # Maximum number of photometric records
    MIN_UNIX_TIME = 0            # Thu Jan  1 01:00:00 1970
    MAX_UNIX_TIME = time.time()  # Minimum and maximum Unix time
    MIN_MAG = 1.47   # Minimum value for random magnitudes (Saturn)
    MAX_MAG = 25     # Maximum value for random magnitudes (Fenrir)
    MIN_SNR = 2      # Minimum value for random signal-to-noise ratios
    MAX_SNR = 10000  # Maximum value for random signal-to-noise ratios

    @classmethod
    def random_data(cls):
        """ Return the information needed to instantiate a random DBStar """

        id_ = random.randint(cls.MIN_ID, cls.MAX_ID)
        pfilter = passband.Passband.random()
        size = random.randint(cls.MIN_SIZE, cls.MAX_SIZE)
        phot_info = numpy.empty((3, size))
        times_indexes = {}
        for index in xrange(size):
            unix_time = random.uniform(cls.MIN_UNIX_TIME, cls.MAX_UNIX_TIME)
            magnitude = random.uniform(cls.MIN_MAG, cls.MAX_MAG)
            snr = random.uniform(cls.MIN_SNR, cls.MAX_SNR)
            phot_info[:,index] = unix_time, magnitude, snr
            times_indexes[unix_time] = index
        return id_, pfilter, phot_info, times_indexes

    @classmethod
    def random(cls):
        """ Return a random DBStar """
        args = cls.random_data()
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
        self.assertRaises(ValueError, DBStar, id_, pfilter, phot_info, times_indexes)

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

        for index in xrange(to_add):
            unix_time = random.uniform(cls.MIN_UNIX_TIME, cls.MAX_UNIX_TIME)
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
        self.assertRaises(KeyError, star._trim_to, nonsubset)

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
        cls = self.__class__
        for _ in xrange(NITERS):
            id_, pfilter, phot_info, times_indexes = cls.random_data()
            # Construct the row of three-element tuples, out of this random
            # NumPy array, and then check that the array in the DBStar returned
            # by DBStar.make_star is equal than this input, original array.
            rows = [phot_info[:, index] for index in xrange(phot_info.shape[1])]
            star = DBStar.make_star(id_, pfilter, rows)
            self.assertEqual(star.id, id_)
            self.assertEqual(star.pfilter, pfilter)
            self.assertTrue(numpy.all(numpy.equal(star._phot_info, phot_info)))
            self.assertEqual(star._time_indexes, times_indexes)


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

    MIN_UNIX_TIME = 0            # Thu Jan  1 01:00:00 1970
    MAX_UNIX_TIME = time.time()  # Minimum and maximum Unix time
    MIN_AIRMASS = 1     # Minimum value for random airmasses (zenith)
    MAX_AIRMASS = 2.92  # Maximum value for random airmasses (~70 dec)
    MIN_GAIN = 1   # Minimum value for random CCD gains
    MAX_GAIN = 10  # Maximum value for random CCD gains
    MIN_XOFFSET = MIN_YOFFSET = 0    # Minimum and maximum values for random
    MAX_XOFFSET = MAX_YOFFSET = 500  # ... offsets in both axes (pixels)
    MIN_XOVERLAP = MIN_YOVERLAP = 0    # Minimum and maximum values for random
    MAX_XOVERLAP = MAX_YOVERLAP = 2000 # ... numbers of stars that overlap

    @classmethod
    def random_data(cls):
        """ Return the information needed to instantiate a random Image """

        fd, path = tempfile.mkstemp(suffix = '.fits')
        os.close(fd)
        os.unlink(path)

        pfilter = passband.Passband.random()
        pparams = PhotometricParametersTest.random()
        unix_time = random.uniform(cls.MIN_UNIX_TIME, cls.MAX_UNIX_TIME)
        airmass = random.uniform(cls.MIN_AIRMASS, cls.MAX_AIRMASS)
        gain = random.uniform(cls.MIN_GAIN, cls.MAX_GAIN)
        xoffset = random.uniform(cls.MIN_XOFFSET, cls.MAX_XOFFSET)
        yoffset = random.uniform(cls.MIN_YOFFSET, cls.MAX_YOFFSET)
        xoverlap = random.randint(cls.MIN_XOVERLAP, cls.MAX_XOVERLAP)
        yoverlap = random.randint(cls.MIN_YOVERLAP, cls.MAX_YOVERLAP)

        return (path, pfilter, pparams, unix_time, airmass,
                gain, xoffset, yoffset, xoverlap, yoverlap)

    @classmethod
    def random(cls):
        """ Return a random Image """
        args = cls.random_data()
        return Image(*args)

    def test_init_(self):
        for _ in xrange(NITERS):
            args = self.random_data()
            img = Image(*args)
            self.assertEqual(img.path, args[0])
            self.assertEqual(img.pfilter, args[1])
            pparams = args[2]
            self.assertEqual(img.pparams.aperture, pparams.aperture)
            self.assertEqual(img.pparams.annulus, pparams.annulus)
            self.assertEqual(img.pparams.dannulus, pparams.dannulus)
            self.assertEqual(img.unix_time, args[3])
            self.assertEqual(img.airmass, args[4])
            self.assertEqual(img.gain, args[5])
            self.assertEqual(img.xoffset, args[6])
            self.assertEqual(img.yoffset, args[7])
            self.assertEqual(img.xoverlap, args[8])
            self.assertEqual(img.yoverlap, args[9])


class ReferenceImageTest(unittest.TestCase):

    @classmethod
    def random_data(cls):
        """ Return the args needed to instantiate a random ReferenceImage"""
        args = ImageTest.random_data()
        # Return the first six elements, except for 'pparams'
        return args[:2] + args[3:6]

    def test_init_(self):
        for _ in xrange(NITERS):
            args = self.random_data()
            path, pfilter, unix_time, airmass, gain = args
            rimage = ReferenceImage(*args)
            self.assertEqual(rimage.path, path)
            self.assertEqual(rimage.pfilter, pfilter)
            self.assertEqual(rimage.unix_time, unix_time)
            self.assertEqual(rimage.airmass, airmass)
            self.assertEqual(rimage.gain, gain)


class LightCurveTest(unittest.TestCase):

    MIN_ID = 1     # Minimum value for random DBStar IDs
    MAX_ID = 9999  # Maximum value for random DBStar IDs
    MIN_UNIX_TIME = 0            # Thu Jan  1 01:00:00 1970
    MAX_UNIX_TIME = time.time()  # Minimum and maximum Unix time
    MIN_MAG = 1.47   # Minimum value for random magnitudes (Saturn)
    MAX_MAG = 25     # Maximum value for random magnitudes (Fenrir)
    MIN_SNR = 2      # Minimum value for random signal-to-noise ratios
    MAX_SNR = 10000  # Maximum value for random signal-to-noise ratios

    @classmethod
    def random_point(cls):
        """ Return a random (time, magnitude, snr) tuple """
        unix_time = random.uniform(cls.MIN_UNIX_TIME, cls.MAX_UNIX_TIME)
        magnitude = random.uniform(cls.MIN_MAG, cls.MAX_MAG)
        snr = random.uniform(cls.MIN_SNR, cls.MAX_SNR)
        return unix_time, magnitude, snr

    @classmethod
    def random_data(cls):
        """ Return the args needed to instantiate a random LightCurve"""
        pfilter = passband.Passband.random()
        size = random.randint(MIN_NSTARS, MAX_NSTARS)
        cstars = [random.randint(cls.MIN_ID, cls.MAX_ID) for x in xrange(size)]
        cweights = Weights.random(size)
        return pfilter, cstars, cweights

    @classmethod
    def random(cls):
        """ Return a random, empty LightCurve"""
        args = cls.random_data()
        return LightCurve(*args)

    def test_init(self):
        for _ in xrange(NITERS):
            pfilter, cstars, cweights = args = self.random_data()
            curve = LightCurve(*args)
            self.assertEqual(curve.pfilter, pfilter)
            self.assertEqual(curve.cstars, cstars)
            self.assertTrue(numpy.all(numpy.equal(curve.cweights, cweights)))

        # ValueError raised if number of weigts != number of comparison stars
        pfilter, cstars, cweights = self.random_data()
        rindex = random.choice(range(len(cweights)))
        cweights = cweights.rescale(rindex) # remove rindex-th coefficient
        assert len(cstars) != len(cweights)
        self.assertRaises(ValueError, LightCurve, pfilter, cstars, cweights)

        # Same exception also raised if there are no comparison stars
        cstars = cweights = []
        self.assertRaises(ValueError, LightCurve, pfilter, cstars, cweights)

    def test_add_len_and_getitem(self):
        for _ in xrange(NITERS):
            curve = self.random()
            size = random.randint(MIN_NSTARS, MAX_NSTARS)
            for index in xrange(size):
                point = self.random_point()
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
        self.assertRaises(StopIteration, lambda: points.next())

        # ... and the random test cases
        for _ in xrange(NITERS):
            curve = self.random()
            size = random.randint(MIN_NSTARS, MAX_NSTARS)
            curve_points = [self.random_point() for x in xrange(size)]
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
            for x in xrange(size):
                unix_time, mag, snr = point = self.random_point()
                curve.add(*point)
                magnitudes.append(mag)
            self.assertAlmostEqual(curve.stdev, numpy.std(magnitudes))

        # ValueError is raised if the LightCuve is empty
        curve = self.random()
        assert not len(curve)
        self.assertRaises(ValueError, lambda: curve.stdev)

    def test_weights(self):
        for _ in xrange(NITERS):
            pfilter, cstars, cweights = self.random_data()
            curve = LightCurve(pfilter, cstars, cweights)
            for index, (cstar_id, cweight) in enumerate(curve.weights()):
                self.assertEqual(cstar_id, cstars[index])
                self.assertEqual(cweight, cweights[index])

