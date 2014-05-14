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

import copy
import math
import numpy
import random

from test import unittest
import passband
import test_database
from database import DBStar
from diffphot import Weights, StarSet

NITERS = 50  # How many times some test cases are run with random data

def assertSequencesAlmostEqual(cls, first, second):
    """ Test that the elements of these two sequences are almost equal.

    The first parameter, 'cls', must be the subclass of unittest.TestCase
    instance in which this test is done. Values are considered to be almost
    equal according to the rules of the TestCase.assertAlmostEqual method.

    """

    if len(first) != len(second):
        cls.fail("arrays must be of the same size")
    for values in zip(first, second):
        cls.assertAlmostEqual(*values)

class WeightsTest(unittest.TestCase):

    MIN_COEFFS = 2   # Minimum number of coefficients for random Weights
    MAX_COEFFS = 99  # Maximum number of coefficients for random Weights

    @staticmethod
    def rcoeffs(n):
        """ Return a list of n [0.0, 1.0) floating point numbers """
        return [random.random() for index in xrange(n)]

    def setUp(self):
        """ Prepare the test fixture.

        Sets the 'rcoefficients' attibute to a list of random size, each of its
        elements being a random list of coefficients in the [0.0, 1.0) range"""

        self.rcoefficients = []
        for x in xrange(NITERS):
            size = random.randint(self.MIN_COEFFS, self.MAX_COEFFS)
            self.rcoefficients.append(self.rcoeffs(size))

    def test_new(self):
        self.assertRaises(ValueError, Weights, [])
        for coeffs in self.rcoefficients:
            w = Weights(coeffs)
            assertSequencesAlmostEqual(self, w, coeffs)

    def test_rescale(self):

        # Basic, specific test cases
        w1 = Weights([0.5, 0.3, 0.2])
        assertSequencesAlmostEqual(self, w1.rescale(0), [0.6, 0.4])
        assertSequencesAlmostEqual(self, w1.rescale(1), [5/7, 2/7])
        assertSequencesAlmostEqual(self, w1.rescale(2), [5/8, 3/8])

        w2 = Weights([0.2, 0.4, 0.4])
        assertSequencesAlmostEqual(self, w2.rescale(0), [0.5, 0.5])
        assertSequencesAlmostEqual(self, w2.rescale(1), [1/3, 2/3])
        assertSequencesAlmostEqual(self, w2.rescale(2), [1/3, 2/3])

        w3 = Weights([3])
        self.assertRaises(ValueError, w3.rescale, 0)

        # Random cases in which one of the coefficients is discarded
        for coeffs in self.rcoefficients:
            weights  = Weights(coeffs)
            oweights = copy.deepcopy(weights)

            index = random.choice(xrange(len(weights)))
            rweights = weights.rescale(index)

            # The original instance must not have been modified
            assertSequencesAlmostEqual(self, weights, oweights)
            self.assertAlmostEqual(rweights.total, 1.0)
            self.assertEqual(len(rweights), len(weights)- 1)

            # Test that normalized coefficients have the expected value
            ecoeffs = coeffs[:index] + coeffs[index+1:]
            eweights = Weights(ecoeffs).normalize()
            assertSequencesAlmostEqual(self, rweights, eweights)

    def test_total(self):
        for coeffs in self.rcoefficients:
            w = Weights(coeffs)
            self.assertAlmostEqual(w.total, math.fsum(coeffs))

    def test_normalize(self):
        for coeffs in self.rcoefficients:
            nweights = Weights(coeffs).normalize()
            self.assertAlmostEqual(nweights.total, 1.0)
            eweights = [x / math.fsum(coeffs) for x in coeffs]
            assertSequencesAlmostEqual(self, nweights, eweights)

    def test_inversely_proportional(self):

        # Basic, specific test cases
        w1 = Weights.inversely_proportional([3])
        assertSequencesAlmostEqual(self, w1, [1])

        w2 = Weights.inversely_proportional([1, 2])
        assertSequencesAlmostEqual(self, w2, [2/3, 1/3])

        w3 = Weights.inversely_proportional([1, 3])
        assertSequencesAlmostEqual(self, w3, [3/4, 1/4])

        w4 = Weights.inversely_proportional([1, 2, 3])
        assertSequencesAlmostEqual(self, w4, [6/11, 3/11, 2/11])

        self.assertRaises(ValueError, Weights.inversely_proportional, [])
        self.assertRaises(ValueError, Weights.inversely_proportional, [0, 1, 2])

        for coeffs in self.rcoefficients:
            iweights = Weights.inversely_proportional(coeffs)
            icoeffs  = [1 / x for x in coeffs]
            eweights = [1 / x / math.fsum(icoeffs) for x in coeffs]
            self.assertAlmostEqual(iweights.total, 1.0)
            assertSequencesAlmostEqual(self, eweights, iweights)

    def test_absolute_percent_change(self):

        # Basic, specific test cases
        w1 = Weights([1, 2, 3])
        w2 = Weights([1, 2, 3])
        self.assertAlmostEqual(w1.absolute_percent_change(w2), 0.0)

        w1 = Weights([1, 4, 1])
        w2 = Weights([1, 2, 1])
        self.assertAlmostEqual(w1.absolute_percent_change(w2), 0.5)

        w1 = Weights([1, 0, 4])  # 0->9 ignored
        w2 = Weights([1, 9, 5])
        self.assertAlmostEqual(w1.absolute_percent_change(w2), 0.25)

        w1 = Weights([1, 4, 4])
        w2 = Weights([1, 0, 5])
        self.assertAlmostEqual(w1.absolute_percent_change(w2), 1)

        # An example with decimals:
        # 1st coefficient, percent change = -0.37125788 <-- maximum (abs)
        # 2nd coefficient, percent change =  0.354215852
        # 3rd coefficient, percent change = -0.006970987
        w1 = Weights([0.295238595044, 0.317072842106, 0.38768856285])
        w2 = Weights([0.185628940097, 0.429385068928, 0.38498599097])
        self.assertAlmostEqual(w1.absolute_percent_change(w2), 0.37125788019233952)

        # ValueError raised if sequences are of a different size
        w1 = Weights([1, 2, 3])
        w2 = Weights([1, 2, 3, 4])
        self.assertRaises(ValueError, w1.absolute_percent_change, w2)

        # We also need to test the 'minimum' keyword:
        # 1st coefficient, percent change =  2.5 <-- maximum
        # 2nd coefficient, percent change = -0.069999999999999951 <-- 2nd highest
        # 3rd coefficient, percent change =  0.047619047619047658
        w1 = Weights([0.5,  1,    2.1])
        w2 = Weights([1.75, 0.93, 2.2])

        # If 'minimum = 0.75' is used, the first coefficient is ignored (as one
        # of the values, 0.5, is smaller than this threshold), so the maximum
        # absolute percent change is found among the two remaining coefficients
        # and the right answer becomes 0.09524.
        pc = w1.absolute_percent_change(w2, minimum = 0.75)
        self.assertAlmostEqual(pc, 0.069999999999999951)

        # If 'minimum' is set to 1.25 the second coefficient is also ignored
        pc = w1.absolute_percent_change(w2, minimum = 1.25)
        self.assertAlmostEqual(pc, 0.047619047619047658)

        # Finally, a too high minimum value discards all the weights, which
        # causes the method to raise the ValueError exception. This also
        # happens in the value of 'minimum' is negative or zero.
        args = w1.absolute_percent_change, w2
        self.assertRaises(ValueError, *args, minimum = 2.5)
        self.assertRaises(ValueError, *args, minimum = -0.132)
        self.assertRaises(ValueError, *args, minimum =  0)

    def test_random(self):
        for _ in xrange(NITERS):
            size = random.randint(self.MIN_COEFFS, self.MAX_COEFFS)
            w = Weights.random(size)
            self.assertEqual(len(w), size)
            self.assertAlmostEqual(w.total, 1.0)


class StarSetTest(unittest.TestCase):

    # These two-element tuples (a, b) define the range from within random
    # integers N such that a <= N <= b will be chosen when random sets of
    # data are used by some tests cases.

    NSTARS_RANGE  = (10, 50)   # Number of DBStars in random StarSets
    NRECORDS_RANGE = (25, 100) # Number of records in DBStars
    IDS_RANGE = (1, 9999) # IDs for random DBStars

    # And these tuples define the range from within random floating point
    # numbers N such that a <= N <= b will be chosen for some random sets
    # of data.

    MAG_RANGE = (1.47, 25) # Instrumental magnitudes (Saturn, Fenrir)
    SNR_RANGE = (2, 10000) # Signal-to-noise ratios

    def rDBStars(self, size = None, nrecords = None, pfilter = None):
        """ Return a random list of DBStars.

        The random DBStars are guaranteed to have the same photometric filter
        and records for the same Unix times, as well as different IDs. If the
        'size', 'nrecords' or 'pfilter' parameters are not given, random
        numbers of stars and photometric records or a random photometric
        filter are used, respectively.

        """

        if size is None:
            size = random.randint(*self.NSTARS_RANGE)
        elif size <= 0:
            raise ValueError("at least one star is needed")

        if nrecords is None:
            nrecords = random.randint(*self.NRECORDS_RANGE)
        if pfilter is None:
            pfilter = passband.Passband.random()

        stars_ids = random.sample(xrange(*self.IDS_RANGE), size)
        unix_times = test_database.runix_times(nrecords)

        stars = []
        for star_id in stars_ids:
            rows = []
            for unix_time in unix_times:
                magnitude = random.uniform(*self.MAG_RANGE)
                snr = random.uniform(*self.SNR_RANGE)
                rows.append((unix_time, magnitude, snr))
            stars.append(DBStar.make_star(star_id, pfilter, rows))

        self.assertEqual(size, len(stars))
        # There must be no duplicate IDs
        self.assertEqual(size, len(set(s.id for s in stars)))
        # All the DBStars must have the same photometric filter...
        self.assertEqual(set([pfilter]), set(s.pfilter for s in stars))
        # ... and records for the same Unix times
        for s in stars:
            self.assertEqual(set(s._unix_times), set(unix_times))

        return stars

    def test_init_and_get(self):

        # DBStars are internally stored in a three-dimensional, private NumPy
        # array. Accessing this array element by element in order to check if
        # it has the expected values would be extremely ugly, unpythonic and
        # would fail the day the inner workings of the class were changed.
        # Instead, we test the __init__ method by feeding some DBStars to it
        # and then verifying that they are equal to those returned by the
        # __get__ method. We do not care about what the method does internally
        # as long as it returns the same DBStars.
        #
        # Note that, implicitly, this test case also verifies the StarSet._add
        # method, since it is internally called by StarSet.__init__ to parse
        # and store in the NumPy array each of the DBStars.

        for _ in xrange(NITERS):
            stars = self.rDBStars()
            set_ = StarSet(stars)

            for istar, ostar in zip(stars, set_):
                _eq_ = test_database.DBStarTest.equal
                self.assertTrue(_eq_(istar, ostar))

        # The sequence of stars cannot be empty
        self.assertRaises(ValueError, StarSet, [])

        # The DBStars cannot be empty either
        stars = self.rDBStars(nrecords = 0)
        self.assertRaises(ValueError, StarSet, stars)

        # All DBStars must have the same photometric filter...
        stars = self.rDBStars()
        rstar = random.choice(stars)
        rstar.pfilter = rstar.pfilter.different()
        self.assertEqual(2, len(set(s.pfilter for s in stars)))
        self.assertRaises(ValueError, StarSet, stars)

        # ... and Unix times. To check for this we first generate a list of
        # random DBStars and then add to it another one with the same number
        # of photometric records and filter, but different Unix times.
        stars = self.rDBStars()
        while True:
            dstar = self.rDBStars(size = len(stars),
                                  nrecords = len(stars[0]),
                                  pfilter = stars[0].pfilter)[0]
            if set(dstar._unix_times) != (stars[0]._unix_times):
                stars.append(dstar)
                break

        self.assertRaises(ValueError, StarSet, stars)

        # There cannot be duplicate IDs among the DBStars
        stars = self.rDBStars()
        rstar1, rstar2 = random.sample(stars, 2)
        rstar1.id = rstar2.id
        self.assertEqual(len(stars) - 1, len(set(s.id for s in stars)))
        self.assertRaises(ValueError, StarSet, stars)

    def test_star_ids_len_and_nimages(self):

        # Generate some random DBStars and take note of which are their IDs and
        # number of photometric records. Then, create a set with them and check
        # that values returned by StarSet.__len__ and nimages are equal to the
        # number of input stars and their number of photometric records,
        # respectively, and that the IDs returned by star_ids are the expected.

        for _ in xrange(NITERS):
            stars = self.rDBStars()
            stars_ids = [s.id for s in stars]
            set_ = StarSet(stars)
            self.assertEqual(set_.star_ids, stars_ids)
            self.assertEqual(len(set_), len(stars))
            self.assertEqual(set_.nimages, len(stars[0]))

    def test_delitem(self):

        # Instantiate a StarSet from a sequence of random DBStars, then
        # randomly choose one of them and delete it from the StarSet. The stars
        # in the set up to the index of the removed one should be the same, and
        # those after it should have been moved one position to the left.

        for _ in xrange(NITERS):
            stars = self.rDBStars()
            set_ = StarSet(stars)
            rindex = random.choice(xrange(len(stars)))
            del set_[rindex]

            self.assertEqual(len(set_), len(stars) - 1)
            _eq_ = test_database.DBStarTest.equal
            for index in xrange(rindex):
                self.assertTrue(_eq_(set_[index], stars[index]))
            for index in xrange(rindex, len(set_)):
                self.assertTrue(_eq_(set_[index], stars[index + 1]))

    def _populate_set(self, magnitudes, snrs = None):
        """ Return a random StarSet whose stars have the specified magnitudes.

        The method returns a StarSet with a random photometric filter, Unix
        times and signal-to-noise ratios, but the magnitudes of whose stars are
        defined by the values given in the 'magnitudes' parameter. It must be a
        sequence of sequences, each one of the latter containing the magnitudes
        of that star in the different images. For example, magnitudes[0][2]
        gives the magnitude of the first star in the third image, while the
        magnitude of the fourth star in the second image will be the value
        stored in magnitudes[3][1].

        It goes without saying that the second-level sequences (i.e., each one
        of the elements in 'magnitudes') must have all the same length. Also,
        at least the magnitudes of one star must be given, as StarSet instances
        cannot be empty.

        These instrumental magnitudes are assigned random signal-to-noise
        ratios unless the 'snrs' keyword argument is given. It must be a second
        sequence of sequences, of the same dimensions that 'magnitudes', giving
        for each star the signal-to-noise ratio of each of its instrumental
        magnitudes. To continue with the above example, snrs[0][2] would be
        the SNR of the first star in the third image.

        """

        if not len(magnitudes):
            msg = "the magnitudes of at least one star are needed"
            self.fail(msg)

        if len(set(len(mags) for mags in magnitudes)) != 1:
            msg = "the same number of magnitudes is needed for all stars"
            self.fail(msg)


        # If the sequence of sequences with the signal-to-noise ratio of
        # each instrumental magnitude is not given, use random values.
        if not snrs:
            snrs = []
            n = len(magnitudes[0]) # number of magnitudes of each star
            for _ in xrange(len(magnitudes)):
                star_snrs = [random.uniform(*self.SNR_RANGE) for x in xrange(n)]
                snrs.append(star_snrs)

        if [len(x) for x in magnitudes] != [len(x) for x in snrs]:
            msg = "'magnitudes' and 'snrs' must have the same dimensions"
            self.fail(msg)

        size = len(magnitudes)
        nrecords = len(magnitudes[0])

        star_ids = random.sample(xrange(*self.IDS_RANGE), size)
        unix_times = sorted(test_database.runix_times(nrecords))
        pfilter = passband.Passband.random()

        stars = []
        for star_id, star_mags, star_snrs in zip(star_ids, magnitudes, snrs):
            rows = []
            args = unix_times, star_mags, star_snrs
            for unix_time, magnitude, snr in zip(*args):
                rows.append((unix_time, magnitude, snr))
            stars.append(DBStar.make_star(star_id, pfilter, rows))

        set_ = StarSet(stars)
        self.assertEqual(len(set_), size)
        self.assertEqual(set_.nimages, nrecords)
        self.assertEqual(set_.star_ids, star_ids)
        for star, star_mags, star_snrs in zip(set_, magnitudes, snrs):
            for index in xrange(len(star)):
                self.assertEqual(star.mag(index), star_mags[index])
                self.assertEqual(star.snr(index), star_snrs[index])
        return set_

    def test_flux_proportional_weights(self):

        star_mags = \
        [[6.4, 6.9],  # star = 0
         [8.5, 8.7]]  # star = 1
        # img1 img2

        # (1) Normalized (for each image)
        # [0.7529411764705882, 0.7931034482758622]
        # [1.                , 1.                ]
        #
        # (2) Median (of each star)
        # [0.77302231237322516, 1.0]
        #
        # (3) Magnitude-inversely proportional weights:
        expected = Weights([0.55207403, 0.44792597])
        weights = self._populate_set(star_mags).flux_proportional_weights()
        assertSequencesAlmostEqual(self, weights, expected)

        star_mags = \
        [[12.4, 11.3, 12.3, 11.8], # star = 0
         [13.4, 12.5, 13.2, 12.7], # star = 1
         [11.3, 10.3, 11.4, 10.7]] # star = 2
        # img1  img2  img3  img4

        # (1) Normalized (for each image)
        # [0.92537313432835822, 0.904, 0.93181818181818188, 0.92913385826771666]
        # [1.                 , 1.   , 1.                 , 1.                 ]
        # [0.84328358208955223, 0.824, 0.86363636363636376, 0.84251968503937003]
        #
        # (2) Median (of each star)
        # [0.92725349629803744, 1. , 0.84290163356446113]
        #
        # (3) Magnitude-inversely proportional weights:
        expected = Weights([0.33156698, 0.31007922, 0.3583538])
        weights = self._populate_set(star_mags).flux_proportional_weights()
        assertSequencesAlmostEqual(self, weights, expected)

        star_mags = \
        [[11.3, 12.3, 14.3, 12.1], # star = 0
         [10.4, 11.8, 13.4, 11.4], # star = 1
         [12.8, 14.2, 15.8,  8.1], # star = 2
         [10.4, 11.4, 13.0, 11.9]] # star = 3
        # img1  img2  img3  img4

        # (1) Normalized (for each image)
        # [0.8828125, 0.8661971830985916, 0.9050632911392406, 1.                ]
        # [0.8125   , 0.8309859154929579, 0.8481012658227848, 0.9421487603305786]
        # [1.       , 1.                , 1.                , 0.6694214876033058]
        # [0.8125   , 0.8028169014084507, 0.8227848101265822, 0.9834710743801653]
        #
        # (2) Median (of each star)
        # [0.89393789556962022, 0.83954359065787132, 1. , 0.81764240506329111]
        #
        # (3) Magnitude-inversely proportional weights:
        expected = Weights([0.24807083, 0.26081552, 0.2249836 , 0.26613004])
        weights = self._populate_set(star_mags).flux_proportional_weights()
        assertSequencesAlmostEqual(self, weights, expected)

    def random_set(self, size = None, nrecords = None):
        """ Return a random StarSet and a DBStar for which it is complete.

        The method returns a two-element tuple: a StarSet of random photometric
        filter, 'size' stars and 'nrecords' photometric records, and a DBStar
        with the same Unix times that the set has. The StarSet is thus said to
        be 'complete' for the star (see method DBStar.complete_for) and can be
        used to generate its light curve (method StarSet.light_curve). If the
        size or number of photometric records of the StarSet are not given,
        random values are used.

        """

        if size is None:
            size = random.randint(*self.NSTARS_RANGE)

        set_ = StarSet(self.rDBStars(size = size + 1, nrecords = nrecords))
        star = set_[-1]
        del set_[-1]
        self.assertEqual(len(set_), size)
        self.assertEqual(list(star._unix_times), list(set_._unix_times))
        self.assertTrue(star.id not in set_.star_ids)
        return set_, star

    def _assert_light_curve(self, cmp_stars, weights, star_mags, star_snrs,
                            expected_mags, expected_snrs, _exclude_index = None):
        """ Assert that the StarSet class generates the expected LightCurve.

        The method receives a StarSet, 'cmp_stars', and creates a random DBStar
        with the same photometric filter and Unix times that the set, and which
        uses the instrumental magnitudes and signal-to-noise ratios stored in
        'star_mags' and 'star_snrs', respectively. Then it computes its light
        curve (method StarSet.light_curve), using the given Weights, and
        compares the differential magnitudes and SNRs to those defined in,
        respectively, 'expected_mags' and 'expected_snrs'. If the values are
        not 'almost' equal, as per the rules of the TestCase.assertAlmostEqual
        method, an unconditional test failure is signaled.

        The value of the '_exclude_index' keyword argument is passed down to
        StarSet.light_curve. The value of 'no_snr', the other keyword argument
        of the method, needs not be set: each light curve is in actuality
        computed twice, once with signal-to-noise ratios and once without
        them. In this manner, we always verify that the differential magnitudes
        are the same regardless of the value of 'no_snr', and that the SNRs are
        set to None when it evaluates to True.

        It is probably sort of obvious that there must be as many weights as
        stars in the StarSet and that 'star_mags', 'star_snrs', 'expected_mags'
        and 'expected_snrs' must have the same length, which must also match
        the number of photometric records of the comparison stars. If any of
        these conditions is not met, an error is raised.

        """

        if cmp_stars.nimages != len(star_mags):
            msg = ("number of instrumental magnitudes of the DBStar must be "
                   "equal to the number of records of each comparison star")
            self.fail(msg)

        if len(star_mags) != len(star_snrs):
            msg = "each magnitude of the DBStar must have a SNR"
            self.fail(msg)

        if len(star_mags) != len(expected_mags):
            msg = ("number of instrumental magnitudes of the DBStar must "
                   "match that of the expected differential magnitudes")
            self.fail(msg)

        if len(expected_mags) != len(expected_snrs):
            msg = "each expected differential magnitude must have a SNR"
            self.fail(msg)

        rows = []
        args = cmp_stars._unix_times, star_mags, star_snrs
        for unix_time, magnitude, snr in zip(*args):
            rows.append((unix_time, magnitude, snr))

        # Make sure not to select an ID already in the StarSet
        candidate_ids = set(xrange(*self.IDS_RANGE)) - set(cmp_stars.star_ids)
        id_ = random.choice(list(candidate_ids))
        star = DBStar.make_star(id_, cmp_stars.pfilter, rows)

        kwargs = dict(_exclude_index = _exclude_index, no_snr = False)
        curve = cmp_stars.light_curve(weights, star, **kwargs)
        returned_mags, returned_snrs = zip(*[x[1:] for x in curve])
        assertSequencesAlmostEqual(self, returned_mags, expected_mags)
        assertSequencesAlmostEqual(self, returned_snrs, expected_snrs)

        # The differental magnitudes must be the same when the 'no_snr' keyword
        # argument is set to True. It is only the signal-to-noise ratios, which
        # are set to None, that should be affected.
        kwargs['no_snr'] = True
        curve = cmp_stars.light_curve(weights, star, **kwargs)
        returned_mags, returned_snrs = zip(*[x[1:] for x in curve])
        assertSequencesAlmostEqual(self, returned_mags, expected_mags)
        self.assertEqual(set(returned_snrs), set([None]))

    def test_light_curve_errors(self):

        nstars = random.randint(*self.NSTARS_RANGE)
        set_, star = self.random_set(nstars)

        # ValueError if the number of coefficients in Weights does not match
        # that of stars in the StarSet.
        weights = Weights.random(nstars - 1)
        args = set_.light_curve, weights, star
        self.assertRaises(ValueError, *args)
        weights = Weights.random(nstars + 1)
        self.assertRaises(ValueError, *args)
        weights = Weights.random(nstars)
        set_.light_curve(weights, star) # this works!

        # The Unix times of the DBStar must exactly match those of the StarSet;
        # in practical terms this means that they will be those returned by the
        # DBStar.complete_for method, which identifies precisely the stars that
        # can be used to create the artificial comparison star.
        while True:
            # Loop until we get a random DBStar with different Unix times
            kwargs = dict(size = 1, nrecords = len(star), pfilter = star.pfilter)
            dstar = self.rDBStars(**kwargs)[0]
            if set(dstar._unix_times) != set(star._unix_times):
                weights = Weights.random(nstars)
                args = set_.light_curve, weights, dstar
                self.assertRaises(ValueError, *args)
                break

    def test_light_curve_basic_case(self):

        cmp_mags = \
        [[14.5, 14.2, 13.4], # (a)
         [13.2, 13.1, 12.8], # (b)
         [10.1, 10.0,  9.7], # (c)
         [11.3, 12.0, 10.3]] # (d)
        # img1  img2  img3

        cmp_snrs = \
        [[100, 250, 115], # (a)
         [120, 280, 110], # (b)
         [150, 400, 170], # (c)
         [140, 310, 160]] # (d)
        # img1  img2  img3

        star_set = self._populate_set(cmp_mags, cmp_snrs)

        #                  (a)  (b)  (c)  (d)
        weights = Weights([0.1, 0.4, 0.3, 0.2])

        # The instrumental magnitude and signal-to-noise ratio of the DBStar
        # whose light curve will be generated in the three images (img1, img2
        # ang img3, respectively)
        star_mags = [12.3, 12.7, 11.1]
        star_snrs = [115, 150, 200]

        # The comparison star is 0.1 x a + 0.4 x b + 0.3 x c + 0.2 x d:
        # img1 = 0.1 * 14.5  + 0.4 * 13.2 + 0.3 * 10.1 + 0.2 * 11.3 = 12.02
        # img2 = 0.1 * 14.2  + 0.4 * 13.1 + 0.3 * 10.0 + 0.2 * 12.0 = 12.06
        # img3 = 0.1 * 13.4  + 0.4 * 12.8 + 0.3 *  9.7 + 0.2 * 10.3 = 11.43

        # To calculate the SNR of the comparison star, we convert from SNR to
        # error in magnitudes, compute the error of the weighted mean of these
        # errors (method snr.mean_snr) and finally convert the resulting error
        # back to the equivalent SNR. The SNR of the differential magnitude is
        # then obtained with the snr.difference_snr method. Please refer to the
        # 'snr' module for further information.
        #
        # cmp img1 = 234.28290196143436 SNR
        # cmp img2 = 560.36980119965619  "
        # cmp img3 = 231.05922831061324  "
        #
        # Finally, calculate the difference between the two magnitudes:
        # img1 = 12.3 - 12.02 = 0.28  (SNR = 103.32946003293335)
        # img2 = 12.7 - 12.06 = 0.64  (SNR = 144.93925214045919)
        # img3 = 11.1 - 11.43 = -0.33 (SNR = 151.36317796943533)

        expected_mags = [0.28, 0.64, -0.33]
        expected_snrs = [103.32946003293335, 144.93925214045919, 151.36317796943533]
        self._assert_light_curve(star_set, weights, star_mags, star_snrs,
                                 expected_mags, expected_snrs)

    def test_light_curve_more_complex_case(self):

        cmp_mags = \
        [[12.1, 13.2, 12.3, 12.4, 11.5], # (a)
         [10.3, 11.2, 10.6, 10.5,  9.7], # (b)
         [ 9.2, 10.5,  9.7,  9.8,  8.6], # (c)
         [11.8, 12.8, 11.9, 11.9, 10.4], # (d)
         [13.5, 14.3, 13.7, 13.8, 13.3]] # (e)
        # img1  img2  img3  img4  img5

        cmp_snrs = \
        [[100,   90,  400,  120,  130], # (a)
         [200,  180,  600,  160,  170], # (b)
         [300,  260,  850,  210,  190], # (c)
         [250,  200,  700,  300,  275], # (d)
         [ 75,  125,  300,  105,  120]] # (e)
        # img1  img2  img3  img4  img5

        star_set = self._populate_set(cmp_mags, cmp_snrs)

        #                  (a)  (b)   (c)   (d)  (e)
        weights = Weights([0.3, 0.15, 0.25, 0.1, 0.2])
        star_mags = [12.3, 14.4, 13.9, 13.4, 11.4]
        star_snrs = [200, 250, 300, 150, 425]

        # The comparison star is 0.3 x a + 0.15 x b + 0.25 x c + 0.1 x d + 0.2 x e:
        # img1 = 0.3 * 12.1 + 0.15 * 10.3 + 0.25 *  9.2 + 0.1 * 11.8 + 0.2 * 13.5 = 11.355
        # img2 = 0.3 * 13.2 + 0.15 * 11.2 + 0.25 * 10.5 + 0.1 * 12.8 + 0.2 * 14.3 = 12.405
        # img3 = 0.3 * 12.3 + 0.15 * 10.6 + 0.25 *  9.7 + 0.1 * 11.9 + 0.2 * 13.7 = 11.635
        # img4 = 0.3 * 12.4 + 0.15 * 10.5 + 0.25 *  9.8 + 0.1 * 11.9 + 0.2 * 13.8 = 11.695
        # img5 = 0.3 * 11.5 + 0.15 *  9.7 + 0.25 *  8.6 + 0.1 * 10.4 + 0.2 * 13.3 = 10.755
        #
        # cmp img1 = 238.05305761543815 (SNR)
        # cmp img2 = 252.92441662287322   "
        # cmp img3 = 921.27421827134822   "
        # cmp img4 = 284.64870682870139   "
        # cmp img5 = 304.45126229797097   "
        #
        # Finally, calculate the difference between the two magnitudes:
        # img1 = 12.3 - 11.355 = 0.945 (SNR = 153.27250023365929)
        # img2 = 14.4 - 12.405 = 1.995 (SNR = 177.94827806386499)
        # img3 = 13.9 - 11.635 = 2.265 (SNR = 285.31216200419124)
        # img4 = 13.4 - 11.695 = 1.705 (SNR = 132.80545638560201)
        # img5 = 11.4 - 10.755 = 0.645 (SNR = 247.63237581648795)

        expected_mags = [0.945, 1.995, 2.265, 1.705, 0.645]
        expected_snrs = [153.27250023365929, 177.94827806386499,
                         285.31216200419124, 132.80545638560201,
                         247.63237581648795]
        self._assert_light_curve(star_set, weights, star_mags, star_snrs,
                                 expected_mags, expected_snrs)

    def test_light_curve_exclude_index(self):

        cmp_mags = \
        [[12.1, 11.4, 10.5], # (a)
         [11.9, 10.3,  9.5], # (b)
         [ 8.2,  9.6,  9.9]] # (c)
        # img1  img2  img3

        cmp_snrs = \
        [[100, 200, 150], # (a)
         [155, 110,  95], # (b)
         [175, 180, 200]] # (c)
        # img1  img2  img3

        star_set = self._populate_set(cmp_mags, cmp_snrs)

        #                  (a)  (b)  (c)
        weights = Weights([0.4, 0.4, 0.2])
        star_mags = [9.5, 8.5, 8.8]
        star_snrs = [200, 225, 300]

        # First case: _exclude_index = 0 (star 'a')
        # Rescaled Weights are [0.66666667, 0.33333333]
        #
        # Comparison star (0.667 x b + 0.333 x c:
        # img1 = 0.667 * 11.9 + 0.333 * 8.2 = 10.666666667 (SNR = 212.41299513302027)
        # img2 = 0.667 * 10.3 + 0.333 * 9.6 = 10.066666667 (SNR = 157.60400333233619)
        # img3 = 0.667 *  9.5 + 0.333 * 9.9 =  9.633333333 (SNR = 138.4335007832733)
        #
        # Finally, calculate the difference between the two magnitudes:
        # img1 = 9.5 - 10.666666667 = -1.1666666669999994 (SNR = 145.75823127087546)
        # img2 = 8.5 - 10.066666667 = -1.5666666669999998 (SNR = 129.2172884922613)
        # img3 = 8.8 -  9.633333333 = -0.8333333329999986 (SNR = 125.78548090920276)

        expected_mags = [-1.1666666669999994, -1.5666666669999998, -0.8333333329999986]
        expected_snrs = [145.75823127087546, 129.2172884922613, 125.78548090920276]
        args = [star_set, weights, star_mags, star_snrs, expected_mags, expected_snrs]
        self._assert_light_curve(*args, _exclude_index = 0)

        # Second case: _exclude_index = 1 (star 'b')
        # Rescaled Weights are (again) [0.66666667, 0.33333333]
        #
        # Comparison star (0.667 x a + 0.333 x c:
        # img1 = 0.667 * 12.1 + 0.333 * 8.2 = 10.8 (SNR = 144.03013882286783)
        # img2 = 0.667 * 11.4 + 0.333 * 9.6 = 10.8 (SNR = 262.07411913056194)
        # img3 = 0.667 * 10.5 + 0.333 * 9.9 = 10.3 (SNR = 210.49309502135651)
        #
        # Finally, calculate the difference between the two magnitudes:
        # img1 = 9.5 - 10.8 = -1.3 (SNR = 117.01037414473119)
        # img2 = 8.5 - 10.8 = -2.3 (SNR = 170.85882319080514)
        # img3 = 8.8 - 10.3 = -1.5 (SNR = 172.44067461856571)

        args[-2] = [-1.3, -2.3, -1.5]
        args[-1] = [117.01037414473119, 170.85882319080514, 172.44067461856571]
        self._assert_light_curve(*args, _exclude_index = 1)

        # Third case: _exclude_index = 2 (star 'c')
        # Rescaled Weights are [0.5, 0.5]
        #
        # Comparison star (0.667 x b + 0.333 x c:
        # img1 = 0.5 * 12.1 + 0.5 * 11.9 = 12.0  (SNR = 167.80583173706879)
        # img2 = 0.5 * 11.4 + 0.5 * 10.3 = 10.85 (SNR = 192.48195714828401)
        # img3 = 0.5 * 10.5 + 0.5 *  9.5 = 10.0  (SNR = 160.25858807928134)
        #
        # Finally, calculate the difference between the two magnitudes:
        # img1 = 9.5 - 12.0  = -2.5  (SNR = 128.69392743319793)
        # img2 = 8.5 - 10.85 = -2.35 (SNR = 146.40717535218673)
        # img3 = 8.8 - 10.0  = -1.2  (SNR = 141.45872827469148)

        args[-2] = [-2.5, -2.35, -1.2]
        args[-1] = [128.69392743319793, 146.40717535218673, 141.45872827469148]
        self._assert_light_curve(*args, _exclude_index = 2)

        # IndexError is raised if _exclude_index is out of range
        nstars = random.randint(*self.NSTARS_RANGE)
        set_, star = self.random_set(nstars)
        weights = Weights.random(nstars)

        args = set_.light_curve, weights, star
        self.assertRaises(IndexError, *args, _exclude_index = -2)
        self.assertRaises(IndexError, *args, _exclude_index = -1)
        self.assertRaises(IndexError, *args, _exclude_index = nstars)
        self.assertRaises(IndexError, *args, _exclude_index = nstars + 1)

        # Nothing goes wrong, of course, for valid indexes
        for index in xrange(nstars):
            set_.light_curve(weights, star, _exclude_index = index)

    def _assert_broeg_weights(self, star_mags, eweights, pct, max_iters):
        """ Assert that the StarSet class returns the expected Broeg weights.

        The method populates a StarSet with a random photometric filter, Unix
        times and signal-to-noise ratios, but the magnitudes of whose stars are
        defined by the values given in the 'star_mags' parameter. It must be a
        sequence of sequences, each one of the latter containing the magnitudes
        of that star in the different images, which is fed into
        StarSetTest._populate_set. Then, the Broeg weights of the stars in the
        set are computed (method StarSet.broeg_weights) and compared to
        'eweights'. If the values are not 'almost' equal, as per the rules of
        the TestCase.assertAlmostEqual method, an unconditional test failure is
        signaled.

        The value of the pct and max_iters parameters is passed down, as
        keyword arguments, to the method StarSet.light_curve.

        Note that there must be as many weights as stars in the StarSet;
        otherwise, an error is raised.

        """

        if len(star_mags) != len(eweights):
            msg = "number of stars in the StarSet must match that of weights"
            self.fail(msg)

        set_ = self._populate_set(star_mags)
        weights = set_.broeg_weights(pct = pct, max_iters = max_iters)
        assertSequencesAlmostEqual(self, weights, eweights)

    def test_broeg_weights_fewer_than_two_stars(self):

        # Zero is also less than two, but that scenario cannot even take place
        # since StarSet.__init__ raises ValueError if it receives no DBStars
        expected_weights = [1]
        for _ in xrange(NITERS):
            set_ = self.random_set(size = 1)[0]
            bweights = set_.broeg_weights()
            self.assertEqual(bweights, expected_weights)

    def test_broeg_weights_two_stars(self):

        expected_weights = [0.5, 0.5]
        for _ in xrange(NITERS):
            set_ = self.random_set(size = 2)[0]
            bweights = set_.broeg_weights()
            self.assertEqual(list(bweights), expected_weights)

    def test_broeg_weights_fewer_than_two_images(self):

        # StarSet.__init__ raises ValueError if it receives empty DBStars,
        # so the number of images of the stars in a set can never be zero
        for index in xrange(NITERS):
            set_ = self.random_set(nrecords = 1)[0]
            self.assertRaises(ValueError, set_.broeg_weights)

    # The small data set used to test StarSet.broeg_weights: the instrumental
    # magnitudes of three stars, observed in three images, whose light curves
    # we have to compute manually to compare them to the output of the method.

    broeg_mags =  \
    [[14.5, 14.2, 13.4], # (a)
     [13.2, 13.1, 12.8], # (b)
     [10.1, 10.0,  9.7]] # (c)
    # img1  img2  img3

    # Flux-propotional weights (initial values in Broeg's algorithm)
    # Weights = [0.295238595044, 0.317072842106, 0.38768856285] (w0)
    #
    # Now, generate the light curves and calculate their standard deviations;
    # for each star, use the other two to compute the artificial comparison
    # star, with the weights previously obtained.
    #
    # --------------------- First iteration ---------------------
    # Star a --> Weights become [0.317072842106, 0.38768856285] =
    #                (rescaled) [0.44990097, 0.55009903]
    # Light curve:
    # img1 = 14.5 - (0.44990097 x 13.2 + 0.55009903 x 10.1) = 3.00530698245
    # img2 = 14.2 - (0.44990097 x 13.1 + 0.55009903 x 10.0) = 2.80530698245
    # img3 = 13.4 - (0.44990097 x 12.8 + 0.55009903 x  9.7) = 2.30530698245
    # Standard deviation: 0.29439202887759480671
    #
    # Star b --> Weights become [0.295238595044, 0.38768856285] =
    #                (rescaled) [0.43231345, 0.56768655]
    # Light curve:
    # img1 = 13.2 - (0.43231345 x 14.5 + 0.56768655 x 10.1) = 1.19782082441
    # img2 = 13.1 - (0.43231345 x 14.2 + 0.56768655 x 10.0) = 1.28428351421
    # img3 = 12.8 - (0.43231345 x 13.4 + 0.56768655 x 9.7 ) = 1.50044023870
    # Standard deviation: 0.12726963336178350794
    #
    # Star c --> Weights become [0.295238595044, 0.317072842106]
    #                (rescaled) [0.48217064, 0.51782936]
    # Light curve:
    # img1 = 10.1 - (0.48217064 x 14.5 + 0.51782936 x 13.2) = -3.72682182672
    # img2 = 10.0 - (0.48217064 x 14.2 + 0.51782936 x 13.1) = -3.63038769953
    # img3 =  9.7 - (0.48217064 x 13.4 + 0.51782936 x 12.8) = -3.38930238156
    # Standard deviation: 0.1419471917795997788
    #
    # Weights inversely proportional to these standard deviations:
    w1 = [0.185628940097, 0.429385068928, 0.384985990976]
    #
    # Absolute percent change between (w0, w1) = 0.3712578801933831773 (pct1)
    #
    # --------------------- Second iteration ---------------------
    # Star a --> Weights become [0.429385068928, 0.384985990976] =
    #                (rescaled) [0.52725973, 0.47274027]
    # Light curve:
    # img1 = 14.5 - (0.52725973 x 13.2 + 0.47274027 x 10.1) = 2.76549482267
    # img2 = 14.2 - (0.52725973 x 13.1 + 0.47274027 x 10.0) = 2.56549482267
    # img3 = 13.4 - (0.52725973 x 12.8 + 0.47274027 x  9.7) = 2.06549482267
    # Standard deviation: 0.29439202887759486837
    #
    # Star b --> Weights become [0.185628940097, 0.384985990976] =
    #                (rescaled) [0.32531385, 0.67468615]
    # Light curve:
    # img1 = 13.2 - (0.32531385 x 14.5 + 0.67468615 x 10.1) = 1.66861905998
    # img2 = 13.1 - (0.32531385 x 14.2 + 0.67468615 x 10.0) = 1.73368182998
    # img3 = 12.8 - (0.32531385 x 13.4 + 0.67468615 x 9.7 ) = 1.89633875499
    # Standard deviation: 0.09576980432459545768
    #
    # Star c --> Weights become [0.185628940097, 0.429385068928] =
    #                (rescaled) [0.3018288, 0.6981712]
    # Light curve:
    # img1 = 10.1 - (0.3018288 x 14.5 + 0.6981712 x 13.2) = -3.49237743951
    # img2 = 10.0 - (0.3018288 x 14.2 + 0.6981712 x 13.1) = -3.43201167959
    # img3 =  9.7 - (0.3018288 x 13.4 + 0.6981712 x 12.8) = -3.28109727978
    # Standard deviation: 0.088855992695495706987
    #
    # Weights inversely proportional to these standard deviations:
    w2 = [0.13537128412, 0.41612517916, 0.44850353672]
    #
    # Absolute percent change between (w1, w2) = 0.27074256823392388096 (pct2)

    def test_broeg_weights_one_iteration(self):

        # StarSet.broeg_mags should, for the above set, stop right after the
        # first iteration if the percent change threshold is set to a value
        # greater than pct1 (0.37), or if the max_iters keyword argument is
        # set to one, or both.

        args = self.broeg_mags, self.w1
        self._assert_broeg_weights(*args, pct = 0.40, max_iters = 1)
        self._assert_broeg_weights(*args, pct = 0.45, max_iters = None)
        self._assert_broeg_weights(*args, pct = None, max_iters = 1)

    def test_broeg_weights_two_iterations(self):

        # For the method to stop after the second iteration, the percent change
        # threshold must be set to a value in the range [pct2, pct1), that is,
        # [0.27, 0.37), or max_iters set to two, or both.

        args = self.broeg_mags, self.w2
        self._assert_broeg_weights(*args, pct = 0.275, max_iters = 2)
        self._assert_broeg_weights(*args, pct = 0.315, max_iters = None)
        self._assert_broeg_weights(*args, pct = None, max_iters = 2)

    def test_worst_fraction_out_of_range(self):

        # # Valid fractions are in the range (0, 1]
        set_, _ = self.random_set()
        self.assertRaises(ValueError, set_.worst, -0.1)
        self.assertRaises(ValueError, set_.worst,  0.0)
        self.assertRaises(ValueError, set_.worst,  1.1)
        self.assertRaises(ValueError, set_.worst,  2.0)

        # Try now with different values, evenly spaced in the interval.
        # The method always returns at least one star, the worst in the
        # set, no matter how low the value of 'fraction' is.
        for fraction in numpy.linspace(1, 0, NITERS, endpoint = False):
            nworst = max(1, round(fraction * len(set_)))
            worst_indexes = set_.worst(fraction)
            self.assertEqual(nworst, len(worst_indexes))

    def test_worst_less_than_three_stars(self):

        # Alternate between StarSets of one and two stars
        for index in xrange(NITERS):
            size = 1 + (index % 2)
            self.assertTrue(1 <= size <= 2)
            set_, _ = self.random_set(size = size)
            fraction = 1.0 - random.random() # (0.0, 1.0]
            self.assertTrue(0 < fraction <= 1)
            self.assertRaises(ValueError, set_.worst, fraction)

    def _assert_worst(self, star_set, fraction, expected_indexes):
        """ Assert that the StarSet correctly identifies the worst stars.

        The method receives a StarSet and calls its method StarSet.worst to
        find the specified fraction of its stars with the highest standard
        deviation in their light curves. If the returned indexes are not
        equal, both in value and order, to those in 'expected_indexes',
        an unconditional test failure is signaled.

        """

        worst_indexes = star_set.worst(fraction)
        self.assertEqual(worst_indexes, expected_indexes)

    # The data set used to test StarSet.worst and best: the instrumental
    # magnitudes of five stars, observed in five images, among which we will
    # manually identify the less and most constant, to compare them to the
    # output of both methods.

    star_mags = \
    [[12.1, 13.2, 12.3, 12.4, 11.5], # index = 0
     [10.3, 11.2, 10.6, 10.5,  9.7], # index = 1
     [ 9.2, 10.5,  9.7,  9.8,  8.6], # index = 2
     [11.8, 12.8, 11.9, 11.9, 10.4], # index = 3
     [13.5, 14.3, 13.7, 13.8, 13.3]] # index = 4
    # img1  img2  img3  img4  img5

    def test_worst(self):

        set_ = self._populate_set(self.star_mags)

        # Use Broeg's algorithm to find the weights of the optimal artificial
        # comparison star for this StarSet and use these (rescaled) values to
        # generate the light curve of each star. The stars, if sorted by the
        # standard deviation of their curves, are the following:
        #
        # Star 3 (stdev = 0.303995955838)
        # Star 4 (  "   = 0.246020144827)
        # Star 2 (  "   = 0.148178935224)
        # Star 1 (  "   = 0.098037225466)
        # Star 0 (  "   = 0.067647973668)
        #
        # Therefore, the method must return [3, 4, 2, 1, 0] for those values of
        # 'fraction' which cause the indexes of all the stars to be returned.
        # For smaller values of the fraction, the returned indexes must be a
        # sublist, always starting with 3, as that is the index of the star
        # with the highest standard deviation.

        indexes = [3, 4, 2, 1, 0]

        # (a) 0.01 x 5 = 0.05, rounded up to one (the minimum)
        # (b) 0.1  x 5 = 0.5,                "
        # (c) 0.25 x 5 = 1.25, rounded down to one
        # (d) 0.29 x 5 = 1.45,            "
        self._assert_worst(set_, 0.01, indexes[0:1])
        self._assert_worst(set_, 0.1,  indexes[0:1])
        self._assert_worst(set_, 0.25, indexes[0:1])
        self._assert_worst(set_, 0.29, indexes[0:1])

        # (e) 0.35 x 5 = 1.75, rounded up to two
        # (f) 0.4  x 5 = 2
        # (g) 0.47 x 5 = 2.35, rounded down to two
        self._assert_worst(set_, 0.35, indexes[0:2])
        self._assert_worst(set_, 0.4,  indexes[0:2])
        self._assert_worst(set_, 0.47, indexes[0:2])

        # (h) 0.5  x 5 = 2.5,  rounded up to three
        # (i) 0.65 x 5 = 3.25, rounded down to three
        self._assert_worst(set_, 0.5,  indexes[0:3])
        self._assert_worst(set_, 0.65, indexes[0:3])

        # (j) 0.8  x 5 = 4
        # (k) 0.89 x 5 = 4.45, rounded down to four
        self._assert_worst(set_, 0.8, indexes[0:4])
        self._assert_worst(set_, 0.8, indexes[0:4])

        # (l) For fractions >= 0.9, must return of the indexes
        self._assert_worst(set_, 0.9,  indexes)
        self._assert_worst(set_, 0.95, indexes)
        self._assert_worst(set_, 0.99, indexes)
        self._assert_worst(set_, 1,    indexes)

        # Test the method now with random data. Not being the right answer
        # known beforehand, we need to determine at runtime which are the
        # indexes that StarSet.worst must return. In order to do that, we
        # compute the Broeg's weights and select the indexes of the stars with
        # the lowest coefficients, as these values are inversely proportional
        # to the standard deviation of the light curves.

        # We could simply use numpy.argsort() here, but we do not fancy the
        # idea of basically replicating the code of StarSet.worst. Although
        # inefficient, at the very least the same functionality should be
        # implemented in a different way when finding the correct answer
        # of random test cases.

        set_ = self.random_set()[0]
        for _ in xrange(NITERS):
            fraction = fraction = 1.0 - random.random() # (0.0, 1.0]
            nworst = max(1, round(fraction * len(set_)))
            self.assertTrue(1 <= nworst <= len(set_))

            # Find, one by one, the indexes of the 'nworst' lowest values
            expected_indexes = []
            bweights = list(set_.broeg_weights())
            while len(expected_indexes) != nworst:
                index = bweights.index(min(bweights))
                expected_indexes.append(index)
                bweights[index] = float('inf')

            worst_indexes = set_.worst(fraction)
            self.assertEqual(worst_indexes, expected_indexes)

    def _assert_best(self, star_set, how_many, fraction, expected_indexes):
        """ Assert that the StarSet correctly identifies the best stars.

        The method receives a StarSet and calls its method StarSet.best to find
        the 'how_many' most constant stars by iteratively discarding the given
        fraction of those with the highest standard deviation in their light
        curves. If the stars in the returned StarSet are not, regardless of
        their order, those whose indexes are given in 'expected_indexes', an
        unconditional test failure is signaled. This method allows us, thus, to
        work with positions instead of with IDs (i.e., 'the first star' vs 'the
        star with ID 345').

        """

        # Find which are the IDs of the stars in the StarSet whose indexes are
        # listed in 'expected_indexes' (e.g., store the ID of the first and
        # third stars in the set) and then, when the best stars have been
        # identified, assert that the IDs match.

        expected_ids = [star_set[index].id for index in expected_indexes]
        best_indexes = star_set.best(how_many, fraction = fraction)
        self.assertEqual(set(best_indexes.star_ids), set(expected_ids))

    def test_best(self):

        set_ = self._populate_set(self.star_mags)

        # In this test case we use a fraction of 0.1 to iteratively discard the
        # less constant stars. As we did in StarSetTest.test_worst, we first
        # use Broeg's algorithm to find the weights of the optimal artificial
        # comparison star for the StarSet and use these (rescaled) values to
        # generate the light curve and standard deviation of each star:
        #
        # Star 0 (stdev = 0.067647973668)
        # Star 1 (  "   = 0.098037225466)
        # Star 2 (  "   = 0.148178935224)
        # Star 3 (  "   = 0.303995955838)
        # Star 4 (  "   = 0.246020144827)
        #
        # First iteration: 0.1 x 5 = 0.5 stars, rounded up to one, are
        # discarded, which means that star 3 (stdev = 0.303995955838) is
        # removed from the set. The remaining four stars have their light
        # curves recomputed, which gives us the following standard deviations:
        #
        # Star 0 (stdev = 0.083712838645)
        # Star 1 (  "   = 0.087727681361)
        # Star 2 (  "   = 0.171624384684)
        # Star 4 (  "   = 0.218442974750)
        #
        # Second iteration: 0.1 x 4 = 0.4 stars, rounded up to the minimum of
        # one, are discarded. This means that star 4 (stdev = 0.218442974750)
        # is removed from the set. Working with the remaining three stars, we
        # have the following standard deviations:
        #
        # Star 0 (stdev = 0.077868981218)
        # Star 1 (  "   = 0.117431394038)
        # Star 2 (  "   = 0.146149524067)
        #
        # The method cannot iterate any further when the minimum of three stars
        # is reached, so these standard deviations must be used to identify the
        # best stars in the set, be it one, two or three. The star 0, that with
        # the lowest standard deviation, is returned first, then follows star 1
        # and the 'less best' one is star 2.

        self._assert_best(set_, 1, 0.1, [0])
        self._assert_best(set_, 2, 0.1, [0, 1])
        self._assert_best(set_, 3, 0.1, [0, 1, 2])

        # Let's try now with a fraction of 0.4. This means that in the first
        # iteration the 0.4 x 5 = 2 worst stars (numbers three and four, as it
        # can be seen from the initial standard deviations we calculated above)
        # are discarded. We are left, again, with stars 0, 1 and 2, among which
        # the best star(s) are directly identifed.

        self._assert_best(set_, 1, 0.4, [0])
        self._assert_best(set_, 2, 0.4, [0, 1])
        self._assert_best(set_, 3, 0.4, [0, 1, 2])

        # StarSet.best always returns at least one star, so at most four stars
        # can be discarded in this test case, independently of the fraction
        # used. We write 'at most' because we may indeed ask the method for
        # more than one star: for example, if the value of the fraction results
        # in four stars being discarded at the first iteration but we want the
        # best three, only two (5 - 2 = 3) are discarded.

        self._assert_best(set_, 1, 0.95, [0])
        self._assert_best(set_, 2, 0.80, [0, 1])
        self._assert_best(set_, 3, 0.90, [0, 1, 2])
        self._assert_best(set_, 4, 0.85, [0, 1, 2, 4])
        self._assert_best(set_, 5, 0.99, [0, 1, 2, 4, 3])

        # A second data set: seven stars, observed in seven images
        other_star_mags = \
        [[11.4, 12.3, 11.5, 12.4, 14.1, 12.6, 11.2], # index = 0
         [12.5, 15.3, 14.8, 14.9, 14.4, 11.4, 12.3], # index = 1
         [10.3, 11.8, 10.6, 12.4, 11.4,  9.2,  8.9], # index = 2
         [14.2, 13.7, 13.8, 12.4, 13.2, 10.2, 12.3], # index = 3
         [12.3, 12.2, 12.6, 12.5, 12.7, 11.3, 10.3], # index = 4
         [ 9.8,  7.3,  8.3,  8.1,  9.3,  9.6, 10.3], # index = 5
         [10.3, 11.4, 10.5, 11.5, 12.9, 11.6, 10.7]] # index = 6
        # img1  img2  img3  img4  img5  img6  img7

        set_ = self._populate_set(other_star_mags)

        # This time we will use a fraction value of 0.4. As always, we start by
        # finding the weights of the optimal articial comparison star for the
        # StarSet, thanks to Broeg's algorithm, and use these (rescaled) values
        # to generate the light curve of each star:
        #
        # Star 0 (stdev = 0.821609782333)
        # Star 1 (  "   = 1.091553614770)
        # Star 2 (  "   = 0.803179601172)
        # Star 3 (  "   = 1.217139839440)
        # Star 4 (  "   = 0.514281645173)
        # Star 5 (  "   = 1.623944041070)
        # Star 6 (  "   = 0.816862488012)
        #
        # First iteration: 0.4 x 7 = 2.8 stars, rounded up to three, are
        # discarded: numbers 1, 3 and 5, those with the highest standard
        # deviations. After their removal, the light curves of the remaining
        # four stars are recomputed, which gives us the following standard
        # deviations:

        # Star 0 (stdev = 0.578335746321)
        # Star 2 (  "   = 0.968634900697)
        # Star 4 (  "   = 0.762522738044)
        # Star 6 (  "   = 0.609073938371)

        self._assert_best(set_, 4, 0.4, [0, 2, 4, 6])
        self._assert_best(set_, 5, 0.4, [0, 2, 4, 6, 1])
        self._assert_best(set_, 6, 0.4, [0, 2, 4, 6, 1, 3])
        self._assert_best(set_, 7, 0.4, [0, 2, 4, 6, 1, 3, 5])

        # Second iteration: 0.4 x 4 = 1.6 stars, rounded up to two, should be
        # discarded: numbers 2 and 4. However, three is the minimum number of
        # stars that there must be in a StarSet for it be able to establish the
        # variability of its stars. Therefore, only one star (number 2, as it
        # has the highest standard deviation) will be discarded; once the
        # minimum of three stars is reached, the light curves of the best
        # three stars in the set are recomputed again:
        #
        # Star 0 (stdev = 0.259713715458)
        # Star 4 (  "   = 0.952603228565)
        # Star 6 (  "   = 0.305675915384)

        self._assert_best(set_, 1, 0.4, [0])
        self._assert_best(set_, 2, 0.4, [0, 6])
        self._assert_best(set_, 3, 0.4, [0, 6, 4])

        # Use the same data set, but now with a fraction value of 0.8.
        # As a reminder, the initial standard deviations of the stars are:
        #
        # Star 0 (stdev = 0.821609782333)
        # Star 1 (  "   = 1.091553614770)
        # Star 2 (  "   = 0.803179601172)
        # Star 3 (  "   = 1.217139839440)
        # Star 4 (  "   = 0.514281645173)
        # Star 5 (  "   = 1.623944041070)
        # Star 6 (  "   = 0.816862488012)
        #
        # First iteration: 0.8 x 7 = 5.6 stars, rounded up to six, should be
        # discarded. However, the StarSet must have a minimum of three stars,
        # so only 7 - 3 = 4 stars can be removed: numbers 0, 1, 3 and 5, the
        # four with the highest standard deviations. After they are discarded,
        # the light curves of the best three stars in the set are recomputed
        # again, giving us the following standard deviations:
        #
        # Star 2 (  "   = 0.800223863178)
        # Star 4 (  "   = 0.631528873069)
        # Star 6 (  "   = 1.004307069400)
        #
        # Note how the best star (number 4) identified with this fraction, 0.8,
        # differs from that (number 0) found earlier with a smaller fraction,
        # 0.4. This is an excellent example of how the value of the 'fraction'
        # keyword parameter determines not only the performance of StarSet.best
        # but also its reliability.

        self._assert_best(set_, 1, 0.8, [4])
        self._assert_best(set_, 2, 0.8, [4, 2])
        self._assert_best(set_, 3, 0.8, [4, 2, 6])

        self._assert_best(set_, 4, 0.8, [4, 2, 6, 0])
        self._assert_best(set_, 5, 0.8, [4, 2, 6, 0, 1])
        self._assert_best(set_, 6, 0.8, [4, 2, 6, 0, 1, 3])
        self._assert_best(set_, 7, 0.8, [4, 2, 6, 0, 1, 3, 5])

    def test_best_fraction_out_of_range(self):

        # Valid fractions are in the range (0, 1]
        set_ = self.random_set()[0]
        for _ in xrange(NITERS):
            fraction = 1.0 - random.random() # (0.0, 1.0]
            self.assertTrue(0 < fraction <= 1)
            set_.best(1, fraction = fraction)

        self.assertRaises(ValueError, set_.best, 1, fraction = -0.1)
        self.assertRaises(ValueError, set_.best, 2, fraction =  0.0)
        self.assertRaises(ValueError, set_.best, 3, fraction =  1.1)
        self.assertRaises(ValueError, set_.best, 4, fraction =  2.0)

    def test_best_less_than_three_stars(self):

        # Alternate between StarSets of one and two stars
        for index in xrange(NITERS):
            size = 1 + (index % 2)
            self.assertTrue(1 <= size <= 2)
            set_ = self.random_set(size = size)[0]
            fraction = 1.0 - random.random() # (0.0, 1.0]
            self.assertTrue(0 < fraction <= 1)
            self.assertRaises(ValueError, set_.best, fraction)

    def test_best_number_of_stars_out_of_range(self):

        set_ = self.random_set()[0]
        for _ in xrange(NITERS):
            how_many = random.randint(1, len(set_))
            set_.best(how_many)

        self.assertRaises(ValueError, set_.best, -5)
        self.assertRaises(ValueError, set_.best, -1)
        self.assertRaises(ValueError, set_.best, 0)
        self.assertRaises(ValueError, set_.best, len(set_) + 1)
        self.assertRaises(ValueError, set_.best, len(set_) + 5)

