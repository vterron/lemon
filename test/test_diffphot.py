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
import random
import unittest

from diffphot import Weights

NITERS = 100  # How many times some test cases are run with random data

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

