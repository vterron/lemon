#! /usr/bin/env python2

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

import random
import numpy
import uncertainties

from test import unittest
from snr import *

NITERS  = 1000   # How many times each test case is run with random data
MIN_SNR = 2      # Minimum value for random signal-to-noise ratios
MAX_SNR = 10000  # Maximum value for random signal-to-noise ratios

MIN_MAG = 1.47   # Minimum value for random magnitudes (Saturn)
MAX_MAG = 25     # Maximum value for random magnitudes (Fenrir)

MIN_NERR = 1     # Minimum number of random errors when combining them
MAX_NERR = 10    # Maximum number of random errors when combining them

NMEANS = 25 # The number of times each mean_* method is run with random data
MIN_WEIGHT = 0.0001  # Minimum value for random coefficients (weighted mean)
MAX_WEIGHT = 99.999  # Maximum value for random coefficients (weighted mean)

class SNRTest(unittest.TestCase):

    def setUp(self):

        # A list of three-element tuples: the signal-to-noise ratio and the
        # positive and negative errors induced by noise. These values have been
        # manually (metaforically speaking, there was a computer involved)
        # calculated, so they can be safely used.
        self.known_snrs = \
            [(3, -0.31234684152074982, 0.44022814763920293),
             (10, -0.10348171289556268, 0.1143937264016878),
             (30, -0.035601097786525711, 0.036808142051765864),
             (500, -0.0021693038280672832, 0.0021736467815722288),
             (1000, -0.0010851936982964822, 0.0010862794350442293)]

        self.random_snrs = \
            (random.uniform(MIN_SNR, MAX_SNR) for x in xrange(NITERS))

    def test_snr_to_error(self):
        for snr, max_error, min_error in self.known_snrs:
            self.assertAlmostEqual(snr_to_error(snr)[0], max_error)
            self.assertAlmostEqual(snr_to_error(snr)[1], min_error)

        # ValueErrror raised id SNR outside of the function's domain
        self.assertRaises(ValueError, snr_to_error, 1)
        self.assertRaises(ValueError, snr_to_error, random.random())
        self.assertRaises(ValueError, snr_to_error, 0)
        self.assertRaises(ValueError, snr_to_error, -1.5)

    def test_error_to_snr(self):
        for snr, max_error, min_error in self.known_snrs:
            self.assertAlmostEqual(error_to_snr(max_error), snr)
            self.assertAlmostEqual(error_to_snr(min_error), snr)

        # Now generate some random signal-to-noise ratio, and check that the
        # conversion to error and then back to SNR gives the original value.
        for snr in self.random_snrs:
            max_error, min_error = snr_to_error(snr)
            self.assertAlmostEqual(error_to_snr(max_error), snr)
            self.assertAlmostEqual(error_to_snr(min_error), snr)


    def _random_magnitude_error(self):
        """ Return a random magnitude with the corresponding error """
        magnitude = random.uniform(MIN_MAG, MAX_MAG)
        snr = random.uniform(MIN_SNR, MAX_SNR)
        merror = snr_to_error(snr)[1] # the positive error
        return magnitude, merror

    def _random_error(self):
        """ Return a random error in magnitudes """
        return self._random_magnitude_error()[1]

    def test_difference_error(self):

        error = 0.05   # The method also accepts a single value
        self.assertAlmostEqual(difference_error(error), error)

        errors = (0.01, 0.01)
        cerror = 0.01414213562373095
        self.assertAlmostEqual(difference_error(*errors), cerror)

        errors = (0.1, 0.05)
        cerror = 0.1118033988749895
        self.assertAlmostEqual(difference_error(*errors), cerror)

        errors = (0.78, 0.56, 0.21)
        cerror = 0.9829038610159185
        self.assertAlmostEqual(difference_error(*errors), cerror)

        # Use some random magnitudes with uncertainty (two-element tuples,
        # called below 'stars'), using the uncertainties module, and make sure
        # that the uncertainty that results from adding or subtracting them
        # matches the error returned by our method.

        for _ in xrange(NITERS):
            how_many = random.randint(MIN_NERR, MAX_NERR)
            stars = [self._random_magnitude_error() for _ in xrange(how_many)]

            # First compute the resulting error returned by our method. We only
            # need the errors in magnitudes (the second element of each tuple)
            our_err = difference_error(*[x[1] for x in stars])

            # And now calculate the same value using somebody else's code: the
            # 'uncertainties' module, that transparently handles calculations
            # with numbers with uncertainties. We could just add all of them,
            # but we are taking it a step further -- the values will be
            # randomly added or subtracted, as snr.diference_error supports
            # subtraction of errors, addition or a combination of them.

            ufloats = [uncertainties.ufloat(*x) for x in stars]

            total = 0
            for number in ufloats:
                if random.choice((True, False)):
                    total += number
                else:
                    total -= number

            self.assertAlmostEqual(total.std_dev, our_err)


    def test_difference_snr(self):

        snr = 115  # The method also accepts a single value
        self.assertAlmostEqual(difference_snr(snr), snr)

        snrs = (100, 100)
        csnr = 70.8577169468
        self.assertAlmostEqual(difference_snr(*snrs), csnr)

        snrs = (100, 200)
        csnr = 89.5403909525
        self.assertAlmostEqual(difference_snr(*snrs), csnr)

        snrs = (90, 170, 350)
        csnr = 77.6904354852
        self.assertAlmostEqual(difference_snr(*snrs), csnr)

        # Now with random data: generate a series of errors in magnitudes and
        # obtain the absolute error that results from calculating their
        # difference. Then, convert them to signal-to-noise ratios to do the
        # analog operation: calculate the SNR that results from calculating
        # their difference and check that it is equal to the value we obtained
        # earlier when working with errors.

        for _ in xrange(NITERS):
            how_many = random.randint(MIN_NERR, MAX_NERR)
            errors = [self._random_error() for i in xrange(how_many)]
            cerror = difference_error(*errors)

            snrs = (error_to_snr(e) for e in errors)
            csnr = difference_snr(*snrs)

            # We compare the second of the two errors (i.e., the positive one)
            # returned by the conversion from signal-to-noise ratio, as this is
            # the one with which the conversion from error to SNR works -- as
            # the addition in quadrature would make the negative error become
            # positive, and as such they would be considered when converting
            # back to SNR.

            back_to_error = snr_to_error(csnr)[1]
            self.assertAlmostEqual(back_to_error, cerror)


    def test_unweighted_mean_error(self):

        error = 0.05   # The method also accepts a single value
        self.assertAlmostEqual(mean_error([error]), error)

        errors = (0.1, 0.05)
        cerror = 0.05590169943749475
        self.assertAlmostEqual(mean_error(errors), cerror)

        errors = (0.1, 0.07, 0.01)
        cerror = 0.040824829046386298
        self.assertAlmostEqual(mean_error(errors), cerror)

        # If we combine a series of observations with the same SNR ratio, we
        # know that the SNR will increase with the square root of the number of
        # observations that were averaged. Although we are working here with
        # errors, we can put this statement to test too: if the average of n
        # times the same error is computed, it should decrease with sqrt(n)
        for _ in xrange(NITERS):
            snr = random.uniform(MIN_SNR, MAX_SNR)
            error = snr_to_error(snr)[1] # the positive error
            for n in xrange(2, MAX_NERR):
                merror = mean_error([error] * n)
                expected_error = error / math.sqrt(n)
                self.assertAlmostEqual(merror, expected_error)

        # Finally, test the unweighted mean using once more the 'uncertainties'
        # module. We generate random numbers with uncertainties (errors) and
        # check that the uncertainty that this package returns for their
        # arithmetic mean matches our value.

        for _ in xrange(NITERS):
            how_many = random.randint(MIN_NERR, MAX_NERR)
            stars = [self._random_magnitude_error() for _ in xrange(how_many)]

            # First compute the resulting error returned by our method. We only
            # need the errors in magnitudes (the second element of each tuple)
            our_err = mean_error([x[1] for x in stars])

            # Now let NumPy calculate the average of the 'ufloat' instances
            # with wich 'uncertainties' represents numbers and their error.
            ufloats = numpy.array([uncertainties.ufloat(*x) for x in stars])
            self.assertAlmostEqual(ufloats.mean().std_dev, our_err)

    def _random_weight(self):
        """ Return a random float in the range [MIN_WEIGHT, MAX_WEIGHT] """
        return random.uniform(MIN_WEIGHT, MAX_WEIGHT)

    def test_weighted_mean_error(self):

        error = 0.05   # The method also accepts a single value
        self.assertAlmostEqual(mean_error([error], weights = [1]), error)

        errors = (0.1, 0.05)
        weights1 = (0.75, 0.25)
        weights2 = (0.56, 0.44)
        cerror1 = 0.076034531628727753
        cerror2 = 0.060166435825965307
        self.assertAlmostEqual(mean_error(errors, weights = weights1), cerror1)
        self.assertAlmostEqual(mean_error(errors, weights = weights2), cerror2)

        errors = (0.1, 0.07, 0.01)
        weights1 = (0.45, 0.33, 0.22)
        weights2 = (0.28, 0.68, 0.04)
        weights3 = (0.28, 0.21, 0.51)
        weights4 = (0.45, 0.55, 0)
        cerror1 = 0.050630524389936947
        cerror2 = 0.055226080795218492
        cerror3 = 0.03203279569441294
        cerror4 = 0.059222039816271117

        self.assertAlmostEqual(mean_error(errors, weights = weights1), cerror1)
        self.assertAlmostEqual(mean_error(errors, weights = weights2), cerror2)
        self.assertAlmostEqual(mean_error(errors, weights = weights3), cerror3)
        self.assertAlmostEqual(mean_error(errors, weights = weights4), cerror4)

        errors = (0.1, 0.25)
        weights = (0.25, 0.50, 0.35)
        with self.assertRaises(ValueError):
            mean_error(errors, weights = weights)

        # It does not matter if the weights do not sum up to one, as they will
        # be internally normalized by the method. This means that the weights
        # [1, 2, 3], for example, are equivalent to [2, 4, 6] or [0.5, 1, 1.5]

        errors = (0.1, 0.06, 0.008)
        weights1 = (1, 2, 3)
        weights2 = (2, 4, 6)
        weights3 = (0.5, 1, 1.5)
        cerror = 0.026339661686851215

        self.assertAlmostEqual(mean_error(errors, weights = weights1), cerror)
        self.assertAlmostEqual(mean_error(errors, weights = weights2), cerror)
        self.assertAlmostEqual(mean_error(errors, weights = weights3), cerror)

        # Now with random data: generate a series of random errors and their
        # corresponding weights and see whether, if multiplied by the same
        # constant, the resulting error stays the same -- it must!
        for _ in xrange(NITERS):
            how_many = random.randint(MIN_NERR, MAX_NERR)
            errors = [self._random_error() for i in xrange(how_many)]
            rweights = [self._random_weight() for i in xrange(how_many)]
            assert len(errors) == len(rweights)
            expected_err = mean_error(errors, weights = rweights)

            for n in xrange(NMEANS):
                constant = random.uniform(MIN_WEIGHT, MAX_WEIGHT)
                cweights = [constant * w for w in rweights]
                returned_err = mean_error(errors, weights = cweights)
                self.assertAlmostEqual(returned_err, expected_err)

        # Using the same weight for all the errors should yield the same
        # result as not using weights at all (the unweighted mean)

        errors = (0.01, 0.08)
        expected_err = mean_error(errors)

        weights1 = (0.5, 0.5)
        weights2 = (1, 1)
        weights3 = (math.pi, math.pi)
        self.assertEqual(mean_error(errors, weights = weights1), expected_err)
        self.assertEqual(mean_error(errors, weights = weights2), expected_err)
        self.assertEqual(mean_error(errors, weights = weights2), expected_err)

        # With random data, one more time: check that, no matter what the
        # values of our errors are, if the coefficients of the weighted mean
        # are all equal the resulting error result is equal to that of the
        # normal, unweighted mean.
        for _ in xrange(NITERS):
            how_many = random.randint(MIN_NERR, MAX_NERR)
            errors = [self._random_error() for _ in xrange(how_many)]
            expected_err = mean_error(errors)

            for n in xrange(NMEANS):
                weights = [random.uniform(MIN_WEIGHT, MAX_WEIGHT)] * len(errors)
                assert len(set(weights)) == 1
                returned_err = mean_error(errors, weights = weights)
                self.assertAlmostEqual(returned_err, expected_err)

        # To be added here in the future: more tests cases making use of the
        # 'uncertainties' package, once NumPy arrays of ufloats can be passed
        # to numpy.average. [As of Fri Mar 23 2012, it raises "AttributeError:
        # 'AffineScalarFunc' object has no attribute 'dtype'"]


    # The unit tests of mean_snr are mostly expressed in terms of those of
    # mean_error: first, a series of random error in magnitudes are generated
    # and the error of their arithmetic mean is calculated. Then, we convert
    # the errors to signal-to-noise ratios and check that the value returned by
    # mean_snr corresponds to the error obtained in the first step. If we are
    # sure that both operations are equivalent there is no need to test deeper,
    # as mean_snr is simply a wrapper around mean_error to handle to error to
    # SNR to error conversion.

    def test_unweighted_mean_snr(self):

        snr = 100   # The method also accepts a single value
        self.assertAlmostEqual(mean_snr([snr]), snr)

        snrs = (100, 200)
        cerror = 178.57937803657543
        self.assertAlmostEqual(mean_snr(snrs), cerror)

        snrs = (100, 150, 200)
        cerror = 229.99838960854942
        self.assertAlmostEqual(mean_snr(snrs), cerror)

        for _ in xrange(NITERS):
            how_many = random.randint(MIN_NERR, MAX_NERR)
            errors = [self._random_error() for i in xrange(how_many)]
            cerror = mean_error(errors)

            snrs = (error_to_snr(e) for e in errors)
            csnr = mean_snr(snrs)

            back_to_error = snr_to_error(csnr)[1]
            self.assertAlmostEqual(back_to_error, cerror)

    def test_weighted_mean_snr(self):

        snrs = (100, 200)
        weights = (0.25, 0.50, 0.35)
        with self.assertRaises(ValueError):
            mean_snr(snrs, weights = weights)

        for _ in xrange(NITERS):
            how_many = random.randint(MIN_NERR, MAX_NERR)
            errors = [self._random_error() for i in xrange(how_many)]
            weights = [self._random_weight() for i in xrange(how_many)]
            assert len(errors) == len(weights)
            cerror = mean_error(errors, weights = weights)

            snrs = (error_to_snr(e) for e in errors)
            csnr = mean_snr(snrs, weights = weights)
            back_to_error = snr_to_error(csnr)[1]
            self.assertAlmostEqual(back_to_error, cerror)

