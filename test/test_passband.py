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
import string
import unittest

from passband import Passband, NonRecognizedPassband, UnknownPassbandLetter

NITERS  = 100     # How many times each test case is run with random data
NPASSBANDS = 100  # Number of elements for sequences of random Passbands

class PassbandTest(unittest.TestCase):

    def test_init(self):
        # Make sure that the constructor works as expected.

        # Improperly-formatted filter names are expected to be rejected
        self.assertRaises(NonRecognizedPassband, Passband, "V(Johnson)")
        self.assertRaises(NonRecognizedPassband, Passband, "Johnson (V)")
        self.assertRaises(NonRecognizedPassband, Passband, "Johnson(V)")
        self.assertRaises(NonRecognizedPassband, Passband, "Johnson")
        self.assertRaises(NonRecognizedPassband, Passband, " ")
        self.assertRaises(NonRecognizedPassband, Passband, '')

        for letter in string.ascii_uppercase:

            # The letter should be correctly extracted from the filter name...
            if letter in Passband.wavelengths.keys():
                self.assertEqual(Passband(letter).letter, letter)
                self.assertEqual(Passband("%s (Johnson)" % letter).letter, letter)
                self.assertEqual(Passband("Johnson %s" % letter).letter, letter)

            else: # ... unless it does not belong to the photometric system
                self.assertRaises(UnknownPassbandLetter, Passband, letter)
                self.assertRaises(UnknownPassbandLetter, Passband, "%s Johnson" % letter)
                self.assertRaises(UnknownPassbandLetter, Passband, "Johnson %s" % letter)

    def test_all(self):
        # Make sure all the photometric letters are present in all()
        wavelengths = set([x.wavelength for x in Passband.all()])
        self.assertEqual(wavelengths, set(Passband.wavelengths.values()))

    def test_wavelength(self):
        for letter in Passband.wavelengths.keys():
            expected_wavelength = Passband.wavelengths[letter]
            self.assertEqual(Passband(letter).wavelength, expected_wavelength)

    def test_repr(self):
        for letter in Passband.wavelengths.keys():
            self.assertEqual(Passband(letter), eval(`Passband(letter)`))

    def test_cmp(self):
        # Make sure that filters are correctly sorted by their wavelength.

        # A specific test case: B (445 nm) < V (551 nm) < I (806 nm)
        self.assertEqual(min(Passband('V'), Passband('B')), Passband('B'))
        self.assertEqual(min(Passband('B'), Passband('V')), Passband('B'))
        self.assertEqual(min(Passband('V'), Passband('I')), Passband('V'))
        self.assertEqual(min(Passband('I'), Passband('V')), Passband('V'))
        self.assertEqual(min(Passband('B'), Passband('I')), Passband('B'))
        self.assertEqual(min(Passband('I'), Passband('B')), Passband('B'))

        # Now sort a series of lists of random passbands
        for _ in xrange(NITERS):
            passbands = [Passband.random() for x in xrange(NPASSBANDS)]
            passbands.sort()
            for index in xrange(0, len(passbands) - 1):
                self.assertTrue(passbands[index].wavelength <= passbands[index+1].wavelength)

    def test_hash(self):
        # The hash must be its effective wavelength midpoint
        for _ in xrange(NITERS):
            passband = Passband.random()
            self.assertEqual(hash(passband), Passband.wavelengths[passband.letter])

    def test_random(self):
        # Make sure the returned filter is a valid one
        for _ in xrange(NITERS):
            passband = Passband.random()
            self.assertTrue(passband.letter in Passband.wavelengths.iterkeys())

    def test_different(self):
        for _ in xrange(NITERS):
            passband = Passband.random()
            self.assertNotEqual(passband, passband.different())

