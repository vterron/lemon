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

import functools
import os.path
import string
import unittest

from passband import Passband, NonRecognizedPassband, InvalidPassbandLetter

NITERS  = 100     # How many times each test case is run with random data
NPASSBANDS = 100  # Number of elements for sequences of random Passbands

class PassbandTest(unittest.TestCase):

    get_data_path = functools.partial(os.path.join, './test/test_data/filters')
    JOHNSON_TEST_DATA = get_data_path('Johnson')
    COUSINS_TEST_DATA = get_data_path('Cousins')

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

    @staticmethod
    def read_filter_data_file(path):
        """ Read the contents of a file in the ./test_data/filters directory.

        Loop over the lines of the file, each one of them expected to contain a
        two-element tuple which is eval()'uated and yielded. The first element
        of the tuple should be the name of the photometric filter (for example,
        'Johnson V'), while the second should be the letter that the Passband
        class must identify (e.g., 'V'). Empty and comment lines are ignored.

        """

        with open(path, 'r') as fd:

            for line in fd:
                line = line.strip()

                # Ignore empty and comment lines
                if not line or line[0] == '#':
                    continue

                name, letter = eval(line)
                yield name, letter

    def test_johnson_filters(self):

        # Read the file with test data for the Johnson photometric system,
        # containing two-element tuples such as "Johnson U", 'U'. For each one
        # of them, the first element is used to instantiate a Passband object
        # (e.g., `Passband("Johnson U")`), and we make sure that the system is
        # identified as 'Johnson' and the letter (e.g. 'U') correctly parsed.

        for name, letter in self.read_filter_data_file(self.JOHNSON_TEST_DATA):
                passband = Passband(name)
                self.assertEqual(passband.system, 'Johnson')
                self.assertEqual(passband.letter, letter)

        # Letters other than UBVRIJHKLMN raise InvalidPassbandLetter
        for letter in string.ascii_uppercase:
            if letter not in Passband.JOHNSON_LETTERS:
                name = "Johnson %s" % letter
                self.assertRaises(InvalidPassbandLetter, Passband, name)

        # There are some rare cases in which the letter of the photometric
        # filter cannot be identified and NonRecognizedPassband is raised.
        # For example, if more than one letter is given.

        for name in 'Johnson', 'Johnson BV', 'BV (Johnson)', 'RI_(John)':
            self.assertRaises(NonRecognizedPassband, Passband, name)

    def test_cousins_filters(self):

        for name, letter in self.read_filter_data_file(self.COUSINS_TEST_DATA):
                passband = Passband(name)
                self.assertEqual(passband.system, 'Cousins')
                self.assertEqual(passband.letter, letter)

        # Letters other than RI raise InvalidPassbandLetter
        for letter in string.ascii_uppercase:
            if letter not in Passband.COUSINS_LETTERS:
                name = "Cousins %s" % letter
                self.assertRaises(InvalidPassbandLetter, Passband, name)

        # NonRecognizedPassband raised if the letter of the filter cannot be
        # identified (for example, if more than one letter is given).
        for name in 'Cousins', 'Cousins RI', 'RI (Cousins)', 'RI_(Cou)':
            self.assertRaises(NonRecognizedPassband, Passband, name)

