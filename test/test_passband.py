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

from __future__ import division

import functools
import os.path
import random
import string

# LEMON modules
import methods
import passband
from test import unittest
from passband import Passband, NonRecognizedPassband, InvalidPassbandLetter, \
                     JOHNSON, COUSINS, HARRIS, GUNN, SDSS, TWOMASS, STROMGREN, \
                     HALPHA, UNKNOWN, CUSTOM

NITERS  = 100     # How many times each test case is run with random data
NPASSBANDS = 100  # Number of elements for sequences of random Passbands

class PassbandTest(unittest.TestCase):

    TEST_DATA_DIR = './test/test_data/filters'

    @classmethod
    def get_data_path(cls, system):
        """ Return the path to the test data file for a photometric system.

        The data files are located in the TEST_DATA_DIR directory, and their
        filename matches the name of the photometric system (e.g., 'Johnson').
        These files contain a series of two-element tuples, such as ('Johnson
        U', 'U'), one per line.

        """

        return os.path.join(cls.TEST_DATA_DIR, system)

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

    def _test_photometric_system(self, system):
        """ Test that Passband parses a photometric system correctly.

        'system' must be the formal name of the photometric system, adequately
        capitalized, such as 'Johnson' or 'Cousins'. The Passband class must
        set the 'system' attribute to this value at instantiation time.

        """

        data_file = self.get_data_path(system)

        # For each two-element tuple, such as 'Johnson U', 'U', the first
        # element is used to instantiate a Passband object and we then make
        # sure that both the system is the letter are correctly parsed.
        for name, letter in self.read_filter_data_file(data_file):
            passband = Passband(name)
            self.assertEqual(passband.system, system)
            self.assertEqual(passband.letter, letter)

        # A sequence of the letters allowed by the photometric system.
        # In the case of Johnson, for example, they are 'UBVRIJHKLMN'.
        valid_letters = Passband.SYSTEM_LETTERS[system]

        # Letters other than the valid ones raise InvalidPassbandLetter
        for letter in string.ascii_uppercase:
            if letter not in valid_letters:
                name = "%s %s" % (system, letter)
                with self.assertRaises(InvalidPassbandLetter):
                    Passband(name)

        # There are some rare cases in which the letter of the photometric
        # filter cannot be identified and NonRecognizedPassband is raised.
        # For example, if more than one letter is given.

        two_letters = ''.join(random.sample(valid_letters, 2))
        values = dict(system = system, letter = two_letters)
        patterns = ("%(system)s",              # e.g., "Johnson"
                    "%(letter)s %(system)s",   # e.g., "BV Johnson"
                    "%(system)s %(letter)s",   # e.g., "Johnson BV"
                    "%(letter)s (%(system)s)") # e.g., "BV (Johnson)"

        for pattern in patterns:
            name = pattern % values
            with self.assertRaises(NonRecognizedPassband):
                Passband(name)

    def test_johnson_filters(self):
        self._test_photometric_system(JOHNSON)

    def test_cousins_filters(self):
        self._test_photometric_system(COUSINS)

    def test_harris(self):
        self._test_photometric_system(HARRIS)

    def test_gunn_filters(self):
        self._test_photometric_system(GUNN)

    def test_sdss_filters(self):
        self._test_photometric_system(SDSS)

    def test_2mass_filters(self):
        self._test_photometric_system(TWOMASS)

    def test_stromgren_filters(self):
        self._test_photometric_system(STROMGREN)

    def test_halpha_filters(self):
        data_file = self.get_data_path(HALPHA)
        for name, wavelength in self.read_filter_data_file(data_file):
            passband = Passband(name)
            self.assertEqual(passband.system, HALPHA)
            self.assertEqual(passband.letter, wavelength)

    def test_custom_filters(self):
        for name, description in passband.CUSTOM_FILTERS.iteritems():
            pfilter = Passband(name)
            self.assertEqual(pfilter.system, CUSTOM)
            self.assertEqual(pfilter.letter, description)

    def test_unknown_filters(self):
        data_file = self.get_data_path(UNKNOWN)
        for name, letter in self.read_filter_data_file(data_file):
            passband = Passband(name)
            self.assertEqual(passband.system, UNKNOWN)
            self.assertEqual(passband.letter, letter)

        for letter in string.ascii_uppercase:
            if letter not in passband.ALL_LETTERS:
                with self.assertRaises(NonRecognizedPassband):
                    Passband(letter)

    def test_all(self):

        # Make sure that, except for H-alpha, all the photometric systems
        # and letters are present in the list returned by Passband.all()
        pfilters = Passband.all()
        # There must not be duplicate Passband objects
        self.assertEqual(len(pfilters), len(set(pfilters)))
        for system, letters in Passband.SYSTEM_LETTERS.iteritems():
            for letter in letters:
                name = "%s %s" % (system, letter)
                pfilter = Passband(name)
                self.assertIn(pfilter, pfilters)

        # No user-defined filter must be missing either
        for name in passband.CUSTOM_FILTERS.iterkeys():
            pfilter = Passband(name)
            self.assertIn(pfilter, pfilters)

    def test_repr(self):
        for _ in xrange(NITERS):
            pfilter = Passband.random()
            self.assertEqual(pfilter, eval(`pfilter`))

    def test_cmp(self):

        def assert_max(index, *names):
            """ Make sure that Passband(names[index]) is the largest item """
            pfilters = [Passband(x) for x in names]
            self.assertEqual(max(pfilters), pfilters[index])

        assert_max(1, "Johnson B", "Cousins R")
        assert_max(0, "SDSS Z", "Gunn U")
        assert_max(1, "Johnson V", "SDSS I")
        assert_max(0, "2MASS KS", "Johnson R")
        assert_max(1, "Cousins I", "2MASS J")
        assert_max(0, "Gunn G", "Stromgren B")

        def letter_index(letter):
            """ Return the position of a letter in Passband.LETTERS_ORDER """
            return Passband.LETTERS_ORDER.index(letter)

        # Now sort a series of lists of random Passband objects
        for _ in xrange(NITERS):

            pfilters = [Passband.random() for x in xrange(NPASSBANDS)]
            pfilters.sort()

            for index in xrange(0, len(pfilters) - 1):
                first  = pfilters[index]
                second = pfilters[index + 1]
                ncustoms = sum(1 for p in [first, second] if p.system == CUSTOM)
                nhalphas = sum(1 for p in [first, second] if p.system == HALPHA)

                if not nhalphas and not ncustoms:
                    first_index  = letter_index(first.letter)
                    second_index = letter_index(second.letter)
                    self.assertTrue(first_index <= second_index)
                    # If letters are equal, sort lexicographically by the system
                    if first.letter == second.letter:
                        self.assertTrue(first.system <= second.system)

                # Custom filters are smaller than others system
                elif ncustoms == 1:
                    first.system  == CUSTOM
                    second.system != CUSTOM

                elif ncustoms == 2:
                    # Two custom filters are compared lexicographically
                    self.assertTrue(first.letter <= second.letter)

                # H-alpha filters are greater than other systems
                elif nhalphas == 1:
                    first.system  != HALPHA
                    second.system == HALPHA

                else:
                    assert nhalphas == 2
                    self.assertTrue(int(first.letter) <= int(second.letter))

    def test_hash(self):
        for _ in xrange(NITERS):
            pfilter = Passband.random()
            self.assertEqual(hash(pfilter), hash(eval(`pfilter`)))
            self.assertNotEqual(hash(pfilter), hash(pfilter.different()))

    def test_random(self):
        # Make sure the returned filter is a valid one
        for _ in xrange(NITERS):
            pfilter = Passband.random()
            self.assertTrue(pfilter.system in Passband.ALL_SYSTEMS)
            # Neither custom nor H-alpha filters have letter
            if pfilter.system not in [CUSTOM, HALPHA]:
                self.assertTrue(pfilter.letter in Passband.SYSTEM_LETTERS[pfilter.system])

    def test_different(self):
        for _ in xrange(NITERS):
            pfilter = Passband.random()
            self.assertNotEqual(pfilter, pfilter.different())

    def test_load_custom_filters(self):

        section_header = "[%s]" % passband.CUSTOM_SECTION

        # Name and description of user-defined filters
        custom1 = "R-EROS", "R (EROS-2 survey)"
        custom2 = "NO", "Blank Filter"
        custom_filters = [custom1, custom2]

        data = '\n'.join([section_header] +
                         ["%s = %s" % x for x in custom_filters])

        with methods.tempinput(data) as path:
            loaded = list(passband.load_custom_filters(path))
            self.assertEqual(len(custom_filters), len(loaded))
            for input, output in zip(custom_filters, loaded):
                # For example, ('R-EROS', 'R (EROS-2 survey)')
                self.assertEqual(input[0], output[0])
                self.assertEqual(input[1], output[1])

        def assert_returns_nothing(path):
            """ Assert that Passband.load_custom_filters() does not return
            anything when the 'path' configuration file is parsed"""
            self.assertEqual([], list(passband.load_custom_filters(path)))

        # Does not return anything if:
        # (1) The configuration file file does not exist
        with methods.tempinput('') as path:
            pass
        self.assertFalse(os.path.exists(path))
        assert_returns_nothing(path)

        # (2) Does not have the CUSTOM_SECTION section
        data = ""
        with methods.tempinput(data) as path:
            assert_returns_nothing(path)

        # (3) The CUSTOM_SECTION section is empty
        data = section_header
        with methods.tempinput(data) as path:
            assert_returns_nothing(path)

