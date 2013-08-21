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


""" The class Passband encapsulates a filter of the photometric system. It is
    mainly (and probably only) intended to be used so that the filters with
    with the stars are observed can be automatically sorted according to their
    position in the electromagnetic spectrum. Otherwise, if the strings with
    the name of the filters were simply sorted lexicographically, the Johnson I
    filter, for example, would go before the V filter, while the right order is
    the other way around (as V has a shorter wavelenght than I, 551 vs 806 nm).

"""

import random
import re

class NonRecognizedPassband(ValueError):
    """ Raised when a the astronomical filter cannot be identified. """
    pass

class InvalidPassbandLetter(NonRecognizedPassband):
    """ Raised if the letter of the filter does not belong to the system.

    For example, this exception should be raised if we come across something
    like 'Johnson Z', as Z is not a filter of the Johnson photometric system
    (UBVRIJHKLMN).

    """
    pass

class Passband(object):
    """ Encapsulates a passband (or filter) of the photometric system. """

    JOHNSON_LETTERS = 'UBVRIJHKLMN'
    COUSINS_LETTERS = 'VRI'

    # The effective wavelength midpoint for each standard filter, in nanometers.
    # Taken from Binney's "Galactic Astronomy", 1998, ch.2.3.2, pp.53, seen at
    # Wikipedia. Note that this is a class (or static) variable, as there is no
    # need to define this information separately for each instance.

    wavelengths = {'U' : 365,
                   'B' : 445,
                   'V' : 551,
                   'R' : 658,
                   'I' : 806,
                   'Z' : 882,
                   'Y' : 1020,
                   'J' : 1220,
                   'H' : 1630,
                   'KS': 2151, # Omega2K (CAHA)
                   'K' : 2190,
                   'L' : 3450,
                   'M' : 4750}

    def __init__(self, name):
        """ Instantiation method for the Passband class.

        The class constructor receives as its only parameter the name of the
        photometric filter being instantiated, from which its letter is
        automatically extracted. This means the the only truly important
        information is, actually, the letter of the filter, a single character
        which may ge given in lower or upper case.

        'Johnson V', 'V Johnson', 'V (Johnson)' and even simply 'V' are
        acceptable filter names, from which the 'V' letter is easily identified.
        Examples of unacceptable filter names include 'V(Johnson)' [there should
        be a whitespace character between the letter and the left parenthesis],
        'Johnson' [there is no letter!] and 'Johnson (V)' [as there is no point
        in enclosing the letter in parentheses].

        """

        # This regular expression matches the start of the string or a
        # whitespace, followed a single (such as 'K') or two letters (e.g.,
        # 'Ks'), either lower or upper case and followed by the end of the
        # string or another whitespace.

        regexp = re.compile("(^|\s)([A-Z]{1,2})(\s|$)")
        match = regexp.findall(name.upper())
        if not 1 <= len(match) <= 2:
            raise NonRecognizedPassband(name)

        self.name = name
        self.letter = match[0][1] # match = [('', 'V', ' ')], e.g.

        if self.letter.upper() not in Passband.wavelengths:
            raise UnknownPassbandLetter(self.letter)

    @classmethod
    def all(cls):
        """ Return a sequence with all the passbands this class encapsulates """
        return sorted(cls(letter) for letter in cls.wavelengths.iterkeys())

    @property
    def wavelength(self):
        """ Return the effective wavelength midpoint, in nanometers. """
        return Passband.wavelengths[self.letter.upper()]

    def __str__(self):
        """ The 'informal' string representation of the filter. """
        return self.name

    def __repr__(self):
        """ The unambiguous string representation of a Passband object """
        return "%s(%r)" % (self.__class__.__name__,  self.name)

    def __cmp__(self, other):
        """ Called by comparison operations if rich comparison is not defined.

        The method returns a negative integer is 'self' has a shorter wavelength
        than 'other', zero if both wavelengths are equal and a positive integer
        if the wavelength of 'self' is longer than that of 'other'. Internally,
        this is done by means of comparing the effective wavelenght of both
        filters.

        """
        return self.wavelength - other.wavelength

    def __hash__(self):
        return self.wavelength

    @classmethod
    def random(cls):
        """ Return a random Passband """
        return cls(random.choice(cls.wavelengths.keys()))

    def different(self):
        """ Return a random filter other than this one """

        while True:
            passband = random.choice(Passband.all())
            if passband != self:
                return passband

