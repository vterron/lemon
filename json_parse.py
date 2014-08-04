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

import collections
import copy
import json
import operator

# LEMON modules
import passband

typename = 'CandidateAnnuli'
field_names = "aperture, annulus, dannulus, stdev"
class CandidateAnnuli(collections.namedtuple(typename, field_names)):
    """ Encapsulate the quality of a set of photometric parameters.

    How do we determine how good a set of parameters for aperture photometry
    is? In order to compare them, we need to identify the most constant stars
    (or, by extension, any other astronomical object) in the field and compute
    their light curves. The better the aperture, annulus and dannulus that we
    use are, the lower the standard deviation of the resulting curves.

    This class simply encapsulates these four values, mapping the parameters
    for aperture photometry (aperture, annulus and dannulus) to the standard
    deviation of the light curves of the most constant astronomical objects.

    Fields:
    aperture - the aperture radius, in pixels.
    annulus - the inner radius of the sky annulus, in pixels.
    dannulus - the width of the sky annulus, in pixels.
    stdev - the median, arithmetic mean or a similar statistical measure of the
            standard deviation of the light curves of the astronomical objects
            when photometry is done using these aperture, annulus and dannulus
            values.

    """

    @staticmethod
    def dump(annuli, path):
        """ Save a series of CadidateAnnuli objects to a JSON file.

        Serialize 'annuli' to a JSON file. It must be a dictionary which maps
        each photometric filter (a Passband object) to a sequence of the
        corresponding CandidateAnnuli objects -- i.e., the different aperture
        photometric parameters that were evaluated for that filter. The output
        file will be mercilessly overwritten if it already exists.

        """

        # Being a subclass of tuple, JSON serializes namedtuples as lists. We
        # need to convert them to ordered dictionaries (namedtuple._asdict())
        # first, so that the field names are not lost in the serialization.
        data = copy.deepcopy(annuli)
        for values in data.itervalues():
            for index in xrange(len(values)):
                values[index] = values[index]._asdict()
            values.sort(key=operator.itemgetter('stdev'))

        # Use strings, not Passband objects, as keys
        for pfilter in data.keys():
            data[str(pfilter)] = data.pop(pfilter)

        with open(path, 'wt') as fd:
            kwargs = dict(indent=2, sort_keys=True)
            json.dump(data, fd, **kwargs)

    @classmethod
    def load(cls, path):
        """ Load a series of CandidateAnnuli objects from a JSON file.

        Deserialize a JSON file created with CandidateAnnuli.dump(), returning
        a dictionary which maps each photometric filter (a Passband object) to
        a list of the corresponding CandidateAnnuli objects. These lists are
        sorted in increasing order by the standard deviation ('stdev' attribute
        of the namedtuples), so that the one with the lowest standard deviation
        (and therefore the optimal for aperture photometry) is returned first.

        """

        with open(path, 'rt') as fd:
            data = json.load(fd)

        # Convert the dictionaries back to namedtuples, and then sort them by
        # their standard deviation, in increasing order.
        for values in data.itervalues():
            for index in xrange(len(values)):
                values[index] = cls(**values[index])
            values.sort(key=operator.attrgetter('stdev'))

        # Use Passband objects as keys
        for pfilter in data.keys():
            data[passband.Passband(pfilter)] = data.pop(pfilter)

        return data
