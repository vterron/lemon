#! /usr/bin/env python

# Copyright (c) 2019 Victor Terron. All rights reserved.
#
# This file is part of LEMON.
#
# LEMON is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
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

import re

def DMS_to_DD(degrees, arcminutes, arcseconds):
    """ Degrees, arcminutes, arcseconds to decimal degrees conversion. """

    decimal = abs(degrees) + float(arcminutes)/60 + float(arcseconds)/3600
    if degrees < 0:
        decimal = -decimal
    return decimal

def DD_to_DMS(decimal_degrees):
    """ Decimal degrees to degrees, arcminutes, arcseconds conversion. """

    # math.modf(x) returns the fractional and integer parts of x
    arcminutes, degrees = math.modf(decimal_degrees)
    arcminutes = abs(arcminutes) * 60  # do not propagate the minus sign, if any
    arcseconds, arcminutes = math.modf(arcminutes)
    arcseconds *= 60
    return int(degrees), int(arcminutes), arcseconds

def HMS_to_DD(hours, minutes, seconds):
    """ Hours, minutes, seconds to decimal degrees conversion. """

    decimal = abs(hours)*15 + float(minutes)/60*15 + float(seconds)/3600*15
    if hours < 0:
        decimal = -decimal
    return decimal

def DD_to_HMS(decimal_degrees):
    """ Decimal degrees to hours, minutes, seconds conversion. """

    # math.modf(x) returns the fractional and integer parts of x
    minutes, hours = math.modf(decimal_degrees / 15.0)
    hours = abs(hours)    # do not propagate the minus sign, if any
    minutes = abs(minutes) * 60
    seconds, minutes = math.modf(minutes)
    seconds *= 60
    return int(hours), int(minutes), seconds

def ra_str(decimal_degrees):
    """ Return the string representation of a right ascension.
    Example: 14h 03m 12.64s """

    ra_hours, ra_min, ra_sec = DD_to_HMS(decimal_degrees)
    ra_hours = '%02d' % ra_hours
    ra_min = '%02d' % ra_min
    ra_sec = '%05.2f' % ra_sec
    return "%sh %sm %ss" % (ra_hours, ra_min, ra_sec)

def dec_str(decimal_degrees):
    """ Return the string representation of a declination.
    Example: +57d 33m 10.75s """

    dec_deg, dec_arcmin, dec_arcsec = DD_to_DMS(decimal_degrees)
    dec_deg = '%+03d' % dec_deg
    dec_arcmin = '%02d' % dec_arcmin
    dec_arcsec = '%05.2f' % dec_arcsec
    return "%sd %sm %ss" % (dec_deg, dec_arcmin, dec_arcsec)

def load_coordinates(path):
    """ Load a list of celestial coordinates from a text file.

    Parse a text file containing the celestial coordinates of a series of
    astronomical objects, one per line, and return a generator that yields (ra,
    dec, pm_ra, pm_dec) tuples. The file must have two columns, with the right
    ascension and declination (in decimal degrees) and, optionally, two other
    columns with the proper motion in right ascension and declination (in
    seconds of arc per year) surrounded by brackets. For example:

      269.456271 4.665281
      269.452075 4.693391 [-0.79858] [10.32812] # Barnard's Star
      269.466450 4.705625 [0.0036] [-.0064]     # TYC 425-262-1

    The four-element tuples contain the right ascension, declination, proper
    motion in right ascension and proper motion in declination, respectively.
    Nones are used in case the proper motion of an astronomical object is not
    specified. Empty lines, as well as comments (which start with the hash
    character, #, and extend to the end of the physical line) are ignored.

    ValueError is raised (a) if in any line there is a number of coordinates
    other than two (right ascension and declination) or the proper motions are
    not surrounded by brackets, (b) if any right ascension or declination is
    out of range or (c) if the proper motion in right ascension is specified
    but not that in declination, or vice versa.

    """

    with open(path, 'rt') as fd:
        for line in fd:

            # Ignore comments
            line = line.split('#')[0]

            # Ignore empty lines
            regexp = "^\s*$"
            if re.match(regexp, line):
                continue

            kwargs = dict(float = "([+-]?\d+(\.\d+)(?:[eE][+\-]?\d+)*)")
            regexp = ("^\s*"
                      "(?P<ra>{float})"
                      "\s+"
                      "(?P<dec>{float})"
                      "("
                        "\s+"
                        "\[\s*(?P<pm_ra>{float})\s*\]"
                        "\s+"
                        "\[\s*(?P<pm_dec>{float})\s*\]"
                      ")?"
                      "\s*$")

            match = re.match(regexp.format(**kwargs), line)

            if not match:
                msg = ("Unable to parse line %r. Astronomical objects must be "
                       "listed one per line with coordinate values in columns one "
                       "(right ascension) and two (declination). Proper motions "
                       "may be optionally specified in columns three (ra) and "
                       "four (dec), surrounded by brackets -- but, in that case, "
                       "both of them are required." % line)
                raise ValueError(msg)

            ra = float(match.group('ra'))
            if not 0 <= ra < 360:
                msg = "Right ascension '%s' not in range [0, 360[ degrees"
                raise ValueError(msg % ra)

            dec = float(match.group('dec'))
            if not -90 <= dec <= 90:
                msg = "Declination '%s' not in range [-90, 90] degrees"
                raise ValueError(msg % dec)

            pm_ra = match.group('pm_ra')
            if pm_ra is not None:
                pm_ra = float(pm_ra)

            pm_dec = match.group('pm_dec')
            if pm_dec is not None:
                pm_dec = float(pm_dec)

            yield ra, dec, pm_ra, pm_dec
