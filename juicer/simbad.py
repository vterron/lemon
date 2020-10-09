#! /usr/bin/env python2

# Copyright (c) 2014 Victor Terron. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
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

import webbrowser


def coordinate_query(ra, dec):
    """Look up an astronomical object in the SIMBAD database.

    Submit a coordinate-query to the SIMBAD database, using the given right
    ascension and declination (ICRS, J2000, 2000.0). Display it using the
    default browser, opening the page in the same window if possible.

    """

    SIMBAD_URL = (
        "http://simbad.u-strasbg.fr/simbad/sim-basic?"
        "Ident={0}+{1}&submit=SIMBAD+search"
    )
    url = SIMBAD_URL.format(ra, dec)
    webbrowser.open(url)
