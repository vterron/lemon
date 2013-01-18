#! /usr/bin/env python

# Copyright (c) 2013 Victor Terron. All rights reserved.
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

""" This module provides access to a series of FITS images that can be used in
the unit tests. In order to keep the size of the repository as small as
possible, these images are not included along with the source code, but instead
downloaded automatically from the STScI Digitized Sky Survey the first time the
module is imported and copied to the ./test_data/fits directory. Any removed
images will be downloaded again the next time the module is imported. In
practical terms, the only thing that we need to know is that the module-level
variable TEST_IMAGES is a set with the paths to several FITS images that were
transparently downloaded from the Digitized Sky Survey.

"""

import functools
import re
import os
import os.path
import sys
import urllib

DATA_DIR = os.path.join(os.path.dirname(__file__), './test_data')
IMAGES_DIR = os.path.join(DATA_DIR, 'fits')

def get_dss_image(path, ra, dec):
    """ Download an image from the STScI Digitized Sky Survey.

    This method uses the DSS CGI script to automatically download a FITS image
    from the STScI Digitized Sky Survey, copying it to 'path'. The version of
    the survey is POSS2/UKSTU Infrared (second generation, 1 arcsex/pixel).
    The rest of the parameters use their default values, which means that the
    coordinates are assumed to be J2000 and the dimensions of the FITS image
    are 15x15 arcmins.

    DSS Web Access Script Interface Control Document:
    http://archive.stsci.edu/dss/script_usage.html

    If any exception is raised by urllib.urlretrieve() (and this includes
    KeyboardException, if the user decides to stop the execution), the file
    being downloaded is silently removed before the exception is re-raised.
    This guarantees that we will not be left with truncated data.

    """

    base_url = "http://archive.stsci.edu/cgi-bin/dss_search?"
    parameters = dict(v = 'poss2ukstu_ir', ra = ra, dec = dec)
    url = base_url + urllib.urlencode(parameters)

    # For example, "Downloading test/test_data/fits/IC_5146.fits: 87 %"
    status = functools.partial(("Downloading %s: {0:d} %%" % path).format)
    sys.stdout.write(status(0))
    sys.stdout.flush()

    def update_status(count, block_size, total_size):
        percent = int(count * block_size / total_size * 100)
        sys.stdout.write('\r' + status(percent))
        sys.stdout.flush()

    try:
        urllib.urlretrieve(url, filename = path, reporthook = update_status)
    except:
        try: os.unlink(path)
        except: pass
        raise
    finally:
        print

# Map each object to its right ascension and declination
TEST_OBJECTS = {'IC 5070' : (312.75, 44.37),
                'IC 5146' : (328.35, 47.267),
                'Messier 92': (259.281, 43.136),
                'NGC 2264' : (100.242, 9.895),
                'RMC 136' : (98.67, 4.03),
                'Serpens' : (277.454, 1.247),
                'Orion' : (83.822, -5.391),
                'Trapezium' : (83.819, -5.387),
                'Trumpler 37' : (324.536, 57.447)}

def get_image_path(name):
    """ Determine the local path to which to download the image of an object.

    The base name of the image is that of the astronomical object, but with any
    whitespace characters replaced with underscores and the '.fits' extension.
    The directory where all the images are downloaded is IMAGES_DIR.

    """
    basename = '%s.fits' % re.sub('\s+', '_', name)
    path = os.path.join(IMAGES_DIR, basename)
    return os.path.normpath(path)

if not os.path.exists(IMAGES_DIR):
    os.makedirs(IMAGES_DIR)

TEST_IMAGES = set()

for name, (ra, dec) in sorted(TEST_OBJECTS.items()):
    path = get_image_path(name)
    if not os.path.exists(path):
        get_dss_image(path, ra, dec)
    TEST_IMAGES.add(path)

