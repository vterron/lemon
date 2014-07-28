#! /usr/bin/env python

# Copyright (c) 2014 Victor Terron. All rights reserved.
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

import os.path
import pyfits
import tempfile

def fix_DSS_image(path):
    """ Fix images downloaded from STScI DSS before we can do photometry.

    The images that we download from the STScI Digitized Sky Survey (DSS) to
    use in our unit tests are multi-extension FITS files. This causes IRAF's
    qphot to complain that 'FXF: must specify which FITS extension' and exit.
    Also, the exposure time is stored in the header in minutes, instead of in
    seconds as our code expects.

    This function accepts as input a DSS image and returns the path to a
    temporary copy of it with only the primary HDU, which is the one we are
    interested in; and the exposure time in seconds. You are responsible for
    deleting this temporary FITS file when done with it.

    """

    kwargs = dict(suffix = os.path.splitext(path)[-1])
    with tempfile.NamedTemporaryFile(**kwargs) as fd:
        output_path = fd.name

    with pyfits.open(path) as hdu:
        # Keep only the primary header
        del hdu[1:]
        header = hdu[0].header
        # From minutes to seconds
        header['EXPOSURE'] *= 60
        hdu.writeto(output_path)

    return output_path

