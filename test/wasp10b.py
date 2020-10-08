#! /usr/bin/env python2
# -*- coding: utf-8 -*-

# Copyright (c) 2020 Victor Terron. All rights reserved.
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

"""
End-to-end integration test using 198 FITS images of WASP-10b taken at OSN
(Sierra Nevada Observatory, Spain) on 2011-08-23. The URL of the 2.3 GiB .xz
file is stored in the WASP10_URL environment variable, stored in Travis CI
via https://docs.travis-ci.com/user/encryption-keys.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections
import contextlib
import functools
import os
import os.path
import subprocess
import sys

from absl.testing import absltest

CHECKSUMS_FILE = "SHA1SUMS"
COORDINATES_FILE = "WASP10b-coordinates.txt"
GOLDEN_LIGHT_CURVE_FILE = "WASP10b-golden-curve.txt"

cmd = functools.partial(subprocess.check_call, shell=True)

# Returns the absolute path of a file in the test/ directory.
test_path = functools.partial(os.path.join, os.environ["TRAVIS_BUILD_DIR"], 'test')


@contextlib.contextmanager
def cd(dst):
    """Temporarily changes the working directory."""
    cwd = os.getcwd()
    os.chdir(dst)
    try:
        yield
    finally:
        os.chdir(cwd)

def download():
    """Downloads and uncompresses the WASP-10b test data."""
    url = os.environ["WASP10_URL"]
    cmd("curl --remote-header-name --remote-name {}".format(url))
    zx_file = os.listdir(".")[0]
    cmd("tar -xvf {}".format(zx_file))
    cmd("sha1sum -c {}".format(CHECKSUMS_FILE))

def copy_coordinates_file():
    """Makes the coordinates .txt file available in the current directory.

    The end-to-end test runs in a temporary directory, so it is necessary to
    copy there the file with the astronomical coordinates of the objects on
    which to do photometry. The function uses the TRAVIS_BUILD_DIR environment
    variable, which on Travis CI contains the absolute path to the directory
    where the repository being built has been copied on the worker:
    https://docs.travis-ci.com/user/environment-variables/

    """
    cmd("cp -v {} .".format(test_path(COORDINATES_FILE)))


"""A single point in the light curve of an astronomical object."""
DataPoint = collections.namedtuple("DataPoint", "jd, mag, snr")

def load_exported_light_curve(path):
    """Loads a light curve written to disk with `lemon export`.

    These files are formatted like this:

    +--------------------------+----------------+-----------+------------+
    |        Date (UTC)        |       JD       |   Δ Mag   |    SNR     |
    +--------------------------+----------------+-----------+------------+
    | Wed Aug  3 20:49:55 2011 | 2455777.367999 | -0.948379 | 914.218013 |
    | Wed Aug  3 20:52:23 2011 | 2455777.369720 | -0.951874 | 930.153940 |
    [...]

    Args:
        path: the path to the text file with the exported light curve.
    Yields:
        DataPoint objects, one per data point in the light curve.

    """

    with open(path, 'rb') as fd:
        for line in fd:
            line = line.decode('utf8')
            chunks = line.lstrip('|').rstrip('|\n').split('|')
            if len(chunks) != 4:
                continue
            try:
                # Ignore the first column with the data in UTC.
                yield DataPoint(*[float(x) for x in chunks[1:]])
            except ValueError as e:
                assert "could not convert string to float" in str(e)


def load_golden_file():
    """Loads the golden file with the WASP-10b light curve.

    This file looks like the following:

    # Golden file with the light curve of a WASP-10b transit on Aug 4th 2011.
    #
    # Julian Date   Δ magnitude     SNR
    2455777.367999  -0.948379       914.218013
    2455777.369720  -0.951874       930.153940
    [...]

    Yields:
        DataPoint objects, one per data point in the light curve.
    """

    with open(test_path(GOLDEN_LIGHT_CURVE_FILE), 'rb') as fd:
        for line in fd:
            line = line.decode('utf8')
            if line.strip().startswith("#"):
                continue
            chunks = line.split('\t')
            if len(chunks) == 3:
                yield DataPoint(*[float(x) for x in chunks])


class WASP10Test(absltest.TestCase):
    """End-to-end test using WASP-10b data."""

    def test_e2e(self):
        with cd(self.create_tempdir().full_path):
            download()
            copy_coordinates_file()

            # TODO(vterron): look up and set the actual gain at OSN.
            cmd("lemon photometry WASP10b-mosaic.fits WASP10b-*a.fits WASP10b-photometry.LEMONdB --coordinates={} --gain=1".format(COORDINATES_FILE))
            cmd("lemon diffphot WASP10b-photometry.LEMONdB WASP10b-diffphot.LEMONdB")

if __name__ == "__main__":
    absltest.main()
