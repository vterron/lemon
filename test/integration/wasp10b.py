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
See README.md.

Expects two environment variables:
  - $NCORES: number of --cores for the `photometry` and `diffphot` commands.
  - $WASP10_DATA: path to the directory with the images of WASP10b's transit.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections
import contextlib
import functools
import itertools
import os
import os.path
import shutil
import subprocess
import sys

from absl.testing import absltest, parameterized

COORDINATES_FILE = "WASP10b-coordinates.txt"
WASP10_RA = 348.9930
WASP10_DEC = 31.4629
CURVE_FILTER = "V"
DECIMAL_PLACES = 6
GOLDEN_LIGHT_CURVE_FILE = "WASP10b-golden-curve.txt"
EXPORTED_LIGHT_CURVE_FILE = "WASP10b-exported_curve.txt"

cmd = functools.partial(subprocess.check_call, shell=True)

# Returns the absolute path of a file in the test/ directory.
test_path = functools.partial(
    os.path.join, os.environ["TRAVIS_BUILD_DIR"], "test/integration"
)


@contextlib.contextmanager
def cd(dst):
    """Temporarily changes the working directory."""
    cwd = os.getcwd()
    os.chdir(dst)
    try:
        yield
    finally:
        os.chdir(cwd)


def link(path):
    """Hard links in the current directory all the files under 'path'."""
    cmd("cp -v --link {}/* .".format(path))


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

    with open(path, "rb") as fd:
        for line in fd:
            line = line.decode("utf8")
            chunks = line.lstrip("|").rstrip("|\n").split("|")
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

    with open(test_path(GOLDEN_LIGHT_CURVE_FILE), "rb") as fd:
        for line in fd:
            line = line.decode("utf8")
            if line.strip().startswith("#"):
                continue
            chunks = line.split("\t")
            if len(chunks) == 3:
                yield DataPoint(*[float(x) for x in chunks])


class WASP10Test(parameterized.TestCase):
    """End-to-end test using WASP-10b data."""

    @parameterized.named_parameters(
        # Use parameterized instead of reading the environment variable
        # directly in the function body so that the resulting test name
        # includes the CPUs count (e.g. "WASP10Test.test_e2e_4_cores").
        #
        # Underscores in the test names as create_tempdir() includes them
        # verbatim in the resulting path, and IRAF raises an error if the
        # path to the FITS images has whitespaces.
        {
            "testcase_name": "{}_cores".format(os.environ["NCORES"]),
            "ncores": int(os.environ["NCORES"]),
        },
    )
    def test_e2e(self, ncores):

        with cd(self.create_tempdir().full_path):
            copy_coordinates_file()
            link(os.environ["WASP10_DATA"])

            cmd(
                "lemon "
                "photometry "
                "WASP10b-mosaic.fits "
                "WASP10b-*a.fits "
                "WASP10b-photometry.LEMONdB "
                "--coordinates={coordinates} "
                "--gain=1.31 "
                "--cores={cores}".format(
                    coordinates=COORDINATES_FILE,
                    cores=ncores,
                )
            )

            cmd(
                "lemon "
                "diffphot "
                "WASP10b-photometry.LEMONdB "
                "WASP10b-diffphot.LEMONdB "
                "--cores={cores}".format(
                    cores=ncores,
                )
            )

            # Export the generated light curve to a text file.
            cmd(
                "lemon "
                "export "
                "WASP10b-diffphot.LEMONdB "
                "{ra} "
                "{dec} "
                "{pfilter} "
                "--decimal_places={nplaces} "
                "--output_file={path}".format(
                    ra=WASP10_RA,
                    dec=WASP10_DEC,
                    pfilter=CURVE_FILTER,
                    nplaces=DECIMAL_PLACES,
                    path=EXPORTED_LIGHT_CURVE_FILE,
                )
            )

            # Now load the generated light curve and compare it to the one in the golden file.
            golden_curve = load_golden_file()
            actual_curve = load_exported_light_curve(EXPORTED_LIGHT_CURVE_FILE)

            # izip_longest() so that we catch regressions where fewer than the
            # expected light curve data points are generated.
            for want, got in itertools.izip_longest(golden_curve, actual_curve):
                print("{}: ".format(want), end="")
                self.assertAlmostEqual(want.jd, got.jd)
                self.assertAlmostEqual(want.mag, got.mag)
                self.assertAlmostEqual(want.snr, got.snr)
                print("OK")


if __name__ == "__main__":
    absltest.main()
