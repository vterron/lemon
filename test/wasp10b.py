#! /usr/bin/env python2

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

import contextlib
import functools
import os
import os.path
import subprocess
import sys

from absl.testing import absltest

CHECKSUMS_FILE="SHA1SUMS"

cmd = functools.partial(subprocess.check_call, shell=True)

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


class WASP10Test(absltest.TestCase):
    """End-to-end test using WASP-10b data."""

    def test_e2e(self):
        with cd(self.create_tempdir().full_path):
            download()
            # TODO(vterron): look up and set the actual gain at OSN.
            cmd("lemon photometry WASP10b-mosaic.fits WASP10b-*a.fits WASP10b-photometry.LEMONdB --gain=1")
            cmd("lemon diffphot WASP10b-photometry.LEMONdB WASP10b-diffphot.LEMONdB")

if __name__ == "__main__":
    absltest.main()