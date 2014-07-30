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
import subprocess
import tempfile

# LEMON module
import methods

def get_git_revision():
    """ Return a human-readable revision number of the LEMON Git repository.

    Return the output of git-describe, the command that returns an identifier
    which tells us how far (number of commits) off a tag we are and the hash of
    the current HEAD. This allows us to precisely pinpoint where we are in the
    Git repository.

    """

    # --long: always output the long format even when it matches a tag
    # --dirty: describe the working tree; append '-dirty' if necessary
    # --tags: use any tag found in refs/tags namespace

    # check_output() is new in 2.7; we need 2.6 compatibility
    args = ['git', 'describe', '--long', '--dirty', '--tags']
    path = os.path.dirname(os.path.abspath(__file__))
    with methods.tmp_chdir(path):
        with tempfile.TemporaryFile() as fd:
            subprocess.check_call(args, stdout = fd)
            fd.seek(0)
            return fd.readline().strip()
