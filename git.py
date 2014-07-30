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
import requests
import subprocess
import tempfile

# LEMON module
import methods

LEMON_DIR = os.path.dirname(os.path.abspath(__file__))
COMMITS_URL = 'https://api.github.com/repos/vterron/lemon/commits?page=1&per_page=1'

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
    with methods.tmp_chdir(LEMON_DIR):
        with tempfile.TemporaryFile() as fd:
            subprocess.check_call(args, stdout = fd)
            fd.seek(0)
            return fd.readline().strip()

def git_update():
    """ Merge upstream changes into the local repository with `git pull` """

    args = ['git', 'pull']
    with methods.tmp_chdir(LEMON_DIR):
        return subprocess.call(args)

def get_last_github_commit():
    """ Return the short SHA1 of the last commit pushed to GitHub.

    Use the GitHub API to get the SHA1 hash of the last commit pushed to the
    LEMON repository, and then obtain its short version with `git rev-parse`.
    Returns a two-element tuple with (a) the short SHA1 and (b) date of the
    last commit.

    """

    # All API requests MUST include a valid User-Agent header [..] We request
    # that you use your GitHub username, or the name of your application, for
    # the User-Agent header value. This allows us to contact you if there are
    # problems. [https://developer.github.com/v3/#user-agent-required]

    headers = {'User-Agent': 'vterron'}
    r = requests.get(COMMITS_URL, headers = headers)
    last_commit = r.json()[0]
    hash_ = last_commit['sha']
    date_ = last_commit['commit']['committer']['date']

    args = ['git', 'rev-parse', '--short', hash_]
    with methods.tmp_chdir(LEMON_DIR):
        with tempfile.TemporaryFile() as fd:
            subprocess.check_call(args, stdout = fd)
            fd.seek(0)
            short_hash = fd.readline().strip()
            return short_hash, date_
