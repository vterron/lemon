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

import calendar
import functools
import json
import os.path
import requests
import subprocess
import tempfile
import time
import warnings

# LEMON module
import methods

LEMON_DIR = os.path.dirname(os.path.abspath(__file__))
COMMITS_URL = 'https://api.github.com/repos/vterron/lemon/commits?page=1&per_page=1'
GITHUB_CACHE_FILE = os.path.join(LEMON_DIR, '.last-github-commit-cache.json')

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

    # Delete the cache file, forcing a request to the GitHub API the next time
    # get_last_github_commit() is called. Otherwise, although up-to-date after
    # running `git pull`, we would be comparing the SHA1 hash of HEAD with the
    # cached one.

    try:
        os.unlink(GITHUB_CACHE_FILE)
    except OSError:
        pass

    args = ['git', 'pull']
    with methods.tmp_chdir(LEMON_DIR):
        return subprocess.call(args)

class FileCache(object):
    """ Interface to cache data to disk.

    This class allows to easily write and read data from a JSON file, via
    its set() and get() methods. The file is opened and closed every time an
    operation is made, so we do not need to worry about closing it when done
    with the FileCache object.

    """

    def __init__(self, path):
        self.path = path

    def up_to_date(self, max_hours = 1):
        """ Determine whether the cache file has expired.

        Return True if the cache file was last modified less than 'max_hours'
        ago; and False otherwise. If the file does not exist, returns False
        too. In this manner, any time that False is returned we know that we
        cannot use the cached value.

        """

        try:
            max_seconds = max_hours * 3600
            cache_mtime = os.path.getmtime(self.path)
            return (time.time() - cache_mtime) <= (max_seconds)
        except OSError:
            return False

    def get(self):
        """ Return the contents of the JSON cache file """

        with open(self.path, 'rt') as fd:
            return json.load(fd)

    def set(self, *args):
        """ Write the received arguments to the JSON cache file """

        with open(self.path, 'wt') as fd:
            json.dump(args, fd)

def github_cache(func):
    """ Decorator to avoid unnecessarily querying the GitHub API too much.

    This decorator uses the GITHUB_CACHE_FILE file to cache the values returned
    by the decorated function. If the cache file was last modified less than
    an hour ago, the function call is skipped and the contents of the file
    returned. Otherwise, the function is called and the result file-cached
    before being returned.

    This function is expected to be used to decorate get_last_github_commit()
    since it would be impolite (and also inefficient) to make too many queries
    to the GitHub API. Anyway, if we did not care about that, the rate limit
    for unauthenticated requests would only allow us to make up to sixty
    requests per hour.

    """

    cache = FileCache(GITHUB_CACHE_FILE)
    @functools.wraps(func)
    def cachedf(*args, **kwargs):
        if not cache.up_to_date(max_hours = 1):
            cache.set(*func(*args, **kwargs))
        return cache.get()
    return cachedf

@github_cache
def get_last_github_commit(timeout = None):
    """ Return the short SHA1 of the last commit pushed to GitHub.

    Use the GitHub API to get the SHA1 hash of the last commit pushed to the
    LEMON repository, and then obtain its short version with `git rev-parse`.
    Returns a two-element tuple with (a) the short SHA1 and (b) date of the
    last commit as a Unix timestamp.

    The 'timeout' keyword argument defines the number of seconds after which
    the requests.exceptions.Timeout exception is raised if the server has not
    issued a response. Note that this is not the same as a time limit on the
    entire response download.

    """

    # All API requests MUST include a valid User-Agent header [..] We request
    # that you use your GitHub username, or the name of your application, for
    # the User-Agent header value. This allows us to contact you if there are
    # problems. [https://developer.github.com/v3/#user-agent-required]

    headers = {'User-Agent': 'vterron'}
    kwargs = dict(headers = headers, timeout = timeout)
    r = requests.get(COMMITS_URL, **kwargs)
    last_commit = r.json()[0]
    hash_ = last_commit['sha']
    date_str = last_commit['commit']['committer']['date']

    # Timestamps are returned in ISO 8601 format: "YYYY-MM-DDTHH:MM:SSZ", where
    # Z is the zone designator for the zero UTC offset (that is, the time is in
    # UTC). Parse the string and convert it to a Unix timestamp value.

    fmt = "%Y-%m-%dT%H:%M:%SZ"
    date_struct = time.strptime(date_str, fmt)
    date_ = calendar.timegm(date_struct)

    args = ['git', 'rev-parse', '--short', hash_]
    with methods.tmp_chdir(LEMON_DIR):
        with tempfile.TemporaryFile() as fd:
            subprocess.check_call(args, stdout = fd)
            fd.seek(0)
            short_hash = fd.readline().strip()
            return short_hash, date_

def check_up_to_date(timeout = None):
    """ Issue a warning if there are unmerged changes on GitHub.

    Compare the SHA1 hash of the last commit in the local LEMON Git repository
    with that pushed to GitHub. If they differ, issue a warning to let the user
    know that there is a newer version available and that `lemon --update` can
    be used to update the installation. Do nothing if we are up to date.

    The 'timeout' parameter is passed to get_last_github_commit().

    """

    current_revision  = get_git_revision()
    short_hash, date_ = get_last_github_commit(timeout = timeout)
    if short_hash not in current_revision:
        msg = ("Your current revision is '%s', but there is a more recent "
               "version (%s, %s) available on GitHub. You may use `lemon "
               "--update` to retrieve these changes.")
        args = (current_revision, short_hash, date_)
        warnings.warn(msg % args)
