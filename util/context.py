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

import contextlib
import os
import tempfile


@contextlib.contextmanager
def tmp_chdir(path):
    """A context manager to temporarily change the working directory.

    This is a rather simple context manager to change the current working
    directory within a with statement, restoring the original one upon exit.

    """

    cwd = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(cwd)


@contextlib.contextmanager
def tempinput(data):
    """A context manager to work with StringIO-like temporary files.

    The open() built-in command only takes filenames, so we cannot feed to it a
    StringIO. This context manager writes 'data' to a temporary file, ensuring
    that it is cleaned up afterwards. In this manner, we can work in a similar
    manner to what we would have done with StringIO, as the data is written to
    disk only temporarily and transparently to us.

    with tempinput('some data\more data') as path:
        open(path)

    Taken from Martijn Pieters's answer on Stack Overflow:
    [URL] https://stackoverflow.com/a/11892712/184363

    """

    fd, path = tempfile.mkstemp()
    os.write(fd, data)
    os.close(fd)
    yield path
    os.unlink(path)
