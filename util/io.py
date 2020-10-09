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

import logging
import os
import os.path
import shutil
import stat

# LEMON modules
import style


def determine_output_dir(output_directory, dir_suffix=None, quiet=False):
    """Ensure that the specified directory does exist.

    The method receives the intented output directory and guarantees that it
    can be used. If None was specified, it creates a temporary directory
    (readable, writable, and searchable only by the creating user ID), while if
    the specified directory does not exists, it creates it. The method always
    return the path to the output directory, whether it is the path to the
    just-created temporary directory or that originally specified.

    Keyword arguments:
    dir_suffix - the string to be appended to the output directory if it is to
                 be a temporary directory (that is, if the specified output
                 directory is None.
    quiet - if True, no messages will be printed to standard output.

    """

    # If no output directory was specified, create a temporary directory
    if not output_directory:
        if not quiet:
            print "%sNo output directory was specified." % style.prefix
        if dir_suffix is not None:
            temporary_directory = tempfile.mkdtemp(suffix=dir_suffix)
        else:
            temporary_directory = tempfile.mkdtemp()
        if not quiet:
            print "%sImages will be saved to temporary directory %s" % (
                style.prefix,
                temporary_directory,
            )
        return temporary_directory

    # If the path to the output directory was specified, check that it exists
    # and that it is writable. If it does not exist, try to create it.
    else:
        if not os.path.exists(output_directory):
            try:
                os.makedirs(output_directory)
            except OSError:
                msg = (
                    "The output directory '%s' could not be created. "
                    "Is the directory writable?" % output_directory
                )
                raise IOError(msg)

            if not quiet:
                print "%sThe output directory '%s' did not exist, so it " "had to be created." % (
                    style.prefix,
                    output_directory,
                )
        else:
            if not os.path.isdir(output_directory):
                raise IOError("%s is not a directory" % output_directory)
            if not os.access(output_directory, os.W_OK):
                raise IOError(
                    "The output directory '" + output_directory + "' is not writable."
                )
            if not quiet:
                print "%sImages will be saved to directory '%s'" % (
                    style.prefix,
                    output_directory,
                )

        return output_directory


def owner_writable(path, add):
    """Make the file owner writeable or unwriteable.

    The method either adds or removes write permission for the file's owner,
    depending on the value of the 'add' parameter. In other words, if 'add' is
    True, the method is equivalent to chmod u+w on the file, while if the
    parameter evaluates to False it is equivalent to a chmod u-w on 'path'.

    Note that doing this is not as straight-forward as you may think, as
    os.chmod(path, mode) _sets_ the mode of path to a numeric mode, which means
    that before doing that we need to determine which the current permissions
    of the file add.

    """

    # stat.ST_MODE returns the protection bits of the file.
    # stat.S_IMODE returns the portion of the file's mode that can be set by
    # os.chmod(), i.e., the file's permission bits, plus the sticky bit,
    # set-group-id, and set-user-id bits (on systems that support them).
    # stat.S_IWUSR = owner has write permission.

    mode = stat.S_IMODE(os.stat(path)[stat.ST_MODE])

    if add:
        mode |= stat.S_IWUSR  # bitwise (inclusive) OR
    else:
        mode ^= stat.S_IWUSR  # bitwise XOR (exclusive OR)

    os.chmod(path, mode)


def which(*names):
    """Search PATH for executable files with the given names.

    Replicate the functionality of Unix 'which', returning a list of the full
    paths to the executables that would be executed in the current environment
    if the arguments were given as commands in a POSIX-conformant shell. This
    is done by searching, in the directories listed in the environment variable
    PATH, for executable files matching the names of the arguments. If all the
    command are nonexistent or not executable, an empty list is returned.

    The code in this function is largely ripped from Twister's repository:
    https://twistedmatrix.com/trac/browser/trunk/twisted/python/procutils.py

    """

    result = []
    for directory in os.environ.get("PATH", "").split(os.pathsep):
        for name in names:
            path = os.path.join(directory, name)
            if os.access(path, os.X_OK) and os.path.isfile(path):
                result.append(path)
    return result


def clean_tmp_files(*paths):
    """Try to remove multiple temporary files and directories.

    Loop over the provided positional arguments, calling os.unlink() on files
    and shutil.rmtree() on directories. Errors never raise an exception, but
    are logged at DEBUG level. These files are considered to be 'temporary' in
    the sense that, being no longer necessary, they must be cleaned up, but
    they are not important enough as to require special handling if they cannot
    be deleted. After all, if they are located in /tmp/, as they are expected,
    they will eventually get cleared.

    """

    for path in paths:

        if os.path.isdir(path):

            msg = "Cleaning up temporary directory '%s'"
            logging.debug(msg % path)

            error_count = [0]

            def log_error(function, path, excinfo):
                """ Error handler for shutil.rmtree() """

                # nonlocal is not available in Python 2.x so, being it outside
                # of the local scope, we cannot use 'error_count' as a counter
                # and rebind it each time we come across an error. But we can
                # make it a list, which being mutable allows us to modify its
                # elements inside the function.

                error_count[0] += 1
                msg = "%s: error deleting '%s' (%s)"
                args = function, path, excinfo[1]
                logging.debug(msg % args)

            try:
                kwargs = dict(ignore_errors=False, onerror=log_error)
                shutil.rmtree(path, **kwargs)

            finally:
                msg = "Temporary directory '%s' deleted"
                if max(error_count) > 0:
                    msg += " (but there were failed removals)"
                logging.debug(msg % path)

        else:

            msg = "Cleaning up temporary file '%s'"
            logging.debug(msg % path)

            try:
                os.unlink(path)
            except OSError, e:
                msg = "Cannot delete '%s' (%s)"
                logging.debug(msg % (path, e))
            else:
                msg = "Temporary file '%s' removed"
                logging.debug(msg % path)
