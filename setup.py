#! /usr/bin/env python2

# Copyright (c) 2012 Victor Terron. All rights reserved.
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

from __future__ import with_statement

""" Just a shadow of what the actual setup.py will look like, this simple
script, when run, creates the IRAF login script (login.cl file) and the user
parameters directory (uparm) in the directory where LEMON is installed. This
directory is automatically detected as it is assumed to be the same in which
this file is located, so it does not depend on where it is executed from.

"""

import os
import os.path
import errno
import shutil
import subprocess

MKIRAF_BIN = 'mkiraf'
TERM_TYPE = 'xgterm'
LOGIN_FILE = 'login.cl'
UPARM_DIR = 'uparm'
PYRAF_CACHE = 'pyraf'

# The LEMON configuration file
CONFIG_FILENAME = '~/.lemonrc'
CONFIG_PATH = os.path.expanduser(CONFIG_FILENAME)

def mkiraf(path):
    """ Create the IRAF login script and the uparm directory.

    This function implements a high-level wrapper around IRAF's mkiraf, which
    initializes the IRAF login script and the user parameters directory in the
    path given as input. Any existing login script or uparm directory is
    silently overwritten, although mkiraf makes a backup of the former by
    appending '.OLD' to its name. In this manner, the original login script
    'login.cl' becomes 'login.cl.OLD'.

    The login script that this function creates chooses 'xgterm' as its
    terminal type. It is usually the best choice, being a xterm-like terminal
    program written specifically to work with IRAF. Note that, although we are
    required to choose a terminal type when mkiraf is run, LEMON never needs
    to use it.

    """

    os.chdir(path)

    # Avoid having to answer the "Initialize uparm? (y|n):" question
    if os.path.exists(UPARM_DIR):
        shutil.rmtree(UPARM_DIR)

    # mkiraf reads the terminal type ("Enter terminal type:") from the user,
    # not accepting it as a command line argument. Thus, we need to use a pipe
    # to send the terminal type to mkiraf's stdin.
    args = [MKIRAF_BIN]
    p = subprocess.Popen(args, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    p.communicate(input = TERM_TYPE)

    if p.returncode:
        msg = "execution of %s failed" % MKIRAF_BIN
        raise subprocess.CalledProcessError(p.returncode, args)

    with open(LOGIN_FILE, 'rt') as fd:
        for line in fd:
            splitted = line.split()
            if len(splitted) == 2:
                if splitted[0] == 'stty':
                    if splitted[1] == TERM_TYPE:
                        break
                    else:
                        msg = "terminal type wasn't set correctly"
                        raise ValueError(msg)
        else:
            msg = "terminal type not defined in %s" % LOGIN_FILE
            raise ValueError(msg)

def mkdir_p(path):
    """ Create a directory, give no error if it already exists.

    This implements the functionality of Unix `mkdir -p`, creating a
    directory but without giving any error if it already exists and
    making parent directories as needed.
    [URL] http://stackoverflow.com/a/600612/184363

    """

    try:
        os.mkdir(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


if __name__ == "__main__":

    lemon_path = os.path.dirname(os.path.realpath(__file__))
    print "Setting up IRAF's %s in %s ..." % (LOGIN_FILE, lemon_path) ,
    mkiraf(lemon_path)
    print 'done.'

    print "Creating pyraf/ directory for cache...",
    pyraf_cache_path = os.path.join(lemon_path, PYRAF_CACHE)
    mkdir_p(pyraf_cache_path)
    print 'done.'

