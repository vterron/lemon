#! /usr/bin/env python

# Copyright (c) 2012 Victor Terron. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
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

import ConfigParser
import os.path

CONFIG_FILENAME = '.juicerc'
CONFIG_PATH = os.path.expanduser('~/%s' % CONFIG_FILENAME)

VIEW_SECTION = 'view'
VIEW_SEXAGESIMAL = 'sexagesimal'
VIEW_DECIMAL = 'decimal'
PERIODS_UNIT = 'periods'
PERIODS_DAYS, PERIODS_HHMMSS, PERIODS_SECONDS = range(3)
PLOT_AIRMASSES = 'airmasses'

class Configuration(ConfigParser.SafeConfigParser):
    """ Just a quite simple wrapper to automatically have the configuration
    file loaded at instantiation and written to disk on deletion"""

    DEFAULT_CONFIG = '\n'.join(
    ["[view]",
     "sexagesimal = 1",
     "periods = 1",
     "decimal = 0",
     "airmasses = 1"])

    def __init__(self, path, update = True):
        """ Parse a configuration file, creating and populating it with
        the default options in case 'path' does not exist """

        ConfigParser.SafeConfigParser.__init__(self)

        if not os.path.exists(path):
            with open(path, 'wt') as fd:
                fd.write(self.DEFAULT_CONFIG)

        self.read([path])
        self.path = path

    def __del__(self):
        with open(self.path, 'wt') as fd:
            self.write(fd)

