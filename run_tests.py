#! /usr/bin/env python

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

"""
This is a convenience script for loading all the unit tests in the test/
directory and running them. Test modules are automatically detected: those
files which start with 'test_' and have the '.py' extension (such as
test/test_passband.py) are treated as such. For Python 2.6 compatibility,
we cannot use the 'failfast' option, which causes the test run to stop on
the first error or failure.

"""

import os
import unittest

TESTS_PACKAGE = 'test'
TESTS_PREFIX = 'test_'
TESTS_EXTENSION = '.py'

if __name__ == "__main__":

    for path in os.listdir(TESTS_PACKAGE):
        if path.startswith(TESTS_PREFIX) and path.endswith(TESTS_EXTENSION):

            # The module name to import, such as 'test.test_passband'
            module = '%s.%s' % (TESTS_PACKAGE, path[:-len(TESTS_EXTENSION)])

            # The 'exit' parameter --so that the result is displayed on
            # standard output without calling sys.exit()-- was not added until
            # Python 2.7. For 2.6 compatibility, we do not use it and manually
            # catch the exception that sys.exit() raises (SystemExit) instead.
            try:
                retcode = unittest.main(module = module)
            except SystemExit:
                pass

