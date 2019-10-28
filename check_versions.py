#! /usr/bin/env python

# Copyright (c) 2013 Victor Terron. All rights reserved.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

""" Use Python import hooks to enforce a minimum version of those modules
defined in the requirements.txt file. Another option would have been to check
the versions manually, but that would have forced us to do it multiple times
and across all the code, every time one of the modules were imported. And, as
more code is added in the future, we could forget to test for this. With this
solution we only need to define the hooks once, here, and this guarantee that
the modules will always have the required version, no matter how, when or
where they are imported.

"""

import imp
import re
import sys

# LEMON module
import astromatic

def version_to_str(version):
    """ From (0, 2, 4) to '0.2.4' for example """
    return '.'.join(str(x) for x in version)

def str_to_version(version):
    """ From '0.2.4' to (0, 2, 4), for example """
    return tuple(int(x) for x in version.split('.'))

class RequireModuleVersionHook(object):
    """ An import hook to enforce minimum versions of Python modules.

    This class implements both the module finder (find_module) and loader
    (load_module). It can be installed in sys.meta_path to intercept every
    attempt to import a new module, checking whether it has the required
    version and raising ImportError otherwise.

    """

    def __init__(self, fullname, min_version, vfunc):
        """ Instantiation method for the RequireModuleVersionHook class.

        The 'fullname' parameter is the a fully qualified name of the module
        that we want to import, such as 'pyfits' or 'scipy'. 'min_version' is a
        tuple of integers specifying the minimum version of the module, such as
        (1, 2, 1). Finally, 'vfunc' is a hook function that will be passed the
        module object, after it is imported, and that must return its version
        as a tuple of integers.

        """

        self.fullname = fullname
        self.min_version = min_version
        self.vfunc = vfunc

    def find_module(self, fullname, path=None):
        """ The module finder.

        Receive the fully qualified name of the module to be imported, along
        with, optionally, the path where it is supposed to be found. Return
        None if the module cannot be found by this particular finder and self
        (since this class also implements the module loader) otherwise.

        """

        if fullname == self.fullname:
            self.path = path
            return self
        return None

    def load_module(self, fullname):
        """ The module loader.

        Receive the fully qualified name of the module that we want to import.
        Then, import the module, call vfunc() to get its version as a tuple of
        integers and compare it to the specified minimum version: if the module
        version is equal to or higher than the required one; return the module
        object; otherwise, raise ImportError.

        """

        # If module already imported, return it
        if fullname in sys.modules:
            return sys.modules[fullname]

        module_info = imp.find_module(fullname, self.path)
        module = imp.load_module(fullname, *module_info)

        version = self.vfunc(module)
        if not version >= self.min_version:
            msg = "%s >= %s is required, found %s"
            args = (fullname,
                    version_to_str(self.min_version),
                    version_to_str(version))
            raise ImportError(msg % args)
        else:
            sys.modules[fullname] = module
            return module


def get__version__(module):
    """ Return module.__version__, as a tuple of integers """

    version = module.__version__
    # Extract '2.1.1' from '2.1.1-r1785' / '3.2' from '3.2.dev'
    regexp = "\d+(\.\d+)+"
    match = re.match(regexp, version)
    if match is None:
        msg = "cannot extract version from '%s' (%s)" % (version, module)
        raise Exception(msg)
    else:
        version = match.group(0)
        return str_to_version(version)

# For each module whose minimum version has been defined in requirements.txt,
# create an import hook and add it to sys.meta_path, which is searched before
# any implicit default finders or sys.path.

for module, version in [
  ('numpy', (1, 7, 1)),
  ('aplpy', (0, 9, 9)),
  ('scipy', (0, 12, 0)),
  ('matplotlib', (1, 2, 1)),
  ('mock', (1, 0, 1)),
  ('pyfits', (3, 1, 2)),
  ('pyraf', (2, 1, 1)),
  ('uncertainties', (2, 4, 1))]:
      hook = RequireModuleVersionHook(module, version, get__version__)
      sys.meta_path.append(hook)

# If the minimum version is not met, raise SExtractorUpgradeRequired as
# soon as possible, instead of waiting until we attempt to run SExtractor.
if astromatic.sextractor_version() < astromatic.SEXTRACTOR_REQUIRED_VERSION:
    raise astromatic.SExtractorUpgradeRequired()
