#! /usr/bin/env python

import imp
import sys

def version_to_str(version):
    """ From (0, 2, 4) to '0.2.4' for example """
    return '.'.join(str(x) for x in version)

def str_to_version(version):
    """ From '0.2.4' to (0, 2, 4), for example """
    return tuple(int(x) for x in version.split('.'))

class RequireModuleVersionHook(object):

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

