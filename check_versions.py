#! /usr/bin/env python

def version_to_str(version):
    """ From (0, 2, 4) to '0.2.4' for example """
    return '.'.join(str(x) for x in version)

def str_to_version(version):
    """ From '0.2.4' to (0, 2, 4), for example """
    return tuple(int(x) for x in version.split('.'))

