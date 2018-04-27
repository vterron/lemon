#! /usr/bin/env python2

import sys

# Several convenient features of unittest, such as assertRaises as a context
# manager and test skipping, are not available until Python 2.7. In previous
# versions, use unittest2, a backport of the new features.

if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest

