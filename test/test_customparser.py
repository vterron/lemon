#! /usr/bin/env python2

# Copyright (c) 2015 Victor Terron. All rights reserved.
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

# LEMON modules
from test import unittest
import customparser


class CustomParserTest(unittest.TestCase):
    def test_additional_options_callback(self):

        parser = customparser.get_parser(None)
        parser.add_option(
            "-o",
            action="callback",
            type="str",
            dest="parsed_options",
            default={},
            callback=customparser.additional_options_callback,
        )

        arguments = [
            "-o",
            "-v",
            "-o",
            "--batch",
            "-o",
            "-r=2",
            "-o",
            "--sigma 3",
            "-o=-h",
            "-o=--invert",
            "-o=-d=cubic",
            "-o=--no-verify = True",
        ]

        (options, args) = parser.parse_args(args=arguments)
        options = options.parsed_options
        self.assertIsNone(options["-v"])
        self.assertIsNone(options["--batch"])
        self.assertEqual(options["-r"], "2")
        self.assertEqual(options["--sigma"], "3")
        self.assertIsNone(options["-h"])
        self.assertIsNone(options["--invert"])
        self.assertEqual(options["-d"], "cubic")
        self.assertEqual(options["--no-verify"], "True")
