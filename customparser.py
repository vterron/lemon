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

import copy
import optparse
import textwrap

# LEMON modules
import passband

class NewlinesFormatter(optparse.IndentedHelpFormatter):
    """ This quick-and-dirty trick prevents optparse from stripping newlines
    (using textwrap) when the description of the module is printed. This should
    be acceptable enough until the transition to argparse is made. """

    def _format_text(self, text):
        text_width = self.width - self.current_indent
        indent = ' ' * self.current_indent
        # Wrap one paragraph at a time, then concatenate them
        formatted_text = ""
        for paragraph in text.split('\n\n'):

            formatted_text += textwrap.fill(paragraph.strip(),
                                            text_width,
                                            initial_indent=indent,
                                            subsequent_indent=indent)
            formatted_text += '\n\n'

        return formatted_text.rstrip()

def check_passband(option, opt, value):
    """ Type-checking function for the 'passband' optparse option type.

    This is the type-checking function for 'passband', a custom optparse type
    that accepts a string with the name of a photometric filter and returns
    it as a passband.Passband object. 'option' is an optpase.Option instance,
    'opt' is the option string (e.g., -f), and 'value' is the string from the
    command line that must be checked and converted to a Passband object.

    In case of doubt, please refer to:
    http://docs.python.org/2.7/library/optparse.html#adding-new-types

    """

    try:
        return passband.Passband(value)
    except ValueError, e:
        msg = "option %s: invalid photometric filter: %r (%s)"
        raise optparse.OptionValueError(msg % (opt, value, e))

class PassbandOption(optparse.Option):
    """ Custom optparse option type encapsulating a photometric filter.

    This subclass of optparse's Option class implements 'passband', a custom
    option type: it receives a string with the name of a photometric filter,
    such as 'Johnson V', and returns it as a passband.Passband object. This
    option supports all the photometric systems allowed by the Passband class.

    """

    # A tuple of type names
    TYPES = optparse.Option.TYPES + ('passband',)
    # A dictionary mapping type names to type-checking functions
    TYPE_CHECKER = copy.copy(optparse.Option.TYPE_CHECKER)
    TYPE_CHECKER['passband'] = check_passband

def get_parser(description):
    """ Return the OptionParser object used in the LEMON modules.

    This function instantiates an optparse.OptionParser object and returns it.
    Its 'description' argument (a paragraph of text giving a brief overview of
    the program) is set to the value of the argument of the same name, while
    the NewlinesFormatter class is used for printing help text ('formatter'
    argument). This parser adds a custom option type, 'passband', which
    receives a string with the name of a photometric filter and converts it to
    a passband.Passband object. Neither the default -h nor --help options are
    included.

    """

    kwargs = dict(description = description,
                  add_help_option = False,
                  formatter = NewlinesFormatter(),
                  option_class = PassbandOption)
    parser = optparse.OptionParser(**kwargs)
    return parser

def clear_metavars(parser):
    """ Set all the meta-variables of an OptionParser to a whitespace.

    This is a hackish convenience function to set the meta-variables of all the
    options of an OptionParser, independently of whether they are contained in
    option groups, to the string ' '. This is not an empty string, but a string
    consisting of a whitespace: this is necessary because if the meta-variable
    evaluates to False the OptionParser converts the destination variable name
    to uppercase and uses that instead.

    Using this function on an OptionParser instance clears the meta-variables,
    which means that where the help message showed '--filename=FILE' now only
    '--filename=' will be displayed, with the equals sign indicating the fact
    that the option takes a value.

    """

    EMPTY_VALUE = ' '
    for option in parser.option_list:
        option.metavar = EMPTY_VALUE
    for group in parser.__dict__['option_groups']:
        for option in group.option_list:
            option.metavar = EMPTY_VALUE

