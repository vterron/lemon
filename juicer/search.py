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

import gtk

# LEMON modules
import glade
import util

def amplitudes_search(parent_window, builder):

    DEFAULT_NUMBER_MIN_MAX_POINTS = 5
    DEFAULT_NUMBER_STDEVS = 10
    DEFAULT_AMPSTDEV_RATIO = 2

    builder.add_from_file(glade.AMPLITUDES_DIALOG)
    object_name = 'amplitudes-search-dialog'
    with util.destroying(builder.get_object(object_name)) as dialog:

        dialog.set_resizable(False)
        dialog.set_title("Select stars by their amplitudes")
        dialog.add_button(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL)
        dialog.add_button(gtk.STOCK_OK, gtk.RESPONSE_OK)

        # Although the value of the adjustment is set in Glade, the spinbuttons
        # are all zero when the window is created, so we need to set them here
        namplitudes = builder.get_object('amplitudes-how-many')
        nstdevs = builder.get_object('comparison-stdevs-how-many')
        ratio = builder.get_object('min-amplitude-stdev-ratio')
        namplitudes.set_value(DEFAULT_NUMBER_MIN_MAX_POINTS)
        nstdevs.set_value(DEFAULT_NUMBER_STDEVS)
        ratio.set_value(DEFAULT_AMPSTDEV_RATIO)

        dialog.set_transient_for(parent_window)
        dialog.set_position(gtk.WIN_POS_CENTER_ON_PARENT)

        response =  dialog.run()
        if response == gtk.RESPONSE_OK:
            pass
        else:
            pass

