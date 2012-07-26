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

import contextlib
import functools
import gtk

@contextlib.contextmanager
def destroying(thing):
    try:
        yield thing
    finally:
        thing.destroy()

def show_message_dialog(parent_window, title, msg, msg_type = gtk.MESSAGE_INFO):

    kwargs = dict(type = msg_type,
                  buttons = gtk.BUTTONS_CLOSE,
                  message_format = msg)
    msg_dlg = gtk.MessageDialog(**kwargs)
    msg_dlg.set_title(title)
    msg_dlg.set_transient_for(parent_window)
    msg_dlg.set_position(gtk.WIN_POS_CENTER_ON_PARENT)
    res = msg_dlg.run()
    msg_dlg.destroy()
    return res

show_error_dialog = functools.partial(show_message_dialog,
                                      msg_type = gtk.MESSAGE_ERROR)
