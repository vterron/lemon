#! /usr/bin/env python2

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


def show_message_dialog(
    parent_window, title, msg, msg_type=gtk.MESSAGE_INFO, buttons=gtk.BUTTONS_CLOSE
):

    kwargs = dict(type=msg_type, buttons=buttons, message_format=msg)
    msg_dlg = gtk.MessageDialog(**kwargs)
    msg_dlg.set_title(title)
    msg_dlg.set_transient_for(parent_window)
    msg_dlg.set_position(gtk.WIN_POS_CENTER_ON_PARENT)
    res = msg_dlg.run()
    msg_dlg.destroy()
    return res


show_error_dialog = functools.partial(show_message_dialog, msg_type=gtk.MESSAGE_ERROR)


@contextlib.contextmanager
def gtk_sync():
    """A context manager to force an update to the GTK application window.

    By design, all GTK events (including window refreshing and updates) are
    handled in the main loop, which cannot handle window update events while
    the application or callback code is running [1]. This means that nothing
    will happen in the application windows unless we explicitly tell GTK to
    process any events that have been left pending. This is what this context
    manager does, running any pending iteration of the main loop immediately
    before and after the 'with' block.

    [1] PyGTK: FAQ Entry 3.7
    http://faq.pygtk.org/index.py?req=show&file=faq03.007.htp

    """

    def sync():
        # Ensure rendering is done immediately
        while gtk.events_pending():
            gtk.main_iteration()

    sync()
    try:
        yield
    finally:
        sync()


@contextlib.contextmanager
def disable_while(widget):
    """A context manager to temporarily disable a GTK widget.

    Use the gtk.Widget.set_sensitive() method to set the 'sensitive' property
    of the widget to False (therefore disabling it, so that it appears 'grayed
    out' and the user cannot interact with it) when the block is entered and to
    True (enabling it again) after the block is exited. In this manner, we can
    disable a GTK widget while in the 'with block'.

    """

    with gtk_sync():
        widget.set_sensitive(False)

    yield

    with gtk_sync():
        widget.set_sensitive(True)
