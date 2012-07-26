#! /usr/bin/env python
# encoding:UTF-8

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

from __future__ import division

import datetime
import pygtk
pygtk.require ('2.0')
import gtk
import os.path
import sys

# http://stackoverflow.com/a/1054293/184363
path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.append(path)

# LEMON modules
import database
import glade
import methods
import util

class LEMONJuicerGUI(object):

    def __init__ (self, *args, **kwds):
        super(LEMONJuicerGUI, self).__init__(*args, **kwds)

        builder = gtk.Builder()
        builder.add_from_file(glade.GUI_MAIN)
        builder.connect_signals(self)

        self._main_window = builder.get_object('main-window')
        self._notebook    = builder.get_object('main-notebook')
        self._box_toolbar = builder.get_object('box-toolbar')
        self._status_bar  = builder.get_object('status-bar')

        self._main_window.show()
        self._builder = builder


    def run(self):
        gtk.main()

    def handle_close(self, obj):
        self._main_window.destroy()

    def handle_show_about(self, obj):
        self._builder.add_from_file(glade.GUI_ABOUT)
        about = self._builder.get_object('about-dialog')
        about.set_transient_for(self._main_window)
        about.run()
        about.destroy()

    def handle_destroy(self, win):
        gtk.main_quit()

    def handle_open(self, window):
        kwargs = dict(title = None,
                      action = gtk.FILE_CHOOSER_ACTION_OPEN,
                      buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                                 gtk.STOCK_OPEN, gtk.RESPONSE_OK))

        with util.destroying(gtk.FileChooserDialog(**kwargs)) as dialog:

            filt = gtk.FileFilter()
            pattern = '*.LEMONdB'
            filt.set_name('LEMON Database (%s)' % pattern)
            filt.add_mime_type('application/octet-stream; charset=binary')
            filt.add_pattern(pattern)
            dialog.add_filter(filt)

            response = dialog.run()
            if response == gtk.RESPONSE_OK:
                path = dialog.get_filename()
                dialog.destroy()
                self.open_db(path)

    def open_db(self, path):

        self._builder.add_from_file(glade.GUI_OVERVIEW)
        overview = self._builder.get_object('database-overview')
        view = self._builder.get_object('table-view')

        # Display a dialog with a progress bar which is updated as all the
        # stars are loaded into memory, as this may take a while. A 'Cancel'
        # button allows the process to be interrupted at any time.
        self._builder.add_from_file(glade.LOADING_DIALOG)
        dialog = self._builder.get_object('loading-dialog')
        dialog.set_transient_for(self._main_window)
        progressbar = self._builder.get_object('loading-progressbar')

        self._aborted = False
        def cancel(*args):
            self._aborted = True
        dialog.connect('response', cancel)
        dialog.show()

        try:

            self.db = database.LEMONdB(path)
            db_pfilters = self.db.pfilters

            # Two columns are used for the right ascension and declination of
            # the stars, with the sexagesimal and decimal coordinates. They
            # will be shown depending on the radio button selected at View -
            # Coordinates, but the decimal coordinates are always needed, as
            # sexagesimal ones are sorted by them, not lexicographically.

            star_attrs = ['ID', 'α', 'α', 'δ', 'δ', 'm']

            for pfilter in db_pfilters:
                star_attrs.append("Period %s" % pfilter.letter)

            args = [str] * len(self.db.pfilters)
            store = gtk.ListStore(int, str, float, str, float, float, *args)

            for index, attribute in enumerate(star_attrs):
                render = gtk.CellRendererText()
                column = gtk.TreeViewColumn(attribute, render, text = index)
                column.props.resizable = False
                view.append_column(column)

            nstars = len(self.db)
            for star_index, star_id in enumerate(self.db.star_ids):

                # Has the user pressed 'Cancel'?
                if self._aborted:
                    break

                x, y, ra, dec, imag = self.db.get_star(star_id)
                ra_str  = methods.ra_str(ra)
                dec_str = methods.dec_str(dec)
                row = [star_id, ra_str, ra, dec_str, dec, imag]

                for pfilter in db_pfilters:
                    # Returns a two-element tuple, with the period of the star,
                    # in seconds, and the step that the string-length method
                    # used. In case the period is unknown, None is returned.
                    star_period = self.db.get_period(star_id, pfilter)
                    if star_period is None:
                        period_str = ''
                    else:
                        period, step = star_period
                        period_str = datetime.timedelta(seconds = period)
                    row.append(period_str)

                store.append(row)

                fraction = star_index / nstars
                progressbar.set_fraction(fraction)
                # Ensure the progress bar is immediately rendered
                while gtk.events_pending():
                    gtk.main_iteration()

            if not self._aborted:
                label = gtk.Label(os.path.basename(path))
                self._notebook.append_page(overview, label)
                self._notebook.set_tab_reorderable(overview, False)
                self._notebook.set_current_page(-1)

                view.set_model(store)

                title = 'Success'
                msg = "%s stars loaded" % nstars
                util.show_message_dialog(self._main_window, title, msg)

        except Exception, err:

            path = os.path.basename(path)
            title = "Error while loading LEMON database"
            msg = "File '%s' could not be loaded: %s" % (path, str(err))
            util.show_error_dialog(self._main_window, title, msg)

        finally:
            dialog.destroy()
