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

import pygtk
pygtk.require ('2.0')
import gtk

import datetime
import os.path
import sys
import time

import matplotlib.figure
from matplotlib.backends.backend_gtkagg \
     import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg \
     import NavigationToolbar2GTKAgg as NavigationToolbar

# http://stackoverflow.com/a/1054293/184363
path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.append(path)

# LEMON modules
import database
import glade
import methods
import plot
import snr
import util

class StarDetailsGUI(object):
    """ A tabs of the notebook with all the details of a star """

    def update_curve(self, curve, show_airmasses):

        if show_airmasses:
            airmasses = self.db.airmasses(curve.pfilter)
        else:
            airmasses = None

        kwargs = dict(airmasses = airmasses,
                      delta = self.airmasses_interval())
        plot.curve_plot(self.figure, curve, **kwargs)
        self.figure.canvas.draw()

    def update_light_curve_points(self, curve):
        """ Update the list of points of the light curve """

        self.curve_store.clear()
        for unix_time, magnitude, noise in curve:
            row = []
            row.append(time.ctime(unix_time))
            row.append(unix_time)
            row.append(magnitude)
            row.append(noise)
            # Returns two errors in mags, positive and negative
            nmerr, pmerr = snr.snr_to_error(noise)
            row.append(pmerr)
            row.append(nmerr)
            self.curve_store.append(row)

        # Automatically change the order after showing the model (without
        # the user having to click on the column header), unless the user
        # has selected a diferent column or order.
        if self.curve_store.get_sort_column_id() == (None, None):
            self.curve_store.set_sort_column_id(1, gtk.SORT_ASCENDING)

    def update_reference_stars(self, curve):
        """ Update the list of reference stars """

        self.refstars_store.clear()
        for weight in curve.weights():
            self.refstars_store.append(weight)

        # Sort automatically, but respect the user preferences
        if self.refstars_store.get_sort_column_id() == (None, None):
            self.refstars_store.set_sort_column_id(1, gtk.SORT_DESCENDING)

    def handle_toggle_airmasses_checkbox(self, widget):

        # Disable (gray out) the spin button, and therefore don't let the user
        # interact with it, if airmasses are unchecked and won't be plotted
        self.airmasses_interval_spinbutton.set_sensitive(self.airmasses_visible())

        # Replot the light curve again, with or without airmasses
        curve = self.db.get_light_curve(self.id, self.shown)
        self.update_curve(curve, self.airmasses_visible())

    def show_pfilter(self, widget, pfilter):
        """ Display the information of the star in this photometric filter """

        # 'show' is initialized to NoneType; so this method will
        # be always executed the first time it is called
        if not self.shown or self.shown != pfilter:

            self.shown = pfilter
            curve = self.db.get_light_curve(self.id, pfilter)
            self.update_curve(curve, self.airmasses_visible())
            self.update_light_curve_points(curve)
            self.update_reference_stars(curve)

    def airmasses_visible(self):
        """ Return the state (active or not) of the airmasses checkbox """
        return self.airmasses_checkbox.get_active()

    def airmasses_interval(self):
        """ Return the value of the airmasses spin button, in seconds """
        # Hours to seconds conversion
        return 3600 * self.airmasses_interval_spinbutton.get_value()

    def __init__(self, db, star_id):

        builder = gtk.Builder()
        builder.add_from_file(glade.STAR_DETAILS)
        self.vbox = builder.get_object('star-details')
        # Also store the ID of the star in the VBox object; so that we can
        # later loop over the tabs of the notebook and determine to which
        # star each tabs corresponds.
        self.vbox.id = star_id
        builder.connect_signals(self.vbox)

        # Use Matplotlib to plot the light curve(s) of the star
        container = builder.get_object('image-container')
        self.figure = matplotlib.figure.Figure()
        canvas = FigureCanvas(self.figure)
        navig = NavigationToolbar(canvas, container)

        container.pack_start(canvas)
        container.pack_start(navig, False)
        container.show_all()

        # GTKTreeView used to display the list of points of the curve; dates
        # are plotted twice: hh:mm:ss and also in Unix time, the latter of
        # which is used to sort the columns by their date.
        attrs = ('Date', 'Date', 'Mag', 'SNR', 'merr (+)', 'merr (-)')
        self.curve_store = gtk.ListStore(str, float, float, float, float, float)
        self.curve_view = builder.get_object('curve-points-view')
        for index, title in enumerate(attrs):
            render = gtk.CellRendererText()
            column = gtk.TreeViewColumn(title, render, text = index)
            column.props.resizable = False
            # The first column (index = 0) is sorted by the second
            column.set_sort_column_id(1 if not index else index)

            # The column with dates in Unix time is not shown
            if index == 1:
                column.set_visible(False)

            self.curve_view.append_column(column)
        self.curve_view.set_model(self.curve_store)

        # GTKTreeView used to display the reference stars and their weight
        self.refstars_store = gtk.ListStore(int, float)
        self.refstars_view = builder.get_object('refstars-view')
        for index, title in enumerate(('Star', 'Weight')):
            render = gtk.CellRendererText()
            column = gtk.TreeViewColumn(title, render, text = index)
            column.props.resizable = False
            column.set_sort_column_id(index)
            self.refstars_view.append_column(column)
        self.refstars_view.set_model(self.refstars_store)

        self.shown = None  # pfilter currently shown
        self.db = db
        self.id = star_id

        x, y, ra, dec, imag = db.get_star(star_id)
        # TODO: Display star information

        self.buttons = {}  # map each pfilter to its button
        self.HButtonBox = builder.get_object('filter-buttons')
        group = None
        for index, pfilter in enumerate(db.pfilters):
            button = gtk.RadioButton(group, label = str(pfilter))
            if not index:
                group = button

            # Make the radio button look like a normal button
            button.set_mode(draw_indicator = False)
            button.connect('clicked', self.show_pfilter, pfilter)

            self.buttons[pfilter] = button
            self.HButtonBox.pack_end(button)
            button.show()

        # The checkbox to enable-disable airmasses in the plots...
        self.airmasses_checkbox = builder.get_object('show-airmasses-checkbox')
        args = 'toggled', self.handle_toggle_airmasses_checkbox
        self.airmasses_checkbox.connect(*args)

        # ... and the spin button which allows us to navigate through a range
        # of values in order to select the maximum distance between points in
        # the light curve if their airmasses are to be plotted connected.
        arg = 'airmass-interval-spinbutton'
        self.airmasses_interval_spinbutton = builder.get_object(arg)
        # https://mail.gnome.org/archives/gtk-app-devel-list/2010-April/msg00127.html
        self.airmasses_interval_spinbutton.set_value(3)
        args = 'output', self.handle_toggle_airmasses_checkbox
        self.airmasses_interval_spinbutton.connect(*args)


        # Activate the button of the first filter
        pfilter = min(self.buttons.keys())
        self.buttons[pfilter].clicked()
        self.shown = pfilter


class LEMONJuicerGUI(object):

    # The label on the tab for those pages with the details of a star
    TABS_LABEL = "Star %d"

    def __init__(self, *args, **kwds):
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

        self.db = None  # instance of the LEMONdB class
        # Map the ID of each open star to its StarDetailsGUI instance
        self.open_stars = {}

        # The length, in characters, that all the tabs will have. It cannot be
        # determined until the database is loaded and the number of stars that
        # it has is known.
        self.tabs_length = -1

    def run(self):
        gtk.main()

    def handle_close(self, obj):
        # If the notebook has no pages, -1 is returned
        index = self._notebook.get_current_page()

        # We're closing the list of stars!
        if index == 0:

            really_close = True

            # Warn the user that all open tabs, if any, will be lost
            if self._notebook.get_n_pages() > 1:
                title = "Are you done with this database?"
                msg = ("If you close the LEMONdB all the open tabs will "
                       "also be closed. Do you really want to proceed?")

                args = self._main_window, title, msg
                kwargs = dict(msg_type = gtk.MESSAGE_QUESTION,
                              buttons = gtk.BUTTONS_OK_CANCEL)
                response = util.show_message_dialog(*args, **kwargs)
                really_close = response == gtk.RESPONSE_OK

            if really_close:
                self.db = None # forget about this LEMONdB
                while self._notebook.get_n_pages():
                    self._notebook.remove_page(-1)

        elif index >= 1:
            self._notebook.remove_page(index)

    def handle_quit(self, obj):
        self._main_window.destroy()

    def handle_show_about(self, obj):
        self._builder.add_from_file(glade.GUI_ABOUT)
        about = self._builder.get_object('about-dialog')
        about.set_transient_for(self._main_window)
        about.run()
        about.destroy()

    def handle_destroy(self, win):
        gtk.main_quit()

    def handle_toggle_view_sexagesimal(self, *args):
        button = self._builder.get_object('radio-view-sexagesimal')
        active = button.get_active()
        try:
            self.view.get_column( self.ra_sex_index).set_visible(active)
            self.view.get_column(self.dec_sex_index).set_visible(active)
        # Raised if no LEMONdB has yet been loaded
        except AttributeError:
            pass

    def handle_toggle_view_decimal(self, *args):
        button = self._builder.get_object('radio-view-decimal')
        active = button.get_active()
        try:
            self.view.get_column( self.ra_dec_index).set_visible(active)
            self.view.get_column(self.dec_dec_index).set_visible(active)
        except AttributeError:
            pass

    def handle_select_period_units(self, button):

        def set_column(index):
            """ Set the visibility of the index-th column of the GtkTreeView,
            depending on whether the button is active or not"""
            self.view.get_column(index).set_visible(button.get_active())

        try:
            if button.get_label() == 'Days':
                for index in self.period_days_indexes:
                    set_column(index)

            elif button.get_label() == 'hh:mm:ss':
                for index in self.period_hhmmss_indexes:
                    set_column(index)

            elif button.get_label() == 'Seconds':
                for index in self.period_seconds_indexes:
                    set_column(index)
            else:
                msg = "unknown button label"
                raise ValueError(msg)

        except AttributeError:
            pass

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

        # If this is not the first LEMONdB that is opened, warn the user that
        # the current one will be closed, unless the operation is aborted.
        if self.db is not None:
            title = "Are you done with this database?"
            msg = ("Unfortunately, data analysis must be done on one LEMONdB "
                   "at a time, so opening a new one will close all your "
                   "current tabs. Do you really want to proceed?")
            args = self._main_window, title, msg
            kwargs = dict(msg_type = gtk.MESSAGE_QUESTION,
                          buttons = gtk.BUTTONS_OK_CANCEL)
            response = util.show_message_dialog(*args, **kwargs)

            # If the user presses OK, close all the tabs of the notebook.
            # Otherwise, the load is cancelled and nothing is modified.
            if response == gtk.RESPONSE_OK:
                while self._notebook.get_n_pages():
                    self._notebook.remove_page(-1)
            else:
                return

        self._builder.add_from_file(glade.GUI_OVERVIEW)
        overview = self._builder.get_object('database-overview')
        self.view = self._builder.get_object('table-view')
        self.view.connect('row-activated', self.handle_row_activated)

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

            db = database.LEMONdB(path)
            db_pfilters = db.pfilters

            # Two columns are used for the right ascension and declination of
            # the stars, with the sexagesimal and decimal coordinates. They
            # will be shown depending on the radio button selected at View -
            # Coordinates, but the decimal coordinates are always needed, as
            # sexagesimal ones are sorted by them, not lexicographically.

            star_attrs = ['ID', 'α', 'α', 'δ', 'δ', 'm']

            assert star_attrs.count('ID') == 1
            assert star_attrs.count('α') == 2
            assert star_attrs.count('δ') == 2
            self.id_index = star_attrs.index('ID')
            self.ra_sex_index = star_attrs.index('α')
            self.ra_dec_index = self.ra_sex_index + 1
            self.dec_sex_index = star_attrs.index('δ')
            self.dec_dec_index = self.dec_sex_index + 1

            # Three columns are used for each period, to show it in days (a
            # real number), hh:mm:ss (str) or seconds (int). Periods will be
            # sorted always by their value in seconds. We need to keep an
            # internal list with the indexes of each type of column, to make
            # them visible or not depending on the option selected at View -
            # Periods.

            self.period_days_indexes = []
            self.period_hhmmss_indexes = []
            self.period_seconds_indexes = []

            for pfilter in db_pfilters:
                label = "Period %s" % pfilter.letter
                star_attrs += [label] * 3
                length = len(star_attrs)
                self.period_days_indexes.append(length - 3)
                self.period_hhmmss_indexes.append(length - 2)
                self.period_seconds_indexes.append(length - 1)

            args = [float, str, int] * len(db.pfilters)
            self.store = gtk.ListStore(int, str, float, str, float, float, *args)

            for index, attribute in enumerate(star_attrs):
                render = gtk.CellRendererText()
                column = gtk.TreeViewColumn(attribute, render, text = index)
                column.props.resizable = False

                # Set the current sort comparison function of the column
                if index == self.ra_sex_index:
                    sort_index = self.ra_dec_index
                elif index == self.dec_sex_index:
                    sort_index = self.dec_dec_index

                # Periods are in the following order: days, hh:mm:ss, seconds.
                # We want to sort hh:mm:ss periods not lexicographically, but
                # by the value of the third column, which stores seconds.
                elif index in self.period_hhmmss_indexes:
                    sort_index = index + 1
                else:
                    sort_index = index
                column.set_sort_column_id(sort_index)

                self.view.append_column(column)

            nstars = len(db)
            for star_index, star_id in enumerate(db.star_ids):

                # Has the user pressed 'Cancel'?
                if self._aborted:
                    break

                x, y, ra, dec, imag = db.get_star(star_id)
                ra_str  = methods.ra_str(ra)
                dec_str = methods.dec_str(dec)
                row = [star_id, ra_str, ra, dec_str, dec, imag]

                for pfilter in db_pfilters:
                    # Returns a two-element tuple, with the period of the star,
                    # in seconds, and the step that the string-length method
                    # used. In case the period is unknown, None is returned.
                    star_period = db.get_period(star_id, pfilter)
                    if star_period is None:
                        # Use -1 for unknown values
                        row += [-1, '-1', -1]
                    else:
                        period, step = star_period
                        row.append(period / 3600 / 24) # hours
                        row.append(str(datetime.timedelta(seconds = period)))
                        row.append(period) # seconds

                self.store.append(row)

                # Update the progress bar only when the percentage varies; if,
                # e.g., it is 0.971 (97%), setting the fraction to 0.972 (still
                # 97%) would only unnecessarily slow down the execution.
                fraction = round(star_index / nstars, 2)
                if fraction != progressbar.get_fraction():
                    progressbar.set_fraction(fraction)
                    # Ensure rendering is done immediately
                    while gtk.events_pending():
                        gtk.main_iteration()

            if not self._aborted:

                # Save the LEMONdB for later use now that the load has been
                # completed; we could not do it earlier as it could have been
                # aborted. Also clear the set of IDs of open stars.
                self.db = db
                self.open_stars.clear()

                # Find the maximum number of characters needed for the label on
                # the tabs; that is, find the length of the tab for the star in
                # the database with the highest ID (these are assumed to be
                # positive integers, because of the AUTOINCREMENT keyword, so
                # this ID must be also that with the maximum number of digits).
                self.tabs_length = len(self.TABS_LABEL % max(db.star_ids))

                # Show sexagesimal, decimal coordinates or both
                self.handle_toggle_view_sexagesimal()
                self.handle_toggle_view_decimal()

                # Hide all the period columns but one
                buttons = [self._builder.get_object('radio-view-period-days'),
                           self._builder.get_object('radio-view-period-hhmmss'),
                           self._builder.get_object('radio-view-period-seconds')]
                for button in buttons:
                    self.handle_select_period_units(button)


                label = gtk.Label(os.path.basename(path))
                self._notebook.append_page(overview, label)
                self._notebook.set_tab_reorderable(overview, False)
                self._notebook.set_current_page(-1)

                self.view.set_model(self.store)

        except Exception, err:
            path = os.path.basename(path)
            title = "Error while loading LEMON database"
            msg = "File '%s' could not be loaded: %s" % (path, str(err))
            util.show_error_dialog(self._main_window, title, msg)

        finally:
            dialog.destroy()

    def show_star(self, star_id):
        """ Append a page with the star whose ID is 'star_id'.
        Returns the StarDetailsGUI instance which encapsulates all the
        informacion of the star"""

        details = StarDetailsGUI(self.db, star_id)
        details.refstars_view.connect('row-activated', self.handle_row_activated)

        tab_label = gtk.Label(self.TABS_LABEL % star_id)
        tab_label.set_width_chars(self.tabs_length)

        self._notebook.append_page(details.vbox, tab_label)
        self._notebook.set_tab_reorderable(details.vbox, False)
        self._notebook.set_current_page(-1)
        return details

    def switch_to_tab(self, star_id):
        """ Switch to the page with the star whose ID is 'star_id' """

        for index, page in enumerate(self._notebook):
            try:
                if page.id == star_id:
                    self._notebook.set_current_page(index)
                    break
            # The main tab (the list of stars) doesn't have an ID
            except AttributeError:
                pass
        else:
            msg = "star %d is not being shown" % star_id
            raise ValueError(msg)

    def handle_row_activated(self, view, row, column, id_index = 0):
        """ Open the details of the star, or switch to them if already open.

        Keyword arguments:

        id_index - the column of the model to which the view is associated
        which contains the ID of the star.

        """

        # determine the ID of the star on which the user has clicked
        star_id = view.get_model()[row][id_index]

        if star_id in self.open_stars.iterkeys():
            self.switch_to_tab(star_id)

        else:
            details = self.show_star(star_id)
            self.open_stars[star_id] = details

