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

import ConfigParser
import datetime
import functools
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
import config
import database
import glade
import methods
import plot
import snr
import search
import util

# The value which will be inserted in the cells of the gtk.ListStore for which
# there is no data. It must evaluate to False (and, therefore, must be zero) as
# it is used by gtk.TreeViewColumn to determine the visibility of the cells:
# instead of showing a zero, nothing is displayed.
UNKNOWN_VALUE = 0

class StarDetailsGUI(object):
    """ A tabs of the notebook with all the details of a star """

    def update_curve(self, curve, show_airmasses):

        if show_airmasses:
            airmasses = self.db.airmasses(curve.pfilter)
        else:
            airmasses = None

        kwargs = dict(airmasses = airmasses, delta = 3 * 3600,
                      color = self.config.color(curve.pfilter.letter))

        plot.curve_plot(self.figure, curve, **kwargs)
        self.figure.canvas.draw()

    def update_light_curve_points(self, curve):
        """ Update the list of points of the light curve """

        self.curve_store.clear()
        for unix_time, magnitude, noise in curve:
            row = []
            row.append(datetime.datetime.utcfromtimestamp(unix_time).ctime())
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
        for star_id, weight in curve.weights():
            imag = self.db.get_star(star_id)[-1]
            stdev = self.db.get_light_curve(star_id, curve.pfilter).stdev
            args = (star_id, weight, imag, stdev)
            self.refstars_store.append(args)

        # Sort automatically, but respect the user preferences
        if self.refstars_store.get_sort_column_id() == (None, None):
            self.refstars_store.set_sort_column_id(1, gtk.SORT_DESCENDING)

    def handle_toggle_airmasses_checkbox(self, widget):

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
            # Keep track of the filter in the TreeView object, so that
            # LEMONJuicerGUI.handle_row_activated can know in which filter we
            # were working when we clicked on one of the reference stars
            self.refstars_view.pfilter = pfilter

    def handle_toggle_view_sexagesimal(self, *args):
        button = self._builder.get_object('radio-view-sexagesimal')
        for index in self.sex_indexes:
            self.starinfo_store[index][-1] = button.get_active()

    def handle_toggle_view_decimal(self, *args):
        button = self._builder.get_object('radio-view-decimal')
        for index in self.dec_indexes:
            self.starinfo_store[index][-1] = button.get_active()

    def handle_select_period_units(self, button):

        def set_row(index):
            """ Set the visibility of the index-th column of the GtkTreeView,
            depending on whether the button is active or not"""
            self.starinfo_store[index][-1] = button.get_active()

        try:
            if button.get_label() == 'Days':
                for index in self.period_days_indexes:
                    set_row(index)

            elif button.get_label() == 'hh:mm:ss':
                for index in self.period_hhmmss_indexes:
                    set_row(index)

            elif button.get_label() == 'Seconds':
                for index in self.period_seconds_indexes:
                    set_row(index)
            else:
                msg = "unknown button label"
                raise ValueError(msg)

        except AttributeError:
            pass

    def airmasses_visible(self):
        """ Return the state (active or not) of the airmasses checkbox """
        return self.airmasses_checkbox.get_active()

    def __init__(self, parent, star_id, init_pfilter = None):
        """ Instantiate a notebook page with all the star data.

        Keyword arguments:
        init_pfilter - the photometric filter in which to display the data for
                       the star. If not given, the first filter, when sorted by
                       their wavelength, is used.

        """

        # Keep reference to parent LEMONJuicerGUI
        self.parent = parent
        self._builder = parent._builder
        self.config = parent.config
        self.db = parent.db

        self._builder.add_from_file(glade.STAR_DETAILS)
        self.vbox = self._builder.get_object('star-details')
        # Also store the ID of the star in the VBox object; so that we can
        # later loop over the tabs of the notebook and determine to which
        # star each tabs corresponds.
        self.vbox.id = star_id
        self._builder.connect_signals(self.vbox)

        # Use Matplotlib to plot the light curve(s) of the star
        container = self._builder.get_object('image-container')
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
        self.curve_view = self._builder.get_object('curve-points-view')
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

        # GTKTreeView used to display the reference stars and their weight,
        # instrumental magnitude and the standard deviation of their curve
        self.refstars_store = gtk.ListStore(int, float, float, float)
        self.refstars_view = self._builder.get_object('refstars-view')
        for index, title in enumerate(('Star', 'Weight', 'Mag', 'Stdev')):
            render = gtk.CellRendererText()
            column = gtk.TreeViewColumn(title, render, text = index)
            column.props.resizable = False
            column.set_sort_column_id(index)
            self.refstars_view.append_column(column)
        self.refstars_view.set_model(self.refstars_store)
        self.refstars_view.pfilter = None # the filter being currently shown
        self.refstars_view.connect('row-activated', parent.handle_row_activated)

        # GTKTreeView which displays the information for the star
        self.starinfo_store = gtk.ListStore(str, str, bool)
        self.starinfo_view = self._builder.get_object('star-info-view')
        self.starinfo_view.set_headers_visible(False)
        # Column titles not visible; used only for internal reference
        for index, title in enumerate(('Attribute', 'Value', 'Visible')):
            render = gtk.CellRendererText()
            column = gtk.TreeViewColumn(title, render, text = index)
            column.props.resizable = False

            # The third column is only used to determine whether the row
            # is displayed or not, and thus does not need be shown
            if index == 2:
                column.set_visible(False)
            self.starinfo_view.append_column(column)

        store = self.starinfo_store
        x, y, ra, dec, imag = self.db.get_star(star_id)
        store.append(('Right ascension', '%s' % methods.ra_str(ra), True))
        store.append(('Right ascension', '%.4f deg' % ra, True))
        store.append(('Declination', '%s' % methods.dec_str(dec), True))
        store.append(('Declination', '%.4f deg' % dec, True))
        store.append(('Magnitude', '%.3f' % imag, True))
        store.append(('x-coordinate', '%.2f' % x, True))
        store.append(('y-coordinate', '%.2f' % y, True))

        # Two rows (sexagesimal and decimal) are used for the coordinates, and
        # three for each period (days, hh:mm:ss and seconds), but only one will
        # be shown at a time. To hide some of the rows of a TreeView, we need
        # to use a TreeModelFilter, which acts as a wrapper for the TreeModel,
        # allowing you to choose which rows are displayed based on the value of
        # a gobject.TYPE_BOOLEAN column, or based on the output of a certain
        # function. [http://faq.pygtk.org/index.py?file=faq13.048.htp&req=show]

        self.sex_indexes = [0, 2]
        self.dec_indexes = [1, 3]
        self.period_days_indexes = []
        self.period_hhmmss_indexes = []
        self.period_seconds_indexes = []

        for pfilter in self.db.pfilters:
            star_period = self.db.get_period(star_id, pfilter)
            if star_period is not None:
                period, step = star_period
                name = 'Period %s' % pfilter.letter

                store.append((name, "%.4f days" % (period / 3600 / 24), True))
                hhmmss = str(datetime.timedelta(seconds = period))
                store.append((name, "%s hours" % hhmmss, True))
                store.append((name, "%d secs" % period, True))

                length = len(store)
                self.period_days_indexes.append(length - 3)
                self.period_hhmmss_indexes.append(length - 2)
                self.period_seconds_indexes.append(length - 1)

        # A row per filter with the standard deviation of the light curve;
        # these rows are inserted even for the photometric filters in which the
        # curve could not be generated (using UNKNOWN_VALUE instead), but they
        # are hidden from the user.

        for pfilter in self.db.pfilters:
            label = "Stdev %s" % pfilter.letter
            curve = self.db.get_light_curve(star_id, pfilter)
            if curve:
                stdev = curve.stdev if curve else UNKNOWN_VALUE
                store.append((label, stdev, bool(stdev)))

        # Creation of the filter, from the model
        self.starinfo_filter = self.starinfo_store.filter_new()
        self.starinfo_filter.set_visible_column(2)
        self.starinfo_view.set_model(self.starinfo_filter)

        # Show sexagesimal, decimal coordinates or both, depending on
        # what is selected in the View - Coordinates submenu
        args = 'toggled', self.handle_toggle_view_sexagesimal
        self._builder.get_object('radio-view-sexagesimal').connect(*args)
        args = 'toggled', self.handle_toggle_view_decimal
        self._builder.get_object('radio-view-decimal').connect(*args)
        self.handle_toggle_view_sexagesimal()
        self.handle_toggle_view_decimal()

        # Hide all the period rows but one
        buttons = [self._builder.get_object('radio-view-period-days'),
                   self._builder.get_object('radio-view-period-hhmmss'),
                   self._builder.get_object('radio-view-period-seconds')]
        args = 'toggled', self.handle_select_period_units
        for button in buttons:
            button.connect(*args)
            self.handle_select_period_units(button)

        self.shown = None  # pfilter currently shown
        self.id = star_id

        self.buttons = {}  # map each non-empty pfilter to its button
        self.HButtonBox = self._builder.get_object('filter-buttons')
        group = None
        for index, pfilter in enumerate(self.db.pfilters):
            button = gtk.RadioButton(group, label = str(pfilter))
            if not index:
                group = button

            # Make the radio button look like a normal button
            button.set_mode(draw_indicator = False)

            # Disable the button if there is no curve in this filter
            if not self.db.get_light_curve(star_id, pfilter):
                button.set_sensitive(False)
            else:
                button.connect('clicked', self.show_pfilter, pfilter)
                self.buttons[pfilter] = button

            self.HButtonBox.pack_end(button)
            button.show()

        # The checkbox to enable/disable airmasses in the plots, located in
        # the View submenu and shared by all the StarDetailsGUI instances
        self.airmasses_checkbox = self._builder.get_object('plot-airmasses-checkbox')
        args = 'toggled', self.handle_toggle_airmasses_checkbox
        self.airmasses_checkbox.connect(*args)

        # Activate the button of the first filter for which there is data
        # (those in the self.buttons dictionary), unless indicated otherwise
        if not init_pfilter:
            init_pfilter = min(self.buttons.keys())
        self.buttons[init_pfilter].clicked()
        self.shown = init_pfilter


class LEMONJuicerGUI(object):

    # The minimum size (width, height) of the GUI, in pixels
    MIN_SIZE = (1024, 768)

    # The label on the tab for those pages with the details of a star
    TABS_LABEL = "Star %d"

    def __init__(self, *args, **kwds):
        super(LEMONJuicerGUI, self).__init__(*args, **kwds)

        builder = gtk.Builder()
        builder.add_from_file(glade.GUI_MAIN)
        builder.connect_signals(self)

        # The configuration file parser, so that the options selected by the
        # user (such as whether coordinates should be shown in sexagesimal,
        # decimal or both) are persistent.
        self.config = config.Configuration(config.CONFIG_PATH)

        self._main_window = builder.get_object('main-window')
        self._notebook    = builder.get_object('main-notebook')
        self._box_toolbar = builder.get_object('box-toolbar')
        self._status_bar  = builder.get_object('status-bar')
        self.close_button = builder.get_object('close-button')
        self.close_menu_item = builder.get_object('close-menu-item')

        self._main_window.set_size_request(*self.MIN_SIZE)
        self._main_window.show()
        self._builder = builder

        self.db = None  # instance of the LEMONdB class
        # Map the ID of each open star to its StarDetailsGUI instance
        self.open_stars = {}

        # The length, in characters, that all the tabs will have. It cannot be
        # determined until the database is loaded and the number of stars that
        # it has is known.
        self.tabs_length = -1

        # Override the Glade definitions, and set the buttons to the values
        # defined in the configuration file:

        # ================== View submenu ==================
        args = self.config.getboolean, config.VIEW_SECTION
        get_view_booloption = functools.partial(*args)

        checkbox = builder.get_object('radio-view-sexagesimal')
        checkbox.set_active(get_view_booloption(config.VIEW_SEXAGESIMAL))

        checkbox = builder.get_object('radio-view-decimal')
        checkbox.set_active(get_view_booloption(config.VIEW_DECIMAL))

        checkbox = builder.get_object('plot-airmasses-checkbox')
        checkbox.set_active(get_view_booloption(config.PLOT_AIRMASSES))

        # Activate one of the radio buttons (periods expressed in days,
        # hh:mm:ss or seconds) depending on the integer value of the option
        args = config.VIEW_SECTION, config.PERIODS_UNIT
        periods_unit = self.config.getint(*args)
        if periods_unit == config.PERIODS_DAYS:
            name = 'radio-view-period-days'
        elif periods_unit == config.PERIODS_HHMMSS:
            name = 'radio-view-period-hhmmss'
        elif periods_unit == config.PERIODS_SECONDS:
            name = 'radio-view-period-seconds'
        else:
            msg = "invalid value for option '%s'" % periods_unit
            raise ConfigParser.ParsingError(msg)
        builder.get_object(name).set_active(True)

    def save_plot_airmasses_checkbox(self, widget):
        """ Airmasses are not plotted here (that is done StarDetailsGUI), but
        we need to update the configuration file with the new value of the
        option every time this checkbox is toggled """

        checkbox = self._builder.get_object('plot-airmasses-checkbox')
        active = checkbox.get_active()
        value = '1' if active else '0'
        args = config.VIEW_SECTION, config.PLOT_AIRMASSES, value
        self.config.set(*args)

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
            # Update list of open stars and close the page
            star_id = self._notebook.get_nth_page(index).id
            del self.open_stars[star_id]
            self._notebook.remove_page(index)

        # Disable the close button / menu item if the notebook becomes empty
        npages_left = self._notebook.get_n_pages()
        self.close_button.set_sensitive(npages_left)
        self.close_menu_item.set_sensitive(npages_left)

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

        # Update the configuration file with the new value of the option (not
        # actually written to disk until the execution of the program finishes)
        args = config.VIEW_SECTION, config.VIEW_SEXAGESIMAL, str(active and 1 or 0)
        self.config.set(*args)

        try:
            self.view.get_column( self.ra_sex_index).set_visible(active)
            self.view.get_column(self.dec_sex_index).set_visible(active)
        # Raised if no LEMONdB has yet been loaded
        except AttributeError:
            pass

    def handle_toggle_view_decimal(self, *args):
        button = self._builder.get_object('radio-view-decimal')
        active = button.get_active()

        # Update the configuration file with the new value of the option
        args = config.VIEW_SECTION, config.VIEW_DECIMAL, str(active and 1 or 0)
        self.config.set(*args)

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

    def save_periods_unit_radio_item(self, button):
        """ Update the configuration file with the new value of the option
         every time a radio item with the unit of the periods is selected"""

        # This function gets called twice every time a new option is selected,
        # as two different buttons are being toggled. We are only interesed in
        # the one which has been activated, though.
        if not button.get_active():
            return

        if button.get_label() == 'Days':
            option = config.PERIODS_DAYS

        elif button.get_label() == 'hh:mm:ss':
            option = config.PERIODS_HHMMSS

        elif button.get_label() == 'Seconds':
            option = config.PERIODS_SECONDS

        else:
            msg = "unknown button label"
            raise ValueError(msg)

        # Update the configuration file with the new value of the option
        args = config.VIEW_SECTION, config.PERIODS_UNIT, str(option)
        self.config.set(*args)

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

            # Two additional columns for each photometric filter, with (a) the
            # standard deviation of the points of each light curve and (b) a
            # boolean value indicating whether (a) must be shown. The latter
            # column, which is never visible, is required because standard
            # deviations are real values and therefore PyGTK cannot set the
            # 'visible' attribute from them, even if they contain
            # UNKNOWN_VALUE.

            self.stdevs_indexes = []
            self.stdevs_visibility_indexes = []

            for pfilter in db_pfilters:
                label = "Stdev %s" % pfilter.letter
                star_attrs.append(label)
                star_attrs.append("Visible (%s)" % label)
                length = len(star_attrs)
                self.stdevs_indexes.append(length - 2)
                self.stdevs_visibility_indexes.append(length - 1)

            args = [int, str, float, str, float, float]
            args += [float, str, int] * len(db.pfilters)
            args += [float, bool] * len(db.pfilters)
            self.store = gtk.ListStore(*args)

            for index, attribute in enumerate(star_attrs):
                render = gtk.CellRendererText()

                # 'text': column from which the text of the column is taken
                # 'visible': column which determines the visibility of the cell
                kwargs = dict(text = index)

                # Determine the column by which this column is sorted
                if index == self.ra_sex_index:
                    sort_index = self.ra_dec_index
                elif index == self.dec_sex_index:
                    sort_index = self.dec_dec_index

                # Periods are in the following order: days, hh:mm:ss, seconds.
                # We want to sort hh:mm:ss periods not lexicographically, but
                # by the value of the third column, which stores seconds.

                # The first of these columns also be sorted by the period in
                # days (a real number) but then, when setting the 'visible'
                # attribute below, we would get the "unable to set property
                # `visible' of type `gboolean' from value of type `gdouble'"
                # error. To avoid this, we sort all periods and determine their
                # visibility using the columns which contain them in seconds.
                elif index in self.period_days_indexes:
                    kwargs['visible'] = sort_index = index + 2

                elif index in self.period_hhmmss_indexes:
                    kwargs['visible'] = sort_index = index + 1

                elif index in self.period_seconds_indexes:
                    kwargs['visible'] = sort_index = index

                # The cells with standard deviations are only displayed if the
                # value in the adjacent column, which stores whether the light
                # curve could be generated, evaluates to True.
                elif index in self.stdevs_indexes:
                    sort_index = index
                    kwargs['visible'] = index + 1

                else:
                    sort_index = index

                column = gtk.TreeViewColumn(attribute, render, **kwargs)
                column.props.resizable = False
                column.set_sort_column_id(sort_index)

                # Hide columns which determine the visibility of the stdevs
                if index in self.stdevs_visibility_indexes:
                    column.set_visible(False)

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
                        row += [UNKNOWN_VALUE, str(UNKNOWN_VALUE), UNKNOWN_VALUE]
                    else:
                        period, step = star_period
                        row.append(period / 3600 / 24) # days
                        row.append(str(datetime.timedelta(seconds = period)))
                        row.append(period) # seconds

                for pfilter in db_pfilters:
                    # None returned if the star doesn't have this light curve
                    curve = db.get_light_curve(star_id, pfilter)
                    if curve:
                        row += [curve.stdev, True]
                    else:
                        row += [UNKNOWN_VALUE, False]

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

                self.view.set_model(self.store)

                label = gtk.Label(os.path.basename(path))
                self._notebook.append_page(overview, label)
                self._notebook.set_tab_reorderable(overview, False)
                self._notebook.set_current_page(-1)

                # Now that there is at least one page in the notebook, make the
                # 'Close' button and menu items sensitive, so that the user can
                # interact with them.
                self.close_button.set_sensitive(True)
                self.close_menu_item.set_sensitive(True)

        except Exception, err:
            path = os.path.basename(path)
            title = "Error while loading LEMON database"
            msg = "File '%s' could not be loaded: %s" % (path, str(err))
            util.show_error_dialog(self._main_window, title, msg)

        finally:
            dialog.destroy()

    def show_star(self, star_id, pfilter = None):
        """ Append a page with the star whose ID is 'star_id'.
        Returns the StarDetailsGUI instance which encapsulates all the
        information of the star.

        Keyword arguments:
        pfilter - the photometric filter in which information must be displayed
                  when the StarDetailsGUI is created. If not given, data will
                  be shown for the first filter (when sorted by wavelength).

        """

        details = StarDetailsGUI(self, star_id, init_pfilter = pfilter)

        tab_label = gtk.Label(self.TABS_LABEL % star_id)
        tab_label.set_width_chars(self.tabs_length)

        self._notebook.append_page(details.vbox, tab_label)
        self._notebook.set_tab_reorderable(details.vbox, False)
        self._notebook.set_current_page(-1)
        return details

    def get_star_index(self, star_id):
        """ Return the page in the notebook with the star whose ID is 'star_id'.
        Returns -1 if none of the pages of the notebook contains that star """

        for index, page in enumerate(self._notebook):
            try:
                if page.id == star_id:
                    return index
            # The main tab (the list of stars) doesn't have an ID
            except AttributeError:
                pass
        else:
            return -1

    def switch_to_tab(self, star_id, pfilter = None):
        """ Switch to the page with the star whose ID is 'star_id'.

        Keyword arguments:
        pfilter - after switching to the page of the star, switch also to this
                  photometric filter. If not given, the filter will not be
                  modified, so it will be the last selected by the user. The
                  KeyError exception is raised if the star has no data for the
                  specified photometric filter.

        """

        index = self.get_star_index(star_id)
        if index != -1:
            self._notebook.set_current_page(index)
            # Do we have to move to a specific filter after switching?
            if pfilter:
                self.open_stars[star_id].buttons[pfilter].clicked()
        else:
            msg = "star %d is not being shown" % star_id
            raise ValueError(msg)

    def handle_row_activated(self, view, row, column, id_index = 0):
        """ Open the details of the star, or switch to them if already open.

        Keyword arguments:

        id_index - the column of the model to which the view is associated
        which contains the ID of the star.

        """

        # Determine the ID of the star on which the user has clicked
        star_id = view.get_model()[row][id_index]

        # If we clicked on the star not in the main window (with all the stars
        # in the database) but in a list of reference stars, we want to open
        # the star (or switch to it) in the same filter from which we come.
        try:
            pfilter = view.pfilter
        except AttributeError:
            pfilter = None

        if star_id in self.open_stars.iterkeys():
            self.switch_to_tab(star_id, pfilter = pfilter)

        else:
            details = self.show_star(star_id, pfilter = pfilter)
            self.open_stars[star_id] = details


    def search_by_amplitudes(self, window):
        """ Identify stars with amplitudes correlated to the wavelength. These
        are listed in a gtk.ScrolledWindow which is appended to the notebook"""

        result = search.amplitudes_search(self._main_window, self._builder, self.db.path)
        if result is not None:
            label = gtk.Label(result.get_label())
            window = result.get_window()
            self._notebook.append_page(window, label)
            self._notebook.set_tab_reorderable(window, False)
            self._notebook.set_current_page(-1)

