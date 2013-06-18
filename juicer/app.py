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
from __future__ import with_statement

import pygtk
pygtk.require ('2.0')
import gtk

import ConfigParser
import datetime
import functools
import operator
import os.path
import re
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
import chart
import config
import glade
import methods
import mining
import plot
import snr
import search
import util

# Make sure that stock icons are shown in buttons, as with the default options
# this is not always the case. I have not been able to find out whether this is
# caused by GTK+ or GNOME (see http://stackoverflow.com/questions/2188659/ and
# http://www.daa.com.au/pipermail/pygtk/2009-November/017826.html), but in any
# case it can be fixed without having to change any system configuration files.
settings = gtk.settings_get_default()
settings.props.gtk_button_images = True

# The value which will be inserted in the cells of the gtk.ListStore for which
# there is no data. It must evaluate to False (and, therefore, must be zero) as
# it is used by gtk.TreeViewColumn to determine the visibility of the cells:
# instead of showing a zero, nothing is displayed.
UNKNOWN_VALUE = 0

class ExportCurveDialog(object):
    """ Encapsulates a dialog window that allows the user to select which
    attributes of a light curve (such as the time of observation in Unix
    seconds, or the signal-to-noise ratio, or the error in magnitudes) are
    dumped to a text file, selected using a gtk.FileChooserDialog """

    def get(self, name):
        """ Access a widget in the interface """
        return self.builder.get_object(name)

    def __init__(self, parent_window, builder, config, id_, pfilter,
                 field_name, curve_store):
        """ Instantiation method for the ExportCurveDialog class.

        The 'parent_window' parameter must be the transient parent of the
        dialog, while 'builder' and 'config' are the gtk.GtkBuilder and
        Configuration instances, respectively, of the parent GTK widget.
        The 'id_' parameter is the ID of the star, 'pfilter' a Passband
        instance encapsulating the photometric filter of the light curve
        and 'field_name' a string containing the name of the observed
        field, used to suggest a filename in a 'Save As...' dialog.

        Lastly, 'curve_store' must be a gtk.ListStore with the data that will
        be dumped to a file, and should contain six columns: (1) a textual
        representation of the date of observation, (2) the date of observation
        in seconds after the Unix epoch, (3) the differential magnitude, (4)
        the signal-to-noise ratio, (5 and 6) the maximum and minimum errors
        induced by the noise, respectively. The columns must be in this exact
        order and of type float, except for the first one, which is a str.

        """

        self.parent_window = parent_window
        self.builder = builder
        self.builder.add_from_file(glade.EXPORT_CURVE_DIALOG)
        self.config = config
        self.id = id_
        self.pfilter = pfilter
        self.field_name = field_name
        self.store = curve_store

        self.dialog = self.get('export-curve-dialog')
        self.dialog.set_resizable(False)
        self.dialog.set_title("Export to text file")
        self.dialog.add_button(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL)
        self.dialog.add_button(gtk.STOCK_SAVE, gtk.RESPONSE_OK)

        self.dialog.set_transient_for(self.parent_window)
        self.dialog.set_position(gtk.WIN_POS_CENTER_ON_PARENT)

        text = \
        "Please select the attributes of the light curve of <b>star %d</b> " \
        "in <b>filter %s</b> that will be saved to disk. The tab character, " \
        "'\\t', will be used as separator." % (self.id, self.pfilter)
        self.get('dialog-description').set_label(text)

        self.date_str_checkbox = self.get('date-str-checkbox')
        self.date_secs_checkbox = self.get('date-secs-checkbox')
        self.mags_checkbox = self.get('mags-checkbox')
        self.snr_checkbox = self.get('snr-checkbox')
        self.merr_pos_checkbox = self.get('merr-pos-checkbox')
        self.merr_neg_checkbox = self.get('merr-neg-checkbox')
        self.spinbutton = self.get('ndecimals-spinbutton')
        self.update()

        # Handlers that update the options in the configuration file every time
        # that the value of a widget (either check or spin button) is modified
        def save_widget_update(option, func = 'get_active'):
            """ Return the function that, when called, updates 'option' in the
            curve export section of the configuration file. 'func' is the
            method used to get the value of the widget, cast to integer """

            def handler(widget):
                self.config.dumpset(option, int(getattr(widget, func)()))
            return handler

        f = save_widget_update
        self.date_str_checkbox.connect('toggled', f('dump_date_text'))
        self.date_secs_checkbox.connect('toggled', f('dump_date_seconds'))
        self.mags_checkbox.connect('toggled', f('dump_magnitude'))
        self.snr_checkbox.connect('toggled', f('dump_snr'))
        self.merr_pos_checkbox.connect('toggled', f('dump_max_merr'))
        self.merr_neg_checkbox.connect('toggled', f('dump_min_merr'))
        self.spinbutton.connect('output', f('decimal_places', 'get_value'))

    def update(self):
        """ Update the buttons to the values given in the configuration file """

        config = self.config
        self.date_str_checkbox.set_active(config.dumpint('dump_date_text'))
        self.date_secs_checkbox.set_active(config.dumpint('dump_date_seconds'))
        self.mags_checkbox.set_active(config.dumpint('dump_magnitude'))
        self.snr_checkbox.set_active(config.dumpint('dump_snr'))
        self.merr_pos_checkbox.set_active(config.dumpint('dump_max_merr'))
        self.merr_neg_checkbox.set_active(config.dumpint('dump_min_merr'))
        self.spinbutton.set_value(config.dumpint('decimal_places'))

    def dump(self, path, separator = '\t'):
        """ Save a light curve to the plain text file 'path'.

        Iterate over the gtk.ListStore (which is expected, although this is not
        enforced, to be sorted chronologically; that is, sorted by the value of
        the second column, which contains the date of observation in Unix time)
        and save a textual representation of its rows to 'path', truncating the
        file if it already exists. Not all the rows of the gtk.ListStore (i.e.,
        the attributes of the light curve, such as the signal-to-noise ratio or
        the maximum and minimum error induced by the noise) are saved to the
        file, but only those for which the corresponding checkbox is active
        (checked).

        """

        def parse_float(value):
            """ Cast value to str; use exactly 'decimals' decimal digits """
            ndecimals = int(self.spinbutton.get_value())
            return '%.*f' % (ndecimals, value)

        with open(path, 'wt') as fd:

            for row in self.store:

                values = []
                if self.date_str_checkbox.get_active():
                    values.append(row[0])
                if self.date_secs_checkbox.get_active():
                    values.append(parse_float(row[1]))
                if self.mags_checkbox.get_active():
                    values.append(parse_float(row[2]))
                if self.snr_checkbox.get_active():
                    values.append(parse_float(row[3]))
                if self.merr_pos_checkbox.get_active():
                    values.append(parse_float(row[4]))
                if self.merr_neg_checkbox.get_active():
                    values.append(parse_float(row[5]))

                fd.write('%s\n' % separator.join(values))

    def run(self):
        """ Run the dialog window in a recursive loop.

        This method shows the dialog window and allows the user to select,
        using checkboxes, which attributes of the light curve of the star will
        be dumped. Then, upon clicking 'Save', a gtk.FileChooserDialog lets the
        user browse for the location where the plain text file will be saved.

        """

        response = self.dialog.run()
        if response == gtk.RESPONSE_OK:

            kwargs = dict(title = "Export light curve to...",
                          parent = self.parent_window,
                          action = gtk.FILE_CHOOSER_ACTION_SAVE,
                          buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                                     gtk.STOCK_SAVE, gtk.RESPONSE_OK))

            with util.destroying(gtk.FileChooserDialog(**kwargs)) as chooser:

                # Ask the user to confirm overwriting an existing file
                chooser.set_do_overwrite_confirmation(True)

                # Suggest a name for the plain text file
                field = re.sub(r'\s', '_', self.field_name.lower())
                args = field, self.id, self.pfilter
                filename = '%s_star_%d_curve_%s' % args
                chooser.set_current_name(filename)
                response = chooser.run()

                if response == gtk.RESPONSE_OK:
                    self.dump(chooser.get_filename())

        self.dialog.destroy()


class StarDetailsGUI(object):
    """ A tabs of the notebook with all the details of a star """

    def set_canvas(self, visible):
        """ Make rows in the 'matplotlib-container' VBox visible/invisible.

        Adjust the visibility of the three rows in the vertical box: if
        'visible' is True, the first two rows (canvas and navigation toolbar)
        are made visible, while the error message (stating that no point is
        above the SNR threshold) is made invisible. Their visibility is
        reversed if it is False: only the error message can be been.

        """

        self.image_box.set_visible(visible)
        self.navigation_box.set_visible(visible)
        self.error_msg.set_visible(not visible)

    def update_curve(self, curve, show_airmasses):

        if show_airmasses:
            airmasses = self.db.airmasses(curve.pfilter)
        else:
            airmasses = None

        threshold = self.config.get_minimum_snr()
        curve = curve.ignore_noisy(threshold)

        if not curve:
            # Display error message: no points of the light curve are above
            # the SNR threshold. Note that 'curve' can only be None for this
            # reason: filters for which there is no light curve are disabled
            # at __init__, so the user cannot plot anything in these filters.
            self.set_canvas(False)

        else:
            self.set_canvas(True)
            kwargs = dict(airmasses = airmasses, delta = 3 * 3600,
                          color = self.config.color(curve.pfilter.letter))

            plot.curve_plot(self.figure, curve, **kwargs)
            self.figure.canvas.draw()

    def update_light_curve_points(self, curve):
        """ Update the list of points of the light curve """

        self.curve_store.clear()
        for unix_time, magnitude, noise in curve:
            row = []
            row.append(methods.utctime(unix_time, suffix = False))
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

    def redraw_light_curve(self, widget):
        """ Replot the light curve """

        curve = self.db.get_light_curve(self.id, self.shown)
        self.update_curve(curve, self.airmasses_visible())

    def show_pfilter(self, button, event, pfilter):
        """ Display the information of the star in this photometric filter """


        # Ignore clicks on already-pressed buttons (i.e., do not draw the same
        # light curve twice in a row). Note that 'show' is initialized to None,
        # so this method will be always executed the first time it is called.
        if self.shown and self.shown == pfilter:
            return

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

        # We use a three-row vertical box to pack the Matplotlib canvas (used
        # to plot the light curve(s) of the star) and its navigation toolbar.
        # The last row contains a label with the error message that is shown
        # (after making invisible the other two rows) when there is no light
        # curve to plot because all the points are below the SNR threshold.

        matplotlib_container = self._builder.get_object('matplotlib-container')
        self.image_box = self._builder.get_object('image-container-box')
        self.figure = matplotlib.figure.Figure()
        canvas = FigureCanvas(self.figure)
        self.image_box.pack_start(canvas)

        self.navigation_box = self._builder.get_object('navigation-toolbar-box')
        navig = NavigationToolbar(canvas, self.image_box)
        self.navigation_box.pack_start(navig)
        matplotlib_container.show_all()

        self.error_msg = self._builder.get_object('error-messages-label')
        self.error_msg.set_visible(False)

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
                button.connect('button-press-event', self.show_pfilter, pfilter)
                self.buttons[pfilter] = button

            self.HButtonBox.pack_end(button)
            button.show()

        # The checkbox to enable/disable airmasses in the plots, located in
        # the View submenu and shared by all the StarDetailsGUI instances
        self.airmasses_checkbox = self._builder.get_object('plot-airmasses-checkbox')
        args = 'toggled', self.redraw_light_curve
        self.airmasses_checkbox.connect(*args)

        # Activate the button of the first filter for which there is data
        # (those in the self.buttons dictionary), unless indicated otherwise
        if not init_pfilter:
            init_pfilter = min(self.buttons.keys())

        self.buttons[init_pfilter].set_active(True)
        event = gtk.gdk.Event(gtk.gdk.NOTHING)
        self.buttons[init_pfilter].emit('button-press-event', event)
        self.shown = init_pfilter

        # The button to export the light curve to a text file
        export_button = self._builder.get_object('save-curve-points-button')
        args = 'clicked', self.save_light_curve_points
        export_button.connect(*args)

    def save_light_curve_points(self, widget):
        """ Dump the points of the light curve to a plain text file """

        args = (self.parent._main_window, self._builder, self.config,
                self.id, self.shown, self.db.field_name, self.curve_store)
        dialog = ExportCurveDialog(*args)
        dialog.run()


class SNRThresholdDialog(object):
    """ A GTK.Dialog to select the minimum SNR required for plots """

    def __init__(self, parent_window, builder, config):
        self.builder = builder
        self.builder.add_from_file(glade.SNR_THRESHOLD_DIALOG)
        self.config = config

        self.dialog = self.builder.get_object('snr-treshold-dialog')
        self.dialog.set_resizable(False)
        self.dialog.set_title("Mininum SNR in plots")
        self.dialog.add_button(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL)
        self.dialog.add_button(gtk.STOCK_OK, gtk.RESPONSE_OK)
        self.dialog.set_default_response(gtk.RESPONSE_OK)

        self.dialog.set_transient_for(parent_window)
        self.dialog.set_position(gtk.WIN_POS_CENTER_ON_PARENT)

        self.spinbutton = self.builder.get_object('snr-threshold-spinbutton')
        self.spinbutton.set_value(self.config.get_minimum_snr())

    def run(self):
        """ Run the dialog window in a recursive loop.

        This method shows the dialog window and allows the user to select,
        using a spin button, the minimum SNR to be used in plots. If the value
        of the spin button is modified and the 'Ok' button is clicked, the new
        SNR threshold is returned. On the other hand, clicking 'Ok' without
        having modified the SNR, clicking 'Cancel' or destroying the window
        return None.

        Note that returning None when the SNR is not modified but 'Ok' clicked
        avoids unnecessary updates of the light curves: without this, we would
        replot them with the same SNR threshold.

        """

        try:
            old_snr = self.spinbutton.get_value()
            response = self.dialog.run()
            if response == gtk.RESPONSE_OK:
                snr = self.spinbutton.get_value()
                if snr != old_snr:
                    self.config.set_minimum_snr(snr)
                    return snr

        finally:
            self.dialog.destroy()


class LEMONJuicerGUI(object):

    # The minimum size (width, height) of the GUI, in pixels
    MIN_SIZE = (1024, 768)

    # The label on the tab for those pages with the details of a star
    TABS_LABEL = "Star %d"

    get_abspath = functools.partial(os.path.join, os.path.dirname(__file__))
    COMPASS_ICON = get_abspath('./gui/img/compass.png')

    def _add_custom_stock_icons(self):
        """ Register our own stock icon names """

        factory = gtk.IconFactory()
        pixbuf = gtk.gdk.pixbuf_new_from_file(self.COMPASS_ICON)
        iconset = gtk.IconSet(pixbuf)
        factory.add('Compass', iconset)
        factory.add_default()

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

        self._add_custom_stock_icons()
        self.finding_chart_menuitem = builder.get_object('finding-chart-item')
        self.finding_chart_button = builder.get_object('finding-chart-button')
        self.finding_chart_button.set_stock_id('Compass')
        self.finding_chart_menuitem.set_sensitive(False)
        self.finding_chart_button.set_sensitive(False)

        # The dialog window with the finding chart is kept in memory, hidden
        # (instead of destroyed) when the user closes the window. The dialog
        # can in this manner be shown repeatedly, and it remembers its size,
        # location and any Matplotlib panning and zooming.
        self.finding_chart_dialog = None

        self.amplitudes_search_button = builder.get_object('amplitudes-search-button')
        self.amplitudes_search_menuitem = builder.get_object('amplitudes-search-item')

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

        # The searches for stars by amplitudes (Find → Amplitudes-wavelength
        # correlation) are numbered sequentially and in Roman numerals: i.e.,
        # the first one is labeled with 'I', the second with 'II', etc.
        self.nampl_searches = 0

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
                self.finding_chart_dialog = None
                while self._notebook.get_n_pages():
                    self._notebook.remove_page(-1)

        elif index >= 1:

            # Pages of the notebook that correspond to a star (StarDetailsGUI)
            # have the 'id' attribute, which we need to update the list of open
            # stars before closing the page. If the attribute does not exist we
            # are dealing with a different kind of page and then we can simply
            # remove it.

            try:
                star_id = self._notebook.get_nth_page(index).id
                del self.open_stars[star_id]
            except AttributeError:
                pass
            self._notebook.remove_page(index)

        # Disable the close button / menu item if the notebook becomes empty
        npages_left = self._notebook.get_n_pages()
        self.close_button.set_sensitive(npages_left)
        self.close_menu_item.set_sensitive(npages_left)
        self.finding_chart_menuitem.set_sensitive(npages_left)
        self.finding_chart_button.set_sensitive(npages_left)
        self.amplitudes_search_button.set_sensitive(npages_left)
        self.amplitudes_search_menuitem.set_sensitive(npages_left)

    def handle_quit(self, obj):
        """ Close the application after removing all the pages of the notebook """

        # Update the configuration file with any changes
        self.config.update()

        # I'm not sure why it takes way too longer to destroy the window if
        # there are pages left in the notebook. Are we doing something wrong
        # in our code, or is it inherent to how (Py)GTK works?
        while self._notebook.get_n_pages():
            self._notebook.remove_page(-1)
        gtk.main_quit()

    def handle_show_about(self, obj):
        self._builder.add_from_file(glade.GUI_ABOUT)
        about = self._builder.get_object('about-dialog')
        about.set_transient_for(self._main_window)

        # Replace 'YYYY' in copyright notice with the current year
        copyright = about.get_copyright()
        now = datetime.datetime.now()
        copyright = copyright.replace('YYYY', str(now.year))
        about.set_copyright(copyright)

        about.run()
        about.destroy()

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

            # By default, both LEMON databases (.LEMONdB extension) and XML
            # files (with a previous search for wavelength-amplitude correlated
            # stars) are displayed by the gtk.FileChooserDialog. Two additional
            # filters are added in case the user is interested in a specific
            # type of file.

            db_extension = '.LEMONdB'
            db_pattern = '*' + db_extension
            xml_extension = '.xml'
            xml_pattern = '*' + xml_extension

            # Return the function that can be used to filter files with the
            # gtk.FileFilter.add_custom method. 'filter_info' is a 4-element
            # tuple: the full pathname of the file, the URI of the file, the
            # display name of the file and the MIME type of the file.
            def ends_with(extension):
                def func(filter_info):
                    return filter_info[0].lower().endswith(extension.lower())
                return func
            db_filter = ends_with(db_extension)
            xml_filter = ends_with(xml_extension)

            filt = gtk.FileFilter()
            filt.set_name("All LEMON Files")
            filt.add_custom(gtk.FILE_FILTER_FILENAME, db_filter)
            filt.add_custom(gtk.FILE_FILTER_FILENAME, xml_filter)
            dialog.add_filter(filt)

            filt = gtk.FileFilter()
            filt.set_name("LEMON Database (%s)" % db_pattern)
            filt.add_custom(gtk.FILE_FILTER_FILENAME, db_filter)
            dialog.add_filter(filt)

            filt = gtk.FileFilter()
            filt.set_name('LEMON XML File (%s)' % xml_pattern)
            filt.add_custom(gtk.FILE_FILTER_FILENAME, xml_filter)
            dialog.add_filter(filt)

            response = dialog.run()
            if response == gtk.RESPONSE_OK:
                path = dialog.get_filename()
                dialog.destroy()

                if path.lower().endswith(db_extension.lower()):
                    self.open_db(path)
                else:
                    assert path.lower().endswith(xml_extension.lower())
                    self.open_amplitudes_xml(path)

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

        # Reset the counter of searches for stars by amplitudes: if we are
        # working with a LEMONdB, have made four searches and open another
        # LEMONdB, we do not want the next search to be considered the
        # fifth — it is the *first* search for this database.
        self.nampl_searches = 0

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

            db = mining.LEMONdBMiner(path)
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

                # Although an extremely rare case, LEMONdB.field_name will be
                # None if the name of the observed object cannot be determined
                # (i.e., if there is no prefix common to the object names of
                # the images stored in the database). In those cases, we use
                # the name of the .LEMONdB file for the tab label for the page.

                text = db.field_name or os.path.basename(db.path)
                label = gtk.Label(text)
                self._notebook.insert_page(overview, label, 0)
                self._notebook.set_tab_reorderable(overview, False)
                self._notebook.set_current_page(0)

                # Now that there is at least one page in the notebook, make the
                # 'Close' button and menu item sensitives, as well as other
                # widgets that need an open database with which interact.
                self.close_button.set_sensitive(True)
                self.close_menu_item.set_sensitive(True)
                self.finding_chart_menuitem.set_sensitive(True)
                self.finding_chart_button.set_sensitive(True)
                self.amplitudes_search_button.set_sensitive(True)
                self.amplitudes_search_menuitem.set_sensitive(True)

        except Exception, err:
            path = os.path.basename(path)
            title = "Error while loading LEMON database"
            msg = "File '%s' could not be loaded: %s" % (path, str(err))
            util.show_error_dialog(self._main_window, title, msg)

        finally:
            dialog.destroy()

    def add_star_page(self, star_id, pfilter = None):
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

    def view_star(self, star_id, view = None):
        """ Open the details of the star, or switch to them if already open. """

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
            details = self.add_star_page(star_id, pfilter = pfilter)
            self.open_stars[star_id] = details

    def handle_row_activated(self, view, row, column, id_index = 0):
        """ Handler for when the user double-clicks on a star.

        Keyword arguments:
        id_index - the column of the model to which the view is associated
        which contains the ID of the star.

        """

        # Determine the ID of the star on which the user has clicked
        star_id = view.get_model()[row][id_index]
        self.view_star(star_id, view)

    def append_amplitudes_search(self, result, connect):
        """ Append an AmplitudesSearchPage to the gtk.Notebook.

        Append the gtk.ScrolledWindow (with the result of the search) of an
        AmplitudesSearchPage object to the gtk.Notebook. The 'row-activated'
        signal is connected to the LEMONJuicerGUI.handle_row_activated method
        depending on the truth value of 'connect'.

        """

        if connect:
            view = result.view # gtk.GtkTreeView
            view.connect('row-activated', self.handle_row_activated)

        self.nampl_searches += 1
        label = gtk.Label(result.get_label(self.nampl_searches))
        window = result.get_window()
        self._notebook.append_page(window, label)
        self._notebook.set_tab_reorderable(window, False)
        self._notebook.set_current_page(-1)

    def search_by_amplitudes(self, window):
        """ Identify stars with amplitudes correlated to the wavelength. These
        are listed in a gtk.ScrolledWindow which is appended to the notebook"""

        args = self._main_window, self._builder, self.db.path, self.config
        result = search.amplitudes_search(*args)
        if result is not None:
            self.append_amplitudes_search(result, True)

    def change_snr_threshold(self, widget):
        """ Select a new SNR threshold and update all the plots with it """

        args = self._main_window, self._builder, self.config
        dialog = SNRThresholdDialog(*args)
        threshold = dialog.run()

        # SNRThresholdDialog.run() returns a value other than None when a new
        # minimum SNR value is input and the 'Ok' button is clicked. If that
        # is the case, we iterate over all the StarDetailsGUI's, updating the
        # plots with the new SNR threshold.
        if threshold is not None:
            for details in self.open_stars.itervalues():
                details.redraw_light_curve(None)

    def open_amplitudes_xml(self, path):
        """ Parse an XML file and deserialize an AmplitudesSearchPage object.

        The gtk.ScrolledWindow returned by the AmplitudesSearchPage instance is
        appended to the gtk.Notebook so that the stars found in the serialized
        search are presented to the user exactly as other search results. The
        only difference is that the 'row-activated' signal is ignored if the
        XML file has been opened before the LEMONdB: without access to the
        database we cannot show the details of any star.

        """

        try:
            result = search.AmplitudesSearchPage.xml_load(path)
        except Exception, err:
            title = "Error while loading XML file"
            msg = "File '%s' could not be loaded. Please make sure that you " \
            "are opening an XML file created by LEMON through 'Save As'. " \
            "The error was: \n\n%s" % (path, str(err))
            util.show_error_dialog(self._main_window, title, msg)
            return

        # The search results stored in the XML file must refer to the current
        # LEMONdB, if any. This can be verified as the 'database_id' attribute
        # of the root element must equal the value of the LEMONdB.id property,
        # as the latter is written to the XML file when results are serialized.

        if self.db is not None:
            ids = result.id, self.db.id
            if operator.ne(*ids):
                title = "XML file does not match LEMONdB"
                msg = "The LEMONdB to which the search results in this XML " \
                "file correspond (ID = %s) is different from the one " \
                "currently open (ID = %s). The file, therefore, cannot be " \
                "be loaded, as it refers to a different LEMONdB." % ids
                args = self._main_window, title, msg
                kwargs = dict(buttons = gtk.BUTTONS_CLOSE)
                util.show_error_dialog(*args, **kwargs)
                return

        connect_double_click = self.db is not None
        self.append_amplitudes_search(result, connect_double_click)

        # Make sure that the 'Close' button and menu item are sensitive (they
        # will be disabled if the XML file is opened before the LEMONdB, but
        # we need a way to close the page when we are done with it).
        self.close_button.set_sensitive(True)
        self.close_menu_item.set_sensitive(True)

        # If no LEMONdB is open, double-clicking on a row will have no effect
        if self.db is None:
            title = "Double-click has been disabled"
            msg = "Please note that double-clicking on a star will not open " \
            "a page with its information: as the XML file has been opened " \
            "before the LEMON database, no other information than what is " \
            "presented here is available."
            args = self._main_window, title, msg
            kwargs = dict(msg_type = gtk.MESSAGE_WARNING,
                          buttons = gtk.BUTTONS_CLOSE)
            util.show_message_dialog(*args, **kwargs)

    def view_finding_chart(self, widget):
        """ Display the reference frame in a new dialog """

        # Create the FindingChartDialog the first time this method is called,
        # hiding the dialog when run() exits. In this manner, we can run() as
        # many times as we want to display the dialog again.
        if self.finding_chart_dialog is None:
            self.finding_chart_dialog = chart.FindingChartDialog(self)
        self.finding_chart_dialog.show()

