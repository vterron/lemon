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

import astropy.time
import ConfigParser
import datetime
import functools
import operator
import os.path
import random
import re
import sys
import time
import warnings

import matplotlib
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
import simbad
import util
from version import __version__

# Workaround for GTK+ bug 632538. Until GTK+ 2.24.19 (the PyGTK version does
# not matter, this is all in the C library), some GtkSettings properties are
# registered by other classes. This leads to the "interesting" issue that
# setting GtkSettings:gtk-button-images requires that the GtkButton class
# is referenced first - or that a GtkButton is created.
# [URL] https://bugzilla.gnome.org/show_bug.cgi?id=632538

if gtk.gtk_version < (2, 24, 19):
    with util.destroying(gtk.Button()):
        pass

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

    def __init__(self, parent_window, builder, config, id_, pfilter, db,
                 curve_store):
        """ Instantiation method for the ExportCurveDialog class.

        The 'parent_window' parameter must be the transient parent of the
        dialog, while 'builder' and 'config' are the gtk.GtkBuilder and
        Configuration instances, respectively, of the parent GTK widget. The
        'id_' parameter is the ID of the star, 'pfilter' a Passband instance
        encapsulating the photometric filter of the light curve, and 'db' a
        LEMONdB handle.

        Lastly, 'curve_store' must be a gtk.ListStore with the data that will
        be dumped to a file, and should contain seven columns: (1) the date of
        observation in Unix time, (2) a textual representation of the date, (3)
        the Julian date, (4) the differential magnitude, (5) the SNR, (6 and 7)
        the maximum and minimum errors induced by the noise, respectively. The
        columns must be in this exact order and of type float, except for the
        second one, which is a str.

        """

        self.parent_window = parent_window
        self.builder = builder
        self.builder.add_from_file(glade.EXPORT_CURVE_DIALOG)
        self.config = config
        self.id = id_
        self.pfilter = pfilter
        self.db = db
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
        self.date_julian_checkbox = self.get('date-julian-checkbox')
        self.date_secs_checkbox = self.get('date-secs-checkbox')
        self.mags_checkbox = self.get('mags-checkbox')
        self.snr_checkbox = self.get('snr-checkbox')
        self.merr_pos_checkbox = self.get('merr-pos-checkbox')
        self.merr_neg_checkbox = self.get('merr-neg-checkbox')
        self.inst_mags_checkbox = self.get('inst-mags-checkbox')
        self.inst_snr_checkbox = self.get('inst-snr-checkbox')
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
        self.date_julian_checkbox.connect('toggled', f('dump_date_julian'))
        self.date_secs_checkbox.connect('toggled', f('dump_date_seconds'))
        self.mags_checkbox.connect('toggled', f('dump_magnitude'))
        self.snr_checkbox.connect('toggled', f('dump_snr'))
        self.merr_pos_checkbox.connect('toggled', f('dump_max_merr'))
        self.merr_neg_checkbox.connect('toggled', f('dump_min_merr'))
        self.inst_mags_checkbox.connect('toggled', f('dump_instrumental_magnitude'))
        self.inst_snr_checkbox.connect('toggled', f('dump_instrumental_snr'))
        self.spinbutton.connect('output', f('decimal_places', 'get_value'))

    def update(self):
        """ Update the buttons to the values given in the configuration file """

        config = self.config
        self.date_str_checkbox.set_active(config.dumpint('dump_date_text'))
        self.date_julian_checkbox.set_active(config.dumpint('dump_date_julian'))
        self.date_secs_checkbox.set_active(config.dumpint('dump_date_seconds'))
        self.mags_checkbox.set_active(config.dumpint('dump_magnitude'))
        self.snr_checkbox.set_active(config.dumpint('dump_snr'))
        self.merr_pos_checkbox.set_active(config.dumpint('dump_max_merr'))
        self.merr_neg_checkbox.set_active(config.dumpint('dump_min_merr'))
        self.inst_mags_checkbox.set_active(config.dumpint('dump_instrumental_magnitude'))
        self.inst_snr_checkbox.set_active(config.dumpint('dump_instrumental_snr'))
        self.spinbutton.set_value(config.dumpint('decimal_places'))

    def dump(self, path, separator = '\t'):
        """ Save a light curve to the plain text file 'path'.

        Iterate over the gtk.ListStore and save a textual representation of its
        rows to 'path', truncating the file if it already exists. The rows are
        written to the file chronologically by the date of observation. Not all
        the rows of the gtk.ListStore (i.e., the attributes of the light curve,
        such as the signal-to-noise ratio or the maximum and minimum error
        induced by the noise) are saved to the file, but only those for which
        the corresponding checkbox is active (checked). Magnitudes,
        signal-to-noise ratios and errors are written with the number of
        decimal places set in the spin button.

        """

        def parse_float(value):
            """ Cast value to str; use exactly 'decimals' decimal digits """
            ndecimals = int(self.spinbutton.get_value())
            return '%.*f' % (ndecimals, value)

        # A dictionary mapping each Unix time to a two-element namedtuple with
        # the instrumental magnitude of the star and the SNR of the measurement.
        if self.inst_mags_checkbox.get_active() or \
           self.inst_snr_checkbox.get_active():
            args = (self.id, self.pfilter)
            instrumental = self.db.get_instrumental_magnitudes(*args)

        with open(path, 'wt') as fd:

            # First element of each row is the date of observation (Unix)
            for row in sorted(self.store, key = operator.itemgetter(0)):
                unix_time = row[0]

                values = []
                assert len(row) == 7
                if self.date_secs_checkbox.get_active():
                    values.append(str(unix_time))
                if self.date_str_checkbox.get_active():
                    values.append(row[1])
                if self.date_julian_checkbox.get_active():
                    values.append(str(row[2]))
                if self.mags_checkbox.get_active():
                    values.append(parse_float(row[3]))
                if self.snr_checkbox.get_active():
                    values.append(parse_float(row[4]))
                if self.merr_pos_checkbox.get_active():
                    values.append(parse_float(row[5]))
                if self.merr_neg_checkbox.get_active():
                    values.append(parse_float(row[6]))
                if self.inst_mags_checkbox.get_active():
                    mag = instrumental[unix_time].magnitude
                    values.append(parse_float(mag))
                if self.inst_snr_checkbox.get_active():
                    snr = instrumental[unix_time].snr
                    values.append(parse_float(snr))

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
                field = re.sub(r'\s', '_', self.db.field_name.lower())
                args = field, self.id, self.pfilter
                filename = '%s_star_%d_curve_%s' % args
                chooser.set_current_name(filename)
                response = chooser.run()

                if response == gtk.RESPONSE_OK:
                    self.dump(chooser.get_filename())

        self.dialog.destroy()


class StarDetailsGUI(object):
    """ A tabs of the notebook with all the details of a star """

    SNR_THRESHOLD_ERROR_MSG = \
    "No points have a signal-to-noise ratio above the current threshold\n " \
    "(to modify it, go to View → Plots → SNR threshold)"

    NO_CURVE_IN_ANY_FILTER_ERROR_MSG = \
    "The light curve of the astronomical object has not been generated\n" \
    "in any photometric filter, so there is nothing to be displayed here."

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

    def set_error_msg(self, msg):
        """ Set the text of the error message label """
        self.error_msg.set_label(msg)

    def update_curve(self, curve, show_airmasses, show_julian_dates):

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
            self.set_error_msg(self.SNR_THRESHOLD_ERROR_MSG)
            self.set_canvas(False)

        else:
            self.set_canvas(True)

            # The color to use for each photometric filter is defined in the
            # [colors] section of the .juicer configuration file. However, the
            # user may have added a custom filter (via the [custom_filters]
            # section in .lemonrc) without necessarily specifying a color for
            # it. When that happens, warn the user and pick a random color.

            try:
                color = self.config.color(curve.pfilter.letter)
            except ConfigParser.NoOptionError:
                cycle = matplotlib.rcParams['axes.color_cycle']
                color = random.choice(cycle)
                msg = ("cannot find a color specified for photometric filter "
                       "'{0}', so we will use a random one instead: {1}. This "
                       "probably means that it is a user-defined filter. If "
                       "that is the case, you may want to add its color to "
                       "the [colors] section of the configuration file.".
                       format(curve.pfilter, color))
                warnings.warn(msg)

            kwargs = dict(airmasses = airmasses,
                          julian = show_julian_dates,
                          delta = 3 * 3600,
                          color = color)

            plot.curve_plot(self.figure, curve, **kwargs)
            self.figure.canvas.draw()

    def update_light_curve_points(self, curve):
        """ Update the list of points of the light curve """

        self.curve_store.clear()
        for unix_time, magnitude, noise in curve:
            row = []
            row.append(unix_time)
            row.append(methods.utctime(unix_time, suffix = False))
            row.append(astropy.time.Time(unix_time, format = 'unix').jd)
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
        for star_id, weight, stdev in curve.weights():
            imag = self.db.get_star(star_id)[-1]
            args = (star_id, weight, imag, stdev)
            self.refstars_store.append(args)

        # Sort automatically, but respect the user preferences
        if self.refstars_store.get_sort_column_id() == (None, None):
            self.refstars_store.set_sort_column_id(1, gtk.SORT_DESCENDING)

    def redraw_light_curve(self, widget):
        """ Replot the light curve """

        curve = self.db.get_light_curve(self.id, self.shown)
        args = self.airmasses_visible(), self.julian_dates_visible()
        self.update_curve(curve, *args)

    def update_file_selector_name(self):
        """ Update the name suggested by the 'Save' button FileChooserDialog.

        When the 'Save' button on the matplotlib navigation toolbar is clicked,
        the gtk.FileChooserDialog suggests as filename the value returned by
        NavigationToolbar.canvas.get_default_filename(). If we take a look at
        its source code, we can see that it returns `self.get_window_title() or
        'image'`. That is: if the title of the window is undefined, 'image' is
        used as the default filename.

        Unfortunately, set_window_title() has no effect here because there is
        no window containing the figure. Work around this by monkey-patching
        the method, making it return what we want to be the default filename.
        The string includes the name of the field, the star ID and the filter
        of the light curve: for example, 'ngc_2264_star_1856_curve_h.png'.

        """

        args = self.db.field_name, self.id, self.shown
        filename = '%s_star_%d_curve_%s' % args
        self.canvas.get_window_title = lambda: filename

    def show_pfilter(self, button, event, pfilter):
        """ Display the information of the star in this photometric filter """

        # Ignore clicks on already-pressed buttons (i.e., do not draw the same
        # light curve twice in a row). Note that 'show' is initialized to None,
        # so this method will be always executed the first time it is called.
        if self.shown and self.shown == pfilter:
            return

        self.shown = pfilter
        curve = self.db.get_light_curve(self.id, pfilter)
        args = self.airmasses_visible(), self.julian_dates_visible()
        self.update_curve(curve, *args)
        self.update_light_curve_points(curve)
        self.update_reference_stars(curve)
        # Keep track of the filter in the TreeView object, so that
        # LEMONJuicerGUI.handle_row_activated can know in which filter we
        # were working when we clicked on one of the reference stars
        self.refstars_view.pfilter = pfilter
        self.update_file_selector_name()

    def handle_toggle_view_sexagesimal(self, *args):
        button = self._builder.get_object('radio-view-sexagesimal')
        for index in self.sex_indexes:
            self.starinfo_store[index][-1] = button.get_active()

    def handle_toggle_view_decimal(self, *args):
        button = self._builder.get_object('radio-view-decimal')
        for index in self.dec_indexes:
            self.starinfo_store[index][-1] = button.get_active()

    def airmasses_visible(self):
        """ Return the state (active or not) of the airmasses checkbox """
        return self.airmasses_checkbox.get_active()

    def julian_dates_visible(self):
        """ Return the state (active or not) of the Julian dates checkbox """
        return self.julian_dates_checkbox.get_active()

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
        self.canvas = FigureCanvas(self.figure)
        self.image_box.pack_start(self.canvas)

        self.navigation_box = self._builder.get_object('navigation-toolbar-box')
        navig = NavigationToolbar(self.canvas, self.image_box.get_window())

        self.navigation_box.pack_start(navig)
        matplotlib_container.show_all()
        self.error_msg = self._builder.get_object('error-messages-label')
        self.error_msg.set_visible(False)

        # Button to open the Finding Chart and mark this star on it. Render a
        # stock button, STOCK_JUMP_TO, but with a different label. Do this by
        # getting the label child of the stock button and setting the text
        # directly, as suggested by Christian Reis in the PyGTK FAQ.
        # http://faq.pygtk.org/index.py?req=show&file=faq09.005.htp

        arg = 'view-star-in-finding-chart-button'
        self.view_in_chart_button = self._builder.get_object(arg)
        alignment = self.view_in_chart_button.get_children()[0]
        hbox = alignment.get_children()[0]
        image, label = hbox.get_children()
        label.set_text('View in Finding Chart')

        args = 'clicked', self.handle_view_star_in_chart
        self.view_in_chart_button.connect(*args)

        # GTKTreeView used to display the list of points of the curve; dates
        # are plotted three times: as (a) Unix times, (b) 24-character strings
        # with the format 'Www Mmm dd hh:mm:ss yyyy' and (b) Julian dates. The
        # date strings are not a floating-point numbers and thus cannot be
        # sorted by themselves, so we use the Unix time (the first column).

        attrs = ('Date', 'Date', 'JD', 'Δ Mag', 'SNR', 'merr (+)', 'merr (-)')
        column_types = float, str, float, float, float, float, float
        self.curve_store = gtk.ListStore(*column_types)
        self.curve_view = self._builder.get_object('curve-points-view')
        for index, title in enumerate(attrs):
            render = gtk.CellRendererText()
            column = gtk.TreeViewColumn(title, render, text = index)
            column.props.resizable = False
            # Date strings are sorted by the Unix time
            column.set_sort_column_id(0 if index == 1 else index)

            # The column with dates in Unix time is not shown
            if not index:
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

            # Instead of the six decimal places that gtk.TreeViewColumn uses by
            # default, limit the instrumental magnitudes to three. This is the
            # maximum precision that we can get using IRAF's qphot, and it does
            # not look like more decimal places will be available any time
            # soon, so there is no point in using additional zeros.
            # http://iraf.net/forum/viewtopic.php?showtopic=1467248

            if index == 2:

                def render_imag(column, cell, model, iter):
                    """ Render the cell with three decimal places """
                    orig_str = cell.get_property('text')
                    new_str = '%.3f' % float(orig_str)
                    cell.set_property('text', new_str)

                # http://faq.pygtk.org/index.py?req=show&file=faq13.024.htp
                column.set_cell_data_func(render, render_imag)
                assert title == 'Mag'

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
        x, y, self.ra, self.dec, _, _, _, imag = self.db.get_star(star_id)
        store.append(('Right ascension', '%s' % methods.ra_str(self.ra), True))
        store.append(('Right ascension', '%.4f deg' % self.ra, True))
        store.append(('Declination', '%s' % methods.dec_str(self.dec), True))
        store.append(('Declination', '%.4f deg' % self.dec, True))
        store.append(('Magnitude', '%.3f' % imag, True))
        store.append(('x-coordinate', '%.2f' % x, True))
        store.append(('y-coordinate', '%.2f' % y, True))

        # Two rows (sexagesimal and decimal) are used for the coordinates, but
        # only one will be shown at a given time. To hide some of the rows of a
        # TreeView, we need to use a TreeModelFilter, which acts as a wrapper
        # for the TreeModel, allowing us to choose which rows are displayed
        # based on the value of a gobject.TYPE_BOOLEAN column, or based on the
        # output of a certain function.
        # [http://faq.pygtk.org/index.py?file=faq13.048.htp&req=show]

        self.sex_indexes = [0, 2]
        self.dec_indexes = [1, 3]

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

        # The checkbox to alternate between strings representing the date and
        # time (for example, 'Jan 02 2012' or '08:15:31', depending on the date
        # range) and Julian dates (JD) such as 2456877.9660300924.
        object_name = 'plot-julian-dates-checkbox'
        self.julian_dates_checkbox = self._builder.get_object(object_name)
        args = 'toggled', self.redraw_light_curve
        self.julian_dates_checkbox.connect(*args)

        # The button to export the light curve to a text file
        self.export_button = self._builder.get_object('save-curve-points-button')
        args = 'clicked', self.save_light_curve_points
        self.export_button.connect(*args)

        # Activate the button of the first filter for which there is data
        # (those in the self.buttons dictionary), unless indicated otherwise

        if self.buttons.keys():

            if not init_pfilter:
                init_pfilter = min(self.buttons.keys())

            self.buttons[init_pfilter].set_active(True)
            event = gtk.gdk.Event(gtk.gdk.NOTHING)
            self.buttons[init_pfilter].emit('button-press-event', event)
            self.shown = init_pfilter

        else:
            # If the light curve of the astronomical object has not been
            # computed in any photometric filter, show the error message
            # explaining that there is nothing that we can display here.
            self.set_error_msg(self.NO_CURVE_IN_ANY_FILTER_ERROR_MSG)
            self.set_canvas(False)
            # There are no light curves whose data we can save to disk
            self.export_button.set_sensitive(False)

        # Button to look up the star in the SIMBAD astronomical database.
        # Render a stock button, STOCK_INFO, but use a different label.

        arg = 'look-up-star-in-simbad-button'
        self.look_up_in_simbad_button = self._builder.get_object(arg)
        alignment = self.look_up_in_simbad_button.get_children()[0]
        hbox = alignment.get_children()[0]
        image, label = hbox.get_children()
        label.set_text("Look up in SIMBAD")

        args = 'clicked', self.handle_look_up_in_simbad
        self.look_up_in_simbad_button.connect(*args)

        # The column with the Julian Dates is only visible when the View ->
        # Plot -> 'Julian dates' checkbox is toggled. This code is here, not
        # where gtk.TreeView self.curve_view is created, because at that point
        # julian_dates_visible() would raise AttributeError ("'StarDetailsGUI'
        # object has no attribute 'julian_dates_checkbox'"), since the object
        # does not yet exist.

        julian_column = self.curve_view.get_column(2)
        assert julian_column.get_title() == 'JD'
        visible = self.julian_dates_visible()
        julian_column.set_visible(visible)

    def save_light_curve_points(self, widget):
        """ Dump the points of the light curve to a plain text file """

        args = (self.parent._main_window, self._builder, self.config,
                self.id, self.shown, self.db, self.curve_store)
        dialog = ExportCurveDialog(*args)
        dialog.run()

    def handle_view_star_in_chart(self, widget):
        """ Show the Finding Chart window and mark the star on it.

        This is the callback function for the 'View in Finding Chart' button
        (button_press_event). It shows the finding chart gtk.Dialog, if it's
        not already visible, and overlays a green marker on this star.

        """

        main_window = self.parent
        main_window.handle_finding_chart(self, set_visibility = True)
        main_window.finding_chart_dialog.mark_star(self.id)

    def handle_look_up_in_simbad(self, widget):
        """ Look up the star in the SIMBAD database.

        This is the callback function for the 'Look up in SIMBAD' button
        (button_press_event). It submits a coordinate-query to the SIMBAD
        database, opening the page in the same window (if possible) of the
        default browser.

        """

        args = (self.ra, self.dec)
        simbad.coordinate_query(*args)


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
    LEMON_ICON = get_abspath('./gui/img/lemon.png')
    COMPASS_ICON = get_abspath('./gui/img/compass.png')

    def _add_custom_stock_icons(self):
        """ Register our own stock icon names """

        factory = gtk.IconFactory()
        pixbuf = gtk.gdk.pixbuf_new_from_file(self.COMPASS_ICON)
        iconset = gtk.IconSet(pixbuf)
        factory.add('Compass', iconset)
        factory.add_default()

    def __init__(self, db_path = None):
        """ Initialize a LEMONJuicerGUI object.

        Keyword arguments:
        db_path - LEMONdB to open when Juicer starts.

        """

        super(LEMONJuicerGUI, self).__init__()

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

        # Append version number to the title
        title = self._main_window.get_title()
        title += " " + __version__
        self._main_window.set_title(title)

        # Manually connect to the 'toggled' signal, storing the handler ID. We
        # need it in order to suppress the signal that is emitted when we set
        # the state of the gtk.ToggleToolButton with set_active(), but there
        # is no way to get the signal handler IDs if they are connected by
        # Gtk.Builder. [http://stackoverflow.com/a/11803778/184363]
        args = 'toggled', self.handle_finding_chart
        self.chart_handler_id = self.finding_chart_button.connect(*args)

        # The dialog window with the finding chart is kept in memory, hidden
        # (instead of destroyed) when the user closes the window. The dialog
        # can in this manner be shown repeatedly, and it remembers its size,
        # location and any Matplotlib panning and zooming.
        self.finding_chart_dialog = None

        # Create an accelerator group and add it to the toplevel window. These
        # accelerators will be always available (except when a modal dialog is
        # active, of course), since they are 'inherited' by the child windows.
        self.global_accelators = gtk.AccelGroup()
        self._main_window.add_accel_group(self.global_accelators)

        # Connect <Ctrl>F to LEMONJuicerGUI.chart_callback()
        key, modifier = gtk.accelerator_parse('<Control>F')
        args = key, modifier, gtk.ACCEL_VISIBLE, self.chart_callback
        self.global_accelators.connect_group(*args)

        # Connect <Ctrl>V to LEMONJuicerGUI.handle_view_in_chart_accelerator()
        key, modifier = gtk.accelerator_parse('<Control>V')
        args = (key, modifier, gtk.ACCEL_VISIBLE,
                self.handle_view_in_chart_accelerator)
        self.global_accelators.connect_group(*args)

        # Use <Ctrl>S to look up the astronomical object in SIMBAD
        key, modifier = gtk.accelerator_parse('<Control>S')
        args = (key, modifier, gtk.ACCEL_VISIBLE,
                self.handle_look_up_in_simbad_accelerator)
        self.global_accelators.connect_group(*args)

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

        checkbox = builder.get_object('plot-julian-dates-checkbox')
        checkbox.set_active(get_view_booloption(config.PLOT_JULIAN))

        if db_path:
            self.open_db(db_path)

    def _activate_StarDetailsGUI_button(self, name):
        """ Activate a button of the current StarDetailsGUI.

        If the current page in the gtk.Notebook corresponds to a StarDetailsGUI
        instance (that is, all pages except for the first one, with index zero,
        which contains the list of stars), call the activate() method of its
        'name' button.

        """

        index = self._notebook.get_current_page()
        if index:
            star_id = self._notebook.get_nth_page(index).id
            star_details = self.open_stars[star_id]
            getattr(star_details, name).activate()

    def handle_view_in_chart_accelerator(self, *args):
        """ Show the star in the Finding Chart.

        Activate the 'Show in Finding Chart' button of the star in the current
        StarDetailsGUI.

        """

        name = 'view_in_chart_button'
        self._activate_StarDetailsGUI_button(name)

    def handle_look_up_in_simbad_accelerator(self, *args):
        """ Look up the star in the SIMBAD database, using the default browser.

        Activate the 'Look up in SIMBAD' button of the star in the current
        StarDetailsGUI.

        """

        name = 'look_up_in_simbad_button'
        self._activate_StarDetailsGUI_button(name)

    def save_widget_state(self, widget, section, option):
        """ Update the specified section and option of the configuration file
        with the state of the gtk.Widget: 1 if its active and 0 otherwise """

        active = widget.get_active()
        value = '1' if active else '0'
        self.config.set(section, option, value)

    # Airmasses and Julian dates (JDs) are not plotted here (that is done in
    # StarDetailsGUI), but we need to update the configuration file with the
    # new values of the options every time their checkboxes are toggled.

    def save_plot_airmasses_checkbox(self, widget):
        args = widget, config.VIEW_SECTION, config.PLOT_AIRMASSES
        self.save_widget_state(*args)

    def save_plot_julian_dates_checkbox(self, widget):
        args = widget, config.VIEW_SECTION, config.PLOT_JULIAN
        self.save_widget_state(*args)

    def run(self):
        gtk.main()

    def close_database(self):
        """ Forget about the LEMONdB to which we are currently connected """

        self.db = None

        # Close all the tabs of the notebook.
        while self._notebook.get_n_pages():
            self._notebook.remove_page(-1)

        # Destroy and forget about the current FindingChartDialog; if we're
        # going to open another LEMONdB we will need a new finding chart.
        if self.finding_chart_dialog is not None:
            self.finding_chart_dialog.destroy()
            self.finding_chart_dialog = None

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
                self.close_database()

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

        # Replace 'x.x.x' with the current version
        version_str = about.get_version()
        assert version_str == 'x.x.x'
        about.set_version(__version__)

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

    def handle_open(self, window):
        kwargs = dict(title = None,
                      action = gtk.FILE_CHOOSER_ACTION_OPEN,
                      buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                                 gtk.STOCK_OPEN, gtk.RESPONSE_OK))

        with util.destroying(gtk.FileChooserDialog(**kwargs)) as dialog:

            # Only LEMON databases (.LEMONdB extension) are displayed by the
            # gtk.FileChooserDialog.

            db_extension = '.LEMONdB'
            db_pattern = '*' + db_extension

            # Return the function that can be used to filter files with the
            # gtk.FileFilter.add_custom method. 'filter_info' is a 4-element
            # tuple: the full pathname of the file, the URI of the file, the
            # display name of the file and the MIME type of the file.
            def ends_with(extension):
                def func(filter_info):
                    return filter_info[0].lower().endswith(extension.lower())
                return func
            db_filter = ends_with(db_extension)

            filt = gtk.FileFilter()
            filt.set_name("LEMON Database (%s)" % db_pattern)
            filt.add_custom(gtk.FILE_FILTER_FILENAME, db_filter)
            dialog.add_filter(filt)

            dialog.set_icon_from_file(self.LEMON_ICON)

            response = dialog.run()
            if response == gtk.RESPONSE_OK:
                path = dialog.get_filename()
                dialog.destroy()

                assert path.lower().endswith(db_extension.lower())
                self.open_db(path)

    def open_db(self, path):

        if not os.path.exists(path):
            msg = "database '%s' does not exist" % path
            raise IOError(msg)

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
                self.close_database()
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

                x, y, ra, dec, _, _, _, imag = db.get_star(star_id)
                ra_str  = methods.ra_str(ra)
                dec_str = methods.dec_str(dec)
                row = [star_id, ra_str, ra, dec_str, dec, imag]

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
                    with util.gtk_sync():
                        progressbar.set_fraction(fraction)

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

    def change_JDs_visibility(self, widget):
        """ Set the visibility of the column with Julian dates in all the
        StarDetailsGUI objects, depending on the state (active or not) of
        the View -> Plots -> 'Julian dates' checkbox. """

        for details in self.open_stars.itervalues():
            julian_column = details.curve_view.get_column(2)
            assert julian_column.get_title() == 'JD'
            visible = widget.get_active()
            julian_column.set_visible(visible)

    def handle_finding_chart(self, widget, set_visibility = None):
        """ Display the reference frame in a new dialog.

        Create a gtk.Dialog() the first time this method is called and show()
        it. Subsequent calls to this method will alternate between hide()'ing
        and show()'ing the dialog.

        The status of the 'Finding Chart' button (a gtk.ToggleToolButton) in
        the toolbar is toggled or untoggled to reflect whether the chart is
        being shown. Note that this happens independently on where the user
        clicked: showing or hiding the chart through the appropriate item in
        the View menu or with the <Control>F accelerator also causes the
        ToggleToolButton to be updated accordingly.

        Keyword arguments:
        set_visibility - if set to a value other than None, the gtk.Dialog() is
                         show()'n or hide()'n depending on its truth value, nor
                         on its visibility. In other words: if this argument is
                         True, the gtk.Dialog is shown if currently hidden, but
                         nothing happens if it is already visible. If False,
                         the gtk.Dialog is hidden if visible, and nothing
                         happens if if is already hidden.

        """

        # Create the FindingChartDialog the first time this method is called,
        # hiding the dialog when run() exits. In this manner, we can run() as
        # many times as we want to display the dialog again.
        if self.finding_chart_dialog is None:
            self.finding_chart_dialog = chart.FindingChartDialog(self)

        def show():
            self.finding_chart_dialog.show()
            self.set_finding_chart_button_active(True)

        def hide():
            self.finding_chart_dialog.hide()
            self.set_finding_chart_button_active(False)

        is_visible = self.finding_chart_dialog.is_visible()
        if set_visibility is True:
            if not is_visible:
                show()
        elif set_visibility is False:
            if is_visible:
                hide()
        elif not is_visible:
            show()
        else:
            hide()

    def set_finding_chart_button_active(self, active):
        """ Update state of 'Finding Chart' button without emitting a signal.

        The 'Finding Chart' button is a gtk.ToggleToolButton, but its status
        can change from other events than just clicking it: (a) the item with
        the same name in the View menu and (b) the <Control>F accelerator. When
        the chart is shown or hidden using these other two options, we want to
        update the status of the button in the toolbar accordingly. Use this
        method to toggle or untoggle the 'Finding Chart' without emitting the
        'toggled' signal, thanks to gobject.GObject.handler_block().

        """

        self.finding_chart_button.handler_block(self.chart_handler_id)
        self.finding_chart_button.set_active(active)
        self.finding_chart_button.handler_unblock(self.chart_handler_id)

    def chart_callback(self, accel_group, acceleratable, keyval, modifier):
        """ Callback function for the <Ctrl>F accelerator.

        This method calls LEMONJuicerGUI.handle_finding_chart(), but only if a
        LEMONdB is open (otherwise, there is no chart to show and, in the same
        way that the 'Finding Chart' button is disabled, <Ctrl>F should not
        allow us to display anything either). It returns True because callbacks
        assigned in connect_group() must return this value if the accelerator
        was handled by the callback; otherwise GTK+ will segfault sooner or
        later.

        """

        if self.db is not None:
            self.handle_finding_chart(acceleratable)
        return True

