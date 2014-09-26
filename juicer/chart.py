#! /usr/bin/env python

# Copyright (c) 2013 Victor Terron. All rights reserved.
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

import astropy.wcs
import atexit
import aplpy
import gtk
import logging
import methods
import numpy
import os
import pyfits

import matplotlib.figure
from matplotlib.backends.backend_gtkagg \
     import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg \
     import NavigationToolbar2GTKAgg as NavigationToolbar

# LEMON modules
import glade
import util

class PreferencesDialog(object):
    """ gtk.Dialog to configure the finding chart normalization parameters.

    This class encapsulates a GTK dialog window that allows the user the adjust
    the value of the 'vmin' and 'vmax' parameters of APLpy's logarithmic
    normalization algorithm (stretch = 'log'). This scale transforms the values
    from the image to the range [0, 1], which then corresponds to the endpoints
    of the grayscale:

    http://aplpy.readthedocs.org/en/latest/normalize.html#log-scale

    """

    def __init__(self, parent):
        """ Initialize a new PreferencesDialog object.

        The gtk.Dialog has two spin buttons, so that the user can select the
        value for both Vmin and Vmax. The lower and upper range values of these
        parameters is determined by the minimum and maximum values of the image
        pixels. There is also the restriction that Vmin <= Vmax. Additionally,
        two non-sensitive (disabled) gtk.Entry fields let the user know what
        the absolute minimum and maximum allowed values are.

        The dialog has three buttons: (a) 'Close', (b) 'Apply' and (c) 'Save'.
        Clicking 'Apply' immediately updates the finding chart image
        (parent.aplpy_plot), using the logarithmic normalization algorithm
        (stretch = 'log') and the Vmin and Vmax values of the spin buttons.
        'Save', on the other hand, stores Vmin and Vmax in the LEMONdB, so
        that they can be reused in the future.

        The initial value of the two spin buttons is read from the LEMONdB (via
        the 'vmin' and 'vmax' properties). If these values have not been stored
        in the database (and, therefore, None is returned), the normalization
        parameters (stretch algorithm, Vmin and Vmax) that we use are those
        stored in the APLpyNormalize object (parent.aplpy_plot.image.norm)
        created by FITSFigure.show_grayscale(), method that must have been
        previously called in FindingChartDialog.__init__(). By default,
        FITSFigure.show_grayscale() normalizes the image linearly: once
        we click 'Apply' we switch to a logarithmic stretch function.

        """

        assert isinstance(parent, FindingChartDialog)

        self.parent = parent
        self.db = self.parent.db

        builder = self.parent.builder
        builder.add_from_file(glade.CHART_PREFERENCES_DIALOG)

        self.dialog = builder.get_object('chart-preferences-dialog')
        self.dialog.set_transient_for(self.parent.dialog)
        self.dialog.set_title("Finding Chart: Preferences")
        self.dialog.set_resizable(False)

        # Note: gtk.RESPONSE_SAVE doesn't exist; we use gtk.RESPONSE_OK
        self.close_button = self.dialog.add_button(gtk.STOCK_CLOSE, gtk.RESPONSE_CLOSE)
        self.apply_button = self.dialog.add_button(gtk.STOCK_APPLY, gtk.RESPONSE_APPLY)
        self.save_button  = self.dialog.add_button(gtk.STOCK_SAVE, gtk.RESPONSE_OK)
        self.dialog.set_default_response(gtk.RESPONSE_CLOSE)
        self.dialog.set_focus(self.close_button)

        text = "Update chart with these parameters"
        self.apply_button.set_tooltip_text(text)

        text = "Store these parameters in the LEMONdB"
        self.save_button.set_tooltip_text(text)

        # Spin buttons to select the value of Vmin / Vmax
        self.vmin_button = builder.get_object('vmin-spinbutton')
        self.vmax_button = builder.get_object('vmax-spinbutton')

        # If the values of both Vmin and Vmax are stored in the LEMONdB, assume
        # a logarithmic scale. Otherwise, use the normalization algorithm and
        # Vmin / Vmax values defined by the APLpyNormalize object of the
        # finding chart (parent.aplpy_plot.image.norm). Note that, by default,
        # FITSFigure.show_grayscale() uses a linear stretch.

        try:
            self.stretch = 'log'
            vmin = self.db.vmin
            vmax = self.db.vmax

            msg1 = "Normalization parameters (vmin and vmax) read from LEMONdB"
            msg2 = "Assuming logarithmic normalization (stretch = 'log')"
            for message in msg1, msg2:
                logging.debug(message)

        except AttributeError:
            normalize = self.parent.aplpy_plot.image.norm
            self.stretch = normalize.stretch
            vmin    = normalize.vmin
            vmax    = normalize.vmax

            msg1 = "Normalization parameters not stored in the LEMONdB"
            msg2 = "Algorithm and values read from APLpyNormalize object"
            for message in msg1, msg2:
                logging.debug(message)

        # The absolute minimum and maximum values that the two spin buttons can
        # take are those of the minimum and maximum pixels of the finding chart
        # image. These are read from the 'data_min' and 'data_max' attributes
        # of the parent FindingChartDialog object.

        data_min = numpy.ceil (min(self.parent.data_min, vmin))
        data_max = numpy.floor(max(self.parent.data_max, vmax))
        assert hasattr(self, 'stretch')

        kwargs = dict(lower = data_min, upper = data_max, step_incr = 1)
        vmin_adjust = gtk.Adjustment(value = vmin, **kwargs)
        vmax_adjust = gtk.Adjustment(value = vmax, **kwargs)
        self.vmin_button.set_adjustment(vmin_adjust)
        self.vmax_button.set_adjustment(vmax_adjust)

        def ndigits(n):
            """ Return the number of digits of an integer """
            return len(str(abs(n)))

        # The desired width of the button, in characters
        self.vmin_button.set_width_chars(ndigits(data_min))
        self.vmax_button.set_width_chars(ndigits(data_max))

        # Show the absolute minimum and maximum allowed values
        data_min_entry = builder.get_object('data-min-entry')
        data_min_entry.set_width_chars(ndigits(data_min))
        data_min_entry.set_text(str(data_min))
        data_min_entry.set_sensitive(False)

        data_max_entry = builder.get_object('data-max-entry')
        data_max_entry.set_width_chars(ndigits(data_max))
        data_max_entry.set_text(str(data_max))
        data_max_entry.set_sensitive(False)

        # Both spin buttons must be in the range [data_min, data_max], but
        # there is a second restriction: Vmin must be at all times <= Vmax.
        # Use the 'value-changed' signal, emitted when any of the settings
        # (i.e. value, digits) that change the display of the spinbutton are
        # changed, to enforce this. Every time that Vmin is changed we make
        # sure that it is <= Vmax; otherwise we set it to Vmax. The same is
        # done with Vmax, ensuring that it is always >= Vmin.

        def vmin_changed_callback(*args):
            upper = self.vmax_button.get_value()
            if self.vmin_button.get_value() > upper:
                self.vmin_button.set_value(upper)

        def vmax_changed_callback(*args):
            lower = self.vmin_button.get_value()
            if self.vmax_button.get_value() < lower:
                self.vmax_button.set_value(lower)

        self.vmin_button.connect('value-changed', vmin_changed_callback)
        self.vmax_button.connect('value-changed', vmax_changed_callback)
        self.dialog.connect('response', self.handle_response)

    def hide(self):
        """ Hide the gtk.Dialog """
        self.dialog.hide()

    def show(self):
        """ Display the gtk.Dialog """
        self.dialog.set_default_response(gtk.RESPONSE_CLOSE)
        self.dialog.show()

    def destroy(self):
        """ Destroy the gtk.Dialog """
        self.dialog.destroy()

    def normalize_plot(self):
        """ Update the finding chart with the current normalization parameters.

        Create an APLpyNormalize object (an APLpy class that allows different
        stretching functions for astronomical images) and set it as the
        normalization instance of the FITSFigure (self.parent.aplpy_plot).
        This method immediately updates the finding chart with the current
        normalization algorithm and the Vmin and Vmax values selected in the
        spin buttons.

        """

        kwargs = dict(stretch = self.stretch,
                      vmin = self.vmin_button.get_value(),
                      vmax = self.vmax_button.get_value())
        norm = aplpy.normalize.APLpyNormalize(**kwargs)
        self.parent.aplpy_plot.image.set_norm(norm)
        self.parent.aplpy_plot.refresh()

    def handle_response(self, widget, response):
        """ Callback function for the 'response' signal.

        This is the callback function that should be registered as the handler
        for the 'response' signal emitted by the gtk.Dialog. If the user clicks
        'Close' (gtk.RESPONSE_CLOSE) or closes the dialog, it is hidden (never
        destroyed, so that the changes to the two spin buttons are not lost).
        'Apply' (gtk.RESPONSE_APPLY) calls PreferencesDialog.normalize_plot(),
        thus immediately updating the finding chart with the new normalization
        parameters. Lastly, 'Save' (we are using gtk.RESPONSE_OK because there
        is no such a thing as gtk.RESPONSE_SAVE) stores the values of Vmin and
        Vmax in the LEMONdB.

        """


        # Don't destroy the gtk.Dialog if we click on the close button
        if response in (gtk.RESPONSE_CLOSE, gtk.RESPONSE_DELETE_EVENT):
            self.hide()

        elif response == gtk.RESPONSE_APPLY:
            with util.disable_while(self.apply_button):
                self.normalize_plot()
        else:
            # For lack of a better response ID
            assert response == gtk.RESPONSE_OK
            with util.disable_while(self.save_button):
                self.db.vmin = self.vmin_button.get_value()
                self.db.vmax = self.vmax_button.get_value()
                self.db.commit()


class FindingChartDialog(object):
    """ GTK.Dialog that displays the reference frame """

    # Width of the window, in pixels
    WIDTH = 650

    # For markers plotted with FITSFigure.show_markers()
    MARK_RADIUS = 60

    # The name of the scatter layer (the 'layer' keyword argument) when
    # FITSFigure.show_markers() is called to overlay markers on a plot.
    MARKERS_LAYER = 'markers'

    def handle_response(self, widget, response):
        """ Handler for the dialog 'response' event """

        if response == gtk.RESPONSE_APPLY:
            self.goto_star()
        # For lack of a better response ID
        elif response == gtk.RESPONSE_OK:
            self.preferences_dialog.dialog.run()
        else:
            # Untoggle gtk.ToggleToolButton in toolbar (main window)
            assert response in (gtk.RESPONSE_CLOSE, gtk.RESPONSE_DELETE_EVENT)
            self.toggle_toolbar_button(False)
            self.hide()

    def __init__(self, parent):

        self.db = parent.db
        parent_window = parent._main_window
        self.view_star = parent.view_star # LEMONJuicerGUI.view_star()
        self.toggle_toolbar_button = parent.set_finding_chart_button_active

        self.builder = gtk.Builder()
        self.builder.add_from_file(glade.FINDING_CHART_DIALOG)
        self.dialog = self.builder.get_object('finding-chart-dialog')
        self.dialog.set_resizable(True)
        self.dialog.set_title("Finding Chart: %s" % self.db.field_name)

        # gtk.RESPONSE_PREFERENCES doesn't exist: use gtk.RESPONSE_OK
        self.dialog.add_button(gtk.STOCK_PREFERENCES, gtk.RESPONSE_OK)
        self.dialog.add_button(gtk.STOCK_CLOSE, gtk.RESPONSE_CLOSE)
        self.dialog.set_default_response(gtk.RESPONSE_CLOSE)

        # Connect to the accelerators of the parent window
        self.dialog.add_accel_group(parent.global_accelators)

        # This private variable stores whether the gtk.Dialog is currently
        # visible or not. It is update to True and False when the show() and
        # hide() methods are called, respectively.
        self._currently_shown = False

        # The matplotlib figure...
        matplotlib_container = self.builder.get_object('matplotlib-container')
        self.image_box = self.builder.get_object('image-container-box')
        self.figure = matplotlib.figure.Figure()
        canvas = FigureCanvas(self.figure)
        self.image_box.pack_start(canvas)

        # ... and the navigation toolbar
        self.navigation_box = self.builder.get_object('navigation-toolbar-box')
        self.navig = NavigationToolbar(canvas, self.image_box.get_window())
        self.navigation_box.pack_start(self.navig)
        matplotlib_container.show_all()

        # Use the field name as the suggested filename of the FileChooserDialog
        # when the user clicks on the 'Save' button on the navigation toolbar.
        # See the docstring of app.StarDetailsGUI.update_file_selector_name()
        # for an explanation on why we are monkey-patching this method.
        canvas.get_window_title = lambda: self.db.field_name

        self.dialog.set_transient_for(parent_window)
        self.dialog.set_position(gtk.WIN_POS_CENTER_ON_PARENT)

        ax1 = self.figure.add_subplot(111)
        ax1.get_xaxis().set_visible(False)
        ax1.get_yaxis().set_visible(False)

        # Temporarily save to disk the FITS file used as a reference frame
        path = self.db.mosaic
        atexit.register(methods.clean_tmp_files, path)
        self.wcs = astropy.wcs.WCS(path)
        with pyfits.open(path) as hdu:
            data = hdu[0].data
            # Ignore any NaN pixels
            self.data_min = numpy.nanmin(data)
            self.data_max = numpy.nanmax(data)

        self.aplpy_plot = aplpy.FITSFigure(path, figure = self.figure)
        self.figure.canvas.mpl_connect('button_press_event', self.mark_closest_star)
        self.aplpy_plot.show_grayscale()

        self.aplpy_plot.add_grid()
        self.aplpy_plot.grid.set_alpha(0.2)

        # Create a PreferencesDialog object, whose __init__() method reads Vmin
        # and Vmax from the LEMONdB or, in case these values are not stored in
        # the database, uses the values computed by FITSFigure.show_grayscale()
        # and stored in the APLpyNormalize object. After the PreferencesDialog
        # is created we can call its normalize_plot() method for the first time
        # in order to apply the normalization parameters to the finding chart.
        #
        # Note: we must create the PreferencesDialog *after* show_grayscale()
        # has been called, because until then the underlying APLpyNormalize
        # object (i.e., aplpy.FITSFigure.image.norm) does not exist, and the
        # __init__() method of the PreferencesDialog class needs to access to
        # it if the values of Vmin and Vmax cannot be read from the LEMONdB.

        assert hasattr(self.aplpy_plot.image, 'norm')
        self.preferences_dialog = PreferencesDialog(self)
        self.preferences_dialog.normalize_plot()

        # The dialog has always the same width; the height is adjusted
        # proportionally depending on the dimensions of the FITS image.
        size = data.shape[::-1]
        size_ratio = size[1] / size[0]
        new_size = self.WIDTH, int(self.WIDTH * size_ratio)
        self.dialog.resize(*new_size)

        # We cannot run() the dialog because it blocks in a recursive main
        # loop, while we want it to be non-modal. Therefore, we need to show()
        # it. But this means that we cannot get the response ID directly from
        # run(), so we have to connect to the dialog 'response' event.
        self.dialog.connect('response', self.handle_response)

        # Don't destroy the gtk.Dialog if we click on the window's close button
        self.dialog.connect('delete-event', self.on_delete_event)

        # Button to, when a star is selected, view its details. We want to
        # render a stock button, STOCK_GO_FORWARD, but with a different label.
        # In order to achieve this, we register our own stock icon, reusing
        # the image from the existing stock item, as seen in the PyGTK FAQ:
        # http://faq.pygtk.org/index.py?req=show&file=faq09.005.htp

        STOCK_GOTO_STAR = 'goto-star-custom'
        gtk.stock_add([(STOCK_GOTO_STAR, '_Go to Star', 0, 0, None)])
        factory = gtk.IconFactory()
        factory.add_default()
        style = self.dialog.get_style()
        icon_set = style.lookup_icon_set(gtk.STOCK_GO_FORWARD)
        factory.add(STOCK_GOTO_STAR, icon_set)

        args = STOCK_GOTO_STAR, gtk.RESPONSE_APPLY
        self.goto_button = self.dialog.add_button(*args)
        self.goto_button.set_sensitive(False)

        # <Ctrl>G also activates the 'Go to Star' button
        accelerators = gtk.AccelGroup()
        self.dialog.add_accel_group(accelerators)
        key, modifier = gtk.accelerator_parse('<Control>G')
        args = 'activate', accelerators, key, modifier, gtk.ACCEL_VISIBLE
        self.goto_button.add_accelerator(*args)


    def mark_closest_star(self, event):
        """ Callback function for 'button_press_event'.

        Find the closest star to the right ascension and declination where the
        user has right-clicked and overlay a red marker of radius MARK_RADIUS
        on the APLpy plot. This marker disappears when the user clicks again,
        so only one marker is displayed at all times. Clicks outside of the
        plot (axes) are ignored. This method must be connected to the
        Matplotlib event manager, which is part of the FigureCanvasBase.

        """

        # The attribute 'button' from event is an integer. A left-click is
        # mapped to 1, a middle-click to 2 and a right-click to 3. If we click
        # outside of the axes: event.xdata and event.ydata hold the None value.

        click = (event.xdata, event.ydata)

        if event.button == 3 and None not in click:
            # Get the alpha and delta for these x- and y-coordinates
            coords = self.wcs.all_pix2world(event.xdata, event.ydata, 1)
            star_id = self.db.star_closest_to_world_coords(*coords)[0]
            # LEMONdB.get_star() returns (x, y, ra, dec, epoch, pm_ra, pm_dec, imag)
            ra, dec = self.db.get_star(star_id)[2:4]
            kwargs = dict(layer = self.MARKERS_LAYER,
                          edgecolor = 'red',
                          s = self.MARK_RADIUS)
            self.aplpy_plot.show_markers(ra, dec, **kwargs)
            self.selected_star_id = star_id
            self.goto_button.set_sensitive(True)
            # Pressing Enter activates 'Go to Star'
            self.dialog.set_default_response(gtk.RESPONSE_APPLY)

    def mark_star(self, star_id):
        """ Mark a star in the finding chart.

        Read from the LEMONdB the right ascension and declination of the star
        whose ID is 'star_id' and overlay a green marker of radius MARK_RADIUS
        on the APLpy plot. Any existing markers are removed. The original view
        of the plot is restored (as if the user had clicked the 'Home' button
        in the navigation toolbar), undoing any zooming and panning and taking
        us to the first, default view of the FITS image.

        """

        ra, dec = self.db.get_star(star_id)[2:4]
        kwargs = dict(layer = self.MARKERS_LAYER,
                      edgecolor = '#24ff29',
                      s = self.MARK_RADIUS)
        self.aplpy_plot.show_markers(ra, dec, **kwargs)
        self.navig.home()

        self.selected_star_id = star_id
        self.goto_button.set_sensitive(True)

    def goto_star(self):
        """ Show the details of the selected star.

        This method calls the parent LEMONdB.view_star() with the ID of the
        star currently selected in the finding chart. It adds a notebook page
        with the details of the star, or switches to it if it already exists.

        """

        self.view_star(self.selected_star_id)
        # Now pressing Enter closes the finding chart window
        self.dialog.set_default_response(gtk.RESPONSE_CLOSE)

    def hide(self):
        """ Hide the GTk.Dialog """
        self.dialog.hide()
        self._currently_shown = False

    def show(self):
        """ Display the GTK.Dialog """

        # Until a star is selected, Enter closes the window
        self.dialog.set_default_response(gtk.RESPONSE_CLOSE)
        self.dialog.show()
        self._currently_shown = True

    def is_visible(self):
        """ Return True if the gtk.Dialog is shown, False if hidden """
        return self._currently_shown

    def on_delete_event(self, widget, event):
        """ Callback handler for the 'delete-event' signal.

        Closing a window using the window manager (i.e., clicking on the
        window's close button), by default, causes it to be destroyed, so after
        that there is nothing left to be redisplayed. Instead of destroying it,
        hide the gtk.Dialog. Returns True in order to indicate that the default
        handler is *not* to be called.

        [http://faq.pygtk.org/index.py?req=show&file=faq10.006.htp]

        """

        self.hide()
        return True

    def destroy(self):
        """ Destroy the gtk.Dialog """
        self.preferences_dialog.destroy()
        self.dialog.destroy()

