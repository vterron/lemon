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

import atexit
import aplpy
import gtk
import methods
import os
import pyfits

import matplotlib.figure
from matplotlib.backends.backend_gtkagg \
     import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg \
     import NavigationToolbar2GTKAgg as NavigationToolbar

# LEMON modules
import glade

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
        elif response in (gtk.RESPONSE_CLOSE, gtk.RESPONSE_DELETE_EVENT):
            # Untoggle gtk.ToggleToolButton in toolbar (main window)
            self.toggle_toolbar_button(False)
            self.hide()
        else:
            raise ValueError("unexpected dialog response")

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
        with pyfits.open(path) as hdu:
            data = hdu[0].data

        self.aplpy_plot = aplpy.FITSFigure(path, figure = self.figure)
        self.figure.canvas.mpl_connect('button_press_event', self.mark_closest_star)
        self.aplpy_plot.show_grayscale()

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
            star_id = self.db.star_closest_to_image_coords(*click)[0]
            # LEMONdB.get_star() returns (x, y, ra, dec, imag)
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

        Read from the LEMONdB the x- and y-image coordinates of the star whose
        ID is 'star_id' and overlay a green marker of radius MARK_RADIUS on the
        APLpy plot. Any existing markers are removed. The original view of the
        plot is restored (as if the user had clicked the 'Home' button in the
        navigation toolbar), undoing any zooming and panning and taking us to
        the first, default view of the FITS image.

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
        self.dialog.destroy()

