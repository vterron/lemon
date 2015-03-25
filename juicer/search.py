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

import functools
import gtk
import json
import operator
import re

# LEMON modules
import glade
import methods
import mining
import passband
import util

class AmplitudesSearchMessageWindow(object):

    def get(self, name):
        """ Access a widget in the interface """
        return self.builder.get_object(name)

    def __init__(self, parent_window, builder, db_path, config):

        self.builder = builder
        self.builder.add_from_file(glade.AMPLITUDES_DIALOG)
        self.config = config
        self.miner = mining.LEMONdBMiner(db_path)

        self.dialog = self.get('amplitudes-search-dialog')
        self.dialog.set_resizable(False)
        self.dialog.set_title("Select stars by their amplitudes")
        self.dialog.add_button(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL)
        self.ok_button = self.dialog.add_button(gtk.STOCK_OK, gtk.RESPONSE_OK)

        self.dialog.set_transient_for(parent_window)
        self.dialog.set_position(gtk.WIN_POS_CENTER_ON_PARENT)

        self.progressbar = self.get('search-progress-bar')
        self.increasing_button = self.get('direct-correlation')
        self.decreasing_button = self.get('inverse-correlation')
        self.namplitudes = self.get('amplitudes-how-many')
        self.amplitudes_median = self.get('amplitudes-median')
        self.amplitudes_mean = self.get('amplitudes-mean')
        self.exclude_checkbox = self.get('filter-out-noisy')

        # The radio and spin buttons below, that adjust how stars with noisy
        # amplitudes are excluded, are grouped in a dictionary, so that they
        # all can be easily enabled or disabled by iterating over the items.
        w = {}
        w['nstdevs'] = self.get('comparison-stdevs-how-many')
        w['stdevs_mean'] = self.get('comparison-stdevs-mean')
        w['stdevs_median'] = self.get('comparison-stdevs-median')
        w['min_stdev_ratio'] = self.get('min-amplitude-stdev-ratio')
        self.exclude_widgets = w

        args = 'toggled', self.handle_toggle_exclude_noisy
        self.exclude_checkbox.connect(*args)

        # Set the labels of the radio buttons that select whether amplitudes
        # must increase or decrease to the photometric filters in the database
        pfilters = sorted(self.miner.pfilters)
        letters = [p.letter for p in pfilters]
        label = functools.partial('-'.join)
        self.increasing_button.set_label(label(letters))
        self.decreasing_button.set_label(label(reversed(letters)))

        self.update()

        # Handlers that update the options in the configuration file every time
        # that the value of a widget (radio, check or spin button) is modified

        def save_widget_update(option, func, type_ = int):
            """ Return the function that, when called, updates 'option' in the
            amplitudes search section of the configuration file. 'func' is the
            method used to get the value of the widget, cast to 'type_' """

            def handler(widget):
                self.config.amplset(option, type_(getattr(widget, func)()))
            return handler

        f = save_widget_update
        self.increasing_button.connect('toggled', f('increasing', 'get_active'))
        self.amplitudes_median.connect('toggled', f('use_median', 'get_active'))
        self.namplitudes.connect('output', f('npoints', 'get_value'))
        self.exclude_checkbox.connect('toggled', f('exclude_noisy', 'get_active'))
        w['nstdevs'].connect('output', f('noisy_nstdevs', 'get_value'))
        w['stdevs_median'].connect('toggled', f('noisy_use_median', 'get_active'))
        w['min_stdev_ratio'].connect('output', f('noisy_min_ratio', 'get_value', type_ = float))

    def update(self):
        """ Update the buttons to the values given in the configuration file"""

        config = self.config
        increasing = config.amplint('increasing')
        self.increasing_button.set_active(increasing)
        self.decreasing_button.set_active(int(not increasing))
        use_median = config.amplint('use_median')
        self.amplitudes_median.set_active(use_median)
        self.amplitudes_mean.set_active(int(not use_median))
        self.namplitudes.set_value(config.amplint('npoints'))
        self.exclude_checkbox.set_active(config.amplint('exclude_noisy'))
        w = self.exclude_widgets
        w['nstdevs'].set_value(config.amplint('noisy_nstdevs'))
        stdevs_median = config.amplint('noisy_use_median')
        w['stdevs_median'].set_active(stdevs_median)
        w['stdevs_mean'].set_active(int(not stdevs_median))
        w['min_stdev_ratio'].set_value(config.amplfloat('noisy_min_ratio'))

    def set_fraction(self, fraction):
        """ Fill in the portion of the bar specified by 'fraction'.

        The progress bar is only updated when the percentage varies; if, for
        example, it is 0.971 (97%), setting the fraction to 0.972 (still 97%)
        would only unnecessarily slow down the execution.

        """

        if fraction != self.progressbar.get_fraction():
            with util.gtk_sync():
                self.progressbar.set_fraction(fraction)

    @property
    def description(self):
        """ Return a text description of the search, explaining exactly what
        was made, including the value of all the parameters."""

        increasing_order = self.increasing_button.get_active()
        using_median = self.amplitudes_median.get_active()
        params = dict(order = 'increase' if increasing_order else 'decrease',
                      pfilters = self.increasing_button.get_label(),
                      mode = 'median' if using_median  else 'mean',
                      how_many = int(self.namplitudes.get_value()))

        # These parameters will be shown in bold in the GtkLabel
        # http://faq.pygtk.org/index.py?req=show&file=faq07.003.htp
        def bold(str_):
            return "<b>%s</b>" % str_
        params = dict((k, bold(v)) for k, v in params.iteritems())
        params['field'] = self.miner.field_name # not in bold

        text = \
        "These are the stars in %(field)s whose amplitudes %(order)s with " \
        "%(pfilters)s when the peak and trough of each light curve are " \
        "obtained by taking the %(mode)s of the %(how_many)s highest and " \
        "lowest magnitudes." % params

        if self.exclude_checkbox.get_active():

            w = self.exclude_widgets
            using_median = w['stdevs_median'].get_active()
            params = dict(mode = 'median' if using_median else 'mean',
                          how_many = int(w['nstdevs'].get_value()),
                          ratio = w['min_stdev_ratio'].get_value())
            params = dict((k, bold(v)) for k, v in params.iteritems())

            text += "\n\n" + \
            "Those stars with one or more noisy amplitudes were " \
            "excluded from the search. To determine if an amplitude was " \
            "<i>noisy</i>, we divided it by the %(mode)s of the standard " \
            "deviation of the light curves of the %(how_many)s stars with " \
            "the most similar brightnesses: if the ratio was smaller than " \
            "%(ratio)s, the amplitude was considered noisy." % params

        return text

    def run(self):
        """ Run the dialog window in a recursive loop.

        This method shows the dialog window and allows the user to adjust the
        parameters that will be used in the search. The search can be cancelled
        at any time by clicking 'Cancel', while clicking it when no search is
        in progress closes the window. The 'Ok' button is disabled to avoid
        running two searches in parallel.

        """

        def get_response(widget, response):
            self.response = response
        self.dialog.connect('response', get_response)

        while True:

            self.response = self.dialog.run()
            if self.response == gtk.RESPONSE_OK:

                # Save the description of the search right away, before even
                # starting the search; otherwise, the user could modify some of
                # the parameters while the search is in progress, and thus the
                # description would not match the values used in actuality.
                description = self.description

                increasing = self.increasing_button.get_active()
                exclude_noisy = self.exclude_checkbox.get_active()
                w = self.exclude_widgets

                args = (increasing,
                        int(self.namplitudes.get_value()),
                        self.amplitudes_median.get_active(),
                        exclude_noisy,
                        int(w['nstdevs'].get_value()),
                        w['stdevs_median'].get_active(),
                        w['min_stdev_ratio'].get_value())

                g = self.miner.amplitudes_by_wavelength(*args)
                self.progressbar.set_text("Please wait...")
                self.ok_button.set_sensitive(False)

                pfilters = sorted(self.miner.pfilters)
                id_ = self.miner.id
                field = self.miner.field_name
                args = pfilters, exclude_noisy, description, id_, field
                result = AmplitudesSearchPage(*args)

                nstars = len(self.miner)
                for star_index, star_data in enumerate(g):

                    # Has the user pressed 'Cancel'?
                    if self.response == gtk.RESPONSE_CANCEL:
                        self.progressbar.set_fraction(0.0)
                        self.progressbar.set_text('')
                        self.ok_button.set_sensitive(True)
                        self.response = None
                        break

                    if star_data:
                        star_id = star_data[0]

                        # Sort the list of three-element tuples (photometric
                        # filter, amplitude and comparison standard deviation)
                        # by the first element, as columns must be sorted by
                        # wavelength in the AmplitudesSearchPage.
                        arg = sorted(star_data[-1], key = operator.itemgetter(0))
                        _ , amplitudes, stdevs = zip(*arg)

                        if exclude_noisy:
                            ratios = [a / s for a, s in zip(amplitudes, stdevs)]
                        else:
                            ratios = None

                        result.add(star_id, amplitudes, ratios)

                    fraction = round(star_index / nstars, 2)
                    self.set_fraction(fraction)

                else:
                    self.dialog.destroy()
                    return result

            if self.response in [gtk.RESPONSE_CANCEL, gtk.RESPONSE_DELETE_EVENT]:
                self.dialog.destroy()
                break

    def handle_toggle_exclude_noisy(self, widget):
        """ Enable / disable widgets, as needed, when the 'Filter out those
        stars with one or more noisy amplitudes' check button is toggled """

        checkbox_enabled = widget.get_active()
        for w in self.exclude_widgets.itervalues():
            w.set_sensitive(checkbox_enabled)
