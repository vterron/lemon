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
import lxml.etree
import operator
import re

# LEMON modules
import glade
import methods
import mining
import passband
import util
import xmlparse

class AmplitudesSearchPage(object):
    """ Encapsulates a gtk.HBox with (1) a description of the search for stars
    whose light curve amplitudes are correlated with the wavelength and (2) a
    gtk.GtkTreeView (and the corresponding tree store) with the found stars.
    Use the 'get_window' method to access the gtk.HBox, meant to be used in
    the main gtk.GtkNotebook of the LEMONJuicerGUI class.

    """

    XML_DTD = [
    "",
    "<!DOCTYPE AmplitudesSearchResult [",
    "<!ELEMENT AmplitudesSearchResult (description, field, filter+, star+)>",
    "<!ATTLIST AmplitudesSearchResult include_ratios CDATA  #REQUIRED>",
    "<!ATTLIST AmplitudesSearchResult database_id CDATA  #REQUIRED>",
    "",
    "<!ELEMENT description (#PCDATA)>",
    "<!ELEMENT field (#PCDATA)>",
    "",
    "<!ELEMENT filter (#PCDATA)>",
    "<!ATTLIST filter order CDATA  #REQUIRED>",
    "",
    "<!ELEMENT star (amplitude+)>",
    "<!ATTLIST star id CDATA  #REQUIRED>",
    "<!ELEMENT amplitude (#PCDATA)>",
    "<!ATTLIST amplitude filter CDATA  #REQUIRED>",
    "<!ATTLIST amplitude ratio CDATA  #IMPLIED>",
    "]>",
    ""]

    def __init__(self, pfilters, include_ratios, description, id_, field_name):
        """ Instantiation method for the AmplitudesSearchPage class.

        The 'pfilters' parameter is a sequence with the passband.Passband
        instances of the photometric filters of the amplitudes that will be
        added later with the 'add' method. If 'include_ratios' is True,
        additional columns will be created to store the ratio between each
        amplitude and its comparison standard deviation. The text of the
        GtkLabel on the left column, expected to contain an account of the
        parameters used in the search, is set to 'description'.

        The 'id_' parameter is the unique identifier of the LEMONdB: when the
        object is serialized it is written to the XML file to map the results
        of the search to the database to which they correspond. It is necessary
        to avoid, when we are working with a LEMONdB, loading search results
        for a different one. Lastly, 'field_name' is a string containing the
        name of the observed field, used to suggest a filename in the 'Save
        As...' dialog.

        """

        self.builder = gtk.Builder()
        self.builder.add_from_file(glade.AMPLITUDES_RESULTS)

        self.pfilters = pfilters
        self.include_ratios = include_ratios

        # ID, amplitudes and (optionally) Δ:stdev ratios
        args = [int] + [float] * len(pfilters)
        if include_ratios:
            args += [float] * len(pfilters)
        self.store = gtk.ListStore(*args)

        attrs = ['ID']
        attrs += ["Δ %s" % p.letter for p in pfilters]
        if include_ratios:
            attrs += ["Ratio %s" % p.letter for p in pfilters]

        self.view = self.builder.get_object('amplitudes-search-found-view')
        for index, title in enumerate(attrs):
            render = gtk.CellRendererText()
            column = gtk.TreeViewColumn(title, render, text = index)
            column.props.resizable = False
            column.set_sort_column_id(index)
            self.view.append_column(column)

        self.description = self.builder.get_object('search-description-label')
        self.description.set_label(description)
        self.id = id_
        self.field_name = field_name

        # The button to export the search results to an XML file
        save_button = self.builder.get_object('save-button')
        save_button.connect('clicked', self.export_to_file)


    def add(self, star_id, amplitudes, ratios = None):
        """ Append a new row to the store.

        The 'star_id' must be an integer, the ID of the star, 'amplitudes' a
        sequence of floats with its amplitudes (in the same photometric filters
        and order that were passed to __init__), and 'ratios' a second sequence
        with the ratio between each amplitude and its comparison standard
        deviations. Note, however, that these ratios are ignored if
        'include_ratios' was not set to True at instantiation time.

        """

        row = [star_id] + list(amplitudes)
        if self.include_ratios:
            row += list(ratios)
        self.store.append(row)

    def get_window(self):
        """ Return a gtk.HBox with the found stars """
        self.view.set_model(self.store)
        return self.builder.get_object('amplitudes-search-result')

    def get_label(self, order):
        """ Return the title (text of the label) of the search.

        The method returns a string, 'Amplitudes', followed by the value given
        in 'order' but expressed in Roman numerals (e.g., 'Amplitudes IV' and
        'Amplitudes VII', for values 4 and 7, respectively). Searches should
        be numbered sequentially to allow the user to refer to them easily.

        The returned title is also stored in the 'label' attribute. Note that
        this attribute is not initialized at __init__, so attempting to access
        it before this method is called for the first time will unavoidably
        raise the AttributeError exception.

        """

        self.label = 'Amplitudes %s' % methods.int_to_roman(order)
        return self.label

    def dump(self, path):
        """ Write JSON serialization to a file, overwrite if it exists """

        data = dict(database_id = self.id,
                    include_ratios = self.include_ratios,
                    description = self.description.get_label(),
                    field = self.field_name,
                    stars = [])

        filters_order = {}
        # Map each photometric filter to its order
        for index, pfilter in enumerate(self.pfilters):
            filters_order[str(pfilter)] = index
        data['filters_order'] = filters_order

        for row in self.store:

            # Needed as gtk.ListStore doesn't allow slices as indexes
            row = list(row)

            id_ = row[0]
            star_data = dict(id = id_, amplitudes = [])

            # When the ratios (between each amplitude and its comparison
            # standard deviation) are not included, the columns in the
            # gtk.ListStore, except for the first one, only contain the
            # different amplitudes of each star. When ratios are present, they
            # are listed after the amplitudes and in the same order: we need to
            # pair them up; i.e., match the amplitude for the first filter with
            # its ratio, the second amplitude with the second ratio, and so on.

            if not self.include_ratios:
                amplitudes = row[1:]
                ratios = [None] * len(amplitudes)
            else:
                assert not len(row[1:]) % 2
                n = int((len(row) - 1) / 2) # how many amplitude-ratio pairs
                amplitudes = row[1:1+n]
                ratios = row[1+n:]
            assert len(self.pfilters) == len(amplitudes) == len(ratios)

            args = self.pfilters, amplitudes, ratios
            for pfilter, amplitude, ratio in zip(*args):
                record = dict(filter = str(pfilter),
                              value = amplitude)
                if self.include_ratios:
                    record['ratio'] = ratio

                star_data['amplitudes'].append(record)

            data['stars'].append(star_data)

        with open(path, 'wt') as fd:
            kwargs = dict(indent=2, sort_keys=True)
            json.dump(data, fd, **kwargs)

    def export_to_file(self, widget):
        """ Let the user choose where the object is saved to an XML file.

        Prompt a gtk.FileChooserDialog to let the user browse for the location
        where the XML file with the representation of the AmplitudesSearchPage
        will be saved, and write the file to disk (see 'xml_dump' method) only
        if 'Save' is clicked. The dialog suggests a filename based on the name
        of the observed field and the title of the search, but only if the
        'get_label' method has been called at least once.

        """

        kwargs = dict(title = "Export light curve to...",
                      action = gtk.FILE_CHOOSER_ACTION_SAVE,
                      buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                                 gtk.STOCK_SAVE, gtk.RESPONSE_OK))

        with util.destroying(gtk.FileChooserDialog(**kwargs)) as chooser:

            # Ask the user to confirm overwriting an existing file
            chooser.set_do_overwrite_confirmation(True)

            # If AmplitudesSearchPage.get_label has been called, use the label
            # that it returned, and the name of the observed field, to suggest
            # a filename for the XML file. Except for the Roman numerals, it is
            # in lower case (e.g., "ngc2264_amplitudes_IV.xml" is suggested for
            # "NGC2264" and "Amplitudes IV"). Do not suggest anything if
            # get_label has not been called.

            try:
                # 'Amplitudes IV' to 'amplitudes-IV'
                title = re.sub(r'\s', '_', self.label)
                title = title[0].lower() + title[1:]
                # 'IC 5146' to 'ic_5146'
                field = re.sub(r'\s', '_', self.field_name.lower())
                args = field, title
                filename = '%s_%s.xml' % args
                chooser.set_current_name(filename)
            except AttributeError:
                pass

            response = chooser.run()
            if response == gtk.RESPONSE_OK:
                self.xml_dump(chooser.get_filename())

    @classmethod
    def xml_load(cls, path):
        """ Load an AmplitudesSearchPage object from an XML file.

        The opposite of toxml, this class method parses 'path', de-serializing
        the XML representation of an AmplitudesSearchPage object and returning
        it. The XML files generated by LEMON are always standalone documents,
        meaning that the Document Type Definitions (DTD), defining the document
        structure with a list of legal elements, are also included. If the XML
        file cannot be validated, the appropriate exception will be raised.

        """

        # Raise a Python built-in exception if the file cannot be read; these
        # are more understandable than those that lxml would throw, with error
        # messages such as "failed to load external entity".
        with open(path, 'r') as _: pass
        root = lxml.etree.parse(path).getroot()

        elements = list(root.iterchildren(tag = 'description'))
        assert len(elements) == 1
        description = elements[0].text

        elements = list(root.iterchildren(tag = 'field'))
        assert len(elements) == 1
        field_name = elements[0].text

        # Encapsulate the photometric filters as Passband instances,
        # and sort them by the value of their 'order' XML attribute.
        filter_elements = root.iterchildren(tag = 'filter')
        orders = dict((e.text, e.get('order')) for e in filter_elements)
        pfilters = [passband.Passband(k) for k in orders.iterkeys()]
        pfilters.sort(key = lambda x: orders[str(x)])

        id_ = root.get('database_id')
        include_ratios = root.get('include_ratios').title() == str(True)
        result = cls(pfilters, include_ratios, description, id_, field_name)

        for star_element in root.iterchildren(tag = 'star'):
            star_id = int(star_element.get('id'))

            # The order of the amplitudes must match that of the filters
            elements = list(star_element.iterchildren(tag = 'amplitude'))
            elements.sort(key = lambda x: orders[x.get('filter')])

            amplitudes = []
            ratios = []

            for e in elements:
                amplitudes.append(float(e.text))
                ratio = float(e.get('ratio')) if include_ratios else None
                ratios.append(ratio)

            result.add(star_id, amplitudes, ratios = ratios)

        return result


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
            self.progressbar.set_fraction(fraction)
            # Ensure rendering is done immediately
            while gtk.events_pending():
                gtk.main_iteration()

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


def amplitudes_search(parent_window, builder, db, config):
    dialog = AmplitudesSearchMessageWindow(parent_window, builder, db, config)
    return dialog.run()

