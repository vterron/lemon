#! /usr/bin/env python

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

description = """
This module implements several data mining and data analysis algorithms to
identify potentially interesting objects, such as, for example, those whose
light curves are strongly correlated in different photometric filters or those
that have the lowest standard deviation in their curves and are therefore the
most constant.

"""

import collections
import datetime
import itertools
import numpy
import operator
import scipy.stats

# LEMON modules
import database

class NoStarsSelectedError(ValueError):
    """ Raised when no stars can be returned by the LEMONdBMiner """
    pass


class LEMONdBMiner(database.LEMONdB):
    """ Interface to identify potentially interesting stars in a LEMONdB """

    @staticmethod
    def _ascii_table(headers, table_rows, sort_index = 1, descending = True,
                    ndecimals = 8, dates_columns = None):
        """ Format data as an ASCII table.

        Returns a string with the data in 'table_rows' formatted as an ASCII
        table, sorted by the values of the sort_index-th column. 'table_rows'
        must be a sequence (table rows) of sequences (table columns), all of
        the latter having the same number of elements; use Nones for missing
        values and therefore empty cells in the table.

        In addition to the data received in 'table_rows', a first column is
        automatically prefixed with the index of the row. There is, thus, no
        need to manually include it. This numbering starts at zero, as we
        belive in spreading the word of Dijkstra's teachings.

        Keyword arguments:
        sort_index - the column by which the rows of the table are sorted.
        descending - sort the columns into decreasing order?
        ndecimals - number of decimals shown for real numbers.
        dates_columns - the indexes of the columns, containing seconds, to be
                        formatted as a string indicating hours, minutes and
                        seconds (e.g., '0:11:06'). The number of days is also
                        included in case the number of seconds is greater than
                        or equal to twenty-four hours ('17 days, 6:12:00').

        """

        if not dates_columns:
            dates_columns = ()

        if len(set(len(row) for row in table_rows)) != 1:
            raise ValueError("number of columns must be the same for all rows")

        if len(headers) != len(table_rows[0]):
            raise ValueError("elements in 'headers' must equal that of columns")

        # Sort the table by the specified column
        if sort_index is not None:
            table_rows.sort(key = operator.itemgetter(sort_index),
                            reverse = descending)

        # Give format to all the data in the table as a string, as (1) a date,
        # (2) a real number with the specified number of decimals or (3) by
        # simply casting it to a string. The formatted data is stored as a
        # two-level dictionary, where the first index maps to the row of the
        # table and the second to the column.

        table_data = collections.defaultdict(dict)
        table_data[0][0] = '' # top-left 'corner' of table left empty

        # Populate the table with the headers...
        for column_index in range(len(headers)):
            table_data[0][column_index + 1] = headers[column_index]

        # ... and with all the data
        for row_index, row in enumerate(table_rows):
            table_data[row_index + 1][0] = str(row_index)
            for column_index, value in enumerate(row):
                if value is None:
                    table_cell = ''
                elif column_index in dates_columns:
                    table_cell = str(datetime.timedelta(seconds = value))
                elif isinstance(value, float):
                    table_cell = '%.*f' % (ndecimals, value)
                else:
                    table_cell = str(value)
                table_data[row_index + 1][column_index + 1] = table_cell

        # Determine how many characters are needed for each column. Two
        # different measures are made: the width of the entire column (i.e, the
        # maximum length among the cells that share the same column) and the
        # width of only the data (the maximum length among the cells that share
        # the same column, excluding the header). These values will be used to
        # right-justify the values and at the same time center them.

        widths = {}
        data_widths = {}
        for column_index in table_data[0].iterkeys():
            column_values = (table_data[row][column_index]
                             for row in table_data.iterkeys())
            column_lengths = [len(str(x))
                              for x in column_values if x is not None]
            # The maximum width is increased by two, so that no value in the
            # table touches the vertical cell borders (the '|' characters)
            widths[column_index] = max(column_lengths) + 2
            data_widths[column_index] = max(column_lengths[1:]) # only the data


        # Finally, populate the ASCII table with the headers...
        header_values = table_data[0].itervalues()
        header_str = '|'.join([x.center(widths[index])
                               for index, x in enumerate(header_values)])
        output = header_str + '\n' + '-' * len(header_str) + '\n'

        # ... and with the rest of the data
        for row_index in sorted(table_data.iterkeys()):
            if not row_index:  # headers (index == 0) already formatted
                continue
            row_values = table_data[row_index].values()
            output += row_values[0].ljust(data_widths[0]).center(widths[0]) + '|'
            output += '|'.join(x.rjust(data_widths[index]).center(widths[index])
                                for index, x in enumerate(row_values) if index)
            output += '\n'

        return output

    def sort_by_period_similarity(self, minimum = 2, normalization = 'max'):
        """ Sort the stars by their period similarity.

        This method sorts the stars in a LEMONdB by the similarity between
        their periods in the different photometric filters in which they were
        observed. How similar the periods of a star are is determined by
        computing their standard deviation after normalizing them. A list of
        two-element tuples, with the ID of the star and the aforementioned
        standard deviation, respectively, is returned, sorted in increasing
        order by the latter.

        The NoStarsSelectedError exception is raised in case no stars can be
        returned, either because 'minimum' is set to a too high value or
        because the database has no information about the star periods.

        Keyword arguments:
        minimum - ignore stars whose periods have been calculated in fewer than
                  this number of photometric filters. As the standard deviation
                  of a single element is zero, 'minimum' must have a value of
                  at least two, since otherwise the first returned stars would
                  almost always be those with a single period, which most of
                  the time are not those of interest to us.
        normalization - defines how the periods of a star are normalized (i.e.,
                        the value by which they are divided) before their
                        standard deviation is calculated. May be any of the
                        following, self-explanatory values: 'max', 'median'
                        or 'mean'.

        """

        if minimum < 2:
            raise ValueError("a minimum of at least two periods is needed")

        if normalization == 'max':
            norm_func = numpy.max
        elif normalization == 'mean':
            norm_func = numpy.average
        elif normalization == 'median':
            norm_func = numpy.median
        else:
            raise ValueError("unrecognized normalization method")

        periods_stdevs = []
        for star_id in self.star_ids:
            star_periods = self.get_periods(star_id)
            if len(star_periods) < minimum:
                continue

            normalized_periods = star_periods / norm_func(star_periods)
            periods_stdev = float(numpy.std(normalized_periods, dtype = self.dtype))
            periods_stdevs.append((star_id, periods_stdev))

        if not periods_stdevs:
            msg = "no stars with at least %d known periods"
            raise NoStarsSelectedError(msg % minimum)

        return sorted(periods_stdevs, key = operator.itemgetter(1))

    def period_similarity(self, how_many, minimum = 2,
                          normalization = 'max', ndecimals = 8):
        """ Return an ASCII table with the stars with the most similar periods.

        This is a wrapper around the LEMONdBMiner.sort_by_period_similarity and
        LEMONdBMiner._ascii_table methods, to avoid having to invoke them
        sequentially. In other words: first, the 'how_many' stars with the most
        similar periods are found, and then, sorted increasingly by the
        standard deviation of their normalized periods, returned in an ASCII
        table. Along with this value, the period of each star in each
        photometric filter is shown, formatted such as '0:11:06' or
        '17 days, 6:12:00'.

        Keyword arguments:
        minimum - ignore stars whose periods have been calculated in fewer than
                  this number of photometric filters. As the standard deviation
                  of a single element is zero, 'minimum' must have a value of
                  at least two, since otherwise the first returned stars would
                  almost always be those with a single period, which most of
                  the time are not those of interest to us.
        normalization - defines how the periods of a star are normalized (i.e.,
                        the value by which they are divided) before their
                        standard deviation is calculated. May be any of the
                        following, self-explanatory values: 'max', 'median'
                        or 'mean'.
        ndecimals - number of decimals shown for real numbers.

        """

        # The 'how_many' stars with the most similar normalized periods
        most_similar = \
            self.sort_by_period_similarity(normalization = normalization,
                                           minimum = minimum)[:how_many]
        header = []
        header.append("Star")
        header.append("Stdev. Norm.")
        for pfilter in self.pfilters:
            header.append("Period %s" % pfilter.letter)

        table_columns = []
        for star_id, star_stdev in most_similar:
            column = [star_id, star_stdev]
            for pfilter in self.pfilters:

                # LEMONdB.get_period returns a two-element tuple (period,
                # step), but a single None if the star period is unknown
                star_period = self.get_period(star_id, pfilter)
                if star_period is None:
                    period_secs = step_secs = None
                else:
                    period_secs, step_secs = star_period

                column.append(period_secs)
            table_columns.append(column)

        return LEMONdBMiner._ascii_table(header, table_columns, sort_index = 1,
                                        descending = False, ndecimals = ndecimals,
                                        dates_columns = (2, 3, 4))

    def match_bands(self, star_id, first_pfilter, second_pfilter,
                    delta = 3600 * 1.5):
        """ Return the paired magnitudes from two light curves of the star.

        The method takes the light curves of a star in two photometric filters
        and 'matches' each point of the first curve to the closest point in
        time in the second. If they are not more than 'delta' seconds apart,
        they are considered to be paired. Returns a list of two-element tuples
        which contain the two matched magnitudes. If the star has no light
        curve in one (or both) of the light curves, None is returned.

        """

        try:
            # Extract the light curves of the star in the two filters, as
            # (Unix time, magnitude, S/N) three-element tuples. NoneType
            # is returned (and therefore the conversion to list raises
            # TypeError) if the star has no light curve in one filter.
            phot1_points = list(self.get_light_curve(star_id, first_pfilter))
            phot2_points = list(self.get_light_curve(star_id, second_pfilter))

            matches = []
            for point1 in phot1_points:
                # Raises ValueError for empty sequences
                closest = min(phot2_points,
                              key = lambda x: abs(point1[0] - x[0]))

                # Second element of the tuple (index = 1) is the magnitude
                if abs(point1[0] - closest[0]) < delta: # we have a match!
                    matches.append((point1[1], closest[1]))

            return matches

        except (TypeError, ValueError):
            return None

    def star_correlation(self, star_id, first_pfilter, second_pfilter,
                         min_matches = 10, delta = 3600 * 1.5):
        """ Return the correlation between two photometric filters of a star.

        Take the output of LEMONdBMiner.match_bands, which pairs the magnitudes
        from two light curves of the star, and compute the least-squares
        regression of the matched elements. Returns a four-element tuple,
        containing (a) the ID of the star, (b) the slope of the regression
        line, (c) the correlation coefficient and (d) the number of 'paired'
        elements that were used for the regression.

        If the star has absolutely no photometric information in one or both
        filters, or if less than 'min_matches' points are paired, NoneType is
        returned.

        """

        matches = self.match_bands(star_id, first_pfilter, second_pfilter,
                                   delta = delta)

        if not matches:  # either None or empty
            return None
        if len(matches) < min_matches:
            return None

        match_x, match_y = zip(*matches)
        slope, intercept, r_value, p_value, std_err = \
            scipy.stats.linregress(match_x, match_y)
        return star_id, slope, r_value, len(matches)

    def band_correlation(self, how_many, sort_index = 0, delta = 3600 * 1.5,
                         min_matches = 10, ndecimals = 8):
        """ Return an ASCII table with the most correlated stars.

        The method takes the combinations of two elements out of all the
        photometric filters in the database and, for each one, computes the
        correlation coefficient of the least-squares regression for the two
        light curves and takes the 'how_many' stars with the highest value.
        Not all the points of the curve are used for the regression, but only
        those which are at most 'delta' seconds apart. Those stars whose light
        curves have fewer than 'min_matches' are ignored. The returned table
        has the stars sorted in decreasing order by the correlation coefficient
        of the 'sort_index'-th combination of photometric filters.

        The NoStarsSelectedError exception is raised in case no stars can be
        returned, either because 'min_matches' is set to a too high value or
        because the database has no information on the curves of the stars.

        """

        pfilter_combinations = list(itertools.combinations(self.pfilters, 2))

        # First, calculate the correlation for all the stars between the
        # 'sort_index'-th combination of photometric filters (the one by which
        # the table will be sorted) and take the 'how_many' most correlated. In
        # this manner, the correlation in the other filters will only have to
        # be calculated for these 'how_many' stars, instead of for all the
        # stars in the database.

        main_combination = pfilter_combinations[sort_index]
        main_correlation = []

        for star_id in self. star_ids:
            kwargs = {'min_matches' : min_matches, 'delta' : delta}
            star_corr = self.star_correlation(star_id, *main_combination, **kwargs)
            if star_corr:
                main_correlation.append(star_corr)

        if not main_correlation:
            msg = "no stars with at least %d correlated points in %s-%s"
            raise NoStarsSelectedError(msg % ((min_matches,) + main_combination))

        # Sort the stars by the correlation between these two photometric
        # filters, decreasing order, and take the best 'how-many' stars. The
        # correlation coefficient is the third element of the tuples, while
        # the ID of each star is the first one (indexes 2 and 0, respectively)
        main_correlation.sort(key = operator.itemgetter(2), reverse = True)
        most_correlated = main_correlation[:how_many]
        most_correlated_ids = [x[0] for x in most_correlated]

        # A two-level dictionary, mapping the combination of filters (first
        # level) to the ID of the star (second level), which in turn maps to a
        # two-element tuple: (a) the correlation coefficient and (b) the number
        # of elements that were paired and used in the linear regression.
        correlations = collections.defaultdict(dict)
        for star_id, slope, r_value, nmatches in most_correlated:
            correlations[main_combination][star_id] = (r_value, nmatches)

        # Find the correlation coefficients in the other filters, but only for
        # these 'how_many', most correlated stars that we have just identified.
        for index, secondary_combination in enumerate(pfilter_combinations):
            if index == sort_index:
                continue # correlations already calculated!

            for star_id in most_correlated_ids:
                star_corr = self.star_correlation(star_id, *secondary_combination,
                                                  min_matches = min_matches,
                                                  delta = delta)
                if star_corr:
                    r_value, nmatches = star_corr[-2:]
                else:
                    # NoneType leaves table cells empty
                    r_value, nmatches = None, None
                correlations[secondary_combination][star_id] = (r_value, nmatches)

        # Now that all the correlations have been calculated, store the
        # information in a list of lists, properly sorted, and create the
        # ASCII table which displays all the data.
        header = []
        header.append('Star')
        for first_pfilter, second_pfilter in pfilter_combinations:
            pfilters_str = '%s-%s' % (first_pfilter.letter, second_pfilter.letter)
            header.append("r-value %s" % pfilters_str)
            header.append('Matches')

        table_columns = []
        for star_id in most_correlated_ids:
            column = [star_id]
            for pfilters_pair in pfilter_combinations:
                r_value, nmatches = correlations[pfilters_pair][star_id]
                column.append(r_value)
                column.append(nmatches)
            table_columns.append(column)

        return LEMONdBMiner._ascii_table(header, table_columns, sort_index = None,
                                         descending = True, ndecimals = ndecimals)

    @staticmethod
    def dump(path, data, decimals = 8):
        """ Dump a sequence of floating point numbers to a text file.

        The method iterates over 'data', each one of whose elements must be a
        sequence of floating point numbers, and saves them to the 'path' text
        file. Each iterable in 'data' is written to a different row, with its
        elements tab-separated. If it already exists, the output file is
        silently overwritten.

        Keyword arguments:
        ndecimals - number of decimals with which real numbers are saved.

        """

        with open(path, 'wt') as fd:
            for row in data:
                fd.write('\t'.join(('%.*f' % (decimals, column)
                                    for column in row)) + '\n')

    def sort_by_curve_stdev(self, pfilter, minimum = 10):
        """ Sort the stars by the standard deviation of their light curves.

        Take the light curves of the stars in a photometric filter and compute
        the standard deviation of their points. Returns a list of two-element
        tuples, with the ID of each star and the said standard deviation,
        respectively, sorted in increasing order by the latter. Those stars
        with fewer than 'minimum' points in their light curve are ignored.

        The NoStarsSelectedError exception is raised in case no stars can be
        returned, either because 'minimum' is set to a too high value or
        because the database has no information on the curves of the stars.

        """

        curves_stdevs = []
        for star_id in self.star_ids:
            star_curve = self.get_light_curve(star_id, pfilter)
            # NoneType returned if the star has no light curve
            if star_curve is None or len(star_curve) < minimum:
                continue
            curves_stdevs.append((star_id, star_curve.stdev))

        if not curves_stdevs:
            msg = "no light curves with at least %d points in %s"
            raise NoStarsSelectedError(msg % (minimum, pfilter))
        return sorted(curves_stdevs, key = operator.itemgetter(1))

    def curve_stdev(self, how_many, sort_index = 0,
                    minimum = 10, ndecimals = 8):
        """ Return an ASCII table with the most constant stars.

        The method computes the standard deviation of the light curve of each
        star in all the photometric filters and takes the 'how_many' stars for
        which this value is lowest (that is, the most constant stars) in the
        'sort_index'-th photometric filter. The order of the photometric
        filters is given by their effective wavelength midpoint. Those stars
        whose light curves have fewer than 'minimum' points are ignored. The
        returned table has the stars sorted in increasing order by the standard
        deviation of the light curve in the 'sort-index'-th filter.

        Keyword arguments:
        ndecimals - number of decimals shown for real numbers.

        """

        # Sort the stars by the standard deviation of their light curve in the
        # 'sort-index'-th photometric filter, ascending order, and take the
        # best 'how-many' stars. The first element of each tuple is the ID of
        # the star; the second is the standard deviation of its light curve.
        main_pfilter = self.pfilters[sort_index]
        most_similar = \
            self.sort_by_curve_stdev(main_pfilter, minimum = minimum)[:how_many]
        most_similar_ids = [x[0] for x in most_similar]

        # A two-level dictionary, mapping the photometric filter (first level)
        # to the ID of the star (second level), which in turn maps to the
        # standard deviation of the corresponding light curve.
        stdevs = collections.defaultdict(dict)
        for star_id, curve_std in most_similar:
            stdevs[main_pfilter][star_id] = curve_std

        # Find the standard deviation of the light curves in the other
        # photometric filters, but only for these 'how_many', most similar
        # stars that we have just identified.
        for index, pfilter in enumerate(self.pfilters):
            if index == sort_index:
                continue # stdevs already calculated!

            for star_id in most_similar_ids:
                star_curve = self.get_light_curve(star_id, pfilter)
                # NoneType returned if the star has no light curve
                if star_curve is None or len(star_curve) < minimum:
                    stdevs[pfilter][star_id] = None
                else:
                    stdevs[pfilter][star_id] = star_curve.stdev

        header = []
        header.append('Star')
        for pfilter in self.pfilters:
            header.append("%s Stdev" % pfilter)

        table_columns = []
        for star_id in most_similar_ids:
            column = [star_id]
            for pfilter in self.pfilters:
                curve_stdev = stdevs[pfilter][star_id]
                column.append(curve_stdev)
            table_columns.append(column)

        # sort_index + 1 because the first column contains the star IDs
        return LEMONdBMiner._ascii_table(header, table_columns,
                                         sort_index = sort_index + 1,
                                         descending = False,
                                         ndecimals = ndecimals)

