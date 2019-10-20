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

import astropy.time
import datetime
import matplotlib
import matplotlib.dates

# LEMON modules
import snr

def split_by_diff(iterable, delta = 3):
    """ Split a sequence by the difference between consecutive elements.

    The method returns an interator over the result of splitting the input
    sequence, stopping each sub-sequence at the element at which its difference
    with the next one is greater than 'delta'. In other words: the difference
    between consecutive elements of each of the returned sublists will be
    smaller than or equal to delta.

    For example, split_by_diff([1, 2, 3, 8, 9, 15], delta = 3) returns an
    iterator over three lists: [1, 2, 3], [8, 9] and [15]

    """

    differences = numpy.diff(iterable)
    sublist_indexes = numpy.where(differences > delta)[0]

    sublists = []
    iterable = list(iterable)  # work on a copy
    for index in reversed(sublist_indexes):
        sublists.append(iterable[index + 1:])
        del iterable[index + 1:]
    sublists.append(iterable)
    return reversed(sublists)

def curve_plot(figure, curve, marker = 'o', color = '',
               airmasses = None, delta = 3600 * 3, julian = False):
    """ Plot the light curve as a subplot of the matplotlib Figure.

    The method removes from 'figure' all the existing axes and adds a new one,
    plotting the differential magnitudes as a function of time, and including
    the error bars which can be derived from the signal-to-noise ratios. The
    original matplotlib Figure instance is modified, so nothing is returned.

    Keyword arguments:
    marker - format string character, which is passed down to matplotlib, that
             controls the line style. The default value, 'o', defines circle
             markers, while '-.', for example, uses the dash-dot line style.
             Please refer to the matplotlib documentation for details.
    color - the string which determines the color of the markers. Abbreviations
            such as 'g', full names ('green') and hexadecimal strings
            ('#008000'), among others, can be given. With an empty string,
            the first element of the matplotlib cycle of colors is used.
    airmasses - a dictionary mapping each Unix time of the light curve to its
                airmass. In case one or more of the Unix times of the curve
                cannot be found among the keys of the dictionary, the KeyError
                exception will be raised. If None is used, the airmasses of
                the light curve will not be plotted.
    delta - the maximum difference, in seconds, between consecutive points of
            the light curve if their airmasses are to be connected with blue
            dashed lines. This argument prevents the airmasses from different
            nights to be connected by the same line, which does not add any
            information and only clutters the plot unnecessarily.
    julian - plot the points of the light curve using Julian dates (JD) (such
             as 2456877.9660300924) instead of strings representing the date
             and time. For the latter, matplotlib automatically figures out
             the best format to use (for example, 'Jan 02 2012' or '08:15:31')
             depending of the date range of the light curve.

    """

    # Remove existing subplots plots, if any
    figure.clear()
    for axes in figure.get_axes():
        figure.delaxes(axes)

    unix_times, magnitudes, snrs = zip(*curve)

    # By default, convert the Unix timestamps to datetime objects, but use
    # Julian days in case the 'julian' keyword argument evaluates to True.
    if julian:
        dates = astropy.time.Time(unix_times, format = 'unix').jd
    else:
        dates = [datetime.datetime.utcfromtimestamp(x) for x in unix_times]

    # For each signal-to-noise ratio we want the equivalent error in mags;
    # that is, both the maximum and minimum values induced by the noise.
    negative_errors = []
    positive_errors = []
    for s in snrs:
        min_error, max_error = snr.snr_to_error(s)
        negative_errors.append(abs(min_error))
        positive_errors.append(max_error)

    ax1 = figure.add_subplot(111)
    ax1.errorbar(dates, magnitudes, color = color, marker = marker,
                 fmt = 's', yerr = (negative_errors, positive_errors))
    ax1.set_ylabel(r"$\Delta$ magnitude") # mathtext
    ax1.grid(True)

    # Use 5% of the height of the plot as the upper and bottom margin of the
    # magnitudes, so that the curves do not touch the border. Also, reverse the
    # y-axis, as we want the light curve to go up when the star is brighter;
    # that is, as the differential magnitude becomes smaller.
    max_mag = max(magnitudes)
    min_mag = min(magnitudes)
    margin_mag = (max_mag - min_mag) * 0.05
    ax1.set_ylim(max_mag + margin_mag, min_mag - margin_mag)

    # We also want a similar margin on the x-axis; the dates of the
    # observations should not touch the border of the plot either.
    elapsed_seconds = (unix_times[-1] - unix_times[0])
    margin_seconds = elapsed_seconds * 0.05

    if julian:
        kwargs = dict(format = 'sec')
        margin_delta = astropy.time.TimeDelta(margin_seconds, **kwargs).jd
    else:
        margin_delta = datetime.timedelta(seconds = margin_seconds)

    ax1.set_xlim(dates[0] - margin_delta, dates[-1] + margin_delta)

    # In interactive navigation, show the full date of the x-locations instead
    # of formatting them the same way the tick labels are (e.g., "Jan 04 2012")
    xdata_formatter = matplotlib.dates.DateFormatter("%a %b %d %H:%M:%S %Y")
    ax1.fmt_xdata = xdata_formatter

    # Datetime instances are often squashed together when used for the x-axis
    # tick labels. Until matplotlib can handle this better, we can fix this by
    # setting the maximum number of ticks desired (AutoDateLocator, 'maxticks'
    # keyword argument, which by default is None; seven seems to work fine with
    # our minimum width of 768 pixels). This keyword was not available until
    # matplotlib version 1.0.0. This only applies, of course, when we use
    # datetime objects.

    if not julian and matplotlib.__version__ >= '1.0.0':
        xticks_locator = matplotlib.dates.AutoDateLocator(maxticks = 8)
        xticks_formatter = matplotlib.dates.AutoDateFormatter(xticks_locator)
        ax1.xaxis.set_major_locator(xticks_locator)
        ax1.xaxis.set_major_formatter(xticks_formatter)

    # Now (optionally) plot the airmasses
    if airmasses:

        ax2 = ax1.twinx()
        ax2.fmt_xdata = xdata_formatter

        if not julian and matplotlib.__version__ >= '1.0.0':
            ax2.xaxis.set_major_locator(xticks_locator)
            ax2.xaxis.set_major_formatter(xticks_formatter)

        ax2.set_ylabel('Airmass', rotation = 270)
        ax2.yaxis.labelpad = 17 # padding between axis and label
        periods = split_by_diff(unix_times, delta = delta)
        for period_unix_times in periods:

            if julian:
                utimes = period_unix_times
                kwargs = dict(format = 'unix')
                period_dates = astropy.time.Time(utimes, **kwargs).jd
            else:
                func = datetime.datetime.utcfromtimestamp
                period_dates = [func(x) for x in period_unix_times]

            period_airmasses = [airmasses[x] for x in period_unix_times]
            ax2.plot(period_dates, period_airmasses,
                     color = color, linestyle = '--')

        # Let the airmasses have a little margin too. We cannot find the
        # maximum and minimum values directly among the keys of the dictionary
        # of airmasses as it may have information for more Unix times than
        # those in the light curve.
        plotted_airmasses = [airmasses[x] for x in unix_times]
        max_airmass = max(plotted_airmasses)
        min_airmass = min(plotted_airmasses)
        margin_airmass = (max_airmass - min_airmass) * 0.05
        ax2.set_ylim(min_airmass - margin_airmass, max_airmass + margin_airmass)

        # Margins on the x-axis must be set again after plotting the right
        # y-axis; otherwise, the margins we set for ax1 will be ignored.
        ax2.set_xlim(dates[0] - margin_delta, dates[-1] + margin_delta)

