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
This module receives as input a LEMON database and finds the period in each
photometric filter of each of the stars using the string-length method proposed
by Dworetsky in 1983. All the periods in the range ]0,n/2], with n being the
time elapsed between the first and last observations of each star, are
evaluated. The reason for this is the Nyquist-Shannon sampling theorem,
according to which at least two samples per period are necessary for detecting
a certain frequency.

"""

import copy
import logging
import optparse
import multiprocessing
import numpy
import os
import shutil
import socket
import sys
import time

# LEMON modules
import database
import defaults
import methods
import style

class StringLengthFinder(object):
    """ The implementation of the method proposed by M. M. Dworetsky in his
    1983 paper 'A period-finding method for sparse randomly spaced observations
    or 'How long is a piece of string'. This is, in short, a force-brute
    approach, in which the optimal period is that which minimizes the sum of
    the lengths of line segments joining successive points in a phase diagram.

    If in doubt, you may refer to the original paper:
    http://adsabs.harvard.edu/abs/1983MNRAS.203..917D

   """

    def __init__(self, size, dtype = numpy.longdouble):
        """ Instantiation method for the StringLengthFinder class.

        The 'size' argument is the number of points of the light curve of the
        star whose period will be computed using the string-length method. It
        must be known beforehand as we use internally a bi-dimensional NumPy
        array to compute the photometric information. In this array, the first
        row (index = zero) stores the Unix times of the observations and the
        second one (index = one) the differential magnitudes.

        """
        self._data = numpy.empty((2, size), dtype = dtype)

    @property
    def unix_times(self):
        """ Return a NumPy array with all the Unix times"""
        return self._data[0]

    @unix_times.setter
    def unix_times(self, value):
        self._data[0] = value

    def set_time(self, index, unix_time):
        """ Set the Unix time of the index-th record """
        self.unix_times[index] = unix_time

    @property
    def magnitudes(self):
        """ Return a NumPy array with all the magnitudes"""
        return self._data[1]

    @magnitudes.setter
    def magnitudes(self, value):
        self._data[1] = value

    def set_mag(self, index, magnitude):
        """ Set the magnitude of the index-th record """
        self.magnitudes[index] = magnitude

    def __len__(self):
        """ Return how many Unix time - magnitude pairs there are """
        return self._data.shape[1] # number of columns

    @classmethod
    def from_curve(cls, light_curve):
        """ Instantiate and return a StringLengthFinder populated with the data
        (Unix times and magnitudes) of a LightCurve instance."""

        finder = cls(len(light_curve))
        for index, (unix_time, magnitude, snr) in enumerate(light_curve):
            finder.unix_times[index] = unix_time
            finder.magnitudes[index] = magnitude
        return finder

    def normalize(self):
        """ Scale the observations to give equal weight to measures and phases
        when periods are detected. Returns a new, normalized instance. For each
        magnitude: m'_{i} = (m_{i} - m_{min}) / 2(m_{max} - m_{min}) - 0.25

        """
        max_mag = self.magnitudes.max()
        min_mag = self.magnitudes.min()
        diff_mag = (max_mag - min_mag)
        normalized = copy.deepcopy(self)
        normalized.magnitudes = \
            (normalized.magnitudes - min_mag) / (2 * diff_mag) - 0.25
        return normalized

    def phase(self, period):
        """ Return the phase diagram as a new StringLengthFinder instance"""

        zero_t = self.unix_times.min()
        phased = copy.deepcopy(self)
        phased.unix_times = numpy.modf((phased.unix_times - zero_t) / period)[0]
        # Sort the phases into increasing order; i.e., find the indexes that
        # sort the phased Unix times and move the magnitudes correspondingly
        indexes = phased.unix_times.argsort()
        phased.unix_times = phased.unix_times[indexes]
        phased.magnitudes = phased.magnitudes[indexes]
        return phased

    def _point(self, index):
        """ Return the index-th Unix time and magnitude as a two-element NumPy
        array. Used to easily calculate the Euclidean distance in length() """
        return self._data[:,index]

    @property
    def length(self):
        """ Return the length of the line ('a piece of string') connecting each
        point with its nearest neighbor in phase, such as that returned by the
        phase method. Raises value error if the Unix times are not sorted."""

        if __debug__:
            for index in xrange(1, len(self)):
                if not self.unix_times[index - 1] <= self.unix_times[index]:
                    msg = "Unix times not sorted (is this a phase diagram?)"
                    raise ValueError(msg)

        length = 0
        for index in xrange(1, len(self)):
            # http://stackoverflow.com/q/1401712/
            length += numpy.linalg.norm(self._point(index) - self._point(index - 1))
        last = numpy.copy(self._point(-1))
        last[-1] += 1  # see formula at original paper
        length += numpy.linalg.norm(self._point(0) - last)
        return length

    @staticmethod
    def _best_period_range(finder, lower, upper, step):
        """ Find the best period in the range [lower, upper]

        The method evaluates all the periods using a step of 'step' seconds,
        returning that which minimizes the length of the line that connects
        each point with its nearest neighbor in phase.

        """
        best_length = float('inf')
        best_period = None

        for period in numpy.arange(lower, upper, step):
            length = finder.phase(period).length
            if length < best_length:
                best_length = length
                best_period = period
        return best_period

    def best_period(self, initial_step, exhaustive_step):
        """ Find by brute force the period of a star.

        Evaluate all the periods in the range ]0, n/2], with 'n' being the time
        elapsed between the first and last observations of the star. Return
        that which minimizes the length of the line that connects each point
        with its nearest neighbor in phase.

        The search is done in two parts: in the first one, a step of
        'initial_step' seconds is used to tentatively explore the search
        space. After the approximate location of the star period is determined,
        the search is intensified and periods are valuated again, centered at
        this preliminary result, plus and minus the value of the initial step,
        with a step of 'exhaustive_step' seconds. This is just a simple
        approach the diversification-intensification pattern.

        The reason why n/2 is the maximum detectable period is the
        Nyquist-Shannon sampling theorem, according to which at least two
        samples per period are necessary for detecting a certain frequency.

        """

        normalized = self.normalize()
        max_period = (self.unix_times.max() - self.unix_times.min()) / 2
        approx_period = \
          StringLengthFinder._best_period_range(normalized, exhaustive_step,
                                                max_period, initial_step)

        # And now the exhaustive search, centered at 'approx_period'. Note how
        # we guarantee, through the use of max(), that the lower bound of the
        # periods range is greater than or equal to the size of the step.
        exhaustive_lower = max(approx_period - initial_step, exhaustive_step)
        exhaustive_upper = approx_period + initial_step
        return StringLengthFinder._best_period_range(normalized, exhaustive_lower,
                                                     exhaustive_upper, exhaustive_step)


parser = optparse.OptionParser(description = description,
                               formatter = style.NewlinesFormatter())

parser.usage = "%prog [OPTION]... INPUT_DB"

parser.add_option('-o', action = 'store', type = 'str',
                  dest = 'output_db', default = None,
                  help = "path to the output LEMON database "
                  "[default: %default]")

parser.add_option('-w', action = 'store_true', dest = 'overwrite',
                  help = "overwrite output database if it already exists")

parser.add_option('--is', action = 'store', type = 'int',
                  dest = 'initial_step', default = 300,
                  help = "first step for the string-length method, in number "
                  "of seconds. This step will be used to determine the "
                  "approximate position of the period which minimizes the "
                  "sum of the lengths of line segments joining successive "
                  "points in a phase diagram [default: %default]")

parser.add_option('--es', action = 'store', type = 'int',
                  dest = 'exhaustive_step', default = 30,
                  help = "secondary step for the string-length method, in "
                  "number of seconds. A secondary, exhaustive search will "
                  "follow the initial one, exploring the search area along "
                  "the found minima in order to determine the exact value "
                  "of the star period [default: %default]")

parser.add_option('--cores', action = 'store', type = 'int',
                  dest = 'ncores', default = defaults.ncores,
                  help = defaults.desc['ncores'])

parser.add_option('-v', '--verbose', action = 'count',
                  dest = 'verbose', default = defaults.verbosity,
                  help = defaults.desc['verbosity'])

# The Queue is global -- this works, but note that we could have
# passed its reference to the function managed by pool.map_async.
# See http://stackoverflow.com/a/3217427/184363
queue = multiprocessing.Queue()

def parallel_periods(args):
    """ Method argument of map_async to compute periods in parallel.

    Functions defined in classes don't pickle, so we have moved this code here
    in order to be able to use it with multiprocessing's map_async. As it
    receives a single argument, values are passed in a tuple which is then
    unpacked.

    """

    star_id, light_curve, options = args

    if light_curve is None:
        logging.debug("Star %d: ignored (no light curve)" % star_id)
        queue.put((star_id, None))
        return

    if not len(light_curve) >=2:
        logging.debug("%Star %d: ignored (minimum of two points in "
                     "light curve not met)" % star_id)
        queue.put((star_id, None))
        return

    period_finder = StringLengthFinder.from_curve(light_curve)
    best_period = period_finder.best_period(options.initial_step,
                                            options.exhaustive_step)
    logging.debug("Star %d: period sucessfully determined (%.2f)" % \
                  (star_id, best_period))
    queue.put((star_id, best_period))


def main(arguments = None):
    """ main() function, encapsulated in a method to allow for easy invokation.

    This method follows Guido van Rossum's suggestions on how to write Python
    main() functions in order to make them more flexible. By encapsulating the
    main code of the script in a function and making it take an optional
    argument the script can be called not only from other modules, but also
    from the interactive Python prompt.

    Guido van van Rossum - Python main() functions:
    http://www.artima.com/weblogs/viewpost.jsp?thread=4829

    Keyword arguments:
    arguments - the list of command line arguments passed to the script.

    """

    if arguments is None:
        arguments = sys.argv[1:] # ignore argv[0], the script name
    (options, args) = parser.parse_args(args = arguments)

    # Adjust the logger level to WARNING, INFO or DEBUG, depending on the
    # given number of -v options (none, one or two or more, respectively)
    logging_level = logging.WARNING
    if options.verbose == 1:
        logging_level = logging.INFO
    elif options.verbose >= 2:
        logging_level = logging.DEBUG
    logging.basicConfig(format = style.LOG_FORMAT, level = logging_level)

    if len(args) != 1:
        parser.print_help()
        return 2 # used for command line syntax errors
    else:
        assert len(args) == 1
        input_db_path = args[0]

    if not os.path.exists(input_db_path):
        print "%sError. Database '%s' does not exist." % \
              (style.prefix, input_db_path)
        print style.error_exit_message
        return 1

    if not options.output_db:
        if not os.access(input_db_path, os.W_OK):
            print "%sError. Input database '%s' cannot be updated (no write " \
                  "access)" % (style.prefix, input_db_path)
            print style.error_exit_message
            return 1
        else:
            output_db_path = input_db_path
    else:
        output_db_path = options.output_db
        if os.path.exists(output_db_path):
            if not options.overwrite:
                print "%sError. The output database '%s' already exists." % \
                      (style.prefix, output_db_path)
                print style.error_exit_message
                return 1
            else:
                os.unlink(output_db_path)

        print "%sMaking a copy of the input database..." % style.prefix ,
        sys.stdout.flush()
        shutil.copy2(input_db_path, output_db_path)
        methods.owner_writable(output_db_path, True) # chmod u+w
        print 'done.'

    db = database.LEMONdB(output_db_path);
    nstars = len(db)
    print "%sThere are %d stars in the database" % (style.prefix, nstars)

    for pfilter in sorted(db.pfilters):

        print style.prefix
        print "%sPeriods for the %s filter will now be calculated." % \
              (style.prefix, pfilter)
        print "%sLoading light curves from the database..." % style.prefix ,
        sys.stdout.flush()
        all_light_curves = [db.get_light_curve(star_id, pfilter)
                            for star_id in db.star_ids]
        print 'done.'

        print "%sUsing the string-length method to determine the " \
              "periods..." % style.prefix
        pool = multiprocessing.Pool(options.ncores)
        map_async_args = ((star_id, light_curve, options)
                          for star_id, light_curve
                          in zip(db.star_ids, all_light_curves))
        result = pool.map_async(parallel_periods, map_async_args)
        #[parallel_periods(x) for x in map_async_args]

        methods.show_progress(0.0)
        while not result.ready():
            time.sleep(1)
            methods.show_progress(queue.qsize() / len(all_light_curves) * 100)
            # Do not update the progress bar when debugging; instead, print it
            # on a new line each time. This prevents the next logging message,
            # if any, from being printed on the same line that the bar.
            if logging_level < logging.WARNING:
                print

        result.get() # reraise exceptions of the remote call, if any
        methods.show_progress(100) # in case the queue was ready too soon
        print

        # The queue contains two-element tuples, (star_id, period)
        print "%sStoring the periods in the database..." % style.prefix
        methods.show_progress(0)
        periods = (queue.get() for x in xrange(queue.qsize()))
        for index, (star_id, star_period) in enumerate(periods):

            # NoneType is returned by parallel_periods when the period could
            # not be determined -- namely because there was no light curve!
            if star_period is None:
                msg = "Nothing for star %d; there was no light curve" % star_id
                logging.debug(msg)
                continue

            logging.debug("Storing period for star %d in database" % star_id)
            db.add_period(star_id, pfilter, star_period, options.exhaustive_step)
            logging.debug("Period for star %d successfully stored" % star_id)

            methods.show_progress(100 * (index + 1) / len(db))
            if logging_level < logging.WARNING:
                print

        else:
            logging.info("Periods for %s calculated" % pfilter)
            logging.debug("Committing database transaction")
            db.commit()
            logging.info("Database transaction commited")

            methods.show_progress(100.0)
            print

    # Update LEMONdB metadata
    db.date = time.time()
    db.author = os.getlogin()
    db.hostname = socket.gethostname()
    db.commit()

    methods.owner_writable(options.output_db, False) # chmod u-w
    print "%sYou're done ^_^" % style.prefix
    return 0

if __name__ == "__main__":
    sys.exit(main())

