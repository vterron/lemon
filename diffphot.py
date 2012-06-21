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
This module computes the light curve of each star in each photometric filter
using a moderately-modified version of Broeg's algorithm for computing an
optimal artificial comparison star. In short, what for each light curve we do
here is (a) to identify the most constant stars in the field, (b) compute their
weights, inversely proportional to its standard deviation and (c) generate the
light curve, subtracting, for each point in time, the magnitude of the star
from that of the comparison star. In this manner, positive values are returned
when the star is *brighter* than the comparison star, while *negative* means
that the star is fainter.

"""

import copy
import logging
import optparse
import os
import multiprocessing
import numpy
import random
import shutil
import sys
import time

# LEMON modules
import database
import methods
import snr
import style

class Weights(numpy.ndarray):
    """ Encapsulate the weights associated with some values """

    def __new__(cls, coefficients, dtype = numpy.float128):
        """ The constructor does not normalize the Weights, so they do not
        necessarily have to add up to 1.0. Use Weights.normalize() for that """
        if not len(coefficients):
            raise ValueError("arg is an empty sequence")
        return numpy.asarray(coefficients, dtype = dtype).view(cls)

    def __str__(self):
        coeffs_str = ", ".join(["%s" % x for x in self])
        return "%s(%s)" % (self.__class__.__name__, coeffs_str)

    def rescale(self, key):
        """ Exclude the key-th coefficient and return the rescaled Weights """
        if len(self) == 1:
            raise ValueError("cannot rescale one-element instance")
        return Weights(numpy.delete(self, key)).normalize()

    @property
    def total(self):
        return numpy.sum(self)

    def normalize(self):
        return Weights(self / self.total)

    @classmethod
    def inversely_proportional(cls, values, dtype = numpy.float128):
        """ Receives a series of values and returns the weights that are
        inversely proportional to them. Note that at least one value is
        required, and none of them can be zero (as we would be dividing by
        zero)"""

        if not len(values):
            raise ValueError("'values' is an empty sequence")
        if not all(values):
            raise ValueError("'values' cannot contain zeros")
        values = numpy.array(values,  dtype = dtype)
        return cls(1 / values).normalize()

    def absolute_percent_change(self, other, minimum = None):
        """ Return the percent change of two Weights of the same size. More
        specifically, the method returns the maximum absolute percent change
        between the coefficients of both instances that have the same
        position. In other words: the percent change between the first
        coefficient in 'self' and the first one in 'other' is calculated, and
        so on, and the maximum value among them is then returned. This method
        can be used in order to check how similar two weights are in percentage
        terms, or whether some Weights have converged to a given solution. Note
        that pairs which contain a zero in self are ignored (i.e., the maximum
        percent change between 0 and 4 cannot be calculated).

        The 'minimum' keyword, which must be a positive number, sets the lowest
        value that a coefficient must have (in both 'self' and 'other') in
        order to be taken into account. In this manner, extremely small
        coefficients, such as 1.09720056734e-15 (true story, we came across
        this exact value while reducing real data) may be excluded from the
        comparison, as at this point we are comparing mostly noise. ValueError
        is raised is 'minimum' causes all the coefficients to be discarded.

        """

        if minimum is not None and minimum <= 0:
            raise ValueError("'minimum' must be a positive number")

        if len(self) != len(other):
            raise ValueError("both instances must be of the same size")

        coefficients = []
        for x, y in zip(self, other):
            if x and (not minimum or min(x, y) >= minimum):
                coefficients.append((x, y))

        if not coefficients:
            msg = "all coefficients were discarded ('minimum' too high?)"
            raise ValueError(msg, self)
        else:
            return max((abs(methods.percentage_change(x, y))
                            for x, y in coefficients))

    @classmethod
    def random(cls, size):
        weights = cls([random.uniform(0, 1) for x in xrange(size)])
        return weights.normalize()


class StarSet(object):
    """ A collection of DBStars *with information for the exact same images*,
    all of them also taken in the same filter. Internally, and for performance
    purposes, all the photometric information is kept in a three-dimensional
    (star, data, image) array. """

    def _add(self, index, star):
        """ Add a DBStar to the set.

        In order not to maintain duplicate information in memory, the DBStar
        instances are parsed and their photometric information stored in the
        three-dimensional, low-level internal array.

        """

        # We cannot do time-series photometry with stars for which no have
        # information at all. It makes no sense, so don't even allow it.
        if not len(star):
            raise ValueError("star cannot be empty")

        if not len(self.star_ids):
            self._star_ids.append(star.id)
            self.pfilter = star.pfilter
            self._unix_times = star._unix_times

            # A cache, mapping each Unix time to its index in phot_info; passed
            # to the constructor of DBStar for O(1) lookups of Unix times
            self._times_indexes = {}
            for time_index, unix_time in enumerate(self._unix_times):
                self._times_indexes[unix_time] = time_index

            # Three-dimensional array: first dimension maps to the star; second
            # to the type of info (index 0 for mag, 1 for SNR) and third to the
            # Unix time (in case we need to know to which time a given index
            # correspond, we can use self._unix_times). For example, to get all
            # the magnitudes of the third star, we will do self.phot_info[2][0]
            self._phot_info[index] = star._phot_info[1:]

        else:
            if star.pfilter != self.pfilter:
                raise ValueError("star has filter '%s', expected '%s'" % \
                                 (star.pfilter, self.pfilter))

            if star.id in self.star_ids:
                raise ValueError("star with ID = %d already in the set" % star.id)

            for unix_time in self._unix_times:
                if unix_time not in star._time_indexes:
                    raise ValueError("stars must have info for the same Unix times")

            self._star_ids.append(star.id)
            self._phot_info[index] = star._phot_info[1:]

    def __init__(self, stars, dtype = numpy.float128):
        """ Instantiation method for the StarSet class.

        The 'stars' argument must be a sequence of DBStars, all of which will
        belong to the set must be given at instantiation time. In other words:
        you are expected not to call StarSet._add by yourself.

        """

        if not stars:
            raise ValueError("at least one star is needed")

        self.dtype = dtype
        self._star_ids = []
        # Do not set the array values to zero (marginally faster)
        self._phot_info = numpy.empty((len(stars), 2, len(stars[0])), dtype = self.dtype)
        for index, star in enumerate(stars):
            self._add(index, star)

    @property
    def star_ids(self):
        """ Return a list with the IDs of the stars.

        Note that this value is updated as the .phot_info attribute of each
        star is copied (see DBStar.add) at instantiation time; so you should
        use it with caution until the instantation of the class finishes

        """
        return self._star_ids

    def __len__(self):
        """ The number of stars in the set """
        return self._phot_info.shape[0]  # 1st dimension = stars

    @property
    def nimages(self):
        """ The number of Unix times for which the stars have info """
        return len(self._unix_times)

    def __delitem__(self, index):
        """ Delete the index-th star """
        self._phot_info = numpy.delete(self._phot_info, index, axis = 0)
        assert self._phot_info.shape[0] == len(self)  # first dimension: star
        del self._star_ids[index]
        assert len(self.star_ids) == len(self), "%d vs %d" % (len(self.star_ids), len(self))

    def __getitem__(self, index):
        """ Return the index-th star as a DBStar instance """

        # The constructor of a DBStar receives a NumPy array with three rows
        # (the first for the time, the second for the magnitude and the last
        # for the SNR), and as many columns as records for which there is
        # photometric information.
        sphot_info = numpy.empty((3, self.nimages), dtype = self.dtype)
        sphot_info[0] = self._unix_times
        sphot_info[1] = self._phot_info[index][0]
        sphot_info[2] = self._phot_info[index][1]

        id_ = self._star_ids[index]
        return database.DBStar(id_, self.pfilter, sphot_info,
                               self._times_indexes, self.dtype)

    def magnitude_inversely_proportional_weights(self):
        """ Returns the Weights which are proportional to the median of the
        normalized flux of each star. The logic behind this is that the
        brightest star in the field, for example, should always be the
        brightest, no matter what the atmospheric conditions are. This only
        applies, of course, for non-variable stars, which are assumed to be the
        great majority of the objects in the field. This approach means that
        all the instrumental magnitudes of all the stars are taken into account
        in order to calculate the weights, which is more robust than relying on
        a single image, as we did before, even if it was that with the best
        seeing.

        Note that we work with fluxes: if the instrumental magnitude of A is
        one magnitude greater than that of B, it means that A is 2.512x fainter
        than B, so its weight must be 2.512 times smaller than that of B."""

        # bi-dimensional matrix (as many rows as stars, as many columns as
        # images) in which the normalized magnitude of the stars are stored
        norm_mags = numpy.empty((len(self), self.nimages), dtype = self.dtype)

        # For each image (third dimension of self.phot_info), normalize the
        # magnitudes (second dimension, x = 0) of each star (firstdimension).
        # That is: the magnitudes for all the stars in the image are divided
        # by the maximum magnitude.
        for image_index in xrange(self.nimages):
            image_mags = self._phot_info[:, 0, image_index]
            norm_mags[:,image_index] = image_mags / image_mags.max()

        # Now, for each star, calculate the median of the normalized magnitudes
        mag_medians = numpy.empty(len(self), dtype = self.dtype)
        for star_index in xrange(len(self)):
            mag_medians[star_index] = numpy.median(norm_mags[star_index])

        pogsonr = 100 ** 0.2  # fifth root of 100 (Pogson's Ratio)
        return Weights.inversely_proportional(pogsonr ** mag_medians)


    def light_curve(self, weights, star, no_snr = False, _exclude_index = None):
        """ Compute the light curve of 'star' using the stars in the set as the
        comparison star and these weights. For each point in time, subtracts
        the magnitude of the comparison star from that of 'star'.

        Note that it is mandatory that 'star' has photometric information for
        the same Unix times for which the stars in 'self' have info.

        The calculation of the signal-to-noise ratio of each differential
        magnitude is currently, being written in pure Python, extremely
        slow. It is barely noticeable when a few conversions from SNR to error
        in magnitudes (or vice versa) are done, but when invoked thousands of
        times it causes an unbearable performance hit. Until the 'snr' module
        is rewritten, 'no_snr', if set to True, will make the method skip the
        calculation of the resulting SNRs, using None instead. This is probably
        only needed by StarSet.broeg_weights(), which does not care about the
        signal-to-noise rations, but only the standard deviation of the light
        curves, and can therefore skip their computation.

        If given, '_exclude_index' determines the index of the star in the set
        that will not be used as comparison star. Note that, if that is the
        case, there must be as many weights as the number of stars in the set
        minus one. This parameter is here (and probably only needed because of)
        the StarSet.broeg_weights method, which needs to compute the light
        curve of each star in the set using as comparison all the others.

        """

        if len(weights) != len(self) - (1 if _exclude_index is not None else 0):
            msg = "number of weights must match that of comparison stars"
            raise ValueError(msg)

        if _exclude_index is not None and \
            (_exclude_index < 0 or _exclude_index >= len(self)):
              raise IndexError("_exclude_index out of range")

        if not numpy.all(numpy.equal(star._unix_times, self._unix_times)):
            raise ValueError("the star must have photometric records for "
                             "the same Unix times for which StarSet has")

        cstars = self.star_ids[:]
        if _exclude_index is not None:
            del cstars[_exclude_index]

        len(weights) == len(weights)
        curve = database.LightCurve(self.pfilter, cstars, weights,
                                    dtype = self.dtype)

        # Instead of, for each image, taking the list of magnitudes and
        # removing that of the star to be ignored, it is simpler and
        # faster to assign to it a weight of zero and proceed as usual.
        if _exclude_index is not None:
            weights = numpy.insert(weights, _exclude_index, 0.0)

        assert len(weights) == len(self)
        for index, unix_time in enumerate(self._unix_times):
            # The magnitude and SNR of the comparison star
            cmags = self._phot_info[:, 0, index]
            csnrs = self._phot_info[:, 1, index]
            cmag = numpy.average(cmags, weights = weights)
            csnr = None if no_snr else snr.mean_snr(csnrs, weights = weights)

            smag  = star.mag(index)
            ssnr  = star.snr(index)

            dmag = smag - cmag
            dsnr = None if no_snr else snr.difference_snr(ssnr, csnr)
            curve.add(unix_time, dmag, dsnr)

        return curve

    def broeg_weights(self, pct = 0.01, max_iters = None, minimum = None):
        """ Determine the weights that give the optimum comparison star.

        This is our moderately modified implementation of Broeg's algorithm,
        used to identify a proper artificial comparison star by means of
        determining the variability of each star and giving to it a weight
        inversely proportional to its standard deviation. This is just a
        cursory explanation; a detailed documentation of these methods will
        follow at some point in the future.

        ValueError is raised if there are less than two stars in the set, as
        we cannot calculate the standard deviation of a single differential
        magnitude (and assuming it to be zero, as NumPy does, is equally
        useless as all the stars would have the same standard deviation)

        Keyword arguments:
        - pct: the convergence threshold, given as a percentage change. The
               optimum weights will be considered to have been found when the
               percentage change (see Weights.absolute_percent_change) between
               the weights computed in the last two iterations is less than or
               equal to this value.
        - max_iters: the maximum number of iterations of the algorithm. Any
                     value which evaluates to False (such as zero or None) will
                     be interpreted as to mean 'use the current value of the
                     recursion limit' - as returned by sys.getrecursionlimit().
        - minimum: the minimum value for a coefficient to be taken into account
                   when calculating the percentage change between two Weights;
                   used in order to prevent scientifically-insignificant values
                   from making the algorithm stop or iterate more than needed.

        """

        if not len(self):
            raise ValueError("cannot work with an empty instance")

        # Special case: only one star, then its weight cannot be other than one
        if len(self) == 1:
            return Weights([1.0])

        # Special case: two stars, both will have the same standard deviation
        # (as their light curves will be generated by comparing each one to the
        # other), so they will have the same weight (0.5)
        if len(self) == 2:
            return Weights([0.5, 0.5])

        # At least two images are needed for a light curve. We cannot calculate
        # the standard deviation of a single differential magnitude.
        if self.nimages < 2:
            raise ValueError("at least two images are needed")

        # Initial weights are inversely proportional to the magnitude of each
        # star. Then, the differential magnitude of each comparison star is
        # calculated with respect to the others using the same (rescaled)
        # weights. The standard deviation of the light curve that results from
        # using these (rescaled) weights is used in order to compute the new
        # weights for each star. We stop when the absolute percent change
        # between the old weights and the new one is below the threshold

        weights = [self.magnitude_inversely_proportional_weights()]
        curves_stdevs = numpy.empty(len(self), dtype = self.dtype)
        for iteration in xrange(max_iters or sys.getrecursionlimit()):
            for star_index in xrange(len(self)):
                rweights = weights[-1].rescale(star_index)
                star_curve = self.light_curve(rweights, self[star_index],
                                              _exclude_index = star_index,
                                              no_snr = True)
                curves_stdevs[star_index] = star_curve.stdev

            # Added Thu 01 Mar 2012 09:24:45 PM CET. If somehow a star ends up
            # having all the weight (I thought this was impossible), avoid the
            # division by zero -- Weights.inversely_proportional would raise
            # ValueError. Instead, stop the computation so that the last valid
            # weights are returned.
            #
            # How this may happen? The only time it's ever happened to us we
            # were computing the weights of a set of stars with had information
            # only in two Unix dates. With that few points in the light curve
            # of each star we may have, if only by chance, two differential
            # magnitudes almost identical and therefore a standard deviation
            # which would converge to zero.

            if not all(curves_stdevs):
                break

            weights.append(Weights.inversely_proportional(curves_stdevs))
            if weights[-2].absolute_percent_change(weights[-1], minimum = minimum) < pct:
                break

        return weights[-1]

    def worst(self, fraction, pct = 0.01, max_iters = None, minimum = None):
        """ Return the indexes of the 'fraction' (range (0, 1]) of the less
        constant stars among those in the set. This is done by taking the Broeg
        weights and identifying the n stars with the lowest ones (identified
        thus as those with the less constant light curves). n is obtained by
        multiplying 'fraction' by the number of stars in the set, rounded to
        the closest integer and having always a minimum value of one. The order
        of the indexes is not casual: instead, they are order by the
        'unconstantness' of the stars -- the index of the most variable star is
        returned first, after it that of the second most variable star, and so
        on.

        Keyword arguments:
        - pct: the convergence threshold, given as a percentage change. The
               optimum weights will be considered to have been found when the
               percentage change (see Weights.absolute_percent_change) between
               the weights computed in the last two iterations is less than or
               equal to this value.
        - max_iters: the maximum number of iterations of the algorithm. Any
                     value which evaluates to False (such as zero or None) will
                     be interpreted as to mean 'use the current value of the
                     recursion limit' - as returned by sys.getrecursionlimit().
        - minimum: the minimum value for a coefficient to be taken into account
                   when calculating the percentage change between two Weights;
                   used in order to prevent scientifically-insignificant values
                   from making the algorithm stop or iterate more than needed.

        """

        if not 0 < fraction <= 1:
            raise ValueError("'fraction' must be in the range (0,1]")

        if len(self) < 3:
            raise ValueError("at least three stars are needed")

        # The number of stars with the lowest weights (and therefore the
        # highest standard deviation in their light curves) to be returned
        nstars = int(round(fraction * len(self)))
        if not nstars:
            nstars = 1

        # Get the indexes of the n minimum values among the Weights, as these
        # will correspond to the stars with the maximum standard deviation - the
        # 'worst'. Code courtesy of 'aix' at: http://stackoverflow.com/q/6910672
        bweights = self.broeg_weights(pct = pct, max_iters = max_iters, minimum = minimum)
        return list(bweights.argsort()[:nstars])

    def best(self, n, fraction = 0.1, pct = 0.01, max_iters = None, minimum = None):
        """ Return the n most constant stars in the set, found by iteratively
        identifying the 'fraction' times the len(self) less constant stars,
        discarding them and recomputing all the light curves until only n stars
        are left.  Ideally stars would be discarded one by one until only n
        remained, but this would be terribly CPU-expensive for large data sets,
        so we have to consider and discard several of them at once.

        ValueError is raised if the method is invoked on a StarSet with fewer
        than three stars, as that is the minimum number of stars needed in
        order to be able to determine their variability reliably, or if the
        number of stars in the set is smaller than the number of best stars
        requested.

        Keyword arguments:
        - pct: the convergence threshold, given as a percentage change. The
               optimum weights will be considered to have been found when the
               percentage change (see Weights.absolute_percent_change) between
               the weights computed in the last two iterations is less than or
               equal to this value.
        - max_iters: the maximum number of iterations of the algorithm. Any
                     value which evaluates to False (such as zero or None) will
                     be interpreted as to mean 'use the current value of the
                     recursion limit' - as returned by sys.getrecursionlimit().
        - minimum: the minimum value for a coefficient to be taken into account
                   when calculating the percentage change between two Weights;
                   used in order to prevent scientifically-insignificant values
                   from making the algorithm stop or iterate more than needed.

        """

        if not 0 < fraction <= 1:
            raise ValueError("'fraction' must be in the range (0,1]")

        if len(self) < 3:
            raise ValueError("at least three stars are needed")

        if len(self) < n:
            msg = "'n' is greater than the number of stars in the set"
            raise ValueError(msg)

        stars = copy.deepcopy(self)

        # We do not discard stars here until only 'n' stars are left, as at
        # least three stars are needed in order to determine their variability
        # reliably. This means that, once we are left with fewer than tree
        # stars (that is, two), we cannot iterate more. The solution, thus, is
        # to discard stars until n or three stars are left, whatever happens
        # first. After that, in we stopped at three when n was lower, we can
        # discard the last batch of stars until only n are left.
        while len(stars) > max(n, 3):
            worst_indexes = \
                stars.worst(fraction, pct = pct,
                            max_iters = max_iters, minimum = minimum)

            # The stars whose indexes have been returned cannot be blindly
            # deleted, as the difference between the number of indexes returned
            # and the number of stars left in 'stars' may be higher than the
            # number of stars that have to be deleted in order to get 'n' stars
            # left. For example, assume there are 100 stars and we want the
            # best 30, with a fraction of 0.5. 50 stars would be deleted in the
            # first loop, while 25 more would be removed in the second.  There
            # would then be 100-50-25 = 25 stars left, when we wanted 30!

            del worst_indexes[(len(stars) - max(n, 3)):]
            # Remove the stars going backwards, from highest to lowest index
            for index in sorted(worst_indexes, reverse = True):
                del stars[index]

        # If there are only three stars left but there are still stars to
        # discard we must identify them all at once, independently of the value
        # of 'fraction'. The reason for this is that a minimum of three stars
        # in needed to realibly determine their variability.

        assert len(stars) >= n
        if len(stars) != n:
            worst_indexes = stars.worst(1.0, pct = pct, max_iters = max_iters)
            del worst_indexes[(len(stars) - n):]
            for index in sorted(worst_indexes, reverse = True):
                del stars[index]

        assert len(stars) == n
        return stars

# The Queue is global -- this works, but note that we could have
# passed its reference to the function managed by pool.map_async.
# See http://stackoverflow.com/a/3217427/184363
queue = multiprocessing.Queue()

def parallel_light_curves(args):
    """ Method argument of map_async to compute light curves in parallel.

    Functions defined in classes don't pickle, so we have moved this code here
    in order to be able to use it with multiprocessing's map_async. As it
    receives a single argument, values are passed in a tuple which is then
    unpacked.

    """

    star, all_stars, options = args
    logging.debug("Star %d: photometry on %d images, enforced minimum of %d" %
                 (star.id, len(star), options.min_images))

    if len(star) < options.min_images:
        logging.debug("Star %d: ignored (minimum of %d images not met)" %
                     (star.id, options.min_images))
        queue.put((star.id, None))
        return

    complete_for = star.complete_for(all_stars)
    logging.debug("Star %d: %d complete stars, enforced minimum = %d" %
                 (star.id, len(complete_for), options.min_cstars))

    ncstars = min(len(complete_for), options.ncstars)
    if ncstars < options.min_cstars:
        logging.debug("Star %d: ignored (minimum of %d comparison stars "
                      "not met)" % (star.id, options.min_cstars))
        queue.put((star.id, None))
        return

    logging.debug("Star %d: will use %d complete stars (out of %d) as "
                  "comparison)" % (star.id, ncstars, len(complete_for)))

    logging.debug("Star %d: fraction of worst stars to discard = %.4f" %
                 (star.id, options.worst_fraction))
    logging.debug("Star %d: weights percentage change threshold = %.4f" %
                 (star.id, options.pct))
    logging.debug("Star %d: minimum Weights coefficients = %.4f" %
                 (star.id, options.wminimum))
    logging.debug("Star %d: maximum Broeg iterations: %.4f" %
                 (star.id, options.max_iters))

    complete_stars = StarSet(complete_for)
    comparison_stars = \
        complete_stars.best(ncstars, fraction = options.worst_fraction,
                            pct = options.pct, minimum = options.wminimum,
                            max_iters = options.max_iters)

    logging.debug("Star %d: identified the %d best stars (out of %d)" %
                 (star.id, len(comparison_stars), len(complete_stars)))
    logging.debug("Star %d: best stars IDs: %s" %
                 (star.id, [x.id for x in comparison_stars]))

    cweights = \
        comparison_stars.broeg_weights(pct = options.pct,
                                       minimum = options.wminimum,
                                       max_iters = options.max_iters)

    logging.debug("Star %d: Broeg weights: %s" % (star.id, str(cweights)))
    light_curve = comparison_stars.light_curve(cweights, star)
    logging.debug("Star %d: light curve sucessfully generated "
                  "(stdev = %.4f)" % (star.id, light_curve.stdev))
    queue.put((star.id, light_curve))


parser = optparse.OptionParser(description = description,
                               formatter = style.NewlinesFormatter())

parser.usage = "%prog [OPTION]... INPUT_DB"
parser.add_option('-o', action = 'store', type = 'str',
                  dest = 'output_db', default = 'diffphotometry.db',
                  help = "path to the output LEMON database "
                  "[default: %default]")

parser.add_option('-w', action = 'store_true', dest = 'overwrite',
                  help = "overwrite output database if it already exists")

parser.add_option('--cores', action = 'store', type = 'int',
                  dest = 'ncores', default = multiprocessing.cpu_count(),
                  help = "the maximum number of cores available to the "
                  "module. This option defaults to the number of CPUs in "
                  "the system, which are automatically detected "
                  "[default: %default]")

parser.add_option('-v', '--verbose', action = 'count',
                  dest = 'verbose', default = 0,
                  help = "increase the amount of information given during "
                  "the execution. A single -v tracks INFO events, while two "
                  "or more enable DEBUG messages. These should only be used "
                  "when debugging the module.")

curves_group = optparse.OptionGroup(parser, "Light Curves", "")
curves_group.add_option('--mimages', action = 'store', type = 'int',
                        dest = 'min_images', default = 10,
                        help = "the minimum number of images in which a star "
                        "must have been observed; the light curve will not be "
                        "calculated for those filters for which the star was "
                        "observed a number of times lower than this value "
                        "[default: %default]")

curves_group.add_option('-n', action = 'store', type = 'int',
                        dest = 'ncstars', default = 20,
                        help = "number of complete stars that will be used as "
                        "the artificial comparison star. For each star, its "
                        "'complete' set are those stars for which there is "
                        "photometric information in at least each of the "
                        "images in which it was observed [default: %default]")

curves_group.add_option('--mincstars', action = 'store', type = 'int',
                        dest = 'min_cstars', default = 8,
                        help = "the minimum number of stars used to compute "
                        "the artificial comparison star, regarless of the "
                        "value of the -n option. The light curve will not be "
                        "generated if this minimum value cannot be reached. "
                        "It follows that if this option is set to a value "
                        "greater than that of -n, no curve is generated, so "
                        "the module exits with an error. Although acceptable, "
                        "a value equal to that of -n, is not recommended. "
                        "[default: %default]")
parser.add_option_group(curves_group)

broeg_group = optparse.OptionGroup(parser, "Broeg's Algorithm", "")
broeg_group.add_option('--pct', action = 'store', type = 'float',
                       dest = 'pct', default = 0.01,
                       help = "the convergence threshold of the weights "
                       "determined by Broeg algorithm, given as a percentage "
                       "change. Iteration will stop when the percentage "
                       "change between the last two weights is less than or "
                       "equal to this value [default: %default]")

broeg_group.add_option('--wminimum', action = 'store', type = 'float',
                       dest = 'wminimum', default = 0.0001,
                       help = "the minimum value for a coefficient to be "
                       "taken into account when calculating the percentage "
                       "change between two Weights; needed to prevent "
                       "scientifically-insignificant values from making the "
                       "algorithm stop earlier or iterate more than needed "
                       "[default: %default]")

broeg_group.add_option('--max-iters', action = 'store', type = 'int',
                       dest = 'max_iters', default = sys.getrecursionlimit(),
                       help = "the maximum number of iterations of the "
                       "algorithm. If exceeded, the last computed weights "
                       "will be taken, regardless of the percentage change. "
                       "This option defaults to the maximum depth of the "
                       "Python interpreter stack [default: %default]")
parser.add_option_group(broeg_group)

best_group = optparse.OptionGroup(parser, "Worst and Best Stars", "")
best_group.add_option('--wfraction', action = 'store', type = 'float',
                      dest = 'worst_fraction', default = 0.75,
                      help = "the fraction of the stars that will be "
                      "discarded at each step when identifying which are the "
                      "most constant stars. The lower this value, the most "
                      "reliable the identification of the constant stars in "
                      "the field will be, but also more CPU-expensive "
                      "[default: %default]")
parser.add_option_group(best_group)

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
        return 2  # used for command line syntax errors
    else:
        assert len(args) == 1
        input_db_path = args[0]

    if options.min_cstars > options.ncstars:
        print "%sError. The value of --mincstars must be <= -n." % style.prefix
        print style.error_exit_message
        return 1

    if not os.path.exists(input_db_path):
        print "%sError. Database '%s' does not exist." % (style.prefix, input_db_path)
        print style.error_exit_message
        return 1

    if os.path.exists(options.output_db):
        if not options.overwrite:
            print "%sError. The output database '%s' already exists." % \
                  (style.prefix, options.output_db)
            print style.error_exit_message
            return 1
        else:
            os.unlink(options.output_db)

    # Note that we do not allow the user to update an existing LEMON database,
    # although it would make sense if we had just done photometry and now we
    # wanted to generate the light curves. However, and to be on the safe side,
    # we prefer to preserve a copy of the original database, as it is not
    # inconceivable that the astronomer may need to recompute the curves more
    # than once, each time with a different set of parameters.

    print "%sMaking a copy of the input database..." % style.prefix ,
    sys.stdout.flush()
    shutil.copy2(input_db_path, options.output_db)
    methods.owner_writable(options.output_db, True) # chmod u+w
    print 'done.'

    db = database.LEMONdB(options.output_db);
    nstars = len(db)
    print "%sThere are %d stars in the database" % (style.prefix, nstars)

    for pfilter in sorted(db.pfilters):

        print style.prefix
        print "%sLight curves for the %s filter will now be generated." % \
              (style.prefix, pfilter)
        print "%sLoading photometric information..." % style.prefix ,
        sys.stdout.flush()
        all_stars = [db.get_photometry(star_id, pfilter) for star_id in db.star_ids]
        print 'done.'

        # The generation of each light curve is a task independent from the
        # others, so we can use a pool of workers and do it in parallel.
        pool = multiprocessing.Pool(options.ncores)
        map_async_args = ((star, all_stars, options) for star in all_stars)
        result = pool.map_async(parallel_light_curves, map_async_args)

        methods.show_progress(0.0)
        while not result.ready():
            time.sleep(1)
            methods.show_progress(queue.qsize() / len(all_stars) * 100)

        result.get() # reraise exceptions of the remote call, if any
        methods.show_progress(100) # in case the queue was ready too soon
        print

        # The multiprocessing queue contains two-element tuples,
        # mapping the ID of each star to its light curve.
        print "%sStoring the light curves in the database..." % style.prefix
        methods.show_progress(0)
        light_curves = (queue.get() for x in xrange(queue.qsize()))
        for index, (star_id, curve) in enumerate(light_curves):

            # NoneType is returned by parallel_light_curves when the light
            # curve could not be calculated -- because it did not meet the
            # minimum number of images or comparison stars.
            if curve is None:
                logging.debug("Nothing for star %d; light curve could not "
                              "be generated" % star_id)
                continue

            logging.debug("Storing light curve for star %d in database" % star_id)
            db.add_light_curve(star_id, curve)
            logging.debug("Light curve for star %d successfully stored" % star_id)

            methods.show_progress(100 * (index + 1) / len(all_stars))

        else:
            logging.info("Light curves for %s generated" % pfilter)
            logging.debug("Committing database transaction")
            db.commit()
            logging.info("Database transaction commited")

            methods.show_progress(100.0)
            print

    methods.owner_writable(options.output_db, False) # chmod u-w
    print "%sYou're done ^_^" % style.prefix
    return 0

if __name__ == "__main__":
    sys.exit(main())

