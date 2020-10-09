#! /usr/bin/env python2
# encoding:UTF-8

# Copyright (c) 2012 Victor Terron. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of LEMON.
#
# LEMON is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division

import math
import operator
import numpy


def snr_to_error(snr):
    """Signal-to-noise ratio to error in magnitudes conversion.

    The method computes the error in magnitudes to which a determined SNR is
    equivalent. This is achieved by comparing the mean number of counts, c, to
    the maximum and minimum values, both of which are returned, induced by the
    noise, i.e, Δm = -2.5 log10((c +- c/SNR)/c) = -2.5 log10(1 ± 1/SNR). The
    two values derived from this formula are returned as a two-element tuple:
    the first one will correspond to the solution which the addition and the
    second one to the subtraction of the 'c/SNR' part of the equation.

    Note that this means that a S/N of 100 implies an observational error or
    0.01. This formula is taken from 'Astronomical Photometry', by Arne A.
    Henden and Ronald H. Kaitchuc [1982, page 79]

    The domain of the formula is the set of all non-negative numbers, as the
    base-10 logarithm of zero or a negative value cannot be calculated. This,
    ValueError is raised if the signal-to-noise ratio is not above one.

    """

    if not numpy.any(snr > 1):
        raise ValueError("SNR cannot be less than or equal to one")

    operators = (operator.add, operator.sub)
    return tuple([-2.5 * numpy.log10(f(1, 1 / snr)) for f in operators])


def error_to_snr(error):
    """Error in magnitudes to signal-to-noise ratio conversion.

    If we find the value of Δm for a given SNR, we can also obtain the
    signal-to-noise ration to which a determined error in magnitudes is
    equivalent. To do this we need to isolate in Henden's equation the
    value of SNR:

    1. Δm = -2.5 log10(1 ± 1/SNR)
    2. Δm / -2.5 = log10(1 ± 1/SNR)
    3. 10 ** (Δm / -2.5) = 10 ** (log10(1 ± 1/SNR))
    4. 10 ** (Δm / -2.5) = 1 ± 1/SNR
    5. 10 ** (Δm / -2.5) -1 = ± 1/SNR
    6. ± 1 / (10 ** (Δm / -2.5) -1) = SNR

    Note that the equation has a plus-minus sign, while the function only
    receives one error in magnitudes. How to handle this? Easily. Examining
    Henden's equation we can notice that when plus sign is used the noise is
    always negative, while a minus makes the returned value to be positive.
    Therefore, which sign the equation must use can be straight-forwardly
    determined: a plus if the error is negative and a minus if it turns out
    to be negative.

    """

    return (1 if error < 0 else -1) / (math.pow(10, error / -2.5) - 1)


def difference_error(*errors):
    """Return the absolute error of the difference of a series of errors.

    The error of a combination (addition / subtraction) of *independent* errors
    is given by their addition in quadrature. A pessimistic approach could
    consist in simply adding them, obtaining the maximum and minimum probable
    errors, but since they are independent quantities the first error is just
    as likely to (partially) cancel the second error as it is to increase it.

    The correct way of combining independent errors, therefore, is to add them
    in quadrature; in other words, the absolute error in a combination of sums
    and differences is obtained by taking the root of the sum of the squares of
    the absolute errors in the original quantities. For example, if we had two
    measures whose errors were 0.10 and 0.05, then the error of the addition or
    difference of these measures would be be sqrt(0.1 ** 2 + 0.05 ** 2) = 0.11.

    Note: the name of this error is so because, in our code, it is expected to
    be used only when calculating the error of the difference in magnitude
    between two or more stars. However, it may be perfectly used for additions
    too, and of course also for a combination of additions and subtractions.

    """

    return math.sqrt(sum(e ** 2 for e in errors))


def difference_snr(*snrs):
    """Return the SNR of the difference of a series of SNRs.

    The method returns the signal-to-noise ratio of the difference of a series
    of signal-to-noise ratios. This is internally done by converting the SNRs
    to errors in magnitudes, computing the absolute error and converting the
    resulting value back to its equivalent signal-to-noise ratio.

    As it is the case with difference_error, this method, despite its name, may
    also be used to compute the resulting signal-to-noise ratio of the addition
    of different SNRs, as well as the combination of additions and subtractions.

    """

    # We take the second of the two errors (i.e., the positive one) returned by
    # the conversion from signal-to-noise ratio. Using the negative value (the
    # minimum noise) would not work, as the addition is quadrature would make
    # the errors become positive, and as such they would be considered when
    # converting back to SNR.

    errors = [snr_to_error(s)[1] for s in snrs]
    assert all(e >= 0 for e in errors)
    error = difference_error(*errors)
    return error_to_snr(error)


def mean_error(errors, weights=None):
    """Return the absolute error of the arithmetic mean of errors.

    If we have a series of independent values (for example, instrumental
    magnitudes) and their corresponding uncertainty (the error in magnitudes),
    what is the uncertainty (error) if we calculate their arithmetic mean? By
    default all the values are assumed to have contributed equally to the mean,
    but the method can also calculate the error resulting from a weighted one.

    The arithmetic mean is just a scaled version of the sum, so the error
    scales as the quantity itself under scaling. Therefore, the error in the
    arithmetic mean of a series on n independent values is given by the
    equation \sqtr{e_{1}^{2} + e_{2}^{2} + ... + e_{n}^{2}} / n

    As for the weighted mean of n independent values with their corresponding
    weights, c_{1}, c_{2}... c{n}, and where \sum{i=0}^n{c_{i}} = 1, the error
    is \sqrt{c_{1}^2 e_{1}^2 + c_{2}^2 e_{2}^2 + ... + c_{n}^2 e_{n}^2). This
    is the equation implemented by the method, indeed, where the coefficients
    default to 1/n if no weights are given.

    Thanks so much to the people at Math Stack Exchange for their help:
    http://math.stackexchange.com/q/123276/

    Keyword arguments:
    weights - the coefficients of the weighted mean. The i-th weight is
              interpreted to correspond to the i-th error received by the
              method. The coefficients do not have to sum up to one; they will
              internally normalized to guarantee this. Therefore, the weights
              [0.5, 0.5], [1.0, 1.0] and [2.6, 2.6], e.g., are equivalent.
    """

    if weights is None:
        # All the values contribute equally (and weights sum up to one)
        weights = [1 / len(errors)] * len(errors)
    elif len(weights) != len(errors):
        raise ValueError("number of weights must equal that of errors")
    else:
        # Normalize the values so that they sum up to one
        weights = [w / math.fsum(weights) for w in weights]
    return math.sqrt(math.fsum(((w ** 2) * (e ** 2) for e, w in zip(errors, weights))))


def mean_snr(snrs, weights=None):
    """Return the SNR of the arithmetic mean of a series of SNSRs.

    The method returns the signal-to-noise ratio of the arithmetic or weighted
    mean of a series of signal-to-noise ratios. This is internally done by
    converting the SNRs to errors in magnitudes, computing the absolute error
    and converting the resulting value back to its equivalent SNR.

    Keyword arguments:
    weights - the coefficients of the weighted mean. The i-th weight is
              interpreted to correspond to the i-th error received by the
              method. The coefficients do not have to sum up to one; they will
              internally normalized to guarantee this. Therefore, the weights
              [0.5, 0.5], [1.0, 1.0] and [2.6, 2.6], e.g., are equivalent.

    """

    # We take the second of the two errors (i.e., the positive one) returned by
    # the conversion from signal-to-noise ratio. Using the negative value (the
    # minimum noise) would not work, as the addition is quadrature would make
    # the errors become positive, and as such they would be considered when
    # converting back to SNR.

    errors = [snr_to_error(s)[1] for s in snrs]
    assert all(e >= 0 for e in errors)
    error = mean_error(errors, weights=weights)
    return error_to_snr(error)


if __name__ == "__main__":
    pass
