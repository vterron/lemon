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
We aim high here, and hopefully will not fail too miserably, as this module
attempts to automatically determine which are the optimal aperture and sky
annuli for aperture photometry in the context of time-series analysis.

This is possible because of the premise that, the better the photometric
parameters are, the most stable the light curve (i.e., the lowest its standard
deviation) of the most constant stars in the field will be. In this manner, how
good the parameters are can be easily measured, as the aforementioned standard
deviation will get lower as we get closer to the optimal aperture and sky
annuli.

Ideally, we would compare two candidate photometric parameters by means of
evaluating the light curves of all the stars, but this would be rather
impractical for large data sets (hundreds of images, thousands of stars, or
even larger). Therefore, what the module does is to initially identify which
are the most constant objects in the field so that only they have photometry
done and their light curves is compared, thus increasing performance many-fold.

The dimensions of the search space are given in terms of the median FWHM of the
images in each photometric filter, using this value a number of times in order
to establish both the lower and upper bounds for the apertures that have to be
evaluated. The search space is explored in a series of steps of a fixed number
of pixels. The sky annulus, although also expressed in FWHMs, remains, however,
the same for all the candidate apertures that are evaluated.

As it is the case with photometry.py, this module receives two parameters, the
first being the reference image, on which sources are detected, and the second
a LEMON XML file, produced by offsets.py, which lists the translation offsets
of all the images with respect to the aforementioned reference image.

The output of the module is a LEMON XML file which lists all the photometric
parameters that, for each photometric filter, were tested during the
exploration of the search space. This file is later parsed by the photometry.py
module in order to automatically use these optimal parameters.

"""

import collections
import logging
import numpy
import operator
import optparse
import os
import style
import sys
import tempfile

# LEMON modules
import customparser
import diffphot
import keywords
import mining
import photometry
import subprocess
import xmlparse

class NotEnoughImages(ValueError):
    pass

class NotEnoughConstantStars(ValueError):
    pass


parser = customparser.get_parser(description)
parser.usage = "%prog [OPTION]... OFFSETS_XML_FILE"
parser.add_option('--output', action = 'store', type = 'str',
                  dest = 'xml_output', default = 'annuli.xml',
                  help = "path of the XML file to which the evaluated "
                  "photometric parameters will be saved. This file can be "
                  "later parsed by photometry.py in order to use the "
                  "optimal parameters for each photometric band "
                  "[default: %default]")

parser.add_option('--overwrite', action = 'store_true', dest = 'overwrite',
                  help = "overwrite output XML file if it already exists")

parser.add_option(photometry.parser.get_option('--margin'))
parser.add_option(photometry.parser.get_option('--gain'))
parser.add_option(photometry.parser.get_option('--cores'))
parser.add_option(photometry.parser.get_option('--verbose'))

qphot_group = optparse.OptionGroup(parser, "Initial Photometry",
              "In what may seem sort of a recursive problem, we need to do an "
              "initial (aperture) photometry in order to detect all the stars "
              "in the field and determine which among them are most constant. "
              "The value of the photometry parameters (aperture and sky "
              "annuli) are defined in terms of the *median* FWHM of the "
              "images in each band.\n\n"

              "And yes, we are aware that all this means that the search for "
              "the optimal photometric aperture parameters is dependent upon "
              "an initial photometry whose parameters have to be, even in "
              "terms of the FWHM, initially specified. Kind of paradoxical, "
              "we know.")

qphot_group.add_option(photometry.parser.get_option('--aperture'))
qphot_group.add_option(photometry.parser.get_option('--annulus'))
qphot_group.add_option(photometry.parser.get_option('--dannulus'))

qphot_group.add_option('--min-sky', action = 'store',
                       type = 'float', dest = 'min',
                       default = photometry.parser.defaults['min'],
                       help ="the minimum width of the sky annulus, in "
                       "pixels, regardless of the value derived from the "
                       "FWHM. This option is intended to prevent small FWHMs "
                       "from resulting in too thin an sky annulus, and "
                       "applies to both the initial photometry and "
                       "subsequent explorations of the search space "
                       "[default = %default]")

parser.add_option_group(qphot_group)

stats_group = optparse.OptionGroup(parser, "Comparison stars",
              "Ideally, two candidate apertures would be compared by means "
              "of evaluating the light curves of all the stars: the premise "
              "here is that, the better the photometric parameters happen to "
              "be, the most stable the light curve of the constant stars "
              "will be. However, doing photometry and generating the light "
              "curves for each star and aperture would be extremely "
              "impractical for large data sets, as it is quite often our "
              "case. Therefore, our approach is to compare the light curves "
              "of only a subset of the stars: those that are the most "
              "constant. Furthermore, and in order to allow for rare "
              "anomalies in the light curve of some stars, not all of these "
              "light curves are used; instead, they are evaluated as a whole "
              "by taking the median of the standard deviation of the most "
              "stable among them. Note that this means that the stars whose "
              "light curves are used to compute the median are not "
              "necessarily always the same.")

stats_group.add_option('--constant', action = 'store', type = 'float',
                       dest = 'nconstant', default = 20,
                       help = "the number of stars used in order to "
                       "compare each set of photometric parameters. For "
                       "each of these stars, its light curve will be "
                       "generated using all the other as comparison "
                       "(-n option at diffphot.py) [default: %default]")

stats_group.add_option('--minimum-constant', action = 'store', type = int,
                       dest = 'pminimum', default = 5,
                       help = "the minimum number of constant stars which "
                       "must have had their light curves for the percentile "
                       "of the standard deviations to be calculated. This "
                       "option is intended to prevent a nonrepresentative "
                       "value (for example, a single light curve would "
                       "result in a standard deviation of zero) from being "
                       "returned if most light curves cannot be calculated "
                       "with a certain set of photometric paramaters. If "
                       "less than this number of light curves are computed, "
                       "the parameters are ignored. [default: %default]")
parser.add_option_group(stats_group)

search_group = optparse.OptionGroup(parser, "Search space",
               "The number of apertures that will be evaluated is determined "
               "by the median FWHM of the images in each photometric filter, "
               "from which the lower and upper bounds of the candidate "
               "apertures are derived. The sky annulus, however, remains the "
               "same for all the candidate apertures that are evaluated.")

search_group.add_option('--lower', action = 'store', type = 'float',
                        dest = 'lower', default = 0.5,
                        help = "the lower bound of the range of aperture "
                        "annuli that will be evaluated, in number of times "
                        "the median FWHM [default = %default]")

search_group.add_option('--upper', action = 'store', type = 'float',
                        dest = 'upper', default = 4.5,
                        help = "the upper bound of the range of aperture "
                        "annuli that will be evaluated, in number of times "
                        "the median FWHM [default = %default]")

search_group.add_option('--step', action = 'store', type = 'float',
                        dest = 'step', default = 1,
                        help = "the number of pixels by which the candidate "
                        "apertures in the range [lower x FWHM, upper x FWHM] "
                        "will be incremented each time [default = %default]")

search_group.add_option('--sky', action = 'store', type = 'float',
                        dest = 'sky', default = 4.6,
                        help = "the inner radius of the sky annulus, in "
                        "number of times the median FWHM. Note that the "
                        "sky annulus remains the same for all the "
                        "evaluated apertures. [default = %default]")

search_group.add_option('--width', action = 'store', type = 'float',
                        dest = 'width', default = 1.0,
                        help = "the width of the sky annulus, in number of "
                        "times each candidate aperture. Note that the sky "
                        "annulus remains the same for all the evaluated "
                        "apertures. [default = %default]")
parser.add_option_group(search_group)

const_group = optparse.OptionGroup(parser, "Stars eligibility",
              "Apart from a stable light curve, a star must not have "
              "saturated in any of the images if it is to be eligible as a "
              "constant star. Note that this means that the candidates to "
              "constant stars are those that are below the saturation level "
              "for all the images")
const_group.add_option(photometry.parser.get_option('--maximum'))
parser.add_option_group(const_group)

diffphot_group = optparse.OptionGroup(parser, "Differential Photometry",
                 "These options, the same that can be found in the "
                 "diffphot.py, determine how light curves are generated.")

diffphot_group.add_option(diffphot.parser.get_option('--minimum-images'))
diffphot_group.add_option(diffphot.parser.get_option('--minimum-stars'))
diffphot_group.add_option(diffphot.parser.get_option('--pct'))
diffphot_group.add_option(diffphot.parser.get_option('--weights-threshold'))
diffphot_group.add_option(diffphot.parser.get_option('--max-iters'))
diffphot_group.add_option(diffphot.parser.get_option('--worst-fraction'))
parser.add_option_group(diffphot_group)

key_group = optparse.OptionGroup(parser, "FITS Keywords",
                                 keywords.group_description)

key_group.add_option(photometry.parser.get_option('--expk'))
key_group.add_option(photometry.parser.get_option('--coaddk'))
key_group.add_option(photometry.parser.get_option('--gaink'))
key_group.add_option(photometry.parser.get_option('--uik'))
parser.add_option_group(key_group)
customparser.clear_metavars(parser)

def check_run(function, *args):
    """ Run the function and raise CalledProcessError for non-zero retcodes.

    This is a very simple convenience function to feed a function with some
    data, check the returned value and raise subprocess.CalledProcressError in
    case it is other than zero, i.e., if the execution of the function was
    unsuccessful. Note that the function may raise its own errors, so a
    different exception may be raised if something goes wrong.

    """

    retcode = function(*args)
    if retcode:
        cmd = '%s%s' % (function.__name__, args)
        raise subprocess.CalledProcessError(retcode, cmd)

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
        arguments = sys.argv[1:]  # ignore argv[0], the script name
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

    assert len(args) == 1
    offsets_xml_path = args[0]

    # The execution of this module, especially when doing long-term monitoring
    # of reasonably crowded fields, may easily take several *days*. The least
    # we can do, in order to spare the end-user from insufferable grief because
    # of the waste of billions of valuable CPU cycles, is to avoid to have the
    # output file accidentally overwritten.

    if os.path.exists(options.xml_output):
        if not options.overwrite:
            print "%sError. The output file '%s' already exists." % \
                  (style.prefix, options.xml_output)
            print style.error_exit_message
            return 1
        else:
            os.unlink(options.xml_output)

    # To begin with, we need to identify the most constant stars, something for
    # which we have to do photometry on all the stars and for all the images of
    # the campaign. But fret not, as this has to be done only this time: once
    # we get the light curves of all the stars and for all the images, we will
    # be able to determine which are the most constant among them and work
    # always with this subset in order to determine which aperture and sky
    # annulus are the optimal.

    # mkstemp() returns a tuple containing an OS-level handle to an open file
    # and its absolute pathname. Thus, we need to close the file right after
    # creating it, and tell the photometry module to overwrite (-w) it.

    phot_db_handle, phot_db_path = \
        tempfile.mkstemp(prefix = 'photometry_', suffix = '.LEMONdB')
    os.close(phot_db_handle)

    args = [offsets_xml_path,
            '--output', phot_db_path, '--overwrite',
            '--maximum', options.maximum,
            '--margin', options.margin,
            '--cores', options.ncores,
            '--aperture', options.aperture,
            '--annulus', options.annulus,
            '--dannulus', options.dannulus,
            '--min-sky', options.min,
            '--expk', options.exptimek,
            '--coaddk', options.coaddk,
            '--gaink', options.gaink,
            '--uik', options.uncimgk]

    # The --gain option defaults to None, so we add it to the list of arguments
    # only if it was given by the user. Otherwise, it would be given a value od
    # 'None', a string, which would result in an error when attempted to be
    # converted to float by optparse.
    if options.gain:
        args += ['--gain', options.gain]

    # Pass as many '-v' options as we have received here
    [args.append('-v') for x in xrange(options.verbose)]

    print "%sDoing initial photometry with FWHM-derived apertures..." % style.prefix
    print style.prefix

    try:
        # Non-zero return codes raise subprocess.CalledProcessError
        check_run(photometry.main, [str(a) for a in args])

        # Now we need to compute the light curves and find those that are most
        # constant. This, of course, has to be done for each filter, as a star
        # identified as constant in Johnson I may be too faint in Johnson B.
        # In other words: we need to calculate the light curve of each star and
        # for each filter, and then determine which are the options.nconstant
        # stars with the lowest standard deviation.

        diffphot_db_handle, diffphot_db_path = \
            tempfile.mkstemp(prefix = 'diffphot_', suffix = '.LEMONdB')
        os.close(diffphot_db_handle)

        args = [phot_db_path,
                '--output', diffphot_db_path, '--overwrite',
                '--cores', options.ncores,
                '--minimum-images', options.min_images,
                '--stars', options.nconstant,
                '--minimum-stars', options.min_cstars,
                '--pct', options.pct,
                '--weights-threshold', options.wminimum,
                '--max-iters', options.max_iters,
                '--worst-fraction', options.worst_fraction]

        [args.append('-v') for x in xrange(options.verbose)]

        print style.prefix
        print "%sGenerating light curves for initial photometry." % style.prefix
        print style.prefix
        check_run(diffphot.main, [str(a) for a in args])
        print style.prefix

        # Map each photometric filter to the path of the temporary file where
        # the x- and y-, coordinates of each constant star, one per line, will
        # be saved. This file is from now on passed to the photometry module so
        # that it is not done on all stars (after automatically detecting them)
        # but instead exclusively on those pixels.
        pixels_files = {}

        miner = mining.LEMONdBMiner(diffphot_db_path)
        for pfilter in miner.pfilters:

            # LEMONdBMiner.sort_by_curve returns a list of two-element tuples,
            # mapping the ID of each star to the standard deviation of its
            # light curve in this filter. The list is sorted in increasing
            # order by the standard deviation. We are only interested in the
            # first 'options.nconstant', needing at least 'options.pminimum'.

            print "%sIdentifying the %d most constant stars for the %s " \
                  "filter..." % (style.prefix, options.nconstant, pfilter) ,
            sys.stdout.flush()

            stars_stdevs = \
              miner.sort_by_curve_stdev(pfilter, minimum = options.min_images)
            cstars = stars_stdevs[:options.nconstant]

            if len(cstars) < options.pminimum:
                msg = "fewer than %d stars identified as constant in the " \
                      "initial photometry for %s filter" % \
                      (options.pminimum, pfilter)
                raise NotEnoughConstantStars(msg)
            print 'done.'

            if len(cstars) < options.nconstant:
                print "%sBut only %d stars were available. Using them all, " \
                      "anyway." % (style.prefix, len(cstars))

            # Replacing whitespaces with underscores is easier than having to
            # quote the path to the --pixels file in case the name of the
            # filter contains them (as optparse would only see up to the
            # first whitespace).
            prefix = '%s_' % pfilter.name.replace(' ', '_')
            pixels_fd, pixels_files[pfilter] = \
                 tempfile.mkstemp(prefix = prefix, suffix = '.pixels')

            # LEMONdBMiner.get_star returns a five-element tuple with
            # the x and y coordinates, right ascension, declination and
            # instrumental magnitude of the star in the reference image.
            for star_id, stdev in cstars:
                star_x, star_y = miner.get_star(star_id)[:2]
                os.write(pixels_fd, "%f\t%f\n" % (star_x, star_y))
            os.close(pixels_fd)

            print "%sStar coordinates for %s temporarily saved to %s" % \
                  (style.prefix, pfilter, pixels_files[pfilter])

        # The constant stars, the only ones to which we will pay attention from
        # now on, have been identified. So far, so good. Now we will generate
        # the light curves of these stars for each candidate set of photometric
        # parameters. We store the evaluated values in a dictionary in which
        # each filter maps to a list of xmlparse.CandidateAnnuli instances.
        evaluated_annuli = collections.defaultdict(list)

        # Needed in-memory to compute the median FWHM of each filter
        print "%sLoading XML offsets into memory..." % style.prefix ,
        sys.stdout.flush()
        xml_offsets = xmlparse.XMLOffsetFile.load(offsets_xml_path)
        print 'done.'

        for pfilter, path in pixels_files.iteritems():

            print style.prefix
            print "%sFinding the optimal photometric parameters for the " \
                  "%s filter." % (style.prefix, pfilter)

            band_offsets = [offset for offset in xml_offsets
                            if offset.filter == pfilter]

            if len(band_offsets) < options.min_images:
                msg = "fewer than %d images (--mimages option) for %s" % \
                      (options.min_images, pfilter)
                raise NotEnoughConstantStars(msg)

            # The median FWHM of the images in this filter is needed in order
            # to calculate the range of aperture annuli to be evaluated.
            print "%sCalculating the median FWHM for this filter..." % style.prefix,
            sys.stdout.flush()
            band_fwhms = [offset.fwhm for offset in band_offsets]
            fwhm = numpy.median(band_fwhms)
            print 'done.'

            # FWHM to range of pixels conversion
            min_aperture = fwhm * options.lower
            max_aperture = fwhm * options.upper
            annulus      = fwhm * options.sky
            dannulus     = fwhm * options.width

            # The dimensions of the sky annulus remain fixed, while the
            # aperture is in the range [lower * FWHM, upper FWHM], with
            # increments of options.step pixels.
            band_apertures = numpy.arange(min_aperture, max_aperture, options.step)
            assert band_apertures[0] == min_aperture

            print "%sFWHM (%s passband) = %.3f pixels, therefore:" % \
                  (style.prefix, pfilter, fwhm)
            print "%sAperture annulus, minimum = %.3f x %.2f = %.3f pixels " % \
                  (style.prefix, fwhm, options.lower, min_aperture)
            print "%sAperture annulus, maximum = %.3f x %.2f = %.3f pixels " % \
                  (style.prefix, fwhm, options.upper, max_aperture)
            print "%sAperture annulus, step = %.2f pixels, which means that:" % \
                  (style.prefix, options.step)
            print "%sAperture annulus, actual maximum = %.3f + %d x %.2f = %.3f pixels" % \
                  (style.prefix, min_aperture, len(band_apertures),
                   options.step, max(band_apertures))
            print "%sSky annulus, inner radius = %.3f x %.2f = %.3f pixels" % \
                  (style.prefix, fwhm, options.sky, annulus)
            print "%sSky annulus, width = %.3f x %.2f = %.3f pixels" % \
                  (style.prefix, fwhm, options.width, dannulus)
            print "%s%d different aperture annuli in the range [%.2f, %.2f] to be evaluated:" % \
                  (style.prefix, len(band_apertures), band_apertures[0], band_apertures[-1])

            # For each candidate aperture, and only with the images taken in
            # this filter, do photometry on the constant stars and compute the
            # median of the standard deviation of their light curves as a means
            # of evaluating the suitability of this combination of parameters.
            for index, aperture in enumerate(band_apertures):

                aper_phot_db_handle, aper_phot_db_path = \
                    tempfile.mkstemp(prefix = 'photometry_', suffix = '.LEMONdB')
                os.close(aper_phot_db_handle)

                args = [offsets_xml_path,
                        '--output', aper_phot_db_path, '--overwrite',
                        '--passband', pfilter.letter,
                        '--pixels', pixels_files[pfilter],
                        '--maximum', options.maximum,
                        '--margin', options.margin,
                        '--cores', options.ncores,
                        '--aperture-pix', aperture,
                        '--annulus-pix', annulus,
                        '--dannulus-pix', dannulus,
                        '--min-sky', options.min,
                        '--expk', options.exptimek,
                        '--coaddk', options.coaddk,
                        '--gaink', options.gaink,
                        '--uik', options.uncimgk]

                if options.gain:
                    args += ['--gain', options.gain]

                [args.append('-v') for x in xrange(options.verbose)]

                try:
                    print style.prefix
                    check_run(photometry.main, [str(a) for a in args])

                    aper_diff_db_handle, aper_diff_db_path = \
                        tempfile.mkstemp(prefix = 'diffphot_', suffix = '.LEMONdB')
                    os.close(aper_diff_db_handle)

                    args = [aper_phot_db_path,
                            '--output', aper_diff_db_path, '--overwrite',
                            '--cores', options.ncores,
                            '--minimum-images', options.min_images,
                            '--stars', options.nconstant,
                            '--minimum-stars', options.min_cstars,
                            '--pct', options.pct,
                            '--weights-threshold', options.wminimum,
                            '--max-iters', options.max_iters,
                            '--worst-fraction', options.worst_fraction]

                    [args.append('-v') for x in xrange(options.verbose)]

                    check_run(diffphot.main, [str(a) for a in args])

                    # Extract the standard deviations and compute their median.
                    miner = mining.LEMONdBMiner(aper_diff_db_path)
                    cstars = miner.sort_by_curve_stdev(pfilter, minimum = options.min_images)

                    # There must be at most 'nconstant' stars, but there may be
                    # fewer if this aperture causes one or more of the constant
                    # stars to be too faint (INDEF) in so many images as to
                    # prevent their lights curve from being computed.
                    assert len(cstars) <= options.nconstant

                    if len(cstars) < options.pminimum:
                        print "%sJust %d constant stars, fewer than the " \
                              "minimum allowed of %d, had their" % \
                              (style.prefix, len(cstars), options.pminimum)
                        print "%slight curves calculated for this " \
                              "aperture. Ignoring it..." % style.prefix

                    else:
                        stdevs_median = numpy.median([x[1] for x in cstars])
                        params = (aperture, annulus, dannulus, stdevs_median)
                        candidate = xmlparse.CandidateAnnuli(*params)
                        evaluated_annuli[pfilter].append(candidate)
                        print "%sAperture = %.3f, median stdev (%d stars) = %.4f" % \
                              (style.prefix, aperture, len(cstars), stdevs_median)

                # Raised by miner.sort_by_curve_stdev if no stars are found.
                # In other words: for this aperture not even a single star has
                # a light curve with at least options.min_images points; thus,
                # much to our sorrow, we cannot evaluate it.
                except mining.NoStarsSelectedError:
                    msg = "%sNo constant stars for this aperture. Ignoring it..."
                    print msg % style.prefix

                finally:

                    percentage = (index + 1) / len(band_apertures) * 100
                    msg = "%s%s progress: %.2f %%"
                    print msg % (style.prefix, pfilter, percentage)

                    # Temporary databases no longer needed; ignore errors in
                    # case the file cannot be removed or if something went
                    # wrong even before these variables could be defined.

                    try: os.unlink(aper_phot_db_path)
                    except (NameError, OSError): pass

                    try: os.unlink(aper_diff_db_path)
                    except (NameError, OSError): pass

            # Let the user know of the best 'annuli', that is, the one for
            # which the standard deviation of the constant stars is minimal
            best_candidate = min(evaluated_annuli[pfilter],
                                 key = operator.attrgetter('stdev'))

            print "%sBest aperture found at %.3f pixels with stdev = %.4f" % \
                  (style.prefix, best_candidate.aperture, best_candidate.stdev)

        print style.prefix
        print "%sSaving the evaluated annuli to the '%s' XML file ..." % \
              (style.prefix, options.xml_output) ,
        sys.stdout.flush()
        xmlparse.CandidateAnnuli.xml_dump(options.xml_output, evaluated_annuli)
        print 'done.'

        print "%sYou're done ^_^" % style.prefix
        return 0

    finally:

        # Last but not least, delete the databases that were initially used to
        # identify the constant stars and the files to which their coordinates
        # were temporarily saved.

        try: os.unlink(phot_db_path)
        except (NameError, OSError): pass

        try: os.unlink(diffphot_db_path)
        except (NameError, OSError): pass

        try:
            for path in pixels_files.itervalues():
                try: os.unlink(path)
                except OSError: pass
        except UnboundLocalError: pass

if __name__ == "__main__":
    sys.exit(main())

