#! /usr/bin/env python

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

description = """
The purpose of this module, the first of the LEMON pipeline, is to
automatically detect all the FITS files belonging to an observation campaign
and save them to another directory. This is particularly needed when reducing
images taken at observatories that enforce specific naming conventions, such as
using a different directory for each night. In situations like that, the
astronomer may easily end up with hundreds of images scattered over dozens of
directories, many times also with duplicate filenames. As if, for example, the
first image of the night must be named 'ferM_0001.fits' and so on, there will
be _many_ files with that name, this script may be used in order to
sequentially rename all the images.

What this script does, in short, is to: (a) automatically detect all the FITS
files within a directory tree, (b) discard those saturated (i.e., whose median
number of ADUs is above a certain threshold) or whose keyword with the name for
the object observed does not match any of the given patterns, such as 'bias',
'skyflat*' or 'ngc2264*', (c) sort them by their date of observation and,
finally, (d) copy them to the output directory, renaming them
sequentially. Some statistics on the imported FITS files are shown to the user
before the script exits. It is enforced that all the imported images have the
same size: in case there are multiple dimensions among the detected FITS
images, only those with the most common size will be imported, while the others
will be ignored.

"""

import collections
import fnmatch
import numpy
import operator
import optparse
import os
import os.path
import pyfits
import re
import shutil
import stat
import sys

# LEMON modules
import customparser
import keywords
import fitsimage
import methods
import style

parser = customparser.get_parser(description)
parser.usage = "%prog [OPTION]... INPUT_DIRS... OUTPUT_DIR"

parser.add_option('--object', action = 'callback', type = 'str',
                  dest = 'objectn', default = ['*'],
                  callback = methods.str_split_callback,
                  help = "list of case-insensitive patterns, according to the "
                  "rules used by the Unix shell, and separated by commas, of "
                  "the object names to import. Those FITS images whose object "
                  "keyword (see --objectk option) matches none of these "
                  "patterns will be ignored [default: %default]")

parser.add_option('--pattern', action = 'store', type = 'str',
                  dest = 'pattern',
                  help = "pattern, according to the rules used by the Unix "
                  "shell, that the filename of a FITS image must match if "
                  "it is to be considered when scanning the directories. "
                  "Non-matching images will be ignored.")

parser.add_option('--counts', action = 'store', type = 'int',
                  dest = 'max_counts', default = None,
                  help = "median number of counts, or ADUs "
                  "(analog-to-digital units) at which saturation arises. "
                  "Images above this value will be ignored. If not set, "
                  "no image is discarded because of its median ADUs "
                  "[default: %default]")

parser.add_option('--filename', action = 'store', type = 'str',
                  dest = 'filename',
                  help = "name shared by all output images, to which the "
                  "sequence number will be appended. If not set, it will be "
                  "automatically detected by finding the most common "
                  "filename among the input images.")

parser.add_option('--follow', action = 'store_true', default = False,
                  dest = 'followlinks',
                  help = "walk down into symbolic links that resolve to "
                  "directories, on systems that support them. This can lead "
                  "to infinite recursion if a link points to a parent "
                  "directory of itself.")

parser.add_option('--exact', action = 'store_true', default = False,
                  dest = 'exact',
                  help = "do not modify the imported files, but just rename "
                  "their exact copies. This means that FITS headers are not "
                  "altered and, in particular, that the --uik keyword is not "
                  "added. The SHA-1 hash is used to verify that the copy of "
                  "the FITS images is identical.")

key_group = optparse.OptionGroup(parser, "FITS Keywords",
                                 keywords.group_description)

key_group.add_option('--datek', action = 'store', type = 'str',
                     dest = 'datek', default = keywords.datek,
                     help = keywords.desc['datek'])

key_group.add_option('--timek', action = 'store', type = 'str',
                     dest = 'timek', default = keywords.timek,
                     help = keywords.desc['timek'])

key_group.add_option('--expk', action = 'store', type = 'str',
                     dest = 'exptimek', default = keywords.exptimek,
                     help = keywords.desc['exptimek'])

key_group.add_option('--objectk', action = 'store', type = 'str',
                     dest = 'objectk', default = keywords.objectk,
                     help = keywords.desc['objectk'])

key_group.add_option('--uik', action = 'store', type = 'str',
                     dest = 'uncimgk', default = keywords.uncimgk,
                     help = "keyword to which the path of each imported image "
                     "is written to its own FITS header. In this manner, when "
                     "later on we do photometry we can use the original FITS "
                     "image, before any possible calibration was performed, "
                     "to check for saturation -- as the overscan, bias and "
                     "(particularly) flat-fielding steps may take a saturated "
                     "pixel below the saturation threshold.")

parser.add_option_group(key_group)
customparser.clear_metavars(parser)

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

    # Print the help message and abort the execution if there are not two
    # positional arguments left after parsing the options, as the user must
    # specify the path to both the input and output directories.

    if len(args) < 2:
        parser.print_help()
        return 2  # 2 is generally used for command line syntax errors
    else:
        input_dirs = args[:-1]
        output_dir = args[-1]

    # Make sure that all the input directories exist, abort otherwise.
    for path in input_dirs:
        if not os.path.exists(path):
            print "%sThe input directory, '%s', does not exist. Exiting." % \
                  (style.prefix, path)
            return 1

    # The input and output directories must be different, as otherwise some
    # files (especially if the filename of the output files is automatically
    # detected) could be overwritten.
    for path in input_dirs:
        if os.path.abspath(path) == os.path.abspath(output_dir):
            print "%s[INPUT_DIRS] and OUTPUT_DIR must be different. " \
                  "Exiting." % style.prefix
            return 1

    # Make sure that the output directory exists, create it otherwise
    methods.determine_output_dir(output_dir)

    # Recursively walk down the input directories, obtaining a list of all the
    # regular files. Then, and while a progress bar is shown to let the user
    # estimate how much longer it is, detect which among them are FITS files.

    print "%sIndexing regular files within directory trees starting at " \
          "INPUT_DIRS..." % style.prefix ,
    files_paths = fitsimage.find_files(input_dirs,
                                       followlinks = options.followlinks,
                                       pattern = options.pattern)
    print 'done.'

    print "%sDetecting FITS images among the %d indexed regular files..." % \
          (style.prefix, len(files_paths))

    images_set = fitsimage.FITSet()
    methods.show_progress(0.0)
    for path_index, path in enumerate(files_paths):
        try:
            images_set.append(fitsimage.FITSImage(path))
            fraction = (path_index + 1) / len(files_paths) * 100
            methods.show_progress(fraction)
        except fitsimage.NonStandardFITS:
            pass
    else:
        methods.show_progress(100)
        print

    if not len(images_set):
        print "%sNo FITS files were found. Exiting." % style.prefix
        return 1
    else:
        print "%s%d FITS files detected." % (style.prefix, len(images_set))

    # All the images must have the same size; otherwise, only those with the
    # most common dimensions will be imported, while the rest will be ignored
    print style.prefix
    print "%sChecking the sizes of the detected images..." % style.prefix,
    img_sizes = collections.defaultdict(int) # dimensions counter
    for img in images_set:
        img_sizes[img.size] += 1
    print 'done.'

    # The most common size is the only one element in case len(img_sizes) == 1
    x_size, y_size = max(img_sizes.iterkeys(), key = img_sizes.get)[:2]

    if len(img_sizes) == 1:
        print "%sAll the FITS images have the same size: %d x %d pixels" % \
              (style.prefix, x_size, y_size)
    else:

        print "%sMultiple sizes were detected among the FITS images." % style.prefix
        print "%sDiscarding images with a size other than %d x %d pixels, " \
              "the most common..." % (style.prefix, x_size, y_size) ,
        old_size = len(images_set)
        selected = [img for img in images_set if img.size == (x_size, y_size)]
        images_set = fitsimage.FITSet(selected)
        print 'done.'

        if not images_set:
            print "%sThere are no FITS files left. Exiting." % style.prefix
            return 1
        else:
            print "%s%d FITS files were discarded because of their size, " \
                  "%s remain." % (style.prefix, old_size - len(images_set),
                                  len(images_set))

    # Those FITS images whose object names do not match any of the given
    # patterns, or which do not even have the keyword which contains the
    # name for the object observed, are discarded.
    print style.prefix
    print "%sImporting only those FITS files whose %s keyword can be found " \
          "and matches" % (style.prefix, options.objectk)
    print "%sone of the following Unix patterns: %s ..." % \
          (style.prefix, options.objectn)

    # We first test that the keyword exists (hence the pass for the KeyError
    # exception, which means that the image is filtered out) and, after that,
    # check whether its value matches one of the regular expressions which
    # define the object names to be imported.
    object_set = fitsimage.FITSet()

    # Keep the track of how many images are ignored for each reason
    saturated_excluded = 0
    non_match_excluded = 0

    for img in images_set:

        try:
            object_name = img.read_keyword(options.objectk)
            for pattern in options.objectn:
                regexp = re.compile(fnmatch.translate(pattern), re.IGNORECASE)
                if regexp.match(object_name):
                    # Even if the object name matchs, the median number of
                    # counts must still be below the threshold, if any. If the
                    # number of ADUs is irrelevant we can avoid having to
                    # unnecessarily compute it.
                    if options.max_counts:
                        with pyfits.open(img.path, readonly = True) as hdu:
                            median_counts = numpy.median(hdu[0].data)
                        if median_counts > options.max_counts:
                            print "%s%s excluded (matched, but saturated " \
                                  "with %d ADUs)" % (style.prefix, img.path,
                                                     median_counts)
                            saturated_excluded += 1
                            break

                    # This point reached if median number of ADUs of image is
                    # above the threshold or irrelevant, so it can be imported.
                    print "%s%s imported (%s matches '%s')" % (style.prefix,
                           img.path, object_name, pattern)

                    object_set.append(img)
                    break

            else: # only executed if for loop exited cleanly
                print "%s%s excluded (%s does not match anything)" % \
                      (style.prefix, img.path, object_name)
                non_match_excluded += 1
        except KeyError:
            pass

    if not saturated_excluded and not non_match_excluded:
        print "%sNo images were filtered out. Hooray!" % style.prefix
    if saturated_excluded:
        print "%s%d files were discarded because they were saturated " \
              "(> %d ADUs)." % (style.prefix, saturated_excluded,
                                options.max_counts)
    if non_match_excluded:
        print "%s%d files were discarded because of their non-matching " \
              "object names." % (style.prefix, non_match_excluded)

    # Abort the execution if all the FITS files were filtered out
    if not object_set:
        print "%sThere are no FITS files left. Exiting." % style.prefix
        return 1

    # Sort the FITS files by their date of observation, according to the header
    print style.prefix
    print "%sSorting the FITS files by their date of observation " \
          "[keyword: %s]..." % (style.prefix, options.datek) ,

    kwargs = dict(date_keyword = options.datek,
                  time_keyword = options.timek,
                  exp_keyword = options.exptimek)
    get_date = operator.methodcaller('date', **kwargs)
    sorted_set = sorted(object_set, key = get_date)

    # Let the user know if one or more images could not be sorted (because of
    # problems when parsing the FITS keywords from which the observation date
    # is derived) and thus discarded.
    difference = len(object_set) - len(sorted_set)
    assert difference >= 0
    if difference:
        print
        print "%s%d files were discarded as the observation date keyword " \
              "was not found or the " % (style.prefix, difference)
        print "%sdate in it represented did not conform to the FITS " \
              "standard." % style.prefix

        # Execution is aborted if all the FITS files were filtered out
        if not sorted_set:
            print "%sThere are no FITS files left. Exiting." % style.prefix
            return 1
    else:
        print 'done.'

    # If no filename for the output images was specified, attempt to
    # automatically detect the most common basename among the FITS files.
    # This is doing by extracting the leftmost non-numeric substring of
    # all the filenames and taking that which repeats the most.

    if not options.filename:
        print style.prefix
        print "%sDetecting the most common name among input files..." % \
              style.prefix ,
        sys.stdout.flush()

        # Use a dictionary in order to keep the track of how many times we
        # have come across each prefix (leftmost non-numeric substring in
        # the filename) and select that with most occurrences.

        prefixes = collections.defaultdict(int)
        for prefix in (img.prefix for img in sorted_set):
            prefixes[prefix] += 1

        # Select the prefix (key) that is repeated the most
        options.filename = max(prefixes, key = prefixes.get)
        print 'done.'

    print "%sImported FITS filenames will start with the string: '%s'" % \
          (style.prefix, options.filename)

    # Now we have to copy the FITS files. The basename of each imported file
    # will be options.filename + its sequence number. Filling zeros will be
    # affixed to each number so that the lenth of all the basenames is the
    # same. Following Dijkstra's teachings, we start numbering at zero.

    assert len(sorted_set)
    ndigits = len(str(len(sorted_set) - 1))
    print "%s%d digits are needed in order to enumerate %d files." % \
          (style.prefix, ndigits, len(sorted_set))

    print style.prefix
    print "%sCopying the FITS files to '%s'..." % \
          (style.prefix, output_dir)

    for index, fits_file in enumerate(sorted_set):

        # i.e., 'ferM_' + '0000' + '.fits' = 'ferM_0000.fits'
        dest_name = '%s%0*d.fits' % (options.filename, ndigits, index)
        dest_path = os.path.join(output_dir, dest_name)

        shutil.copy2(fits_file.path, dest_path)

        # The permission bits have been copied, but we need to make sure
        # that the copy of the FITS file is always writable, no matter what
        # the original permissions were. This is equivalent to `chmod u+w`
        mode = os.stat(dest_path)[stat.ST_MODE] | stat.S_IWUSR
        os.chmod(dest_path, mode)

        dest_img = fitsimage.FITSImage(dest_path)

        # Add some information to the FITS header...
        if not options.exact:

            msg1 = "File imported by LEMON on %s" % methods.utctime()
            dest_img.add_history(msg1)

            # If the --uik option is given, store in this keyword the absolute
            # path to the image of which we made a copy. This allows other
            # LEMON commands, if necessary, to access the original FITS files
            # in case the imported images are modified (e.g., bias subtraction
            # or flat-fielding) before these other commands are executed.

            if options.uncimgk:

                comment = "before any calibration task"
                dest_img.update_keyword(options.uncimgk,
                                        os.path.abspath(dest_img.path),
                                        comment = comment)

                msg2 = "[Import] Original image: %s"
                dest_img.add_history(msg2 % os.path.abspath(fits_file.path))

        # ... unless we want an exact copy of the images. If that is the case,
        # verify that the SHA-1 checksum of the original and the copy matches
        elif fits_file.sha1sum != dest_img.sha1sum:
            msg = "copy of %s not identical (SHA-1 differs)" % fits_file.path
            raise IOError(msg)

        # Show which file has been copied, using the format of the
        # 'cp -v' command: `./ultra2/ferM_11.fits' -> `imported/img_01.fits'
        print  "%s`%s' -> `%s'" % (style.prefix, fits_file.path, dest_path)

    # Finally, let the user know how many FITS images, and the fraction of
    # the total, that were imported, as well as their size in megabytes.
    print style.prefix
    ifraction = len(sorted_set) / len(images_set) * 100
    print "%sFITS files detected: %d" % (style.prefix, len(images_set))
    print "%sFITS files successfully imported: %d (%.2f%%)" % \
          (style.prefix, len(sorted_set), ifraction)

    total_size = 0.0
    for fits_file in sorted_set:
        total_size += os.path.getsize(fits_file.path) # in bytes

    print "%sTotal size of imported files: %.2f MB" % \
          (style.prefix, total_size / (1024.0 ** 2))
    print "%sYou're done ^_^" % style.prefix
    return 0

if __name__ == "__main__":
    sys.exit(main())

