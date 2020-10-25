#! /usr/bin/env python2
# encoding: UTF-8

# Author: Victor Terron (c) 2020
# Email: `echo vt2rron1iaa32s | tr 132 @.e`
# License: GNU GPLv3

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

_DESCRIPTION = """
Print the light curve of an object stored in a LEMONdB.

This command takes as input the right ascension and declination of an object,
and finds the one stored in the LEMONdB that's close to these coordinates. It
then prints to standard output the (a) time, (b) differential magnitude and
(c) signal-to-noise ratio of all the points in the light curve of the object
in the specified photometric filter.
"""

import argparse
import astropy.time
import os.path
import prettytable
import re
import sys
import time

# LEMON modules
import database
import passband
import util.coords
import util

parser = argparse.ArgumentParser(description=_DESCRIPTION)
parser.add_argument(
    "db_path",
    metavar="LEMON_DB",
    type=str,
    help="the LEMON database with the light curves",
)
parser.add_argument(
    "ra",
    metavar="<right ascension>",
    type=float,
    help="Right adcension of the astronomical object, " "in decimal degrees.",
)
parser.add_argument(
    "dec",
    metavar="<declination>",
    type=float,
    help="Declination of the astronomical object, in " "decimal degrees.",
)
parser.add_argument(
    "filter",
    metavar="<photometric filter>",
    type=passband.Passband,
    help="The name of the photometric filter.",
)
parser.add_argument(
    "--decimal_places",
    dest="places",
    type=int,
    default=3,
    help="Round floating-point numbers to this many decimal places.",
)
parser.add_argument(
    "--output_file",
    dest="output",
    type=argparse.FileType("w"),
    default=sys.stdout,
    help="File to which to write the light curve data points",
)


def main(arguments=None):

    if arguments is None:
        arguments = sys.argv[1:]
    args = parser.parse_args(args=arguments)

    with database.LEMONdB(args.db_path) as db:
        print("Input coordinates:")
        print("α: {} ({})".format(args.ra, util.coords.ra_str(args.ra)))
        print("δ: {} ({})".format(args.dec, util.coords.dec_str(args.dec)))

        star_id, distance = db.star_closest_to_world_coords(args.ra, args.dec)
        star = db.get_star(star_id)

        print()
        print("Selected star:")
        print("ID: {}".format(star_id))
        print("α: {} ({})".format(star.ra, util.coords.ra_str(star.ra)))
        print("δ: {} ({})".format(star.dec, util.coords.dec_str(star.dec)))
        print("Distance to input coordinates: {} deg".format(distance))
        print()

        if args.output == sys.stdout:
            print("Light curve in {!r} photometric filter:".format(args.filter))

        star_diff = db.get_light_curve(star_id, args.filter)
        if star_diff is None:
            raise ValueError(
                "no light curve for {!r} photometric filter".format(args.filter)
            )

        table = prettytable.PrettyTable()
        table.field_names = ["Date (UTC)", "JD", "Δ Mag", "SNR"]

        def format_float(f):
            """Returns f as a string rounded to parser.places decimal places."""
            return "{:.{places}f}".format(f, places=args.places)

        for unix_time, magnitude, snr in star_diff:
            jd = astropy.time.Time(unix_time, format="unix").jd
            table.add_row(
                [
                    util.utctime(unix_time, suffix=False),
                    format_float(jd),
                    format_float(magnitude),
                    format_float(snr),
                ]
            )

        args.output.write(str(table))
        args.output.write("\n")

        if args.output != sys.stdout:
            print(
                "Wrote light curve in {!r} photometric filter to {!r}.".format(
                    args.filter, args.output.name
                )
            )


if __name__ == "__main__":
    sys.exit(main())
