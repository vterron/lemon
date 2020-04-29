#! /usr/bin/env python2
# encoding: UTF-8

# Author: Victor Terron (c) 2020
# Email: `echo vt2rron1iaa32s | tr 132 @.e`
# License: GNU GPLv3

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

_DESCRIPTION = ("""
Print the light curve of an object stored in a LEMONdB.

This command takes as input the right ascension and declination of an object,
and finds the one stored in the LEMONdB that's close to these coordinates. It
then prints to standard output the (a) time, (b) differential magnitude and
(c) signal-to-noise ratio of all the points in the light curve of the object
in the specified photometric filter.
""")

import argparse
import astropy.time
import os.path
import re
import sys

# LEMON modules
import database
import passband
import util.coords

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=_DESCRIPTION)
    parser.add_argument('db_path', metavar='LEMON_DB', type=str,
                        help="the LEMON database with the light curves")
    parser.add_argument('ra', metavar='RA', type=float,
                        help="Right adcension of the astronomical object, "
                             "in decimal degrees.")
    parser.add_argument('dec', metavar='DEC', type=float,
                        help="Declination of the astronomical object, in "
                             "decimal degrees.")
    parser.add_argument('filter', metavar='FILTER', type=str,
                        help="The name of the photometric filter.")

    args = parser.parse_args()

    # TODO(vterron): create the PhotometricFilter object directly with argparse.
    pfilter = passband.Passband(args.filter)

    with database.LEMONdB(args.db_path) as db:
        print("Input coordinates:")
        print("α: {} ({})".format(args.ra, util.coords.ra_str(args.ra)))
        print("δ: {} ({})".format(args.dec, util.coords.dec_str(args.dec)))

        star_id, distance = db.star_closest_to_world_coords(args.ra, args.dec)
        star_info = db.get_star(star_id)
        # TODO(vterron): make get_star return a namedtuple.
        ra, dec = star_info[2:4]

        print()
        print("Selected star:")
        print("ID: {}".format(star_id))
        print("α: {} ({})".format(ra, util.coords.ra_str(ra)))
        print("δ: {} ({})".format(dec, util.coords.dec_str(dec)))
        print("Distance to input coordinates: {} deg".format(distance))

        print()
        print("Light curve in {!r} photometric filter:".format(pfilter))
        star_phot = db.get_photometry(star_id, pfilter)
        if star_phot is None:
            raise ValueError("no light curve for {!r} photometric filter".format(pfilter))

        for index in range(len(star_phot)):
            utime = star_phot.time(index)
            jd    = astropy.time.Time(utime, format='unix').jd
            mag   = star_phot.mag(index)
            snr   = star_phot.snr(index)
            # TODO(vterron): consider using PrettyTable.
            print("{}\t{:.3f}\t{}".format(jd, mag, snr))
