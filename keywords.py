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

# Description of the optparse.OptionGroup
group_description = \
"These options customize the FITS keywords from which some information " \
"needed by the module is extracted. The default values are those of the " \
"images taken with the CCD SITE#2b instrument at the 1,23m CAHA " \
"telescope. Setting one of these options to an incorrect value will " \
"almost certainly result in apocalyptic consequences."


desc = {} # option descriptions (for optparse)
filterk = 'FILTER'
desc['filterk'] = \
"keyword for the name of the filter of the observation [default: %default]"

rak = 'RA'
desc['rak'] = \
"keyword for the right ascension of the astronomical object(s), " \
"in decimal degrees [default: %default]"

deck = 'DEC'
desc['deck'] = \
"keyword for the declination of the astronomical object(s), in " \
"decimal degrees [default: %default]"

datek = 'DATE-OBS'
desc['datek'] = \
"keyword for the date of the observation, in the new Y2K compliant " \
"date format: 'yyyy-mm-dd' or 'yyyy-mm-ddTHH:MM:SS[.sss] " \
"[default: %default]"

timek = 'TIME-OBS'
desc['timek'] = \
"keyword for the time at which the observation started, in the format " \
"HH:MM:SS[.sss]. This keyword is used in conjunction with --datek to " \
"determine the starting time of the observation: --datek gives the " \
"starting calendar date and this keyword the time within that day. This " \
"keyword is not necessary (and thus this option ignored) if the time is " \
"included directly as part of the --datek keyword value with the format " \
"yyyy-mm-ddTHH:MM:SS[.sss] [default: %default]"

exptimek = 'EXPTIME'
desc['exptimek'] = \
"keyword for the exposure time [default: %default]"

airmassk = 'AIRMASS'
desc['airmassk'] = \
"keyword for the airmass [default: %default]"

gaink = 'GAIN'
desc['gaink'] = \
"keyword for the gain of the CCD, in e-/ADU. Needed in order to " \
"accurately calculate the SNR of each measurement [default: %default]"

uncimgk = 'UNCIMG'
desc['uncimgk'] = \
"keyword for the path to the uncalibrated image. This will be one " \
"used to check whether pixels are saturated -- as the overscan, bias and " \
"(particularly) flat-fielding steps may make a pixel be below the " \
"saturation threshold, when in actuality, before the calibration was " \
"performed, they were above. This value may be set to an empty string " \
"('') if saturation is to be checked for on the same image in which " \
"photometry is done [default: %default]"

fwhmk = 'LEMON FWHM'
desc['fwhmk'] = \
"keyword for the Full Width at Half Maximum (FWHM) of the image, which is " \
"written to the FITS header by the 'seeing' command [default: %default]"

objectk = 'OBJECT'
desc['objectk'] = \
"keyword for the name of the object observed [default: %default]"

typek = 'IMAGETYP'
desc['typek'] = \
"keyword that identifies the type of image, with values such as 'dark', " \
"'flat' or 'object', to cite some of the most common [default: %default]"

# Used by seeing.FITSeeingImage to 'cache' the SExtractor catalog
sex_catalog = 'SEX-CAT'
sex_md5sum  = 'SEX-MD5'

coaddk = 'NCOADDS'
desc['coaddk'] = \
"keyword for the number of effective coadds. This value is essential to " \
"determine the number of counts at which saturation arises in coadded " \
"observations. If the keyword is missing, we assume a value of one (that " \
"is, that the observation consisted of a single exposure) [default: %default]"
