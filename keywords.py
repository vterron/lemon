#! /usr/bin/env python2

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
group_description = (
"These options customize the FITS keywords from which some required "
"information is extracted. The default values are expected to work "
"well with any file that conforms to the FITS standard. If that is "
"not your case, you should make sure to set these options to the "
"correct values, or otherwise face apocalyptic consequences.")

desc = {} # option descriptions (for optparse)
filterk = 'FILTER'
desc['filterk'] = \
"keyword for the name of the filter of the observation [default: %default]"

rak = 'RA'
desc['rak'] = \
"keyword for the right ascension of the astronomical object(s), expressed " \
"either as a floating point number in decimal degrees, or as a string in " \
"the 'hh:mm:ss[.sss]' format [default: %default]"

deck = 'DEC'
desc['deck'] = \
"keyword for the declination of the astronomical object(s), expressed " \
"either as a floating point number in decimal degrees, or as a string " \
"in the 'dd:mm:ss[.sss]' format [default: %default]"

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

uncimgk = None
desc['uncimgk'] = \
"keyword that stores the path to the uncalibrated image used to check for " \
"saturation -- as the overscan, bias and (particularly) flat-fielding steps " \
"may take a saturated pixel below the saturation threshold. If (as by " \
"default) this option is not set, saturation is checked for on the same " \
"image on which we do photometry."

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

bjdk = 'BJD_TDB'
desc['bjdk'] = (
"keyword for the barycentric JD (TDB, Barycentric Dynamical Time) at "
"mid-exposure. Internally, LEMON stores all times as UTC Unix timestamps, so "
"BJD_TDB is converted to JD_UTC using Jason Eastman's web applet for time "
"correction (htttp://astroutils.astronomy.ohio-state.edu/time/bjd2utc.html) "
"[default: %default]")
