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

import functools
import os.path

GLADE_DIR = os.path.join(os.path.dirname(__file__), './gui/')
get = functools.partial(os.path.join, GLADE_DIR)

GUI_MAIN  = get('main.glade')
GUI_ABOUT = get('about.glade')
GUI_OVERVIEW = get('overview.glade')
LOADING_DIALOG = get('loading-dialog.glade')
STAR_DETAILS = get('star-details.glade')
AMPLITUDES_DIALOG = get('amplitudes-search-dialog.glade')
AMPLITUDES_RESULTS = get('amplitudes-search-results.glade')
SNR_THRESHOLD_DIALOG = get('snr-threshold-dialog.glade')

