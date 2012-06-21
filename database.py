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

"""
This module implements LEMONdB, the interface to the databases to which
photometric information (photometry.py), light curves (diffphot.py) and star
periods (periods.py) are saved. These databases contain all the information
relative to the campaign that may be needed for the data analysis process.

"""

import math
import numpy
import operator
import random
import string
import sqlite3
import time

# LEMON modules
import passband

class DBStar(object):
    """ Encapsulates the instrumental photometric information for a star.

    This class is used as a container for the instrumental photometry of a
    star, observed in a specific photometric filter. It implements both
    high-level and low-level routines, the latter of which are fundamental for
    a scalable implementation of the differential photometry algorithms.

    """

    def __init__(self, id_, pfilter, phot_info, times_indexes,
                 dtype = numpy.float128):
        """ Instantiation method for the DBStar class.

        This is an abrupt descent in the abstraction ladder, but needed in
        order to compute the light curves fast and minimize the memory usage.

        Arguments:
        id_ - the ID of the star in the LEMONdB.
        pfilter - the photometric filter of the information being stored.
        phot_info - a two-dimensional NumPy array with the photometric
                    information. It *must* have three rows (the first for the
                    time, the second for the magnitude and the last for the
                    SNR) and as many columns as records for which there is
                    photometric information. For example, in order to get the
                    magnitude of the third image, we would do phot_info[1][2]
        times_indexes - a dictionary mapping each Unix time for which the star
                        was observed to its index in phot_info; this gives us
                        O(1) lookups when 'trimming' an instance. See the
                        BDStar.issubset and complete_for for further
                        information. Note that the values in this dictionary
                        are trusted blindly, so they better have the correct
                        values for phot_info!

        """

        self.id = id_
        self.pfilter = pfilter
        if phot_info.shape[0] != 3: # number of rows
            raise ValueError("'phot_info' must have exactly three rows")
        self._phot_info = phot_info
        self._time_indexes = times_indexes
        self.dtype = dtype

    def __str__(self):
        """ The 'informal' string representation """
        return "%s(ID = %d, filter = %s, %d records)" % \
               (self.__class__.__name__, self.id, self.pfilter, len(self))

    def __len__(self):
        """ Return the number of records for the star """
        return self._phot_info.shape[1] # number of columns

    def time(self, index):
        """ Return the Unix time of the index-th record """
        return self._phot_info[0][index]

    def _time_index(self, unix_time):
        """ Return the index of the Unix time in '_phot_info' """
        return self._time_indexes[unix_time]

    def mag(self, index):
        """ Return the magnitude of the index-th record """
        return self._phot_info[1][index]

    def snr(self, index):
        """ Return the SNR of the index-th record """
        return self._phot_info[2][index]

    @property
    def _unix_times(self):
        """ Return the Unix times at which the star was observed """
        return self._phot_info[0]

    def issubset(self, other):
        """ Return True if for each Unix time at which 'self' was observed,
        there is also an observation for 'other'; False otherwise """

        for unix_time in self._unix_times:
            if unix_time not in other._time_indexes:
                return False
        return True

    def _trim_to(self, other):
        """ Return a new DBStar which contains the records of 'self' that were
        observed at the Unix times that can be found in 'other'. KeyError will
        be raised if self if not a subset of other -- so you should check for
        that before trimming anything"""

        phot_info = numpy.empty((3, len(other)), dtype = self.dtype)
        for oindex, unix_time in enumerate(other._unix_times):
            sindex = self._time_index(unix_time)
            phot_info[0][oindex] = self.time(sindex)
            phot_info[1][oindex] = self.mag(sindex)
            phot_info[2][oindex] = self.snr(sindex)
        return DBStar(self.id, self.pfilter, phot_info,
                      other._time_indexes, dtype = self.dtype)

    def complete_for(self, iterable):
        """ Iterate over the supplied DBStars and trim them.

        The method returns a list with the 'trimmed' version of those DBStars
        which are different than 'self' (i.e., a star instance will not be
        considered to be a subset of itself) and of which it it is a subset.

        """

        complete_stars = []
        for star in iterable:
            if self is not star and self.issubset(star):
                complete_stars.append(star._trim_to(self))
        return complete_stars


class PhotometricParameters(object):
    """ Encapsulates the parameters used for photometry """
    def __init__(self, aperture, annulus, dannulus):
        self.aperture = aperture
        self.annulus  = annulus
        self.dannulus = dannulus


class ReferenceImage(object):
    """ Encapculates the image used for the offsets calculation """
    def __init__(self, path, pfilter, unix_time, airmass, gain):
        self.path = path
        self.pfilter  = pfilter
        self.unix_time = unix_time
        self.airmass = airmass
        self.gain = gain


class Image(ReferenceImage):
    """ Encapsulates the images shifted from the reference one """
    def __init__(self, path, pfilter, pparams, unix_time, airmass, gain,
                 xoffset, yoffset, xoverlap, yoverlap):
        super(Image, self).__init__(path, pfilter, unix_time, airmass, gain)
        self.pparams = pparams  # photometric parameters
        self.xoffset = xoffset
        self.yoffset = yoffset
        self.xoverlap = xoverlap
        self.yoverlap = yoverlap


class LightCurve(object):
    """ The data points of a graph of light intensity of a celestial object.

    Encapsulates a series of Unix times linked to a differential magnitude with
    a signal-to-noise ratio. Internally stored as a list of three-element
    tuples, but we are implementing the add method so that we can interact with
    it as if it were a set, moving us up one level in the abstraction ladder.

    """

    def __init__(self, pfilter, cstars, cweights, dtype = numpy.float128):
        """ 'cstars' is a sequence or iterable of the IDs in the LEMONdB of the
        stars that were used as comparison stars when the light curve was
        computed, while 'cweights' is another sequence or iterable with the
        corresponding weights. The i-th comparison star (cstars) is assigned
        the i-th weight (cweights). The sum of all weights should equal one.

        """

        if len(cstars) != len(cweights):
            msg = "number of weights must equal that of comparison stars"
            raise ValueError(msg)
        if not cstars:
            msg = "at least one comparison star is needed"
            raise ValueError(msg)

        self._data = []
        self.pfilter = pfilter
        self.cstars = cstars
        self.cweights = cweights
        self.dtype = dtype

    def add(self, unix_time, magnitude, snr):
        """ Add a data point to the light curve """
        self._data.append((unix_time, magnitude, snr))

    def __len__(self):
        return len(self._data)

    def __getitem__(self, index):
        return self._data[index]

    def __iter__(self):
        """ Return a copy of the (unix_time, magnitude, snr) tuples,
        chronologically sorted"""
        return iter(sorted(self._data, key = operator.itemgetter(0)))

    @property
    def stdev(self):
        if not self:
            raise ValueError("light curve is empty")
        magnitudes = [mag for unix_time, mag, snr in self._data]
        return numpy.std(numpy.array(magnitudes, dtype = self.dtype))

    def weights(self):
        """ Return a generator over the pairs of comparison stars and their
            corresponding weights """
        for cstar_id, cweight in zip(self.cstars, self.cweights):
            yield cstar_id, cweight


class DuplicateImageError(KeyError):
    """ Raised if two Images with the same Unix time are added to a LEMONdB """
    pass

class DuplicateStarError(KeyError):
    """ Raised if tho stars with the same ID are added to a LEMONdB """
    pass


class LEMONdB(object):
    """ Interface to the SQLite database used to store our results """

    def __init__(self, path, dtype = numpy.float128):

        self.path = path
        self.dtype = dtype
        self.connection = sqlite3.connect(self.path, isolation_level = None)
        self._cursor = self.connection.cursor()

        # Enable foreign key support (SQLite >= 3.6.19)
        self._execute("PRAGMA foreign_keys = ON")
        self._execute("PRAGMA foreign_keys")
        if not self._rows.fetchone()[0]:
            raise sqlite3.NotSupportedError("foreign key support is not enabled")

        self._start()
        self._create_tables()
        self.commit()

    def __del__(self):
        self._cursor.close()
        self.connection.close()

    def _execute(self, query, t = ()):
        """ Execute SQL query; returns nothing """
        self._cursor.execute(query, t)

    @property
    def _rows(self):
        """ Return an iterator over the rows returned by the last query """
        return self._cursor

    def _start(self):
        """ Start a new transaction """
        self._execute("BEGIN TRANSACTION")

    def _end(self):
        """ End the current transaction """
        self._execute("END TRANSACTION")

    def commit(self):
        """ Make the changes of the current transaction permanent.
        Automatically starts a new transaction """
        self._end()
        self._start()

    def _savepoint(self, name = None):
        """ Start a new savepoint, use a random name if not given any.
        Returns the name of the savepoint that was started. """

        if not name:
            name = ''.join(random.sample(string.letters, 12))
        self._execute("SAVEPOINT %s" % name)
        return name

    def _rollback_to(self, name):
        """ Revert the state of the database to a savepoint """
        self._execute("ROLLBACK TO %s" % name)

    def _release(self, name):
        """ Remove from the transaction stack all savepoints back to and
        including the most recent savepoint with this name """
        self._execute("RELEASE %s" % name)

    def _create_tables(self):
        """ Create, if needed, the tables used by the database """

        self._execute('''
        CREATE TABLE IF NOT EXISTS stars (
            id   INTEGER PRIMARY KEY,
            x    REAL NOT NULL,
            y    REAL NOT NULL,
            ra   REAL,
            dec  REAL,
            imag REAL NOT NULL)
        ''')

        self._execute('''
        CREATE TABLE IF NOT EXISTS photometric_filters (
            wavelength INTEGER PRIMARY KEY,
            name       TEXT UNIQUE NOT NULL)
        ''')

        self._execute('''
        CREATE TABLE IF NOT EXISTS photometric_parameters (
            id       INTEGER PRIMARY KEY,
            aperture INTEGER NOT NULL,
            annulus  INTEGER NOT NULL,
            dannulus INTEGER NOT NULL)
        ''')

        self._execute('''
        CREATE TABLE IF NOT EXISTS rimage (
           path       TEXT NOT NULL,
           wavelength INTEGER NOT NULL,
           unix_time  REAL NOT NULL,
           airmass    REAL NOT NULL,
           gain       REAL NOT NULL,
           FOREIGN KEY (wavelength) REFERENCES photometric_filters(wavelength))
        ''')

        self._execute('''
        CREATE TABLE IF NOT EXISTS images (
            unix_time  REAL PRIMARY KEY,
            path       TEXT NOT NULL,
            wavelength INTEGER NOT NULL,
            pparams_id INTEGER NOT NULL,
            airmass    REAL NOT NULL,
            gain       REAL NOT NULL,
            xoffset    REAL NOT NULL,
            xoverlap   INTEGER NOT NULL,
            yoffset    REAL NOT NULL,
            yoverlap   INTEGER NOT NULL,
            FOREIGN KEY (wavelength) REFERENCES photometric_filters(wavelength),
            FOREIGN KEY (pparams_id) REFERENCES photometric_parameters(id))
        ''')

        self._execute('''
        CREATE TABLE IF NOT EXISTS photometry (
            id         INTEGER PRIMARY KEY,
            star_id    INTEGER NOT NULL,
            unix_time  INTEGER NOT NULL,
            magnitude  REAL NOT NULL,
            snr        REAL NOT NULL,
            FOREIGN KEY (star_id)   REFERENCES stars(id),
            FOREIGN KEY (unix_time) REFERENCES images(unix_time),
            UNIQUE (star_id, unix_time))
        ''')

        self._execute('''
        CREATE TABLE IF NOT EXISTS light_curves (
            id         INTEGER PRIMARY KEY,
            star_id    INTEGER NOT NULL,
            unix_time  INTEGER NOT NULL,
            magnitude  REAL NOT NULL,
            snr        REAL,
            FOREIGN KEY (star_id)   REFERENCES stars(id),
            FOREIGN KEY (unix_time) REFERENCES images(unix_time),
            UNIQUE (star_id, unix_time))

        ''')

        self._execute('''
        CREATE TABLE IF NOT EXISTS cmp_stars (
            id         INTEGER PRIMARY KEY,
            star_id    INTEGER NOT NULL,
            wavelength INTEGER NOT NULL,
            cstar_id   INTEGER NOT NULL,
            weight     REAL NOT NULL,
            FOREIGN KEY (star_id)    REFERENCES stars(id),
            FOREIGN KEY (wavelength) REFERENCES photometric_filters(wavelength),
            FOREIGN KEY (cstar_id)   REFERENCES stars(id))
        ''')

        self._execute('''
        CREATE TABLE IF NOT EXISTS periods (
            id         INTEGER PRIMARY KEY,
            star_id    INTEGER NOT NULL,
            wavelength INTEGER NOT NULL,
            step       REAL NOT NULL,
            period     REAL NOT NULL,
            FOREIGN KEY (star_id)    REFERENCES stars(id),
            FOREIGN KEY (wavelength) REFERENCES photometric_filters(wavelength),
            UNIQUE (star_id, wavelength))
        ''')

        self._execute("CREATE INDEX IF NOT EXISTS phot_by_star "
                      "ON photometry(star_id)")
        self._execute("CREATE INDEX IF NOT EXISTS phot_by_unix_time "
                      "ON photometry(unix_time)")
        self._execute("CREATE INDEX IF NOT EXISTS img_by_unix_time "
                      "ON images(unix_time)")
        self._execute("CREATE INDEX IF NOT EXISTS img_by_wavelength "
                      "ON images(wavelength)")
        self._execute("CREATE INDEX IF NOT EXISTS cstars_by_star_and_wavelength "
                      "ON cmp_stars(star_id, wavelength)")
        self._execute("CREATE INDEX IF NOT EXISTS curve_by_star "
                      "ON light_curves(star_id)")
        self._execute("CREATE INDEX IF NOT EXISTS period_by_star "
                      "ON periods(star_id)")
        self._execute("CREATE INDEX IF NOT EXISTS period_by_star_and_wavelength ON "
                      "periods(star_id, wavelength)")

    def _table_count(self, table):
        """ Return the number of rows in 'table' """
        self._execute("SELECT COUNT(*) FROM %s" % table)
        rows = list(self._rows) # from iterator to list
        assert len(rows) == 1
        return rows[0][0]

    def _add_pfilter(self, pfilter):
        """ Store a photometric filter in the database. The primary
        key of the filters in their table is their wavelength """

        t = (pfilter.wavelength, pfilter.name)
        self._execute("INSERT OR IGNORE INTO photometric_filters VALUES (?, ?)", t)

    @property
    def _pparams_ids(self):
        """ Return the ID of the photometric parameters, in ascending order"""
        self._execute("SELECT id "
                      "FROM photometric_parameters "
                      "ORDER BY id ASC")
        return list(x[0] for x in self._rows)

    def _get_pparams(self, id_):
        """ Return the PhotometricParamaters with this ID.
        Raises KeyError if the database has nothing for this ID """

        self._execute("SELECT aperture, annulus, dannulus "
                      "FROM photometric_parameters "
                      "WHERE id = ?", (id_,))
        rows = list(self._rows)
        if not rows:
            raise KeyError('%d' % id_)
        else:
            assert len(rows) == 1
            args = rows[0]
            return PhotometricParameters(*args)

    def _add_pparams(self, pparams):
        """ Add a PhotometricParameters instance and return its ID or do
        nothing and simply return the ID if already present in the database"""

        t = [pparams.aperture, pparams.annulus, pparams.dannulus]
        self._execute("SELECT id "
                      "FROM photometric_parameters "
                      "WHERE aperture = ? "
                      "  AND annulus  = ? "
                      "  AND dannulus = ?", t)
        try:
            return list(self._rows)[0][0]
        except IndexError:
            t.insert(0, None)
            self._execute("INSERT INTO photometric_parameters VALUES (?, ?, ?, ?)", t)
            return self._cursor.lastrowid

    @property
    def rimage(self):
        """ Return a ReferenceImage instance, or None if there isn't any"""
        self._execute("SELECT r.path, p.name, r.unix_time, r.airmass, r.gain "
                      "FROM rimage AS r, photometric_filters AS p "
                      "ON r.wavelength = p.wavelength")
        rows = list(self._rows)
        if not rows:
            return None
        else:
            assert len(rows) == 1
            args = list(rows[0])
            args[1] = passband.Passband(args[1])
            return ReferenceImage(*args)

    @rimage.setter
    def rimage(self, rimage):
        """ Set the reference image that was used to compute the offsets """
        self._add_pfilter(rimage.pfilter)
        t = (rimage.path, rimage.pfilter.wavelength, rimage.unix_time,
             rimage.airmass, rimage.gain)
        self._execute("DELETE FROM rimage")
        self._execute("INSERT INTO rimage VALUES (?, ?, ?, ?, ?)", t)

    def add_image(self, image):
        """ Store an image into the database.

        Raises DuplicateImageError if the Image has the same Unix time than
        another instance already stored in the database (as the primary key
        of the table of images is the date of observation).

        """

        # Use a SAVEPOINT to, if the insertion of the Image fails, be able to
        # roll back the insertion of the filter and photometric parameters

        mark = self._savepoint()
        self._add_pfilter(image.pfilter)
        pparams_id = self._add_pparams(image.pparams)
        t = (image.unix_time, image.path, image.pfilter.wavelength, pparams_id,
             image.airmass, image.gain, image.xoffset, image.xoverlap,
             image.yoffset, image.yoverlap)
        try:
            self._execute("INSERT INTO images "
                          "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", t)
            self._release(mark)

        except sqlite3.IntegrityError:
            self._rollback_to(mark)
            unix_time = image.unix_time
            if __debug__:
                self._execute("SELECT unix_time FROM images")
                assert (unix_time,) in self._rows
            msg = "Unix time %.4f (%s) already in database" % \
                  (unix_time, time.ctime(unix_time))
            raise DuplicateImageError(msg)

    def get_image(self, unix_time):
        """ Return the Image observed at a Unix time. Raises
        KeyError if there is no image for this observation date."""

        self._execute("SELECT i.path, p.name, i.unix_time, i.airmass, i.gain, "
                      "       i.xoffset, i.yoffset, i.xoverlap, i.yoverlap, "
                      "       i.pparams_id "
                      "FROM images AS i, photometric_filters AS p "
                      "ON i.wavelength = p.wavelength "
                      "WHERE i.unix_time = ?", (unix_time,))
        rows = list(self._rows)
        if not rows:
            raise KeyError("%.4f (%s)" % (unix_time, time.ctime(unix_time)))
        else:
            assert len(rows) == 1
            row = rows[0]

            pparams_id = row[-1]
            pparams = self._get_pparams(pparams_id)

            args = list(row[:-1])
            args[1] = passband.Passband(args[1])
            args.insert(2, pparams)
            return Image(*args)

    def add_star(self, star_id, x, y, ra, dec, imag):
        """ Add an star to th database.

        This method only stores the 'description' of the star, that is, its
        image and celestial coordinates, as well as its instrumental magnitude
        in the reference image. To add the photometric records and the light
        curves, use LEMONdB.add_photometry and add_light_curve, respectively.
        Raises DuplicateStarError if the specified ID was already used for
        another star in the database.

        """

        t = (star_id, x, y, ra, dec, imag)
        try:
            self._execute("INSERT INTO stars VALUES (?, ?, ?, ?, ?, ?)", t)
        except sqlite3.IntegrityError:
            if __debug__:
                self._execute("SELECT id FROM stars")
                assert (star_id,) in self._rows
            msg = "star with ID = %d already in database" % star_id
            raise DuplicateStarError(msg)

    def add_photometry(self, star_id, unix_time, magnitude, snr):
        """ Store a photometric record for an star at a given time """
        t = (None, star_id, unix_time, magnitude, snr)
        self._execute("INSERT INTO photometry VALUES (?, ?, ?, ?, ?)", t)

    def get_star(self, star_id, pfilter):
        """ Return the photometric information of the star in a photometric
        filter, encapsulated as a DBStar instance """

        t = (star_id, pfilter.wavelength)
        self._execute("SELECT phot.unix_time, phot.magnitude, phot.snr "
                      "FROM photometry AS phot, images AS img "
                      "ON phot.unix_time = img.unix_time "
                      "WHERE phot.star_id = ? "
                      "  AND img.wavelength = ? "
                      "ORDER BY phot.unix_time ASC", t)

        # NumPy arrays are stored in contiguous blocks of memory, so
        # adding rows or columns to an existing one would require to
        # copy the entire array to a new block of memory. It is much
        # better to first create an array as big as will be needed --
        # numpy.empty(), unlike zeros, does not initializes its
        # entries and may therefore be marginally faster

        rows = list(self._rows)
        phot_info = numpy.empty((3, len(rows)), dtype = self.dtype)

        # A cache, mapping each Unix time to its index in phot_info; passed
        # to the constructor of DBStar for O(1) lookups of Unix times
        times_indexes = {}

        for index, row in enumerate(rows):
            unix_time, magnitude, snr = row
            phot_info[0][index] = unix_time
            phot_info[1][index] = magnitude
            phot_info[2][index] = snr
            times_indexes[unix_time] = index
        return DBStar(star_id, pfilter, phot_info,
                      times_indexes, dtype = self.dtype)

    def get_star_info(self, star_id):
        """ Return the coordinates and magnitude of a star.

        The method returns a five-element tuple with, in this order: the x- and
        y- coordinates of the star in the reference image, the right ascension
        and declination and its instrumental magnitude in the reference image.
        Raises KeyError is no star in the database has this ID.

        """

        t = (star_id, )
        self._execute("SELECT x, y, ra, dec, imag "
                      "FROM stars "
                      "WHERE id = ?", t)
        try:
            return self._rows.next()
        except StopIteration:
            msg = "star with ID = %d not in database" % star_id
            raise KeyError(msg)

    def __len__(self):
        """Return the number of stars in the database """
        return self._table_count('STARS')

    @property
    def star_ids(self):
        """ Return a list with the ID of the stars, in ascending order """
        self._execute("SELECT id FROM stars ORDER BY id ASC")
        return list(x[0] for x in self._rows)

    @property
    def pfilters(self):
        """ Return a sorted list (by wavelength) of the photometric filters (as
        Passband instances) for which there is information stored. The filter
        of the reference image is ignored, so if, say, it was observed in the
        Johnson I filter while the rest of the images of the campaign were
        taken in Johnson B, only the latter will be returned"""

        self._execute("""SELECT DISTINCT f.name
                         FROM photometry AS phot
                         INNER JOIN images AS img
                         ON phot.unix_time = img.unix_time
                         INNER JOIN photometric_filters AS f
                         ON img.wavelength = f.wavelength
                         ORDER BY f.wavelength ASC """)

        return [passband.Passband(x[0]) for x in self._rows]

    def _add_curve_point(self, star_id, unix_time, magnitude, snr):
        """ Store a point of the light curve of a star """

        # Note the casts to Python's built-in float otherwise, if the method
        # receives a NumPy float, SQLite will raise "sqlite3.InterfaceError:
        # Error binding parameter - probably unsupported type"
        t = (None, star_id, float(unix_time), float(magnitude), float(snr))
        self._execute("INSERT INTO light_curves "
                      "VALUES (?, ?, ?, ?, ?)", t)

    def _add_cmp_star(self, star_id, pfilter, cstar_id, cweight):
        """ Add a comparison star to the light curve of a star.

        The method stores 'cstar_id' as the ID of one of the comparison stars,
        with a weight of 'cweight', that were used to compute the light curve
        of the star with ID 'star_id' in the 'pfilter' photometric filter.

        """

        # Note the cast to Python's built-in float otherwise, if the method
        # receives a NumPy float, SQLite will raise "sqlite3.InterfaceError:
        # Error binding parameter - probably unsupported type"
        t = (None, star_id, pfilter.wavelength, cstar_id, float(cweight))
        self._execute("INSERT INTO cmp_stars "
                      "VALUES (?, ?, ?, ?, ?)", t)

    def add_light_curve(self, star_id, light_curve):
        """ Store a light curve in the database """
        for point in light_curve:
            self._add_curve_point(star_id, *point)
        for cstar_id, cweight in light_curve.weights():
            self._add_cmp_star(star_id, light_curve.pfilter, cstar_id, cweight)

    def get_light_curve(self, star_id, pfilter):
        """ Return the light curve of a star in a filter.

        None is returned if for this star and filter the database has no
        information about the comparison stars or points of the light curve.
        Both cases mean, assuming there were no errors during the database
        population, that there is no light curve for this star"""

        # The instantiation method of the LightCurve class receives the two
        # sequences with the IDs of the comparison stars and their weights
        t = (star_id, pfilter.wavelength)
        self._execute("SELECT cstar_id, weight "
                      "FROM cmp_stars "
                      "WHERE star_id = ? "
                      "  AND wavelength = ? ", t)
        cstars   = []
        cweights = []
        for cstar, cweight in self._rows:
            cstars.append(cstar)
            cweights.append(cweight)

        if not cstars:
            assert not cweights
            return None

        assert len(cstars) == len(cweights)
        curve = LightCurve(pfilter, cstars, cweights, dtype = self.dtype)
        self._execute("SELECT curve.unix_time, curve.magnitude, curve.snr "
                      "FROM light_curves AS curve, images AS img "
                      "ON curve.unix_time = img.unix_time "
                      "WHERE curve.star_id = ? "
                      "  AND img.wavelength = ?", t)

        for row in self._rows:
            curve.add(*row)

        if not curve:
            return None
        else:
            return curve

    def add_period(self, star_id, pfilter, step, period):
        """ Store the string-length period of a star """
        t = (None, star_id, pfilter.wavelength, step, float(period))
        self._execute("INSERT INTO periods "
                      "VALUES (?, ?, ?, ?, ?)", t)

    def get_period(self, star_id, pfilter):
        """ Return a two-element tuple with the string-length period of a star
        and the step that was used when the candidate periods where evaluated
        (in other words, the precision with which the period was determined) """

        t = (star_id, pfilter.wavelength)
        self._execute("SELECT period, step "
                      "FROM periods "
                      "WHERE star_id = ? "
                      "  AND wavelength = ?", t)
        try:
            return list(self._rows)[0]
        except IndexError:
            return None, None

    def get_periods(self, star_id):
        """ Return a NumPy array with the periods, in the different filters, of
        the star. This is only a convenience function to retrieve all the
        periods of the star -- if, for example, we need to examine how similar
        they are. Note that they may be returned in any order, so there is no
        way of knowing to which photometric filter each one belongs. For that,
        use the LEMONdB.get_period method instead."""

        t = (star_id,)
        self._execute("SELECT period "
                      "FROM periods "
                      "WHERE star_id = ? ", t)
        return numpy.array(list(x[0] for x in self._rows))

    def airmasses(self, pfilter):
        """ Return a dictionary which maps each Unix time for which there
        is information in this filter to the corresponding airmass"""

        t = (pfilter.wavelength,)
        self._execute("SELECT unix_time, airmass "
                      "FROM images "
                      "WHERE wavelength = ? ", t)
        return dict(self._rows)

    def get_phase_diagram(self, star_id, pfilter, period, repeat = 1):
        """ Return the folded light curve of a star.

        The method returns a LightCurve instance with the phase diagram of the
        star in a photometric filter: 'Phase diagrams (also known as 'folded
        light curves') are a useful tool for studying the behavior of periodic
        stars such as Cepheid variables and eclipsing binaries. In a phase
        diagram, multiple cycles of brightness variation are superimposed on
        each other. Instead of plotting magnitude versus Julian date as with a
        regular light curve, each observation is plotted as a function of 'how
        far into the cycle' it is' [http://www.aavso.org/data/lcg/curve.shtml].

        The 'repeat' keyword argument determines how many times the cycle is
        repeated, in order help us more easily appreciate what the shape of the
        period is. A 'phased Unix time' of 0,05, for example, becomes 1.05 the
        first time the phase diagram is repeated, 2.05 the second time, etc.

        """

        curve = self.get_light_curve(star_id, pfilter)
        if curve is None:
            return None

        phase = LightCurve(pfilter, curve.cstars,
                           curve.cweights, dtype = curve.dtype)
        unix_times, magnitudes, snrs = zip(*curve)
        zero_t = min(unix_times)

        phased_x = []
        for utime, mag, snr in zip(unix_times, magnitudes, snrs):
            # How far into the cycle is this Unix time?
            fractional_part, integer_part = math.modf((utime - zero_t) / period)
            phased_x.append(fractional_part)
        assert len(phased_x) == len(unix_times)

        x_max = 1;
        phased_unix_times = phased_x[:]
        for _ in xrange(repeat - 1): # -1 as there is already one (phased_x)
            phased_unix_times += [utime + x_max for utime in phased_x]
            x_max += 1;

        assert len(phased_unix_times) == len(unix_times) * repeat
        phased_magnitudes = magnitudes * repeat
        phased_snr = snrs * repeat

        for utime, mag, snr in \
            zip(phased_unix_times, phased_magnitudes, phased_snr):
                phase.add(utime, mag, snr)

        assert len(phase) == len(curve) * repeat
        return phase

