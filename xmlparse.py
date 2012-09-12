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

import calendar
import collections
import lxml.etree
import operator
import os.path
import re
import time
import xml.dom.minidom

# LEMON modules
import passband

def validate_dtd(path):
    """ Validate a XML file against a DTD.

    The method validates an Extensible Markup Language (XML) against a Document
    Type Definition (DTD) referenced by the document, raising the appropiate
    exception if an error is encountered, and doing nothing otherwise.

    """

    dtd_parser = lxml.etree.XMLParser(dtd_validation = True)
    lxml.etree.parse(path, dtd_parser)

def xml_header(version = 1.0, encoding = 'utf-8', standalone = True):
    """ Return a XML declaration.

    Return a stromg with a XML declaration, the processing instruction that
    identifies the document as being XML and with which all XML documents
    should begin, situated at the first position of the first line. The
    declaration includes the 'version', 'encoding' and 'standalone'
    attributes. The returned string includes the newline character.

    """
    return """<?xml version="%.1f" encoding="%s" standalone="%s" ?>\n""" % \
           (version, encoding, standalone and 'yes' or 'no')

def toreallyprettyxml(xml_code):
    """ Fix minidom.toprettyxml's silly whitespace.

    This is a not-very-Pythonic hack to solve the problem that toprettyxml adds
    an extra white space when printing the contents of text nodes, which is not
    only ugly but also problematic in case we need to re-parse the 'pretty' XML
    code. This fix uses a regular expression to match a newline character
    followed by tabs characters if they come between '>' and other text, or
    between text and '>'. In other words: it allows us to strip out the
    unwanted whitespace added by toprettyxml().

    The credit goes to BrendanM, who posted it at:
    http://ronrothman.com/public/leftbraned/xml-dom-minidom-toprettyxml-and-silly-whitespace/

    """

    fix = re.compile(r'((?<=>)(\n[\t]*)(?=[^<\t]))|(?<=[^>\t])(\n[\t]*)(?=<)')
    return re.sub(fix, '', xml_code)

class XMLOffset(object):
    """ The translation offset between two FITS images.

    This class encapsulates the translation offset between two images, the
    'reference' one and the 'shifted' one, which displacement is given, in
    pixels, in terms of the reference image.

    """

    def __init__(self, reference, shifted, shifted_filter, shifted_date,
                 x_offset, y_offset, x_overlap, y_overlap):
        """ Instantiation method for the XMLOffset class.

        reference - the path to the reference, 'unmoved' FITS image.
        shifted - the path to the displaced FITS image.
        shifted_filter - the photometric filter of the shifted image,
                         encapsulated as a passband.Passband instance.
        shifted_date - the data of observation of the shifted image, in
                       seconds after the Unix epoch (aka Unix time)
        x_offset - the offset, in pixels, in the x-axis.
        y_offset - the offset, in pixels, in the y-axis.
        x_overlap - the number of stars that overlapped in the x-axis when
                    the offset was calculated and the images were aligned.
        y_overlap - the number of stars that overlapped in the y-axis when
                    the offset was calculated and the images were aligned.

        """

        self.reference = reference
        self.shifted   = shifted
        self.filter    = shifted_filter
        self.date      = shifted_date
        self.x         = x_offset
        self.y         = y_offset
        self.x_overlap = x_overlap
        self.y_overlap = y_overlap

    def __cmp__(self, other):
        """ Comparison operation, sorts XMLOffset instances by their date """
        return self.date - other.date


class XMLOffsetFile(object):
    """ Interface to write and read XMLOffsets to a standalone XML file.

    This class is used as a container of XMLOffset instances, whose XML
    representation can be written to a XML file and also read from it. These
    XML files are self-contained, which means that the DTD declaration is at
    the top of the document, right after the XML declaration.

    Note that, although kept in memory in seconds after the Unix epoch, dates
    are written to disk, for user's convenience, as a string of the following
    form: 'Sun Jun 20 23:21:05 1993 UTC' (Coordinated Universal Time).

    """

    STRPTIME_FORMAT = "%a %b %d %H:%M:%S %Y UTC" # default format + 'UTC'

    XML_DTD = [
    "",
    "<!DOCTYPE offsets [",
    "<!ELEMENT offsets (offset+)>",
    "<!ATTLIST offsets size CDATA  #REQUIRED>",
    "",
    "<!ELEMENT offset (reference, shifted, x_offset, y_offset)>",
    "<!ELEMENT reference (#PCDATA)>",
    "<!ELEMENT shifted   (#PCDATA)>",
    "<!ATTLIST shifted date   CDATA  #REQUIRED>",
    "<!ATTLIST shifted filter CDATA  #REQUIRED>",
    "<!ELEMENT x_offset  (#PCDATA)>",
    "<!ATTLIST x_offset overlap CDATA #REQUIRED>",
    "<!ELEMENT y_offset  (#PCDATA)>",
    "<!ATTLIST y_offset overlap CDATA #REQUIRED>",
    "]>",
    ""]

    def __init__(self, path = None):
        """ Load a XML file into memory, or create a new one """

        if not path:
            self.root = lxml.etree.Element('offsets', size = '0')
        else:
            # There is no need to open the file to parse it with lxml, but we
            # do it here to check that it exists and is readable by the user.
            # This gives us a clearer error message than what lxml would raise
            # if it failed because the file did not exist (would be something
            # like "Error reading file: failed to load external entity)
            with open(path, 'r') as _: pass
            self.root = lxml.etree.parse(path).getroot()

    def __len__(self):
        """ Return the number of XMLOffsets saved in the XML file """
        return len(self.root)

    def add(self, offset):
        """ Add a XMLOffset to the XML file.

        This updates the in-memory contents of the XML file, adding to it
        another XMLOffset. These changes are not written to disk until the
        XMLOffsetFile.dump method is called.

        """

        element = lxml.etree.Element('offset')
        reference = lxml.etree.Element('reference')
        reference.text = offset.reference
        element.append(reference)

        shifted_date_str = "%s UTC" % time.asctime(time.gmtime(offset.date))
        kwargs = {'date' : shifted_date_str, 'filter' : str(offset.filter)}
        shifted = lxml.etree.Element('shifted', **kwargs)
        shifted.text = offset.shifted
        element.append(shifted)

        x = lxml.etree.Element('x_offset', overlap = str(offset.x_overlap))
        x.text = str(offset.x)
        element.append(x)

        y = lxml.etree.Element('y_offset', overlap = str(offset.y_overlap))
        y.text = str(offset.y)
        element.append(y)

        self.root.append(element)

        # Update the 'nimages' attribute of the root element
        self.root.set('size', str(len(self)))

    def _toxml(self, encoding = 'utf-8'):
        """ Return the XML representation of the XMLOffsets.

        The method returns a string with the standalone XML representation of
        the XMLOffsets that have been added to the instance, ready to be
        written to disk. Includes the XML header and the DTD declaration.

        """

        kwargs = {'encoding' : encoding, 'xml_declaration': True,
                  'pretty_print' : True, 'standalone' : True}
        xml_content = lxml.etree.tostring(self.root, **kwargs)

        # Insert between the XML header and the content a comment with the
        # creation time of the file and the XML Document Type Definition.
        lines = xml_content.split('\n')
        comment = "<!-- File generated by LEMON on %s UTC -->"
        lines.insert(1, comment % time.asctime(time.gmtime()))
        lines = lines[:2] + self.XML_DTD + lines[2:]
        return '\n'.join(lines)

    def dump(self, path, encoding = 'utf-8'):
        """ Write the XMLOffsets to a XML file.

        The method saves the XMLOffsets to a standalone XML file, silently
        overwriting it if it already exists. The output document includes
        the XML header and the DTD declaration.

        """

        with open(path, 'wt') as fd:
            fd.write(self._toxml(encoding = encoding))
        validate_dtd(path)

    def __getitem__(self, index):
        """ Return the index-th XMLOffset in the XML file """

        offset = self.root[index]
        reference_path = offset[0].text

        element = offset[1]
        shifted_path = element.text
        # From string to struct_time in UTC to Unix seconds
        shifted_date_str = element.get('date')
        args = shifted_date_str, self.STRPTIME_FORMAT
        shifted_date = calendar.timegm(time.strptime(*args))
        shifted_pfilter = passband.Passband(element.get('filter'))

        element = offset[2]
        x_offset = float(element.text)
        x_overlap = int(element.get('overlap'))

        element = offset[3]
        y_offset = float(element.text)
        y_overlap = int(element.get('overlap'))

        args = (reference_path, shifted_path, shifted_pfilter, shifted_date,
                x_offset, y_offset, x_overlap, y_overlap)
        return XMLOffset(*args)


class CandidateAnnuli(object):
    """ Encapsulates the quality of a set of photometric parameters.

    How do we determine how 'good' a set of aperture, annulus and dannulus
    values are for photometry? What we do is to look at the median (or even the
    arithmetic mean, for this matter both approaches are statistically sound)
    standard deviation of the light curves of the most constant stars. It
    follows that the better (i.e., most appropiate for the images being
    reduced) the parameters, the lower this standard deviation will be.

    This class simply encapsulates these four values. You may think of it as a
    surjective function (as two different sets of parameters may result in the
    same values) which links a three-element tuple with the parameters used for
    photometry (aperture, annulus, dannulus) to the standard deviation of the
    light curves of the most constant stars.

    """

    def __init__(self, aperture, annulus, dannulus, stdev):
        """ Instantiation method.

        aperture - the aperture radius, in pixels.
        annulus - the inner radius of the sky annulus, in pixels.
        dannulus - the width of the sky annulus, in pixels.
        stdev - the median, arithmetic mean or a similar statistical measure
                of the standard deviation of the light curves of the evaluated
                stars when photometry is done using these aperture, annulus
                and dannulus values.

        """

        self.aperture = aperture
        self.annulus  = annulus
        self.dannulus = dannulus
        self.stdev    = stdev

    def __repr__(self):
        return "%s(%f, %f, %f, %f)" % (self.__class__.__name__, self.aperture,
                                       self.annulus, self.dannulus, self.stdev)

    @staticmethod
    def xml_dump(xml_path, annuli):
        """ Save multiple CadidateAnnuli instances to a XML file.

        This method dumps to a file the XML representation of a dictionary
        which maps each photometric filter to a list of the CandidateInstances
        that for it were evaluated. This offers a functionality similar to that
        of the pickle module, with the additional advantages of being
        human-readable, easily understood and parseable virtually everywhere.

        The generated XML file is a standalone document, which means that the
        Document Type Definitions (DTD), defining the document structure with a
        list of legal elements, is also included. This information is used by
        the XML processor in order to validate the code.

        xml_path - the path to which to save the XML file. Any existing file
                   will be mercilessly overwritten without warning.
        annuli - a dictionary mapping each photometric filter to a list of
                 CandidateInstances, encapsulating the quality of a set of
                 photometric parameters.

        """

        header = xml_header(standalone = True)

        xml_dtd = [
        "",
        "<!DOCTYPE annuli [",
        "<!ELEMENT annuli (band*)>",
        "",
        "<!ELEMENT band (candidate*)>",
        "<!ATTLIST band name     CDATA #REQUIRED>",
        "<!ATTLIST band aperture CDATA #REQUIRED>",
        "<!ATTLIST band annulus  CDATA #REQUIRED>",
        "<!ATTLIST band dannulus CDATA #REQUIRED>",
        "<!ATTLIST band stdev    CDATA #REQUIRED>",
        "",
        "<!ELEMENT candidate EMPTY>",
        "<!ATTLIST candidate aperture CDATA #REQUIRED>",
        "<!ATTLIST candidate annulus  CDATA #REQUIRED>",
        "<!ATTLIST candidate dannulus CDATA #REQUIRED>",
        "<!ATTLIST candidate stdev    CDATA #REQUIRED>",
        "]>",
        ""]

        # This approach is more elegant than having to manually add the
        # trailing newline to each line of the XML Document Type Definition.
        xml_dtd = '\n'.join(xml_dtd)

        impl = xml.dom.minidom.getDOMImplementation()
        dom  = impl.createDocument(None, 'annuli', None)
        root = dom.documentElement

        for pfilter in sorted(annuli.iterkeys()):

            # Identify the candidate parameters for which the standard
            # deviation of the curves of the most constant stars is minimal.
            best = min(annuli[pfilter], key = operator.attrgetter('stdev'))

            # Three attributes (i.e., values within the start-tag) identify the
            # optimal photometric parameters for this photometric filters. In
            # this manner, the aperture and sky annuli that have to be used for
            # this photometric filter can be directly and easily extracted from
            # the XML tree. The median or arithmetic mean of the standard
            # deviation of the light curves of the constant stars when
            # photometry was done in order to evaluate these parameters
            # is also stored for debugging purposes.

            band_element = dom.createElement('band')
            band_element.setAttribute('name', pfilter.name)
            band_element.setAttribute('aperture', '%.5f' % best.aperture)
            band_element.setAttribute('annulus',  '%.5f' % best.annulus)
            band_element.setAttribute('dannulus', '%.5f' % best.dannulus)
            band_element.setAttribute('stdev',    '%.8f' % best.stdev)

            # Although most of the time only the optimal photometric parameters
            # will be of interest, it is also worth saving all the aperture and
            # sky annuli that were evaluated and the median or mean stardard
            # deviation that resulted from using them. Note that the optimal
            # parameters are included here again, as the purpose of this
            # listing is to provide a compendium of all the photometric
            # parameters that were taken into consideration. The photometric
            # parameters will be listed sorted on two keys: the aperture
            # annulus itself (primary) and the sky annulus (secondary).

            annuli[pfilter].sort(key = operator.attrgetter('annulus', 'aperture'))
            for candidate in annuli[pfilter]:
                cand_element = dom.createElement('candidate')
                cand_element.setAttribute('aperture', '%.5f' % candidate.aperture)
                cand_element.setAttribute('annulus',  '%.5f' % candidate.annulus)
                cand_element.setAttribute('dannulus', '%.5f' % candidate.dannulus)
                cand_element.setAttribute('stdev',    '%.8f' % candidate.stdev)
                band_element.appendChild(cand_element)

            root.appendChild(band_element)

        fout = open(xml_path, 'wt')
        fout.write(header)
        fout.write("<!-- File generated by LEMON on %s UTC -->\n" % \
                   time.asctime(time.gmtime()))
        fout.write(xml_dtd + '\n')
        fout.write(toreallyprettyxml(root.toprettyxml()))
        fout.close()

        # Make sure that the written XML is valid
        validate_dtd(xml_path)

    @staticmethod
    def xml_load(xml_path, best_only = False):
        """ Load a series of CandidateAnnuli instances from a XML file.

        This method reverses the functionality of xml_dump(), reading a XML
        file and returning a dictionary which maps each photometric filter to
        a list of the CandidateAnnuli instances that were saved to it.

        By default, all the photometric parameters (encapsulated as
        CandidateAnnuli) that were evaluated for each passband are returned. In
        case only the best (and thus the ones that in theory should be used for
        photometry) are of interest, 'best_only' should be set to True so that
        only the optimal CandidateAnnuli for each passband is parsed.

        The XML files generated by LEMON are always standalone documents,
        meaning that the Document Type Definitions (DTD), defining the document
        structure with a list of legal elements, are also included. If the XML
        file cannot be validated, the appropiate exception will be raised.

        xml_path - the path to the XML file to which the CandidateAnnuli
                   instances were saved and from which they will be now loaded.

        Keyword arguments:
        best_only - if True, only the optimal CandidateAnnuli for for each
                    photometric filter is parsed. Otherwise the entire XML
                    file (and thus all the Cxinstances) are loaded.

        """

        # validate_dtd() will also return False (instead of an IOError
        # exception) if the file does not exists, so we raise it manually
        if not os.path.exists(xml_path):
            raise IOError("'%s' does not exist" % xml_path)

        validate_dtd(xml_path)

        dom = xml.dom.minidom.parse(xml_path)
        root = dom.getElementsByTagName('annuli')[0]

        # For each passband, the optimal aperture and sky annuli are stored
        # as attributes of the <band> entity, so that they can be directly
        # extracted in case we are only interested in the optimal paramaters
        # for photometry. Note that, when all the candidate annuli are to be
        # extracted, these attributes can be safely ignored, as the best
        # parameters are also listed as <candidate> entities for each
        # photometric passband.

        annuli = collections.defaultdict(list)
        for band in root.getElementsByTagName('band'):
            pfilter = passband.Passband(band.getAttribute('name'))
            if best_only:
                annuli[pfilter].append(
                    CandidateAnnuli(
                        *[float(x) for x in band.getAttribute('aperture'),
                                            band.getAttribute('annulus'),
                                            band.getAttribute('dannulus'),
                                            band.getAttribute('stdev')]))
            else:
                for candidate in band.getElementsByTagName('candidate'):
                    annuli[pfilter].append(
                        CandidateAnnuli(
                            *[float(x) for x in candidate.getAttribute('aperture'),
                                                candidate.getAttribute('annulus'),
                                                candidate.getAttribute('dannulus'),
                                                candidate.getAttribute('stdev')]))
        return annuli

