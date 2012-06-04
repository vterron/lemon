LEMON
=====

LEMON is a CCD differential-photometry pipeline, written in Python, developed at the `Institute of Astrophysics of Andalusia (CSIC) <http://www.iaa.es/>`_ and originally designed for its use at the `1.23m CAHA telescope <http://www.caha.es/telescopes-overview-and-instruments-manuals.html/>`_ for automated variable stars detection and analysis. The aim of this tool is to make it possible to *completely reduce thousands of images of time series* in a matter of hours and with minimal user interaction, if not none at all, automatically detecting variable stars and presenting the results to the astronomer.

A first overview of LEMON was presented some time ago at `<http://adsabs.harvard.edu/abs/2011hsa6.conf..755T>`_.


Modules
=======

The pipeline currently consists of **ten stages**, which are usually run sequentially, although depending on your needs only a specific subset of them may be used. Proper and detailed documentation, the bane of many a programmer, is still pending and will probably be not written until there is definitive proof that somebody other than us, the developers, is using it. In the meantime, you should be able to get by with the scripts help, which can be accessed by running them without arguments.

1. import.py — group all the images of a campaign
#. seeing.py — find image with best seeing
#. offsets.py — determine translation offsets
#. mosaic.py — create master frame
#. astrometry.py — astrometric calibration
#. annuli.py — find best parameters for photometry
#. photometry.py — simply do photometry
#. diffphot.py — generate light curves
#. periods.py — string-length method
#. mining.py — identify *interesting* stars
