.. _commands-index:

##############
LEMON commands
##############

The LEMON pipeline consists of ten commands, implementing different atronomical
data reduction and analysis steps, that are usually run sequentially, although
depending on your needs only a specific subset of them may be used. You may,
for example, be interested in just aligning your FITS images or update their
headers with astrometric information. In this sense, and following the Unix
tools philosophy (*write programs that do one thing and do it well*), LEMON
can be viewed as a set of tasks that *may* be used as a pipeline.

.. toctree::
   :numbered:

   import.rst
