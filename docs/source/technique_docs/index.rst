.. _technqiues:

Techniques: ``ixdat``'s measurement subclasses
==============================================

TechniqueMeasurement classes (interchangeable with Techniques or Measurement subclasses)
inherit from the ``Measurement`` class (see :ref:`measurement`)

A full list of the techniques and there names is in the ``TECHNIQUE_CLASSES`` dictionary:

>>> from ixdat.techniques import TECHNIQUE_CLASSES
>>> TECHNIQUE_CLASSES


.. toctree::
    :maxdepth: 2

    electrochemistry
    mass_spec
    ec_ms
