.. _techniques:

Techniques: ``ixdat``'s measurement subclasses
==============================================
Source: https://github.com/ixdat/ixdat/tree/user_ready/src/ixdat/techniques

TechniqueMeasurement classes (interchangeable with Techniques or Measurement subclasses)
inherit from the ``Measurement`` class (see :ref:`measurement`)

A full list of the techniques and there names is in the ``TECHNIQUE_CLASSES`` dictionary::

    >>> from ixdat.techniques import TECHNIQUE_CLASSES
    >>> TECHNIQUE_CLASSES  # note, more techniques may have been added since!
    {
        'simple': <class 'ixdat.measurements.Measurement'>,
        'EC': <class 'ixdat.techniques.ec.ECMeasurement'>,
        'CV': <class 'ixdat.techniques.cv.CyclicVoltammagram'>,
        'MS': <class 'ixdat.techniques.ms.MSMeasurement'>,
        'EC-MS': <class 'ixdat.techniques.ec_ms.ECMSMeasurement'>
    }

.. toctree::
    :maxdepth: 2

    electrochemistry
    mass_spec
    ec_ms
    sec
