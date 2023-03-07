.. _mass-spec:

Mass Spectrometry
=================

Mass spectrometry is commonly used in catalysis and electrocatalysis for two different
types of data - spectra, where intensity is taken while scanning over m/z, and
multiple ion detection (MID) where the intensity of a small set of m/z values are
tracked in time.

The main TechniqueMeasurement class for MID data is the ``MSMeasurement``.

The position of spectra in ixdat is not completely set yet. Currently, they exist as ``Spectrum`` base
class alongside ``Measurement``. For more details see :ref:`spectra`

The ``ms`` module
.................

.. automodule:: ixdat.techniques.ms
    :members: