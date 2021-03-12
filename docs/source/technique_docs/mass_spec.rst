.. _mass-spec:

Mass Spectrometry
=================
Source: https://github.com/ixdat/ixdat/tree/user_ready/src/ixdat/techniques/ms

Mass spectrometry is commonly used in catalysis and electrocatalysis for two different
types of data - spectra, where intensity is taken while scanning over m/z, and
mass intensity detection (MID) where the intensity of a small set of m/z values are
tracked in time.

The main TechniqueMeasurement class for MID data is the :ref:`MSMeasurement`.

Classes dealing with spectra are under development.

The ``ms`` module
........................

.. automodule:: ixdat.techniques.ms
    :members: