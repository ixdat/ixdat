.. _spectra:

Spectrum
========

The position of spectra is not yet completely set in ixdat.

A spectrum is in essence just a 1-D field where a response variable (e.g. counts,
detector current, adsorption) lives on a space defined by a scanning variable (e.g. wavelength,
mass-to-charge ratio, two-theta).
As such it could be a DataSeries on par with ValueSeries (a 1-D field with a value living in
a space defined by a TimeSeries).
A Spectrum is however a stand-alone output of an experimental technique. It could also be a type of Measurement.

As it is now, ``Spectrum`` is its own base class on par with ``Measurement``, with its own table (i.e. the ``Spectrum`` class
inherits directly from ``Saveable``). It has properties which give quick access to the scanning
variable and response as ``x`` and ``y``, respectively. It also has its
own :ref:`plotter <spectrum-plotter>` and :ref:`exporter <spectrum-exporter>`.

Similar questions can be raised about a sequence of spectra - whether it is a Measurement or a 2-D field.
As it is now, sequences of spectra are represented by ``SpectrumSeries``, which inherits from
``Spectrum``.

The ``spectra`` module
......................

.. automodule:: ixdat.spectra
    :members: