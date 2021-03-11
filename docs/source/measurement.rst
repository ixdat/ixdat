.. _measurement:

The measurement structure
=========================

The **measurement** (``meas``) is the central object in the pluggable structure of ixdat, and the
main interface for user interaction. A measurement is an object of the generalized class
main interface for user interaction. A measurement is an object of the generalized class
``Measurement``, defined in the ``measurements`` module, or an inheriting
***TechniqueMeasurement*** class defined in a module of the ``techniques`` folder
(see :ref:`techniques`_).

The general pluggable structure is defined by ``Measurement``, connecting every
measurement to a *reader* for importing from text, a *backend* for saving and loading in
``ixdat``, a *plotter* for visualization, and an *exporter* for saving outside of ``ixdat``.
Each TechniqueMeasurement class will likely hav its own default reader, plotter, and
exporter, while an ``ixdat`` session will typically work with one backend handled by the
``db`` model.

.. image:: figures/pluggable.svg
  :width: 400
  :alt: Design: pluggability

Classes for measurement techniques
----------------------------------

Inheritance in TechniqueMeasurement classes makes it so that related techniques
can share functionality. Here is an illustration of the role of inheritence, using
EC, MS, and EC-MS as an example:

.. image:: figures/inheritance.svg
  :width: 400
  :alt: Design: inheritance

A full list of TechniqueMeasurements is in :ref:`techniques`_.

Initiating a measurement
------------------------

A typical workflow is to start by reading a file. For convenience, most readers are
accessible directly from ``Measurement``. So, for example, to read a .mpt file exported
by Biologic's EC-Lab, one can type:

>>> from ixdat import Measurement
>>> ec_meas = Measurement.read("my_file.mpt", reader="biologic")

See :ref:`readers`_ for a description of the available readers.

The biologic reader (``ixdat.readers.biologic.BiologicMPTReader``) ensures that the
object returned, ``ec_meas``, is of type ``ECMeasurement``.

Another workflow starts with loading a measurement from the active ``ixdat`` backend.
This can also be done straight from ``Measurement``, as follows:

>>> from ixdat import Measurement
>>> ec_meas = Measurement.get(3)

Where the row with id=3 of the measurements table represents an electrochemistry
measurement. Here the column "technique" in the measurements table specifies which
TechniqueMeasurement class is returned. Here, it for row three of the measurements
table, the entry "technique" is "EC", ensuring ``ec_meas`` is an object of type
``ECMeasurement``.

What's in a measurement
-----------------------
A measurement is basically a wrapper around a collection of ``data_series`` (see
:ref:`data_series`).

There are several ways of interracting with a measurement's ``data_series``:

- ``meas.grab()`` is the canonical way of getting numerical data out of a
  measurement. Given the name of a ``ValueSeries``, it returns two numpy arrays, ``t`` and ``v``
  where ``t`` is the time (wrt ``meas.tstamp``) and ``v`` is the value as a function of that
  time vector. ``grab`` takes a series name as its first argument and can also take a ``tspan``
  argument in which case it cuts the vectors to return data for the specific timespan of
  the measurement.
- Indexing a measurement with the name of a data series returns that data series, with
  any time values tstamp'd at ``meas.tstamp``
- Most TechniqueMeasureemnts provide attribute-style access to essential DataSeries and
  data. For example, ``ECMeasurement`` has properties for ``potential`` and ``current`` series,
  as well as ``t``, ``v``, and ``j`` for data.
- The names of the series are available in ``meas.series_names``.
- The raw series are available in ``meas.series_list``.


The ``measurements`` module
---------------------------
Here is the full in-line documentation of the ``measurements`` module containing the
``Measurement`` class.

.. automodule:: ixdat.measurements
    :members:
