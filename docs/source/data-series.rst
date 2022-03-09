.. _`data_series`:

The data series structure
=========================

A **data series** is an object of the ``DataSeries`` class in the ``data_series`` module
or an inheriting class. It is basically a wrapper around a numpy array, which is its
``data`` attribute, together with a name and a unit. Most data series also contain some
additional metadata and/or references to other data series. The most important function
of these is to keep track of everything in time, as described below.

data series and time
--------------------

(Copied from text in design workshop 2, in December 2020):

Time is special!
In some deeper way, time is just another dimension…
but for hyphenated laboratory measurements, as well as multi-technique experimental projects in general (so long as samples and equipment are moving slow compared to the speed of light), time is special because it is the one measurable quantity that is always shared between all detectors.

Absolute time (epoch timestamp) exists in two places:

- ``Measurement.tstamp``: This timestamp is a bit decorative – it tells the measurement’s plotter and data selection methods what to use as t=0
- ``TSeries.tstamp``: This timestamp is truth. It defines the t=0 for the primary time data of any measurement.

Data carriers:

- ``DataSeries``: The ``TimeSeries`` is a special case of the ``DataSeries``. All data
  carried by ixdat will be as a numpy array in a Series. All Series share a primary key
  (id in table series in the db diagram on the left), and in addition to the data have a ``name`` (think “column header”) and a ``unit``. Series is a table in the ``ixdat`` database structure, with helper tables for special cases.
- ``TimeSeries``: The only additional row for ``TimeSeries`` (table tseries) is ``tstamp``, as described above.
- ``Field``: Some series consist of values spanning a space defined by other series. Such a series is called a ``Field``, and defined by a list of references to the series which define their space. In the database, this is represented in a field_axis table, of which n rows (with axis_number from 0 to n-1) will be dedicated to representing each n-D ``Field``.
- ``ValueSeries``: Finally, a very common series type is a scalar value over time. This is called a ``ValueSeries``, and must have a corresponding ``TimeSeries``. A ``ValueSeries`` is actually a special case of a ``Field``, spanning a 1-d space, and so doesn’t need a new table in the db.

**Immutability!** All of the data carriers above will be immutable! This means that,
even though truth is preserved by adding ``dt`` to a ``tsteries.tstamp`` and subtracting
``dt`` from ``tseries.data``, we will never do this! This is a cheap calculation that
``ixdat`` can do on demand. Same with appending corresponding
series from consecutive measurements. Performing these operations on every series in
a measurement set is referred to as building a combined measurement, and is only done
when explicitly asked for (f.ex. to export or save the combined measurement). Building
makes new Series rather than mutating existing ones. A possible exception to immutability
may be appending data to use ``ixdat`` on an ongoing measurement.


The ``data_series`` module
--------------------------

.. automodule:: ixdat.data_series
    :members: