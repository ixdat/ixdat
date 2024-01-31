Instructions
============

Use this as a staging ground for CHANGES.rst. In other words, describe the
changes and additions to ixdat's API associated with your contribution. The idea is
that what you write here informs other developers what is on its way and then will be
copied to CHANGES.rst when the next version of ixdat is distributed. Please include
links to relevant Issues, Discussions, and PR's on github with the following format
(replace XX):

`Issue #XX <https://github.com/ixdat/ixdat/issues/XX>`_

`PR #XX <https://github.com/ixdat/ixdat/pulls/XX>`_


ixdat 0.2.9
===========


API changes
-----------

readers
^^^^^^^
- ``BiologicReader`` can now also read biologic .mpr files using an external package.
  When reading a file ending in ".mpr", it first tries the ``galvani`` package, which
  seems to work for LSV, CA, and CVA files. If that fails, it tries the ``eclabfiles``
  package, which seems to work for OCV and CP files. See:
  - https://github.com/echemdata/galvani
  - https://github.com/vetschn/eclabfiles
  These packages are not added as a requirement, but instead imported dynamically.
  If the user tries to read a .mpr file without the needed package installed, they are
  pointed to the package but also encouraged to export and read ".mpt" files instead.
  ".mpr" files are recognized as biologic, and ``reader="biologic"`` works for both types.
  Resolves `Issue #132 <https://github.com/ixdat/ixdat/issues/132`_


- Zilien MS spectrum reader fixed.
  Resolves `Issue #117 <https://github.com/ixdat/ixdat/issues/117>`_

- ``Spectrum.read(reader="zilien")`` rather than ``reader="zilien_spec"`` as
  before for reading in a zilien spectrum. This is accomplished by different
  groups of reader for ``Spectrum.read`` and ``Measurement.read``
  Also, zilien-read spectra now know their duration.

- ``Measurement.read(..., reader="zilien")`` now returns a ``SpectroMSMeasurement``
  when the reader can find zilien mass scans taken during the measurement. It
  looks for the mass scans folder as zilien names it.
  The default plotter is a ``SpectroMSPlotter`` which includes the MS spectra
  data in a separate panel. The spectra are accessible by:

    meas = Measurement.read("my_MID_with_spectra.tsv", reader="zilien")
    meas.spectrum_series[0].plot()

  which plots the first MS spectrum.
  To leave out the mass scan data, include the argument ``include_spectra=False``
  in the call to ``read``.
  This finishes `Issue #117 <https://github.com/ixdat/ixdat/issues/117`_

- If a series name is present in the raw data *and* in in a measurement's ``aliases``,
  the raw data series matching the name and the aliased series are appended. (Before,
  only the raw data series matching the name would be returned.)

techniques
^^^^^^^^^^
- ECOpticalMeasurement.get_spectrum() now has an option not to interpolate.
  Passing the argument `interpolate=False` gets it to choose nearest spectrum instead.

- Indexing a ``SpectroMeasurement`` with an integer returns a ``Spectrum``.
  For example, ``zilien_meas_with_spectra[0].plot()``  plots the first mass scan

plotters
^^^^^^^^
- ``SpectrumPlotter.heat_plot()`` and methods that rely on it can now plot discrete heat plots, with
  each spectrum only plotted for its duration, if available. If spectrum durations are not available,
  it plots each spectrum until the start of the next spectrum, i.e. like the previous (continuous)
  behaviour but without interpolation.
  Discrete heat plotting is now the default behavior for ``MSSpectroMeasurement``s read by "zilien"
  (in which case durations are available).
  ``ECOpticalMeasurement``s read by "msrh_sec" have an unchanged (continuous) default plot.
  resolves `Issue #140 <https://github.com/ixdat/ixdat/issues/140

General
^^^^^^^

- The string representation, which is what is printed if an object is printed, has been
  changed for ``TimeSeries``, ``ValueSeries`` and ``Measurement``. The data series have
  changed, so that they will return a helpful summary, displaying the name, min, max and
  the timestamp for a ``TimeSeries`` as opposed to the class name and ``__init__``
  argument form, which is normally inherited from ``__repr__``. In short::

    Before: TimeSeries(id=1, name='Potential time [s]')
    After: TimeSeries: 'Potential time [s]'. Min, max: 699, 1481s @ 21B01 17:44:12

    Before: ValueSeries(id=2, name='Voltage [V]')
    After: ValueSeries: 'Voltage [V]'. Min, max: 1.4e+00, 5.4e+00 [V]

  These new string representations should be helpful on their own, but the main goal of
  changing them, was to make them useful in the new string representation of
  ``Measurement``, which is inherited by all measurements. It will now display a summary
  of all data series in the measurement, along with information about their internal
  connections, like which ``ValueSeries`` depends on which ``TimeSeries``. For an ECMS
  measurement the form is::

    ECMSMeasurement '2021-02-01 17_44_12' with 48 series

    Series list:
    ┏ TimeSeries: 'Potential time [s]'. Min, max: 699, 1481 [s] @ 21B01 17:44:12
    ┣━ ValueSeries: 'Voltage [V]'. Min, max: 1.4e+00, 5.4e+00 [V]
    ┣━ ValueSeries: 'Current [mA]'. Min, max: -2.5e-02, 2.5e-02 [mA]
    ┗━ ValueSeries: 'Cycle [n]'. Min, max: 1.0e+00, 1.0e+00 [n]
    ┏ TimeSeries: 'Iongauge value time [s]'. Min, max: 1, 3042 [s] @ 21B01 17:44:12
    ┗━ ValueSeries: 'Iongauge value [mbar]'. Min, max: 6.6e-09, 6.9e-07 [mbar]
    << SNIP MORE SYSTEM CHANNELS >>
    ┏ TimeSeries: 'C0M2 time [s]'. Min, max: 1, 3041 [s] @ 21B01 17:44:12
    ┗━ ValueSeries: 'M2 [A]'. Min, max: 3.4e-12, 2.0e-11 [A]
    ┏ TimeSeries: 'C1M4 time [s]'. Min, max: 1, 3041 [s] @ 21B01 17:44:12
    ┗━ ValueSeries: 'M4 [A]'. Min, max: 1.2e-17, 2.7e-10 [A]
    << SNIP MORE MASS CHANNELS>>

- Reading measurement without the need to specify ``technique`` keyword argument.
  The technique is determined from dataset's metadata. The ``MSMeasurement`` is used
  when it is a Mass Spec measurement. And when it includes an electrochemistry
  data, then ``ECMSMeasurement`` is used. The default/safe case is ``MSMeasurement``.
  Resolves `Issue #159 <https://github.com/ixdat/ixdat/pull/159>`_

dev
^^^

- Enable running external tests in CI
