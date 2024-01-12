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

readers
^^^^^^^

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
