Instructions
============

Use this as a staging ground for CHANGES.rst. In other words, describe the
changes and additions to ixdat's API associated with your contribution. The idea is
that what you write here informs other developers what is on its way and then will be
copied to CHANGES.rst when the next version of ixdat is distributed. Please include
links to relevant Issues, Discussions, and PR's on github with the following format
(replace XX):

`Issue #XX <https://github.com/ixdat/ixdat/issues/XX>`_
`PR #XX <https://github.com/ixdat/ixdat/pull/XX>`_

For ixdat 0.4.1
===============

readers
^^^^^^^

- The ``OceanViewTimeSeriesReader`` (reader="oceanview") has been added for reading
  time series of spectra exported by Ocean Insight's OceanView software. Read a file
  with ``Spectrum.read(path, reader="oceanview")`` to get an ``OpticalSpectrumSeries``.
  Use ``spectra_type`` to say whether the file has intensity (default), transmission,
  or absorption data.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

- The ``AndorKineticsCSVReader`` (reader="andor") has been added for reading kinetic
  series of spectra exported as .csv by Andor Solis. It also returns an
  ``OpticalSpectrumSeries``.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

techniques
^^^^^^^^^^

- Adding an ``ECMeasurement`` or ``CyclicVoltammogram`` to an ``OpticalSpectrumSeries``
  now gives an ``ECOpticalMeasurement``.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

- ``ECOpticalMeasurement`` calculates dOD according to the ``spectra_type`` of its
  spectra. For transmission and absorption data, a reference spectrum is optional.
  The dOD is available as ``sec["dOD"]``, and ``select`` and ``cut`` keep the reference
  spectrum, so that e.g. ``sec.select(cycle=1)["dOD"]`` works.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

- ``CyclicVoltammogram.redefine_cycle(turning_point=True)`` starts a new cycle at the
  start of each potential sweep, so that each cycle is one sweep. With ``redox=True`` or
  ``False``, only anodic or cathodic sweeps start a new cycle. Sweeps are found with the
  same tools as ``get_timed_sweeps``, and set with its ``v_scan_res`` and
  ``res_points``. The new ``merge_short_sections`` in
  ``ixdat.calculators.scan_rate_tools`` keeps noise during a hold from starting new
  cycles.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

calculators
^^^^^^^^^^^

- A calculator can now respond to the end of a series name. For example, the
  calculator handling ``">smooth"`` returns a smoothed version of whatever comes
  before it, so ``sec["dOD>diff>smooth"]`` is the smoothed ``sec["dOD>diff"]``.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

- Three calculators for spectroelectrochemistry have been added in
  ``ixdat.calculators.sec_calculators``. The ``DeltaODCalculator`` comes with every
  ``ECOpticalMeasurement`` and gives ``sec["dOD"]`` and difference spectra
  ``sec["dOD>diff"]``. The ``Denoiser`` gives denoised spectra, e.g.
  ``sec["dOD>diff>smooth"]``. The ``SpectrumDeconvolutor`` finds the converged spectra
  of a measurement and fits the dOD as a sum of them.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

plotters
^^^^^^^^

- ``ECOpticalPlotter`` has new methods for plotting difference spectra, converged
  spectra, the fit and its residuals, and the change between two cycles.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

- ``plot_wavelengths_vs_potential(coplot=True)`` plots the CV on the same axis. The SEC
  heat plots take ``zspan=[min, max]``.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

API changes
^^^^^^^^^^^

- ``ECOpticalMeasurement.plot_wavelengths`` and
  ``ECOpticalMeasurement.plot_wavelengths_vs_potential`` are deprecated, to be removed in
  0.5.0. Use the plotter instead, e.g. ``sec.plotter.plot_wavelengths()``.
  ``plot_waterfall`` stays on the measurement.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_

bug fixes
^^^^^^^^^

- Indexing a ``SpectrumSeries`` returns a single spectrum again. Since 0.4.0 it could
  return a series. `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_
- ``SpectrumSeries.cut`` works when the series has no ``durations``.
  `PR #197 <https://github.com/ixdat/ixdat/pull/197>`_
