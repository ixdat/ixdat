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


ixdat 0.2.7
===========

API changes
-----------

techniques
^^^^^^^^^^

- ``MSInlet.gas_flux_calibration_curve()`` has the additional option of passing
  a boolean ``axes_measurement_raw``. Set to True if the axes passed to 
  ``axes_measurement`` are raw signal data (i.e. not background subtracted)
  Mentioned in `Issue #94 <https://github.com/ixdat/ixdat/issues/94`_

- ``ECMSMeasurement.ecms_calibration_curve()`` has the additional option of
  forcing the line of best fit through zero with ``force_through_zero``
  Resolves `Issue #109 <https://github.com/ixdat/ixdat/issues/109`_

- ``ECMSMeasurement.ecms_calibration_curve()`` has the additional option of 
  passing a J_name to be used for highlighting the integrated current passed to
  ``axes_measurement``. This does not affect the calculation of sensitivity factors, 
  only plotting.
  Resolves `Issue #118 <https://github.com/ixdat/ixdat/issues/118`_

readers
^^^^^^^
- ``Spectrum.read(reader="zilien")`` rather than ``reader="zilien_spec"`` as 
  before for reading in a zilien spectrum. This is accomplished by different 
  groups of reader for ``Spectrum.read`` and ``Measurement.read``

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

- biologic readers now recognize both "Ece/V" and "<Ece>/V" as "raw_CE_potential".
  Resolves `Issue #110 <https://github.com/ixdat/ixdat/issues/110`_


Debugging
---------

general
^^^^^^^
- ``EC_MS`` is no longer a dependency
  Resolves `Issue #128 <https://github.com/ixdat/ixdat/issues/124>`_

measurement
^^^^^^^^^^^
- ``cut`` no longer crashes when one of the component measurements is empty.
  Resolves `Issue #93 <https://github.com/ixdat/ixdat/issues/93>`_

techniques
^^^^^^^^^^^

- ``MSInlet.gas_flux_calibration_curve()`` now works also when passing an
  axes_measurement 
  Resolves `Issue #94 <https://github.com/ixdat/ixdat/issues/94>`_

- ``CyclicVoltammogram.calibrate()`` now works, passing arguments on to a new ``ECCalibration``
  Resolves `Issue #111 <https://github.com/ixdat/ixdat/issues/111>`_

readers
^^^^^^^
- Zilien MS spectrum reader fixed.
  Resolves `Issue #117 <https://github.com/ixdat/ixdat/issues/117>`_

plotters
^^^^^^^^
- ``ECOpticalPlotter.plot_wavelengths_vs_potential()`` now returns a list of axes.
  Resolves `Issue #121 <https://github.com/ixdat/ixdat/issues/121>`_

exporters
^^^^^^^^^
- Fixed exporting and re-importing of ``ECOpticalMeasurment``s for new pandas version
  Resolves `Issue #124 <https://github.com/ixdat/ixdat/issues/124>`_

constants
^^^^^^^^^
- ``BOLTZMAN_CONSTANT`` renamed ``BOLTZMANN_CONSTANT``
  Resolves `Issue #125 <https://github.com/ixdat/ixdat/issues/125>`_
