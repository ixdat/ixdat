ixdat 0.2.2 (2022-04-08)
========================

API Changes
-----------

measurements
^^^^^^^^^^^^

- Ability to change the selector increment of a measurement using
  ``Measurement.rebuild_selector``. This returns a ``ValueSeries`` which counts cumulative
  changes in specified columns (common ones include "file_number", "cycle_number", etc).
  Especially useful for compound techniques with biologic potentiostats.

techniques.ms
^^^^^^^^^^^^^
- ``ECMSMeasurement.ecms_calibration_curve()`` no longer returns axes by default.
  Instead it only returns the calculated ``MSCalResult``. The axes on which the result
  is plotted can be returned by setting the keyword argument ``return_axes=True``

- ``MSMeasurement.grab_flux()`` can take as its first argument the name of a molecule
  for which the measurement has a calibration, but can also take a ``MSCalResult`` directly.
  As a result, MS plotting functions also accept ``MSCalResults`` in the requested
  ``mol_list``.  For example::

    cal_H2_M2 = my_ecms_meas.ecms_calibration(
        mol="H2", mass="M2", n_el=-2, tspan=[100, 200), tspan_bg=[0, 20]
    )
    t, n_dot = my_ms_meas.grab_flux(cal_H2_M2)
    my_ms_meas.plot(mol_list=[cal_H2_M2])

readers
^^^^^^^

- An ixdat native reader for Spectro Inlets - Zilien files has been
  implemented, as a replacement for leveraging the one in the legacy
  EC_MS project. While re-implementing it, support was also added for
  all auxiliary data series like MFC, PC, ion gauge and pirani values
  etc. Also, the Zilien tsv files can now be loaded as both EC-MS, MS
  and EC measurements, which solves some issues pertaining to
  plotting. To load a Zilien tsv file as something else than the
  default ``ECMSMeasurement``, do either:

    >>> ms_measurement = MSMeasurement.load(path_to_file="...")

  or:

    >>> ms_measurement = Measurement.load(path_to_file="...", technique="MS")

  This added functionality solves
  `Issue #31 <https://github.com/ixdat/ixdat/issues/31>`_

Bug Fixes
---------

- ``ECMSMeasurement.ecms_calibration_curve()`` now plots on the correct axes when given
  ``axes_measurement``.

- Adding calibration data to a measurement with ``calibrate()`` now clears the cache.
  This ensures that calibrated series are generated when needed rather than reusing cached
  uncalibrated series.

- ``Measurement.read(... reader="zilien", technique="MS")`` now returns an
  ``MSMeasurement``. Before, the user had to import ``MSMeasurement`` and use its read
  function or they'd get the reader complaining because it couldn't find the EC data.

- ``ECMeasurement["selector"]`` now correctly counts cumulative changes in file number,
  loop number, and cycle/step number. Previously there had been bugs of ignoring "cycle number"
  or "Ns" from biologic file sets if that variable wasn't included in all files, and
  including it even in file types where it is less useful as a counter (for example,
  Biologic seems to increment cycle number every time current changes sign during a
  chronoamperometric potential hold).



ixdat 0.2.1 (2022-04-01)
========================

API changes
-----------

techniques
^^^^^^^^^^

- ``ECMSCyclicVoltammogram.diff_with()`` raises a ``NotImplementedError`` rather than
  an obscure error message.

- ``ECMSPlotter.plot_vs_potential`` can accept ``color`` as a keyword argument. This
  colors both the U-J curve in the lower panel and all the signals in the top panel,
  so best to use together with a ``mass_list`` or ``mol_list`` of length 1.

exporters
^^^^^^^^^

- ``export()`` functions now take a ``time_step`` argument. It must be used together with
  a ``tspan``. It is effective in reducing the file size of exported .csv's. `PR #60 <https://github.com/ixdat/ixdat/pull/60>`_

- Renamed keyword argument: ``columns`` replaces ``v_list`` in exporters (``ECExporter`` in v0.2.0, ``ECMSExporter`` in 0.2.1)

  The old name is deprecated.

plotters
^^^^^^^^

- Added interactive range selection functionality to all plotters,
  `PR #61 <https://github.com/ixdat/ixdat/pull/61>`_. Left mouse click will set the
  left marker, right mouse click the right marker, and double clicking with any of the
  buttons will remove that marker. When both left and right markers are in place, the
  selected span will be printed out.

constants
^^^^^^^^^
- Dynamic viscosities are temperature-dependent, `PR #55 <https://github.com/ixdat/ixdat/pull/55>`_
  This enables accurate gas flux MS calibration (for the gases included) accross a range of temperatures.

Bug Fixes
---------

- ``ECMSExporter`` works as of 0.2.1 (it had been broken in 0.2.0).

  This solves `Issue #56 <https://github.com/ixdat/ixdat/issues/56>`_

- Other exporters were also sometimes prone to having extra ``aliases`` leading to
  circular lookups in the measurement when they were loaded. This has been fixed (as of 0.2.1)

- ``RE_vs_RHE=0`` counts as calibrated potential!



ixdat 0.2.0 (2022-03-25)
========================

API Changes
-----------

ixdat.measurement
^^^^^^^^^^^^^^^^^

- Generalization of multiple calibrations:

  ``Measurement.calibration`` is deprecated in favor of ``Measurement.calibration_list``.
  A dummy calibration property is put in with a getter and setter which raise
  ``DepreciationWarning``'s. Use, e.g.:

  - ``meas.calibration_list[0]`` (replaces ``meas.calibration``).

  - ``meas.add_calibration(my_cal)`` (replaces ``meas.calibration = my_cal``).

- New method ``Measurement.calibrate(*args, **kwargs)``

  ``meas.calibrate(...)`` is a shortcut for ``meas.add_calibration(calibration_class(...))``
  where calibration_class is found by looking up ``meas.techniques`` in
  ``ixdat.technqiues.CALIBRATION_CLASSES``. For esample, if ``meas`` is an
  ``ECMSMeasurement``, calibration_class is ``ECMSCalibration``.

  Together with the generalization of multiple calibrations, this enables very flexible
  calibration. All of the following code examples work.

  1. When measurements are appended or hyphenated, all their calibrations carry over.
  This example results in ``my_ecms_meas.calibration_list`` having an
  ``ECCalibration`` and an ``MSCalibration``, both of which are accessible to ``grab``
  and plotting functions.

    >>> my_ec_meas.calibrate(RE_vs_RHE=0.715, A_el=0.196)
    >>> my_ms_meas.calibrate(ms_cal_results=[my_H2_at_M2, my_O2_at_M32])
    >>> my_ecms_meas = my_ec_meas + my_ms_meas
    >>> my_ecms_meas.plot(mol_list=["H2", "O2"])

  2. You can calibrate one measurement multiple times. When two calibrations contain the
  same parameter, the last one added is used:

    >>> my_ecms_meas = my_ec_meas + my_ms_meas
    >>> my_ecms_meas.calibrate(RE_vs_RHE=0.715, A_el=0.196)
    >>> my_ecms_meas.calibrate(ms_cal_results=[my_H2_at_M2, my_O2_at_M32])
    >>> my_ecms_meas.calibrate(RE_vs_RHE=0.656)   # overshadows the first RE_vs_RHE
    >>> my_ecms_meas.plot(mol_list=["H2", "O2"])

  3. You can calibrate all at once.

    >>> my_ecms_meas = my_ec_meas + my_ms_meas
    >>> my_ecms_meas.calibrate(
    >>>     RE_vs_RHE=0.715, A_el=0.196, ms_cal_results=[my_H2_at_M2, my_O2_at_M32]
    >>> )   # note that all of these are keyword arguments to ECMSCalibration.
    >>> my_ecms_meas.plot(mol_list=["H2", "O2"])


ixdat.techniques
^^^^^^^^^^^^^^^^

- Renamed measurement technique class: ``CyclicVoltammogram`` (replaces ``CyclicVoltammagram``).
  The old name is deprecated.

- Renamed properties in ``ECMeasurement`` and inheriting classes:

  - ``U_name`` (replaces ``V_str``)
  - ``J_name`` (replaces ``J_str``)
  - ``E_name`` (replaces ``E_str``)
  - ``I_name`` (replaces ``I_str``)

  The old property names are deprecated.

- Lookups instead of properties

  - ``my_ec_meas["raw_potential"]`` replaces ``ECMeasurement.raw_potential``
  - ``my_ec_meas["raw_current"]`` replaces ``ECMeasurement.raw_current``
  - ``my_cv["scan_rate"]`` replaces ``CyclicVoltammogram.scan_rate``

  The old properties are deprecated.

- Renamed keyword in ``MSMeasurement.grab_flux()`` and related methods:

  - ``remove_background`` (replaces ``removebackground``)

  The old argument name is deprecated

- ``MSCalResult.name`` is by default set to ``{mol}@{mass}``, e.g. "H2@M2" instead of
  ``{mol}_{mass}``, e.g. "H2_M2". However, you can
  use just the mol name in ``grab`` and plotting functions if a measurement has an ``MSCalResult``
  with the molecule in the ``ms_cal_list`` of one of the calibrations in its ``calibration_list``:

    >>> cal = my_ms_inlet.gas_flux_calibration(
    >>>     my_ms_measurement, mol="H2", mass="M2", tspan=[0, 20]
    >>> )
    >>> cal.name
    'H2@M2'
    >>> cal.mol
    'H2'
    >>> my_ms_measurement.calibrate(ms_cal_results=[cal])
    >>> meas.grab("n_dot_H2")
    numpy.Array([...]), numpy.Array([...])
    >>> meas.plot(mol_list=["H2"])


ixdat.plotters
^^^^^^^^^^^^^^

- Axes order: Two-panel figures with shared x-axis always return a list of axes in the order
  ``axes=[top_left, bottom_left]`` in the case of one y-axis of each panel, or
  ``[top_left, bottom_left, top_right, bottom_right]`` in the case of an additional y-axis
  in either panel. Note that ``axes[2]`` or ``axes[3]`` might be ``None``. For example,
  ``axes = ECMSPlotter.plot_measurement()`` by default results in ``axes[0]`` being the
  MS data, ``axes[1]`` electrode potential, ``axes[2]=None`` and ``axes[3]`` being the
  electrode current.

- Renamed keyword arguments in EC, EC-MS, and SEC plotting functions:

  - ``U_name`` (replaces ``V_str``)
  - ``J_name`` (replaces ``J_str``)
  - ``U_color`` (replaces ``V_color``)

  The old keyword argument names are deprecated.

- Renamed keyword argument in MS and EC-MS plotting functions:

  - ``remove_background`` (replaces ``removebackground``)

  The old argument name is deprecated

ixdat.readers
^^^^^^^^^^^^^
- ixdat-exported .csv back compatability

  The "ixdat" reader is no longer automatically compatible with .csv files exported by
  the ``IxdatCSVExporter`` of ixdat v0.1.x. You can, however, get ixdat v0.2.0 to read
  your v0.1.x .csv's by giving it a little help, in the form of aliases. For example:

     >>> meas = Measurement.read_url(
     >>>     "https://raw.githubusercontent.com/ixdat/tutorials/"
     >>>     + "ixdat_v0p1/loading_appending_and_saving/co_strip.csv",
     >>>     reader="ixdat",
     >>>     aliases={
     >>>         "t": ["time/s"],
     >>>         "raw_current": ["raw current / [mA]"],
     >>>         "raw_potential": ["raw potential / [V]"]
     >>>     }
     >>> )

  Can read `this file <https://raw.githubusercontent.com/ixdat/tutorials/ixdat_v0p1/loading_appending_and_saving/co_strip.csv>`_
  into ixdat.

  To know which aliases to use, you should check the file and the ``essential_series_names``
  of the technique class. For example:

    >>> form ixdat.techniques import ECMeasurement
    >>> ECMeasurement.essential_series_names
    {'t', 'raw_potential', 'raw_current'}

  For even earlier .csv files (exported by ixdat version <0.1.5), you will also need to
  specify the technique.

  Starting in ixdat v0.2.0, ixdat-exported files will have the ixdat version in the header.

Bug Fixes
---------

Ixdat 0.2.0 has a number of deep improvements in data series handling, which help
fix the following bugs:


- Reliable construction of `selector` and `cycle` counters. This is done via the
  ``selector_name`` and ``_series_constructors`` class attributes, which can be
  customized for every ``Measurement`` subclass.

  This solves `#14 <https://github.com/ixdat/ixdat/issues/14>`_

- Carrying calibrations over through measurement combination (hyphenation or appending
  with the `+` operator) or transformation (through copying and technique-changing
  functions like `ECMeasurement.as_cv()`. This is done via replacing the use of a single
  ``Measurement.calibration`` with a ``calibration_list`` which can be appended,
  treating the ``calibration_list`` as a measurement's linked objects, and a
  ``MemoryBackend`` which stores such linked objects while the main object is being
  copied via its dictionary representation.

  This solves `#20 <https://github.com/ixdat/ixdat/issues/20>`_

  This also solves `#22 <https://github.com/ixdat/ixdat/issues/22>`_

  This also solves `#29 <https://github.com/ixdat/ixdat/issues/29>`_





