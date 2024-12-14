

ixdat 0.3.0 (2024-09-25)
========================

API Changes
-----------
This release marks an important change in ixdat's organization resulting from the
disscusions under the `calculators label. <https://github.com/ixdat/ixdat/issues?q=is%3Aissue+label%3Acalculators>`_
and `PR #182 <https://github.com/ixdat/ixdat/pull/182>`_.

- Nearly all calculations for processing and analyzing raw data are moved into a
  new type of Saveable class, ``Calculators``. These take over from the
  ``Calibration`` classes of ixdat 0.1.x and 0.2.x but the name is more general in
  order to accomodate e.g. ``Filters`` and ``Backgrounds`` which have been long
  anticipated (see `Issue #147 <https://github.com/ixdat/ixdat/issues/147>`_), as well as what were formerly known as
  ``series_constructor``s

- To use language that we've employed for quantitative mass spectrometry,
  ``Calibration`` classes are in general for both *calibration*, whereby values
  needed to perform a certain type of data analysis, like sensitivity factors,
  are extracted from one measurement (representing a calibration experiment);
  and *quantification* where these values are applied to analyze data from another
  measurement (representing the experiment of interest).
  - Calibration is generally done by constructor class methods, e.g.::

    calc = MSCalibration.gas_flux_calibration(measurement_1, ...)

  where ``MSCalibration`` is a class that inherits from ``Calculator``.

  - Quantification is generally done by normal object methods, e.g.::

      n_dot_H2 = calc.calc_flux(measurement_2, mol="H2", ...)

    All ``Calculator`` classes have a special quantification method called
    ``calculate_series(measurement, series_name)`` which, if possible,
    returns a ``ValueSeries`` for ``series_name`` derived from data in
    ``measurement`` and the ``Calculator``'s own internal data.

  Ideas behind this are described in `Issue #164 <https://github.com/ixdat/ixdat/issues/164>`_.

- All of the ``series_name``s that can be used in a ``Calculator``'s
  ``calculate_series`` method are included in the set
  ``Calculator.available_series_names``. In general, this is a dynamic property, as
  it may depend on which internal data the calculator object has. For example,
  an ``MSCalculator`` object which has a sensitivity factor for "H2" at "M2" will
  include ``"n_dot_H2"`` in its ``available_series_names`` whereas one which
  instead only has a sensitivity factor for "O2" at "M32" will not.

- A ``Calculator`` object (``calc``) can be attached to a ``Measurement`` object
  (``meas``) by calling ``meas.add_calculator(calc)``. This makes it possible to
  directly ``grab`` or look up any of the series in that calculator's
  ``available_series_names``. A few challenges arise here with regards to calculator
  chaining, (an example is when we want to grab a molecular flux that is both quantified from
  the raw mass signals and also background subtracted, or not). Thise issues
  are described in `issue #183 <https://github.com/ixdat/ixdat/issues/183>`_.
  The solution implemented is described in the first comment to that issue, and
  demonstrated in **demo_calculators.py**.

- A ``Measurement`` class can have a ``default_calculator`` class. An object of this
  class is initiated (using ``__init__``, not a constructor method) and automatically
  attached to the measurement object (``meas``) by calling the method
  ``meas.calibrate(...)``. A ``Measurement`` class also has a list of
  ``built_in_calculator_types`` for which a calculator is initiated with default
  arguments and appended to ``meas`` upon its initation, and can have a list of
  ``background_calculator_types`` which are ignored when the argument
  ``remove_background=False`` is passed to ``grab`` and methods that use ``grab``.

Calculators
^^^^^^^^^^^

The code for ``Calculator`` classes is in the modules in **src/ixdat/calculators/**.
As of now, the following are included:

- ``indexer.Indexer`` is in the built-in
  calculators of the ``Measurement`` base class. Every ``Measurement`` object thus
  has an indexer. It is used for easily choosing a subset of the measurement.
  - available_series_names: "file_number" and "selector". ``calculate_series``
    uses the following methods, respectively, for those two series names:
  - ``_build_file_number_series``, formerly the series constructor
    ``Measurement._build_file_number_series``
  - ``_build_selector``, formerly the series constructor
    ``Measurement._build_selector_series``

- ``ec_calculators.ECCalibration`` is the default_calculator of ``ECMeasurement``.
  Otherwise, it is basically just moved from **techniques/ec.py**.
  - available_series_names: "potential" and "current"
  - "potential" is the raw_potential (i) adjusted to the RHE scale if ``RE_vs_RHE``
    is included and (ii) corrected for ohmic drop if ``R_Ohm`` is included
  - "current" is the raw_current (i) normalized to electrode area if ``A_el``
    is included.
  - Note that looking up "potential" or "current" in an ``ECMeasurement`` without
    an ``ECCalibration`` returns "raw_potential" or "raw_current", respectively, as
    in ixdat v0.2.x

- ``ec_calculators.ScanRateCalculator`` is a built-in calculator for ``CyclicVoltammogram``,
  previously implemented as a function and methods of ``CyclicVoltammogram``
  - available_series_names: "scan_rate"

- ``ms_calculators.MSCalibration`` is, like in ixdat 0.2.x, a wrapper around a list
  of sensitivity factors (themselves saveable as ``ms_calculators.MSCalResult``s).
  - available_series_names: "n_dot_{mol}" for all the mol in its ``mol_list``
  Constructor methods:
  - ``gas_flux_calibration``, formerly a method of ``MSMeasurement``.
  - ``gas_flux_calibration_curve``, formerly a method of ``MSMeasurement``.

- ``ms_calculators.MSBackgroundSet`` is a new calculator that replaces the poor
  incomplete implementation of backgrounds previously in ``MSCalibration``. The
  structure is similar to ``MSCalibration`` in that a ``MSBackgroundSet`` contains
  a set of saveable ``MSBackground`` objects, each for a single m/z. So far, only
  one type of background is implemented, the ``MSConstantBackground``.
  - available_series_names: mass_list.
  - Note that The available series names have the same names as the corresponding
    raw data series before background subtraction. To get the raw series, grab or
    look up ``f"{mass}-raw"`` or ``grab`` with ``remove_background=False``.

- ``ecms_calculators.ECMSCalibration`` is not a real calculator in the sense that
  it doesn't do *quantification*. Instead, it does *calibration*, and its calibration
  methods, listed below, all return ``MSCalibration`` objects. This is consistent
  with the fact that an EC-MS calibration experiment can be used to obtain sensitivity
  factors for a setup which is then used without electrochemistry (e.g. for thermal
  catalysis measurements). The calibration methods are
  - ``ecms_calibration``, formerly a method of ``ECMSMeasurement``
  - ``ecms_calibration_curve``, also formerly a method of ``ECMSMeasurement``

- ``ecms_calculators.ECMSImpulseResponse``, moved from **deconvolution.py**, is the
  deconvolution calculator. An ``ECMSImpulseResponse`` object (``imp_resp``)
  describes the response of one ``mol``. It's demonstrated in
  **deconvolution_demo.py**.
  - available_series_names: "n_dot_{mol}-deconvoluted". Getting this from an
    ``MSMeasurement`` with ``imp_resp`` attached requires that there is also another
    calculator which provides "n_dot_{mol}".
  Constructor methods:
  - ``from_measurement`` takes the shape of the impulse response from a measurement
    representing an impulse experiment, i.e. one where a short burst of product (e.g.
    "H2" form hydrogen evolution) is produced at the electrode and its mass signal,
    after broadening by mass transport between the electrode and the inlet, recorded.
  - ``from_model`` calculates the shape of the impulse response according to a mass
    transport model.

- ``xrf_calculators.TRXRFCalculator`` is a built-in calculator of ``TRXRFMeasurement``
  which replaces a simple series_constructor.
  - available_series_names: "FF_over_I0"


One ``Calculator``, the ``siqCalculator`` for advanced MS and EC-MS calibration,
is implemented as a plugin. At present, this is in **src/ixdat/plugins/siq_plugin.py**.

- ``siqCalculator`` implements all the calibration and quantification methods that
  make use of the external ``spectro_inlets_quantification`` package. The calibration
  methods were previously methods of ``MSMeasurement`` and ``ECMSMeasurement`` prefixed
  "siq_", and the quantification methods were accessed through an overloading of
  ``MSMeasurement.grab()``.

  Constructor methods:
  - ``gas_flux_calibration``, formerly ``MSMeasurement.siq_gas_flux_calibration``.
  - ``gas_flux_calibration_curve``, formerly ``MSMeasurement.siq_gas_flux_calibration_curve``.
  - ``multicomp_gas_flux_calibration``, formerly ``MSMeasurement.siq_multicomp_gas_flux_calibration``
  - ``ecms_calibration``, formerly ``ECMSMeasurement.siq_ecms_calibration``
  - ``ecms_calibration_curve``, formerly ``ECMSMeasurement.siq_ecms_calibration_curve``
  Useage:
  - A ``siqCalculator`` object (``siqcalc``) inherits from
    ``spectro_inlets_quantification.Calibration`` (as well as ``ixdat.Calculator``),
    which implements addition, visualization, and sensitivity factor prediction
  - Before use, a ``siqCalculator`` object must be given a ``mol_list`` and ``mass_list``,
    which are used to define a ``SensitivityMatrix``, as well as a ``carrier`` gas
    (typically "He") as needed by ``siq.Quantifier``. These parameters are given
    by the method ``siqcalc.set_quantifier(mol_list=..., mass_list=..., carrier=...)``.
  - available_series_names: "n_dot_{mol}" for all the mol in its ``mol_list``, but
    only once ``mol_ist`` has been given via the ``set_quantifier`` method

  Use of ``siqCalculator`` is demonstrated in **demo_siq_integration.py** and
  **demo_ecms_calibration_curve.py**.



ixdat 0.2.13 (2024-09-03)
=========================

Debugging
---------
- Fixed timestamp form in ``QexafsDATReader`` to correctly parse timezone all year.


API changes
-----------

- Time-resolved x-ray flouresence (``technique = "TRXRF"``) implemented in `PR #168 <https://github.com/ixdat/ixdat/pull/168>`_:

  - ``B18TRXRFReader`` (reader="b18_trxrf") implemented for reading TRXRF data from the Diamond lightsource beamline B18TRXRFReader

  - ``TRXRFMeasurement`` with a series constructor method for the value series of interest, "FF_over_I0", and ``TRXRFPlotter`` for plotting the TRXRF data.

  - Hyphenation of TRXRF with EC (``technique = "EC-TRXRF"``) implemented (syntax: ``ec_txrf = ec + trxrf``) in ``ECTRXRFMeasurement`` and ``ECTRXRFPlotter``


- Deconvolution module based on Krempl et al. 2021 https://pubs.acs.org/doi/abs/10.1021/acs.analchem.1c00110
  is revived. ``ECMSImpulseResponse`` is a class for calculating an impulse response
  as input for deconvolution. It can generated either from a measured impulse response using class method
  ``.from_measurement()`` or from mass transport parameters using class method ``.from_parameters``.
  Several methods of ECMSMeasurement class use this new class: ``grab_deconvoluted_signal()`` allows to grab
  a an tuple of time and value arrays (similar to other ``grab()`` methods). ``deconvolute_for_tspans()`` loops
  through a number of tspans for which to deconvolute data with options to plot and export the original + decon-
  voluted data. For examples see deconvolution_demo.py in development_scripts




ixdat 0.2.12 (2024-05-17)
=========================

API changes
-----------

techniques
^^^^^^^^^^

- Added ``faradaic_efficiency`` to correct the molecule flux as calculated from EC if not 100% FE,
	also works for isotopes

- Added ``integrate_flux()`` method to techniques\ms.py (analogous to existing integrate_signal(),
  allows for integration of the calibrated signal)


plotters
^^^^^^^^

- Added default colors for "H2@M2", "H2@M3", "H2@M4" to ms_plotter.py


ixdat 0.2.11 (2024-04-16)
========================

debugging
---------
- Fixed three bugs in 0.2.10 on spectrum series heat plotting:
  - The x axis limit is now such that the full duration
    of all the spectra is included when ``spectrum_series.continous`` is ``False`` (not just the last one).
  - A universal ``vmin`` and ``vmax`` are determined for all the spectra to be plotted. This avoids effective
    normalization of each individual spectrum when ``spectrum_series.continuous`` is ``False`` (which made it
    impossible to see changes in overall intensity).
  - The ``vs`` variable is now handled correctly when ``continuous`` is ``False``, enabling plotting vs potential
    of non-continuous spectrum series data.

API changes
-----------

- ``OpusFTIRReader`` and ``MsrhSECDecayReader`` both set the ``continuous`` attribute of the ``SpectrumSeries``
  object to ``True``, since these are readers for fast and continuous-scanning instruments.

- The ``continuous`` attribute is now exported and read from ixdat-exported .csv's.

- The ``continuous`` property of ``SpectroMeasurement`` now has a setter. To change it, e.g.::

    ftir = Spectrum.read(...)
    ec = Measurement.read(...)
    ecftir = ec + ftir
    ecftir.continous = False
    ecftir.plot_measurement()

ixdat 0.2.10 (2024-04-15)
========================

debugging
---------
- Fixed bug in 0.2.9 on exporter initiation that had caused a crash on initation of some Spectrum
  and SpectroMeasurement objects due to removed ``delim`` argument.

- Fixed bug in 0.2.9 on spectrum reading that had caused a crash when making SpectrumSeries objects
  from auxiliary files representing simple Spectrum objects.


ixdat 0.2.9 (2024-04-12)
========================

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
  With `PR #134 <https://github.com/ixdat/ixdat/pull/134>`_

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
  To leave out the mass scan data, include the argument ``include_mass_scans=False``
  in the call to ``read``. To leave out EC data but include the mass scans, include
  the argument ``technique="MS-MS_spectra"`` in the call to ``read``.
  This finishes `Issue #117 <https://github.com/ixdat/ixdat/issues/117`_

- If a series name is present in the raw data *and* in in a measurement's ``aliases``,
  the raw data series matching the name and the aliased series are appended. (Before,
  only the raw data series matching the name would be returned.)

- Added "opus_ftir" reader for text-exported files from Opus FTIR spectrumeter.

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

techniques
^^^^^^^^^^

- Added FTIR and EC-FTIR. The latter for now just inherits from SpectroECMeasurement

- ECOpticalMeasurement.get_spectrum() now has an option not to interpolate.
  Passing the argument `interpolate=False` gets it to choose nearest spectrum instead.

- Indexing a ``SpectroMeasurement`` with an integer returns a ``Spectrum``.
  For example, ``zilien_meas_with_spectra[0].plot()``  plots the first mass scan

- ``SpectrumSeries`` and inheriting classes now have a ``cut()`` method which, like
  ``Measurement.cut``, returns an object containing the subset of the data falling in
  a specified timespan.

- ``SpectrumSeries`` and ``SpectroMeasurement`` objects can be added to each other,
  given that their ``xseries`` have the same shape and the technique resulting from
  ``get_combined_technique(meas_or_spec_1.technique, meas_or_spec_2.technique)``
  is present in ``TECHNIQUE_CLASSES``. For example:
      ms_with_spectra (technique="MS") + ms_without_spectra (technique="MS-MS_spectra")
          --> ms_with_spectra (technique="MS-MS_spectra")

  - `Issue #158 <https://github.com/ixdat/ixdat/issues/158>`_
  - `PR #166 <https://github.com/ixdat/ixdat/pull/166>`_

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

- Added a ``plot_stacked_spectra`` method to ``SpectrumSeriesPlotter``, ``SpectroMeasurementPlotter``,
  and ``SECPlotter``. This is the default plotting method for FTIR.

exporters
^^^^^^^^^
- ``Spectrum`` and ``SpectrumSeries`` objects now have working ``export`` functions, as does
  ``MSSpectroMeasurement``. In each case, the exported file can be read (``reader="ixdat"``)

- Default delimiter is now just ",", not "\t,", as the latter doesn't open well in excel.

- ``time_step`` argument to ``export`` can now be used without specifying a ``tspan``

  - `Issue #153 <https://github.com/ixdat/ixdat/issues/153>`_
  - `PR #166 <https://github.com/ixdat/ixdat/pull/166>`_

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


ixdat 0.2.8 (2023-12-05)
========================

Documentation
-------------
Jupyter notebook tutorials are now compiled to .rst with nbsphinx.
This solves `Issue #115 <https://github.com/ixdat/ixdat/issues/115>`_


API changes
-----------

techniques
^^^^^^^^^^
- Native and ``spectro_inlets_quantification`` (siq) calibration objects can be interconverted:
  - ``cal.to_siq()`` for ixdat ``MSCalResult`` --> siq ``CalPoint``
  - ``MSCalResult.from_siq(cal)`` for siq ``CalPoint`` --> ixdat ``MSCalResult``
  - ``cals.to_siq()`` for ixdat ``MSCalibration`` --> siq ``Calibration``
  - ``MSCalibration.from_siq(cals)`` for siq ``Calibration`` --> ixdat ``MSCalibration``

- ``MSInlet.gas_flux_calibration_curve()`` has the additional option of passing
  a boolean ``axes_measurement_raw``. Set to True if the axes passed to
  ``axes_measurement`` are raw signal data (i.e. not background subtracted)
  Mentioned in `Issue #94 <https://github.com/ixdat/ixdat/issues/94>`_

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
- biologic readers now recognize both "Ece/V" and "<Ece>/V" as "raw_CE_potential".
  Resolves `Issue #110 <https://github.com/ixdat/ixdat/issues/110`_

- Native MS calibration methods are moved from ``MSInlet`` to ``MSMeasurement``.
  Meanwhile, those using ``spectro_inlets_quantification`` are prefixed with ``siq_``,
  e.g. use ``ms_meas.gas_flux_calibration(...)`` for native ixdat MS calibration and
  ``ms_meas.siq_gas_flux_calibration(...)`` to use the more powerful siq equivalents.
  Resolves `Issue #122 <https://github.com/ixdat/ixdat/issues/122>`_

- ``spectro_inlets_quantification`` is consistently abbreviated ``siq``.


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


ixdat 0.2.7 (2023-12-05)
========================


ixdat 0.2.6 (2023-10-04)
========================

API changes
-----------

Readers
^^^^^^^
- A function ``fix_WE_potential`` has been added to the ``biologic`` module. This
  corrects the mistake in some .mpt files that are exported with zeros for "<Ewe>/V".

Techniques
^^^^^^^^^^
- MS measurements now have access to the powerful ``spectro_inlets_quantification`` (siq)
  package as a plugin. See `siq's documentation <https://spectroinlets-spectro-inlets-quantification.readthedocs-hosted.com/en/v1.1/theory/quant_ms.html>`_
  for details.

  To activate the plugin, use::

    import ixdat
    ixdat.plugins.activate_si_quant()

  after activation, the following calibration and quantification methods will use siq:
  - ``MSMeasurement.gas_flux_calibration`` now returns a siq ``CalPoint``. Note that this
    method is only available with siq activated - native ixdat gas flux calibration
    is only available as a method of ``MSInlet``. Otherwise they work basically the same.
  - ``MSMeasurement.multicomp_gas_flux_calibration`` returns a ``Calibration``. Note that it
    solves a matrix equation to deconvolute overlapping peaks in a multi-analyte
    calibration gas.
  - ``ECMSMeasurement.ecms_calibration`` and ``ECMSMeasuremnt.ecms_calibration_curve`` each
    return an object of ``CalPoint``.
  - ``MSMeasurement.set_quantifier`` is used to attache a siq quantifier to the measurement.
    Usage: ``ms.set_quantifier(calibration=my_siq_calibration)``
  - ``MSMeasurement.grab_fluxes`` is a new method which uses the attached quantifier to
    quantify the fluxes of all molecules covered therein. In general it solves a matrix
    equation to deconvolute overlapping signals.
  - ``MSMeasurement.grab_flux`` will, when siq is activated, run ``MSMeasurement.grab_fluxes``
    and return the requested flux vector.

  SIQ comes with data for a small batch of molecules. To supplement this collection of
  yaml-formatted molecule data, place new files in the QUANT_DIRECTORY. This can be set
  as follows (change to the folder where you have your data)::

    ixdat.plugins.si_quant.QUANT_DIRECTORY = "~/projects/batteries/quantification_data"


debugging
---------

readers
^^^^^^^

- ``cinfdata_db`` is a new direct db reader for cinfdatabase at DTU SurfCat

plotters
^^^^^^^^

- ``MSPlotter.plot_measurment`` now create a new axis if no initial axis is given
  before initiating right axes in the case of dual plotting on left and right axes.
  Fixes `Issue #97 <https://github.com/ixdat/ixdat/issues/97>`_


- ``SpectrumSeriesPlotter.heat_plot`` now accept max_threshold and min_threshold and
   scanning_mask to include or exclude specific values from scanning variable

- ``SpectroMSPlotter`` new plotter for ``SpectroMSMeasurment`` now create a new axis if no initial axis is given
  before initiating right axes in the case of dual plotting on left and right axes.

techniques
^^^^^^^^^^

- ``ReactorMeasurement`` class, technique="reactor", with a ``TPMSPlotter``. This
  technique is analogous to EC-MS with temperature replacing potential and
  pressure replacing current.

- ``SpectroMSMeasurement`` class set ``SpectroMSPlotter`` as default plotter

dev
^^^
- Renamed development scripts which are not software tests "demo" instead of "test".

- Skip py36 because github is having a problem building it. See, for example, here:
https://github.com/ixdat/ixdat/actions/runs/3876991446/jobs/6611480640#step:3:7

- Do black test before software tests in github CI to save time


ixdat 0.2.5 (2022-11-11)
========================

Quant integration
-----------------

general
^^^^^^^
An object called ``plugins`` can be imported from the ``config`` module. This gives access
to variables controlling plugin behaviour.``plugins.USE_QUANT`` determines whether an
external quantification package is used for quantification by ``MSMeasurement``.
``plugins.QUANT_DIRECTORY`` determines where that external quantification package looks
for reference data.

A module ``quant_patch`` is included with functions that act on objects of the external
quantification package, with the hope that the functionality migrates to that package in
the near future.

techniques
^^^^^^^^^^
``MSMeasurement`` has new methods which can be used if ``plugins.USE_QUANT`` is True:

- ``MSMeasurement.gas_flux_calibration`` for sensitivity factor determination by flowing
  a pure gas through an inlet chip
- ``MSMeasurement.multicomp_gas_flux_calibration`` for sensitivity matrix determination
  by flowing a gas with multiple analyte components. This uses reference spectra and
  results in sensitivity factors that are as inaccurate as the reference spectra. Its use
  is therefore discouraged (but sometimes necessary).
- ``MSMeasurement.set_quantifier`` sets the *quantifier*, which then determines how
  ``MSMeasurement.grab_flux`` and ``MSMeasurement.grab_fluxes`` calculate fluxes of
  molecules to the vacuum chamber.
- ``MSMeasurement.grab_fluxes`` uses the measurement's quantifier to calculate the fluxes
  of all the molecules in the quantifier's ``mol_list`` with the signals at all the
  masses in the quantifier's ``mass_list`` as inputs. It takes the tspan and background
  arguments familiar in ixdat from other ``grab`` methods.

The workings of ``MSMeasurement.grab_flux`` are changed if ``plugins.USE_QUANT`` is True.
In that case, it invokes the *quantifier* via ``MSMeasurement.grab_fluxes`` and returns
just the flux of the requested molecule.

New Zilien reader
-----------------

A new Zilien reader for the SpectroInlets' new Zilien dataset file format
version. The new dataset version is able to integrate the Biologic EC-lab
dataset during a measurement. Such dataset contains a new series with a header
name ``EC-lab`` and two meta columns ``experiment_number`` and
``technique_number`` in the series. The new reader is using the columns during
the process of creating Ixdat series objects. The objects, created from the
Zilien dataset, match exactly the objects created from the Biologic MPT files,
like when read one by one, **except for** the timestamps from the Biologic
series.  The Biologic series timestamps are incremented by a time offset when
the Biologic EC-lab measurement was triggered.

E.g. when you start a Zilien measurement and then trigger an EC-lab measurement
after five seconds, the timestamps in the series of the Biologic dataset part
will be higher by five, comparing to the timestamps in the MPT files.


Debugging
---------

- xrdml reader can now import files where the data is labeled "counts" rather than
  "intensities", as the text exports from the Royce Institute XRD


ixdat 0.2.4 (2022-07-25)
========================

API changes
-----------
techniques
^^^^^^^^^^^
- Improved docstring for ``ECMSMeasurement.ecms_calibration_curve()`` to include the new additions from previous release.

- Added MSInlet ``gas_flux_calibration_curve`` to enable multiple point calibration using calculated gas flux
  either with different concentrations in carrier gas or at different inlet pressures. Note, concentration needs to be given in ppm, as the flux calculation uses various constants from the carrier gas molecule instead of a mixture, which will lead to significant inaccuracy for high concentrations.

- ``Measurement.select`` is now even more versatile. A user can specify a ``selector_name``
  for args to work on. This enables selection based on columns with a space in them, like
  "cycle number".
  Resolves `Issue #77 <https://github.com/ixdat/ixdat/issues/77>`_

readers
^^^^^^^
- ``Mesurement.read_set`` can now take a folder as the first argument, in which case  all
  files in that folder (with the specified suffix) are appended.
  Resolves `Issue #88 <https://github.com/ixdat/ixdat/issues/88>`_

- ``Measurement.read_set`` now also raises a ``ReadError`` rather than returning ``None`` in
  the case of no matching files.

- ``Measurement.read`` (and by extension ``Measurement.read_set``) can now be called
  without a specified reader for certain known file types. To see which file types, use::

    from ixdat.readers.reading_tools import DEFAULT_READER_NAMES
    print(DEFAULT_READER_NAMES)

plotters
^^^^^^^^

- ``ECMSPlotter.plot_measurment`` and ``plot_potential`` can take string options "all",
  "ms", and "ec" to specify the time span of the data to plot.
  Resolves `Issue #91 <https://github.com/ixdat/ixdat/issues/91>`_

Debugging
---------

readers
^^^^^^^
- The biologic reader now checks for "Ns" and "cycle number" rather than assuming it
  knows which EC techniques have which of these selector-defining series.
  Resolves `Issue #87 <https://github.com/ixdat/ixdat/issues/87>`_

techniques
^^^^^^^^^^^
- ``Measurement.cut`` now skips empty component measurements rather than crashing.
  Resolves `Issue #93 <https://github.com/ixdat/ixdat/issues/93>`_

- ``MSMeasurement.reset_bg`` works again! It now adds a new calibration with bg=0 for
  masses that had previously had a bg set.
  Resolves `Issue #82 <https://github.com/ixdat/ixdat/issues/82>`_

- ``_get_tspan_list`` in ``ECMSMeasurement`` now defaults ``t_steady_pulse`` to ``None``
  instead of ``0``, which simplifies explanation in docstring and makes it more clear what
  it does (i.e. if now a pulse time of 0 is given it will actually use 0s instead of the
  entire pulse)


Other changes
-------------

- Tests for EC-lab parser using a dataset with multiple techniques and a dataset with looping techniques.


ixdat 0.2.3 (2022-06-10)
========================

API Changes
-----------

measurement
^^^^^^^^^^^
- Added example usage to docstring of ``Measurement.select_values``

readers
^^^^^^^
- Added an ``XRDMLReader`` (reader="xrdml") for xml-formatted XRD spectrum files from,
  for example, Empyrion's software. Usage::

    from ixdat import Spectrum
    my_xrd = Spectrum.read("my_file.xrdml", reader="xrdml")
    my_xrd.plot()

- Added an ``AvantageAVGReader`` (reader="avantage") for avantage-exported spectrum files from,
  for example, Thermo Fisher's K-alpha equipment. Usage::

    from ixdat import Spectrum
    my_xps = Spectrum.read("my_file.avg", reader="avantage")
    my_xps.plot()

- Added a ``QEXAFSDATReader`` (reader="qexafs") for .dat files exported by Diamond
  Synchrotron's B18-Core beamline. These data files include X-Ray Absorption across
  multiple detector elements, X-ray fluorescence, and diagnostic information such as primary
  beam intensity. All the data is read together as a new object called a ``MultiSpectrum``.
  If the ``read`` function is called with ``technique="XAS"``, a single spectrum with the
  XAS data (the "FFIO" column) is returned::

    from ixdat import Spectrum
    multi_spectrum = Spectrum.read("my_file.dat", reader="qexafs")
    my_xas = multi_spectrum["QexafsFFI0"]   # index retrieves a Spectrum from a MultiSpectrum
    # OR
    my_xas = Spectrum.read("my_file.dat", reader="qexafs", technique="XAS")
    my_xas.plot()

- A new reading method ``Spectrum.read_set()`` imports and appends multiple spectrum
  files. It is thus similar to ``Measurement.read_set()``. The appended spectra are
  returned as a ``SpectrumSeries``, which has a heat plot as its default ``plot()``.

- Both ``Measurement.read_set()`` and the new ``Spectrum.read_set()`` allow the user
  to specify a string that appears in the middle of the names of the files to be read
  and appended, rather than just the beginning::

    from ixdat import Spectrum
    spectrum_series = Spectrum.read_set(
       part="data/XAS_spectrum", suffix=".dat", technique="XAS"
    )  # reads .dat files in the directory "data" with "XAS_spectrum" in their name
    spectrum_series.plot()  # heat plot with time on x-axis and energy on y-axis

techniques
^^^^^^^^^^

- ``CyclicVoltammogram`` now has a ``plot_cycles`` function. This plots all
  the cycles in the cv color coded by cycle number, with a scale bar.

- Refactoring of ``SpectroECMeasurement`` in order to generalize aspects of combining
  time-resolved measurements with spectra taken simultaneously. This results in the
  following classes:

  - ``SpectroMeasurement(Measurement)`` for any time-resolved measurement data with spectra.
  - ``SpectroECMeasurement(SpectroMeasurement, ECMeasurement)``, ``technique="SEC"``
    for an EC measurement with spectra. By default plots the EC data on a lower panel and
    the spectral data as a heat plot on the upper panel.
  - ``ECXASMeasurement(SpectroECMeasurement)``, ``technique="EC-XAS"`` for an EC measurement
    with XAS spectra.
  - ``ECOpticalMeasurement(SpectroECMeasurement)``, ``technique=EC-Optical`` for an EC
    measurement with optical spectroscopy. A reference spectrum is required, and can also
    be set using a reference time or reference potential after the data is loaded. By
    default, optical density = **-log(y/y0)**, with ``y=spectrum_series.y`` and
    ``y0=reference_spectrum.y``, is plotted instead of raw data.

  Before the refactor, ``ECOpticalMeasurement`` had been called ``SpectroECMeasurement``.
  This is all discussed in `PR #73 <https://github.com/ixdat/ixdat/pull/73>`_

- Addition of a ``Measurement`` and a ``SpectrumSeries`` gives a ``SpectroMeasurement``
  or a subclass thereof determined by hyphenating the technique. For example::

    from ixdat import Measurement, Spectrum

    my_ec = Measurement.read("my_ec_data.mpt", reader="biologic")
    print(my_ec, ",", my_ec.technique)
    # ECMeasurement(...) , EC

    my_xas = Spectrum.read_set("data/my_xas", suffix=".dat", technique="XAS")
    print(my_xas, ",", my_ec.technique)
    # SpectrumSeries(...) , XAS

    my_ec_xas = my_ec + my_xas
    print(my_ec_xas, ",", my_ec_xas.technique)
    # ECXASMeasurement(...) , EC-XAS

    my_ec_xas.plot()  # XAS data heat plot in top panel and EC data in bottom panel.

- ``ECMSMeasurement.ecms_calibration_curve`` now supports data specification using a
  a selector. To do so, specify the section to use as numbers in the argument ``selector_list``,
  the counter defining the sections (if different from the default selector) in ``selector_name``,
  and the length of the steady-state period at the end of the pulse in ``t_steady_pulse``.
  This can be much more convenient than manually specifying ``tspans``.
  Implemented in #74.

plotters
^^^^^^^^

- Plotters for spectroelectrochemistry were refactored in a similar way to the technique
  classes (see above, under "techniques"). All ``SpectroECMeasurements`` come with plotters
  with ``plot_measurement()`` (heat plot and EC data vs time) and ``plot_vs_potential()``
  (heat plot and EC data vs potential). For ``ECOpticalMeasurements``, these plot
  optical density rather than raw data on the heat plots.

exporters
^^^^^^^^^
- The main .csv file exported from an ``ECOpticalMeasurement`` refers to its auxilliary
  files as::

    'spectrum_series' in file: 'exported_sec_spectra.csv'
    'reference_spectrum' in file: 'exported_sec_reference.csv'

  Files exported by ``SpectroECMeasurement`` in ixdat<0.2.3 will need these two lines modified, i.e.
  "spectra"->"spectrum_series" and "reference"->"reference_spectrum", before they can
  be imported. After this modification they can be read as before by
  ``Measurement.read(path_to_main_csv_file, reader="ixdat")``.

Debugging
---------

- ``PVMassSpecReader`` (reader="pfeiffer") now identifies columns as masses when they have
  names assigned during the data acquisition (e.g. "32_amu_(Oxygen)" as a column name),
  whereas before it would only have identified as masses columns which ended with "_amu".

- ``ECPlotter.plot_measurement`` and ``ECMSPlotter.plot_measuremnt`` raise warnings instead
  of ``SeriesNotFoundError``\ s if they can't find the requested
  potential (``U_name``) or current (``J_name``).

- ``SpectrumSeriesPlotter.heat_plot()`` now correctly orients its data. Previously it
  had been flipped in the vertical direction.


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





