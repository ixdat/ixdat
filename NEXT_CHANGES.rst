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


For ixdat 0.3.0
===============

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
  it doesn't do *quantification*. Instead, it does *calibrtion*, and its calibration
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
is implemented as a plugin. At present, this is in **src/ixdat/config.py** but it
will likely move.

- ``siqCalculator`` implements all the calibration and quantification methods that
  make use of the external ``spectro_inlets_quantification`` package. The calibration
  methods were previously methods of ``MSMeasurement`` and ``ECMSMeasurement`` prefixed
  "siq_", and the quantification methods were accessed through an overloading of
  ``MSMeasurement.grab()``.

  Constructor methods:
  - ``gas_flux_calibration``, formerly ``MSMeasurement.siq_gas_flux_calibration``.
  - ``gas_flux_calibration_curve``, formerly ``MSMeasurement.siq_gas_flux_calibration_curve``.
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


For ixdat 0.2.13
===============

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

