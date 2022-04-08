"""Module for representation and analysis of EC measurements"""

from ..measurements import Measurement, Calibration
from ..data_series import ValueSeries
from ..exporters.ec_exporter import ECExporter
from ..plotters.ec_plotter import ECPlotter
from ..tools import deprecate

EC_FANCY_NAMES = {
    "t": "time / [s]",
    "raw_potential": "raw potential / [V]",
    "potential": "$U_{RHE}$ / [V]",
    "raw_current": "raw current / [mA]",
    "current": "J / [mA cm$^{-2}$]",
}


class ECMeasurement(Measurement):
    """Class implementing electrochemistry measurements

    TODO: Implement a unit library for current and potential, A_el and RE_vs_RHE
      so that e.g. current can be seamlessly normalized to mass OR area.

    The main job of this class is making sure that the ValueSeries most essential for
    visualizing and normal electrochemistry measurements (i.e. excluding impedance
    spec., RRDE, etc, which would need new classes) are always available in the
    correct form as the measurement is added with others, reduced to a selection,
    calibrated and normalized, etc. These most important ValueSeries are:

    - `potential`: The working-electrode potential typically in [V].
      If `ec_meas` is an `ECMeasurement`, then `ec_meas["potential"]` always returns a
      `ValueSeries` characterized by:

        - calibrated and/or corrected, if the measurement has been calibrated with the
          reference electrode potential (`RE_vs_RHE`, see `calibrate`) and/or corrected
          for ohmic drop (`R_Ohm`, see `correct_ohmic_drop`).
        - A name that makes clear any calibration and/or correction
        - Data which spans the entire timespan of the measurement - i.e. whenever EC
          data is being recorded, `potential` is there, even if the name of the raw
          `ValueSeries` (what the acquisition software calls it) changes. Indeed
          `ec_meas["potential"].tseries` is the measurement's definitive time variable.

    - `current`: The working-electrode current typically in [mA] or [mA/cm^2].
      `ec_meas["current"]` always returns a `ValueSeries` characterized by:

        - normalized if the measurement has been normalized with the electrode area
          (`A_el`, see `normalize`)
        - A name that makes clear whether it is normalized
        - Data which spans the entire timespan of the measurement

    - `selector`: A counter series distinguishing sections of the measurement program.
      This is essential for analysis of complex measurements as it allows for
      corresponding parts of experiments to be isolated and treated identically.
      `selector` in `ECMeasurement` is defined to increment each time one or more of
      the following changes:

        - `loop_number`: A parameter saved by some potentiostats (e.g. BioLogic) which
          allow complex looped electrochemistry programs.
        - `file_number`: The id of the component measurement from which each section of
          the data (the origin of each `ValueSeries` concatenated to `potential`)
        - `cycle_number`: An incrementer within a file saved by a potentiostat.

    The names of these ValueSeries, which can also be used to index the measurement, are
    conveniently available as properties:

    - `ec_meas.t_name` is the name of the definitive time, i.e. that of the potential
    - `ec_meas.E_name` is the name of the raw potential
    - `ec_meas.U_name` is the name of the calibrated and/or corrected potential
    - `ec_meas.I_name` is the name of the raw current
    - `ec_meas.J_name` is the name of the normalized current
    - `ec_meas.selector_name` is the name of the default selector, i.e. "selector"

    Numpy arrays from important `DataSeries` are directly accessible via attributes:

    - `ec_meas.t` for `ec_meas["potential"].t`
    - `ec_meas.U` for `ec_meas["potential"].data`
    - `ec_meas.J` for `ec_meas["current"].data`

    `ECMeasurement` comes with an `ECPlotter` which either plots `potential` and
    `current` against time (`ec_meas.plot_measurement()`) or plots `current` against
    `potential (`ec_meas.plot_vs_potential()`).

    It turns out that keeping track of current, potential, and selector when combining
    datasets is enough of a job to fill a class. Thus, the more exciting
    electrochemistry-related functionality should be implemented in inheriting classes
    such as `CyclicVoltammogram`.
    """

    extra_column_attrs = {
        "ec_meaurements": {
            "ec_technique",
        }
    }
    control_series_name = "raw_potential"
    essential_series_names = ("t", "raw_potential", "raw_current")
    selection_series_names = ("file_number", "loop_number", "cycle number", "Ns")
    default_exporter = ECExporter
    default_plotter = ECPlotter

    def __init__(
        self,
        name,
        *,
        ec_technique=None,
        RE_vs_RHE=None,
        R_Ohm=None,
        A_el=None,
        **kwargs,
    ):
        """initialize an electrochemistry measurement

        Args:
            name (str): The name of the measurement
            ec_technique (str): The electrochemistry sub-technique
            RE_vs_RHE (float): The reference electrode potential on the RHE scale in [V]
            R_Ohm (float): The ohmic drop resistance in [Ohm]
            A_el (float): The electrode area in [cm^2]

        Kwargs, passed on to `Measurement.__init__` (see :class:`.Measurement`):
            metadata (dict): Free-form measurement metadata. Must be json-compatible.
            technique (str): The measurement technique
            s_ids (list of int): The id's of the measurement's DataSeries, if
                to be loaded (instead of given directly in series_list)
            series_list (list of DataSeries): The measurement's DataSeries
            m_ids (list of int): The id's of the component measurements, if to be
                loaded. None unless this is a combined measurement (typically
                corresponding to more than one file).
            component_measurements (list of Measurements): The measurements of which
                this measurement is a combination
            reader (Reader): The file reader (None unless read from a file)
            plotter (Plotter): The visualization tool for the measurement
            exporter (Exporter): The exporting tool for the measurement
            sample (Sample or str): The sample being measured
            lablog (LabLog): The log entry with e.g. notes taken during the measurement
            tstamp (float): The nominal starting time of the measurement, used for
                data selection, visualization, and exporting.
        """
        super().__init__(name, **kwargs)

        self.ec_technique = ec_technique
        if RE_vs_RHE is not None or A_el is not None or R_Ohm is not None:
            self.calibrate(RE_vs_RHE, A_el, R_Ohm)
        self.plot_vs_potential = self.plotter.plot_vs_potential

    @property
    def E_name(self):
        return self["raw_potential"].name

    @property
    def I_name(self):
        return self["raw_current"].name

    @property
    def U_name(self):
        return self.potential.name

    @property
    def J_name(self):
        return self.current.name

    @property
    @deprecate("0.1", "Use `E_name` instead.", "0.3")
    def E_str(self):
        return self.E_name

    @property
    @deprecate("0.1", "Use `I_name` instead.", "0.3")
    def I_str(self):
        return self.I_name

    @property
    @deprecate("0.1", "Use `U_name` instead.", "0.3")
    def V_str(self):
        return self.U_name

    @property
    @deprecate("0.1", "Use `J_name` instead.", "0.3")
    def J_str(self):
        return self.J_name

    @property
    def aliases(self):
        """A dictionary with the names of other data series a given name can refer to"""
        a = super().aliases.copy()
        return a

    @property
    def calibrations(self):
        """The list of calibrations of the measurement.

        The following is necessary to ensure that all EC Calibration parameters are
        joined in a single calibration when processing. So that "potential" is both
        calibrated to RHE and ohmic drop corrected, even if the two calibration
        parameters were added separately.
        """
        full_calibration_list = self.calibration_list
        good_calibration_list = [self.ec_calibration]
        for calibration in full_calibration_list:
            if calibration.__class__ is ECCalibration:
                # Then we have all we need from it
                continue
            good_calibration_list.append(calibration)
        return good_calibration_list

    @property
    def ec_calibration(self):
        """A calibration joining the first RE_vs_RHE, A_el, and R_Ohm"""
        return ECCalibration(RE_vs_RHE=self.RE_vs_RHE, A_el=self.A_el, R_Ohm=self.R_Ohm)

    @property
    def RE_vs_RHE(self):
        """The refernce electrode potential on the RHE scale in [V]"""
        for calibration in self.calibration_list:
            if getattr(calibration, "RE_vs_RHE", None) is not None:
                return calibration.RE_vs_RHE

    @property
    def A_el(self):
        """The electrode area in [cm^2]"""
        for calibration in self.calibration_list:
            if getattr(calibration, "A_el", None) is not None:
                return calibration.A_el

    @property
    def R_Ohm(self):
        """The ohmic drop resistance in [Ohm]"""
        for calibration in self.calibration_list:
            if getattr(calibration, "R_Ohm", None) is not None:
                return calibration.R_Ohm

    def calibrate_RE(self, RE_vs_RHE):
        """Calibrate the reference electrode by providing `RE_vs_RHE` in [V]."""
        new_calibration = ECCalibration(
            RE_vs_RHE=RE_vs_RHE,
            measurement=self,
        )
        self.add_calibration(new_calibration)

    def normalize_current(self, A_el):
        """Normalize current to electrode surface area by providing `A_el` in [cm^2]."""
        new_calibration = ECCalibration(
            A_el=A_el,
            measurement=self,
        )
        self.add_calibration(new_calibration)

    def correct_ohmic_drop(self, R_Ohm):
        """Correct for ohmic drop by providing `R_Ohm` in [Ohm]."""
        new_calibration = ECCalibration(
            R_Ohm=R_Ohm,
            measurement=self,
        )
        self.add_calibration(new_calibration)

    @property
    def potential(self):
        return self["potential"]

    @property
    def current(self):
        return self["current"]

    @property
    @deprecate("0.1", "Use a look-up, i.e. `ec_meas['raw_potential']`, instead.", "0.3")
    def raw_potential(self):
        return self["raw_potential"]

    @property
    @deprecate("0.1", "Use a look-up, i.e. `ec_meas['raw_current']`, instead.", "0.3")
    def raw_current(self):
        return self["raw_current"]

    @property
    def U(self):
        """The potential [V] numpy array of the measurement"""
        return self.potential.data.copy()

    @property
    def J(self):
        """The current ([mA] or [mA/cm^2]) numpy array of the measurement"""
        return self.current.data.copy()

    @property
    @deprecate("0.1", "Use `U` instead.", "0.3")
    def v(self):
        """The potential [V] numpy array of the measurement"""
        return self.potential.data.copy()

    @property
    @deprecate("0.1", "Use `J` instead.", "0.3")
    def j(self):
        """The current ([mA] or [mA/cm^2]) numpy array of the measurement"""
        return self.current.data.copy()

    def as_cv(self):
        """Convert self to a CyclicVoltammogram"""
        from .cv import CyclicVoltammogram

        cv_as_dict = self.as_dict()
        cv_as_dict["technique"] = "CV"
        # Note, this works perfectly! All needed information is in self_as_dict :)
        return CyclicVoltammogram.from_dict(cv_as_dict)


class ECCalibration(Calibration):
    """An electrochemical calibration with RE_vs_RHE, A_el, and/or R_Ohm"""

    extra_column_attrs = {"ec_calibration": {"RE_vs_RHE", "A_el", "R_Ohm"}}
    # TODO: https://github.com/ixdat/ixdat/pull/11#discussion_r677552828

    def __init__(
        self,
        name=None,
        technique="EC",
        tstamp=None,
        measurement=None,
        RE_vs_RHE=None,
        A_el=None,
        R_Ohm=None,
    ):
        """Initiate a Calibration

        Args:
            name (str): The name of the calibration
            technique (str): The technique of the calibration
            tstamp (float): The time at which the calibration took place or is valid
            measurement (ECMeasurement): Optional. A measurement to calibrate by default
            RE_vs_RHE (float): The reference electrode potential on the RHE scale in [V]
            A_el (float): The electrode area in [cm^2]
            R_Ohm (float): The ohmic drop resistance in [Ohm]
        """
        super().__init__(
            name=name, technique=technique, tstamp=tstamp, measurement=measurement
        )
        self.RE_vs_RHE = RE_vs_RHE
        self.A_el = A_el
        self.R_Ohm = R_Ohm

    def __repr__(self):
        # TODO: make __repr__'s consistent throught ixdat
        return (
            f"{self.__class__.__name__}"
            f"(RE_vs_RHE={self.RE_vs_RHE}, A_el={self.A_el}, R_Ohm={self.R_Ohm})"
        )

    def calibrate_series(self, key, measurement=None):
        """Return a calibrated series for key based on the raw data in the measurement.

        Key should be "potential" or "current". Anything else will return None.

        - "potential": the calibration looks up "raw_potential" in the measurement,
        shifts it to the RHE potential if RE_vs_RHE is available, corrects it for
        Ohmic drop if R_Ohm is available, and then returns a calibrated potential
        series with a name indicative of the corrections done.
        - "current": The calibration looks up "raw_current" in the measurement,
        normalizes it to the electrode area if A_el is available, and returns a
        calibrated current series with a name indicative of whether the normalization
        was done.
        """
        measurement = measurement or self.measurement
        if key == "potential":
            raw_potential = measurement["raw_potential"]
            name = raw_potential.name
            U = raw_potential.data
            if self.RE_vs_RHE is not None:
                U = U + self.RE_vs_RHE
                name = EC_FANCY_NAMES["potential"]
            if self.R_Ohm is not None:
                I_mA = measurement.grab_for_t("raw_current", t=raw_potential.t)
                U = U - self.R_Ohm * I_mA * 1e-3  # [V] = [Ohm*mA*(A/mA)]
                name = name + " $_{ohm. corr.}$"
            return ValueSeries(
                name=name,
                unit_name=raw_potential.unit_name,
                data=U,
                tseries=raw_potential.tseries,
            )

        if key == "current":
            raw_current = measurement["raw_current"]
            name = raw_current.name
            J = raw_current.data
            if self.A_el is not None:
                J = J / self.A_el
                name = EC_FANCY_NAMES["current"]
            return ValueSeries(
                name=name,
                unit_name=raw_current.unit_name,
                data=J,
                tseries=raw_current.tseries,
            )
