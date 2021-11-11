"""Module for representation and analysis of EC measurements"""

import numpy as np

from ..measurements import Measurement, Calibration
from ..data_series import ValueSeries
from ..exporters.ec_exporter import ECExporter
from ..plotters.ec_plotter import ECPlotter


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
    TODO:   so that e.g. current can be seamlessly normalized to mass OR area.
    The main job of this class is making sure that the ValueSeries most essential for
    visualizing and normal electrochemistry measurements (i.e. excluding impedance spec.,
    RRDE, etc, which would need new classes) are always available in the correct form as
    the measurement is added with others, reduced to a selection, calibrated and
    normalized, etc. These most important ValueSeries are:

    - `potential`: The working-electrode potential typically in [V].
      If `ec_meas` is an `ECMeasurement`, then `ec_meas["potential"]` always returns a
      `ValueSeries` characterized by:

        - calibrated and/or corrected, if the measurement has been calibrated with the
          reference electrode potential (`RE_vs_RHE`, see `calibrate`) and/or corrected
          for ohmic drop (`R_Ohm`, see `correct_ohmic_drop`).
        - A name that makes clear any calibration and/or correction
        - Data which spans the entire timespan of the measurement - i.e. whenever EC data
          is being recorded, `potential` is there, even if the name of the raw
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
    - `ec_meas.v_name` is the name of the calibrated and/or corrected potential
    - `ec_meas.I_name` is the name of the raw current
    - `ec_meas.j_name` is the name of the normalized current
    - `ec_meas.selector_name` is the name of the default selector, i.e. "selector"

    Numpy arrays from important `DataSeries` are also directly accessible via attributes:

    - `ec_meas.t` for `ec_meas["potential"].t`
    - `ec_meas.v` for `ec_meas["potential"].data`
    - `ec_meas.j` for `ec_meas["current"].data`

    `ECMeasurement` comes with an `ECPlotter` which either plots `potential` and
    `current` against time (`ec_meas.plot_measurement()`) or plots `current` against
    `potential (`ec_meas.plot_vs_potential()`).

    It turns out that keeping track of current, potential, and selector when combining
    datasets is enough of a job to fill a class. Thus, the more exciting
    electrochemistry-related functionality should be implemented in inheriting classes
    such as `CyclicVoltammagram`.
    """

    extra_column_attrs = {
        "ec_meaurements": {
            "ec_technique",
        }
    }
    control_series_name = "raw_potential"
    essential_series_names = ("t", "raw_potential", "raw_current", "cycle")
    selection_series_names = ("file_number", "loop_number", "cycle")
    default_exporter = ECExporter
    default_plotter = ECPlotter
    v_name = EC_FANCY_NAMES["potential"]
    j_name = EC_FANCY_NAMES["current"]
    E_name = EC_FANCY_NAMES["raw_potential"]
    I_name = EC_FANCY_NAMES["raw_current"]

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
        if RE_vs_RHE or A_el or R_Ohm:
            self.calibrate(RE_vs_RHE, A_el, R_Ohm)
        self.plot_vs_potential = self.plotter.plot_vs_potential

    @property
    def aliases(self):
        """A dictionary with the names of other data series a given name can refer to"""
        a = super().aliases.copy()
        a.update({value: [key] for (key, value) in EC_FANCY_NAMES.items()})
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

    def calibrate(
        self,
        RE_vs_RHE=None,
        A_el=None,
        R_Ohm=None,
        tstamp=None,
        cal_name=None,
    ):
        """Calibrate the EC measurement (all args optional)

        Args:
            RE_vs_RHE (float): reference electode potential on RHE scale in [V]
            A_el (float): electrode area in [cm^2]
            R_Ohm (float): ohmic drop resistance in [Ohm]
            tstamp (flaot): The timestamp at which the calibration was done (defaults
                to now)
            cal_name (str): The name of the calibration.
        """
        if not (RE_vs_RHE or A_el or R_Ohm):
            print("Warning! Ignoring attempt to calibrate without any parameters.")
            return
        new_calibration = ECCalibration(
            RE_vs_RHE=RE_vs_RHE or self.RE_vs_RHE,
            A_el=A_el or self.A_el,
            R_Ohm=R_Ohm or self.R_Ohm,
            tstamp=tstamp,
            name=cal_name,
            measurement=self,
        )
        self.add_calibration(new_calibration)
        self.clear_cache()

    def calibrate_RE(self, RE_vs_RHE):
        """Calibrate the reference electrode by providing `RE_vs_RHE` in [V]."""
        new_calibration = ECCalibration(
            RE_vs_RHE=RE_vs_RHE,
            measurement=self,
        )
        self.add_calibration(new_calibration)
        self.clear_cache()

    def normalize_current(self, A_el):
        """Normalize current to electrod surface area by providing `A_el` in [cm^2]."""
        new_calibration = ECCalibration(
            A_el=A_el,
            measurement=self,
        )
        self.add_calibration(new_calibration)
        self.clear_cache()

    def correct_ohmic_drop(self, R_Ohm):
        """Correct for ohmic drop by providing `R_Ohm` in [Ohm]."""
        new_calibration = ECCalibration(
            R_Ohm=R_Ohm,
            measurement=self,
        )
        self.add_calibration(new_calibration)
        self.clear_cache()

    @property
    def potential(self):
        return self["potential"]

    @property
    def current(self):
        return self["current"]

    def grab_potential(self, tspan=None):
        """Return the time [s] and potential [V] vectors cut by tspan

        TODO: I think this is identical, now that __getitem__ finds potential, to
            self.grab("potential", tspan=tspan)
        """
        t = self.potential.t.copy()
        v = self.potential.data.copy()
        if tspan:
            mask = np.logical_and(tspan[0] < t, t < tspan[-1])
            t = t[mask]
            v = v[mask]
        return t, v

    def grab_current(self, tspan=None):
        """Return the time [s] and current ([mA] or [mA/cm^2]) vectors cut by tspan"""
        t = self.current.t.copy()
        j = self.current.data.copy()
        if tspan:
            mask = np.logical_and(tspan[0] < t, t < tspan[-1])
            t = t[mask]
            j = j[mask]
        return t, j

    @property
    def v(self):
        """The potential [V] numpy array of the measurement"""
        return self.potential.data.copy()

    @property
    def j(self):
        """The current ([mA] or [mA/cm^2]) numpy array of the measurement"""
        return self.current.data.copy()

    def as_cv(self):
        """Convert self to a CyclicVoltammagram"""
        from .cv import CyclicVoltammagram

        cv_as_dict = self.as_dict()
        cv_as_dict["technique"] = "CV"
        # Note, this works perfectly! All needed information is in self_as_dict :)
        return CyclicVoltammagram.from_dict(cv_as_dict)


class ECCalibration(Calibration):
    """An electrochemical calibration with RE_vs_RHE, A_el, and/or R_Ohm"""

    extra_column_attrs = {"ec_calibration": {"RE_vs_RHE", "A_el", "R_Ohm"}}
    # TODO: https://github.com/ixdat/ixdat/pull/11#discussion_r677552828

    def __init__(
        self,
        technique="EC",
        tstamp=None,
        name=None,
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
            measurement (ECMeasurement): Optional. A measurement to calibrate by default.
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

        - potential: the calibration looks up "raw_potential" in the measurement, shifts
        it to the RHE potential if RE_vs_RHE is available, corrects it for Ohmic drop if
        R_Ohm is available, and then returns a calibrated potential series with a name
        indicative of the corrections done.
        - current: The calibration looks up "raw_current" in the measurement, normalizes
        it to the electrode area if A_el is available, and returns a calibrated current
        series with a name indicative of whether the normalization was done.
        """
        measurement = measurement or self.measurement
        if key == "potential":
            raw_potential = measurement["raw_potential"]
            name = raw_potential.name
            v = raw_potential.data
            if self.RE_vs_RHE:
                v = v + self.RE_vs_RHE
                name = measurement.v_name or EC_FANCY_NAMES["potential"]
            if self.R_Ohm:
                i_mA = measurement.grab_for_t("raw_current", t=raw_potential.t)
                v = v - self.R_Ohm * i_mA * 1e-3  # [V] = [Ohm*mA*(A/mA)]
                name = name + " $_{ohm. corr.}$"
            return ValueSeries(
                name=name,
                unit_name=raw_potential.unit_name,
                data=v,
                tseries=raw_potential.tseries,
            )

        if key == "current":
            raw_current = measurement["raw_current"]
            name = raw_current.name
            v = raw_current.data
            if self.A_el:
                v = v / self.A_el
                name = measurement.j_name or EC_FANCY_NAMES["current"]
            return ValueSeries(
                name=name,
                unit_name=raw_current.unit_name,
                data=v,
                tseries=raw_current.tseries,
            )
