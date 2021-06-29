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

    TODO:
        Implement a unit library for current and potential, A_el and RE_vs_RHE
        so that e.g. current can be seamlessly normalized to mass OR area.

    The main job of this class is making sure that the ValueSeries most essential for
    visualizing and normal electrochemistry measurements (i.e. excluding impedance spec,
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
            is being recorded, `potential` is there, even the name of the raw
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
        `selector` in `ECMeasurement` is defined to incriment each time one or more of
        the following changes:
        - `loop_number`: A parameter saved by some potentiostats (e.g. BioLogic) which
            allow complex looped electrochemistry programs.
        - `file_number`: The id of the component measurement from which each section of
            the data (the origin of each `ValueSeries` concatenated to `potential`)
        - `cycle_number`: An incrementer within a file saved by a potentiostat.

    The names of these ValueSeries, which can also be used to index the measurement, are
    conveniently available as properties:
    - `ec_meas.t_str` is the name of the definitive time, which corresponds to potential.
    - `ec_meas.E_str` is the name of the raw potential
    - `ec_meas.V_str` is the name to the calibrated and/or corrected potential
    - `ec_meas.I_str` is the name of the raw current
    - `ec_meas.J_str` is the name of the normalized current.
    - `ec_meas.sel_str` is the name of the default selector, i.e. "selector"

    Numpy arrays from important `DataSeries` are also directly accessible via attributes:
    - `ec_meas.t` for `ec_meas["potential"].t`
    - `ec_meas.v` for `ec_meas["potential"].data`
    - `ec_meas.j` for `ec_meas["current"].data

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
    control_str = "raw_potential"
    essential_series = ("t", "raw_potential", "raw_current", "cycle")
    select_on = ("file_number", "loop_number", "cycle")
    default_exporter_class = ECExporter
    default_plotter_class = ECPlotter
    V_str = EC_FANCY_NAMES["potential"]
    J_str = EC_FANCY_NAMES["current"]
    E_str = EC_FANCY_NAMES["raw_potential"]
    I_str = EC_FANCY_NAMES["raw_current"]

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
                TODO: Decide if metadata needs the json string option.
                    See: https://github.com/ixdat/ixdat/pull/1#discussion_r546436991
            metadata (dict): Free-form measurement metadata
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
            sample (Sample): The (already loaded) sample being measured
            lablog (LabLog): The log entry with e.g. notes taken during the measurement
            tstamp (float): The nominal starting time of the measurement, used for
                data selection, visualization, and exporting.
            t_str (str): Name of the main time variable (corresponding to potential)
            E_str (str): Name of raw potential (so called because potential is saved
                as "Ewe/V" in biologic .mpt files)
            V_str (str): Name of calibrated potential
            RE_vs_RHE (float): Reference electrode potential in [V] on the RHE scale.
                If RE_vs_RHE is not None, the measurement is considered *calibrated*,
                and will use the calibrated potential `self[self.V_str]` by default
                TODO: Unit
            R_Ohm (float): Ohmic drop in [Ohm]. If R_Ohm is not None, the ohmic drop
                is corrected for when returning potential.
                TODO: Unit
            raw_potential_names (tuple of str): The names of the VSeries which
                represent raw working electrode current. This is typically how the data
                acquisition software saves potential
            I_str (str): Name of raw current
            J_str (str): Name of normalized current
            A_el (float): Area of electrode in [cm^2].
                If A_el is not None, the measurement is considered *normalized*,
                and will use the calibrated potential `self[self.V_str]` by default
                TODO: Unit
            raw_current_names (tuple of str): The names of the VSeries which represent
                raw working electrode current. This is typically how the data
                acquisition software saves current.
        """
        super().__init__(name, **kwargs)

        self.ec_technique = ec_technique
        self.calibrate(RE_vs_RHE, A_el, R_Ohm)
        self.plot_vs_potential = self.plotter.plot_vs_potential

    @property
    def aliases(self):
        a = self._aliases.copy()
        a.update({value: [key] for (key, value) in EC_FANCY_NAMES.items()})
        return a

    @property
    def calibration_list(self):
        """The list of calibrations of the measurement"""
        full_calibration_list = super().calibration_list
        # The following is necessary to ensure that all EC Calibration parameters are
        # joined in a single calibration when processing. So that "potential" is both
        # calibrated to RHE and ohmic drop corrected, even if the two calibration
        # parameters were added separately.
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
        for calibration in self._calibration_list:
            if hasattr(calibration, "RE_vs_RHE") and calibration.RE_vs_RHE is not None:
                return calibration.RE_vs_RHE

    @property
    def A_el(self):
        for calibration in self._calibration_list:
            if hasattr(calibration, "A_el") and calibration.A_el is not None:
                return calibration.A_el

    @property
    def R_Ohm(self):
        for calibration in self._calibration_list:
            if hasattr(calibration, "R_Ohm") and calibration.R_Ohm is not None:
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
        new_calibration = ECCalibration(
            RE_vs_RHE=RE_vs_RHE or self.RE_vs_RHE,
            A_el=A_el or self.A_el,
            R_Ohm=R_Ohm or self.R_Ohm,
            tstamp=tstamp,
            name=cal_name,
            measurement=self,
        )
        self._calibration_list = [new_calibration] + self._calibration_list
        self.clear_cache()

    def calibrate_RE(self, RE_vs_RHE):
        """Calibrate the reference electrode by providing `RE_vs_RHE` in [V]."""
        new_calibration = ECCalibration(
            RE_vs_RHE=RE_vs_RHE,
            measurement=self,
        )
        self._calibration_list = [new_calibration] + self._calibration_list
        self.clear_cache()

    def normalize_current(self, A_el):
        """Normalize current to electrod surface area by providing `A_el` in [cm^2]."""
        new_calibration = ECCalibration(
            A_el=A_el,
            measurement=self,
        )
        self._calibration_list = [new_calibration] + self._calibration_list
        self.clear_cache()

    def correct_ohmic_drop(self, R_Ohm):
        """Correct for ohmic drop by providing `R_Ohm` in [Ohm]."""
        new_calibration = ECCalibration(
            R_Ohm=R_Ohm,
            measurement=self,
        )
        self._calibration_list = [new_calibration] + self._calibration_list
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

        self_as_dict = self.as_dict()
        self_as_dict["series_list"] = self.series_list
        self_as_dict["calibration_list"] = self._calibration_list
        self_as_dict["technique"] = "CV"
        del self_as_dict["s_ids"]
        # Note, this works perfectly! All needed information is in self_as_dict :)
        return CyclicVoltammagram.from_dict(self_as_dict)


class ECCalibration(Calibration):
    extra_column_attrs = {"ec_calibration": {"RE_vs_RHE", "A_el", "R_Ohm"}}

    def __init__(
        self,
        technique="EC",
        RE_vs_RHE=None,
        A_el=None,
        R_Ohm=None,
        tstamp=None,
        name=None,
        measurement=None,
    ):
        super().__init__(
            name=name, technique=technique, tstamp=tstamp, measurement=measurement
        )
        self.RE_vs_RHE = RE_vs_RHE
        self.A_el = A_el
        self.R_Ohm = R_Ohm

    def __repr__(self):
        return (
            f"{self.__class__.__name__}"
            f"(RE_vs_RHE={self.RE_vs_RHE}, A_el={self.A_el}, R_Ohm={self.R_Ohm})"
        )

    def calibrate_series(self, key, measurement=None):
        measurement = measurement or self.measurement
        if key == "potential":
            raw_potential = measurement["raw_potential"]
            name = raw_potential.name
            v = raw_potential.data
            if self.RE_vs_RHE:
                v = v + self.RE_vs_RHE
                name = measurement.V_str or EC_FANCY_NAMES["potential"]
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
                name = measurement.J_str or EC_FANCY_NAMES["current"]
            return ValueSeries(
                name=name,
                unit_name=raw_current.unit_name,
                data=v,
                tseries=raw_current.tseries,
            )
