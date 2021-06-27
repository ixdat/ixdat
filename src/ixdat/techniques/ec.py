"""Module for representation and analysis of EC measurements"""

import numpy as np

from ..measurements import Measurement, append_series
from ..data_series import ValueSeries, ConstantValue, time_shifted
from ..exceptions import SeriesNotFoundError
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
            "RE_vs_RHE",
            "R_Ohm",
            "A_el",
        }
    }
    control_str = "raw_potential"
    default_exporter_class = ECExporter
    default_plotter_class = ECPlotter

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

        self.RE_vs_RHE = RE_vs_RHE
        self.R_Ohm = R_Ohm
        self.A_el = A_el

        self.V_str = EC_FANCY_NAMES["potential"]
        self.J_str = EC_FANCY_NAMES["current"]
        self.E_str = EC_FANCY_NAMES["raw_potential"]
        self.I_str = EC_FANCY_NAMES["raw_current"]

        self.plot_vs_potential = self.plotter.plot_vs_potential

        self._raw_potential = None
        self._raw_current = None
        self._selector = None
        self._file_number = None

    def _populate_constants(self):
        """Replace any ConstantValues with ValueSeries on potential's tseries

        TODO: This function flagrantly violates laziness. Not only does it fill up all
            the ConstantValue's with long vectors before they're needed, it also forces
            raw_potential to be built before it is needed.
            A lazier solution is needed.
        """
        for (i, s) in enumerate(self.series_list):
            if isinstance(s, ConstantValue):
                tseries = self.potential.tseries
                series = s.get_vseries(tseries=tseries)
                self.series_list[i] = series

    @property
    def aliases(self):
        a = self._aliases.copy()
        a.update({value: [key] for (key, value) in EC_FANCY_NAMES.items()})
        return a

    def calibrate(self, RE_vs_RHE=None, A_el=None, R_Ohm=None):
        """Calibrate the EC measurement (all args optional)

        Args:
            RE_vs_RHE (float): reference electode potential on RHE scale in [V]
            A_el (float): electrode area in [cm^2]
            R_Ohm (float): ohmic drop resistance in [Ohm]
        """
        if RE_vs_RHE:
            self.calibrate_RE(RE_vs_RHE=RE_vs_RHE)
        if A_el:
            self.normalize_current(A_el=A_el)
        if R_Ohm:
            self.correct_ohmic_drop(R_Ohm=R_Ohm)

    def calibrate_RE(self, RE_vs_RHE):
        """Calibrate the reference electrode by providing `RE_vs_RHE` in [V]."""
        self.RE_vs_RHE = RE_vs_RHE

    def normalize_current(self, A_el):
        """Normalize current to electrod surface area by providing `A_el` in [cm^2]."""
        self.A_el = A_el

    def correct_ohmic_drop(self, R_Ohm):
        """Correct for ohmic drop by providing `R_Ohm` in [Ohm]."""
        self.R_Ohm = R_Ohm

    @property
    def potential(self):
        """The ValueSeries with the ECMeasurement's potential.

        This is result of the following:
        - Starts with `self.raw_potential`
        - if the measurement is "calibrated" i.e. `RE_vs_RHE` is not None: add
            `RE_vs_RHE` to the potential data and change its name from `E_str` to `V_str`
        - if the measurement is "corrected" i.e. `R_Ohm` is not None: subtract
            `R_Ohm` times the raw current from the potential and add " (corrected)" to
            its name.
        """
        raw_potential = self["raw_potential"]
        if self.RE_vs_RHE is None and self.R_Ohm is None:
            return raw_potential
        fixed_V_str = raw_potential.name
        fixed_potential_data = raw_potential.data
        fixed_unit_name = raw_potential.unit_name
        if self.RE_vs_RHE:
            fixed_V_str = self.V_str
            fixed_potential_data = fixed_potential_data + self.RE_vs_RHE
            fixed_unit_name = "V <RHE>"
        if self.R_Ohm:
            fixed_V_str += " (corrected)"
            fixed_potential_data = (
                fixed_potential_data - self.R_Ohm * self["raw_current"].data * 1e-3
            )  # TODO: Units. The 1e-3 here is to bring raw_current.data from [mA] to [A]
        return ValueSeries(
            name=fixed_V_str,
            data=fixed_potential_data,
            unit_name=fixed_unit_name,
            tseries=raw_potential.tseries,
        )  # TODO: Better cache'ing. This is not cached at all.

    @property
    def current(self):
        """The ValueSeries with the ECMeasurement's current.

        This is result of the following:
        - Starts with `self.raw_current`
        - if the measurement is "normalized" i.e. `A_el` is not None: divide the current
            data by `A_el`, change its name from `I_str` to `J_str`, and add `/cm^2` to
            its unit.
        """
        raw_current = self["raw_current"]
        if self.A_el is None:
            return raw_current
        else:
            return ValueSeries(
                name=self.J_str,
                data=raw_current.data / self.A_el,
                unit_name=raw_current.unit_name + "/cm^2",
                tseries=raw_current.tseries,
            )

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
        self_as_dict["technique"] = "CV"
        del self_as_dict["s_ids"]
        # Note, this works perfectly! All needed information is in self_as_dict :)
        return CyclicVoltammagram.from_dict(self_as_dict)
