"""Module for representation and analysis of EC measurements"""

from ..measurement_base import Measurement
from ..exporters import ECExporter
from ..plotters.ec_plotter import ECPlotter
from ..tools import deprecate
from ..calculators.ec_calculators import ECCalibration


class ECMeasurement(Measurement):
    """Class implementing electrochemistry measurements

    TODO: Implement a unit library for current and potential, A_el and RE_vs_RHE
      so that e.g. current can be seamlessly normalized to mass OR area.

    The main job of this class is making sure that the ValueSeries most essential for
    visualizing any normal electrochemistry measurements (i.e. excluding impedance
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
    default_calibration = ECCalibration

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
        if "potential" not in self.aliases:
            self._aliases.update({"potential": ["raw_potential"]})
        if "current" not in self.aliases:
            self._aliases.update({"current": ["raw_current"]})

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
    @deprecate("0.1", "Use `E_name` instead.", "0.3.1")
    def E_str(self):
        return self.E_name

    @property
    @deprecate("0.1", "Use `I_name` instead.", "0.3.1")
    def I_str(self):
        return self.I_name

    @property
    @deprecate("0.1", "Use `U_name` instead.", "0.3.1")
    def V_str(self):
        return self.U_name

    @property
    @deprecate("0.1", "Use `J_name` instead.", "0.3.1")
    def J_str(self):
        return self.J_name

    @property
    def aliases(self):
        """A dictionary with the names of other data series a given name can refer to"""
        a = super().aliases.copy()
        return a

    @property
    def ec_calibration(self):
        """A calibration joining the first RE_vs_RHE, A_el, and R_Ohm"""
        return self.calculators["EC calibration"]

    @property
    def RE_vs_RHE(self):
        """The refernce electrode potential on the RHE scale in [V]"""
        for calculator in self.calculator_list:
            if getattr(calculator, "RE_vs_RHE", None) is not None:
                return calculator.RE_vs_RHE

    @property
    def A_el(self):
        """The electrode area in [cm^2]"""
        for calculator in self.calculator_list:
            if getattr(calculator, "A_el", None) is not None:
                return calculator.A_el

    @property
    def R_Ohm(self):
        """The ohmic drop resistance in [Ohm]"""
        for calculator in self.calculator_list:
            if getattr(calculator, "R_Ohm", None) is not None:
                return calculator.R_Ohm

    def calibrate_RE(self, RE_vs_RHE):
        """Calibrate the reference electrode by providing `RE_vs_RHE` in [V]."""
        new_calibration = ECCalibration(
            RE_vs_RHE=RE_vs_RHE,
            measurement=self,
        )
        self.add_calculator(new_calibration)

    def normalize_current(self, A_el):
        """Normalize current to electrode surface area by providing `A_el` in [cm^2]."""
        new_calibration = ECCalibration(
            A_el=A_el,
            measurement=self,
        )
        self.add_calculator(new_calibration)

    def correct_ohmic_drop(self, R_Ohm):
        """Correct for ohmic drop by providing `R_Ohm` in [Ohm]."""
        new_calibration = ECCalibration(
            R_Ohm=R_Ohm,
            measurement=self,
        )
        self.add_calculator(new_calibration)

    @property
    def potential(self):
        return self["potential"]

    @property
    def current(self):
        return self["current"]

    @property
    @deprecate(
        "0.1", "Use a look-up, i.e. `ec_meas['raw_potential']`, instead.", "0.3.1"
    )
    def raw_potential(self):
        return self["raw_potential"]

    @property
    @deprecate("0.1", "Use a look-up, i.e. `ec_meas['raw_current']`, instead.", "0.3.1")
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
    @deprecate("0.1", "Use `U` instead.", "0.3.1")
    def v(self):
        """The potential [V] numpy array of the measurement"""
        return self.potential.data.copy()

    @property
    @deprecate("0.1", "Use `J` instead.", "0.3.1")
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
