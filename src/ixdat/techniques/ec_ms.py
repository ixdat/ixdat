"""Module for representation and analysis of EC-MS measurements"""
import numpy as np
import warnings
from scipy.optimize import minimize
from ..constants import FARADAY_CONSTANT
from .ec import ECMeasurement, ECCalibration
from .ms import MSMeasurement, MSSpectroMeasurement, MSCalResult, MSCalculator
from .cv import CyclicVoltammogram
from ..exceptions import QuantificationError
from ..exporters import ECMSExporter
from ..plotters import ECMSPlotter
from ..plotters.ms_plotter import STANDARD_COLORS
from ..config import plugins
from ..tools import deprecate


class ECMSMeasurement(ECMeasurement, MSMeasurement):
    """Class for raw EC-MS functionality. Parents: ECMeasurement and MSMeasurement"""

    extra_column_attrs = {
        "ecms_meaurements": {"ec_technique", "tspan_bg"},
    }
    # FIXME: It would be much more elegant if this carried over automatically from
    #  *both* parents, by appending the table columns...
    #  We'll see how the problem changes with the metaprogramming work.

    default_plotter = ECMSPlotter
    default_exporter = ECMSExporter

    @property
    def ec_plotter(self):
        """A plotter for just plotting the ec data"""
        return self.plotter.ec_plotter  # the ECPlotter of the measurement's ECMSPlotter

    @property
    def ms_plotter(self):
        """A plotter for just plotting the ms data"""
        return self.plotter.ms_plotter  # the MSPlotter of the measurement's ECMSPlotter

    @classmethod
    def from_dict(cls, obj_as_dict):
        """Initiate an ECMSMeasurement from a dictionary representation.

        This unpacks the ECMSCalibration from its own nested dictionary
        TODO: Figure out a way for that to happen automatically.
        """

        if "calibration" in obj_as_dict:
            if isinstance(obj_as_dict["calibration"], dict):
                # FIXME: This is a mess
                obj_as_dict["calibration"] = ECMSCalculator.from_dict(
                    obj_as_dict["calibration"]
                )
        obj = super(ECMSMeasurement, cls).from_dict(obj_as_dict)
        return obj

    @property
    def tspan(self):
        """The tspan of an MS measurement is the tspan of its potential data"""
        return [self.t[0], self.t[-1]]

    def as_cv(self):
        self_as_dict = self.as_dict()

        # FIXME: The following lines are only necessary because
        #  PlaceHolderObject.get_object isn't able to find things in the MemoryBackend
        del self_as_dict["s_ids"]
        self_as_dict["series_list"] = self.series_list

        ecms_cv = ECMSCyclicVoltammogram.from_dict(self_as_dict)

        return ecms_cv

    # --- METHODS WHICH HAVE BEEN MOVED TO `Calculator` CLASSES ---- #

    def _get_tspan_list(
        self,
        selector_list,
        selector_name=None,
        t_steady_pulse=None,
    ):
        """
        Generate a t_span list from input of selectors.

        This is useful for e.g. calibration curves

        Args:
            selector_list (list of selector): selector numbers that define the
                                            tspans over which data should be integrated
            selector_name (str): name of selector that will be used to determine sections
                                of data. Will refer to data['selector'] by default.
                                selector_name cannot contain a space character due to
                                limitations of self.select_values().
            t_steady_pulse (float): length of steady state pulse period to integrate
                                    (will choose the last x seconds of the period).
                                    Defaults to None: uses entire steady state pulse

        Returns tspan_list(list of tspan)
        """
        selector_name = selector_name or "selector"
        t_idx = -1
        if not t_steady_pulse:
            t_idx = 0
            t_steady_pulse = 0
        tspan_list = [
            [
                self.select_values(**{selector_name: selector_value}).grab("t")[0][t_idx]
                - t_steady_pulse,
                self.select_values(**{selector_name: selector_value}).grab("t")[0][-1],
            ]
            for selector_value in selector_list
        ]
        print("Following tspans were selected for calibration: " + str(tspan_list))
        return tspan_list

    @deprecate(
        "0.2.13",
        "Use `ECMSCalulator.ecms_calibration` instead.",
        "0.3",
    )
    def ecms_calibration(self, mol, mass, n_el, tspan, tspan_bg=None):
        return ECMSCalculator.ecms_calibration(
            measurement=self,
            mol=mol,
            mass=mass,
            n_el=n_el,
            tspan=tspan,
            tspan_bg=tspan_bg,
        )

    @deprecate(
        "0.2.13",
        "Use `ECMSCalulator.ecms_calibration_curve` instead.",
        "0.3",
    )
    def ecms_calibration_curve(self, mol, mass, n_el, *args, **kwargs):
        return ECMSCalculator.ecms_calibration_curve(
            measurement=self, mol=mol, mass=mass, n_el=n_el, *args, **kwargs
        )

    @deprecate(
        "0.2.13",
        "Use `plugins.siq.Calculator.ecms_calibration` instead.",
        "0.3",
    )
    def siq_ecms_calibration(self, mol, mass, n_el, tspan, tspan_bg=None):
        return plugins.siq.Calculator.ecms_calibration(
            measurement=self,
            mol=mol,
            mass=mass,
            n_el=n_el,
            tspan=tspan,
            tspan_bg=tspan_bg,
        )

    @deprecate(
        "0.2.13",
        "Use `plugins.siq.Calculator.ecms_calibration_curve` instead.",
        "0.3",
    )
    def siq_ecms_calibration_curve(
        self,
        mol,
        mass,
        n_el,
        *args,
        **kwargs,
    ):
        return plugins.siq.Calculator.ecms_calibration_curve(
            measurement=self, mol=mol, mass=mass, n_el=n_el, *args, **kwargs
        )


class ECMSCyclicVoltammogram(CyclicVoltammogram, ECMSMeasurement):
    """Class for raw EC-MS functionality. Parents: CyclicVoltammogram, ECMSMeasurement"""


class ECMSCalculator(ECCalibration, MSCalculator):
    """Class for calibrations useful for ECMSMeasurements"""

    extra_column_attrs = {
        "ecms_calibrations": {"date", "setup", "RE_vs_RHE", "A_el", "L"}
    }
    # FIXME: The above should be covered by the parent classes. Needs metaprogramming!
    # NOTE: technique, name, and tstamp in column_attrs are inherited from Calibration
    # NOTE: ms_results_ids in extra_linkers is inherited from MSCalibration.
    # NOTE: signal_bgs is left out

    def __init__(
        self,
        name=None,
        measurement=None,
        date=None,
        tstamp=None,
        setup=None,
        ms_cal_results=None,
        signal_bgs=None,
        RE_vs_RHE=None,
        A_el=None,
        R_Ohm=None,
        L=None,
        technique="EC-MS",
    ):
        """
        Args:
            name (str): Name of the ms_calibration
            date (str): Date of the ms_calibration
            setup (str): Name of the setup where the ms_calibration is made
            ms_cal_results (list of MSCalResult): The mass spec calibrations
            RE_vs_RHE (float): the RE potential in [V]
            A_el (float): The geometric electrode area in [cm^2]
            R_Ohm (float): The Ohmic drop in [Ohm]
            L (float): The working distance in [m]
        """
        ECCalibration.__init__(
            self,
            A_el=A_el,
            RE_vs_RHE=RE_vs_RHE,
            R_Ohm=R_Ohm,
        )
        MSCalculator.__init__(
            self,
            name=name,
            date=date,
            tstamp=tstamp,
            measurement=measurement,
            setup=setup,
            ms_cal_results=ms_cal_results,
            signal_bgs=signal_bgs,
        )
        self.technique = technique
        self.L = L

    @classmethod
    def ecms_calibration(cls, measurement, mol, mass, n_el, tspan, tspan_bg=None):
        """Calibrate for mol and mass based on one period of steady electrolysis

        Args:
            measurement (ECMSMeasurement): measurement with calibration data
            mol (str): Name of the molecule to calibrate
            mass (str): Name of the mass at which to calibrate
            n_el (str): Number of electrons passed per molecule produced (remember the
                sign! e.g. +4 for O2 by OER and -2 for H2 by HER)
            tspan (tspan): The timespan of steady electrolysis
            tspan_bg (tspan): The time to use as a background

        Return MSCalResult: The result of the ecms_calibration
        """
        if plugins.use_siq:
            warnings.warn(
                "spectro_inlets_quantification is active but you are using the native "
                "ixdat version of `ECMSMeasurement.ecms_calibration`"
            )
        Y = measurement.integrate_signal(mass, tspan=tspan, tspan_bg=tspan_bg)
        Q = measurement.integrate("raw_current", tspan=tspan) * 1e-3
        n = Q / (n_el * FARADAY_CONSTANT)
        F = Y / n
        cal = MSCalResult(
            name=f"{mol}@{mass}",
            mol=mol,
            mass=mass,
            cal_type="ecms_calibration",
            F=F,
        )
        return cls(ms_cal_results=[cal], measurement=measurement)

    @classmethod
    def ecms_calibration_curve(
        cls,
        measurement,
        mol,
        mass,
        n_el,
        tspan_list=None,
        selector_name=None,
        selector_list=None,
        t_steady_pulse=None,
        tspan_bg=None,
        force_through_zero=False,
        ax="new",
        axes_measurement=None,
        axes_measurement_J_name="raw_current",
        return_ax=False,
    ):
        """Fit mol's sensitivity at mass based on steady periods of EC production.

        Args:
            mol (str): Name of the molecule to calibrate
            mass (str): Name of the mass at which to calibrate
            n_el (str): Number of electrons passed per molecule produced (remember the
                sign! e.g. +4 for O2 by OER and -2 for H2 by HER)
            tspan_list (list of tspan): The timespans of steady electrolysis
            selector_name (str): Name of selector which identifies the periods
                of steady electrolysis for automatic selection of timespans of steady
                electrolysis. E.g. "selector" or "Ns" for biologic EC data
            selector_list (list): List of values for selector_name for automatic
                selection of timespans of steady electrolysis
            t_steady_pulse (float): Length of steady electrolysis for each segment
                given by selector_list. Defaults to None = entire length of segment
            tspan_bg (tspan): The time to use as a background
            force_through_zero (boolean): Whether to force the calibration curve through
                zero. This can be done when confident in the background subtraction.
            ax (Axis): The axis on which to plot the ms_calibration curve result.
                Defaults to a new axis.
            axes_measurement (list of Axes): The EC-MS plot axes to highlight the
                ms_calibration on. Defaults to None. These axes are not returned.
            axes_measurement_J_name (str): The J_name used in the axis passed
                to axes_measurement. Must be passed manually as the axis does not "know"
                its J_name. Defaults to "raw_current". IMPORTANT: the method still uses
                "raw_current" to calculate the sensitivity factor, this J_name is only
                used for plotting.
            return_ax (bool): Whether to return the axis on which the calibration curve
                is plotted together with the MSCalResult. Defaults to False.

        Return MSCalResult(, Axis): The result of the ms_calibration (and calibration
            curve axis if requested) based on integration of selected time periods.
        """

        axis_ms = axes_measurement[0] if axes_measurement else None
        axis_current = axes_measurement[3] if axes_measurement else None
        Y_list = []
        n_list = []
        if not tspan_list:
            tspan_list = measurement._get_tspan_list(
                selector_list, selector_name, t_steady_pulse
            )
        for tspan in tspan_list:
            Y = measurement.integrate_signal(
                mass, tspan=tspan, tspan_bg=tspan_bg, ax=axis_ms
            )
            # FIXME: plotting current by giving integrate() an axis doesn't work great.
            if (
                axes_measurement
            ):  # FIXME: need to run twice, once to plot, once to calculate Q
                measurement.integrate(
                    axes_measurement_J_name, tspan=tspan, ax=axis_current
                )
            Q = measurement.integrate("raw_current", tspan=tspan)
            Q *= 1e-3  # mC --> [C]
            n = Q / (n_el * FARADAY_CONSTANT)
            Y_list.append(Y)
            n_list.append(n)
        n_vec = np.array(n_list)
        Y_vec = np.array(Y_list)
        n_fit = np.array([0, max(n_vec)])
        if force_through_zero:

            def rms_error(F_guess):
                return np.mean((Y_vec - F_guess * n_vec) ** 2)

            F_guess_0 = np.sum(Y_vec) / np.sum(n_vec)
            res = minimize(rms_error, F_guess_0)
            F = res.x[0]
            Y_fit = n_fit * F
        else:
            pfit = np.polyfit(n_vec, Y_vec, deg=1)
            F = pfit[0]
            Y_fit = n_fit * pfit[0] + pfit[1]

        if ax:
            color = STANDARD_COLORS[mass]
            if ax == "new":
                ax = measurement.plotter.new_ax()
                ax.set_xlabel("amount produced / [nmol]")
                ax.set_ylabel("integrated signal / [nC]")
            ax.plot(n_vec * 1e9, Y_vec * 1e9, "o", color=color)
            ax.plot(n_fit * 1e9, Y_fit * 1e9, "--", color=color)

        cal = MSCalResult(
            name=f"{mol}@{mass}",
            mol=mol,
            mass=mass,
            cal_type="ecms_calibration_curve",
            F=F,
        )
        return cls(ms_cal_results=[cal], measurement=measurement)

    @property
    def available_series_names(self):
        names = set([f"n_dot_{mol}" for mol in self.mol_list])
        names.add("potential")
        names.add("current")
        return names

    def calculate_series(self, key, measurement=None):
        measurement = measurement or self.measurement
        try_1 = ECCalibration.calculate_series(self, key, measurement)
        if try_1:
            return try_1
        try_2 = MSCalculator.calculate_series(self, key, measurement)
        if try_2:
            return try_2


class ECMSSpectroMeasurement(ECMSMeasurement, MSSpectroMeasurement):
    pass
