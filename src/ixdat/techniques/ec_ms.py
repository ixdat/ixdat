"""Module for representation and analysis of EC-MS measurements"""
import warnings
from .ec import ECMeasurement
from .ms import MSMeasurement, MSSpectroMeasurement
from .cv import CyclicVoltammogram
from ..exporters import ECMSExporter
from ..plotters import ECMSPlotter
from ..plugins import plugins
from ..tools import deprecate
from ..calculators.ecms_calculators import ECMSCalibration


class ECMSMeasurement(ECMeasurement, MSMeasurement):
    """Class for raw EC-MS functionality. Parents: ECMeasurement and MSMeasurement"""

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
        "0.3.1",
        kwarg_name="ms_cal_results",
    )
    def calibrate(self, *args, **kwargs):
        ms_cal_results = kwargs.pop("ms_cal_results", None)
        if ms_cal_results:
            from ..calculators.ms_calculators import MSCalibration

            warnings.warn(
                "Giving `ms_cal_results` to `ECMSMeasurement.calibrate` is DEPRECATED"
                "and will give an error in 0.3.1. Use instead:\n"
                "`ecms.add_calculator(MSCalculator(ms_cal_results=ms_cal_results))`"
            )
            cal = MSCalibration(ms_cal_results=ms_cal_results)
            self.add_calculator(cal)
            return cal

        return super().calibrate(*args, **kwargs)

    @deprecate(
        "0.2.13",
        "Use `ECMSCalulator.ecms_calibration` instead.",
        "0.3.1",
    )
    def ecms_calibration(self, mol, mass, n_el, tspan, tspan_bg=None):
        return ECMSCalibration.ecms_calibration(
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
        "0.3.1",
    )
    def ecms_calibration_curve(self, mol, mass, n_el, *args, **kwargs):
        return ECMSCalibration.ecms_calibration_curve(
            measurement=self, mol=mol, mass=mass, n_el=n_el, *args, **kwargs
        )

    @deprecate(
        "0.2.13",
        "Use `plugins.siq.Calculator.ecms_calibration` instead.",
        "0.3.1",
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
        "0.3.1",
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


class ECMSSpectroMeasurement(ECMSMeasurement, MSSpectroMeasurement):
    pass
