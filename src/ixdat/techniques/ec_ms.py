"""Module for representation and analysis of EC-MS measurements"""
import numpy as np
from ..constants import FARADAY_CONSTANT
from .ec import ECMeasurement
from .ms import MSMeasurement, MSCalResult
from .cv import CyclicVoltammagram
from ..exporters.ecms_exporter import ECMSExporter
from ..plotters.ms_plotter import STANDARD_COLORS


class ECMSMeasurement(ECMeasurement, MSMeasurement):
    """Class for raw EC-MS functionality. Parents: ECMeasurement and MSMeasurement"""

    extra_column_attrs = {
        # FIXME: It would be more elegant if this carried over from both parents
        #   That might require some custom inheritance definition...
        "ecms_meaurements": {
            "mass_aliases",
            "signal_bgs",
            "ec_technique",
            "RE_vs_RHE",
            "R_Ohm",
            "raw_potential_names",
            "A_el",
            "raw_current_names",
        },
    }

    def __init__(self, **kwargs):
        """FIXME: Passing the right key-word arguments on is a mess"""
        ec_kwargs = {
            k: v for k, v in kwargs.items() if k in ECMeasurement.get_all_column_attrs()
        }
        ms_kwargs = {
            k: v for k, v in kwargs.items() if k in MSMeasurement.get_all_column_attrs()
        }
        # FIXME: I think the lines below could be avoided with a PlaceHolderObject that
        #  works together with MemoryBackend
        if "series_list" in kwargs:
            ec_kwargs.update(series_list=kwargs["series_list"])
            ms_kwargs.update(series_list=kwargs["series_list"])
        if "component_measurements" in kwargs:
            ec_kwargs.update(component_measurements=kwargs["component_measurements"])
            ms_kwargs.update(component_measurements=kwargs["component_measurements"])
        ECMeasurement.__init__(self, **ec_kwargs)
        MSMeasurement.__init__(self, **ms_kwargs)

    @property
    def plotter(self):
        """The default plotter for ECMSMeasurement is ECMSPlotter"""
        if not self._plotter:
            from ..plotters.ecms_plotter import ECMSPlotter

            self._plotter = ECMSPlotter(measurement=self)

        return self._plotter

    @property
    def exporter(self):
        """The default plotter for ECMSMeasurement is ECMSExporter"""
        if not self._exporter:
            self._exporter = ECMSExporter(measurement=self)
        return self._exporter

    def as_cv(self):
        self_as_dict = self.as_dict()

        # FIXME: The following lines are only necessary because
        #  PlaceHolderObject.get_object isn't able to find things in the MemoryBackend
        del self_as_dict["s_ids"]
        self_as_dict["series_list"] = self.series_list

        return ECMSCyclicVoltammogram.from_dict(self_as_dict)

    def ecms_calibration_curve(
        self,
        mol,
        mass,
        n_el,
        tspan_list=None,
        tspan_bg=None,
        ax="new",
        axes_measurement=None,
    ):
        """Fit mol's sensitivity at mass based on steady periods of EC production

        Args:
            mol (str): Name of the molecule to calibrate
            mass (str): Name of the mass at which to calibrate
            n_el (str): Number of electrons passed per molecule produced (remember the
                sign! e.g. +4 for O2 by OER and -2 for H2 by HER)
            tspan_list (list of tspan): THe timespans of steady electrolysis
            tspan_bg (tspan): The time to use as a background
            ax (Axis): The axis on which to plot the calibration curve result. Defaults
                to a new axis.
            axes_measurement (list of Axes): The EC-MS plot axes to highlight the
                calibration on. Defaults to None.

        Return MSCalResult: The result of the calibration
        """
        axis_ms = axes_measurement[0] if axes_measurement else None
        axis_current = axes_measurement[0] if axes_measurement else None
        Y_list = []
        n_list = []
        for tspan in tspan_list:
            Y = self.integrate_signal(mass, tspan=tspan, tspan_bg=tspan_bg, ax=axis_ms)
            # FIXME: plotting current by giving integrate() an axis doesn't work great.
            I = self.integrate("raw current / [mA]", tspan=tspan) * 1e-3
            n = I / (n_el * FARADAY_CONSTANT)
            Y_list.append(Y)
            n_list.append(n)
        n_vec = np.array(n_list)
        Y_vec = np.array(Y_list)
        pfit = np.polyfit(n_vec, Y_vec, deg=1)
        F = pfit[0]
        if ax:
            color = STANDARD_COLORS[mass]
            if ax == "new":
                ax = self.plotter.new_ax()
                ax.set_xlabel("amount produced / [nmol]")
                ax.set_ylabel("integrated signal / [nC]")
            ax.plot(n_vec * 1e9, Y_vec * 1e9, "o", color=color)
            n_fit = np.array([0, max(n_vec)])
            Y_fit = n_fit * pfit[0] + pfit[1]
            ax.plot(n_fit * 1e9, Y_fit * 1e9, "--", color=color)
        cal = MSCalResult(
            name=f"{mol}_{mass}",
            mol=mol,
            mass=mass,
            cal_type="ecms_calibration_curve",
            F=F,
        )
        if ax:
            if axes_measurement:
                return cal, ax, axes_measurement
            return cal, ax
        return cal


class ECMSCyclicVoltammogram(CyclicVoltammagram, MSMeasurement):
    """Class for raw EC-MS functionality. Parents: CyclicVoltammogram, MSMeasurement

    FIXME: Maybe this class should instead inherit from ECMSMeasurement and
        just add the CyclicVoltammogram functionality?
    """

    extra_column_attrs = {
        # FIXME: It would be more elegant if this carried over from both parents
        #   That might require some custom inheritance definition...
        "ecms_meaurements": {
            "mass_aliases",
            "signal_bgs",
            "ec_technique",
            "RE_vs_RHE",
            "R_Ohm",
            "raw_potential_names",
            "A_el",
            "raw_current_names",
        },
    }

    def __init__(self, **kwargs):
        """FIXME: Passing the right key-word arguments on is a mess"""
        ec_kwargs = {
            k: v for k, v in kwargs.items() if k in ECMeasurement.get_all_column_attrs()
        }
        ec_kwargs.update(series_list=kwargs["series_list"])
        ECMeasurement.__init__(self, **ec_kwargs)
        ms_kwargs = {
            k: v for k, v in kwargs.items() if k in MSMeasurement.get_all_column_attrs()
        }
        ms_kwargs.update(series_list=kwargs["series_list"])
        MSMeasurement.__init__(self, **ms_kwargs)
        self.plot = self.plotter.plot_vs_potential

    @property
    def plotter(self):
        """The default plotter for ECMSCyclicVoltammogram is ECMSPlotter"""
        if not self._plotter:
            from ..plotters.ecms_plotter import ECMSPlotter

            self._plotter = ECMSPlotter(measurement=self)

        return self._plotter

    @property
    def exporter(self):
        """The default plotter for ECMSCyclicVoltammogram is ECMSExporter"""
        if not self._exporter:
            self._exporter = ECMSExporter(measurement=self)
        return self._exporter


class ECMSCalibration:
    pass
