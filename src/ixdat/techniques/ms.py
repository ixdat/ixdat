"""Module for representation and analysis of MS measurements"""

from ..measurements import Measurement
from ..plotters.ms_plotter import MSPlotter
from ..exceptions import SeriesNotFoundError
import re
import numpy as np


class MSMeasurement(Measurement):
    """Class implementing raw MS functionality"""

    extra_column_attrs = {
        "ms_meaurements": {
            "mass_aliases",
        }
    }

    def __init__(self, name, mass_aliases=None, **kwargs):
        """Initializes a MS Measurement

        Args:
            name (str): The name of the measurement
            calibration (dict): calibration constants whereby the key
                corresponds to the respective signal name.
            mass_aliases (dict): {mass: mass_name} for any masses that
                do not have the standard 'M<x>' format used by ixdat.
        """
        super().__init__(name, **kwargs)
        self.calibration = None  # TODO: Not final implementation
        self.mass_aliases = mass_aliases or {}

    def __getitem__(self, item):
        """Adds to Measurement's lookup to check if item is an alias for a mass"""
        try:
            return super().__getitem__(item)
        except SeriesNotFoundError:
            if item in self.mass_aliases:
                return self[self.mass_aliases[item]]
            else:
                raise

    def grab_signal(self, signal_name, tspan=None, t_bg=None):
        """Returns raw signal for a given signal name

        Args:
            signal_name (str): Name of the signal.
            tspan (list): Timespan for which the signal is returned.
            t_bg (list): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
        """
        time, value = self.grab(signal_name, tspan=tspan)

        if t_bg is None:
            return time, value

        else:
            _, bg = self.grab(signal_name, tspan=t_bg)
            return time, value - np.average(bg)

    def grab_cal_signal(self, signal_name, tspan=None, t_bg=None):
        """Returns a calibrated signal for a given signal name. Only works if
        calibration dict is not None.

        Args:
            signal_name (str): Name of the signal.
            tspan (list): Timespan for which the signal is returned.
            t_bg (list): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
        """
        # TODO: Not final implementation
        if self.calibration is None:
            print("No calibration dict found.")
            return

        time, value = self.grab_signal(signal_name, tspan=tspan, t_bg=t_bg)

        return time, value * self.calibration[signal_name]

    @property
    def mass_list(self):
        """List of the masses for which ValueSeries are contained in the measurement"""
        return [self.as_mass(col) for col in self.series_names if self.is_mass(col)]

    def is_mass(self, item):
        if re.search("^M[0-9]+$", item):
            return True
        if item in self.mass_aliases.values():
            return True
        return False

    def as_mass(self, item):
        if re.search("^M[0-9]+$", item):
            return item
        else:
            try:
                return next(k for k, v in self.mass_aliases.items() if v == item)
            except StopIteration:
                raise TypeError(f"{self} does not recognize '{item}' as a mass.")

    @property
    def plotter(self):
        """The default plotter for ECMeasurement is ECPlotter"""
        if not self._plotter:
            from ..plotters.ec_plotter import ECPlotter

            self._plotter = MSPlotter(measurement=self)

        return self._plotter
