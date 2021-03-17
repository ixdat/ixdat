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
            "signal_bgs",
        },
    }

    def __init__(
        self, name, mass_aliases=None, signal_bgs=None, tspan_bg=None, **kwargs
    ):
        """Initializes a MS Measurement

        Args:
            name (str): The name of the measurement
            calibration (dict): calibration constants whereby the key
                corresponds to the respective signal name.
            mass_aliases (dict): {mass: mass_name} for any masses that
                do not have the standard 'M<x>' format used by ixdat.
            signal_bgs (dict): {mass: S_bg} where S_bg is the background signal
                in [A] for the mass (typically set with a timespan by `set_bg()`)
            tspan_bg (timespan): background time used to set masses
        """
        super().__init__(name, **kwargs)
        self.calibration = None  # TODO: Not final implementation
        self.mass_aliases = mass_aliases or {}
        self.signal_bgs = signal_bgs or {}
        self.tspan_bg = tspan_bg

    def __getitem__(self, item):
        """Adds to Measurement's lookup to check if item is an alias for a mass"""
        try:
            return super().__getitem__(item)
        except SeriesNotFoundError:
            if item in self.mass_aliases:
                return self[self.mass_aliases[item]]
            else:
                raise

    def set_bg(self, tspan_bg=None, mass_list=None):
        """Set background values for mass_list to the average signal during tspan_bg."""
        mass_list = mass_list or self.mass_list
        tspan_bg = tspan_bg or self.tspan_bg
        for mass in mass_list:
            t, v = self.grab(mass, tspan_bg)
            self.signal_bgs[mass] = np.mean(v)

    def reset_bg(self, mass_list=None):
        """Reset background values for the masses in mass_list"""
        mass_list = mass_list or self.mass_list
        for mass in mass_list:
            if mass in self.signal_bgs:
                del self.signal_bgs[mass]

    def grab_signal(
        self,
        signal_name,
        tspan=None,
        t_bg=None,
        removebackground=False,
        include_endpoints=False,
    ):
        """Returns raw signal for a given signal name

        Args:
            signal_name (str): Name of the signal.
            tspan (list): Timespan for which the signal is returned.
            t_bg (list): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
            removebackground (bool): Whether to remove a pre-set background if available
        """
        time, value = self.grab(
            signal_name, tspan=tspan, include_endpoints=include_endpoints
        )

        if t_bg is None:
            if removebackground and signal_name in self.signal_bgs:
                return time, value - self.signal_bgs[signal_name]
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
        """The default plotter for MSMeasurement is MSPlotter"""
        if not self._plotter:
            self._plotter = MSPlotter(measurement=self)
        return self._plotter
