"""Module for representation and analysis of MS measurements"""

from ..measurements import Measurement
import numpy as np


class MSMeasurement(Measurement):
    """Class implementing raw MS functionality"""

    def __init__(self, name, **kwargs):
        """Initializes a MS Measurement

        Args:
            name (str): The name of the measurement
            calibration (dict): calibration constants whereby the key
                corresponds to the respective signal name."""
        super().__init__(name, **kwargs)
        self.calibration = None  # TODO: Not final implementation

    def get_signal(self, signal_name, tspan=None, t_bg=None):
        """Returns raw signal for a given signal name

        Args:
            signal_name (str): Name of the signal.
            tspan (list): Timespan for which the signal is returned.
            t_bg (list): Timespan that corresponds to the background signal.
                If not given, no background is subtracted.
        """
        time, value = self.get_t_and_v(signal_name, tspan=tspan)

        if t_bg is None:
            return time, value

        else:
            _, bg = self.get_t_and_v(signal_name, tspan=t_bg)
            return time, value - np.average(bg)

    def get_calib_signal(self, signal_name, tspan=None, t_bg=None):
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

        time, value = self.get_signal(signal_name, tspan=tspan, t_bg=t_bg)

        return time, value * self.calibration[signal_name]
