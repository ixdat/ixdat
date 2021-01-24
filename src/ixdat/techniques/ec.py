"""Module for representation and analysis of EC measurements"""

from ..measurements import Measurement


class ECMeasurement(Measurement):
    """Class implementing raw electrochemistry measurements"""

    def __init__(self, name, ec_technique=None, **kwargs):
        super().__init__(name, **kwargs)
        self.ec_technique = ec_technique
        self.E_str = "Ewe/V"
        self.I_str = "I/mA"

    def get_potential(self, tspan=None):
        """Returns measured electrochemical potential

        Args:
            tspan (list): Timespan for which the potential is returned.
        """
        return self.get_t_and_v(self.E_str, tspan=tspan)

    def get_current(self, tspan=None):
        """Returns measured electrochemical current

        Args:
            tspan (list): Timespan for which the current is returned.
        """
        return self.get_t_and_v(self.I_str, tspan=tspan)
