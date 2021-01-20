"""Module for representation and analysis of EC measurements"""

from ..measurements import Measurement


class ECMeasurement(Measurement):
    """Class implementing raw electrochemistry measurements"""

    extra_column_attrs = {"ec_techniques": {"ec_technique"}}
    technique = "EC"

    def __init__(self, ec_technique=None, **kwargs):
        super().__init__(**kwargs)
        self.ec_technique = ec_technique
        self.E_str = "Ewe/V"
        self.I_str = "I/mA"

    def get_potential(self, tspan=None):
        return self.get_t_and_v(self.E_str, tspan=tspan)
