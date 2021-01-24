"""Module for representation and analysis of EC-MS measurements"""

from .ec import ECMeasurement
from .ms import MSMeasurement


class ECMSMeasurement(ECMeasurement, MSMeasurement):
    """Class implementing raw EC-MS functionality"""

    def __init__(self, name, **kwargs):
        super().__init__(name, **kwargs)
