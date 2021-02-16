"""Module for representation and analysis of EC-MS measurements"""
from .ec import ECMeasurement
from .ms import MSMeasurement


class ECMSMeasurement(ECMeasurement, MSMeasurement):
    """Class implementing raw EC-MS functionality"""

    def __init__(self, name, **kwargs):
        super().__init__(name, **kwargs)

    @property
    def plotter(self):
        """The default plotter for ECMeasurement is ECPlotter"""
        if not self._plotter:
            from ..plotters.ecms_plotter import ECMSPlotter

            self._plotter = ECMSPlotter(measurement=self)

        return self._plotter
