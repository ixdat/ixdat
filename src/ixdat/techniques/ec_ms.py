"""Module for representation and analysis of EC-MS measurements"""
from .ec import ECMeasurement
from .ms import MSMeasurement


class ECMSMeasurement(ECMeasurement, MSMeasurement):
    """Class for raw EC-MS functionality. Parents: ECMeasurement and MSMeasurement"""

    def __init__(self, name, **kwargs):
        super().__init__(name, **kwargs)  # FIXME: This seems to just be ECMeasurement.

    @property
    def plotter(self):
        """The default plotter for ECMSMeasurement is ECMSPlotter"""
        if not self._plotter:
            from ..plotters.ecms_plotter import ECMSPlotter

            self._plotter = ECMSPlotter(measurement=self)

        return self._plotter
