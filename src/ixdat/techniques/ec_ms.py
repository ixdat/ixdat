"""Module for representation and analysis of EC-MS measurements"""
import re

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

    @property
    def mass_list(self):
        return [col for col in self.series_names if re.search("^M[0-9]+$", col)]
