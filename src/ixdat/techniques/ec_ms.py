"""Module for representation and analysis of EC-MS measurements"""
from .ec import ECMeasurement
from .ms import MSMeasurement
from .cv import CyclicVoltammagram


class ECMSMeasurement(ECMeasurement, MSMeasurement):
    """Class for raw EC-MS functionality. Parents: ECMeasurement and MSMeasurement"""

    @property
    def plotter(self):
        """The default plotter for ECMSMeasurement is ECMSPlotter"""
        if not self._plotter:
            from ..plotters.ecms_plotter import ECMSPlotter

            self._plotter = ECMSPlotter(measurement=self)

        return self._plotter

    def as_cv(self):
        self_as_dict = self.as_dict()

        # FIXME: The following lines are only necessary because
        #  PlaceHolderObject.get_object isn't able to find things in the MemoryBackend
        del self_as_dict["s_ids"]
        self_as_dict["series_list"] = self.series_list

        return ECMSCyclicVoltammogram.from_dict(self_as_dict)


class ECMSCyclicVoltammogram(CyclicVoltammagram, MSMeasurement):
    """Class for raw EC-MS functionality. Parents: CyclicVoltammogram, MSMeasurement"""

    @property
    def plotter(self):
        """The default plotter for ECMSCyclicVoltammogram is ECMSPlotter"""
        if not self._plotter:
            from ..plotters.ecms_plotter import ECMSPlotter

            self._plotter = ECMSPlotter(measurement=self)

        return self._plotter
