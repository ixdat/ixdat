from .ms import MSMeasurement
from ..plotters.tpms_plotter import TPMSPlotter


class ReactorMeasurement(MSMeasurement):

    default_plotter = TPMSPlotter
    essential_series_names = ("temperature", "pressure")

    @property
    def T_name(self):
        return self["temperature"].name

    @property
    def P_name(self):
        return self["pressure"].name

    @property
    def t_name(self):
        return self["temperature"].tseries.name

    @property
    def T(self):
        return self["temperature"].data

    @property
    def P(self):
        return self["pressure"].data

    @property
    def t(self):
        return self["temperature"].t
