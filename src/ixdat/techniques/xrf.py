# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 16:58:47 2024

@author: cl2120
"""

from ..measurements import Measurement
from ..data_series import ValueSeries
from ..plotters.xrf_plotter import TRXRFPlotter, ECTRXRFPlotter
from .ec import ECMeasurement


class TRXRFMeasurement(Measurement):

    series_constructors = Measurement.series_constructors
    series_constructors.update({"FF_over_I0": "_build_FF_over_I0"})
    default_plotter = TRXRFPlotter

    def _build_FF_over_I0(self):
        return ValueSeries(
            name="FF_over_I0",
            unit_name="",
            data=self["FF"].data / self["I0"].data,
            tseries=self["FF"].tseries,
        )


class ECTRXRFMeasurement(ECMeasurement, TRXRFMeasurement):
    default_plotter = ECTRXRFPlotter
