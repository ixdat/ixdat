# -*- coding: utf-8 -*-
"""Module for representation and analysis of TRXRF measurements"""

from ..measurements import Measurement
from ..data_series import ValueSeries
from ..plotters.xrf_plotter import TRXRFPlotter, ECTRXRFPlotter
from .ec import ECMeasurement


class TRXRFMeasurement(Measurement):
    """Time-resolved X-ray Fluoresence."""
    series_constructors = Measurement.series_constructors
    series_constructors.update({"FF_over_I0": "_build_FF_over_I0"})
    default_plotter = TRXRFPlotter

    def _build_FF_over_I0(self):
        """Construct and return a new series of FF/I0, using FF and I0 data series in raw data file"""
        
        return ValueSeries(
            name="FF_over_I0",
            unit_name="",
            data=self["FF"].data / self["I0"].data,
            tseries=self["FF"].tseries,
        )
        

class ECTRXRFMeasurement(ECMeasurement, TRXRFMeasurement):
    """Electrochemistry with time-resolved X-ray Fluoresence."""
    
    default_plotter = ECTRXRFPlotter
