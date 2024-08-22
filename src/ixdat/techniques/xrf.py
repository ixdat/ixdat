# -*- coding: utf-8 -*-
"""Module for representation and analysis of TRXRF measurements"""

from ..measurements import Measurement
from ..calculators.xrf_calculators import TRXRFCalculator
from ..plotters.xrf_plotter import TRXRFPlotter, ECTRXRFPlotter
from .ec import ECMeasurement


class TRXRFMeasurement(Measurement):
    """Time-resolved X-ray Fluoresence."""

    default_plotter = TRXRFPlotter
    default_calibration = TRXRFCalculator


class ECTRXRFMeasurement(ECMeasurement, TRXRFMeasurement):
    """Electrochemistry with time-resolved X-ray Fluoresence."""

    default_plotter = ECTRXRFPlotter
