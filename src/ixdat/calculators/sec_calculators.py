# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 14:51:49 2026

@author: Olivia
"""

import warnings
import numpy as np
import pandas as pd
import json
from scipy.optimize import curve_fit, minimize  # noqa
from scipy.interpolate import interp1d  # noqa
from scipy import signal  # noqa
from mpmath import invertlaplace, sinh, cosh, sqrt, exp, erfc, pi, tanh, coth  # noqa
from numpy.fft import fft, ifft, ifftshift, fftfreq  # noqa
from ..plotters.ms_plotter import STANDARD_COLORS
from ..constants import FARADAY_CONSTANT, R, STANDARD_TEMPERATURE, STANDARD_PRESSURE
from .ms_calculators import MSCalibration, MSCalResult
from .ec_calculators import ECCalibration
from ..measurement_base import Calculator
from ..data_series import ValueSeries, TimeSeries
from ..exceptions import TechniqueError, QuantificationError
from ..plugins import plugins
from ..plotters.plotting_tools import calc_linear_background
from ..tools import deprecate


class Spectroelectrochemistry(Calculator):
    """Class for calculations done using ECOpticalMeasurements

    Its classmethods return SECCalibration objects.
    """