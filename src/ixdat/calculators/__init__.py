# -*- coding: utf-8 -*-

from .indexer import Indexer
from .ec_calculators import ECCalibration, ScanRateCalculator
from .ms_calculators import MSBackgroundSet, MSCalibration
from .ecms_calculators import ECMSCalibration
from .spectral_calculators import Surfer, DeltaODCalculator


CALCULATOR_CLASSES = {
    cls.calculator_type: cls
    for cls in [
        Indexer,
        ECCalibration,
        ScanRateCalculator,
        MSBackgroundSet,
        MSCalibration,
        DeltaODCalc,
        Surfer
    ]
}
