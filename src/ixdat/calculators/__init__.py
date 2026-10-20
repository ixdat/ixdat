# -*- coding: utf-8 -*-

from .indexer import Indexer
from .ec_calculators import ECCalibration, ScanRateCalculator
from .ms_calculators import MSBackgroundSet, MSCalibration
from .ecms_calculators import ECMSCalibration
from .sec_calculators import DeltaODCalculator, Denoiser, SpectrumDeconvolutor


CALCULATOR_CLASSES = {
    cls.calculator_type: cls
    for cls in [
        Indexer,
        ECCalibration,
        ScanRateCalculator,
        MSBackgroundSet,
        MSCalibration,
        DeltaODCalculator,
        Denoiser,
        SpectrumDeconvolutor,
    ]
}
