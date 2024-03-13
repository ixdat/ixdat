"""Import techniques and build the technique_classes dictionary for direct import

Constants:
    technique_classes (dict): Dictionary of {technique_name: technique_class} where
        technique_name is the name of the technique (like "EC") and technique_class
        is the technique class (inheriting from Measurement) which implements the
        technique-specific functionality.
"""

from .ec import ECMeasurement, ECCalibration
from .cv import CyclicVoltammogram, CyclicVoltammagram  # The latter is deprecated.
from .ms import MSMeasurement, MSCalibration, MSSpectrum, MSSpectroMeasurement
from .ec_ms import ECMSMeasurement, ECMSCalibration, ECMSSpectroMeasurement
from .spectroelectrochemistry import (
    SpectroECMeasurement,
    ECXASMeasurement,
    ECOpticalMeasurement,
)
from .reactor import ReactorMeasurement, ReactorSpectroMeasurement, ReactorCalibration

from ..spectra import Spectrum
from ..measurements import Measurement  # for importing in the technique modules

# TODO: Is something like DecoMeasurement a Measurement or something else?

TECHNIQUE_CLASSES = {
    "simple": Measurement,
    "EC": ECMeasurement,
    "CV": CyclicVoltammogram,
    "MS": MSMeasurement,
    "EC-MS": ECMSMeasurement,
    "XRD": Spectrum,
    "XPS": Spectrum,
    "XAS": Spectrum,
    "MS_spectra": MSSpectrum,
    "SEC": SpectroECMeasurement,
    "EC-Optical": ECOpticalMeasurement,
    "EC-XAS": ECXASMeasurement,
    "MS-MS_spectra": MSSpectroMeasurement,
    "reactor": ReactorMeasurement,
    "reactor-MS_spectra": ReactorSpectroMeasurement,
    "S-EC": SpectroECMeasurement,
    "EC-MS-MS_spectra": ECMSSpectroMeasurement,
}

CALIBRATION_CLASSES = {
    "EC": ECCalibration,
    "CV": ECCalibration,
    "MS": MSCalibration,
    "EC-MS": ECMSCalibration,
    "reactor": ReactorCalibration,
}
