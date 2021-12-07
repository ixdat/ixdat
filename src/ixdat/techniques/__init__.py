"""Import techniques and build the technique_classes dictionary for direct import

Constants:
    technique_classes (dict): Dictionary of {technique_name: technique_class} where
        technique_name is the name of the technique (like "EC") and technique_class
        is the technique class (inheriting from Measurement) which implements the
        technique-specific functionality.
"""

from .ec import ECMeasurement, ECCalibration
from .cv import CyclicVoltammagram
from .ms import MSMeasurement, MSCalibration
from .ec_ms import ECMSMeasurement, ECMSCalibration
from .spectroelectrochemistry import SpectroECMeasurement
from ..measurements import Measurement  # for importing in the technique modules

# TODO: Is something like DecoMeasurement a Measurement or something else?

TECHNIQUE_CLASSES = {
    "simple": Measurement,
    "EC": ECMeasurement,
    "CV": CyclicVoltammagram,
    "MS": MSMeasurement,
    "EC-MS": ECMSMeasurement,
    "S-EC": SpectroECMeasurement,
}

CALIBRATION_CLASSES = {
    "EC": ECCalibration,
    "MS": MSCalibration,
    "EC-MS": ECMSCalibration,
}
