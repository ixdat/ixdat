"""Import techniques and build the technique_classes dictionary for direct import

Constants:
    technique_classes (dict): Dictionary of {technique_name: technique_class} where
        technique_name is the name of the technique (like "EC") and technique_class
        is the technique class (inheriting from Measurement) which implements the
        technique-specific functionality.
"""

from .ec import ECMeasurement
from .cv import CyclicVoltammagram
from .ms import MSMeasurement
from .ec_ms import ECMSMeasurement
from ..measurements import Measurement  # for importing in the technique modules

# TODO: Is something like DecoMeasurement a Measurement or something else?

TECHNIQUE_CLASSES = {
    "simple": Measurement,
    "EC": ECMeasurement,
    "CV": CyclicVoltammagram,
    "MS": MSMeasurement,
    "EC-MS": ECMSMeasurement,
}
