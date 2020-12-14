from .ec import ECMeasurement
from .ms import MSMeasurement
from .ec_ms import ECMSMeasurement

technique_classes = {
    "EC": ECMeasurement,
    "MS": MSMeasurement,
    "EC-MS": ECMSMeasurement,
}
