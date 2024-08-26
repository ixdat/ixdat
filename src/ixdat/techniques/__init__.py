"""Import techniques and build the technique_classes dictionary for direct import

Constants:
    TECHNIQUE_CLASSES (dict): Dictionary of {technique_name: technique_class} where
        technique_name is the name of the technique (like "EC") and technique_class
        is the technique class (inheriting from Measurement) which implements the
        technique-specific functionality.
"""

from .ec import ECMeasurement
from .cv import CyclicVoltammogram, CyclicVoltammagram  # The latter is deprecated.
from .ms import MSMeasurement, MSSpectrum, MSSpectrumSeries, MSSpectroMeasurement
from .ec_ms import ECMSMeasurement, ECMSSpectroMeasurement
from .spectroelectrochemistry import (
    SpectroECMeasurement,
    ECXASMeasurement,
    ECOpticalMeasurement,
)

from .xrf import TRXRFMeasurement, ECTRXRFMeasurement
from .reactor import ReactorMeasurement, ReactorSpectroMeasurement, ReactorCalibration
from .ftir import FTIRSpectrum, ECFTIRMeasurement
from ..spectra import Spectrum
from ..measurements import Measurement


TECHNIQUE_CLASSES = {
    "simple": Measurement,
    "EC": ECMeasurement,
    "CV": CyclicVoltammogram,
    "MS": MSMeasurement,
    "EC-MS": ECMSMeasurement,
    "XRD": Spectrum,
    "XPS": Spectrum,
    "XAS": Spectrum,
    "MS_spectrum": MSSpectrum,
    "MS_spectra": MSSpectrumSeries,
    "SEC": SpectroECMeasurement,
    "EC-Optical": ECOpticalMeasurement,
    "EC-XAS": ECXASMeasurement,
    "MS-MS_spectra": MSSpectroMeasurement,
    "reactor": ReactorMeasurement,
    "reactor-MS_spectra": ReactorSpectroMeasurement,
    "S-EC": SpectroECMeasurement,
    "EC-MS-MS_spectra": ECMSSpectroMeasurement,
    "FTIR": FTIRSpectrum,
    "EC-FTIR": ECFTIRMeasurement,
}
