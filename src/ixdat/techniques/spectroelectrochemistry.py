from .ec import ECMeasurement
from ..spectra import Spectrum


class SpectroECMeasurement(ECMeasurement):
    @property
    def whitelight_spectrum(self):
        return Spectrum.from_field(self["white_light"])
