from ..measurement_base import Calculator


class Surfer(Calculator):
    pass


class ODCalculator:
    def __init__(self, spectrum_series=None, reference_spectrum=None):
        self.spectrum_series = spectrum_series
        self.reference_spectrum = reference_spectrum


class SECCalculator:
    def __init__(self, spectrum_series=None, reference_spectrum=None):
        self.spectrum_series = spectrum_series
        self.reference_spectrum = reference_spectrum

