"""Unit tests for ixdat.spectra"""

import numpy as np

from ixdat.spectra import Spectrum, SpectrumSeries
from ixdat.techniques.ms import MSSpectrumSeries


def test_indexing_ms_spectrum_series_returns_a_spectrum():
    """Indexing a series whose technique is registered to a SpectrumSeries class
    returns a single Spectrum."""
    x = np.linspace(1, 50, 50)
    spectrum_list = [
        Spectrum.from_data(x, x * i, tstamp=1e9 + i, x_name="m/z", y_name="signal")
        for i in range(3)
    ]
    series = MSSpectrumSeries.from_spectrum_list(spectrum_list, technique="MS_spectra")

    spectrum = series[1]

    assert not isinstance(spectrum, SpectrumSeries)
    assert np.array_equal(spectrum.y, x * 1)
