"""NMR technique classes."""

from ..spectra import Spectrum, SpectrumSeries
from ..plotters.spectrum_plotter import SpectrumPlotter, SpectrumSeriesPlotter


class NMRSpectrum(Spectrum):
    """A single NMR spectrum (1D).

    The x-axis is chemical shift in ppm and the y-axis is signal intensity.
    Acquisition and processing parameters are stored in :attr:`metadata`.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.plotter = SpectrumPlotter(spectrum=self)
        self.plot = self.plotter.plot


class NMRSpectrumSeries(SpectrumSeries):
    """A series of NMR spectra sharing a common chemical-shift axis."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.plotter = SpectrumSeriesPlotter(spectrum_series=self)
        self.plot = self.plotter.plot_stacked_spectra
        self.plot_stacked_spectra = self.plotter.plot_stacked_spectra
        self.heat_plot = self.plotter.heat_plot
        self.plot_waterfall = self.plotter.plot_waterfall
