"""NMR technique classes."""

from ..spectra import Spectrum, SpectrumSeries
from ..plotters.nmr_plotter import NMRPlotter
from ..plotters.spectrum_plotter import SpectrumSeriesPlotter


class NMRSpectrum(Spectrum):
    """A single NMR spectrum (1D).

    The x-axis is chemical shift in ppm and the y-axis is signal intensity.
    Acquisition and processing parameters are stored in :attr:`metadata`.

    Chemical shift (ppm) expresses a nucleus's resonance frequency relative to a
    reference compound: delta = (freq_sample - freq_ref) / freq_ref * 1e6. Using
    parts per million makes peak positions independent of the spectrometer field
    strength, so the same compound looks the same whether measured at 400 MHz or
    700 MHz. By NMR convention the ppm axis runs right-to-left (high ppm on the left).
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.plotter = NMRPlotter(spectrum=self)
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
