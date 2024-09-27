# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:59:53 2023

@author: Søren
"""

from ..spectra import Spectrum, SpectrumSeries
from .spectroelectrochemistry import SpectroECMeasurement
from ..plotters.spectrum_plotter import SpectrumSeriesPlotter
from ..plotters.sec_plotter import SECPlotter


class FTIRSpectrum(Spectrum):
    """FTIR Spectrum"""

    pass  # Just a generic Spectrum for now.


class FTIRSpectrumSeries(SpectrumSeries):
    """FTIR Spectrum Series"""

    # A generic SpectrumSeries but with a custom plotter, set in the __init__.
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.plotter = SpectrumSeriesPlotter(spectrum_series=self)
        self.plot = self.plotter.plot_stacked_spectra
        self.plot_stacked_spectra = self.plotter.plot_stacked_spectra
        self.heat_plot = self.plotter.heat_plot
        self.plot_waterfall = self.plotter.plot_waterfall


class ECFTIRMeasurement(SpectroECMeasurement):
    """FTIR Spectrum Series"""

    # A generic SpectroECMeasurement but with a custom plotter, set in the __init__.
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.plotter = SECPlotter(measurement=self)
        self.plot = self.plotter.plot_stacked_spectra
        self.plot_stacked_spectra = self.plotter.plot_stacked_spectra
        self.plot_measurement = self.plotter.plot_measurement
