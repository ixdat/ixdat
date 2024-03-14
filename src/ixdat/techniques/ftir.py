# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:59:53 2023

@author: Søren
"""

from ..spectra import Spectrum, SpectrumSeries
from .spectroelectrochemistry import SpectroECMeasurement
from ..plotters.ftir_plotter import FTIRPlotter


class FTIRSpectrum(Spectrum):
    pass


class FTIRSpectrumSeries(SpectrumSeries):
    default_plotter = FTIRPlotter

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.plotter = FTIRPlotter(spectrum_series=self)
        self.plot = self.plotter.plot_stacked_spectra
        self.heat_plot = self.plotter.heat_plot
        self.plot_waterfall = self.plotter.plot_waterfall


class ECFTIRMeasurement(SpectroECMeasurement):
    pass
