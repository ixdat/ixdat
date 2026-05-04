"""Module for representation and analysis of XRD spectra."""

from ..spectra import MultiSpectrum
from ..plotters.xrd_plotter import XRDSpectrumPlotter


class XRDSpectrum(MultiSpectrum):
    """XRD diffraction pattern with optional per-point intensity error.

    Returned by XRDXYReader for .xye files (three columns: x, intensity, error).
    The first field is the intensity, the second is the per-point error.
    Index by field name to get individual Spectrum objects, e.g. xrd["intensity"].
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.plotter = XRDSpectrumPlotter(spectrum=self)
        self.plot = self.plotter.plot

    @property
    def x_name(self):
        return self.xseries.name

    @property
    def y(self):
        return self.fields[0].data

    @property
    def y_name(self):
        return self.fields[0].name

    @property
    def y_err(self):
        if len(self.fields) >= 2:
            return self.fields[1].data
        return None
