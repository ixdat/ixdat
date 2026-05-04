"""Plotter for NMR spectra."""

from .spectrum_plotter import SpectrumPlotter


class NMRPlotter(SpectrumPlotter):
    """A plotter for NMR spectra. Inverts the x-axis so high ppm is on the left,
    matching the standard NMR convention."""

    def plot(self, *, spectrum=None, ax=None, **kwargs):
        ax = super().plot(spectrum=spectrum, ax=ax, **kwargs)
        ax.invert_xaxis()
        return ax
