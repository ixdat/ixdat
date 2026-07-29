"""Plotter for NMR spectra."""

from .spectrum_plotter import SpectrumPlotter
from .renderers import _is_matplotlib_backend


class NMRPlotter(SpectrumPlotter):
    """A plotter for NMR spectra. Inverts the x-axis so high ppm is on the left,
    matching the standard NMR convention."""

    def plot(self, *, spectrum=None, ax=None, backend=None, figure=None, **kwargs):
        if _is_matplotlib_backend(backend):
            ax = super().plot(spectrum=spectrum, ax=ax, **kwargs)
            ax.invert_xaxis()
            return ax
        return super().plot(
            spectrum=spectrum,
            ax=ax,
            backend=backend,
            figure=figure,
            inverted_x=True,
            **kwargs,
        )
