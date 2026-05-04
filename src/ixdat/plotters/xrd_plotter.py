"""Plotter for XRD spectra."""

from .base_mpl_plotter import MPLPlotter


class XRDSpectrumPlotter(MPLPlotter):
    """Plots XRD intensity vs x, with a shaded error band for .xye data."""

    def __init__(self, spectrum=None):
        super().__init__()
        self.spectrum = spectrum

    def plot(self, *, spectrum=None, ax=None, **kwargs):
        """Plot intensity vs x with a shaded y+/-e band if error data is present.

        Args:
            spectrum (XRDSpectrum): Defaults to self.spectrum.
            ax (mpl.Axis): Axis to plot on. A new one is made by default.
            kwargs: Passed to ax.plot() for the intensity line.
        """
        spectrum = spectrum or self.spectrum
        if not ax:
            ax = self.new_ax()
        x = spectrum.x
        y = spectrum.y
        ax.plot(x, y, **kwargs)
        ax.set_xlabel(spectrum.x_name)
        ax.set_ylabel(spectrum.y_name)
        if spectrum.y_err is not None:
            e = spectrum.y_err
            color = kwargs.get("color", None)
            ax.fill_between(x, y - e, y + e, alpha=0.3, color=color)
        return ax
