from .base_mpl_plotter import MPLPlotter


class SpectrumPlotter(MPLPlotter):
    def __init__(self, spectrum=None):
        self.spectrum = spectrum

    def plot(self, spectrum=None, ax=None, **kwargs):
        spectrum = spectrum or self.spectrum
        if not ax:
            ax = self.new_ax()
        ax.plot(spectrum.x, spectrum.y, **kwargs)
        ax.set_xlabel(spectrum.x_name)
        ax.set_ylabel(spectrum.y_name)
        return ax
