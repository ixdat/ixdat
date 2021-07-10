import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
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


class SpectrumSeriesPlotter(MPLPlotter):
    def __init__(self, spectrum_series=None):
        self.spectrum_series = spectrum_series

    @property
    def plot(self):
        return self.heat_plot

    def plot_average(self, spectrum_series=None, ax=None, **kwargs):
        spectrum_series = spectrum_series or self.spectrum_series
        if not ax:
            ax = self.new_ax()
        ax.plot(spectrum_series.x, spectrum_series.y_average, **kwargs)
        ax.set_xlabel(spectrum_series.x_name)
        ax.set_ylabel(spectrum_series.y_name + " (average)")
        return ax

    def heat_plot(
        self,
        spectrum_series=None,
        tspan=None,
        xspan=None,
        ax=None,
        cmap_name="inferno",
        make_colorbar=False,
    ):
        field = spectrum_series.field
        tseries = field.axes_series[0]
        t = tseries.data
        xseries = field.axes_series[1]
        x = xseries.data
        data = field.data

        if tspan:
            t_mask = np.logical_and(tspan[0] < t, t < tspan[-1])
            t = t[t_mask]
            data = data[t_mask, :]
        if xspan:
            wl_mask = np.logical_and(xspan[0] < x, x < xspan[-1])
            x = x[wl_mask]
            data = data[:, wl_mask]

        ax.imshow(
            data.swapaxes(0, 1),
            cmap=cmap_name,
            aspect="auto",
            extent=(t[0], t[-1], x[0], x[-1]),
        )

        ax.set_xlabel(
            (spectrum_series.t_str if hasattr(spectrum_series, "t_str") else None)
            or tseries.name
        )
        ax.set_ylabel(xseries.name)
        if make_colorbar:
            cmap = mpl.cm.get_cmap(cmap_name)
            norm = mpl.colors.Normalize(vmin=np.min(data), vmax=np.max(data))
            cb = plt.colorbar(
                mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                ax=ax,
                use_gridspec=True,
                anchor=(0.75, 0),
            )
            cb.set_label("$\Delta$ opdical density")
        return ax

    def plot_waterfall(
        self, spectrum_series=None, cmap_name="jet", make_colorbar=True, ax=None
    ):
        spectrum_series = spectrum_series or self.spectrum_series
        if not ax:
            ax = self.new_ax()
        field = spectrum_series.field
        data = field.data
        t = spectrum_series.t
        wl = spectrum_series.wl

        cmap = mpl.cm.get_cmap(cmap_name)
        norm = mpl.colors.Normalize(vmin=np.min(t), vmax=np.max(t))

        for i, v_i in enumerate(t):
            spec = data[i]
            color = cmap(norm(v_i))
            ax.plot(wl, spec, color=color)

        ax.set_xlabel(spectrum_series.x_name)
        ax.set_ylabel()

        if make_colorbar:
            cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
            cb.set_label(spectrum_series.potential.name)

        return ax
