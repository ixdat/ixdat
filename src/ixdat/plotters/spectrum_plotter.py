"""Plotters for spectra and spectrumseries."""

import numpy as np
import matplotlib as mpl
from ixdat.plotters.plotting_tools import add_colorbar
from matplotlib import pyplot as plt
from .base_mpl_plotter import MPLPlotter


class SpectrumPlotter(MPLPlotter):
    """A plotter for spectrum_list"""

    def __init__(self, spectrum=None):
        super().__init__()
        self.spectrum = spectrum

    def plot(self, spectrum=None, ax=None, **kwargs):
        """Plot a spectrum as y (signal) vs x (scanning variable)

        Args:
            spectrum (Spectrum): The spectrum to plot if different from self.spectrum
            ax (mpl.Axis): The axis to plot on. A new one is made by default.
            kwargs: additional key-word arguments are given to ax.plot()
        """
        spectrum = spectrum or self.spectrum
        if not ax:
            ax = self.new_ax()
        ax.plot(spectrum.x, spectrum.y, **kwargs)
        ax.set_xlabel(spectrum.x_name)
        ax.set_ylabel(spectrum.y_name)
        return ax


class SpectrumSeriesPlotter(MPLPlotter):
    """A plotter for spectrum series, f.ex. spectra taken continuously over time"""

    def __init__(self, spectrum_series=None):
        super().__init__()
        self.spectrum_series = spectrum_series

    @property
    def plot(self):
        """The default plot of a SpectrumSeries is heat_plot"""
        return self.heat_plot

    def plot_average(self, spectrum_series=None, ax=None, **kwargs):
        """Take an average of the spectra and plot that."""
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
        field=None,
        tspan=None,
        xspan=None,
        ax=None,
        cmap_name="inferno",
        make_colorbar=False,
    ):
        """Plot with time as x, the scanning variable as y, and color as signal

        See SpectrumSeriesPlotter.heat_plot_vs(). This function calls it with vs="t".
        """
        return self.heat_plot_vs(
            spectrum_series=spectrum_series,
            field=field,
            vspan=tspan,
            xspan=xspan,
            ax=ax,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            vs="t",
        )

    def heat_plot_vs(
        self,
        spectrum_series=None,
        field=None,
        vspan=None,
        xspan=None,
        ax=None,
        cmap_name="inferno",
        make_colorbar=False,
        vs=None,
    ):
        """Plot an SECMeasurement in two panels with time as x-asis.

        The top panel is a heat plot with wavelength on y-axis and color representing
        spectrum. At most one of V_ref and t_ref should be given, and if neither are
        given the measurement's default reference_spectrum is used to calculate the
        optical density.

        Args:
            spectrum_series (SpectrumSeries): The spectrum series to be plotted, if
                different from self.spectrum_series.
                FIXME: spectrum_series needs to actually be a Measurement to have other
                FIXME:  series to plot against if vs isn't in field.series_axes
            field (Field): The field to be plotted, if different from
                spectrum_series.field
            xspan (iterable): The span of the spectral data to plot
            ax (mpl.Axis): The axes to plot on. A new one is made by default
            cmap_name (str): The name of the colormap to use. Defaults to "inferno",
                which ranges from black through red and orange to yellow-white. "jet"
                is also good.
            make_colorbar (bool): Whether to make a colorbar.
                FIXME: colorbar at present mis-alignes axes
            vs (str): The ValueSeries or TimeSeries to plot against.
        """
        spectrum_series = spectrum_series or self.spectrum_series
        field = field or spectrum_series.field

        xseries = field.axes_series[1]
        x = xseries.data
        tseries = field.axes_series[0]

        v_name = vs
        if vs in ("t", tseries.tseries.name):
            v = tseries.t
            if hasattr(spectrum_series, "t_str") and spectrum_series.t_str:
                v_name = spectrum_series.t_str
        else:
            v = spectrum_series.grab_for_t(vs, t=tseries.t)

        data = field.data
        # ^ FIXME: The heat plot will be distorted if spectra are not taken at even
        #     spacing on the "vs" variable. They will be especially meaningless if
        #     the U variable itself is not always increasing or decreasing.

        if vspan:
            v_mask = np.logical_and(vspan[0] < v, v < vspan[-1])
            v = v[v_mask]
            data = data[v_mask, :]
            if (v[0] < v[-1]) != (vspan[0] < vspan[-1]):  # this is an XOR.
                # Then we need to plot the data against U in the reverse direction:
                v = np.flip(v, axis=0)
                data = np.flip(data, axis=0)
        if xspan:
            wl_mask = np.logical_and(xspan[0] < x, x < xspan[-1])
            x = x[wl_mask]
            data = data[:, wl_mask]

        ax.imshow(
            data.swapaxes(0, 1),
            cmap=cmap_name,
            aspect="auto",
            extent=(v[0], v[-1], x[0], x[-1]),
        )
        ax.set_xlabel(v_name)
        ax.set_ylabel(xseries.name)
        if make_colorbar:
            add_colorbar(ax, cmap_name, vmin=np.min(data), vmax=np.max(data))
        return ax

    def plot_waterfall(
        self,
        spectrum_series=None,
        field=None,
        cmap_name="jet",
        make_colorbar=True,
        vs=None,
        ax=None,
    ):
        """Plot a SpectrumSeries as spectra colored by the value at which they are taken

        Args:
            spectrum_series (SpectrumSeries): The spectrum series to be plotted, if
                different from self.spectrum_series.
                FIXME: spectrum_series needs to actually be a Measurement to have other
                FIXME: ...series to plot against if vs isn't in field.series_axes
            field (Field): The field to be plotted, if different from
                spectrum_series.field
            ax (matplotlib Axis): The axes to plot on. A new one is made by default.
            cmap_name (str): The name of the colormap to use. Defaults to "inferno",
                which ranges from black through red and orange to yellow-white. "jet"
                is also good.
            make_colorbar (bool): Whether to make a colorbar.
            vs (str): The name of the value to use for the color scale. Defaults to time
        """
        spectrum_series = spectrum_series or self.spectrum_series
        field = field or spectrum_series.field

        data = field.data
        t = field.axes_series[0].t
        x = field.axes_series[1].data

        if vs:
            v = spectrum_series.grab_for_t(vs, t=t)
        else:
            v = t

        cmap = mpl.cm.get_cmap(cmap_name)
        norm = mpl.colors.Normalize(vmin=np.min(v), vmax=np.max(v))

        if not ax:
            ax = self.new_ax()

        for i, v_i in enumerate(v):
            spec = data[i]
            color = cmap(norm(v_i))
            ax.plot(x, spec, color=color)

        ax.set_xlabel(field.axes_series[1].name)
        ax.set_ylabel(field.name)

        if make_colorbar:
            cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
            cb.set_label(vs)

        return ax
