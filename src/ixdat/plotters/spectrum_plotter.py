"""Plotters for spectra and spectrumseries."""

import numpy as np
import matplotlib as mpl
from ixdat.plotters.plotting_tools import add_colorbar
from matplotlib import pyplot as plt
from .base_mpl_plotter import MPLPlotter


class SpectrumPlotter(MPLPlotter):
    """A plotter for a spectrum"""

    def __init__(self, spectrum=None):
        super().__init__()
        self.spectrum = spectrum

    def plot(self, *, spectrum=None, ax=None, **kwargs):
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

    def plot_average(self, *, spectrum_series=None, ax=None, **kwargs):
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
        *,
        spectrum_series=None,
        field=None,
        tspan=None,
        xspan=None,
        ax=None,
        cmap_name="inferno",
        make_colorbar=False,
        t=None,
        t_name=None,
        max_threshold=None,
        min_threshold=None,
        scanning_mask=None,
        vmin=None,
        vmax=None,
        continuous=None,
    ):
        """
        Plot a spectrum series with `t` on the horizontal axis, `x` on the vertical axis,
        and color representing `y`.

        Args:
            spectrum_series (SpectrumSeries): The spectrum series to be plotted, if
                different from self.spectrum_series.
            field (Field): The field to be plotted, if different from
                spectrum_series.field
            tspan (iterable): The span of the time data to plot
            xspan (iterable): The span of the spectral data to plot
            ax (mpl.Axis): The axes to plot on. A new one is made by default

            cmap_name (str): The name of the colormap to use. Defaults to "inferno", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html#sequential
            make_colorbar (bool): Whether to make a colorbar.
                FIXME: colorbar at present mis-alignes axes
            t (numpy array): Time data to use if not the data in spectrum_series
            t_name (str): Name of time variable if not the one in spectrum_series
            max_threshold (float): Maximum value to display.
                Values above are set to zero.
            min_threshold (float): Minimum value to display.
                Values below are set to 0.
            scanning_mask (list): List of booleans to exclude from scanning variable
                before plotting data by setting y values to 0 (zero).
            vmin (float): minimum value to represent in colours.
            vmax (float): maximum value to represent in colours.
            continuous (bool): Optional. Whether to make a continuous heat plot (True) or
                a discrete heat plot for each spectrum (False). In the discrete case,
                each heat plot is a rectangle with the spectrum's duration as its width,
                if available. If the duration is not available, each spectrum heat plot
                extends to the start of the next one.
                Defaults to the `spectrum_series.continuous`.
        """
        spectrum_series = spectrum_series or self.spectrum_series
        field = field or spectrum_series.field
        if continuous is None:
            continuous = spectrum_series.continuous

        xseries = field.axes_series[1]
        x = xseries.data
        t = t if t is not None else field.axes_series[0].t
        t_name = t_name or field.axes_series[0].name

        data = field.data

        if max_threshold:
            # data = np.minimum(max_threshold, data)
            data[data > max_threshold] = 0
        if min_threshold:
            data[data < min_threshold] = 0

        if np.any(scanning_mask):
            data[:, scanning_mask] = 0

        if tspan:
            t_mask = np.logical_and(tspan[0] < t, t < tspan[-1])
            t = t[t_mask]
            data = data[t_mask, :]
            if (t[0] < t[-1]) != (tspan[0] < tspan[-1]):  # this is an XOR.
                # Then we need to plot the data against U in the reverse direction:
                t = np.flip(t, axis=0)
                data = np.flip(data, axis=0)

        if xspan:
            x_mask = np.logical_and(xspan[0] < x, x < xspan[-1])
            x = x[x_mask]
            data = data[:, x_mask]

        if not ax:
            ax = self.new_ax()

        if continuous:
            ax.imshow(
                np.flip(data.swapaxes(0, 1), axis=0),
                cmap=cmap_name,
                aspect="auto",
                extent=(t[0], t[-1], x[0], x[-1]),
                vmin=vmin,
                vmax=vmax,
            )
        else:
            for i, t_i in enumerate(spectrum_series.t):
                if tspan and (t_i < min(tspan) or t_i > max(tspan)):
                    continue
                try:
                    duration = spectrum_series.durations[i]
                    # ^ raises TypeError if durations is None.
                    t_f = t_i + duration  # raises TypeError if durations[i] is None.
                except TypeError:
                    if i < len(t) - 1:
                        t_f = t[i + 1]
                    else:
                        # If its duration is unknown, we don't plot the last spectrum.
                        break
                y = data[i]
                yy = np.stack([y, y])
                ax.imshow(
                    np.flip(yy.swapaxes(0, 1), axis=0),
                    cmap=cmap_name,
                    aspect="auto",
                    extent=(t_i, t_f, x[0], x[-1]),
                    vmin=vmin,
                    vmax=vmax,
                )

        ax.set_xlabel(t_name)
        ax.set_ylabel(xseries.name)

        if make_colorbar:
            add_colorbar(
                ax,
                cmap_name,
                vmin=(vmin if vmin else np.min(data)),
                vmax=(vmax if vmax else np.max(data)),
            )
        return ax

    def plot_waterfall(
        self,
        *,
        spectrum_series=None,
        field=None,
        ax=None,
        cmap_name="jet",
        make_colorbar=True,
        t=None,
        t_name=None,
    ):
        """Plot a SpectrumSeries as spectra colored by the time at which they are taken

        Args:
            spectrum_series (SpectrumSeries): The spectrum series to be plotted, if
                different from self.spectrum_series.
            field (Field): The field to be plotted, if different from
                spectrum_series.field
            ax (matplotlib Axis): The axes to plot on. A new one is made by default.

            cmap_name (str): The name of the colormap to use. Defaults to "inferno", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html#sequential
            make_colorbar (bool): Whether to make a colorbar.
            t (numpy array): Time data to use if not the data in spectrum_series
            t_name (str): Name of time variable if not the one in spectrum_series
        """
        spectrum_series = spectrum_series or self.spectrum_series
        field = field or spectrum_series.field

        data = field.data
        x = field.axes_series[1].data
        t = t if t is not None else field.axes_series[0].t
        t_name = t_name or field.axes_series[0].name

        cmap = mpl.cm.get_cmap(cmap_name)
        norm = mpl.colors.Normalize(vmin=np.min(t), vmax=np.max(t))

        if not ax:
            ax = self.new_ax()

        for i, t_i in enumerate(t):
            spec = data[i]
            color = cmap(norm(t_i))
            ax.plot(x, spec, color=color)

        ax.set_xlabel(field.axes_series[1].name)
        ax.set_ylabel(field.name)

        if make_colorbar:
            cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
            cb.set_label(t_name)

        return ax


class SpectroMeasurementPlotter(MPLPlotter):
    """Plotter for measurements with spectrum_series

    This makes use of the methods in `SpectrumSeriesPlotter`, but allows a second
    scanned variable (such as `potential` in the case of spectroelectrochemistry) to be
    used from the measurement instead of time data.
    """

    def __init__(self, measurement=None):
        super().__init__()
        self.measurement = measurement
        self.spectrum_series_plotter = SpectrumSeriesPlotter()

    def heat_plot_vs(
        self,
        *,
        measurement=None,
        field=None,
        vs=None,
        vspan=None,
        xspan=None,
        ax=None,
        cmap_name="inferno",
        make_colorbar=False,
    ):
        """Plot a SpectroMeasurement in two panels with time as x-asis.

        The top panel is a heat plot with wavelength on y-axis and color representing
        spectrum. At most one of V_ref and t_ref should be given, and if neither are
        given the measurement's default reference_spectrum is used to calculate the
        optical density.

        Args:
            measurement (SpectrumSeries): The spectrum series to be plotted, if
                different from self.spectrum_series.
            field (Field): The field to be plotted, if different from
                spectrum_series.field
            vs (str): The name of the value series or time series to plot against.
            vspan (iterable): The span of the value series or time series to include
            xspan (iterable): The span of the spectral data to plot
            ax (mpl.Axis): The axes to plot on. A new one is made by default

            cmap_name (str): The name of the colormap to use. Defaults to "inferno", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html#sequential
            make_colorbar (bool): Whether to make a colorbar.
                FIXME: colorbar at present mis-alignes axes
        """
        measurement = measurement or self.measurement

        tseries = field.axes_series[0]
        v_name = vs
        if vs in ("t", tseries.tseries.name):
            v = tseries.t
            if hasattr(measurement, "t_str") and measurement.t_str:
                v_name = measurement.t_str
        else:
            v = measurement.grab_for_t(vs, t=tseries.t)

        return self.spectrum_series_plotter.heat_plot(
            spectrum_series=measurement.spectrum_series,
            field=field,
            tspan=vspan,
            xspan=xspan,
            ax=ax,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            t=v,
            t_name=v_name,
        )

    def plot_waterfall_vs(
        self,
        *,
        measurement=None,
        field=None,
        vs=None,
        ax=None,
        cmap_name="jet",
        make_colorbar=True,
    ):
        """Plot a SpectrumSeries as spectra colored by the value at which they are taken

        Args:
            measurement (SpectroMeasurement): The measurement to be plotted if different
                from self.measurement.
            field (Field): The field to be plotted, if different from
                spectrum_series.field
            vs (str): The name of the value to use for the color scale. Defaults to time
            ax (matplotlib Axis): The axes to plot on. A new one is made by default.
            cmap_name (str): The name of the colormap to use. Defaults to "inferno", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html#sequential
            make_colorbar (bool): Whether to make a colorbar.
        """
        measurement = measurement or self.measurement

        tseries = field.axes_series[0]
        v_name = vs
        if vs in ("t", tseries.tseries.name):
            v = tseries.t
            if hasattr(measurement, "t_str") and measurement.t_str:
                v_name = measurement.t_str
        else:
            v = measurement.grab_for_t(vs, t=tseries.t)

        return self.spectrum_series_plotter.plot_waterfall(
            spectrum_series=measurement.spectrum_series,
            field=field,
            ax=ax,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            t=v,
            t_name=v_name,
        )
