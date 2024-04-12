"""Plotters for spectra and spectrumseries."""

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from . import MPLPlotter, add_colorbar
from .plotting_tools import get_indeces_and_times


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

        This is commonly used for e.g. in-situ UV-Vis spectrometry

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

    def plot_stacked_spectra(
        self,
        spectrum_series=None,
        dt=None,
        t_list=None,
        dn=None,
        index_list=None,
        average=False,
        xspan=None,
        xspan_bg=None,
        scale_mode="auto",
        scale_factor=1,
        y_values="time",
        ax=None,
        color="k",
        **kwargs,
    ):
        """Plot a selection of spectra, stacked

        This is commonly used for e.g. FTIR.
        Specify which spectra to plot by one of four ways: dt, t_list,
        dn, or n_list. See descriptions below.

        Args:
            spectrum_series (SpectrumSeries): What to plot from, if
                different from self.spectrum_series
            dt (float): time interval between spectra to plot, [s]. The
                first spectrum and those taken at times closest to each
                integer multiple of dt after are plotted.
            t_list (list of float): List of times for which to plot the
                spectrum, [s]. The closest spectrum to each time in the
                list is plotted.
            dn (int): number of spectra between plotted spectra
            index_list (list of int): List of indeces of spectra to plot
            average (bool or int): Whether and how to average spectra for
                plotting. False means no averaging. True means average all
                the spectra in the interval between spectra. An integer
                `n` means average the `n/2` spectra before and `n/2` spectra
                after the spectra at the given time or index
            xspan (list of float): Range of x-axis variable to include.
            xspan_bg (list of float): Range of x-axis variable for which to
                consider the signal at background. For each spectrum, the
                average y value in this range is subtracted.
            scale_mode (str): The way to initially scale the spectra.
                Options:
                - "auto": scale uniformly such that all spectra fit in
                    their given interval. The raw y-values are scaled by
                    min(interval) / max(y range)
                - [no other scale_mode options yet]
            scale_factor: A factor to apply on top of the initial scaling
            y_values (str): What to plot on the y-axis. Options: "time", "n"
            ax (Axis): axis to plot on, if not a new axis
            color (str): color of traces. Defaults to black.
            **kwargs: Additional key-word args are passed on to ax.plot()
        """

        spectrum_series = spectrum_series or self.spectrum_series

        # whichever of the four were specified, we need t_list and index_list:
        t_vec = spectrum_series.t
        index_list, t_list = get_indeces_and_times(
            t_vec, dt=dt, t_list=t_list, dn=dn, index_list=index_list
        )

        # Now we get the y-values of the spectra to plot
        y_vec_list = []
        for i, index in enumerate(index_list):
            if average:
                # The challenge with averaging is figuring out the range of ideces over
                # which to average.
                if type(average) is int:
                    range_start = max(0, index - average)
                    range_end = min(index + average, len(t_vec))
                else:
                    # If average is set to "True", we need to figure out the median
                    # points between the spectra, as follows.
                    if i == 0:
                        range_start = index
                    else:
                        range_start = int((index + index_list[i - 1]) / 2)
                    if i + 1 == len(index_list):
                        range_end = index
                    else:
                        range_end = int((index + index_list[i + 1]) / 2)
                # Once we have the range of indeces, the SpectrumSeries makes it easy:
                y_vec = spectrum_series[range_start:range_end].y_average
            else:
                # The simple case where we're not averaging.
                y_vec = spectrum_series[index].y
            y_vec_list.append(y_vec)

        # Establish the axis
        if not ax:
            ax = self.new_ax()
            ax.set_xlabel(spectrum_series.xseries.name)

        # Whatever y's we are plotting, the x is the same:
        x = spectrum_series.x

        # If an xspan_bg is given, zero the corresponding part of each y vector:
        if xspan_bg:
            mask_bg = np.logical_and(xspan_bg[0] < x, x < xspan_bg[-1])
            for i, y_vec in enumerate(y_vec_list):
                y_vec_list[i] = y_vec_list[i] - np.mean(y_vec_list[i][mask_bg])
        # If an xspan (for plotting) is given, cut all the vectors accordingly:
        if xspan:
            mask = np.logical_and(xspan[0] < x, x < xspan[-1])
            x = x[mask]
            for i, y_vec in enumerate(y_vec_list):
                y_vec_list[i] = y_vec[mask]

        # a list of y ranges is useful in all cases for figuring out how to scale the
        # y-vectors onto the y axis, which a bit confusingly represents something else.
        y_ranges = np.max(y_vec_list, axis=1) - np.min(y_vec_list, axis=1)

        # No matter what, we loop through the y-vectors to plot, but the code is a bit
        # different depending on what the y-axis represents, so consider those two cases
        # ("time" and "n", which means index) seperately.
        if y_values == "time":
            ax.set_ylabel(spectrum_series.t_name)
            if scale_mode == "auto":
                # time intervals:
                dts = [t_list[i + 1] - t_list[i] for i in range(len(t_list) - 1)]
                # The scale is determined by the minimum time interval divided by
                # the maximum y range, so that by default (scale_factor=1) all spectra
                # fit completely in thier own time interval
                t_per_y = min(dts) / max(y_ranges) * scale_factor
                scaled_y_vec_list = [y_vec * t_per_y for y_vec in y_vec_list]
            else:
                raise ValueError(f"scale_mode={scale_mode} not implemented.")
            for t, scaled_y_vec in zip(t_list, scaled_y_vec_list):
                ax.plot(x, t + scaled_y_vec, label=t, color=color, **kwargs)
        elif y_values == "n":
            ax.set_ylabel("spectrum number")
            if scale_mode == "auto":
                # index intervals:
                dns = [
                    index_list[i + 1] - index_list[i] for i in range(len(index_list) - 1)
                ]
                # The scale is determined by the minimum index interval divided by
                # the maximum y range, so that by default (scale_factor=1) all spectra
                # fit completely in thier own index interval
                t_per_n = min(dns) / max(y_ranges) * scale_factor
                scaled_y_vec_list = [y_vec * t_per_n for y_vec in y_vec_list]
            else:
                raise ValueError(f"scale_mode='{scale_mode}' not implemented.")
            for n, scaled_y_vec in zip(index_list, scaled_y_vec_list):
                ax.plot(x, n + scaled_y_vec, label=n, color=color, **kwargs)
        else:
            raise ValueError(f"y_values='{y_values}' not implemented.")

        # Make the spectra go all the way to the plot edges, with ticks on both sides:
        ax.set_xlim(min(x), max(x))
        ax.tick_params(axis="y", right=True)

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

    def plot_stacked_spectra_vs(
        self,
        *,
        measurement=None,
        vs=None,
        dt=None,
        t_list=None,
        dn=None,
        index_list=None,
        average=False,
        xspan=None,
        xspan_bg=None,
        scale_mode="auto",
        scale_factor=1,
        y_values="time",
        ax=None,
        color="k",
        **kwargs,
    ):
        """Plot stacked spectra with a time-dependent value on the y-axis

        This is commonly used for e.g. FTIR.
        Specify which spectra to plot by one of four ways: dt, t_list,
        dn, or n_list. See descriptions below.

        Args:
            measurement (SpectroMeasurement): The spectromeasurement to
                plot data form, if different from self.measurement
            vs (str): The name of the value series to stack spectra
                according to.
            dt (float): time interval between spectra to plot, [s]. The
                first spectrum and those taken at times closest to each
                integer multiple of dt after are plotted.
            t_list (list of float): List of times for which to plot the
                spectrum, [s]. The closest spectrum to each time in the
                list is plotted.
            dn (int): number of spectra between plotted spectra
            index_list (list of int): List of indeces of spectra to plot
            average (bool or int): Whether and how to average spectra for
                plotting. False means no averaging. True means average all
                the spectra in the interval between spectra. An integer
                `n` means average the `n/2` spectra before and `n/2` spectra
                after the spectra at the given time or index
            xspan (list of float): Range of x-axis variable to include.
            xspan_bg (list of float): Range of x-axis variable for which to
                consider the signal at background. For each spectrum, the
                average y value in this range is subtracted.
            scale_mode (str): The way to initially scale the spectra.
                Options:
                - "auto": scale uniformly such that all spectra fit in
                    their given interval. The raw y-values are scaled by
                    min(interval) / max(y range)
                - [no other scale_mode options yet]
            scale_factor: A factor to apply on top of the initial scaling
            y_values (str): What to plot on the y-axis. Options: "time", "n"
            ax (Axis): axis to plot on, if not a new axis
            color (str): color of traces. Defaults to black.
            **kwargs: Additional key-word args are passed on to ax.plot()
        """
        measurement = measurement or self.measurement

        t_vec = measurement.spectrum_series.t

        ax = self.spectrum_series_plotter.plot_stacked_spectra(
            spectrum_series=measurement.spectrum_series,
            dt=dt,
            t_list=t_list,
            dn=dn,
            index_list=index_list,
            average=average,
            xspan=xspan,
            xspan_bg=xspan_bg,
            scale_mode=scale_mode,
            scale_factor=scale_factor,
            y_values=y_values,
            ax=ax,
            color=color,
        )
        index_list, t_list = get_indeces_and_times(t_vec, dt, t_list, dn, index_list)

        v_list = measurement.grab_for_t(vs, np.array(t_list))

        ax.set_yticks(t_list)
        ax.set_yticklabels([np.round(v, 2) for v in v_list])

        ax.set_ylabel(vs)
