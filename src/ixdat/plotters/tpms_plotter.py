import warnings
from .base_mpl_plotter import MPLPlotter
from .ms_plotter import MSPlotter
from .plotting_tools import color_axis
import numpy as np


class TPMSPlotter(MPLPlotter):
    """A matplotlib plotter for TP-MS measurements."""

    def __init__(self, measurement=None):
        """Initiate the TPMSPlotter with its default Measurement to plot"""
        super().__init__()
        self.measurement = measurement
        self.ms_plotter = MSPlotter(measurement=measurement)

    def plot_measurement(
        self,
        *,
        measurement=None,
        axes=None,
        mass_list=None,
        mass_lists=None,
        mol_list=None,
        mol_lists=None,
        tspan=None,
        tspan_bg=None,
        remove_background=None,
        unit=None,
        x_unit=None,
        meta_units=None,
        meta_list=None,
        meta_lists=None,
        T_name=None,
        P_name=None,
        T_color=None,
        P_color=None,
        logplot=None,
        logdata=None,
        legend=True,
        emphasis="top",
        **kwargs,
    ):
        """Make a two panel plot with mass spec data on top panel and meta data on bottom
        panel. Default Temperature and Pressure. TP-MS plot return the axis handles.

        Allocates some tasks to MSPlotter.plot_measurement()

        Args:
            measurement (TPMSMeasurement): Defaults to the measurement to which the
                plotter is bound (self.measurement)
            axes (list of three matplotlib axes): axes[0] plots the MID data,
                axes[1] the variable given by T_name (temperature), and axes[3] the
                variable given by P_name (reactor pressure). By default three axes are
                made with axes[0] a top panel with 3/5 the area, and axes[1] and axes[3]
                are the left and right y-axes of the lower panel with 2/5 the area.
            mass_list (list of str): The names of the m/z values, eg. ["M2", ...] to
                plot. Defaults to all of them (measurement.mass_list)
            mass_lists (list of list of str): Alternately, two lists can be given for
                masses in which case one list is plotted on the left y-axis and the other
                on the right y-axis of the top panel.
            mol_list (list of str): The names of the molecules, eg. ["H2", ...] to
                plot. Defaults to all of them if quantified (measurement.mass_list)
            mol_lists (list of list of str): Alternately, two lists can be given for
                molecules in which case one list is plotted on the left y-axis and the
                other on the right y-axis of the top panel.
            tspan (iter of float): The time interval to plot, wrt measurement.tstamp
            tspan_bg (timespan): A timespan for which to assume the signal is at its
                background. The average signals during this timespan are subtracted.
                If `mass_lists` are given rather than a single `mass_list`, `tspan_bg`
                must also be two timespans - one for each axis. Default is `None` for no
                background subtraction.
            remove_background (bool): Whether otherwise to subtract pre-determined
                background signals if available. Defaults to (not logplot)
            unit (str): the unit for the MS data. Defaults to "A" for Ampere
            x_unit (str): unit for x axis defaults to 's' for seconds
            meta_units (dict): dictionary with new units for y axis on bottom panel.
                Default to ValueSeries.unit.name for last plotted meta_name
            meta_list (list of str): The names of the data_series plotted on axes[1].
                Defulats None.
            meta_lists (list of list of str): The names of the data series in two list to
                be plotted on axes[1] and axes[3] respectively. Defaults None.
            T_name (str): The name of the value to plot on the lower left y-axis.
                Defaults to the name of the series `measurement.temperature`
            P_name (str): The name of the value to plot on the lower right y-axis.
                Defaults to the name of the series `measurement.pressure`
            T_color (str): Overwrite standard color to plot the variable given by T_name
            P_color (str): Overwrite standard color to plot the variable given by P_name
            logplot (bool): Whether to plot the MS data on a log scale (default True
                unless mass_lists are given)
            logdata (bool): Whether to plot the MS data on a log scale (default True
                unless mass_lists are given)
            legend (bool): Whether to use a legend for the MS data (default True)
            emphasis (str or None): "top" for bigger top panel, "bottom" for bigger
                bottom panel, None for equal-sized panels, "one figure" to plot all in
                one figure
            kwargs (dict): Additional kwargs go to all calls of matplotlib's plot()

        Returns:
            list of Axes: (top_left, bottom_left, top_right, bottom_right) where:
                axes[0] is top_left is MS data;
                axes[1] is bottom_left is temperature;
                axes[2] is top_right is additional MS data if left and right mass_lists
                    or mol_lists were plotted (otherwise axes[2] is None); and
                axes[3] is bottom_right is pressure.
        """

        measurement = measurement or self.measurement

        axes = self.new_two_panel_axes(
            n_bottom=2,
            n_top=(2 if (mass_lists or mol_lists) else 1),
            emphasis=emphasis,
        )

        if (
            mass_list
            or mass_lists
            or mol_list
            or mol_lists
            or hasattr(measurement, "mass_list")
        ):
            # then we have MS data!
            self.ms_plotter.plot_measurement(
                measurement=measurement,
                axes=[axes[0], axes[2]] if (mass_lists or mol_lists) else [axes[0]],
                tspan=tspan,
                tspan_bg=tspan_bg,
                remove_background=remove_background,
                mass_list=mass_list,
                mass_lists=mass_lists,
                mol_list=mol_list,
                mol_lists=mol_lists,
                unit=unit,
                x_unit=x_unit,
                logplot=logplot,
                logdata=logdata,
                legend=legend,
                **kwargs,
            )

        T_name = T_name or measurement.T_name
        P_name = P_name or measurement.P_name

        if not meta_lists:
            meta_lists = [[T_name], [P_name]] if not meta_list else [meta_list]
        for i, meta_list in enumerate(meta_lists):
            # plot on dataseries on correct axis
            ax = [axes[1], axes[3]][i]
            for n, meta_name in enumerate(meta_list):
                if meta_name == T_name and T_color:
                    color = T_color
                elif meta_name == P_name and P_color:
                    color = P_color
                else:
                    color = STANDARD_COLORS.get(meta_name, "k")
                t, v = measurement.grab_signal(
                    meta_name,
                    tspan=tspan,
                    include_endpoints=False,
                )

                y_label, y_unit, y_unit_factor = _get_y_unit_and_label(
                                                        measurement[meta_name],
                                                        meta_units=meta_units
                                                        )

                # expect always to plot against time on x axis
                x_unit_factor, x_unit = _get_unit_factor_and_name(
                                                        new_unit_name=x_unit,
                                                        from_unit_name="s"
                                                        )

                ax.plot(
                    t * x_unit_factor,
                    v * y_unit_factor,
                    color=color,
                    label=meta_name,
                    **kwargs,
                )

            if logplot and y_unit == "mbar":
                ax.set_yscale("log")
            if legend:
                ax.legend()

            if not n:
                color_axis(ax, color=color, lr=['left', 'right'][i])

            ax.set_ylabel(f"{y_label} / [{y_unit}]")
            ax.set_xlabel(f"time / [{x_unit}]")

        if not i: # remove last axis if nothing is plotted on it
            axes[3].remove()

        return axes


    def plot_arrhenius(
        self,
        *,
        T_name,
        measurement=None,
        ax=None,
        axes=None,
        mass_list=None,
        mol_list=None,
        tspan=None,
        tspan_bg=None,
        remove_background=None,
        unit=None,
        x_unit=None,
        arrh_color="r",
        logplot=None,
        logdata=None,
        legend=True,
        **kwargs,
    ):

        measurement = measurement or self.measurement
        T_name = T_name or measurement.T_name
        if (
            mass_list
            or mol_list
            or hasattr(measurement, "mass_list")
        ):
            # then we have MS data!

            self.ms_plotter.plot_vs(
                measurement=measurement,
                ax=ax,
                mass_list=mass_list,
                mol_list=mol_list,
                tspan=tspan,
                tspan_bg=tspan_bg,
                remove_background=remove_background,
                unit=unit,
                x_unit=x_unit,
                x_name=T_name,
                logplot=logplot,
                logdata=logdata,
                legend=legend,
                **kwargs,
            )

        return axes

    def plot_measurement_in_same_figure(
        self,
        *,
        measurement=None,
        ax=None,
        axes=None,
        mass_list=None,
        mol_list=None,
        meta_list=None,
        tspan=None,
        tspan_bg=None,
        remove_background=None,
        unit=None,
        x_unit=None,
        meta_units=None,
        T_name=None,
        T_color=None,
        P_name=None,
        P_color=None,
        logplot=None,
        logdata=None,
        legend=True,
        hightlighted=None,
        **kwargs,
    ):
        """Make a one panel plot with mass spec data on left y-axis and meta data right
        y-axis. Defaultis to all masses and Temperature and Pressure.
        TP-MS plot return the axis handles.

        Allocates some tasks to MSPlotter.plot_measurement()

        Args:
            measurement (TPMSMeasurement): Defaults to the measurement to which the
                plotter is bound (self.measurement)
            ax (matplotlib ax): If only ax is given, this is used to plot MID data on and
                create a second twonx() ax to plot meta_list on.
            axes (list of two matplotlib axes): axes[0] plots the MID data,
                axes[1] the variable given by T_name (temperature) and the variable given
                by P_name (reactor pressure). By default three axes are
                made with axes[0] a top panel with 3/5 the area, and axes[1] and axes[3]
                are the left and right y-axes of the lower panel with 2/5 the area.
            mass_list (list of str): The names of the m/z values, eg. ["M2", ...] to
                plot. Defaults to all of them (measurement.mass_list)
            mol_list (list of str): The names of the molecules, eg. ["H2", ...] to
                plot. Defaults to all of them if quantified (measurement.mass_list)
            tspan (iter of float): The time interval to plot, wrt measurement.tstamp
            tspan_bg (timespan): A timespan for which to assume the signal is at its
                background. The average signals during this timespan are subtracted.
                If `mass_lists` are given rather than a single `mass_list`, `tspan_bg`
                must also be two timespans - one for each axis. Default is `None` for no
                background subtraction.
            remove_background (bool): Whether otherwise to subtract pre-determined
                background signals if available. Defaults to (not logplot)
            unit (str): the unit for the MS data. Defaults to "A" for Ampere
            x_unit (str): unit for x axis defaults to 's' for seconds
            meta_units (dict): Dictionary with new units meta series to be plotted.
                Default to ValueSeries.unit.name for last plotted ValueSeries
            meta_list (list of str): The names of the data_series plotted on axes[1].
                Defulats None.
            T_name (str): The name of the value to plot on the right y-axis.
                Defaults to the name of the series `measurement.temperature`
            P_name (str): The name of the value to plot on the right y-axis.
                Defaults to the name of the series `measurement.pressure`
            T_color (str): Overwrite standard color to plot the variable given by T_name
            P_color (str): Overwrite standard color to plot the variable given by P_name
            logplot (bool): Whether to plot the MS data on a log scale (default True
                unless mass_lists are given)
            logdata (bool): Whether to plot the MS data on a log scale (default True
                unless mass_lists are given)
            legend (bool): Whether to use a legend for the MS data (default True)
            emphasis (str or None): "top" for bigger top panel, "bottom" for bigger
                bottom panel, None for equal-sized panels, "one figure" to plot all in
                one figure
            kwargs (dict): Additional kwargs go to all calls of matplotlib's plot()

        Returns:
                axes[0] is MS data;
                axes[1] is meta data (default temperature and pressure);
        """

        measurement = measurement or self.measurement

        if not ax:
            ax = (
                axes[0]
                if axes
                else self.new_ax(ylabel=f"signal / [{unit}]", xlabel="time / [s]")
            )
        axes = axes or [ax, ax.twinx()]  # prepare an axis unless we were given two.

        # then we have MS data!
        self.ms_plotter.plot_measurement(
                measurement=measurement,
                ax=axes[0],
                tspan=tspan,
                tspan_bg=tspan_bg,
                remove_background=remove_background,
                mass_list=mass_list,
                mol_list=mol_list,
                unit=unit,
                x_unit=x_unit,
                logplot=logplot,
                logdata=logdata,
                legend=legend,
                **kwargs,
        )

        T_name = T_name or measurement.T_name
        P_name = P_name or measurement.P_name

        # figure out if one ot two axes is to be plotted in bottom panel
        if not meta_list:
            meta_list = [T_name, P_name]

        for n, meta_name in enumerate(meta_list):
            if meta_name == T_name and T_color:
                color = T_color
            elif meta_name == P_name and P_color:
                color = P_color
            else:
                color = STANDARD_COLORS.get(meta_name, "k")

            t, v = measurement.grab_signal(
                meta_name,
                tspan=tspan,
                include_endpoints=False,
            )

            y_label, y_unit, y_unit_factor = _get_y_unit_and_label(
                                                    measurement[meta_name],
                                                    meta_units=meta_units
                                                    )

            # expect always to plot against time
            x_unit_factor, x_unit = _get_unit_factor_and_name(
                                                    new_unit_name=x_unit,
                                                    from_unit_name="s"
                                                    )

            axes[1].plot(
                t * x_unit_factor,
                v * y_unit_factor,
                color=color,
                label=meta_name,
                alpha=1,
                **kwargs,
            )
        if logplot and y_unit == "mbar":
            axes[1].set_yscale("log")
        if legend:
            axes[1].legend()


        if n > 1:
            color = 'k'
            y_label = "signal"
            y_unit = "mixed"

        color_axis(axes[1], color=color, lr='right')
        axes[1].set_ylabel(f"{y_label} / [{y_unit}]")
        axes[1].set_xlabel(f"time / [{x_unit}]")
        #if hightlighted:
        #    for ax in axes:
        #        for line in ax.get_lines:
        #            line.set_alpha(0.3)
        #            if line.get_label() in highlighted:
        #                line.set_alpha(1)

        return axes

class SpectroTPMSPlotter(MPLPlotter):
    def __init__(self, measurement=None):
        """Initiate the Spectro-TPMSPlotter with its default Measurement to plot"""
        super().__init__()
        self.measurement = measurement
        self.tpms_plotter = TPMSPlotter(measurement=measurement)

    def plot_measurement(
        self,
        *,
        measurement=None,
        axes=None,
        mass_list=None,
        mass_lists=None,
        mol_list=None,
        mol_lists=None,
        tspan=None,
        tspan_bg=None,
        remove_background=None,
        unit=None,
        x_unit=None,
        meta_units=None,
        meta_lists=None,
        meta_list=None,
        T_name=None,
        P_name=None,
        T_color="k",
        P_color="r",
        logplot=None,
        legend=True,
        xspan=None,
        cmap_name="inferno",
        make_colorbar=False,
        aspect=1.25,
        max_threshold=None,
        min_threshold=None,
        scanning_mask=None,
        _sort_indicies=None,
        vmin=None,
        vmax=None,
        emphasis="middle",
        **kwargs,
    ):
        """Make a spectro TP-MS plot vs time and return the axis handles.

        Allocates some tasks to TPMSPlotter.plot_measurement()

        Args:
            measurement (SpectroReactorMeasurement): Defaults to the measurement to which
                the plotter is bound (self.measurement)
            axes (list of four matplotlib axes): axes[0] plots the spectral, axes[1] MS,
                axes[2] the variable given by T_name (temperature), and axes[4] the
                variable given by P_name (reactor pressure). By default four axes are made
                with axes[0] a top panel, axes[1] a middle panel, axes[2] and axes[4]
                the left and right yaxes of the bottom panel
            mass_list (list of str): The names of the m/z values, eg. ["M2", ...] to
                plot. Defaults to all of them (measurement.mass_list)
            mass_lists (list of list of str): Alternately, two lists can be given for
                masses in which case one list is plotted on the left y-axis and the other
                on the right y-axis of the top panel.
            mol_list (list of str): The names of the molecules, eg. ["H2", ...] to
                plot. Defaults to all of them (measurement.mass_list)
            mol_lists (list of list of str): Alternately, two lists can be given for
                molecules in which case one list is plotted on the left y-axis and the
                other on the right y-axis of the top panel.
            tspan (iter of float): The time interval to plot, wrt measurement.tstamp
            tspan_bg (timespan): A timespan for which to assume the signal is at its
                background. The average signals during this timespan are subtracted.
                If `mass_lists` are given rather than a single `mass_list`, `tspan_bg`
                must also be two timespans - one for each axis. Default is `None` for no
                background subtraction.
            remove_background (bool): Whether otherwise to subtract pre-determined
                background signals if available. Defaults to (not logplot)
            unit (str): the unit for the MS data. Defaults to "A" for Ampere
            T_name (str): The name of the value to plot on the lower left y-axis.
                Defaults to the name of the series `measurement.temperature`
            P_name (str): The name of the value to plot on the lower right y-axis.
                Defaults to the name of the series `measurement.pressure`
            T_color (str): The color to plot the variable given by 'T_name'
            P_color (str): The color to plot the variable given by 'P_name'
            logplot (bool): Whether to plot the MS data on a log scale (default True
                unless mass_lists are given)
            legend (bool): Whether to use a legend for the MID data (default True)
            xspan (iterable): The span of the spectral data to plot
            cmap_name (str): The name of the colormap to use. Defaults to "inferno", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html#sequential
            make_colorbar (bool): Whether to make a colorbar.
                FIXME: colorbar at present mis-alignes axes
            aspect (float): aspect ratio. Defaults to 1.25 times taller than wide.
            max_threshold (float): Set maximum value in scanning data and set to zero if
                data is above.
            min_threshold (float): Set minimum value in scanning data and set to zero if
                data is below.
            scanning_mask (boolean list): List of booleans to include/ exclude specfic
                data in scanning plot (specific masses that are monitored otherwise)
            _sort_indicies (list): list of floats to sort data. Defaults low to high.
            vmin (float): Value used to shift colours in the colorbar to lower values
            vmax (float): Value to shift colours in the colorbar to higher values
            kwargs (dict): Additional kwargs go to all calls of matplotlib's plot()

        Returns:
            list of Axes: (top, mid_left, bottom_left, mid_right, bottom_right) where:
                axes[0] is top is MS spectral data
                axes[1] is mid_left is Mass ID data;
                axes[2] is bottom_left is temperature;
                axes[3] is mid_right is additional MS data if left and right mass_lists
                    or mol_lists were plotted (otherwise axes[3] is None); and
                axes[4] is bottom_right is pressure.
        """
        measurement = measurement or self.measurement

        if not axes:
            axes = self.new_three_panel_axes(
                n_top=1, n_middle=(2 if (mass_lists or mol_lists) else 1), n_bottom=2
            )

        measurement.spectrum_series.heat_plot(
            ax=axes[0],
            tspan=tspan,
            xspan=xspan,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            max_threshold=max_threshold,
            min_threshold=min_threshold,
            scanning_mask=scanning_mask,
            _sort_indicies=_sort_indicies,
            vmin=vmin,
            vmax=vmax,
        )

        self.tpms_plotter.plot_measurement(
            measurement=measurement,
            axes=[axes[1], axes[2], axes[4], axes[5]],
            tspan=tspan,
            tspan_bg=tspan_bg,
            remove_background=remove_background,
            mass_list=mass_list,
            mass_lists=mass_lists,
            mol_list=mol_list,
            mol_lists=mol_lists,
            unit=unit,
            logplot=logplot,
            legend=legend,
            T_name=T_name,
            P_name=P_name,
            T_color=T_color,
            P_color=P_color,
            **kwargs,
        )

        axes[0].set_xlim(axes[1].get_xlim())

        fig = axes[0].get_figure()
        fig.set_figheight(fig.get_figwidth() * aspect)

        return axes

    def plot_measurement_vs(
        self,
        *,
        x_name,
        measurement=None,
        axes=None,
        mass_list=None,
        mass_lists=None,
        mol_list=None,
        mol_lists=None,
        tspan=None,
        tspan_bg=None,
        remove_background=None,
        unit=None,
        logplot=True,
        legend=True,
        xspan=None,
        cmap_name="inferno",
        make_colorbar=False,
        emphasis="top",
        ms_data="top",
        max_threshold=None,
        min_threshold=None,
        scanning_mask=None,
        _sort_indicies=None,
        **kwargs,
    ):

        if logplot is None:
            logplot = not mol_lists and not mass_lists

        if not axes:
            if ms_data == "top":
                n_bottom = 1
                n_top = 2 if (mass_lists or mol_lists) else 1
                ms_axes = 0
                ms_spec_axes = 1
            else:
                n_top = 1
                n_bottom = 2 if (mass_lists or mol_lists) else 1
                ms_axes = 1
                ms_spec_axes = 0

            axes = self.new_two_panel_axes(
                n_bottom=n_bottom,
                n_top=n_top,
                emphasis=emphasis,
            )

        measurement = measurement or self.measurement

        if (
            mass_list
            or mass_lists
            or mol_list
            or mol_lists
            or hasattr(measurement, "mass_list")
        ):
            # then we have MS data!
            self.tpms_plotter.ms_plotter.plot_vs(
                x_name=x_name,
                measurement=measurement,
                axes=[axes[ms_axes], axes[2]]
                if (mass_lists or mol_lists)
                else [axes[ms_axes]],
                tspan=tspan,
                tspan_bg=tspan_bg,
                remove_background=remove_background,
                mass_list=mass_list,
                mass_lists=mass_lists,
                mol_list=mol_list,
                mol_lists=mol_lists,
                unit=unit,
                logplot=logplot,
                legend=legend,
                **kwargs,
            )

        _tseries = measurement.spectrum_series.field.axes_series[0]
        _v = measurement.grab_for_t(item=x_name, t=_tseries.t)
        if not _sort_indicies:
            _sort_indicies = np.argsort(_v)

        measurement.spectrum_series.heat_plot(
            ax=axes[ms_spec_axes],
            t=_v[_sort_indicies],
            tspan=tspan,
            xspan=xspan,
            t_name=x_name,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            max_threshold=max_threshold,
            min_threshold=min_threshold,
            scanning_mask=scanning_mask,
            _sort_indicies=_sort_indicies,
            **kwargs,
        )

        return axes


def _get_y_unit_and_label(data_series, meta_units):
    name = data_series.name
    if "temperatur" in name.lower():
        ylabel = "temperature"
    elif "flow" in name.lower():
        ylabel = "flow rate"
    elif "pressure" in name.lower():
        ylabel = "pressure"
    elif "current" in name.lower():
        ylabel = "current"
    elif "voltage" in name.lower():
        ylabel = "voltage"
    elif "power" in name.lower():
        ylabel = "power"
    else:
        ylabel = " "

    new_unit_name = data_series.unit.name
    if meta_units:
        try:
            new_unit_name = meta_units[name]
        except KeyError:
            pass

    y_unit_factor, y_unit_name = _get_unit_factor_and_name(
            new_unit_name = new_unit_name,
            from_unit_name = data_series.unit.name,
        )

    return ylabel, y_unit_name, y_unit_factor


def _get_unit_factor_and_name(
    new_unit_name,
    from_unit_name,
):
    try:
        if from_unit_name == "celsius" or from_unit_name == "C":
            warnings.warn(
                "Temperature is not factorial converted and should be done in technique"
                " method prior to plotting with this plotter. ",
                stacklevel=2,
            )
            unit_factor = 1
            new_unit_name = ' C' #from_unit_name

        elif from_unit_name == "kelvin" or from_unit_name == "K":
            warnings.warn(
                "Temperature is not factorial converted and should be done in technique"
                " method prior to plotting with this plotter. ",
                stacklevel=2,
            )
            unit_factor = 1
            new_unit_name = ' K'#from_unit_name

        elif from_unit_name == "mbar":
            unit_factor = {
                "mbar": 1,
                "bar": 1e-3,
                "hPa": 1,
                "kPa": 0.1e-3,
            }[new_unit_name]

        elif from_unit_name == "bar":
            unit_factor = {
                "mbar": 1e3,
                "bar": 1,
                "hPa": 1e3,
                "kPa": 0.1e3,
            }[new_unit_name]

        else:
            unit_factor = {
                # Time conversion
                "s": 1,
                "min": 1 / 60,
                "minutes": 1 / 60,
                "h": 1 / 3600,
                "hr": 1 / 3600,
                "hour": 1 / 3600,
                "hours": 1 / 3600,
                "d": 1 / (3600 * 24),
                "days": 1 / (3600 * 24),
                # Pressure conversion
                "mbar": 1,
                "bar": 1000,
                "hPa": 1,
                "kPa": 0.1,
            }[new_unit_name]
    except KeyError:
        warnings.warn(
            f"Can't convert original unit '{from_unit_name}' to new unit"
            f"'{new_unit_name}'. Plotting using original unit!",
            stacklevel=2,
        )
        unit_factor = 1
        new_unit_name = from_unit_name
    return unit_factor, new_unit_name


#  ----- These are the standard colors for TP-MS plots! ------- #

MIN_SIGNAL = 1e-14  # So that the bottom half of the plot isn't wasted on log(noise)
# TODO: This should probably be customizeable from a settings file.

STANDARD_COLORS = {
    # Inset of meta channels #
    "TC temperature": "#000075",#"#808000",
    "RTD temperature": "#4363d8",
    "Reactor pressure": "#808000",
    "Baratron pressure": "#808000",
    "Containment pressure": "#808000",
    "Flow1": "k",
    "Flow2": "brown",
    "Flow3": "c",
    "Flow4": "b",
    "Flow5": "r",
    "Flow6": "0.5",
    # Inset for anodic bonding #
    "TC anodic bonding (top)": "#000075",  # "#808000",
    "TC anodic bonding (bottom)": "#4363d8",  # "#9A6324",
    "Total power": "#800000",
    "Heater voltage 1": "#fabed4",
    "Heater current 1": "#ffd8b1",
}
