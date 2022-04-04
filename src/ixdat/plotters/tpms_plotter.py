from .base_mpl_plotter import MPLPlotter
from .ms_plotter import MSPlotter
from .plotting_tools import color_axis


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
        T_name=None,
        P_name=None,
        T_color="k",
        P_color="r",
        logplot=None,
        legend=True,
        emphasis="top",
        **kwargs,
    ):
        """Make a TP-MS plot vs time and return the axis handles.

        Allocates some tasks to MSPlotter.plot_measurement()

        Args:
            measurement (TPMSMeasurement): Defaults to the measurement to which the
                plotter is bound (self.measurement)
            axes (list of three matplotlib axes): axes[0] plots the MID data,
                axes[1] the variable given by T_str (temperature), and axes[2] the
                variable given by P_str (reactor pressure). By default three axes are made with
                axes[0] a top panel with 3/5 the area, and axes[1] and axes[3] are
                the left and right y-axes of the lower panel with 2/5 the area.
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
            T_color (str): The color to plot the variable given by 'T_str'
            P_color (str): The color to plot the variable given by 'P_str'
            logplot (bool): Whether to plot the MS data on a log scale (default True
                unless mass_lists are given)
            legend (bool): Whether to use a legend for the MS data (default True)
            emphasis (str or None): "top" for bigger top panel, "bottom" for bigger
                bottom panel, None for equal-sized panels
            kwargs (dict): Additional kwargs go to all calls of matplotlib's plot()

        Returns:
            list of Axes: (top_left, bottom_left, top_right, bottom_right) where:
                axes[0] is top_left is MS data;
                axes[1] is bottom_left is temperature;
                axes[2] is top_right is additional MS data if left and right mass_lists
                    or mol_lists were plotted (otherwise axes[2] is None); and
                axes[3] is bottom_right is pressure.
        """

        if logplot is None:
            logplot = not mol_lists and not mass_lists

        if not axes:
            axes = self.new_two_panel_axes(
                n_bottom=2,
                n_top=(2 if (mass_lists or mol_lists) else 1),
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
                logplot=logplot,
                legend=legend,
                **kwargs,
            )

        T_name = T_name or measurement.T_name
        P_name = P_name or measurement.P_name

        t_T, T = measurement.grab(T_name, tspan=tspan)
        axes[1].plot(t_T, T, color=T_color)
        axes[1].set_xlabel("time / [s]")
        axes[1].set_ylabel(T_name)

        t_P, P = measurement.grab(P_name, tspan=tspan)
        axes[3].plot(t_P, P, color=P_color)
        axes[3].set_xlabel("time / [s]")
        axes[3].set_ylabel(P_name)
        color_axis(axes[3], P_color)

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
        **kwargs,
    ):
        """Make a spectro TP-MS plot vs time and return the axis handles.

        Allocates some tasks to TPMSPlotter.plot_measurement()

        Args:
            measurement (SpecroReactorMeasurement): Defaults to the measurement to which
                the plotter is bound (self.measurement)
            axes (list of four matplotlib axes): axes[0] plots the spectral, axes[1] MS,
                axes[2] the variable given by T_str (temperature), and axes[4] the
                variable given by P_str (reactor pressure). By default four axes are made
                with axes[0] a top panel, axes[1] a middle panel, axes[2] and axes[4]
                the left and right yaxes of the bottom panel
                are the left and right y-axes of the lower panel with 2/5 the area.
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
            T_color (str): The color to plot the variable given by 'T_str'
            P_color (str): The color to plot the variable given by 'P_str'
            logplot (bool): Whether to plot the MS data on a log scale (default True
                unless mass_lists are given)
            legend (bool): Whether to use a legend for the MID data (default True)
            xspan (iterable): The span of the spectral data to plot
            cmap_name (str): The name of the colormap to use. Defaults to "inferno", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html#sequential
            make_colorbar (bool): Whether to make a colorbar.
                FIXME: colorbar at present mis-alignes axes
            aspect (float): aspect ratio. Defaults to 1.25 times taller than wide.
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
