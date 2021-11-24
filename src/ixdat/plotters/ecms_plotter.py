from .base_mpl_plotter import MPLPlotter
from .ec_plotter import ECPlotter
from .ms_plotter import MSPlotter


class ECMSPlotter(MPLPlotter):
    """A matplotlib plotter for EC-MS measurements."""

    def __init__(self, measurement=None):
        """Initiate the ECMSPlotter with its default Meausurement to plot"""
        self.measurement = measurement
        self.ec_plotter = ECPlotter(measurement=measurement)
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
        removebackground=None,
        unit=None,
        V_str=None,
        J_str=None,
        V_color="k",
        J_color="r",
        logplot=None,
        legend=True,
        emphasis="top",
        **kwargs,
    ):
        """Make an EC-MS plot vs time and return the axis handles.

        Allocates tasks to ECPlotter.plot_measurement() and MSPlotter.plot_measurement()

        TODO: add all functionality in the legendary plot_experiment() in EC_MS.Plotting
            - variable subplot sizing (emphasizing EC or MS)
            - plotting of calibrated data (mol_list instead of mass_list)
            - units!

        Args:
            measurement (ECMSMeasurement): defaults to the measurement to which the
                plotter is bound (self.measurement)
            axes (list of three matplotlib axes): axes[0] plots the MID data,
                axes[1] the variable given by V_str (potential), and axes[2] the
                variable given by J_str (current). By default three axes are made with
                axes[0] a top panel with 3/5 the area, and axes[1] and axes[2] are
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
            removebackground (bool): Whether otherwise to subtract pre-determined
                background signals if available. Defaults to (not logplot)
            unit (str): the unit for the MS data. Defaults to "A" for Ampere
            V_str (str): The name of the value to plot on the lower left y-axis.
                Defaults to the name of the series `measurement.potential`
            J_str (str): The name of the value to plot on the lower right y-axis.
                Defaults to the name of the series `measurement.current`
            V_color (str): The color to plot the variable given by 'V_str'
            J_color (str): The color to plot the variable given by 'J_str'
            logplot (bool): Whether to plot the MS data on a log scale (default True
                unless mass_lists are given)
            legend (bool): Whether to use a legend for the MS data (default True)
            emphasis (str or None): "top" for bigger top panel, "bottom" for bigger
                bottom panel, None for equal-sized panels
            kwargs (dict): Additional kwargs go to all calls of matplotlib's plot()
        """
        measurement = measurement or self.measurement

        logplot = (not mass_lists) if logplot is None else logplot

        if not axes:
            axes = self.new_two_panel_axes(
                n_bottom=2,
                n_top=(2 if (mass_lists or mol_lists) else 1),
                emphasis=emphasis,
            )

        if not tspan:
            if hasattr(measurement, "potential") and measurement.potential:
                t, _ = measurement.grab("potential")
                tspan = [t[0], t[-1]]
            else:
                tspan = measurement.tspan

        if hasattr(measurement, "potential") and measurement.potential:
            # then we have EC data!
            self.ec_plotter.plot_measurement(
                measurement=measurement,
                axes=[axes[1], axes[2]],
                tspan=tspan,
                V_str=V_str,
                J_str=J_str,
                V_color=V_color,
                J_color=J_color,
                **kwargs,
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
                ax=axes[0],
                axes=[axes[0], axes[3]] if (mass_lists or mol_lists) else axes[0],
                tspan=tspan,
                tspan_bg=tspan_bg,
                removebackground=removebackground,
                mass_list=mass_list,
                mass_lists=mass_lists,
                mol_list=mol_list,
                mol_lists=mol_lists,
                unit=unit,
                logplot=logplot,
                legend=legend,
                **kwargs,
            )
        axes[1].set_xlim(axes[0].get_xlim())
        return axes

    def plot_vs_potential(
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
        removebackground=None,
        unit=None,
        logplot=False,
        legend=True,
        emphasis="top",
        **kwargs,
    ):
        """ "Make an EC-MS plot vs time and return the axis handles.

        Allocates tasks to ECPlotter.plot_measurement() and MSPlotter.plot_measurement()

        TODO: add all functionality in the legendary plot_experiment() in EC_MS.Plotting
            - variable subplot sizing (emphasizing EC or MS)
            - plotting of calibrated data (mol_list instead of mass_list)
            - units!
        
        Args:
            measurement (ECMSMeasurement): defaults to the measurement to which the
                plotter is bound (self.measurement)
            axes (list of three matplotlib axes): axes[0] plots the MID data,
                axes[1] the current vs potential. By default three axes are made with
                axes[0] a top panel with 3/5 the area.
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
            removebackground (bool): Whether otherwise to subtract pre-determined
                background signals if available. Defaults to (not logplot)
            unit (str): the unit for the MS data. Defaults to "A" for Ampere
            logplot (bool): Whether to plot the MS data on a log scale (default False)
            legend (bool): Whether to use a legend for the MS data (default True)
            emphasis (str or None): "top" for bigger top panel, "bottom" for bigger
                bottom panel, None for equal-sized panels
            kwargs (dict): Additional kwargs go to all calls of matplotlib's plot()
        """

        if not axes:
            axes = self.new_two_panel_axes(
                n_bottom=1,
                n_top=(2 if (mass_lists or mol_lists) else 1),
                emphasis=emphasis,
            )

        self.ec_plotter.plot_vs_potential(
            measurement=measurement, tspan=tspan, ax=axes[1], **kwargs
        )
        self.ms_plotter.plot_vs(
            x_name="potential",
            measurement=measurement,
            ax=axes[0],
            axes=[axes[0], axes[2]] if (mass_lists or mol_lists) else axes[0],
            tspan=tspan,
            tspan_bg=tspan_bg,
            removebackground=removebackground,
            mass_list=mass_list,
            mass_lists=mass_lists,
            mol_list=mol_list,
            mol_lists=mol_lists,
            unit=unit,
            logplot=logplot,
            legend=legend,
            **kwargs,
        )
        axes[1].set_xlim(axes[0].get_xlim())
        return axes
