"""Plotter for Electrochemistry - Mass Spectrometry"""

from .base_mpl_plotter import MPLPlotter
from .ec_plotter import ECPlotter
from .ms_plotter import MSPlotter
from ..tools import deprecate


class ECMSPlotter(MPLPlotter):
    """A matplotlib plotter for EC-MS measurements."""

    def __init__(self, measurement=None):
        """Initiate the ECMSPlotter with its default Measurement to plot"""
        super().__init__()
        self.measurement = measurement
        self.ec_plotter = ECPlotter(measurement=measurement)
        self.ms_plotter = MSPlotter(measurement=measurement)

    @deprecate("0.1", "Use `U_name` instead.", "0.3", kwarg_name="V_str")
    @deprecate("0.1", "Use `J_name` instead.", "0.3", kwarg_name="J_str")
    @deprecate("0.1", "Use `U_color` instead.", "0.3", kwarg_name="V_color")
    @deprecate(
        "0.1", "Use `remove_background` instead.", "0.3", kwarg_name="removebackground"
    )
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
        removebackground=None,
        unit=None,
        U_name=None,
        J_name=None,
        U_color="k",
        J_color="r",
        V_str=None,
        J_str=None,
        V_color=None,
        logplot=None,
        legend=True,
        emphasis="top",
        **kwargs,
    ):
        """Make an EC-MS plot vs time and return the axis handles.

        Allocates tasks to ECPlotter.plot_measurement() and MSPlotter.plot_measurement()

        Args:
            measurement (ECMSMeasurement): Defaults to the measurement to which the
                plotter is bound (self.measurement)
            axes (list of matplotlib axes): axes[0] plots the MID data,
                axes[1] the variable given by `J_name` (potential), and axes[3] the
                variable given by `J_name` (current). By default three axes are made
                with axes[0] a top panel with 3/5 the area, and axes[1] and axes[3] are
                the left and right y-axes of the lower panel with 2/5 the area.
                axes[2], typically the top right panel, will only be used if two MS
                axes are requested (see `mass_lists` and `mol_lists`).
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
            tspan (str or iter of float): The time interval `(t_start, t_end)` to plot,
                wrt measurement.tstamp. Use `tspan="all"` to plot all the data.
                Defaults to EC data tspan. See :fun:`determine_tspan`_ for more details.
            tspan_bg (timespan): A timespan for which to assume the signal is at its
                background. The average signals during this timespan are subtracted.
                If `mass_lists` are given rather than a single `mass_list`, `tspan_bg`
                must also be two timespans - one for each axis. Default is `None` for no
                background subtraction.
            remove_background (bool): Whether otherwise to subtract pre-determined
                background signals if available. Defaults to (not logplot)
            removebackground (bool): DEPRECATED. Use `remove_background`
            unit (str): the unit for the MS data. Defaults to "A" for Ampere
            U_name (str): The name of the value to plot on the lower left y-axis.
                Defaults to the name of the series `measurement.potential`
            J_name (str): The name of the value to plot on the lower right y-axis.
                Defaults to the name of the series `measurement.current`
            U_color (str): The color to plot the variable given by 'V_str'
            J_color (str): The color to plot the variable given by 'J_str'
            V_str (str): DEPRECATED. Use `U_name`.
            J_str (str): DEPRECATED. Use `J_name`.
            V_color (str): DEPRECATED. Use `U_color`.
            logplot (bool): Whether to plot the MS data on a log scale (default True
                unless mass_lists are given)
            legend (bool): Whether to use a legend for the MS data (default True)
            emphasis (str or None): "top" for bigger top panel, "bottom" for bigger
                bottom panel, None for equal-sized panels
            kwargs (dict): Additional kwargs go to all calls of matplotlib's plot()

        Returns:
            list of Axes: (top_left, bottom_left, top_right, bottom_right) where:
                axes[0] is top_left is MS data;
                axes[1] is bottom_left is potential;
                axes[2] is top_right is additional MS data if left and right mass_lists
                    or mol_lists were plotted (otherwise axes[2] is None); and
                axes[3] is bottom_right is current.
        """
        measurement = measurement or self.measurement
        tspan = determine_tspan(tspan, measurement)

        logplot = (not mass_lists) if logplot is None else logplot

        # apply deprecated arguments (the user will get a warning):
        U_name = U_name or V_str
        J_name = J_name or J_str
        U_color = U_color or V_color
        if removebackground is not None:
            # note removebackground can be set to `False`
            remove_background = removebackground

        if not axes:
            axes = self.new_two_panel_axes(
                n_bottom=2,
                n_top=(2 if (mass_lists or mol_lists) else 1),
                emphasis=emphasis,
            )

        # plot the EC data (note that the EC plotter will skip current and/or potential
        #   in the case that the data is missing).
        self.ec_plotter.plot_measurement(
            measurement=measurement,
            axes=[axes[1], axes[3]],
            tspan=tspan,
            U_name=U_name,
            J_name=J_name,
            U_color=U_color,
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
                axes=[axes[0], axes[2]] if (mass_lists or mol_lists) else axes[0],
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
        axes[1].set_xlim(axes[0].get_xlim())
        return axes

    @deprecate(
        "0.1", "Use `remove_background` instead.", "0.3", kwarg_name="removebackground"
    )
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
        remove_background=None,
        removebackground=None,
        unit=None,
        logplot=False,
        legend=True,
        emphasis="top",
        **kwargs,
    ):
        """Make an EC-MS plot vs time and return the axis handles.

        Allocates tasks to ECPlotter.plot_measurement() and MSPlotter.plot_measurement()

        Args:
            measurement (ECMSMeasurement): Defaults to the measurement to which the
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
            tspan (str or iter of float): The time interval `(t_start, t_end)` to plot,
                wrt measurement.tstamp. Use `tspan="all"` to plot all the data.
                Defaults to EC data tspan. See :fun:`determine_tspan`_ for more details.
            tspan_bg (timespan): A timespan for which to assume the signal is at its
                background. The average signals during this timespan are subtracted.
                If `mass_lists` are given rather than a single `mass_list`, `tspan_bg`
                must also be two timespans - one for each axis. Default is `None` for no
                background subtraction.
            remove_background (bool): Whether otherwise to subtract pre-determined
                background signals if available. Defaults to (not logplot)
            removebackground: DEPRECATED. Use `remove_background`
            unit (str): the unit for the MS data. Defaults to "A" for Ampere
            logplot (bool): Whether to plot the MS data on a log scale (default False)
            legend (bool): Whether to use a legend for the MS data (default True)
            emphasis (str or None): "top" for bigger top panel, "bottom" for bigger
                bottom panel, None for equal-sized panels
            kwargs (dict): Additional kwargs go to all calls of matplotlib's plot()
        """
        measurement = measurement or self.measurement
        tspan = determine_tspan(tspan, measurement)
        if removebackground is not None:
            remove_background = removebackground
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
        axes[1].set_xlim(axes[0].get_xlim())
        return axes


def determine_tspan(tspan=None, measurement=None):
    """Return a timespan based on a requested timespan and the EC-MS measurement

    `tspan` is most directly given as two numbers, which are taken to be
    start time and end time (wrt to measurement.tstamp), e.g. `[800, 1000]`.
    `tspan` can also be one of these strings:
        "all": The full timespan from the first data point to the last
        "ec": The full timespan of the EC data.
        "ms": The full timespan of the MS data.
    If tspan is not used, the measurement's default timespan `measurement.tspan` is
    used, which for EC-MS measurements will be "ec".

    Args:
         tspan (str or iter of float): The timespan specification
         measurement (ECMSMeasurement): The measurement for which to determine a timespan
    """
    if isinstance(tspan, str):
        ec_start = measurement.t[0]
        ec_end = measurement.t[-1]
        ms_start = None
        ms_end = None
        for mass in measurement.mass_list:
            t_M, y = measurement.grab(mass)
            ms_start = t_M[0] if ms_start is None else min(ms_start, t_M[0])
            ms_end = t_M[-1] if ms_end is None else max(ms_end, t_M[-1])
        t_start = min(ec_start, ms_start)
        t_end = max(ec_end, ms_end)
        if tspan == "all":
            return [t_start, t_end]
        elif tspan == "ec":
            return [ec_start, ec_end]
        elif tspan == "ms":
            return [ms_start, ms_end]
        else:
            raise ValueError(
                f"Invalid tspan: '{tspan}'. String options are 'all', 'ec', and 'ms'."
            )
    elif tspan:
        return tspan
    else:
        return measurement.tspan
