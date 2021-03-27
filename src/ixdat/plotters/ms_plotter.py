import numpy as np
from .base_mpl_plotter import MPLPlotter


class MSPlotter(MPLPlotter):
    """A matplotlib plotter specialized in mass spectrometry MID measurements."""

    def __init__(self, measurement=None):
        """Initiate the ECMSPlotter with its default Meausurement to plot"""
        self.measurement = measurement

    def plot_measurement(
        self,
        *,
        measurement=None,
        ax=None,
        axes=None,
        mass_list=None,
        mass_lists=None,
        tspan=None,
        tspan_bg=None,
        removebackground=None,
        unit="A",
        logplot=True,
        legend=True,
        **kwargs,
    ):
        """Plot m/z signal vs time (MID) data and return the axis.

        Args:
            measurement (MSMeasurement): defaults to the one that initiated the plotter
            ax (matplotlib axis): Defaults to a new axis
            axes (list of matplotlib axis): Left and right y-axes if mass_lists are given
            mass_list (list of str): The names of the m/z values, eg. ["M2", ...] to
                plot. Defaults to all of them (measurement.mass_list)
            mass_lists (list of list of str): Alternately, two lists can be given for
                masses in which case one list is plotted on the left y-axis and the other
                on the right y-axis of the top panel.
            tspan (iter of float): The time interval to plot, wrt measurement.tstamp
            tspan_bg (timespan): A timespan for which to assume the signal is at its
                background. The average signals during this timespan are subtracted.
                If `mass_lists` are given rather than a single `mass_list`, `tspan_bg`
                must also be two timespans - one for each axis. Default is `None` for no
                background subtraction.
            removebackground (bool): Whether otherwise to subtract pre-determined
                background signals if available. Defaults to (not logplot)
            logplot (bool): Whether to plot the MS data on a log scale (default True)
            legend (bool): Whether to use a legend for the MS data (default True)
            kwargs: extra key-word args are passed on to matplotlib's plot()
        """
        measurement = measurement or self.measurement
        if removebackground is None:
            removebackground = not logplot
        if not ax:
            ax = (
                axes[0]
                if axes
                else self.new_ax(ylabel=f"signal / [{unit}]", xlabel="time / [s]")
            )
        tspan_bg_right = None
        if mass_lists:
            axes = axes or [ax, ax.twinx()]
            ax = axes[0]
            mass_list = mass_lists[0]
            try:
                tspan_bg_right = tspan_bg[1]
                if isinstance(tspan_bg_right, (float, int)):
                    raise TypeError
            except (KeyError, TypeError):
                tspan_bg_right = None
            else:
                tspan_bg = tspan_bg[0]
        unit_factor = {"pA": 1e12, "nA": 1e9, "uA": 1e6, "A": 1}[unit]
        # TODO: Real units with a unit module! This should even be able to figure out the
        #  unit prefix to put stuff in a nice 1-to-1e3 range

        mass_list = mass_list or measurement.mass_list
        for mass in mass_list:
            t, v = measurement.grab_signal(
                mass,
                tspan=tspan,
                t_bg=tspan_bg,
                removebackground=removebackground,
                include_endpoints=False,
            )
            if logplot:
                v[v < MIN_SIGNAL] = MIN_SIGNAL
            ax.plot(
                t,
                v * unit_factor,
                color=STANDARD_COLORS.get(mass, "k"),
                label=mass,
                **kwargs,
            )
        if mass_lists:
            self.plot_measurement(
                measurement=measurement,
                ax=axes[1],
                mass_list=mass_lists[1],
                unit=unit,
                tspan=tspan,
                tspan_bg=tspan_bg_right,
                logplot=logplot,
                legend=legend,
                **kwargs,
            )

        if logplot:
            ax.set_yscale("log")
        if legend:
            ax.legend()

        return axes if axes else ax

    def plot_vs(
        self,
        *,
        x_name,
        measurement=None,
        ax=None,
        axes=None,
        mass_list=None,
        mass_lists=None,
        tspan=None,
        tspan_bg=None,
        removebackground=None,
        unit="A",
        logplot=True,
        legend=True,
        **kwargs,
    ):
        """Plot m/z signal (MID) data against a specified variable and return the axis.

        Args:
            x_name (str): Name of the variable to plot on the x-axis
            measurement (MSMeasurement): defaults to the one that initiated the plotter
            ax (matplotlib axis): Defaults to a new axis
            axes (list of matplotlib axis): Left and right y-axes if mass_lists are given
            mass_list (list of str): The names of the m/z values, eg. ["M2", ...] to
                plot. Defaults to all of them (measurement.mass_list)
            mass_lists (list of list of str): Alternately, two lists can be given for
                masses in which case one list is plotted on the left y-axis and the other
                on the right y-axis of the top panel.
            tspan (iter of float): The time interval to plot, wrt measurement.tstamp
            tspan_bg (timespan): A timespan for which to assume the signal is at its
                background. The average signals during this timespan are subtracted.
                If `mass_lists` are given rather than a single `mass_list`, `tspan_bg`
                must also be two timespans - one for each axis. Default is `None` for no
                background subtraction.
            removebackground (bool): Whether otherwise to subtract pre-determined
                background signals if available
            logplot (bool): Whether to plot the MS data on a log scale (default True)
            legend (bool): Whether to use a legend for the MS data (default True)
            kwargs: key-word args are passed on to matplotlib's plot()
        """
        measurement = measurement or self.measurement
        if removebackground is None:
            removebackground = not logplot
        if not ax:
            ax = (
                axes[0]
                if axes
                else self.new_ax(ylabel=f"signal / [{unit}]", xlabel=x_name)
            )
        tspan_bg_right = None
        if mass_lists:
            axes = axes or [ax, ax.twinx()]
            ax = axes[0]
            mass_list = mass_lists[0]
            try:
                tspan_bg_right = tspan_bg[1]
                if isinstance(tspan_bg_right, (float, int)):
                    raise TypeError
            except (KeyError, TypeError):
                tspan_bg_right = None
            else:
                tspan_bg = tspan_bg[0]
        unit_factor = {"pA": 1e12, "nA": 1e9, "uA": 1e6, "A": 1}[unit]
        # TODO: Real units with a unit module! This should even be able to figure out the
        #  unit prefix to put stuff in a nice 1-to-1e3 ranges
        t, x = measurement.grab(x_name, tspan=tspan, include_endpoints=True)
        mass_list = mass_list or measurement.mass_list
        for mass in mass_list:
            t_mass, v = measurement.grab_signal(
                mass,
                tspan=tspan,
                t_bg=tspan_bg,
                removebackground=removebackground,
                include_endpoints=False,
            )
            if logplot:
                v[v < MIN_SIGNAL] = MIN_SIGNAL
            x_mass = np.interp(t_mass, t, x)
            ax.plot(
                x_mass,
                v * unit_factor,
                color=STANDARD_COLORS.get(mass, "k"),
                label=mass,
                **kwargs,
            )
        if mass_lists:
            self.plot_vs(
                x_name=x_name,
                measurement=measurement,
                ax=axes[1],
                mass_list=mass_lists[1],
                unit=unit,
                tspan=tspan,
                tspan_bg=tspan_bg_right,
                logplot=logplot,
                legend=legend,
                **kwargs,
            )

        if logplot:
            ax.set_yscale("log")
        if legend:
            ax.legend()

        return axes if axes else ax


MIN_SIGNAL = 1e-14  # So that the bottom half of the plot isn't wasted on log(noise)

#  ----- These are the standard colors for EC-MS plots! ------- #

STANDARD_COLORS = {
    "M2": "b",
    "M4": "m",
    "M18": "y",
    "M28": "0.5",
    "M32": "k",
    "M40": "c",
    "M44": "brown",
    "M15": "r",
    "M26": "g",
    "M27": "limegreen",
    "M30": "darkorange",
    "M31": "yellowgreen",
    "M43": "tan",
    "M45": "darkgreen",
    "M34": "r",
    "M36": "g",
    "M46": "purple",
    "M48": "darkslategray",
    "M20": "slateblue",
    "M16": "steelblue",
    "M19": "teal",
    "M17": "chocolate",
    "M41": "#FF2E2E",
    "M42": "olive",
    "M29": "#001146",
    "M70": "purple",
    "M3": "orange",
    "M73": "crimson",
    "M74": "r",
    "M60": "g",
    "M58": "darkcyan",
    "M88": "darkred",
    "M89": "darkmagenta",
    "M130": "purple",
    "M132": "purple",
}
