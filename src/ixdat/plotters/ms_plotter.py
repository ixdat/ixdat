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
        mol_list=None,
        mol_lists=None,
        tspan=None,
        tspan_bg=None,
        removebackground=None,
        unit=None,
        logplot=True,
        legend=True,
        **kwargs,
    ):
        """Plot m/z signal vs time (MID) data and return the axis.

        There are four ways to specify what to plot. Only specify one of these::
            mass_list: Uncalibrated signals in [(u/n/p)A] on on axis
            mass_lists: Uncalibrated signals in [(u/n/p)A] on two axes
            mol_list: Calibrated signals in [(u/n/p)mol/s] on on axis
            mol_lists: Calibrated signals in [(u/n/p)mol/s] on two axes

        Two axes refers to separate left and right y-axes. Default is to use all
        available masses as mass_list.

        Args:
            measurement (MSMeasurement): defaults to the one that initiated the plotter
            ax (matplotlib axis): Defaults to a new axis
            axes (list of matplotlib axis): Left and right y-axes if mass_lists are given
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
            unit (str): defaults to "A" or "mol/s"
            logplot (bool): Whether to plot the MS data on a log scale (default True)
            legend (bool): Whether to use a legend for the MS data (default True)
            kwargs: extra key-word args are passed on to matplotlib's plot()
        """
        measurement = measurement or self.measurement
        if removebackground is None:
            removebackground = not logplot

        quantified, specs_this_axis, specs_next_axis = self._parse_overloaded_inputs(
            mass_list,
            mass_lists,
            mol_list,
            mol_lists,
            unit,
            tspan_bg,
            ax,
            axes,
            measurement,
        )
        ax = specs_this_axis["ax"]
        v_list = specs_this_axis["v_list"]
        tspan_bg = specs_this_axis["tspan_bg"]
        unit = specs_this_axis["unit"]
        unit_factor = specs_this_axis["unit_factor"]
        for v_or_v_name in v_list:
            if isinstance(v_or_v_name, str):
                v_name = v_or_v_name
                color = STANDARD_COLORS.get(v_name, "k")
            else:
                v_name = v_or_v_name.name
                color = v_or_v_name.color
            if quantified:
                t, v = measurement.grab_flux(
                    v_or_v_name,
                    tspan=tspan,
                    tspan_bg=tspan_bg,
                    removebackground=removebackground,
                    include_endpoints=False,
                )
            else:
                t, v = measurement.grab_signal(
                    v_or_v_name,
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
                color=color,
                label=v_name,
                **kwargs,
            )
        ax.set_ylabel(f"signal / [{unit}]")
        ax.set_xlabel("time / [s]")
        if specs_next_axis:
            self.plot_measurement(
                measurement=measurement,
                ax=specs_next_axis["ax"],
                mass_list=specs_next_axis["mass_list"],
                mol_list=specs_next_axis["mol_list"],
                unit=specs_next_axis["unit"],
                tspan=tspan,
                tspan_bg=specs_next_axis["tspan_bg"],
                logplot=logplot,
                legend=legend,
                **kwargs,
            )
            axes = [ax, specs_next_axis["ax"]]
        else:
            axes = None

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
        mol_list=None,
        mol_lists=None,
        tspan=None,
        tspan_bg=None,
        removebackground=None,
        unit=None,
        logplot=True,
        legend=True,
        **kwargs,
    ):
        """Plot m/z signal (MID) data against a specified variable and return the axis.

        There are four ways to specify what to plot. Only specify one of these::
            mass_list: Uncalibrated signals in [(u/n/p)A] on on axis
            mass_lists: Uncalibrated signals in [(u/n/p)A] on two axes
            mol_list: Calibrated signals in [(u/n/p)mol/s] on on axis
            mol_lists: Calibrated signals in [(u/n/p)mol/s] on two axes

        Two axes refers to seperate left and right y-axes. Default is to use all
        available masses as mass_list.

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
                background signals if available
            logplot (bool): Whether to plot the MS data on a log scale (default True)
            legend (bool): Whether to use a legend for the MS data (default True)
            kwargs: key-word args are passed on to matplotlib's plot()
        """
        measurement = measurement or self.measurement
        if removebackground is None:
            removebackground = not logplot

        # The overloaded inputs are a pain in the ass. This function helps:
        quantified, specs_this_axis, specs_next_axis = self._parse_overloaded_inputs(
            mass_list,
            mass_lists,
            mol_list,
            mol_lists,
            unit,
            tspan_bg,
            ax,
            axes,
            measurement,
        )
        ax = specs_this_axis["ax"]
        v_list = specs_this_axis["v_list"]
        tspan_bg = specs_this_axis["tspan_bg"]
        unit = specs_this_axis["unit"]
        unit_factor = specs_this_axis["unit_factor"]

        t, x = measurement.grab(x_name, tspan=tspan, include_endpoints=True)
        for v_name in v_list:
            if quantified:
                t_v, v = measurement.grab_flux(
                    v_name,
                    tspan=tspan,
                    tspan_bg=tspan_bg,
                    removebackground=removebackground,
                    include_endpoints=False,
                )
            else:
                t_v, v = measurement.grab_signal(
                    v_name,
                    tspan=tspan,
                    t_bg=tspan_bg,
                    removebackground=removebackground,
                    include_endpoints=False,
                )
            if logplot:
                v[v < MIN_SIGNAL] = MIN_SIGNAL
            x_mass = np.interp(t_v, t, x)
            ax.plot(
                x_mass,
                v * unit_factor,
                color=STANDARD_COLORS.get(v_name, "k"),
                label=v_name,
                **kwargs,
            )
        ax.set_ylabel(f"signal / [{unit}]")
        ax.set_xlabel(x_name)
        if specs_next_axis:
            self.plot_vs(
                x_name=x_name,
                measurement=measurement,
                ax=specs_next_axis["ax"],
                mass_list=specs_next_axis["mass_list"],
                mol_list=specs_next_axis["mol_list"],
                unit=specs_next_axis["unit"],
                tspan=tspan,
                tspan_bg=specs_next_axis["tspan_bg"],
                logplot=logplot,
                legend=legend,
                **kwargs,
            )
            axes = [ax, specs_next_axis["ax"]]
        else:
            axes = None

        if logplot:
            ax.set_yscale("log")
        if legend:
            ax.legend()

        return axes if axes else ax

    def _parse_overloaded_inputs(
        self,
        mass_list,
        mass_lists,
        mol_list,
        mol_lists,
        unit,
        tspan_bg,
        ax,
        axes,
        measurement,
    ):
        """From the overloaded function inputs, figure out what the user wants to do.

        This includes:
        1. determine if we're doing quantifed results (mols) or raw (masses)
        2. figure out if there's one or two axes (remaining) and what goes on them.
        3. figure out what to multiply numbers by when plotting to match the unit.
        """
        # TODO: Maybe there's a way to do this function as a decorator?
        # So this function is overloaded in the sense that the user can give
        #   exactly one of mol_list, mol_lists, mass_list, mass_lists.
        # To manage that complexity, first we reduce it to two options, that down to
        #   either v_list or v_lists and a boolean "quantified":
        quantified = False  # default, if they give nothing
        v_lists = None  # default, if they give nothing
        v_list = measurement.mass_list  # default, if they give nothing
        if mol_list:
            quantified = True
            v_list = mol_list
        elif mol_lists:
            quantified = True
            v_lists = mol_lists
        elif mass_list:
            quantified = False
            v_list = mass_list
        elif mass_lists:
            quantified = False
            v_lists = mass_lists

        # as the next simplification, if they give two things (v_lists), we pretend we
        #   got one (v_list) but prepare an axis for a recursive call of this function.
        if v_lists:
            axes = axes or [ax, ax.twinx()]  # prepare an axis unless we were given two.
            ax_right = axes[-1]
            ax = axes[0]
            v_list = v_lists[0]
            v_list_right = v_lists[1]
            # ah, and to enable different background tspans for the two axes:
            try:
                tspan_bg_right = tspan_bg[1]
                if isinstance(tspan_bg_right, (float, int)):
                    raise TypeError
            except (KeyError, TypeError):
                tspan_bg_right = None
            else:
                tspan_bg = tspan_bg[0]
            if isinstance(unit, str) or not unit:
                unit_right = unit
            else:
                unit_right = unit[1]
                unit = unit[0]
            specs_next_axis = {
                "ax": ax_right,
                "unit": unit_right,
                "mass_list": None if quantified else v_list_right,
                "mol_list": v_list_right if quantified else None,
                "tspan_bg": tspan_bg_right,
            }
        else:
            specs_next_axis = None

        if quantified:
            unit = unit or "mol/s"
            unit_factor = {
                "pmol/s": 1e12,
                "nmol/s": 1e9,
                "umol/s": 1e6,
                "mol/s": 1,  # noqa
                "pmol/s/cm^2": 1e12,
                "nmol/s/cm^2": 1e9,
                "umol/s/cm^2": 1e6,
                "mol/s/cm^2": 1,  # noqa
            }[unit]
            if "/cm^2" in unit:
                unit_factor = unit_factor / measurement.A_el
        else:
            unit = unit or "A"
            unit_factor = {"pA": 1e12, "nA": 1e9, "uA": 1e6, "A": 1}[unit]
        # TODO: Real units with a unit module! This should even be able to figure out the
        #  unit prefix to put stuff in a nice 1-to-1e3 range

        if not ax:
            ax = (
                axes[0]
                if axes
                else self.new_ax(ylabel=f"signal / [{unit}]", xlabel="time / [s]")
            )
        specs_this_axis = {
            "ax": ax,
            "v_list": v_list,
            "unit": unit,
            "unit_factor": unit_factor,
            "tspan_bg": tspan_bg,
        }

        return quantified, specs_this_axis, specs_next_axis


#  ----- These are the standard colors for EC-MS plots! ------- #

MIN_SIGNAL = 1e-14  # So that the bottom half of the plot isn't wasted on log(noise)

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
    # and now, molecules:
    "H2": "b",
    "He": "m",
    "H2O": "y",
    "CO": "0.5",
    "N2": "0.5",
    "O2": "k",
    "Ar": "c",
    "CO2": "brown",
    "CH4": "r",
    "C2H4": "g",
    "O2_M32": "k",
    "O2_M34": "r",
    "O2_M36": "g",
    "CO2_M44": "brown",
    "CO2_M46": "purple",
    "CO2_M48": "darkslategray",
}
