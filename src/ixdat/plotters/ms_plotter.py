"""Plotter for Mass Spectrometry"""
import warnings
import numpy as np
from .base_mpl_plotter import MPLPlotter


class MSPlotter(MPLPlotter):
    """A matplotlib plotter specialized in mass spectrometry MID measurements."""

    def __init__(self, measurement=None):
        """Initiate the ECMSPlotter with its default Meausurement to plot"""
        super().__init__()
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
        remove_background=None,
        unit=None,
        x_unit=None,
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
            measurement (MSMeasurement): Defaults to the one that initiated the plotter
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
            remove_background (bool): Whether otherwise to subtract pre-determined
                background signals if available. Defaults to (not logplot)
            unit (str): defaults to "A" or "mol/s"
            x_unit (str): defaults to "s"
            logplot (bool): Whether to plot the MS data on a log scale (default True)
            legend (bool): Whether to use a legend for the MS data (default True)
            kwargs: extra key-word args are passed on to matplotlib's plot()
        """
        measurement = measurement or self.measurement
        if remove_background is None:
            remove_background = not logplot

        # Figure out, based on the inputs, whether or not to plot calibrated results
        # (`quantified`), specifications for the axis to plot on now (`specs_this_axis`)
        # and specifications for the next axis to plot on, if any (`specs_next_axis`):
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
                print("ender jeg nogensinde her")
                v_name = v_or_v_name.name
                color = v_or_v_name.color
            if quantified:
                t, v = measurement.grab_flux(
                    v_or_v_name,
                    tspan=tspan,
                    tspan_bg=tspan_bg,
                    remove_background=remove_background,
                    include_endpoints=False,
                )
            else:
                t, v = measurement.grab_signal(
                    v_or_v_name,
                    tspan=tspan,
                    tspan_bg=tspan_bg,
                    remove_background=remove_background,
                    include_endpoints=False,
                )
            if logplot:
                v[v < MIN_SIGNAL] = MIN_SIGNAL

            # expect always to plot against time
            x_unit_factor, x_unit = self._get_x_unit_factor(x_unit, "s")
            ax.plot(
                t * x_unit_factor,
                v * unit_factor,
                color=color,
                label=v_name,
                **kwargs,
            )
        ax.set_ylabel(f"signal / [{unit}]")
        ax.set_xlabel(f"time / [{x_unit}]")
        if specs_next_axis:
            self.plot_measurement(
                measurement=measurement,
                ax=specs_next_axis["ax"],
                mass_list=specs_next_axis["mass_list"],
                mol_list=specs_next_axis["mol_list"],
                unit=specs_next_axis["unit"],
                x_unit=x_unit,
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
        remove_background=None,
        unit=None,
        x_unit=None,
        logplot=True,
        legend=True,
        x_inverse=False,
        **plot_kwargs,
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
            measurement (MSMeasurement): Defaults to the one that initiated the plotter
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
            remove_background (bool): Whether otherwise to subtract pre-determined
                background signals if available
            unit (str): defaults to "A" or "mol/s"
            x_unit (str): defaults to x_name.unit.name
            logplot (bool): Whether to plot the MS data on a log scale (default True)
            legend (bool): Whether to use a legend for the MS data (default True)
            x_inverse (bool): Plot against inverse of x_name
            plot_kwargs: additional key-word args are passed on to matplotlib's plot()
        """
        measurement = measurement or self.measurement
        if remove_background is None:
            remove_background = not logplot

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
                    remove_background=remove_background,
                    include_endpoints=False,
                )
            else:
                t_v, v = measurement.grab_signal(
                    v_name,
                    tspan=tspan,
                    tspan_bg=tspan_bg,
                    remove_background=remove_background,
                    include_endpoints=False,
                )
            if logplot:
                v[v < MIN_SIGNAL] = MIN_SIGNAL
            x_mass = np.interp(t_v, t, x)
            plot_kwargs_this_mass = plot_kwargs.copy()
            if "color" not in plot_kwargs:
                plot_kwargs_this_mass["color"] = STANDARD_COLORS.get(v_name, "k")
            if "label" not in plot_kwargs:
                plot_kwargs_this_mass["label"] = v_name


            x_unit_factor, x_unit = self._get_x_unit_factor(
                                         x_unit, measurement[x_name].unit.name
                                         )
            if 'emperatur' in x_name:
                x_mass += x_unit_factor
            else:
                x_mass *= x_unit_factor

            if x_inverse:
                x_mass = 1 / x_mass
                if not logplot:
                    v = np.log(v * unit_factor) * 1 / unit_factor

            ax.plot(
                x_mass,
                v * unit_factor,
                **plot_kwargs_this_mass,
            )

        ax.set_ylabel(f"signal / [{unit}]")

        if x_unit:
            if x_inverse:
                ax.set_xlabel(f"{x_name} / [1/{x_unit}]")
                if not logplot:
                    ax.set_ylabel(f"signal / [ln({unit})]")
            else:
                ax.set_xlabel(f"{x_name} / [{x_unit}]")
        else:
            ax.set_xlabel(f"{x_name}")

        if specs_next_axis:
            self.plot_vs(
                x_name=x_name,
                measurement=measurement,
                ax=specs_next_axis["ax"],
                mass_list=specs_next_axis["mass_list"],
                mol_list=specs_next_axis["mol_list"],
                unit=specs_next_axis["unit"],
                x_unit=x_unit,
                tspan=tspan,
                tspan_bg=specs_next_axis["tspan_bg"],
                logplot=logplot,
                legend=legend,
                x_inverse=x_inverse,
                **plot_kwargs,
            )
            axes = [ax, specs_next_axis["ax"]]
        else:
            axes = None

        if logplot:
            ax.set_yscale("log")
        if legend:
            ax.legend()

        return axes if axes else ax

    def _get_x_unit_factor(
        self,
        x_unit,
        x_unit_name,
    ):
        try:
            if x_unit_name == "celcius" or x_unit_name == "C":
                x_unit_factor = {
                    "C": 1,
                    "celcius": 1,
                    "K": 273.15,
                    "kelvin": 273.15,
                }[x_unit]

            elif x_unit_name == "kelvin" or x_unit_name == "K":
                x_unit_factor = {
                    "K": 1,
                    "kelvin": 1,
                    "celcius": -273.15,
                    "C": -273.15,
                }[x_unit]

            elif x_unit_name == "mbar":
                x_unit_factor = {
                    "mbar": 1,
                    "bar": 1000,
                    "hPa": 1,
                    "kPa": 0.1,
                }[x_unit]

            elif x_unit_name == "bar":
                x_unit_factor = {
                    "mbar": 1e-3,
                    "bar": 1,
                    "hPa": 1e-3,
                    "kPa": 0.1e-3,
                }[x_unit]

            else:
                x_unit_factor = {
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
                    # Temperature conversion
                    "K": 273.15,
                    "kelvin": 273.15,
                }[x_unit]
        except KeyError:
            warnings.warn(
                f"Can't convert original unit '{x_unit_name}' to new unit"
                f"'{x_unit}'. Plotting using original unit!",
                stacklevel=2,
            )
            x_unit_factor = 1
            x_unit = x_unit_name
        return x_unit_factor, x_unit

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

        if not ax:
            ax = (
                axes[0]
                if axes
                else self.new_ax(ylabel=f"signal / [{unit}]", xlabel="time / [s]")
            )
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

        specs_this_axis = {
            "ax": ax,
            "v_list": v_list,
            "unit": unit,
            "unit_factor": unit_factor,
            "tspan_bg": tspan_bg,
        }

        return quantified, specs_this_axis, specs_next_axis


class SpectroMSPlotter(MPLPlotter):
    """A matplotlib plotter specialized in mass spectrometry MID measurements."""

    def __init__(self, measurement=None):
        """Initiate the SpectroMSPlotter with its default Meausurement to plot"""
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
            self.ms_plotter.plot_measurement(
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

        measurement.spectrum_series.heat_plot(
            ax=axes[ms_spec_axes],
            tspan=tspan,
            xspan=xspan,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            max_threshold=max_threshold,
            min_threshold=min_threshold,
            scanning_mask=scanning_mask,
        )

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
            self.ms_plotter.plot_vs(
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
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            max_threshold=max_threshold,
            min_threshold=min_threshold,
            scanning_mask=scanning_mask,
            _sort_indicies=_sort_indicies,
        )

        return axes


#  ----- These are the standard colors for EC-MS plots! ------- #

MIN_SIGNAL = 1e-14  # So that the bottom half of the plot isn't wasted on log(noise)
# TODO: This should probably be customizeable from a settings file.

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
    "N2": "8f8fffff",  # light blue-ish-purple
    "O2": "k",
    "Ar": "c",
    "CO2": "brown",
    "CH4": "r",
    "C2H4": "g",
    "O2@M32": "k",
    "O2@M34": "r",
    "O2@M36": "g",
    "CO2@M44": "brown",
    "CO2@M46": "purple",
    "CO2@M48": "darkslategray",
    # Inset of meta channels #
    "TC temperature": "#808000",
    "RTD temperature": "#808000",
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
