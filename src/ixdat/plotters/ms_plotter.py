"""Plotter for Mass Spectrometry"""
import warnings
from ..data_series import Field
import numpy as np
from .base_mpl_plotter import MPLPlotter, ConversionError
from ..units import ureg, DimensionalityError, Quantity


class MSPlotter(MPLPlotter):
    """A matplotlib plotter specialized in mass spectrometry MID measurements."""

    def __init__(self, measurement=None):
        """Initiate the ECMSPlotter with its default Meausurement to plot"""
        super().__init__()
        self.measurement = measurement
        self.ureg = ureg

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
        logdata=False,
        legend=True,
        use_quantity=True,
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
            unit (str): unit of the y axis. defaults to "A" or "mol/s"
            x_unit (str): unit of the x axis variable (usually time). defaults to "s"
            logplot (bool): Whether to plot the MS data on a log scale (default True)
            logdata (bool): Whether to plot the natural logarithm of MS data on a
                linear scale (default False)
            legend (bool): Whether to use a legend for the MS data (default True)
            use_quantity (bool): boolean if plotter should plot data according to the
                units of the dataseries or just the magnitude of the dataseries
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
            use_quantity,
        )
        ax = specs_this_axis["ax"]

        v_list = specs_this_axis["v_list"]
        tspan_bg = specs_this_axis["tspan_bg"]
        unit = specs_this_axis["unit"]

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
                    remove_background=remove_background,
                    include_endpoints=False,
                    return_quantity=use_quantity,
                )
            else:
                t, v = measurement.grab_signal(
                    v_or_v_name,
                    tspan=tspan,
                    tspan_bg=tspan_bg,
                    remove_background=remove_background,
                    include_endpoints=False,
                    return_quantity=use_quantity,
                )

            if ureg("mol / s / cm ** 2").check(unit):
                if not hasattr(measurement, "A_el"):
                    warnings.warn(
                        "Measurement does not have an attribute A_el and cannot "
                        "calibrate to A_el"
                    )
                    continue
                if isinstance(measurement.A_el, type(Quantity)):
                    v = v / measurement.A_el
                else:
                    v = v / (measurement.A_el * ureg.cm**2)
                    warnings.warn(
                        f"You should explicit set the unit of A_el: {measurement.A_el} "
                        f"before plotting.Falling back to default unit {ureg.cm**2:~P}"
                    )

            if logplot:
                try:
                    v.m[v.m < MIN_SIGNAL] = MIN_SIGNAL
                except (DimensionalityError, AttributeError):
                    v[v < MIN_SIGNAL] = MIN_SIGNAL

            if logdata:
                logplot = False
                vu = v.to(unit).u if use_quantity else unit
                v = (
                    np.log(v.to(unit).m) * ureg.dimensionless
                    if use_quantity
                    else np.log(v / (1 * unit).to_base_units().m) * ureg.dimensionless
                )
                ylabel = f"{vu:~ixdat_log}"

            if use_quantity and not v.check(unit) and not logdata:
                vn = f"n_dot_{v_name}" if quantified else f"{v_name}"
                error_message = (
                    f"cannot plot {vn} with units {v.u:~P} on axis with units {unit}."
                    f"If the unit of '{vn}' is wrong you can set the correct unit with "
                    f"measurement['{vn}'].unit.set_unit('correct_unit_name')"
                    "Otherwise choose a compatible unit to plot on axis. List of "
                    f"compatible units to dataseries: {ureg.get_compatible_units(v.u)}"
                )
                print(error_message)
                continue

            ax.plot(
                t if use_quantity else t * ureg.s,  # ixdat internal time is seconds
                v
                if use_quantity
                else (
                    v * ureg.mol / ureg.s if quantified else v * ureg.A
                ),  # ixdat internal default is A and mol/s for un/calibrated signals
                color=color,
                label=v_name,
                **kwargs,
            )

        if logdata or (use_quantity and v.check(ureg.dimensionless)):
            ax.set_ylabel(ylabel) if logdata else ax.set_ylabel("signal / [a.u.]")
            ax.yaxis.isDefault_label = True

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
                logdata=logdata,
                legend=legend,
                use_quantity=use_quantity,
                **kwargs,
            )
            axes = [ax, specs_next_axis["ax"]]
        else:
            axes = None

        if logplot:
            ax.set_yscale("log")
        if legend:
            ax.legend()

        if x_unit:
            x_unit = (
                x_unit
                if isinstance(x_unit, type(ureg.s))
                else (ureg(x_unit).u if isinstance(x_unit, str) else x_unit.u)
            )
            ax.xaxis.set_units(x_unit)
            

        if unit and not logdata:
            unit = (
                unit
                if isinstance(unit, type(ureg.A))
                else (ureg(unit).u if isinstance(unit, str) else unit.u)
            )
            ax.yaxis.set_units(unit)
                    

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
        logdata=False,
        legend=True,
        use_quantity=True,
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
            logdata (bool): Whether to plot the natural logarithm of MS data on a
                linear scale (default False)
            legend (bool): Whether to use a legend for the MS data (default True)
            use_quantity (bool): boolean if plotter should plot data according to the
                units of the dataseries or just the magnitude of the dataseries
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
            use_quantity,
        )
        ax = specs_this_axis["ax"]
        v_list = specs_this_axis["v_list"]
        tspan_bg = specs_this_axis["tspan_bg"]
        unit = specs_this_axis["unit"]
        # unit_factor = specs_this_axis["unit_factor"]

        t, x = measurement.grab(
            x_name, tspan=tspan, include_endpoints=True, return_quantity=use_quantity
        )

        if use_quantity and x_unit and not x.check(x_unit):
            error_message = (
                f"cannot plot {x_name} with units {x.u:~P} on axis with units {x_unit}."
                f"If the unit of '{x_name}' is wrong you can set the correct unit with "
                f"measurement['{x_name}'].unit.set_unit('correct_unit_name')"
                "Otherwise choose a compatible unit to plot on axis. List of "
                f"compatible units to dataseries: {ureg.get_compatible_units(x.u)}"
            )

            return print(error_message)

        for i, v_name in enumerate(v_list):
            if quantified:
                t_v, v = measurement.grab_flux(
                    v_name,
                    tspan=tspan,
                    tspan_bg=tspan_bg,
                    remove_background=remove_background,
                    include_endpoints=False,
                    return_quantity=use_quantity,
                )
            else:
                t_v, v = measurement.grab_signal(
                    v_name,
                    tspan=tspan,
                    tspan_bg=tspan_bg,
                    remove_background=remove_background,
                    include_endpoints=False,
                    return_quantity=use_quantity,
                )

            if ureg("mol / s / cm ** 2").check(unit):
                if not hasattr(measurement, "A_el"):
                    warnings.warn(
                        "Measurement does not have an attribute A_el and cannot "
                        "calibrate to A_el"
                    )
                    continue
                if isinstance(measurement.A_el, type(Quantity)):
                    v = v / measurement.A_el
                else:
                    v = v / (measurement.A_el * ureg.cm**2)
                    warnings.warn(
                        f"You should explicit set the unit of A_el: {measurement.A_el} "
                        f"before plotting. Falling back to default unit {ureg.cm**2:~P}"
                    )

            if use_quantity and not v.check(unit) and not logdata:
                vn = f"n_dot_{v_name}" if quantified else f"{v_name}"
                error_message = (
                    f"cannot plot {vn} with units {v.u:~P} on axis with units {unit}."
                    f"If the unit of '{vn}' is wrong you can set the correct unit with "
                    f"measurement['{vn}'].unit.set_unit('correct_unit_name')"
                    "Otherwise choose a compatible unit to plot on axis. List of "
                    f"compatible units to dataseries: {ureg.get_compatible_units(v.u)}"
                )
                print(error_message)
                continue

            if logplot:
                try:
                    v.m[v.m < MIN_SIGNAL] = MIN_SIGNAL
                except (DimensionalityError, AttributeError):
                    v[v < MIN_SIGNAL] = MIN_SIGNAL

            if logdata:
                logplot = False
                vu = v.to(unit).u if use_quantity else unit
                v = (
                    np.log(v.to(unit).m) * ureg.dimensionless
                    if use_quantity
                    else np.log(v / (1 * unit).to_base_units().m) * ureg.dimensionless
                )
                ylabel = f"{vu:~ixdat_log}"

            x_mass = np.interp(t_v, t, x)
            plot_kwargs_this_mass = plot_kwargs.copy()
            if "color" not in plot_kwargs:
                plot_kwargs_this_mass["color"] = STANDARD_COLORS.get(v_name, "k")
            if "label" not in plot_kwargs:
                plot_kwargs_this_mass["label"] = v_name

            ax.plot(
                x_mass if use_quantity else x_mass * ureg.dimensionless,
                v
                if use_quantity
                else (
                    v * ureg.mol / ureg.s if quantified else v * ureg.A
                ),  # ixdat internal default is A and mol/s for un/calibrated signals
                **plot_kwargs_this_mass,
            )

        if logdata or (use_quantity and v.check(ureg.dimensionless)):
            ax.set_ylabel(ylabel) if logdata else ax.set_ylabel("signal / [a.u.]")
            ax.yaxis.isDefault_label = True

        if ax.xaxis.get_units() == ureg.dimensionless:
            ax.set_xlabel(f"{x_name}")
            ax.xaxis.isDefault_label = True

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
                logdata=logdata,
                **plot_kwargs,
            )
            axes = [ax, specs_next_axis["ax"]]
        else:
            axes = None

        if logplot:
            ax.set_yscale("log")
        if legend:
            ax.legend()

        if x_unit:
            x_unit = (
                x_unit
                if isinstance(x_unit, type(ureg.s))
                else (ureg(x_unit).u if isinstance(x_unit, str) else x_unit.u)
            )
            ax.xaxis.set_units(x_unit)

        if unit and not logdata:
            unit = (
                unit
                if isinstance(unit, type(ureg.A))
                else (ureg(unit).u if isinstance(unit, str) else unit.u)
            )
            ax.yaxis.set_units(unit)

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
        use_quantity,
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
                else self.new_ax()  # ylabel=f"signal / [{DEFAULT_UNIT['signal']}]",
                # xlabel=f"time / [{DEFAULT_UNIT['time']}]",
                #             use_quantity=use_quantity)
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
            if isinstance(unit, (list, tuple)):
                unit_right = unit[1]
                unit = unit[0]
            else:
                unit_right = unit

            specs_next_axis = {
                "ax": ax_right,
                "unit": unit_right,
                "mass_list": None if quantified else v_list_right,
                "mol_list": v_list_right if quantified else None,
                "tspan_bg": tspan_bg_right,
            }
        else:
            specs_next_axis = None

        if isinstance(unit, (ureg.Unit, ureg.Quantity)):
            unit = unit.u if hasattr(unit, "units") else unit
            
        elif isinstance(unit, str):
            unit = ureg(unit).u
        
        else:
            unit = ureg.mol / ureg.s if quantified else ureg.A

        specs_this_axis = {
            "ax": ax,
            "v_list": v_list,
            "unit": unit,
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
        x_unit=None,
        logplot=True,
        logdata=False,
        legend=True,
        xspan=None,
        cmap_name="inferno",
        make_colorbar=False,
        emphasis="top",
        ms_data="top",
        max_threshold=None,
        min_threshold=None,
        scanning_mask=None,
        vmin=None,
        vmax=None,
        use_quantity=True,
        **kwargs,
    ):
        """Plot m/z signal, mass spectra vs time (MID) data and return the axes of a two
        panel figure.

        There are four ways to specify what to plot. Only specify one of these::
            mass_list: Uncalibrated signals in [(u/n/p)A] on on axis
            mass_lists: Uncalibrated signals in [(u/n/p)A] on two axes
            mol_list: Calibrated signals in [(u/n/p)mol/s] on on axis
            mol_lists: Calibrated signals in [(u/n/p)mol/s] on two axes

        Two axes refers to separate left and right y-axes. Default is to use all
        available masses as mass_list.

        Args:
            measurement (SpectroMSMeasurement): Defaults to the one that initiated the
                plotter
            axes (list of matplotlib axis): Left and right y-axes if mass_lists are given
                default to axes[0] as left and axes[2] as right axis for plotting masses
                and axes[1] for plotting MSSpectra
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
            logdata (bool): Whether to plot the natural logarithm of MS data on a
                linear scale (default False)
            legend (bool): Whether to use a legend for the MS data (default True)
            xspan (iter of float): The physical span for spectra to plot
            cmap_name (str): Colour map to pass to heat_plot method
            make_colorbar (bool): Include a colour bar. Misalignes time axis with other
                panels in same figure
            emphasis (str): Whether to emphasise top or bottom panel 3/5 fig size or eq.
            ms_data (str): Whether to plot ms_data on "top" or "bottom" panel
            max_threshold (int or float): Only applies to spectra plotted with heat_plot.
                All values above max threshold is set to 0 (zero).
            min_threshold (int or float): Only applies to spectra plotted with heat_plot.
                All values below threshold is set to 0 (zero).
            scanning_mask (boolean list): Only applies to spectra plotted with heat_plot.
                List of booleans of same shape as SpectrumSeries.data to exclude specific
                data prior to plotting
            vmin (int or float): Shift minimum value in color bar. Default lowest value
                in measurement.spectrum_series.
            vmax (int or float): Shift maximum value in color bar. Default highest value
                in measurement.spectrum_series.
            use_quantity (bool): boolean if plotter should plot data according to the
                units of the dataseries or just the magnitude of the dataseries
            kwargs: extra key-word args are passed on to matplotlib's plot()
        """

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
                x_unit=x_unit,
                logplot=logplot,
                logdata=logdata,
                legend=legend,
                use_quantity=use_quantity,
                **kwargs,
            )

        # then we have the spectrum series to plot
        measurement.spectrum_series.heat_plot(
            ax=axes[ms_spec_axes],
            tspan=tspan,
            xspan=xspan,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            max_threshold=max_threshold,
            min_threshold=min_threshold,
            scanning_mask=scanning_mask,
            vmin=vmin,
            vmax=vmax,
        )

        return axes

    def plot_measurement_vs(
        self,
        *,
        vs_name,
        measurement=None,
        axes=None,
        mass_list=None,
        mass_lists=None,
        mol_list=None,
        mol_lists=None,
        vspan=None,
        tspan=None,
        tspan_bg=None,
        remove_background=None,
        unit=None,
        vs_unit=None,
        logplot=True,
        logdata=False,
        legend=True,
        xspan=None,
        cmap_name="inferno",
        make_colorbar=False,
        vmin=None,
        vmax=None,
        emphasis="top",
        ms_data="top",
        max_threshold=None,
        min_threshold=None,
        scanning_mask=None,
        sort_spectra="linear",
        use_quantity=True,
        **kwargs,
    ):
        """Plot m/z signal and MSSpectra data in a two panel subfigure vs a specified
        variable and return the axes.

        There are four ways to specify which (MID) signals to plot in panel.
        Only specify one of these:
            mass_list: Uncalibrated signals in [(u/n/p)A] on on axis
            mass_lists: Uncalibrated signals in [(u/n/p)A] on two axes
            mol_list: Calibrated signals in [(u/n/p)mol/s] on on axis
            mol_lists: Calibrated signals in [(u/n/p)mol/s] on two axes

        Two axes refers to separate left and right y-axes. Default is to use all
        available masses as mass_list.


        Args:
            measurement (SpectroMSMeasurement): Defaults to the one that initiated the
                plotter
            vs_name (str): Name of the series to plot versus.
            axes (list of matplotlib axis): Defaults to axes[0], axes[2] for left and
                right axis for plotting masses and axes[1] for MSSpectra data.
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
            vspan (iter of float): The value interval to plot on x_axis
            tspan (iter of float): The time interval to plot, wrt measurement.tstamp
            tspan_bg (timespan): A timespan for which to assume the signal is at its
                background. The average signals during this timespan are subtracted.
                If `mass_lists` are given rather than a single `mass_list`, `tspan_bg`
                must also be two timespans - one for each axis. Default is `None` for no
                background subtraction.
            remove_background (bool): Whether otherwise to subtract pre-determined
                background signals if available. Defaults to (not logplot)
            unit (str): defaults to "A" or "mol/s"
            vs_unit (str): defaults to v_name.unit.name
            logplot (bool): Whether to plot the MS data on a log scale (default True)
            logdata (bool): Whether to plot the natural logarithm of MS data on a
                linear scale (default False)
            legend (bool): Whether to use a legend for the MS data (default True)
            xspan (iter of float): The physical span for spectra to plot
            cmap_name (str): Colour map to pass to heat_plot method
            make_colorbar (bool): Include a colour bar. Misalignes time axis with other
                panels in same figure
            vmin (int or float): Shift minimum value in color bar. Default lowest value
                in measurement.spectrum_series.
            vmax (int or float): Shift maximum value in color bar. Default highest value
                in measurement.spectrum_series.
            emphasis (str): Whether to emphasise top or bottom panel 3/5 fig size or eq.
                Default 'top'
            ms_data (str): Whether to plot ms_data on "top" or "bottom" panel.
                Default 'top'
            max_threshold (int or float): Only applies to spectra plotted with heat_plot.
                All values above max threshold is set to 0 (zero).
            min_threshold (int or float): Only applies to spectra plotted with heat_plot.
                All values below threshold is set to 0 (zero).
            scanning_mask (boolean list): Only applies to spectra plotted with heat_plot.
                List of booleans of same shape as SpectrumSeries.data to exclude specific
                data prior to plotting
            sort_spectra (list or str): Whether or not to sort the spectra data prior to
                plotting.
                There is three specifers:
                    'none':
                        This gives no new sorting. Effectively the spectras are sorted
                        by time. (tstamp for each of the spectrum).
                    'linear' (default):
                        the spectras are sorted linear to v_name from low to high.
                    a list of same shape as field.data to be sorted:
                        This list is passed directly as the indices to sort the spectras.
                        Defaults to sort lowest to highest value. Example
                Note: If tspan spans a time span of a measurement with up and down
                cycles in v_name, this might yield funny looking heat_plots.

                    Example:
                    Scanning up and down in temperature from T_low to T_high the spectras
                    obtained wil be plotted from [T_low_start ..  T_high .. T_low_end].

                    - If 'none' sorting is specified leads to heat_plot xaxis linearly
                    from T_low_start to T_low_end missing representation of the high
                    values in the middle of the axis.

                    - If 'linear' sorting is specified all spectras obtained
                        are sorted linearly from lowest v_name_value to highest v_name.
                        When data is assymetric from scanning up or down in v_name this
                        leads to abrupt looking figures since two non similair spectras
                        are obtained at similair v_name value and hence plotted next to
                        eachother.
            use_quantity (bool): boolean if plotter should plot data according to the
                units of the dataseries or just the magnitude of the dataseries
            kwargs: extra key-word args are passed on to matplotlib's plot()
        """

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
            # define where one or two axes is plotted for MS data!
            ms_plot_axes = (
                [axes[ms_axes], axes[2]]
                if (mass_lists or mol_lists)
                else [axes[ms_axes]]
            )

            # then we have MS data!
            self.ms_plotter.plot_vs(
                x_name=vs_name,
                measurement=measurement,
                axes=ms_plot_axes,
                tspan=tspan,
                tspan_bg=tspan_bg,
                remove_background=remove_background,
                mass_list=mass_list,
                mass_lists=mass_lists,
                mol_list=mol_list,
                mol_lists=mol_lists,
                unit=unit,
                x_unit=vs_unit,
                logplot=logplot,
                logdata=logdata,
                legend=legend,
                use_quantity=use_quantity,
                **kwargs,
            )

        # To plot heat plot.
        # First get all values for v_name at all spectrum times
        field = measurement.spectrum_series.field
        _data = field.data.copy()

        _t = field.axes_series[0].t
        _v = measurement.grab_for_t(item=vs_name, t=_t)

        if tspan:
            # create t_mask from tspan
            t_mask = np.logical_and(tspan[0] < _t, _t < tspan[-1])
            # apply t_mask to field.data and vs_name.data
            _data = _data[t_mask]
            _v = _v[t_mask]

        if isinstance(sort_spectra, str):
            if sort_spectra == "linear":
                sorted_indicies = np.argsort(_v)

            elif sort_spectra == "none":
                sorted_indicies = np.array([])

            else:
                warnings.warn(
                    f"Recieved {sort_spectra} for sort_spectra."
                    "sort_spectra has to be 'linear',  'none' or "
                    "of type list with same length as spectrum_series.",
                    stacklevel=2,
                )

        elif isinstance(sort_spectra, list):
            if len(_t) == len(sort_spectra):
                sorted_indicies = sort_spectra
            else:
                warnings.warn(
                    f"length [{len(sort_spectra)}] of 'sort_spectra' has to be equal"
                    f"to [{len(_t)}]."
                    "sort_spectra can be 'linear', 'none' or "
                    "of type list with same length as spectrum_series.",
                    stacklevel=2,
                )
        else:
            warnings.warn(
                "sort_spectra has to be 'linear', 'none' or "
                "of type list with same length as spectrum_series."
                "Recived {sort_spectra}.",
                stacklevel=2,
            )
        # sort field data and x_axis equal
        new_field_data = _data[sorted_indicies, :] if len(sorted_indicies) > 0 else _data
        new_x_axis = _v[sorted_indicies] if len(sorted_indicies) > 0 else _v

        new_field = Field(
            name=field.name + f"_sorted_vs_{vs_name}_for_{tspan}",
            unit_name=field.unit_name,
            axes_series=field.axes_series,
            data=new_field_data,
        )
        # Now we can plot the heat plot
        measurement.spectrum_series.heat_plot(
            ax=axes[ms_spec_axes],
            t=new_x_axis,
            field=new_field,
            t_name=vs_name,
            tspan=vspan,
            xspan=xspan,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            vmin=vmin,
            vmax=vmax,
            max_threshold=max_threshold,
            min_threshold=min_threshold,
            scanning_mask=scanning_mask,
        )

        axes[ms_spec_axes].set_xlim(axes[ms_axes].get_xlim())

        if vspan:
            axes[ms_axes].set_xlim([vspan[0], vspan[-1]])
            axes[ms_spec_axes].set_xlim([vspan[0], vspan[-1]])

        return axes


#  ----- These are the standard colors for EC-MS plots! ------- #

MIN_SIGNAL = 1e-14  # So that the bottom half of the plot isn't wasted on log(noise)
# TODO: This should probably be customizeable from a settings file.

DEFAULT_UNIT = {"time": "s", "signal": "A"}

DEFAULT_YLABELS = {
    "signal": {"[current]": 1},
    "cal. sig": {"[substance]": 1, "[time]": -1},
    "time": {"[time]": 1},
}

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
    # FIXME: Upgrade to include user defined colours from file or other module
    # https://github.com/ixdat/ixdat/pull/101/files#r1088739480
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
