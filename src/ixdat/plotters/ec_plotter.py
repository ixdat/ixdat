"""Plotter for Electrochemistry"""

import warnings
import numpy as np
from .base_mpl_plotter import MPLPlotter
from .plotting_tools import color_axis
from ..tools import deprecate
from ..exceptions import SeriesNotFoundError


class ECPlotter(MPLPlotter):
    """A matplotlib plotter specialized in electrochemistry measurements."""

    def __init__(self, measurement=None):
        """Initiate the ECPlotter with its default Meausurement to plot"""
        super().__init__()
        self.measurement = measurement

    @deprecate("0.1", "Use `U_name` instead.", "0.3", kwarg_name="V_str")
    @deprecate("0.1", "Use `J_name` instead.", "0.3", kwarg_name="J_str")
    @deprecate("0.1", "Use `U_color` instead.", "0.3", kwarg_name="V_color")
    def plot_measurement(
        self,
        *,
        measurement=None,
        tspan=None,
        U_name=None,
        J_name=None,
        U_color=None,
        J_color=None,
        V_str=None,
        J_str=None,
        V_color=None,
        axes=None,
        **plot_kwargs,
    ):
        """Plot two variables on two y-axes vs time

        All arguments are optional. By default it plots potential in black on the left
        y-axis and current in red on the right y-axis, using data from its entire
        measurement. The axes are colored to match the traces and labeled with the
        respective series names.

        Args:
            measurement (Measurement): The measurement to plot, if not the one the
                plotter was initiated with.
            tspan (iter of float): The timespan (wrt to measurement.tstamp) to plot.
            axes (list of matplotlib.Axis): Two axes to plot on, if not the default
                new twinx()'d axes. axes[0] is for `U_name` and axes[1] for `J_name`.
            U_name (string): The name of the ValueSeries to plot on the left y-axis.
                Defaults to measurement.V_str, which for an ECMeasurement is the name
                of its most calibrated/correct potential.
            J_name (string): The name of the ValueSeries to plot on the right y-axis.
                Defaults to measurement.J_str, which for an ECMeasurement is the name
                of its most normalized/correct current.
            U_color (str): The color to plot U_name. Defaults to black.
            J_color (str): The color to plot J_name. Defaults to red.
            V_str (str): DEPRECATED. Use `U_name`.
            J_str (str): DEPRECATED. Use `J_name`.
            V_color (str): DEPRECATED. Use `U_color`.
            **plot_kwargs (dict): Additional key-word arguments are passed to
                matplotlib's plot() function. See below for a few examples

        Keyword Arguments:
            linestle (str): Type of line, e.g. "-" for solid or "--" for dotted

        Returns list of matplotlib.pyplot.Axis: The axes plotted on.
        """
        measurement = measurement or self.measurement

        # apply deprecated arguments (the user will get a warning):
        U_name = U_name or V_str
        J_name = J_name or J_str
        U_color = U_color or V_color

        # apply defaults.
        U_name = U_name or measurement.U_name
        J_name = J_name or measurement.J_name
        U_color = U_color or "k"
        J_color = J_color or "r"

        if axes:
            ax1, ax2 = axes
        else:
            ax1 = self.new_ax()
            ax2 = ax1.twinx()
            axes = [ax1, ax2]
        ax1.set_xlabel("time / [s]")
        ax1.set_ylabel(U_name)
        ax2.set_ylabel(J_name)
        color_axis(ax1, U_color, lr="left")
        color_axis(ax2, J_color, lr="right")

        try:
            t_v, v = measurement.grab(U_name, tspan=tspan)
        except SeriesNotFoundError:
            warnings.warn(f"No '{U_name}' found in {measurement}")
        else:
            ax1.plot(t_v, v, "-", color=U_color, label=U_name, **plot_kwargs)

        try:
            t_j, j = measurement.grab(J_name, tspan=tspan)
        except SeriesNotFoundError:
            warnings.warn(f"No '{J_name}' found in {measurement}")
        else:
            ax2.plot(t_j, j, "-", color=J_color, label=J_name, **plot_kwargs)

        return axes

    def plot_vs_potential(
        self,
        measurement=None,
        tspan=None,
        U_name=None,
        J_name=None,
        ax=None,
        **plot_kwargs,
    ):
        """Plot an ECMeasurement with electrode potential on the x-axis.

        This can actually plot with anything on the x-axis, by specifying what you want
        on the x-axis using V_str. The y-axis variable, which can be specified by J_str,
        is interpolated onto the time corresponding to the x-axis variable.
        .. TODO::
            This is a special case of the not-yet-implemented generalized
            `plot_vs`. Consider an inheritance structure to reduce redundancy in
            future plotters.
            sub-TODO: hide or fix TODO's using sphix boxes.
        All arguments are optional. By default it will plot current vs potential in
        black on a single axis for the whole experiment.
            TODO: color gradient (cmap=inferno) from first to last cycle.

        Args:
            measurement (Measurement): What to plot. Defaults to the measurement the
                plotter was initialized with
            tspan (iter of float): The timespan, relative to vs measurement.tstamp, on
                which to plot.
            U_name (str): Name of the x-axis variable. Defaults to calibrated potential
            J_name (str): Name of the y-axis variable. Defaults to normalized current.
            ax (matplotlib.pyplot.Axis): The axis to plot on, if not a new one.
            **plot_kwargs (dict): Additional key-word arguments are passed to
                matplotlib's plot() function. See below for a few examples

        Keyword Arguments:
            color (color): Color of the trace, e.g. "r", "blue", or RGB like [0, 0, 1]
            linestle (str): Type of line, e.g. "-" for solid or "--" for dotted

        Returns matplotlib.pyplot.axis: The axis plotted on.
        """

        measurement = measurement or self.measurement
        U_name = U_name or measurement.U_name
        J_name = J_name or measurement.J_name
        t_v, v = measurement.grab(U_name, tspan=tspan)
        t_j, j = measurement.grab(J_name, tspan=tspan)

        j_v = np.interp(t_v, t_j, j)
        if not ax:
            ax = self.new_ax()

        if "color" not in plot_kwargs:
            plot_kwargs["color"] = "k"
        ax.plot(v, j_v, **plot_kwargs)
        ax.set_xlabel(U_name)
        ax.set_ylabel(J_name)
        return ax


class CVDiffPlotter(MPLPlotter):
    """A matplotlib plotter for highlighting the difference between two cv's."""

    def __init__(self, measurement=None):
        """Initiate the ECPlotter with its default CyclicVoltammagramDiff to plot"""
        super().__init__()
        self.measurement = measurement

    def plot(self, measurement=None, ax=None):
        """Plot the two cycles of the CVDiff measurement and fill in the areas between

        example: https://ixdat.readthedocs.io/en/latest/_images/cv_diff.svg
        """
        measurement = measurement or self.measurement
        # FIXME: This is probably the wrong use of plotter functions.
        #    see https://github.com/ixdat/ixdat/pull/30/files#r810926968
        ax = ECPlotter.plot_vs_potential(
            self, measurement=measurement.cv_compare_1, ax=ax, color="g"
        )
        ax = ECPlotter.plot_vs_potential(
            self, measurement=measurement.cv_compare_2, ax=ax, color="k", linestyle="--"
        )
        t1, U1 = measurement.cv_compare_1.grab("potential")
        J1 = measurement.cv_compare_1.grab_for_t("current", t=t1)
        J_diff = measurement.grab_for_t("current", t=t1)
        # a mask which is true when cv_1 had bigger current than cv_2:
        v_scan = measurement.grab_for_t("scan_rate", t=t1)
        mask = np.logical_xor(0 < J_diff, v_scan < 0)

        ax.fill_between(U1, J1 - J_diff, J1, where=mask, alpha=0.2, color="g")
        ax.fill_between(
            U1,
            J1 - J_diff,
            J1,
            where=np.logical_not(mask),
            alpha=0.1,
            hatch="//",
            color="g",
        )

        return ax

    def plot_measurement(self, measurement=None, axes=None, **kwargs):
        """Plot the difference between the two cv's vs time"""
        measurement = measurement or self.measurement
        # FIXME: not correct useage of
        return ECPlotter.plot_measurement(
            self, measurement=measurement, axes=axes, **kwargs
        )

    def plot_diff(self, measurement=None, tspan=None, ax=None):
        """Plot the difference between the two cv's vs potential.

        The trace is solid where the current in cv_2 is greater than cv_1 in the anodic
        scan or the current cv_2 is more negative than cv_1 in the cathodic scan.
        """
        measurement = measurement or self.measurement
        t, U = measurement.grab("potential", tspan=tspan, include_endpoints=False)
        j_diff = measurement.grab_for_t("current", t)
        v_scan = measurement.grab_for_t("scan_rate", t)
        # a mask which is true when cv_1 had bigger current than cv_2:
        mask = np.logical_xor(0 < j_diff, v_scan < 0)

        if not ax:
            ax = self.new_ax()

        ax.plot(U[mask], j_diff[mask], "k-", label="cv1 > cv2")
        ax.plot(
            U[np.logical_not(mask)],
            j_diff[np.logical_not(mask)],
            "k--",
            label="cv1 < cv2",
        )
        return ax

    def plot_vs_potential(self):
        """FIXME: This is needed to satisfy ECMeasurement.__init__"""
        pass
