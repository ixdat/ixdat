"""Plotter for Electrochemistry"""

import numpy as np
from .base_mpl_plotter import MPLPlotter
from .plotting_tools import color_axis


class ECPlotter(MPLPlotter):
    """A matplotlib plotter specialized in electrochemistry measurements."""

    def __init__(self, measurement=None):
        """Initiate the ECPlotter with its default Meausurement to plot"""
        super().__init__()
        self.measurement = measurement

    def plot_measurement(
        self,
        *,
        measurement=None,
        tspan=None,
        v_name=None,
        j_name=None,
        axes=None,
        v_color="k",
        j_color="r",
        V_str=None,
        J_str=None,
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
            v_name (string): The name of the ValueSeries to plot on the left y-axis.
                Defaults to measurement.V_str, which for an ECMeasurement is the name
                of its most calibrated/correct potential.
            j_name (string): The name of the ValueSeries to plot on the right y-axis.
                Defaults to measurement.J_str, which for an ECMeasurement is the name
                of its most normalized/correct current.
            axes (list of matplotlib.Axis): Two axes to plot on, if not the default
                new twinx()'d axes. axes[0] is for V_str and axes[1] for J_str.
            V_str (str): DEPRECIATED. now v_name
            J_str (str): DEPRECIATED. now j_name
            v_color (str): The color to plot v_name. Defaults to black.
            j_color (str): The color to plot j_name. Defaults to red.
            **plot_kwargs (dict): Additional key-word arguments are passed to
                matplotlib's plot() function. See below for a few examples

        Keyword Arguments:
            linestle (str): Type of line, e.g. "-" for solid or "--" for dotted

        Returns list of matplotlib.pyplot.Axis: The axes plotted on.
        """
        measurement = measurement or self.measurement
        if V_str or J_str:
            print(
                "DEPRECIATION WARNING! V_str has been renamed v_name and J_str has "
                "been renamed j_name. Get it right next time."
            )
        v_name = v_name or V_str or measurement.v_name
        # FIXME: We need a better solution for V_str and J_str that involves the
        #   Calibration and is generalizable. see:
        #   https://github.com/ixdat/ixdat/pull/11#discussion_r679290123
        j_name = j_name or J_str or measurement.j_name
        t_v, v = measurement.grab(v_name, tspan=tspan)
        t_j, j = measurement.grab(j_name, tspan=tspan)
        if axes:
            ax1, ax2 = axes
        else:
            ax1 = self.new_ax()
            ax2 = ax1.twinx()
            axes = [ax1, ax2]
        ax1.plot(t_v, v, "-", color=v_color, label=v_name, **plot_kwargs)
        ax2.plot(t_j, j, "-", color=j_color, label=j_name, **plot_kwargs)
        ax1.set_xlabel("time / [s]")
        ax1.set_ylabel(v_name)
        ax2.set_ylabel(j_name)
        color_axis(ax1, v_color, lr="left")
        color_axis(ax2, j_color, lr="right")
        return axes

    def plot_vs_potential(
        self,
        measurement=None,
        tspan=None,
        v_name=None,
        j_name=None,
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
            v_name (str): Name of the x-axis variable. Defaults to calibrated potential
            j_name (str): Name of the y-axis variable. Defaults to normalized current.
            ax (matplotlib.pyplot.Axis): The axis to plot on, if not a new one.
            **plot_kwargs (dict): Additional key-word arguments are passed to
                matplotlib's plot() function. See below for a few examples

        Keyword Arguments:
            color (color): Color of the trace, e.g. "r", "blue", or RGB like [0, 0, 1]
            linestle (str): Type of line, e.g. "-" for solid or "--" for dotted

        Returns matplotlib.pyplot.axis: The axis plotted on.
        """

        measurement = measurement or self.measurement
        v_name = v_name or (
            measurement.v_name
            if measurement.RE_vs_RHE is not None
            else measurement.E_name
        )
        j_name = j_name or (
            measurement.j_name if measurement.A_el is not None else measurement.I_name
        )
        t_v, v = measurement.grab(v_name, tspan=tspan)
        t_j, j = measurement.grab(j_name, tspan=tspan)

        j_v = np.interp(t_v, t_j, j)
        if not ax:
            ax = self.new_ax()

        if "color" not in plot_kwargs:
            plot_kwargs["color"] = "k"
        ax.plot(v, j_v, **plot_kwargs)
        ax.set_xlabel(v_name)
        ax.set_ylabel(j_name)
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
            self, measurement=measurement.cv_compare_1, axes=ax, color="g"
        )
        ax = ECPlotter.plot_vs_potential(
            self, measurement=measurement.cv_compare_2, ax=ax, color="k", linestyle="--"
        )
        t1, v1 = measurement.cv_compare_1.grab("potential")
        j1 = measurement.cv_compare_1.grab_for_t("current", t=t1)
        j_diff = measurement.grab_for_t("current", t=t1)
        # a mask which is true when cv_1 had bigger current than cv_2:
        v_scan = measurement.scan_rate.data
        mask = np.logical_xor(0 < j_diff, v_scan < 0)

        ax.fill_between(v1, j1 - j_diff, j1, where=mask, alpha=0.2, color="g")
        ax.fill_between(
            v1,
            j1 - j_diff,
            j1,
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
        t, v = measurement.grab("potential", tspan=tspan, include_endpoints=False)
        j_diff = measurement.grab_for_t("current", t)
        v_scan = measurement.scan_rate.data
        # a mask which is true when cv_1 had bigger current than cv_2:
        mask = np.logical_xor(0 < j_diff, v_scan < 0)

        if not ax:
            ax = self.new_ax()

        ax.plot(v[mask], j_diff[mask], "k-", label="cv1 > cv2")
        ax.plot(
            v[np.logical_not(mask)],
            j_diff[np.logical_not(mask)],
            "k--",
            label="cv1 < cv2",
        )
        return ax

    def plot_vs_potential(self):
        """FIXME: This is needed to satisfy ECMeasurement.__init__"""
        pass
