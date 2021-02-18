import numpy as np
from matplotlib import pyplot as plt
from .plotting_tools import color_axis


class ECPlotter:
    """A matplotlib plotter specialized in electrochemistry measurements."""

    def __init__(self, measurement=None):
        """Initiate the ECPlotter with its default Meausurement to plot"""
        self.measurement = measurement

    def plot_measurement(
        self,
        *,
        measurement=None,
        tspan=None,
        V_str=None,
        J_str=None,
        axes=None,
        V_color="k",
        J_color="r",
        **kwargs,
    ):
        """Plot two variables on two y-axes vs time

        All arguments are optional. By default it plots potential in black on the left
        y-axis and current in red on the right y-axis, using data from its entire
        measurement. The axes are colored to match the traces and labled witht the
        respective series names.

        Args:
            measurement (Measurement): The measurement to plot, if not the one the
                plotter was initiated with.
            tspan (iter of float): The timespan (wrt to measurement.tstamp) to plot.
            V_str (string): The name of the ValueSeries to plot on the left y-axis.
                Defaults to measurement.V_str, which for an ECMeasurement is the name
                of its most calibrated/correct potential.
            J_str (string): The name of the ValueSeries to plot on the right y-axis.
                Defaults to measurement.J_str, which for an ECMeasurement is the name
                of its most normalized/correct current.
            axes (list of matplotlib.Axis): Two axes to plot on, if not the default
                new twinx()'d axes. axes[0] is for V_str and axes[1] for J_str.
            V_color (str): The color to plot the series called V_str. Defaults to black.
            J_color (str): The color to plot the series called J_str. Defaults to red.
            kwargs (dict): Additional key-word arguments are passed to matplotlib's
                plot() function, for both potential and current.

        Returns list of matplotlib.pyplot.Axis: The axes plotted on.
        """
        measurement = measurement or self.measurement
        V_str = V_str or (
            measurement.V_str
            if measurement.RE_vs_RHE is not None
            else measurement.E_str
        )
        J_str = J_str or (
            measurement.J_str if measurement.A_el is not None else measurement.I_str
        )
        t_v, v = measurement.grab(V_str, tspan=tspan, include_endpoints=False)
        t_j, j = measurement.grab(J_str, tspan=tspan, include_endpoints=False)
        if axes:
            ax1, ax2 = axes
        else:
            fig, ax1 = plt.subplots()
            ax2 = ax1.twinx()
            axes = [ax1, ax2]
        ax1.plot(t_v, v, "-", color=V_color, label=V_str, **kwargs)
        ax2.plot(t_j, j, "-", color=J_color, label=J_str, **kwargs)
        ax1.set_xlabel("time / [s]")
        ax1.set_ylabel(V_str)
        ax2.set_ylabel(J_str)
        color_axis(ax1, V_color, lr="left")
        color_axis(ax2, J_color, lr="right")
        return axes

    def plot_vs_potential(
        self, measurement=None, tspan=None, V_str=None, J_str=None, ax=None, **kwargs
    ):
        """Plot an ECMeasurement with electrode potential on the x-axis.

        This can actually plot with anything on the x-axis, by specifying what you want
        on the x-axis using V_str. The y-axis variable, which can be specified by J_str,
        is interpolated onto the time corresponding to the x-axis variable.
            TODO: This is a special case of the not-yet-implemented generalized
                `plot_vs`. Consider an inheritance structure to reduce redundancy in
                future plotters.
        All arguments are optional. By default it will plot current vs potential in
        black on a single axis for the whole experiment.
            TODO: color gradient (cmap=inferno) from first to last cycle.

        Args:
            measurement (Measurement): What to plot. Defaults to the measurement the
                plotter was initialized with
            tspan (iter of float): The timespan, relative to vs measurement.tstamp, on
                which to plot.
            V_str (str): Name of the x-axis ValueSeries. Defaults to calibrated potential
            J_str (str): Name of the y-axis ValueSeries. Defaults to normalized current.
            ax (matplotlib.pyplot.Axis): The axis to plot on, if not a new one.
            kwargs: Additional key-word arguments are passed to matplotlib's
                plot() function, including `color`.

        Returns matplotlib.pyplot.axis: The axis plotted on.
        """
        measurement = measurement or self.measurement
        V_str = V_str or (
            measurement.V_str
            if measurement.RE_vs_RHE is not None
            else measurement.E_str
        )
        J_str = J_str or (
            measurement.J_str if measurement.A_el is not None else measurement.I_str
        )
        t_v, v = measurement.grab(V_str, tspan=tspan, include_endpoints=False)
        t_j, j = measurement.grab(J_str, tspan=tspan, include_endpoints=False)

        j_v = np.interp(t_v, t_j, j)
        if not ax:
            fig, ax = plt.subplots()

        if "color" not in kwargs:
            kwargs["color"] = "k"
        ax.plot(v, j_v, **kwargs)
        ax.set_xlabel(V_str)
        ax.set_ylabel(J_str)
        return ax


class CVDiffPlotter:
    """A matplotlib plotter for highlighting the difference between two cv's."""

    def __init__(self, measurement=None):
        """Initiate the ECPlotter with its default CyclicVoltammagramDiff to plot"""
        self.measurement = measurement

    def plot(self, measurement=None, ax=None):
        measurement = measurement or self.measurement
        ax = ECPlotter.plot_vs_potential(
            self, measurement=measurement.cv_1, axes=ax, color="g"
        )
        ax = ECPlotter.plot_vs_potential(
            self, measurement=measurement.cv_2, ax=ax, color="k", linestyle="--"
        )
        t1, v1 = measurement.cv_1.grab("potential")
        j1 = measurement.cv_1.grab_for_t("current", t=t1)
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
        measurement = measurement or self.measurement
        return ECPlotter.plot_measurement(
            self, measurement=measurement, axes=axes, **kwargs
        )

    def plot_diff(self, measurement=None, tspan=None, ax=None):
        measurement = measurement or self.measurement
        t, v = measurement.grab("potential", tspan=tspan, include_endpoints=False)
        j_diff = measurement.grab_for_t("current", t)
        v_scan = measurement.scan_rate.data
        # a mask which is true when cv_1 had bigger current than cv_2:
        mask = np.logical_xor(0 < j_diff, v_scan < 0)

        if not ax:
            fig, ax = plt.subplots()

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
