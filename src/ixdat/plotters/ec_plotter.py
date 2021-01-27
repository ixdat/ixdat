import numpy as np
from matplotlib import pyplot as plt
from .plotting_tools import color_axis


class ECPlotter:
    def __init__(self, measurement=None):
        self.measurement = measurement

    def plot_measurement(
        self,
        measurement=None,
        tspan=None,
        V_str=None,
        J_str=None,
        axes=None,
        V_color="k",
        J_color="r",
        **kwargs,
    ):
        measurement = measurement or self.measurement
        V_str = V_str or measurement.V_str
        J_str = J_str or measurement.J_str
        t_v, v = measurement.get_t_and_v(V_str, tspan=tspan)
        t_j, j = measurement.get_t_and_v(J_str, tspan=tspan)
        if axes:
            ax1, ax2 = axes
        else:
            fig, ax1 = plt.subplots()
            ax2 = ax1.twinx()
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
        measurement = measurement or self.measurement
        V_str = V_str or measurement.V_str
        J_str = J_str or measurement.J_str
        t_v, v = measurement.get_t_and_v(V_str, tspan=tspan)
        t_j, j = measurement.get_t_and_v(J_str, tspan=tspan)

        j_v = np.interp(t_v, t_j, j)
        if not ax:
            fig, ax = plt.subplots()

        if "color" not in kwargs:
            kwargs["color"] = "k"
        ax.plot(v, j_v, **kwargs)
        ax.set_xlabel(V_str)
        ax.set_ylabel(J_str)
        return ax
