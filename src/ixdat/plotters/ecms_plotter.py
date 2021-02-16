from matplotlib import pyplot as plt
from matplotlib import gridspec
from .ec_plotter import ECPlotter


class ECMSPlotter:
    """A matplotlib plotter specialized in electrochemistry measurements."""

    def __init__(self, measurement=None):
        """Initiate the ECMSPlotter with its default Meausurement to plot"""
        self.measurement = measurement

    def plot_measurement(
        self,
        *,
        measurement=None,
        mass_list=None,
        tspan=None,
        V_str=None,
        J_str=None,
        axes=None,
        V_color="k",
        J_color="r",
        logplot=True,
        **kwargs,
    ):
        measurement = measurement or self.measurement

        gs = gridspec.GridSpec(5, 1)
        # gs.update(hspace=0.025)
        if not axes:
            axes = [plt.subplot(gs[0:3, 0])]
            axes += [plt.subplot(gs[3:5, 0])]
            axes += [axes[1].twinx()]
        mass_list = mass_list or measurement.mass_list
        if mass_list:
            for mass in mass_list:
                t, v = measurement.grab(mass, tspan=tspan)
                v[v < MIN_SIGNAL] = MIN_SIGNAL
                axes[0].plot(t, v, color=STANDARD_COLORS.get(mass, "k"), label=mass)
        if measurement.potential:
            ECPlotter.plot_measurement(
                self,
                measurement=measurement,
                axes=[axes[1], axes[2]],
                V_str=V_str,
                J_str=J_str,
                V_color=V_color,
                J_color=J_color,
                **kwargs,
            )
        axes[0].xaxis.set_label_position("top")
        axes[0].tick_params(
            axis="x", top=True, bottom=False, labeltop=True, labelbottom=False
        )
        axes[0].set_xlabel("time / [s]")
        axes[1].set_xlabel("time / [s]")
        axes[0].set_ylabel("signal / [A]")
        axes[1].set_xlim(axes[0].get_xlim())
        if logplot:
            axes[0].set_yscale("log")

    def plot_vs_potential(self):
        pass


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

MIN_SIGNAL = 1e-14
