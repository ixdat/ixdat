from .base_mpl_plotter import MPLPlotter


class MSPlotter(MPLPlotter):
    """A matplotlib plotter specialized in mass spectrometry MID measurements."""

    def __init__(self, measurement=None):
        """Initiate the ECMSPlotter with its default Meausurement to plot"""
        self.measurement = measurement

    def plot_measurement(
        self,
        measurement=None,
        ax=None,
        mass_list=None,
        tspan=None,
        logplot=True,
        legend=True,
    ):
        """Plot m/z signal vs time (MID) data and return the axis handle.

        Args:
            measurement (MSMeasurement): defaults to the one that initiated the plotter
            ax (matplotlib axis): Defaults to a new axis
            mass_list (list of str): The names of the m/z values, eg. ["M2", ...] to
                plot. Defaults to all of them (measurement.mass_list)
            tspan (iter of float): The timespan, wrt measurement.tstamp, on which to plot
            logplot (bool): Whether to plot the MS data on a log scale (default True)
            legend (bool): Whether to use a legend for the MS data (default True)
        """
        if not ax:
            ax = self.new_ax(ylabel="signal / [A]", xlabel="time / [s]")
        measurement = measurement or self.measurement
        mass_list = mass_list or measurement.mass_list
        for mass in mass_list:
            t, v = measurement.grab(mass, tspan=tspan, include_endpoints=False)
            v[v < MIN_SIGNAL] = MIN_SIGNAL
            ax.plot(t, v, color=STANDARD_COLORS.get(mass, "k"), label=mass)
        if logplot:
            ax.set_yscale("log")
        if legend:
            ax.legend()


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
