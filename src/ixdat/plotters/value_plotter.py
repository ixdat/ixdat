"""Classes for plotting measurement data"""

from .base_mpl_plotter import MPLPlotter
from ..exceptions import SeriesNotFoundError


class ValuePlotter(MPLPlotter):
    """Default plotter. By default plots all of the VSeries vs time on a single axis"""

    def __init__(self, measurement=None):
        super().__init__()
        self.measurement = measurement

    def plot(self, *args, **kwargs):
        """Plot the exporter's measurement via plotter.plot_measurement()"""
        return self.plot_measurement(measurement=self.measurement, *args, **kwargs)

    def plot_measurement(
        self,
        measurement=None,
        v_list=None,
        tspan=None,
        ax=None,
        legend=True,
        logscale=False,
    ):
        """Plot a measurement's values vs time

        Args:
            measurement (Measurement): The measurement to plot
            v_list (list of str): The names of the data series to include. Defaults to
                names of all VSeries in the measurement.
            tspan (timespan): The timespan to include in the file, defaults to all of it
            legend (bool): Whether to include a legend. Defaults to True.
            logscale (bool): Whether to use a log-scaled y-axis. Defaults to False.
        """
        measurement = measurement or self.measurement
        if not ax:
            ax = self.new_ax()
        v_list = v_list or measurement.value_names

        for v_name in v_list:
            try:
                v, t = measurement.grab(v_name, tspan=tspan)
            except SeriesNotFoundError as e:
                print(f"WARNING!!! {e}")
                continue
            ax.plot(v, t, label=v_name)

        if legend:
            ax.legend()
        if logscale:
            ax.set_yscale("log")

        return ax
