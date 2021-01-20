"""Classes for plotting measurement data"""

from matplotlib import pyplot as plt


class ValuePlotter:
    """Default plotter. By default plots all of the VSeries vs time on a single axis"""

    def __init__(self, measurement=None):
        self.measurement = measurement

    def plot(self, *args, **kwargs):
        """Plot the exporter's measurement via plotter.plot_measurement()"""
        return self.plot_measurement(measurement=self.measurement, *args, **kwargs)

    def plot_measurement(
        self, measurement, v_list=None, tspan=None, ax=None, legend=True, logscale=False
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
        if not ax:
            fig, ax = plt.subplots()
        v_list = v_list or measurement.value_names

        for v_name in v_list:
            v, t = measurement.get_t_and_v(v_name, tspan=tspan)
            ax.plot(v, t, label=v_name)

        if legend:
            ax.legend()

        return ax
