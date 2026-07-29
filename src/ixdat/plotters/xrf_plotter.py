"""Plotter for time-resolved X-ray fluoresence"""

from .base_mpl_plotter import MPLPlotter
from .ec_plotter import ECPlotter
from .plot_spec import combine_plot_specs, ec_measurement_spec, time_series_spec
from .renderers import _is_matplotlib_backend, get_renderer


class TRXRFPlotter(MPLPlotter):
    """A time-resolved X-ray fluoresence matplotlib plotter."""

    def __init__(self, measurement=None):
        """Initiate the TRXRFPlotter with its default Meausurement to plot"""
        super().__init__()
        self.measurement = measurement

    def plot_measurement(
        self,
        measurement=None,
        ax=None,
        tspan=None,
        y_name="FF_over_I0",
        backend=None,
        figure=None,
        **kwargs,
    ):
        """Plot FF_over_I0 signal vs time (MID) data and return the axis.

        Args:
            measurement (MSMeasurement): Defaults to the one that initiated the plotter
            ax (matplotlib axis): Defaults to a new axis
            tspan (iter of float): The time interval to plot, wrt measurement.tstamp
            y_name （list of str): The names of siganl in X-ray data, eg. "FF", "I0", "It"
                to plot. Default to FF_over_I0, corresponding to a newly build data from
                original lsit of data: FF/I0.
            backend (str): ``"matplotlib"`` or ``"plotly"``.
            figure: Plotly figure to add the plot to.
            kwargs: extra key-word args are passed on to matplotlib's plot()
        """
        measurement = measurement or self.measurement
        if not _is_matplotlib_backend(backend):
            plot_spec = time_series_spec(
                measurement,
                left_names=[y_name],
                tspan=tspan,
                left_label=y_name,
                line_style=kwargs.get("linestyle", kwargs.get("ls")),
                line_width=kwargs.get("linewidth", kwargs.get("lw")),
                show_legend=False,
            )
            return get_renderer(backend).render(plot_spec, figure=figure)

        t, y = measurement.grab(y_name, tspan=tspan)

        time_name = measurement["FF_over_I0"].tseries.name

        if not ax:
            ax = self.new_ax(xlabel=time_name, ylabel=y_name)

        ax.plot(t, y, **kwargs)

        return ax


class ECTRXRFPlotter(MPLPlotter):
    """A matplotlib plotter for EC-TRXRF measurements."""

    def __init__(self, measurement=None):
        """Initiate the ECTRXRFPlotter with its default Meausurement to plot"""
        super().__init__()
        self.measurement = measurement
        self.ec_plotter = ECPlotter(measurement=self.measurement)

    def plot_measurement(
        self,
        measurement=None,
        axes=None,
        tspan=None,
        y_name="FF_over_I0",
        U_name=None,
        J_name=None,
        U_color=None,
        J_color=None,
        backend=None,
        figure=None,
        **kwargs,
    ):
        """Make an EC-TRXRF plot vs time and return the axis handles.

        Args:
            measurement (MSMeasurement): Defaults to the one that initiated the plotter
            ax (matplotlib axis): Defaults to a new axis
            tspan (iter of float): The time interval to plot, wrt measurement.tstamp
            y_name （str): The names of siganl in X-ray data, eg. "FF", "I0", "It"...
                to plot. Default to FF_over_I0, corresponding to a newly build data from
                original lsit of data: FF/I0.

            U_name (str): The name of the value to plot on the lower left y-axis.
                Defaults to the name of the series `measurement.potential`
            J_name (str): The name of the value to plot on the lower right y-axis.
                Defaults to the name of the series `measurement.current`
            U_color (str): The color to plot the variable given by 'V_str'
            J_color (str): The color to plot the variable given by 'J_str'
            backend (str): ``"matplotlib"`` or ``"plotly"``.
            figure: Plotly figure to add the plot to.
            kwargs (dict): Additional kwargs go to all calls of matplotlib's plot()

        Returns:
            list of Axes: (top_left, bottom_left, top_right, bottom_right) where:
                axes[0] is top_left is MS data;
                axes[1] is bottom_left is potential;
                axes[2] is top_right is additional TRXRF data if left and right
                    X-ray siganl were plotted (otherwise axes[2] is None); and
                axes[3] is bottom_right is current.
        """

        measurement = measurement or self.measurement
        if not _is_matplotlib_backend(backend):
            xrf_spec = time_series_spec(
                measurement,
                left_names=[y_name],
                tspan=tspan,
                left_label=y_name,
                line_style=kwargs.get("linestyle", kwargs.get("ls")),
                line_width=kwargs.get("linewidth", kwargs.get("lw")),
                show_legend=False,
            )
            ec_spec = ec_measurement_spec(
                measurement,
                tspan=tspan,
                U_name=U_name,
                J_name=J_name,
                U_color=U_color,
                J_color=J_color,
                line_style=kwargs.get("linestyle", kwargs.get("ls")),
                line_width=kwargs.get("linewidth", kwargs.get("lw")),
            )
            return get_renderer(backend).render(
                combine_plot_specs(xrf_spec, ec_spec),
                figure=figure,
            )

        if not axes:
            axes = self.new_two_panel_axes(n_bottom=2, n_top=1, emphasis=None)

        t, y = measurement.grab(y_name, tspan=tspan)
        time_name = measurement["FF_over_I0"].tseries.name

        axes[0].set_ylabel(y_name)
        axes[0].set_xlabel(time_name)
        axes[0].plot(t, y, **kwargs)

        self.ec_plotter.plot_measurement(
            measurement=measurement,
            axes=[axes[1], axes[3]],
            tspan=tspan,
            U_name=U_name,
            J_name=J_name,
            U_color=U_color,
            J_color=J_color,
            **kwargs,
        )
        axes[1].set_xlim(axes[0].get_xlim())

        return axes

    plot_vs_potential = ECPlotter.plot_vs_potential
