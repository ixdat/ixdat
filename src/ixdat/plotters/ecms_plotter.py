from matplotlib import pyplot as plt
from matplotlib import gridspec
from .ec_plotter import ECPlotter
from .ms_plotter import MSPlotter


class ECMSPlotter:
    """A matplotlib plotter for EC-MS measurements."""

    def __init__(self, measurement=None):
        """Initiate the ECMSPlotter with its default Meausurement to plot"""
        self.measurement = measurement

    def plot_measurement(
        self,
        *,
        measurement=None,
        axes=None,
        mass_list=None,
        tspan=None,
        V_str=None,
        J_str=None,
        V_color="k",
        J_color="r",
        logplot=True,
        legend=True,
        **kwargs,
    ):
        """Make an EC-MS plot vs time and return the axis handles.

        Allocates tasks to ECPlotter.plot_measurement() and MSPlotter.plot_measurement()

        TODO: add all functionality in the legendary plot_experiment() in EC_MS.Plotting
            - variable subplot sizing (emphasizing EC or MS)
            - plotting of calibrated data (mol_list instead of mass_list)
            - units!
            - optionally two y-axes in the upper panel
        Args:
            measurement (ECMSMeasurement): defaults to the measurement to which the
                plotter is bound (self.measurement)
            axes (list of three matplotlib axes): axes[0] plots the MID data,
                axes[1] the variable given by V_str (potential), and axes[2] the
                variable given by J_str (current). By default three axes are made with
                axes[0] a top panel with 3/5 the area, and axes[1] and axes[2] are
                the left and right y-axes of the lower panel with 2/5 the area.
            mass_list (list of str): The names of the m/z values, eg. ["M2", ...] to
                plot. Defaults to all of them (measurement.mass_list)
            tspan (iter of float): The time interval to plot, wrt measurement.tstamp
            V_str (str): The name of the value to plot on the lower left y-axis.
                Defaults to the name of the series `measurement.potential`
            J_str (str): The name of the value to plot on the lower right y-axis.
                Defaults to the name of the series `measurement.current`
            V_color (str): The color to plot the variable given by 'V_str'
            J_color (str): The color to plot the variable given by 'J_str'
            logplot (bool): Whether to plot the MS data on a log scale (default True)
            legend (bool): Whether to use a legend for the MS data (default True)
            kwargs (dict): Additional kwargs go to all calls of matplotlib's plot()
        """
        measurement = measurement or self.measurement

        if not axes:
            gs = gridspec.GridSpec(5, 1)
            # gs.update(hspace=0.025)
            axes = [plt.subplot(gs[0:3, 0])]
            axes += [plt.subplot(gs[3:5, 0])]
            axes += [axes[1].twinx()]
        if hasattr(measurement, "potential") and measurement.potential:
            # then we have EC data!
            ECPlotter.plot_measurement(
                self,
                measurement=measurement,
                axes=[axes[1], axes[2]],
                tspan=tspan,
                V_str=V_str,
                J_str=J_str,
                V_color=V_color,
                J_color=J_color,
                **kwargs,
            )
        if mass_list or hasattr(measurement, "mass_list"):
            # then we have MS data!
            MSPlotter.plot_measurement(
                self,
                ax=axes[0],
                tspan=tspan,
                mass_list=mass_list,
                logplot=logplot,
                legend=legend,
            )
        axes[0].xaxis.set_label_position("top")
        axes[0].tick_params(
            axis="x", top=True, bottom=False, labeltop=True, labelbottom=False
        )
        axes[1].set_xlim(axes[0].get_xlim())
        return axes

    def plot_vs_potential(self):
        """FIXME: This is needed due to assignment in ECMeasurement.__init__"""
        pass
