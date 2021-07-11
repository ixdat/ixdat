import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from .base_mpl_plotter import MPLPlotter
from .ec_plotter import ECPlotter
from .spectrum_plotter import SpectrumSeriesPlotter


class SECPlotter(MPLPlotter):
    """An spectroelectrochemsitry (SEC) matplotlib plotter.

    FIXME: This should make use of the code in spectrum_plotter.SpectrumSeriesPlotter
    """

    def __init__(self, measurement=None):
        """Initiate the plotter with its default Meausurement to plot"""
        self.measurement = measurement
        self.ec_plotter = ECPlotter(measurement=measurement)
        self.spectrum_series_plotter = SpectrumSeriesPlotter(
            spectrum_series=self.measurement
            # FIXME: Maybe SpectrumSeries should inherit from Measurement?
        )

    def plot_measurement(
        self,
        measurement=None,
        tspan=None,
        wlspan=None,
        axes=None,
        V_ref=None,
        cmap_name="inferno",
        make_colorbar=False,
        **kwargs,
    ):
        measurement = measurement or self.measurement

        if not axes:
            axes = self.new_two_panel_axes(
                n_bottom=2,
                n_top=1,
                emphasis="top",
            )
        self.ec_plotter.plot_measurement(
            measurement=measurement,
            axes=[axes[1], axes[2]],
            tspan=tspan,
            **kwargs,
        )

        dOD_series = measurement.calc_dOD(V_ref=V_ref)
        axes[0] = self.spectrum_series_plotter.heat_plot(
            field=dOD_series,
            tspan=tspan,
            xspan=wlspan,
            ax=axes[0],
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
        )
        if make_colorbar:
            pass  # TODO: adjust EC plot to be same width as heat plot despite colorbar.

        axes[1].set_xlim(axes[0].get_xlim())

        return axes

    def plot_waterfall(
        self, measurement=None, cmap_name="jet", make_colorbar=True, V_ref=None, ax=None
    ):
        measurement = measurement or self.measurement
        dOD = measurement.calc_dOD(V_ref=V_ref)

        return self.spectrum_series_plotter.plot_waterfall(
            field=dOD,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.V_str,
        )

    def plot_vs_potential(
        self,
        measurement=None,
        tspan=None,
        vspan=None,
        V_str=None,
        J_str=None,
        axes=None,
        wlspan=None,
        V_ref=None,
        cmap_name="inferno",
        make_colorbar=False,
        **kwargs,
    ):
        measurement = measurement or self.measurement

        if not axes:
            axes = self.new_two_panel_axes(
                n_bottom=1,
                n_top=1,
                emphasis="top",
            )

        self.ec_plotter.plot_vs_potential(
            measurement=measurement,
            tspan=tspan,
            V_str=V_str,
            J_str=J_str,
            ax=axes[1],
            **kwargs,
        )

        dOD_series = measurement.calc_dOD(V_ref=V_ref)
        axes[0] = self.spectrum_series_plotter.heat_plot_vs(
            field=dOD_series,
            vspan=vspan,
            xspan=wlspan,
            ax=axes[0],
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            vs=V_str or measurement.V_str,
        )
        axes[1].set_xlim(axes[0].get_xlim())
        return axes
