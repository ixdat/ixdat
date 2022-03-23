"""Plotters for spectroelectrochemistry. Makes use of those in spectrum_plotter.py"""

import matplotlib as mpl

from .base_mpl_plotter import MPLPlotter
from .ec_plotter import ECPlotter
from .spectrum_plotter import SpectrumSeriesPlotter
from ..exceptions import SeriesNotFoundError


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
        t_ref=None,
        cmap_name="inferno",
        make_colorbar=False,
        **kwargs,
    ):
        """Plot an SECMeasurement in two panels with time as x-asis.

        The top panel is a heat plot with wavelength on y-axis and color representing
        spectrum. At most one of V_ref and t_ref should be given, and if neither are
        given the measurement's default reference_spectrum is used to calculate the
        optical density.

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            tspan (timespan): The timespan of data to keep for the measurement.
            wlspan (iterable): The wavelength span of spectral data to plot
            axes (list of mpl.Axis): The axes to plot on. axes[0] is for the heat
                plot, axes[1] for potential, and axes[2] for current. The axes are
                optional and a new set of axes, where axes[1] and axes[2] are twinned on
                x, are generated if not provided.
            V_ref (float): Potential to use as reference for calculating optical density
            t_ref (float): Time to use as a reference for calculating optical density
            cmap_name (str): The name of the colormap to use. Defaults to "inferno",
                which ranges from black through red and orange to yellow-white. "jet"
                is also good.
            make_colorbar (bool): Whether to make a colorbar.
                FIXME: colorbar at present mis-alignes axes
            kwargs: Additional key-word arguments are passed on to
                ECPlotter.plot_measurement().

        Returns:
            list of Axes: axes=[spectra, potential, None, current]
                axes[0] is the top axis with the heat map of the spectra
                axes[1] is the bottom left axis with electrochemical potential
                axes[2] is None (this is where a top right axis would go)
                axes[3] is the bottom right axis with electrode current
        """
        measurement = measurement or self.measurement

        if not axes:
            axes = self.new_two_panel_axes(
                n_bottom=2,
                n_top=1,
                emphasis="top",
            )
        self.ec_plotter.plot_measurement(
            measurement=measurement,
            axes=[axes[1], axes[3]],
            tspan=tspan,
            **kwargs,
        )

        dOD_series = measurement.calc_dOD(V_ref=V_ref, t_ref=t_ref)
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
        self,
        measurement=None,
        ax=None,
        V_ref=None,
        t_ref=None,
        cmap_name="jet",
        make_colorbar=True,
    ):
        """Plot an SECMeasurement as spectra colored based on potential.

        The top panel is a heat plot with wavelength on y-axis and color representing
        spectrum. At most one of V_ref and t_ref should be given, and if neither are
        given the measurement's default reference_spectrum is used to calculate the
        optical density.

        This uses SpectrumSeriesPlotter.plot_waterfall()

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            tspan (timespan): The timespan of data to keep for the measurement.
            wlspan (iterable): The wavelength span of spectral data to plot
            ax (matplotlib Axis): The axes to plot on. A new one is made by default.
            V_ref (float): potential to use as reference for calculating optical density
            t_ref (float): time to use as a reference for calculating optical density
            cmap_name (str): The name of the colormap to use. Defaults to "inferno",
                which ranges from black through red and orange to yellow-white. "jet"
                is also good.
            make_colorbar (bool): Whether to make a colorbar.
        """
        measurement = measurement or self.measurement
        dOD = measurement.calc_dOD(V_ref=V_ref, t_ref=t_ref)

        return self.spectrum_series_plotter.plot_waterfall(
            field=dOD,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.U_name,
        )

    def plot_vs_potential(
        self,
        measurement=None,
        tspan=None,
        vspan=None,
        v_name=None,
        j_name=None,
        axes=None,
        wlspan=None,
        V_ref=None,
        cmap_name="inferno",
        make_colorbar=False,
        **kwargs,
    ):
        """Plot an SECMeasurement in two panels with potential as x-asis.

        The top panel is a heat plot with wavelength on y-axis and color representing
        spectrum. At most one of V_ref and t_ref should be given, and if neither are
        given the measurement's default reference_spectrum is used to calculate the
        optical density.

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            tspan (timespan): The timespan of data to keep for the measurement.
            vspan (timespan): The potential span of data to keep for the measurement.
            v_name (str): Optional. The name of the data series to use as potential.
            j_name (str): Optional. The name of the data series to use as current.
            wlspan (iterable): The wavelength span of spectral data to plot
            axes (list of numpy Axes): The axes to plot on. axes[0] is for the heat
                plot and axes[1] for potential. New are made by default.
            V_ref (float): potential to use as reference for calculating optical density
            t_ref (float): time to use as a reference for calculating optical density
            cmap_name (str): The name of the colormap to use. Defaults to "inferno",
                which ranges from black through red and orange to yellow-white. "jet"
                is also good.
            make_colorbar (bool): Whether to make a colorbar.
            kwargs: Additional key-word arguments are passed on to
                ECPlotter.plot_vs_potential().
        """
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
            U_name=v_name,
            J_name=j_name,
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
            vs=v_name or measurement.U_name,
        )
        axes[1].set_xlim(axes[0].get_xlim())
        return axes

    def plot_wavelengths(
        self,
        measurement=None,
        wavelengths=None,
        axes=None,
        cmap_name="jet",
        tspan=None,
        **kwargs,
    ):
        """Plot the dO.D. for specific wavelength in the top panel and EC in bottom

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            wavelengths (list of str): The names of the wavelengths to track as strings,
                e.g. "w400" for 400 nm
            axes (list of Ax): The axes to plot on, defaults to new matplotlib axes
            cmap_name (str): Name of the colormap. Defaults to "jet"
            tspan (timespan): The timespan to plot
            **kwargs: Additional key-word arguments are passed on to
                ECPlotter.plot_measurement
        """
        measurement = measurement or self.measurement
        wavelengths = wavelengths or measurement.tracked_wavelengths

        cmap = mpl.cm.get_cmap(cmap_name)
        norm = mpl.colors.Normalize(vmin=min(measurement.wl), vmax=max(measurement.wl))

        if not axes:
            axes = self.new_two_panel_axes(n_bottom=2)
        for wl_str in wavelengths:
            x = float(wl_str[1:])
            try:
                t, y = measurement.grab(wl_str, tspan=tspan)
            except SeriesNotFoundError:
                measurement.track_wavelength(x)
                t, y = measurement.grab(wl_str, tspan=tspan)
            axes[0].plot(t, y, color=cmap(norm(x)), label=wl_str)
        axes[0].legend()
        axes[0].set_ylabel(r"$\Delta$O.D.")

        self.ec_plotter.plot_measurement(
            measurement=measurement, axes=[axes[1], axes[3]], tspan=tspan, **kwargs
        )

    def plot_wavelengths_vs_potential(
        self,
        measurement=None,
        wavelengths=None,
        axes=None,
        cmap_name="jet",
        tspan=None,
        **kwargs,
    ):
        """Plot the dO.D. for specific wavelength in the top panel vs potential

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            wavelengths (list of str): The names of the wavelengths to track as strings,
                e.g. "w400" for 400 nm
            axes (list of Ax): The axes to plot on, defaults to new matplotlib axes
            cmap_name (str): Name of the colormap. Defaults to "jet"
            tspan (timespan): The timespan to plot
            **kwargs: Additional key-word arguments are passed on to
                ECPlotter.plot_vs_potential
        """
        measurement = measurement or self.measurement
        wavelengths = wavelengths or measurement.tracked_wavelengths

        cmap = mpl.cm.get_cmap(cmap_name)
        norm = mpl.colors.Normalize(vmin=min(measurement.wl), vmax=max(measurement.wl))

        if not axes:
            axes = self.new_two_panel_axes()
        for wl_str in wavelengths:
            x = float(wl_str[1:])
            try:
                t, y = measurement.grab(wl_str, tspan=tspan)
            except SeriesNotFoundError:
                measurement.track_wavelength(x)
                t, y = measurement.grab(wl_str, tspan=tspan)
            v = measurement.U
            axes[0].plot(v, y, color=cmap(norm(x)), label=wl_str)
        axes[0].legend()
        axes[0].set_ylabel(r"$\Delta$O.D.")

        self.ec_plotter.plot_vs_potential(
            measurement=measurement, ax=axes[1], tspan=tspan, **kwargs
        )
