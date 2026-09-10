"""Plotters for spectroelectrochemistry. Makes use of those in spectrum_plotter.py"""

import matplotlib as mpl
from . import ECPlotter, SpectrumSeriesPlotter, SpectroMeasurementPlotter
from ..exceptions import SeriesNotFoundError
from matplotlib import pyplot as plt
import numpy as np


class SECPlotter(SpectroMeasurementPlotter):
    """A spectroelectrochemistry (SEC) matplotlib plotter."""

    def __init__(self, measurement=None):
        """Initiate the plotter with its default Meausurement to plot"""
        super().__init__()
        self.measurement = measurement
        self.ec_plotter = ECPlotter(measurement=measurement)
        self.spectrum_series_plotter = SpectrumSeriesPlotter()

    def plot_measurement(
        self,
        *,
        measurement=None,
        field=None,
        tspan=None,
        xspan=None,
        axes=None,
        cmap_name="inferno",
        make_colorbar=False,
        continuous=None,
        max_threshold=None,
        min_threshold=None,
        **kwargs,
    ):
        """Plot an SECMeasurement in two panels with time as x-asis.

        The top panel is a heat plot with the spectral scanning variable (x) on y-axis
        and color representing the value of the spectral data.
        The bottom panel contains electrochemistry data.

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            field (Field): The field with the spectral data to plot. Defaults to
                `measurement.spectra`
            tspan (timespan): The timespan of data to keep for the measurement.
            xspan (iterable): The span of spectral data to plot
            axes (list of mpl.Axis): The axes to plot on. axes[0] is for the heat
                plot, axes[1] for potential, and axes[2] for current. The axes are
                optional and a new set of axes, where axes[1] and axes[2] are twinned on
                x, are generated if not provided.
            cmap_name (str): The name of the colormap to use. Defaults to "inferno", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html#sequential
            make_colorbar (bool): Whether to make a colorbar.
                FIXME: colorbar at present mis-alignes axes
            kwargs: Additional key-word arguments are passed on to
                ECPlotter.plot_measurement().
            continuous (bool): Optional. Whether to make a continuous heat plot (True) or
                a discrete heat plot for each spectrum (False). In the discrete case,
                each heat plot is a rectangle with the spectrum's duration as its width,
                if available. If the duration is not available, each spectrum heat plot
                extends to the start of the next one.
                Defaults to `measurement.spectrum_series.continuous`.
            max_threshold (float): Maximum value to display.
                Values above are set to max_threshold.
            min_threshold (float): Minimum value to display.
                Values below are set to min_threshold.

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
        axes[0] = self.spectrum_series_plotter.heat_plot(
            spectrum_series=measurement.spectrum_series,
            field=field,
            tspan=tspan,
            xspan=xspan,
            ax=axes[0],
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            continuous=continuous,
            max_threshold=max_threshold,
            min_threshold=min_threshold,
        )
        if make_colorbar:
            pass  # TODO: adjust EC plot to be same width as heat plot despite colorbar.

        axes[1].set_xlim(axes[0].get_xlim())

        return axes

    def plot_vs_potential(
        self,
        *,
        measurement=None,
        field=None,
        tspan=None,
        vspan=None,
        U_name=None,
        J_name=None,
        xspan=None,
        axes=None,
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
            field (Field): The field with the spectral data to plot. Defaults to
                `measurement.spectra`
            tspan (timespan): The timespan of data to keep for the measurement.
            vspan (timespan): The potential span of data to keep for the measurement.
            U_name (str): Optional. The name of the data series to use as potential.
            J_name (str): Optional. The name of the data series to use as current.
            xspan (iterable): The span of spectral data to plot
            axes (list of numpy Axes): The axes to plot on. axes[0] is for the heat
                plot and axes[1] for potential. New are made by default.
            cmap_name (str): The name of the colormap to use. Defaults to "inferno", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html#sequential
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
            U_name=U_name,
            J_name=J_name,
            ax=axes[1],
            **kwargs,
        )

        super().heat_plot_vs(
            measurement=measurement,
            field=field or measurement.spectra,
            vspan=vspan,
            xspan=xspan,
            ax=axes[0],
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            vs=U_name or measurement.U_name,
        )
        axes[1].set_xlim(axes[0].get_xlim())
        return axes

    def plot_stacked_spectra(self, **kwargs):
        kwargs.update(vs="potential")
        return super().plot_stacked_spectra_vs(**kwargs)


class ECOpticalPlotter(SECPlotter):
    def plot_measurement(
        self,
        *,
        measurement=None,
        tspan=None,
        wlspan=None,
        axes=None,
        V_ref=None,
        t_ref=None,
        cmap_name="inferno",
        make_colorbar=False,
        max_threshold=None,
        min_threshold=None,
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
            cmap_name (str): The name of the colormap to use. Defaults to "inferno", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html#sequential
            make_colorbar (bool): Whether to make a colorbar.
                FIXME: colorbar at present misaligns axes
            max_threshold (float): Maximum value to display.
                Values above are set to max_threshold.
            min_threshold (float): Minimum value to display.
                Values below are set to min_threshold.
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
        dOD_series = measurement.calc_dOD(V_ref=V_ref, t_ref=t_ref)

        return super().plot_measurement(
            measurement=measurement,
            tspan=tspan,
            xspan=wlspan,
            axes=axes,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            field=dOD_series,
            max_threshold=max_threshold,
            min_threshold=min_threshold,
            **kwargs,
        )

    def plot_vs_potential(
        self,
        *,
        measurement=None,
        tspan=None,
        vspan=None,
        U_name=None,
        J_name=None,
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
            vspan (timespan): The potential span of data to keep for the measurement.
            U_name (str): Optional. The name of the data series to use as potential.
            J_name (str): Optional. The name of the data series to use as current.
            wlspan (iterable): The wavelength span of spectral data to plot
            axes (list of mpl.Axis): The axes to plot on. axes[0] is for the heat
                plot, axes[1] for potential, and axes[2] for current. The axes are
                optional and a new set of axes, where axes[1] and axes[2] are twinned on
                x, are generated if not provided.
            V_ref (float): Potential to use as reference for calculating optical density
            t_ref (float): Time to use as a reference for calculating optical density
            cmap_name (str): The name of the colormap to use. Defaults to "inferno", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html#sequential
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

        dOD_series = measurement.calc_dOD(V_ref=V_ref, t_ref=t_ref)

        return super().plot_vs_potential(
            measurement=measurement,
            tspan=tspan,
            vspan=vspan,
            U_name=U_name,
            J_name=J_name,
            xspan=wlspan,
            axes=axes,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            field=dOD_series,
            **kwargs,
        )

    def plot_waterfall(
        self,
        *,
        measurement=None,
        ax=None,
        V_ref=None,
        t_ref=None,
        cmap_name="jet",
        make_colorbar=True,
        xlim=None,
        ylim=None,
        tspan=None,
        boxcar=None,
    ):
        """Plot an SECMeasurement as spectra colored based on potential.

        At most one of V_ref and t_ref should be given, and if neither are
        given the measurement's default reference_spectrum is used to calculate the
        optical density.

        This uses :func:`~spectrum_plotter.SpectrumSeriesPlotter.plot_waterfall()`

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            ax (matplotlib Axis): The axes to plot on. A new one is made by default.
            V_ref (float): potential to use as reference for calculating optical density
            t_ref (float): time to use as a reference for calculating optical density
            cmap_name (str): The name of the colormap to use. Defaults to "jet", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html
            make_colorbar (bool): Whether to make a colorbar.
            boxcar = boxcar filter size. By default it is one.
        """

        measurement = measurement or self.measurement
        dOD = measurement.calc_dOD(V_ref=V_ref, t_ref=t_ref)

        return super().plot_waterfall_vs(
            measurement=self.measurement,
            field=dOD,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.U_name,
            xlim=xlim,
            ylim=ylim,
            tspan=tspan,
        )

    def plot_wavelengths(
        self,
        *,
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
            cmap_name (str): The name of the colormap to use. Defaults to "jet", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html
            tspan (timespan): The timespan to plot
            **kwargs: Additional key-word arguments are passed on to
                ECPlotter.plot_measurement
        """
        measurement = measurement or self.measurement
        wavelengths = wavelengths or measurement.tracked_wavelengths

        cmap = plt.get_cmap(cmap_name)
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

        axes[1].set_xlim(axes[0].get_xlim())

        self.ec_plotter.plot_measurement(
            measurement=measurement, axes=[axes[1], axes[3]], tspan=tspan, **kwargs
        )
        return axes

    def plot_wavelengths_vs_cv(
        self,
        *,
        measurement=None,
        wavelengths=None,
        axes=None,
        cmap_name="jet",
        tspan=None,
        cycle=None,
        **kwargs,
    ):
        """Plot the dO.D. for specific wavelength vs potential on the left x axis, with
        current on the right y axis.

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            wavelengths (list of str): The names of the wavelengths to track as strings,
                e.g. "w400" for 400 nm
            axes (list of Ax): The axes to plot on, defaults to new matplotlib axes
            cmap_name (str): The name of the colormap to use. Defaults to "jet", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html
            tspan (timespan): The timespan to plot
            **kwargs: Additional key-word arguments are passed on to
                ECPlotter.plot_measurement
        """

        measurement = measurement or self.measurement
        wavelengths = wavelengths or measurement.tracked_wavelengths

        if not axes:
            fig, ax = plt.subplots()
            ax_current = ax.twinx()
            axes = ax, ax_current
        else:
            ax, ax_current = axes

        cmap = plt.get_cmap(cmap_name)
        norm = mpl.colors.Normalize(vmin=min(measurement.wl), vmax=max(measurement.wl))
        t_v, v = measurement.grab("potential")

        # Plot CV on right-hand axis
        self.ec_plotter.plot_vs_potential(
            measurement=measurement, ax=ax_current, tspan=tspan, **kwargs
        )

        for wl_str in wavelengths:
            x = float(wl_str[1:])
            try:
                t, y = measurement.grab(wl_str, tspan=tspan)
            except SeriesNotFoundError:
                measurement.track_wavelength(x)
                t, y = measurement.grab(wl_str, tspan=tspan)
            v_interp = np.interp(t, t_v, v)
            ax.plot(v_interp, y, color=cmap(norm(x)), label=wl_str)
        ax.legend()
        ax.set_xlabel(measurement.U_name)
        ax.set_ylabel(r"$\Delta$O.D.")

        return axes

    def plot_wavelengths_vs_potential(
        self,
        *,
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
            cmap_name (str): The name of the colormap to use. Defaults to "jet", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html
            tspan (timespan): The timespan to plot
            **kwargs: Additional key-word arguments are passed on to
                ECPlotter.plot_vs_potential
        """
        measurement = measurement or self.measurement
        wavelengths = wavelengths or measurement.tracked_wavelengths

        cmap = plt.get_cmap(cmap_name)
        norm = mpl.colors.Normalize(vmin=min(measurement.wl), vmax=max(measurement.wl))
        t_v, v = measurement.grab(measurement.U_name)

        if not axes:
            axes = self.new_two_panel_axes()
        for wl_str in wavelengths:
            x = float(wl_str[1:])
            try:
                t, y = measurement.grab(wl_str, tspan=tspan)
            except SeriesNotFoundError:
                measurement.track_wavelength(x)
                t, y = measurement.grab(wl_str, tspan=tspan)
            v_interp = np.interp(t, t_v, v)
            axes[0].plot(v_interp, y, color=cmap(norm(x)), label=wl_str)
        axes[0].legend()
        axes[0].set_ylabel(r"$\Delta$O.D.")

        self.ec_plotter.plot_vs_potential(
            measurement=measurement, ax=axes[1], tspan=tspan, **kwargs
        )
        return axes

    def plot_waterfall_cycle(
        self,
        *,
        measurement=None,
        J_name="cycle",
        cycle_number=None,
        V_ref=None,
        t_ref=None,
        index_ref=None,
        cmap_name="jet",
        make_colorbar=True,
        ax=None,
        xlim=None,
        ylim=None,
        tspan=None,
        direction=None,
        N_points=10,
        **kwargs,
    ):
        """Plot one cycle of an SECMeasurement as spectra colored based on potential.

        At most one of V_ref and t_ref should be given, and if neither are
        given the measurement's default reference_spectrum is used to calculate the
        optical density.

        This uses :func:`~spectrum_plotter.SpectrumSeriesPlotter.plot_waterfall()`

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            J_name (str): Name of the cycle series
            cycle_number (int) : the cycle to be plotted
            V_ref (float) : reference spectrum potential
            t_ref (float) : reference spectrum time
            **kwargs: Additional key-word argumetns passed on to plot_waterfall
        """
        measurement = measurement or self.measurement
        # get cycle values
        dOD_one_cycle = measurement.get_dOD_cycle(
            J_name=J_name,
            V_ref=V_ref,
            t_ref=t_ref,
            index_ref=index_ref,
            cycle_number=cycle_number,
            direction=direction,
            N_points=N_points,
        )

        return super().plot_waterfall_vs(
            measurement=measurement,
            field=dOD_one_cycle,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.U_name,
            xlim=xlim,
            ylim=ylim,
            tspan=tspan,
        )

    def plot_dOD_cycle_diff(
        self,
        *,
        measurement=None,
        J_name="cycle",
        cycle_1=None,
        cycle_2=None,
        V_ref1=None,
        t_ref1=None,
        index_ref1=None,
        V_ref2=None,
        t_ref2=None,
        index_ref2=None,
        cmap_name="jet",
        make_colorbar=True,
        ax=None,
        xlim=None,
        ylim=None,
        tspan=None,
        upper_noise_bound=None,
        lower_noise_bound=None,
        direction=None,
        N_points=10,
        **kwargs,
    ):
        """Plot the difference between spectra generated by two cycles of an
        SECMeasurement as spectra colored based on potential. This is to determine if
        there
        is any irreversibility in the system.

        At most one of V_ref and t_ref should be given, and if neither are
        given the measurement's default reference_spectrum is used to calculate the
        optical density.

        This uses :func:`~spectrum_plotter.SpectrumSeriesPlotter.plot_waterfall()`

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            J_name (str): Name of the cycle series
            cycle_1 (int) : the cycle number for the first cycle of interest.
            cycle_2 (int) : the cycle number for the second cycle of interest.
            V_ref (float) : reference spectrum potential
            t_ref (float) : reference spectrum time
            upper_noise_bound (float) : upper noise bound of the system. Found by
                measuring blank FTO substrate.
            lower_noise_bound (float): lower noise bound of the systme. Found by
                measuring blank FTO substrate.
            **kwargs: Additional key-word argumetns passed on to plot_waterfall
        """
        measurement = measurement or self.measurement
        # get difference values
        dOD_diff = measurement.get_dOD_cycle_noise(
            J_name=J_name,
            cycle_1=cycle_1,
            cycle_2=cycle_2,
            V_ref1=V_ref1,
            t_ref1=t_ref1,
            index_ref1=index_ref1,
            V_ref2=V_ref2,
            t_ref2=t_ref2,
            index_ref2=index_ref2,
            direction=direction,
            N_points=N_points,
        )

        axes = super().plot_waterfall_vs(
            measurement=measurement,
            field=dOD_diff,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.U_name,
            xlim=xlim,
            ylim=ylim,
            tspan=tspan,
        )

        # add upper and lower noise bounds
        # Add noise region
        if upper_noise_bound is not None and lower_noise_bound is not None:
            axes.axhspan(
                lower_noise_bound,
                upper_noise_bound,
                facecolor="gainsboro",
                edgecolor="black",
                linestyle="--",
            )

        return axes

    def plot_dOD_difference_spectra(
        self,
        *,
        measurement=None,
        diff_spectra=None,
        J_name="cycle",
        cycle_number=None,
        V_ref=None,
        t_ref=None,
        index_ref=None,
        direction=None,
        normalise=True,
        cmap_name="jet",
        make_colorbar=True,
        ax=None,
        xlim=None,
        ylim=None,
        tspan=None,
        wlmin=None,
        wlmax=None,
        step=None,
        denoise=False,
        denoise_method=None,
        PCA_explained_variance=None,
        sg_window=None,
        sg_poly_order=None,
    ):
        """Plot the finite difference spectra between U and delta U (delta U defined by
        the time resolution of the data collection)

        At most one of V_ref and t_ref should be given, and if neither are
        given the measurement's default reference_spectrum is used to calculate the
        optical density.

        This uses :func:`~spectrum_plotter.SpectrumSeriesPlotter.plot_waterfall()`

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            J_name (str): Name of the cycle series
            cycle_1 (int) : the cycle number for the first cycle of interest.
            cycle_2 (int) : the cycle number for the second cycle of interest.
            V_ref (float) : reference spectrum potential
            t_ref (float) : reference spectrum time
            upper_noise_bound (float) : upper noise bound of the system. Found by
                measuring blank FTO substrate.
            lower_noise_bound (float): lower noise bound of the systme. Found by
                measuring blank FTO substrate.
            **kwargs: Additional key-word argumetns passed on to plot_waterfall
        """
        measurement = measurement or self.measurement
        # get finite difference spectra
        if diff_spectra is not None:
            dOD_finite_diff = diff_spectra
        else:
            dOD_finite_diff = measurement.get_dOD_difference_spectra(
                J_name=J_name,
                cycle_number=cycle_number,
                V_ref=V_ref,
                t_ref=t_ref,
                index_ref=index_ref,
                direction=direction,
                normalise=normalise,
                wlmin=wlmin,
                wlmax=wlmax,
                step=step,
            )

        if denoise:
            dOD_finite_diff = measurement.denoise_spectra(
                spectra_field=dOD_finite_diff,
                denoise_method=denoise_method,
                PCA_explained_variance=PCA_explained_variance,
                sg_window=sg_window,
                sg_poly_order=sg_poly_order,
            )

        return super().plot_waterfall_vs(
            measurement=measurement,
            field=dOD_finite_diff,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.U_name,
            xlim=xlim,
            ylim=ylim,
            tspan=tspan,
        )

    def plot_convergent_spectra(
        self,
        measurement=None,
        converge_output=None,
        J_name="cycle",
        cycle_number=None,
        V_ref=None,
        t_ref=None,
        index_ref=None,
        direction=None,
        normalise=True,
        wlmin=400,
        wlmax=900,
        step=1,
        conv_limit=0.01,
        convergence_metric="correlation",
        smooth_distances=None,
        window_length=4,
        polyorder=3,
        prominence=None,
        threshold=None,
        distance=None,
        height=None,
        min_region_width=2,
        min_region_separation=2,
        converged_spectrum_form="average",
        cmap_name="jet",
        make_colorbar=True,
        ax=None,
        vs=None,
        xlim=None,
        ylim=None,
        tspan=None,
        vspan=None,
    ):
        """Plots convergent dOD spectra.
        If V_ref, t_ref, or index_ref are provided, they specify what to reference dOD
            to. Otherwise, dOD is referenced to the SECMeasurement's reference_spectrum.
            Note: If data is too noisy, it may not meet convergence criterion.

        Args:
            J_name (str): Name of the cycle series. Defaults to "cycle"
            cycle_number (int) : the cycle of interest. Defaults to 1.
            V_ref (float): The potential at which to get the reference spectrum
            t_ref (float): The time at which to get the reference spectrum
            index_ref (int): The index of the reference spectrum
            direction (0, 1): the scan direction. 0 = anodic, 1 = cathodic, nothing is
                full cycle. Full cycle by default.
            normalise (bool): Whether or not to normalise spectra. Defaults to true
            step (int): the step size for difference spectra. 1 by default
            wlmin (float): minimum wavelength to consider. 400nm by default.
            wlmax (float): maximum wavelength to consider. 900 nm by default
            conv_limt (float): the limit for convergence of spectra. Set to 0.01 by
                default.
            convergence_metric (str) : the choice of distance algorithm to pass to pdist.
                Correlation by default.
            smooth_distances (bool) : whether to smooth the distance array (used to find
                similar spectra) with a Savitzky-Golay filter before peak-finding.
            window_length (int) : the window length to pass to the Savitzky-Golay
                filter when smooth_distances is True. 4 by default.
            polyorder (int) : the polynomial order to pass to the Savitzky-Golay
                filter when smooth_distances is True. 3 by default.
            prominence (float): To be passed to find_peaks. How far the peak protrudes
                from the rest
            threshold (float) : To be passed to find_peaks. The required threshold of
                peaks, the vertical distance to its neighboring samples
            distance (float) : To be passed to find_peaks. The required minimum
                horizontal distance in samples between neighbouring peaks. Smaller peaks
                    are
                removed first until the condition is met.
            height (float) : To be passed to find_peaks. The required height of peaks.
            min_region_width (int) : Minimum number of points in a convergent region. 2
                by default.
            min_region_separation (int) : The minimum number of spectra that need to be
                between two convergent regions.
            converged_spectrum_form (str) : The way the converged spectra are handled.
                "average" returns an average spectrum over the convergent region.

        """
        measurement = measurement or self.measurement

        # wl = measurement.axes_series[1].data

        if converge_output is None:
            converge = measurement.get_convergent_spectra(
                J_name=J_name,
                cycle_number=cycle_number,
                V_ref=V_ref,
                t_ref=t_ref,
                index_ref=index_ref,
                direction=direction,
                normalise=normalise,
                wlmin=wlmin,
                wlmax=wlmax,
                step=step,
                conv_limit=conv_limit,
                convergence_metric=convergence_metric,
                smooth_distances=smooth_distances,
                window_length=window_length,
                polyorder=polyorder,
                prominence=prominence,
                threshold=threshold,
                distance=distance,
                height=height,
                min_region_width=min_region_width,
                min_region_separation=min_region_separation,
                converged_spectrum_form=converged_spectrum_form,
            )
        else:
            converge = converge_output

        converged_spectra = converge[0]
        converged_potentials = converge[1]
        if vspan is None:
            full_potentials = measurement.grab(measurement.U_name)[1]
        else:
            full_potentials = vspan
        norm = mpl.colors.Normalize(
            vmin=np.min(full_potentials),
            vmax=np.max(full_potentials),
        )

        return super().plot_waterfall_vs(
            measurement=measurement,
            field=converged_spectra,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.U_name,
            t=converged_potentials,
            xlim=xlim,
            ylim=ylim,
            tspan=tspan,
            norm=norm,
        )

    def plot_fit_and_residuals(
        self,
        *,
        measurement=None,
        cycle_data=None,
        converged_data=None,
        wlmin=400,
        wlmax=900,
        bounds_bool=None,
        noise_bound=0,
        cmap_name="jet",
        make_colorbar=True,
        axes=None,
        xlim=None,
        ylim=None,
        tspan=None,
        direction=None,
        vspan=None,
        xspan=None,
        **kwargs,
    ):
        """Plot the fitting coefficients of convergent spectra for one cycle of an
        SECMeasurement vs potential on the upper figure
        and plot the residuals of the fit on the lower figure.

        This uses :func:`~spectrum_plotter.SpectrumSeriesPlotter.plot_waterfall()`

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            cycle_data: output of get_dOD_cycle
            converged_data: output of get_convergent_spectra
            wlmin (float): minimum wavelength to consider. 400nm by default.
            wlmax (float): maximum wavelength to consider. 900 nm by default
            bounds_bool (bool) : switch to apply a physicality constraint such that the
                components must increase as the spectral intensity is always increasing
                                The strictness of this constraint removes some
                                    flexibility so in real data one might relax it to say
                                    that the next component can only be the noise
                                        threshold
                                    below the previous component.
                                False by default
            noise_bound (float) : the noise bound of the system. 0 by default
            **kwargs: Additional key-word argumetns passed on to plot_waterfall
        """
        from scipy.interpolate import interp1d

        measurement = measurement or self.measurement

        fit_result = measurement.fit_convergent_spectra(
            cycle_data=cycle_data,
            converged_data=converged_data,
            bounds_bool=bounds_bool,
            wlmin=wlmin,
            wlmax=wlmax,
            noise_bound=noise_bound,
        )

        coefficients_field = fit_result[0]

        redox_potentials = converged_data[1]

        residuals_field = fit_result[2]

        if not axes:
            axes = self.new_two_panel_axes(
                n_bottom=1,
                n_top=1,
                emphasis="top",
            )

        coefficients = coefficients_field.data
        coefficients_time = coefficients_field.axes_series[0].data

        times, potentials = measurement.grab(
            "potential", tspan=(coefficients_time.min(), coefficients_time.max())
        )
        norm = mpl.colors.Normalize(
            vmin=np.min(potentials),
            vmax=np.max(potentials),
        )
        cmap = mpl.colormaps[cmap_name]
        end_potentials = (np.min(potentials), np.max(potentials))
        potential_interp = interp1d(
            times,
            potentials,
            kind="linear",
            bounds_error=False,
            fill_value=end_potentials,
        )
        coeff_potentials = potential_interp(coefficients_time)

        axes[0]

        for i in range(coefficients.shape[1]):
            axes[0].plot(
                coeff_potentials,
                coefficients[:, i],
                color=cmap(norm(redox_potentials[i])),
            )

        axes[0].set_ylabel("ΔA (OD)")
        axes[0].set_xlabel(measurement.U_name)

        super().plot_waterfall_vs(
            measurement=measurement,
            field=residuals_field,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=axes[1],
            vs=measurement.U_name,
            xlim=xlim,
            ylim=ylim,
            tspan=tspan,
        )

        # add upper and lower noise bounds
        # Add noise region
        if noise_bound is not None:
            axes[1].axhspan(
                -noise_bound,
                noise_bound,
                facecolor="gainsboro",
                edgecolor="black",
                linestyle="--",
            )

        # axes[1].set_xlim(axes[0].get_xlim())

        return axes

    def plot_fit_reconstruction(
        self,
        *,
        measurement=None,
        cycle_data=None,
        converged_data=None,
        wlmin=400,
        wlmax=900,
        bounds_bool=False,
        noise_bound=0,
        cmap_name="jet",
        make_colorbar=True,
        ax=None,
        xlim=None,
        ylim=None,
        tspan=None,
        direction=None,
        **kwargs,
    ):

        measurement = measurement or self.measurement

        fit_result = measurement.fit_convergent_spectra(
            cycle_data=cycle_data,
            converged_data=converged_data,
            bounds_bool=bounds_bool,
            wlmin=wlmin,
            wlmax=wlmax,
            noise_bound=noise_bound,
        )

        reconstruction = fit_result[1]

        return super().plot_waterfall_vs(
            measurement=self.measurement,
            field=reconstruction,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.U_name,
            xlim=xlim,
            ylim=ylim,
            tspan=tspan,
        )
