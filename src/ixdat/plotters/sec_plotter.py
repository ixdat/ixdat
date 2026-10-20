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
        zspan=None,
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
            zspan (iterable): [min, max] of the values to show in the heat plot.
                Values outside are shown as min or max.

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
            min_threshold=zspan[0] if zspan is not None else None,
            max_threshold=zspan[1] if zspan is not None else None,
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
        ax = super().plot_stacked_spectra_vs(**kwargs)

        return ax


class ECFTIRPlotter(SECPlotter):
    def plot_measurement(
        self,
        *,
        measurement=None,
        tspan=None,
        wavenumber_span=None,
        axes=None,
        index_ref=None,
        V_ref=None,
        t_ref=None,
        cmap_name="inferno",
        make_colorbar=False,
        zspan=None,
        **kwargs,
    ):
        """Plot an ECFTIRMeasurement in two panels with time as x-asis.

        The top panel is a heat plot with wavenumber on y-axis and color representing
        spectrum. At most one of index_ref, V_ref and t_ref should be given, and if
        neither are given the measurement's default reference_spectrum is used for the
        relative spectra

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            tspan (timespan): The timespan of data to keep for the measurement.
            wlspan (iterable): The wavelength span of spectral data to plot
            axes (list of mpl.Axis): The axes to plot on. axes[0] is for the heat
                plot, axes[1] for potential, and axes[2] for current. The axes are
                optional and a new set of axes, where axes[1] and axes[2] are twinned on
                x, are generated if not provided.
            index_ref (int): Spectrum index to use as reference for calculating
                relative spectra
            V_ref (float): Potential to use as reference for calculating relative spectra
            t_ref (float): Time to use as a reference for calculating relative spectra
            cmap_name (str): The name of the colormap to use. Defaults to "inferno", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html#sequential
            make_colorbar (bool): Whether to make a colorbar.
                FIXME: colorbar at present misaligns axes
            zspan (iterable): [min, max] of the values to show in the heat plot.
                Values outside are shown as min or max.
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
        relative_spectra = measurement.apply_reference(
            V_ref=V_ref, t_ref=t_ref, index_ref=index_ref
        )

        axes = super().plot_measurement(
            measurement=measurement,
            tspan=tspan,
            xspan=wavenumber_span,
            axes=axes,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            field=relative_spectra,
            zspan=zspan,
            **kwargs,
        )

        # axes[0].invert_xaxis()

        return axes

    def plot_waterfall(
        self,
        *,
        measurement=None,
        ax=None,
        V_ref=None,
        t_ref=None,
        index_ref=None,
        cmap_name="jet",
        make_colorbar=True,
        tspan=None,
    ):
        """Plot an ECFTIRMeasurement as spectra colored based on potential.

        At most one of index_ref, V_ref and t_ref should be given, and if none are
        given the measurement's default reference_spectrum is used to calculate the
        relative spectra.

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
        """

        measurement = measurement or self.measurement
        relative_spectra = measurement.apply_reference(
            V_ref=V_ref, t_ref=t_ref, index_ref=index_ref
        )

        ax = super().plot_waterfall_vs(
            measurement=measurement,
            field=relative_spectra,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.U_name,
            tspan=tspan,
        )

        ax.invert_xaxis

        return ax

    def plot_stacked_spectra(
        self, measurement=None, V_ref=None, t_ref=None, index_ref=None, vs=None, **kwargs
    ):
        kwargs.update(vs=vs)
        measurement = measurement or self.measurement
        relative_spectra = measurement.apply_reference(
            V_ref=V_ref, t_ref=t_ref, index_ref=index_ref
        )
        ax = super().plot_stacked_spectra_vs(field=relative_spectra, **kwargs)

        ax.invert_xaxis()

        return ax


class ECOpticalPlotter(SECPlotter):
    def plot_measurement(
        self,
        *,
        measurement=None,
        tspan=None,
        wlspan=None,
        axes=None,
        index_ref=None,
        V_ref=None,
        t_ref=None,
        cmap_name="inferno",
        make_colorbar=False,
        zspan=None,
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
            zspan (iterable): [min, max] of the values to show in the heat plot.
                Values outside are shown as min or max.
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
            zspan=zspan,
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
        index_ref=None,
        cmap_name="jet",
        make_colorbar=True,
        tspan=None,
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
        """

        measurement = measurement or self.measurement
        dOD = measurement.calc_dOD(V_ref=V_ref, t_ref=t_ref, index_ref=index_ref)

        return super().plot_waterfall_vs(
            measurement=measurement,
            field=dOD,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.U_name,
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

    def plot_wavelengths_vs_potential(
        self,
        *,
        measurement=None,
        wavelengths=None,
        axes=None,
        cmap_name="jet",
        tspan=None,
        coplot=False,
        **kwargs,
    ):
        """Plot the dO.D. for specific wavelengths vs potential, with the CV

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            wavelengths (list of str): The names of the wavelengths to track as strings,
                e.g. "w400" for 400 nm
            axes (list of Ax): The axes to plot on, defaults to new matplotlib axes.
                axes[0] is for the wavelengths and axes[1] for the CV.
            cmap_name (str): The name of the colormap to use. Defaults to "jet", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html
            tspan (timespan): The timespan to plot
            coplot (bool): If True, the new axes are one plot with the CV on a right
                y-axis. By default the wavelengths are in a top panel and the CV in a
                bottom panel.
            **kwargs: Additional key-word arguments are passed on to
                ECPlotter.plot_vs_potential
        """
        measurement = measurement or self.measurement
        wavelengths = wavelengths or measurement.tracked_wavelengths

        cmap = plt.get_cmap(cmap_name)
        norm = mpl.colors.Normalize(vmin=min(measurement.wl), vmax=max(measurement.wl))
        t_v, v = measurement.grab(measurement.U_name)

        if not axes:
            if coplot:
                fig, ax = plt.subplots()
                axes = [ax, ax.twinx()]
            else:
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
        if coplot:
            axes[0].set_xlabel(measurement.U_name)

        self.ec_plotter.plot_vs_potential(
            measurement=measurement, ax=axes[1], tspan=tspan, **kwargs
        )
        return axes

    def plot_dOD_cycle_diff(
        self,
        *,
        measurement=None,
        field=None,
        cycle_1=1,
        cycle_2=2,
        c_name="cycle",
        cmap_name="jet",
        make_colorbar=True,
        ax=None,
        tspan=None,
        upper_noise_bound=None,
        lower_noise_bound=None,
    ):
        """Plot the change in dOD from one cycle to another, colored by potential

        This shows irreversible changes in the system. See
        :func:`~calculators.sec_calculators.calc_cycle_difference`.

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            field (Field): The cycle difference to plot. Defaults to
                `calc_cycle_difference(measurement, cycle_1, cycle_2, c_name)`.
            cycle_1 (int): The first cycle. Defaults to 1.
            cycle_2 (int): The second cycle. Defaults to 2.
            c_name (str): The name of the cycle series. Defaults to "cycle".
            cmap_name (str): The name of the colormap to use. Defaults to "jet".
            make_colorbar (bool): Whether to make a colorbar.
            ax (matplotlib Axis): The axes to plot on. A new one is made by default.
            tspan (iterable): The timespan of the spectra to plot
            upper_noise_bound (float): Upper noise bound of the system in [mOD],
                e.g. found by measuring a blank FTO substrate.
            lower_noise_bound (float): Lower noise bound of the system in [mOD]
        """
        from ..calculators.sec_calculators import calc_cycle_difference

        measurement = measurement or self.measurement
        if field is None:
            field = calc_cycle_difference(
                measurement, cycle_1=cycle_1, cycle_2=cycle_2, c_name=c_name
            )
        ax = super().plot_waterfall_vs(
            measurement=measurement,
            field=field,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.U_name,
            tspan=tspan,
        )
        if upper_noise_bound is not None and lower_noise_bound is not None:
            ax.axhspan(
                lower_noise_bound,
                upper_noise_bound,
                facecolor="gainsboro",
                edgecolor="black",
                linestyle="--",
            )
        return ax

    def plot_dOD_difference_spectra(
        self,
        *,
        measurement=None,
        field=None,
        cmap_name="jet",
        make_colorbar=True,
        ax=None,
        tspan=None,
    ):
        """Plot difference spectra colored by potential

        The difference spectra are calculated by the measurement's `dOD_calculator`
        (see :class:`~calculators.sec_calculators.DeltaODCalculator`). Set its step
        first, e.g. `measurement.dOD_calculator.Vstep = 0.005`.

        Args:
            measurement (Measurement): The measurement to be plotted, if different from
                self.measurement
            field (Field): The difference spectra to plot. Defaults to
                `measurement["dOD>diff"]`. Give e.g. `measurement["dOD>diff>smooth"]`
                to plot denoised difference spectra.
            cmap_name (str): The name of the colormap to use. Defaults to "jet".
            make_colorbar (bool): Whether to make a colorbar.
            ax (matplotlib Axis): The axes to plot on. A new one is made by default.
            tspan (iterable): The timespan of the spectra to plot
        """
        measurement = measurement or self.measurement
        field = field if field is not None else measurement["dOD>diff"]
        return super().plot_waterfall_vs(
            measurement=measurement,
            field=field,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.U_name,
            tspan=tspan,
        )

    @staticmethod
    def _get_deconvolutor(measurement):
        """Return the SpectrumDeconvolutor added to the measurement"""
        try:
            return measurement.calculators["spectrum_deconvolutor"]
        except KeyError:
            raise ValueError(
                f"{measurement!r} has no SpectrumDeconvolutor. Make one with"
                " `SpectrumDeconvolutor.from_spectromeasurement(...)` and add it with"
                " `measurement.add_calculator(...)`."
            )

    def plot_converged_spectra(
        self,
        *,
        measurement=None,
        cmap_name="jet",
        make_colorbar=True,
        ax=None,
    ):
        """Plot the converged spectra of the measurement's SpectrumDeconvolutor

        Each converged spectrum is colored by the average potential of its
        convergent region, on the scale of the measurement's full potential range.

        Args:
            measurement (Measurement): The measurement, if different from
                self.measurement. A SpectrumDeconvolutor must be added to it.
            cmap_name (str): The name of the colormap to use. Defaults to "jet".
            make_colorbar (bool): Whether to make a colorbar.
            ax (matplotlib Axis): The axes to plot on. A new one is made by default.
        """
        measurement = measurement or self.measurement
        decon = self._get_deconvolutor(measurement)
        potentials = [np.mean(span) for span in decon.potentialspans]
        full_potentials = measurement.grab(measurement.U_name)[1]
        norm = mpl.colors.Normalize(
            vmin=np.min(full_potentials), vmax=np.max(full_potentials)
        )
        return super().plot_waterfall_vs(
            measurement=measurement,
            field=decon.converged_spectra,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.U_name,
            t=np.array(potentials),
            norm=norm,
        )

    def plot_fit_and_residuals(
        self,
        *,
        measurement=None,
        cmap_name="jet",
        make_colorbar=True,
        axes=None,
        tspan=None,
    ):
        """Plot the fit coefficients vs potential on top and the residuals below

        Uses the SpectrumDeconvolutor added to the measurement. Each coefficient is
        colored by the potential of its converged spectrum, and the residuals by the
        potential at which they are taken. The grey band shows the deconvolutor's
        `noise_bound`.

        Args:
            measurement (Measurement): The measurement, if different from
                self.measurement. A SpectrumDeconvolutor must be added to it.
            cmap_name (str): The name of the colormap to use. Defaults to "jet".
            make_colorbar (bool): Whether to make a colorbar for the residuals.
            axes (list of Axis): [coefficients axis, residuals axis]. New ones are
                made by default.
            tspan (iterable): The timespan of the residuals to plot
        """
        measurement = measurement or self.measurement
        decon = self._get_deconvolutor(measurement)
        n_converged = len(decon.converged_data)
        coefficients = [measurement[f"coefficient_{i}"] for i in range(n_converged)]
        residuals = measurement["residuals"]
        redox_potentials = [np.mean(span) for span in decon.potentialspans]

        if not axes:
            axes = self.new_two_panel_axes(n_bottom=1, n_top=1, emphasis="top")

        full_potentials = measurement.grab(measurement.U_name)[1]
        norm = mpl.colors.Normalize(
            vmin=np.min(full_potentials), vmax=np.max(full_potentials)
        )
        cmap = mpl.colormaps[cmap_name]
        for coefficient, redox_potential in zip(coefficients, redox_potentials):
            t = coefficient.tseries.t + coefficient.tseries.tstamp - measurement.tstamp
            axes[0].plot(
                measurement.grab_for_t(measurement.U_name, t=t),
                coefficient.data,
                color=cmap(norm(redox_potential)),
            )
        axes[0].set_ylabel("ΔA (OD)")
        axes[0].set_xlabel(measurement.U_name)

        super().plot_waterfall_vs(
            measurement=measurement,
            field=residuals,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=axes[1],
            vs=measurement.U_name,
            tspan=tspan,
        )
        if decon.noise_bound:
            axes[1].axhspan(
                -decon.noise_bound,
                decon.noise_bound,
                facecolor="gainsboro",
                edgecolor="black",
                linestyle="--",
            )
        return axes

    def plot_fit_reconstruction(
        self,
        *,
        measurement=None,
        cmap_name="jet",
        make_colorbar=True,
        ax=None,
        tspan=None,
    ):
        """Plot the spectra reconstructed by the SpectrumDeconvolutor's fit

        Args:
            measurement (Measurement): The measurement, if different from
                self.measurement. A SpectrumDeconvolutor must be added to it.
            cmap_name (str): The name of the colormap to use. Defaults to "jet".
            make_colorbar (bool): Whether to make a colorbar.
            ax (matplotlib Axis): The axes to plot on. A new one is made by default.
            tspan (iterable): The timespan of the spectra to plot
        """
        measurement = measurement or self.measurement
        self._get_deconvolutor(measurement)
        return super().plot_waterfall_vs(
            measurement=measurement,
            field=measurement["reconstructed"],
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.U_name,
            tspan=tspan,
        )
