"""Plotters for spectroelectrochemistry. Makes use of those in spectrum_plotter.py"""

import matplotlib as mpl
from . import ECPlotter, SpectrumSeriesPlotter, SpectroMeasurementPlotter,MPLPlotter
from ..exceptions import SeriesNotFoundError
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

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
    ):
        """Plot an SECMeasurement as spectra colored based on potential.

        The top panel is a heat plot with wavelength on y-axis and color representing
        spectrum. At most one of V_ref and t_ref should be given, and if neither are
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
        dOD = measurement.calc_dOD(V_ref=V_ref, t_ref=t_ref)

        return super().plot_waterfall_vs(
            measurement=self.measurement,
            field=dOD,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            ax=ax,
            vs=measurement.U_name,
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
        return axes

class StaircaseSECPlotter(MPLPlotter):
    def __init__(self, measurement=None):

        self.measurement = measurement
        self.smooth = self.measurement.smooth
        self.with_voltage_interval = self.measurement.with_voltage_interval
        self.with_wavelength_interval = self.measurement.with_wavelength_interval
        self.data_config = self.measurement.data_config
        self.around_WL = self.measurement.around_WL
        self.around_V = self.measurement.around_V

    def plot_waterfall(
        self,
        fig=None,
        ax=None,
        dataFrame=None,
        make_colorbar=True,  # GUI
        cmap_name="jet",  # GUI
        t_ref=None,
        WL_interval=None,
        V_interval=None,
        smooth_param=None,
        reduce_factor = 1,
        **kwargs
    ):  # GUI
        """waterfall plot

        Args:
            fig : the matplotlib figure object
            ax : the matplotlib ax object which is on the fig
            dataFrame : the data to plot
            make_colorbar : whether add a colorbar to the fig
            cmap_name : color map name
            t_ref : time reference (not in used currently)
            WL_interval : wavelength interval
            V_interval : voltage interval
            smooth_param : smooth parameters in the form of [window length, polyorder]

        Returns:
            the fig
        """

        # return an fig of waterfall plot
        if dataFrame is None:
            dataFrame = self.measurement.SEC_dataFrame.copy() * 1000
        dataFrame = self.data_config(
            dataFrame=dataFrame, WL_interval=WL_interval, V_interval=V_interval
        )

        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)

        cmap = plt.get_cmap(cmap_name, len(dataFrame.columns))
        sm = plt.cm.ScalarMappable(
            cmap=cmap,
            norm=plt.Normalize(vmin=dataFrame.columns[0], vmax=dataFrame.columns[-1]),
        )

        for i in range(0, len(dataFrame.columns),reduce_factor):
            if smooth_param != None:
                y_data = savgol_filter(
                    dataFrame[dataFrame.columns[i]], smooth_param[0], smooth_param[1]
                )
            else:
                y_data = dataFrame[dataFrame.columns[i]]
            ax.plot(dataFrame.index, y_data, c=cmap(i),**kwargs)

        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Absorption (m $\Delta$ O.D.)")

        if make_colorbar:
            cbar = ax.figure.colorbar(sm, ax=ax)
            cbar.set_label("U vs RHE (V)")

        return fig,ax, cbar
    


    
    def plot_differential(
        self,
        fig=None,
        ax=None,
        dataFrame=None,
        differential_interval=0.1,
        make_colorbar=True,
        cmap_name="jet",
        normalize=True,
        normalize_method="max",
        skip_negative=True,
        WL_interval=None,
        V_interval=None,
        smooth_param=None,
        capacitance_func = None
    ):
        """differential plot

        Args:
            fig : the matplotlib figure object
            ax : the matplotlib ax object which is on the fig
            dataFrame : the data to plot
            differential_interval : differential interval in Volts
            make_colorbar : whether add a colorbar to the fig
            cmap_name (str) : color map name
            normalize (boolean) : normalize differential spectra with respect to the maximum
            normalize_method (str) : max or norm
            skip_negative (boolean) : skip spectra with a negative average
            WL_interval (list): wavelength interval
            V_interval (list) : voltage interval
            smooth_param (list) : smooth parameters in the form of [window length, polyorder]

        Returns:
            the fig
        """
        if dataFrame is None:

            dataFrame = self.measurement.SEC_dataFrame.copy() * 1000
        dataFrame = self.data_config(
            dataFrame=dataFrame, WL_interval=WL_interval, V_interval=V_interval
        )
        if V_interval is None:
            V_interval = [dataFrame.columns[0], dataFrame.columns[-1]]
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)
        fig.suptitle(
            f"Voltage range: {V_interval[0]} V - {V_interval[1]} V, differentiate interval = {differential_interval} V",
            fontsize=10,
        )

        self.differential_dict = {}

        differential_interval = int(
            round(
                differential_interval / (dataFrame.columns[1] - dataFrame.columns[0]), 0
            )
        )

        cmap = plt.get_cmap(cmap_name, len(dataFrame.columns))
        # the indexes for the dataframe.columns are the first i in the following loop and the last one
        sm = plt.cm.ScalarMappable(
            cmap=cmap,
            norm=plt.Normalize(
                vmin=dataFrame.columns[
                    [
                        i
                        for i in range(
                            differential_interval,
                            len(dataFrame.columns),
                            differential_interval,
                        )
                    ][0]
                ],
                vmax=dataFrame.columns[
                    [
                        i
                        for i in range(
                            differential_interval,
                            len(dataFrame.columns),
                            differential_interval,
                        )
                    ][-1]
                ],
            ),
        )
        for i in range(
            differential_interval, len(dataFrame.columns), differential_interval
        ):

            original_spectrum = (
                dataFrame[dataFrame.columns[i]]
                - dataFrame[dataFrame.columns[i - differential_interval]]
            )

            if sum(original_spectrum) < 0 and skip_negative:
                continue

            if smooth_param is not None:
                spectrum_to_plot = pd.Series(
                    savgol_filter(original_spectrum, smooth_param[0], smooth_param[1]),
                    index=original_spectrum.index,
                )
            else:
                spectrum_to_plot = original_spectrum

            if normalize:
                if normalize_method == "max":
                    spectrum_to_plot = spectrum_to_plot / np.abs(np.max(spectrum_to_plot))
                elif normalize_method == "norm":
                    spectrum_to_plot = spectrum_to_plot / np.linalg.norm(spectrum_to_plot)
                elif normalize_method == 'capacitance':
                    if capacitance_func is not None:
                        spectrum_to_plot = spectrum_to_plot/capacitance_func(dataFrame.columns[i])

            ax.plot(spectrum_to_plot, c=cmap(i), linewidth=0.8, alpha=0.9)

            self.differential_dict[dataFrame.columns[i]] = spectrum_to_plot

        ax.set_xlabel("Wavelength (nm)")
        if normalize:
            ax.set_ylabel("Normalized\n Differential Absorption")
        else:
            ax.set_ylabel("Differential Absorption (m $\Delta$ O.D.)")
        if make_colorbar:
            cbar = ax.figure.colorbar(sm, ax=ax)
            cbar.set_label("U vs RHE (V)")
            return fig,ax,cbar
        else:
            return fig,ax

        
        

    def plot_spectrum(
        self,
        fig=None,
        ax=None,
        V_range=[0, 1],
        normalize=True,
        WL_interval=None,
        V_interval=None,
        smooth_param=None,
        normalize_method = "max",
        **kwargs
    ):
        """single spectrum plot

        Args:
            fig : the matplotlib figure object
            ax : the matplotlib ax object which is on the fig
            dataFrame : the data to plot
            V_range : the voltage range that defines the spectrum based on the SEC data
            normalize (boolean) : normalize differential spectra with respect to the maximum
            WL_interval (list): wavelength interval
            V_interval (list) : voltage interval
            smooth_param (list) : smooth parameters in the form of [window length, polyorder]

        """
        dataFrame = self.measurement.SEC_dataFrame.copy() * 1000
        dataFrame = self.data_config(
            dataFrame=dataFrame, WL_interval=WL_interval, V_interval=V_interval
        )
        spectrum_to_plot = (
            dataFrame.loc[:, self.around_V(V_range[1])]
            - dataFrame.loc[:, self.around_V(V_range[0])]
        )

        if smooth_param is not None:
            spectrum_to_plot = pd.Series(
                savgol_filter(spectrum_to_plot, smooth_param[0], smooth_param[1]),
                index=spectrum_to_plot.index,
            )
        if normalize:
            if normalize_method == "max":
                spectrum_to_plot = spectrum_to_plot / np.abs(np.max(spectrum_to_plot))
            elif normalize_method == "norm":
                spectrum_to_plot = spectrum_to_plot / np.linalg.norm(spectrum_to_plot)

        ax.plot(spectrum_to_plot,**kwargs)






    def plot_stack_differential(
        self,
        fig=None,
        ax=None,
        dataFrame=None,
        differential_interval=0.1,
        make_colorbar=True,
        cmap_name="jet",
        normalize=True,
        normalize_method="max",
        skip_negative=True,
        WL_interval=None,
        V_interval=None,
        smooth_param=None,
        capacitance_func = None,
        spacing = 0.018,
    ):
        """differential plot

        Args:
            fig : the matplotlib figure object
            ax : the matplotlib ax object which is on the fig
            dataFrame : the data to plot
            differential_interval : differential interval in Volts
            make_colorbar : whether add a colorbar to the fig
            cmap_name (str) : color map name
            normalize (boolean) : normalize differential spectra with respect to the maximum
            normalize_method (str) : max or norm
            skip_negative (boolean) : skip spectra with a negative average
            WL_interval (list): wavelength interval
            V_interval (list) : voltage interval
            smooth_param (list) : smooth parameters in the form of [window length, polyorder]

        Returns:
            the fig
        """
        if dataFrame is None:

            dataFrame = self.measurement.SEC_dataFrame.copy() * 1000
        dataFrame = self.data_config(
            dataFrame=dataFrame, WL_interval=WL_interval, V_interval=V_interval
        )
        if V_interval is None:
            V_interval = [dataFrame.columns[0], dataFrame.columns[-1]]
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)
        # fig.suptitle(
        #     f"Voltage range: {V_interval[0]} V - {V_interval[1]} V, differentiate interval = {differential_interval} V",
        #     fontsize=10,
        # )

        self.differential_dict = {}

        differential_interval = int(
            round(
                differential_interval / (dataFrame.columns[1] - dataFrame.columns[0]), 0
            )
        )

        cmap = plt.get_cmap(cmap_name, len(dataFrame.columns)+50)
        # the indexes for the dataframe.columns are the first i in the following loop and the last one
        sm = plt.cm.ScalarMappable(
            cmap=cmap,
            norm=plt.Normalize(
                vmin=dataFrame.columns[
                    [
                        i
                        for i in range(
                            differential_interval,
                            len(dataFrame.columns),
                            differential_interval,
                        )
                    ][0]
                ],
                vmax=dataFrame.columns[
                    [
                        i
                        for i in range(
                            differential_interval,
                            len(dataFrame.columns),
                            differential_interval,
                        )
                    ][-1]
                ],
            ),
        )
        n = 0
        for i in range(
            differential_interval, len(dataFrame.columns), differential_interval
        ):

            original_spectrum = (
                dataFrame[dataFrame.columns[i]]
                - dataFrame[dataFrame.columns[i - differential_interval]]
            )

            if sum(original_spectrum) < 0 and skip_negative:
                continue

            if smooth_param is not None:
                spectrum_to_plot = pd.Series(
                    savgol_filter(original_spectrum, smooth_param[0], smooth_param[1]),
                    index=original_spectrum.index,
                )
            else:
                spectrum_to_plot = original_spectrum

            if normalize:
                if normalize_method == "max":
                    spectrum_to_plot = spectrum_to_plot / np.abs(np.max(spectrum_to_plot))
                elif normalize_method == "norm":
                    spectrum_to_plot = spectrum_to_plot / np.linalg.norm(spectrum_to_plot)
                elif normalize_method == "zero_norm":
                    spectrum_to_plot = spectrum_to_plot - np.min(spectrum_to_plot)
                    spectrum_to_plot = spectrum_to_plot / np.linalg.norm(spectrum_to_plot)
                    
                elif normalize_method == 'capacitance':
                    if capacitance_func is not None:
                        spectrum_to_plot = spectrum_to_plot/capacitance_func(dataFrame.columns[i])
            line_lable = str(round(dataFrame.columns[i - differential_interval],2))+"V - "+str(round(dataFrame.columns[i],2))+"V"
            ax.plot(spectrum_to_plot+n*spacing, c=cmap(i), linewidth=1, alpha=1)
            ax.annotate(text=line_lable,xy = [735,spectrum_to_plot.iloc[-1]+n*spacing+0.005],c=cmap(i),fontsize = 7)
            

            self.differential_dict[dataFrame.columns[i]] = original_spectrum
            n+=1
        
        ax.set_xlabel("Wavelength (nm)",fontsize = 12)
        ax.set_yticks([])

        if normalize:
            ax.set_ylabel("Normalized\n Differential Absorption",fontsize = 12)
        else:
            ax.set_ylabel("Differential Absorption (m $\Delta$ O.D.)",fontsize = 12)
        if make_colorbar:
            cbar = ax.figure.colorbar(sm, ax=ax)
            cbar.set_label("U vs RHE (V)")
            return fig,ax,cbar
        else:
            return fig,ax