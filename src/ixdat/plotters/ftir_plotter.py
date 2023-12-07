# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 18:28:47 2023

@author: Søren
"""

from .spectrum_plotter import SpectrumSeriesPlotter


class FTIRPlotter(SpectrumSeriesPlotter):

    def heat_plot(
        self,
        *,
        spectrum_series=None,
        field=None,
        tspan=None,
        xspan=None,
        ax=None,
        cmap_name="inferno",
        make_colorbar=False,
        t=None,
        t_name=None,
        max_threshold=None,
        min_threshold=None,
        scanning_mask=None,
        vmin=None,
        vmax=None,
    ):
        """
        Plot a spectrum series with `t` on the horizontal axis, `x` on the vertical axis,
        and color representing `y`.

        Args:
            spectrum_series (SpectrumSeries): The spectrum series to be plotted, if
                different from self.spectrum_series.
            field (Field): The field to be plotted, if different from
                spectrum_series.field
            tspan (iterable): The span of the time data to plot
            xspan (iterable): The span of the spectral data to plot
            ax (mpl.Axis): The axes to plot on. A new one is made by default

            cmap_name (str): The name of the colormap to use. Defaults to "inferno", see
                https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html#sequential
            make_colorbar (bool): Whether to make a colorbar.
                FIXME: colorbar at present mis-alignes axes
            t (numpy array): Time data to use if not the data in spectrum_series
            t_name (str): Name of time variable if not the one in spectrum_series
            max_threshold (float): Maximum value to display.
                Values above are set to zero.
            min_threshold (float): Minimum value to display.
                Values below are set to 0.
            scanning_mask (list): List of booleans to exclude from scanning variable
                before plotting data by setting y values to 0 (zero).
            vmin (float): minimum value to represent in colours.
            vmax (float): maximum value to represent in colours.
        """

        return super().heat_plot(
            spectrum_series=spectrum_series,
            field=field,
            tspan=tspan,
            xspan=xspan,
            ax=ax,
            cmap_name=cmap_name,
            make_colorbar=make_colorbar,
            t=t,
            t_name=t_name,
            max_threshold=max_threshold,
            min_threshold=min_threshold,
            scanning_mask=scanning_mask,
            vmin=vmin,
            vmax=vmax,
        )

    def waterfall_plot(
        self,
        spectrum_series=None,
        ax=None,
        offset=0.2
    ):
        print(offset)
        spectrum_series = spectrum_series or self.spectrum_series
        if not ax:
            ax = self.new_ax()
            ax.set_xlabel(spectrum_series.xseries.name)
            ax.set_ylabel(spectrum_series.field.name)

        for (n, spectrum) in enumerate(spectrum_series):
            ax.plot(spectrum.x, spectrum.y + n*offset)

        return ax