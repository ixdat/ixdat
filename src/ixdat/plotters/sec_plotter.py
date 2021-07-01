import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from .base_mpl_plotter import MPLPlotter
from .ec_plotter import ECPlotter


class SECPlotter(MPLPlotter):
    """An spectroelectrochemsitry (SEC) matplotlib plotter."""

    def __init__(self, measurement=None):
        """Initiate the plotter with its default Meausurement to plot"""
        self.measurement = measurement
        self.ec_plotter = ECPlotter(measurement=measurement)

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
        tseries = dOD_series.axes_series[0]
        t = tseries.data
        wlseries = dOD_series.axes_series[1]
        wl = wlseries.data
        dOD_data = dOD_series.data

        if tspan:
            t_mask = np.logical_and(tspan[0] < t, t < tspan[-1])
            t = t[t_mask]
            dOD_data = dOD_data[t_mask, :]
        if wlspan:
            wl_mask = np.logical_and(wlspan[0] < wl, wl < wlspan[-1])
            wl = wl[wl_mask]
            dOD_data = dOD_data[:, wl_mask]

        axes[0].imshow(
            dOD_data.swapaxes(0, 1),
            cmap=cmap_name,
            aspect="auto",
            extent=(t[0], t[-1], wl[0], wl[-1]),
        )

        axes[0].set_xlabel(
            (measurement.t_str if hasattr(measurement, "t_str") else None)
            or tseries.name
        )
        axes[0].set_ylabel(wlseries.name)
        if make_colorbar:
            cmap = mpl.cm.get_cmap(cmap_name)
            norm = mpl.colors.Normalize(vmin=np.min(dOD_data), vmax=np.max(dOD_data))
            cb = plt.colorbar(
                mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                ax=[axes[0], axes[1]],
                use_gridspec=True,
                anchor=(0.75, 0),
            )
            cb.set_label("$\Delta$ opdical density")
        return axes

    def plot_waterfall(
        self, measurement=None, cmap_name="jet", make_colorbar=True, V_ref=0.66, ax=None
    ):
        measurement = measurement or self.measurement
        if not ax:
            ax = self.new_ax()
        dOD = measurement.calc_dOD(V_ref=V_ref)
        dOD_data = dOD.data
        v = measurement.v
        wl = measurement.wl

        cmap = mpl.cm.get_cmap(cmap_name)
        norm = mpl.colors.Normalize(vmin=np.min(v), vmax=np.max(v))

        for i, v_i in enumerate(v):
            spec = dOD_data[i]
            color = cmap(norm(v_i))
            ax.plot(wl, spec, color=color)

        ax.set_xlabel(measurement.wavelength.name)
        ax.set_ylabel(dOD.name)

        if make_colorbar:
            cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
            cb.set_label(measurement.potential.name)

        return ax

    def plot_vs_potential(self, *args, **kwargs):
        return self.ec_plotter.plot_vs_potential(*args, **kwargs)
