from .plotting_tools import (
    color_axis,
    add_colorbar,
    get_color_from_cmap,
    # FIXME: the following should be Calculators.
    #   see https://github.com/ixdat/ixdat/issues/164.
    smooth_vector,
    calc_linear_background,
)
from .base_mpl_plotter import MPLPlotter
from .backends import (
    PlotterBackendWarning,
    available_plotter_backends,
    bind_plotter,
    get_renderer,
    register_plotter_adapter,
    register_plotter_backend,
    unregister_plotter_backend,
)
from .plot_spec import (
    AxisSpec,
    ColorScaleSpec,
    ErrorBand,
    HeatmapTrace,
    LineTrace,
    PanelSpec,
    PlotSpec,
)
from .renderers import MatplotlibRenderer, PlotlyRenderer
from .value_plotter import ValuePlotter
from .spectrum_plotter import (
    SpectrumPlotter,
    SpectrumSeriesPlotter,
    SpectroMeasurementPlotter,
)
from .ec_plotter import ECPlotter, CVDiffPlotter
from .ms_plotter import MSPlotter, MSSpectroPlotter
from .ecms_plotter import ECMSPlotter
from .sec_plotter import SECPlotter, ECOpticalPlotter
from .tpms_plotter import TPMSPlotter, TPMSSpectroPlotter

# Importing this module registers ixdat's built-in plot-description adapters after
# their plotter classes are available.
from . import plot_adapters  # noqa: F401, E402
