"""Tests for backend-neutral plot descriptions and renderers."""

import inspect
from pathlib import Path

from matplotlib import pyplot as plt
import numpy as np
import pytest

from ixdat import Measurement, Spectrum
from ixdat.data_series import TimeSeries, ValueSeries
from ixdat.plotters import (
    AxisSpec,
    LineTrace,
    MatplotlibRenderer,
    MPLPlotter,
    PanelSpec,
    PlotSpec,
    PlotterBackendWarning,
    available_plotter_backends,
    register_plotter_adapter,
    register_plotter_backend,
    unregister_plotter_backend,
)
from ixdat.plotters.value_plotter import ValuePlotter
from ixdat.spectra import SpectrumSeries
from ixdat.techniques.ec import ECMeasurement
from ixdat.techniques.ec_ms import ECMSMeasurement
from ixdat.techniques.ms import MSMeasurement, MSSpectroMeasurement
from ixdat.techniques.nmr import NMRSpectrum
from ixdat.techniques.reactor import ReactorSpectroMeasurement
from ixdat.techniques.spectroelectrochemistry import SpectroECMeasurement
from ixdat.techniques.xrf import ECTRXRFMeasurement, TRXRFMeasurement


class ExampleRenderer:
    """Return the plot description unchanged."""

    def render(self, plot_spec, **kwargs):
        return plot_spec


class MatplotlibOnlyPlotter(MPLPlotter):
    """Provide one plot method without a registered plot-description adapter."""

    def __init__(self, measurement=None):
        super().__init__()
        self.measurement = measurement

    def plot_measurement(self, measurement=None, ax=None):
        measurement = measurement or self.measurement
        if ax is None:
            ax = self.new_ax()
        time, value = measurement.grab("value")
        ax.plot(time, value)
        return ax


@pytest.fixture
def measurement():
    """Return a generic measurement with one value series."""
    time = TimeSeries(
        name="time",
        unit_name="s",
        data=np.array([0.0, 1.0]),
        tstamp=0,
    )
    value = ValueSeries(
        name="value",
        unit_name="A",
        data=np.array([2.0, 3.0]),
        tseries=time,
    )
    return Measurement(name="generic", series_list=[time, value], tstamp=0)


@pytest.fixture
def ec_measurement():
    """Return a small EC measurement with angle brackets in the current name."""
    time = TimeSeries(
        name="time",
        unit_name="s",
        data=np.array([0.0, 1.0]),
        tstamp=0,
    )
    potential = ValueSeries(
        name="Ewe/V",
        unit_name="V",
        data=np.array([0.1, 0.2]),
        tseries=time,
    )
    current = ValueSeries(
        name="<I>/mA",
        unit_name="mA",
        data=np.array([1.0, 2.0]),
        tseries=time,
    )
    return ECMeasurement(
        name="ec",
        series_list=[time, potential, current],
        aliases={
            "t": ["time"],
            "raw_potential": ["Ewe/V"],
            "raw_current": ["<I>/mA"],
        },
        tstamp=0,
    )


@pytest.fixture
def ms_measurement():
    """Return a small mass-spectrometry measurement."""
    time = TimeSeries(
        name="time",
        unit_name="s",
        data=np.array([0.0, 1.0]),
        tstamp=0,
    )
    signal = ValueSeries(
        name="M2",
        unit_name="A",
        data=np.array([1e-9, 2e-9]),
        tseries=time,
    )
    return MSMeasurement(
        name="ms",
        series_list=[time, signal],
        tstamp=0,
    )


@pytest.fixture
def ecms_measurement():
    """Return a small EC-MS measurement."""
    time = TimeSeries(
        name="time",
        unit_name="s",
        data=np.array([0.0, 1.0]),
        tstamp=0,
    )
    series = [
        ValueSeries(
            name=name,
            unit_name=unit,
            data=np.asarray(data),
            tseries=time,
        )
        for name, unit, data in [
            ("Ewe/V", "V", [0.1, 0.2]),
            ("<I>/mA", "mA", [1.0, 2.0]),
            ("M2", "A", [1e-9, 2e-9]),
        ]
    ]
    return ECMSMeasurement(
        name="ecms",
        series_list=[time] + series,
        aliases={
            "t": ["time"],
            "raw_potential": ["Ewe/V"],
            "raw_current": ["<I>/mA"],
        },
        tstamp=0,
    )


@pytest.fixture
def spectrum_series():
    """Return two spectra sampled on the same x-axis."""
    first = Spectrum.from_data(
        x=np.array([1.0, 2.0]),
        y=np.array([3.0, 4.0]),
        name="first",
        tstamp=0,
    )
    second = Spectrum.from_data(
        x=np.array([1.0, 2.0]),
        y=np.array([5.0, 6.0]),
        name="second",
        tstamp=1,
    )
    return SpectrumSeries.from_spectrum_list([first, second])


@pytest.fixture
def sec_measurement(spectrum_series):
    """Return a spectroelectrochemistry measurement."""
    time = TimeSeries(
        name="time",
        unit_name="s",
        data=np.array([0.0, 1.0]),
        tstamp=0,
    )
    series = [
        ValueSeries(
            name=name,
            unit_name=unit,
            data=np.asarray(data),
            tseries=time,
        )
        for name, unit, data in [
            ("Ewe/V", "V", [0.1, 0.2]),
            ("<I>/mA", "mA", [1.0, 2.0]),
        ]
    ]
    return SpectroECMeasurement(
        name="sec",
        series_list=[time] + series,
        spectrum_series=spectrum_series,
        aliases={
            "t": ["time"],
            "raw_potential": ["Ewe/V"],
            "raw_current": ["<I>/mA"],
        },
        tstamp=0,
    )


@pytest.fixture
def ms_spectro_measurement(spectrum_series):
    """Return an MS measurement with a spectrum series."""
    time = TimeSeries(
        name="time",
        unit_name="s",
        data=np.array([0.0, 1.0]),
        tstamp=0,
    )
    signal = ValueSeries(
        name="M2",
        unit_name="A",
        data=np.array([1e-9, 2e-9]),
        tseries=time,
    )
    return MSSpectroMeasurement(
        name="ms-spectro",
        series_list=[time, signal],
        spectrum_series=spectrum_series,
        tstamp=0,
    )


@pytest.fixture
def reactor_spectro_measurement(spectrum_series):
    """Return a TP-MS measurement with a spectrum series."""
    time = TimeSeries(
        name="time",
        unit_name="s",
        data=np.array([0.0, 1.0]),
        tstamp=0,
    )
    value_series = [
        ValueSeries(
            name=name,
            unit_name=unit,
            data=np.asarray(data),
            tseries=time,
        )
        for name, unit, data in [
            ("M2", "A", [1e-9, 2e-9]),
            ("M32", "A", [4e-10, 8e-10]),
            ("temperature", "K", [300.0, 310.0]),
            ("pressure", "bar", [1.0, 1.1]),
        ]
    ]
    return ReactorSpectroMeasurement(
        name="reactor-spectro",
        series_list=[time] + value_series,
        spectrum_series=spectrum_series,
        tstamp=0,
    )


@pytest.fixture
def xrf_measurement():
    """Return a small time-resolved XRF measurement."""
    time = TimeSeries("time", "s", np.array([0.0, 1.0]), tstamp=0)
    signal = ValueSeries(
        "FF_over_I0",
        "",
        np.array([0.2, 0.3]),
        tseries=time,
    )
    return TRXRFMeasurement(
        name="xrf",
        series_list=[time, signal],
        tstamp=0,
    )


@pytest.fixture
def ec_xrf_measurement():
    """Return a small EC and time-resolved XRF measurement."""
    time = TimeSeries("time", "s", np.array([0.0, 1.0]), tstamp=0)
    series = [
        ValueSeries(name, unit, np.asarray(data), tseries=time)
        for name, unit, data in [
            ("Ewe/V", "V", [0.1, 0.2]),
            ("<I>/mA", "mA", [1.0, 2.0]),
            ("FF_over_I0", "", [0.2, 0.3]),
        ]
    ]
    return ECTRXRFMeasurement(
        name="ec-xrf",
        series_list=[time] + series,
        aliases={
            "t": ["time"],
            "raw_potential": ["Ewe/V"],
            "raw_current": ["<I>/mA"],
        },
        tstamp=0,
    )


@pytest.fixture
def backend_name():
    """Provide a backend name and remove it after each test."""
    name = "example"
    yield name
    if name in available_plotter_backends():
        unregister_plotter_backend(name)


def test_default_plotter_is_unchanged(measurement):
    """Measurements keep their technique-specific Matplotlib plotter."""
    assert isinstance(measurement.plotter, ValuePlotter)


def test_backend_routing_stays_outside_matplotlib_plotter_classes(measurement):
    """Bound data methods expose renderers while plotter classes stay Matplotlib-only."""
    assert "backend" not in inspect.signature(ValuePlotter.plot_measurement).parameters
    assert "backend" in inspect.signature(measurement.plot).parameters


def test_missing_adapter_warns_and_uses_matplotlib(measurement):
    """A plotter without an adapter remains available through Matplotlib."""
    measurement = Measurement(
        name=measurement.name,
        series_list=measurement.series_list,
        tstamp=measurement.tstamp,
        plotter=MatplotlibOnlyPlotter(measurement=measurement),
    )

    with pytest.warns(PlotterBackendWarning, match="Using Matplotlib"):
        axis = measurement.plot(backend="plotly")

    assert list(axis.lines[0].get_xdata()) == [0.0, 1.0]
    assert list(axis.lines[0].get_ydata()) == [2.0, 3.0]


def test_adapter_can_add_renderer_support_later(measurement, backend_name):
    """A plotter gains renderer coverage through an external adapter."""

    class AdaptedPlotter(MatplotlibOnlyPlotter):
        pass

    @register_plotter_adapter(AdaptedPlotter, "plot_measurement")
    def build_plot_spec(owner, arguments):
        assert arguments["measurement"] is None
        assert arguments["ax"] is None
        time, value = owner.grab("value")
        return PlotSpec(
            [
                PanelSpec(
                    traces=[LineTrace(time, value, name="value")],
                    x_axis=AxisSpec("time / [s]"),
                    left_axis=AxisSpec("value / [A]"),
                )
            ]
        )

    adapted_measurement = Measurement(
        name=measurement.name,
        series_list=measurement.series_list,
        tstamp=measurement.tstamp,
        plotter=AdaptedPlotter(measurement=measurement),
    )
    register_plotter_backend(backend_name, ExampleRenderer)

    result = adapted_measurement.plot(backend=backend_name)

    assert isinstance(result, PlotSpec)
    assert list(result.panels[0].traces[0].y) == [2.0, 3.0]


def test_generic_matplotlib_path_keeps_its_existing_output(measurement):
    """The generic Matplotlib path keeps its line data and axis settings."""
    axis = measurement.plot()

    assert list(axis.lines[0].get_xdata()) == [0.0, 1.0]
    assert list(axis.lines[0].get_ydata()) == [2.0, 3.0]
    assert axis.get_xlabel() == ""
    assert axis.get_ylabel() == ""


def test_normalized_matplotlib_name_uses_the_direct_path(measurement):
    """Whitespace and capitalization preserve the Matplotlib result."""
    axis = measurement.plot(backend=" Matplotlib ")

    assert list(axis.lines[0].get_xdata()) == [0.0, 1.0]
    assert axis.get_xlabel() == ""
    assert axis.get_ylabel() == ""


def test_builtin_backends_are_available():
    """Matplotlib and Plotly ship as renderers."""
    assert available_plotter_backends() == ("matplotlib", "plotly")


def test_generic_measurement_can_plot_with_plotly(measurement):
    """A generic measurement renders the same values with Plotly."""
    pytest.importorskip("plotly")

    figure = measurement.plot(backend="plotly")

    assert list(figure.data[0].x) == [0.0, 1.0]
    assert list(figure.data[0].y) == [2.0, 3.0]
    assert figure.data[0].name == "value"


def test_plotly_uses_shared_legend_setting(measurement):
    """The plot specification carries the legend setting."""
    pytest.importorskip("plotly")

    figure = measurement.plot(backend="plotly", legend=False)

    assert not figure.layout.showlegend


def test_plotly_follows_matplotlib_figure_dimensions(measurement):
    """Plotly uses the canvas dimensions configured for Matplotlib."""
    pytest.importorskip("plotly")

    with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": 100}):
        figure = measurement.plot(backend="plotly")

    assert figure.layout.width == 500
    assert figure.layout.height == 500


def test_ec_plotly_labels_show_literal_angle_brackets(ec_measurement):
    """Plotly escapes HTML-sensitive current names."""
    pytest.importorskip("plotly")

    figure = ec_measurement.plot(backend="plotly")

    assert figure.data[1].name == "&lt;I&gt;/mA"
    assert figure.layout.yaxis2.title.text == "&lt;I&gt;/mA"
    assert not figure.data[0].showlegend
    assert not figure.data[1].showlegend
    assert figure.layout.xaxis.ticks == "outside"
    assert figure.layout.xaxis.ticklen == 5
    assert figure.layout.yaxis.tickcolor == figure.data[0].line.color
    assert figure.layout.yaxis.tickfont.color == figure.data[0].line.color
    assert figure.layout.yaxis.title.font.color == figure.data[0].line.color
    assert figure.layout.yaxis2.tickcolor == figure.data[1].line.color
    assert figure.layout.yaxis2.tickfont.color == figure.data[1].line.color
    assert figure.layout.yaxis2.title.font.color == figure.data[1].line.color


def test_plotly_adds_an_ec_plot_to_a_plain_figure(ec_measurement):
    """The figure argument accepts a Plotly figure without a subplot grid."""
    go = pytest.importorskip("plotly.graph_objects")
    figure = go.Figure(layout={"width": 800, "height": 500})

    result = ec_measurement.plot(backend="plotly", figure=figure)

    assert result is figure
    assert [trace.name for trace in figure.data] == ["Ewe/V", "&lt;I&gt;/mA"]
    assert figure.data[1].yaxis == "y2"
    assert figure.layout.width == 800
    assert figure.layout.height == 500


def test_plotly_rejects_an_incompatible_subplot_grid(
    measurement,
    ec_measurement,
):
    """A supplied subplot grid produces a clear compatibility error."""
    pytest.importorskip("plotly")
    figure = measurement.plot(backend="plotly")

    with pytest.raises(ValueError, match="incompatible subplot grid"):
        ec_measurement.plot(backend="plotly", figure=figure)


def test_matplotlib_renderer_applies_axis_colors_and_trace_legends():
    """The shared Matplotlib renderer follows axis and trace descriptions."""
    plot_spec = PlotSpec(
        [
            PanelSpec(
                traces=[
                    LineTrace([0, 1], [1, 2], "left", "black", show_legend=False),
                    LineTrace(
                        [0, 1],
                        [2, 3],
                        "right",
                        "red",
                        y_axis="right",
                        show_legend=False,
                    ),
                ],
                x_axis=AxisSpec("time / [s]"),
                left_axis=AxisSpec("left", color="black"),
                right_axis=AxisSpec("right", color="red"),
            )
        ],
        show_legend=True,
    )

    left_axis, right_axis = MatplotlibRenderer().render(plot_spec)

    assert left_axis.yaxis.label.get_color() == "black"
    assert right_axis.yaxis.label.get_color() == "red"
    assert left_axis.get_legend() is None
    assert right_axis.get_legend() is None


def test_cv_style_plot_can_use_plotly(ec_measurement):
    """Current against potential uses the shared line description."""
    pytest.importorskip("plotly")

    figure = ec_measurement.plot_vs_potential(backend="plotly")

    assert list(figure.data[0].x) == [0.1, 0.2]
    assert list(figure.data[0].y) == [1.0, 2.0]


def test_spectrum_can_plot_with_plotly():
    """A spectrum carries common line settings into Plotly."""
    pytest.importorskip("plotly")
    spectrum = Spectrum.from_data(
        x=np.array([1.0, 2.0]),
        y=np.array([3.0, 4.0]),
        name="spectrum",
        tstamp=0,
    )

    figure = spectrum.plot(
        backend="plotly",
        color="green",
        linewidth=2,
        linestyle="--",
    )

    assert list(figure.data[0].x) == [1.0, 2.0]
    assert list(figure.data[0].y) == [3.0, 4.0]
    assert figure.data[0].line.color == "#008000"
    assert figure.data[0].line.width == 2
    assert figure.data[0].line.dash == "dash"
    assert not figure.data[0].showlegend


def test_spectrum_matplotlib_path_keeps_its_existing_output():
    """A spectrum continues to use its direct Matplotlib path."""
    spectrum = Spectrum.from_data(
        x=np.array([1.0, 2.0]),
        y=np.array([3.0, 4.0]),
        x_name="x",
        y_name="y",
        name="spectrum",
        tstamp=0,
    )

    axis = spectrum.plot(color="green", linewidth=2)

    assert list(axis.lines[0].get_xdata()) == [1.0, 2.0]
    assert list(axis.lines[0].get_ydata()) == [3.0, 4.0]
    assert axis.lines[0].get_color() == "green"
    assert axis.lines[0].get_linewidth() == 2
    assert axis.get_xlabel() == "x"
    assert axis.get_ylabel() == "y"


def test_ms_measurement_can_plot_with_plotly(ms_measurement):
    """Mass traces use the same backend-neutral line description."""
    pytest.importorskip("plotly")

    figure = ms_measurement.plot(backend="plotly", logplot=False)

    assert figure.data[0].name == "M2"
    assert list(figure.data[0].y) == [1e-9, 2e-9]


def test_ms_plotly_uses_the_matplotlib_log_floor():
    """The two backends clip low positive MS values at the same floor."""
    pytest.importorskip("plotly")
    time = TimeSeries("time", "s", np.array([0.0, 1.0]), tstamp=0)
    signal = ValueSeries(
        "M2",
        "A",
        np.array([1e-16, 1e-12]),
        tseries=time,
    )
    measurement = MSMeasurement(
        name="low-signal-ms",
        series_list=[time, signal],
        tstamp=0,
    )

    figure = measurement.plot(backend="plotly")

    assert list(figure.data[0].y) == [1e-14, 1e-12]


def test_ecms_measurement_composes_plotly_panels(ecms_measurement):
    """EC-MS combines MS, potential, and current specifications."""
    pytest.importorskip("plotly")

    figure = ecms_measurement.plot(backend="plotly", logplot=False)

    assert [trace.name for trace in figure.data] == [
        "M2",
        "Ewe/V",
        "&lt;I&gt;/mA",
    ]
    assert figure.layout.xaxis.matches == "x2"
    assert not figure.layout.xaxis.showticklabels
    assert figure.layout.xaxis.title.text is None
    assert figure.layout.xaxis2.title.text == "time / [s]"
    assert figure.layout.hovermode == "x unified"
    assert figure.layout.hoversubplots == "axis"
    assert {trace.xaxis for trace in figure.data} == {"x2"}
    assert not figure.layout.xaxis.showgrid
    assert not figure.layout.xaxis2.showgrid
    assert not figure.layout.yaxis.showgrid
    assert not figure.layout.yaxis2.showgrid
    assert len(figure.layout.shapes) == 2
    assert all(
        shape.type == "rect"
        and shape.xref == "paper"
        and shape.yref == "paper"
        and shape.line.color == "black"
        and shape.line.width == 1
        for shape in figure.layout.shapes
    )
    top_height = figure.layout.yaxis.domain[1] - figure.layout.yaxis.domain[0]
    bottom_height = figure.layout.yaxis2.domain[1] - figure.layout.yaxis2.domain[0]
    assert top_height / bottom_height == pytest.approx(1.5)
    assert figure.layout.width == 640
    assert figure.layout.height == 480


def test_sec_measurement_composes_plotly_panels(sec_measurement):
    """Spectroelectrochemistry combines a heatmap with EC traces."""
    pytest.importorskip("plotly")

    figure = sec_measurement.plot(backend="plotly")

    assert [trace.type for trace in figure.data] == [
        "heatmap",
        "scatter",
        "scatter",
    ]


def test_time_resolved_xrf_measurement_can_plot_with_plotly(xrf_measurement):
    """A time-resolved XRF measurement renders its selected signal."""
    pytest.importorskip("plotly")

    figure = xrf_measurement.plot(backend="plotly")

    assert [trace.name for trace in figure.data] == ["FF_over_I0"]
    assert list(figure.data[0].y) == [0.2, 0.3]
    assert not figure.data[0].showlegend
    assert figure.layout.yaxis.tickcolor == "#000000"


def test_ec_xrf_measurement_composes_plotly_panels(ec_xrf_measurement):
    """EC and time-resolved XRF values share one Plotly time axis."""
    pytest.importorskip("plotly")

    figure = ec_xrf_measurement.plot(backend="plotly")

    assert [trace.name for trace in figure.data] == [
        "FF_over_I0",
        "Ewe/V",
        "&lt;I&gt;/mA",
    ]
    assert {trace.xaxis for trace in figure.data} == {"x2"}
    assert figure.layout.hovermode == "x unified"


def test_nmr_spectrum_reverses_plotly_x_axis():
    """The shared axis description keeps the NMR direction."""
    pytest.importorskip("plotly")
    spectrum = NMRSpectrum.from_data(
        x=np.array([1.0, 2.0]),
        y=np.array([3.0, 4.0]),
        name="nmr",
        tstamp=0,
    )

    figure = spectrum.plot(backend="plotly")

    assert figure.layout.xaxis.autorange == "reversed"


def test_nmr_spectrum_keeps_matplotlib_axis_direction():
    """NMR keeps its direct Matplotlib path and reversed x-axis."""
    spectrum = NMRSpectrum.from_data(
        x=np.array([1.0, 2.0]),
        y=np.array([3.0, 4.0]),
        name="nmr",
        tstamp=0,
    )

    axis = spectrum.plot()

    assert axis.xaxis_inverted()
    assert list(axis.lines[0].get_xdata()) == [1.0, 2.0]
    assert list(axis.lines[0].get_ydata()) == [3.0, 4.0]


def test_xrd_error_band_can_plot_with_plotly():
    """XRD error values render as a filled Plotly trace."""
    pytest.importorskip("plotly")
    data_file = Path(__file__).parents[2] / "test_data/xrd/no_header.xye"
    spectrum = Spectrum.read(data_file, reader="xrdxy")

    figure = spectrum.plot(backend="plotly")

    assert [trace.type for trace in figure.data] == ["scatter", "scatter"]
    assert figure.data[1].fill == "toself"
    assert figure.data[1].fillcolor.startswith("rgba(31, 119, 180,")
    assert figure.data[0].line.color == "#1f77b4"
    assert not figure.data[0].showlegend


def test_spectrum_series_can_make_continuous_plotly_heatmap(spectrum_series):
    """A continuous spectrum series renders one Plotly heatmap."""
    pytest.importorskip("plotly")

    figure = spectrum_series.plot(backend="plotly", continuous=True)

    assert figure.data[0].type == "heatmap"
    assert np.asarray(figure.data[0].z).shape == (2, 2)


def test_spectrum_series_discrete_heatmap_uses_spectrum_durations(spectrum_series):
    """A discrete heatmap gives each spectrum its recorded time interval."""
    pytest.importorskip("plotly")
    spectrum_series.durations = [0.25, 0.5]

    figure = spectrum_series.plot(backend="plotly", continuous=False)

    assert [trace.type for trace in figure.data] == ["heatmap", "heatmap"]
    assert list(figure.data[0].x) == [0.0, 0.25]
    assert list(figure.data[1].x) == [1.0, 1.5]
    assert all(np.asarray(trace.z).shape == (2, 1) for trace in figure.data)


def test_spectrum_series_inferred_intervals_share_one_heatmap(spectrum_series):
    """Contiguous inferred intervals render as one heatmap without seams."""
    pytest.importorskip("plotly")

    figure = spectrum_series.plot(backend="plotly", continuous=False)

    assert len(figure.data) == 1
    assert list(figure.data[0].x) == [0.0, 1.0]
    assert np.asarray(figure.data[0].z).shape == (2, 1)


def test_spectrum_series_can_make_plotly_waterfall(spectrum_series):
    """Waterfall lines use a continuous Plotly color scale."""
    pytest.importorskip("plotly")

    figure = spectrum_series.plotter.plot_waterfall(backend="plotly")

    assert len(figure.data) == 3
    assert figure.data[0].line.color != figure.data[1].line.color
    assert figure.data[2].marker.showscale
    assert figure.data[2].marker.colorbar.title.text == spectrum_series.t_name


def test_spectrum_series_can_make_plotly_stacked_plot(spectrum_series):
    """Stacked Plotly spectra keep their offsets and use distinct colors."""
    pytest.importorskip("plotly")

    figure = spectrum_series.plotter.plot_stacked_spectra(
        backend="plotly",
        index_list=[0, 1],
    )

    assert len(figure.data) == 2
    assert list(figure.data[0].y) == [3.0, 4.0]
    assert list(figure.data[1].y) == [6.0, 7.0]
    assert figure.data[0].line.color != figure.data[1].line.color


def test_spectrum_series_stacked_matplotlib_default_stays_black(spectrum_series):
    """Stacked Matplotlib spectra retain their black default."""
    axis = spectrum_series.plotter.plot_stacked_spectra(index_list=[0, 1])

    assert {line.get_color() for line in axis.lines} == {"k"}


def test_ms_spectro_measurement_composes_plotly_panels(ms_spectro_measurement):
    """Spectro-MS combines line and heatmap specifications."""
    pytest.importorskip("plotly")

    figure = ms_spectro_measurement.plot(backend="plotly", logplot=False)

    assert [trace.type for trace in figure.data] == ["scatter", "heatmap"]


def test_reactor_spectro_measurement_composes_plotly_panels(
    reactor_spectro_measurement,
):
    """Spectro-TP-MS combines spectra, MS, temperature, and pressure."""
    pytest.importorskip("plotly")

    figure = reactor_spectro_measurement.plot(backend="plotly", logplot=False)

    assert [trace.type for trace in figure.data] == [
        "heatmap",
        "scatter",
        "scatter",
        "scatter",
        "scatter",
    ]
    mass_colors = {
        trace.name: trace.line.color
        for trace in figure.data
        if trace.name in {"M2", "M32"}
    }
    assert set(mass_colors) == {"M2", "M32"}
    assert mass_colors["M2"] != mass_colors["M32"]
    for trace in figure.data[-2:]:
        axis_name = trace.yaxis.replace("y", "yaxis", 1)
        axis = figure.layout[axis_name]
        assert not trace.showlegend
        assert axis.tickcolor == trace.line.color
        assert axis.tickfont.color == trace.line.color
        assert axis.title.font.color == trace.line.color
    assert figure.layout.width == 640
    assert figure.layout.height == 480
    assert figure.layout.hoversubplots == "axis"
    assert {trace.xaxis for trace in figure.data} == {"x3"}
    assert figure.layout.xaxis.matches == "x3"
    assert figure.layout.xaxis2.matches == "x3"
    assert not figure.layout.xaxis.showticklabels
    assert not figure.layout.xaxis2.showticklabels
    assert figure.layout.xaxis3.title.text == "time / [s]"
    assert figure.layout.yaxis3.title.text == "temperature / [K]"
    assert figure.layout.yaxis4.title.text == "pressure / [bar]"
    assert figure.layout.yaxis.automargin
    assert figure.layout.yaxis.title.standoff == 14
    assert all(
        not figure.layout[axis_name].showgrid
        for axis_name in figure.layout
        if axis_name.startswith(("xaxis", "yaxis"))
    )
    assert len(figure.layout.shapes) == 3
    assert all(
        shape.type == "rect" and shape.line.color == "black" and shape.line.width == 1
        for shape in figure.layout.shapes
    )


def test_spectro_tpms_converts_every_shared_time_axis(
    reactor_spectro_measurement,
):
    """MS, heatmap, and metadata traces use one requested time unit."""
    pytest.importorskip("plotly")

    figure = reactor_spectro_measurement.plot(
        backend="plotly",
        logplot=False,
        x_unit="min",
    )

    assert figure.layout.xaxis3.title.text == "time / [min]"
    assert all(np.allclose(trace.x, [0.0, 1.0 / 60]) for trace in figure.data)


def test_paired_ms_axes_name_and_color_each_single_trace(
    reactor_spectro_measurement,
):
    """A pair of MS axes identifies its traces through labels and colors."""
    pytest.importorskip("plotly")

    figure = (
        reactor_spectro_measurement.plotter.tpms_plotter.ms_plotter.plot_measurement(
            measurement=reactor_spectro_measurement,
            backend="plotly",
            mass_lists=[["M2"], ["M32"]],
            logplot=False,
        )
    )

    assert figure.layout.yaxis.title.text == "M2 signal / [A]"
    assert figure.layout.yaxis2.title.text == "M32 signal / [A]"
    for trace in figure.data:
        axis_name = trace.yaxis.replace("y", "yaxis", 1)
        assert not trace.showlegend
        assert figure.layout[axis_name].tickcolor == trace.line.color


def test_tpms_plotly_labels_and_converts_metadata_units(
    reactor_spectro_measurement,
):
    """TP-MS labels metadata units and scales requested pressure units."""
    pytest.importorskip("plotly")

    figure = reactor_spectro_measurement.plotter.tpms_plotter.plot_measurement(
        measurement=reactor_spectro_measurement,
        backend="plotly",
        logplot=False,
        TP_units={"pressure": "mbar"},
    )
    traces = {trace.name: trace for trace in figure.data}

    assert figure.layout.yaxis2.title.text == "temperature / [K]"
    assert figure.layout.yaxis3.title.text == "pressure / [mbar]"
    assert list(traces["temperature"].y) == [300.0, 310.0]
    assert list(traces["pressure"].y) == [1000.0, 1100.0]


def test_custom_renderer_receives_plot_specs_across_data_types(
    measurement,
    ms_measurement,
    backend_name,
):
    """One renderer handles generic, technique, and spectrum plotters."""
    register_plotter_backend(backend_name, ExampleRenderer)
    spectrum = Spectrum.from_data(
        x=np.array([1.0, 2.0]),
        y=np.array([3.0, 4.0]),
        name="spectrum",
        tstamp=0,
    )

    assert isinstance(measurement.plot(backend=backend_name), PlotSpec)
    assert isinstance(
        ms_measurement.plot(backend=backend_name, logplot=False),
        PlotSpec,
    )
    assert isinstance(spectrum.plot(backend=backend_name), PlotSpec)


def test_reader_has_no_backend_responsibility(measurement, backend_name):
    """Reader inputs stay independent of the selected renderer."""

    class ExampleReader:
        def read(self, path_to_file, cls=None, **kwargs):
            assert kwargs == {}
            return measurement

    register_plotter_backend(backend_name, ExampleRenderer)
    read_measurement = Measurement.read("unused", reader=ExampleReader())

    assert isinstance(read_measurement.plot(backend=backend_name), PlotSpec)


def test_duplicate_registration_requires_overwrite(backend_name):
    """A backend name changes ownership only with explicit overwrite."""
    register_plotter_backend(backend_name, ExampleRenderer)

    with pytest.raises(ValueError, match="overwrite=True"):
        register_plotter_backend(backend_name, ExampleRenderer)

    register_plotter_backend(backend_name, ExampleRenderer, overwrite=True)


def test_builtin_backend_cannot_be_removed():
    """Built-in renderer names remain available."""
    with pytest.raises(ValueError, match="built-in"):
        unregister_plotter_backend("plotly")


def test_builtin_backend_cannot_be_replaced():
    """Built-in renderer names keep their shipped implementation."""
    with pytest.raises(ValueError, match="built-in"):
        register_plotter_backend("plotly", ExampleRenderer, overwrite=True)


def test_invalid_renderer_is_rejected():
    """Registration requires a renderer class with render()."""
    with pytest.raises(TypeError, match="must be a class"):
        register_plotter_backend("invalid", object())


def test_unknown_backend_lists_available_names(measurement):
    """An unknown name produces an actionable error."""
    with pytest.raises(ValueError, match="matplotlib, plotly"):
        measurement.plot(backend="missing")
