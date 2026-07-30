"""Built-in adapters from Matplotlib plotter calls to ``PlotSpec`` objects."""

from .backends import register_plotter_adapter
from .ec_plotter import ECPlotter
from .ecms_plotter import ECMSPlotter, determine_tspan
from .ms_plotter import MSPlotter, MSSpectroPlotter, STANDARD_COLORS as MS_COLORS
from .nmr_plotter import NMRPlotter
from .plot_spec import (
    combine_plot_specs,
    ec_measurement_spec,
    ec_vs_potential_spec,
    ms_measurement_spec,
    spectrum_series_heatmap_spec,
    spectrum_series_stacked_spec,
    spectrum_series_waterfall_spec,
    spectrum_spec,
    time_series_spec,
    value_measurement_spec,
    xrd_spectrum_spec,
)
from .sec_plotter import ECOpticalPlotter, SECPlotter
from .spectrum_plotter import SpectrumPlotter, SpectrumSeriesPlotter
from .tpms_plotter import (
    STANDARD_COLORS as TP_COLORS,
    TPMSSpectroPlotter,
    TPMSPlotter,
    _get_y_unit_and_label,
)
from .value_plotter import ValuePlotter
from .xrd_plotter import XRDSpectrumPlotter
from .xrf_plotter import ECTRXRFPlotter, TRXRFPlotter


def _keyword_arguments(arguments):
    """Return ordinary and variadic keyword arguments in one dictionary."""
    parameters = dict(arguments)
    parameters.pop("args", None)
    for name in ("kwargs", "plot_kwargs"):
        parameters.update(parameters.pop(name, {}))
    return parameters


def _line_settings(parameters):
    """Return the common line settings from plot keyword arguments."""
    return {
        "line_style": parameters.get("linestyle", parameters.get("ls")),
        "line_width": parameters.get("linewidth", parameters.get("lw")),
    }


def _measurement(owner, parameters):
    """Return the explicitly supplied measurement or the bound owner."""
    return parameters.get("measurement") or owner


def _spectrum(owner, parameters):
    """Return the explicitly supplied spectrum or the bound owner."""
    return parameters.get("spectrum") or owner


def _spectrum_series(owner, parameters):
    """Return the explicitly supplied spectrum series or the bound owner."""
    return parameters.get("spectrum_series") or owner


def _value_measurement_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    measurement = _measurement(owner, parameters)
    plot_spec = value_measurement_spec(
        measurement,
        v_list=parameters.get("v_list"),
        tspan=parameters.get("tspan"),
        logscale=parameters.get("logscale", False),
    )
    plot_spec.show_legend = parameters.get("legend", True)
    return plot_spec


def _value_plot_adapter(owner, arguments):
    parameters = dict(arguments.get("kwargs", {}))
    positional = arguments.get("args", ())
    names = ("measurement", "v_list", "tspan", "ax", "legend", "logscale")
    parameters.update(zip(names, positional))
    return _value_measurement_adapter(owner, parameters)


def _ec_measurement_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    measurement = _measurement(owner, parameters)
    return ec_measurement_spec(
        measurement,
        tspan=parameters.get("tspan"),
        U_name=parameters.get("U_name") or parameters.get("V_str"),
        J_name=parameters.get("J_name") or parameters.get("J_str"),
        U_color=parameters.get("U_color") or parameters.get("V_color"),
        J_color=parameters.get("J_color"),
        **_line_settings(parameters),
    )


def _ec_vs_potential_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    measurement = _measurement(owner, parameters)
    return ec_vs_potential_spec(
        measurement,
        tspan=parameters.get("tspan"),
        U_name=parameters.get("U_name"),
        J_name=parameters.get("J_name"),
        color=parameters.get("color"),
        **_line_settings(parameters),
    )


def _ms_spec(measurement, parameters, color_map=None):
    """Build an MS plot description from bound method arguments."""
    return ms_measurement_spec(
        measurement,
        mass_list=parameters.get("mass_list"),
        mass_lists=parameters.get("mass_lists"),
        mol_list=parameters.get("mol_list"),
        mol_lists=parameters.get("mol_lists"),
        tspan=parameters.get("tspan"),
        tspan_bg=parameters.get("tspan_bg"),
        remove_background=parameters.get("remove_background"),
        unit=parameters.get("unit"),
        x_unit=parameters.get("x_unit"),
        logplot=parameters.get("logplot"),
        logdata=parameters.get("logdata", False),
        color_map=color_map or MS_COLORS,
        **_line_settings(parameters),
    )


def _ms_measurement_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    plot_spec = _ms_spec(_measurement(owner, parameters), parameters)
    plot_spec.show_legend = parameters.get("legend", True)
    return plot_spec


def _ms_spectro_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    measurement = _measurement(owner, parameters)
    ms_spec = _ms_spec(measurement, parameters)
    heatmap_spec = spectrum_series_heatmap_spec(
        measurement.spectrum_series,
        tspan=parameters.get("tspan"),
        xspan=parameters.get("xspan"),
        cmap_name=parameters.get("cmap_name", "inferno"),
        make_colorbar=parameters.get("make_colorbar", False),
        max_threshold=parameters.get("max_threshold"),
        min_threshold=parameters.get("min_threshold"),
        scanning_mask=parameters.get("scanning_mask"),
        vmin=parameters.get("vmin"),
        vmax=parameters.get("vmax"),
        x_unit=parameters.get("x_unit"),
    )
    plot_specs = (
        (ms_spec, heatmap_spec)
        if parameters.get("ms_data", "top") == "top"
        else (heatmap_spec, ms_spec)
    )
    plot_spec = combine_plot_specs(*plot_specs)
    plot_spec.show_legend = parameters.get("legend", True)
    plot_spec.row_heights = {
        "top": [3, 2],
        "bottom": [2, 3],
    }.get(parameters.get("emphasis"), [1, 1])
    return plot_spec


def _ecms_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    measurement = _measurement(owner, parameters)
    parameters["tspan"] = determine_tspan(parameters.get("tspan"), measurement)
    parameters["logplot"] = (
        not parameters.get("mass_lists")
        if parameters.get("logplot") is None
        else parameters["logplot"]
    )
    if parameters.get("removebackground") is not None:
        parameters["remove_background"] = parameters["removebackground"]
    ms_spec = _ms_spec(measurement, parameters)
    ec_spec = ec_measurement_spec(
        measurement,
        tspan=parameters["tspan"],
        U_name=parameters.get("U_name") or parameters.get("V_str"),
        J_name=parameters.get("J_name") or parameters.get("J_str"),
        U_color=parameters.get("U_color") or parameters.get("V_color"),
        J_color=parameters.get("J_color"),
        **_line_settings(parameters),
    )
    plot_spec = combine_plot_specs(ms_spec, ec_spec)
    plot_spec.show_legend = parameters.get("legend", True)
    plot_spec.row_heights = {
        "top": [3, 2],
        "bottom": [2, 3],
    }.get(parameters.get("emphasis"), [1, 1])
    return plot_spec


def _sec_spec(measurement, parameters, field=None):
    """Build a spectral heatmap and EC time plot."""
    heatmap_spec = spectrum_series_heatmap_spec(
        measurement.spectrum_series,
        field=field or parameters.get("field"),
        tspan=parameters.get("tspan"),
        xspan=parameters.get("xspan"),
        cmap_name=parameters.get("cmap_name", "inferno"),
        make_colorbar=parameters.get("make_colorbar", False),
        continuous=parameters.get("continuous"),
    )
    ec_spec = ec_measurement_spec(
        measurement,
        tspan=parameters.get("tspan"),
        U_name=parameters.get("U_name"),
        J_name=parameters.get("J_name"),
        U_color=parameters.get("U_color"),
        J_color=parameters.get("J_color"),
        **_line_settings(parameters),
    )
    plot_spec = combine_plot_specs(heatmap_spec, ec_spec)
    plot_spec.row_heights = [3, 2]
    return plot_spec


def _sec_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    return _sec_spec(_measurement(owner, parameters), parameters)


def _ec_optical_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    measurement = _measurement(owner, parameters)
    parameters["xspan"] = parameters.get("wlspan")
    field = measurement.calc_dOD(
        V_ref=parameters.get("V_ref"),
        t_ref=parameters.get("t_ref"),
    )
    return _sec_spec(measurement, parameters, field=field)


def _metadata_label_and_factors(measurement, names, meta_units=None):
    """Return an axis label and scale factors for metadata series."""
    factors = {}
    y_label = "signal"
    y_unit = "mixed"
    for name in names:
        data_series = measurement[name]
        y_label, y_unit, factors[name] = _get_y_unit_and_label(
            data_series,
            meta_units=meta_units,
        )
    return f"{y_label} / [{y_unit.strip()}]", factors


def _tp_spec(measurement, parameters, meta_units_name="TP_units"):
    """Build temperature and pressure traces against time."""
    temperature_name = parameters.get("T_name") or measurement.T_name
    pressure_name = parameters.get("P_name") or measurement.P_name
    temperature_names = parameters.get("T_names") or [temperature_name]
    pressure_names = parameters.get("P_names") or [pressure_name]
    meta_units = parameters.get(meta_units_name)
    temperature_label, temperature_factors = _metadata_label_and_factors(
        measurement,
        temperature_names,
        meta_units=meta_units,
    )
    pressure_label, pressure_factors = _metadata_label_and_factors(
        measurement,
        pressure_names,
        meta_units=meta_units,
    )
    colors = dict(TP_COLORS)
    if parameters.get("T_color"):
        colors[temperature_name] = parameters["T_color"]
    if parameters.get("P_color"):
        colors[pressure_name] = parameters["P_color"]
    return time_series_spec(
        measurement,
        left_names=temperature_names,
        right_names=pressure_names,
        tspan=parameters.get("tspan"),
        left_label=temperature_label,
        right_label=pressure_label,
        colors=colors,
        value_factors={**temperature_factors, **pressure_factors},
        x_unit=parameters.get("x_unit"),
        **_line_settings(parameters),
    )


def _tpms_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    measurement = _measurement(owner, parameters)
    ms_spec = _ms_spec(
        measurement,
        parameters,
        color_map={**MS_COLORS, **TP_COLORS},
    )
    plot_spec = combine_plot_specs(ms_spec, _tp_spec(measurement, parameters))
    plot_spec.show_legend = parameters.get("legend", True)
    plot_spec.row_heights = {
        "top": [3, 2],
        "bottom": [2, 3],
    }.get(parameters.get("emphasis"), [1, 1])
    return plot_spec


def _tpms_spectro_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    measurement = _measurement(owner, parameters)
    heatmap_spec = spectrum_series_heatmap_spec(
        measurement.spectrum_series,
        tspan=parameters.get("tspan"),
        xspan=parameters.get("xspan"),
        cmap_name=parameters.get("cmap_name", "inferno"),
        make_colorbar=parameters.get("make_colorbar", False),
        max_threshold=parameters.get("max_threshold"),
        min_threshold=parameters.get("min_threshold"),
        scanning_mask=parameters.get("scanning_mask"),
        vmin=parameters.get("vmin"),
        vmax=parameters.get("vmax"),
        x_unit=parameters.get("x_unit"),
    )
    ms_spec = _ms_spec(
        measurement,
        parameters,
        color_map={**MS_COLORS, **TP_COLORS},
    )
    plot_spec = combine_plot_specs(
        heatmap_spec,
        ms_spec,
        _tp_spec(measurement, parameters, meta_units_name="meta_units"),
    )
    plot_spec.show_legend = parameters.get("legend", True)
    return plot_spec


def _trxrf_spec(measurement, parameters):
    """Build a selected XRF signal against time."""
    y_name = parameters.get("y_name", "FF_over_I0")
    return time_series_spec(
        measurement,
        left_names=[y_name],
        tspan=parameters.get("tspan"),
        left_label=y_name,
        show_legend=False,
        **_line_settings(parameters),
    )


def _trxrf_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    return _trxrf_spec(_measurement(owner, parameters), parameters)


def _ec_trxrf_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    measurement = _measurement(owner, parameters)
    ec_spec = ec_measurement_spec(
        measurement,
        tspan=parameters.get("tspan"),
        U_name=parameters.get("U_name"),
        J_name=parameters.get("J_name"),
        U_color=parameters.get("U_color"),
        J_color=parameters.get("J_color"),
        **_line_settings(parameters),
    )
    return combine_plot_specs(_trxrf_spec(measurement, parameters), ec_spec)


def _spectrum_adapter(owner, arguments, inverted_x=False):
    parameters = _keyword_arguments(arguments)
    return spectrum_spec(
        _spectrum(owner, parameters),
        color=parameters.get("color"),
        inverted_x=inverted_x,
        **_line_settings(parameters),
    )


def _nmr_adapter(owner, arguments):
    return _spectrum_adapter(owner, arguments, inverted_x=True)


def _xrd_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    return xrd_spectrum_spec(
        _spectrum(owner, parameters),
        color=parameters.get("color"),
        **_line_settings(parameters),
    )


def _heatmap_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    return spectrum_series_heatmap_spec(
        _spectrum_series(owner, parameters),
        field=parameters.get("field"),
        tspan=parameters.get("tspan"),
        xspan=parameters.get("xspan"),
        cmap_name=parameters.get("cmap_name", "inferno"),
        make_colorbar=parameters.get("make_colorbar", False),
        t=parameters.get("t"),
        t_name=parameters.get("t_name"),
        max_threshold=parameters.get("max_threshold"),
        min_threshold=parameters.get("min_threshold"),
        vmin=parameters.get("vmin"),
        vmax=parameters.get("vmax"),
        scanning_mask=parameters.get("scanning_mask"),
        continuous=parameters.get("continuous"),
    )


def _waterfall_adapter(owner, arguments):
    parameters = _keyword_arguments(arguments)
    return spectrum_series_waterfall_spec(
        _spectrum_series(owner, parameters),
        field=parameters.get("field"),
        cmap_name=parameters.get("cmap_name", "jet"),
        make_colorbar=parameters.get("make_colorbar", True),
        t=parameters.get("t"),
        t_name=parameters.get("t_name"),
    )


def _stacked_adapter(owner, arguments):
    supplied_arguments = getattr(arguments, "supplied", ())
    parameters = _keyword_arguments(arguments)
    color = parameters.get("color") if "color" in supplied_arguments else None
    return spectrum_series_stacked_spec(
        _spectrum_series(owner, parameters),
        dt=parameters.get("dt"),
        t_list=parameters.get("t_list"),
        dn=parameters.get("dn"),
        index_list=parameters.get("index_list"),
        average=parameters.get("average", False),
        xspan=parameters.get("xspan"),
        xspan_bg=parameters.get("xspan_bg"),
        scale_mode=parameters.get("scale_mode", "auto"),
        scale_factor=parameters.get("scale_factor", 1),
        y_values=parameters.get("y_values", "time"),
        color=color,
        **_line_settings(parameters),
    )


def _register_builtin_adapters():
    """Register the plot methods covered by ixdat's shared plot descriptions."""
    registrations = (
        (ValuePlotter, "plot", _value_plot_adapter),
        (ValuePlotter, "plot_measurement", _value_measurement_adapter),
        (ECPlotter, "plot_measurement", _ec_measurement_adapter),
        (ECPlotter, "plot_vs_potential", _ec_vs_potential_adapter),
        (MSPlotter, "plot_measurement", _ms_measurement_adapter),
        (MSSpectroPlotter, "plot_measurement", _ms_spectro_adapter),
        (ECMSPlotter, "plot_measurement", _ecms_adapter),
        (SECPlotter, "plot_measurement", _sec_adapter),
        (ECOpticalPlotter, "plot_measurement", _ec_optical_adapter),
        (TPMSPlotter, "plot_measurement", _tpms_adapter),
        (TPMSSpectroPlotter, "plot_measurement", _tpms_spectro_adapter),
        (TRXRFPlotter, "plot_measurement", _trxrf_adapter),
        (ECTRXRFPlotter, "plot_measurement", _ec_trxrf_adapter),
        (ECTRXRFPlotter, "plot_vs_potential", _ec_vs_potential_adapter),
        (SpectrumPlotter, "plot", _spectrum_adapter),
        (NMRPlotter, "plot", _nmr_adapter),
        (XRDSpectrumPlotter, "plot", _xrd_adapter),
        (SpectrumSeriesPlotter, "heat_plot", _heatmap_adapter),
        (SpectrumSeriesPlotter, "plot_waterfall", _waterfall_adapter),
        (SpectrumSeriesPlotter, "plot_stacked_spectra", _stacked_adapter),
    )
    for plotter_class, method_name, adapter in registrations:
        register_plotter_adapter(plotter_class, method_name, adapter)


_register_builtin_adapters()
