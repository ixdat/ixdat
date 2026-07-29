"""Backend-neutral descriptions of plots."""

import numpy as np

from ..exceptions import SeriesNotFoundError
from .plotting_tools import get_indeces_and_times


DEFAULT_LINE_COLORS = (
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
)
MS_MIN_SIGNAL = 1e-14


class AxisSpec:
    """Describe one axis without referring to a plotting library."""

    def __init__(self, label=None, scale="linear", inverted=False, color=None):
        self.label = label
        self.scale = scale
        self.inverted = inverted
        self.color = color


class LineTrace:
    """Describe one line."""

    def __init__(
        self,
        x,
        y,
        name=None,
        color=None,
        y_axis="left",
        line_style=None,
        line_width=None,
        color_value=None,
        color_range=None,
        color_scale=None,
        show_legend=True,
    ):
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.name = name
        self.color = color
        self.y_axis = y_axis
        self.line_style = line_style
        self.line_width = line_width
        self.color_value = color_value
        self.color_range = color_range
        self.color_scale = color_scale
        self.show_legend = show_legend


class ErrorBand:
    """Describe a filled interval around a line."""

    def __init__(self, x, lower, upper, color=None, alpha=0.3, y_axis="left"):
        self.x = np.asarray(x)
        self.lower = np.asarray(lower)
        self.upper = np.asarray(upper)
        self.color = color
        self.alpha = alpha
        self.y_axis = y_axis


class HeatmapTrace:
    """Describe values shown as colors across two axes."""

    def __init__(
        self,
        x,
        y,
        z,
        colorscale="inferno",
        colorbar=False,
        colorbar_label=None,
        zmin=None,
        zmax=None,
    ):
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.z = np.asarray(z)
        self.colorscale = colorscale
        self.colorbar = colorbar
        self.colorbar_label = colorbar_label
        self.zmin = zmin
        self.zmax = zmax


class ColorScaleSpec:
    """Describe a continuous color scale shared by line traces."""

    def __init__(self, label, value_range, colorscale, visible=True):
        self.label = label
        self.value_range = value_range
        self.colorscale = colorscale
        self.visible = visible


class PanelSpec:
    """Describe one panel and the traces drawn in it."""

    def __init__(
        self,
        traces=None,
        x_axis=None,
        left_axis=None,
        right_axis=None,
        color_scale=None,
    ):
        self.traces = traces or []
        self.x_axis = x_axis or AxisSpec()
        self.left_axis = left_axis or AxisSpec()
        self.right_axis = right_axis
        self.color_scale = color_scale


class PlotSpec:
    """Describe a complete plot as a list of panels."""

    def __init__(
        self,
        panels,
        title=None,
        show_legend=None,
        shared_x_axes=False,
        row_heights=None,
    ):
        self.panels = panels
        self.title = title
        self.show_legend = show_legend
        self.shared_x_axes = shared_x_axes
        self.row_heights = row_heights


def combine_plot_specs(*plot_specs):
    """Join complete plot descriptions into one multi-panel description."""
    return PlotSpec(
        panels=[
            panel
            for plot_spec in plot_specs
            if plot_spec is not None
            for panel in plot_spec.panels
        ],
        shared_x_axes=True,
    )


def time_series_spec(
    measurement,
    left_names,
    right_names=None,
    tspan=None,
    left_label=None,
    right_label=None,
    colors=None,
    value_factors=None,
    x_unit=None,
    line_style=None,
    line_width=None,
    show_legend=None,
):
    """Describe named measurement values against time on one or two y-axes."""
    traces = []
    colors = colors or {}
    value_factors = value_factors or {}
    paired_single_traces = len(left_names) == 1 and len(right_names or []) == 1
    trace_show_legend = not paired_single_traces if show_legend is None else show_legend
    left_color = colors.get(left_names[0]) or DEFAULT_LINE_COLORS[0]
    right_color = (
        (colors.get(right_names[0]) or DEFAULT_LINE_COLORS[1]) if right_names else None
    )
    x_factor, shown_x_unit = _time_unit_factor(x_unit)
    for index, name in enumerate(left_names):
        t, value = measurement.grab(name, tspan=tspan)
        traces.append(
            LineTrace(
                t * x_factor,
                value * value_factors.get(name, 1),
                name,
                colors.get(name)
                or DEFAULT_LINE_COLORS[index % len(DEFAULT_LINE_COLORS)],
                line_style=line_style,
                line_width=line_width,
                show_legend=trace_show_legend,
            )
        )
    for index, name in enumerate(right_names or []):
        t, value = measurement.grab(name, tspan=tspan)
        traces.append(
            LineTrace(
                t * x_factor,
                value * value_factors.get(name, 1),
                name,
                colors.get(name)
                or DEFAULT_LINE_COLORS[
                    (index + len(left_names)) % len(DEFAULT_LINE_COLORS)
                ],
                y_axis="right",
                line_style=line_style,
                line_width=line_width,
                show_legend=trace_show_legend,
            )
        )
    return PlotSpec(
        [
            PanelSpec(
                traces=traces,
                x_axis=AxisSpec(f"time / [{shown_x_unit}]"),
                left_axis=AxisSpec(
                    left_label or ", ".join(left_names),
                    color=left_color if paired_single_traces else None,
                ),
                right_axis=(
                    AxisSpec(
                        right_label or ", ".join(right_names),
                        color=right_color if paired_single_traces else None,
                    )
                    if right_names
                    else None
                ),
            )
        ]
    )


def value_measurement_spec(
    measurement,
    v_list=None,
    tspan=None,
    logscale=False,
):
    """Describe a generic measurement's values against time."""
    traces = []
    for v_name in v_list or measurement.value_names:
        try:
            t, value = measurement.grab(v_name, tspan=tspan)
        except SeriesNotFoundError as error:
            print(f"WARNING!!! {error}")
            continue
        traces.append(LineTrace(t, value, name=v_name))
    return PlotSpec(
        [
            PanelSpec(
                traces=traces,
                x_axis=AxisSpec("time / [s]"),
                left_axis=AxisSpec(scale="log" if logscale else "linear"),
            )
        ]
    )


def ec_measurement_spec(
    measurement,
    tspan=None,
    U_name=None,
    J_name=None,
    U_color=None,
    J_color=None,
    line_style=None,
    line_width=None,
):
    """Describe potential and current against time."""
    U_name = U_name or measurement.U_name
    J_name = J_name or measurement.J_name
    U_color = U_color or "black"
    J_color = J_color or "red"
    t_u, potential = measurement.grab(U_name, tspan=tspan)
    t_j, current = measurement.grab(J_name, tspan=tspan)
    return PlotSpec(
        [
            PanelSpec(
                traces=[
                    LineTrace(
                        t_u,
                        potential,
                        U_name,
                        U_color,
                        line_style=line_style,
                        line_width=line_width,
                        show_legend=False,
                    ),
                    LineTrace(
                        t_j,
                        current,
                        J_name,
                        J_color,
                        y_axis="right",
                        line_style=line_style,
                        line_width=line_width,
                        show_legend=False,
                    ),
                ],
                x_axis=AxisSpec("time / [s]"),
                left_axis=AxisSpec(U_name, color=U_color),
                right_axis=AxisSpec(J_name, color=J_color),
            )
        ]
    )


def ec_vs_potential_spec(
    measurement,
    tspan=None,
    U_name=None,
    J_name=None,
    color=None,
    line_style=None,
    line_width=None,
):
    """Describe current against potential."""
    U_name = U_name or measurement.U_name
    J_name = J_name or measurement.J_name
    t_u, potential = measurement.grab(U_name, tspan=tspan)
    t_j, current = measurement.grab(J_name, tspan=tspan)
    return PlotSpec(
        [
            PanelSpec(
                traces=[
                    LineTrace(
                        potential,
                        np.interp(t_u, t_j, current),
                        J_name,
                        color or "black",
                        line_style=line_style,
                        line_width=line_width,
                    )
                ],
                x_axis=AxisSpec(U_name),
                left_axis=AxisSpec(J_name),
            )
        ]
    )


def ms_measurement_spec(
    measurement,
    mass_list=None,
    mass_lists=None,
    mol_list=None,
    mol_lists=None,
    tspan=None,
    tspan_bg=None,
    remove_background=None,
    unit=None,
    x_unit=None,
    logplot=True,
    logdata=False,
    color_map=None,
    line_style=None,
    line_width=None,
):
    """Describe mass signals or calibrated molecular fluxes against time."""
    quantified, value_groups = _ms_value_groups(
        measurement,
        mass_list,
        mass_lists,
        mol_list,
        mol_lists,
    )
    units = _as_group_values(unit, len(value_groups))
    background_spans = _as_group_values(tspan_bg, len(value_groups), span=True)
    traces = []
    axis_specs = []
    paired_single_traces = len(value_groups) == 2 and all(
        len(group) == 1 for group in value_groups
    )
    remove_background = not logplot if remove_background is None else remove_background
    x_factor, shown_x_unit = _time_unit_factor(x_unit)

    for group_index, values in enumerate(value_groups):
        group_unit = units[group_index] or ("mol/s" if quantified else "A")
        unit_factor = _ms_unit_factor(group_unit, quantified, measurement)
        display_unit = f"ln({group_unit})" if logdata else group_unit
        axis_color = None
        axis_value_names = []
        for value_name in values:
            value_name_string = getattr(value_name, "name", value_name)
            axis_value_names.append(value_name_string)
            trace_color = (color_map or {}).get(
                value_name_string
            ) or DEFAULT_LINE_COLORS[len(traces) % len(DEFAULT_LINE_COLORS)]
            if paired_single_traces:
                axis_color = trace_color
            if quantified:
                t, value = measurement.grab_flux(
                    value_name,
                    tspan=tspan,
                    tspan_bg=background_spans[group_index],
                    remove_background=remove_background,
                    include_endpoints=False,
                )
            else:
                t, value = measurement.grab_signal(
                    value_name,
                    tspan=tspan,
                    tspan_bg=background_spans[group_index],
                    remove_background=remove_background,
                    include_endpoints=False,
                )
            value = np.array(value, copy=True)
            if logplot:
                value[value < MS_MIN_SIGNAL] = MS_MIN_SIGNAL
            if logdata:
                value = np.log(value * unit_factor) / unit_factor
            traces.append(
                LineTrace(
                    t * x_factor,
                    value * unit_factor,
                    name=value_name_string,
                    color=trace_color,
                    y_axis="right" if group_index else "left",
                    line_style=line_style,
                    line_width=line_width,
                    show_legend=not paired_single_traces,
                )
            )
        axis_specs.append(
            AxisSpec(
                (
                    f"{axis_value_names[0]} signal / [{display_unit}]"
                    if paired_single_traces
                    else f"signal / [{display_unit}]"
                ),
                scale="linear" if logdata or not logplot else "log",
                color=axis_color,
            )
        )
    return PlotSpec(
        [
            PanelSpec(
                traces=traces,
                x_axis=AxisSpec(f"time / [{shown_x_unit}]"),
                left_axis=axis_specs[0],
                right_axis=axis_specs[1] if len(axis_specs) > 1 else None,
            )
        ]
    )


def spectrum_spec(
    spectrum,
    color=None,
    inverted_x=False,
    line_style=None,
    line_width=None,
):
    """Describe a spectrum as one line."""
    return PlotSpec(
        [
            PanelSpec(
                traces=[
                    LineTrace(
                        spectrum.x,
                        spectrum.y,
                        name=spectrum.y_name,
                        color=color,
                        line_style=line_style,
                        line_width=line_width,
                        show_legend=False,
                    )
                ],
                x_axis=AxisSpec(spectrum.x_name, inverted=inverted_x),
                left_axis=AxisSpec(spectrum.y_name),
            )
        ]
    )


def xrd_spectrum_spec(spectrum, color=None, line_style=None, line_width=None):
    """Describe an XRD spectrum with its error interval when present."""
    color = color or DEFAULT_LINE_COLORS[0]
    traces = [
        LineTrace(
            spectrum.x,
            spectrum.y,
            spectrum.y_name,
            color,
            line_style=line_style,
            line_width=line_width,
            show_legend=False,
        )
    ]
    if spectrum.y_err is not None:
        traces.append(
            ErrorBand(
                spectrum.x,
                spectrum.y - spectrum.y_err,
                spectrum.y + spectrum.y_err,
                color=color,
            )
        )
    return PlotSpec(
        [
            PanelSpec(
                traces=traces,
                x_axis=AxisSpec(spectrum.x_name),
                left_axis=AxisSpec(spectrum.y_name),
            )
        ]
    )


def spectrum_series_heatmap_spec(
    spectrum_series,
    field=None,
    tspan=None,
    xspan=None,
    cmap_name="inferno",
    make_colorbar=False,
    t=None,
    t_name=None,
    max_threshold=None,
    min_threshold=None,
    vmin=None,
    vmax=None,
    scanning_mask=None,
    x_unit=None,
    continuous=None,
):
    """Describe a spectrum series as a heatmap."""
    field = field or spectrum_series.field
    data = np.array(field.data, copy=True)
    xseries = field.axes_series[1]
    x = np.asarray(xseries.data)
    t = np.asarray(t if t is not None else field.axes_series[0].t)
    t_name = t_name or field.axes_series[0].name

    if max_threshold is not None:
        data[data > max_threshold] = max_threshold
    if min_threshold is not None:
        data[data < min_threshold] = min_threshold
    if scanning_mask is not None and np.any(scanning_mask):
        data[:, scanning_mask] = 0
    if xspan:
        x_mask = np.logical_and(xspan[0] < x, x < xspan[-1])
        x = x[x_mask]
        data = data[:, x_mask]
    shown_zmin = np.min(data) if vmin is None else vmin
    shown_zmax = np.max(data) if vmax is None else vmax
    if x_unit:
        x_factor, shown_x_unit = _time_unit_factor(x_unit)
        t_name = f"time / [{shown_x_unit}]"
    else:
        x_factor = 1
    continuous = spectrum_series.continuous if continuous is None else continuous

    if continuous:
        if tspan:
            t_mask = np.logical_and(min(tspan) < t, t < max(tspan))
            t = t[t_mask]
            data = data[t_mask, :]
            if tspan[0] > tspan[-1]:
                t = np.flip(t)
                data = np.flip(data, axis=0)
        traces = [
            HeatmapTrace(
                x=t * x_factor,
                y=x,
                z=data.T,
                colorscale=cmap_name,
                colorbar=make_colorbar,
                colorbar_label=field.name,
                zmin=shown_zmin,
                zmax=shown_zmax,
            )
        ]
    else:
        traces = []
        durations = spectrum_series.durations
        inferred_durations = durations is None or all(
            duration is None for duration in durations
        )
        if inferred_durations:
            interval_indices = [
                index
                for index, t_i in enumerate(t[:-1])
                if not tspan or min(tspan) <= t_i <= max(tspan)
            ]
            if interval_indices:
                edges = [t[index] for index in interval_indices]
                edges.append(t[interval_indices[-1] + 1])
                traces.append(
                    HeatmapTrace(
                        x=np.asarray(edges) * x_factor,
                        y=x,
                        z=data[interval_indices].T,
                        colorscale=cmap_name,
                        colorbar=make_colorbar,
                        colorbar_label=field.name,
                        zmin=shown_zmin,
                        zmax=shown_zmax,
                    )
                )
        else:
            for index, t_i in enumerate(t):
                if tspan and (t_i < min(tspan) or t_i > max(tspan)):
                    continue
                duration = durations[index]
                if duration is None:
                    if index + 1 == len(t):
                        break
                    t_f = t[index + 1]
                else:
                    t_f = t_i + duration
                traces.append(
                    HeatmapTrace(
                        x=np.array([t_i, t_f]) * x_factor,
                        y=x,
                        z=data[index][:, np.newaxis],
                        colorscale=cmap_name,
                        colorbar=make_colorbar and not traces,
                        colorbar_label=field.name,
                        zmin=shown_zmin,
                        zmax=shown_zmax,
                    )
                )

    return PlotSpec(
        [
            PanelSpec(
                traces=traces,
                x_axis=AxisSpec(t_name),
                left_axis=AxisSpec(xseries.name),
            )
        ]
    )


def spectrum_series_waterfall_spec(
    spectrum_series,
    field=None,
    cmap_name="jet",
    make_colorbar=True,
    t=None,
    t_name=None,
):
    """Describe every spectrum as a line colored by its acquisition time."""
    field = field or spectrum_series.field
    x = field.axes_series[1].data
    t = np.asarray(t if t is not None else field.axes_series[0].t)
    t_name = t_name or field.axes_series[0].name
    color_range = (np.min(t), np.max(t))
    traces = [
        LineTrace(
            x,
            field.data[index],
            name=str(t_i),
            color_value=t_i,
            color_range=color_range,
            color_scale=cmap_name,
            show_legend=False,
        )
        for index, t_i in enumerate(t)
    ]
    return PlotSpec(
        [
            PanelSpec(
                traces=traces,
                x_axis=AxisSpec(field.axes_series[1].name),
                left_axis=AxisSpec(field.name),
                color_scale=ColorScaleSpec(
                    label=t_name,
                    value_range=color_range,
                    colorscale=cmap_name,
                    visible=make_colorbar,
                ),
            )
        ]
    )


def spectrum_series_stacked_spec(
    spectrum_series,
    dt=None,
    t_list=None,
    dn=None,
    index_list=None,
    average=False,
    xspan=None,
    xspan_bg=None,
    scale_mode="auto",
    scale_factor=1,
    y_values="time",
    color=None,
    line_style=None,
    line_width=None,
):
    """Describe selected spectra offset along time or spectrum number."""
    t_vec = spectrum_series.t
    index_list, t_list = get_indeces_and_times(
        t_vec,
        dt=dt,
        t_list=t_list,
        dn=dn,
        index_list=index_list,
    )
    y_vectors = []
    for list_index, spectrum_index in enumerate(index_list):
        if average:
            if type(average) is int:
                range_start = max(0, spectrum_index - average)
                range_end = min(spectrum_index + average, len(t_vec))
            else:
                range_start = (
                    spectrum_index
                    if list_index == 0
                    else int((spectrum_index + index_list[list_index - 1]) / 2)
                )
                range_end = (
                    spectrum_index
                    if list_index + 1 == len(index_list)
                    else int((spectrum_index + index_list[list_index + 1]) / 2)
                )
            y_vector = spectrum_series[range_start:range_end].y_average
        else:
            y_vector = spectrum_series[spectrum_index].y
        y_vectors.append(y_vector)

    x = spectrum_series.x
    if xspan_bg:
        background_mask = np.logical_and(xspan_bg[0] < x, x < xspan_bg[-1])
        y_vectors = [
            y_vector - np.mean(y_vector[background_mask]) for y_vector in y_vectors
        ]
    if xspan:
        x_mask = np.logical_and(xspan[0] < x, x < xspan[-1])
        x = x[x_mask]
        y_vectors = [y_vector[x_mask] for y_vector in y_vectors]
    if scale_mode != "auto":
        raise ValueError(f"scale_mode='{scale_mode}' not implemented.")

    y_ranges = np.max(y_vectors, axis=1) - np.min(y_vectors, axis=1)
    if y_values == "time":
        offsets = t_list
        intervals = np.diff(t_list)
        y_label = spectrum_series.t_name
    elif y_values == "n":
        offsets = index_list
        intervals = np.diff(index_list)
        y_label = "spectrum number"
    else:
        raise ValueError(f"y_values='{y_values}' not implemented.")
    scale = min(intervals) / max(y_ranges) * scale_factor
    traces = [
        LineTrace(
            x,
            offset + y_vector * scale,
            name=str(offset),
            color=color or DEFAULT_LINE_COLORS[index % len(DEFAULT_LINE_COLORS)],
            line_style=line_style,
            line_width=line_width,
        )
        for index, (offset, y_vector) in enumerate(zip(offsets, y_vectors))
    ]
    return PlotSpec(
        [
            PanelSpec(
                traces=traces,
                x_axis=AxisSpec(spectrum_series.xseries.name),
                left_axis=AxisSpec(y_label),
            )
        ]
    )


def _ms_value_groups(
    measurement,
    mass_list,
    mass_lists,
    mol_list,
    mol_lists,
):
    """Return whether MS values are calibrated and group them by y-axis."""
    if mol_list:
        return True, [mol_list]
    if mol_lists:
        return True, list(mol_lists)
    if mass_lists:
        return False, list(mass_lists)
    return False, [mass_list or measurement.mass_list]


def _as_group_values(value, group_count, span=False):
    """Repeat one setting or split settings supplied for two groups."""
    if group_count == 1:
        return [value]
    if (
        span
        and value
        and len(value) == 2
        and all(isinstance(item, (int, float)) for item in value)
    ):
        return [value, None]
    if isinstance(value, (list, tuple)) and len(value) == group_count:
        return list(value)
    return [value] * group_count


def _ms_unit_factor(unit, quantified, measurement):
    """Return the numerical factor for a supported MS display unit."""
    if quantified:
        factor = {
            "pmol/s": 1e12,
            "nmol/s": 1e9,
            "umol/s": 1e6,
            "mmol/s": 1e3,
            "mol/s": 1,
            "pmol/s/cm^2": 1e12,
            "nmol/s/cm^2": 1e9,
            "umol/s/cm^2": 1e6,
            "mmol/s/cm^2": 1e3,
            "mol/s/cm^2": 1,
        }[unit]
        return factor / measurement.A_el if "/cm^2" in unit else factor
    return {"pA": 1e12, "nA": 1e9, "uA": 1e6, "mA": 1e3, "A": 1}[unit]


def _time_unit_factor(unit):
    """Return the factor and label for a supported time unit."""
    unit = unit or "s"
    return {
        "s": 1,
        "min": 1 / 60,
        "minutes": 1 / 60,
        "h": 1 / 3600,
        "hr": 1 / 3600,
        "hour": 1 / 3600,
        "hours": 1 / 3600,
        "d": 1 / (3600 * 24),
        "days": 1 / (3600 * 24),
    }[unit], unit
