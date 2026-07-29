"""Render backend-neutral plot descriptions."""

from html import escape

from matplotlib import pyplot as plt

from .plot_spec import ErrorBand, HeatmapTrace, LineTrace


_RENDERER_CLASSES = {}


class MatplotlibRenderer:
    """Turn a plot description into Matplotlib axes."""

    def render(self, plot_spec, axes=None):
        """Render ``plot_spec`` and return one axis or a list of axes."""
        supplied_right_axes = []
        if axes is None:
            gridspec_kw = (
                {"height_ratios": plot_spec.row_heights}
                if plot_spec.row_heights
                else None
            )
            _, panel_axes = plt.subplots(
                len(plot_spec.panels),
                1,
                squeeze=False,
                sharex=plot_spec.shared_x_axes,
                gridspec_kw=gridspec_kw,
            )
            left_axes = [row[0] for row in panel_axes]
        elif isinstance(axes, (list, tuple)):
            if len(plot_spec.panels) == 1 and plot_spec.panels[0].right_axis:
                left_axes = [axes[0]]
                supplied_right_axes = [axes[1]] if len(axes) > 1 else []
            else:
                left_axes = list(axes)
        else:
            left_axes = [axes]

        rendered_axes = []
        for panel_index, (panel, left_ax) in enumerate(zip(plot_spec.panels, left_axes)):
            if panel.right_axis and panel_index < len(supplied_right_axes):
                right_ax = supplied_right_axes[panel_index]
            else:
                right_ax = left_ax.twinx() if panel.right_axis else None
            for trace in panel.traces:
                ax = right_ax if getattr(trace, "y_axis", "left") == "right" else left_ax
                if isinstance(trace, LineTrace):
                    kwargs = {}
                    color = trace.color or _matplotlib_scaled_color(trace)
                    if color is not None:
                        kwargs["color"] = color
                    if trace.line_style is not None:
                        kwargs["linestyle"] = trace.line_style
                    if trace.line_width is not None:
                        kwargs["linewidth"] = trace.line_width
                    label = trace.name if trace.show_legend else "_nolegend_"
                    ax.plot(trace.x, trace.y, label=label, **kwargs)
                elif isinstance(trace, ErrorBand):
                    ax.fill_between(
                        trace.x,
                        trace.lower,
                        trace.upper,
                        color=trace.color,
                        alpha=trace.alpha,
                    )
                elif isinstance(trace, HeatmapTrace):
                    image = ax.pcolormesh(
                        trace.x,
                        trace.y,
                        trace.z,
                        shading="auto",
                        cmap=trace.colorscale,
                        vmin=trace.zmin,
                        vmax=trace.zmax,
                    )
                    if trace.colorbar:
                        colorbar = left_ax.figure.colorbar(image, ax=left_ax)
                        if trace.colorbar_label:
                            colorbar.set_label(trace.colorbar_label)
            _set_matplotlib_axis(left_ax, panel.x_axis, "x")
            _set_matplotlib_axis(left_ax, panel.left_axis, "y")
            if plot_spec.shared_x_axes and panel_index + 1 < len(plot_spec.panels):
                left_ax.set_xlabel("")
                left_ax.tick_params(axis="x", labelbottom=False)
            if right_ax:
                _set_matplotlib_axis(right_ax, panel.right_axis, "y")
                rendered_axes.extend([left_ax, right_ax])
            else:
                rendered_axes.append(left_ax)
            if panel.color_scale and panel.color_scale.visible:
                _add_matplotlib_colorbar(left_ax, panel.color_scale)
            if plot_spec.show_legend:
                for ax in [left_ax, right_ax]:
                    if ax and ax.get_legend_handles_labels()[0]:
                        ax.legend()
        if plot_spec.title:
            rendered_axes[0].figure.suptitle(plot_spec.title)
        return rendered_axes[0] if len(rendered_axes) == 1 else rendered_axes


class PlotlyRenderer:
    """Turn a plot description into an interactive Plotly figure."""

    def render(self, plot_spec, figure=None):
        """Render ``plot_spec`` and return a Plotly figure."""
        go, make_subplots = _import_plotly()
        panel_count = len(plot_spec.panels)
        secondary_y = [panel.right_axis is not None for panel in plot_spec.panels]
        subplot_settings = {
            "rows": panel_count,
            "cols": 1,
            "specs": [[{"secondary_y": value}] for value in secondary_y],
            "shared_xaxes": plot_spec.shared_x_axes,
            "vertical_spacing": 0.035 if panel_count > 1 else 0,
            "row_heights": plot_spec.row_heights,
        }
        if figure is None:
            figure = make_subplots(**subplot_settings)
        elif getattr(figure, "_grid_ref", None) is None:
            figure.set_subplots(**subplot_settings)
        else:
            _validate_plotly_subplot_grid(figure, secondary_y)
        first_trace_index = len(figure.data)

        for row, panel in enumerate(plot_spec.panels, start=1):
            for trace in panel.traces:
                use_secondary_y = getattr(trace, "y_axis", "left") == "right"
                if isinstance(trace, LineTrace):
                    line = {}
                    color = trace.color or _plotly_scaled_color(trace)
                    if color is not None:
                        line["color"] = _plotly_color(color)
                    if trace.line_style is not None:
                        line["dash"] = _plotly_dash(trace.line_style)
                    if trace.line_width is not None:
                        line["width"] = trace.line_width
                    plotly_trace = go.Scatter(
                        x=trace.x,
                        y=trace.y,
                        name=_escape_plotly_label(trace.name),
                        mode="lines",
                        line=line or None,
                        showlegend=trace.show_legend,
                    )
                    figure.add_trace(
                        plotly_trace,
                        row=row,
                        col=1,
                        secondary_y=use_secondary_y,
                    )
                elif isinstance(trace, ErrorBand):
                    figure.add_trace(
                        go.Scatter(
                            x=list(trace.x) + list(trace.x[::-1]),
                            y=list(trace.upper) + list(trace.lower[::-1]),
                            fill="toself",
                            fillcolor=_rgba(trace.color, trace.alpha),
                            line={"color": "rgba(0,0,0,0)"},
                            hoverinfo="skip",
                            showlegend=False,
                        ),
                        row=row,
                        col=1,
                        secondary_y=use_secondary_y,
                    )
                elif isinstance(trace, HeatmapTrace):
                    colorbar = None
                    if trace.colorbar and trace.colorbar_label:
                        colorbar = {"title": _escape_plotly_label(trace.colorbar_label)}
                    figure.add_trace(
                        go.Heatmap(
                            x=trace.x,
                            y=trace.y,
                            z=trace.z,
                            colorscale=trace.colorscale,
                            showscale=trace.colorbar,
                            colorbar=colorbar,
                            zmin=trace.zmin,
                            zmax=trace.zmax,
                        ),
                        row=row,
                        col=1,
                    )
            if panel.color_scale and panel.color_scale.visible:
                _add_plotly_colorbar(
                    figure,
                    go,
                    panel.color_scale,
                    row=row,
                )
            show_x_axis = not plot_spec.shared_x_axes or row == panel_count
            left_axis_color = _plotly_color(panel.left_axis.color or "black")
            figure.update_xaxes(
                title_text=(
                    _escape_plotly_label(panel.x_axis.label) if show_x_axis else None
                ),
                autorange="reversed" if panel.x_axis.inverted else True,
                type=_plotly_axis_type(panel.x_axis.scale),
                showticklabels=show_x_axis,
                showgrid=False,
                zeroline=False,
                ticks="outside",
                ticklen=5,
                tickwidth=1,
                tickcolor="black",
                tickfont={"color": "black"},
                ticklabelposition="outside",
                automargin=True,
                title_standoff=14,
                title_font={"color": "black"},
                row=row,
                col=1,
            )
            figure.update_yaxes(
                title_text=_escape_plotly_label(panel.left_axis.label),
                type=_plotly_axis_type(panel.left_axis.scale),
                showgrid=False,
                zeroline=False,
                ticks="outside",
                ticklen=5,
                tickwidth=1,
                tickcolor=left_axis_color,
                tickfont={"color": left_axis_color},
                ticklabelposition="outside",
                automargin=True,
                title_standoff=14,
                title_font={"color": left_axis_color},
                row=row,
                col=1,
                secondary_y=False,
            )
            if panel.right_axis:
                right_axis_color = _plotly_color(panel.right_axis.color or "black")
                figure.update_yaxes(
                    title_text=_escape_plotly_label(panel.right_axis.label),
                    type=_plotly_axis_type(panel.right_axis.scale),
                    showgrid=False,
                    zeroline=False,
                    ticks="outside",
                    ticklen=5,
                    tickwidth=1,
                    tickcolor=right_axis_color,
                    tickfont={"color": right_axis_color},
                    ticklabelposition="outside",
                    automargin=True,
                    title_standoff=14,
                    title_font={"color": right_axis_color},
                    row=row,
                    col=1,
                    secondary_y=True,
                )
        for row in range(1, panel_count + 1):
            subplot = figure.get_subplot(row, 1)
            figure.add_shape(
                type="rect",
                xref="paper",
                yref="paper",
                x0=subplot.xaxis.domain[0],
                x1=subplot.xaxis.domain[1],
                y0=subplot.yaxis.domain[0],
                y1=subplot.yaxis.domain[1],
                line={"color": "black", "width": 1},
                fillcolor="rgba(0,0,0,0)",
                layer="above",
            )
        if plot_spec.shared_x_axes:
            new_traces = figure.data[first_trace_index:]
            if new_traces:
                shared_xaxis = new_traces[-1].xaxis or "x"
                for trace in new_traces:
                    trace.xaxis = shared_xaxis
        if plot_spec.title:
            figure.update_layout(title_text=_escape_plotly_label(plot_spec.title))
        if plot_spec.show_legend is not None:
            figure.update_layout(showlegend=plot_spec.show_legend)
        matplotlib_width, matplotlib_height = plt.rcParams["figure.figsize"]
        matplotlib_dpi = plt.rcParams["figure.dpi"]
        width = figure.layout.width or round(matplotlib_width * matplotlib_dpi)
        height = figure.layout.height or round(matplotlib_height * matplotlib_dpi)
        layout = {
            "width": width,
            "height": height,
            "hovermode": "x unified" if plot_spec.shared_x_axes else "closest",
            "legend": {
                "orientation": "h",
                "x": 0,
                "xanchor": "left",
                "y": 1.03,
                "yanchor": "bottom",
            },
            "margin": {"l": 85, "r": 85, "t": 85, "b": 70},
        }
        if plot_spec.shared_x_axes:
            layout["hoversubplots"] = "axis"
        figure.update_layout(**layout)
        return figure


def register_plotter_backend(name, renderer_class, overwrite=False):
    """Register a class that renders :class:`~ixdat.plotters.plot_spec.PlotSpec`.

    Args:
        name (str): Name accepted by plotting methods as ``backend``.
        renderer_class (type): Class with a ``render(plot_spec, ...)`` method.
        overwrite (bool): Replace a renderer registered under the same name.
    """
    name = _normalize_backend_name(name)
    if not isinstance(renderer_class, type):
        raise TypeError("A plotter backend renderer must be a class.")
    if not callable(getattr(renderer_class, "render", None)):
        raise TypeError("A plotter backend renderer must define render().")
    if name in ("matplotlib", "plotly"):
        raise ValueError(f"'{name}' is a built-in plotter backend.")
    if name in _RENDERER_CLASSES and not overwrite:
        raise ValueError(
            f"Plotter backend '{name}' is registered. "
            "Pass overwrite=True to replace it."
        )
    _RENDERER_CLASSES[name] = renderer_class
    return renderer_class


def unregister_plotter_backend(name):
    """Remove and return a registered renderer class."""
    name = _normalize_backend_name(name)
    if name in ("matplotlib", "plotly"):
        raise ValueError(f"'{name}' is a built-in plotter backend.")
    return _RENDERER_CLASSES.pop(name)


def available_plotter_backends():
    """Return backend names accepted by plotting methods."""
    return tuple(sorted(_RENDERER_CLASSES))


def get_renderer(backend=None):
    """Return the renderer named by ``backend``."""
    name = "matplotlib" if backend is None else _normalize_backend_name(backend)
    try:
        renderer_class = _RENDERER_CLASSES[name]
    except KeyError:
        available = ", ".join(available_plotter_backends())
        raise ValueError(
            f"Unknown plot backend '{backend}'. Available backends: {available}."
        ) from None
    return renderer_class()


def _is_matplotlib_backend(backend):
    """Return whether a backend selection uses ixdat's direct Matplotlib path."""
    return backend is None or _normalize_backend_name(backend) == "matplotlib"


def _set_matplotlib_axis(ax, axis_spec, direction):
    """Apply an axis description to one Matplotlib axis."""
    if direction == "x":
        if axis_spec.label:
            ax.set_xlabel(axis_spec.label)
        if axis_spec.color:
            ax.xaxis.label.set_color(axis_spec.color)
            ax.tick_params(axis="x", colors=axis_spec.color)
        ax.set_xscale(axis_spec.scale)
        if axis_spec.inverted:
            ax.invert_xaxis()
    else:
        if axis_spec.label:
            ax.set_ylabel(axis_spec.label)
        if axis_spec.color:
            ax.yaxis.label.set_color(axis_spec.color)
            ax.tick_params(axis="y", colors=axis_spec.color)
        ax.set_yscale(axis_spec.scale)


def _plotly_axis_type(scale):
    """Translate a common scale name to Plotly."""
    return "log" if scale == "log" else "linear"


def _validate_plotly_subplot_grid(figure, secondary_y):
    """Check that a supplied Plotly subplot grid can hold the plot description."""
    grid = figure._grid_ref
    compatible_rows = len(grid) == len(secondary_y) and all(
        len(row) == 1 for row in grid
    )
    compatible_axes = compatible_rows and all(
        grid[row_index][0] is not None
        and (not needs_secondary or len(grid[row_index][0]) > 1)
        for row_index, needs_secondary in enumerate(secondary_y)
    )
    if not compatible_axes:
        raise ValueError(
            "The supplied Plotly figure has an incompatible subplot grid. "
            "Pass a plain plotly.graph_objects.Figure or a figure returned by "
            "the same ixdat plot method."
        )


def _plotly_dash(line_style):
    """Translate common Matplotlib line styles to Plotly."""
    return {
        "-": "solid",
        "--": "dash",
        ":": "dot",
        "-.": "dashdot",
        "solid": "solid",
        "dashed": "dash",
        "dotted": "dot",
        "dashdot": "dashdot",
    }.get(line_style, line_style)


def _plotly_color(color):
    """Translate Matplotlib's one-letter colors to CSS color names."""
    color = {
        "k": "black",
        "r": "red",
        "b": "blue",
        "g": "green",
        "c": "cyan",
        "m": "magenta",
        "y": "yellow",
        "w": "white",
    }.get(color, color)
    try:
        from matplotlib.colors import to_hex

        return to_hex(color)
    except ValueError:
        return color


def _matplotlib_scaled_color(trace):
    """Return a Matplotlib color sampled from a trace's color scale."""
    if trace.color_value is None or not trace.color_range:
        return None
    from matplotlib import cm, colors

    color_map = cm.get_cmap(trace.color_scale)
    normalizer = colors.Normalize(*trace.color_range)
    return color_map(normalizer(trace.color_value))


def _plotly_scaled_color(trace):
    """Return a Plotly color sampled from a trace's color scale."""
    if trace.color_value is None or not trace.color_range:
        return None
    from plotly.colors import sample_colorscale

    lower, upper = trace.color_range
    position = (trace.color_value - lower) / (upper - lower) if upper != lower else 0
    return sample_colorscale(trace.color_scale, [position])[0]


def _add_matplotlib_colorbar(ax, color_scale):
    """Add the continuous line-color scale to a Matplotlib axis."""
    from matplotlib import cm, colors

    normalizer = colors.Normalize(*color_scale.value_range)
    colorbar = ax.figure.colorbar(
        cm.ScalarMappable(norm=normalizer, cmap=color_scale.colorscale),
        ax=ax,
    )
    colorbar.set_label(color_scale.label)


def _add_plotly_colorbar(figure, go, color_scale, row):
    """Add the continuous line-color scale to a Plotly panel."""
    lower, upper = color_scale.value_range
    figure.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker={
                "color": [lower],
                "cmin": lower,
                "cmax": upper,
                "colorscale": color_scale.colorscale,
                "showscale": True,
                "colorbar": {"title": _escape_plotly_label(color_scale.label)},
            },
            hoverinfo="skip",
            showlegend=False,
        ),
        row=row,
        col=1,
    )


def _rgba(color, alpha):
    """Return a Plotly fill color."""
    color = _plotly_color(color or "black")
    try:
        from matplotlib.colors import to_rgb

        red, green, blue = to_rgb(color)
    except ValueError:
        return color
    return (
        f"rgba({round(red * 255)}, {round(green * 255)}, "
        f"{round(blue * 255)}, {alpha})"
    )


def _import_plotly():
    """Import Plotly or raise an installation-focused error."""
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError as e:
        raise ImportError(
            "The Plotly backend requires Plotly. "
            'Install it with `pip install "ixdat[plotly]"`.'
        ) from e
    return go, make_subplots


def _escape_plotly_label(label):
    """Escape label text that Plotly interprets as HTML."""
    return None if label is None else escape(str(label))


def _normalize_backend_name(name):
    """Return the normalized form of a backend name."""
    if not isinstance(name, str):
        raise TypeError("A plotter backend name must be a string.")
    name = name.strip().lower()
    if not name:
        raise ValueError("A plotter backend name cannot be empty.")
    return name


_RENDERER_CLASSES.update(
    {
        "matplotlib": MatplotlibRenderer,
        "plotly": PlotlyRenderer,
    }
)
