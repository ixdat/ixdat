"""Select renderers through adapters around ixdat's Matplotlib plotters."""

from functools import wraps
import inspect
import warnings

from .renderers import (
    _is_matplotlib_backend,
    available_plotter_backends,
    get_renderer,
    register_plotter_backend,
    unregister_plotter_backend,
)

_PLOTTER_ADAPTERS = {}


class _AdapterArguments(dict):
    """Carry bound method arguments and the names supplied by the caller."""

    def __init__(self, arguments, supplied):
        super().__init__(arguments)
        self.supplied = supplied


class PlotterBackendWarning(UserWarning):
    """Warn that a plot uses Matplotlib after an unsupported backend request."""


def register_plotter_adapter(plotter_class, method_name, adapter, overwrite=False):
    """Register a plot-description adapter for one plotter method.

    The adapter receives the plotted object and the arguments bound to the original
    Matplotlib method. It returns a backend-neutral ``PlotSpec``.
    """
    if not isinstance(plotter_class, type):
        raise TypeError("A plotter adapter owner must be a class.")
    if not isinstance(method_name, str) or not method_name:
        raise TypeError("A plotter adapter method name must be a non-empty string.")
    if not callable(adapter):
        raise TypeError("A plotter adapter must be callable.")
    key = (plotter_class, method_name)
    if key in _PLOTTER_ADAPTERS and not overwrite:
        raise ValueError(
            f"An adapter for {plotter_class.__name__}.{method_name}() is registered. "
            "Pass overwrite=True to replace it."
        )
    _PLOTTER_ADAPTERS[key] = adapter
    return adapter


def bind_plotter(plotter, owner):
    """Add backend dispatch to the plotting methods of one bound plotter tree."""
    _bind_plotter(plotter, owner, visited=set())
    return plotter


def _bind_plotter(plotter, owner, visited):
    """Bind one plotter and any plotters to which it delegates."""
    if id(plotter) in visited:
        return
    visited.add(id(plotter))

    for attribute_name, value in vars(plotter).items():
        if attribute_name.endswith("_plotter") and value is not None:
            _bind_plotter(value, owner, visited)

    for method_name in dir(plotter):
        if not _is_plot_method_name(method_name):
            continue
        method = getattr(plotter, method_name)
        if not callable(method):
            continue
        static_attribute = inspect.getattr_static(type(plotter), method_name, None)
        if isinstance(static_attribute, property):
            continue
        original_method = getattr(method, "_ixdat_original_plot_method", method)
        setattr(
            plotter,
            method_name,
            _backend_dispatcher(
                plotter=plotter,
                owner=owner,
                method_name=method_name,
                original_method=original_method,
            ),
        )


def _backend_dispatcher(plotter, owner, method_name, original_method):
    """Return a bound plotting function with renderer selection."""

    @wraps(original_method)
    def dispatch(*args, backend=None, figure=None, **kwargs):
        if _is_matplotlib_backend(backend):
            return original_method(*args, **kwargs)

        renderer = get_renderer(backend)
        adapter = _get_plotter_adapter(type(plotter), method_name)
        if adapter is None:
            warnings.warn(
                f"No adapter for backend '{backend}' is registered for "
                f"{type(plotter).__name__}.{method_name}(). Using Matplotlib.",
                PlotterBackendWarning,
                stacklevel=2,
            )
            return original_method(*args, **kwargs)

        arguments = inspect.signature(original_method).bind(*args, **kwargs)
        supplied_arguments = frozenset(arguments.arguments)
        arguments.apply_defaults()
        adapter_arguments = _AdapterArguments(
            arguments.arguments,
            supplied=supplied_arguments,
        )
        plot_spec = adapter(owner, adapter_arguments)
        return renderer.render(plot_spec, figure=figure)

    dispatch._ixdat_original_plot_method = original_method
    dispatch.__signature__ = _backend_signature(original_method)
    dispatch.__doc__ = (
        (original_method.__doc__ or "")
        + "\n\n"
        + "Renderer options:\n"
        + '    backend (str): Renderer name, such as "matplotlib" or "plotly".\n'
        + "    figure: Plotly figure to extend when using the Plotly backend.\n"
    )
    return dispatch


def _get_plotter_adapter(plotter_class, method_name):
    """Return the closest adapter registered for a plotter class."""
    for candidate in plotter_class.__mro__:
        adapter = _PLOTTER_ADAPTERS.get((candidate, method_name))
        if adapter is not None:
            return adapter
    return None


def _backend_signature(method):
    """Expose renderer options on a bound plotting method."""
    signature = inspect.signature(method)
    parameters = list(signature.parameters.values())
    insertion_index = next(
        (
            index
            for index, parameter in enumerate(parameters)
            if parameter.kind == inspect.Parameter.VAR_KEYWORD
        ),
        len(parameters),
    )
    renderer_parameters = [
        inspect.Parameter(name, inspect.Parameter.KEYWORD_ONLY, default=None)
        for name in ("backend", "figure")
        if name not in signature.parameters
    ]
    parameters[insertion_index:insertion_index] = renderer_parameters
    return signature.replace(parameters=parameters)


def _is_plot_method_name(name):
    """Return whether an attribute is a public plotting entry point."""
    return not name.startswith("_") and (name.startswith("plot") or name == "heat_plot")


__all__ = [
    "PlotterBackendWarning",
    "available_plotter_backends",
    "bind_plotter",
    "get_renderer",
    "register_plotter_adapter",
    "register_plotter_backend",
    "unregister_plotter_backend",
]
