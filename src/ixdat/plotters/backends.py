"""Connect ixdat plot methods to optional figure libraries.

An ixdat data object keeps its normal Matplotlib plotter. :func:`bind_plotter`
adds ``backend`` and ``figure`` options to that plotter's public methods on the
specific plotter instance. A call for another figure library follows a registered
adapter, which describes the requested plot as a :class:`~ixdat.plotters.PlotSpec`.
The selected renderer then draws that description.

The Matplotlib plotter class keeps its original methods and signature. A new reader
can assign its Matplotlib plotter first and add an adapter in a later change.
"""

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
    """Map each original method parameter to the value used for this plot call.

    ``supplied`` records which parameter names the caller provided. The mapping also
    contains default values from the Matplotlib method signature.
    """

    def __init__(self, arguments, supplied):
        super().__init__(arguments)
        self.supplied = supplied


class PlotterBackendWarning(UserWarning):
    """Tell the user that ixdat used Matplotlib because an adapter was unavailable."""


def register_plotter_adapter(
    plotter_class,
    method_name,
    adapter=None,
    overwrite=False,
):
    """Connect one Matplotlib plot method to ``PlotSpec`` creation.

    A plotter adapter is a function that describes one plot call without creating
    Matplotlib or Plotly objects. ixdat passes it the measurement or spectrum being
    plotted and a mapping of the original method's arguments. The adapter returns a
    :class:`~ixdat.plotters.PlotSpec`, which any registered renderer can draw.

    The function can be used as a decorator::

        from ixdat.plotters import register_plotter_adapter
        from ixdat.plotters.plot_spec import value_measurement_spec
        from my_ixdat_extension import MyPlotter

        @register_plotter_adapter(MyPlotter, "plot_measurement")
        def my_plot_adapter(owner, args):
            measurement = args["measurement"] or owner
            return value_measurement_spec(
                measurement,
                v_list=args["v_list"],
                tspan=args["tspan"],
            )

    Import the module containing this registration when ixdat starts. Built-in
    adapters live in :mod:`ixdat.plotters.plot_adapters`.

    Args:
        plotter_class (type): Matplotlib plotter class that owns the method.
        method_name (str): Public plot method, such as ``"plot_measurement"``.
        adapter (callable): Function that returns a ``PlotSpec``. Omitting it returns
            a decorator.
        overwrite (bool): Replace a registration for the same class and method.
    """
    if not isinstance(plotter_class, type):
        raise TypeError("A plotter adapter owner must be a class.")
    if not isinstance(method_name, str) or not method_name:
        raise TypeError("A plotter adapter method name must be a non-empty string.")
    if adapter is None:
        return lambda adapter_function: register_plotter_adapter(
            plotter_class,
            method_name,
            adapter_function,
            overwrite=overwrite,
        )
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
    """Connect a data object's plot methods to registered renderers.

    ``owner`` is the measurement or spectrum that supplies the data. ``plotter`` is
    the object's Matplotlib plotter. ixdat wraps public methods such as
    ``plot_measurement()`` on this plotter instance and follows these routes:

    - no ``backend``, or ``backend="matplotlib"``: call the original method;
    - a registered backend with an adapter: build a ``PlotSpec`` and render it;
    - a registered backend without an adapter: warn and call the original method.

    Composite plotters often contain helpers such as ``ms_plotter``. This function
    connects those child plotters to the same owner as well. Plotter classes and
    reader classes need no renderer-specific methods.
    """
    _bind_plotter(plotter, owner, visited=set())
    return plotter


def _bind_plotter(plotter, owner, visited):
    """Connect one plotter and each child plotter to the same data object."""
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
            _renderer_routing_method(
                plotter=plotter,
                owner=owner,
                method_name=method_name,
                original_method=original_method,
            ),
        )


def _renderer_routing_method(plotter, owner, method_name, original_method):
    """Wrap one plot method so each backend request follows the correct route."""

    @wraps(original_method)
    def route_plot_call(*args, backend=None, figure=None, **kwargs):
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

    route_plot_call._ixdat_original_plot_method = original_method
    route_plot_call.__signature__ = _backend_signature(original_method)
    route_plot_call.__doc__ = (
        (original_method.__doc__ or "")
        + "\n\n"
        + "Renderer options:\n"
        + '    backend (str): Renderer name, such as "matplotlib" or "plotly".\n'
        + "    figure: Plotly figure to extend when using the Plotly backend.\n"
    )
    return route_plot_call


def _get_plotter_adapter(plotter_class, method_name):
    """Find the adapter registered for this plotter or its nearest parent class."""
    for candidate in plotter_class.__mro__:
        adapter = _PLOTTER_ADAPTERS.get((candidate, method_name))
        if adapter is not None:
            return adapter
    return None


def _backend_signature(method):
    """Add ``backend`` and ``figure`` to the displayed method signature."""
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
    """Identify public methods that a user can call to create a plot."""
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
