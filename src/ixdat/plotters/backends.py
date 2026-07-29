"""Registration of renderers for backend-neutral plot descriptions."""

from .renderers import (
    available_plotter_backends,
    get_renderer,
    register_plotter_backend,
    unregister_plotter_backend,
)


__all__ = [
    "available_plotter_backends",
    "get_renderer",
    "register_plotter_backend",
    "unregister_plotter_backend",
]
