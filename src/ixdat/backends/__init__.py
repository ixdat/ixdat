"""Import techniques and build the technique_classes dictionary for direct import

Constants:
    DATABASE_BACKENDS (dict): Dictionary of {backend_name: BackendClass} where
        backend_name is the name of the backend (like "directory") and BackendClass
        is the backend class (inheriting from Backend) for saving and loading things.
"""

from .backend_base import BackendBase
from .memory_backend import MemoryBackend
from .directory_backend import DirBackend


BACKEND_CLASSES = {
    "none": BackendBase,
    "memory": MemoryBackend,
    "directory": DirBackend,
}

database_backends = {
    "none": BackendBase(),
    "memory": MemoryBackend(),
}  # This will store the initiated
