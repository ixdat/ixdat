"""Import techniques and build the technique_classes dictionary for direct import

Constants:
    DATABASE_BACKENDS (dict): Dictionary of {backend_name: BackendClass} where
        backend_name is the name of the backend (like "directory") and BackendClass
        is the backend class (inheriting from Backend) for saving and loading things.
"""

from .memory_backend import BackendBase, MemoryBackend
from .directory_backend import DirBackend


DATABASE_BACKENDS = {
    "None": BackendBase,
    "memory": MemoryBackend,
    "directory": DirBackend,
}
