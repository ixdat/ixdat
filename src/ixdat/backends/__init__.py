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
# FIXME: should automate that all initiated backends get added to database_backends?
database_backends = {
    "none": BackendBase(),  # Just assigns id's but doesn't keep track.
    "memory": MemoryBackend(),  # Keeps track so child objects can be passed around
    "directory": DirBackend(),  # Saves json files, stand-in for SQL, mainly for testing
}
