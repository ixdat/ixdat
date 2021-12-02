"""This module implements the simplest backend, which just counts objects."""


class BackendBase:
    """Base class listing the functions that must be implemented in a database backend.

    A backend defines where and how all Savable objects (inheriting from the Savable
    class) save their data. This is a key part of the seamless interoperability of
    ixdat classes and experimental databases. Each Savable class roughly corresponds
    to a table, and the save() and get() functions correspond to inserting and
    selecting from a database table.

    Backends inheriting from this base class are in modules of the ixdat.backends folder

    Savable objects will by default be initiated with this base (backend="none")
    backend. Objects which ixdat should keep track of for use in multiple places in
    a workflow should instead be initiated with MemoryBackend (backend="memory").
    Other backends imply an object is saved to or loaded form a database including
    ixdat's directory backend (backend="directory")
    """

    backend_type = "none"
    address = "none"
    """A location uniquely identifying the backend together with its type"""

    def __init__(self):
        """Initialize the backend with dict for {table_name (str): id_counter (int)}"""
        self.next_available_ids = {}

    def get_next_available_id(self, table_name, obj=None):
        """Return the id counter for table_name, starting with 1."""
        if table_name in self.next_available_ids:
            i = self.next_available_ids[table_name]
            self.next_available_ids[table_name] += 1
        else:
            i = 1
            self.next_available_ids[table_name] = 2
        return i

    def save(self, obj):
        """Save a Savable object and return its id. Must be implemented."""
        raise NotImplementedError

    def get(self, cls, i):
        """Load the object with id=i of a Savable class. Must be implemented."""
        raise NotImplementedError

    def load_obj_data(self, obj):
        """Load and return the 'data' of a saveable object. Must be implemented."""
        raise NotImplementedError

    @property
    def name(self):
        return f"{self.__class__.__name__}({self.backend_type}, address={self.address})"

    def __repr__(self):
        return self.name
