"""This module implements the simplest backend, which just counts rows in memory"""


class BackendBase:
    """Base class listing the functions that must be implemented in a database backend.

    A backend defines where and how all Saveable objects (inheriting from the Saveable
    class) save their data. This is a key part of the seamless interoperability of
    ixdat classes and experimental databases. Each Saveable class roughly corresponds
    to a table, and the save() and get() functions correspond to inserting and
    selecting from a database table.

    Backends inheriting from this base class are in modules of the ixdat.backends folder

    Saveable objects will by default be initiated with this base (backend="none")
    backend. Objects which ixdat should keep track of for use in multiple places in
    a workflow should instead be initiated with MemoryBackend (backend="memory").
    Other backends are imply an object is saved to or loaded form a database including
    ixdat's directory backend (backend="directory")
    """

    backend_type = "none"

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
        """Save a Saveable object and return its id. Must be implemented."""
        raise NotImplementedError

    def get(self, cls, i):
        """Load the object with id=i of a Saveable class. Must be implemented."""
        raise NotImplementedError

    def load_obj_data(self, obj):
        """Load and return the 'data' of a saveable object. Must be implemented."""
        raise NotImplementedError

    @property
    def name(self):
        return self.backend_type
