from .backend_base import BackendBase


class MemoryBackend(BackendBase):
    """A backend that assigns id's and keeps track of objects for later access.

    This means that a Saveable object can make a serializable representation of itself
    that includes the identity in the memory backend of objects that it references, and
    a new Saveable object made from this representation can find the original objects.
    See db.Saveable.as_dict(), which ensures that "child objects" are in memory.
    """

    backend_type = "memory"
    address = "here"

    def __init__(self):
        """Initialize the backend with dict for {table_name (str): id_counter (int)}"""
        super().__init__()
        self.objects = {}

    def get_next_available_id(self, table_name, obj=None):
        """Return the id counter for table_name, starting with 1."""
        i = super().get_next_available_id(table_name)
        if table_name not in self.objects:
            self.objects[table_name] = {}
        self.objects[table_name][i] = obj
        return i

    def get(self, cls, i):
        """Return an object of a specified class and id by looking up in memory"""
        return self.objects[cls.table_name][i]

    def remember(self, obj):
        """Save the object into memory, and change its backend to this memory backend.

        FIXME: This should probably just be called save() since it behaves the same.
        """
        if obj.backend is self:
            return obj.id
        table_name = obj.table_name
        i = self.get_next_available_id(table_name, obj)
        obj.set_id(i)
        obj.set_backend(self)
        return i
