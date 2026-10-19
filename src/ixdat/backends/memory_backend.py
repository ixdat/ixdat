from .backend_base import BackendBase
from ..exceptions import DataBaseError


class MemoryBackend(BackendBase):
    """A backend that assigns id's and keeps track of objects for later access.

    This means that a Savable object can make a serializable representation of itself
    that includes the short_identity in the memory backend of objects that it
    references, and a new Savable object made from this representation can find the
    original objects. See db.Saveable.as_dict(), which ensures that related objects
    are in memory.
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
        return i

    def get(self, cls, i):
        """Return an object of a specified class and id by looking up in memory"""
        return self.objects[cls.table_name][i]

    def load(self, cls, name):
        """Return the most recently saved object of cls with the given name"""
        objects_by_id = self.objects.get(cls.table_name, {})
        # id's count up as objects are saved, so the highest one is the newest:
        for i in sorted(objects_by_id, reverse=True):
            if objects_by_id[i].name == name:
                return objects_by_id[i]
        raise DataBaseError(
            f"{self} has no object named '{name}' in table '{cls.table_name}'"
        )

    def save(self, obj):
        """Save the object into memory, and change its backend to this backend."""
        if obj.backend is self:
            return obj.id
        table_name = obj.table_name
        i = self.get_next_available_id(table_name, obj)
        self.objects[obj.table_name][i] = obj
        obj.set_id(i)
        obj.set_backend(self)
        return i
