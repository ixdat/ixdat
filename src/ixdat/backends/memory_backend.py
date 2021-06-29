from .backend_base import BackendBase


class MemoryBackend(BackendBase):
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

    def get(self, table_name, i):
        return self.objects[table_name][i]
