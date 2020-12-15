"""This module implements the simplest backend, which just counts rows in memory"""
from . import BackendBase


class MemoryBackend(BackendBase):
    """Simplest possible backend. No saving or loading, just counting."""

    def __init__(self):
        """Initialize the backend with dict for {table_name (str): id_counter (int)}"""
        self.next_available_ids = {}

    def get_next_available_id(self, table_name):
        """Return the id counter for table_name, starting with 1."""
        if table_name in self.next_available_ids:
            i = self.next_available_ids[table_name]
            self.next_available_ids[table_name] += 1
        else:
            i = 1
            self.next_available_ids[table_name] = 2
        return i
