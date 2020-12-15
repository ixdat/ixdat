from . import BackendBase


class MemoryBackend(BackendBase):
    def __init__(self):
        self.next_available_ids = {}

    def get_next_available_id(self, table_name):
        if table_name in self.next_available_ids:
            i = self.next_available_ids[table_name]
            self.next_available_ids[table_name] += 1
        else:
            i = 1
            self.next_available_ids[table_name] = 2
        return i
