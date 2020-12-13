from .db import Saveable


class LabLog(Saveable):
    table_name = "lablog"
    column_attrs = {"name": "name"}

    def __init__(self, name):
        super().__init__()
        self.name = name

    @classmethod
    def load_or_make(cls, name):
        return cls(name)
