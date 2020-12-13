from .db import Saveable


class Sample(Saveable):
    table_name = "sample"
    column_attrs = {"name": "name"}

    def __init__(self, name):
        super().__init__()
        self.name = name

    @classmethod
    def load_or_make(cls, name):
        return cls(name)
