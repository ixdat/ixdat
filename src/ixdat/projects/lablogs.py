from ixdat.db import Saveable


class LabLog(Saveable):
    """TODO: flush out this class"""

    table_name = "lablog"
    column_attrs = {"name": "name"}

    def __init__(self, name, metadata=None, notes=None):
        """Initiate the lablog with its name, metadata, and notes"""
        super().__init__()
        self.name = name
        self.metadata = metadata or {}
        self.notes = notes or ""

    @classmethod
    def load_or_make(cls, name):
        """Load a lab log or make a new one if it hasn't been saved before"""
        return cls(name)
