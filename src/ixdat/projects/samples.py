"""The module implements the sample class"""

from ixdat.db import Saveable


class Sample(Saveable):
    """TODO: flush out this class"""

    table_name = "sample"
    column_attrs = {"name": "name"}

    def __init__(self, name):
        """Initate the sample with its name"""
        super().__init__()
        self.name = name

    @classmethod
    def load_or_make(cls, name):
        """Load the sample or make a new one if it's not already saved."""
        return cls(name)
