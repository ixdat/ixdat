"""This module defines the Dataset class, the central data structure of ixdat

An ixdat Dataset is a collection of references to DataSeries with the metadata required
to combine them, i.e. "build" the combined dataset. It has a number of general methods
to visualize and analyze the combined dataset. Dataset is also the base class for a
number of technique-specific Dataset-derived classes.
"""
from .db import Saveable
from .data_series import DataSeries


class Measurement(Saveable):
    """The Measurement class"""

    table_name = "measurement"
    column_attrs = {
        "id": "i", "name": "name", "metadata": "metadata_json", "s_ids": "s_ids",
    }

    def __init__(self, i, name, metadata, s_ids=None, series=None):
        """initialize a measurement"""
        super().__init__()
        self.id = i
        self.name = name
        self.metadata = metadata
        self._s_ids = s_ids
        self._series = series

