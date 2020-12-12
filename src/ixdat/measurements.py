"""This module defines the Dataset class, the central data structure of ixdat

An ixdat Dataset is a collection of references to DataSeries with the metadata required
to combine them, i.e. "build" the combined dataset. It has a number of general methods
to visualize and analyze the combined dataset. Dataset is also the base class for a
number of technique-specific Dataset-derived classes.
"""
import json
from .db import Saveable
from .data_series import PlaceHolderSeries, DataSeries, TimeSeries, ValueSeries, Field
from .samples import Sample
from .lablogs import LabLog
from .plotters import ValuePlotter
from .exporters import CSVExporter


class Measurement(Saveable):
    """The Measurement class"""

    table_name = "measurement"
    column_attrs = {
        "id": "i",
        "name": "name",
        "metadata": "metadata_json",
        "sample": "sample_name",
    }
    extra_linkers = {"measurement_series": ("data_series", {"s_ids": "s_ids"})}

    def __init__(
        self,
        i,
        name,
        metadata,
        s_ids=None,
        series_list=None,
        plotter=None,
        exporter=None,
        sample=None,
        lablog=None,
    ):
        """initialize a measurement"""
        super().__init__()
        self.id = i
        self.name = name
        if isinstance(metadata, str):
            metadata = json.loads(metadata)
        self.metadata = metadata
        self.plotter = plotter or ValuePlotter(measurement=self)
        self.exporter = exporter or CSVExporter(measurement=self)
        if isinstance(sample, str):
            sample = Sample.open_or_make(sample)
        self.sample = sample
        if isinstance(lablog, str):
            lablog = LabLog.open_or_make(lablog)
        self.lablog = lablog
        self._series_list = fill_series_list(series_list, s_ids)

    @property
    def metadata_json(self):
        return json.dumps(self.metadata)

    @property
    def sample_name(self):
        return self.sample.name

    @property
    def series_list(self):
        for i, s in enumerate(self._series_list):
            if isinstance(s, PlaceHolderSeries):
                self._series_list[i] = s.get_series()
        return self._series_list

    @property
    def s_ids(self):
        return [series.id for series in self._series_list]

    @property
    def series_dict(self):
        return {s.id: s for s in self.series_list}

    @property
    def series_names(self):
        return [series.name for series in self.series_names]

    @property
    def value_series(self):
        return [
            series.name
            for series in self.series_list
            if isinstance(series, ValueSeries)
        ]

    @property
    def time_series(self):
        return [
            series.name for series in self.series_list if isinstance(series, TimeSeries)
        ]

    def __getitem__(self, name):
        ss = [s for s in self.series_list if s.name == name]
        if len(ss) == 1:
            return ss[0]

    @property
    def data_cols(self):
        return set([s.name for s in (self.value_series + self.time_series)])

    def plot(self, plotter=None, *args, **kwargs):
        if plotter:
            return plotter.plot_measurement(self, *args, **kwargs)
        return self.plotter.plot(*args, **kwargs)

    def export(self, exporter=None, *args, **kwargs):
        if exporter:
            return exporter.export_measurement(self, *args, **kwargs)
        return self.exporter.export(*args, **kwargs)


def fill_series_list(series_list, s_ids):
    series_list = series_list or []
    provided_series_ids = [s.id for s in series_list]
    for i in s_ids:
        if i not in provided_series_ids:
            series_list.append(PlaceHolderSeries(i))
    return series_list
