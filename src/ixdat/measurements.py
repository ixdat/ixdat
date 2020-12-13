"""This module defines the Dataset class, the central data structure of ixdat

An ixdat Dataset is a collection of references to DataSeries with the metadata required
to combine them, i.e. "build" the combined dataset. It has a number of general methods
to visualize and analyze the combined dataset. Dataset is also the base class for a
number of technique-specific Dataset-derived classes.
"""
import json
import numpy as np
from .db import Saveable
from .data_series import PlaceHolderSeries, TimeSeries, ValueSeries
from .samples import Sample
from .lablogs import LabLog
from .plotters import ValuePlotter
from .exporters import CSVExporter
from .exceptions import BuildError, SeriesNotFoundError


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
            sample = Sample.load_or_make(sample)
        self.sample = sample
        if isinstance(lablog, str):
            lablog = LabLog.load_or_make(lablog)
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
        return [series.name for series in self.series_list]

    @property
    def value_names(self):
        return [vseries.name for vseries in self.value_series]

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

    def __getitem__(self, item):
        ss = [s for s in self.series_list if s.name == item]
        if len(ss) == 1:
            s = ss[0]
        elif len(ss) > 1:
            s = append_vseries_by_time(ss)
        elif item[-2:] in ["-t", "-x", "-v", "-y"]:
            ss = [s for s in self.series_list if s.name == item[:-2]]
            if len(ss) == 1:
                s = ss[0]
            else:
                s = append_vseries_by_time(ss)
            if item[-2:] in ["-t", "-x"]:
                s = s.tseries
        else:
            raise SeriesNotFoundError
        return s

    def get_t_and_v(self, item, tspan=None):
        vseries = self[item]
        tseries = vseries.tseries
        v = vseries.data
        t = tseries.data
        mask = np.logical_and(tspan[0] < t < tspan[-1])
        return t[mask], v[mask]

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


def append_vseries_by_time(series_list):
    name = series_list[0].name
    cls = series_list[0].__class__
    unit = series_list[0].unit
    data = np.array([])

    for s in series_list:
        if not (s.name == name and s.unit == unit and s.__class__ == cls):
            raise BuildError(f"can't append {series_list}")
        data = np.append(data, s.data)

    tseries = append_tseries([s.tseries for s in series_list])
    vseries = cls(name=name, unit=unit, data=data, tseries=tseries)
    return vseries


def append_tseries(series_list, tstamp=None):
    name = series_list[0].name
    cls = series_list[0].__class__
    unit = series_list[0].unit
    tstamp = tstamp or series_list[0].tstamp
    data = np.array([])

    for s in series_list:
        if not (s.name == name and s.unit == unit and s.__class__ == cls):
            raise BuildError(f"can't append {series_list}")
        offset = s.tstamp - tstamp
        data = np.append(data, s.data - offset)

    tseries = cls(name=name, unit=unit, data=data, tstamp=tstamp)
    return tseries


def fill_series_list(series_list, s_ids):
    series_list = series_list or []
    provided_series_ids = [s.id for s in series_list]
    for i in s_ids:
        if i not in provided_series_ids:
            series_list.append(PlaceHolderSeries(i))
    return series_list
