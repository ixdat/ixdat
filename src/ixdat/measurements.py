"""This module defines the Dataset class, the central data structure of ixdat

An ixdat Dataset is a collection of references to DataSeries with the metadata required
to combine them, i.e. "build" the combined dataset. It has a number of general methods
to visualize and analyze the combined dataset. Dataset is also the base class for a
number of technique-specific Dataset-derived classes.
"""
import json
import numpy as np
from .db import Saveable, PlaceHolderObject
from .data_series import DataSeries, TimeSeries, ValueSeries
from .samples import Sample
from .lablogs import LabLog
from .exporters import CSVExporter
from .exceptions import BuildError, SeriesNotFoundError


class Measurement(Saveable):
    """The Measurement class"""

    table_name = "measurement"
    column_attrs = {
        "id": "id",
        "name": "name",
        "technique": "technique",
        "metadata": "metadata_json_string",
        "sample": "sample_name",
        "tstamp": "tstamp",
    }
    extra_linkers = {
        "measurement_series": ("data_series", {"s_ids": "s_ids"}),
        "component_measurements": ("measurements", {"m_ids": "m_ids"}),
    }

    def __init__(
        self,
        name,
        technique=None,
        metadata=None,
        s_ids=None,
        series_list=None,
        m_ids=None,
        component_measurements=None,
        reader=None,
        plotter=None,
        exporter=None,
        sample=None,
        lablog=None,
        tstamp=None,
    ):
        """initialize a measurement

        Args:
            name (str): The name of the measurement
            TODO: Decide if metadata needs the json string option.
                See: https://github.com/ixdat/ixdat/pull/1#discussion_r546436991
            metadata (dict or json string): Free-form measurement metadata
            technique (str): The measurement technique
            s_ids (list of int): The id's of the measurement's DataSeries, if
                to be loaded (instead of given directly in series_list)
            series_list (list of DataSeries): The measurement's DataSeries
            m_ids (list of int): The id's of the component measurements, if to be
                loaded. None unless this is a combined measurement (typically
                corresponding to more than one file).
            component_measurements (list of Measurements): The measurements of which
                this measurement is a combination
            reader (Reader): The file reader (None unless read from a file)
            plotter (Plotter): The visualization tool for the measurement
            exporter (Exporter): The exporting tool for the measurement
            sample (Sample): The sample being measured
            lablog (LabLog): The log entry with e.g. notes taken during the measurement
            tstamp (float): The nominal starting time of the measurement, used for
                data selection, visualization, and exporting.
        """
        super().__init__()
        self.name = name
        self.technique = technique
        if isinstance(metadata, str):
            metadata = json.loads(metadata)
        self.metadata = metadata or {}
        self.reader = reader
        self.plotter = plotter
        self.exporter = exporter or CSVExporter(measurement=self)
        if isinstance(sample, str):
            sample = Sample.load_or_make(sample)
        self.sample = sample
        if isinstance(lablog, str):
            lablog = LabLog.load_or_make(lablog)
        self.lablog = lablog
        self._series_list = fill_object_list(series_list, s_ids, cls=DataSeries)
        self._component_measurements = fill_object_list(
            component_measurements, m_ids, cls=Measurement
        )
        self.tstamp = tstamp

    @classmethod
    def from_dict(cls, obj_as_dict):
        """Return an object of the measurement class of the right technique

        Args:
              obj_as_dict (dict): The full serializaiton (rows from table and aux
                tables) of the measurement. obj_as_dict["technique"] specifies the
                technique class to use, from TECHNIQUE_CLASSES
        """
        # TODO: see if there isn't a way to put the import at the top of the module.
        #    see: https://github.com/ixdat/ixdat/pull/1#discussion_r546437410
        from .techniques import TECHNIQUE_CLASSES

        if obj_as_dict["technique"] in TECHNIQUE_CLASSES:
            technique_class = TECHNIQUE_CLASSES[obj_as_dict["technique"]]
        else:
            technique_class = cls
        return technique_class(**obj_as_dict)

    @classmethod
    def read(cls, path_to_file, reader, **kwargs):
        """Return a Measurement object from parsing a file with the specified reader"""
        if isinstance(reader, str):
            # TODO: see if there isn't a way to put the import at the top of the module.
            #    see: https://github.com/ixdat/ixdat/pull/1#discussion_r546437471
            from .readers import READER_CLASSES

            reader = READER_CLASSES[reader]()
        return reader.read(path_to_file, **kwargs)

    @property
    def metadata_json_string(self):
        """Measurement metadata as a JSON-formatted string"""
        return json.dumps(self.metadata, indent=4)

    @property
    def sample_name(self):
        """Name of the sample on which the measurement was conducted"""
        if self.sample:
            return self.sample.name

    @property
    def series_list(self):
        """List of the DataSeries containing the measurement's data"""
        for i, s in enumerate(self._series_list):
            if isinstance(s, PlaceHolderObject):
                self._series_list[i] = s.get_object()
        return self._series_list

    @property
    def data_objects(self):
        """This is what the DB backend knows to save separately, here the series"""
        return self.series_list

    @property
    def component_measurements(self):
        """List of the component measurements of which this measurement is a combination

        For a pure measurement (not a measurement set), this is None.
        """
        if not self._component_measurements:
            return None
        for i, m in enumerate(self._component_measurements):
            if isinstance(m, PlaceHolderObject):
                self._component_measurements[i] = m.get_object()
        return self._component_measurements

    @property
    def s_ids(self):
        """List of the id's of the measurement's DataSeries"""
        return [series.id for series in self._series_list]

    @property
    def m_ids(self):
        """List of the id's of a combined measurement's component measurements"""
        if not self._component_measurements:
            return None
        return [m.id for m in self._component_measurements]

    @property
    def series_dict(self):
        """Dictionary mapping the id's of the measurement's series to the DataSeries"""
        return {s.id: s for s in self.series_list}

    @property
    def series_names(self):
        """List of the names of the series in the measurement"""
        return [series.name for series in self.series_list]

    @property
    def value_names(self):
        """List of the names of the VSeries in the measurement's DataSeries"""
        return [vseries.name for vseries in self.value_series]

    @property
    def value_series(self):
        """List of the VSeries in the measurement's DataSeries"""
        return [
            series for series in self.series_list if isinstance(series, ValueSeries)
        ]

    @property
    def time_series(self):
        """List of the TSeries in the measurement's DataSeries. NOT timeshifted!"""
        return [
            series.name for series in self.series_list if isinstance(series, TimeSeries)
        ]

    def __getitem__(self, item):
        """Return the built measurement DataSeries with its name specified by item

        The item is interpreted as the name of a series. VSeries names can have "-v"
        or "-y" as a suffix. The suffix "-t" or "-x" to a VSeries name can be used to
        get instead its corresponding TSeries. In any case, if there are more than one
        series with the name specified by item, they are appended. The timestamp is
        always shifted to the measurement's tstamp

        Args:
            item (str): The name of a DataSeries (see above)
        """
        ss = [s for s in self.series_list if s.name == item]
        if len(ss) == 1:
            s = ss[0]
        elif len(ss) > 1:
            s = append_series(ss)
        elif item[-2:] in ["-t", "-x", "-v", "-y"]:
            ss = [s for s in self.series_list if s.name == item[:-2]]
            if len(ss) == 1:
                s = ss[0]
            else:
                s = append_series(ss)
        else:
            raise SeriesNotFoundError
        return time_shifted(s, self.tstamp)

    def get_t_and_v(self, item, tspan=None):
        """Return the time and value vectors for a given VSeries name cut by tspan"""
        vseries = self[item]
        tseries = vseries.tseries
        v = vseries.data
        t = tseries.data + tseries.tstamp - self.tstamp
        if tspan:
            mask = np.logical_and(tspan[0] < t, t < tspan[-1])
            t, v = t[mask], v[mask]
        return t, v

    @property
    def data_cols(self):
        """Return a set of the names of all of the measurement's VSeries and TSeries"""
        return set([s.name for s in (self.value_series + self.time_series)])

    def plot(self, plotter=None, *args, **kwargs):
        """Plot the measurement using its plotter (see its Plotter for details)"""
        if plotter:
            return plotter.plot_measurement(self, *args, **kwargs)
        if not self.plotter:
            from .plotters import ValuePlotter

            self.plotter = ValuePlotter(measurement=self)
        return self.plotter.plot(*args, **kwargs)

    def export(self, exporter=None, *args, **kwargs):
        """Export the measurement using its exporter (see its Exporter for details)"""
        if exporter:
            return exporter.export_measurement(self, *args, **kwargs)
        return self.exporter.export(*args, **kwargs)

    def __add__(self, other):
        """Addition of measurements appends the series and component measurements lists.

        Adding results in a new Measurement. If the combination of the two measurements'
        techniques is a recognized hyphenated technique, it returns an object of that
        technique's measurement class. Otherwise it returns an object of Measurement.
        metadata, sample, and logentry come from the first measurement.

        An important point about addition is that it is almost but not quite associative
        and commutative i.e.
        A + (B + C) == (A + B) + C == C + B + A   is not quite true
        Each one results in the same series and component measurements. They will even
        appear in the same order in A + (B + C) and (A + B) + C. However, the technique
        might be different, as a new technique might be determined each time.

        Note also that there is no difference between hyphenating (simultaneous EC and
        MS datasets, for example) and appending (sequential EC datasets). Either way,
        all the raw series (or their placeholders) are just stored in the lists.
        TODO: Make sure with tests this is okay, differentiate using | operator if not.
        """
        obj_as_dict = self.as_dict()
        new_name = self.name + " AND " + other.name
        new_technique = self.technique + " AND " + other.technique

        # TODO: see if there isn't a way to put the import at the top of the module.
        #    see: https://github.com/ixdat/ixdat/pull/1#discussion_r546437410
        from .techniques import TECHNIQUE_CLASSES

        if new_technique in TECHNIQUE_CLASSES:
            cls = TECHNIQUE_CLASSES[new_technique]
        elif self.__class__ is other.__class__:
            cls = self.__class__
        else:
            cls = Measurement

        new_series_list = self.series_list
        new_series_list.append(other.series_list)
        new_component_measurements = self.component_measurements or [self]
        new_component_measurements.append(other.component_measurements or [other])
        obj_as_dict.update(
            name=new_name,
            series_list=new_series_list,
            component_measurements=new_component_measurements,
        )
        return cls.from_dict(obj_as_dict)


#  ------- Now come a few module-level functions for series manipulation ---------
#  TODO: move to a .build module or similar


def append_series(series_list):
    """Return series appending series_list relative to series_list[0].tseries.tstamp"""
    s0 = series_list[0]
    if isinstance(s0, TimeSeries):
        return append_tseries(series_list)
    elif isinstance(s0, ValueSeries):
        return append_vseries_by_time(series_list)
    raise BuildError(
        f"An algorithm of append_series for series like {s0} is not yet implemented"
    )


def append_vseries_by_time(series_list):
    """Return new series with the data in series_list appended"""
    name = series_list[0].name
    cls = series_list[0].__class__
    unit = series_list[0].unit
    data = np.array([])

    for s in series_list:
        if not (s.name == name and s.unit == unit and s.__class__ == cls):
            raise BuildError(f"can't append {series_list}")
        data = np.append(data, s.data)

    tseries = append_tseries([s.tseries for s in series_list])
    return cls(name=name, unit=unit, data=data, tseries=tseries)


def append_tseries(series_list, tstamp=None):
    """Return new TSeries with the data appended relative to series_list[0].tstamp"""
    name = series_list[0].name
    cls = series_list[0].__class__
    unit = series_list[0].unit
    tstamp = tstamp or series_list[0].tstamp
    data = np.array([])

    for s in series_list:
        if not (s.name == name and s.unit == unit and s.__class__ == cls):
            raise BuildError(f"can't append {series_list}")
        data = np.append(data, s.data + s.tstamp - tstamp)

    tseries = cls(name=name, unit=unit, data=data, tstamp=tstamp)
    return tseries


def fill_object_list(object_list, obj_ids, cls=None):
    """Add PlaceHolderObjects to object_list for any unrepresented obj_ids.

    Args:
        object_list (list of objects or None): The objects already known,
            in a list. This is the list to be appended to. If None, an empty
            list will be appended to.
        obj_ids (list of ints or None): The id's of objects to ensure are in
            the list. Any id in obj_ids not already represented in object_list
            is added to the list as a PlaceHolderObject
        cls (Saveable class): the class remembered by any PlaceHolderObjects
            added to the object_list, so that eventually the right object will
            be loaded.
    """
    cls = cls or object_list[0].__class__
    object_list = object_list or []
    provided_series_ids = [s.id for s in object_list]
    if not obj_ids:
        return object_list
    for i in obj_ids:
        if i not in provided_series_ids:
            object_list.append(PlaceHolderObject(i=i, cls=cls))
    return object_list


def time_shifted(series, tstamp=None):
    """Return a series with the time shifted to be relative to tstamp"""
    if not tstamp:
        return series
    cls = series.__class__
    if isinstance(series, TimeSeries):
        return cls(
            name=series.name,
            unit_name=series.unit.name,
            data=series.data + series.tstamp - tstamp,
            tstamp=tstamp,
        )
    elif isinstance(series, ValueSeries):
        series = cls(
            name=series.name,
            unit_name=series.unit.name,
            data=series.data,
            tseries=time_shifted(series.tseries, tstamp=tstamp),
        )
    return series
