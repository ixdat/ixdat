"""This module defines the Dataset class, the central data structure of ixdat

An ixdat Dataset is a collection of references to DataSeries with the metadata required
to combine them, i.e. "build" the combined dataset. It has a number of general methods
to visualize and analyze the combined dataset. Dataset is also the base class for a
number of technique-specific Dataset-derived classes.
"""
from pathlib import Path
import json
import numpy as np
from .db import Saveable, PlaceHolderObject
from .data_series import DataSeries, TimeSeries, ValueSeries
from .projects.samples import Sample
from .projects.lablogs import LabLog
from .exporters.csv_exporter import CSVExporter
from .exceptions import BuildError, SeriesNotFoundError  # , TechniqueError


class Measurement(Saveable):
    """The Measurement class"""

    table_name = "measurement"
    column_attrs = {
        "name",
        "technique",
        "metadata",
        "sample_name",
        "tstamp",
    }
    extra_linkers = {
        "measurement_series": ("data_series", "s_ids"),
        "component_measurements": ("measurements", "m_ids"),
    }

    sel_str = None  # the default thing to select on.
    #  FIXME: this is here because otherwise MSMeasurement.__init__ overwrites what it
    #   gets set to by ECMeasurement.__init__ in ECMSMeasurement.__init__

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
            metadata (dict): Free-form measurement metadata. Must be json-compatible.
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
            sample (Sample or str): The sample being measured
            lablog (LabLog): The log entry with e.g. notes taken during the measurement
            tstamp (float): The nominal starting time of the measurement, used for
                data selection, visualization, and exporting.
        """
        super().__init__()
        self.name = name
        self.technique = technique
        self.metadata = metadata or {}
        self.reader = reader
        self._plotter = plotter
        self._exporter = exporter
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

        # defining these methods here gets them the right docstrings :D
        self.plot_measurement = self.plotter.plot_measurement
        self.plot = self.plotter.plot_measurement

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

        # certain objects stored in the Measurement, but only saved as their names.
        #   __init__() will get the object from the name, but the argument is
        #   called like the object either way. For example __init__() takes an argument
        #   called `sample` which can be an ixdat.Sample or a string interpreted as the
        #   name of the sample to load. Subsequently, the sample name is accessible as
        #   the property `sample_name`. But in the database is only saved the sample's
        #   name as a string with the key/column "sample_name". So
        #   obj_as_dict["sample_name"] needs to be renamed obj_as_dict["sample"] before
        #   obj_as_dict can be passed to __init__.
        #   TODO: This is a rather general problem (see, e.g. DataSeries.unit vs
        #       DataSeries.unit_name) and as such should be moved to db.Saveable
        #       see: https://github.com/ixdat/ixdat/pull/5#discussion_r565090372
        objects_saved_as_their_name = [
            "sample",
        ]
        for object_type_str in objects_saved_as_their_name:
            object_name_str = object_type_str + "_name"
            if object_name_str in obj_as_dict:
                obj_as_dict[object_type_str] = obj_as_dict[object_name_str]
                del obj_as_dict[object_name_str]

        if obj_as_dict["technique"] in TECHNIQUE_CLASSES:
            # This makes it so that from_dict() can be used to initiate for any more
            # derived technique, so long as obj_as_dict specifies the technique name!
            technique_class = TECHNIQUE_CLASSES[obj_as_dict["technique"]]
            if not issubclass(technique_class, cls):
                # But we never want obj_as_dict["technique"] to take us to a *less*
                # specific technique, if the user has been intentional about which
                # class they call `as_dict` from (e.g. via a Reader)!
                technique_class = cls
        else:
            # Normally, we're going to want to make sure that we're in
            technique_class = cls

        if technique_class is cls:
            return cls(**obj_as_dict)
        else:  # Then its from_dict() might have more than ours:
            return technique_class.from_dict(obj_as_dict)

    @classmethod
    def read(cls, path_to_file, reader, **kwargs):
        """Return a Measurement object from parsing a file with the specified reader

        Args:
            path_to_file (Path or str): The path to the file to read
            reader (str or Reader class): The (name of the) reader to read the file with.
            kwargs: key-word arguments are passed on to the reader's read() method.
        """
        if isinstance(reader, str):
            # TODO: see if there isn't a way to put the import at the top of the module.
            #    see: https://github.com/ixdat/ixdat/pull/1#discussion_r546437471
            from .readers import READER_CLASSES

            reader = READER_CLASSES[reader]()
        # print(f"{__name__}. cls={cls}")  # debugging
        return reader.read(path_to_file, cls=cls, **kwargs)

    @classmethod
    def read_url(cls, url, reader, **kwargs):
        """Read a url (via a temporary file) using the specified reader"""
        from .readers.reading_tools import url_to_file

        path_to_temp_file = url_to_file(url)
        measurement = cls.read(path_to_temp_file, reader=reader, **kwargs)
        path_to_temp_file.unlink()
        return measurement

    @classmethod
    def read_set(
        cls, path_to_file_start, reader, suffix=None, file_list=None, **kwargs
    ):
        """Read and append a set of files.

        Args:
            path_to_file_start (Path or str): The path to the files to read including
                the shared start of the file name: `Path(path_to_file).parent` is
                interpreted as the folder where the file are.
                `Path(path_to_file).name` is interpreted as the shared start of the files
                to be appended.
            reader (str or Reader class): The (name of the) reader to read the files with
            file_list (list of Path): As an alternative to path_to_file_start, the
                exact files to append can be specified in a list
            suffix (str): If a suffix is given, only files with the specified ending are
                added to the file list
            kwargs: Key-word arguments are passed via cls.read() to the reader's read()
                method, AND to cls.from_component_measurements()
        """
        base_name = None
        if not file_list:
            folder = Path(path_to_file_start).parent
            base_name = Path(path_to_file_start).name
            file_list = [f for f in folder.iterdir() if f.name.startswith(base_name)]
            if suffix:
                file_list = [f for f in file_list if f.suffix == suffix]

        component_measurements = [
            cls.read(f, reader=reader, **kwargs) for f in file_list
        ]

        if base_name and "name" not in kwargs:
            kwargs["name"] = base_name
        measurement = cls.from_component_measurements(component_measurements, **kwargs)
        return measurement

    @classmethod
    def from_component_measurements(
        cls, component_measurements, keep_originals=True, sort=True, **kwargs
    ):
        """Return a measurement with the data contained in the component measurements

        TODO: This function "builds" the resulting measurement, i.e. it appends series
            of the same name rather than keeping all the original copies. This should be
            made more explicit, and a `build()` method should take over some of the work.

        Args:
            component_measurements (list of Measurement)
            keep_originals: Whether to keep a list of component_measurements referenced.
                This may result in redundant numpy arrays in RAM.
            sort (bool): Whether to sort the series according to time
            kwargs: key-word arguments are added to the dictionary for cls.from_dict()

        Returns cls: a Measurement object of the
        """

        # First prepare everything but the series_list in the object dictionary
        obj_as_dict = component_measurements[0].as_dict()
        obj_as_dict.update(kwargs)
        del obj_as_dict["m_ids"], obj_as_dict["s_ids"]
        if keep_originals:
            obj_as_dict["component_measurements"] = component_measurements

        # Now, prepare the built series. First, we loop through the component
        # measurements and get all the data and metadata organized in a dictionary:
        series_as_dicts = {}
        tstamp = component_measurements[0].tstamp
        for meas in component_measurements:
            tstamp_i = meas.tstamp  # save this for later.
            meas.tstamp = tstamp  # so that the time vectors share a t=0
            for s_name in meas.series_names:
                series = meas[s_name]
                if s_name in series_as_dicts:
                    series_as_dicts[s_name]["data"] = np.append(
                        series_as_dicts[s_name]["data"], series.data
                    )
                else:
                    series_as_dicts[s_name] = series.as_dict()
                    series_as_dicts[s_name]["data"] = series.data
                    if isinstance(series, ValueSeries):
                        # This will serve to match it to a TimeSeries later:
                        series_as_dicts[s_name]["t_name"] = series.tseries.name
            meas.tstamp = tstamp_i  # so it's not changed in the outer scope

        # Now we make DataSeries, starting with all the TimeSeries
        tseries_dict = {}
        sort_indeces = {}
        for name, s_as_dict in series_as_dicts.items():
            if "tstamp" in s_as_dict:
                if sort:
                    sort_indeces[name] = np.argsort(s_as_dict["data"])
                    s_as_dict["data"] = s_as_dict["data"][sort_indeces[name]]
                tseries_dict[name] = TimeSeries.from_dict(s_as_dict)
        # And then ValueSeries, and put both in with the TimeSeries
        series_list = []
        for name, s_as_dict in series_as_dicts.items():
            if name in tseries_dict:
                series_list.append(tseries_dict[name])
            elif "t_name" in s_as_dict:
                tseries = tseries_dict[s_as_dict["t_name"]]
                if s_as_dict["data"].shape == tseries.shape:
                    # Then we assume that the time and value data have lined up
                    # successfully! :D
                    if sort:
                        s_as_dict["data"] = s_as_dict["data"][
                            sort_indeces[tseries.name]
                        ]
                    vseries = ValueSeries(
                        name=name,
                        data=s_as_dict["data"],
                        unit_name=s_as_dict["unit_name"],
                        tseries=tseries,
                    )
                else:
                    # this will be the case if vseries sharing the same tseries
                    # are not present in the same subset of component_measurements.
                    # In that case just append the vseries even though some tdata gets
                    # duplicated.
                    vseries = append_series(
                        [
                            s
                            for m in component_measurements
                            for s in m.series_list
                            if s.name == name
                        ],
                        sort=sort,
                    )
                series_list.append(vseries)

        # Finally, add this series to the dictionary representation and return the object
        obj_as_dict["series_list"] = series_list
        return cls.from_dict(obj_as_dict)

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
        # TimeSeries have to go first, so that ValueSeries are saved with the right t_id!
        data_object_list = self.time_series
        for s in self.series_list:
            if s not in data_object_list:
                if s.tseries not in data_object_list:
                    # FIXME: some tseries, likely with duplicate data, seem to not
                    #  make it into series_list
                    data_object_list.append(s.tseries)
                data_object_list.append(s)
        return data_object_list

    @property
    def component_measurements(self):
        """List of the component measurements of which this measurement is a combination

        For a pure measurement (not a measurement set), this is itself in a list.
        """
        if not self._component_measurements:
            return [
                self,
            ]
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
        return {(s.id, s.backend_name): s for s in self.series_list}

    @property
    def series_names(self):
        """List of the names of the series in the measurement"""
        return set([series.name for series in self.series_list])

    @property
    def value_names(self):
        """List of the names of the VSeries in the measurement's DataSeries"""
        return set([vseries.name for vseries in self.value_series])

    @property
    def value_series(self):
        """List of the VSeries in the measurement's DataSeries"""
        return [
            series for series in self.series_list if isinstance(series, ValueSeries)
        ]

    @property
    def time_names(self):
        """List of the names of the VSeries in the measurement's DataSeries"""
        return set([tseries.name for tseries in self.time_series])

    @property
    def time_series(self):
        """List of the TSeries in the measurement's DataSeries. NOT timeshifted!"""
        return [series for series in self.series_list if isinstance(series, TimeSeries)]

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
            raise SeriesNotFoundError(f"{self} has no series called {item}")
        if hasattr(s, "tstamp") and not s.tstamp == self.tstamp:
            s = time_shifted(s, self.tstamp)
        return s

    def __setitem__(self, series_name, series):
        """Append `series` with name=`series_name` to `series_list` and remove others."""
        if not series.name == series_name:
            raise SeriesNotFoundError(
                f"Can't set {self}[{series_name}] = {series}. Series names don't agree."
            )
        del self[series_name]
        self.series_list.append(series)

    def __delitem__(self, series_name):
        """Remove all series which have `series_name` as their name from series_list"""
        new_series_list = []
        for s in self.series_list:
            if not s.name == series_name:
                new_series_list.append(s)
        self._series_list = new_series_list

    def correct_data(self, value_name, new_data):
        """Replace the old data for ´value_name´ (str) with ´new_data` (np array)"""
        old_vseries = self[value_name]
        new_vseries = ValueSeries(
            name=value_name,
            unit_name=old_vseries.unit_name,
            data=new_data,
            tseries=old_vseries.tseries,
        )
        self[value_name] = new_vseries

    def grab(self, item, tspan=None, include_endpoints=False, tspan_bg=None):
        """Return a value vector with the corresponding time vector

        Grab is the *canonical* way to retrieve numerical time-dependent data from a
        measurement in ixdat. The first argument is always the name of the value to get
        time-resolved data for (the name of a ValueSeries). The second, optional,
        argument is a timespan to select the data for.
        Two vectors are returned: first time (t), then value (v). They are of the same
        length so that `v` can be plotted against `t`, integrated over `t`, interpolated
        via `t`, etc. `t` and `v` are returned in the units of their DataSeries.
        TODO: option to specifiy desired units

        Typical usage::
            t, v = measurement.grab(potential, tspan=[0, 100])

        Args:
            item (str): The name of the DataSeries to grab data for
            tspan (iter of float): Defines the timespan with its first and last values.
                Optional. By default the entire time of the measurement is included.
            include_endpoints (bool): Whether to add a points at t = tspan[0] and
                t = tspan[-1] to the data returned. This makes trapezoidal integration
                less dependent on the time resolution. Default is False.
            tspan_bg (iterable): Optional. A timespan defining when `item` is at its
                baseline level. The average value of `item` in this interval will be
                subtracted from the values returned.
        """
        vseries = self[item]
        tseries = vseries.tseries
        v = vseries.data
        t = tseries.data + tseries.tstamp - self.tstamp
        if tspan is not None:  # np arrays don't boolean well :(
            if include_endpoints:
                if t[0] < tspan[0]:  # then add a point to include tspan[0]
                    v_0 = np.interp(tspan[0], t, v)
                    t = np.append(tspan[0], t)
                    v = np.append(v_0, v)
                if tspan[-1] < t[-1]:  # then add a point to include tspan[-1]
                    v_end = np.interp(tspan[-1], t, v)
                    t = np.append(t, tspan[-1])
                    v = np.append(v, v_end)
            mask = np.logical_and(tspan[0] <= t, t <= tspan[-1])
            t, v = t[mask], v[mask]
        if tspan_bg:
            t_bg, v_bg = self.grab(item, tspan=tspan_bg)
            v = v - np.mean(v_bg)
        return t, v

    def grab_for_t(self, item, t, tspan_bg=None):
        """Return a numpy array with the value of item interpolated to time t

        Args:
            item (str): The name of the value to grab
            t (np array): The time vector to grab the value for
            tspan_bg (iterable): Optional. A timespan defining when `item` is at its
                baseline level. The average value of `item` in this interval will be
                subtracted from what is returned.
        """
        vseries = self[item]
        tseries = vseries.tseries
        v_0 = vseries.data
        t_0 = tseries.data + tseries.tstamp - self.tstamp
        v = np.interp(t, t_0, v_0)
        if tspan_bg:
            t_bg, v_bg = self.grab(item, tspan=tspan_bg)
            v = v - np.mean(v_bg)
        return v

    def integrate(self, item, tspan=None, ax=None):
        """Return the time integral of item in the specified timespan"""
        t, v = self.grab(item, tspan, include_endpoints=True)
        if ax:
            if ax == "new":
                ax = self.plotter.new_ax(ylabel=item)
                # FIXME: xlabel=self[item].tseries.name gives a problem :(
            ax.plot(t, v, color="k", label=item)
            ax.fill_between(t, v, np.zeros(t.shape), where=v > 0, color="g", alpha=0.3)
            ax.fill_between(
                t, v, np.zeros(t.shape), where=v < 0, color="g", alpha=0.1, hatch="//"
            )

        return np.trapz(v, t)

    @property
    def data_cols(self):
        """Return a set of the names of all of the measurement's VSeries and TSeries"""
        return set([s.name for s in (self.value_series + self.time_series)])

    @property
    def plotter(self):
        """The default plotter for Measurement is ValuePlotter."""
        if not self._plotter:
            from .plotters import ValuePlotter

            # FIXME: I had to import here to avoid running into circular import issues

            self._plotter = ValuePlotter(measurement=self)
        return self._plotter

    @property
    def exporter(self):
        """The default exporter for Measurement is CSVExporter."""
        if not self._exporter:
            self._exporter = CSVExporter(measurement=self)
        return self._exporter

    def export(self, *args, exporter=None, **kwargs):
        """Export the measurement using its exporter (see its Exporter for details)"""
        if exporter:
            return exporter.export_measurement(self, *args, **kwargs)
        return self.exporter.export(*args, **kwargs)

    def get_original_m_id_of_series(self, series):
        """Return the id(s) of component measurements to which `series` belongs."""
        m_id_list = []
        for m in self.component_measurements:
            if series.id in m.s_ids:
                m_id_list.append(m.id)
        if len(m_id_list) == 1:
            return m_id_list[0]
        return m_id_list

    def cut(self, tspan, t_zero=None):
        """Return a new measurement with the data in the given time interval

        Args:
            tspan (iter of float): The time interval to use, relative to self.tstamp
                tspan[0] is the start time of the interval, and tspan[-1] is the end
                time of the interval. Using tspan[-1] means you can directly use a
                long time vector that you have at hand to describe the time interval
                you're looking for.
            t_zero (float or str): Where to put the tstamp of the returned measurement.
                Default is to keep it the same as the present tstamp. If instead it is
                a float, this adds the float to the present tstamp. If t_zero is "start",
                tspan[0] is added to the present tstamp.
        """
        new_series_list = []
        obj_as_dict = self.as_dict()
        time_cutting_stuff = {}  # {tseries_id: (mask, new_tseries)}
        for series in self.series_list:
            try:
                tseries = series.tseries
                if tseries is None:
                    raise AttributeError
            except AttributeError:  # series independent of time are uneffected by cut
                new_series_list.append(series)
            else:
                t_id = (tseries.id, tseries.backend_name)
                # FIXME: Beautiful, met my first id clash here. Local memory and loaded
                #    each had a timeseries with id=1, but different length. Previously
                #    the above line of code was just t_id = tseries.id as you'd expect,
                #    meaning that time_cutting_stuff appeared to already have the needed
                #    tseries but didn't!
                #    Note that the id together with the backend works but should be
                #    replaced by a single Universal Unique Identifier, or perhaps just
                #    a property `Saveable.uid`, returning `(self.id, self.backend_name)`

                if t_id in time_cutting_stuff:
                    mask, new_tseries = time_cutting_stuff[t_id]
                else:
                    t = tseries.t + tseries.tstamp - self.tstamp
                    mask = np.logical_and(tspan[0] <= t, t <= tspan[-1])
                    new_tseries = TimeSeries(
                        name=tseries.name,
                        unit_name=tseries.unit_name,
                        tstamp=tseries.tstamp,
                        data=tseries.data[mask],
                    )
                    time_cutting_stuff[t_id] = (mask, new_tseries)
                if True not in mask:
                    continue
                if False not in mask:
                    new_series_list.append(series)
                elif (series.id, series.backend_name) == t_id:
                    new_series_list.append(new_tseries)
                else:
                    new_series = series.__class__(
                        name=series.name,
                        unit_name=series.unit_name,
                        data=series.data[mask],
                        tseries=new_tseries,
                    )
                    new_series_list.append(new_series)
        obj_as_dict["series_list"] = new_series_list
        del obj_as_dict["s_ids"]
        if t_zero:
            if t_zero == "start":
                t_zero = tspan[0]
            obj_as_dict["tstamp"] += t_zero
        new_measurement = self.__class__.from_dict(obj_as_dict)
        return new_measurement

    def select_value(self, *args, **kwargs):
        """Return a new Measurement with the time(s) meeting criteria.

        Can only take one arg or kwarg!
        The `series_name` is `self.sel_str` if given an arg, kw if given a kwarg.
        Either way the argument is the `value` to be selected for.

        The method finds all time intervals for which `self[series_name] == value`
        It then cuts the measurement according to each time interval and adds these
        segments together. TODO: This can be done better, i.e. without chopping series.

        TODO: greater-than and less-than kwargs.
            Ideally you should be able to say e.g., `select(cycle=1, 0.5<potential<1)`
        """
        if len(args) >= 1:
            if not self.sel_str:
                raise BuildError(
                    f"{self} does not have a default selection string "
                    f"(Measurement.sel_str), and so selection only works with kwargs."
                )
            kwargs[self.sel_str] = args
        if len(kwargs) > 1:
            raise BuildError(
                f"select_value got kwargs={kwargs} but can only be used for one value "
                f"at a time. Use select_values for more."
            )
        new_measurement = self
        ((series_name, value),) = kwargs.items()

        t, v = self.grab(series_name)
        mask = v == value  # linter doesn't realize this is a np array
        mask_prev = np.append(False, mask[:-1])
        mask_next = np.append(mask[1:], False)
        interval_starts_here = np.logical_and(
            np.logical_not(mask_prev), mask
        )  # True at [0] if mask[0] is True.
        interval_ends_here = np.logical_and(
            mask, np.logical_not(mask_next)
        )  # True at [-1] if mask[-1] is True.
        t_starts = list(t[interval_starts_here])
        t_ends = list(t[interval_ends_here])
        tspans = zip(t_starts, t_ends)
        meas = None
        for tspan in tspans:
            if meas:
                meas = meas + new_measurement.cut(tspan)
            else:
                meas = new_measurement.cut(tspan)
        new_measurement = meas
        return new_measurement

    def select_values(self, *args, **kwargs):
        """Return a new Measurement with the time(s) in the measurement meeting criteria

        Any series can be selected for using the series name as a key-word. Arguments
        can be single acceptable values or lists of acceptable values. In the latter
        case, each acceptable value is selected for on its own and the resulting
        measurements added together.
        FIXME: That is sloppy because it multiplies the number of DataSeries
        FIXME:  containing the same amount of data.
        If no key-word is given, the series name is assumed to
        be the default selector, which is named by self.sel_str. Multiple criteria are
        applied sequentially, i.e. you get the intersection of satisfying parts.

        Args:
            args (tuple): Argument(s) given without key-word are understood as acceptable
                value(s) for the default selector (that named by self.sel_str)
            kwargs (dict): Each key-word arguments is understood as the name
                of a series and its acceptable value(s).
        """

        if len(args) >= 1:
            if not self.sel_str:
                raise BuildError(
                    f"{self} does not have a default selection string "
                    f"(Measurement.sel_str), and so selection only works with kwargs."
                )
            if len(args) == 1:
                args = args[0]
            kwargs[self.sel_str] = args
        new_measurement = self
        for series_name, allowed_values in kwargs.items():
            if not hasattr(allowed_values, "__iter__"):
                allowed_values = [allowed_values]
            meas = None
            for value in allowed_values:
                m = new_measurement.select_value(**{series_name: value})
                if meas:
                    meas = meas + m
                else:
                    meas = m
            new_measurement = meas
        return new_measurement

    def select(self, *args, tspan=None, **kwargs):
        """`cut` (with tspan) and `select_values` (with args and/or kwargs)."""
        new_measurement = self
        if tspan:
            new_measurement = new_measurement.cut(tspan=tspan)
        if args or kwargs:
            new_measurement = new_measurement.select_values(*args, **kwargs)
        return new_measurement

    @property
    def tspan(self):
        """Return `(t_start, t_finish)` interval including all data in the measurement"""
        t_start = None
        t_finish = None
        for tcol in self.time_names:
            t = self[tcol].data
            t_start = min(t_start, t[0]) if t_start else t[0]
            t_finish = max(t_finish, t[-1]) if t_finish else t[-1]
        return t_start, t_finish

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

        # First we prepare a dictionary for all but the series_list.
        # This has both dicts, but prioritizes self's dict for all that appears twice.
        obj_as_dict = self.as_dict()
        other_as_dict = other.as_dict()
        for k, v in other_as_dict.items():
            # Looking forward to the "|" operator!
            if k not in obj_as_dict:
                obj_as_dict[k] = v
        new_name = self.name + " AND " + other.name
        new_technique = get_combined_technique(self.technique, other.technique)

        # TODO: see if there isn't a way to put the import at the top of the module.
        #    see: https://github.com/ixdat/ixdat/pull/1#discussion_r546437410
        from .techniques import TECHNIQUE_CLASSES

        if new_technique in TECHNIQUE_CLASSES:
            cls = TECHNIQUE_CLASSES[new_technique]
        elif self.__class__ is other.__class__:
            cls = self.__class__
        else:
            cls = Measurement

        new_series_list = self.series_list + other.series_list
        new_component_measurements = (
            self.component_measurements + other.component_measurements
        )
        obj_as_dict.update(
            name=new_name,
            technique=new_technique,
            series_list=new_series_list,
            component_measurements=new_component_measurements,
        )
        return cls.from_dict(obj_as_dict)

    def join(self, other, join_on=None):
        """Join two measurements based on a shared data series

        This involves projecting all timeseries from other's data series so that the
        variable named by `join_on` is shared between all data series.
        This is analogous to an explicit inner join.

        Args:
            other (Measurement): a second measurement to join to self
            join_on (str or tuple): Either a string, if the value to join on is called
                the same thing in both measurements, or a tuple of two strings if it is
                not.
                The variable described by join_on must be monotonically increasing in
                both measurements.
        """


#  ------- Now come a few module-level functions for series manipulation ---------
# TODO: move to an `ixdat.build` module or similar.
#   There's a lot of stuff that should go there. Basically anything in ECMeasurement
#   that can be reasonably converted to a module level function to decrease the
#   awkwardness there.


def append_series(series_list, sort=True, tstamp=None):
    """Return series appending series_list relative to series_list[0].tseries.tstamp

    Args:
        series_list (list of Series): The series to append (must all be of same type)
        sort (bool): Whether to sort the data so that time only goes forward
        tstamp (unix tstamp): The t=0 of the returned series or its TimeSeries.
    """
    s0 = series_list[0]
    if isinstance(s0, TimeSeries):
        return append_tseries(series_list, sort=sort, tstamp=tstamp)
    elif isinstance(s0, ValueSeries):
        return append_vseries_by_time(series_list, sort=sort, tstamp=tstamp)
    raise BuildError(
        f"An algorithm of append_series for series like {s0} is not yet implemented"
    )


def append_vseries_by_time(series_list, sort=True, tstamp=None):
    """Return new ValueSeries with the data in series_list appended

    Args:
        series_list (list of ValueSeries): The value series to append
        sort (bool): Whether to sort the data so that time only goes forward
        tstamp (unix tstamp): The t=0 of the returned ValueSeries' TimeSeries.
    """
    name = series_list[0].name
    cls = series_list[0].__class__
    unit = series_list[0].unit
    data = np.array([])
    tseries_list = [s.tseries for s in series_list]
    tseries, sort_indeces = append_tseries(
        tseries_list, sort=sort, return_sort_indeces=True, tstamp=tstamp
    )

    for s in series_list:
        if not (s.unit == unit and s.__class__ == cls):
            raise BuildError(f"can't append {series_list}")
        data = np.append(data, s.data)
    if sort:
        data = data[sort_indeces]

    return cls(name=name, unit_name=unit.name, data=data, tseries=tseries)


def append_tseries(series_list, sort=True, return_sort_indeces=False, tstamp=None):
    """Return new TimeSeries with the data appended.

    Args:
        series_list (list of TimeSeries): The time series to append
        sort (bool): Whether to sort the data so that time only goes forward
        return_sort_indeces (bool): Whether to return the indeces that sort the data
        tstamp (unix tstamp): The t=0 of the returned TimeSeries.
    """
    name = series_list[0].name
    cls = series_list[0].__class__
    unit = series_list[0].unit
    tstamp = tstamp or series_list[0].tstamp
    data = np.array([])

    for s in series_list:
        if not (s.unit == unit and s.__class__ == cls):
            raise BuildError(f"can't append {series_list}")
        data = np.append(data, s.data + s.tstamp - tstamp)

    if sort:
        sort_indices = np.argsort(data)
        data = data[sort_indices]
    else:
        sort_indices = None

    tseries = cls(name=name, unit_name=unit.name, data=data, tstamp=tstamp)
    if return_sort_indeces:
        return tseries, sort_indices
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
    if tstamp is None or not series:
        return series
    if tstamp == series.tstamp:
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


def get_combined_technique(technique_1, technique_2):
    """Return the name of the technique resulting from adding two techniques"""
    # TODO: see if there isn't a way to put the import at the top of the module.
    #    see: https://github.com/ixdat/ixdat/pull/1#discussion_r546437410
    if technique_1 == technique_2:
        return technique_1

    # if we're a component technique of a hyphenated technique to that hyphenated
    # technique, the result is still the hyphenated technique. e.g. EC-MS + MS = EC-MS
    if "-" in technique_1 and technique_2 in technique_1.split("-"):
        return technique_1
    elif "-" in technique_2 and technique_1 in technique_2.split("-"):
        return technique_2

    # if we're adding two independent technique which are components of a hyphenated
    # technique, then we want that hyphenated technique. e.g. EC + MS = EC-MS
    from .techniques import TECHNIQUE_CLASSES

    for hyphenated in [
        technique_1 + "-" + technique_2,
        technique_2 + "-" + technique_1,
    ]:
        if hyphenated in TECHNIQUE_CLASSES:
            return hyphenated

    # if all else fails, we just join them with " AND ". e.g. MS + XRD = MS AND XRD
    return technique_1 + " AND " + technique_2
