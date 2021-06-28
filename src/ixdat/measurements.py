"""This module defines the Measurement class, the central data structure of ixdat

An ixdat Measurement is a collection of references to DataSeries with the metadata needed
to combine them, i.e. "build" the combined dataset. It has a number of general methods
to visualize and analyze the combined dataset. Measurement is also the base class for a
number of technique-specific Measurement-derived classes.

A Measurement will typically be accompanied by one or more Calibration. This module
also defines the base class for Calibration, while technique-specific Calibration
classes will be defined in the corresponding module in ./techniques/
"""
import json
import numpy as np
from .db import Saveable, PlaceHolderObject, fill_object_list
from .data_series import (
    DataSeries,
    TimeSeries,
    ValueSeries,
    ConstantValue,
    append_series,
)
from .samples import Sample
from .lablogs import LabLog
from .exporters.csv_exporter import CSVExporter
from .plotters.value_plotter import ValuePlotter
from .exceptions import BuildError, SeriesNotFoundError


class Measurement(Saveable):
    """The Measurement class"""

    # ------ table description class attributes --------
    table_name = "measurement"
    column_attrs = {
        "name",
        "technique",
        "metadata",
        "aliases",
        "sample_name",
        "tstamp",
    }
    extra_linkers = {
        "component_measurements": ("measurements", "m_ids"),
        "measurement_calibrations": ("calibration", "c_ids"),
        "measurement_series": ("data_series", "s_ids"),
    }
    child_attrs = ["_component_measurements", "_calibration_list", "series_list"]

    # ---- measurement class attributes, can be overwritten in inheriting classes ---- #
    control_technique = None
    """Name of technique primarily used to control the experiment"""
    control_str = None
    """Name (or alias) for main time variable or main time-dependent value variable,
    typically of the control_technique"""
    sel_str = "selector"
    """Name of the default selector"""
    select_on = ("file_number",)
    """Name of the default things to use to construct the selector"""
    series_constructors = {
        "file_number": "_build_file_number",
        "selector": "_build_selector",
    }
    """Series which should be constructed from other series by the specified method
    and cached the first time they are looked up"""
    default_plotter_class = ValuePlotter
    default_exporter_class = CSVExporter

    def __init__(
        self,
        name,
        technique=None,
        metadata=None,
        s_ids=None,
        series_list=None,
        c_ids=None,
        calibration_list=None,
        m_ids=None,
        component_measurements=None,
        aliases=None,
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
        self._calibration_list = fill_object_list(
            calibration_list, c_ids, cls=Calibration
        )
        self.tstamp = tstamp

        self._cached_series = {}
        self._aliases = aliases or {}

        # defining these methods here gets them the right docstrings :D
        self.plot_measurement = self.plotter.plot_measurement
        self.plot = self.plotter.plot_measurement
        self.export = self.exporter.export

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
            technique_class = TECHNIQUE_CLASSES[obj_as_dict["technique"]]
        else:
            technique_class = cls
        try:
            measurement = technique_class(**obj_as_dict)
        except Exception:
            raise
        return measurement

    @classmethod
    def read(cls, path_to_file, reader, **kwargs):
        """Return a Measurement object from parsing a file with the specified reader"""
        if isinstance(reader, str):
            # TODO: see if there isn't a way to put the import at the top of the module.
            #    see: https://github.com/ixdat/ixdat/pull/1#discussion_r546437471
            from .readers import READER_CLASSES

            reader = READER_CLASSES[reader]()
        return reader.read(path_to_file, **kwargs)  # TODO: take cls as kwarg

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
    def component_measurements(self):
        """List of the component measurements of which this measurement is a combination

        For a pure measurement (not a measurement set), this is itself in a list.
        """
        if not self._component_measurements:
            return [self]
        for i, m in enumerate(self._component_measurements):
            if isinstance(m, PlaceHolderObject):
                self._component_measurements[i] = m.get_object()
        return self._component_measurements

    @property
    def m_ids(self):
        """List of the id's of a combined measurement's component measurements"""
        if not self._component_measurements:
            return None
        return [m.id for m in self._component_measurements]

    @property
    def calibration_list(self):
        for i, c in enumerate(self._calibration_list):
            if isinstance(c, PlaceHolderObject):
                self._calibration_list[i] = c.get_object()
        return self._calibration_list

    @property
    def c_ids(self):
        """List of the id's of the measurement's Calibrations"""
        return [c.id for c in self._calibration_list]

    @property
    def series_list(self):
        """List of the DataSeries containing the measurement's data"""
        for i, s in enumerate(self._series_list):
            if isinstance(s, PlaceHolderObject):
                self._series_list[i] = s.get_object()
        return self._series_list

    @property
    def s_ids(self):
        """List of the id's of the measurement's DataSeries"""
        return [series.id for series in self._series_list]

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
    def time_names(self):
        """List of the names of the VSeries in the measurement's DataSeries"""
        return set([tseries.name for tseries in self.time_series])

    @property
    def value_series(self):
        """List of the VSeries in the measurement's DataSeries"""
        return [
            series for series in self.series_list if isinstance(series, ValueSeries)
        ]

    @property
    def time_series(self):
        """List of the TSeries in the measurement's DataSeries. NOT timeshifted!"""
        return [series for series in self.series_list if isinstance(series, TimeSeries)]

    @property
    def aliases(self):
        """Dictionary of {key: series_names} pointing to where desired raw data is"""
        return self._aliases

    def get_series_names(self, key):
        """Return list: series names for key found by (recursive) lookup in aliases"""
        keys = [key] if key in self.series_names else []
        for k in self.aliases.get(key, []):
            keys += self.get_series_names(k)
        return keys

    def __getitem__(self, key):
        """Return the built measurement DataSeries with its name specified by key

        The item is interpreted as the name of a series. VSeries names can have "-v"
        or "-y" as a suffix. The suffix "-t" or "-x" to a VSeries name can be used to
        get instead its corresponding TSeries.
        The measurement also checks its `aliases` for other name(s) which may refer to
        what is being looked up.
        If there are more than one series with the name(s) specified, they are appended.
        The timestamp is always shifted to the measurement's tstamp.

        Args:
            key (str): The name of a DataSeries (see above)
        """
        if key in self._cached_series:
            return self._cached_series[key]
        if key in self.series_constructors:
            series = getattr(self, self.series_constructors[key])()
        else:
            for calibration in self.calibration_list:
                series = calibration.calibrate_series(key, measurement=self)
                # ^ the calibration will call this __getitem__ with the name of the
                #   corresponding raw data and return a new series with calibrated data
                #   if possible. Otherwise it will return None.
                if series:
                    break
            else:
                # only if the requested series name is neither cached nor the name of
                #   a calibrated series do we go into raw data:
                series_to_append = []
                if key in self.series_names:
                    # Then we'll append any series matching the desired name
                    series_to_append += [s for s in self.series_list if s.name == key]
                elif key in self.aliases:
                    # Then we'll look up the aliases instead and append them
                    for k in self.aliases[key]:
                        try:
                            series_to_append.append(self[k])
                        except SeriesNotFoundError:
                            continue
                # If the key is something in the data, by now we have series to append.
                if not series_to_append:  # If not...
                    # check if it's because key uses a suffix:
                    if key.endswith("-t") or key.endswith("-x"):
                        return self[key[:-2]].tseries
                    if key.endswith("-v") or key.endswith("-y"):
                        return self[key[:-2]]
                    # if not, we've exhausted possibilities! Raise an error.
                    raise SeriesNotFoundError(f"{self} has no series named {key}.")
                # Now, append the series
                name = key if len(set(s.name for s in series_to_append)) > 1 else None
                series = append_series(series_to_append, name=name, tstamp=self.tstamp)
        # Finally, wherever we found the series, cache it and return it.
        self._cached_series[key] = series
        return series

    def __setitem__(self, series_name, series):
        """Append `series` with name=`series_name` to `series_list` and remove others"""
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

    def clear_cache(self):
        self._cached_series = {}

    def grab(self, item, tspan=None):
        """Return the time and value vectors for a given VSeries name cut by tspan"""
        vseries = self[item]
        tseries = vseries.tseries
        v = vseries.data
        t = tseries.data + tseries.tstamp - self.tstamp
        if tspan:
            mask = np.logical_and(tspan[0] < t, t < tspan[-1])
            t, v = t[mask], v[mask]
        return t, v

    def grab_for_t(self, item, t):
        """Return a numpy array with the value of item interpolated to time t"""
        vseries = self[item]
        tseries = vseries.tseries
        v_0 = vseries.data
        t_0 = tseries.data + tseries.tstamp - self.tstamp
        v = np.interp(t, t_0, v_0)
        return v

    @property
    def t(self):
        return self[self.control_str].t

    @property
    def t_str(self):
        return self[self.control_str].tseries.name

    def _build_file_number(self):
        series_to_append = []
        for i, m in enumerate(self.component_measurements):
            if self.control_technique and not m.technique == self.control_technique:
                continue
            if not self.control_str:
                tseries = m.time_series[0]
            else:
                try:
                    tseries = m[self.control_str].tseries
                except KeyError:
                    continue
            series_to_append.append(
                ConstantValue(name="file_number", unit_name="", data=i, tseries=tseries)
            )
        return append_series(series_to_append, name="file_number", tstamp=self.tstamp)

    def _build_selector(self, sel_str=None, col_list=None):
        """Build `selector` from col_list

        See the class docstring for details.
        """
        sel_str = sel_str or self.sel_str
        changes = np.tile(False, self.t.shape)
        col_list = col_list or self.select_on
        for col in col_list:
            try:
                vseries = self[col]
            except SeriesNotFoundError:
                continue
            values = vseries.data
            if len(values) == 0:
                print("WARNING: " + col + " is empty")
                continue
            elif not len(values) == len(changes):
                print("WARNING: " + col + " has an unexpected length")
                continue
            n_down = np.append(
                values[0], values[:-1]
            )  # comparing with n_up instead puts selector a point ahead
            changes = np.logical_or(changes, n_down < values)
        selector = np.cumsum(changes)
        selector_series = ValueSeries(
            name=sel_str,
            unit_name="",
            data=selector,
            tseries=self[self.control_str].tseries,
        )
        return selector_series

    @property
    def selector(self):
        return self[self.sel_str]

    @property
    def data_cols(self):
        """Return a set of the names of all of the measurement's VSeries and TSeries"""
        return set([s.name for s in (self.value_series + self.time_series)])

    @property
    def plotter(self):
        """The default plotter for Measurement is ValuePlotter."""
        if not self._plotter:
            self._plotter = self.default_plotter_class(measurement=self)

        return self._plotter

    @property
    def exporter(self):
        """The default exporter for Measurement is CSVExporter."""
        if not self._exporter:
            self._exporter = self.default_exporter_class(measurement=self)
        return self._exporter

    def get_original_m_id_of_series(self, series):
        """Return the id(s) of component measurements to which `series` belongs."""
        m_id_list = []
        for m in self.component_measurements:
            if series.id in m.s_ids:
                m_id_list.append(m.id)
        if len(m_id_list) == 1:
            return m_id_list[0]
        return m_id_list

    @property
    def tspan(self):
        """The minimum timespan (with respect to self.tstamp) containing all the data"""
        t_start = None
        t_finish = None
        for t_name in self.time_names:
            t = self[t_name].data
            t_start = min(t_start, t[0]) if t_start else t[0]
            t_finish = max(t_finish, t.data[-1]) if t_finish else t[-1]
        return [t_start, t_finish]

    def cut(self, tspan):
        """Return a new measurement with the data in the given time interval

        Args:
            tspan (iter of float): The time interval to use, relative to self.tstamp
                tspan[0] is the start time of the interval, and tspan[-1] is the end
                time of the interval. Using tspan[-1] means you can directly use a
                long time vector that you have at hand to describe the time interval
                you're looking for.
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
        obj_as_dict["calibration_list"] = self.calibration_list
        del obj_as_dict["c_ids"]
        obj_as_dict["component_measurements"] = []
        for m in self._component_measurements:
            dt = m.tstamp - self.tstamp
            tspan_m = [tspan[0] - dt, tspan[1] - dt]
            if m.tspan[-1] < tspan_m[0] or tspan_m[-1] < m.tspan[0]:
                continue
            obj_as_dict["component_measurements"].append(m.cut(tspan_m))
        del obj_as_dict["m_ids"]
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
        # FIXME: That is sloppy because it mutliplies the number of DataSeries
            containing the same amount of data.
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
        """`cut` (with tspan) and `select_values` (with *args and/or **kwargs)."""
        new_measurement = self
        if tspan:
            new_measurement = new_measurement.cut(tspan=tspan)
        if args or kwargs:
            new_measurement = new_measurement.select_values(*args, **kwargs)
        return new_measurement

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

        new_series_list = list(set(self.series_list + other.series_list))
        new_component_measurements = list(
            set(self.component_measurements + other.component_measurements)
        )
        new_calibration_list = list(
            set(self._calibration_list + other._calibration_list)
        )
        new_aliases = self.aliases.copy()
        for key, names in other.aliases.items():
            if key in new_aliases:
                new_aliases[key] = list(set(new_aliases[key] + other.aliases[key]))
            else:
                new_aliases[key] = other.aliases[key]
        obj_as_dict.update(
            name=new_name,
            series_list=new_series_list,
            component_measurements=new_component_measurements,
            calibration_list=new_calibration_list,
            aliases=new_aliases,
        )
        del obj_as_dict["c_ids"]
        del obj_as_dict["m_ids"]
        del obj_as_dict["s_ids"]
        return cls.from_dict(obj_as_dict)


class Calibration(Saveable):
    table_name = "calibration"
    column_attrs = {
        "name",
        "technique",
        "tstamp",
    }

    def __init__(self, name=None, technique=None, tstamp=None, measurement=None):
        super().__init__()
        self.name = name or f"{self.__class__.__name__}({measurement})"
        self.technique = technique
        self.tstamp = tstamp or (measurement.tstamp if measurement else None)
        self.measurement = measurement

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
        from .techniques import CALIBRATION_CLASSES

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

        if obj_as_dict["technique"] in CALIBRATION_CLASSES:
            calibration_class = CALIBRATION_CLASSES[obj_as_dict["technique"]]
        else:
            calibration_class = cls
        try:
            measurement = calibration_class(**obj_as_dict)
        except Exception:
            raise
        return measurement
