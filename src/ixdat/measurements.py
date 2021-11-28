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
from .db import Savable, PlaceHolderObject, fill_object_list
from .data_series import (
    DataSeries,
    TimeSeries,
    ValueSeries,
    ConstantValue,
    append_series,
    time_shifted,
    get_tspans_from_mask,
)
from .samples import Sample
from .lablogs import LabLog
from .exporters.csv_exporter import CSVExporter
from .plotters.value_plotter import ValuePlotter
from .exceptions import BuildError, SeriesNotFoundError
from .tools import dict_is_close


class Measurement(Savable):
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
    child_attrs = ["component_measurements", "calibration_list", "series_list"]
    # TODO: child_attrs should be derivable from extra_linkers?

    # ---- measurement class attributes, can be overwritten in inheriting classes ---- #
    control_technique_name = None
    """Name of the control technique primarily used to control the experiment"""
    control_series_name = None
    """Name (or alias) for main time variable or main time-dependent value variable,
    typically of the control technique"""
    selector_name = "selector"
    """Name of the default selector"""
    selection_series_names = ("file_number",)
    """Name of the default things to use to construct the selector"""
    series_constructors = {
        "file_number": "_build_file_number_series",
        "selector": "_build_selector_series",
    }
    """Series which should be constructed from other series by the specified method
    and cached the first time they are looked up"""
    essential_series_names = None
    """Series which should always be present"""
    default_plotter = ValuePlotter
    default_exporter = CSVExporter

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
            c_ids (list of int): The id's of the measurement's Calibrations, if
                to be loaded (instead of given directly in calibration_list)
            calibration_list: The measurement's Calibrations
            m_ids (list of int): The id's of the component measurements, if to be
                loaded. None unless this is a combined measurement (typically
                corresponding to more than one file).
            component_measurements (list of Measurements): The measurements of which
                this measurement is a combination
            aliases (dict): Alternative names for DataSeries for versatile access
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

        self.plotter = plotter or self.__class__.default_plotter(measurement=self)
        self.exporter = exporter or self.__class__.default_exporter(measurement=self)
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
        #       DataSeries.unit_name) and as such should be moved to db.Savable
        #       see: https://github.com/ixdat/ixdat/pull/5#discussion_r565090372.
        #       Will be fixed with the table definition PR.
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
        measurement = technique_class(**obj_as_dict)
        return measurement

    @classmethod
    def read(cls, path_to_file, reader, **kwargs):
        """Return a Measurement object from parsing a file with the specified reader"""
        if isinstance(reader, str):
            # TODO: see if there isn't a way to put the import at the top of the module.
            #    see: https://github.com/ixdat/ixdat/pull/1#discussion_r546437471
            from .readers import READER_CLASSES

            reader = READER_CLASSES[reader]()
        obj = reader.read(path_to_file, **kwargs)  # TODO: take cls as kwarg

        if obj.__class__.essential_series_names:
            for series_name in obj.__class__.essential_series_names:
                try:
                    _ = obj[series_name]  # this also caches it.
                except SeriesNotFoundError:
                    raise SeriesNotFoundError(
                        f"{reader} loaded without {obj.__class__.__name__} "
                        f"essential series '{series_name}'"
                    )
        return obj

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
        for i, m in enumerate(self._component_measurements):
            if isinstance(m, PlaceHolderObject):
                # This is where we find objects from a Backend including MemoryBackend:
                self._component_measurements[i] = m.get_object()
        return self._component_measurements

    @property
    def m_ids(self):
        """List of the id's of a combined measurement's component measurements
        FIXME: m.id can be (backend, id) if it's not on the active backend.
            This is as of now necessary to find it if you're only given self.as_dict()
            see https://github.com/ixdat/ixdat/pull/11#discussion_r746632897
        """
        if not self._component_measurements:
            return None
        return [m.short_identity for m in self.component_measurements]

    @property
    def calibration_list(self):
        """List of calibrations (with placeholders filled)"""
        for i, c in enumerate(self._calibration_list):
            if isinstance(c, PlaceHolderObject):
                # This is where we find objects from a Backend including MemoryBackend:
                self._calibration_list[i] = c.get_object()
        return self._calibration_list

    @property
    def calibrations(self):
        """For overriding: List of calibrations with any needed manipulation done."""
        return self.calibration_list

    @property
    def c_ids(self):
        """List of the id's of the measurement's Calibrations
        FIXME: c.id can be (backend, id) if it's not on the active backend.
            This is as of now necessary to find it if you're only given self.as_dict()
             see https://github.com/ixdat/ixdat/pull/11#discussion_r746632897
        """
        return [c.short_identity for c in self.calibration_list]

    def add_calibration(self, calibration):
        self._calibration_list = [calibration] + self._calibration_list

    @property
    def series_list(self):
        """List of the DataSeries containing the measurement's data"""
        for i, s in enumerate(self._series_list):
            if isinstance(s, PlaceHolderObject):
                # This is where we find objects from a Backend including MemoryBackend:
                self._series_list[i] = s.get_object()
        return self._series_list

    @property
    def s_ids(self):
        """List of the id's of the measurement's DataSeries
        FIXME: m.id can be (backend, id) if it's not on the active backend.
            This is as of now necessary to find it if you're only given self.as_dict()
            see https://github.com/ixdat/ixdat/pull/11#discussion_r746632897
        """
        return [series.short_identity for series in self._series_list]

    @property
    def series_names(self):
        """Set of the names of the series in the measurement"""
        return set([series.name for series in self.series_list])

    @property
    def value_names(self):
        """Set of the names of the VSeries in the measurement's DataSeries"""
        return set([vseries.name for vseries in self.value_series])

    @property
    def time_names(self):
        """Set of the names of the VSeries in the measurement's DataSeries"""
        return set([tseries.name for tseries in self.time_series])

    @property
    def value_series(self):
        """Set of the VSeries in the measurement's DataSeries"""
        return [
            series for series in self.series_list if isinstance(series, ValueSeries)
        ]

    @property
    def time_series(self):
        """List of the TSeries in the measurement's DataSeries. NOT timeshifted!"""
        return [series for series in self.series_list if isinstance(series, TimeSeries)]

    @property
    def aliases(self):
        """Dictionary of {key: series_names} pointing to where desired raw data is

        TODO: get the possible aliases based on calibrations, etc, in here.
        """
        return self._aliases.copy()

    def get_series_names(self, key):
        """Return list: series names for key found by (recursive) lookup in aliases"""
        keys = [key] if key in self.series_names else []
        for k in self.aliases.get(key, []):
            keys += self.get_series_names(k)
        return keys

    def __getitem__(self, key):
        """Return the built measurement DataSeries with its name specified by key

        This method does the following:
        1. check if `key` is in in the cache. If so return the cached data series
        2. find or build the desired data series by the first possible of:
            A. Check if `key` corresponds to a method in `series_constructors`. If
                so, build the data series with that method.
            B. Check if the `calibration`'s `calibrate_series` returns a data series
                for `key` given the data in this measurement. (Note that the
                `calibration` will typically start with raw data looked C, below.)
            C. Generate a list of data series and append them:
                i. Check if `key` is in `aliases`. If so, append all the data series
                    returned for each key in `aliases[key]`.
                ii. Otherwise, check if there are data series in `series_list` that
                    have `key` as their `name`. If so, append them.
            D. Finally, check if the user is using a suffix.
                i. If `key` ends with "-y" or "-v", look it up with the suffix removed.
                ii. If `key` ends with "-x" or "-t", look up `key` with the suffix
                    removed and use instead the corresponding `tseries`.
        3. Cache and return the data series found or built in (2).

        Step (2) above, the searching step, is outsourced to the method
        `get_series(key)`.
        Notice that some calls of `__getitem__` can be recursive. For example, we
        suppose that a new `ECMeasurement` is read from a source that calls raw
        potential `Ewe/V`, and that this measurement is then calibrated:

        >>> ec_meas = Measurement.read(...)
        >>> ec_meas.aliases
        {..., 'raw_potential': ['Ewe/V'], ...}
        >>> ec_meas["raw_potential"]  # first lookup, explained below
        ValueSeries("Ewe/V", ...)
        >>> ec_meas.calibrate_RE(RE_vs_RHE=0.7)
        >>> ec_meas["potential"]   # second lookup, explained below
        ValueSeries("U_{RHE} / [V]", ...)

        - The first lookup, with `key="raw_potential"`, (1) checks for
        "raw_potential" in the cache, doesn't find it; then (2A) checks in
        `series_constructors`, doesn't find it; (2B) asks the calibration for
        "raw_potential" and doesn't get anything back; and finally (2Ci) checks
        `aliases` for raw potential where it finds that "raw_potential" is called
        "Ewe/V". Then it looks up again, this time with `key="Ewe/V"`, which it doesn't
        find in (1) the cache, (2A) `series_consturctors`, (2B) the calibration, or
        (2Ci) `aliases`, but does find in (2Cii) `series_list`. There is only one
        data series named "Ewe/V" so no appending is necessary, but it does ensure that
        the series has the measurement's `tstamp` before cache'ing and returning it.
        Now we're back in the original lookup, from which __getitem__ (3) caches
        the data series (which still has the name "Ewe/V") as "raw_potential" and
        returns it.
        - The second lookup, with `key="potential"`, (1) checks for "potential" in the
        cache, doesn't find it; then (2A) checks in `series_constructors`, doesn't find
        it; and then (2B) asks the calibration for "potential". The calibration knows
        that when asked for "potential" it should look for "raw_potential" and add
        `RE_vs_RHE`. So it does a lookup with `key="raw_potential"` and (1) finds it
        in the cache. The calibration does the math and returns a new data series for
        the calibrated potential, bringing us back to the original lookup. The data
        series returned by the calibration is then (3) cached and returned to the user.

        Note that, if the user had not looked up "raw_potential" before looking up
        "potential", "raw_potential" would not have been in the cache and the first
        lookup above would have been nested in the second.

        Args:
            key (str): The name of a DataSeries (see above)
        Raises:
            SeriesNotFoundError if none of the above lookups find the key.
        Side-effects:
            if key is not already in the cache, it gets added
        Returns:
            The (calibrated) (appended) dataseries for key with the right t=0.
        """
        # step 1
        if key in self._cached_series:
            return self._cached_series[key]
        # step 2
        series = self.get_series(key)
        # Finally, wherever we found the series, cache it and return it.
        # step 3.
        self._cached_series[key] = series
        return series

    def get_series(self, key):
        """Find or build the data series corresponding to key without direct cache'ing

        See more detailed documentation under `__getitem__`, for which this is a
        helper method. This method (A) looks for a method for `key` in the measurement's
        `series_constructors`; (B) requests its `calibration` for `key`; and if those
        fails appends the data series that either (Ci) are returned by looking up the
        key's `aliases` or (Cii) have `key` as their name; and finally (D) check if the
        user was using a key with a suffix.

        Args:
            key (str): The key to look up

        Returns DataSeries: the data series corresponding to key
        Raises SeriesNotFoundError if no series found for key
        """
        # A
        if key in self.series_constructors:
            return getattr(self, self.series_constructors[key])()
        # B
        for calibration in self.calibrations:
            series = calibration.calibrate_series(key, measurement=self)
            # ^ the calibration will call __getitem__ with the name of the
            #   corresponding raw data and return a new series with calibrated data
            #   if possible. Otherwise it will return None.
            if series:
                return series
        # C
        series_to_append = []
        if key in self.aliases:  # i
            # Then we'll look up the aliases instead and append them
            for k in self.aliases[key]:
                try:
                    series_to_append.append(self[k])
                except SeriesNotFoundError:
                    continue
        elif key in self.series_names:  # ii
            # Then we'll append any series matching the desired name
            series_to_append += [s for s in self.series_list if s.name == key]
        # If the key is something in the data, by now we have series to append.
        if series_to_append:
            # the following if's are to do as little extra manipulation as possible:
            if len(series_to_append) == 1:  # no appending needed
                if series_to_append[0].tstamp == self.tstamp:  # no time-shifting needed
                    return series_to_append[0]
                return time_shifted(series_to_append[0], tstamp=self.tstamp)
            return append_series(series_to_append, name=key, tstamp=self.tstamp)
        # D
        if key.endswith("-t") or key.endswith("-x"):
            return self[key[:-2]].tseries
        if key.endswith("-v") or key.endswith("-y"):
            return self[key[:-2]]

        raise SeriesNotFoundError(f"{self} does not contain '{key}'")

    def clear_cache(self):
        """Clear the cache so derived series are constructed again with updated info"""
        self._cached_series = {}

    def grab(self, item, tspan=None):
        """Return the time and value vectors for a given VSeries name cut by tspan"""
        series = self[item]
        v = series.v
        t = series.t
        if tspan:
            mask = np.logical_and(tspan[0] < t, t < tspan[-1])
            t, v = t[mask], v[mask]
        return t, v

    def grab_for_t(self, item, t):
        """Return a numpy array with the value of item interpolated to time t"""
        series = self[item]
        v_0 = series.v
        t_0 = series.t
        v = np.interp(t, t_0, v_0)
        return v

    @property
    def t(self):
        return self[self.control_series_name].t

    @property
    def t_name(self):
        return self[self.control_series_name].tseries.name

    def _build_file_number_series(self):
        """Build a `file_number` series based on component measurements times."""
        series_to_append = []
        for i, m in enumerate(self.component_measurements):
            if (
                self.control_technique_name
                and not m.technique == self.control_technique_name
            ):
                continue
            if not self.control_series_name:
                tseries = m.time_series[0]
            else:
                try:
                    tseries = m[self.control_series_name].tseries
                except KeyError:
                    continue
            series_to_append.append(
                ConstantValue(name="file_number", unit_name="", data=i, tseries=tseries)
            )
        return append_series(series_to_append, name="file_number", tstamp=self.tstamp)

    def _build_selector_series(
        self, selector_string=None, col_list=None, extra_col_list=None
    ):
        """Build a `selector` series which demarcates the data.

        The `selector` is a series which can be used to conveniently and powerfully
        grab sections of the data. It is built up from less powerful demarcation series
        in the raw data (like `cycle_number`, `step_number`, `loop_number`, etc) and
        `file_number` by counting the cumulative changes in those series.
        See slide 3 of:
        https://www.dropbox.com/s/sjxzr52fw8yml5k/21E18_DWS3_cont.pptx?dl=0

        Args:
            selector_string (str): The name to use for the selector series
            col_list (list): The list of demarcation series. The demarcation series have
                to have the same tseries, which should be the one pointed to by the
                meausrement's `control_series_name`.
            extra_col_list (list): Extra demarcation series to include if needed.
        """
        # the name of the selector series:
        selector_string = selector_string or self.selector_name
        # a vector that will be True at the points where a series changes:
        changes = np.tile(False, self.t.shape)
        # the names of the series which help demarcate the data
        col_list = col_list or self.selection_series_names
        if extra_col_list:
            col_list += extra_col_list
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
            # a vector which is shifted one.
            last_value = np.append(values[0], values[:-1])
            # comparing value and last_value shows where in the vector changes occur:
            changes = np.logical_or(changes, last_value != values)
        # taking the cumsum makes a vector that increases 1 each time one of the
        #   original demarcation vector changes
        selector_data = np.cumsum(changes)
        selector_series = ValueSeries(
            name=selector_string,
            unit_name="",
            data=selector_data,
            tseries=self[self.control_series_name].tseries,
        )
        return selector_series

    @property
    def selector(self):
        return self[self.selector_name]

    @property
    def data_cols(self):
        """Return a set of the names of all of the measurement's VSeries and TSeries"""
        return set([s.name for s in (self.value_series + self.time_series)])

    def get_original_m_ids_of_series(self, series):
        """Return a list of id's of component measurements to which `series` belongs."""
        m_id_list = []
        for m in self.component_measurements:
            if series.short_identity in m.s_ids:
                # FIXME: the whole id vs short_identity issue
                #   see https://github.com/ixdat/ixdat/pull/11#discussion_r746632897
                m_id_list.append(m.id)
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

    def cut(self, tspan, t_zero=None):
        """Return a new measurement with the data in the given time interval

        Args:
            tspan (iter of float): The time interval to use, relative to self.tstamp
                tspan[0] is the start time of the interval, and tspan[-1] is the end
                time of the interval. Using tspan[-1] means you can directly use a
                long time vector that you have at hand to describe the time interval
                you're looking for.
            t_zero (float or str): The time in the measurement to set to t=0. If a
                float, it is interpreted as wrt the original tstamp. String options
                include "start", which puts t=0 at the start of the cut interval.
        """
        # Start with self's dictionary representation, but
        # we don't want original series (s_ids) or component_measurements (m_ids):
        obj_as_dict = self.as_dict(exclude=["s_ids", "m_ids"])

        # first, cut the series list:
        new_series_list = []
        time_cutting_stuff = {}  # {tseries_id: (mask, new_tseries)}
        for series in self.series_list:
            try:
                tseries = series.tseries
                if tseries is None:
                    raise AttributeError
            except AttributeError:  # series independent of time are uneffected by cut
                new_series_list.append(series)
            else:
                t_identity = tseries.full_identity

                if t_identity in time_cutting_stuff:
                    mask, new_tseries = time_cutting_stuff[t_identity]
                else:
                    t = tseries.t + tseries.tstamp - self.tstamp
                    mask = np.logical_and(tspan[0] <= t, t <= tspan[-1])
                    new_tseries = TimeSeries(
                        name=tseries.name,
                        unit_name=tseries.unit_name,
                        tstamp=tseries.tstamp,
                        data=tseries.data[mask],
                    )
                    time_cutting_stuff[t_identity] = (mask, new_tseries)
                if True not in mask:
                    continue
                if False not in mask:
                    new_series_list.append(series)
                elif series.full_identity == t_identity:
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

        # then cut the component measurements.
        new_component_measurements = []
        for m in self._component_measurements:
            # FIXME: This is perhaps overkill, to make new cut component measurements,
            #    as it duplicates data (a big no)... especially bad because
            #    new_measurement.save() saves them.
            #    The step is here in order for file_number to get built correctly.
            dt = m.tstamp - self.tstamp
            tspan_m = [tspan[0] - dt, tspan[1] - dt]
            if m.tspan[-1] < tspan_m[0] or tspan_m[-1] < m.tspan[0]:
                continue
            new_component_measurements.append(m.cut(tspan_m))
        obj_as_dict["component_measurements"] = new_component_measurements

        new_measurement = self.__class__.from_dict(obj_as_dict)
        if t_zero:
            if t_zero == "start":
                new_measurement.tstamp += tspan[0]
            else:
                new_measurement.tstamp += t_zero
        return new_measurement

    def multicut(self, tspans):
        """Return a selection of the measurement including each of the given tspans"""
        # go through the tspans, cuting the measurement and appending the results
        new_measurement = None
        for tspan in tspans:
            if new_measurement:
                new_measurement = new_measurement + self.cut(tspan)
            else:
                new_measurement = self.cut(tspan)
        return new_measurement

    def select_value(self, *args, **kwargs):
        """Return a selection of the measurement where a criterion is matched.

        Specifically, this method returns a new Measurement where the time(s) returned
        are those where the values match the provided criteria, i.e. the part of the
        measurement where `self[series_name] == value`

        Can only take one arg or kwarg!
        The `series_name` is `self.selector_name` if given an argument without keyword.
        If given a keyword argument, the kyword is the name of the series to select on.
        Either way the argument is the `value` to be selected for.

        The method finds all time intervals for which `self[series_name] == value`
        It then cuts the measurement according to each time interval and adds these
        segments together.
        TODO: This can maybe be done better, i.e. without chopping series.
        TODO: Some way of less than and greater than kwargs.
            Ideally you should be able to say e.g., `select(cycle=1, 0.5<potential<1)`
            But this is hard,
            see: https://github.com/ixdat/ixdat/pull/11#discussion_r677272239
        """
        if len(args) + len(kwargs) != 1:
            raise BuildError("Need exactly 1 arg. Use `select_values` for more.")
        if args:
            if not self.selector_name:
                raise BuildError(
                    f"{self} does not have a default selection string "
                    f"(Measurement.sel_str), and so selection only works with kwargs."
                )
            kwargs[self.selector_name] = args[0]

        ((series_name, value),) = kwargs.items()

        # The time and values of the series to be selected on:
        t, v = self.grab(series_name)
        # This mask is true everywhere on `t` that the condition is met:
        mask = v == value  # linter doesn't realize this is a np array

        # Now we have to convert that to timespans on which `t` is met. This means
        #  finding the start and finish times of the intervals on which mask is True.
        #  this is done with a helper function:
        tspans = get_tspans_from_mask(t, mask)

        # now we go through the tspans, cuting the measurement and appending the results:
        return self.multicut(tspans)

    def select_values(self, *args, **kwargs):
        """Return a selection of the measurement based on one or several criteria

        Specifically, this method returns a new Measurement where the time(s) returned
        are those where the values match the provided criteria, i.e. the part of the
        measurement where `self[series_name] == value`

        TODO: Testing and documentation with examples like those suggested here:
            https://github.com/ixdat/ixdat/pull/11#discussion_r677324246

        Any series can be selected for using the series name as a key-word. Arguments
        can be single acceptable values or lists of acceptable values. In the latter
        case, each acceptable value is selected for on its own and the resulting
        measurements added together.
        # FIXME: That is sloppy because it mutliplies the number of DataSeries
            containing the same amount of data.
        Arguments without key-word are considered valid values of the default selector,
        which is named by `self.selelector_name`. Multiple criteria are
        applied sequentially, i.e. you get the intersection of satisfying parts.

        Args:
            args (tuple): Argument(s) given without keyword are understood as acceptable
                value(s) for the default selector (that named by self.sel_str)
            kwargs (dict): Each key-word arguments is understood as the name
                of a series and its acceptable value(s).
        """
        if args:
            if not self.selector_name:
                raise BuildError(
                    f"{self} does not have a default selection string "
                    f"(Measurement.sel_str), and so selection only works with kwargs."
                )
            flat_args = []
            for arg in args:
                if hasattr(arg, "__iter__"):
                    flat_args += list(arg)
                else:
                    flat_args.append(arg)
            if self.selector_name in kwargs:
                raise BuildError(
                    "Don't call select values with both arguments and "
                    "'{self.selector_name}' as a key-word argument"
                )
            kwargs[self.selector_name] = flat_args

        t = self.t
        mask = np.tile(np.array([True]), t.shape)
        for series_name, allowed_values in kwargs.items():
            if not hasattr(allowed_values, "__iter__"):
                allowed_values = [allowed_values]
            v = self.grab_for_t(series_name, t)
            submask = np.tile(np.array([False]), t.shape)
            for allowed_value in allowed_values:
                submask = np.logical_or(submask, v == allowed_value)
            mask = np.logical_and(mask, submask)

        tspans = get_tspans_from_mask(t, mask)

        return self.multicut(tspans)

    def select(self, *args, tspan=None, **kwargs):
        """`cut` (with tspan) and `select_values` (with *args and/or **kwargs).

        These all work:
        - `meas.select_values(1, 2)`
        - `meas.select_values(tspan=[200, 300])`
        - `meas.select_values(range(10))`
        - `meas.select_values(cycle=4)`
        - `meas.select_values(1, range(5, 20), file_number=1, tspan=[1000, 2000])`
        """
        new_measurement = self
        if tspan:
            new_measurement = new_measurement.cut(tspan=tspan)
        if args or kwargs:
            new_measurement = new_measurement.select_values(*args, **kwargs)
        return new_measurement

    def copy(self):
        """Make a copy of the Measurement via its dictionary representation"""
        return self.__class__.from_dict(self.as_dict())

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

        # breakpoint()
        new_series_list = list(set(self.series_list + other.series_list))
        new_component_measurements = list(
            set(
                (self.component_measurements or [self])
                + (other.component_measurements or [other])
            )
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
        obj_as_dict = self.as_dict()
        obj_as_dict.update(
            name=new_name,
            series_list=new_series_list,
            component_measurements=new_component_measurements,
            calibration_list=new_calibration_list,
            aliases=new_aliases,
        )
        # don't want the original calibrations, component measurements, or series:
        del obj_as_dict["c_ids"]
        del obj_as_dict["m_ids"]
        del obj_as_dict["s_ids"]
        return cls.from_dict(obj_as_dict)

    def __eq__(self, other):
        """Return whether this object is equal to `other`"""
        # IMPORTANT, the order of the checks, in general, and in particular of the
        # property names further down, is intentional to keep cheap result determining
        # comparisons first and expensive ones last, for performance reasons
        if self is other:
            return True

        if not np.isclose(self.tstamp, other.tstamp):
            return False

        if not dict_is_close(self.metadata, other.metadata):
            return False

        # TODO This list could probably be replaced completely or in part with some of
        # the DB-related class attributes like column_attrs, but this shouldn't be done
        # before those attributes have stabilized
        for property_name in (
            "__class__",
            "name",
            "technique",
            "sample",
            "lablog",
            "aliases",
            "calibration_list",
            "series_list",
            "component_measurements",  # Assuming these can't reference each other
        ):
            if getattr(self, property_name) != getattr(other, property_name):
                return False
        return True

    # This is necessary, because overriding __eq__ means that __hash__ is set to None
    # https://docs.python.org/3/reference/datamodel.html#object.__hash__
    __hash__ = DataSeries.__hash__


class Calibration(Savable):
    """Base class for calibrations."""

    table_name = "calibration"
    column_attrs = {
        "name",
        "technique",
        "tstamp",
    }

    def __init__(self, name=None, technique=None, tstamp=None, measurement=None):
        """Initiate a Calibration

        Args:
            name (str): The name of the calibration
            technique (str): The technique of the calibration
            tstamp (float): The time at which the calibration took place or is valid
            measurement (Measurement): Optional. A measurement to calibrate by default.
        """
        super().__init__()
        self.name = name or f"{self.__class__.__name__}({measurement})"
        self.technique = technique
        self.tstamp = tstamp or (measurement.tstamp if measurement else None)
        self.measurement = measurement

    def __eq__(self, _):
        """Return whether the calibration is equal to another

        MUST be overwritten in subclasses

        """
        raise NotImplementedError

    @classmethod
    def from_dict(cls, obj_as_dict):
        """Return an object of the Calibration class of the right technique

        Args:
              obj_as_dict (dict): The full serializaiton (rows from table and aux
                tables) of the measurement. obj_as_dict["technique"] specifies the
                technique class to use, from TECHNIQUE_CLASSES
        """
        # TODO: see if there isn't a way to put the import at the top of the module.
        #    see: https://github.com/ixdat/ixdat/pull/1#discussion_r546437410
        from .techniques import CALIBRATION_CLASSES

        if obj_as_dict["technique"] in CALIBRATION_CLASSES:
            calibration_class = CALIBRATION_CLASSES[obj_as_dict["technique"]]
        else:
            calibration_class = cls
        try:
            measurement = calibration_class(**obj_as_dict)
        except Exception:
            raise
        return measurement

    def calibrate_series(self, key, measurement=None):
        """This should be overwritten in real calibration classes.

        FIXME: Add more documentation about how to write this in inheriting classes.
        """
        raise NotImplementedError
