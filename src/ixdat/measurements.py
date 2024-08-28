"""This module defines the Measurement class, the central data structure of ixdat

An ixdat Measurement is a collection of references to DataSeries and the metadata needed
to combine them, i.e. "build" the combined dataset. It has a number of general methods
to visualize and analyze the combined dataset. Measurement is also the base class for a
number of technique-specific Measurement-derived classes.

A Measurement will typically be accompanied by one or more Calculator. This module
also defines the base class for Calculator, while technique-specific Calculator
classes will be defined in the corresponding module in ./techniques/
"""

import warnings
import json
import numpy as np
from .db import Saveable, PlaceHolderObject, fill_object_list
from .data_series import (
    DataSeries,
    TimeSeries,
    ValueSeries,
    append_series,
    time_shifted,
    get_tspans_from_mask,
)
from .projects.samples import Sample
from .projects.lablogs import LabLog
from .exporters.csv_exporter import CSVExporter
from .plotters.value_plotter import ValuePlotter
from .exceptions import BuildError, SeriesNotFoundError, TechniqueError, ReadError
from .tools import tstamp_to_string


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
        "measurement_calculators": ("calculators", "c_ids"),
        "measurement_series": ("data_series", "s_ids"),
    }
    child_attrs = ["component_measurements", "calculator_list", "series_list"]
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
    essential_series_names = None
    """Series which should always be present"""
    default_plotter = ValuePlotter
    default_exporter = CSVExporter
    default_calibration = None
    built_in_calculator_types = ["indexer"]
    background_calculator_types = []

    def __init__(
        self,
        name,
        technique=None,
        metadata=None,
        s_ids=None,
        series_list=None,
        c_ids=None,
        calculator_list=None,
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
            c_ids (list of int): The id's of the measurement's calculators, if
                to be loaded (instead of given directly in calculator_list)
            calculator_list: The measurement's calculators
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
        self._calculator_list = fill_object_list(calculator_list, c_ids, cls=Calculator)

        self._temp_calculator_list = None
        self._calculator_dict = None
        self._built_in_calculators = None

        self._tstamp = tstamp

        self._cached_series = {}
        self._aliases = aliases or {}

        self.plotter = plotter or self.__class__.default_plotter(measurement=self)
        self.exporter = exporter or self.__class__.default_exporter(measurement=self)
        # defining these methods here gets them the right docstrings :D
        self.plot_measurement = self.plotter.plot_measurement
        self.plot = self.plotter.plot_measurement
        self.export = self.exporter.export
        # TODO: ... but we need to think a bit more about how to most elegantly and
        #    dynamically choose plotters (Nice idea from Anna:
        #    https://github.com/ixdat/ixdat/issues/32)

    def __str__(self):
        """Return string representation"""
        tseries_to_valueseries = {}
        for series in self.series_list:
            if isinstance(series, TimeSeries):
                if series not in tseries_to_valueseries:
                    tseries_to_valueseries[series] = []
            else:
                if series.tseries in tseries_to_valueseries:
                    tseries_to_valueseries[series.tseries].append(series)
                else:
                    tseries_to_valueseries[series.tseries] = [series]

        out = []
        for tseries, value_serieses in tseries_to_valueseries.items():
            out.append("┏ " + str(tseries))
            for n, value_series in enumerate(value_serieses):
                if n == len(value_serieses) - 1:
                    out.append("┗━ " + str(value_series))
                else:
                    out.append("┣━ " + str(value_series))

        calc_out = []
        for (n, (name, calculator)) in enumerate(self.calculator_dict.items()):
            if n == len(self.calculator_dict) - 1:
                calc_out.append("┗━ " + str(calculator))
            else:
                calc_out.append("┣━ " + str(calculator))

        return (
            f"{self.__class__.__name__} '{self.name}' with {len(self.series_list)} "
            "series\n\n"
            "Series list:\n"
            + "\n".join(out)
            + "\n\nCalculators:\n"
            + "\n".join(calc_out)
        )

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
        #       see: https://github.com/ixdat/ixdat/pull/5#discussion_r565090372.
        #       Will be fixed with the table definition PR.
        objects_saved_as_their_name = ["sample"]
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
            technique_class = cls
        try:
            measurement = technique_class(**obj_as_dict)
        except TypeError as e:
            raise TechniqueError(
                "ixdat ran into an error while trying to set up an object of type "
                f"{technique_class}. This usually happens when ixdat wasn't able"
                "to correctly determine the measurement technique. Consider"
                "passing the `technique` argument into the read() function. \n"
                "For a list of available techniques use: \n "
                ">>> from ixdat.techniques import TECHNIQUE_CLASSES\n"
                ">>> print(TECHNIQUE_CLASSES.keys())\n"
                f"{e}"
            )
            raise

        return measurement

    @classmethod
    def read(cls, path_to_file, reader=None, **kwargs):
        """Return a Measurement object from parsing a file with the specified reader

        Args:
            path_to_file (Path or str): The path to the file to read
            reader (str or Reader class): The (name of the) reader to read the file
                with. If not specified, ixdat will try to determine the reader from the
                file suffix.
            kwargs: key-word arguments are passed on to the reader's read() method.
        """
        if not reader:
            # Check if there is a default reader based on the file's suffix
            from .readers.reading_tools import get_default_reader_name

            reader = get_default_reader_name(path_to_file)
            if not reader:
                raise ValueError(
                    f"There is no default reader for files of the type {path_to_file}. "
                    "Please specify a reader to read this file."
                )

        if isinstance(reader, str):
            # TODO: see if there isn't a way to put the import at the top of the module.
            #    see: https://github.com/ixdat/ixdat/pull/1#discussion_r546437471
            from .readers import READER_CLASSES

            reader = READER_CLASSES[reader]()
        obj = reader.read(path_to_file, cls=cls, **kwargs)

        if getattr(obj.__class__, "essential_series_names", None):
            for series_name in obj.__class__.essential_series_names:
                try:
                    _ = obj[series_name]  # this also caches it.
                except SeriesNotFoundError:
                    raise SeriesNotFoundError(
                        f"{reader} loaded without {obj.__class__.__name__} "
                        f"essential series '{series_name}'"
                    )
        return obj

    @classmethod
    def read_url(cls, url, reader=None, **kwargs):
        """Read a url (via a temporary file) using the specified reader"""
        from .readers.reading_tools import url_to_file

        path_to_temp_file = url_to_file(url)
        measurement = cls.read(path_to_temp_file, reader=reader, **kwargs)
        path_to_temp_file.unlink()
        return measurement

    @classmethod
    def read_set(
        cls,
        path_to_file_start=True,
        part=None,
        suffix=None,
        file_list=None,
        reader=None,
        **kwargs,
    ):
        """Read and append a set of files.

        Args:
            path_to_file_start (Path or str): The path to the files to read including
                the shared start of the file name: `Path(path_to_file).parent` is
                interpreted as the folder where the file are.
                `Path(path_to_file).name` is interpreted as the shared start of the files
                to be appended.
                Alternatively, path_to_file_start can be a folder, in which case all
                files in that folder (with the specified suffix) are included.
            part (Path or str): A path where the folder is the folder containing data
                and the name is a part of the name of each of the files to be read and
                combined.
            suffix (str): If a suffix is given, only files with the specified ending are
                added to the file list
            file_list (list of Path): As an alternative to path_to_file_start, the
                exact files to append can be specified in a list
            reader (str or Reader class): The (name of the) reader to read the files with
            kwargs: Key-word arguments are passed via cls.read() to the reader's read()
                method, AND to cls.from_component_measurements()
        """
        from .readers.reading_tools import get_file_list

        file_list = file_list or get_file_list(path_to_file_start, part, suffix)
        if not file_list:
            raise ReadError(
                "No files found! Please check that there are files satisfying:\n"
                f"path_to_file_start={path_to_file_start}, part={part}, suffix={suffix}"
            )
        component_measurements = [
            cls.read(f, reader=reader, **kwargs) for f in file_list
        ]
        measurement = None
        for meas in component_measurements:
            measurement = measurement + meas if measurement else meas
        return measurement

    @classmethod
    def from_component_measurements(
        cls, component_measurements, keep_originals=True, sorted=True, **kwargs
    ):
        """Return a measurement with the data contained in the component measurements

        TODO: This function "builds" the resulting measurement, i.e. it appends series
            of the same name rather than keeping all the original copies. This should be
            made more explicit, and a `build()` method should take over some of the work.

        Args:
            component_measurements (list of Measurement)
            keep_originals: Whether to keep a list of component_measurements referenced.
                This may result in redundant numpy arrays in RAM.
            sorted (bool): Whether to sort the series according to time
            kwargs: key-word arguments are added to the dictionary for cls.from_dict()

        Returns cls: the combined measurement.
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
                if sorted:
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
                    if sorted:
                        s_as_dict["data"] = s_as_dict["data"][sort_indeces[tseries.name]]
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
                        sorted=sorted,
                    )
                series_list.append(vseries)

        # Finally, add the series to the dictionary representation and return the object
        obj_as_dict["series_list"] = series_list
        return cls.from_dict(obj_as_dict)

    @property
    def tstamp(self):
        """Float: The unix epoch time used by the measurement as t=0"""
        return self._tstamp

    @tstamp.setter
    def tstamp(self, tstamp):
        # Resetting the tstamp needs to clear the cache, so series are returned wrt the
        # new timestamp.
        self.clear_cache()
        self._tstamp = tstamp

    @property
    def yyMdd(self):
        return tstamp_to_string(self.tstamp, string_format="native_date")

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
    def calculator_list(self):
        """List of calculators (with placeholders filled)"""
        for i, c in enumerate(self._calculator_list):
            if isinstance(c, PlaceHolderObject):
                # This is where we find objects from a Backend including MemoryBackend:
                self._calculator_list[i] = c.get_object()
        if self._temp_calculator_list:
            return self._temp_calculator_list
        return self._calculator_list

    @property
    def calculator_dict(self):
        if self._calculator_dict is None:
            self.consolidate_calculators()
        return self._calculator_dict

    def consolidate_calculators(self):
        """Dictionary of calculators, consolidated if needed to one per type."""
        calculators = self.built_in_calculators.copy()  # <-- result of a tough debug :)
        for cal in self.calculator_list:
            name = cal.calculator_type
            if name in calculators:
                calculators[name] = calculators[name] + cal
            else:
                calculators[name] = cal
        self._calculator_dict = calculators
        return self._calculator_dict

    @property
    def built_in_calculators(self):
        if not self._built_in_calculators:
            from .calculators import CALCULATOR_CLASSES

            self._built_in_calculators = {
                key: CALCULATOR_CLASSES[key]() for key in self.built_in_calculator_types
            }
        return self._built_in_calculators

    @property
    def c_ids(self):
        """List of the id's of the measurement's calculators
        FIXME: c.id can be (backend, id) if it's not on the active backend.
            This is as of now necessary to find it if you're only given self.as_dict()
             see https://github.com/ixdat/ixdat/pull/11#discussion_r746632897
        """
        return [c.short_identity for c in self.calculator_list]

    def add_calculator(self, calculator):
        redundants = calculator.available_series_names.intersection(
            self.available_calculated_series
        )
        if redundants:
            warnings.warn(
                f"adding {calculator} to {self!r} "
                f"results in redundant calculators for: {redundants}"
            )

        self._calculator_list = [calculator] + self._calculator_list
        ctype = calculator.calculator_type
        if not self._calculator_dict:
            self.consolidate_calculators()
        elif ctype in self._calculator_dict:
            self._calculator_dict[ctype] = self.calculator_dict[ctype] + calculator
        else:
            self._calculator_dict[ctype] = calculator
        self.clear_cache()

    @property
    def available_calculated_series(self):
        calculated_series_names = set([])
        for ctype, cal in self.calculator_dict.items():
            calculated_series_names = calculated_series_names.union(
                cal.available_series_names
            )
        return calculated_series_names

    def calibrate(self, *args, **kwargs):
        """Add a calculator of the Measurement's default calculator type

        The calculator class is determined by the measurement's `technique`.
        *args and **kwargs are passed to the calculator class's `__init__`.

        Raises:
            TechniqueError if no calculator class for the measurement's technique
        """
        new_calculator = self.default_calibration(*args, **kwargs)
        self.add_calculator(new_calculator)
        self.clear_cache()
        return new_calculator

    def rebuild_selector(self, selector_name=None, columns=None, extra_columns=None):
        from .calculators import CALCULATOR_CLASSES

        indexer = CALCULATOR_CLASSES["indexer"](
            measurement=self,
            selector_name=selector_name,
            columns=columns,
            extra_columns=extra_columns,
        )
        self.add_calculator(indexer)

    # Note: Not all ´Calculator´s are ´Calibraiton´s There are also, e.g., `Filter`s and
    # `Background`s. One could, for completeness, also implement the methods
    # `Measurement.background_correct` and `Measurement.filter` which, like the above,
    # just add an object of some default calculator class to the measurement's
    # `calculator_list`. However, perhaps one global default way to add a calculator to
    # a measurement's calculator list is enough, and `calibrate` is a more natural
    # and broad English verb than the others.

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
        """List of the VSeries in the measurement's DataSeries"""
        return [series for series in self.series_list if isinstance(series, ValueSeries)]

    @property
    def time_series(self):
        """List of the TSeries in the measurement's DataSeries. NOT timeshifted!"""
        return [series for series in self.series_list if isinstance(series, TimeSeries)]

    @property
    def aliases(self):
        """Dictionary of {key: series_names} pointing to where desired raw data is

        TODO: get the possible aliases based on calculators, etc, in here.
        """
        return self._aliases.copy()

    @property
    def reverse_aliases(self):
        """{series_name: standard_names} indicating how raw data can be accessed"""
        rev_aliases = {}
        for name, other_names in self.aliases.items():
            for other_name in other_names:
                if other_name in rev_aliases:
                    rev_aliases[other_name].append(name)
                else:
                    rev_aliases[other_name] = [name]
        return rev_aliases

    def get_series_names(self, key):
        """Return list: series names for key found by (recursive) lookup in aliases"""
        keys = [key] if key in self.series_names else []
        for k in self.aliases.get(key, []):
            keys += self.get_series_names(k)
        return keys

    def __getitem__(self, key):
        """Return the built measurement DataSeries with its name specified by key

        This method does the following:
        0. Check that the key is a string. If a technique supports lookup of other
           types, the technique class should implement that in its `__getitem__`
           before calling `super().__getitem__`.
        1. check if `key` is in in the cache. If so return the cached data series
        2. find or build the desired data series by the first possible of:
            A. Check if `key` corresponds to a method in `series_constructors`. If
                so, build the data series with that method.
            B. Check if the `calculator`'s `calculate_series` returns a data series
                for `key` given the data in this measurement. (Note that the
                `calculator` will typically start with raw data looked C, below.)
               If the key has the suffix "-raw", skip the calculator lookup; instead
               remove the suffix and continue to C.
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
        >>> ec_meas["potential"]      # second lookup, explained below
        ValueSeries("U_{RHE} / [V]", ...)

        - The first lookup, with `key="raw_potential"`, (1) checks for
        "raw_potential" in the cache, doesn't find it; then (2A) checks in
        `series_constructors`, doesn't find it; (2B) asks the calculator for
        "raw_potential" and doesn't get anything back; and finally (2Ci) checks
        `aliases` for raw potential where it finds that "raw_potential" is called
        "Ewe/V". Then it looks up again, this time with `key="Ewe/V"`, which it doesn't
        find in (1) the cache, (2A) `series_consturctors`, (2B) the calculator, or
        (2Ci) `aliases`, but does find in (2Cii) `series_list`. There is only one
        data series named "Ewe/V" so no appending is necessary, but it does ensure that
        the series has the measurement's `tstamp` before cache'ing and returning it.
        Now we're back in the original lookup, from which __getitem__ (3) caches
        the data series (which still has the name "Ewe/V") as "raw_potential" and
        returns it.
        - The second lookup, with `key="potential"`, (1) checks for "potential" in
        the cache, doesn't find it; then (2A) checks in `series_constructors`,
        doesn't find it; and then (2B) asks the calculator for "potential". The
        calculator knows that when asked for "potential" it should look for
        "raw_potential" and add `RE_vs_RHE`. So it does a lookup with
        `key="raw_potential"` and (1) finds it in the cache. The calculator does
        the math and returns a new data series for the calculated potential, bringing
        us back to the original lookup. The data series returned by the
        calculator is then (3) cached and returned to the user.

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
            The (calculated) (appended) dataseries for key with the right t=0.
        """
        # step 0
        if not isinstance(key, str):
            message = f"Invalid lookup for {type(self)} object: {key}."
            message += f" The key type was {type(key)}. Expected a string."
            if isinstance(key, int):
                message += (
                    " Note: Integer lookup is possible for SpectroMeasurement and"
                    " CyclicVoltammogram objects. If you expected a measurement"
                    " containing spectra or index-able cycles,"
                    " please check your file reading."
                )
            raise TypeError(message)

        # step 1
        if key in self._cached_series:
            return self._cached_series[key]
        # step 2
        series = self.get_series(key)
        # Finally, wherever we found the series, cache it and return it.
        # step 3.
        self._cache_series(key, series)
        return series

    def _cache_series(self, key, series):
        """Cache `series` such that it can be looked up with its name or with `key`."""
        self._cached_series[key] = series  # now it can be looked up with by `key`
        # If the name of the series is not `key`, we can get in a situation where
        # looking up the series name raises a SeriesNotFoundError. To avoid this
        # problematic situation, we check if it can be looked up, and if not,
        # add it also to aliases, now under `series.name`
        try:
            _ = self[series.name]
        except SeriesNotFoundError:
            self._aliases[series.name] = [key]

    def get_series(self, key):
        """Find or build the data series corresponding to key without direct cache'ing

        See more detailed documentation under `__getitem__`, for which this is a
        helper method. This method (A) looks for a method for `key` in the measurement's
        `series_constructors`; (B) requests its `calculator` for `key`; and if those
        fail appends the data series that either (Ci) are returned by looking up the
        key's `aliases` or (Cii) have `key` as their name; and finally (D) check if the
        user was using a key with a suffix.

        Args:
            key (str): The key to look up

        Returns DataSeries: the data series corresponding to key
        Raises SeriesNotFoundError if no series found for key
        """
        # B
        if key.endswith("-raw"):
            # A "-raw" suffix means to skip the calculators
            key = key.removesuffix("-raw")
        else:
            for calculator in self.calculator_dict.values():
                if key in calculator.available_series_names:
                    series = calculator.calculate_series(key, measurement=self)
                    if series:
                        return series
                    warnings.warn(
                        f"The Calulator {calculator} inscludes {key} in its"
                        f" `available_series_names`, yet returns `None` when asked to"
                        " calculate it."
                    )
        # C
        series_to_append = []
        if key in self.series_names:  # ii
            # Then we'll append any series matching the desired name
            series_to_append += [s for s in self.series_list if s.name == key]
        if key in self.aliases:  # i
            # Then we'll look up the aliases instead and append them
            for k in self.aliases[key]:
                if k == key:  # this would result in infinite recursion.
                    warnings.warn(
                        f"{self!r} has {key} in its aliases for {key}:\n"
                        f"self.aliases['{key}'] = {self.aliases[key]}"
                    )
                    continue
                try:
                    series_to_append.append(self[k])
                except SeriesNotFoundError:
                    continue
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

        raise SeriesNotFoundError(f"{self!r} does not contain '{key}'")

    def replace_series(self, series_name, new_series=None):
        """Remove an existing series, add a series to the measurement, or both.

        FIXME: This will not appear to change the series for the user if the
            measurement's calculator returns something for ´series_name´, since
            __getitem__ asks the calculator before looking in series_list.

        Args:
            series_name (str): The name of a series. If the measurement has (raw) data
                series with this name, cached series with this name, and/or aliases for
                this name, they will be removed.
            new_series (DataSeries): Optional new series to append to the measurement's
                series_list. To sanity check, it must have ´series_name´ as its ´name´.
        """
        if new_series and not series_name == new_series.name:
            raise TypeError(
                f"Cannot replace {series_name} in {self!r} with {new_series}. "
                f"Names must agree."
            )
        if series_name in self._cached_series:
            del self._cached_series[series_name]
        if series_name in self._aliases:
            del self._aliases[series_name]
        new_series_list = [s for s in self.series_list if not s.name == series_name]
        if new_series:
            new_series_list.append(new_series)
        self._series_list = new_series_list

    def clear_cache(self):
        """Clear the cache so derived series are constructed again with updated info"""
        self._cached_series = {}

    def correct_data(self, value_name, new_data):
        """Replace the old data for ´value_name´ (str) with ´new_data` (np array)"""
        old_vseries = self[value_name]
        new_vseries = ValueSeries(
            name=value_name,
            unit_name=old_vseries.unit_name,
            data=new_data,
            tseries=old_vseries.tseries,
        )
        self.replace_series(value_name, new_vseries)

    def grab(
        self,
        item,
        tspan=None,
        include_endpoints=False,
        tspan_bg=None,
        remove_background=None,
        calculator_list=None,
    ):
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
            t, v = measurement.grab("potential", tspan=[0, 100])

        Args:
            item (str): The name of the DataSeries to grab data for
                TODO: Should this be called "name" or "key" instead? And/or, should
                   the argument to __getitem__ be called "item" instead of "key"?
            tspan (iter of float): Defines the timespan with its first and last values.
                Optional. By default the entire time of the measurement is included.
            include_endpoints (bool): Whether to add a points at t = tspan[0] and
                t = tspan[-1] to the data returned. This makes trapezoidal integration
                less dependent on the time resolution. Default is False.
            tspan_bg (iterable): Optional. A timespan defining when `item` is at its
                baseline level. The average value of `item` in this interval will be
                subtracted from the values returned.
            remove_background (boolean): Whether to subtract a pre-set background, if
                available. This is True by default. If set to False, it suppresses
                background calculators, taken to be those specified by the class
                attribute `background_calculator_types`
            calculator_list (list of Calculators): A list of ixdat.Calculator instances
                to apply in place of the measurement's existing calculator_list. These
                calculators are given a chance, starting from the front of the list,
                to calculate a `DataSeries` for `item`.

        Returns: tuple of np.Array. The first array is time, and the second array is the
            corresponding (calculated) values of the requested item.
        """
        if (calculator_list is not None) or (remove_background is not None):
            self.clear_cache()

        if remove_background is False:
            if calculator_list is None:
                calculator_list = self._calculator_list
            calculator_list = [
                cal
                for cal in calculator_list
                if not type(cal) in self.background_calculator_types
            ]

        if calculator_list is not None:
            # this will now be the case if either (i) a calculator_list was given, or
            # (ii) remove_background was set to False.
            self._temp_calculator_list = calculator_list
            self.consolidate_calculators()

        vseries = self[item]
        self._temp_calculator_list = None
        self.consolidate_calculators()

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
            t_bg, v_bg = self.grab(
                item, tspan=tspan_bg, remove_background=remove_background
            )
            v = v - np.mean(v_bg)
        return t, v

    def grab_for_t(
        self, item, t, tspan_bg=None, remove_background=None, calculator_list=None
    ):
        """Return a numpy array with the value of item interpolated to time t

        Args:
            item (str): The name of the value to grab
            t (np array): The time vector to grab the value for
            tspan_bg (iterable): Optional. A timespan defining when `item` is at its
                baseline level. The average value of `item` in this interval will be
                subtracted from what is returned.
            remove_background (boolean): Whether to subtract a pre-set background, if
                available. This is True by default. If set to False, it suppresses
                background calculators, taken to be those specified by the class
                attribute `background_calculator_types`
        """
        t_0, v_0 = self.grab(
            item,
            tspan=[t[0], t[-1]],
            include_endpoints=True,
            tspan_bg=tspan_bg,
            remove_background=remove_background,
            calculator_list=calculator_list,
        )
        v = np.interp(t, t_0, v_0)
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
    def t(self):
        return self[self.control_series_name].t

    @property
    def t_name(self):
        return self[self.control_series_name].tseries.name

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
        if not self.time_names:  # No TimeSeries in the measurement means no tspan.
            return None
        for t_name in self.time_names:
            t = self[t_name].data
            if len(t) == 0:
                return None
            t_start = t[0] if t_start is None else min(t_start, t[0])
            t_finish = t[-1] if t_finish is None else max(t_finish, t[-1])
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
            if not m.tspan:
                # if it has no TimeSeries it must be a "constant". Best to include:
                new_component_measurements.append(m)
                continue
            # Otherwise we have to cut it according to the present tspan.
            dt = m.tstamp - self.tstamp
            try:
                tspan_m = [tspan[0] - dt, tspan[1] - dt]
            except IndexError:  # Apparently this can happen for empty files. See:
                continue  # https://github.com/ixdat/ixdat/issues/93
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
                    f"{self!r} does not have a default selection string "
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

    def select_values(self, *args, selector_name=None, **kwargs):
        """Return a selection of the measurement based on one or several criteria

        Specifically, this method returns a new Measurement where the time(s) returned
        are those where the values match the provided criteria, i.e. the part of the
        measurement where `self[series_name] == value`

        Any series can be selected for using the series name as a key-word. Arguments
        can be single acceptable values or lists of acceptable values.
        You can select for one or more series without valid python variable names by
        providing the kwargs using ** notation (see last example below).

        Arguments without key-word are considered valid values of the default
        selector, which is normally `self.selector_name` but can also be specified
        here using the key-word argument `selector_name`. Multiple criteria are
        applied sequentially, i.e. you get the intersection of satisfying parts.

        Examples of valid calls given a measurement `meas`:
        ```
        # to select where the default selector is 3, use:
        meas.select_values(3)
        # to select for where the default selector is 4 or 5:
        meas.select_values(4, 5)
        # to select for where "cycle" (i.e. the value of meas["cycle"].data) is 4:
        meas.select_values(cycle=4)
        # to select for where "loop_number" is 1 AND "cycle" is 3, 4, or 5:
        meas.select_values(loop_number=1, cycle=[3, 4, 5])
        # to select for where "cycle number" (notice the space) is 2 or 3:
        meas.select_values([2, 3], selector_name="cycle number")
        # which is equivalent to:
        meas.select_values(**{"cycle number": [2, 3]})

        Args:
            args (tuple): Argument(s) given without keyword are understood as acceptable
                value(s) for the selector (that named by selector_name or
                self.selector_name).
            selector_name: The name of the selector to which the args specify
            kwargs (dict): Each key-word arguments is understood as the name
                of a series and its acceptable value(s).
        """
        if args:
            # Then we must interpret the arguments as allowed values of a selector,
            #  either specified in the kwargs or the Measurement's default selector:
            selector_name = selector_name or self.selector_name
            if not selector_name:
                raise BuildError(
                    f"{self:r} does not have a default selector_name "
                    f"(Measurement.selector_name), and so selection only works "
                    f"with a selector_name specified "
                    f"(see `help(Measurement.select_values)`)"
                )
            # Get the args into a simple list:
            flat_args = []
            for arg in args:
                if hasattr(arg, "__iter__"):
                    flat_args += list(arg)
                else:
                    flat_args.append(arg)
            if selector_name in kwargs:
                raise ValueError(
                    "Don't call select_values with both arguments and "
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

        These all work for measurements that have a default selector and/or the
        indicated columns:
        - `meas.select(1, 2)`
        - `meas.select(tspan=[200, 300])`
        - `meas.select(range(10))`
        - `meas.select(cycle=4)`
        - `meas.select(**{"cycle number": [20, 21]})
        - `meas.select(loop_number=1, tspan=[1000, 2000])
        - `meas.select(1, range(5, 20), file_number=1, tspan=[1000, 2000])`
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
        from .spectra import SpectrumSeries, add_spectrum_series_to_measurement

        if isinstance(other, SpectrumSeries):
            return add_spectrum_series_to_measurement(self, other)

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

        new_series_list = list(set(self.series_list + other.series_list))
        new_component_measurements = list(
            set(
                (self.component_measurements or [self])
                + (other.component_measurements or [other])
            )
        )
        new_calculator_list = list(set(self._calculator_list + other._calculator_list))
        new_aliases = self.aliases.copy()
        for key, names in other.aliases.items():
            if key in new_aliases:
                new_aliases[key] = list(set(new_aliases[key] + other.aliases[key]))
            else:
                new_aliases[key] = other.aliases[key]
        obj_as_dict = self.as_dict()
        other_as_dict = other.as_dict()
        for k, v in other_as_dict.items():
            # Looking forward to the "|" operator!
            if k not in obj_as_dict:
                obj_as_dict[k] = v
        obj_as_dict.update(
            name=new_name,
            technique=new_technique,
            series_list=new_series_list,
            component_measurements=new_component_measurements,
            calculator_list=new_calculator_list,
            aliases=new_aliases,
        )
        # don't want the original calculators, component measurements, or series:
        del obj_as_dict["c_ids"]
        del obj_as_dict["m_ids"]
        del obj_as_dict["s_ids"]
        return cls.from_dict(obj_as_dict)

    def join(self, other, join_on=None):
        """Join two measurements based on a shared data series

        This involves projecting all timeseries from other's data series so that the
        variable named by `join_on` is shared between all data series.
        This is analogous to an explicit inner join.

        Args:
            other (Measurement): a second measurement to join to self
            join_on (str or tuple): Either a string, if the value to join on is called
                the same thing in both measurements, or a tuple of two strings where
                the first is the name of the variable in self and the second in other.
                The variable described by join_on must be monotonically increasing in
                both measurements.
        """
        raise NotImplementedError


class Calculator(Saveable):
    """Base class for calculators."""

    table_name = "calculator"
    calculator_type = None  # to be overwritten
    column_attrs = {"name", "technique", "tstamp", "calculator_type"}

    def __init__(self, *, name=None, technique=None, tstamp=None, measurement=None):
        """Initiate a Calculator

        Args:
            name (str): The name of the calculator
            technique (str): The technique of the calculator
            tstamp (float): The time at which the calculator took place or is valid
            measurement (Measurement): Optional. The measurement to use by default for
                raw data form which to calculate new series.
        """
        super().__init__()
        # NOTE: The :r syntax in f-strings doesn't work on None
        self.name = name or f"{self.__class__.__name__}()"
        self.technique = technique
        self.tstamp = tstamp or (measurement.tstamp if measurement else None)
        self.measurement = measurement

    def __str__(self):
        rep = (
            f"{self.__class__.__name__} providing: " + f"{self.available_series_names} "
        )
        return rep

    @classmethod
    def from_dict(cls, obj_as_dict):
        """Return an object of the Calculator class of the right technique

        Args:
              obj_as_dict (dict): The full serializaiton (rows from table and aux
                tables) of the measurement. obj_as_dict["technique"] specifies the
                technique class to use, from TECHNIQUE_CLASSES
        """
        # TODO: see if there isn't a way to put the import at the top of the module.
        #    see: https://github.com/ixdat/ixdat/pull/1#discussion_r546437410
        from .calculators import CALCULATOR_CLASSES

        calculator_type = obj_as_dict.pop("calculator_type")
        if calculator_type in CALCULATOR_CLASSES:
            calculator_class = CALCULATOR_CLASSES[calculator_type]
        else:
            calculator_class = cls
        try:
            calculator = calculator_class(**obj_as_dict)
        except Exception:
            raise
        return calculator

    def export(self, path_to_file=None):
        """Export a Calculator as a json-formatted text file"""
        path_to_file = path_to_file or (self.name + ".ix")
        self_as_dict = self.as_dict()
        with open(path_to_file, "w") as f:
            json.dump(self_as_dict, f, indent=4)

    @classmethod
    def read(cls, path_to_file):
        """Read a Calculator from a json-formatted text file"""
        with open(path_to_file) as f:
            obj_as_dict = json.load(f)
        return cls.from_dict(obj_as_dict)

    @property
    def available_series_names(self):
        """The set of the names of the series the Calculater can provide

        FIXME: Add more documentation about how to write this in inheriting classes.
        """
        raise NotImplementedError(
            f"class {type(self)} does not implement available_series_names"
        )

    def calculate_series(self, key, measurement=None):
        """This should be overwritten in real Calculator classes.

        FIXME: Add more documentation about how to write this in inheriting classes.
        """
        raise NotImplementedError

    def __add__(self, other):
        warnings.warn(
            f"Addition is not implemented for {type(self)}. Adding a second to a"
            " measurement causes the first to be ignored."
        )
        return other


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

    # if all else fails, we just join them with " and ". e.g. MS + XRD = MS and XRD
    return technique_1 + " and " + technique_2
