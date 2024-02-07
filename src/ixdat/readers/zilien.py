"""Readers for files produces by the Zilien software from Spectro Inlets.

Zilien tsv files have two data header lines to define each of the data columns.
The first one is referred to as "series header" and explains what the data describes,
and the second one is called "column header" and specifies the specific column.
It is done in order to keep the columns headers more readable.
Typically, a series header will specify the measuring device (e.g. "iongauge value")
or MS channel (e.g. "C0M2") and will apply for two or more column headers
where the first is time ("Time [s]", "time/s") and the subsequent are the corresponding
value(s) ("Pressure [mbar]" or "M2-H2 [A]" etc.).
Zilien files version 2 and higher may also include all the data from an integrated
Biologic dataset. These are grouped under the series header "EC-lab".
"""

import re
import sys
import time
from collections import defaultdict
from itertools import groupby, zip_longest
from pathlib import Path

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from ..data_series import DataSeries, TimeSeries, ValueSeries, Field
from ..techniques import (
    ECMSMeasurement,
    MSMeasurement,
    ECMeasurement,
    Measurement,
    TECHNIQUE_CLASSES,
)
from ..techniques.ms import MSSpectrum, MSSpectrumSeries, MSSpectroMeasurement
from .reading_tools import timestamp_string_to_tstamp
from ..exceptions import ReadError, TechniqueError


ZILIEN_TIMESTAMP_FORM = "%Y-%m-%d %H_%M_%S"  # like 2021-03-15 18_50_10
ZILIEN_MASS_COLUMN_NAMES = ["Mass  [AMU]", "Mass [AMU]"]
ZILIEN_EC_ALIASES = {
    "t": ["Potential time [s]"],
    "raw_potential": ["Voltage [V]"],
    "raw_current": ["Current [mA]"],
    "cycle": ["Cycle [n]"],
}
# The Zilien .tsv files can be loaded as three different experiment types. These are the
# aliases for each of them
ZILIEN_ALIASES = {
    ECMSMeasurement: ZILIEN_EC_ALIASES,
    MSMeasurement: {},
    ECMeasurement: ZILIEN_EC_ALIASES,
    MSSpectroMeasurement: {},
}

BIOLOGIC_SERIES_NAME = "EC-lab"

# TODO: When, in the future, Zilien files include the whole EC dataset, remove the
#    unflattering example presently in the docs.
#    https://github.com/ixdat/ixdat/pull/30/files#r810087496


def parse_metadata_line(line):
    """Parse a single metadata line and return the name, value"""
    # The metadata format is a 5 column format:
    name, comment, attach_to_series, type_as_str, value = line.strip("\n").split("\t")

    # Since, as yet, ixdat doesn't support per-series metadata, we prefix the per-series
    # metadata item names with the name of the series, to avoid name clashes while still
    # preserving the data
    if attach_to_series:
        full_name = f"{attach_to_series}_{name}"
    else:
        full_name = name

    # Type convert the metadata (the specification for version 1 also has a color type,
    # but it is not used yet)
    if type_as_str == "string":
        return full_name, value
    elif type_as_str == "int":
        # Python does not seem able to directly evaluate a string like "0.0" as an
        # integer. Therefore, we evaluate it as a float first and convert to int:
        return full_name, int(float(value))
    elif type_as_str == "double":
        return full_name, float(value)
    elif type_as_str == "bool":
        return full_name, value == "true"
    else:
        raise TypeError(f"Unknown metadata type {type_as_str} for {name}")


def to_snake_case(string):
    """Turn a space separated string into a snake_case string"""
    return string.lower().replace(" ", "_")


# Matches: "{name} [{unit}]"
ZILIEN_COLUMN_HEADER_RE = re.compile(r"^(.+?) \[(.+?)\]$")
# Matches: "{name}/{unit}" and "{name1/name2}/{unit}"
BIOLOGIC_COLUMN_HEADER_RE = re.compile(r"^(.+)/(.+)$")
# Matches: "C??M{mass}"
MASS_SERIES_RE = re.compile(r"^C[0-9]+M([0-9]+)$")


def to_mass(string):
    """Return mass (i.e. "18") if `string` matches the C0M18 mass series form or None"""
    possible_match = MASS_SERIES_RE.match(string)
    if possible_match:
        return possible_match.group(1)
    return None


def determine_class(technique):
    """Choose appropriate measurement class according to a given technique."""

    if technique in ("EC-MS", "EC", "MS", "MS-MS_spectra"):
        if technique == "MS-MS_spectra":
            # We read it as a MSMeasurement and then add the SpectrumSeries
            return MSMeasurement
        return TECHNIQUE_CLASSES[technique]
    else:
        raise TechniqueError(
            f'Unknown technique given: "{technique}". '
            "Use one of the following (in upper-case): "
            '"EC-MS", "EC, "MS", "MS-MS_spectra".'
        )


class ZilienTSVReader:
    """Class for reading files saved by Spectro Inlets' Zilien software"""

    def __init__(self):
        self._path_to_file = None
        self._cls = None
        self._measurement = None

        # start time of the Zilien measurement
        self._timestamp = None
        # a dictionary with metadata general information about the Zilien measurement
        self._metadata = None
        # a list with the Zilien TSV series headers,
        # such as "Iongauge value" (see module docstring)
        self._series_headers = None
        # a list with the Zilien TSV columns headers,
        # such as "Time [s]" and "Pressure [mbar]"
        self._column_headers = None
        # a numpy array of the parsed Zilien data (a big rectangle with NaN filling)
        self._data = None

    def read(
        self,
        path_to_file,
        cls=None,
        name=None,
        include_mass_scans=None,
        **kwargs,
    ):
        """Read a Zilien file

        Args:
            path_to_file (Path or str): The path of the file to read
            cls (Measurement): The measurement class to read the file as. Zilien tsv
                files can be read both as an EC-MS measurement, an MS measurement (which
                will exclude the EC series from the measurement) and as an EC measurement
                (which will exclude the MS series from the measurement). To avoid
                importing classes, this behavior can also be controlled by setting the
                `technique` argument to either 'EC-MS', 'MS' or 'EC'. The default is
                determined according to what is parsed from the dataset.
            name (str): The name of the measurement. Will default to the part of the
                filename before the '.tsv' extension
            include_mass_scans (bool): Whether to include mass scans (if available) and
                thereby return a `SpectroMSMeasurement` which can be indexed to give
                the spectrum objects. (Defaults to True if technique not specified.)
            kwargs: All remaining keyword-arguments will be passed onto the `__init__`
                of the Measurement
        """
        if self._path_to_file:
            print(
                f"This {self.__class__.__name__} has already read {self._path_to_file}. "
                "Returning the measurement resulting from the original read. "
                "Use a new Reader if you want to read another file."
            )
            return self._measurement

        self._path_to_file = Path(path_to_file)

        # Parse metadata items
        with open(self._path_to_file, encoding="utf-8") as file_handle:
            (
                self._metadata,
                self._series_headers,
                self._column_headers,
            ) = self._read_metadata(file_handle)
            file_position = file_handle.tell()

        # Read raw data
        with open(self._path_to_file, "rb") as file_handle:
            file_handle.seek(file_position)
            self._data = np.genfromtxt(file_handle, delimiter="\t")

        if "technique" in kwargs:
            technique = kwargs["technique"]
        else:
            # We don't put "technique directly into kwargs because that would prevent
            # reading of mass scans by default.
            technique = self._get_technique()

        if cls is Measurement or cls is None:
            self._cls = determine_class(technique)
        else:
            self._cls = cls

        if include_mass_scans is None:
            # This becomes True if neither a class nor a technique was specified,
            # or if a technique or class including spectra was specified.
            include_mass_scans = (
                not cls or cls is Measurement or issubclass(cls, MSSpectroMeasurement)
            ) and ("MS-MS_spectra" in kwargs.get("technique", "MS-MS_spectra"))

        if "technique" not in kwargs:
            # So that technique gets passed on to the Measurement's __init__:
            kwargs["technique"] = technique

        # Part of filename before the extension
        file_stem = self._path_to_file.stem

        if "start_time_unix" in self._metadata:
            # Extract unix timestamp from metadata
            self._timestamp = float(self._metadata["start_time_unix"])
        else:
            # Extract timestamp from filename on form:
            # 2021-04-20 11_16_18 Measurement name
            self._timestamp = timestamp_string_to_tstamp(
                timestamp_string=" ".join(file_stem.split(" ")[:2]),
                form=ZILIEN_TIMESTAMP_FORM,
            )

        # Extract series data and form series
        series, aliases = self._form_series()
        for standard_name, general_aliases in ZILIEN_ALIASES[self._cls].items():
            aliases[standard_name] += general_aliases
        aliases = dict(aliases)  # Convert from defaultdict to normal dict

        measurement_kwargs = {
            "name": name or file_stem,
            "series_list": series,
            "aliases": aliases,
            "tstamp": self._timestamp,
            "metadata": self._metadata,
        }
        measurement_kwargs.update(kwargs)
        self._measurement = self._cls(**measurement_kwargs)

        if include_mass_scans:
            # Check if there are MS spectra.
            # If the file name is "YYYY-MM-DD hh_mm_ss my_file_name.tsv", then
            # the spectra are the .tsv files in the folder "my_file_name mass scans"
            spectra_folder = self._path_to_file.parent / (
                self._path_to_file.stem[20:] + " mass scans"
            )
            if spectra_folder.exists():  # Then we have a spectra folder!
                spectrum_series = MSSpectrumSeries.read_set(
                    spectra_folder,
                    suffix=".tsv",
                    reader="zilien",
                    t_zero=self._timestamp,
                )
                self._measurement = self._measurement + spectrum_series

        return self._measurement

    @staticmethod
    def _read_metadata(file_handle):
        """Read metadata from `file_handle`.

        Description of the returned variables:
        - metadata: A dictionary of meta dataset information
        - series_headers: A list with string series headers (the first header line).
          A series header is a group of column headers. They have possible one
          or more empty cells (empty strings) inbetween them, so extra columns
          in the group can be aligned. E.g. "Iongauge value", "MFC1 setpoint",
          "MFC1 value", etc.
        - column_headers: A list with string column headers (the second header line).
          A column header is a name of one data column in a series. E.g. "Time [s]",
          "Pressure [mbar]", "Flow [ml/min]", etc.

        The length of the "series_headers" and the "column_headers" is the same.

        Returns:
            Tuple[Dict[str, str | int | float | bool], List[String], List[String]]:
            Three variables with metadata, series and column headers described above.

        """
        # The first 4 lines always include the file version, number of header lines,
        # number of data header lines and data start line in this order.
        # Backwards compatibility is ensured, because the one extra line read will be
        # just added to the 'metadata' dict and ignored later.
        metadata = {}
        fixed_metadata_lines_amount = 4
        for _ in range(fixed_metadata_lines_amount):
            key, value = parse_metadata_line(file_handle.readline())
            metadata[key] = value

        # read the rest when the total amount is known
        for _ in range(metadata["num_header_lines"] - fixed_metadata_lines_amount):
            key, value = parse_metadata_line(file_handle.readline())
            metadata[key] = value

        # version 1 of the file format is sometimes missing this value
        if "file_format_version" not in metadata:
            metadata["file_format_version"] = 1

        series_headers = file_handle.readline().strip("\n").split("\t")
        column_headers = file_handle.readline().strip("\n").split("\t")

        return metadata, series_headers, column_headers

    def _get_technique(self):
        """Get technique according to parsed series headers and column headers.

        This method covers the following cases:
          - "pot" and "EC-lab" is in the metadata and in the headers.
            There is "EC-lab" data.
          - "pot" and "EC-lab" is in the metadata and in the headers.
            There is no "EC-lab" data. (bug)
          - "pot" is not in the metadata, but it is in the headers and the data
            part is filled with NaN. (bug)
          - "pot" is neither in the metadata, nor in the headers.

        """

        # Should stay "pot_pot_count", because the series name is prepended in
        # order to keep unique series names. See `parse_metadata_line()`.
        pot_in_metadata = "pot_pot_count" in self._metadata
        pot_in_series_headers = "pot" in self._series_headers

        result = ""
        # there was EC measurement
        if pot_in_metadata and pot_in_series_headers:
            result = "EC-MS"

            # BUG CASE: EC-lab series header gets included, but .mpt data don't
            missing_column_headers = len(self._column_headers) < len(
                self._series_headers
            )
            empty_eclab_series = self._series_headers[-1] == "EC-lab"
            if missing_column_headers and empty_eclab_series:
                result = "MS"
                sys.stderr.write(
                    "EC-lab data is missing in the dataset. That means Zilien did not"
                    "find .mpt files when creating it. You can convert .mpr files into "
                    ".mpt files in EC-lab (or read the .mpr files directly) and then "
                    "connect the two Measurement objects.\n"
                    "Reading only the Mass Scan part of the dataset.\n\n"
                )
        # BUG CASE: Zilien internal bug (https://github.com/ixdat/ixdat/issues/114)
        elif not pot_in_metadata and pot_in_series_headers:
            result = "MS"
        # no EC did run and EC is not included in the settings dialog
        elif not pot_in_metadata and not pot_in_series_headers:
            result = "MS"
        else:
            # unexpected case, so go with the safest option
            result = "MS"

        return result

    def _form_series(self):
        """Form the series and series aliases

        Returns:
            List[Series], DefaultDict(str, List[str]): List of series and dict of aliases
        """
        aliases = defaultdict(list)
        series = []

        # Get non-empty series headers and their indices
        # in order to process whole series chunks
        series_split_indices, nonempty_headers = self._get_series_splits(
            self._series_headers
        )

        for series_header, (begin, end) in zip(nonempty_headers, series_split_indices):
            # Skip series not relevant for the type of measurement
            if not issubclass(self._cls, ECMeasurement) and series_header in (
                "pot",
                BIOLOGIC_SERIES_NAME,
            ):
                continue
            elif not issubclass(self._cls, MSMeasurement) and to_mass(series_header):
                continue

            column_headers_split = self._column_headers[begin:end]
            data_columns_split = self._data[:, begin:end]

            if series_header == BIOLOGIC_SERIES_NAME:
                column_series = self._biologic_dataset_part(
                    column_headers_split, data_columns_split
                )
            else:
                column_series, aliases_part = self._zilien_dataset_part(
                    series_header, column_headers_split, data_columns_split
                )
                # update aliases
                for standard_name, series_name in aliases_part.items():
                    aliases[standard_name] += series_name

            series += column_series

        return series, aliases

    def _zilien_dataset_part(self, series_header, column_headers, data_columns_split):
        """Process necessary data for a Zilien dataset part.

        Args:
            series_header (str): The current series name.
            column_headers (list): A list with column names from the current series.
            data_columns_split (np.array): A columns taken out from the parsed dataset,
                that represents the current series.

        Returns:
            list, DefaultDict[List]: A list of Ixdat series objects and
            a default dict with standard names for a Mass series.
        """

        count = self._metadata[f"{series_header}_{series_header}_count"]
        names_and_units = [
            self._form_names_and_unit(series_header, column_header)
            for column_header in column_headers
        ]

        # Fill Mass aliases
        aliases = defaultdict(list)
        for series_name, _, standard_name in names_and_units:
            if standard_name:
                aliases[standard_name].append(series_name)

        # Create Ixdat series
        column_series = self._create_series_objects(
            column_headers, names_and_units, data_columns_split[:count, :]
        )

        return column_series, aliases

    def _biologic_dataset_part(self, column_headers, data_columns_split):
        """Process necessary data for a Biologic dataset part.

        The `experiment_number` and the `technique_number` columns are used only
        to create an information how to split the given data part into rows.
        After that they are not used anymore and no Ixdat series objects
        are created from them.

        Args:
            column_headers (list): A list with column names from the current series.
            data_columns_split (np.array): A columns taken out from the parsed dataset,
                that represents the current series.

        Returns:
            list: A list of Ixdat series objects.
        """

        count = self._metadata[f"{BIOLOGIC_SERIES_NAME}_{BIOLOGIC_SERIES_NAME}_count"]
        names_and_units = [
            self._form_names_and_unit(BIOLOGIC_SERIES_NAME, column_header)
            for column_header in column_headers
        ]

        # Split rows according to techniques used according to experiments
        # Combine the experiment numbers and the technique numbers
        # in order to create unique identifiers for a successful split
        exp_nums = data_columns_split[:count, column_headers.index("experiment_number")]
        tech_nums = data_columns_split[:count, column_headers.index("technique_number")]
        split_vector = exp_nums * 1000 + tech_nums
        splits = self._get_biologic_splits(split_vector)

        # Create Ixdat series
        column_series = []
        for begin, end in splits:
            column_series += self._create_series_objects(
                column_headers, names_and_units, data_columns_split[begin:end, :]
            )

        return column_series

    def _create_series_objects(self, column_headers, names_and_units, data_rows_split):
        """Create an Ixdat series objects from a given portion of a dataset.

        The `experiment_number` and the `technique_number` columns are skipped,
        because they are used only to split rows in the Biologic dataset by the
        Zilien reader, in order to create the same series as Ixdat would.

        Args:
            column_headers (list): A list with column names from the current series.
            names_and_units (list): A tuple with three elements. A series name,
                a unit and a standard name for every given column. (The series name
                is same for all columns here.)
            data_rows_split (np.array): A rows taken out from the columns part.
                In Zilien part it represents the Zilien measurement.
                In Biologic part it represents the current technique.

        Returns:
            list: A list of Ixdat series objects.
        """

        series_objects = []
        time_series = None

        for column_number, column_header in enumerate(column_headers):
            # Skip meta columns in the EC-lab dataset
            if column_header in ("experiment_number", "technique_number"):
                continue

            column_data = data_rows_split[:, column_number]

            # Skip holes in the EC-lab dataset
            if np.isnan(column_data).all():
                continue

            # Form series kwargs
            series_name, unit, standard_name = names_and_units[column_number]
            series_kwargs = {
                "name": series_name,
                "unit_name": unit,
                "data": column_data,
            }

            # Create the series
            if column_header in ("Time [s]", "time/s"):
                series_object = TimeSeries(**series_kwargs, tstamp=self._timestamp)
                time_series = series_object
            else:
                if time_series is None:
                    raise ValueError("Time column must be first in a dataset series.")

                series_object = ValueSeries(**series_kwargs, tseries=time_series)

            series_objects.append(series_object)

        return series_objects

    # --- UTILS ---

    @staticmethod
    def _get_series_splits(series_headers):
        """Create series names and their index pairs in the whole read dataset.

        Args:
            series_headers (list): Series name headers (even empty).

        Returns:
            list, list: List of tuples with index pairs
            and a list with non-empty series name headers.
        """
        series_split_indices = []
        nonempty_headers = []

        for index, series_header in enumerate(series_headers):
            if series_header != "":
                series_split_indices.append(index)
                nonempty_headers.append(series_header)

        # create split pairs
        # zip_longest to have the last pair with the index of the last
        # series header to the end
        series_split_indices = [
            (begin, end)
            for begin, end in zip_longest(series_split_indices, series_split_indices[1:])
        ]

        return series_split_indices, nonempty_headers

    @staticmethod
    def _get_biologic_splits(technique_numbers):
        """Create index pairs of row splits in the biologic dataset columns.

        Args:
            technique_numbers (np.array): Numbers of techniques during
                an EC-lab experiment(s).

        Returns:
            list: List of tuples with index pairs.
        """
        index = 0
        splits = []

        for _, group in groupby(technique_numbers):
            group_length = len(list(group))
            splits.append((index, index + group_length))
            index += group_length

        return splits

    @staticmethod
    def _form_names_and_unit(series_header, column_header):
        """Form names and unit from headers.

        Args:
            series_header (str): Something like "Iongauge value" or "C0M18"
            column_header (str): Something like "Time [s]" or "Flow [ml/min]"

        Returns:
            str, str, Optional[str]: Return series_name, unit, standard_name
        """
        standard_name = None
        if column_header in ("Time [s]", "time/s"):  # Form TimeSeries
            unit = "s"
            if series_header == "pot":
                name = f"Potential {column_header.lower()}"
            elif series_header == BIOLOGIC_SERIES_NAME:
                name = f"Biologic {column_header.lower()}"
            else:
                name = f"{series_header} {column_header.lower()}"
        else:  # ValueSeries
            # Perform a bit of reasonable name adaption, first break name and unit out
            # from the column header on the form: Pressure [mbar]
            zilien_components_match = ZILIEN_COLUMN_HEADER_RE.match(column_header)
            biologic_components_match = BIOLOGIC_COLUMN_HEADER_RE.match(column_header)

            if zilien_components_match:
                _, unit = zilien_components_match.groups()
            elif biologic_components_match:
                _, unit = biologic_components_match.groups()
            else:
                _, unit = column_header, ""

            # Is the column a "setpoint" or "value" type
            setpoint_or_value = None
            for option in ("setpoint", "value"):
                if series_header.endswith(option):
                    setpoint_or_value = option

            mass = to_mass(series_header)
            if setpoint_or_value:
                # In that case, the column header is something like "Flow [ml/min]" where
                # "Flow" is unnecessary, because that is apparent from the unit
                # The name will look for example like this "MFC setpoint [ml/min]"
                name = f"{series_header} [{unit}]"
            elif mass is not None:
                # e.g. from series header "C1M4" and column header "M4-He [A]"
                # the name will be "M4 [A]" and standard name will be "M4"
                name = f"M{mass} [{unit}]"
                standard_name = f"M{mass}"
            else:
                name = column_header

        return name, unit, standard_name


class ZilienTMPReader:
    """A class for stitching the files in a Zilien tmp directory to an ECMSMeasurement

    This is necessary because Zilien often crashes, leaving only the tmp directory.
    This is less advanced but more readable than the Spectro Inlets stitching solution.
    """

    def __init__(self, path_to_tmp_dir=None):
        self.path_to_tmp_dir = Path(path_to_tmp_dir) if path_to_tmp_dir else None

    def read(self, path_to_tmp_dir, cls=None, **kwargs):
        """Make a measurement from all the single-value .tsv files in a Zilien tmp dir

        Args:
            path_to_tmp_dir (Path or str): The path to the tmp dir
            cls (Measurement class): Defaults to ECMSMeasurement
        """
        if path_to_tmp_dir:
            self.path_to_tmp_dir = Path(path_to_tmp_dir)
        cls = cls or ECMSMeasurement
        name = self.path_to_tmp_dir.parent.name
        timestamp_string = name[:19]  # the zilien timestamp is the first 19 chars
        tstamp = timestamp_string_to_tstamp(timestamp_string, form=ZILIEN_TIMESTAMP_FORM)
        series_list = []
        for tmp_file in self.path_to_tmp_dir.iterdir():
            series_list += series_list_from_tmp(tmp_file)
        obj_as_dict = {
            "name": name,
            "tstamp": tstamp,
            "series_list": series_list,
            "technique": "EC-MS",
            "reader": self,
        }
        obj_as_dict.update(kwargs)
        return cls.from_dict(obj_as_dict)


def series_list_from_tmp(path_to_file):
    """Return [ValueSeries, TimeSeries] with the data in a zilien tmp .tsv file"""
    file_name = Path(path_to_file).name
    timestamp_string = file_name[:19]  # the zilien timestamp form is 19 chars long
    tstamp = timestamp_string_to_tstamp(timestamp_string, form=ZILIEN_TIMESTAMP_FORM)
    column_match = re.search(r"\.([^\.]+)\.data", file_name)
    if not column_match:
        print(f"could not find column name in {path_to_file}")
        return []
    v_name = column_match.group(1)
    mass_match = re.search("M[0-9]+", v_name)
    if mass_match:
        v_name = mass_match.group()
        unit = "A"
    else:
        unit = None
    t_name = v_name + "-x"
    df = pd.read_csv(path_to_file, delimiter="\t", names=[t_name, v_name], header=0)
    t_data, v_data = df[t_name].to_numpy(), df[v_name].to_numpy()
    tseries = TimeSeries(name=t_name, unit_name="s", data=t_data, tstamp=tstamp)
    vseries = ValueSeries(name=v_name, unit_name=unit, data=v_data, tseries=tseries)
    return [tseries, vseries]


class ZilienSpectrumReader:
    """A reader for individual Zilien spectra
    TODO: A Zilien reader which loads all spectra at once in a SpectrumSeries object
    """

    def __init__(self, path_to_spectrum=None):
        self.path_to_spectrum = Path(path_to_spectrum) if path_to_spectrum else None

    def read(self, path_to_spectrum, cls=None, t_zero=None, **kwargs):
        """Read a Zilien spectrum.
        FIXME: This reader was written hastily and could be designed better.

        Args:
            path_to_spectrum(Path or str): the path to the spectrum file
            cls (Spectrum class): Defaults to MSSpectrum
            t_zero (float): The unix timestamp which the mass scan start time
                is referenced to. Should be the tstamp of the corresponding
                Zilien measurement. If the Spectrum is read individually, it needs
                to be input or defaults to `time.time()`, i.e., now.
            kwargs: Key-word arguments are passed on ultimately to cls.__init__
        """
        if path_to_spectrum:
            self.path_to_spectrum = Path(path_to_spectrum)
        cls = cls or MSSpectrum
        t_zero = t_zero or time.time()

        with open(self.path_to_spectrum, "r") as f:
            in_header = True
            num_line = 1
            metadata_dict = {}
            while in_header:
                name, value = parse_metadata_line(f.readline())
                metadata_dict[name] = value
                if num_line == metadata_dict.get("num_header_lines", 0):
                    in_header = False
                num_line += 1

        tstamp = t_zero + metadata_dict["mass_scan_started_at"]
        duration = (
            (metadata_dict["stop_mass"] - metadata_dict["start_mass"])
            * metadata_dict["points_per_amu"]
            * metadata_dict["dwell_time"]
            * 1e-3
        )

        df = pd.read_csv(
            self.path_to_spectrum,
            header=metadata_dict["data_start"] - 1,
            delimiter="\t",
        )
        y_name = "Current [A]"

        for x_name in ZILIEN_MASS_COLUMN_NAMES:
            try:
                x = df[x_name].to_numpy()
            except KeyError:
                continue
            break
        else:
            raise ReadError(
                f"Can't find a mass column in {self.path_to_spectrum}. "
                f"Looked for one of {ZILIEN_MASS_COLUMN_NAMES}"
            )
        y = df[y_name].to_numpy()

        xseries = DataSeries(data=x, name=x_name, unit_name="m/z")
        field = Field(
            data=np.array(y),
            name=y_name,
            unit_name="A",
            axes_series=[
                xseries,
            ],
        )
        obj_as_dict = {
            "name": self.path_to_spectrum.name,
            "technique": "MS_spectrum",
            "field": field,
            "reader": self,
            "tstamp": tstamp,
            "duration": duration,
        }
        obj_as_dict.update(kwargs)
        return cls.from_dict(obj_as_dict)


def module_demo():
    """Module demo here.

    To run this module in PyCharm, open Run Configuration and set
        Module name = ixdat.readers.zilien,
    and *not*
        Script path = ...
    """
    path_to_test_file = (
        Path(__file__).parent.resolve().parent.parent.parent
        / "test_data"
        / "Zilien version 1"
        / "2022-04-06 16_17_23 full set.tsv"
    )

    ecms_measurement = Measurement.read(
        reader="zilien",
        path_to_file=path_to_test_file,
    )

    ecms_measurement.plot_measurement()
    plt.show()


if __name__ == "__main__":
    module_demo()
