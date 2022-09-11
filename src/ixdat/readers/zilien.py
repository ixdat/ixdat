"""Readers for files produces by the Zilien software from Spectro Inlets"""

import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from ..data_series import ConstantValue, DataSeries, TimeSeries, ValueSeries, Field
from ..techniques import ECMSMeasurement, MSMeasurement, ECMeasurement, Measurement
from ..techniques.ms import MSSpectrum
from .reading_tools import timestamp_string_to_tstamp, FLOAT_MATCH


ZILIEN_TIMESTAMP_FORM = "%Y-%m-%d %H_%M_%S"  # like 2021-03-15 18_50_10

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
}

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
    # but as of yet it is not used)
    if type_as_str == "string":
        return full_name, value
    elif type_as_str == "int":
        return full_name, int(value)
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


class ZilienTSVReader:
    """Class for reading files saved by Spectro Inlets' Zilien software"""

    def __init__(self):
        self._path_to_file = None
        self._cls = None
        self._measurement = None

    def read(self, path_to_file, cls=ECMSMeasurement, name=None, **kwargs):
        """Read a Zilien file

        Args:
            path_to_file (Path or str): The path of the file to read
            cls (Measurement): The measurement class to read the file as. Zilien tsv
                files can be read both as an ECMS measurement, a MS measurement (which
                will exclude the EC series from the meaurement) and as a ECMeasurement
                (which will exclude the MS series from the measurement). To avoid
                importing classes, this behavior can also be controlled by setting the
                `technique` argument to either 'EC-MS', 'MS' or 'EC'. The deafult is a
                ECMSMeasurement.
            name (str): The name of the measurement. Will default to the part of the
                filename before the '.tsv' extension
            kwargs: All remaining keywor-arguments will be passed onto the `__init__` of
                the Meaurement

        """
        if self._path_to_file:
            print(
                f"This {self.__class__.__name__} has already read {self._path_to_file}. "
                "Returning the measurement resulting from the original read. "
                "Use a new Reader if you want to read another file."
            )
            return self._measurement

        if "technique" in kwargs:
            if kwargs["technique"] == "EC-MS":
                cls = ECMSMeasurement
            if kwargs["technique"] == "EC":
                cls = ECMeasurement
            if kwargs["technique"] == "MS":
                cls = MSMeasurement
        else:
            if cls is Measurement:
                cls = ECMSMeasurement
            if issubclass(cls, ECMSMeasurement):
                kwargs["technique"] = "EC-MS"
            elif issubclass(cls, ECMeasurement):
                kwargs["technique"] = "EC"
            elif issubclass(cls, MSMeasurement):
                kwargs["technique"] = "MS"
        self._cls = cls

        self._path_to_file = Path(path_to_file)

        # Extract timestamp from filename on form: 2021-04-20 11_16_18 Measurement name
        file_stem = self._path_to_file.stem  # Part of filename before the extension
        timestamp = timestamp_string_to_tstamp(
            timestamp_string=" ".join(file_stem.split(" ")[:2]),
            form=ZILIEN_TIMESTAMP_FORM,
        )

        # Parse metadata items
        with open(self._path_to_file, encoding="utf-8") as file_handle:
            metadata, series_headers, column_headers = self._read_metadata(file_handle)
            file_position = file_handle.tell()

        # Read raw data
        with open(self._path_to_file, "rb") as file_handle:
            file_handle.seek(file_position)
            data = np.genfromtxt(file_handle, delimiter="\t")

        # Extract series data and form series
        series, aliases = self._form_series(
            data, metadata, timestamp, series_headers, column_headers
        )
        for standard_name, general_aliases in ZILIEN_ALIASES[self._cls].items():
            aliases[standard_name] += general_aliases
        aliases = dict(aliases)  # Convert from defaultdict to normal dict

        measurement_kwargs = {
            "name": name or file_stem,
            "series_list": series,
            "aliases": aliases,
            "tstamp": timestamp,
            "metadata": metadata,
        }
        measurement_kwargs.update(kwargs)
        self._measurement = cls(**measurement_kwargs)
        return self._measurement

    def _form_series(self, data, metadata, timestamp, series_headers, column_headers):
        """Form the series and series aliases

        Args:
            data (numpy.array): The data block of the tsv file as an array
            metadata (dict): Extracted metadata
            timestamp (float): The timestamp of the measurement
            series_headers (List[str]): List of series headers, slots with empty strings
                means the same as the last non-empty one
            column_headers (List[str]): List of column headers

        Returns:
            List[Series], DefaultDict(str, List[str]): List of series and dict of aliases
        """
        series = []
        last_time_series = None
        aliases = defaultdict(list)
        last_series_header = ""

        # Iterate over numbered series and column headers
        for column_number, (series_header, column_header) in enumerate(
            zip(series_headers, column_headers)
        ):
            last_series_header = series_header or last_series_header

            # Skip series not relevant for the type of measurement
            if not issubclass(self._cls, ECMeasurement) and last_series_header in (
                "pot",
                "EC-lab",
            ):
                continue
            elif not issubclass(self._cls, MSMeasurement) and to_mass(
                last_series_header
            ):
                continue

            # Pluck column of the correct length out from the data block and form series
            count = metadata[f"{last_series_header}_{last_series_header}_count"]
            column_data = data[:count, column_number]

            # Form the series_name, unit, aliases and update aliases
            series_name, unit, standard_name = self._form_names_and_unit(
                last_series_header, column_header
            )
            if standard_name:
                aliases[standard_name].append(series_name)

            # Form series kwargs and the series
            series_kwargs = {
                "name": series_name,
                "unit_name": unit,
                "data": column_data,
            }
            if column_header in ("Time [s]", "time/s"):  # Form TimeSeries
                column_series = TimeSeries(**series_kwargs, tstamp=timestamp)
                last_time_series = column_series
            elif np.all(column_data == column_data[0]):
                column_series = ConstantValue(**series_kwargs, tseries=last_time_series)
            # TODO-O remove `column_header` for this one (?)
            elif all(np.isnan(column_data)):
                continue
            else:
                column_series = ValueSeries(**series_kwargs, tseries=last_time_series)

            series.append(column_series)

        return series, aliases

    @staticmethod
    def _form_names_and_unit(series_header, column_header):
        """Return names and unit

        Args:
            series_header (str): Something like "Iongauge value" or "C0M18"
            column_header (str): Something like "Time [s]" or "Time [s]"

        Returns:
            str, str, Optional[str]: Return series_name, unit, standard_name

        """
        standard_name = None
        if column_header in ("Time [s]", "time/s"):  # Form TimeSeries
            unit = "s"
            if series_header == "pot":
                name = f"Potential {column_header.lower()}"
            elif series_header == "EC-lab":
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

            if setpoint_or_value:
                # In that case, the column header is something like "Flow [ml/min]" where
                # "Flow" is unnecessary, because that is apparent from the unit
                name = f"{series_header} [{unit}]"
            elif to_mass(series_header) is not None:
                mass = to_mass(series_header)
                name = f"M{mass} [{unit}]"
                standard_name = f"M{mass}"
            else:
                name = column_header

        return name, unit, standard_name

    @staticmethod
    def _read_metadata(file_handle):
        """Read metadata from `file_handle`"""

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

    def read(self, path_to_spectrum, cls=None, **kwargs):
        """Reat a Zilien spectrum.
        FIXME: This reader was written hastily and could be designed better.

        Args:
            path_to_tmp_dir (Path or str): the path to the tmp dir
            cls (Spectrum class): Defaults to MSSpectrum
            kwargs: Key-word arguments are passed on ultimately to cls.__init__
        """
        if path_to_spectrum:
            self.path_to_spectrum = Path(path_to_spectrum)
        cls = cls or MSSpectrum
        df = pd.read_csv(
            path_to_spectrum,
            header=9,
            delimiter="\t",
        )
        x_name = "Mass  [AMU]"
        y_name = "Current [A]"
        x = df[x_name].to_numpy()
        y = df[y_name].to_numpy()
        with open(self.path_to_spectrum, "r") as f:
            for i in range(10):
                line = f.readline()
                if "Mass scan started at [s]" in line:
                    tstamp_match = re.search(FLOAT_MATCH, line)
                    tstamp = float(tstamp_match.group())
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
            "name": path_to_spectrum.name,
            "technique": "MS",
            "field": field,
            "reader": self,
            "tstamp": tstamp,
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
