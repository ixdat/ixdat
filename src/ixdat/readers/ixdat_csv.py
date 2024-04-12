"""Module defining the ixdat csv reader, so ixdat can read the files it exports."""

from pathlib import Path
import json
import numpy as np
import re
import pandas as pd
from ..exceptions import ReadError
from ..data_series import ValueSeries, TimeSeries, DataSeries, Field
from ..measurements import Measurement
from ..spectra import Spectrum, SpectrumSeries
from ..techniques import TECHNIQUE_CLASSES

regular_expressions = {
    "name": r"^name = (.+)\n",
    "tstamp": r"tstamp = ([0-9\.]+)\n",
    "technique": r"technique = (.+)\n",
    "N_header_lines": r"N_header_lines = ([0-9]+)\n",
    "backend_name": r"backend_name = (.+)\n",
    "id": r"id = ([0-9]+)",
    "timecol": r"timecol '(.+)' for: (?:'(.+)')$",
    "unit": r"/ \[(.+)\]",
    "aux_file": r"'(.*)' in file: '(.*)'",
}
bad_keys = ("time_step",)


class IxdatCSVReader:
    """A class that reads the csv's made by ixdat.exporters.csv_exporter.CSVExporter

    read() is the important method - it takes the path to the mpt file as argument
    and returns an ECMeasurement object (ec_measurement) representing that file.
    The ECMeasurement contains a reference to the BiologicMPTReader object, as
    ec_measurement.reader. This makes available all the following stuff, likely
    useful for debugging.

    Attributes:
        path_to_file (Path): the location and name of the file read by the reader
        n_line (int): the number of the last line read by the reader
        place_in_file (str): The last location in the file read by the reader. This
            is used internally to tell the reader how to parse each line. Options are:
            "header", "column names", and "data".
        header_lines (list of str): a list of the header lines of the files. This
            includes the column name line. The header can be nicely viewed with the
            print_header() function.
        tstamp (str): The unix time corresponding to t=0
        technique (str): The name of the technique
        N_header_lines (int): The number of lines in the header of the file
        column_names (list of str): The names of the data columns in the file
        column_data (dict of str: np.array): The data in the file as a dict.
            Note that the np arrays are the same ones as in the measurement's DataSeries,
            so this does not waste memory.
        file_has_been_read (bool): This is used to make sure read() is only successfully
            called once by the Reader. False until read() is called, then True.
        measurement (Measurement): The measurement returned by read() when the file is
            read. self.measureemnt is None before read() is called.
    """

    delim = ","

    def __init__(self):
        """Initialize a Reader for ixdat-exported .csv files. See class docstring."""
        self.name = None
        self.path_to_file = None
        self.n_line = 0  # TODO: decide if this is part of API.
        # as per https://github.com/ixdat/ixdat/pull/30/files#r816204939
        self.place_in_file = "header"
        self.header_lines = []
        self.tstamp = None
        self.N_header_lines = None
        self.timecols = {}
        self.column_names = []
        self.column_data = {}
        self.technique = None
        self.aux_file_objects = {}
        self.measurement_class = Measurement
        self.file_has_been_read = False
        self.measurement = None
        self.meas_as_dict = {}

    def read(self, path_to_file, name=None, cls=None, **kwargs):
        """Return a Measurement with the data and metadata recorded in path_to_file

        This loops through the lines of the file, processing one at a time. For header
        lines, this involves searching for metadata. For the column name line, this
        involves creating empty arrays for each data series. For the data lines, this
        involves appending to these arrays. After going through all the lines, it
        converts the arrays to DataSeries.
        The technique is specified in the header, and used to pick the
        TechniqueMeasurement class.
        Finally, the method returns a TechniqueMeasurement object `measurement`
        with these DataSeries. All attributes of this reader can be accessed from the
        measurement as `measurement.reader.attribute_name`.

        Args:
            path_to_file (Path): The full abs or rel path including the ".mpt" extension
            name (str): The name of the measurement to return (defaults to path_to_file)
            cls (Measurement subclass): The class of measurement to return. By default,
                cls will be determined from the technique specified in the header of
                path_to_file.
            **kwargs (dict): Key-word arguments are passed to ECMeasurement.__init__

        Returns cls: a Measurement of type cls
        """
        path_to_file = Path(path_to_file) if path_to_file else self.path_to_file
        if self.file_has_been_read:
            print(
                f"This {self.__class__.__name__} has already read {self.path_to_file}."
                " Returning the measurement resulting from the original read. "
                "Use a new Reader if you want to read another file."
            )
            return self.measurement
        self.name = name or path_to_file.name
        self.path_to_file = path_to_file

        with open(self.path_to_file, "r") as f:
            for line in f:
                self.process_line(line)

        for name in self.column_names:
            self.column_data[name] = np.array(self.column_data[name])

        data_series_dict = {}

        for tcol_name in self.timecols:  # then it's time!
            data_series_dict[tcol_name] = TimeSeries(
                name=tcol_name,
                unit_name=get_column_unit(tcol_name) or "s",
                data=self.column_data[tcol_name],
                tstamp=self.tstamp,
            )

        for column_name, data in self.column_data.items():
            if column_name in self.timecols:
                continue
            try:
                tcol_name = next(
                    tcol_name
                    for tcol_name in self.timecols
                    if column_name in self.timecols[tcol_name]
                )
            except StopIteration:  # debugging
                raise ReadError(
                    f"can't find tcol for {column_name}. timecols={self.timecols}"
                )

            tseries = data_series_dict[tcol_name]
            vseries = ValueSeries(
                name=column_name,
                data=data,
                tseries=tseries,
                unit_name=get_column_unit(column_name),
            )
            data_series_dict[column_name] = vseries

        data_series_list = list(data_series_dict.values())
        self.meas_as_dict.update(
            name=self.name,
            technique=self.technique,
            reader=self,
            series_list=data_series_list,
            tstamp=self.tstamp,
        )
        self.meas_as_dict.update(self.aux_file_objects)
        self.meas_as_dict.update(kwargs)

        if issubclass(cls, self.measurement_class):
            self.measurement_class = cls

        self.measurement = self.measurement_class.from_dict(self.meas_as_dict)
        self.file_has_been_read = True
        return self.measurement

    def process_line(self, line):
        """Call the correct line processing method depending on self.place_in_file"""
        if self.place_in_file == "header":
            self.process_header_line(line)
        elif self.place_in_file == "column names":
            self.process_column_line(line)
        elif self.place_in_file == "data":
            self.process_data_line(line)
        else:  # just for debugging
            raise ReadError(f"place_in_file = {self.place_in_file}")
        self.n_line += 1

    def process_header_line(self, line):
        """Search line for important metadata and set the relevant attribute of self"""
        self.header_lines.append(line)
        name_match = re.search(regular_expressions["name"], line)
        if name_match:
            self.name = name_match.group(1)
            return
        N_head_match = re.search(regular_expressions["N_header_lines"], line)
        if N_head_match:
            self.N_header_lines = int(N_head_match.group(1))
            return
        timestamp_match = re.search(regular_expressions["tstamp"], line)
        if timestamp_match:
            self.tstamp = float(timestamp_match.group(1))
            return
        technique_match = re.search(regular_expressions["technique"], line)
        if technique_match:
            self.technique = technique_match.group(1)
            if self.technique in TECHNIQUE_CLASSES:
                if issubclass(TECHNIQUE_CLASSES[self.technique], self.measurement_class):
                    self.measurement_class = TECHNIQUE_CLASSES[self.technique]
            return
        timecol_match = re.search(regular_expressions["timecol"], line)
        if timecol_match:
            tcol = timecol_match.group(1)
            self.timecols[tcol] = []
            for vcol in timecol_match.group(2).split("' and '"):
                self.timecols[tcol].append(vcol)
            return
        aux_file_match = re.search(regular_expressions["aux_file"], line)
        if aux_file_match:
            aux_file_name = aux_file_match.group(1)
            aux_file = self.path_to_file.parent / aux_file_match.group(2)
            self.read_aux_file(aux_file, name=aux_file_name)
            return
        if " = " in line:
            key, value = line.strip().split(" = ")
            if key in bad_keys:
                return
            if key in ("name", "id"):
                return
            try:
                self.meas_as_dict[key] = json.loads(value)
            except json.decoder.JSONDecodeError:
                print(f"skipping the following line:\n{line}")
            return

        if self.N_header_lines and self.n_line >= self.N_header_lines - 2:
            self.place_in_file = "column names"

    def process_column_line(self, line):
        """Split the line to get the names of the file's data columns"""
        self.header_lines.append(line)
        self.column_names = [name.strip() for name in line.split(self.delim)]
        self.column_data.update({name: [] for name in self.column_names})
        self.place_in_file = "data"

    def process_data_line(self, line):
        """Split the line and append the numbers the corresponding data column arrays"""
        data_strings_from_line = line.strip().split(self.delim)
        for name, value_string in zip(self.column_names, data_strings_from_line):
            if value_string:
                try:
                    value = float(value_string)
                except ValueError:
                    # That is probably because different columns are different length.
                    #  so we just skip it!
                    continue
                    # raise ReadError(f"can't parse value string '{value_string}'")
                self.column_data[name].append(value)

    def read_aux_file(self, path_to_aux_file, name):
        """Read an auxiliary file and include its series list in the measurement"""
        spec = IxdatSpectrumReader().read(path_to_aux_file, name=name)
        self.aux_file_objects[name] = spec

    def print_header(self):
        """Print the file header including column names. read() must be called first."""
        header = "".join(self.header_lines)
        print(header)


def get_column_unit(column_name):
    """Return the unit name of an ixdat column, i.e the part of the name after the '/'"""
    unit_match = re.search(regular_expressions["unit"], column_name)
    if unit_match:
        unit_name = unit_match.group(1)
    else:
        unit_name = None
    return unit_name


class IxdatSpectrumReader(IxdatCSVReader):
    """A reader for ixdat spectra."""

    def read(self, path_to_file, name=None, cls=SpectrumSeries, **kwargs):
        """Read an ixdat spectrum.

        This reads the header with the process_line() function inherited from
        IxdatCSVReader. Then it uses pandas to read the data.

        Args:
            path_to_file (Path): The full absolute or relative path including extension
            name (str): The name of the measurement to return (defaults to path_to_file)
            cls (Spectrum subclass): The class of measurement to return. By default,
                cls will be determined from the technique specified in the header of
                path_to_file.
            **kwargs (dict): Key-word arguments are passed to ECMeasurement.__init__

        Returns cls: a Spectrum of type cls
        """
        path_to_file = Path(path_to_file)
        self.name = name or path_to_file.name
        with open(path_to_file, "r") as f:
            for line in f:
                if self.place_in_file == "header":
                    self.process_line(line)
                else:
                    break
        df = pd.read_csv(path_to_file, sep=",", header=self.N_header_lines - 2)
        if self.technique.endswith("spectrum"):
            # FIXME: in the future, this needs to cover all spectrum classes
            x_name, y_name = tuple(df.keys())
            x = df[x_name].to_numpy()
            y = df[y_name].to_numpy()
            cls = cls if issubclass(cls, Spectrum) else Spectrum
            return cls.from_data(  # see Spectrum.from_data()
                x,
                y,
                self.tstamp,
                x_name,
                y_name,
                name=self.name,
                technique=self.technique,
                reader=self,
            )

        elif self.technique.endswith("spectra"):
            # FIXME: in the future, this needs to cover all spectrum series classes
            names = {}
            units = {}
            swap_axes = False
            for line in self.header_lines:
                for line_start in ("values", "first row", "first column"):
                    if line.startswith(line_start):
                        t_x_or_y = re.search("([yxt])=", line).group(1)
                        names[t_x_or_y] = re.search(r"\'(.*)\'", line).group(1)
                        units[t_x_or_y] = re.search(r"\[(.*)\]", line).group(1)
                        if "row" in line_start and t_x_or_y == "t":  # check!
                            swap_axes = True
            z1 = np.array([float(key) for key in list(df.keys())[1:]])
            z1_and_y = df.to_numpy()
            z0 = z1_and_y[:, 0]
            y = z1_and_y[:, 1:]
            if swap_axes:
                # This is the case if the file was export with spectra_as_rows = False.
                t = z1
                x = z0
                y = y.swapaxes(0, 1)
            else:
                t = z0
                x = z1
            tseries = TimeSeries(
                name=names["t"], unit_name=units["t"], data=t, tstamp=self.tstamp
            )
            xseries = DataSeries(name=names["x"], unit_name=units["x"], data=x)
            field = Field(
                name=names["y"],
                unit_name=units["y"],
                data=y,
                axes_series=[tseries, xseries],
            )
            cls = cls if issubclass(cls, SpectrumSeries) else SpectrumSeries
            return cls.from_field(  # see SpectrumSeries.from_field()
                field, name=self.name, technique=self.technique, tstamp=self.tstamp
            )
