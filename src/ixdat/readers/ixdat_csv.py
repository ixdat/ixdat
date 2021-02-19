"""Module defining the ixdat csv reader, so ixdat can read the files it exports."""

from pathlib import Path
import numpy as np
import re
from ..exceptions import ReadError
from ..data_series import ValueSeries, TimeSeries
from ..measurements import Measurement
from ..techniques import TECHNIQUE_CLASSES

regular_expressions = {
    "tstamp": r"tstamp = ([0-9\.]+)",
    "technique": r"technique = (\w+)",
    "N_header_lines": r"N_header_lines = ([0-9]+)",
    "backend_name": r"backend_name = (\w+)",
    "id": r"id = ([0-9]+)",
    "timecol": r"timecol '(.+)' for: (?:'(.+)')$",
    "unit": r"/ [(.+)]",
}


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
        """Initialize a Reader for .mpt files. See class docstring."""
        self.name = None
        self.path_to_file = None
        self.n_line = 0
        self.place_in_file = "header"
        self.header_lines = []
        self.tstamp = None
        self.N_header_lines = None
        self.timecols = {}
        self.column_names = []
        self.column_data = {}
        self.technique = None
        self.measurement_class = Measurement
        self.file_has_been_read = False
        self.measurement = None

    def read(self, path_to_file, name=None, cls=None, **kwargs):
        """Return an ECMeasurement with the data and metadata recorded in path_to_file

        This loops through the lines of the file, processing one at a time. For header
        lines, this involves searching for metadata. For the column name line, this
        involves creating empty arrays for each data series. For the data lines, this
        involves appending to these arrays. After going through all the lines, it
        converts the arrays to DataSeries.
        For .mpt files, there is one TimeSeries, with name "time/s", and all other data
        series are ValueSeries sharing this TimeSeries.
        Finally, the method returns an ECMeasurement with these DataSeries. The
        ECMeasurement contains a reference to the reader.

        Args:
            path_to_file (Path): The full abs or rel path including the ".mpt" extension
            **kwargs (dict): Key-word arguments are passed to ECMeasurement.__init__
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
        obj_as_dict = dict(
            name=self.name,
            technique=self.technique,
            reader=self,
            series_list=data_series_list,
            tstamp=self.tstamp,
        )
        obj_as_dict.update(kwargs)

        if issubclass(cls, self.measurement_class):
            self.measurement_class = cls

        if issubclass(self.measurement_class, TECHNIQUE_CLASSES["EC"]):
            # this is how ECExporter exports current and potential:
            obj_as_dict["raw_potential_names"] = ("raw potential / [V]",)
            obj_as_dict["raw_current_names"] = ("raw current / [mA]",)

        self.measurement = self.measurement_class.from_dict(obj_as_dict)
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
                if issubclass(
                    TECHNIQUE_CLASSES[self.technique], self.measurement_class
                ):
                    self.measurement_class = TECHNIQUE_CLASSES[self.technique]
            return
        timecol_match = re.search(regular_expressions["timecol"], line)
        if timecol_match:
            tcol = timecol_match.group(1)
            self.timecols[tcol] = []
            for vcol in timecol_match.group(2).split("' and '"):
                self.timecols[tcol].append(vcol)
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
                    raise ReadError(f"can't parse value string '{value_string}'")
                self.column_data[name].append(value)

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
