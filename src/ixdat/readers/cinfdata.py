"""Module defining readers for DTU Surfcat's legendary cinfdata system"""

from pathlib import Path
import numpy as np
from ..exceptions import ReadError
from ..data_series import ValueSeries, TimeSeries
from ..techniques import MSMeasurement
from .reading_tools import timestamp_string_to_tstamp


class CinfdataTXTReader:
    """A class that reads the text exported by cinfdata's text export functionality

    TODO: We should also have a reader class that downloads the data from cinfdata like
        `EC_MS`'s `download_cinfdata_set`:
        https://github.com/ScottSoren/EC_MS/blob/master/src/EC_MS/Data_Importing.py#L711

    Attributes:
        path_to_file (Path): the location and name of the file read by the reader
        n_line (int): the number of the last line read by the reader
        place_in_file (str): The last location in the file read by the reader. This
            is used internally to tell the reader how to parse each line. Options are:
            "header", "column names", and "data".
        header_lines (list of str): a list of the header lines of the files. This
            includes the column name line. The header can be nicely viewed with the
            print_header() function.
        tstamp (str): The unix time corresponding to t=0 for the measurement
        tstamp_list (list of float): list of epoch tstamps in the file's timestamp line
        column_tstamps (dict): The unix time corresponding to t=0 for each time column
        technique (str): The name of the technique
        column_names (list of str): The names of the data columns in the file
        t_and_v_cols (dict): {name: (tcol, vcol)} where name is the name of the
            ValueSeries (e.g. "M2"), tcol is the name of the corresponding time column
            in the file (e.g. "M2-x"), and vcol is the the name of the value column in
            the file (e.g. "M2-y).
        column_data (dict of str: np.array): The data in the file as a dict.
            Note that the np arrays are the same ones as in the measurement's DataSeries,
            so this does not waste memory.
        file_has_been_read (bool): This is used to make sure read() is only successfully
            called once by the Reader. False until read() is called, then True.
        measurement (Measurement): The measurement returned by read() when the file is
            read. self.measureemnt is None before read() is called.
    """

    delim = "\t"

    def __init__(self):
        """Initialize a Reader for cinfdata-exported text files. See class docstring."""
        self.name = None
        self.path_to_file = None
        self.n_line = 0
        self.place_in_file = "header"
        self.header_lines = []
        self.tstamp = None
        self.tstamp_list = []
        self.column_tstamps = {}
        self.column_names = []
        self.t_and_v_cols = {}
        self.column_data = {}
        self.technique = "MS"  # TODO: Figure out how to tell if it's something else
        self.measurement_class = MSMeasurement
        self.file_has_been_read = False
        self.measurement = None

    def read(self, path_to_file, name=None, cls=None, **kwargs):
        """Return an MSMeasurement with the data and metadata recorded in path_to_file

        This loops through the lines of the file, processing one at a time. For header
        lines, this involves searching for metadata. For the column name line, this
        involves creating empty arrays for each data series. For the data lines, this
        involves appending to these arrays. After going through all the lines, it
        converts the arrays to DataSeries.
        For cinfdata text files, each value column has its own timecolumn, and they are
        not necessarily all the same length.
        Finally, the method returns an ECMeasurement with these DataSeries. The
        ECMeasurement contains a reference to the reader.
        All attributes of this reader can be accessed from the
        measurement as `measurement.reader.attribute_name`.

        Args:
            path_to_file (Path): The full abs or rel path including the ".txt" extension
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

        data_series_list = []
        for name, (tcol, vcol) in self.t_and_v_cols.items():
            tseries = TimeSeries(
                name=tcol,
                unit_name=get_column_unit(tcol) or "s",
                data=self.column_data[tcol],
                tstamp=self.column_tstamps[tcol],
            )
            vseries = ValueSeries(
                name=name,
                data=self.column_data[vcol],
                tseries=tseries,
                unit_name=get_column_unit(vcol),
            )
            data_series_list.append(tseries)
            data_series_list.append(vseries)

        obj_as_dict = dict(
            name=self.name,
            technique=self.technique,
            reader=self,
            series_list=data_series_list,
            tstamp=self.tstamp,
        )
        # normally MSMeasurement requires mass aliases, but not cinfdata since it uses
        # the ixdat convention (actually, ixdat uses the cinfdata convention) of M<x>
        obj_as_dict.update(kwargs)

        if issubclass(cls, self.measurement_class):
            self.measurement_class = cls

        self.measurement = self.measurement_class.from_dict(obj_as_dict)
        self.file_has_been_read = True
        return self.measurement

    def process_line(self, line):
        """Call the correct line processing method depending on self.place_in_file"""
        if self.place_in_file == "header":
            self.process_header_line(line)
        elif self.place_in_file == "post_header":
            if line.strip():  # then we're in the column headers!
                self.process_column_line(line)
        elif self.place_in_file == "data":
            self.process_data_line(line)
        else:  # just for debugging
            raise ReadError(f"place_in_file = {self.place_in_file}")
        self.n_line += 1

    def process_header_line(self, line):
        """Search line for important metadata and set the relevant attribute of self"""
        self.header_lines.append(line)
        if not line.strip():  # the blank lines between the header and the column names
            self.place_in_file = "post_header"
        elif "Recorded at" in line:
            for s in line.split(self.delim):
                if "Recorded at" not in s:
                    self.tstamp_list.append(
                        timestamp_string_to_tstamp(
                            s.strip()[1:-1],  # remove edge whitespace and quotes.
                            form="%Y-%m-%d %H:%M:%S",  # like "2017-09-20 13:06:00"
                        )
                    )
            self.tstamp = self.tstamp_list[0]

    def process_column_line(self, line):
        """Split the line to get the names of the file's data columns"""
        self.header_lines.append(line)
        self.column_names = [name.strip() for name in line.split(self.delim)]
        self.column_data.update({name: [] for name in self.column_names})
        i = 0  # need a counter to map tstamps to timecols.
        for col in self.column_names:
            if col.endswith("-y"):
                name = col[:-2]
                tcol = f"{name}-x"
                if tcol not in self.column_names:
                    print(f"Warning! No timecol for {col}. Expected {tcol}. Ignoring.")
                    continue
                self.t_and_v_cols[name] = (tcol, col)
                self.column_tstamps[tcol] = self.tstamp_list[i]
                i += 1

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
    if column_name.startswith("M") and column_name.endswith("-y"):
        unit_name = "A"
    elif column_name.startswith("M") and column_name.endswith("-x"):
        unit_name = "s"
    else:
        # TODO: Figure out how cinfdata represents units for other stuff.
        #    see https://github.com/ixdat/ixdat/pull/30/files#r811432543, and
        #    https://github.com/CINF/cinfdata/blob/master/sym-files2/export_data.py#L125
        unit_name = None
    return unit_name
