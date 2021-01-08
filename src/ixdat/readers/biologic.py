import re
import time
import numpy as np

from . import TECHNIQUE_CLASSES
from ..data_series import TimeSeries, ValueSeries
from ..exceptions import ReadError

ECMeasurement = TECHNIQUE_CLASSES["EC"]
delim = "\t"
t_str = "time/s"
time_string_form = "%m/%d/%Y %H:%M:%S"  # like 07/29/2020 10:31:03
regular_expressions = {
    "N_header_lines": "Nb header lines : (.+)\n",
    "timestamp_string": "Acquisition started on : (.+)\n",
    "loop": "Loop ([0-9]+) from point number ([0-9]+) to ([0-9]+)",
}


class BiologicMPTReader:
    def __init__(self):
        self.path_to_file = None
        self.n_line = 0
        self.place_in_file = "header"
        self.header_lines = []
        self.timestamp_string = None
        self.tstamp = None
        self.ec_technique = None
        self.N_header_lines = None
        self.column_names = []
        self.column_data = {}

    def read(self, path_to_file):
        self.path_to_file = path_to_file
        with open(path_to_file) as f:
            for line in f:
                self.process_line(line)

        if not t_str in self.column_data:
            return self  # debugging
            raise ReadError(
                f"{self} did not find any data for t_str='{t_str}'. "
                f"This reader only works for files with a '{t_str}' column"
            )
        tseries = TimeSeries(
            name=t_str,
            data=self.column_data[t_str],
            tstamp=self.tstamp,
            unit_name="s",
        )
        data_series_list = [tseries]
        for column_name, data in self.column_data.items():
            if column_name == t_str:
                continue
            vseries = ValueSeries(
                name=column_name,
                data=data,
                tseries=tseries,
                unit_name=self.get_column_unit(column_name),
            )
            data_series_list.append(vseries)

        return ECMeasurement(
            name=str(self.path_to_file),
            reader=self,
            series_list=data_series_list,
            tstamp=self.tstamp,
            ec_technique=self.ec_technique,
        )

    def process_line(self, line):
        if self.place_in_file == "header":
            self.process_header_line(line)
        elif self.place_in_file == "column names":
            self.process_column_line(line)
        elif self.place_in_file == "data":
            self.process_data_line(line)
        else:
            raise ReadError(f"place_in_file = {self.place_in_file}")
        self.n_line += 1

    def process_header_line(self, line):
        self.header_lines.append(line)
        if not self.N_header_lines:
            N_head_match = re.search(regular_expressions["N_header_lines"], line)
            if N_head_match:
                self.N_header_lines = int(N_head_match.group(1))
                return
        if self.n_line == 3:
            self.ec_technique = line.strip()
            return
        if not self.timestamp_string:
            timestamp_match = re.search(regular_expressions["timestamp_string"], line)
            if timestamp_match:
                self.timestamp_string = timestamp_match.group(1)
                self.tstamp = timestamp_string_to_tstamp(self.timestamp_string)
            return
        loop_match = re.search(regular_expressions["loop"], line)
        if loop_match:
            n = int(loop_match.group(1))
            start = int(loop_match.group(2))
            finish = int(loop_match.group(3))
            if not "loop_number" in self.column_data:
                self.column_data["loop_number"] = np.array([])
            self.column_data["loop_number"] = np.append(
                self.column_data["loop_number"], np.array([n] * (finish - start + 1))
            )
            return

        if self.N_header_lines and self.n_line >= self.N_header_lines - 2:
            self.place_in_file = "column names"

    def process_column_line(self, line):
        self.header_lines.append(line)
        self.column_names = line.strip().split(delim)
        self.column_data = {name: np.array([]) for name in self.column_names}
        self.place_in_file = "data"

    def process_data_line(self, line):
        data_strings_from_line = line.strip().split()
        for name, value_string in zip(self.column_names, data_strings_from_line):
            try:
                value = float(value_string)
            except ValueError:
                if "," in value_string:  # oh my god, why?!
                    value_string = value_string.replace(",", ".")
                if "E" in value_string:  # Biologic uses capital E for sci. notation.
                    value_string = value_string.replace(",", ".")
                    try:
                        value = float(value_string)
                    except ValueError:
                        raise ReadError(f"can't parse value string '{value_string}'")
            self.column_data[name] = np.append(self.column_data[name], value)

    def get_column_unit(self, column_name):
        if "/" in column_name:
            unit_name = column_name.split("/")[-1]
        else:
            unit_name = None
        return unit_name


def timestamp_string_to_tstamp(timestamp_string, form=time_string_form):
    struct = time.strptime(timestamp_string, form)
    tstamp = time.mktime(struct)
    return tstamp
