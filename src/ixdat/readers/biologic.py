"""This module implements the Reader for .mpt files made by BioLogic's EC-Lab software

Demonstrated/tested at the bottom under `if __name__ == "__main__":`
"""

import re
from pathlib import Path
import warnings
import numpy as np
import pandas as pd

from . import TECHNIQUE_CLASSES
from .reading_tools import timestamp_string_to_tstamp, series_list_from_dataframe
from ..data_series import TimeSeries, ValueSeries, ConstantValue
from ..exceptions import ReadError, SeriesNotFoundError

ECMeasurement = TECHNIQUE_CLASSES["EC"]
delim = "\t"
t_str = "time/s"
timestamp_form_strings = [
    "%m/%d/%Y %H:%M:%S",  # like 07/29/2020 10:31:03
    "%m-%d-%Y %H:%M:%S",  # like 01-31-2020 10:32:02
]
regular_expressions = {
    "N_header_lines": "Nb header lines : (.+)\n",
    "timestamp_string": "Acquisition started on : (.+)\n",
    "loop": "Loop ([0-9]+) from point number ([0-9]+) to ([0-9]+)",
}

BIOLOGIC_ALIASES = {
    "t": ["time/s"],
    "raw_potential": ["Ewe/V", "<Ewe>/V"],
    "raw_CE_potential": ["Ece/V", "<Ece>/V"],
    "raw_current": ["I/mA", "<I>/mA"],
    "cycle": ["cycle number"],
}


def fix_WE_potential(measurement):
    """Fix column of zeros in "<Ewe>/V" sometimes exported by EC Lab for CP measurements.

    Some Biologic potentiostats / EC-Lab versions sometimes export a column of zeros for
    "<Ewe>/V" in the .mpt files in chronopotentiometry measurements. This function
    replaces the series of zeros with the correct potential by adding the counter
    electrode potential ("<Ece>/V") and cell potential ("Ewe-Ece/V").

    This function is not called automatically - it needs to be called manually on the
    measurements loaded from the afflicted files. It requires that the counter electrode
    potential was recorded.

    Args:
        measurement(ECMeasurement): The measurement with the column to be replaced
    """
    WE_series = measurement["<Ewe>/V"]
    try:
        CE_data = measurement.grab_for_t("<Ece>/V", WE_series.t)
    except SeriesNotFoundError:
        print(
            "The function `fix_WE_potential` requires that the counter electrode "
            "potential was recorded, and is in the file as '<Ece>/V."
        )
        raise

    cell_potential_data = measurement.grab_for_t("Ewe-Ece/V", WE_series.t)

    WE_potential = cell_potential_data + CE_data
    WE_series = ValueSeries(
        name="<Ewe>/V",
        unit_name="V",
        data=WE_potential,
        tseries=WE_series.tseries,
    )
    measurement.replace_series("<Ewe>/V", WE_series)
    measurement.clear_cache()


class BiologicReader:
    """A class to read .mpt files written by Biologic's EC-Lab.

    read() is the important method - it takes the path to the mpt file as argument
    and returns an ECMeasurement object (ec_measurement) representing that file.
    The ECMeasurement contains a reference to the BiologicMPTReader object, as
    ec_measurement.reader. This makes available all the following stuff, likely
    useful for debugging.

    Attributes:
        file_has_been_read (bool): This is used to make sure read() is only successfully
            called once by the Reader. False until read() is called, then True.
        measurement (Measurement): The measurement returned by read() when the file is
            read. self.measurement is None before read() is called.
        measurement_name (str): The name of the measurement being read
        path_to_file (Path): the location and name of the file read by the reader
        tstamp (float): The unix time corresponding to t=0, parsed from timestamp_string
        measurement_class (class): Type of measurement to return
        data_series_list (list of DataSeries): Data series of the measurement being read
        aliases (dict): Aliases for data series in the measurement being read
        tseries (TimeSeries): Time series for the returned measurement (biologic files
            have one shared time variable)
        ec_technique (str): The name of the electrochemical sub-technique, i.e.
            "Cyclic Voltammetry Advanced", etc.
        n_line (int): the number of the last line read by the reader
        place_in_file (str): The last location in the file read by the reader. This
            is used internally to tell the reader how to parse each line. Options are:
            "header", "column names", and "data".
        header_lines (list of str): a list of the header lines of the files. This
            includes the column name line. The header can be nicely viewed with the
            print_header() function.
        timestamp_string (str): The string identified to represent the t=0 time of the
            measurement recorded in the file.
        N_header_lines (int): The number of lines in the header of the file
        column_names (list of str): The names of the data columns in the file
        column_data (dict of str: np.array): The data in the file as a dict.
            Note that the np arrays are the same ones as in the measurement's DataSeries,
            so this does not waste memory.
        df (Pandas DataFrame): The data from a .mpr as read by an external package
    """

    def __init__(self):
        """Initialize a Reader for .mpt files. See class docstring."""
        super().__init__()

        # These are general for readers
        self.file_has_been_read = False
        self.measurement = None
        self.measurement_name = None
        self.path_to_file = None
        self.tstamp = None

        # These are common to .mpt and .mpr reading
        self.tseries = None
        self.data_series_list = []
        self.aliases = {}
        self.measurement_class = ECMeasurement
        self.ec_technique = None

        # These are special to .mpt reading
        self.n_line = 0
        self.place_in_file = "header"
        self.header_lines = []
        self.timestamp_string = None
        self.N_header_lines = None
        self.column_names = []
        self.column_data = {}

        # this is special to .mpr reading:
        self.df = None

    def read(self, path_to_file, name=None, cls=ECMeasurement, **kwargs):
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
            path_to_file (Path): The full absolute or relative path including the suffix
                (".mpt" or ".mpr")
            **kwargs (dict): Key-word arguments are passed to ECMeasurement.__init__
            name (str): The name to use if not the file name
            cls (Measurement subclass): The Measurement class to return an object of.
                Defaults to `ECMeasurement` and should probably be a subclass thereof in
                any case.
            **kwargs (dict): Key-word arguments are passed to cls.__init__
        """
        self.path_to_file = Path(path_to_file)
        self.measurement_name = name or self.path_to_file.name
        if issubclass(ECMeasurement, cls):
            cls = ECMeasurement
        self.measurement_class = cls

        if self.file_has_been_read:
            print(
                f"This {self.__class__.__name__} has already read {self.path_to_file}."
                " Returning the measurement resulting from the original read. "
                "Use a new Reader if you want to read another file."
            )
            return self.measurement

        if self.path_to_file.suffix == ".mpr":
            self.series_list_from_mpr()
        else:  # if the suffix is not ".mpr", assume it's a ".mpt"
            # read the .mpt file to get the series list using the methods of this class.
            self.series_list_from_mpt()
        self._update_aliases_and_ensure_essential_series()

        obj_as_dict = dict(
            name=self.measurement_name,
            technique="EC",
            reader=self,
            series_list=self.data_series_list,
            tstamp=self.tstamp,
            ec_technique=self.ec_technique,
            aliases=self.aliases,
        )
        obj_as_dict.update(kwargs)

        self.measurement = cls.from_dict(obj_as_dict)
        self.file_has_been_read = True

        return self.measurement

    def series_list_from_mpt(self, path_to_file=None):
        """Read a .mpt file to generate the reader's `data_series_list`"""
        path_to_file = Path(path_to_file or self.path_to_file)

        with open(path_to_file, "r", encoding="ISO-8859-1") as f:
            for line in f:
                self._process_line(line)

        for name in self.column_names:
            self.column_data[name] = np.array(self.column_data[name])

        if t_str not in self.column_data:
            raise ReadError(
                f"{self} did not find any data for t_str='{t_str}'. "
                f"This reader only works for files with a '{t_str}' column"
            )
        self.tseries = TimeSeries(
            name=t_str,
            data=self.column_data[t_str],
            tstamp=self.tstamp,
            unit_name="s",
        )
        self.data_series_list = [self.tseries]
        for column_name, data in self.column_data.items():
            if column_name == t_str:
                continue
            vseries = ValueSeries(
                name=column_name,
                data=data,
                tseries=self.tseries,
                unit_name=get_column_unit_name(column_name),
            )
            self.data_series_list.append(vseries)

    def _process_line(self, line):
        """Call the correct line processing method depending on self.place_in_file"""
        if self.place_in_file == "header":
            self._process_header_line(line)
        elif self.place_in_file == "column names":
            self._process_column_line(line)
        elif self.place_in_file == "data":
            self._process_data_line(line)
        else:  # just for debugging
            raise ReadError(f"place_in_file = {self.place_in_file}")
        self.n_line += 1

    def _process_header_line(self, line):
        """Search line for important metadata and set the relevant attribute of self"""
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
                self.tstamp = timestamp_string_to_tstamp(
                    self.timestamp_string, forms=BIOLOGIC_TIMESTAMP_FORMS
                )
            return
        loop_match = re.search(regular_expressions["loop"], line)
        if loop_match:
            # print(f"loop specified on line='{line}'")  # debugging
            n = int(loop_match.group(1))
            start = int(loop_match.group(2))
            finish = int(loop_match.group(3))
            if "loop_number" not in self.column_data:
                self.column_data["loop_number"] = np.array([])
            self.column_data["loop_number"] = np.append(
                self.column_data["loop_number"], np.array([n] * (finish - start + 1))
            )
            return

        if self.N_header_lines and self.n_line >= self.N_header_lines - 2:
            self.place_in_file = "column names"

    def _process_column_line(self, line):
        """Split the line to get the names of the file's data columns"""
        self.header_lines.append(line)
        self.column_names = line.strip().split(delim)
        self.column_data.update({name: [] for name in self.column_names})
        self.place_in_file = "data"

    def _process_data_line(self, line):
        """Split the line and append the numbers the corresponding data column arrays"""
        data_strings_from_line = line.strip().split()
        for name, value_string in zip(self.column_names, data_strings_from_line):
            try:
                value = float(value_string)
            except ValueError:
                if "," in value_string:  # oh my god, why?!
                    value_string = value_string.replace(",", ".")
                try:
                    value = float(value_string)
                except ValueError:
                    # Some old .mpt files have the value "XXXX" in some columns. To
                    # keep the reader from crashing, we warn and interpret it as 0:
                    warnings.warn(
                        f"Can't parse value string '{value_string}' in "
                        f"column '{name}'. Using Value=0"
                    )
                    value = 0
            self.column_data[name].append(value)

    def series_list_from_mpr(self, path_to_file=None):
        """Read a .mpr file to generate the reader's `data_series_list`

        Raises: ReadError, if no external package succeeds in reading the file. The
            error message includes instructions to install each external package or the
            error that is raised when attempting to use it.
        """
        # read the .mpr file to get the series list using a function which calls
        # an external package.
        warnings.warn(
            "Reading .mpr files is discouraged.\n"
            "We suggest to use the .mpt file if you can.\n"
            "You can set EC-Lab to export .mpt files automatically under "
            "Advanced Settings.\n"
            "ixdat is able to extract more "
            "useful data and metadata (e.g. loop numbers) from .mpt's\n"
        )

        # First, try with `galvani`
        try:
            self.series_list_from_mpr_galvani()
        except ImportError:
            # Nudge the user towards .mpt since we are better at those:
            e_galvani = (
                "To read biologic binary files (.mpr), ixdat first tries to use the "
                "`galvani` package, which could not be imported. "
                "See https://pypi.org/project/galvani/ \n"
                "Install `galvani` and try again."
            )
        except Exception as e:
            e_galvani = e
        else:  # we're good!
            print(f"read with `galvani`: {self.path_to_file}")
            return

        # Then, try with `eclabfiles`
        try:
            self.series_list_from_mpr_eclabfiles()
        except ImportError:
            # Nudge the user towards .mpt since we are better at those:
            e_eclabfiles = (
                "To read biologic binary files (.mpr), ixdat also tries to use the "
                "`eclabfiles` package, which could not be imported. "
                "See https://pypi.org/project/eclabfiles/ \n"
                "Install `eclabfiles` and try again."
            )
        except Exception as e:
            e_eclabfiles = e
        else:  # we're good!
            print(f"read with `eclabfiles`: {self.path_to_file}")
            return

        raise ReadError(
            f"couldn't read {path_to_file} with `galvani` or `eclabfiles`\n"
            f"Consider reading the .mpt file instead. You can set EC-Lab to export "
            ".mpt files automatically under Advanced Settings.\n"
            f"--- `galvani` error:\n{e_galvani}\n"
            f"--- `eclabfiles` error:\n{e_eclabfiles}"
        )

    def series_list_from_mpr_galvani(self, path_to_file=None):
        """Read a biologic .mpr file

        This makes use of the package `galvani`.
        See: https://github.com/echemdata/galvani.
        The dataframe read in by `eclabfiles.to_df()` is stored in the returned
        measurement `meas` as `meas.reader.df`.

        Args:
            path_to_file (Path): The full abs or rel path including the ".mpr" extension
        """
        from galvani import BioLogic

        path_to_file = Path(path_to_file or self.path_to_file)

        with open(path_to_file, "rb") as mpr_binary_file:
            mpr_file = BioLogic.MPRfile(mpr_binary_file)

        self.tstamp = mpr_file.timestamp.timestamp()
        self.df = pd.DataFrame(mpr_file.data)

        # Build the time series from the dataframe
        self.data_series_list = series_list_from_dataframe(
            self.df,
            time_name=t_str,
            tstamp=self.tstamp,
            unit_finding_function=get_column_unit_name,
        )
        self.tseries = self.data_series_list[0]

    def series_list_from_mpr_eclabfiles(self, path_to_file=None):
        """Read a biologic .mpr file

        This makes use of the package `eclabfiles`.
        See: https://github.com/vetschn/eclabfiles.
        The dataframe read in by `eclabfiles.to_df()` is stored in the returned
        measurement `meas` as `meas.reader.df`.

        Args:
            path_to_file (Path): The full abs or rel path including the ".mpr" extension
        """
        import eclabfiles  # not a requirement of ixdat, so we import it here.

        self.df = eclabfiles.to_df(str(path_to_file or self.path_to_file))

        units = self.df.attrs["units"]
        self.tstamp = self.df.attrs["log"]["posix_timestamp"]
        self.ec_technique = self.df.attrs["settings"]["technique"]

        # Build the time series from the dataframe
        self.data_series_list = series_list_from_dataframe(
            self.df,
            # the unit-free name of the time column, i.e. "time":
            time_name=t_str.split("/")[0],
            tstamp=self.tstamp,
            # to find a column's unit, look it up in the dictionary `units`:
            unit_finding_function=units.get,
        )
        self.tseries = self.data_series_list[0]

    def _update_aliases_and_ensure_essential_series(self):
        """A helper function completing the data_series_list for biologic readers."""
        # First, check if any of the biologic aliases are
        for data_series in self.data_series_list:
            for name, alias_list in BIOLOGIC_ALIASES.items():
                # We need to make sure the .mpr (unit-less) names are find-able:
                name_slash_unit = data_series.name + "/" + data_series.unit_name
                # In that case, this adds, e.g. "time" to the aliases for essential
                # series "t", since BIOLOGIC_ALIASES["t"] includes "time/s"
                if data_series.name in alias_list or name_slash_unit in alias_list:
                    if name in self.aliases:
                        self.aliases[name].append(data_series.name)
                    else:
                        self.aliases[name] = [data_series.name]
                    break

        # The following series need to be there:
        essential_series = set(self.measurement_class.essential_series_names).union(
            {"cycle number", "Ns"}
        )
        # To ensure this, we put a ConstantValue of zero in for the series
        series_names = [s.name for s in self.data_series_list]
        for series_name in essential_series:
            if series_name not in series_names and series_name not in self.aliases:
                name_0 = series_name + "=0"  # A series name to indicate it is set to 0.
                self.data_series_list.append(
                    ConstantValue(
                        name=name_0, unit_name="", data=0, tseries=self.tseries
                    )
                )
                # Add an alias to make it find-able with the essential series name:
                self.aliases[series_name] = [name_0]

    def print_header(self):
        """Print the file header including column names. read() must be called first."""
        header = "".join(self.header_lines)
        print(header)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.path_to_file})"


def get_column_unit_name(column_name):
    """Return the unit name of a .mpt column, i.e the part of the name after the '/'"""
    if "/" in column_name:
        unit_name = column_name.split("/")[-1]
    else:
        unit_name = None
    return unit_name


# Formats by which timestamps are saved in various EC-Labs # with example encountered
BIOLOGIC_TIMESTAMP_FORMS = (
    "%m-%d-%Y %H:%M:%S",  # like 01-31-2020 10:32:02
    "%m/%d/%Y %H:%M:%S",  # like 07/29/2020 10:31:03
    "%m-%d-%Y %H:%M:%S.%f",  # (anticipated)
    "%m/%d/%Y %H:%M:%S.%f",  # like 04/27/2021 11:35:39.227 (EC-Lab v11.34)
    "%m/%d/%Y %H.%M.%S",  # like 01/31/2022 11.19.17
    "%m/%d/%Y %H.%M.%S.%f",  # like 09/08/2022 13.08.17.338 (EC-Lab v11.43)
)

# This tuple contains variable names encountered in .mpt files. The tuple can be used by
#   other modules to tell which data is from biologic.
BIOLOGIC_COLUMN_NAMES = (
    "mode",
    "ox/red",
    "error",
    "control changes",
    "time/s",
    "control/V",
    "Ewe/V",
    "<I>/mA",
    "(Q-Qo)/C",
    "P/W",
    "loop number",
    "I/mA",
    "control/mA",
    "Ns changes",
    "counter inc.",
    "cycle number",
    "Ns",
    "(Q-Qo)/mA.h",
    "dQ/C",
    "Q charge/discharge/mA.h",
    "half cycle",
    "Capacitance charge/µF",
    "Capacitance discharge/µF",
    "dq/mA.h",
    "Q discharge/mA.h",
    "Q charge/mA.h",
    "Capacity/mA.h",
    "file number",
    "file_number",
    "Ece/V",
    "Ewe-Ece/V",
    "<Ece>/V",
    "<Ewe>/V",
    "Energy charge/W.h",
    "Energy discharge/W.h",
    "Efficiency/%",
    "Rcmp/Ohm",
)


if __name__ == "__main__":
    """Module demo here.

    To run this module in PyCharm, open Run Configuration and set
        Module name = ixdat.readers.biologic,
    and *not*
        Script path = ...
    """

    from matplotlib import pyplot as plt
    from ixdat.measurements import Measurement

    test_data_dir = Path(__file__).parent.parent.parent.parent / "test_data/biologic"

    path_to_test_file = test_data_dir / "Pt_poly_cv.mpt"

    ec_measurement = Measurement.read(
        reader="biologic",
        path_to_file=path_to_test_file,
    )

    t, v = ec_measurement.grab_potential(tspan=[0, 100])

    ec_measurement.tstamp -= 20
    t_shift, v_shift = ec_measurement.grab_potential(tspan=[0, 100])

    fig, ax = plt.subplots()
    ax.plot(t, v, "k", label="original tstamp")
    ax.plot(t_shift, v_shift, "r", label="shifted tstamp")
    ax.legend()
