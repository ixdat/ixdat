"""Module with possibly general-use tools for readers"""

from pathlib import Path
import time
import urllib.request
from ..config import config
from ..exceptions import ReadError
from ..measurements import TimeSeries, ValueSeries, ConstantValue


STANDARD_TIMESTAMP_FORM = "%d/%m/%Y %H:%M:%S"  # like '31/12/2020 23:59:59'
USA_TIMESTAMP_FORM = "%m/%d/%Y %H:%M:%S"  # like '12/31/2020 23:59:59'
FLOAT_MATCH = "[-]?\\d+[\\.]?\\d*(?:e[-+]?\\d+)?"  # matches floats like '5' or '-2.3e+5'
DEFAULT_READER_NAMES = {
    ".mpt": "biologic",
    ".mpr": "biologic",
    ".tsv": "zilien",
    ".xrdml": "xrdml",
    ".avg": "avantage",
}


def get_default_reader_name(path_to_file):
    """Return a default reader if available given a file's full name with suffix"""
    return DEFAULT_READER_NAMES.get(Path(path_to_file).suffix)


def timestamp_string_to_tstamp(
    timestamp_string,
    form=None,
    forms=(STANDARD_TIMESTAMP_FORM,),
):
    """Return the unix timestamp as a float by parsing timestamp_string

    Args:
        timestamp_string (str): The timestamp as read in the .mpt file
        form (str): The format string used by time.strptime (string-parse time). This is
            optional and overrides `forms` if given.
        forms (iter of str): The formats you want to try for time.strptime, defaults to
            the standard timestamp form.
    """
    if form:
        forms = (form,)
    for form in forms:
        try:
            return time.mktime(time.strptime(timestamp_string, form))
        except ValueError:
            continue

    raise ReadError(f"couldn't parse timestamp_string='{timestamp_string}')")


def prompt_for_tstamp(path_to_file, default="creation", form=STANDARD_TIMESTAMP_FORM):
    """Return the tstamp resulting from a prompt to enter a timestamp, or a default

    Args:
        path_to_file (Path): The file of the measurement that we're getting a tstamp for
        default (str or float): What to use as the tstamp if the user does not enter one.
            This can be a tstamp as a float, or "creation" to use the file creation time,
            or "now" to use `time.time()`.
        form (str): The specification string for the timestamp format. Defaults to
            `ixdat.readers.reading_tools.STANDARD_TIMESTAMP_FORM`
    """
    path_to_file = Path(path_to_file)

    if default == "creation":
        default_tstamp = path_to_file.stat().st_mtime
    elif default == "now":
        default_tstamp = time.time()
    elif type(default) in (int, float):
        default_tstamp = default
    else:
        raise TypeError("`default` must be a number or 'creation' or 'now'.")
    default_timestring = time.strftime(form, time.localtime(default_tstamp))

    tstamp = None
    timestamp_string = "Try at least once."
    while timestamp_string:
        timestamp_string = input(
            f"Please input the timestamp for the measurement at {path_to_file}.\n"
            f"Please use the format {form}.\n"
            "Enter nothing to use the default default,"
            f" '{default}', which is '{default_timestring}'."
        )
        if timestamp_string:
            try:
                tstamp = time.mktime(time.strptime(timestamp_string, form))
            except ValueError:
                print(
                    f"Could not parse '{timestamp_string}' according as '{form}'.\n"
                    f"Try again or enter nothing to use the default."
                )
            else:
                break
    return tstamp or default_tstamp


def series_list_from_dataframe(
    dataframe, time_name, tstamp, unit_finding_function, **kwargs
):
    """Return a list of DataSeries with the data in a pandas dataframe.

    The first series in the returned list is the one shared TimeSeries.

    Args:
        dataframe (pandas dataframe): The dataframe. Column names are used as series
            names, data is taken with series.to_numpy(). The dataframe can only have one
            TimeSeries (if there are more than one, pandas is probably not the right
            format anyway, since it requires columns be the same length).
        time_name (str): The name of the column to use as the TimeSeries
        tstamp (float): The timestamp
        unit_finding_function (function): A function which takes a column name as a
            string and returns its unit.
        kwargs: Additional key-word arguments are interpreted as constants to include
            in the data series list as `ConstantValue`s.
    """
    t = dataframe[time_name].to_numpy()
    tseries = TimeSeries(name=time_name, unit_name="s", data=t, tstamp=tstamp)
    data_series_list = [tseries]
    for column_name, series in dataframe.items():
        if column_name == time_name:
            continue
        data_series_list.append(
            ValueSeries(
                name=column_name,
                unit_name=unit_finding_function(column_name),
                data=series.to_numpy(),
                tseries=tseries,
            )
        )
    for key, value in kwargs.items():
        data_series_list.append(
            ConstantValue(name=key, unit_name="", data=value, tseries=tseries)
        )
    return data_series_list


def url_to_file(url, file_name="temp", directory=None):
    """Copy the contents of the url to a temporary file and return that file's Path."""
    directory = directory or config.ixdat_temp_dir
    suffix = "." + str(url).split(".")[-1]
    path_to_file = (directory / file_name).with_suffix(suffix)
    urllib.request.urlretrieve(url, path_to_file)
    return path_to_file


def get_file_list(path_to_file_start=None, part=None, suffix=None):
    """Get a list of files given their shared start of part.

    Use either `path_to_file_start` OR `part`.

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
            combined. Not to be used together with `path_to_file_start`.
        suffix (str): If a suffix is given, only files with the specified ending are
            added to the file list
    """
    file_list = []
    if path_to_file_start:
        path_to_file_start = Path(path_to_file_start)
        if path_to_file_start.is_dir():
            file_list = [f for f in path_to_file_start.iterdir() if f.is_file()]
        else:
            folder = path_to_file_start.parent
            base_name = path_to_file_start.name
            file_list = [f for f in folder.iterdir() if f.name.startswith(base_name)]
    elif part:
        folder = Path(part).parent
        part_name = Path(part).name
        file_list = [f for f in folder.iterdir() if part_name in f.name]
    if suffix:
        if not suffix.startswith("."):
            # So that the user can type e.g. `suffix="mpt"` as well as `suffix=".mpt"`
            suffix = "." + suffix
        file_list = [f for f in file_list if f.suffix == suffix]
    return file_list
