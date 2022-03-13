"""This module implements the reader for ascii exports from autolab's Nova software"""

import re
from pathlib import Path
import pandas as pd
from .reading_tools import (
    prompt_for_tstamp,
    series_list_from_dataframe,
    STANDARD_TIMESTAMP_FORM,
    timestamp_string_to_tstamp,
)

AUTOLAB_ALIASES = {
    "raw_potential": ("WE(1).Potential (V)",),
    "raw_current": ("WE(1).Current (A)",),
    "t": ("Time (s)",),
}


class NovaASCIIReader:
    """A reader for ascii files exported by Autolab's Nova software"""

    def read(
        self,
        path_to_file,
        cls=None,
        name=None,
        tstamp=None,
        timestring=None,
        timestring_form=STANDARD_TIMESTAMP_FORM,
        **kwargs
    ):
        """Read the ASCII export from Autolab's Nova software

        Args:
            path_to_file (Path): The full absolute or relative path including the suffix
            name (str): The name to use if not the file name
            cls (Measurement subclass): The Measurement class to return an object of.
                Defaults to `ECMeasurement` and should probably be a subclass thereof in
                any case.
            tstamp (float): timestamp of the measurement, if known
            timestring (str): timestring describing the timestamp of the measurement
            timestring_form (str): form of the timestring. Default is "%d/%m/%Y %H:%M:%S"
            **kwargs (dict): Key-word arguments are passed to cls.__init__
        """
        self.path_to_file = Path(path_to_file)
        name = name or self.path_to_file.name
        if not tstamp:
            if timestring:
                tstamp = timestamp_string_to_tstamp(timestring, form=timestring_form)
            else:
                tstamp = prompt_for_tstamp(self.path_to_file)

        dataframe = pd.read_csv(self.path_to_file, delimiter=";")

        data_series_list = series_list_from_dataframe(
            dataframe, "Time (s)", tstamp, get_column_unit
        )
        obj_as_dict = dict(
            name=name,
            technique="EC",
            reader=self,
            aliases=AUTOLAB_ALIASES,
            series_list=data_series_list,
            tstamp=tstamp,
        )
        obj_as_dict.update(kwargs)

        if not cls:
            from ..techniques.ec import ECMeasurement

            cls = ECMeasurement
        return cls.from_dict(obj_as_dict)


def get_column_unit(column_name):
    """Return the unit name of an autolab column, i.e the last part of the name in ()"""
    unit_match = re.search(r"\((.+)\)$", column_name)
    if unit_match:
        unit_name = unit_match.group(1)
    else:
        unit_name = None
    return unit_name
