"""This module implements the autolab txt reader"""
import re
from pathlib import Path
import pandas as pd
from ..data_series import TimeSeries, ValueSeries
from .reading_tools import prompt_for_tstamp


class AutolabTXTReader:
    """A reader for ascii files exported by Autolab's Nova software"""

    def read(self, path_to_file, cls=None, name=None, **kwargs):
        """read the ascii export from Autolab's Nova software

        Args:
            path_to_file (Path): The full abs or rel path including the suffix (.txt)
            name (str): The name to use if not the file name
            cls (Measurement subclass): The Measurement class to return an object of.
                Defaults to `ECMeasurement` and should probably be a subclass thereof in
                any case.
            **kwargs (dict): Key-word arguments are passed to cls.__init__
        """
        self.path_to_file = Path(path_to_file)
        name = name or self.path_to_file.name
        tstamp = prompt_for_tstamp(self.path_to_file)

        dataframe = pd.read_csv(self.path_to_file, delimiter=";")

        t_str = "Time (s)"
        tseries = TimeSeries(
            name=t_str, unit_name="s", data=dataframe[t_str].to_numpy(), tstamp=tstamp
        )
        data_series_list = [tseries]
        for column_name, series in dataframe.items():
            if column_name == t_str:
                continue
            data_series_list.append(
                ValueSeries(
                    name=column_name,
                    unit_name=get_column_unit(column_name),
                    data=series.to_numpy(),
                    tseries=tseries,
                )
            )
        obj_as_dict = dict(
            name=name,
            technique="EC",
            reader=self,
            raw_potential_names=("WE(1).Potential (V)",),
            raw_current_names=("WE(1).Current (A)",),
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
