"""This module implements the reader for the text export of Ivium's software"""

import re
from pathlib import Path
import pandas as pd
from ..techniques.ec import ECMeasurement
from .reading_tools import timestamp_string_to_tstamp, series_list_from_dataframe

IVIUM_ALIASES = {
    "raw_potential": ("E/V",),
    "raw_current": ("I/A",),
    "t": ("time/s",),
}


class IviumDataReader:
    """Class for reading single ivium files"""

    def read(self, path_to_file, cls=None, name=None, cycle_number=0, **kwargs):
        """Read the ASCII export from the Ivium software

        Args:
            path_to_file (Path): The full abs or rel path including the suffix (.txt)
            cls (Measurement subclass): The Measurement class to return an object of.
                Defaults to `ECMeasurement`.
            name (str): The name to use if not the file name
            cycle_number (int): The cycle number of the data in the file (default is 0)

            **kwargs (dict): Key-word arguments are passed to cls.__init__

        Returns:
            cls: technique measurement object with the ivium data
        """
        self.path_to_file = Path(path_to_file)
        name = name or self.path_to_file.name

        with open(self.path_to_file, "r") as f:
            timestring_line = f.readline()  # we need this for tstamp
            columns_line = f.readline()  # we need this to get the column names
            first_data_line = f.readline()  # we need this to check the column names
        tstamp = timestamp_string_to_tstamp(
            timestring_line.strip(),
            form="%d/%m/%Y %H:%M:%S",  # like '04/03/2021 19:42:30'
        )

        # ivium files do something really dumb. They add an extra column of data, which
        # looks like the measured potential (to complement 'E/V' which is presumably the
        # setpoint), but don't add the name of this column in the column name line.
        # So in order for pandas' csv reader to read it, we need assign a name to this
        # extra column (it becomes 'Unlabeled_1') and specify the column names.
        # Here we prepare the thus-corrected column name list, `column_names`:
        column_names = [col.strip() for col in columns_line.split(" ") if col.strip()]
        first_dat = [dat.strip() for dat in first_data_line.split(" ") if dat.strip()]
        if len(first_dat) > len(column_names):
            for i in range(len(first_dat) - len(column_names)):
                column_names.append(f"unlabeled_{i}")

        # And now we can read the data. Notice also the variable whitespace delimiter.
        dataframe = pd.read_csv(
            self.path_to_file, delimiter=r"\s+", header=1, names=column_names
        )

        # All that's left is getting the data from the dataframe into DataSeries and
        # into the Measurement, starting with the TimeSeries:

        data_series_list = series_list_from_dataframe(
            dataframe,
            time_name="time/s",
            tstamp=tstamp,
            unit_finding_function=get_column_unit,
            cycle=cycle_number,
        )
        # With the `series_list` ready, we prepare the Measurement dictionary and
        # return the Measurement object:
        obj_as_dict = dict(
            name=name,
            technique="EC",
            reader=self,
            aliases=IVIUM_ALIASES,
            series_list=data_series_list,
            tstamp=tstamp,
        )
        obj_as_dict.update(kwargs)

        if not issubclass(ECMeasurement, cls):
            cls = ECMeasurement
        return cls.from_dict(obj_as_dict)


class IviumDatasetReader:
    """Class for reading sets of ivium files exported together"""

    def read(self, path_to_file, cls=None, name=None, **kwargs):
        """Return a measurement containing the data of an ivium dataset,

        An ivium dataset is a group of ivium files exported together. They share a
        folder and a base name, and are suffixed "_1", "_2", etc.

        Args:
            path_to_file (Path or str): `Path(path_to_file).parent` is interpreted as the
                folder where the files of the ivium dataset is. `Path(path_to_file).name`
                up to the first "_" is interpreted as the shared start of the files in
                the dataset. You can thus use the base name of the exported files or
                the full path of any one of them.
            cls (Measurement class): The measurement class. Defaults to ECMeasurement.
            name (str): The name of the dataset. Defaults to the base name of the dataset
            kwargs: key-word arguments are included in the dictionary for cls.from_dict()

        Returns cls or ECMeasurement: A measurement object with the ivium data
        """
        self.path_to_file = Path(path_to_file)

        folder = self.path_to_file.parent
        base_name = self.path_to_file.name
        if re.search(r"_[0-9]", base_name):
            base_name = base_name.rpartition("_")[0]
        name = name or base_name

        # With two list comprehensions, we get the Measurement object for each file
        # in the folder who's name starts with base_name:
        all_file_paths = [f for f in folder.iterdir() if f.name.startswith(base_name)]
        component_measurements = [
            IviumDataReader().read(f, cls=cls, cycle_number=i)
            for i, f in enumerate(all_file_paths)
        ]

        # Now we append these using the from_component_measurements class method of the
        # right TechniqueMeasurement class, and return the result.
        if not cls:
            from ..techniques.ec import ECMeasurement

            cls = ECMeasurement
        measurement = cls.from_component_measurements(
            component_measurements, name=name, **kwargs
        )
        return measurement


def get_column_unit(column_name):
    """Return the unit name of an ivium column, i.e what follows the first '/'."""
    if "/" in column_name:
        return column_name.split("/", 1)[1]
