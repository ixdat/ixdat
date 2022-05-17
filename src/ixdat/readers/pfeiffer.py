"""This module implements the reader for Pfeiffer Vacuum's PV Mass Spec software"""

import re
from pathlib import Path
import pandas as pd
from .reading_tools import timestamp_string_to_tstamp, series_list_from_dataframe
from ..techniques import MSMeasurement


class PVMassSpecReader:
    """A reader for (advanced) MID files exported from PVMassSpec ('... - Bin.dat')"""

    def read(self, path_to_file, cls=None, name=None, **kwargs):
        """Return a Measurement with the (advanced) MID data in the PVMassSpec file

        Args:
            path_to_file (Path or str): a path to the file exported by PVMassSpec with
                (advanced) MID data. This file is typically exported with a name that
                ends in '- Bin.dat', and with the timestamp in the file name. Note
                that the file can be renamed, as the original name is in the file,
                and the timestamp is read from there.
            cls (Measurement subclass): The technique class of which to return an object.
                Defaults to MSMeasurement.
            name (str): The name of the measurement. Defaults to Path(path_to_file).name
            kwargs: key-word args are used to initiate the measurement via cls.as_dict()

        Return cls: The measurement object
        """
        self.path_to_file = Path(path_to_file)
        name = name or self.path_to_file.name
        with open(path_to_file, "r") as f:
            # timestamp is on the the third line, which we select here:
            tstamp_line = [f.readline() for _ in range(3)][-1]
        tstamp = timestamp_string_to_tstamp(
            tstamp_line.split(".")[-2][-19:],  # last 19 characters before the last '.'
            form="%m-%d-%Y %H'%M'%S",  # like "03-02-2021 12'58'40"
        )
        df = pd.read_csv(self.path_to_file, header=6, delimiter="\t")
        # PV MassSpec calls masses <x>_amu, information we need to pass on to
        # MSMeasurement, so that the data will be accessible by the 'M<x>' mass string.
        aliases = {
            mass_from_column_name(key): [key] for key in df.keys() if "_amu" in key
        }
        series_list = series_list_from_dataframe(
            df,
            tstamp=tstamp,
            time_name="Time Relative (sec)",
            unit_finding_function=get_column_unit,
        )
        meas_as_dict = {
            "name": name,
            "tstamp": tstamp,
            "series_list": series_list,
            "aliases": aliases,
            "technique": "MS",
        }
        meas_as_dict.update(kwargs)
        cls = cls or MSMeasurement
        return cls.from_dict(meas_as_dict)


class PVMassSpecScanReader:
    """A reader for mass spectra files exported from PVMassSpec ('... - Scan.dat')"""

    pass


def mass_from_column_name(mass):
    """Return the PVMassSpec mass 'M<x>' given the column name '<x>_amu' as string"""
    return f"M{mass.split('_')[0]}"


def get_column_unit(column_name):
    """Return the unit name of an ivium column, i.e what follows the first '/'."""
    unit_match = re.search(r"\((.*)\)$", column_name)
    if unit_match:
        return unit_match.group(1)
    elif "amu" in column_name:
        return "A"
