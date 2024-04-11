# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 16:42:03 2023

@author: Søren
"""
from pathlib import Path
import numpy as np
import pandas as pd
from ..techniques.ftir import FTIRSpectrumSeries
from ..data_series import DataSeries, TimeSeries, Field
from .reading_tools import timestamp_string_to_tstamp


class OpusFTIRReader:
    def read(
        self,
        path_to_file,
        name=None,
        cls=FTIRSpectrumSeries,
        time_first=None,
        time_last=None,
        suffix=".dpt",
    ):
        """Read a set of Opus FTIR spectra

        Args:
            path_to_file (Path): The full path shared by all the spectra-containing files
                Should include folders but not the number suffix or extention. For
                example "path/to/my/data/file_name" instead of
                "path/to/my_data/file_name_1.dpt". All files starting with the specified
                name in the specefied folder will be read and joined
            name (str): The name to use if not the file name
            cls (Spectrum subclass): The class to return an object of. Defaults to
                SpectrumSeries
            time_first (str): Timestamp of first spectrum, formatted
                like
                "05/12/2023 15:20:33.696 (GMT+0)"
            time_last (str): Timestamp of last spectrum, formatted like
                "05/12/2023 17:59:57.725 (GMT+0)"
            suffix (str): Suffix of the files containing text exported spectra. Defaults
                to ".dpt"

        """
        path_to_file = Path(path_to_file)
        name = name or path_to_file.stem

        if not issubclass(cls, FTIRSpectrumSeries):
            cls = FTIRSpectrumSeries

        files = []
        indeces = []
        for file in path_to_file.parent.iterdir():
            file_stem, number = file.stem.rsplit("_", maxsplit=1)
            if not file.suffix == suffix:
                continue
            if not file_stem == path_to_file.stem:
                continue
            indeces.append(int(number))
            files.append(file)

        sorted_indeces = np.argsort(indeces)
        sorted_files = [files[i] for i in sorted_indeces]

        x = None
        ys = []

        for file in sorted_files:
            df = pd.read_csv(file, names=["x", "y"])
            if x is None:
                x = df["x"].to_numpy()
            ys.append(df["y"])

        y_matrix = np.stack(ys)

        tstamp_first = timestamp_string_to_tstamp(
            time_first[:23], form="%d/%m/%Y %H:%M:%S.%f"
        )
        tstamp_last = timestamp_string_to_tstamp(
            time_last[:23], form="%d/%m/%Y %H:%M:%S.%f"
        )
        t = np.linspace(0, tstamp_last - tstamp_first, num=len(ys))

        xseries = DataSeries(name="wavenumber", unit_name="cm^-1", data=x)
        tseries = TimeSeries(name="time", unit_name="s", data=t, tstamp=tstamp_first)
        field = Field(
            name="intensity",
            unit_name=None,
            data=y_matrix,
            axes_series=[tseries, xseries],
        )

        ftir_series = cls(
            name=name,
            reader=self,
            technique="FTIR",
            tstamp=tstamp_first,
            field=field,
        )
        return ftir_series
