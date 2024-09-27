"""Readers for 'qexafs' and 'TRXRF' files exported by Diamond's B18-Core"""

import pandas as pd
from .reading_tools import timestamp_string_to_tstamp
from .. import Spectrum
from ..spectra import MultiSpectrum
from ..data_series import DataSeries, Field, TimeSeries, ValueSeries
from ..techniques.xrf import TRXRFMeasurement


class QexafsDATReader:
    def read(
        self,
        path_to_file,
        cls=Spectrum,
        technique=None,
        y_name=None,
        ref_name=None,
        **kwargs
    ):
        """Read the .dat file exported by Diamond B18-Core.

        Args:
            path_to_file (str or Path): The path to the .dat file
            cls (Spectrum subclass): The class to return an object of (if you use this
                reader via the read() method of a Spectrum sub-class, cls will
                automatically be that subclass. Defaults to `None`, meaning that the
                class will be determined by `technique` and will be `MultiSpectrum` if
                no technique is provided.
            technique (str): The technique to read. This may determine which data is
                incorporated into a `Spectrum` object. Options:
                    i. "XAS". Returns a `Spectrum` object with `x` from the first column
                       (e.g. "qexafs_energy") and `y` based on `y_name` and `ref_name`
                If no technique is requested and cls is None or `Spectrum`, a
                `MultiSpectrum` will be returned including all the data of the .dat file.
            y_name (str): The name of the column to use as the y data. Defaults to
                "QexafsFFI0" for "xas".
            ref_name (str): The name of the column by which the y data should be
                normalized. By default it becomes "I0" unless `y_name` ends with "I0".
                Set `ref_name` to "none" (as a string) to suppress normalization to "I0".
        """
        tstamp = None
        with open(path_to_file, "r") as f:
            for i, line in enumerate(f):
                if not line.startswith("#"):
                    # we're out of the header.
                    break
                if "Date:" in line:
                    timestamp_string = line.split("Date:")[1].strip()
                    tstamp = timestamp_string_to_tstamp(
                        timestamp_string,
                        forms=(
                            "%a, %d %b %Y %H:%M:%S BST",
                            "%a, %d %b %Y %H:%M:%S %Z"
                            # like "Fri, 13 May 2022 19:21:24 BST"
                            # or 'Mon, 6 Feb 2023 05:38:21 GMT'
                        ),
                    )

        df = pd.read_csv(path_to_file, sep="\t", header=i - 1)

        # these lines remove "#" and whitespace from the beginning of column names:
        # see https://pandas.pydata.org/docs/user_guide/text.html#string-methods
        df.columns = df.columns.str.strip("#")
        df.columns = df.columns.str.strip()

        x_name = df.keys()[0]
        x = df[x_name].to_numpy()
        x_unit = "eV"
        xseries = DataSeries(name=x_name, unit_name=x_unit, data=x)
        if technique == "XAS":
            y_name = y_name or "QexafsFFI0"
            if not y_name.endswith("I0"):
                ref_name = ref_name or "I0"
            y = df[y_name].to_numpy()
            if ref_name and not ref_name == "none":
                y = y / df[ref_name].to_numpy()
                unit_name = ""
            else:
                unit_name = "counts"
            yseries = DataSeries(name=y_name, unit_name=unit_name, data=y)
            return cls.from_series(
                xseries,
                yseries,
                tstamp=tstamp,
                technique=technique,
                name=str(path_to_file),
                **kwargs
            )

        fields = []
        for y_name in df.keys()[1:]:
            y = df[y_name].to_numpy()
            fields.append(
                Field(name=y_name, unit_name="", data=y, axes_series=[xseries])
                # TODO: a unit-finding function
            )
        return MultiSpectrum(
            name=str(path_to_file), fields=fields, tstamp=tstamp, technique=technique
        )


class B18TRXRFReader:
    def read(self, path_to_file, cls, seconds_per_x, **kwargs):
        """Read the .dat file exported by Diamond B18-Core, for time-resolved X-ray data
                (or fixed exciting enenrgy data).

        Args:
            path_to_file (str or Path): The path to the .dat file
            cls (Spectrum subclass): The class to return an object of (if you use this
                reader via the read() method of a Spectrum sub-class, cls will
                automatically be that subclass. Defaults to `None`, meaning that the
                class will be determined by `technique`.
            seconds_per_x (float): The time that corrresponds to 1 dummy energy
                (dummy energy is the x of spec, it could be at an interval of 1, 0.1,
                 0.001...) e.g. if the dummy energy is 0.001, 0.002.. but actual time is
                1s, 2s..., then seconds_per_x value should be 1000.
            **kwargs (dict): Key-word arguments are passed to cls.__init__
        """
        if issubclass(TRXRFMeasurement, cls):
            cls = TRXRFMeasurement
        qxafs_reader = QexafsDATReader()
        multi_spec = qxafs_reader.read(path_to_file)

        tstamp = multi_spec.tstamp
        t = multi_spec.x * seconds_per_x

        tseries = TimeSeries(
            name="time from dummy variable",
            unit_name="s",
            data=t,
            tstamp=tstamp,
        )
        series_list = [tseries]
        for spectrum in multi_spec.spectrum_list:
            vseries = ValueSeries(
                name=spectrum.name,
                unit_name=spectrum.yseries.unit_name,
                data=spectrum.y,
                tseries=tseries,
            )
            series_list.append(vseries)

        obj_as_dict = dict(
            name=path_to_file.name,
            series_list=series_list,
            tstamp=tstamp,
            technique="TRXRF",
        )
        obj_as_dict.update(kwargs)
        return cls.from_dict(obj_as_dict)
