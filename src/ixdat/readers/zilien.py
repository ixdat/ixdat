import re
import pandas as pd
import numpy as np
from pathlib import Path
from ..data_series import DataSeries, TimeSeries, ValueSeries, Field
from ..techniques import ECMSMeasurement, MSMeasurement, ECMeasurement, Measurement
from ..techniques.ms import MSSpectrum
from .reading_tools import timestamp_string_to_tstamp, FLOAT_MATCH
from .ec_ms_pkl import measurement_from_ec_ms_dataset

ZILIEN_TIMESTAMP_FORM = "%Y-%m-%d %H_%M_%S"  # like 2021-03-15 18_50_10

ZILIEN_LEGACY_ALIASES = {
    # TODO: These should change to what Zilien calls them. Right now the alias's
    #   reflect the way the lagacy EC_MS code renames essential series
    "t": ["time/s"],
    "raw_potential": ["Ewe/V", "<Ewe>/V"],
    "raw_current": ["I/mA", "<I>/mA"],
    "cycle": ["cycle number"],
}

# TODO: When, in the future, Zilien files include the whole EC dataset, remove the
#    unflattering example presently in the docs.
#    https://github.com/ixdat/ixdat/pull/30/files#r810087496


class ZilienTSVReader:
    """Class for reading files saved by Spectro Inlets' Zilien software"""

    def read(self, path_to_file, cls=None, name=None, **kwargs):
        """Read a zilien file

        TODO: This is a hack using EC_MS to read the .tsv. Will be replaced.
        """

        from EC_MS import Zilien_Dataset

        if cls is Measurement:
            cls = ECMSMeasurement

        if "technique" not in kwargs:
            if issubclass(cls, ECMSMeasurement):
                kwargs["technique"] = "EC-MS"
            elif issubclass(cls, ECMeasurement):
                kwargs["technique"] = "EC"
            elif issubclass(cls, MSMeasurement):
                kwargs["technique"] = "MS"

        ec_ms_dataset = Zilien_Dataset(path_to_file)

        return measurement_from_ec_ms_dataset(
            ec_ms_dataset.data,
            cls=cls,
            name=name,
            reader=self,
            aliases=ZILIEN_LEGACY_ALIASES,
            **kwargs,
        )


class ZilienTMPReader:
    """A class for stitching the files in a Zilien tmp directory to an ECMSMeasurement

    This is necessary because Zilien often crashes, leaving only the tmp directory.
    This is less advanced but more readable than the Spectro Inlets stitching solution.
    """

    def __init__(self, path_to_tmp_dir=None):
        self.path_to_tmp_dir = Path(path_to_tmp_dir) if path_to_tmp_dir else None

    def read(self, path_to_tmp_dir, cls=None, **kwargs):
        """Make a measurement from all the single-value .tsv files in a Zilien tmp dir

        Args:
            path_to_tmp_dir (Path or str): The path to the tmp dir
            cls (Measurement class): Defaults to ECMSMeasurement
        """
        if path_to_tmp_dir:
            self.path_to_tmp_dir = Path(path_to_tmp_dir)
        cls = cls or ECMSMeasurement
        name = self.path_to_tmp_dir.parent.name
        timestamp_string = name[:19]  # the zilien timestamp is the first 19 chars
        tstamp = timestamp_string_to_tstamp(timestamp_string, form=ZILIEN_TIMESTAMP_FORM)
        series_list = []
        for tmp_file in self.path_to_tmp_dir.iterdir():
            series_list += series_list_from_tmp(tmp_file)
        obj_as_dict = {
            "name": name,
            "tstamp": tstamp,
            "series_list": series_list,
            "technique": "EC-MS",
            "reader": self,
        }
        obj_as_dict.update(kwargs)
        return cls.from_dict(obj_as_dict)


def series_list_from_tmp(path_to_file):
    """Return [ValueSeries, TimeSeries] with the data in a zilien tmp .tsv file"""
    file_name = Path(path_to_file).name
    timestamp_string = file_name[:19]  # the zilien timestamp form is 19 chars long
    tstamp = timestamp_string_to_tstamp(timestamp_string, form=ZILIEN_TIMESTAMP_FORM)
    column_match = re.search(r"\.([^\.]+)\.data", file_name)
    if not column_match:
        print(f"could not find column name in {path_to_file}")
        return []
    v_name = column_match.group(1)
    mass_match = re.search("M[0-9]+", v_name)
    if mass_match:
        v_name = mass_match.group()
        unit = "A"
    else:
        unit = None
    t_name = v_name + "-x"
    df = pd.read_csv(path_to_file, delimiter="\t", names=[t_name, v_name], header=0)
    t_data, v_data = df[t_name].to_numpy(), df[v_name].to_numpy()
    tseries = TimeSeries(name=t_name, unit_name="s", data=t_data, tstamp=tstamp)
    vseries = ValueSeries(name=v_name, unit_name=unit, data=v_data, tseries=tseries)
    return [tseries, vseries]


class ZilienSpectrumReader:
    """A reader for individual Zilien spectra
    TODO: A Zilien reader which loads all spectra at once in a SpectrumSeries object
    """

    def __init__(self, path_to_spectrum=None):
        self.path_to_spectrum = Path(path_to_spectrum) if path_to_spectrum else None

    def read(self, path_to_spectrum, cls=None, **kwargs):
        """Reat a Zilien spectrum.
        FIXME: This reader was written hastily and could be designed better.

        Args:
            path_to_tmp_dir (Path or str): the path to the tmp dir
            cls (Spectrum class): Defaults to MSSpectrum
            kwargs: Key-word arguments are passed on ultimately to cls.__init__
        """
        if path_to_spectrum:
            self.path_to_spectrum = Path(path_to_spectrum)
        cls = cls or MSSpectrum
        df = pd.read_csv(
            path_to_spectrum,
            header=9,
            delimiter="\t",
        )
        x_name = "Mass  [AMU]"
        y_name = "Current [A]"
        x = df[x_name].to_numpy()
        y = df[y_name].to_numpy()
        with open(self.path_to_spectrum, "r") as f:
            for i in range(10):
                line = f.readline()
                if "Mass scan started at [s]" in line:
                    tstamp_match = re.search(FLOAT_MATCH, line)
                    tstamp = float(tstamp_match.group())
        xseries = DataSeries(data=x, name=x_name, unit_name="m/z")
        field = Field(
            data=np.array(y),
            name=y_name,
            unit_name="A",
            axes_series=[
                xseries,
            ],
        )
        obj_as_dict = {
            "name": path_to_spectrum.name,
            "technique": "MS",
            "field": field,
            "reader": self,
            "tstamp": tstamp,
        }
        obj_as_dict.update(kwargs)
        return cls.from_dict(obj_as_dict)


if __name__ == "__main__":
    """Module demo here.

    To run this module in PyCharm, open Run Configuration and set
        Module name = ixdat.readers.zilien,
    and *not*
        Script path = ...
    """

    from pathlib import Path  # noqa
    from ixdat.measurements import Measurement  # noqa

    path_to_test_file = Path.home() / (
        "Dropbox/ixdat_resources/test_data/"
        # "zilien_with_spectra/2021-02-01 14_50_40.tsv"
        "zilien_with_ec/2021-02-01 17_44_12.tsv"
    )

    ecms_measurement = Measurement.read(
        reader="zilien",
        path_to_file=path_to_test_file,
    )

    ecms_measurement.plot_measurement()
