import warnings
from pathlib import Path
from ..exceptions import ReadError
from ..techniques import ECMeasurement
from ..data_series import TimeSeries, ValueSeries
import numpy as np


def _parse_ec_macro(path):
    """Parse a Nordic .EC_Macro file into a list of step dicts.

    Each line is a tab-separated sequence: step_type, key, value, key, value, ...
    Returns a list of dicts with a "type" key plus any key-value pairs.
    """
    steps = []
    with open(path) as f:
        for line in f:
            parts = [p for p in line.strip().split("\t") if p]
            if not parts:
                continue
            rest = parts[1:]
            if len(rest) == 1:
                step = {"type": parts[0], "value": rest[0]}
            else:
                step = {"type": parts[0]}
                for key, val in zip(rest[0::2], rest[1::2]):
                    step[key] = val
            steps.append(step)
    return steps


class NordicTDMSReader:
    """Reader for the .tdms files produced by Nordic Electrochemistry potentiostat"""

    def __init__(self):
        self.tdms_file = None  # The raw TdmsFile as imported by npTDMS package

    def read(self, path_to_file, *, cls=ECMeasurement, **kwargs):
        """Read the .tdms file and return the resulting measurement object

        Requires that the npTDMS package has been installed.
        See: https://nptdms.readthedocs.io/en/stable/quickstart.html

        Args:
            path_to_file (Path or str): Path to the .tdms file to read
            cls (Measurement subclass): the class to return an object of. Defaults to
                ECMeasurement
            kwargs: Additional key-word arguments are included in the dictionary used
                to initiate the measurement by cls.from_dict().
        """
        try:
            from nptdms import TdmsFile
        except ImportError as e:
            raise ReadError(
                "To read Nordic Electrochemistry .tdms files, ixdat uses "
                "the npTDMS package. Please install npTDMS and try again.\n"
                "see: https://nptdms.readthedocs.io/en/stable/quickstart.html "
                f"\noriginal error: \n{e}"
            )

        tdms_file = TdmsFile.read(path_to_file)
        self.tdms_file = tdms_file

        name = tdms_file.properties["name"]
        tstamp = tdms_file.properties["dateTime"].astype("datetime64[s]").astype("int")

        ec = tdms_file["EC"]

        t = ec["Time"][:]  # time in [s]
        V = ec["E"][:]  # raw potential
        V_unit = ec["E"].properties["unit_string"]
        i = ec["i"][:]  # raw current; convert A -> mA below
        i_unit = ec["i"].properties["unit_string"]
        if i_unit == "A":
            i = i * 1e3
            i_unit = "mA"

        tseries = TimeSeries(name="Time", unit_name="s", data=t, tstamp=tstamp)
        Vseries = ValueSeries(name="E", unit_name=V_unit, data=V, tseries=tseries)
        Iseries = ValueSeries(name="i", unit_name=i_unit, data=i, tseries=tseries)

        series_list = [tseries, Vseries, Iseries]

        tseries_EIS = None
        try:
            Z = ec["Z_E"][:]
        except KeyError:
            warnings.warn("No Z_E channel found in TDMS file.")
        else:
            # Z and phase have different lengths than t. Assume they are evenly spaced
            # over the same timespan as potential and current:
            t_EIS = np.linspace(t[0], t[-1], num=len(Z))
            tseries_EIS = TimeSeries(
                name="Time EIS", unit_name="s", data=t_EIS, tstamp=tstamp
            )
            Zseries = ValueSeries(
                name="Z_E",
                unit_name=ec["Z_E"].properties["unit_string"],
                data=Z,
                tseries=tseries_EIS,
            )
            series_list.append(Zseries)
        try:
            phase = ec["Phase_E"][:]
        except KeyError:
            warnings.warn("No Phase_E channel found in TDMS file.")
        else:
            if tseries_EIS is None:
                warnings.warn("Phase_E found without Z_E; cannot assign time series.")
            else:
                pseries = ValueSeries(
                    name="phase_E",
                    unit_name=ec["Phase_E"].properties["unit_string"],
                    data=phase,
                    tseries=tseries_EIS,
                )
                series_list.append(pseries)

        aliases = {"t": ["Time"], "raw_potential": ["E"], "raw_current": ["i"]}

        metadata = {}
        metadata["datetime"] = str(tdms_file.properties["dateTime"])
        macro_files = list(Path(path_to_file).parent.glob("*.EC_Macro"))
        if not macro_files:
            warnings.warn(
                "No .EC_Macro file found;"
                "experiment sequence not included in metadata."
            )
        else:
            if len(macro_files) > 1:
                warnings.warn(
                    f"Multiple .EC_Macro files found; using {macro_files[0].name}."
                )
            metadata["macro"] = _parse_ec_macro(macro_files[0])

        obj_as_dict = dict(
            name=name,
            technique="EC",
            tstamp=tstamp,
            series_list=series_list,
            aliases=aliases,
            metadata=metadata,
            reader=self,
        )
        obj_as_dict.update(kwargs)
        return cls.from_dict(obj_as_dict)
