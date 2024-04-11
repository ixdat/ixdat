from pathlib import Path
from . import TECHNIQUE_CLASSES
import pickle
from ..data_series import TimeSeries, ValueSeries
from ..measurements import Measurement
from .biologic import BIOLOGIC_COLUMN_NAMES, get_column_unit_name


ECMSMeasurement = TECHNIQUE_CLASSES["EC-MS"]


class EC_MS_CONVERTER:
    """Imports old .pkl files obtained from the legacy EC-MS package"""

    def __init__(self):
        print("Reader of old ECMS .pkl files")

    def read(self, file_path, cls=None, **kwargs):
        """Return an ECMSMeasurement with the data recorded in path_to_file
        Most of the work is done by module-level function measurement_from_ec_ms_dataset

        Args:
            path_to_file (Path): The full abs or rel path including the
                ".pkl" extension.
        """
        with open(file_path, "rb") as f:
            ec_ms_dict = pickle.load(f)

        return measurement_from_ec_ms_dataset(
            ec_ms_dict,
            name=Path(file_path).name,
            cls=cls,
            reader=self,
            technique="EC-MS",
            **kwargs,
        )


def measurement_from_ec_ms_dataset(
    ec_ms_dict,
    name=None,
    cls=ECMSMeasurement,
    reader=None,
    technique=None,
    **kwargs,
):
    """Return an ixdat Measurement with the data from an EC_MS data dictionary.

    This loops through the keys of the EC-MS dict and searches for MS and
    EC data. Names the dataseries according to their names in the original
    dict. Omits any other data as well as metadata.

    Args:
        ec_ms_dict (dict): The EC_MS data dictionary
        name (str): Name of the measurement
        cls (Measurement class): The class to return a measurement of
        reader (Reader object): The class which read ec_ms_dataset from file
        technique (str): The name of the technique
    """

    if "Ewe/V" in ec_ms_dict and "<Ewe>/V" in ec_ms_dict:
        # EC_MS duplicates the latter as the former, so here we delete it:
        del ec_ms_dict["<Ewe>/V"]
    if "I/mA" in ec_ms_dict and "<I>/mA" in ec_ms_dict:
        # EC_MS duplicates the latter as the former, so here we delete it:
        del ec_ms_dict["<I>/mA"]

    cols_str = ec_ms_dict["data_cols"]
    cols_list = []

    name = name or ec_ms_dict.get("title", None)

    for col in cols_str:
        if col.endswith("-x"):
            cols_list.append(TimeSeries(col, "s", ec_ms_dict[col], ec_ms_dict["tstamp"]))

    if "time/s" in ec_ms_dict:
        cols_list.append(
            TimeSeries("time/s", "s", ec_ms_dict["time/s"], ec_ms_dict["tstamp"])
        )

    tseries_meas = Measurement("tseries_ms", technique="EC_MS", series_list=cols_list)

    for col in cols_str:
        if col not in ec_ms_dict or col in tseries_meas.series_names:
            continue
        if col.endswith("-y"):
            v_name = col[:-2]
            tseries = tseries_meas[col[:-1] + "x"]
            unit_name = "A" if col.startswith("M") else ""
        elif col in BIOLOGIC_COLUMN_NAMES and col not in tseries_meas.series_names:
            v_name = col
            tseries = tseries_meas["time/s"]
            unit_name = get_column_unit_name(col)
        else:
            print(f"Not including '{col}' as I don't know what it is.")
            continue
        data = ec_ms_dict[col]
        if not tseries.data.size == data.size:
            print(f"Not including '{col}' due to mismatch size with {tseries}")
            continue
        cols_list.append(
            ValueSeries(
                name=v_name,
                data=data,
                unit_name=unit_name,
                tseries=tseries,
            )
        )

    aliases = {"t": ["time/s"], "raw_potential": ["Ewe/V"], "raw_current": ["I/mA"]}
    obj_as_dict = dict(
        name=name,
        technique=technique or "EC_MS",
        series_list=cols_list,
        reader=reader,
        tstamp=ec_ms_dict["tstamp"],
        aliases=aliases,
    )
    obj_as_dict.update(kwargs)

    measurement = cls.from_dict(obj_as_dict)
    return measurement
