from . import TECHNIQUE_CLASSES
import pickle
from ..data_series import TimeSeries, ValueSeries
from ..measurements import Measurement
from ..techniques.ec_ms import ECMSMeasurement


ECMSMeasruement = TECHNIQUE_CLASSES["EC-MS"]


class EC_MS_CONVERTER:
    """Imports old .pkl files obtained from the legacy EC-MS package"""

    def __init__(self):
        print("Reader of old ECMS .pkl files")

    def read(self, file_path):
        """Return an ECMSMeasurement with the data recorded in path_to_file

        This loops through the keys of the EC-MS dict and searches for MS and
        EC data. Names the dataseries according to their names in the original
        dict. Omitts any other data as well as metadata.

        Args:
            path_to_file (Path): The full abs or rel path including the
            ".pkl" extension.
        """
        with open(file_path, "rb") as f:
            data = pickle.load(f)

        cols_str = data["data_cols"]
        cols_list = []

        for col in cols_str:
            if col.endswith("-x"):
                cols_list.append(TimeSeries(col, "s", data[col], data["tstamp"]))

            if col == "time/s":
                cols_list.append(TimeSeries(col, "s", data[col], data["tstamp"]))

        measurement = Measurement(
            "tseries_ms", technique="EC_MS", series_list=cols_list
        )

        for col in cols_str:
            if col.startswith("M") and col.endswith("-y"):
                cols_list.append(
                    ValueSeries(
                        col[:-2], "A", data[col], tseries=measurement[col[:-1] + "x"]
                    )
                )

            # TODO: Import all EC data.
            if col == "Ewe/V" or col == "I/mA":
                cols_list.append(
                    ValueSeries(col, "A", data[col], tseries=measurement["time/s"])
                )

        measurement = ECMSMeasurement(
            file_path,
            technique="EC_MS",
            series_list=cols_list,
            reader=self,
            tstamp=data["tstamp"],
        )

        return measurement
