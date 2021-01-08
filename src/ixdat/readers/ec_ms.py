from . import TECHNIQUE_CLASSES
import pickle
from .data_series import TimeSeries, ValueSeries
from .measurements import Measurement

ECMSMeasruement = TECHNIQUE_CLASSES["EC-MS"]


class EC_MS_CONVERTER:
    def __init__(self):
        print('Reader of old ECMS .pkl files')

    def read(self, file_path):

        with open(file_path, "rb") as f:
            data = pickle.load(f)

        cols_str = data['data_cols']
        data = Measurement("test",
            technique="EC_MS",
            )
        for col in cols:

            if col[0] == "M" and col[-1] == "x":
                cols_list.append(TimeSeries(col, "s", data[col], data["tstamp"]))

            if col == "time/s":
                cols_list.append(TimeSeries(col, "s", data[col], data["tstamp"]))

        for col in cols:

            if col[0] == "M" and col[-1] == "y":
                cols_list.append(ValueSeries(col, "A", data[col], data["tstamp"], tseries=))

            if col == "I/mA"
                cols_list.append(ValueSeries(col, "mA", data[col], data["tstamp"]))

            if col == "Ewe/V"
                cols_list.append(ValueSeries(col, "V", data[col], data["tstamp"]))


        pass
