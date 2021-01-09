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
        cols_list = []

        for col in cols_str:
            if col[0] == "M" and col[-1] == "x":
                cols_list.append(TimeSeries(
                    col, "s",
                    data[col], data["tstamp"]
                    ))

            if col == "time/s":
                cols_list.append(TimeSeries(
                    col, "s",
                    data[col], data["tstamp"]
                    ))

        measurement = Measurement('tseries_ms',
            technique = 'EC_MS',
            series_list=cols_list)


        for col in cols_str:
            if col[0] == "M" and col[-1] == "y":
                cols_list.append(ValueSeries(
                col, "A",
                data[col],
                tseries=measurement[col[:-1] + 'x']
                ))
            if col == 'Ewe/V' or col == 'I/mA':
                cols_list.append(ValueSeries(
                col, "A",
                data[col],
                tseries=measurement['time/s']
                ))

        measurement = Measurement('tseries_ms',
            technique = 'EC_MS',
            series_list=cols_list)

        return measurement
