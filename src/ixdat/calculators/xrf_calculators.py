from ..exceptions import QuantificationError
from ..measurement_base import Calculator
from ..data_series import ValueSeries


class TRXRFCalculator(Calculator):
    available_series_names = {"FF_over_I0"}

    def calculate_series(self, key, measurement=None):
        """Construct and return a new series of FF/I0"""
        if key != "FF_over_I0":
            raise QuantificationError(f"{self} cannot calculate {key}")
        return ValueSeries(
            name="FF_over_I0",
            unit_name="",
            data=measurement["FF"].data / measurement["I0"].data,
            tseries=measurement["FF"].tseries,
        )
