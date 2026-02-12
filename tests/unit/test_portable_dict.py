import numpy as np

from ixdat import Measurement
from ixdat.data_series import TimeSeries, ValueSeries


def test_measurement_to_and_from_portable_dict_round_trip():
    time_series = TimeSeries(
        name="time/s",
        unit_name="s",
        data=np.array([0.0, 1.0, 2.0]),
        tstamp=123.0,
    )
    value_series = ValueSeries(
        name="I/mA",
        unit_name="mA",
        data=np.array([1.0, 2.0, 3.0]),
        tseries=time_series,
    )
    measurement = Measurement(
        name="portable-demo",
        technique="simple",
        tstamp=123.0,
        series_list=[time_series, value_series],
        metadata={"numbers": np.array([1, 2, 3])},
        aliases={"current": ["I/mA"]},
    )

    portable = measurement.to_portable_dict()
    assert portable["metadata"]["numbers"] == [1, 2, 3]
    assert portable["series_list"][1]["tseries_name"] == "time/s"

    restored = Measurement.from_portable_dict(portable)
    restored_time = next(s for s in restored.series_list if isinstance(s, TimeSeries))
    restored_value = next(s for s in restored.series_list if isinstance(s, ValueSeries))

    np.testing.assert_allclose(restored_time.data, [0.0, 1.0, 2.0])
    np.testing.assert_allclose(restored_value.data, [1.0, 2.0, 3.0])
    assert restored_value.tseries is restored_time
    assert restored.aliases["current"] == ["I/mA"]
