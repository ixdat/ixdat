import numpy as np

from ixdat import Measurement
from ixdat.data_series import DataSeries, TimeSeries, ValueSeries, Field
from ixdat.spectra import Spectrum, SpectrumSeries


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


def test_spectrum_to_and_from_portable_dict_round_trip():
    x = np.linspace(0, 10, 50)
    y = np.sin(x)
    spectrum = Spectrum.from_data(
        x,
        y,
        tstamp=456.0,
        x_name="energy",
        y_name="intensity",
        x_unit_name="eV",
        y_unit_name="counts",
        name="test-spectrum",
        technique="XPS",
        metadata={"sample": "Cu"},
    )

    portable = spectrum.to_portable_dict()
    assert portable["object_type"] == "spectrum"
    assert portable["field"]["series_type"] == "field"
    assert len(portable["field"]["axes_series"]) == 1

    restored = Spectrum.from_portable_dict(portable)
    np.testing.assert_allclose(restored.x, x)
    np.testing.assert_allclose(restored.y, y)
    assert restored.name == "test-spectrum"
    assert restored.technique == "XPS"
    assert restored.metadata["sample"] == "Cu"


def test_spectrum_series_to_and_from_portable_dict_round_trip():
    x = np.linspace(0, 10, 50)
    spectra = [
        Spectrum.from_data(
            x,
            np.sin(x + i),
            tstamp=100.0 + i,
            x_name="energy",
            y_name="intensity",
            name=f"spec-{i}",
            technique="XPS spectrum",
        )
        for i in range(3)
    ]
    ss = SpectrumSeries.from_spectrum_list(spectra)

    portable = ss.to_portable_dict()
    assert portable["object_type"] == "spectrum_series"

    restored = SpectrumSeries.from_portable_dict(portable)
    np.testing.assert_allclose(restored.x, x)
    assert restored.y.shape == (3, 50)
    np.testing.assert_allclose(restored.y[0], np.sin(x))


def test_field_portable_dict_round_trip():
    axis = DataSeries(name="energy", unit_name="eV", data=np.arange(10.0))
    field = Field(
        name="counts",
        unit_name="arb",
        data=np.arange(10.0) * 2,
        axes_series=[axis],
    )

    portable = field.to_portable_dict()
    assert portable["series_type"] == "field"
    assert len(portable["axes_series"]) == 1

    restored = DataSeries.from_portable_dict(portable)
    assert isinstance(restored, Field)
    np.testing.assert_allclose(restored.data, field.data)
    np.testing.assert_allclose(restored.axes_series[0].data, axis.data)
