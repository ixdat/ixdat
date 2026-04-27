import numpy as np
import pytest
import requests

from ixdat.data_series import TimeSeries, ValueSeries
from ixdat.measurement_base import Measurement
from ixdat.readers.asimov import AsimovConfig, AsimovReader
from ixdat.spectra import Spectrum, SpectrumSeries


class FakeSession:
    """Minimal requests.Session stand-in for AsimovReader tests."""

    def __init__(self, response):
        self.response = response
        self.adapters = {}

    def mount(self, prefix, adapter):
        self.adapters[prefix] = adapter

    def request(self, *args, **kwargs):
        return self.response


def make_response(status_code, body, reason=""):
    response = requests.Response()
    response.status_code = status_code
    response.reason = reason
    response.url = "https://asimov.example.test/api/datasets/dataset-id"
    response._content = body.encode()
    response.headers["content-type"] = "application/json"
    return response


def test_asimov_reader_explains_unauthorized_client():
    response = make_response(
        401, '{"detail": "Unauthorized client"}', reason="Unauthorized"
    )
    reader = AsimovReader(
        token="token",
        config=AsimovConfig(base_url="https://asimov.example.test/api"),
    )
    reader._session = FakeSession(response)

    with pytest.raises(RuntimeError) as exception:
        reader._get_json("datasets/dataset-id", headers={"Authorization": "Bearer t"})

    message = str(exception.value)
    assert "Unauthorized client" in message
    assert "rejected the OAuth client" in message
    assert "KEYCLOAK_CLIENT_ID" in message
    assert "ASIMOV_ACCESS_TOKEN" in message


# ----------------------------------------------------------------------
# Payload-parsing tests. These build Asimov-style JSON payloads inline
# and feed them to the reader's private parse helpers, so they don't
# touch the network.
# ----------------------------------------------------------------------


def _reader():
    # A dummy token short-circuits Keycloak setup in __init__.
    return AsimovReader(token="dummy")


def test_parse_measurement_payload():
    payload = {
        "object_type": "measurement",
        "name": "demo",
        "technique": "simple",
        "tstamp": 123.0,
        "metadata": {"numbers": [1, 2, 3]},
        "aliases": {"current": ["I/mA"]},
        "series_list": [
            {
                "key": "s0",
                "series_type": "tseries",
                "name": "time/s",
                "unit_name": "s",
                "data": [0.0, 1.0, 2.0],
                "tstamp": 123.0,
            },
            {
                "key": "s1",
                "series_type": "vseries",
                "name": "I/mA",
                "unit_name": "mA",
                "data": [1.0, 2.0, 3.0],
                "tseries_key": "s0",
            },
        ],
    }

    measurement = Measurement.from_dict(_reader()._build_kwargs(payload))

    assert isinstance(measurement, Measurement)
    assert measurement.name == "demo"
    assert measurement.aliases["current"] == ["I/mA"]

    tseries = next(s for s in measurement.series_list if isinstance(s, TimeSeries))
    vseries = next(s for s in measurement.series_list if isinstance(s, ValueSeries))
    np.testing.assert_allclose(tseries.data, [0.0, 1.0, 2.0])
    np.testing.assert_allclose(vseries.data, [1.0, 2.0, 3.0])
    # The ValueSeries should be linked to the TimeSeries from the same payload.
    assert vseries.tseries is tseries


def test_parse_spectrum_payload():
    x = np.linspace(0, 10, 50)
    y = np.sin(x)
    payload = {
        "object_type": "spectrum",
        "name": "test-spectrum",
        "technique": "XPS",
        "tstamp": 456.0,
        "metadata": {"sample": "Cu"},
        "field": {
            "series_type": "field",
            "name": "intensity",
            "unit_name": "counts",
            "data": y.tolist(),
            "axes_series": [
                {
                    "series_type": "series",
                    "name": "energy",
                    "unit_name": "eV",
                    "data": x.tolist(),
                }
            ],
        },
    }

    spectrum = Spectrum.from_dict(_reader()._build_kwargs(payload))

    assert isinstance(spectrum, Spectrum)
    assert spectrum.name == "test-spectrum"
    assert spectrum.technique == "XPS"
    assert spectrum.metadata["sample"] == "Cu"
    np.testing.assert_allclose(spectrum.x, x)
    np.testing.assert_allclose(spectrum.y, y)


def test_parse_spectrum_series_payload():
    x = np.linspace(0, 10, 50)
    spectra = np.stack([np.sin(x + i) for i in range(3)])
    payload = {
        "object_type": "spectrum_series",
        "name": "demo-series",
        "technique": "XPS spectra",
        "tstamp": 100.0,
        "field": {
            "series_type": "field",
            "name": "intensity",
            "unit_name": "counts",
            "data": spectra.tolist(),
            "axes_series": [
                {
                    "series_type": "tseries",
                    "name": "Spectrum Time",
                    "unit_name": "s",
                    "data": [0.0, 1.0, 2.0],
                    "tstamp": 100.0,
                },
                {
                    "series_type": "series",
                    "name": "energy",
                    "unit_name": "eV",
                    "data": x.tolist(),
                },
            ],
        },
    }

    series = SpectrumSeries.from_dict(_reader()._build_kwargs(payload))

    assert isinstance(series, SpectrumSeries)
    np.testing.assert_allclose(series.x, x)
    assert series.y.shape == (3, 50)
    np.testing.assert_allclose(series.y[0], np.sin(x))


def test_appended_measurement_same_named_tseries_resolved_by_key():
    """Two TimeSeries with the same name, common in appended measurements,
    are disambiguated via per-entry payload keys rather than name strings."""
    payload = {
        "object_type": "measurement",
        "name": "appended",
        "tstamp": 1.0,
        "series_list": [
            {
                "key": "s0",
                "series_type": "tseries",
                "name": "time/s",
                "unit_name": "s",
                "data": [0.0, 1.0],
                "tstamp": 1.0,
            },
            {
                "key": "s1",
                "series_type": "vseries",
                "name": "I",
                "unit_name": "A",
                "data": [10.0, 11.0],
                "tseries_key": "s0",
            },
            {
                "key": "s2",
                "series_type": "tseries",
                "name": "time/s",
                "unit_name": "s",
                "data": [100.0, 101.0],
                "tstamp": 101.0,
            },
            {
                "key": "s3",
                "series_type": "vseries",
                "name": "I",
                "unit_name": "A",
                "data": [20.0, 21.0],
                "tseries_key": "s2",
            },
        ],
    }

    obj = _reader()._build_kwargs(payload)
    tseriess = [s for s in obj["series_list"] if isinstance(s, TimeSeries)]
    vseriess = [s for s in obj["series_list"] if isinstance(s, ValueSeries)]
    t0, t1 = tseriess
    v0, v1 = vseriess
    assert t0 is not t1  # distinct objects despite shared name
    np.testing.assert_allclose(t0.data, [0.0, 1.0])
    np.testing.assert_allclose(t1.data, [100.0, 101.0])
    assert v0.tseries is t0
    assert v1.tseries is t1


def test_field_axis_reuses_top_level_tseries_via_axes_keys():
    """A Field whose axis is a top-level series references it via
    axes_keys, so reader and writer agree on a single shared object."""
    payload = {
        "object_type": "measurement",
        "name": "demo",
        "tstamp": 1.0,
        "series_list": [
            {
                "key": "s0",
                "series_type": "tseries",
                "name": "time/s",
                "unit_name": "s",
                "data": [0.0, 1.0, 2.0],
                "tstamp": 1.0,
            },
            {
                "key": "s1",
                "series_type": "field",
                "name": "spectrum",
                "unit_name": "counts",
                "data": [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                "axes_keys": ["s0"],
            },
        ],
    }

    obj = _reader()._build_kwargs(payload)
    top_level_t = next(s for s in obj["series_list"] if isinstance(s, TimeSeries))
    field = next(
        s for s in obj["series_list"]
        if not isinstance(s, (TimeSeries, ValueSeries))
    )
    assert field.axes_series[0] is top_level_t
