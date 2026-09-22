"""Reader tests using small payload examples, without contacting Asimov."""

import json

import numpy as np
import pytest
import requests

from ixdat.data_series import TimeSeries, ValueSeries
from ixdat.measurement_base import Measurement
from ixdat.readers.asimov import AsimovConfig, AsimovReader
from ixdat.spectra import Spectrum, SpectrumSeries
from ixdat.techniques.nmr import NMRSpectrum


class FakeSession:
    """Return queued responses and record the requested URLs."""

    def __init__(self, responses):
        self.responses = list(responses)
        self.requests = []
        self.adapters = {}
        self.trust_env = False

    def mount(self, prefix, adapter):
        self.adapters[prefix] = adapter

    def request(self, *args, **kwargs):
        self.requests.append((args, kwargs))
        return self.responses.pop(0)


def make_response(status_code, body, reason=""):
    response = requests.Response()
    response.status_code = status_code
    response.reason = reason
    response.url = "https://asimov.example.test/api/project-files/file-id/ixdat-payload"
    response._content = body.encode()
    response.headers["content-type"] = "application/json"
    return response


def make_payload_response(payload, **envelope_fields):
    envelope = {
        "project_file_id": "file-id",
        "filename": "example.dat",
        "payload_json": payload,
        **envelope_fields,
    }
    return make_response(200, json.dumps(envelope))


def mock_session(monkeypatch, *responses):
    session = FakeSession(responses)
    monkeypatch.setenv("ASIMOV_ACCESS_TOKEN", "dummy")
    monkeypatch.setattr("ixdat.readers.asimov.requests.Session", lambda: session)
    return session


def reader():
    return AsimovReader(token="dummy")


def time_payload(key="time"):
    return {
        "key": key,
        "series_type": "tseries",
        "name": "time/s",
        "unit_name": "s",
        "data": [0.0, 1.0],
        "tstamp": 100.0,
    }


def measurement_payload():
    return {
        "object_type": "measurement",
        "name": "example measurement",
        "technique": "simple",
        "tstamp": 100.0,
        "series_list": [
            time_payload(),
            {
                "key": "current",
                "series_type": "vseries",
                "name": "I/mA",
                "unit_name": "mA",
                "data": [1.0, 2.0],
                "tseries_key": "time",
            },
        ],
    }


def spectrum_payload():
    return {
        "object_type": "spectrum",
        "name": "example NMR",
        "technique": "NMR",
        "field": {
            "series_type": "field",
            "name": "intensity",
            "unit_name": None,
            "data": [1.0, 2.0],
            "axes_series": [
                {
                    "series_type": "series",
                    "name": "chemical shift",
                    "unit_name": "ppm",
                    "data": [1.0, 0.0],
                }
            ],
        },
    }


def spectrum_series_payload():
    payload = spectrum_payload()
    payload.update(object_type="spectrum_series", technique="NMR_spectra", tstamp=100.0)
    payload["field"]["data"] = [[1.0, 2.0], [3.0, 4.0]]
    payload["field"]["axes_series"].insert(0, time_payload())
    return payload


def test_asimov_reader_explains_unauthorized_client():
    reader_obj = AsimovReader(
        token="token",
        config=AsimovConfig(base_url="https://asimov.example.test/api"),
    )
    reader_obj._session = FakeSession(
        [make_response(401, '{"detail": "Unauthorized client"}', "Unauthorized")]
    )

    with pytest.raises(RuntimeError, match="Unauthorized client") as exception:
        reader_obj._get_json(
            "project-files/file-id/ixdat-payload",
            headers={"Authorization": "Bearer t"},
        )

    assert "KEYCLOAK_CLIENT_ID" in str(exception.value)


def test_spectrum_read_uses_file_endpoint(monkeypatch):
    session = mock_session(monkeypatch, make_payload_response(spectrum_payload()))

    spectrum = Spectrum.read("file-id", reader="asimov")

    assert isinstance(spectrum, NMRSpectrum)
    assert spectrum.plotter.__class__.__name__ == "NMRPlotter"
    np.testing.assert_allclose(spectrum.x, [1.0, 0.0])
    np.testing.assert_allclose(spectrum.y, [1.0, 2.0])
    assert spectrum.metadata["asimov"]["project_file_id"] == "file-id"
    assert "project-files/file-id/ixdat-payload" in session.requests[0][1]["url"]


def test_measurement_read_falls_back_to_bundle_endpoint(monkeypatch):
    session = mock_session(
        monkeypatch,
        make_response(404, '{"detail": "Not found"}', "Not Found"),
        make_payload_response(
            measurement_payload(), project_file_id=None, bundle_id="bundle-id"
        ),
    )

    measurement = Measurement.read("bundle-id", reader="asimov")

    tseries, vseries = measurement.series_list
    assert isinstance(tseries, TimeSeries)
    assert isinstance(vseries, ValueSeries)
    assert vseries.tseries is tseries
    np.testing.assert_allclose(vseries.data, [1.0, 2.0])
    assert measurement.metadata["asimov"]["bundle_id"] == "bundle-id"
    assert (
        "project-file-bundles/bundle-id/ixdat-payload" in session.requests[1][1]["url"]
    )


@pytest.mark.parametrize(
    "requested_class,payload_type,recommended_read",
    [
        (Spectrum, "measurement", "Measurement.read"),
        (Measurement, "spectrum", "Spectrum.read"),
        (SpectrumSeries, "spectrum", "Spectrum.read"),
        (NMRSpectrum, "spectrum_series", "SpectrumSeries.read"),
    ],
)
def test_wrong_read_entrypoint_explains_expected_class(
    monkeypatch, requested_class, payload_type, recommended_read
):
    mock_session(monkeypatch, make_payload_response({"object_type": payload_type}))

    with pytest.raises(ValueError, match="incompatible ixdat object type") as exception:
        requested_class.read("file-id", reader="asimov")

    assert recommended_read in str(exception.value)


def test_spectrum_read_accepts_spectrum_series(monkeypatch):
    mock_session(monkeypatch, make_payload_response(spectrum_series_payload()))

    series = Spectrum.read("file-id", reader="asimov")

    assert isinstance(series, SpectrumSeries)
    np.testing.assert_allclose(series.y, [[1.0, 2.0], [3.0, 4.0]])


def test_nested_object_is_rebuilt_by_object_type():
    payload = measurement_payload()
    payload["reference_spectrum"] = spectrum_payload()

    kwargs = reader()._build_kwargs(payload)

    assert isinstance(kwargs["reference_spectrum"], NMRSpectrum)


def test_nested_object_without_object_type_fails():
    payload = measurement_payload()
    payload["reference_spectrum"] = {"field": {}}

    with pytest.raises(ValueError, match="missing 'object_type'"):
        reader()._build_kwargs(payload)


def test_extra_series_metadata_does_not_reach_constructor():
    payload = time_payload()
    payload["class"] = "TimeSeries"

    series = reader()._build_series(payload)

    assert isinstance(series, TimeSeries)


def test_appended_measurement_series_links_use_keys():
    payload = measurement_payload()
    second_time = time_payload(key="second_time")
    second_time["name"] = "time/s"
    second_time["tstamp"] = 200.0
    second_current = {
        "key": "second_current",
        "series_type": "vseries",
        "name": "I/mA",
        "unit_name": "mA",
        "data": [3.0, 4.0],
        "tseries_key": "second_time",
    }
    payload["series_list"].extend([second_time, second_current])

    series = reader()._build_kwargs(payload)["series_list"]
    tseries = [item for item in series if isinstance(item, TimeSeries)]
    vseries = [item for item in series if isinstance(item, ValueSeries)]

    assert vseries[0].tseries is tseries[0]
    assert vseries[1].tseries is tseries[1]


def test_field_axis_reuses_top_level_time_series():
    payload = measurement_payload()
    payload["series_list"].append(
        {
            "key": "field",
            "series_type": "field",
            "name": "intensity",
            "unit_name": "counts",
            "data": [1.0, 2.0],
            "axes_keys": ["time"],
        }
    )

    series = reader()._build_kwargs(payload)["series_list"]

    assert series[2].axes_series[0] is series[0]
