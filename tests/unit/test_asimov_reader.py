import json

import numpy as np
import pytest
import requests

from ixdat.data_series import TimeSeries, ValueSeries
from ixdat.measurement_base import Measurement
from ixdat.readers.asimov import AsimovConfig, AsimovReader
from ixdat.spectra import Spectrum, SpectrumSeries
from ixdat.techniques.spectroelectrochemistry import ECOpticalMeasurement
from ixdat.techniques.nmr import NMRSpectrum


class FakeSession:
    """Minimal requests.Session stand-in for AsimovReader tests."""

    def __init__(self, response):
        self.response = response
        self.adapters = {}

    def mount(self, prefix, adapter):
        self.adapters[prefix] = adapter

    def request(self, *args, **kwargs):
        return self.response


class FakeQueuedSession:
    """requests.Session stand-in that returns queued responses in order."""

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


def make_json_response(body):
    response = requests.Response()
    response.status_code = 200
    response.url = "https://asimov.example.test/api"
    response._content = json.dumps(body).encode()
    response.headers["content-type"] = "application/json"
    return response


def make_payload_envelope(payload, **overrides):
    envelope = {
        "project_file_id": "dataset-id",
        "filename": "remote.dat",
        "parser_name": "parser",
        "created_at": "2026-06-19T12:00:00Z",
        "payload_json": payload,
    }
    envelope.update(overrides)
    return envelope


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
        reader._get_json(
            "project-files/file-id/ixdat-payload",
            headers={"Authorization": "Bearer t"},
        )

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


def test_spectrum_read_uses_asimov_reader_registration(monkeypatch):
    x = np.linspace(0, 10, 50)
    y = np.sin(x)
    payload = {
        "object_type": "spectrum",
        "name": "remote-spectrum",
        "technique": "XPS",
        "tstamp": 456.0,
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
    session = FakeQueuedSession(
        [
            make_json_response(
                make_payload_envelope(payload, filename="remote-spectrum.dat")
            )
        ]
    )
    monkeypatch.setenv("ASIMOV_ACCESS_TOKEN", "dummy")
    monkeypatch.setattr("ixdat.readers.asimov.requests.Session", lambda: session)

    spectrum = Spectrum.read("dataset-id", reader="asimov")

    assert isinstance(spectrum, Spectrum)
    assert spectrum.name == "remote-spectrum"
    np.testing.assert_allclose(spectrum.x, x)
    np.testing.assert_allclose(spectrum.y, y)
    assert spectrum.metadata["asimov"]["project_file_id"] == "dataset-id"
    assert "dataset_id" not in spectrum.metadata["asimov"]
    assert "version" not in spectrum.metadata["asimov"]
    assert "project-files/dataset-id/ixdat-payload" in session.requests[0][1]["url"]
    assert "dataset-versions" not in session.requests[0][1]["url"]


def test_asimov_reader_falls_back_to_bundle_payload_endpoint(monkeypatch):
    payload = {
        "object_type": "measurement",
        "name": "bundle-measurement",
        "technique": "simple",
        "tstamp": 123.0,
        "series_list": [
            {
                "key": "s0",
                "series_type": "tseries",
                "name": "time/s",
                "unit_name": "s",
                "data": [0.0, 1.0],
                "tstamp": 123.0,
            }
        ],
    }
    session = FakeQueuedSession(
        [
            make_response(404, '{"detail": "Not found"}', reason="Not Found"),
            make_json_response(
                make_payload_envelope(
                    payload,
                    project_file_id=None,
                    bundle_id="bundle-id",
                    name="bundle",
                    filename=None,
                )
            ),
        ]
    )
    monkeypatch.setenv("ASIMOV_ACCESS_TOKEN", "dummy")
    monkeypatch.setattr("ixdat.readers.asimov.requests.Session", lambda: session)

    measurement = Measurement.read("bundle-id", reader="asimov")

    assert isinstance(measurement, Measurement)
    assert measurement.metadata["asimov"]["bundle_id"] == "bundle-id"
    assert "project-files/bundle-id/ixdat-payload" in session.requests[0][1]["url"]
    assert "project-file-bundles/bundle-id/ixdat-payload" in session.requests[1][1]["url"]


def test_spectrum_read_explains_measurement_payload_mismatch(monkeypatch):
    payload = {
        "object_type": "measurement",
        "name": "remote-measurement",
        "technique": "simple",
        "tstamp": 123.0,
        "series_list": [
            {
                "key": "s0",
                "series_type": "tseries",
                "name": "time/s",
                "unit_name": "s",
                "data": [0.0, 1.0],
                "tstamp": 123.0,
            },
            {
                "key": "s1",
                "series_type": "vseries",
                "name": "I/mA",
                "unit_name": "mA",
                "data": [1.0, 2.0],
                "tseries_key": "s0",
            },
        ],
    }
    session = FakeQueuedSession(
        [
            make_json_response(make_payload_envelope(payload))
        ]
    )
    monkeypatch.setenv("ASIMOV_ACCESS_TOKEN", "dummy")
    monkeypatch.setattr("ixdat.readers.asimov.requests.Session", lambda: session)

    with pytest.raises(ValueError) as exception:
        Spectrum.read("dataset-id", reader="asimov")

    message = str(exception.value)
    assert "contains a Measurement payload" in message
    assert "incompatible ixdat object type" in message
    assert "Measurement.read" in message


def test_measurement_read_explains_spectrum_payload_mismatch(monkeypatch):
    x = np.linspace(0, 10, 50)
    y = np.sin(x)
    payload = {
        "object_type": "spectrum",
        "name": "remote-spectrum",
        "technique": "XPS",
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
    session = FakeQueuedSession(
        [
            make_json_response(
                make_payload_envelope(payload, filename="remote-spectrum.dat")
            )
        ]
    )
    monkeypatch.setenv("ASIMOV_ACCESS_TOKEN", "dummy")
    monkeypatch.setattr("ixdat.readers.asimov.requests.Session", lambda: session)

    with pytest.raises(ValueError) as exception:
        Measurement.read("dataset-id", reader="asimov")

    message = str(exception.value)
    assert "contains a Spectrum payload" in message
    assert "incompatible ixdat object type" in message
    assert "Spectrum.read" in message


def test_spectrum_series_read_explains_spectrum_payload_mismatch(monkeypatch):
    x = np.linspace(0, 10, 50)
    y = np.sin(x)
    payload = {
        "object_type": "spectrum",
        "name": "remote-spectrum",
        "technique": "XPS",
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
    session = FakeQueuedSession(
        [
            make_json_response(
                make_payload_envelope(payload, filename="remote-spectrum.dat")
            )
        ]
    )
    monkeypatch.setenv("ASIMOV_ACCESS_TOKEN", "dummy")
    monkeypatch.setattr("ixdat.readers.asimov.requests.Session", lambda: session)

    with pytest.raises(ValueError) as exception:
        SpectrumSeries.read("dataset-id", reader="asimov")

    message = str(exception.value)
    assert "contains a Spectrum payload" in message
    assert "incompatible ixdat object type" in message
    assert "Spectrum.read" in message


def test_asimov_nmr_spectrum_uses_nmr_subclass_and_plotter(monkeypatch):
    x = np.linspace(15, -5, 50)
    y = np.sin(x)
    payload = {
        "object_type": "spectrum",
        "name": "remote-nmr",
        "technique": "NMR",
        "tstamp": 456.0,
        "field": {
            "series_type": "field",
            "name": "intensity",
            "unit_name": None,
            "data": y.tolist(),
            "axes_series": [
                {
                    "series_type": "series",
                    "name": "chemical shift",
                    "unit_name": "ppm",
                    "data": x.tolist(),
                }
            ],
        },
    }
    session = FakeQueuedSession(
        [
            make_json_response(make_payload_envelope(payload, filename="nmr.dat"))
        ]
    )
    monkeypatch.setenv("ASIMOV_ACCESS_TOKEN", "dummy")
    monkeypatch.setattr("ixdat.readers.asimov.requests.Session", lambda: session)

    spectrum = Spectrum.read("dataset-id", reader="asimov")

    assert isinstance(spectrum, NMRSpectrum)
    assert spectrum.plotter.__class__.__name__ == "NMRPlotter"
    np.testing.assert_allclose(spectrum.x, x)
    np.testing.assert_allclose(spectrum.y, y)


def test_asimov_nmr_spectrum_can_be_read_from_nmr_class(monkeypatch):
    x = np.linspace(15, -5, 50)
    y = np.sin(x)
    payload = {
        "object_type": "spectrum",
        "name": "remote-nmr",
        "technique": "NMR",
        "tstamp": 456.0,
        "field": {
            "series_type": "field",
            "name": "intensity",
            "unit_name": None,
            "data": y.tolist(),
            "axes_series": [
                {
                    "series_type": "series",
                    "name": "chemical shift",
                    "unit_name": "ppm",
                    "data": x.tolist(),
                }
            ],
        },
    }
    session = FakeQueuedSession(
        [
            make_json_response(make_payload_envelope(payload, filename="nmr.dat"))
        ]
    )
    monkeypatch.setenv("ASIMOV_ACCESS_TOKEN", "dummy")
    monkeypatch.setattr("ixdat.readers.asimov.requests.Session", lambda: session)

    spectrum = NMRSpectrum.read("dataset-id", reader="asimov")

    assert isinstance(spectrum, NMRSpectrum)
    assert spectrum.plotter.__class__.__name__ == "NMRPlotter"
    np.testing.assert_allclose(spectrum.x, x)
    np.testing.assert_allclose(spectrum.y, y)


def test_spectrum_from_dict_dispatches_by_technique():
    x = np.linspace(15, -5, 50)
    y = np.sin(x)
    spectrum = Spectrum.from_dict(
        {
            "name": "nmr",
            "technique": "NMR",
            "field": _reader()._build_series(
                {
                    "series_type": "field",
                    "name": "intensity",
                    "unit_name": None,
                    "data": y.tolist(),
                    "axes_series": [
                        {
                            "series_type": "series",
                            "name": "chemical shift",
                            "unit_name": "ppm",
                            "data": x.tolist(),
                        }
                    ],
                }
            ),
        }
    )

    assert isinstance(spectrum, NMRSpectrum)
    assert spectrum.plotter.__class__.__name__ == "NMRPlotter"
    np.testing.assert_allclose(spectrum.x, x)


def test_asimov_spectrum_read_can_return_spectrum_series(monkeypatch):
    x = np.linspace(0, 10, 50)
    spectra = np.stack([np.sin(x), np.cos(x)])
    payload = {
        "object_type": "spectrum_series",
        "name": "remote-series",
        "technique": "spectra",
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
                    "data": [0.0, 1.0],
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
    session = FakeQueuedSession(
        [
            make_json_response(make_payload_envelope(payload, filename="series.dat"))
        ]
    )
    monkeypatch.setenv("ASIMOV_ACCESS_TOKEN", "dummy")
    monkeypatch.setattr("ixdat.readers.asimov.requests.Session", lambda: session)

    spectrum_series = Spectrum.read("dataset-id", reader="asimov")

    assert isinstance(spectrum_series, SpectrumSeries)
    np.testing.assert_allclose(spectrum_series.x, x)
    np.testing.assert_allclose(spectrum_series.y, spectra)


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


def test_build_kwargs_hydrates_generic_nested_ixdat_objects():
    x = np.linspace(0, 10, 5)
    y = np.sin(x)
    payload = {
        "object_type": "measurement",
        "name": "demo",
        "technique": "simple",
        "tstamp": 1.0,
        "series_list": [
            {
                "key": "s0",
                "series_type": "tseries",
                "name": "time/s",
                "unit_name": "s",
                "data": [0.0, 1.0],
                "tstamp": 1.0,
            }
        ],
        "calibration_spectrum": {
            "object_type": "spectrum",
            "name": "calibration",
            "technique": "spectrum",
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
        },
    }

    obj = _reader()._build_kwargs(payload)

    assert isinstance(obj["calibration_spectrum"], Spectrum)
    np.testing.assert_allclose(obj["calibration_spectrum"].x, x)
    np.testing.assert_allclose(obj["calibration_spectrum"].y, y)


def test_build_kwargs_hydrates_generic_nested_ixdat_object_lists():
    payload = {
        "object_type": "measurement",
        "name": "combined",
        "technique": "simple",
        "tstamp": 1.0,
        "series_list": [
            {
                "key": "s0",
                "series_type": "tseries",
                "name": "time/s",
                "unit_name": "s",
                "data": [0.0, 1.0],
                "tstamp": 1.0,
            }
        ],
        "component_measurements": [
            {
                "object_type": "measurement",
                "name": "part",
                "technique": "simple",
                "tstamp": 1.0,
                "series_list": [
                    {
                        "key": "s0",
                        "series_type": "tseries",
                        "name": "time/s",
                        "unit_name": "s",
                        "data": [0.0],
                        "tstamp": 1.0,
                    }
                ],
            }
        ],
    }

    obj = _reader()._build_kwargs(payload)

    assert len(obj["component_measurements"]) == 1
    assert isinstance(obj["component_measurements"][0], Measurement)
    assert obj["component_measurements"][0].name == "part"


def test_build_object_supports_ixdat_container_object_types():
    reader = _reader()
    measurement = reader._build_object(
        {
            "object_type": "measurement",
            "name": "demo-measurement",
            "technique": "simple",
            "tstamp": 1.0,
            "series_list": [
                {
                    "key": "s0",
                    "series_type": "tseries",
                    "name": "time/s",
                    "unit_name": "s",
                    "data": [0.0, 1.0],
                    "tstamp": 1.0,
                }
            ],
        }
    )
    assert isinstance(measurement, Measurement)

    spectrum = reader._build_object(
        {
            "object_type": "spectrum",
            "name": "demo-spectrum",
            "technique": "spectrum",
            "field": {
                "series_type": "field",
                "name": "intensity",
                "unit_name": "counts",
                "data": [1.0, 2.0],
                "axes_series": [
                    {
                        "series_type": "series",
                        "name": "energy",
                        "unit_name": "eV",
                        "data": [0.0, 1.0],
                    }
                ],
            },
        }
    )
    assert isinstance(spectrum, Spectrum)

    spectrum_series = reader._build_object(
        {
            "object_type": "spectrum_series",
            "name": "demo-series",
            "technique": "spectra",
            "field": {
                "series_type": "field",
                "name": "intensity",
                "unit_name": "counts",
                "data": [[1.0, 2.0], [3.0, 4.0]],
                "axes_series": [
                    {
                        "series_type": "tseries",
                        "name": "Spectrum Time",
                        "unit_name": "s",
                        "data": [0.0, 1.0],
                        "tstamp": 1.0,
                    },
                    {
                        "series_type": "series",
                        "name": "energy",
                        "unit_name": "eV",
                        "data": [0.0, 1.0],
                    },
                ],
            },
        },
        cls=Spectrum,
    )
    assert isinstance(spectrum_series, SpectrumSeries)


def test_build_object_rejects_unknown_container_object_type():
    with pytest.raises(ValueError) as exception:
        _reader()._build_object(
            {
                "object_type": "future_container",
                "name": "demo",
                "technique": "simple",
            }
        )

    assert "unsupported object_type='future_container'" in str(exception.value)
    assert "measurement" in str(exception.value)
    assert "spectrum" in str(exception.value)
    assert "spectrum_series" in str(exception.value)


def test_parse_spectro_measurement_with_nested_spectrum_objects():
    x = np.linspace(400, 700, 4)
    spectra = np.stack([np.sin(x), np.cos(x)])
    payload = {
        "object_type": "measurement",
        "name": "sec-demo",
        "technique": "EC-Optical",
        "tstamp": 100.0,
        "series_list": [
            {
                "key": "s0",
                "series_type": "tseries",
                "name": "time/s",
                "unit_name": "s",
                "data": [0.0, 1.0],
                "tstamp": 100.0,
            },
            {
                "key": "s1",
                "series_type": "vseries",
                "name": "Ewe/V",
                "unit_name": "V",
                "data": [0.1, 0.2],
                "tseries_key": "s0",
            },
        ],
        "spectrum_series": {
            "object_type": "spectrum_series",
            "name": "uv-vis",
            "technique": "spectra",
            "tstamp": 100.0,
            "metadata": {"spectrometer": "demo"},
            "field": {
                "series_type": "field",
                "name": "absorbance",
                "unit_name": "a.u.",
                "data": spectra.tolist(),
                "axes_series": [
                    {
                        "series_type": "tseries",
                        "name": "Spectrum Time",
                        "unit_name": "s",
                        "data": [0.0, 1.0],
                        "tstamp": 100.0,
                    },
                    {
                        "series_type": "series",
                        "name": "wavelength",
                        "unit_name": "nm",
                        "data": x.tolist(),
                    },
                ],
            },
            "durations": [0.5, 0.5],
            "continuous": True,
        },
        "reference_spectrum": {
            "object_type": "spectrum",
            "name": "reference",
            "technique": "spectrum",
            "tstamp": 99.0,
            "metadata": {"dark_corrected": True},
            "field": {
                "series_type": "field",
                "name": "reference intensity",
                "unit_name": "counts",
                "data": [1.0, 2.0, 3.0, 4.0],
                "axes_series": [
                    {
                        "series_type": "series",
                        "name": "wavelength",
                        "unit_name": "nm",
                        "data": x.tolist(),
                    }
                ],
            },
        },
    }

    measurement = Measurement.from_dict(_reader()._build_kwargs(payload))

    assert isinstance(measurement.spectrum_series, SpectrumSeries)
    assert measurement.spectrum_series.metadata == {"spectrometer": "demo"}
    assert measurement.spectrum_series.continuous is True
    np.testing.assert_allclose(measurement.spectrum_series.durations, [0.5, 0.5])
    np.testing.assert_allclose(measurement.spectrum_series.x, x)
    np.testing.assert_allclose(measurement.spectrum_series.y, spectra)

    assert isinstance(measurement.reference_spectrum, Spectrum)
    assert measurement.reference_spectrum.metadata == {"dark_corrected": True}
    np.testing.assert_allclose(measurement.reference_spectrum.x, x)
    np.testing.assert_allclose(measurement.reference_spectrum.y, [1.0, 2.0, 3.0, 4.0])


def test_nested_ixdat_object_without_object_type_fails_clearly():
    payload = {
        "object_type": "measurement",
        "name": "sec-demo",
        "technique": "EC-Optical",
        "tstamp": 100.0,
        "series_list": [
            {
                "key": "s0",
                "series_type": "tseries",
                "name": "time/s",
                "unit_name": "s",
                "data": [0.0, 1.0],
                "tstamp": 100.0,
            }
        ],
        "reference_spectrum": {
            "name": "reference",
            "technique": "spectrum",
            "field": {
                "series_type": "field",
                "name": "reference intensity",
                "unit_name": "counts",
                "data": [1.0, 2.0],
                "axes_series": [
                    {
                        "series_type": "series",
                        "name": "wavelength",
                        "unit_name": "nm",
                        "data": [400.0, 500.0],
                    }
                ],
            },
        },
    }

    with pytest.raises(ValueError) as exception:
        _reader()._build_kwargs(payload)

    message = str(exception.value)
    assert "reference_spectrum" in message
    assert "missing 'object_type'" in message


def test_parse_ec_optical_measurement_preserves_linker_ids():
    payload = {
        "object_type": "measurement",
        "name": "sec-demo",
        "technique": "EC-Optical",
        "tstamp": 100.0,
        "series_list": [
            {
                "key": "s0",
                "series_type": "tseries",
                "name": "time/s",
                "unit_name": "s",
                "data": [0.0, 1.0],
                "tstamp": 100.0,
            },
            {
                "key": "s1",
                "series_type": "vseries",
                "name": "Ewe/V",
                "unit_name": "V",
                "data": [0.1, 0.2],
                "tseries_key": "s0",
            },
        ],
        "spectrum_id": 12,
        "ref_id": 13,
    }

    measurement = _reader()._build_object(payload)

    assert isinstance(measurement, ECOpticalMeasurement)
    assert measurement._spectrum_series.id == 12
    assert measurement._reference_spectrum.id == 13


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

    # At the Measurement level, series with the same name are concatenated and
    # aligned to the measurement's tstamp (1.0). The second segment has
    # tstamp=101.0, so its data is offset by (101.0 - 1.0) = 100 s.
    meas = Measurement.from_dict(obj)
    np.testing.assert_allclose(meas["time/s"].data, [0.0, 1.0, 200.0, 201.0])


def test_series_payload_ignores_non_ixdat_constructor_metadata():
    series = _reader()._build_series(
        {
            "series_type": "tseries",
            "class": "TimeSeries",
            "name": "time/s",
            "unit_name": "s",
            "data": [0.0, 1.0],
            "tstamp": 1.0,
        }
    )

    assert isinstance(series, TimeSeries)
    np.testing.assert_allclose(series.data, [0.0, 1.0])


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
        s for s in obj["series_list"] if not isinstance(s, (TimeSeries, ValueSeries))
    )
    assert field.axes_series[0] is top_level_t
