from ixdat import Measurement
from ixdat.readers.asimov import AsimovReader


def _build_reader_with_stubbed_api(monkeypatch, dataset_id, dataset, versions):
    reader = AsimovReader(base_url="https://asimov.enci.dk/api", token="test-token")

    def fake_get_json(endpoint, headers, params=None):
        assert headers["Authorization"] == "Bearer test-token"
        if endpoint == f"datasets/{dataset_id}":
            return dataset
        if endpoint == "dataset-versions":
            assert params == {"dataset_id": dataset_id}
            return versions
        raise AssertionError(f"Unexpected endpoint: {endpoint}")

    monkeypatch.setattr(reader, "_get_json", fake_get_json)
    return reader


def test_asimov_reader_reads_latest_series_payload(monkeypatch):
    dataset_id = "a47ac10b-58cc-4372-a567-0e02b2c3d479"
    dataset = {
        "id": dataset_id,
        "measurement_run_id": "550e8400-e29b-41d4-a716-446655440000",
        "kind": "eis",
        "label": "EIS at 0.5V",
        "created_at": "2024-02-11T16:00:00Z",
    }
    older_payload = {
        "series": {
            "frequency": [1.0, 0.5],
            "Z_real": [7.0, 8.0],
            "Z_imag": [-0.7, -0.8],
        },
        "metadata": {"electrode_area_cm2": 0.5},
    }
    latest_payload = {
        "series": {
            "frequency": [10.0, 1.0, 0.1],
            "Z_real": [50.2, 250.3, 450.8],
            "Z_imag": [-5.1, -210.5, -350.2],
        },
        "metadata": {"electrode_area_cm2": 0.5},
    }
    versions = [
        {
            "id": "version-old",
            "dataset_id": dataset_id,
            "version": 1,
            "created_at": "2024-02-11T16:15:30Z",
            "summary": {"technique": "EIS"},
            "payload_json": older_payload,
        },
        {
            "id": "version-latest",
            "dataset_id": dataset_id,
            "version": 2,
            "created_at": "2024-02-12T16:15:30Z",
            "summary": {"technique": "EIS"},
            "payload_json": latest_payload,
        },
    ]

    reader = _build_reader_with_stubbed_api(monkeypatch, dataset_id, dataset, versions)
    measurement = reader.read(dataset_id, cls=Measurement)

    assert measurement.name == "EIS at 0.5V"
    assert "Z_real" in measurement.series_names
    assert measurement["Z_real"].data.tolist() == [50.2, 250.3, 450.8]
    assert measurement.metadata["asimov"]["dataset_id"] == dataset_id
    assert measurement.metadata["asimov"]["dataset_version"] == 2


def test_asimov_reader_can_select_explicit_version(monkeypatch):
    dataset_id = "d47ac10b-58cc-4372-a567-0e02b2c3d123"
    dataset = {
        "id": dataset_id,
        "kind": "eis",
        "label": "EIS at OCV",
    }
    versions = [
        {
            "id": "v1-id",
            "dataset_id": dataset_id,
            "version": 1,
            "created_at": "2024-02-11T16:15:30Z",
            "summary": {"technique": "EIS"},
            "payload_json": {
                "series": {
                    "frequency": [1000.0, 100.0],
                    "Z_real": [10.0, 20.0],
                },
                "metadata": {},
            },
        },
        {
            "id": "v2-id",
            "dataset_id": dataset_id,
            "version": 2,
            "created_at": "2024-02-12T16:15:30Z",
            "summary": {"technique": "EIS"},
            "payload_json": {
                "series": {
                    "frequency": [1000.0, 100.0],
                    "Z_real": [30.0, 40.0],
                },
                "metadata": {},
            },
        },
    ]

    reader = _build_reader_with_stubbed_api(monkeypatch, dataset_id, dataset, versions)
    measurement = reader.read(dataset_id, cls=Measurement, version=1)

    assert measurement["Z_real"].data.tolist() == [10.0, 20.0]
    assert measurement.metadata["asimov"]["dataset_version"] == 1


def test_asimov_reader_default_keycloak_settings(monkeypatch):
    monkeypatch.delenv("ASIMOV_BASE_URL", raising=False)
    monkeypatch.delenv("ASIMOV_ACCESS_TOKEN", raising=False)
    monkeypatch.delenv("KEYCLOAK_SERVER_URL", raising=False)
    monkeypatch.delenv("KEYCLOAK_REALM", raising=False)
    monkeypatch.delenv("KEYCLOAK_CLIENT_ID", raising=False)
    monkeypatch.delenv("KEYCLOAK_CLIENT_SECRET", raising=False)

    reader = AsimovReader()
    provider = reader.token_provider
    assert provider is not None
    assert provider.server_url == "https://auth.enci.dk"
    assert provider.realm == "master"
    assert provider.client_id == "ixdat-cli"


def test_asimov_reader_adds_ec_aliases_to_pass_essential_series_check(monkeypatch):
    dataset_id = "82b539f0-6973-4ba3-9cec-b857595a6c8e"
    dataset = {"id": dataset_id, "kind": "ec", "label": "Parsed: 13_PdAg_C01.mpt"}
    versions = [
        {
            "id": "v1-id",
            "dataset_id": dataset_id,
            "version": 1,
            "created_at": "2024-02-12T16:15:30Z",
            "summary": {"technique": "EC"},
            "payload_json": {
                "series": {
                    "time/s": [0.0, 1.0, 2.0],
                    "Ewe/V": [0.1, 0.2, 0.3],
                    "I/mA": [1.0, 1.1, 1.2],
                },
                "metadata": {},
            },
        }
    ]

    reader = _build_reader_with_stubbed_api(monkeypatch, dataset_id, dataset, versions)
    measurement = Measurement.read(dataset_id, reader=reader)

    assert measurement.__class__.__name__ == "ECMeasurement"
    assert measurement.t[0] == 0.0
    assert measurement["raw_potential"].data.tolist() == [0.1, 0.2, 0.3]
    assert measurement["raw_current"].data.tolist() == [1.0, 1.1, 1.2]
