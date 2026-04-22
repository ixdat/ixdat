import pytest
import requests

from ixdat.readers.asimov import AsimovReader


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
    reader = AsimovReader(base_url="https://asimov.example.test/api", token="token")
    reader._session = FakeSession(response)

    with pytest.raises(RuntimeError) as exception:
        reader._get_json("datasets/dataset-id", headers={"Authorization": "Bearer t"})

    message = str(exception.value)
    assert "Unauthorized client" in message
    assert "rejected the OAuth client" in message
    assert "KEYCLOAK_CLIENT_ID" in message
    assert "ASIMOV_ACCESS_TOKEN" in message
