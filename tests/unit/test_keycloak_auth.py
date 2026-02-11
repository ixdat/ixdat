import json
import time

from ixdat.auth.keycloak import KeycloakDeviceTokenProvider


def _build_provider(
    *,
    server_url="https://keycloak.example.org",
    realm="example",
    client_id="ixdat-client",
    cache_path=None,
):
    return KeycloakDeviceTokenProvider(
        server_url=server_url,
        realm=realm,
        client_id=client_id,
        cache_path=cache_path,
        open_browser=False,
    )


def _unexpected_device_flow():
    raise AssertionError("should not trigger device flow")


class _PostResponse:
    def __init__(self, status_code, payload=None, text=""):
        self.status_code = status_code
        self._payload = payload or {}
        self.text = text

    def json(self):
        return self._payload


def test_keycloak_provider_reuses_valid_cached_access_token(tmp_path, monkeypatch):
    cache_path = tmp_path / "token_cache.json"
    token_data = {
        "access_token": "cached-token",
        "refresh_token": "cached-refresh",
        "expires_in": 3600,
        "refresh_expires_in": 7200,
        "obtained_at": time.time(),
    }
    with open(cache_path, "w") as f:
        json.dump(token_data, f)

    provider = _build_provider(cache_path=cache_path)

    monkeypatch.setattr(provider, "_device_authorization_flow", _unexpected_device_flow)

    assert provider.get_access_token() == "cached-token"


def test_keycloak_provider_uses_refresh_before_device_flow(tmp_path, monkeypatch):
    cache_path = tmp_path / "token_cache.json"
    token_data = {
        "access_token": "expired-token",
        "refresh_token": "cached-refresh",
        "expires_in": 10,
        "refresh_expires_in": 7200,
        "obtained_at": time.time() - 300,
    }
    with open(cache_path, "w") as f:
        json.dump(token_data, f)

    provider = _build_provider(cache_path=cache_path)

    refreshed_tokens = {
        "access_token": "refreshed-token",
        "refresh_token": "refreshed-refresh-token",
        "expires_in": 3600,
        "refresh_expires_in": 7200,
        "obtained_at": time.time(),
    }

    monkeypatch.setattr(provider, "_try_refresh_tokens", lambda tokens: refreshed_tokens)
    monkeypatch.setattr(provider, "_device_authorization_flow", _unexpected_device_flow)

    assert provider.get_access_token() == "refreshed-token"


def test_keycloak_provider_uses_realm_based_endpoints():
    provider = _build_provider(
        server_url="https://auth.enci.dk",
        realm="master",
        client_id="ixdat-cli",
    )

    assert provider.token_endpoint.endswith("/realms/master/protocol/openid-connect/token")
    assert provider.device_endpoint.endswith(
        "/realms/master/protocol/openid-connect/auth/device"
    )


def test_keycloak_provider_device_flow_uses_expected_endpoints(tmp_path, monkeypatch):
    provider = _build_provider(
        server_url="https://auth.enci.dk",
        realm="master",
        client_id="ixdat-cli",
        cache_path=tmp_path / "token_cache.json",
    )

    device_calls = []
    token_calls = []

    def fake_post(url, data, timeout):
        if "openid-connect/auth/device" in url:
            device_calls.append(url)
            return _PostResponse(
                200,
                {
                    "device_code": "abc123",
                    "verification_uri_complete": "https://login.example.org/device?user_code=XYZ",
                    "interval": 1,
                    "expires_in": 600,
                },
            )
        if "openid-connect/token" in url:
            token_calls.append(url)
            return _PostResponse(
                200,
                {
                    "access_token": "fresh-access-token",
                    "refresh_token": "fresh-refresh-token",
                    "expires_in": 3600,
                    "refresh_expires_in": 7200,
                },
            )
        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr(provider._session, "post", fake_post)

    token = provider.get_access_token(force_login=True)

    assert token == "fresh-access-token"
    assert device_calls[0].startswith("https://auth.enci.dk/realms/")
    assert token_calls[0].startswith("https://auth.enci.dk/realms/")
