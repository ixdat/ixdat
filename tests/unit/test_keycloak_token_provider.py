from unittest.mock import patch

from ixdat.auth.keycloak import KeycloakDeviceTokenProvider


def test_keycloak_provider_uses_cached_valid_access_token(tmp_path):
    cache_path = tmp_path / "token-cache.json"
    provider = KeycloakDeviceTokenProvider(
        server_url="https://auth.example.test",
        realm="master",
        client_id="ixdat-cli",
        cache_path=cache_path,
        open_browser=False,
    )
    provider._store_tokens(
        {
            "access_token": "cached-access-token",
            "refresh_token": "cached-refresh-token",
            "expires_in": 3600,
            "refresh_expires_in": 7200,
        }
    )

    with patch.object(
        provider,
        "_device_authorization_flow",
        side_effect=AssertionError("device flow should not be called"),
    ):
        token = provider.get_access_token()

    assert token == "cached-access-token"
