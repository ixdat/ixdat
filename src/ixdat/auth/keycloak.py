import json
import os
import subprocess
import sys
import time
from pathlib import Path

import requests

from ..tools import get_default_cache_dir


# Holds references to browser processes we launch so Python does not discard them
# before the browser has finished starting up.
_BROWSER_PROCESSES = []


def _open_in_browser(url):
    """Open ``url`` in the user's default browser without leaking handles."""
    if sys.platform == "win32":
        # ``os.startfile`` hands the URL to Windows shell, which detaches.
        os.startfile(url)  # type: ignore[attr-defined]
        return
    cmd = ["open", url] if sys.platform == "darwin" else ["xdg-open", url]
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
    )
    _BROWSER_PROCESSES.append(proc)


class KeycloakDeviceTokenProvider:
    """Obtain and refresh Keycloak access tokens via the browser-login flow.

    On the first call the user is asked to log in via a browser tab (or a printed
    URL). Two tokens are then saved to a local file:

    - **access token**: a short-lived credential (typically minutes) sent with each
      data request to prove the user is authorised.
    - **refresh token**: a long-lived credential (typically hours/days) used to
      obtain a fresh access token silently, without asking the user to log in again.

    Subsequent calls reuse the cached tokens and only open the browser again when
    both tokens have expired.
    """

    def __init__(
        self,
        server_url,
        realm,
        client_id,
        *,
        client_secret=None,
        scope="openid offline_access",
        cache_path=None,
        open_browser=True,
        connect_timeout=3.0,
        read_timeout=10.0,
    ):
        self.server_url = server_url.rstrip("/")
        self.realm = realm
        self.client_id = client_id
        self.client_secret = client_secret
        self.scope = scope
        self.open_browser = open_browser
        self.connect_timeout = connect_timeout
        self.read_timeout = read_timeout

        if cache_path:
            self.cache_path = Path(cache_path)
        else:
            cache_dir = get_default_cache_dir("ixdat") / "auth" / "keycloak"
            cache_name = f"{self.realm}_{self.client_id}.json"
            self.cache_path = cache_dir / cache_name

        self._session = requests.Session()
        self._session.trust_env = False

    @property
    def token_endpoint(self):
        return f"{self.server_url}/realms/{self.realm}/protocol/openid-connect/token"

    @property
    def device_endpoint(self):
        return (
            f"{self.server_url}/realms/{self.realm}/protocol/openid-connect/auth/device"
        )

    def get_access_token(self, force_login=False):
        """Return a valid access token, refreshing or initiating login if needed."""
        if not force_login:
            tokens = self._load_tokens()
            if tokens:
                if self._access_token_is_valid(tokens):
                    return tokens["access_token"]
                refreshed_tokens = self._try_refresh_tokens(tokens)
                if refreshed_tokens and self._access_token_is_valid(refreshed_tokens):
                    return refreshed_tokens["access_token"]

        new_tokens = self._device_authorization_flow()
        return new_tokens["access_token"]

    def clear_cached_tokens(self):
        """Delete local token cache if present."""
        try:
            self.cache_path.unlink()
        except FileNotFoundError:
            return

    def _access_token_is_valid(self, tokens, leeway=30):
        # Treat the token as expired `leeway` seconds early so it is never sent
        # on a request that starts just before the real expiry time.
        expires_in = tokens.get("expires_in")
        obtained_at = tokens.get("obtained_at")
        access_token = tokens.get("access_token")
        if not access_token or not expires_in or not obtained_at:
            return False
        return time.time() < (obtained_at + expires_in - leeway)

    def _try_refresh_tokens(self, tokens):
        refresh_token = tokens.get("refresh_token")
        if not refresh_token:
            return None

        now = time.time()
        refresh_expires_in = tokens.get("refresh_expires_in")
        obtained_at = tokens.get("obtained_at")
        if (
            refresh_expires_in
            and obtained_at
            and now >= obtained_at + refresh_expires_in
        ):
            return None

        data = {
            "grant_type": "refresh_token",
            "client_id": self.client_id,
            "refresh_token": refresh_token,
        }
        if self.client_secret:
            data["client_secret"] = self.client_secret
        try:
            refreshed = self._post_form(self.token_endpoint, data)
        except RuntimeError:
            return None
        return self._store_tokens(refreshed)

    def _device_authorization_flow(self):
        device_data = {"client_id": self.client_id, "scope": self.scope}
        if self.client_secret:
            device_data["client_secret"] = self.client_secret
        details = self._post_form(self.device_endpoint, device_data)

        verification_uri_complete = details.get("verification_uri_complete")
        verification_uri = details.get("verification_uri")
        user_code = details.get("user_code")
        interval = int(details.get("interval", 5))
        expires_in = int(details.get("expires_in", 600))
        device_code = details.get("device_code")

        if not device_code:
            raise RuntimeError(
                "Keycloak device authorization response did not include `device_code`."
            )

        if verification_uri_complete and self.open_browser:
            _open_in_browser(verification_uri_complete)

        if verification_uri_complete:
            print(
                "Authenticate in your browser to continue:\n"
                f"  {verification_uri_complete}"
            )
        elif verification_uri and user_code:
            print(
                "Authenticate in your browser to continue:\n"
                f"  URL: {verification_uri}\n"
                f"  Code: {user_code}"
            )
        else:
            raise RuntimeError(
                "Keycloak device authorization response did not include a usable "
                "verification URI."
            )

        # Polling: repeatedly ask the server "has the user finished logging in yet?"
        # every `interval` seconds. The server replies with a status code rather
        # than making us wait on an open connection, so we sleep between attempts.
        # Possible outcomes on each check:
        # - 200 OK:
        #      user approved, return tokens and exit
        # - authorization_pending:
        #      user not done yet, wait and retry
        # - slow_down:
        #     server asks us to back off, increase interval and retry
        # - access_denied:
        #     user rejected, raise immediately
        # - expired_token:
        #     login window closed, break and fall through to timeout error
        # Any other error raises immediately. Loop is bounded by `deadline`.
        deadline = time.time() + expires_in
        while time.time() < deadline:
            token_data = {
                "grant_type": "urn:ietf:params:oauth:grant-type:device_code",
                "device_code": device_code,
                "client_id": self.client_id,
            }
            if self.client_secret:
                token_data["client_secret"] = self.client_secret

            response = self._session.post(
                self.token_endpoint,
                data=token_data,
                timeout=(self.connect_timeout, self.read_timeout),
            )
            if response.status_code == 200:
                return self._store_tokens(response.json())

            try:
                error_data = response.json()
            except ValueError:
                error_data = {
                    "error": "unknown_error",
                    "error_description": response.text,
                }

            error = error_data.get("error")
            if error == "authorization_pending":
                time.sleep(interval)
                continue
            if error == "slow_down":
                interval += 5
                time.sleep(interval)
                continue
            if error == "access_denied":
                raise RuntimeError("Keycloak device authentication was denied.")
            if error == "expired_token":
                break

            description = error_data.get("error_description", "No details provided.")
            raise RuntimeError(f"Keycloak device flow failed: {error}: {description}")

        raise RuntimeError("Keycloak device flow timed out before authentication.")

    def _post_form(self, url, data):
        response = self._session.post(
            url,
            data=data,
            timeout=(self.connect_timeout, self.read_timeout),
        )
        if response.status_code >= 400:
            try:
                payload = response.json()
                detail = payload.get("error_description") or payload.get("error")
            except ValueError:
                detail = response.text
            raise RuntimeError(
                f"Token request failed at {url} "
                f"with HTTP {response.status_code}: {detail}"
            )
        try:
            return response.json()
        except ValueError as exc:
            raise RuntimeError("Token request did not return valid JSON.") from exc

    def _load_tokens(self):
        if not self.cache_path.exists():
            return None
        try:
            with open(self.cache_path, "r") as f:
                return json.load(f)
        except (ValueError, OSError):
            return None

    def _store_tokens(self, tokens):
        tokens = tokens.copy()
        tokens["obtained_at"] = time.time()
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.cache_path, "w") as f:
            json.dump(tokens, f)
        try:
            self.cache_path.chmod(0o600)  # owner read/write only. Tokens are credentials
        except OSError:
            pass  # best effort; some platforms (e.g. Windows) ignore Unix permissions
        return tokens
