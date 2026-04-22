"""Reader for ixdat measurements stored in Asimov dataset versions."""

import os
from urllib.parse import urljoin

import requests

from ..auth import KeycloakDeviceTokenProvider
from ..measurement_base import Measurement
from ..spectra import Spectrum, SpectrumSeries
from ..tools import request_with_retries

OBJECT_TYPE_CLASSES = {
    "measurement": Measurement,
    "spectrum": Spectrum,
    "spectrum_series": SpectrumSeries,
}

DEFAULT_ASIMOV_BASE_URL = "https://asimov.enci.dk/api"
DEFAULT_KEYCLOAK_SERVER_URL = "https://auth.enci.dk"
DEFAULT_KEYCLOAK_REALM = "master"
DEFAULT_KEYCLOAK_CLIENT_ID = "ixdat-cli"


class AsimovReader:
    """Read a remote dataset version from Asimov into an ixdat Measurement.

    Authentication is done via bearer token:

    - Direct token: pass ``token=...`` or set ``ASIMOV_ACCESS_TOKEN``.
    - Keycloak provider: pass ``token_provider=...``.
    - Auto-initialise Keycloak device flow with Asimov defaults.
      Optional env-var overrides:
      ``KEYCLOAK_SERVER_URL``, ``KEYCLOAK_REALM``, ``KEYCLOAK_CLIENT_ID``.
    """

    def __init__(
        self,
        base_url=None,
        token=None,
        token_provider=None,
        keycloak_server_url=None,
        keycloak_realm=None,
        keycloak_client_id=None,
        keycloak_client_secret=None,
        keycloak_scope="openid offline_access",
        keycloak_cache_path=None,
        keycloak_open_browser=True,
        connect_timeout=3.0,
        read_timeout=10.0,
        total_timeout=30.0,
    ):
        self.base_url = (
            base_url or os.environ.get("ASIMOV_BASE_URL") or DEFAULT_ASIMOV_BASE_URL
        ).rstrip("/")

        self._token = token or os.environ.get("ASIMOV_ACCESS_TOKEN")
        self.token_provider = token_provider
        self.connect_timeout = connect_timeout
        self.read_timeout = read_timeout
        self.total_timeout = total_timeout

        if self._token is None and self.token_provider is None:
            keycloak_server_url = (
                keycloak_server_url
                or os.environ.get("KEYCLOAK_SERVER_URL")
                or DEFAULT_KEYCLOAK_SERVER_URL
            )
            keycloak_realm = (
                keycloak_realm
                or os.environ.get("KEYCLOAK_REALM")
                or DEFAULT_KEYCLOAK_REALM
            )
            keycloak_client_id = (
                keycloak_client_id
                or os.environ.get("KEYCLOAK_CLIENT_ID")
                or DEFAULT_KEYCLOAK_CLIENT_ID
            )
            keycloak_client_secret = keycloak_client_secret or os.environ.get(
                "KEYCLOAK_CLIENT_SECRET"
            )

            self.token_provider = KeycloakDeviceTokenProvider(
                server_url=keycloak_server_url.rstrip("/"),
                realm=keycloak_realm,
                client_id=keycloak_client_id,
                client_secret=keycloak_client_secret,
                scope=keycloak_scope,
                cache_path=keycloak_cache_path,
                open_browser=keycloak_open_browser,
                connect_timeout=connect_timeout,
                read_timeout=read_timeout,
            )

        self.dataset = None
        self.dataset_version = None
        self._session = requests.Session()
        self._session.trust_env = False

    def read(
        self,
        id,
        cls=None,
        version=None,
        version_id=None,
        force_login=False,
        **kwargs,
    ):
        """Read a dataset version from Asimov as an ixdat object.

        Args:
            id (str): Asimov dataset id (UUID string).
            cls (class, optional): Class to instantiate (Measurement, Spectrum, etc.).
                If None, auto-detected from the payload's ``object_type`` field.
            version (int, optional): Version number to select.
            version_id (str, optional): Dataset-version UUID to select.
            force_login (bool): Force fresh Keycloak device login.
            kwargs: Extra key/value pairs merged into the object dict.
        """
        dataset_id = str(id)
        headers = self._build_auth_headers(force_login=force_login)

        dataset = self._get_json(f"datasets/{dataset_id}", headers=headers)
        versions = self._get_json(
            "dataset-versions",
            headers=headers,
            params={"dataset_id": dataset_id},
        )
        dataset_version = self._select_dataset_version(
            versions, version=version, version_id=version_id
        )

        payload = dataset_version.get("payload_json")
        if not payload and dataset_version.get("payload_uri"):
            payload = self._load_payload_uri(
                dataset_version["payload_uri"], headers=headers
            )
        if not payload:
            raise ValueError(
                f"Dataset version {dataset_version.get('id')} has no "
                "payload_json or payload_uri."
            )

        self.dataset = dataset
        self.dataset_version = dataset_version

        # Inject Asimov provenance into metadata
        meta = dict(payload.get("metadata") or {})
        meta["asimov"] = {
            "dataset_id": dataset.get("id"),
            "dataset_kind": dataset.get("kind"),
            "dataset_label": dataset.get("label"),
            "dataset_version_id": dataset_version.get("id"),
            "dataset_version": dataset_version.get("version"),
            "dataset_version_created_at": dataset_version.get("created_at"),
            "parser_name": dataset_version.get("parser_name"),
            "ixdat_version": dataset_version.get("ixdat_version"),
        }

        d = {**payload, "metadata": meta}
        if not d.get("name"):
            d["name"] = dataset.get("label") or str(dataset.get("id"))

        if cls is None:
            object_type = d.get("object_type", "measurement")
            cls = OBJECT_TYPE_CLASSES.get(object_type, Measurement)

        return cls.from_portable_dict(d, reader=self, **kwargs)

    def _select_dataset_version(self, versions, version=None, version_id=None):
        if not isinstance(versions, list) or not versions:
            raise ValueError("No dataset versions returned from Asimov API.")
        if version_id is not None:
            for v in versions:
                if v.get("id") == version_id:
                    return v
            raise ValueError(f"No dataset version with id={version_id}.")
        if version is not None:
            for v in versions:
                if v.get("version") == version:
                    return v
            raise ValueError(f"No dataset version with version={version}.")
        return sorted(versions, key=lambda v: v.get("created_at", ""), reverse=True)[0]

    def _build_auth_headers(self, force_login=False):
        if self._token:
            return {"Authorization": f"Bearer {self._token}"}
        if self.token_provider:
            token = self.token_provider.get_access_token(force_login=force_login)
            return {"Authorization": f"Bearer {token}"}
        raise RuntimeError(
            "No authentication configured. Provide `token` or `token_provider`."
        )

    def _load_payload_uri(self, payload_uri, headers):
        if payload_uri.startswith(("http://", "https://")):
            url = payload_uri
        else:
            url = urljoin(self.base_url + "/", payload_uri.lstrip("/"))
        response = request_with_retries(
            self._session,
            "GET",
            url,
            headers=headers,
            connect_timeout=self.connect_timeout,
            read_timeout=self.read_timeout,
            total_timeout=self.total_timeout,
            retries=2,
        )
        return response.json()

    def _get_json(self, endpoint, headers, params=None):
        url = urljoin(self.base_url + "/", endpoint.lstrip("/"))
        try:
            response = request_with_retries(
                self._session,
                "GET",
                url,
                headers=headers,
                params=params,
                connect_timeout=self.connect_timeout,
                read_timeout=self.read_timeout,
                total_timeout=self.total_timeout,
                retries=2,
            )
        except RuntimeError as exc:
            if "HTTP 401" in str(exc) and "Unauthorized client" in str(exc):
                raise RuntimeError(
                    str(exc)
                    + "\n\nAsimov rejected the OAuth client used for this token. "
                    f"The default ixdat client is {DEFAULT_KEYCLOAK_CLIENT_ID!r}; "
                    "check that this client is allowed by the Asimov API, or set "
                    "KEYCLOAK_CLIENT_ID/ASIMOV_ACCESS_TOKEN to an API-authorized "
                    "client/token."
                ) from exc
            raise
        return response.json()
