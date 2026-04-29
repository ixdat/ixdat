"""Reader for ixdat measurements stored in Asimov dataset versions."""

import os
from dataclasses import dataclass, fields
from urllib.parse import urljoin

import numpy as np
import requests

from ..auth import KeycloakDeviceTokenProvider
from ..data_series import DataSeries, TimeSeries, ValueSeries, Field
from ..measurement_base import Measurement
from ..spectra import Spectrum, SpectrumSeries
from ..tools import request_with_retries

OBJECT_TYPE_CLASSES = {
    "measurement": Measurement,
    "spectrum": Spectrum,
    "spectrum_series": SpectrumSeries,
}


@dataclass
class AsimovConfig:
    """Connection and authentication settings for the Asimov reader.

    You usually do not need to touch any of this. The defaults point at the
    public Asimov service and use a normal browser login. Override a setting
    only if your group runs a private Asimov deployment, or if you are
    automating ixdat in a place where a browser cannot be opened.

    Asimov uses an external login server (Keycloak) to verify who you
    are, instead of asking for a username and password directly. The
    first time you read from Asimov, ixdat opens a browser tab where
    you log in. A short-lived access token is then cached on disk and
    reused on subsequent reads, including silent renewal, until the
    refresh window expires.

    Attributes:
        base_url (str): Base URL of the Asimov API.
            Env var: ``ASIMOV_BASE_URL``.
        keycloak_server_url (str): URL of the Keycloak login server.
            Env var: ``KEYCLOAK_SERVER_URL``.
        keycloak_realm (str): Keycloak realm. A realm is Keycloak's term
            for an isolated user database
            Env var: ``KEYCLOAK_REALM``.
        keycloak_client_id (str): The "application name" Asimov knows
            ixdat as. Leave at the default unless your administrator
            tells you otherwise.
            Env var: ``KEYCLOAK_CLIENT_ID``.
        keycloak_client_secret (str, optional): Only required for
            "confidential" clients (e.g. server-to-server automation).
            Env var: ``KEYCLOAK_CLIENT_SECRET``.
        keycloak_scope (str): OAuth scope string requested at login.
            ``offline_access`` is included so that ixdat can silently
            refresh expired tokens without prompting again.
        keycloak_cache_path (str or Path, optional): Where to store the
            cached login tokens on disk. Defaults to ixdat's user cache
            directory under ``auth/keycloak/``.
        keycloak_open_browser (bool): If ``True`` (default), the login URL
            is opened in a browser automatically. Set to ``False`` in
            headless environments; ixdat will print the URL instead so
            you can open it manually on another machine.
        connect_timeout (float): TCP-connect timeout for Asimov API calls,
            in seconds.
        read_timeout (float): Read timeout for Asimov API calls, in
            seconds.
        total_timeout (float): Overall time budget for a single API call
            including retries, in seconds.
    """

    base_url: str = "https://asimov.enci.dk/api"
    keycloak_server_url: str = "https://auth.enci.dk"
    keycloak_realm: str = "master"
    keycloak_client_id: str = "ixdat-cli"
    keycloak_client_secret: "str | None" = None
    keycloak_scope: str = "openid offline_access"
    keycloak_cache_path: "str | None" = None
    keycloak_open_browser: bool = True
    connect_timeout: float = 3.0
    read_timeout: float = 10.0
    total_timeout: float = 30.0

    # Mapping of field name -> environment variable name.
    _ENV_VARS = {
        "base_url": "ASIMOV_BASE_URL",
        "keycloak_server_url": "KEYCLOAK_SERVER_URL",
        "keycloak_realm": "KEYCLOAK_REALM",
        "keycloak_client_id": "KEYCLOAK_CLIENT_ID",
        "keycloak_client_secret": "KEYCLOAK_CLIENT_SECRET",
    }

    def __post_init__(self):
        self.base_url = self.base_url.rstrip("/")
        self.keycloak_server_url = self.keycloak_server_url.rstrip("/")

    @classmethod
    def from_env(cls):
        """Build a config, layering environment variables on top of defaults."""
        kwargs = {}
        for f in fields(cls):
            env_name = cls._ENV_VARS.get(f.name)
            if env_name and os.environ.get(env_name):
                kwargs[f.name] = os.environ[env_name]
        return cls(**kwargs)


# Module-level config, which is populated from environment variables at import.
# Edit attributes on this object before the first reader call to override.
ASIMOV_CONFIG = AsimovConfig.from_env()


class AsimovReader:
    """Read a remote dataset version from Asimov into an ixdat object.

    Typical use only needs the dataset id: ``Measurement.read(<id>,
    reader="asimov")``. The first read opens a browser tab where you log
    in once; later reads reuse a cached token.

    Power users can pass a pre-obtained bearer ``token`` (or set
    ``ASIMOV_ACCESS_TOKEN``), or supply a custom ``token_provider`` for
    tests and CI. Connection and Keycloak settings live on
    :class:`AsimovConfig`; see its docstring for what each field means and
    which environment variables it reads.
    """

    def __init__(self, *, token=None, token_provider=None, config=None):
        """Build an AsimovReader.

        Args:
            token (str, optional): A pre-obtained bearer access token. If
                set (or if ``ASIMOV_ACCESS_TOKEN`` is in the environment),
                no Keycloak login is performed.
            token_provider (object, optional): Custom provider exposing
                ``get_access_token(force_login: bool) -> str``. Mostly for
                tests and alternative auth backends.
            config (AsimovConfig, optional): Connection and Keycloak
                settings. Defaults to the module-level ``ASIMOV_CONFIG``.
        """
        self.config = config or ASIMOV_CONFIG
        self.base_url = self.config.base_url
        self.connect_timeout = self.config.connect_timeout
        self.read_timeout = self.config.read_timeout
        self.total_timeout = self.config.total_timeout

        self._token = token or os.environ.get("ASIMOV_ACCESS_TOKEN")
        self.token_provider = token_provider

        if self._token is None and self.token_provider is None:
            self.token_provider = KeycloakDeviceTokenProvider(
                server_url=self.config.keycloak_server_url,
                realm=self.config.keycloak_realm,
                client_id=self.config.keycloak_client_id,
                client_secret=self.config.keycloak_client_secret,
                scope=self.config.keycloak_scope,
                cache_path=self.config.keycloak_cache_path,
                open_browser=self.config.keycloak_open_browser,
                connect_timeout=self.config.connect_timeout,
                read_timeout=self.config.read_timeout,
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
                If None, auto-detected from the payload's object_type field.
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

        # Inject Asimov provenance (dataset versioning) into metadata
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

        return cls.from_dict(self._build_kwargs(d, reader=self, **kwargs))

    def _build_kwargs(self, dct, **kwargs):
        """Translate an Asimov payload into kwargs for cls.from_dict().
        Measurements need an absolute tstamp; Spectrum / SpectrumSeries
        accept tstamp=None natively, so we don't enforce it for them.
        """
        if "series_list" in dct and dct.get("tstamp") is None:
            raise ValueError(
                "Asimov measurement payload is missing 'tstamp'. ixdat requires "
                "an absolute timestamp on every Measurement."
            )
        obj = {
            "name": dct.get("name"),
            "technique": dct.get("technique"),
            "metadata": dct.get("metadata") or {},
            "tstamp": dct.get("tstamp"),
            **kwargs,
        }
        if dct.get("sample_name") is not None:
            obj["sample_name"] = dct["sample_name"]

        key_map = None
        if "series_list" in dct:
            # Build TimeSeries first so vseries / fields can reference them.
            key_map = {}
            series_list = []
            for s in dct["series_list"]:
                if s.get("series_type") == "tseries":
                    ts = self._build_series(s)
                    key_map[s["key"]] = ts
                    series_list.append(ts)
            for s in dct["series_list"]:
                if s.get("series_type") != "tseries":
                    built = self._build_series(s, key_map)
                    key_map[s["key"]] = built
                    series_list.append(built)
            obj["series_list"] = series_list
            obj["aliases"] = dct.get("aliases") or {}

        # this will be the case when retreiving a Spectrum, SpectrumSeries or SpectroMeasurement
        if "field" in dct:
            obj["field"] = self._build_series(dct["field"], key_map)
            obj["duration"] = dct.get("duration")
            if dct.get("object_type") == "spectrum_series":
                obj["durations"] = dct.get("durations")
                obj["continuous"] = dct.get("continuous", False)

        return obj

    @staticmethod
    def _build_series(dct, key_map=None):
        """Build a DataSeries from one payload entry.

        key_map resolves tseries_key / axes_keys references when
        the entry comes from a Measurement's series_list.
        """
        kind = dct.get("series_type", "series")
        name = dct["name"]
        unit_name = dct["unit_name"]
        data = np.asarray(dct["data"])

        if kind == "tseries":
            if dct.get("tstamp") is None:
                raise ValueError(
                    f"Asimov tseries '{name}' is missing 'tstamp'. ixdat requires "
                    "an absolute timestamp on every TimeSeries."
                )
            return TimeSeries(
                name=name, unit_name=unit_name, data=data, tstamp=dct["tstamp"]
            )
        if kind in ("vseries", "constantvalue"):
            ts = None
            if key_map and "tseries_key" in dct:
                ts = key_map.get(dct["tseries_key"])
            return ValueSeries(name=name, unit_name=unit_name, data=data, tseries=ts)
        if kind == "field":
            # axes_keys references axes already built from the shared series_list;
            # axes_series is the inline form where axes are embedded in the payload
            # directly (standalone Spectrum / SpectrumSeries with no series_list).
            if "axes_keys" in dct and key_map:
                axes = [key_map[k] for k in dct["axes_keys"]]
            else:
                # No key_map means axes are not in a shared series_list, so they
                # must be embedded inline under axes_series.
                axes = [
                    AsimovReader._build_series(a, key_map)
                    for a in dct.get("axes_series", [])
                ]
            return Field(name=name, unit_name=unit_name, data=data, axes_series=axes)

        return DataSeries(name=name, unit_name=unit_name, data=data)

    def _select_dataset_version(self, versions, version=None, version_id=None):
        """Return the dataset version dict to load.

        A dataset in Asimov can have multiple versions: each time data is
        re-processed or re-uploaded, a new version is created while older
        ones are kept for provenance (so you can always trace back to exactly
        what data was used in a given analysis). In practice you almost always
        want the latest version, which is the default when neither ``version``
        nor ``version_id`` is given.

        Args:
            versions (list): List of version dicts returned by the API.
            version (int, optional): Version number to select (1-based).
            version_id (str, optional): Exact version UUID to select.
        """
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
        # Asimov returns created_at as an ISO 8601 string (e.g. "2024-06-01T12:00:00Z").
        # Lexicographic sorting of ISO 8601 strings is equivalent to chronological
        # sorting, so reverse=True gives the most recently created version first.
        return sorted(versions, key=lambda v: v.get("created_at", ""), reverse=True)[0]

    def _build_auth_headers(self, force_login=False):
        """Return the HTTP Authorization header dict for an API request.

        Uses a pre-set static token if available, otherwise asks the
        token provider (Keycloak) for a valid access token, triggering
        a browser login if the cached token has expired.
        """
        if self._token:
            return {"Authorization": f"Bearer {self._token}"}
        if self.token_provider:
            token = self.token_provider.get_access_token(force_login=force_login)
            return {"Authorization": f"Bearer {token}"}
        raise RuntimeError(
            "No authentication configured. Provide `token` or `token_provider`."
        )

    def _load_payload_uri(self, payload_uri, headers):
        """Fetch and return a payload stored at a URI rather than inline.

        Large datasets are sometimes stored outside the main API response
        and referenced by a URI. This method resolves relative URIs against
        the base URL and downloads the payload JSON.
        """
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
        """GET an Asimov API endpoint and return the parsed JSON response.

        Constructs the full URL from the base URL and endpoint path, then
        delegates to ``request_with_retries``. Provides a clearer error
        message if the server rejects the OAuth client (HTTP 401).
        """
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
                    "The configured ixdat client is "
                    f"{self.config.keycloak_client_id!r}; "
                    "check that this client is allowed by the Asimov API, or set "
                    "KEYCLOAK_CLIENT_ID/ASIMOV_ACCESS_TOKEN to an API-authorized "
                    "client/token."
                ) from exc
            raise
        return response.json()
