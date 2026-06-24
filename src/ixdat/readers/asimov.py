"""Reader for ixdat objects stored as parsed Asimov project-file payloads."""

import os
from dataclasses import dataclass, fields
from urllib.parse import urljoin

import numpy as np
import requests

from ..auth import KeycloakDeviceTokenProvider
from ..data_series import DataSeries, SERIES_CLASSES
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
    """Read a remote project file or bundle parse result into an ixdat object.

    Typical use only needs the project file or bundle id: ``Measurement.read(<id>,
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

        self.payload_envelope = None
        self._session = requests.Session()
        self._session.trust_env = False

    def read(
        self,
        id,
        cls=None,
        force_login=False,
        **kwargs,
    ):
        """Read an ixdat payload from Asimov.

        Args:
            id (str): Asimov project-file or project-file-bundle id (UUID string).
            cls (class, optional): Class to instantiate (Measurement, Spectrum, etc.).
                If None, auto-detected from the payload's object_type field.
            force_login (bool): Force fresh Keycloak device login.
            kwargs: Extra key/value pairs merged into the object dict.
        """
        asimov_id = str(id)
        headers = self._build_auth_headers(force_login=force_login)

        payload_envelope = self._get_ixdat_payload_envelope(asimov_id, headers=headers)
        payload = payload_envelope.get("payload_json")
        if not payload and payload_envelope.get("payload_uri"):
            payload = self._load_payload_uri(
                payload_envelope["payload_uri"], headers=headers
            )
        if not payload:
            raise ValueError(
                f"Asimov parsed payload response for {asimov_id} has no "
                "payload_json or payload_uri."
            )

        self.payload_envelope = payload_envelope

        # Inject Asimov provenance from the native project-file/bundle envelope.
        meta = dict(payload.get("metadata") or {})
        meta["asimov"] = {
            "project_file_id": payload_envelope.get("project_file_id"),
            "bundle_id": payload_envelope.get("bundle_id"),
            "filename": payload_envelope.get("filename"),
            "name": payload_envelope.get("name"),
            "parser_name": payload_envelope.get("parser_name"),
            "created_at": payload_envelope.get("created_at"),
            "ixdat_version": payload_envelope.get("ixdat_version"),
        }

        d = {**payload, "metadata": meta}
        if not d.get("name"):
            d["name"] = (
                payload_envelope.get("filename")
                or payload_envelope.get("name")
                or payload_envelope.get("project_file_id")
                or payload_envelope.get("bundle_id")
                or asimov_id
            )

        object_type = d.get("object_type", "measurement")
        self._check_requested_class_matches_object_type(
            requested_class=cls, object_type=object_type, asimov_id=asimov_id
        )

        return self._build_object(d, cls=cls, reader=self, **kwargs)

    @staticmethod
    def _check_requested_class_matches_object_type(
        requested_class, object_type, asimov_id
    ):
        """Raise a clear error for the wrong read entrypoint."""
        if requested_class is None:
            return

        if object_type == "measurement":
            expected_cls = Measurement
            recommended_read = "Measurement.read"
            compatible = issubclass(requested_class, Measurement)
        elif object_type == "spectrum":
            expected_cls = Spectrum
            recommended_read = "Spectrum.read"
            compatible = issubclass(requested_class, Spectrum) and not issubclass(
                requested_class, SpectrumSeries
            )
        elif object_type == "spectrum_series":
            expected_cls = SpectrumSeries
            recommended_read = "Spectrum.read or SpectrumSeries.read"
            compatible = issubclass(requested_class, Spectrum)
        else:
            raise ValueError(
                f"Asimov payload {asimov_id} has unsupported "
                f"object_type={object_type!r}. "
                f"Supported object types are {sorted(OBJECT_TYPE_CLASSES)}."
            )

        if compatible:
            return

        raise ValueError(
            f"Asimov payload {asimov_id} contains a {expected_cls.__name__} payload "
            f"(object_type={object_type!r}), but {requested_class.__name__}.read(..., "
            "reader='asimov') requested an incompatible ixdat object type. "
            f"Use {recommended_read}(..., reader='asimov') for this payload instead."
        )

    def _get_ixdat_payload_envelope(self, asimov_id, headers):
        """Return the native Asimov ixdat-payload envelope for a file or bundle."""
        try:
            return self._get_json(
                f"project-files/{asimov_id}/ixdat-payload", headers=headers
            )
        except RuntimeError as exc:
            if not self._is_http_status(exc, 404):
                raise

        try:
            return self._get_json(
                f"project-file-bundles/{asimov_id}/ixdat-payload", headers=headers
            )
        except RuntimeError as exc:
            if self._is_http_status(exc, 404):
                raise ValueError(
                    f"No Asimov project file or project file bundle with id={asimov_id}."
                ) from exc
            raise

    @staticmethod
    def _is_http_status(exc, status_code):
        """Return True if a request RuntimeError mentions an HTTP status code."""
        return f"HTTP {status_code}" in str(exc)

    def _build_kwargs(self, dct, object_class=None, **kwargs):
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
            # Each series entry carries a "key" assigned by Asimov, a unique
            # identifier used to cross-reference series within the same payload
            # (e.g. a ValueSeries points to its TimeSeries via "tseries_key").
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

        # this will be the case when retreiving a Spectrum, SpectrumSeries or
        # SpectroMeasurement
        if "field" in dct:
            obj["field"] = self._build_series(dct["field"], key_map)
            obj["duration"] = dct.get("duration")
            if dct.get("object_type") == "spectrum_series":
                obj["durations"] = dct.get("durations")
                obj["continuous"] = dct.get("continuous", False)

        for key, value in dct.items():
            # A nested dict with object_type is a complete ixdat object payload,
            # such as a reference_spectrum or spectrum_series inside a Measurement.
            # Build it recursively without caring what the attribute is called.
            if key not in obj and isinstance(value, dict) and "object_type" in value:
                obj[key] = self._build_object(value, reader=self)
            # Some ixdat attributes are lists of child objects, for example
            # component_measurements. Hydrate only the list entries that declare
            # object_type and leave ordinary scalar/list metadata untouched.
            elif (
                key not in obj
                and isinstance(value, list)
                and any(isinstance(v, dict) and "object_type" in v for v in value)
            ):
                obj[key] = [
                    self._build_object(v, reader=self)
                    if isinstance(v, dict) and "object_type" in v
                    else v
                    for v in value
                ]
            # If a nested dict has the shape of an ixdat object but omits
            # object_type, the payload is ambiguous. ASIMOV should add object_type
            # rather than ixdat guessing from field names or data shape.
            elif (
                key not in obj
                and isinstance(value, dict)
                and "object_type" not in value
                and ("series_list" in value or "field" in value)
            ):
                raise ValueError(
                    f"Asimov payload field {key!r} is missing 'object_type'. "
                    "Nested ixdat objects must declare object_type so ixdat can "
                    "reconstruct them without key-specific reader logic."
                )
            else:
                # Everything else is ordinary payload metadata or an unsupported
                # auxiliary shape. It is ignored here unless ixdat declares it as a
                # serializable constructor field in the block below.
                pass

        for key in self._get_serializable_attrs(dct, object_class):
            if key not in obj and key in dct:
                obj[key] = dct[key]

        return obj

    @staticmethod
    def _get_serializable_attrs(dct, object_class=None):
        """Return ixdat-declared constructor fields relevant to a payload."""
        if object_class is None:
            return set()
        serializable_attrs = set(object_class.get_all_column_attrs())

        technique = dct.get("technique")
        if technique:
            from ..techniques import TECHNIQUE_CLASSES

            technique_class = TECHNIQUE_CLASSES.get(technique)
            if technique_class and issubclass(technique_class, object_class):
                serializable_attrs.update(technique_class.get_all_column_attrs())

        return serializable_attrs

    def _build_object(self, dct, cls=None, **kwargs):
        """Build an ixdat Measurement, Spectrum, or SpectrumSeries from a payload."""
        object_type = dct.get("object_type", "measurement")
        if cls is None:
            cls = OBJECT_TYPE_CLASSES.get(object_type)
            if cls is None:
                raise ValueError(
                    f"Asimov payload has unsupported object_type={object_type!r}. "
                    f"Supported object types are {sorted(OBJECT_TYPE_CLASSES)}."
                )
        if cls is Spectrum and object_type == "spectrum_series":
            cls = SpectrumSeries
        return cls.from_dict(self._build_kwargs(dct, object_class=cls, **kwargs))

    @staticmethod
    def _build_series(dct, key_map=None):
        """Build a DataSeries from one payload entry.

        key_map resolves tseries_key / axes_keys references when
        the entry comes from a Measurement's series_list.
        """
        kind = dct.get("series_type", "series")
        series_cls = SERIES_CLASSES.get(kind)
        if series_cls is None:
            raise ValueError(
                f"Asimov series '{dct.get('name')}' has unsupported "
                f"series_type={kind!r}. Supported series types are "
                f"{sorted(SERIES_CLASSES)}."
            )
        serializable_attrs = set(series_cls.get_all_column_attrs())
        series_dict = {
            key: value for key, value in dct.items() if key in serializable_attrs
        }
        series_dict["series_type"] = kind
        series_dict["data"] = np.asarray(series_dict["data"])

        if kind == "tseries" and dct.get("tstamp") is None:
            raise ValueError(
                f"Asimov tseries '{dct.get('name')}' is missing 'tstamp'. "
                "ixdat requires an absolute timestamp on every TimeSeries."
            )
        if "tseries_key" in dct and key_map:
            series_dict["tseries"] = key_map[dct["tseries_key"]]
        if "axes_keys" in dct and key_map:
            series_dict["axes_series"] = [key_map[k] for k in dct["axes_keys"]]
        elif "axes_series" in dct:
            series_dict["axes_series"] = [
                AsimovReader._build_series(axis, key_map) for axis in dct["axes_series"]
            ]

        return DataSeries.from_dict(series_dict)

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

        Large payloads are sometimes stored outside the main API response
        and referenced by a URI. This method resolves relative URIs against
        the base URL and downloads the payload JSON.
        """
        if payload_uri.startswith(("http://", "https://")):
            url = payload_uri
        else:
            url = urljoin(self.base_url + "/", payload_uri.lstrip("/"))
        # response is a requests.Response object, i.e. the raw HTTP reply from the
        # server. It contains the status code (200 = OK), headers, and a body with
        # the actual data. Calling .json() on it reads that body and converts it
        # from text into a Python dictionary.
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
