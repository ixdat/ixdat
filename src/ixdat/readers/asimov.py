"""Reader for ixdat measurements stored in Asimov dataset versions."""                                                                                                                                                                                      
                                                                                                                                                                                                                                                            
import os                                                                                                                                                                                                                                                   
from urllib.parse import urljoin                                                                                                                                                                                                                            

import numpy as np                                                                                                                                                                                                                                          
import requests                                                                                                                                                                                                                                             

from ..auth import KeycloakDeviceTokenProvider
from ..data_series import DataSeries, TimeSeries, ValueSeries
from ..measurement_base import Calculator, Measurement
from ..tools import request_with_retries

DEFAULT_ASIMOV_BASE_URL = "https://asimov.enci.dk/api"
DEFAULT_KEYCLOAK_SERVER_URL = "https://auth.enci.dk"
DEFAULT_KEYCLOAK_REALM = "master"
DEFAULT_KEYCLOAK_CLIENT_ID = "ixdat-cli"

# Candidate time-axis column names, checked in priority order.
_TIME_COLUMN_CANDIDATES = [
    "time/s",
    "t",
    "time",
    "elapsed time/s",
    "elapsed_time/s",
    "Time [s]",
    "Time/s",
]


def _parse_unit_from_name(name):
    """Extract a unit string from an ixdat-style column name.

    Supports ``name/unit`` (e.g. ``"Ewe/V"``), ``name [unit]``
    (e.g. ``"Time [s]"``), and ``name (unit)`` conventions.

    Returns (base_name, unit_string).  *unit_string* is empty when no unit
    can be parsed.
    """
    # "Ewe/V", "time/s", "<I>/mA" — last "/" separates name from unit
    if "/" in name:
        base, unit = name.rsplit("/", 1)
        if unit and " " not in unit:
            return base.strip(), unit.strip()
    # "Time [s]", "Voltage [V]"
    if name.endswith("]") and "[" in name:
        idx = name.index("[")
        return name[:idx].strip(), name[idx + 1 : -1].strip()
    # "Time (s)"
    if name.endswith(")") and "(" in name:
        idx = name.rindex("(")
        return name[:idx].strip(), name[idx + 1 : -1].strip()
    return name, ""


def _find_time_column(series_keys):
    """Return the key that represents the time axis, or *None*.

    Checks known candidate names first, then falls back to a heuristic
    that matches columns whose lower-cased name contains ``"time"``.
    """
    keys_set = set(series_keys)
    for candidate in _TIME_COLUMN_CANDIDATES:
        if candidate in keys_set:
            return candidate
    # Heuristic: first column with "time" in the name
    for key in series_keys:
        if "time" in key.lower():
            return key
    return None


def _extract_tstamp(metadata):
    """Return *tstamp* from *metadata* as a float (unix epoch seconds).

    Handles float, int, numeric strings, and ISO-format datetime strings.
    Falls back to ``0.0`` when nothing usable is found.
    """
    raw = metadata.get("tstamp", 0.0)
    if isinstance(raw, (int, float)):
        return float(raw)
    if isinstance(raw, str):
        try:
            return float(raw)
        except ValueError:
            pass
        # ISO datetime string
        from datetime import datetime

        try:
            return datetime.fromisoformat(raw).timestamp()
        except (ValueError, TypeError):
            pass
    return 0.0


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

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def read(
        self,
        path_to_file,
        cls=Measurement,
        version=None,
        version_id=None,
        force_login=False,
        **kwargs,
    ):
        """Read a dataset version from Asimov as an ixdat Measurement.

        Args:
            path_to_file (str): Asimov dataset id (UUID string).
            cls (Measurement class): Measurement class to instantiate.
            version (int, optional): Version number to select.
            version_id (str, optional): Dataset-version UUID to select.
            force_login (bool): Force fresh Keycloak device login.
            kwargs: Extra key/value pairs merged into the measurement dict.
        """
        dataset_id = str(path_to_file)
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
        return self._measurement_from_payload(
            payload=payload,
            dataset=dataset,
            dataset_version=dataset_version,
            cls=cls,
            kwargs=kwargs,
        )

    # ------------------------------------------------------------------
    # Payload → Measurement conversion
    # ------------------------------------------------------------------

    def _measurement_from_payload(self, payload, dataset, dataset_version, cls, kwargs):
        if not isinstance(payload, dict):
            raise ValueError("payload_json must be a dict.")

        # Primary path: Asimov's canonical format from serialize_measurement()
        #   {technique, measurement_class, series: {col: [...]}, columns, metadata}
        series_map = payload.get("series")
        if isinstance(series_map, dict) and series_map:
            return self._from_series_payload(
                payload, dataset, dataset_version, cls, kwargs
            )

        # Fallback: native ixdat dict with series_list (e.g. manual upload)
        if isinstance(payload.get("series_list"), list):
            return self._from_native_ixdat_payload(
                payload, dataset, dataset_version, cls, kwargs
            )

        raise ValueError(
            "payload_json contains neither a 'series' dict nor a 'series_list'. "
            "Cannot reconstruct a Measurement."
        )

    def _from_series_payload(self, payload, dataset, dataset_version, cls, kwargs):
        """Reconstruct an ixdat Measurement from Asimov's flat series format."""
        series_map = payload["series"]
        metadata = (
            payload["metadata"]
            if isinstance(payload.get("metadata"), dict)
            else {}
        )

        tstamp = _extract_tstamp(metadata)

        # --- Identify the time axis ---
        time_key = _find_time_column(series_map)
        if time_key is None:
            # Last resort: use the first column
            time_key = next(iter(series_map))

        time_data = np.asarray(series_map[time_key], dtype=float)
        _, time_unit = _parse_unit_from_name(time_key)
        if not time_unit:
            time_unit = "s"  # assume seconds when unparseable

        tseries = TimeSeries(
            name=time_key,
            unit_name=time_unit,
            data=time_data,
            tstamp=tstamp,
        )

        # --- Build ValueSeries for every other column ---
        series_list = [tseries]
        for col_name, col_data in series_map.items():
            if col_name == time_key:
                continue
            values = np.asarray(col_data, dtype=float)
            if values.shape != time_data.shape:
                # Different-length columns (e.g. from combined techniques)
                # cannot be linked to this TimeSeries — skip them.
                continue
            _, unit = _parse_unit_from_name(col_name)
            series_list.append(
                ValueSeries(
                    name=col_name,
                    unit_name=unit or "a.u.",
                    data=values,
                    tseries=tseries,
                )
            )

        # --- Technique ---
        technique = (
            payload.get("technique")
            or (dataset_version.get("summary") or {}).get("technique")
            or dataset.get("kind")
            or "simple"
        )

        # --- Aliases ---
        series_names = {s.name for s in series_list}
        aliases = self._build_aliases(technique, series_names)

        # --- Metadata ---
        asimov_meta = self._make_asimov_metadata(dataset, dataset_version)
        merged_metadata = metadata.copy()
        merged_metadata["asimov"] = asimov_meta

        measurement_dict = {
            "name": (
                metadata.get("name")
                or dataset.get("label")
                or str(dataset.get("id"))
            ),
            "technique": technique,
            "series_list": series_list,
            "tstamp": tstamp,
            "aliases": aliases,
            "metadata": merged_metadata,
        }
        if metadata.get("sample_name"):
            measurement_dict["sample_name"] = metadata["sample_name"]
        measurement_dict.update(kwargs)
        return cls.from_dict(measurement_dict)

    def _from_native_ixdat_payload(
        self, payload, dataset, dataset_version, cls, kwargs
    ):
        """Reconstruct from a native ixdat-style dict (series_list of dicts)."""
        measurement_dict = payload.copy()

        for key, obj_cls in [
            ("series_list", DataSeries),
            ("component_measurements", Measurement),
            ("calculator_list", Calculator),
        ]:
            items = measurement_dict.get(key)
            if not isinstance(items, list):
                continue
            measurement_dict[key] = [
                obj_cls.from_dict(item) if isinstance(item, dict) else item
                for item in items
            ]

        if "name" not in measurement_dict:
            measurement_dict["name"] = (
                dataset.get("label") or str(dataset.get("id"))
            )

        meta = measurement_dict.get("metadata") or {}
        if not isinstance(meta, dict):
            meta = {"original_metadata": meta}
        meta["asimov"] = self._make_asimov_metadata(dataset, dataset_version)
        measurement_dict["metadata"] = meta
        measurement_dict.update(kwargs)
        return cls.from_dict(measurement_dict)

    # ------------------------------------------------------------------
    # Alias construction
    # ------------------------------------------------------------------

    @staticmethod
    def _build_aliases(technique, series_names):
        """Build ixdat aliases mapping canonical names to actual column names.

        Constructs aliases for ``t``, ``raw_potential``, ``raw_current``, etc.
        based on the *technique* and which *series_names* are present.
        """
        aliases = {}
        technique_upper = str(technique).upper()

        # Time alias — needed for every technique
        _add_alias(
            aliases,
            "t",
            [
                "time/s", "t", "time", "elapsed time/s",
                "elapsed_time/s", "Time [s]", "Time/s",
            ],
            series_names,
        )

        # EC / CV aliases
        if any(kw in technique_upper for kw in ("EC", "CV")):
            _add_alias(
                aliases,
                "raw_potential",
                [
                    "Ewe/V", "<Ewe>/V", "E/V", "potential/V",
                    "Voltage [V]", "U/V", "WE(1).Potential (V)",
                ],
                series_names,
            )
            _add_alias(
                aliases,
                "raw_current",
                [
                    "I/mA", "<I>/mA", "I/A", "Current/A",
                    "current/mA", "WE(1).Current (A)",
                ],
                series_names,
            )
            _add_alias(
                aliases,
                "raw_CE_potential",
                ["Ece/V", "<Ece>/V"],
                series_names,
            )
            _add_alias(
                aliases,
                "cycle",
                ["cycle number", "cycle", "Ns"],
                series_names,
            )

        return aliases

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _make_asimov_metadata(dataset, dataset_version):
        return {
            "dataset_id": dataset.get("id"),
            "dataset_kind": dataset.get("kind"),
            "dataset_label": dataset.get("label"),
            "dataset_version_id": dataset_version.get("id"),
            "dataset_version": dataset_version.get("version"),
            "dataset_version_created_at": dataset_version.get("created_at"),
            "parser_name": dataset_version.get("parser_name"),
            "ixdat_version": dataset_version.get("ixdat_version"),
        }

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

        return sorted(
            versions,
            key=lambda v: v.get("created_at", ""),
            reverse=True,
        )[0]

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
        return response.json()


# ------------------------------------------------------------------
# Module-level helpers for alias construction
# ------------------------------------------------------------------


def _add_alias(aliases, canonical_name, candidates, series_names):
    """Add *canonical_name* → [matched columns] to *aliases* if any match."""
    if canonical_name in series_names:
        return  # Direct name exists, no alias needed
    matched = [c for c in candidates if c in series_names]
    if not matched:
        matched = _heuristic_alias_match(canonical_name, series_names)
    if matched:
        aliases[canonical_name] = matched


def _heuristic_alias_match(canonical_name, series_names):
    """Last-resort case-insensitive matching for standard series names."""
    if canonical_name == "t":
        for name in series_names:
            low = name.lower()
            if "time" in low and ("/s" in low or "sec" in low or "[s]" in low):
                return [name]
    elif canonical_name == "raw_potential":
        for name in series_names:
            low = name.lower()
            if "ewe" in low or "potential" in low or "voltage" in low:
                return [name]
    elif canonical_name == "raw_current":
        for name in series_names:
            low = name.lower()
            if low.startswith("i/") or low.startswith("<i>/") or "current" in low:
                return [name]
    return []
