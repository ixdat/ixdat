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


class AsimovReader:
    """Read a remote dataset version from Asimov into an ixdat Measurement.

    Authentication is done via bearer token:
    - Direct token: pass ``token=...`` or set ``ASIMOV_ACCESS_TOKEN``
    - Keycloak provider: pass ``token_provider=...``
    - Auto-initialize Keycloak device flow with Asimov defaults.
      Optional env var overrides:
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
        path_to_file,
        cls=Measurement,
        version=None,
        version_id=None,
        force_login=False,
        **kwargs,
    ):
        """Read dataset version identified by ``path_to_file`` as dataset_id.

        Args:
            path_to_file (str): Asimov dataset id (UUID string).
            cls (Measurement class): Measurement class to instantiate.
            version (int, optional): Explicit version number to select.
            version_id (str, optional): Explicit dataset-version id to select.
            force_login (bool): Force fresh Keycloak device login if provider is used.
            kwargs: Extra key/value pairs passed into the resulting measurement dict.
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
            payload = self._load_payload_uri(dataset_version["payload_uri"], headers=headers)
        if not payload:
            raise ValueError(
                f"Asimov dataset version {dataset_version.get('id')} does not include "
                "`payload_json` or `payload_uri`."
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

    def _measurement_from_payload(self, payload, dataset, dataset_version, cls, kwargs):
        if self._looks_like_ixdat_measurement(payload):
            measurement_dict = self._deserialize_measurement_dict(payload)
            if "name" not in measurement_dict:
                measurement_dict["name"] = dataset.get("label") or str(dataset.get("id"))
            metadata = measurement_dict.get("metadata") or {}
            if not isinstance(metadata, dict):
                metadata = {"original_metadata": metadata}
            metadata["asimov"] = self._make_asimov_metadata(dataset, dataset_version)
            measurement_dict["metadata"] = metadata
            measurement_dict = self._normalize_measurement_aliases(measurement_dict)
            measurement_dict.update(kwargs)
            return cls.from_dict(measurement_dict)
        return self._fallback_measurement_from_series_payload(
            payload=payload,
            dataset=dataset,
            dataset_version=dataset_version,
            cls=cls,
            kwargs=kwargs,
        )

    def _fallback_measurement_from_series_payload(
        self, payload, dataset, dataset_version, cls, kwargs
    ):
        if not isinstance(payload, dict):
            raise ValueError("Asimov payload_json must be a dict.")

        series_map = payload.get("series")
        if not isinstance(series_map, dict) or not series_map:
            raise ValueError(
                "Asimov payload_json could not be interpreted as an ixdat "
                "measurement, and does not contain `series` data."
            )

        axis_name = (
            "time"
            if "time" in series_map
            else ("frequency" if "frequency" in series_map else next(iter(series_map)))
        )
        axis_data = np.asarray(series_map[axis_name])
        tstamp = (
            (payload.get("metadata") or {}).get("tstamp", 0.0)
            if isinstance(payload.get("metadata"), dict)
            else 0.0
        )
        axis_unit = "s" if axis_name == "time" else "a.u."
        tseries = TimeSeries(
            name=axis_name,
            unit_name=axis_unit,
            data=axis_data,
            tstamp=tstamp,
        )
        series_list = [tseries]
        for series_name, series_data in series_map.items():
            if series_name == axis_name:
                continue
            values = np.asarray(series_data)
            if values.shape != axis_data.shape:
                raise ValueError(
                    f"Series '{series_name}' shape {values.shape} did not match "
                    f"axis '{axis_name}' shape {axis_data.shape}."
                )
            series_list.append(
                ValueSeries(
                    name=series_name,
                    unit_name="a.u.",
                    data=values,
                    tseries=tseries,
                )
            )

        technique = (
            (dataset_version.get("summary") or {}).get("technique")
            or dataset.get("kind")
            or "generic"
        )
        metadata = payload.get("metadata") if isinstance(payload.get("metadata"), dict) else {}
        metadata = metadata.copy()
        metadata["asimov"] = self._make_asimov_metadata(dataset, dataset_version)
        measurement_dict = {
            "name": dataset.get("label") or str(dataset.get("id")),
            "technique": technique,
            "series_list": series_list,
            "tstamp": tstamp,
            "metadata": metadata,
        }
        measurement_dict = self._normalize_measurement_aliases(measurement_dict)
        measurement_dict.update(kwargs)
        return cls.from_dict(measurement_dict)

    def _normalize_measurement_aliases(self, measurement_dict):
        """Add common aliases for known techniques when payload lacks them."""
        aliases = measurement_dict.get("aliases") or {}
        if not isinstance(aliases, dict):
            aliases = {}

        normalized_aliases = {}
        for key, value in aliases.items():
            if isinstance(value, list):
                normalized_aliases[key] = value.copy()
            elif isinstance(value, tuple):
                normalized_aliases[key] = list(value)
            else:
                normalized_aliases[key] = [value]

        series_names = set(self._series_names_from_measurement_dict(measurement_dict))
        technique = str(measurement_dict.get("technique", "")).upper()

        alias_candidates = {
            "t": ["time/s", "time", "Time [s]", "elapsed time/s", "elapsed_time/s"],
        }
        if "EC" in technique:
            alias_candidates.update(
                {
                    "raw_potential": [
                        "Ewe/V",
                        "<Ewe>/V",
                        "E/V",
                        "potential",
                        "Voltage [V]",
                        "U/V",
                    ],
                    "raw_current": [
                        "I/mA",
                        "<I>/mA",
                        "I/A",
                        "current",
                        "Current/A",
                    ],
                }
            )

        for canonical_name, candidates in alias_candidates.items():
            if canonical_name in series_names:
                continue
            existing = normalized_aliases.get(canonical_name, [])
            matched = [name for name in candidates if name in series_names]
            if not matched:
                matched = self._heuristic_alias_match(
                    canonical_name=canonical_name,
                    series_names=series_names,
                )
            for name in matched:
                if name not in existing:
                    existing.append(name)
            if existing:
                normalized_aliases[canonical_name] = existing

        measurement_dict["aliases"] = normalized_aliases
        return measurement_dict

    @staticmethod
    def _series_names_from_measurement_dict(measurement_dict):
        series_names = []
        for series in measurement_dict.get("series_list", []):
            if hasattr(series, "name"):
                series_names.append(series.name)
            elif isinstance(series, dict) and "name" in series:
                series_names.append(series["name"])
        return series_names

    @staticmethod
    def _heuristic_alias_match(canonical_name, series_names):
        if canonical_name == "t":
            for name in series_names:
                lowered = name.lower()
                if "time" in lowered and ("/s" in lowered or "sec" in lowered):
                    return [name]
        if canonical_name == "raw_potential":
            for name in series_names:
                lowered = name.lower()
                if "ewe" in lowered or "potential" in lowered or "voltage" in lowered:
                    return [name]
        if canonical_name == "raw_current":
            for name in series_names:
                lowered = name.lower()
                if lowered.startswith("i/") or "current" in lowered:
                    return [name]
        return []

    def _looks_like_ixdat_measurement(self, payload):
        if not isinstance(payload, dict):
            return False
        has_technique = "technique" in payload
        has_series = "series_list" in payload
        return has_technique and has_series

    def _deserialize_measurement_dict(self, payload):
        payload_dict = payload.copy()

        list_casts = [
            ("series_list", DataSeries),
            ("component_measurements", Measurement),
            ("calculator_list", Calculator),
        ]
        for key, cls in list_casts:
            if key not in payload_dict:
                continue
            items = payload_dict[key]
            if not isinstance(items, list):
                continue
            casted_items = []
            for item in items:
                if isinstance(item, dict):
                    casted_items.append(cls.from_dict(item))
                else:
                    casted_items.append(item)
            payload_dict[key] = casted_items
        return payload_dict

    def _make_asimov_metadata(self, dataset, dataset_version):
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
            for dataset_version in versions:
                if dataset_version.get("id") == version_id:
                    return dataset_version
            raise ValueError(f"No dataset version found with id={version_id}.")

        if version is not None:
            for dataset_version in versions:
                if dataset_version.get("version") == version:
                    return dataset_version
            raise ValueError(f"No dataset version found with version={version}.")

        versions_sorted = sorted(
            versions,
            key=lambda dataset_version: dataset_version.get("created_at", ""),
            reverse=True,
        )
        return versions_sorted[0]

    def _build_auth_headers(self, force_login=False):
        if self._token:
            return {"Authorization": f"Bearer {self._token}"}
        if self.token_provider:
            token = self.token_provider.get_access_token(force_login=force_login)
            return {"Authorization": f"Bearer {token}"}
        raise RuntimeError(
            "AsimovReader has no authentication configured. Provide `token` "
            "or `token_provider`."
        )

    def _load_payload_uri(self, payload_uri, headers):
        if payload_uri.startswith("http://") or payload_uri.startswith("https://"):
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
