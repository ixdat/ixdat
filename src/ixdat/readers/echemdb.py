# -*- coding: utf-8 -*-
import requests
import zipfile
import io
import os
import json
from pathlib import Path
from typing import NamedTuple, TypeVar, Type, Optional, Union, List

import pandas as pd

from ixdat.measurement_base import Measurement
from ixdat.techniques.cv import CyclicVoltammogram
from ixdat.data_series import TimeSeries, ValueSeries, ConstantValue
from ixdat.tools import (
    get_default_cache_dir,
    request_with_retries,
)
from ixdat.exceptions import BuildError


M = TypeVar("M", bound=Measurement)


class EChemDBReader:
    """
    Reader for datasets hosted on EChemDB (via their github release zip).

    Usage:
        reader = EChemDBReader() # uses latest release
        reader = EChemDBReader(version='0.4.1') # pins to v0.4.1
        meas = Measurement.read("alves_2011_electrochemistry_6010", reader='echemdb')

    Notes:
        All network calls use ixdat.tools.request_with_retries.
        This enforces connect, read, and total timeouts,
        allows Ctrl-C to abort immediately,
        and raises clear RuntimeError messages if retries are exhausted.
    """

    # fallback aliases are disabled to avoid circular lookups.
    # Keep empty and handle explicit fallbacks where needed.
    ALIASES = {}

    class _Paths(NamedTuple):
        csv: Path
        json: Path
        bib: Optional[Path]

    def __init__(
        self,
        version: str = "latest",
        connect_timeout: float = 3.0,
        read_timeout: float = 10.0,
        total_timeout: Optional[float] = None,
    ):
        self._version_arg = version
        self._connect_timeout = connect_timeout
        self._read_timeout = read_timeout
        self._total_timeout = total_timeout
        self._session = self._build_session()

    def read(
        self,
        echemdb_identifier: str,
        cls: Type[M] = Measurement,
        version: Optional[str] = None,
        connect_timeout: float = 1.0,
        read_timeout: float = 1.0,
        total_timeout: float = 20.0,
    ) -> M:
        """
        Download (if not in cache), extract, parse CSV+JSON for each measurement

        Args:
            echemdb_identifier: the name on https://www.echemdb.org/
            cls: class to return. If cls is Measurement, return a
                CyclicVoltammogram object
            version: optional override of which release version to use

        Returns an object of type cls
        """

        # allow version override, and resolve "latest" via PyPI
        ver = self._resolve_version(
            version or self._version_arg,
            connect_timeout or self._connect_timeout,
            read_timeout or self._read_timeout,
            total_timeout or self._total_timeout,
        )

        self.zip_url = (
            f"https://github.com/echemdb/electrochemistry-data/"
            f"releases/download/{ver}/data-{ver}.zip"
        )
        self.cache_root = get_default_cache_dir("echemdb") / f"data-{ver}"

        if cls is Measurement:
            cls = CyclicVoltammogram  # type: ignore[arg-type]

        # Extract only this entry if not already cached
        paths = self._find_cached_paths(echemdb_identifier)
        if not paths:
            self._extract_entry(
                echemdb_identifier,
                connect_timeout or self._connect_timeout,
                read_timeout or self._read_timeout,
                total_timeout or self._total_timeout,
            )
            paths = self._find_cached_paths(echemdb_identifier)
            if not paths:
                raise FileNotFoundError(
                    f"Could not find files for entry"
                    f"{echemdb_identifier} in version {ver} @ {self.zip_url}"
                )

        df = self._load_dataframe(paths.csv)
        meta = self._load_metadata(paths.json)
        citation = self._load_citation(paths.bib)
        series_list = self._make_series_list(df, meta)
        measurement = self._assemble_measurement(
            echemdb_identifier, cls, series_list, meta, citation
        )

        return measurement

    def _build_session(self) -> requests.Session:
        s = requests.Session()
        s.trust_env = False  # Ignore HTTP(S)_PROXY for now.
        return s

    def _resolve_version(
        self,
        ver: str,
        read_timeout: float,
        connect_timeout: float,
        total_timeout: Optional[float],
    ) -> str:
        if ver != "latest":
            return ver
        # Fetch with retries and explicit timeouts
        # (raises RuntimeError if retries exhausted: no internet/DNS/TCP issues)
        resp = request_with_retries(
            self._session,
            "GET",
            "https://pypi.org/pypi/echemdb-ecdata/json",
            read_timeout=read_timeout,
            connect_timeout=connect_timeout,
            total_timeout=total_timeout,
            retries=3,
        )
        ver = resp.json()["info"]["version"]
        print(f"[EChemDBReader] using latest release: {ver}")
        return resp.json()["info"]["version"]

    def _find_cached_paths(self, echemdb_identifier: str) -> Optional[_Paths]:
        """
        Search the cache for existing CSV, JSON,
        and optional bib files matching echemdb_identifier.

        Returns:
            _Paths or None if not found.
        """
        root = self.cache_root / "data/generated/svgdigitizer"
        for folder in root.glob("*"):
            if not folder.is_dir():
                continue
            csv_path = folder / f"{echemdb_identifier}.csv"
            json_path = folder / f"{echemdb_identifier}.json"
            if csv_path.exists() and json_path.exists():
                bib_files = list(folder.glob("*.bib"))
                bib_path = bib_files[0] if bib_files else None
                return self._Paths(csv=csv_path, json=json_path, bib=bib_path)

    def _extract_entry(
        self,
        echemdb_identifier: str,
        read_timeout: float,
        connect_timeout: float,
        total_timeout: Optional[float],
    ) -> None:
        """
        Download the release ZIP and extract only the files under the given identifier
        """
        # Fetch with retries and explicit timeouts
        # (raises RuntimeError if retries exhausted: no internet/DNS/TCP issues)
        resp = request_with_retries(
            self._session,
            "GET",
            self.zip_url,
            read_timeout=read_timeout,
            connect_timeout=connect_timeout,
            total_timeout=total_timeout,
            retries=0,
            stream=True,
        )

        with zipfile.ZipFile(io.BytesIO(resp.content)) as z:
            namelist = z.namelist()

            # find the single CSV member
            csv_members = [
                m for m in namelist if m.endswith(f"{echemdb_identifier}.csv")
            ]
            if len(csv_members) != 1:
                raise FileNotFoundError(
                    f"Expected exactly one CSV"
                    f"for '{echemdb_identifier}', found {csv_members}"
                )
            csv_member = csv_members[0]

            # derive the JSON member name
            json_member = csv_member[:-4] + ".json"
            if json_member not in namelist:
                raise FileNotFoundError(f"Missing JSON {json_member} in ZIP")

            # find any .bib in the same directory
            prefix = os.path.dirname(csv_member) + "/"
            bib_members = [
                m
                for m in namelist
                if m.startswith(prefix) and m.lower().endswith(".bib")
            ]
            if not bib_members:
                bib_member = None
            elif len(bib_members) == 1:
                bib_member = bib_members[0]
            else:
                # disambiguate or error
                exact = [
                    m for m in bib_members if m.endswith(f"{echemdb_identifier}.bib")
                ]
                if len(exact) == 1:
                    bib_member = exact[0]
                else:
                    raise FileExistsError(
                        f"Multiple .bib files for '{echemdb_identifier}': {bib_members}"
                    )

            # extract just those files
            to_extract = [csv_member, json_member]
            if bib_member:
                to_extract.append(bib_member)

            for member in to_extract:
                target = self.cache_root / member
                target.parent.mkdir(parents=True, exist_ok=True)
                with z.open(member) as src, open(target, "wb") as dst:
                    dst.write(src.read())

    def _load_dataframe(self, path: Path) -> pd.DataFrame:
        """Read the CSV file into a pandas DataFrame"""
        return pd.read_csv(path)

    def _load_metadata(self, path: Path) -> dict:
        """Read the JSON file into a dict"""
        with path.open("r", encoding="utf-8") as f:
            return json.load(f)

    def _load_citation(self, path: Optional[Path]) -> Optional[str]:
        """Read the BibTex citation, if it exists"""
        return path.read_text(encoding="utf-8") if path is not None else None

    def _make_series_list(
        self, df: pd.DataFrame, full_meta: dict
    ) -> List[Union[TimeSeries, ValueSeries]]:
        """
        Build DataSeries from the DataFrame using the JSON schema.

        Args:
            df (pd.DataFrame): the CSV-loaded data
            full_meta (dict): the entire JSON metadata for this entry

        Returns:
            List of TimeSeries and ValueSeries in schema order.
        """
        # grab the single table resource and its schema
        resource = full_meta["resources"][0]
        schema_fields = resource["schema"]["fields"]

        # determine t=0 timestamp (if present at top level) or default to 0
        tstamp = full_meta.get("tstamp", 0.0)

        series_list: List[Union[TimeSeries, ValueSeries]] = []
        tseries = None

        # iterate the schema in order
        for idx, field in enumerate(schema_fields):
            col_name = field["name"]
            unit = field.get("unit", "")

            # skip columns missing from the CSV
            if col_name not in df.columns:
                continue

            data = df[col_name].to_numpy()

            if idx == 0:
                # first field is TimeSeries
                tseries = TimeSeries(
                    name=col_name,
                    unit_name=unit,
                    data=data,
                    tstamp=tstamp,
                )
                series_list.append(tseries)
            else:
                # subsequent fields are ValueSeries linked to tseries
                series_list.append(
                    ValueSeries(
                        name=col_name,
                        unit_name=unit,
                        data=data,
                        tseries=tseries,
                    )
                )

        # append cycle (singular in the case of EChemDB)
        series_list.append(
            ConstantValue(name="cycle", unit_name="", data=0, tseries=tseries)
        )

        if tseries is None:
            raise BuildError("No time column found according to JSON schema!")

        return series_list

    def _assemble_measurement(
        self,
        echemdb_identifier: str,
        cls: Type[M],
        series_list: list,
        meta: dict,
        citation: Optional[str],
    ) -> M:
        """
        Assemble the cls.from_dict payload and return the instance.
        Attaches citation as `.citation` attribute if provided.
        """

        def get_meta_value(d, *keys, default=None):
            for key in keys:
                value = d.get(key)
                if value is not None:
                    return value
            return default

        resource = meta["resources"][0]
        schema_fields = resource["schema"]["fields"]
        echem_meta = resource["metadata"]["echemdb"]

        # experimental tags
        tags = get_meta_value(echem_meta.get("experimental", {}), "tags")

        # electrode/electrolyte/system info
        system = get_meta_value(echem_meta, "system")
        electrodes = get_meta_value(system or {}, "electrodes")
        electrolyte = get_meta_value(system or {}, "electrolyte")

        # source / citation info
        source = get_meta_value(echem_meta, "source")
        citation_key = get_meta_value(source or {}, "citation key", "citationKey")
        url = get_meta_value(source or {}, "url")
        figure = get_meta_value(source or {}, "figure")
        curve = get_meta_value(source or {}, "curve")
        bibdata = get_meta_value(source or {}, "bibdata")

        # figure description, measurement type, scan rate
        fig_desc = get_meta_value(echem_meta, "figure description", "figureDescription")
        meas_type = get_meta_value(fig_desc or {}, "measurement type", "measurementType")
        scan_rate = get_meta_value(fig_desc or {}, "scan rate", "scanRate")

        # build the alias map ixdat’s CV code expects
        scan_fields = get_meta_value(fig_desc or {}, "fields", default=[])
        horiz = None
        vert = None
        for fld in scan_fields:
            name = fld["name"]
            orient = fld.get("orientation", "").lower()
            if orient == "horizontal":
                horiz = name
            elif orient == "vertical":
                vert = name

        # get reference electrode if available
        reference_electrode = scan_fields[0].get(
            "reference", "[unknown reference electrode]"
        )

        # make aliases with metadata-derived values if available
        aliases = {}
        if horiz or vert:
            if horiz:
                aliases["potential"] = [horiz]
                aliases["raw_potential"] = [horiz]
            if vert:
                aliases["current"] = [vert]
                aliases["raw_current"] = [vert]
        else:
            # fallback only if no orientation info is provided at all
            aliases = {
                "potential": ["E"],
                "current": ["j"],
                "raw_potential": ["E"],
                "raw_current": ["j"],
            }

        # assemble metadata dict
        md: dict = {
            "entry_id": echemdb_identifier,
            "tags": tags,
            "system": system,
            "electrodes": electrodes,
            "electrolyte": electrolyte,
            "citation_key": citation_key,
            "url": url,
            "figure": figure,
            "curve": curve,
            "measurement_type": meas_type,
            "scan_rate": scan_rate,
            "reference_electrode": reference_electrode,
        }

        # preserve the raw schema-fields mapping
        md["schema_fields"] = {fld["name"]: fld for fld in schema_fields}

        # expose full metadata tree in case it is useful for the curious user
        md["_raw_echemdb_meta"] = echem_meta

        # build the payload
        payload = {
            "name": echemdb_identifier,
            "technique": "EC",
            "reader": self,
            "series_list": series_list,
            "tstamp": md["schema_fields"].get("t", {}).get("unit") and 0.0,
            "metadata": md,
            "aliases": aliases,
        }

        meas = cls.from_dict(payload)
        if citation:
            setattr(meas, "citation", citation or bibdata)

        return meas
