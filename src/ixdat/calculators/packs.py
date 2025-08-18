# CalculatorPack class
# @Author: Frederik Lizak Johansen
# @Created: 18/08/2025

from __future__ import annotations

import json
import hashlib
import importlib
from datetime import datetime, timezone
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union

from ixdat.db import Saveable
from ixdat.calculators.indexer import Indexer


JsonDict = Dict[str, Any]
Selector = Union[Iterable[str], Callable[[Any], bool], None]


def _now_utc_iso() -> str:
    """Return current UTC timestamp as ISO 8601 string with 'Z'."""
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _json_canonical(obj: Any) -> str:
    """Canonical JSON representation (sorted keys, no whitespace)."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def _sha256_of_jsonable(obj: Any) -> str:
    """SHA256 hex digest of a JSON-serializable object (canonicalized)."""
    s = _json_canonical(obj).encode("utf-8")
    return hashlib.sha256(s).hexdigest()


def _import_from_class_path(class_path: str) -> type:
    """Import and return a class given its 'module.ClassName' path."""
    mod_name, cls_name = class_path.rsplit(".", 1)
    module = importlib.import_module(mod_name)
    return getattr(module, cls_name)


def _matches_selector(obj: Any, selector: Selector) -> bool:
    """
    Return True if `obj` matches the selector.
    """
    if selector is None:
        return True

    if callable(selector):
        return bool(selector(obj))

    # Normalize selector to a set of strings
    if isinstance(selector, str):
        names = {selector}
    elif isinstance(selector, (set, list, tuple)):
        names = set(selector)
    else:
        # Fallback: try to coerce to set with one string representation
        names = {str(selector)}

    cls = obj.__class__
    class_name = cls.__name__
    class_path = f"{cls.__module__}.{class_name}"
    return (class_name in names) or (class_path in names)


def _ixdat_version() -> Optional[str]:
    """Best-effort retrieval of ixdat version string; returns None if unavailable."""
    try:
        from .. import __version__

        return __version__
    except Exception:
        return None


class CalculatorPack(Saveable):
    """
    A container that serializes multiple calculators together so they can be
    saved/loaded and re-applied to another Measurement in a single operation.

    Overview:
    - Captures a set of calculators (MS, EC, backgrounds, siq, ...)
      with minimal metadata and per-calculator payloads.
    - Supports:
        `CalculatorPack.from_measurement(...)` to capture calculators
        `pack.to_file(path)` and `CalculatorPack.read(path=...)` for file I/O
        `pack.calculators()` for lazy materialization
        `pack.attach_to(measurement, ...)` convenience to add all calculators

    Storage:
    - File format: JSON (`.ixpack.json` recommended)
    - Database: inherits Saveable; `json` column stores `to_json()` blob

    Integrity:
    - Each calculator payload is hashed (SHA256) for integrity verification

    """

    # for `Saveable`
    table_name = "calculator_pack"
    column_attrs = {"name", "schema_version", "ixdat_version", "created_utc", "json"}
    SCHEMA_VERSION = "1.0"  # Keep track of versioning for the pack

    def __init__(
        self,
        *,
        name: Optional[str] = None,
        calculators: Optional[Sequence[Any]] = None,
        metadata: Optional[JsonDict] = None,
    ):
        """
        Create an empty pack or with a given sequence of calculator objects.

        Parameters:
            name: Human-readable name for the pack.
            calculators: Sequence of calculator instances to include.
            If omitted, the pack starts empty and calculators can be
            loaded later via `from_json`.
            metadata: Arbitrary metadata stored under `metadata` in the JSON.
        """
        super().__init__()
        self.name: str = name or "CalculatorPack()"
        self.schema_version: str = self.SCHEMA_VERSION
        self.ixdat_version: Optional[str] = _ixdat_version()
        self.created_utc: str = (
            _now_utc_iso()
        )  # Centralized; can be used anywhere on earth
        self._metadata: JsonDict = metadata or {}

        # Internal storage:
        self._calculators: List[Any] = list(calculators) if calculators else []
        self._items: Optional[List[JsonDict]] = None  # filled by from_json or to_json

        # Backing JSON for Saveable (db column). Keep it None until first to_json()
        self.json: Optional[str] = None

    @staticmethod
    def _serialize_calculator(cal: Any) -> dict:
        """
        Serialize a calculator to a JSON-safe dict.

        Preference order:
          1) Calculator.as_dict()  (ixdat-native format used by Calculator.export())
          2) lightweight __dict__ snapshot (filtered)
        """
        try:
            from ..measurement_base import (
                Calculator as _BaseCalc,
            )  # local import to avoid import cycle errors

            if isinstance(cal, _BaseCalc) and hasattr(cal, "as_dict"):
                payload = cal.as_dict()
                # Assuming for sake of future implementation.
                # Ensure we don't persist a heavy/circular measurement reference:
                payload.pop("measurement", None)
                return {
                    "__format__": "ixdat.calculator.as_dict",
                    "payload": payload,
                    "calculator_type": getattr(cal, "calculator_type", None),
                    "technique": getattr(cal, "technique", None),
                }
        except Exception:
            pass

        # Fallback: filtered __dict__
        safe: dict[str, Any] = {"__format__": "light.dict"}
        for k, v in vars(cal).items():
            if k == "measurement":
                continue  # don't persist the measurement object
            if isinstance(v, (int, float, str, list, dict, type(None))):
                safe[k] = v
        return safe

    @staticmethod
    def _deserialize_calculator(cls: type, payload: dict) -> Any:
        """
        Recreate calculator instance from stored payload.

        Preference order:
          1) Calculator.from_dict() if we stored 'as_dict' format
          2) lightweight __dict__
        """
        fmt = payload.get("__format__")
        if fmt == "ixdat.calculator.as_dict":
            # Try the native method first
            calc_payload = dict(payload.get("payload") or {})
            # from_dict expects calculator metadata in the dict:
            if "calculator_type" not in calc_payload and payload.get("calculator_type"):
                calc_payload["calculator_type"] = payload["calculator_type"]
            if "technique" not in calc_payload and payload.get("technique"):
                calc_payload["technique"] = payload["technique"]
            try:
                from ..measurement_base import (
                    Calculator as _BaseCalc,
                )  # avoid circular import, by loading locally

                if issubclass(cls, _BaseCalc) and hasattr(cls, "from_dict"):
                    return cls.from_dict(calc_payload)
            except Exception:
                pass

        # Fallback: blank instance + shallow dict update
        cal = object.__new__(cls)
        try:
            # Strip helper key if present
            payload = dict(payload)
            payload.pop("__format__", None)
            cal.__dict__.update(payload)
        except Exception:
            pass
        return cal

    @classmethod
    def from_measurement(
        cls,
        measurement: Any,
        *,
        include: Selector = None,
        exclude: Selector = None,
        name: Optional[str] = None,
        notes: Optional[str] = None,
        extra_metadata: Optional[JsonDict] = None,
    ) -> "CalculatorPack":
        """
        Capture calculators from a Measurement into a pack.

        Parameters:
            measurement: Source Measurement providing `.calculators`.
            include: Selector to include some calculators.
            exclude: Selector to exclude some calculators.
            name: Optional name for the created pack.
            notes: Optional human-readable note stored in metadata.history.notes.
            extra_metadata: Extra metadata merged into the pack metadata.

        Returns:
            CalculatorPack
        """
        # Collect calculators
        if hasattr(measurement, "calculators"):
            if isinstance(measurement.calculators, dict):
                calculators = list(measurement.calculators.values())
            else:
                calculators = list(measurement.calculators)

        else:
            raise AttributeError("Measurement has no 'calculators' attribute")

        # Do not include Indexes
        calculators = [c for c in calculators if not isinstance(c, Indexer)]

        # Apply include/exclude filters
        if include is not None:
            calculators = [c for c in calculators if _matches_selector(c, include)]
        if exclude is not None:
            calculators = [c for c in calculators if not _matches_selector(c, exclude)]

        # Build metadata
        meta: JsonDict = {
            "history": {
                "source_measurement_uuid": getattr(measurement, "uuid", None),
                "source_measurement_name": getattr(measurement, "name", None),
                "notes": notes,
            }
        }
        if extra_metadata:
            meta.update(extra_metadata)

        return cls(
            name=name
            or f"CalculatorPack({getattr(measurement, 'name', 'measurement')})",
            calculators=calculators,
            metadata=meta,
        )

    def _ensure_items(self) -> List[JsonDict]:
        """
        Ensure internal `_items` list exists,
        synthesizing it from calculators if needed.
        """
        if self._items is not None:
            return self._items

        items: List[JsonDict] = []
        for cal in self._calculators:
            cls = cal.__class__
            class_path = f"{cls.__module__}.{cls.__name__}"

            payload = self._serialize_calculator(cal)

            item: JsonDict = {
                "class_path": class_path,
                "payload": payload,
                "hash": _sha256_of_jsonable(payload),
            }

            # Optional metadata
            cal_ver = getattr(cal, "version", None)
            if cal_ver:
                item["calculator_version"] = str(cal_ver)

            items.append(item)

        self._items = items
        return items

    def to_json(self) -> JsonDict:
        """
        Serialize the pack (including calculators) to a JSON-serializable dict.

        Returns:
            dict: The full, self-describing pack.
        """
        items = self._ensure_items()
        data: JsonDict = {
            "type": "ixdat.CalculatorPack",
            "name": self.name,
            "schema_version": self.schema_version,
            "ixdat_version": self.ixdat_version,
            "created_utc": self.created_utc,
            "metadata": self._metadata,
            "calculators": items,
        }
        # Also keep a string form for Saveable db column
        self.json = json.dumps(data, indent=2)
        return data

    @classmethod
    def from_json(cls, data: JsonDict) -> "CalculatorPack":
        """
        Construct a pack from a JSON dictionary (inverse of `to_json`).
        """
        if data.get("type") != "ixdat.CalculatorPack":
            raise ValueError(f"Unexpected type: {data.get('type')}")

        pack = cls(name=data.get("name"))
        pack.schema_version = str(data.get("schema_version") or cls.SCHEMA_VERSION)
        pack.ixdat_version = data.get("ixdat_version")
        pack.created_utc = data.get("created_utc") or _now_utc_iso()
        pack._metadata = data.get("metadata", {})
        pack._items = list(data.get("calculators") or [])
        # Keep canonical db json string if caller wants to save to DB
        pack.json = json.dumps(data, indent=2)
        return pack

    def to_file(self, path: str) -> None:
        """
        Write the pack to a JSON file.

        Parameters:
            path: Output file path (e.g., 'my_setup.ixpack.json').
        """
        data = self.to_json()
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)

    @classmethod
    def read(
        cls,
        path: str,
    ) -> "CalculatorPack":
        """
        Read a pack from a file (`path`)

        Returns:
            CalculatorPack
        """
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return cls.from_json(data)

    def calculators(self) -> List[Any]:
        """
        Materialize and return the list of calculator objects contained in the pack.

        Returns:
            list: List of calculator instances reconstructed from their stored payloads.
        """
        if self._calculators:
            return list(self._calculators)

        if not self._items:
            return []

        calculators: List[Any] = []
        for item in self._items:
            class_path = item["class_path"]
            payload = item["payload"]

            # Integrity check (best-effort)
            h_expected = item.get("hash")
            h_actual = _sha256_of_jsonable(payload)
            if h_expected and h_expected != h_actual:
                raise ValueError(
                    f"Integrity check failed for {class_path}: hash mismatch"
                )

            CalCls = _import_from_class_path(class_path)

            # Always use CalculatorPack’s deserializer
            cal = self._deserialize_calculator(CalCls, payload)
            calculators.append(cal)

        self._calculators = calculators

        return list(self._calculators)

    def validate(self, *, strict: bool = True) -> Tuple[bool, List[str]]:
        """
        Validate the internal structure and per-calculator integrity hashes.

        Parameters:
            strict: If True, treat any discrepancy as an error; otherwise aggregate warnings.

        Returns:
            (ok, messages) ok=True if everything looks good (or if not strict and only warnings).
        """
        msgs: List[str] = []
        ok = True

        if self.schema_version != self.SCHEMA_VERSION:
            ok = False
            msgs.append(
                f"Schema version mismatch: pack={self.schema_version}, "
                f"supported={self.SCHEMA_VERSION}"
            )

        if not self._items:
            msgs.append("Pack contains no calculator entries")
            return (ok and not strict, msgs)

        for item in self._items:
            cp = item.get("class_path")
            payload = item.get("payload")
            if not cp or not isinstance(cp, str):
                ok = False
                msgs.append("Entry missing valid 'class_path'")
                if strict:
                    continue
            if not isinstance(payload, dict):
                ok = False
                msgs.append(f"{cp}: payload is not a dict")
                if strict:
                    continue
            h_expected = item.get("hash")
            if h_expected:
                h_actual = _sha256_of_jsonable(payload)
                if h_expected != h_actual:
                    ok = False
                    msgs.append(f"{cp}: integrity hash mismatch")

        if strict and not ok:
            return (False, msgs)
        return (ok, msgs)

    def attach_to(
        self,
        measurement: Any,
        *,
        on_conflict: str = "replace",
    ) -> List[Any]:
        """
        Attach all calculators in the pack to a Measurement.

        Parameters:
            measurement: Target Measurement (must provide `.add_calculator(calculator)`).
            on_conflict:
                "replace" (default): replace existing calculator of same class+name
                "skip": keep existing one, skip new
                "duplicate": add new one and allow duplicates (caller handles names)
        Returns:
            list: The list of calculators actually attached to the measurement.
        """
        attached: List[Any] = []
        calculators = self.calculators()

        # Acquire current calculators on measurement (for conflict handling)
        existing: List[Any]
        if hasattr(measurement, "calculators"):
            if isinstance(measurement.calculators, dict):
                existing = list(measurement.calculators.values())
            else:
                existing = list(measurement.calculators)
        else:
            existing = []

        def _same_kind(a: Any, b: Any) -> bool:
            return a.__class__ is b.__class__ and getattr(a, "name", None) == getattr(
                b, "name", None
            )

        for cal in calculators:
            # Conflict resolution
            conflict = next((c for c in existing if _same_kind(c, cal)), None)
            if conflict:
                if on_conflict == "replace":
                    # remove then add (if API supports remove)
                    if hasattr(measurement, "remove_calculator"):
                        try:
                            measurement.remove_calculator(conflict)
                        except Exception:
                            pass
                    # fall through to add
                elif on_conflict == "skip":
                    continue
                elif on_conflict == "duplicate":
                    # Let both exist; optionally rename to avoid ambiguity
                    if getattr(cal, "name", None) and hasattr(cal, "name"):
                        cal.name = f"{cal.name} (2)"
                else:
                    raise ValueError(f"Unknown on_conflict policy: {on_conflict}")

            # Attach
            if not hasattr(measurement, "add_calculator") or not callable(
                measurement.add_calculator
            ):
                raise AttributeError("Measurement lacks add_calculator(calculator)")

            measurement.add_calculator(cal)
            attached.append(cal)
            existing.append(cal)

        return attached

    def summary(self) -> str:
        """Human-readable one-line-per-entry summary of the pack contents."""
        parts: List[str] = [
            f"CalculatorPack(name={self.name!r}, schema={self.schema_version},"
            f"ixdat={self.ixdat_version}, created={self.created_utc})"
        ]
        items = self._items or []
        if not items and self._calculators:
            # Synthesize items to summarize
            _ = self._ensure_items()
            items = self._items or []

        for i, item in enumerate(items, 1):
            cp = item.get("class_path", "?")
            ver = item.get("calculator_version")
            suffix = []
            if ver:
                suffix.append(f"ver={ver}")
            parts.append(
                f"  {i:02d}. {cp}" + (f" ({', '.join(suffix)})" if suffix else "")
            )
        return "\n".join(parts)

    def __str__(self) -> str:
        return self.summary()

    def __repr__(self) -> str:
        repr_name = (
            f"<CalculatorPack name={self.name!r} "
            f"items={len(self._items or self._calculators)}>"
        )
        return repr_name
