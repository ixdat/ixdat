"""Inspect current Asimov payload structure without printing full data arrays.

Usage:
    python inspect_asimov_payloads.py ID [ID ...]

Set ASIMOV_ACCESS_TOKEN for non-interactive use, or let AsimovReader use the
normal Keycloak login flow.
"""

import argparse
from collections import Counter

from ixdat.readers.asimov import AsimovReader


def _shape(value):
    """Return the nested list shape if value looks like array data."""
    shape = []
    cursor = value
    while isinstance(cursor, list):
        shape.append(len(cursor))
        cursor = cursor[0] if cursor else None
    return tuple(shape)


def _short(value):
    """Return a compact scalar/list summary suitable for terminal inspection."""
    if isinstance(value, list):
        return f"list shape={_shape(value)}"
    if isinstance(value, dict):
        return f"dict keys={sorted(value)}"
    return repr(value)


def _print_series(entry, indent):
    prefix = " " * indent
    print(
        f"{prefix}- {entry.get('series_type', 'series')}: "
        f"{entry.get('name')!r} [{entry.get('unit_name')!r}]"
    )
    if "key" in entry:
        print(f"{prefix}  key={entry['key']!r}")
    if "data" in entry:
        print(f"{prefix}  data={_short(entry['data'])}")
    for ref_key in ("tseries_key", "axes_keys"):
        if ref_key in entry:
            print(f"{prefix}  {ref_key}={entry[ref_key]!r}")
    if entry.get("axes_series"):
        print(f"{prefix}  axes_series:")
        for axis in entry["axes_series"]:
            _print_series(axis, indent + 4)


def _print_object(dct, indent=0, label="payload"):
    prefix = " " * indent
    object_type = dct.get("object_type", "<missing>")
    technique = dct.get("technique", "<missing>")
    print(f"{prefix}{label}: object_type={object_type!r}, technique={technique!r}")
    for key in ("name", "sample_name", "tstamp", "duration", "durations", "continuous"):
        if key in dct:
            print(f"{prefix}  {key}={_short(dct[key])}")
    if isinstance(dct.get("metadata"), dict):
        print(f"{prefix}  metadata keys={sorted(dct['metadata'])}")
    if dct.get("series_list"):
        counts = Counter(s.get("series_type", "series") for s in dct["series_list"])
        print(
            f"{prefix}  series_list len={len(dct['series_list'])} types={dict(counts)}"
        )
        for series in dct["series_list"][:6]:
            _print_series(series, indent + 4)
        if len(dct["series_list"]) > 6:
            print(f"{prefix}    ... {len(dct['series_list']) - 6} more series")
    if isinstance(dct.get("field"), dict):
        print(f"{prefix}  field:")
        _print_series(dct["field"], indent + 4)
    nested = [
        (key, value)
        for key, value in dct.items()
        if isinstance(value, dict) and "object_type" in value
    ]
    for key, value in nested:
        _print_object(value, indent + 2, label=key)


def _get_payload(reader, asimov_id, force_login=False):
    headers = reader._build_auth_headers(force_login=force_login)
    envelope = reader._get_ixdat_payload_envelope(asimov_id, headers=headers)
    payload = envelope.get("payload_json")
    if not payload and envelope.get("payload_uri"):
        payload = reader._load_payload_uri(envelope["payload_uri"], headers=headers)
    return envelope, payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("asimov_ids", nargs="+")
    parser.add_argument("--force-login", action="store_true")
    args = parser.parse_args()

    reader = AsimovReader()
    for asimov_id in args.asimov_ids:
        envelope, payload = _get_payload(
            reader,
            asimov_id,
            force_login=args.force_login,
        )
        print("=" * 88)
        print(
            f"project_file_id={envelope.get('project_file_id')} "
            f"bundle_id={envelope.get('bundle_id')} "
            f"filename={envelope.get('filename')!r} "
            f"name={envelope.get('name')!r}"
        )
        print(
            f"created_at={envelope.get('created_at')} "
            f"parser={envelope.get('parser_name')!r} "
            f"ixdat={envelope.get('ixdat_version')!r}"
        )
        if payload:
            _print_object(payload)
        else:
            print("payload: <missing>")


if __name__ == "__main__":
    main()
