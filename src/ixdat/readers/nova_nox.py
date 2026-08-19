"""ixdat reader for Metrohm Autolab NOVA's native binary ".nox" procedure
files -- fully self-contained (no external nox_reader project dependency,
no .NET runtime, NOVA installation, or Autolab SDK needed).

    from ixdat.techniques.ec import ECMeasurement
    ec = ECMeasurement.read("Cyclic voltammetry potentiostatic(3).nox", reader="nova_nox")

    # or directly:
    from ixdat.readers.nova_nox import NovaNoxReader
    ec = NovaNoxReader().read("Cyclic voltammetry potentiostatic(3).nox")

Unlike `NovaASCIIReader` (`autolab.py`), which reads NOVA's manually
exported ";"-delimited ASCII CSV, this reads the ".nox" file NOVA itself
saves directly -- no export step needed -- and it recovers a real
per-datapoint wall-clock timestamp from the procedure's own serialized
start time (`Procedure.start_time`, decoded from the embedded .NET
DateTime), not just an elapsed-seconds axis with a prompted/assumed t=0.

This module has three layers, in one file for easy distribution:

  1. A generic parser for MS-NRBF (.NET Remoting Binary Format, what
     System.Runtime.Serialization.Formatters.Binary.BinaryFormatter
     produces) -- NOVA saves ".nox" files with it, so this lets us
     deserialize them without needing the original .NET assemblies: class
     field names and array shapes are encoded inline in the stream itself.
     Spec reference: [MS-NRBF] .NET Remoting: Binary Format Data Structure.
     Only the subset of record types actually needed to read .nox files is
     implemented; unsupported records raise NotImplementedError with the
     byte offset so gaps are easy to find and patch.
  2. NOVA/.nox-specific object-graph extraction on top of that: walking the
     deserialized graph for calculated data channels (`load_nox` and
     friends below), grouping them into per-technique "curves" (pandas
     DataFrames), and decoding procedure metadata (name, start time, ...).
  3. `NovaNoxReader`, the actual ixdat reader class, which turns a parsed
     Procedure/Curve into an ECMeasurement (or a `cls=` subclass, e.g.
     CyclicVoltammogram).

Limitations (inherited from the underlying parser):
  - Only *calculated* channels (plain `List<Double>` values) are decoded;
    raw (uncalibrated) ADC buffer channels are skipped, since decoding them
    needs internal gain/offset calibration tables not exposed at the
    object-graph layer.
  - A .nox file with multiple curves (e.g. several sequential CV scan
    segments run back-to-back in one procedure) is concatenated into one
    continuous measurement by default (`concatenate=True`), assuming they
    share the same columns and a continuous elapsed-time axis -- true for
    sequential segments of one technique. Pass `concatenate=False` and
    `curve_index=` to get a single segment instead.
  - Potential is only exposed as `raw_potential` if NOVA recorded it as a
    per-point channel. Techniques that hold potential at a fixed setpoint
    instead (e.g. some chronoamperometry procedures) won't have one -- use
    `cls=ECMeasurement` (the default), not `CyclicVoltammogram`, for those.
    Going through `ECMeasurement.read(...)` (rather than instantiating
    `NovaNoxReader` directly) enforces `raw_potential` as an "essential
    series" and will raise `SeriesNotFoundError` for such files.
  - If a .nox file's object graph uses a record type not seen in NOVA's CV
    and chronoamperometry exports (the two techniques this was developed
    and tested against), the NRBF layer raises NotImplementedError with the
    byte offset -- straightforward to extend against the [MS-NRBF] spec.
"""
from __future__ import annotations

import struct
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from ..data_series import TimeSeries, ValueSeries


# =============================================================================
# 1. Generic MS-NRBF (.NET Remoting Binary Format) parser
# =============================================================================

# ---- record type / enum constants -----------------------------------------

RT_SERIALIZED_STREAM_HEADER = 0
RT_CLASS_WITH_ID = 1
RT_SYSTEM_CLASS_WITH_MEMBERS = 2
RT_CLASS_WITH_MEMBERS = 3
RT_SYSTEM_CLASS_WITH_MEMBERS_AND_TYPES = 4
RT_CLASS_WITH_MEMBERS_AND_TYPES = 5
RT_BINARY_OBJECT_STRING = 6
RT_BINARY_ARRAY = 7
RT_MEMBER_PRIMITIVE_TYPED = 8
RT_MEMBER_REFERENCE = 9
RT_OBJECT_NULL = 10
RT_MESSAGE_END = 11
RT_BINARY_LIBRARY = 12
RT_OBJECT_NULL_MULTIPLE_256 = 13
RT_OBJECT_NULL_MULTIPLE = 14
RT_ARRAY_SINGLE_PRIMITIVE = 15
RT_ARRAY_SINGLE_OBJECT = 16
RT_ARRAY_SINGLE_STRING = 17
RT_METHOD_CALL = 21
RT_METHOD_RETURN = 22

BT_PRIMITIVE = 0
BT_STRING = 1
BT_OBJECT = 2
BT_SYSTEM_CLASS = 3
BT_CLASS = 4
BT_OBJECT_ARRAY = 5
BT_STRING_ARRAY = 6
BT_PRIMITIVE_ARRAY = 7

PT_BOOLEAN = 1
PT_BYTE = 2
PT_CHAR = 3
PT_DECIMAL = 5
PT_DOUBLE = 6
PT_INT16 = 7
PT_INT32 = 8
PT_INT64 = 9
PT_SBYTE = 10
PT_SINGLE = 11
PT_TIMESPAN = 12
PT_DATETIME = 13
PT_UINT16 = 14
PT_UINT32 = 15
PT_UINT64 = 16
PT_NULL = 17
PT_STRING = 18

_PRIMITIVE_STRUCT = {
    PT_BOOLEAN: "<?",
    PT_BYTE: "<B",
    PT_SBYTE: "<b",
    PT_CHAR: None,  # handled specially (UTF-8 char, variable width)
    PT_DOUBLE: "<d",
    PT_INT16: "<h",
    PT_INT32: "<i",
    PT_INT64: "<q",
    PT_SINGLE: "<f",
    PT_UINT16: "<H",
    PT_UINT32: "<I",
    PT_UINT64: "<Q",
    PT_TIMESPAN: "<q",
    PT_DATETIME: "<Q",
}

PRIMITIVE_NAMES = {
    PT_BOOLEAN: "Boolean", PT_BYTE: "Byte", PT_CHAR: "Char",
    PT_DECIMAL: "Decimal", PT_DOUBLE: "Double", PT_INT16: "Int16",
    PT_INT32: "Int32", PT_INT64: "Int64", PT_SBYTE: "SByte",
    PT_SINGLE: "Single", PT_TIMESPAN: "TimeSpan", PT_DATETIME: "DateTime",
    PT_UINT16: "UInt16", PT_UINT32: "UInt32", PT_UINT64: "UInt64",
    PT_NULL: "Null", PT_STRING: "String",
}


class NrbfObject:
    """A deserialized .NET object: class name + ordered member dict."""

    __slots__ = ("class_name", "members", "object_id")

    def __init__(self, class_name: str, object_id: int):
        self.class_name = class_name
        self.object_id = object_id
        self.members: dict[str, Any] = {}

    def __repr__(self):
        return f"<NrbfObject {self.class_name} id={self.object_id} keys={list(self.members)}>"

    def __getitem__(self, key):
        return self.members[key]

    def get(self, key, default=None):
        return self.members.get(key, default)


class NrbfArray(list):
    """A deserialized .NET array/list; subclasses list, tags element type."""

    def __init__(self, *args, element_type=None, object_id=None):
        super().__init__(*args)
        self.element_type = element_type
        self.object_id = object_id


@dataclass
class ClassInfo:
    object_id: int
    name: str
    member_names: list[str]
    # parallel to member_names once type info is read:
    binary_types: list[int] = field(default_factory=list)
    additional_infos: list[Any] = field(default_factory=list)


class NrbfReader:
    def __init__(self, data: bytes):
        self.data = data
        self.pos = 0
        self.objects: dict[int, Any] = {}
        self.libraries: dict[int, str] = {}
        self.class_templates: dict[int, ClassInfo] = {}  # by object_id
        self.root_id = None
        self.header_id = None

    # -- primitive readers ----------------------------------------------
    def _read(self, n: int) -> bytes:
        b = self.data[self.pos:self.pos + n]
        if len(b) != n:
            raise EOFError(f"Unexpected EOF at offset {self.pos}, wanted {n} bytes")
        self.pos += n
        return b

    def read_byte(self) -> int:
        return self._read(1)[0]

    def read_int32(self) -> int:
        return struct.unpack_from("<i", self._read(4))[0]

    def read_7bit_length(self) -> int:
        result = 0
        shift = 0
        while True:
            b = self.read_byte()
            result |= (b & 0x7F) << shift
            if not (b & 0x80):
                break
            shift += 7
        return result

    def read_string(self) -> str:
        length = self.read_7bit_length()
        raw = self._read(length)
        return raw.decode("utf-8")

    def read_primitive(self, ptype: int):
        if ptype == PT_STRING:
            return self.read_string()
        if ptype == PT_DECIMAL:
            # encoded as a length-prefixed ASCII string
            return self.read_string()
        if ptype == PT_CHAR:
            # UTF-8 encoded char, 1-4 bytes; peek to find width
            first = self.data[self.pos]
            if first < 0x80:
                n = 1
            elif first < 0xE0:
                n = 2
            elif first < 0xF0:
                n = 3
            else:
                n = 4
            return self._read(n).decode("utf-8")
        fmt = _PRIMITIVE_STRUCT.get(ptype)
        if fmt is None:
            raise NotImplementedError(f"Unsupported primitive type {ptype} at {self.pos}")
        size = struct.calcsize(fmt)
        return struct.unpack_from(fmt, self._read(size))[0]

    # -- class info / member type info -----------------------------------
    def read_class_info(self) -> ClassInfo:
        object_id = self.read_int32()
        name = self.read_string()
        count = self.read_int32()
        member_names = [self.read_string() for _ in range(count)]
        return ClassInfo(object_id, name, member_names)

    def read_member_type_info(self, ci: ClassInfo):
        n = len(ci.member_names)
        ci.binary_types = [self.read_byte() for _ in range(n)]
        infos = []
        for bt in ci.binary_types:
            if bt in (BT_PRIMITIVE, BT_PRIMITIVE_ARRAY):
                infos.append(self.read_byte())  # PrimitiveTypeEnum
            elif bt == BT_SYSTEM_CLASS:
                infos.append(self.read_string())
            elif bt == BT_CLASS:
                type_name = self.read_string()
                library_id = self.read_int32()
                infos.append((type_name, library_id))
            else:
                # String, Object, StringArray, ObjectArray: no extra info
                infos.append(None)
        ci.additional_infos = infos

    # -- top level ---------------------------------------------------------
    def parse(self):
        rt = self.read_byte()
        if rt != RT_SERIALIZED_STREAM_HEADER:
            raise ValueError(f"Expected header record, got type {rt}")
        self.root_id = self.read_int32()
        self.header_id = self.read_int32()
        major = self.read_int32()
        minor = self.read_int32()
        if (major, minor) != (1, 0):
            raise ValueError(f"Unsupported NRBF version {major}.{minor}")

        while True:
            rt = self.read_byte()
            if rt == RT_MESSAGE_END:
                break
            self._dispatch_top(rt)

        return self.objects.get(self.root_id)

    def _dispatch_top(self, rt: int):
        # Top-level records that don't themselves produce a "value" we need
        # to return anywhere (library refs), plus the same object records
        # that can appear nested. We just reuse read_value_record.
        self._read_record(rt)

    # -- generic record reader ---------------------------------------------
    def _read_record(self, rt: int):
        if rt == RT_BINARY_LIBRARY:
            lib_id = self.read_int32()
            name = self.read_string()
            self.libraries[lib_id] = name
            return None
        if rt == RT_CLASS_WITH_MEMBERS_AND_TYPES:
            return self._read_class_with_members_and_types(system=False)
        if rt == RT_SYSTEM_CLASS_WITH_MEMBERS_AND_TYPES:
            return self._read_class_with_members_and_types(system=True)
        if rt == RT_CLASS_WITH_MEMBERS:
            return self._read_class_with_members(system=False)
        if rt == RT_SYSTEM_CLASS_WITH_MEMBERS:
            return self._read_class_with_members(system=True)
        if rt == RT_CLASS_WITH_ID:
            return self._read_class_with_id()
        if rt == RT_BINARY_OBJECT_STRING:
            return self._read_binary_object_string()
        if rt == RT_MEMBER_REFERENCE:
            idref = self.read_int32()
            return _Ref(idref)
        if rt == RT_OBJECT_NULL:
            return None
        if rt == RT_OBJECT_NULL_MULTIPLE_256:
            count = self.read_byte()
            return _NullRun(count)
        if rt == RT_OBJECT_NULL_MULTIPLE:
            count = self.read_int32()
            return _NullRun(count)
        if rt == RT_ARRAY_SINGLE_PRIMITIVE:
            return self._read_array_single_primitive()
        if rt == RT_ARRAY_SINGLE_OBJECT:
            return self._read_array_single_object()
        if rt == RT_ARRAY_SINGLE_STRING:
            return self._read_array_single_string()
        if rt == RT_BINARY_ARRAY:
            return self._read_binary_array()
        if rt == RT_MEMBER_PRIMITIVE_TYPED:
            ptype = self.read_byte()
            return self.read_primitive(ptype)
        raise NotImplementedError(f"Unsupported record type {rt} at offset {self.pos}")

    def _read_value_member(self, bt: int, add_info):
        """Read one member value given its declared BinaryType/AdditionalInfo
        (used inside ClassWithMembersAndTypes only)."""
        if bt == BT_PRIMITIVE:
            return self.read_primitive(add_info)
        # String, Object, SystemClass, Class, arrays: all show up as a
        # regular record (string/null/ref/nested class/array record).
        rt = self.read_byte()
        val = self._read_record(rt)
        return val

    # -- class records -------------------------------------------------
    def _read_class_with_members_and_types(self, system: bool):
        ci = self.read_class_info()
        self.read_member_type_info(ci)
        library_id = None
        if not system:
            library_id = self.read_int32()
        obj = NrbfObject(ci.name, ci.object_id)
        self.objects[ci.object_id] = obj  # register before reading members (cycles)
        self.class_templates[ci.object_id] = ci
        for name, bt, info in zip(ci.member_names, ci.binary_types, ci.additional_infos):
            obj.members[name] = self._resolve(self._read_value_member(bt, info))
        return obj

    def _read_class_with_members(self, system: bool):
        # No type info given inline; values are plain records. Rare in
        # practice for .nox files, but handle it for completeness.
        ci = self.read_class_info()
        if not system:
            self.read_int32()  # library id
        obj = NrbfObject(ci.name, ci.object_id)
        self.objects[ci.object_id] = obj
        self.class_templates[ci.object_id] = ci
        for name in ci.member_names:
            rt = self.read_byte()
            obj.members[name] = self._resolve(self._read_record(rt))
        return obj

    def _read_class_with_id(self):
        object_id = self.read_int32()
        metadata_id = self.read_int32()
        ci = self.class_templates[metadata_id]
        obj = NrbfObject(ci.name, object_id)
        self.objects[object_id] = obj
        for name, bt, info in zip(ci.member_names, ci.binary_types, ci.additional_infos):
            obj.members[name] = self._resolve(self._read_value_member(bt, info))
        return obj

    def _read_binary_object_string(self):
        object_id = self.read_int32()
        s = self.read_string()
        self.objects[object_id] = s
        return s

    # -- arrays ----------------------------------------------------------
    def _read_array_single_primitive(self):
        object_id = self.read_int32()
        length = self.read_int32()
        ptype = self.read_byte()
        arr = NrbfArray(element_type=PRIMITIVE_NAMES.get(ptype, ptype), object_id=object_id)
        self.objects[object_id] = arr
        if ptype == PT_STRING or ptype == PT_DECIMAL or ptype == PT_CHAR:
            for _ in range(length):
                arr.append(self.read_primitive(ptype))
        else:
            fmt = _PRIMITIVE_STRUCT[ptype]
            size = struct.calcsize(fmt)
            raw = self._read(size * length)
            arr.extend(struct.unpack_from(f"<{length}{fmt[1]}", raw))
        return arr

    def _read_array_single_object(self):
        object_id = self.read_int32()
        length = self.read_int32()
        arr = NrbfArray(element_type="Object", object_id=object_id)
        self.objects[object_id] = arr
        arr.extend([None] * length)
        i = 0
        while i < length:
            rt = self.read_byte()
            val = self._read_record(rt)
            if isinstance(val, _NullRun):
                i += val.count
                continue
            arr[i] = self._resolve(val)
            i += 1
        return arr

    def _read_array_single_string(self):
        object_id = self.read_int32()
        length = self.read_int32()
        arr = NrbfArray(element_type="String", object_id=object_id)
        self.objects[object_id] = arr
        arr.extend([None] * length)
        i = 0
        while i < length:
            rt = self.read_byte()
            val = self._read_record(rt)
            if isinstance(val, _NullRun):
                i += val.count
                continue
            arr[i] = self._resolve(val)
            i += 1
        return arr

    def _read_binary_array(self):
        object_id = self.read_int32()
        array_type = self.read_byte()  # 0..5
        rank = self.read_int32()
        lengths = [self.read_int32() for _ in range(rank)]
        if array_type in (3, 4, 5):  # *Offset variants
            _lower_bounds = [self.read_int32() for _ in range(rank)]
        bt = self.read_byte()
        if bt in (BT_PRIMITIVE, BT_PRIMITIVE_ARRAY):
            add_info = self.read_byte()
        elif bt == BT_SYSTEM_CLASS:
            add_info = self.read_string()
        elif bt == BT_CLASS:
            add_info = (self.read_string(), self.read_int32())
        else:
            add_info = None

        total = 1
        for l in lengths:
            total *= l
        arr = NrbfArray(element_type=bt, object_id=object_id)
        self.objects[object_id] = arr

        if bt == BT_PRIMITIVE and add_info not in (PT_STRING, PT_DECIMAL, PT_CHAR):
            fmt = _PRIMITIVE_STRUCT[add_info]
            size = struct.calcsize(fmt)
            raw = self._read(size * total)
            arr.extend(struct.unpack_from(f"<{total}{fmt[1]}", raw))
        else:
            arr.extend([None] * total)
            i = 0
            while i < total:
                rt = self.read_byte()
                val = self._read_record(rt)
                if isinstance(val, _NullRun):
                    i += val.count
                    continue
                arr[i] = self._resolve(val)
                i += 1
        arr.dims = lengths
        return arr

    # -- reference resolution ---------------------------------------------
    def _resolve(self, val):
        """MemberReference values are resolved lazily via a fixup list if
        the target isn't registered yet (shouldn't normally happen given
        NRBF ordering guarantees, but handled defensively)."""
        if isinstance(val, _Ref):
            if val.id in self.objects:
                return self.objects[val.id]
            return val  # left unresolved; fixed up in a final pass
        return val

    def resolve_all_refs(self, root):
        """Best-effort second pass to replace any remaining _Ref placeholders
        (forward references) now that the whole graph is parsed."""
        seen = set()

        def fix(container, key_or_index):
            v = container[key_or_index]
            if isinstance(v, _Ref):
                if v.id in self.objects:
                    container[key_or_index] = self.objects[v.id]

        def walk(node):
            oid = id(node)
            if oid in seen:
                return
            seen.add(oid)
            if isinstance(node, NrbfObject):
                for k in list(node.members.keys()):
                    fix(node.members, k)
                    walk(node.members[k])
            elif isinstance(node, list):
                for i in range(len(node)):
                    fix(node, i)
                    walk(node[i])
            elif isinstance(node, dict):
                for k in list(node.keys()):
                    walk(node[k])

        walk(root)
        return root


class _Ref:
    __slots__ = ("id",)

    def __init__(self, id_):
        self.id = id_


class _NullRun:
    __slots__ = ("count",)

    def __init__(self, count):
        self.count = count


def parse_nrbf(data: bytes):
    reader = NrbfReader(data)
    root = reader.parse()
    reader.resolve_all_refs(root)
    return root, reader


# =============================================================================
# 2. NOVA/.nox-specific object-graph extraction (Procedure, Curve, load_nox)
# =============================================================================

_HEADER_SIG = bytes([0x00]) + struct.pack("<i", 1) + struct.pack("<i", -1) + struct.pack("<i", 1) + struct.pack("<i", 0)

# .NET DateTime is serialized as a single UInt64: the low 62 bits are Ticks
# (100ns units since 0001-01-01), the top 2 bits are DateTimeKind. Masking
# off the Kind bits and converting the remaining ticks reproduces the exact
# wall-clock moment NOVA recorded -- verified against both sample files'
# on-disk modification times (matched to within ~1s, the OS flush delay).
_DOTNET_TICKS_MASK = 0x3FFFFFFFFFFFFFFF
_DOTNET_EPOCH = datetime(1, 1, 1)


def _dotnet_datetime(raw) -> Optional[datetime]:
    if raw is None:
        return None
    ticks = raw & _DOTNET_TICKS_MASK
    try:
        return _DOTNET_EPOCH + timedelta(microseconds=ticks / 10)
    except OverflowError:
        return None


def _split_streams(data: bytes) -> list[bytes]:
    """Split a .nox file into its concatenated NRBF streams."""
    positions = []
    idx = 0
    while True:
        i = data.find(_HEADER_SIG, idx)
        if i == -1:
            break
        positions.append(i)
        idx = i + 1
    if not positions:
        raise ValueError("No NRBF stream header found -- not a .nox file?")
    positions.append(len(data))
    return [data[positions[i]:positions[i + 1]] for i in range(len(positions) - 1)]


@dataclass
class Curve:
    """A set of named, equal-length data channels sharing a parent command
    (e.g. Time / Potential applied / WE(1).Current for one CV scan)."""
    name: str
    length: int
    columns: dict  # column label -> unit string
    dataframe: pd.DataFrame

    def to_csv(self, path: str, **kwargs):
        self.dataframe.to_csv(path, index=False, **kwargs)


@dataclass
class Procedure:
    name: Optional[str]
    text: Optional[str]
    instrument: Optional[str]
    remarks: list
    estimated_duration: Optional[float]
    rating: Optional[int]
    start_time: Optional[datetime] = None
    modified_time: Optional[datetime] = None
    curves: list = field(default_factory=list)
    root: object = None       # raw NrbfObject graph, for advanced digging
    reader: object = None     # the NrbfReader (object table, libraries, ...)
    n_streams: int = 1
    setpoint_events: list = field(default_factory=list)
    """Every applied-setpoint command found ANYWHERE in the procedure's
    command tree, in chronological order, as {"t_rel", "value", "unit",
    "label", "command_class"} dicts. "Setpoint" here means: NOVA's own
    Autolab.Parameters.*SetpointParameter family -- the parameter class it
    uses for a command that tells the potentiostat/galvanostat to *apply* a
    value, as opposed to e.g. a duration or a channel selector. This is
    deliberately not tied to any one command class (like the top-level
    FHSetSetpoint "Set potential" command): a NOVA "Levels" technique nests
    its own per-level "Step" sub-commands (class LevelShortSetpoint) many
    layers deep in the tree, each with its own real applied value and its
    own timestamp, and other techniques may nest setpoint commands
    differently still. Keying off the underlying parameter class instead of
    the owning command's class means this keeps working (without being
    taught about each new command class by name) for whatever nesting
    structure a given technique happens to use -- chronoamperometry,
    staircase/differential-pulse sub-steps, galvanostatic current setpoints
    (unit "A" instead of "V"), etc.

    `unit` is whatever NOVA's own CommandParameter+_unit says for that
    command ("V" for a plain "Set potential" command's parameter; the
    NOVA "Levels" Step commands observed so far report unit "" even though
    the values are volts -- see `setpoint_values_at`, which treats "V" and
    "" as the same potential-like bucket by default).
    """

    def setpoint_values_at(self, t_rel, unit=("V", None)):
        """Piecewise-constant setpoint value(s) at elapsed time(s) `t_rel`
        (seconds since `start_time`), built from `setpoint_events` filtered
        to `unit` (default: potential-like events -- NOVA reports unit "V"
        for most, but "" i.e. None for some nested ones, like the "Levels"
        Step sub-commands; both are treated as potential by default).
        `t_rel` may be a scalar or a numpy array. Returns NaN before the
        first matching event (nothing had been applied yet).

        Pass e.g. `unit=("A",)` to reconstruct a *current* setpoint trace
        instead, for a galvanostatic technique whose setpoint commands
        report unit "A" -- the same generic mechanism covers it, no new
        code needed."""
        import numpy as np

        events = [e for e in self.setpoint_events if e["unit"] in unit]
        if not events:
            raise ValueError(
                f"no setpoint events with unit in {unit!r} found in this procedure "
                f"(available units: {sorted({e['unit'] for e in self.setpoint_events})!r})"
            )
        xs = np.asarray([e["t_rel"] for e in events], dtype=float)
        vs = np.asarray([e["value"] for e in events], dtype=float)
        t_rel = np.asarray(t_rel, dtype=float)
        idx = np.searchsorted(xs, t_rel, side="right") - 1
        out = np.where(idx >= 0, vs[np.clip(idx, 0, len(vs) - 1)], np.nan)
        return out if out.ndim else float(out)

    def potential_at(self, t_rel):
        """Piecewise-constant applied potential (V) at elapsed time(s)
        `t_rel` -- convenience alias for
        `setpoint_values_at(t_rel, unit=("V", None))`."""
        return self.setpoint_values_at(t_rel, unit=("V", None))

    def summary(self) -> str:
        lines = [
            f"name: {self.name}",
            f"text: {self.text}",
            f"instrument: {self.instrument}",
            f"remarks: {self.remarks}",
            f"estimated_duration_s: {self.estimated_duration}",
            f"start_time: {self.start_time}",
            f"modified_time: {self.modified_time}",
        ]
        for c in self.curves:
            lines.append(f"curve '{c.name}': {c.length} points, columns={list(c.columns)}")
        return "\n".join(lines)

    def plot_preview_png_bytes(self) -> Optional[bytes]:
        """Return the embedded thumbnail PNG NOVA shows in its file browser,
        if present (EcoChemie.Shared.PlotPreview._pngImageBytes)."""
        preview = _get(self.root, "<PlotPreview>k__BackingField")
        raw = _get(preview, "_pngImageBytes")
        if isinstance(raw, NrbfArray):
            return bytes(raw)
        return None

    def save_plot_preview(self, path: str) -> bool:
        png = self.plot_preview_png_bytes()
        if png is None:
            return False
        with open(path, "wb") as f:
            f.write(png)
        return True


def _get(obj, key, default=None):
    if isinstance(obj, NrbfObject):
        return obj.members.get(key, default)
    return default


def _is_double_list(val) -> bool:
    return (
        isinstance(val, NrbfObject)
        and val.class_name.startswith("System.Collections.Generic.List`1[[System.Double")
    )


def _find_data_array_channels(root):
    """Walk the whole object graph and collect every
    EcoChemie.Utils.Sequencer.CommandParameterDataArray that wraps a plain
    List<Double> (the processed/calculated channels -- raw ADC buffer
    channels are skipped since decoding them needs gain/offset calibration
    tables not exposed at this layer)."""
    channels = []
    seen = set()

    def walk(node):
        oid = id(node)
        if oid in seen:
            return
        seen.add(oid)
        if isinstance(node, NrbfObject):
            if node.class_name == "EcoChemie.Utils.Sequencer.CommandParameterDataArray":
                param = node.members.get("_parameter")
                if isinstance(param, NrbfObject) and param.class_name == "EcoChemie.Utils.Sequencer.ParameterObject":
                    val = param.members.get("_value")
                    if _is_double_list(val):
                        items = val.members.get("_items")
                        size = val.members.get("_size", len(items) if items is not None else 0)
                        parent = node.members.get("CommandParameter+_parent")
                        channels.append({
                            "key": node.members.get("CommandParameter+_name"),
                            "label": node.members.get("CommandParameter+_text"),
                            "unit": node.members.get("CommandParameter+_unit"),
                            "data": list(items[:size]) if items is not None else [],
                            "length": size,
                            "parent_id": id(parent) if parent is not None else None,
                            "parent_name": _get(parent, "FunctionHandler+_text") or _get(parent, "CommandBase+_text"),
                        })
            for v in node.members.values():
                if isinstance(v, (NrbfObject, NrbfArray)):
                    walk(v)
        elif isinstance(node, NrbfArray):
            for v in node:
                if isinstance(v, (NrbfObject, NrbfArray)):
                    walk(v)

    walk(root)
    return channels


def _group_into_curves(channels) -> list[Curve]:
    # Group by parent command identity.
    cohorts = {}
    order = []
    for ch in channels:
        key = ch["parent_id"]
        if key not in cohorts:
            cohorts[key] = {"parent_name": ch["parent_name"], "items": []}
            order.append(key)
        cohorts[key]["items"].append(ch)

    # De-duplicate cohorts whose set of underlying data (by value, cheaply
    # approximated with (key, length, first value, last value)) is a subset
    # of an already-accepted, larger cohort -- this happens when a wrapper
    # command re-exposes its child's channels under the same names.
    #
    # NOTE: this used to compare only (key, length), which is a *shape*
    # check, not a content check. That silently collapsed distinct sibling
    # commands that happen to produce same-shaped data (e.g. N repeated
    # "RecordLevelsContainer" bursts of a chronoamperometry loop, each with
    # its own ~30s of real, different data) into a single kept cohort,
    # discarding the rest as if they were duplicates. Including the first/
    # last data values in the fingerprint fixes that while still catching
    # genuine wrapper-reexposes-child-data duplicates (which really do have
    # identical values, not just identical shape).
    def signature(items):
        return frozenset(
            (
                it["key"],
                it["length"],
                it["data"][0] if it["data"] else None,
                it["data"][-1] if it["data"] else None,
            )
            for it in items
            if it["length"] > 1
        )

    accepted = []  # list of (signature, cohort)
    for key in order:
        cohort = cohorts[key]
        sig = signature(cohort["items"])
        if not sig:
            continue
        if any(sig <= existing_sig for existing_sig, _ in accepted):
            continue
        # drop any previously accepted cohort now subsumed by this one
        accepted = [(s, c) for s, c in accepted if not (s <= sig)]
        accepted.append((sig, cohort))

    curves = []
    from collections import Counter

    raw_curves = []
    for i, (sig, cohort) in enumerate(accepted):
        items = [it for it in cohort["items"] if it["length"] > 1]
        if not items:
            continue
        # keep channels that share the majority length (the actual curve
        # table -- a cohort can contain a few stray scalars/short arrays)
        common_len, _ = Counter(it["length"] for it in items).most_common(1)[0]
        rows = [it for it in items if it["length"] == common_len]
        if len(rows) < 2:
            continue  # a single lone column isn't a usable curve

        cols = {}
        data = {}
        used_labels = set()
        for it in rows:
            label = it["label"] or it["key"] or "value"
            if label in used_labels:
                label = f"{label} ({it['key']})"
            used_labels.add(label)
            cols[label] = it["unit"] or ""
            data[label] = it["data"]
        df = pd.DataFrame(data)
        name = cohort["parent_name"] or f"curve_{i+1}"
        raw_curves.append(Curve(name=name, length=common_len, columns=cols, dataframe=df))

    # disambiguate curves that share a parent command name (e.g. repeated
    # scans / segments of the same technique) by appending point counts
    name_counts = Counter(c.name for c in raw_curves)
    seen_counts = Counter()
    for c in raw_curves:
        if name_counts[c.name] > 1:
            seen_counts[c.name] += 1
            c.name = f"{c.name} #{seen_counts[c.name]} ({c.length} pts)"

    return raw_curves


def _add_absolute_timestamps(curves: list[Curve], start_time: Optional[datetime]):
    """Every curve carries an elapsed-time column ('Time', unit 's') that is
    seconds since the *procedure* started (continuous across curves/segments,
    confirmed by consecutive curves' Time ranges abutting each other). Turn
    that into a real, per-datapoint wall-clock timestamp using the
    procedure's decoded start time."""
    if start_time is None:
        return
    for c in curves:
        time_col = None
        for label, unit in c.columns.items():
            if unit == "s" and label.strip().lower() == "time":
                time_col = label
                break
        if time_col is None:
            continue
        ts = pd.Timestamp(start_time) + pd.to_timedelta(c.dataframe[time_col], unit="s")
        c.dataframe.insert(0, "Timestamp", ts)
        c.columns["Timestamp"] = ""


def _command_scalar_parameters(cmd):
    """Return [(text, unit, value, hparam_class), ...] for every scalar
    (numeric) command parameter attached to a command-tree node, e.g. a
    FHSetSetpoint node's "Potential (V)" = 0.85 parameter. Parameters whose
    value is itself a nested object (like the HObjectParameter "channel
    selector" ones NOVA also stores per command) are skipped -- only plain
    numbers come out. `hparam_class` is the short class name of the
    `_hParameter` wrapper itself (e.g. "HSetpointParameter",
    "HDoubleParameter") -- callers use this to tell *what kind* of
    parameter a scalar is (an applied setpoint vs. e.g. a plain duration)
    without depending on which specific command it came from.

    Parameters live under `FunctionHandler+_commandParameters13` for
    FunctionHandler-subclass nodes (FHSetSetpoint, FHWait, LevelShortSetpoint,
    ...) or under `CommandBase+_commandParameters13` for other CommandBase
    nodes -- NOVA uses the FunctionHandler-prefixed field for the nodes that
    matter here. The collection itself is double-wrapped:
    `_commandParameters13.[[Collection`1+items]]._items[:._size]`, each
    item a CommandHParameter whose `_hParameter` -> (`_innerIParameter` /
    `HParameter+_innerIParameter` / `H*Parameter+_innerIParameter`) ->
    ParameterObject._value is the actual scalar.
    """
    p13 = _get(cmd, "FunctionHandler+_commandParameters13") or _get(
        cmd, "CommandBase+_commandParameters13"
    )
    if not isinstance(p13, NrbfObject):
        return []
    wrap = p13.members.get("Collection`1+items")
    if not isinstance(wrap, NrbfObject):
        return []
    arr = wrap.members.get("_items")
    size = wrap.members.get("_size", len(arr) if arr is not None else 0)
    if not isinstance(arr, NrbfArray):
        return []

    out = []
    for hp in list(arr)[:size]:
        if not isinstance(hp, NrbfObject):
            continue
        text = _get(hp, "CommandParameter+_text")
        unit = _get(hp, "CommandParameter+_unit")
        hparam = _get(hp, "_hParameter")
        hparam_class = hparam.class_name if isinstance(hparam, NrbfObject) else None
        inner = None
        for key in (
            "_innerIParameter",
            "HParameter+_innerIParameter",
            "HDoubleParameter+_innerIParameter",
        ):
            inner = _get(hparam, key)
            if inner is not None:
                break
        value = _get(inner, "_value")
        if isinstance(value, (int, float)):
            out.append((text, unit, value, hparam_class))
    return out


def _resolve_command_t_rel(node, start_time: Optional[datetime]):
    """Best-effort elapsed-time-since-procedure-start (seconds) for when a
    command-tree node actually executed. Tries, in order:

    1. A real .NET DateTime in `FunctionHandler+_timeStamp` /
       `CommandBase+_timeStamp` / `ExecCommandBase+_timeStamp` (what
       top-level procedure commands like FHSetSetpoint carry).
    2. `FunctionHandler+_timeStampEmbedded2`: some deeply-nested,
       embedded-controller-generated sub-commands (observed so far on a
       NOVA "Levels" technique's per-level "Step"/LevelShortSetpoint
       commands) don't get a real DateTime timestamp (it reads as the .NET
       default, tick 0) but do carry this integer field -- empirically an
       elapsed **microsecond** count since the procedure's start_time
       (cross-checked against the owning data segment's own elapsed-time
       values; matched to sub-millisecond precision). This convention is
       observed, not documented by Metrohm, so treat it as best-effort.

    Returns None if neither is available/usable.
    """
    ts_raw = (
        _get(node, "FunctionHandler+_timeStamp")
        or _get(node, "CommandBase+_timeStamp")
        or _get(node, "ExecCommandBase+_timeStamp")
    )
    if ts_raw and start_time is not None:
        ts = _dotnet_datetime(ts_raw)
        if ts is not None:
            return (ts - start_time).total_seconds()
    te2 = _get(node, "FunctionHandler+_timeStampEmbedded2")
    if isinstance(te2, (int, float)) and te2:
        return te2 / 1e6
    return None


def _find_setpoint_events(root, start_time: Optional[datetime]) -> list:
    """Walk the *entire* command tree (no class-name filtering) and collect
    every command parameter whose underlying Autolab parameter class is one
    of the `*SetpointParameter` family (NOVA's own marker for "this
    parameter is an applied setpoint", as opposed to e.g. a duration or a
    channel selector), each with the time it actually executed (see
    `_resolve_command_t_rel`). Deliberately not restricted to any specific
    command class (like FHSetSetpoint) so that setpoint structure nested
    inside other technique-specific commands -- e.g. a "Levels" technique's
    per-level LevelShortSetpoint "Step" sub-commands, found many layers
    below their FHSetSetpoint ancestor -- is picked up the same way, without
    this module needing to be taught about every command class NOVA has.
    Returns a chronologically sorted list of {"t_rel", "value", "unit",
    "label", "command_class"} dicts."""
    events = []
    seen = set()

    def walk(node):
        oid = id(node)
        if oid in seen:
            return
        seen.add(oid)
        if isinstance(node, NrbfObject):
            for text, unit, value, hparam_class in _command_scalar_parameters(node):
                if hparam_class and hparam_class.split(".")[-1].endswith("SetpointParameter"):
                    t_rel = _resolve_command_t_rel(node, start_time)
                    if t_rel is not None:
                        events.append(
                            {
                                "t_rel": t_rel,
                                "value": value,
                                "unit": unit or None,
                                "label": text or "Setpoint",
                                "command_class": node.class_name.split(".")[-1],
                            }
                        )
            for v in node.members.values():
                if isinstance(v, (NrbfObject, NrbfArray)):
                    walk(v)
        elif isinstance(node, NrbfArray):
            for v in node:
                if isinstance(v, (NrbfObject, NrbfArray)):
                    walk(v)

    walk(root)
    events.sort(key=lambda e: e["t_rel"])
    return events


def load_nox(path: str) -> Procedure:
    with open(path, "rb") as f:
        data = f.read()
    streams = _split_streams(data)
    # The last stream is the fully-populated object graph.
    root, reader = parse_nrbf(streams[-1])

    channels = _find_data_array_channels(root)
    curves = _group_into_curves(channels)

    start_time = _dotnet_datetime(_get(root, "_timeStamp"))
    modified_time = _dotnet_datetime(_get(root, "_modifiedTimeStamp"))
    _add_absolute_timestamps(curves, start_time)

    setpoint_events = _find_setpoint_events(root, start_time)

    remarks = _get(root, "_remarks")
    remarks_list = list(remarks) if isinstance(remarks, NrbfArray) else ([] if remarks is None else [remarks])

    return Procedure(
        name=_get(root, "_name"),
        text=_get(root, "_text"),
        instrument=_get(root, "_instrument"),
        remarks=remarks_list,
        estimated_duration=_get(root, "_estimatedDuration"),
        rating=_get(root, "_rating"),
        start_time=start_time,
        modified_time=modified_time,
        curves=curves,
        root=root,
        reader=reader,
        n_streams=len(streams),
        setpoint_events=setpoint_events,
    )


# =============================================================================
# 3. The ixdat reader itself
# =============================================================================

def _find_time_column(curve: Curve) -> str:
    for label, unit in curve.columns.items():
        if unit == "s" and label.strip().lower() == "time":
            return label
    raise ValueError(
        f"curve {curve.name!r} has no elapsed-time ('Time', unit 's') column"
    )


def _guess_potential_current_columns(curve: Curve):
    """Best-effort pick of the measured (not setpoint) potential and
    current columns, based on NOVA's own unit + label conventions."""
    potential_col = current_col = None
    for label, unit in curve.columns.items():
        low = label.lower()
        if unit == "V" and "potential" in low and "applied" not in low and potential_col is None:
            potential_col = label
        if unit == "A" and "current" in low and current_col is None:
            current_col = label
    if potential_col is None:
        for label, unit in curve.columns.items():
            if unit == "V":
                potential_col = label
                break
    if current_col is None:
        for label, unit in curve.columns.items():
            if unit == "A":
                current_col = label
                break
    return potential_col, current_col


def _concat_curves(curves: list, name: str = None) -> Curve:
    """Concatenate several Curves that share the same columns and a
    continuous elapsed-Time axis (e.g. several sequential scans of the same
    technique in one .nox procedure) into a single Curve."""
    time_cols = {_find_time_column(c) for c in curves}
    if len(time_cols) != 1:
        raise ValueError(
            "the .nox file's curves don't share a common Time column, "
            "can't concatenate -- pass concatenate=False and curve_index= "
            "to read a single curve instead"
        )
    time_col = time_cols.pop()
    ordered = sorted(curves, key=lambda c: c.dataframe[time_col].iloc[0])
    df = pd.concat([c.dataframe for c in ordered], ignore_index=True)
    columns = {}
    for c in ordered:
        columns.update(c.columns)
    return Curve(name=name or ordered[0].name, length=len(df), columns=columns, dataframe=df)


class NovaNoxReader:
    """Reader for Metrohm Autolab NOVA's native binary ".nox" procedure
    files. Returns an ECMeasurement (or a `cls=` subclass, e.g.
    CyclicVoltammogram) directly -- see module docstring for usage.
    """

    def read(
        self,
        path_to_file,
        cls=None,
        name=None,
        curve_index=None,
        concatenate=True,
        potential_col=None,
        current_col=None,
        use_setpoint_potential=True,
        **kwargs,
    ):
        """Read a NOVA .nox file into an ixdat ECMeasurement (or subclass).

        Args:
            path_to_file (Path or str): path to the .nox file.
            cls (Measurement subclass): defaults to ECMeasurement. Pass
                CyclicVoltammogram for a CV technique that recorded a
                per-point potential channel (auto-detects cycles).
            name (str): defaults to the procedure's own title
                (`Procedure.text`, NOVA's "Chrono amperometry fast"-style
                name), falling back to the curve's name.
            curve_index (int): which curve to read if the procedure has
                more than one and `concatenate` is False (default: 0).
            concatenate (bool): if True (default) and the procedure has
                more than one curve, concatenate them all into a single
                continuous measurement (see `_concat_curves`). If False,
                use `curve_index` to select a single curve/segment.
            potential_col / current_col (str): column labels to use as
                ixdat's "raw_potential" / "raw_current" aliases.
                Auto-detected from column units/labels if not given.
            use_setpoint_potential (bool): Some techniques (e.g. a fixed-
                potential chronoamperometry hold) record no per-point
                measured-potential channel at all -- NOVA only stores the
                *setpoint* (target) potential as a scalar command parameter
                on the FHSetSetpoint ("Set potential") command(s) in the
                procedure's command tree, once per potential step, not once
                per data point. If True (default) and no potential_col was
                found/given, this reconstructs a real, piecewise-constant
                "raw_potential" series from those FHSetSetpoint commands'
                values and their own execution timestamps (so e.g. a 3-step
                0.85 V / 1.0 V / 1.45 V staircase experiment gets a genuine
                step-function raw_potential instead of nothing). Set to
                False to skip this and leave "raw_potential" unset (which
                will make `Measurement.read(..., reader="nova_nox")` /
                `ECMeasurement.read(...)` raise SeriesNotFoundError, since
                raw_potential is an essential series for ECMeasurement --
                use `NovaNoxReader().read(...)` directly in that case).
            **kwargs: additional key-word arguments, passed through to
                `cls.from_dict` (e.g. to override "technique").
        """
        path_to_file = Path(path_to_file)
        proc = load_nox(path_to_file)
        if not proc.curves:
            raise ValueError(
                f"{path_to_file}: .nox file has no curves (no calculated "
                "data channels found in the object graph)"
            )
        if proc.start_time is None:
            raise ValueError(
                f"{path_to_file}: couldn't recover the procedure's start "
                "time (._timeStamp) from the .nox file"
            )
        tstamp = proc.start_time.timestamp()

        if len(proc.curves) > 1 and concatenate:
            curve = _concat_curves(proc.curves, name=name or proc.text)
        else:
            curve = proc.curves[curve_index or 0]

        time_col = _find_time_column(curve)
        tseries = TimeSeries(
            name=time_col,
            unit_name="s",
            data=curve.dataframe[time_col].to_numpy(),
            tstamp=tstamp,
        )
        series_list = [tseries]
        for label, unit in curve.columns.items():
            if label in (time_col, "Timestamp"):
                continue
            series_list.append(
                ValueSeries(
                    name=label,
                    unit_name=unit or "",
                    data=curve.dataframe[label].to_numpy(),
                    tseries=tseries,
                )
            )

        if potential_col is None or current_col is None:
            guessed_p, guessed_c = _guess_potential_current_columns(curve)
            potential_col = potential_col or guessed_p
            current_col = current_col or guessed_c

        has_potential_setpoints = any(
            e["unit"] in ("V", None) for e in proc.setpoint_events
        )
        if potential_col is None and use_setpoint_potential and has_potential_setpoints:
            setpoint_label = "Potential setpoint applied"
            v_setpoint = proc.potential_at(curve.dataframe[time_col].to_numpy())
            series_list.append(
                ValueSeries(
                    name=setpoint_label,
                    unit_name="V",
                    data=v_setpoint,
                    tseries=tseries,
                )
            )
            potential_col = setpoint_label

        aliases = {"t": [time_col]}
        if potential_col:
            aliases["raw_potential"] = [potential_col]
        if current_col:
            aliases["raw_current"] = [current_col]

        if cls is None:
            from ..techniques.ec import ECMeasurement

            cls = ECMeasurement

        technique = kwargs.pop("technique", None)
        if technique is None:
            from ..techniques.cv import CyclicVoltammogram

            technique = "CV" if issubclass(cls, CyclicVoltammogram) else "EC"

        obj_as_dict = dict(
            name=name or proc.text or curve.name,
            technique=technique,
            reader=self,
            aliases=aliases,
            series_list=series_list,
            tstamp=tstamp,
        )
        obj_as_dict.update(kwargs)
        return cls.from_dict(obj_as_dict)