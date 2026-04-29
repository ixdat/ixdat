import re
from pathlib import Path
import numpy as np
from .. import Spectrum
from ..spectra import MultiSpectrum
from ..data_series import DataSeries, Field
from .reading_tools import timestamp_string_to_tstamp

SPEC_DATE_FORM = "%a %b %d %H:%M:%S %Y"
FLUORESCENCE_COLUMNS = ["Det_1", "Det_2", "Det_3", "Det_4", "Det_5", "Det_6"]


class SpecDATReader:
    """Reader for SPEC-format .dat files from synchrotron beamlines.

    SPEC (Certified Scientific Software) is used at ESRF (e.g. BM31), APS, SOLEIL,
    and many other facilities. The test data for this reader is Fe foil EXAFS from
    ESRF BM31.

    Do not confuse with QexafsDATReader (reader="qexafs"), which reads Diamond B18-Core
    files. Those use a simple `# Key: Value` header, tab-separated columns, and energy
    already in eV. SPEC files use a structured header (#F, #E, #S, #L, #O, #P, ...),
    space-separated columns, multiple scans per file delimited by #S, and energy in keV.

    One SPEC file can contain multiple scans (#S 1, #S 2, ...). This reader can return
    a single scan, a subset, or an average over all scans. Column names are
    beamline-specific; the caller must identify the I0 and It columns explicitly for XAS.

    Metadata stored on the returned object includes the scan number, scan command,
    count time, motor positions (#O/#P), UMI lines, original beamline file path (#F),
    and experiment comment (#C).
    """

    def read(
        self,
        path_to_file,
        cls=Spectrum,
        technique=None,
        x_name="ZapEnergy",
        y_name=None,
        ref_name=None,
        scan_numbers=None,
        average_scans=False,
        **kwargs,
    ):
        """Read a SPEC .dat file.

        Args:
            path_to_file (Path or str): Path to the .dat file.
            cls: Class to return. Defaults to Spectrum.
            technique (str): "XAS" for transmission EXAFS (returns -log(It/I0)),
                "XAS_fluorescence" for fluorescence EXAFS (returns sum(Det)/I0).
                None returns a MultiSpectrum with every column as a Field.
            x_name (str): Column used as x axis. Defaults to "ZapEnergy".
                Values are converted from keV to eV.
            y_name (str): Column for y. Required for "XAS" and "XAS_fluorescence"
                (column names are beamline-specific). For "XAS_fluorescence" defaults
                to the sum of Det_1..Det_6 if omitted.
            ref_name (str): Column to normalise y by (I0). Required for "XAS" and
                "XAS_fluorescence". Pass "none" (string) to skip normalisation.
                Column name is beamline-specific.
            scan_numbers (int or list of int): Scan(s) to read. None reads all.
            average_scans (bool): Average y over all selected scans (XAS only).
            **kwargs: Passed to the returned object's initialiser. Explicit values
                override reader defaults (name, metadata, duration, reader).
        """
        path_to_file = Path(path_to_file)
        file_meta, scans = _parse_spec_file(path_to_file)

        if scan_numbers is not None:
            if isinstance(scan_numbers, int):
                scan_numbers = [scan_numbers]
            scans = [s for s in scans if s["number"] in scan_numbers]

        if not scans:
            raise ValueError(f"No matching scans found in {path_to_file}")

        kwargs.setdefault("reader", self)

        if technique in ("XAS", "XAS_fluorescence"):
            if y_name is None and technique == "XAS":
                cols = list(scans[0]["data"].keys())
                raise ValueError(
                    f"technique='XAS' requires explicit y_name. "
                    f"Available columns: {cols}"
                )
            if ref_name is None:
                cols = list(scans[0]["data"].keys())
                raise ValueError(
                    f"technique='{technique}' requires explicit ref_name (I0 column). "
                    f"Available columns: {cols}"
                )
            return self._read_xas(
                scans,
                cls,
                technique,
                x_name,
                y_name,
                ref_name,
                average_scans,
                path_to_file,
                file_meta,
                **kwargs,
            )

        return _multispectrum_from_scan(
            scans[0], x_name, path_to_file, file_meta, technique=technique, **kwargs
        )

    def _read_xas(
        self,
        scans,
        cls,
        technique,
        x_name,
        y_name,
        ref_name,
        average_scans,
        path_to_file,
        file_meta,
        **kwargs,
    ):
        def get_xy(scan):
            data = scan["data"]
            x = data[x_name] * 1000.0  # keV -> eV
            if y_name is None:
                y = sum(data[c] for c in FLUORESCENCE_COLUMNS if c in data)
            else:
                y = data[y_name].copy()
            if ref_name != "none":
                y = y / data[ref_name]
            if technique == "XAS":
                with np.errstate(divide="ignore", invalid="ignore"):
                    y = -np.log(y)
            return x, y

        x, y = get_xy(scans[0])
        if average_scans and len(scans) > 1:
            for scan in scans[1:]:
                xi, yi = get_xy(scan)
                y = y + np.interp(x, xi, yi)
            y = y / len(scans)

        y_label = "µt" if technique == "XAS" else "fluorescence/I0"
        xseries = DataSeries(name="energy", unit_name="eV", data=x)
        yseries = DataSeries(name=y_label, unit_name="", data=y)

        count_time_ms = scans[0]["count_time_ms"]
        duration = count_time_ms * len(x) / 1000.0 if count_time_ms is not None else None

        kwargs.setdefault("name", path_to_file.name)
        kwargs.setdefault("technique", technique)
        kwargs.setdefault("duration", duration)
        kwargs.setdefault("metadata", _build_metadata(scans, file_meta, average_scans))

        return cls.from_series(xseries, yseries, tstamp=scans[0]["tstamp"], **kwargs)


def _multispectrum_from_scan(scan, x_name, path_to_file, file_meta, technique, **kwargs):
    data = scan["data"]
    x = data[x_name] * 1000.0  # keV -> eV
    xseries = DataSeries(name="energy", unit_name="eV", data=x)
    fields = [
        Field(name=col, unit_name="", data=arr, axes_series=[xseries])
        for col, arr in data.items()
        if col != x_name
    ]
    kwargs.setdefault("name", path_to_file.name)
    kwargs.setdefault("technique", technique)
    kwargs.setdefault(
        "metadata", _build_metadata([scan], file_meta, average_scans=False)
    )
    # MultiSpectrum does not accept reader or duration, so remove them
    kwargs.pop("reader", None)
    kwargs.pop("duration", None)
    return MultiSpectrum(fields=fields, tstamp=scan["tstamp"], **kwargs)


def _build_metadata(scans, file_meta, average_scans):
    """Build a JSON-serialisable metadata dict from parsed SPEC scan info."""
    if average_scans and len(scans) > 1:
        meta = {
            "scan_numbers": [s["number"] for s in scans],
            "scan_command": scans[0]["command"],
        }
        if scans[0]["count_time_ms"] is not None:
            meta["count_time_ms"] = scans[0]["count_time_ms"]
    else:
        s = scans[0]
        meta = {
            "scan_number": s["number"],
            "scan_command": s["command"],
        }
        if s["count_time_ms"] is not None:
            meta["count_time_ms"] = s["count_time_ms"]
        if s["umi"]:
            meta["umi"] = s["umi"]
        if s["motor_positions"]:
            meta["motor_positions"] = s["motor_positions"]

    if file_meta.get("beamline_file"):
        meta["beamline_file"] = file_meta["beamline_file"]
    if file_meta.get("comment"):
        meta["comment"] = file_meta["comment"]

    return meta


def _parse_spec_file(path_to_file):
    """Parse a SPEC .dat file.

    Returns:
        file_meta (dict): File-level metadata from #F, #E, #D, #C lines.
        scans (list of dict): One dict per scan with keys: number, command,
            tstamp, count_time_ms, motor_positions, umi, columns, data.
    """
    with open(path_to_file, "r") as f:
        lines = f.readlines()

    file_meta: dict = {}
    motor_names: list = []
    motor_name_lines: dict = {}  # {index: [name, ...]} from #O lines
    scans = []
    current = None
    data_lines: list = []
    file_tstamp = None
    in_file_header = True

    for line in lines:
        line = line.rstrip("\n")

        # ---- file-level header (before the first #S) ----
        if in_file_header:
            if line.startswith("#F "):
                file_meta["beamline_file"] = line[3:].strip()

            elif line.startswith("#E "):
                file_tstamp = float(line[3:].strip())

            elif line.startswith("#D ") and file_tstamp is None:
                try:
                    file_tstamp = timestamp_string_to_tstamp(
                        line[3:].strip(), forms=(SPEC_DATE_FORM,)
                    )
                except Exception:
                    pass

            elif line.startswith("#C "):
                # collect all comment lines, joining with newline
                existing = file_meta.get("comment", "")
                new = line[3:].strip()
                file_meta["comment"] = (
                    (existing + "\n" + new).strip() if existing else new
                )

            else:
                m = re.match(r"#O(\d+) ", line)
                if m:
                    idx = int(m.group(1))
                    names = re.split(r"  +", line[line.index(" ") + 1 :].strip())
                    motor_name_lines[idx] = names

        # ---- start of a new scan block ----
        if line.startswith("#S "):
            in_file_header = False
            if current is not None:
                current["data"] = _arrays_from_lines(current["columns"], data_lines)
                scans.append(current)
            if not motor_names and motor_name_lines:
                for i in sorted(motor_name_lines):
                    motor_names.extend(motor_name_lines[i])
            parts = line[3:].split(None, 1)
            current = {
                "number": int(parts[0]),
                "command": parts[1].strip() if len(parts) > 1 else "",
                "tstamp": file_tstamp,
                "count_time_ms": None,
                "motor_positions": {},
                "umi": [],
                "columns": [],
                "_p_values": {},
            }
            data_lines = []

        elif current is None:
            continue

        elif line.startswith("#D "):
            try:
                current["tstamp"] = timestamp_string_to_tstamp(
                    line[3:].strip(), forms=(SPEC_DATE_FORM,)
                )
            except Exception:
                pass

        elif line.startswith("#T "):
            try:
                current["count_time_ms"] = float(line[3:].split()[0])
            except (ValueError, IndexError):
                pass

        elif line.startswith("#UMI"):
            text = line[line.index(" ") + 1 :].strip() if " " in line else ""
            if text:
                current["umi"].append(text)

        elif line.startswith("#L "):
            current["columns"] = re.split(r"  +", line[3:].strip())

        elif not line.startswith("#") and line.strip():
            data_lines.append(line)

        else:
            m = re.match(r"#P(\d+) ", line)
            if m:
                idx = int(m.group(1))
                try:
                    vals = [float(v) for v in line[line.index(" ") + 1 :].split()]
                    current["_p_values"][idx] = vals
                except ValueError:
                    pass

    if current is not None:
        current["data"] = _arrays_from_lines(current["columns"], data_lines)
        scans.append(current)

    # resolve motor positions using the flat motor_names list
    for scan in scans:
        p = scan.pop("_p_values")
        if motor_names and p:
            flat_vals: list = []
            for i in sorted(p):
                flat_vals.extend(p[i])
            scan["motor_positions"] = {
                name: flat_vals[i]
                for i, name in enumerate(motor_names)
                if i < len(flat_vals)
            }

    return file_meta, scans


def _arrays_from_lines(columns, data_lines):
    if not columns or not data_lines:
        return {}
    arr = np.array([row.split() for row in data_lines], dtype=float)
    return {col: arr[:, i] for i, col in enumerate(columns) if i < arr.shape[1]}
