from pathlib import Path
import numpy as np
import re
from datetime import datetime, timedelta, timezone

from ..data_series import DataSeries, TimeSeries, Field
from ..spectra import SpectrumSeries
from ..techniques.spectroelectrochemistry import OpticalSpectrumSeries


class OceanViewTimeSeriesReader:
    """Reader for Ocean Insight/OceanView 'Data from ...txt Node' exports.

    Produces an OpticalSpectrumSeries with axes = (time, wavelength).

    Parameters (via read):
        path_to_file: str | Path
        name: optional name, defaults to stem
        cls: target SpectrumSeries subclass (defaults to OpticalSpectrumSeries)
        anchor: "header" (default) | "mtime"
            - Which absolute timestamp to use if both exist.
        assume_tz: e.g., "UTC" (default). Only used if header lacks a TZ;
            we anchor the parsed naive datetime to this zone instead of UTC.
        filename_ms: bool (default True)
            - If True, extract milliseconds from pattern '__HH-MM-SS-mmm' and
              inject into the header time (if available).
    """

    def read(
        self,
        path_to_file,
        name=None,
        cls=OpticalSpectrumSeries,
    ):
        path_to_file = Path(path_to_file)
        name = name or path_to_file.stem

        if not issubclass(cls, SpectrumSeries):
            cls = OpticalSpectrumSeries

        with open(path_to_file, encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        # ---- Parse header for Date ----
        dt_header = None
        for ln in lines[:40]:  # only scan the top of the file
            if ln.lower().startswith("date:"):
                date_str = ln.split(":", 1)[1].strip()
                dt_header = self._parse_header_date(date_str)
                break

        # ---- Refine with filename milliseconds if available ----
        ms = self._parse_filename_time(path_to_file)  # int 毫秒
        if ms is not None and dt_header is not None:
            dt_header = dt_header.replace(microsecond=ms * 1000)

        tstamp_first = dt_header.timestamp()

        # ---- Find Begin Spectral Data ----
        start_idx = None
        for i, ln in enumerate(lines):
            if re.search(r"begin\s+spectral\s+data", ln, flags=re.I):
                if re.match(r"^\s*\d", lines[i + 1]):  # next line starts with a number
                    start_idx = i
                    break
        if start_idx is None:
            raise ValueError("No spectral data section found!")

        # Validate/Locate wavelength
        wl_line = lines[start_idx + 1].strip()
        wavelengths = self._parse_float_row(wl_line)
        if wavelengths.size == 0:
            raise ValueError("OceanView: wavelength line is empty or malformed")

        data_lines = [ln for ln in lines[start_idx + 2 :] if ln.strip()]

        # ---- Parse spectra and relative times ----
        spectra = []
        rel_times = []
        for ln in data_lines:

            # Robust handling of seperators
            stamp_str, vals = self._split_stamp(ln)

            rel_sec = self._parse_row_time(stamp_str)
            rel_times.append(rel_sec)
            vals = [float(v.replace(",", ".")) for v in vals.split()]

            if len(vals)>=len(wavelengths):
                spectra.append(vals[len(vals)-len(wavelengths):])
            else:
                spectra.append(vals)

        y_matrix = np.stack(spectra)
        rel_times = np.array(rel_times) - rel_times[0]  # start at 0 s

        # ---- Wrap into ixdat objects ----
        xseries = DataSeries(name="wavelength", unit_name="nm", data=wavelengths)
        tseries = TimeSeries(
            name="time", unit_name="s", data=rel_times, tstamp=tstamp_first
        )

        field = Field(
            name="intensity",
            unit_name="a.u.",
            data=y_matrix,
            axes_series=[tseries, xseries],
        )

        uvvis_series = cls(
            name=name,
            reader=self,
            technique="Optical",
            tstamp=tstamp_first,
            field=field,
            continuous=True,
        )
        return uvvis_series

    # -------- Helpers --------
    @staticmethod
    def _parse_float_row(s):
        parts = re.split(r"[\s\t]+", s.strip())
        return np.array([float(p.replace(",", ".")) for p in parts if p], dtype=float)

    @staticmethod
    def _parse_header_date(date_str: str) -> datetime:

        s = date_str.strip()
        if s.lower().startswith("date:"):
            s = s.split(":", 1)[1].strip()
        s = re.sub(r"\s+", " ", s)


        s = re.sub(r"\b[A-Z]{2,5}\b(?=\s+\d{4}$)", "", s)
        s = re.sub(r"\bGMT[+-]\d{1,2}(?::\d{2})?\b", "", s).strip()

        fmts = [
            "%a %b %d %H:%M:%S %Y",
            "%b %d %H:%M:%S %Y",
            "%a %b %d %H:%M:%S.%f %Y",
            "%b %d %H:%M:%S.%f %Y",
        ]
        for fmt in fmts:
            try:
                return datetime.strptime(s, fmt)  
            except ValueError:
                continue
        raise ValueError(f"Could not parse header date: {date_str!r}")

    def _split_stamp(self, line):
        """Split a data line into (timestamp, rest). Prefer tab if present; else first whitespace run."""
        if "\t" in line:
            stamp, rest = line.split("\t", 1)
        else:
            parts = re.split(r"[\s\t]", line.strip())
            if not parts:
                raise ValueError("empty line")
            stamp = parts[0]
            rest = parts[1] if len(parts) > 1 else ""
        return stamp.strip(), rest.strip()

    @staticmethod
    def _parse_filename_time(path):
        """Extract only millisecond from '__HH-MM-SS-mmm'."""
        m = re.search(r"__(\d{2})-(\d{2})-(\d{2})-(\d{3})", str(path))
        if not m:
            return None
        ms = int(m.group(4))
        return ms

    @staticmethod
    def _parse_row_time(stamp):
        # Row times like '1970-01-01 01:24:29.367452' or '13-42-52-946'
        stamp = stamp.strip()
        fmts = [
            "%Y-%m-%d %H:%M:%S.%f",
            "%H:%M:%S.%f",
            "%Y-%m-%d %H:%M:%S",
            "%H:%M:%S",
            "%H-%M-%S-%f",  # support '13-42-52-946'
        ]
        for fmt in fmts:
            try:
                dt = datetime.strptime(stamp, fmt)
                return dt.hour * 3600 + dt.minute * 60 + dt.second + dt.microsecond / 1e6
            except ValueError:
                continue
        try:
            return float(stamp.replace(",", "."))
        except Exception:
            return np.nan
