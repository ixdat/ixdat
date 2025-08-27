
from pathlib import Path
import numpy as np
import re
from datetime import datetime, timedelta, timezone

from ..data_series import DataSeries, TimeSeries, Field
from ..spectra import SpectrumSeries
from ..techniques.spectroelectrochemistry import OpticalSpectrumSeries

class OceanViewTimeSeriesReader:
    """Reader for Ocean Insight/OceanView 'Data from ...txt Node' exports"""

    def read(
        self,
        path_to_file,
        name=None,
        cls=OpticalSpectrumSeries,
        suffix=".txt",
    ):
        path_to_file = Path(path_to_file)
        name = name or path_to_file.stem

        if not issubclass(cls, SpectrumSeries):
            cls = OpticalSpectrumSeries

        with open(path_to_file, encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        # ---- Parse header for Date ----
        tstamp_first = None
        for ln in lines[:40]:  # only scan the top of the file
            if ln.lower().startswith("date:"):
                date_str = ln.split(":", 1)[1].strip()
                print(date_str)
                tstamp_first = self._parse_header_date(date_str)
                break

        # ---- Find Begin Spectral Data ----
        start_idx = None
        for i, ln in enumerate(lines):
            if "Begin Spectral Data" in ln:
                start_idx = i
                break
        if start_idx is None:
            raise ValueError("No spectral data section found!")

        wl_line = lines[start_idx + 1].strip()
        wavelengths = np.array(
            [float(re.sub(",", ".", w)) for w in wl_line.split() if w]
        )

        data_lines = [ln for ln in lines[start_idx + 2 :] if ln.strip()]

        # ---- Parse spectra and relative times ----
        spectra = []
        rel_times = []
        for j, ln in enumerate(data_lines):
            try:
                stamp_str, vals = ln.split("\t", 1)
            except ValueError:
                parts = ln.split()
                stamp_str, vals = parts[0], " ".join(parts[1:])

            rel_sec = self._parse_row_time(stamp_str)
            rel_times.append(rel_sec)
            vals = [float(v.replace(",", ".")) for v in vals.split()]
            spectra.append(vals)

        y_matrix = np.stack(spectra)
        rel_times = np.array(rel_times) - rel_times[0]  # start at 0 s
        print(tstamp_first)
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
    def _parse_header_date(date_str: str) -> float:
        """Parse OceanView-style header dates like
        'Mon Aug 18 15:23:28 CEST 2025' or 'Mon Aug 18 15:23:28 GMT+2 2025'
        and return a UNIX timestamp (float seconds).
        """
        s = date_str.strip()
        # If the caller passed the whole line ('Date: ...'), trim the label:
        if s.lower().startswith("date:"):
            s = s.split(":", 1)[1].strip()

        # Normalize all whitespace to single spaces:
        s = re.sub(r"\s+", " ", s)

        # Known TZ offsets (hours). Add more if you need them.
        tz_offsets = {
            "UTC": 0, "GMT": 0,
            "CET": 1, "CEST": 2,
            "WET": 0, "WEST": 1,
            "EET": 2, "EEST": 3,
            "PST": -8, "PDT": -7,
            "MST": -7, "MDT": -6,
            "CST": -6, "CDT": -5,
            "EST": -5, "EDT": -4,
            "BST": 1,   # UK summer
            "IST": 5.5, # India (note: ambiguous name globally)
        }

        tz_offset_hours = None

        # 1) GMT±H[:MM] style, e.g. 'GMT+2' or 'GMT-05:30'
        m = re.search(r"\bGMT([+-])(\d{1,2})(?::(\d{2}))?\b", s)
        if m:
            sign = 1 if m.group(1) == "+" else -1
            hours = int(m.group(2))
            minutes = int(m.group(3) or 0)
            tz_offset_hours = sign * (hours + minutes / 60.0)
            # Remove the GMT token so strptime can match:
            s = re.sub(r"\s*GMT[+-]\d{1,2}(?::\d{2})?\s*", " ", s).strip()

        # 2) Abbrev before year, e.g. '... CEST 2025'
        if tz_offset_hours is None:
            m = re.search(r"\b([A-Z]{2,5})\b(?=\s+\d{4}$)", s)
            if m and m.group(1) in tz_offsets:
                tz_offset_hours = tz_offsets[m.group(1)]
                # remove the TZ token so strptime formats will match:
                s = s.replace(" " + m.group(1), "")

        # Try a few likely layouts (with/without weekday, with/without micros)
        fmts = [
            "%a %b %d %H:%M:%S %Y",
            "%b %d %H:%M:%S %Y",
            "%a %b %d %H:%M:%S.%f %Y",
            "%b %d %H:%M:%S.%f %Y",
        ]
        dt = None
        for fmt in fmts:
            try:
                dt = datetime.strptime(s, fmt)
                break
            except ValueError:
                pass

        if dt is None:
            raise ValueError(f"Could not parse header date: {date_str!r}")

        # Apply timezone (default UTC if none found)
        if tz_offset_hours is None:
            tz = timezone.utc
        else:
            tz = timezone(timedelta(hours=tz_offset_hours))

        return dt.replace(tzinfo=tz).timestamp()

    @staticmethod
    def _parse_row_time(stamp):
        # Row times like '1970-01-01 01:24:29.367452'
        stamp = stamp.strip()
        for fmt in ("%Y-%m-%d %H:%M:%S.%f", "%H:%M:%S.%f",
                    "%Y-%m-%d %H:%M:%S", "%H:%M:%S"):
            try:
                dt = datetime.strptime(stamp, fmt)
                return dt.hour * 3600 + dt.minute * 60 + dt.second + dt.microsecond/1e6
            except ValueError:
                continue
        # fallback: try float seconds
        try:
            return float(stamp.replace(",", "."))
        except Exception:
            return np.nan
