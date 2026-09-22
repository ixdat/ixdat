"""Reader for Ocean Insight / OceanView time-series spectral exports."""

from pathlib import Path
import re
from datetime import datetime, timedelta, timezone, tzinfo as TzInfo

import numpy as np
from scipy.ndimage import uniform_filter1d
import warnings

from ..data_series import DataSeries, TimeSeries, Field
from ..spectra import SpectrumSeries
from ..techniques.spectroelectrochemistry import OpticalSpectrumSeries

TZ_ABBREVIATIONS = {
    "UTC": timezone.utc,
    "GMT": timezone.utc,
    "WET": timezone.utc,
    "BST": timezone(timedelta(hours=1)),  # British Summer Time
    "WEST": timezone(timedelta(hours=1)),
    "CET": timezone(timedelta(hours=1)),
    "CEST": timezone(timedelta(hours=2)),
    "EET": timezone(timedelta(hours=2)),
    "EEST": timezone(timedelta(hours=3)),
}

# A clock time alone cannot anchor an absolute timestamp, and would jump
# backwards across midnight, so a full date is required.
ROW_DATETIME_FORMATS = ("%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S")

# OceanView stamps the data rows with 1970-01-01 when it has no real date.
# It parses as a valid datetime, so it has to be rejected explicitly.
PLACEHOLDER_YEAR = 1970


class OceanViewTimeSeriesReader:
    """Reader for Ocean Insight/OceanView 'Data from ...txt Node' exports.

    Produces an OpticalSpectrumSeries with axes = (time, wavelength).

    Parameters (via read):
        path_to_file: str | Path
        name: optional name, defaults to stem
        spectra_type: "Transmission", "Absorption", or None (Intensity)
        boxcar_width: boxcar width for data smoothing. Width is 1 by default.
        tstamp_source: where to take the series' starting timestamp from,
            either "header" (default) or "spectrum". See below.
        assume_timezone: the timezone to use when the file does not state one
            this reader can interpret. See below.
        average_every: average this many consecutive rows into one spectrum
            as they are read, to reduce the size of the resulting series.
            1 (default) keeps every row unaveraged. A group of fewer than
            `average_every` rows left over at the end of the file is still
            averaged and kept, not dropped.
        cls: target SpectrumSeries subclass (defaults to OpticalSpectrumSeries)

    On ``tstamp_source``:

        "header" (default) anchors the time axis to the ``Date:`` line,
        refined with the millisecond field of the filename if one is present.
        A file with no ``Date:`` line has no header timestamp, and says so.

        "spectrum" anchors it to the first data row's own timestamp. OceanView
        writes the header 1-3 s before acquisition starts, so the header value
        runs early; that matters when aligning against another instrument.
        Rows carrying the 1970-01-01 placeholder cannot be used, and say so.

    On timezones:

        Neither the header date nor the data rows state a timezone, so one is
        resolved for both, in this order:

        1. a numeric offset in the header, e.g. ``GMT+02:00``
        2. a recognised abbreviation in the header, e.g. ``CEST``
        3. ``assume_timezone``, passed by the caller
        4. otherwise a ValueError, rather than a guess

        ``assume_timezone`` takes an offset (``"+02:00"``, ``"-05:00"``,
        ``"Z"``) or any ``datetime.tzinfo``.
    """

    def read(
        self,
        path_to_file,
        name=None,
        spectra_type="Intensity",  # spectra type is Intensity by default
        boxcar_width=1,  # boxcar width for data smoothing. Width is 1 by default.
        tstamp_source="header",
        assume_timezone=None,
        average_every=1,
        cls=OpticalSpectrumSeries,
    ):
        path_to_file = Path(path_to_file)
        name = name or path_to_file.stem

        if tstamp_source not in ("header", "spectrum"):
            raise ValueError(
                f"tstamp_source must be 'header' or 'spectrum', not {tstamp_source!r}"
            )

        if not issubclass(cls, SpectrumSeries):
            cls = OpticalSpectrumSeries
            
        dt_header, tz_text = None, None
        start_idx=None
        spectra = []
        rel_times = []
        row_datetimes = []
        rel_time_sum = 0.0
        count = 0

        with open(path_to_file, encoding="utf-8", errors="ignore") as f:
           
            for i,ln in enumerate(f):  
            # ---- Parse header for Date ----
            # Read whatever tstamp_source says, so a bad header is reported in
            # either mode. The date stays naive; the timezone is resolved below.
                if i<40 and not start_idx and ln.lower().startswith("date:"):
                    date_str = ln.split(":", 1)[1].strip()
                    dt_header, tz_text = self._parse_header_date(date_str)
                    # The data rows state no timezone either, so one is used for both.
                    tzinfo = self._resolve_timezone(tz_text, assume_timezone, path_to_file)
                    # ---- Refine with filename milliseconds if available ----
                    # OceanView truncates the header to whole seconds.
                    ms = self._parse_filename_time(path_to_file)  # int milliseconds
                    if ms is not None and dt_header is not None:
                        dt_header = dt_header.replace(microsecond=ms * 1000)
                    continue
                
                # ---- Find Begin Spectral Data and parse wavelengths ----
                if i<40 and re.search(r"begin\s+spectral\s+data", ln, flags=re.I):
                    start_idx = i
                    wl_line=next(f).strip()
                    if not re.match(r"^\s*\d", wl_line):  # next line after begin spectral data should start with a number
                        raise ValueError(
                            f"Expected wavelength line after Begin Spectral Data, "
                            f"but got: {wl_line!r}"
                            )
                    wavelengths = self._parse_float_row(wl_line)
                    if wavelengths.size == 0:
                        raise ValueError("OceanView: wavelength line is empty or malformed")
                    n_wavelengths = len(wavelengths)
                    spectra_sum = np.zeros(n_wavelengths, dtype=np.float64)
                    continue
                if i >=40 and start_idx is None:
                    raise ValueError("No spectral data section found!")
                    
               
                 # ---- Parse spectra and relative times ----
                 # Rows are averaged in groups of `average_every` as they are parsed,
                 # so the full-resolution matrix never has to be held in memory. The
                 # relative time of an averaged row is the mean of its rows' offsets
                 # from the first row -- computed from the parsed datetimes, not from
                 # a bare time-of-day, so it stays correct across midnight.
                
                if start_idx and i >= start_idx+2 and ln.strip():
                    # Robust handling of separators
                    stamp_str, vals = self._split_stamp(ln)

                    row_datetime = self._parse_row_datetime(stamp_str)
                    
                    if row_datetime is None:
                        raise ValueError(
                            f"Row {i} of the spectral data in "
                            f"{path_to_file.name} is stamped {stamp_str!r}, which is "
                            "not a full date and time. This reader needs a stamp like "
                            "'2026-08-05 18:51:45.681484' on every row, so that "
                            "relative times stay correct across midnight."
                        )
                    row_datetimes.append(row_datetime)
                    
                    vals = [float(v.replace(",", ".")) for v in vals.split()]

                    if len(vals) < n_wavelengths:
                        raise ValueError(
                            f"Row {i} of the spectral data in "
                            f"{path_to_file.name} has {len(vals)} values, but the "
                            f"wavelength axis has {n_wavelengths}."
                        )
                    # Take the last N, so an extra leading column cannot shift the
                    # spectrum against the wavelength axis.
                    vals = vals[-n_wavelengths:]

                    # Naive difference: no timezone needed, correct across midnight.
                    rel_sec = (row_datetime - row_datetimes[0]).total_seconds()

                    spectra_sum += vals
                    rel_time_sum += rel_sec
                    count += 1

                    if count == average_every:
                        spectra.append(spectra_sum / average_every)
                        rel_times.append(rel_time_sum / average_every)
                        # Reset
                        spectra_sum.fill(0)
                        rel_time_sum = 0.0
                        count = 0
                if count > 0: # This averages the last (<average every) lines of data
                    spectra.append(spectra_sum / count)
                    rel_times.append(rel_time_sum / count)
         
        if not spectra:
            raise ValueError(
                f"OceanView: {path_to_file.name} has a spectral data section "
                "but no spectra in it"
            )
                                                   
        y_matrix = np.stack(spectra)
        
        rel_times = np.array(rel_times)
        
        # ---- Apply smoothing ----
        
        y_matrix_smoothed = uniform_filter1d(
            y_matrix,
            size=boxcar_width,
            axis=1,
        )      

        # ---- Resolve the starting timestamp ----
        tstamp_header = dt_header.replace(tzinfo=tzinfo).timestamp()   
        tstamp_first = row_datetimes[0].replace(tzinfo=tzinfo).timestamp()
        if tstamp_first < 24 * 60 * 60:    # epoch in any timezone
            if tstamp_source == "spectrum":
                raise ValueError("Timestamp is incorrect")  # the demo script tests for this.
                warnings.warn(
                    f"{path_to_file} saved without correct absolute times. Using the header"
                    " as an approximation, but alignment may be a few seconds off."
                )
            t_zero = tstamp_header
        elif tstamp_source == "header":
            t_zero = tstamp_header
            rel_times += tstamp_first - tstamp_header
        else:
            t_zero = tstamp_first

        # ---- Wrap into ixdat objects ----
        xseries = DataSeries(name="wavelength", unit_name="nm", data=wavelengths)
        # The same tstamp goes on the TimeSeries and on the SpectrumSeries
        # below; if they disagree, ixdat shifts the relative times to match.
        tseries = TimeSeries(
            name="time", unit_name="s", data=rel_times, tstamp=t_zero
        )

        if spectra_type == "Transmission":
            field = Field(
                name="transmission",
                unit_name="a.u",
                data=y_matrix_smoothed,
                axes_series=[tseries, xseries],
            )
        elif spectra_type == "Absorption":
            field = Field(
                name="absorption",
                unit_name="a.u.",
                data=y_matrix_smoothed,
                axes_series=[tseries, xseries],
            )
        else:
            field = Field(
                name="intensity",
                unit_name="a.u.",
                data=y_matrix_smoothed,
                axes_series=[tseries, xseries],
            )

        uvvis_series = cls(
            name=name,
            reader=self,
            technique="Optical",
            tstamp=tstamp_first,
            field=field,
            continuous=True,
            spectra_type=spectra_type,
        )
        return uvvis_series

   
    # -------- Timezone --------
    def _resolve_timezone(self, tz_text, assume_timezone, path_to_file):
        """Return the tzinfo for this file's naive times, or None if unknown.

        The file wins; `assume_timezone` fills the gap. None is not fatal
        here -- only once a naive time actually has to be localised.
        """
        assumed = self._parse_assume_timezone(assume_timezone)

        if tz_text is None:
            return assumed

        from_header = self._timezone_from_text(tz_text)
        if from_header is not None:
            return from_header
        if assumed is not None:
            return assumed
        raise ValueError(
            f"The header of {path_to_file.name} gives its timezone as "
            f"{tz_text!r}, which this reader cannot interpret unambiguously. "
            'Pass the timezone explicitly, for example assume_timezone="+02:00".'
        )

    @staticmethod
    def _timezone_from_text(tz_text):
        """Return a tzinfo for a header timezone token, or None if unknown."""
        m = re.fullmatch(r"GMT([+-])(\d{1,2})(?::?(\d{2}))?", tz_text)
        if m:
            sign = -1 if m.group(1) == "-" else 1
            hours, minutes = int(m.group(2)), int(m.group(3) or 0)
            return timezone(sign * timedelta(hours=hours, minutes=minutes))

        if tz_text in TZ_ABBREVIATIONS:
            return TZ_ABBREVIATIONS[tz_text]
        return None

    @staticmethod
    def _parse_assume_timezone(value):
        """Return a tzinfo for the caller's assume_timezone, or None."""
        if value is None:
            return None
        if isinstance(value, TzInfo):
            return value
        if isinstance(value, str):
            text = value.strip()
            if text.upper() in ("Z", "UTC", "GMT"):
                return timezone.utc
            m = re.fullmatch(r"([+-])(\d{1,2}):?(\d{2})?", text)
            if m:
                sign = -1 if m.group(1) == "-" else 1
                hours, minutes = int(m.group(2)), int(m.group(3) or 0)
                return timezone(sign * timedelta(hours=hours, minutes=minutes))
        raise ValueError(
            f"assume_timezone={value!r} is not a timezone this reader "
            'understands. Pass an offset such as "+02:00", "-05:00" or "Z", '
            "or a datetime.tzinfo such as zoneinfo.ZoneInfo('Europe/Copenhagen')."
        )

    @staticmethod
    def _localise(dt_naive, tzinfo, what, path_to_file):
        """Return the unix timestamp of a naive datetime, or explain why not."""
        if tzinfo is None:
            raise ValueError(
                f"{what} of {path_to_file.name} has no timezone: the file does "
                "not state one this reader can interpret, and none was passed. "
                "The timestamp cannot be placed on an absolute time axis. Pass "
                'the timezone explicitly, for example assume_timezone="+02:00".'
            )
        return dt_naive.replace(tzinfo=tzinfo).timestamp()

    # -------- Helpers --------
    @staticmethod
    def _parse_float_row(s):
        parts = re.split(r"[\s\t]+", s.strip())
        return np.array([float(p.replace(",", ".")) for p in parts if p], dtype=float)

    @staticmethod
    def _parse_header_date(date_str):
        """Parse an OceanView header date into (naive datetime, timezone text).

        Handles 'Mon Aug 18 15:23:28 CEST 2025' and
        'Mon Aug 18 15:23:28 GMT+02:00 2025'. The timezone is left as text so
        that "none found" stays distinguishable from "found and applied".
        """
        s = date_str.strip()
        if s.lower().startswith("date:"):
            s = s.split(":", 1)[1].strip()
        s = re.sub(r"\s+", " ", s)

        # A numeric offset is explicit, so it wins over an abbreviation.
        # OceanView is Java: 'GMT+02:00' means UTC+2, the opposite sign to
        # the POSIX TZ convention.
        tz_text = None
        tz_match = re.search(r"\bGMT[+-]\d{1,2}(?::?\d{2})?\b", s)
        if tz_match:
            tz_text = tz_match.group(0)
        else:
            # An abbreviation sits between the clock time and the year.
            tz_match = re.search(r"\b([A-Z]{2,5})\b(?=\s+\d{4}$)", s)
            if tz_match:
                tz_text = tz_match.group(1)

        if tz_match:
            tz_start, tz_end = tz_match.span()
            s = re.sub(r"\s+", " ", (s[:tz_start] + s[tz_end:]).strip())

        fmts = [
            "%a %b %d %H:%M:%S %Y",
            "%b %d %H:%M:%S %Y",
            "%a %b %d %H:%M:%S.%f %Y",
            "%b %d %H:%M:%S.%f %Y",
        ]
        for fmt in fmts:
            try:
                return datetime.strptime(s, fmt), tz_text
            except ValueError:
                continue
        raise ValueError(f"Could not parse header date: {date_str!r}")

    def _split_stamp(self, line):
        """Split a data line into (timestamp, rest).

        Prefer tab if present; else first whitespace run.
        """
        if "\t" in line:
            stamp, rest = line.split("\t", 1)
        else:
            parts = re.split(r"\s+", line.strip(), maxsplit=1)
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
    def _parse_row_datetime(stamp):
        """Parse a naive datetime from a data-row stamp, or None if it has no date."""
        if not stamp:
            return None
        for fmt in ROW_DATETIME_FORMATS:
            try:
                return datetime.strptime(stamp.strip(), fmt)
            except ValueError:
                continue
        return None
