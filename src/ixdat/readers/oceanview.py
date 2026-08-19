from pathlib import Path
import numpy as np
import re
from datetime import datetime, timedelta, timezone
from scipy.ndimage import uniform_filter1d

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

    """

    def read(
        self,
        path_to_file,
        name=None,
        spectra_type=None,  # spectra type is Intensity by default
        boxcar_width=None,  # boxcar width for data smoothing. Width is 1 by default.
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
        ms = self._parse_filename_time(path_to_file)  # int milliseconds
        if ms is not None and dt_header is not None:
            dt_header = dt_header.replace(microsecond=ms * 1000)

        tstamp_first = (
            dt_header.timestamp() if dt_header else path_to_file.stat().st_mtime
        )

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

            # Robust handling of separators
            stamp_str, vals = self._split_stamp(ln)

            rel_sec = self._parse_row_time(stamp_str)
            rel_times.append(rel_sec)
            vals = [float(v.replace(",", ".")) for v in vals.split()]

            if len(vals) >= len(wavelengths):
                spectra.append(vals[len(vals) - len(wavelengths) :])
            else:
                spectra.append(vals)

        y_matrix = np.stack(spectra)

        # ---- Apply smoothing ----
        if boxcar_width is None:
            boxcar_width = 1

        y_matrix_smoothed = uniform_filter1d(
            y_matrix,
            size=boxcar_width,
            axis=1,
        )

        rel_times = np.array(rel_times) - rel_times[0]  # start at 0 s

        # ---- Wrap into ixdat objects ----
        xseries = DataSeries(name="wavelength", unit_name="nm", data=wavelengths)
        tseries = TimeSeries(
            name="time", unit_name="s", data=rel_times, tstamp=tstamp_first
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

        tzinfo = None
        tz_match = re.search(r"\b([A-Z]{2,5})\b(?=\s+\d{4}$)", s)
        if tz_match:
            tzinfo = {
                "UTC": timezone.utc,
                "GMT": timezone.utc,
                "CET": timezone(timedelta(hours=1)),
                "CEST": timezone(timedelta(hours=2)),
            }.get(tz_match.group(1))
            s = s[: tz_match.start()] + s[tz_match.end() :]
        s = re.sub(r"\bGMT[+-]\d{1,2}(?::\d{2})?\b", "", s).strip()

        fmts = [
            "%a %b %d %H:%M:%S %Y",
            "%b %d %H:%M:%S %Y",
            "%a %b %d %H:%M:%S.%f %Y",
            "%b %d %H:%M:%S.%f %Y",
        ]
        for fmt in fmts:
            try:
                dt = datetime.strptime(s, fmt)
                return dt.replace(tzinfo=tzinfo or timezone.utc)
            except ValueError:
                continue
        raise ValueError(f"Could not parse header date: {date_str!r}")

    def _split_stamp(self, line):
        """Split a data line into (timestamp, rest). Prefer tab if present; else first whitespace run."""
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

class OceanViewTimeSeriesReader_per_spectrum_t(OceanViewTimeSeriesReader):
    """Same as OceanViewTimeSeriesReader, but anchors the returned
    SpectrumSeries' time axis to the spectrometer's own recorded
    acquisition time of the *first spectrum*, instead of the file header's
    declared `Date:` time.

    Why this matters: OceanView writes the file header a beat before the
    first spectrum is actually acquired -- typically 1-3 seconds, the time
    it takes to open the file and start capturing -- so
    `OceanViewTimeSeriesReader`'s header-based tstamp systematically
    understates when the data was actually taken (confirmed here: this
    file's header says 18:51:43, but the first spectral data row's own
    timestamp is 18:51:45.68, a 2.68s gap). That's usually harmless on its
    own, but it matters when aligning against another instrument's
    independently-recorded timestamps -- e.g. a potentiostat's `.nox` file
    via `NovaNoxReader` -- since the resulting offset can be large enough,
    relative to a fast transient, to make an optical response appear to
    *precede* the electrochemical event that caused it.

    Requires the data rows to carry a full absolute datetime (e.g.
    "2026-08-05 18:51:45.681484"), not just a time-of-day -- true for the
    QEP/OceanView exports this is designed for. If the first row's stamp
    can't be parsed as a full datetime, this warns and falls back to the
    header-based tstamp (i.e. behaves exactly like the parent class).
    """

    def read(self, path_to_file, name=None, cls=OpticalSpectrumSeries):
        spectrum_series = super().read(path_to_file, name=name, cls=cls)

        true_tstamp = self._true_first_spectrum_tstamp(path_to_file)
        if true_tstamp is None:
            warnings.warn(
                f"{type(self).__name__}: could not parse a full datetime "
                "('YYYY-MM-DD HH:MM:SS.ffffff') from the first spectral "
                f"data row of {path_to_file} -- falling back to the header "
                "Date: timestamp, same as OceanViewTimeSeriesReader."
            )
            return spectrum_series

        # Mutate the *nested* TimeSeries inside the field, not the outer
        # `spectrum_series.tstamp` attribute. `SpectrumSeries.field` is a
        # lazy property (see ixdat.spectra.Spectrum.field) that re-anchors
        # the field's time axis to `self.tstamp` whenever the two disagree
        # by more than `t_tolerance`, via `time_shifted(...)`, which
        # preserves absolute time (tstamp + data is invariant under the
        # shift). So it doesn't matter that the outer `.tstamp` still holds
        # the (wrong) header value: any later reconciliation just shifts
        # `data` to compensate, carrying this fix's correct absolute times
        # through. Setting the outer attribute instead would do nothing
        # until the next reconciliation, and even then would only relabel
        # t=0 without correcting the underlying absolute times.
        spectrum_series.field.axes_series[0].tstamp = true_tstamp
        return spectrum_series

    @staticmethod
    def _true_first_spectrum_tstamp(path_to_file):
        """Return the Unix timestamp parsed directly from the first
        spectral data row's own absolute timestamp column
        ('YYYY-MM-DD HH:MM:SS.ffffff'), or None if the rows don't carry a
        full datetime (e.g. time-of-day only, with no date part)."""
        with open(path_to_file, encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        start_idx = None
        for i, ln in enumerate(lines):
            if re.search(r"begin\s+spectral\s+data", ln, flags=re.I):
                start_idx = i
                break
        if start_idx is None:
            return None

        for ln in lines[start_idx + 2:]:
            if not ln.strip():
                continue
            stamp_str = ln.split("\t", 1)[0].strip()
            try:
                return datetime.strptime(stamp_str, "%Y-%m-%d %H:%M:%S.%f").timestamp()
            except ValueError:
                return None
        return None