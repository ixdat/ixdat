from pathlib import Path
import numpy as np
import re
from datetime import datetime

from ..data_series import DataSeries, TimeSeries, Field
from ..spectra import SpectrumSeries
from ..techniques.spectroelectrochemistry import OpticalSpectrumSeries


class AndorKineticsCSVReader:
    """
    Reader for Andor/Solis/Kymera Kinetics CSV exports.

    Time axis:
        rel_time = start_offset
                   + AccumulateCycleTime * NumAccum
                   + KineticCycleTime * frame_index

        start_offset = 0  (we subtract rel_time[0] to make time start at 0)

    """

    HEADER_SCAN_LINES = 120

    def read(self, path_to_file, name=None, cls=OpticalSpectrumSeries):
        path = Path(path_to_file)
        name = name or path.stem

        if not issubclass(cls, SpectrumSeries):
            cls = OpticalSpectrumSeries

        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        # ---------------- HEADER FIELDS ----------------
        accum_cycle = self._get_float(lines, "Accumulate Cycle Time")
        kinetic_cycle = self._get_float(lines, "Kinetic Cycle Time")
        num_accum = self._get_int(lines, "Number of Accumulations")
        n_frames = self._get_int(lines, "Number in Kinetics Series")

        if accum_cycle is None or kinetic_cycle is None or num_accum is None:
            raise ValueError("Missing required timing metadata in header")

        # datetime (for absolute time anchor)
        dt_header = self._parse_datetime_header(lines)
        tstamp_first = dt_header.timestamp() if dt_header else None

        # ---------------- FIND DATA START ----------------
        data_start = None
        for i, ln in enumerate(lines[:self.HEADER_SCAN_LINES]):
            if self._row_starts_with_float(ln):
                data_start = i
                break
        if data_start is None:
            raise ValueError("Cannot locate data block.")

        # ---------------- READ MATRIX ----------------
        data_block = []
        for ln in lines[data_start:]:
            if not ln.strip():
                continue
            if not self._row_starts_with_float(ln):
                continue

            parts = [float(x) for x in ln.strip().split(",") if x != ""]
            data_block.append(parts)

        data_block = np.array(data_block)
        wavelengths = data_block[:, 0]
        intens = data_block[:, 1:]

        if n_frames is None:
            n_frames = intens.shape[1]

        # ---------------- CORRECT TIME AXIS ----------------
        # start_offset = AccumulateCycleTime * NumAccum
        start_offset = accum_cycle * num_accum  # constant shift

        # frame increment = Kinetic Cycle Time
        rel_times = start_offset + kinetic_cycle * np.arange(n_frames)

        # normalize to zero (ixdat requires time start at 0)
        rel_times = rel_times - rel_times[0]

        # ---------------- BUILD ixdat OBJECTS ----------------
        xseries = DataSeries("wavelength", "nm", wavelengths)

        tseries = TimeSeries(
            name="time",
            unit_name="s",
            data=rel_times,
            tstamp=tstamp_first,
        )

        field = Field(
            name="intensity",
            unit_name="counts",
            data=intens.T,  # (time, wavelength)
            axes_series=[tseries, xseries],
        )

        spectra = cls(
            name=name,
            reader=self,
            technique="Optical",
            tstamp=tstamp_first,
            field=field,
            continuous=True,
        )
        return spectra

    # ===================================================
    # Helpers
    # ===================================================
    @staticmethod
    def _row_starts_with_float(line):
        try:
            float(line.strip().split(",")[0])
            return True
        except Exception:
            return False

    @staticmethod
    def _get_float(lines, key):
        for ln in lines:
            if key.lower() in ln.lower():
                m = re.search(r"([-+]?\d*\.\d+|\d+)", ln)
                if m:
                    return float(m.group(0))
        return None

    @staticmethod
    def _get_int(lines, key):
        for ln in lines:
            if key.lower() in ln.lower():
                m = re.search(r"(\d+)", ln)
                if m:
                    return int(m.group(1))
        return None

    @staticmethod
    def _parse_datetime_header(lines):
        for ln in lines[:50]:
            if ln.lower().startswith("date and time"):
                s = ln.split(":", 1)[1].strip()

                # REMOVE ALL TRAILING NON-DATE GARBAGE (commas, spaces)
                s = re.sub(r"[, ]+$", "", s)

                # Normalize day (1 → 01)
                s = re.sub(
                    r"(^[A-Za-z]{3}\s+[A-Za-z]{3}\s+)(\d)(\s+)",
                    r"\g<1>0\2\3",
                    s,
                )

                fmts = [
                    "%a %b %d %H:%M:%S.%f %Y",
                    "%a %b %d %H:%M:%S %Y",
                ]
                for fmt in fmts:
                    try:
                        return datetime.strptime(s, fmt)
                    except Exception:
                        continue
        return None
