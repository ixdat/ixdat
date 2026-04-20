"""Reader for .xy and .xye powder diffraction files"""

from pathlib import Path
import numpy as np
from ..spectra import Spectrum
from ..data_series import DataSeries, Field

# Known header prefixes used by common XRD software
_COMMENT_CHARS = ("#", "!", ";", "'")

# Heuristics for x-axis label/unit from header tokens
_Q_TOKENS = {"q", "q(a^-1)", "q(nm^-1)", "q/a^-1", "q/nm^-1"}
_TWOTHETA_TOKENS = {"2theta", "2-theta", "2th", "angle"}


def _parse_header_and_data(path_to_file):
    """Return (x_name, x_unit, y_name, data_array) parsed from an .xy/.xye file.

    Header lines starting with a comment character, or bare non-numeric first lines
    (e.g. ``Q,I(Q)`` from Debye calculator output), are scanned for column labels.
    Comma-separated and whitespace-separated data are both supported.
    """
    header_lines = []
    data_lines = []
    with open(path_to_file, "r") as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped[0] in _COMMENT_CHARS:
                header_lines.append(stripped.lstrip("".join(_COMMENT_CHARS)).strip())
            elif not data_lines and not _is_numeric(
                stripped.replace(",", " ").split()[0]
            ):
                # Bare column-label line before any data (e.g. "Q,I(Q)")
                header_lines.append(stripped)
            else:
                data_lines.append(stripped)

    x_name, x_unit, y_name = _infer_axis_labels(header_lines)

    # Detect delimiter from first data line
    delimiter = "," if data_lines and "," in data_lines[0] else None
    data = np.loadtxt(data_lines, dtype=float, delimiter=delimiter)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return x_name, x_unit, y_name, data


def _infer_axis_labels(header_lines):
    """Return (x_name, x_unit, y_name) inferred from header comment lines."""
    x_name = "two theta"
    x_unit = "degree"
    y_name = "intensity"

    for line in header_lines:
        # Look for a line that looks like column headers (non-numeric tokens)
        tokens = line.replace(",", " ").split()
        if not tokens or all(_is_numeric(t) for t in tokens):
            continue
        lower = [t.lower() for t in tokens]
        # Check if first token suggests Q-space
        if lower[0] in _Q_TOKENS or lower[0].startswith("q"):
            x_name = tokens[0]
            x_unit = _guess_q_unit(lower[0])
        elif any(t in _TWOTHETA_TOKENS for t in lower):
            x_name = tokens[0]
            x_unit = "degree"
        # Second token, if present, may name the y axis
        if len(tokens) >= 2 and not _is_numeric(tokens[1]):
            y_name = tokens[1]

    return x_name, x_unit, y_name


def _guess_q_unit(token):
    if "nm" in token:
        return "1/nm"
    return "1/angstrom"


def _is_numeric(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


class XRDXYReader:
    """Reader for .xy and .xye powder diffraction files.

    These are simple columnar text files with no mandatory metadata:
    - .xy  : two columns  (x, intensity)
    - .xye : three columns (x, intensity, error)

    The x axis is assumed to be two-theta in degrees unless a comment header
    indicates Q-space (e.g. ``# Q  Intensity``).
    """

    def read(self, path_to_file, cls=None, **kwargs):
        """Read an .xy or .xye file and return a Spectrum.

        Args:
            path_to_file (str or Path): Path to the file.
            cls (Spectrum subclass): Class to instantiate. Defaults to Spectrum.
            kwargs: Passed to cls.from_field (e.g. name, sample_name).

        Returns:
            Spectrum with x (angle or Q) and y (intensity) data. For .xye files
            the error column is stored in metadata.
        """
        path_to_file = Path(path_to_file)
        cls = cls or Spectrum

        x_name, x_unit, y_name, data = _parse_header_and_data(path_to_file)

        x_vec = data[:, 0]
        y_vec = data[:, 1]

        metadata = kwargs.pop("metadata", {}) or {}
        if data.shape[1] >= 3:
            metadata["intensity_error"] = data[:, 2].tolist()

        xseries = DataSeries(name=x_name, unit_name=x_unit, data=x_vec)
        field = Field(name=y_name, unit_name="counts", data=y_vec, axes_series=[xseries])

        if "name" not in kwargs:
            kwargs["name"] = path_to_file.stem

        return cls.from_field(field, metadata=metadata or None, **kwargs)
