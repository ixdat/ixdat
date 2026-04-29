from pathlib import Path
import numpy as np
from ..spectra import Spectrum
from ..data_series import DataSeries, Field

# Characters that mark a line as a comment in common XRD text formats.
_COMMENT_CHARS = ("#", "!", ";", "'")

# Exact lowercased column-header tokens that unambiguously indicate Q-space.
# Written by software such as GSAS-II, FullProf, and DebyeCalculator.
_Q_TOKENS = {"q", "q(a^-1)", "q(nm^-1)", "q/a^-1", "q/nm^-1"}

# Tokens that indicate a 2-theta x axis (angle in degrees).
_TWOTHETA_TOKENS = {"2theta", "2-theta", "2th", "angle"}


def _parse_header_and_data(path_to_file):
    """Return (x_name, x_unit, y_name, y_unit, data_array) from an .xy/.xye file."""
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

    x_name, x_unit, y_name, y_unit = _infer_axis_labels(header_lines)

    # Detect delimiter from first data line
    delimiter = "," if data_lines and "," in data_lines[0] else None
    data = np.loadtxt(data_lines, dtype=float, delimiter=delimiter)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return x_name, x_unit, y_name, y_unit, data


def _infer_axis_labels(header_lines):
    """Return (x_name, x_unit, y_name, y_unit) inferred from header lines.

    Defaults to two-theta in degrees with intensity in counts when no
    column-label line is found.
    """
    x_name = "two theta"
    x_unit = "degree"
    y_name = "intensity"
    y_unit = "counts"

    for line in header_lines:
        tokens = line.replace(",", " ").split()
        if not tokens or all(_is_numeric(t) for t in tokens):
            continue
        lower = [t.lower() for t in tokens]
        first = lower[0]

        if first in _Q_TOKENS:
            x_name = tokens[0]
            x_unit = _guess_q_unit(first)
            y_unit = "a.u."
        elif first.startswith("q") and first not in _TWOTHETA_TOKENS:
            # Looser fallback for non-standard exporters (e.g. "Q[1/A]").
            x_name = tokens[0]
            x_unit = _guess_q_unit(first)
            y_unit = "a.u."
        elif any(t in _TWOTHETA_TOKENS for t in lower):
            x_name = tokens[0]
            x_unit = "degree"

        # Second token, if non-numeric, names the y axis
        if len(tokens) >= 2 and not _is_numeric(tokens[1]):
            y_name = tokens[1]

    return x_name, x_unit, y_name, y_unit


def _guess_q_unit(token):
    """Return the Q unit string from a lowercased column-header token.

    1/angstrom is the crystallography convention; 1/nm is common in
    small-angle scattering and simulation tools.
    """
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
    """Reader for .xy and .xye powder diffraction text files.

    These are simple columnar text files exported by many XRD programs
    (GSAS-II, TOPAS, FullProf, DebyeCalculator, etc.):

    - .xy  : two columns  (x, intensity)
    - .xye : three columns (x, intensity, per-point error)

    The x axis is either 2-theta (the raw detector angle in degrees) or Q
    (momentum transfer, Q = 4*pi*sin(theta)/lambda, in 1/angstrom or 1/nm).
    Q is wavelength-independent and preferred when comparing data across
    instruments or with simulations. The reader infers which is present from
    optional column-label lines in the file header, defaulting to 2-theta in
    degrees when the header gives no information.

    The y-axis unit is set to "counts" for 2-theta data and "a.u." for
    Q-space data, which is typically scaled or simulated rather than a raw
    photon/neutron count.

    The per-point intensity error from .xye files (Poisson counting
    statistics for experiment, propagated error for simulation) is stored
    in ``metadata["intensity_error"]`` for use by Rietveld refinement or
    PDF analysis software.
    """

    def read(self, path_to_file, cls=None, **kwargs):
        """Read an .xy or .xye file and return a Spectrum.

        Args:
            path_to_file (str or Path): Path to the file.
            cls (Spectrum subclass): Class to instantiate. Defaults to Spectrum.
            kwargs: Passed to cls.from_field (e.g. name, sample_name).

        Returns:
            Spectrum with x (2-theta or Q) and y (intensity) series. For .xye
            files the error column is stored in ``metadata["intensity_error"]``.
        """
        path_to_file = Path(path_to_file)
        cls = cls or Spectrum

        x_name, x_unit, y_name, y_unit, data = _parse_header_and_data(path_to_file)

        x_vec = data[:, 0]
        y_vec = data[:, 1]

        metadata = kwargs.pop("metadata", {})
        if data.shape[1] >= 3:
            metadata["intensity_error"] = data[:, 2].tolist()

        xseries = DataSeries(name=x_name, unit_name=x_unit, data=x_vec)
        field = Field(name=y_name, unit_name=y_unit, data=y_vec, axes_series=[xseries])

        if "name" not in kwargs:
            kwargs["name"] = path_to_file.stem

        return cls.from_field(field, metadata=metadata or None, **kwargs)
