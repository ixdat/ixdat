"""Reader for Bruker TopSpin NMR experiment folders.

A Bruker TopSpin experiment is a directory containing:

* The raw FID file (binary, ``fid`` for 1D / ``ser`` for nD).
* The acquisition parameter file ``acqus`` (and ``acqu2s`` etc. for nD).
* Optional pulse-program file ``pulseprogram``.
* A processed-data subfolder ``pdata/<n>/`` with parameter files ``procs``
  (``proc2s`` ...) and processed binary spectrum files (``1r``, ``1i`` for 1D).

This reader uses :mod:`nmrglue` to parse the parameter files and the binary
data and returns an :class:`NMRSpectrum` (1D) populated with the chemical-shift
axis (in ppm), the processed real spectrum, and a metadata dictionary holding
the most relevant acquisition + processing parameters.
"""

import datetime
from pathlib import Path

import numpy as np

from ..data_series import DataSeries, Field
from ..exceptions import ReadError
from ..techniques.nmr import NMRSpectrum


# Acquisition parameters lifted from `acqus` into a flat metadata dict.
# Anything not in this list is still preserved in metadata["acqus"].
#
# Definitions follow the Bruker TopSpin Acquisition Commands and Parameters Reference Manual.
# A copy is available at:
# https://www.nmr.ucdavis.edu/sites/g/files/dgvnsk4156/files/inline-files/acqu_commands_parameters.pdf
#
# Key Bruker acqus abbreviations:
#   SOLVENT  - NMR solvent (e.g. D2O, DMSO)
#   TE       - sample temperature in Kelvin
#   BF1      - basic transmitter frequency for channel 1 in MHz (e.g. 700 for a 700 MHz instrument;
#              set automatically by TopSpin for each nucleus)
#   SFO1     - exact irradiation frequency of channel 1 in MHz: SFO1 = BF1 + O1/1e6; this is the
#              frequency at the center of the recorded spectral window
#   O1       - carrier frequency offset from BF1 in Hz; sets the center of the spectral region to
#              be acquired; adjust to place the window over the region of interest
#   SW       - total spectral width in ppm (the chemical-shift range recorded)
#   SW_h     - spectral width in Hz (the equivalent of SW in Hz units)
#   NS       - number of scans: each scan is one RF pulse + FID acquisition; signals add
#              coherently while random noise adds as sqrt(NS), so signal-to-noise ratio grows
#              as sqrt(NS); NS should be a multiple of the phase-cycle length (often 8)
#   DS       - dummy scans: loops of the pulse program executed without digitizing or accumulating
#              data, run before the NS counted scans to reach a repeatable steady state
#   TD       - time-domain size: number of complex data points digitised in the FID
#   PULPROG  - pulse program name (defines the RF pulse sequence, e.g. zg, noesypr1d)
#   AQ_mod   - acquisition mode (e.g. DQD = digital quadrature detection)
#   NUC1     - nucleus assigned to frequency channel 1 (e.g. 1H, 13C, 31P)
#   P        - array of rectangular pulse lengths in microseconds (P[1] is typically the
#              90-degree hard pulse that tips magnetisation from the z-axis into the xy-plane)
#   D        - array of inter-pulse delay times in seconds
_ACQUS_KEYS = (
    "SOLVENT",
    "TE",
    "BF1",
    "SFO1",
    "O1",
    "SW",
    "SW_h",
    "NS",
    "DS",
    "TD",
    "PULPROG",
    "AQ_mod",
    "DATE",
    "NUC1",
    "P",
    "D",
)

# Processing parameters lifted from `procs`.
#
# Definitions follow the Bruker TopSpin Processing Commands and Parameters Reference Manual.
# A copy is available at:
# https://www.nmr.ucdavis.edu/sites/g/files/dgvnsk4156/files/inline-files/proc_commands_references.pdf
#
#   SI     - size of the frequency-domain spectrum in points (often 2 times TD after zero-filling;
#             zero-filling pads the FID with zeros before Fourier transform to interpolate the axis)
#   OFFSET - ppm value of the leftmost (highest-frequency) point of the processed spectrum
#   SF     - spectrometer reference frequency in MHz used to convert Hz to ppm for axis calibration
#   SW_p   - spectral width of the processed spectrum in ppm
#   LB     - line-broadening in Hz: an exponential decay multiplied into the FID before Fourier
#             transform; positive LB reduces noise at the cost of broader (lower-resolution) peaks
#   WDW    - window/apodization function code applied to the FID before Fourier transform
#             (0 = none, 1 = exponential, 2 = Gaussian, etc.)
#   PHC0   - zero-order (constant) phase correction angle in degrees applied after Fourier transform
#   PHC1   - first-order phase correction angle in degrees; varies linearly across the spectrum to
#             correct the phase roll caused by a delayed FID start or finite pulse width
_PROCS_KEYS = ("SI", "OFFSET", "SF", "SW_p", "LB", "WDW", "PHC0", "PHC1")


class BrukerNMRReader:
    """Reader for a Bruker TopSpin NMR experiment folder."""

    def __init__(self):
        self.path_to_folder = None
        self.dic = None  # full nmrglue parameter dictionary

    def read(
        self,
        path_to_file,
        name=None,
        cls=NMRSpectrum,
        procno=1,
        prefer_processed=True,
        **kwargs,
    ):
        """Read a Bruker TopSpin experiment folder.

        Args:
            path_to_file (Path or str): Path to the experiment folder (the one
                containing ``acqus``). May also be the path to a file inside
                that folder, in which case the parent folder is used.
            name (str): Optional name override. Defaults to the folder name.
            cls (Spectrum subclass): Class to instantiate. Defaults to
                :class:`NMRSpectrum`. Anything that is not a subclass falls
                back to :class:`NMRSpectrum`.
            procno (int): Which processed-data subfolder to read from
                (``pdata/<procno>``). Defaults to 1.
            prefer_processed (bool): If True (default) and a processed
                spectrum is available, return it. If False, or if no processed
                data is found, return the magnitude of the raw FID instead.
            kwargs: Forwarded to ``cls``.
        """
        try:
            import nmrglue as ng
        except ImportError as e:
            raise ReadError(
                "To read Bruker TopSpin NMR data, ixdat uses the nmrglue "
                "package. Install it with `pip install nmrglue` and try again."
                f"\noriginal error:\n{e}"
            )

        folder = Path(path_to_file)
        if folder.is_file():
            folder = folder.parent
        if not (folder / "acqus").is_file():
            raise ReadError(
                f"No 'acqus' file found in {folder}. Point BrukerNMRReader at a "
                "Bruker TopSpin experiment folder (the one containing acqus)."
            )

        self.path_to_folder = folder
        name = name or folder.name

        if not (isinstance(cls, type) and issubclass(cls, NMRSpectrum)):
            cls = NMRSpectrum

        pdata_dir = folder / "pdata" / str(procno)
        have_processed = prefer_processed and (pdata_dir / "procs").is_file()

        if have_processed:
            dic, data = ng.bruker.read_pdata(str(pdata_dir))
            # nmrglue returns the processed real spectrum (1r) for 1D pdata.
            y = np.real(np.asarray(data))
            udic = ng.bruker.guess_udic(dic, data)
            # uc (unit-conversion object) maps data-point indices to ppm using the
            # spectrometer frequency, spectral width, and offset stored in the procs file.
            uc = ng.fileiobase.uc_from_udic(udic, dim=0)
            # ppm (parts per million): chemical shift = (freq_sample - freq_ref) / freq_ref * 1e6.
            # Expressing frequencies as ppm makes peak positions instrument-independent.
            x_ppm = uc.ppm_scale()
        else:
            dic, data = ng.bruker.read(str(folder))
            # Bruker ADCs apply a digital filter that shifts the FID by a fixed number of
            # points; removing it realigns the time-domain signal so point 0 is t = 0.
            # Without this step the Fourier-transformed baseline is severely distorted.
            data = ng.bruker.remove_digital_filter(dic, data)
            # The raw FID (Free Induction Decay) is the time-domain signal recorded after
            # an RF pulse: exponentially decaying oscillations from precessing nuclear spins.
            # It is complex (two receiver channels shifted 90 degrees apart). Without a
            # Fourier transform + phase correction + chemical-shift referencing we cannot
            # assign ppm positions, so we return the magnitude (envelope) against a plain
            # point index instead.
            y = np.abs(np.asarray(data))
            x_ppm = np.arange(y.shape[-1], dtype=float)

        self.dic = dic
        metadata = self._build_metadata(dic, have_processed=have_processed)
        tstamp = _extract_tstamp(dic)

        xseries = DataSeries(
            name="chemical shift" if have_processed else "point",
            unit_name="ppm" if have_processed else None,
            data=np.asarray(x_ppm),
        )
        field = Field(
            name="intensity",
            unit_name=None,
            data=np.asarray(y),
            axes_series=[xseries],
        )

        return cls(
            name=name,
            technique="NMR",
            reader=self,
            tstamp=tstamp,
            metadata=metadata,
            field=field,
            **kwargs,
        )

    def _build_metadata(self, dic, *, have_processed):
        """Flatten the most-used parameters into a json-friendly dict."""
        acqus = dic.get("acqus", {}) if isinstance(dic, dict) else {}
        procs = {}
        if have_processed:
            procs = dic.get("procs", {}) if isinstance(dic, dict) else {}

        metadata = {
            "source": "bruker_topspin",
            "folder": str(self.path_to_folder),
            "processed": have_processed,
        }
        for key in _ACQUS_KEYS:
            if key in acqus:
                metadata[key] = _to_jsonable(acqus[key])
        for key in _PROCS_KEYS:
            if key in procs:
                metadata[f"proc_{key}"] = _to_jsonable(procs[key])

        # Keep the full raw dictionaries too, normalised for json.
        metadata["acqus"] = _to_jsonable(acqus)
        if procs:
            metadata["procs"] = _to_jsonable(procs)
        return metadata


def _extract_tstamp(dic):
    """Best-effort unix timestamp from acqus DATE (seconds since epoch)."""
    acqus = dic.get("acqus", {}) if isinstance(dic, dict) else {}
    date = acqus.get("DATE")
    if date is None:
        return None
    try:
        return float(date)
    except (TypeError, ValueError):
        pass
    try:
        dt = datetime.datetime.fromisoformat(str(date))
        return dt.timestamp()
    except ValueError:
        return None


def _to_jsonable(obj):
    """Recursively convert numpy / nmrglue values into json-safe primitives."""
    if isinstance(obj, dict):
        return {str(k): _to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (bytes, bytearray)):
        try:
            return obj.decode("utf-8")
        except UnicodeDecodeError:
            return obj.decode("latin-1", errors="replace")
    return obj
