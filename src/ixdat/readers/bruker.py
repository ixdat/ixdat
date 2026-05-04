import datetime
from pathlib import Path

import numpy as np

from ..data_series import DataSeries, Field
from ..exceptions import ReadError
from ..techniques.nmr import FIDSpectrum, NMRSpectrum
from ..tools import to_jsonable


# Acquisition parameters lifted from `acqus` into a flat metadata dict.
# Anything not in this list is still preserved in metadata["acqus"].
#
# Definitions follow the Bruker TopSpin Acquisition Commands and Parameters
# Reference Manual. A copy is available at:
# https://www.nmr.ucdavis.edu/sites/g/files/dgvnsk4156/files/inline-files/acqu_commands_parameters.pdf  # noqa: E501
#
# Key Bruker acqus abbreviations:
#   SOLVENT  - NMR solvent (e.g. D2O, DMSO)
#   TE       - requested sample temperature in Kelvin (the set-point, not a
#              measured value)
#   BF1      - basic transmitter frequency for channel 1 in MHz (e.g. 700 MHz);
#              set automatically by TopSpin for each nucleus
#   SFO1     - exact irradiation frequency of channel 1 in MHz; SFO1 = BF1 + O1/1e6,
#              centered in the recorded spectral window
#   O1       - carrier frequency offset from BF1 in Hz; centers the spectral
#              window and determines which chemical-shift region is recorded
#   SW       - total spectral width in ppm (the chemical-shift range recorded)
#   SW_h     - spectral width in Hz; SW_h = SW * SFO1
#   NS       - number of scans; signals add coherently while noise adds only as
#              sqrt(NS), so SNR grows as sqrt(NS). Should be a multiple of the
#              phase-cycle length (often 8)
#   DS       - dummy scans: loops of the pulse program executed without digitizing
#              or accumulating data, run before the NS counted scans to reach a
#              repeatable steady state
#   TD       - number of raw data points in the FID; for quadrature acquisition
#              (the standard) TD is twice the number of complex points, because
#              real and imaginary channels are stored interleaved
#   PULPROG  - pulse program name (defines the RF pulse sequence, e.g. zg, noesypr1d)
#   AQ_mod   - acquisition mode (e.g. DQD = digital quadrature detection)
#   NUC1     - nucleus assigned to frequency channel 1 (e.g. 1H, 13C, 31P)
#   P        - array of rectangular pulse lengths in microseconds (P[1] is typically
#              the 90-degree hard pulse that tips magnetisation from z to the xy-plane)
#   D        - array of inter-pulse delay times in seconds
ACQUS_KEYS = (
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
# Definitions follow the Bruker TopSpin Processing Commands and Parameters
# Reference Manual. A copy is available at:
# https://www.nmr.ucdavis.edu/sites/g/files/dgvnsk4156/files/inline-files/proc_commands_references.pdf  # noqa: E501
#
#   SI     - size of the real processed spectrum in points; standard Bruker
#             setups use SI = TD/2 (the stored 1r/1i files each hold SI points,
#             so total stored data is 2*SI). Zero-filling pads the FID with
#             zeros before Fourier transform to interpolate the frequency axis.
#   OFFSET - ppm of the leftmost (highest-frequency) point in the processed spectrum
#   SF     - spectrometer reference frequency in MHz used to convert Hz to ppm
#   SW_p   - spectral width of the processed spectrum in ppm; typically inherited
#             from the acquisition SW unless changed during processing
#   LB     - line-broadening in Hz: an exponential decay multiplied into the FID
#             before Fourier transform; positive LB reduces noise at the cost of
#             broader (lower-resolution) peaks
#   WDW    - apodization function applied to the FID before Fourier transform;
#             Bruker values: no, em, gm, sine, qsine, trap, user, sinc, qsinc,
#             traf, trafs
#   PHC0   - zero-order phase correction angle in degrees (constant across the spectrum)
#   PHC1   - first-order phase correction angle in degrees; varies linearly across
#             the spectrum to correct the phase roll from a delayed FID start or
#             finite pulse width
PROCS_KEYS = ("SI", "OFFSET", "SF", "SW_p", "LB", "WDW", "PHC0", "PHC1")


class BrukerNMRReader:
    """Reader for a Bruker TopSpin NMR experiment folder.

    A Bruker TopSpin experiment folder contains:

    * The raw FID file (binary, ``fid`` for 1D / ``ser`` for nD).
    * The acquisition parameter file ``acqus`` (and ``acqu2s`` etc. for nD).
    * Optional pulse-program file ``pulseprogram``.
    * A processed-data subfolder ``pdata/<n>/`` with parameter files ``procs``
      (``proc2s`` ...) and binary spectrum files (``1r``, ``1i`` for 1D).

    Uses :mod:`nmrglue` to parse the parameter and data files.
    """

    def __init__(self):
        self.path_to_folder = None
        self.dic = None  # full nmrglue parameter dictionary

    def read(
        self,
        path_to_file,
        name=None,
        cls=None,
        procno=1,
        processed=True,
        **kwargs,
    ):
        """Read a Bruker TopSpin experiment folder.

        Args:
            path_to_file (Path or str): Path to the experiment folder (the one
                containing ``acqus``). May also be the path to a file inside
                that folder, in which case the parent folder is used.
            name (str): Optional name override. Defaults to the folder name.
            cls (Spectrum subclass): Class to instantiate. Defaults to
                :class:`NMRSpectrum` when ``processed=True`` and
                :class:`FIDSpectrum` when ``processed=False``.
            procno (int): Which processed-data subfolder to read from
                (``pdata/<procno>``). Defaults to 1.
            processed (bool): If True (default), read the processed spectrum
                from ``pdata/<procno>/``. If False, read the raw FID instead.
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

        default_cls = NMRSpectrum if processed else FIDSpectrum
        try:
            if not issubclass(cls, default_cls):
                cls = default_cls
        except TypeError:
            cls = default_cls

        pdata_dir = folder / "pdata" / str(procno)
        if processed and not (pdata_dir / "procs").is_file():
            raise ReadError(
                f"processed=True but no 'procs' file found in {pdata_dir}. "
                "Either run TopSpin processing first or pass processed=False to "
                "read the raw FID."
            )

        if processed:
            dic, data = ng.bruker.read_pdata(str(pdata_dir))
            # nmrglue returns the processed real spectrum (1r) for 1D pdata.
            y = np.real(np.asarray(data))
            udic = ng.bruker.guess_udic(dic, data)
            # uc (unit-conversion object) maps data-point indices to ppm using
            # the spectrometer frequency, spectral width, and offset in procs.
            uc = ng.fileiobase.uc_from_udic(udic, dim=0)
            # ppm (parts per million): shift = (f_sample - f_ref) / f_ref * 1e6.
            # Using ppm makes peak positions instrument-independent.
            x_ppm = uc.ppm_scale()
        else:
            dic, data = ng.bruker.read(str(folder))
            # Bruker ADCs apply a digital filter that shifts the FID by a fixed number of
            # points; removing it realigns the time-domain signal so point 0 is t = 0.
            # Without this step the Fourier-transformed baseline is severely distorted.
            try:
                data = ng.bruker.remove_digital_filter(dic, data)
            except Exception:
                # Older or non-standard acqus files may lack the DECIM/DSPFVS/GRPDLY
                # parameters needed to determine the filter delay; proceed with the raw
                # (unshifted) FID in that case.
                pass
            # The raw FID (Free Induction Decay) is the time-domain signal: the
            # exponentially decaying oscillations from precessing nuclear spins,
            # recorded after an RF pulse. It is complex (two receiver channels
            # shifted 90 degrees apart). Without Fourier transform + phase
            # correction + chemical-shift referencing we cannot assign ppm
            # positions, so we return the magnitude (envelope) against a time axis.
            y = np.abs(np.asarray(data))
            # Dwell time (time between samples) = 1 / (2 * SW_h): the factor of 2
            # comes from Nyquist — the ADC samples at twice the spectral width.
            sw_h = dic["acqus"]["SW_h"]
            x_ppm = np.arange(y.shape[-1]) / (2.0 * sw_h)

        self.dic = dic
        metadata = self._build_metadata(dic, processed=processed)
        tstamp = _extract_tstamp(dic)

        xseries = DataSeries(
            name="chemical shift" if processed else "time / s",
            unit_name="ppm" if processed else "s",
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

    def _build_metadata(self, dic, *, processed):
        """Flatten the most-used parameters into a json-friendly dict."""
        acqus = dic.get("acqus", {}) if isinstance(dic, dict) else {}
        procs = {}
        if processed:
            procs = dic.get("procs", {}) if isinstance(dic, dict) else {}

        metadata = {
            "source": "bruker_topspin",
            "folder": str(self.path_to_folder),
            "processed": processed,
        }
        for key in ACQUS_KEYS:
            if key in acqus:
                metadata[key] = to_jsonable(acqus[key])
        for key in PROCS_KEYS:
            if key in procs:
                metadata[f"proc_{key}"] = to_jsonable(procs[key])

        # Keep the full raw dictionaries too, normalised for json.
        metadata["acqus"] = to_jsonable(acqus)
        if procs:
            metadata["procs"] = to_jsonable(procs)
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
