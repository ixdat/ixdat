"""Functional tests for the Bruker TopSpin NMR reader.

The test fixture is a real 1D 1H NMR experiment vendored under
``test_data/bruker/MTBLS1_ADG19007u_162_10`` (see CREDITS.txt in that
folder for attribution; data originates from MetaboLights study MTBLS1
and was redistributed via the MIT-licensed nmrML repository).

Tests skip cleanly when the optional ``nmrglue`` dependency is missing.
"""

from pathlib import Path

import numpy as np
import pytest

nmrglue = pytest.importorskip("nmrglue")

from ixdat import Spectrum  # noqa: E402
from ixdat.exceptions import ReadError  # noqa: E402
from ixdat.readers.bruker import BrukerNMRReader  # noqa: E402
from ixdat.techniques.nmr import NMRSpectrum  # noqa: E402


DATA_DIR = (
    Path(__file__).parent.parent.parent
    / "test_data"
    / "bruker"
    / "MTBLS1_ADG19007u_162_10"
)


@pytest.fixture(scope="module")
def processed_spectrum():
    return BrukerNMRReader().read(DATA_DIR)


@pytest.fixture(scope="module")
def fid_spectrum():
    return BrukerNMRReader().read(DATA_DIR, prefer_processed=False)


class TestBrukerNMRReaderProcessed:
    """Read the processed (pdata/1) 1D spectrum."""

    def test_returns_nmr_spectrum(self, processed_spectrum):
        assert isinstance(processed_spectrum, NMRSpectrum)
        assert processed_spectrum.technique == "NMR"
        assert processed_spectrum.name == DATA_DIR.name

    def test_processed_flag_set(self, processed_spectrum):
        assert processed_spectrum.metadata["processed"] is True

    def test_chemical_shift_axis(self, processed_spectrum):
        assert processed_spectrum.x_name == "chemical shift"
        assert processed_spectrum.xseries.unit_name == "ppm"
        x = processed_spectrum.x
        assert x.ndim == 1
        # 1H NMR with SI=64k -> 65536 points
        assert x.shape == (65536,)
        # ppm scale spans roughly -5 to +15 ppm for this experiment
        assert -10 < x.min() < 0
        assert 10 < x.max() < 20

    def test_intensity_data(self, processed_spectrum):
        y = processed_spectrum.y
        assert y.shape == processed_spectrum.x.shape
        assert np.isrealobj(y)
        # The processed spectrum should have a clearly positive max peak.
        assert y.max() > 0
        assert y.max() > abs(y.min())

    def test_acquisition_metadata_lifted(self, processed_spectrum):
        md = processed_spectrum.metadata
        # Lifted acquisition keys: pulse program, scans, solvent, field
        assert md["PULPROG"] == "noesypr1d"
        assert md["NS"] == 128
        assert md["SOLVENT"] == "Urine"
        # 700 MHz spectrometer, allow rounding either way
        assert 690 < float(md["BF1"]) < 710
        # 300 K (the dataset's TE)
        assert float(md["TE"]) == pytest.approx(300, abs=1)

    def test_processing_metadata_lifted(self, processed_spectrum):
        md = processed_spectrum.metadata
        assert md["proc_SI"] == 65536
        assert "proc_SF" in md
        assert "proc_OFFSET" in md

    def test_full_dictionaries_preserved(self, processed_spectrum):
        md = processed_spectrum.metadata
        assert isinstance(md["acqus"], dict)
        assert isinstance(md["procs"], dict)
        # Folder path round-tripped
        assert Path(md["folder"]) == DATA_DIR

    def test_metadata_is_json_serialisable(self, processed_spectrum):
        import json

        json.dumps(processed_spectrum.metadata)

    def test_tstamp_is_unix_seconds(self, processed_spectrum):
        ts = processed_spectrum.tstamp
        assert ts is not None
        # Some sensible bracket: between 1990 and 2100
        assert 631152000 < ts < 4102444800


class TestBrukerNMRReaderFidFallback:
    """Read the raw FID when no processed data is requested."""

    def test_fid_path_marked_unprocessed(self, fid_spectrum):
        assert fid_spectrum.metadata["processed"] is False

    def test_fid_x_axis_is_index(self, fid_spectrum):
        # Without an FT/reference we don't expose a ppm axis.
        assert fid_spectrum.x_name == "point"
        # DataSeries may store None as the empty string; accept either.
        assert not fid_spectrum.xseries.unit_name

    def test_fid_data_real_and_nonempty(self, fid_spectrum):
        y = fid_spectrum.y
        assert y.size > 0
        assert np.isrealobj(y)
        # Magnitude is non-negative
        assert (y >= 0).all()

    def test_fid_metadata_still_has_pulprog(self, fid_spectrum):
        assert fid_spectrum.metadata["PULPROG"] == "noesypr1d"


class TestBrukerNMRReaderInterface:
    """Reader integration with the high-level Spectrum.read() entrypoint."""

    def test_spectrum_read_dispatch(self):
        spec = Spectrum.read(DATA_DIR, reader="bruker")
        assert isinstance(spec, NMRSpectrum)
        assert spec.metadata["PULPROG"] == "noesypr1d"

    def test_spectrum_read_through_nmrspectrum_subclass(self):
        spec = NMRSpectrum.read(DATA_DIR, reader="bruker")
        assert isinstance(spec, NMRSpectrum)

    def test_read_accepts_path_inside_folder(self):
        spec = BrukerNMRReader().read(DATA_DIR / "acqus")
        assert isinstance(spec, NMRSpectrum)
        assert spec.metadata["PULPROG"] == "noesypr1d"

    def test_read_missing_folder_raises(self, tmp_path):
        with pytest.raises(ReadError, match="acqus"):
            BrukerNMRReader().read(tmp_path)

    def test_custom_name(self):
        spec = BrukerNMRReader().read(DATA_DIR, name="my_nmr")
        assert spec.name == "my_nmr"
