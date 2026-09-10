""""Tests that an ECMeasurement read from test data behaves as it should"""

import numpy as np
import pytest

from ixdat import Measurement
from ixdat.calculators.ms_calculators import (
    MSBackgroundSet,
    MSCalibration,
    MSCalResult,
    MSConstantBackground,
)
from ixdat.measurement_base import Calculator


#  If tox crashes when trying to import matplotlib, see:
#    https://github.com/ixdat/ixdat/issues/10

# NOTE The `ec_measurement` and `fresh_backend` arguments are provided by
# shared fixtures in conftest.py in this directory. The `fresh_backend` fixture is
# parametrized, so each test here runs once per database backend.


class TestBackends:
    """Functional tests for the data backends"""

    def test_measurement_round_trips(
        self, ec_measurement, composed_measurement, fresh_backend
    ):
        """Test id/name loading and a composed-measurement round trip."""
        measurement_id = ec_measurement.save()
        assert Measurement.get(measurement_id) == ec_measurement
        assert Measurement.load(ec_measurement.name) == ec_measurement

        composed_id = composed_measurement.save()
        assert Measurement.get(composed_id) == composed_measurement

    def test_round_trip_of_ms_calibration(self, fresh_backend):
        """Test that an MS calibration's sensitivity factors survive a round trip"""
        calibration = MSCalibration(
            name="test ms calibration",
            tstamp=1.6e9,
            ms_cal_results=[
                MSCalResult(name="O2 at M32", mol="O2", mass="M32", F=1.5),
                MSCalResult(name="H2 at M2", mol="H2", mass="M2", F=3.0),
            ],
        )
        loaded = Calculator.get(calibration.save())
        assert isinstance(loaded, MSCalibration)
        assert sorted((cal.mol, cal.mass, cal.F) for cal in loaded.ms_cal_results) == [
            ("H2", "M2", 3.0),
            ("O2", "M32", 1.5),
        ]
        assert loaded.get_F("O2", "M32") == 1.5

    def test_round_trip_of_ms_background_set(self, fresh_backend):
        """Test that a set of MS backgrounds survives a round trip"""
        backgrounds = MSBackgroundSet(
            name="test ms backgrounds",
            tstamp=1.6e9,
            bg_list=[
                MSConstantBackground(name="M32 background", mass="M32", bg=1e-12),
                MSConstantBackground(name="M2 background", mass="M2", bg=2e-12),
            ],
        )
        loaded = Calculator.get(backgrounds.save())
        assert isinstance(loaded, MSBackgroundSet)
        assert sorted((bg.mass, bg.bg) for bg in loaded.bg_list) == [
            ("M2", 2e-12),
            ("M32", 1e-12),
        ]

    @pytest.mark.parametrize("with_reference", [True, False])
    def test_round_trip_of_ec_optical(
        self, with_reference, ec_optical_measurement_factory, fresh_backend
    ):
        """Test EC-Optical spectrum links with and without a reference."""
        measurement = ec_optical_measurement_factory(with_reference)
        loaded = Measurement.get(measurement.save())
        assert np.allclose(loaded.spectra.data, measurement.spectra.data)
        if with_reference:
            assert loaded.reference_spectrum.name == "ref spectrum"
            assert np.allclose(loaded.reference_spectrum.y, np.ones(4))
        else:
            assert loaded.reference_spectrum is None

    def test_force_updates_scalar_relationship(
        self, ec_optical_measurement_factory, fresh_backend
    ):
        """A forced save writes changes through a one-object relationship."""
        measurement = ec_optical_measurement_factory(with_reference=True)
        measurement_id = measurement.save()

        measurement.reference_spectrum.name = "updated reference"
        measurement.reference_spectrum.field._data = np.full(4, 2.0)

        # A normal repeated save follows ixdat's established no-update rule.
        assert measurement.save() is None
        loaded = Measurement.get(measurement_id)
        assert loaded.reference_spectrum.name == "ref spectrum"
        assert np.allclose(loaded.reference_spectrum.y, np.ones(4))

        fresh_backend.save(measurement, force=True)
        loaded = Measurement.get(measurement_id)
        assert loaded.reference_spectrum.name == "updated reference"
        assert np.allclose(loaded.reference_spectrum.y, np.full(4, 2.0))
        if fresh_backend.backend_type == "directory":
            spectrum_files = list(
                (fresh_backend.project_directory / "spectrums").glob(
                    f"{measurement.reference_spectrum.id}_*"
                    f"{fresh_backend.metadata_suffix}"
                )
            )
            assert len(spectrum_files) == 1

        replacement = ec_optical_measurement_factory(
            with_reference=True
        ).reference_spectrum
        replacement.name = "replacement reference"
        replacement.field._data = np.full(4, 3.0)
        measurement.set_reference_spectrum(spectrum=replacement)

        fresh_backend.save(measurement, force=True)
        loaded = Measurement.get(measurement_id)
        assert loaded.reference_spectrum.name == "replacement reference"
        assert np.allclose(loaded.reference_spectrum.y, np.full(4, 3.0))
