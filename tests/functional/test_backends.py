""""Tests that an ECMeasurement read from test data behaves as it should"""

import numpy as np

from ixdat import Measurement
from ixdat.calculators.ms_calculators import MSBackgroundSet, MSCalibration
from ixdat.measurement_base import Calculator


#  If tox crashes when trying to import matplotlib, see:
#    https://github.com/ixdat/ixdat/issues/10

# NOTE The `ec_measurement` and `fresh_backend` arguments are provided by
# shared fixtures in conftest.py in this directory. The `fresh_backend` fixture is
# parametrized, so each test here runs once per database backend.


class TestBackends:
    """Functional tests for the data backends"""

    def test_round_trip(self, ec_measurement, fresh_backend):
        """Test load/save round trip of a measurement"""
        id_ = ec_measurement.save()
        reloaded_ec_measurement = Measurement.get(id_)
        assert ec_measurement == reloaded_ec_measurement

    def test_round_trip_of_composed(self, composed_measurement, fresh_backend):
        """Test a load/save round trip of a composed measurement"""
        id_ = composed_measurement.save()
        loaded = Measurement.get(id_)
        assert composed_measurement == loaded

    def test_load_by_name(self, ec_measurement, fresh_backend):
        """Test that a saved measurement can be loaded by its name"""
        ec_measurement.save()
        loaded = Measurement.load(ec_measurement.name)
        assert ec_measurement == loaded

    def test_round_trip_of_ms_calibration(self, ms_calibration, fresh_backend):
        """Test that an MS calibration's sensitivity factors survive a round trip"""
        loaded = Calculator.get(ms_calibration.save())
        assert isinstance(loaded, MSCalibration)
        assert sorted((cal.mol, cal.mass, cal.F) for cal in loaded.ms_cal_results) == [
            ("H2", "M2", 3.0),
            ("O2", "M32", 1.5),
        ]
        assert loaded.get_F("O2", "M32") == 1.5

    def test_round_trip_of_ms_background_set(self, ms_background_set, fresh_backend):
        """Test that a set of MS backgrounds survives a round trip"""
        loaded = Calculator.get(ms_background_set.save())
        assert isinstance(loaded, MSBackgroundSet)
        assert sorted((bg.mass, bg.bg) for bg in loaded.bg_list) == [
            ("M2", 2e-12),
            ("M32", 1e-12),
        ]

    def test_round_trip_of_ec_optical(self, ec_optical_measurement, fresh_backend):
        """Test that both of an EC-Optical measurement's spectrum references survive"""
        id_ = ec_optical_measurement.save()
        loaded = Measurement.get(id_)
        assert np.allclose(loaded.spectra.data, ec_optical_measurement.spectra.data)
        assert loaded.reference_spectrum.name == "ref spectrum"
        assert np.allclose(
            loaded.reference_spectrum.y, ec_optical_measurement.reference_spectrum.y
        )
