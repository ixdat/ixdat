""""Tests that an ECMeasurement read from test data behaves as it should"""

import numpy as np

from ixdat import Measurement


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

    def test_round_trip_of_ec_optical(self, ec_optical_measurement, fresh_backend):
        """Test that both of an EC-Optical measurement's spectrum references survive"""
        id_ = ec_optical_measurement.save()
        loaded = Measurement.get(id_)
        assert np.allclose(loaded.spectra.data, ec_optical_measurement.spectra.data)
        assert loaded.reference_spectrum.name == "ref spectrum"
        assert np.allclose(
            loaded.reference_spectrum.y, ec_optical_measurement.reference_spectrum.y
        )
