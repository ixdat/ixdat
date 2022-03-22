""""Tests that an ECMeasurement read from test data behaves as it should"""

from ixdat import Measurement


#  If tox crashes when trying to import matplotlib, see:
#    https://github.com/ixdat/ixdat/issues/10

# NOTE The `ec_measurement` and `fresh_directory_backend` arguments are provided by
# shared fixtures in conftest.py in this directory

# NOTE Ideally, we should add more functional tests here and make it possible to
# parametrize over different database backends


class TestBackends:
    """Functional tests for the data backends"""

    def test_round_trip(self, ec_measurement, fresh_directory_backend):
        """Test load/save round trip of a measurement"""
        id_ = ec_measurement.save()
        reloaded_ec_measurement = Measurement.get(id_)
        assert ec_measurement == reloaded_ec_measurement

    def test_round_trip_of_composed(self, composed_measurement, fresh_directory_backend):
        """Test a load/save round trip of a composed measurement"""
        id_ = composed_measurement.save()
        loaded = Measurement.get(id_)
        assert composed_measurement == loaded
