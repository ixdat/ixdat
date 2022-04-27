"""Tests that an ECMeasurement read from test data behaves as it should"""

from pytest import approx

from ixdat import Measurement


#  If tox crashes when trying to import matplotlib, see:
#    https://github.com/ixdat/ixdat/issues/10

# NOTE The `ec_measurement` and `fresh_directory_backend` arguments are provided by
# shared fixtures in conftest.py in this directory


class TestBiologicEC:
    """This class tests the Biologic EC measurement"""

    def test_basic_data(self, ec_measurement):
        """Test basic data properties"""
        assert ec_measurement.technique == "EC"
        assert len(ec_measurement["potential"].data) == len(
            ec_measurement["time/s"].data
        )

    def test_calibrate_and_append(self, ec_measurement):
        """Test that measurement ms_calibration works"""
        ec_measurement.calibrate_RE(RE_vs_RHE=1)
        assert ec_measurement.U[0] - ec_measurement["raw_potential"].data[0] == approx(
            ec_measurement.RE_vs_RHE
        )

        # To make it complex, we first select a couple cycles, this time by converting
        # to a cyclic voltammagram and indexing, and append them:
        cv = ec_measurement.as_cv()
        cvs_1_plus_2 = cv[1] + cv[2]

        # Check that the ms_calibration survived all that:
        assert cvs_1_plus_2.RE_vs_RHE == ec_measurement.RE_vs_RHE
        # Check that the main time variable, that of potential, wasn't corrupted:
        assert len(cvs_1_plus_2.grab("potential")[0]) == len(cvs_1_plus_2["time/s"].data)
        # Check that the selector is still available and works with the main time var:
        assert len(cvs_1_plus_2.selector.data) == len(cvs_1_plus_2.t)

    def test_calibration_over_save_load(self, composed_measurement):
        """Test save and load functionality"""
        # Now let's try calibrating just the loaded one:
        composed_measurement_copy = Measurement.from_dict(composed_measurement.as_dict())

        RE_vs_RHE = 1
        A_el = 2
        composed_measurement_copy.calibrate(RE_vs_RHE=RE_vs_RHE, A_el=A_el)

        # First we check that its U attribute is calibrated by comparing to the grabbed
        # uncalibrated potential:
        assert (
            composed_measurement_copy.U[0]
            == composed_measurement.grab("raw_potential")[1][0] + RE_vs_RHE
        )
        # ^ note this only builds the potential for cv12_loaded

        # First we check that its J attribute is calibrated by comparing to the same
        # property of the original, still uncalibrated, measurement:
        assert composed_measurement_copy.J[0] == composed_measurement.J[0] / A_el
        # ^ note this builds the current for both cv12 and cv12_loaded

        # As a result: there should be more cached series in cv12_loaded.
        # However, the number of series in the series_list should be conserved:
        assert len(composed_measurement.series_list) == len(
            composed_measurement_copy.series_list
        )

        # Now, try copying the calibrated measurement by as_dict() and from_dict():
        meas12_copied = Measurement.from_dict(composed_measurement_copy.as_dict())
        # And check if it still has the ms_calibration:
        assert meas12_copied.A_el == A_el
        # And that it can still apply it:
        assert meas12_copied.grab("potential")[1][0] == approx(
            composed_measurement.grab("potential")[1][0] + RE_vs_RHE
        )

        # Finally, save and load the copied measurement and check that it's still there:
        # FIXME the load save test add to the total test time, this one is a primary
        # candidate to be removed if we properly trust the backend tests
        i2 = meas12_copied.save()
        meas12_copied_loaded = Measurement.get(i2)
        assert meas12_copied_loaded.RE_vs_RHE == RE_vs_RHE
        assert len(meas12_copied_loaded.U) == len(composed_measurement.U)

    def test_to_and_as_dict(self, composed_measurement):
        """Test a to_dict/from_dict round trip yields the same result"""
        round_trip = Measurement.from_dict(composed_measurement.as_dict())
        assert round_trip == composed_measurement
