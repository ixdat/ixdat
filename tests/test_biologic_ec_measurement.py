"""Tests that an ECMeasurement read from test data behaves as it should"""
from pathlib import Path
import time
from ixdat import Measurement
from ixdat.db import change_database

#  If tox crashes when trying to import matplotlib, see:
#    https://github.com/ixdat/ixdat/issues/10


path_to_file = Path(__file__).parent.parent / "test_data/biologic/Pt_poly_cv.mpt"
change_database("directory", project_name="test_biologic_ec_measurement")


class TestLog:
    def __init__(self, path_to_file="~/ixdat_test_output.txt"):
        self.path_to_file = Path(path_to_file).expanduser()
        with open(self.path_to_file, "w") as f:
            f.write("\n" + "-" * 80 + "\n")
            f.write(f"RUNNING MODULE {__file__} AT t={time.time()}\n")
            f.write("-" * 80 + "\n")

    def print(self, s):
        with open(self.path_to_file, "a") as f:
            f.write(s + "\n")


log = TestLog()
log.print(f"cwd = {Path('.').absolute().resolve()}")


def test_append_essential_series():
    meas = Measurement.read(path_to_file, reader="biologic")
    log.print(f"successfully read {meas}")

    assert meas.technique == "EC"
    assert len(meas["potential"].data) == len(meas["time/s"].data)

    # Calibrate the potential and check that it works.
    meas.calibrate_RE(RE_vs_RHE=1)
    assert meas.v[0] - meas["raw_potential"].data[0] == meas.RE_vs_RHE

    # To make it complex, we first select a couple cycles, this time by converting to
    # a cyclic voltammagram and indexing, and append them:
    cv = meas.as_cv()
    cvs_1_plus_2 = cv[1] + cv[2]

    # Check that the calibration survived all that:
    assert cvs_1_plus_2.RE_vs_RHE == meas.RE_vs_RHE
    # Check that the main time variable, corresponding to potential, wasn't corrupted:
    assert len(cvs_1_plus_2.grab("potential")[0]) == len(cvs_1_plus_2["time/s"].data)
    # Check that the selector is still available and works with the main time variable:
    assert len(cvs_1_plus_2.selector.data) == len(cvs_1_plus_2.t)
    # ^ perfect! Tests that don't work now but should work after this restructure :D


def test_save_load():
    meas = Measurement.read(path_to_file, reader="biologic")
    log.print(f"successfully read {meas}")

    # To make it complex, we first select a couple cycles, this time with select(),
    # and append them:
    meas1 = meas.select(cycle=1)
    meas2 = meas.select(cycle=3)
    meas12 = meas1 + meas2

    # Save and load the series:
    i = meas12.save()
    meas12_loaded = Measurement.get(i)

    # First, check that we have the full time vector and no more:
    assert len(meas12.t) == len(meas12_loaded.t)

    # Now let's try calibrating just the loaded one:
    RE_vs_RHE = 1
    A_el = 2
    meas12_loaded.calibrate(RE_vs_RHE=RE_vs_RHE, A_el=A_el)

    # First we check that its v attribute is calibrated by comparing to the grabbed
    #   uncalibrated potential:
    assert meas12_loaded.v[0] == meas12.grab("raw_potential")[1][0] + RE_vs_RHE
    # ^ note this only builds the potential for cv12_loaded
    # First we check that its j attribute is calibrated by comparing to the same property
    #   of the original, still uncalibrated, measurement:
    assert meas12_loaded.j[0] == meas12.j[0] / A_el
    # ^ note this builds the current for both cv12 and cv12_loaded

    # As a result: there should be more cached series in cv12_loaded.
    # However, the number of series in the series_list should be conserved:
    assert len(meas12.series_list) == len(meas12_loaded.series_list)

    # Now, try copying the calibrated meaurement by as_dict() and from_dict():
    meas12_copied = Measurement.from_dict(meas12_loaded.as_dict())
    # And check if it still has the calibration:
    assert meas12_copied.A_el == A_el
    # And that it can still apply it:
    assert (
        meas12_copied.grab("potential")[1][0]
        == meas12.grab("potential")[1][0] + RE_vs_RHE
    )

    # Finally, save and load the copied measurement and check that it's still there:
    i2 = meas12_copied.save()
    meas12_copied_loaded = Measurement.get(i2)
    assert meas12_copied_loaded.RE_vs_RHE == RE_vs_RHE
    assert len(meas12_copied_loaded.v) == len(meas12.v)


if __name__ == "__main__":
    meas = Measurement.read(path_to_file, reader="biologic")
    log.print(f"successfully read {meas}")

    # To make it complex, we first select a couple cycles, this time with select(),
    # and append them:
    meas1 = meas.select(cycle=1)
    meas2 = meas.select(cycle=3)
    meas12 = meas1 + meas2

    # Save and load the series:
    i = meas12.save()
    meas12_loaded = Measurement.get(i)

    # First, check that we have the full time vector and no more:
    assert len(meas12.t) == len(meas12_loaded.t)

    # Now let's try calibrating just the loaded one:
    RE_vs_RHE = 1
    A_el = 2
    meas12_loaded.calibrate(RE_vs_RHE=RE_vs_RHE, A_el=A_el)

    # First we check that its v attribute is calibrated by comparing to the grabbed
    #   uncalibrated potential:
    assert meas12_loaded.v[0] == meas12.grab("raw_potential")[1][0] + RE_vs_RHE
    # ^ note this only builds the potential for cv12_loaded
    # First we check that its j attribute is calibrated by comparing to the same property
    #   of the original, still uncalibrated, measurement:
    assert meas12_loaded.j[0] == meas12.j[0] / A_el
    # ^ note this builds the current for both cv12 and cv12_loaded

    # As a result: there should be more cached series in cv12_loaded.
    # However, the number of series in the series_list should be conserved:
    assert len(meas12.series_list) == len(meas12_loaded.series_list)

    # Now, try copying the calibrated meaurement by as_dict() and from_dict():
    meas12_copied = Measurement.from_dict(meas12_loaded.as_dict())
    # And check if it still has the calibration:
    assert meas12_copied.A_el == A_el
    # And that it can still apply it:
    assert (
        meas12_copied.grab("potential")[1][0]
        == meas12.grab("potential")[1][0] + RE_vs_RHE
    )

    # Finally, save and load the copied measurement and check that it's still there:
    i2 = meas12_copied.save()
    meas12_copied_loaded = Measurement.get(i2)
    assert meas12_copied_loaded.RE_vs_RHE == RE_vs_RHE
    assert len(meas12_copied_loaded.v) == len(meas12.v)
