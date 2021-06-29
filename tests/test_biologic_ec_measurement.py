"""Tests that an ECMeasurement read from test data behaves as it should"""
from pathlib import Path
import time
from ixdat import Measurement

#  If tox crashes when trying to import matplotlib, see:
#    https://github.com/ixdat/ixdat/issues/10


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
    path_to_file = Path("./test_data/biologic/Pt_poly_cv.mpt")
    meas = Measurement.read(path_to_file, reader="biologic")
    log.print(f"successfully read {meas}")

    assert meas.technique == "EC"
    assert len(meas["potential"].data) == len(meas["time/s"].data)

    meas.calibrate_RE(RE_vs_RHE=1)
    assert meas.v[0] - meas["raw_potential"].data[0] == meas.RE_vs_RHE

    cv = meas.as_cv()
    cvs_1_plus_2 = cv[1] + cv[2]

    assert cvs_1_plus_2.RE_vs_RHE == meas.RE_vs_RHE
    assert len(cvs_1_plus_2.selector.data) == len(cvs_1_plus_2.t)
    assert len(cvs_1_plus_2["potential"].data) == len(cvs_1_plus_2["time/s"].data)
    # ^ perfect! Tests that don't work now but should work after this restructure :D


if __name__ == "__main__":
    path_to_file = Path(__file__).parent.parent / "test_data/biologic/Pt_poly_cv.mpt"
    meas = Measurement.read(path_to_file, reader="biologic")
    log.print(f"successfully read {meas}")

    assert meas.technique == "EC"
    assert len(meas["potential"].data) == len(meas["time/s"].data)

    meas.calibrate_RE(RE_vs_RHE=1)
    assert meas.v[0] - meas["raw_potential"].data[0] == meas.RE_vs_RHE

    cv = meas.as_cv()
    cvs_1_plus_2 = cv[1] + cv[2]

    assert cvs_1_plus_2.RE_vs_RHE == meas.RE_vs_RHE
    assert len(cvs_1_plus_2.selector.data) == len(cvs_1_plus_2.t)
    assert len(cvs_1_plus_2["potential"].data) == len(cvs_1_plus_2["time/s"].data)
    # ^ perfect! Tests that don't work now but should work after this restructure :D