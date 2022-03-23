"""Fixtures used across the functional tests"""

from pathlib import Path
import numpy as np
from tempfile import TemporaryDirectory

from ixdat import Measurement
from ixdat.db import change_database

temporary_directory = TemporaryDirectory()
# Set the directory backup up to work of a temporary directory
change_database(
    "directory",
    directory=Path(temporary_directory.name),
    project_name="test_biologic_ec_measurement",
)

PATH_TO_DATAFILE = Path(__file__).parent / "../test_data/biologic/Pt_poly_cv_CUT.mpt"

ec_measurement = Measurement.read(PATH_TO_DATAFILE, reader="biologic")

measurement1 = ec_measurement.select(cycle=1)
measurement2 = ec_measurement.select(cycle=3)

composed_measurement = measurement1 + measurement2

ec_measurement.calibrate_RE(RE_vs_RHE=1)

assert np.isclose(
    ec_measurement.U[0] - ec_measurement["raw_potential"].data[0],
    ec_measurement.RE_vs_RHE,
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

cv_copy = cv.copy()

assert cv_copy == cv

"""Test a load/save round trip of a composed measurement"""
id_ = composed_measurement.save()
loaded = Measurement.get(id_)

assert composed_measurement == loaded
