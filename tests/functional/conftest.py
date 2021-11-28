"""Fixtures used across the functional tests"""

from pathlib import Path
from tempfile import TemporaryDirectory

from pytest import fixture

from ixdat import Measurement
from ixdat.db import change_database

# FIXME The size of this data file is at present one of the largest contributors to the
# time it takes to run the test suite. We should consider cutting it down or replacing
# it for test purposes
PATH_TO_DATAFILE = Path(__file__).parent / "../../test_data/biologic/Pt_poly_cv_CUT.mpt"


@fixture(scope="function")
def fresh_directory_backend():
    """Fixture that provides a directory backend in a fresh temporary directory"""
    temporary_directory = TemporaryDirectory()
    # Set the directory backup up to work of a temporary directory
    change_database(
        "directory",
        directory=Path(temporary_directory.name),
        project_name="test_biologic_ec_measurement",
    )
    return temporary_directory


@fixture(scope="function")
def ec_measurement():
    """Fixture that sets up an Biologic EC measurement"""
    return Measurement.read(PATH_TO_DATAFILE, reader="biologic")


@fixture(scope="function")
def composed_measurement(ec_measurement):
    """Fixture that returns a composed measurement"""
    measurement1 = ec_measurement.select(cycle=1)
    measurement2 = ec_measurement.select(cycle=3)
    return measurement1 + measurement2
