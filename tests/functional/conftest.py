"""Fixtures used across the functional tests"""

from pathlib import Path
from tempfile import TemporaryDirectory

from pytest import fixture

from ixdat import Measurement
from ixdat.db import DB, change_database

# FIXME The size of this data file is at present one of the largest contributors to the
# time it takes to run the test suite. We should consider cutting it down or replacing
# it for test purposes
PATH_TO_DATAFILE = Path(__file__).parent / "../../test_data/biologic/Pt_poly_cv_CUT.mpt"


@fixture(scope="function", params=["directory", "sqlite"])
def fresh_backend(request):
    """Fixture that provides each database backend in a fresh temporary directory"""
    temporary_directory = TemporaryDirectory()
    original_backend = DB.backend
    if request.param == "directory":
        change_database(
            "directory",
            directory=Path(temporary_directory.name),
            project_name="test_biologic_ec_measurement",
        )
    else:
        db_path = Path(temporary_directory.name) / "test_biologic_ec_measurement.sqlite"
        change_database("sqlite", db_path=db_path)
    backend = DB.backend
    yield temporary_directory
    if hasattr(backend, "close"):
        backend.close()
    DB.set_backend(original_backend)
    temporary_directory.cleanup()


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
