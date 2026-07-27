"""Fixtures used across the functional tests"""

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from pytest import fixture

from ixdat import Measurement, Spectrum
from ixdat.calculators.ms_calculators import (
    MSBackgroundSet,
    MSCalibration,
    MSCalResult,
    MSConstantBackground,
)
from ixdat.data_series import DataSeries, Field, TimeSeries
from ixdat.db import DB, change_database
from ixdat.spectra import SpectrumSeries
from ixdat.techniques.spectroelectrochemistry import ECOpticalMeasurement

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


@fixture(scope="function")
def ms_calibration():
    """Fixture that sets up an MS calibration with two sensitivity factors"""
    return MSCalibration(
        name="test ms calibration",
        tstamp=1.6e9,
        ms_cal_results=[
            MSCalResult(name="O2 at M32", mol="O2", mass="M32", F=1.5),
            MSCalResult(name="H2 at M2", mol="H2", mass="M2", F=3.0),
        ],
    )


@fixture(scope="function")
def ms_background_set():
    """Fixture that sets up a set of constant MS backgrounds"""
    return MSBackgroundSet(
        name="test ms backgrounds",
        tstamp=1.6e9,
        bg_list=[
            MSConstantBackground(name="M32 background", mass="M32", bg=1e-12),
            MSConstantBackground(name="M2 background", mass="M2", bg=2e-12),
        ],
    )


@fixture(scope="function")
def ec_optical_measurement():
    """Fixture that sets up a synthetic EC-Optical measurement

    Unlike the other measurements here this one is built rather than read, since
    ixdat ships no EC-Optical data file. What makes it worth testing is that it
    refers to two different spectrum objects: a `SpectrumSeries` with the spectra
    and a single `Spectrum` used as the reference for optical density.
    """
    tseries = TimeSeries(
        name="t", unit_name="s", data=np.array([0.0, 1.0]), tstamp=1.6e9
    )
    wavelength = DataSeries(
        name="wavelength / nm", unit_name="nm", data=np.linspace(400, 700, 4)
    )
    spectrum_series = SpectrumSeries(
        name="optical spectra",
        technique="EC-Optical",
        tstamp=1.6e9,
        field=Field(
            name="spectra",
            unit_name="counts",
            data=np.array([[1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 4.0, 5.0]]),
            axes_series=[tseries, wavelength],
        ),
    )
    reference_spectrum = Spectrum(
        name="ref spectrum",
        technique="optical",
        tstamp=1.6e9,
        field=Field(
            name="reference",
            unit_name="counts",
            data=np.array([1.0, 1.0, 1.0, 1.0]),
            axes_series=[wavelength],
        ),
    )
    return ECOpticalMeasurement(
        name="synthetic EC-Optical",
        technique="EC-Optical",
        tstamp=1.6e9,
        series_list=[tseries],
        spectrum_series=spectrum_series,
        reference_spectrum=reference_spectrum,
    )
