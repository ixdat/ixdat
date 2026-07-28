"""Fixtures used across the functional tests"""

from pathlib import Path

import numpy as np
from pytest import fixture

from ixdat import Measurement, Spectrum
from ixdat.data_series import DataSeries, Field, TimeSeries
from ixdat.db import DB, change_database
from ixdat.spectra import SpectrumSeries
from ixdat.techniques.spectroelectrochemistry import ECOpticalMeasurement

# FIXME The size of this data file is at present one of the largest contributors to the
# time it takes to run the test suite. We should consider cutting it down or replacing
# it for test purposes
PATH_TO_DATAFILE = Path(__file__).parent / "../../test_data/biologic/Pt_poly_cv_CUT.mpt"


@fixture(scope="function", params=["directory", "sqlite"])
def fresh_backend(request, tmp_path):
    """Fixture that provides each database backend in a fresh temporary directory"""
    original_backend = DB.backend
    if request.param == "directory":
        change_database(
            "directory",
            directory=tmp_path,
            project_name="test_biologic_ec_measurement",
        )
    else:
        change_database(
            "sqlite", db_path=tmp_path / "test_biologic_ec_measurement.sqlite"
        )
    backend = DB.backend
    try:
        yield backend
    finally:
        if hasattr(backend, "close"):
            backend.close()
        DB.set_backend(original_backend)


@fixture(scope="function")
def ec_measurement():
    """Fixture that sets up an Biologic EC measurement"""
    return Measurement.read(PATH_TO_DATAFILE, reader="biologic")


@fixture(scope="function")
def composed_measurement(ec_measurement):
    """Fixture that returns a composed measurement"""
    return ec_measurement.select(cycle=1) + ec_measurement.select(cycle=3)


@fixture(scope="function")
def ec_optical_measurement_factory():
    """Build a synthetic EC-Optical measurement with an optional reference."""

    def make_measurement(with_reference):
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
        reference_spectrum = None
        if with_reference:
            reference_spectrum = Spectrum(
                name="ref spectrum",
                technique="optical",
                tstamp=1.6e9,
                field=Field(
                    name="reference",
                    unit_name="counts",
                    data=np.ones(4),
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

    return make_measurement
