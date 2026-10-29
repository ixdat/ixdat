"""Tests for calculator lookup of spectral data in a SpectroMeasurement."""

import numpy as np
import pytest

from ixdat.data_series import DataSeries, Field, TimeSeries, ValueSeries
from ixdat.measurement_base import Calculator
from ixdat.spectra import SpectroMeasurement, SpectrumSeries


TSTAMP = 1756386155.0


class DoublingCalculator(Calculator):
    """Doubles the field looked up under the key with ">double" removed."""

    calculator_type = "doubler"
    available_series_names = set()
    available_suffixes = {">double"}

    def calculate_series(self, key, measurement=None):
        field = measurement[key[: -len(">double")]]
        return Field(
            name=key,
            unit_name=field.unit_name,
            data=field.data * 2,
            axes_series=field.axes_series,
        )


def make_spectro_measurement():
    """Return a SpectroMeasurement with 4 spectra of 3 wavelengths and one value."""
    tseries = TimeSeries(name="time", unit_name="s", data=np.arange(4.0), tstamp=TSTAMP)
    wavelength = DataSeries(
        name="wavelength", unit_name="nm", data=np.array([400.0, 500.0, 600.0])
    )
    field = Field(
        name="intensity",
        unit_name="counts",
        data=np.arange(12.0).reshape(4, 3),
        axes_series=[tseries, wavelength],
    )
    spectrum_series = SpectrumSeries(
        name="spectra", tstamp=TSTAMP, field=field, continuous=True
    )
    value = ValueSeries(name="value", unit_name="", data=np.arange(4.0), tseries=tseries)
    return SpectroMeasurement(
        name="test",
        tstamp=TSTAMP,
        series_list=[tseries, value],
        spectrum_series=spectrum_series,
    )


@pytest.mark.parametrize("key", ["spectra", "intensity"])
def test_spectra_field_is_found_by_lookup(key):
    meas = make_spectro_measurement()
    assert meas[key] is meas.spectra


def test_other_keys_use_the_normal_lookup():
    meas = make_spectro_measurement()
    np.testing.assert_array_equal(meas["value"].data, np.arange(4.0))


def test_calculator_claims_key_by_suffix():
    meas = make_spectro_measurement()
    meas.add_calculator(DoublingCalculator())
    doubled = meas["spectra>double"]
    np.testing.assert_array_equal(doubled.data, meas.spectra.data * 2)


def test_suffix_lookups_chain():
    meas = make_spectro_measurement()
    meas.add_calculator(DoublingCalculator())
    quadrupled = meas["spectra>double>double"]
    np.testing.assert_array_equal(quadrupled.data, meas.spectra.data * 4)


def test_suffix_result_is_cached():
    meas = make_spectro_measurement()
    meas.add_calculator(DoublingCalculator())
    assert meas["spectra>double"] is meas["spectra>double"]


def test_cut_measurement_looks_up_its_own_spectra():
    meas = make_spectro_measurement()
    meas.add_calculator(DoublingCalculator())
    _ = meas["spectra>double"]  # populate the cache of the original
    cut = meas.cut(tspan=[0.5, 2.5])
    assert cut["spectra"].data.shape[0] == cut.spectra.data.shape[0] == 2
    np.testing.assert_array_equal(cut["spectra>double"].data, cut.spectra.data * 2)


def test_resetting_the_spectrum_series_clears_the_cache():
    meas = make_spectro_measurement()
    old = meas["spectra"]
    meas.set_spectrum_series(meas.spectrum_series.cut(tspan=[0.5, 2.5]))
    assert meas["spectra"] is not old
    assert meas["spectra"] is meas.spectra
