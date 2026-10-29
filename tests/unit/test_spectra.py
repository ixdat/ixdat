"""Tests for spectrum indexing and averaging."""

import numpy as np
import pytest

from ixdat.data_series import DataSeries, Field, TimeSeries
from ixdat.spectra import Spectrum, SpectrumSeries
from ixdat.techniques.ms import MSSpectrumSeries
from ixdat.techniques.spectroelectrochemistry import OpticalSpectrumSeries


def make_spectrum_series(cls=SpectrumSeries, n=4, **kwargs):
    tseries = TimeSeries("time", "s", 2 * np.arange(float(n)), tstamp=1000.0)
    field = Field(
        name="intensity",
        unit_name="counts",
        data=np.arange(float(n * 2)).reshape(n, 2),
        axes_series=[tseries, DataSeries("wavelength", "nm", np.array([400, 500]))],
    )
    return cls(name="spectra", field=field, tstamp=1000.0, **kwargs)


@pytest.mark.parametrize(
    "average_number, expected_y, expected_t",
    [
        (1, [[0, 1], [2, 3], [4, 5], [6, 7]], [0, 2, 4, 6]),
        (2, [[1, 2], [5, 6]], [1, 5]),
        (np.int64(4), [[3, 4]], [3]),
    ],
)
def test_average_spectra_includes_every_complete_group(
    average_number, expected_y, expected_t
):
    averaged = make_spectrum_series().calc_average_spectra(average_number)
    np.testing.assert_allclose(averaged.y, expected_y)
    np.testing.assert_allclose(averaged.t, expected_t)


def test_average_spectra_drops_an_incomplete_trailing_group():
    averaged = make_spectrum_series(n=5).calc_average_spectra(2)
    np.testing.assert_allclose(averaged.y, [[1, 2], [5, 6]])
    np.testing.assert_allclose(averaged.t, [1, 5])


@pytest.mark.parametrize("cls", [SpectrumSeries, OpticalSpectrumSeries])
@pytest.mark.parametrize("continuous", [False, True])
def test_average_spectra_preserves_type_and_metadata(cls, continuous):
    kwargs = {"spectra_type": "Absorption"} if cls is OpticalSpectrumSeries else {}
    reader = object()
    series = make_spectrum_series(
        cls=cls,
        technique="Optical" if cls is OpticalSpectrumSeries else "spectra",
        metadata={"instrument": "test"},
        sample_name="sample",
        continuous=continuous,
        reader=reader,
        **kwargs,
    )
    averaged = series.calc_average_spectra(2)
    assert type(averaged) is cls
    assert averaged.technique == series.technique
    assert averaged.metadata == series.metadata
    assert averaged.sample_name == "sample"
    assert averaged.reader is reader
    assert averaged.continuous is continuous
    assert averaged.field.unit_name == "counts"
    assert averaged.tseries.unit_name == "s"
    assert averaged.tstamp == 1000.0
    np.testing.assert_array_equal(averaged.x, series.x)
    if cls is OpticalSpectrumSeries:
        assert averaged.spectra_type == "Absorption"


def test_average_spectra_uses_the_supplied_field():
    series = make_spectrum_series()
    series.durations = [1, 2, 3, 4]
    tseries = TimeSeries("processed time", "s", np.array([10, 12, 14, 16]), 2000.0)
    field = Field(
        name="processed spectra",
        unit_name="OD",
        data=series.y * 2,
        axes_series=[tseries, series.xseries],
    )
    averaged = series.calc_average_spectra(2, field=field)
    np.testing.assert_allclose(averaged.y, [[2, 4], [10, 12]])
    np.testing.assert_allclose(averaged.t, [11, 15])
    assert averaged.tstamp == 2000.0
    assert averaged.tseries.name == "processed time"
    assert averaged.field.unit_name == "OD"
    assert averaged.durations is None
    np.testing.assert_allclose(series.y, [[0, 1], [2, 3], [4, 5], [6, 7]])


def test_average_spectra_sums_recorded_durations():
    series = make_spectrum_series(n=5, durations=[1, 2, 3, 4, 100])
    averaged = series.calc_average_spectra(2)
    np.testing.assert_allclose(averaged.durations, [3, 7])


def test_average_spectra_preserves_masked_values():
    series = make_spectrum_series()
    field = Field(
        name="masked spectra",
        unit_name="counts",
        data=np.ma.array(
            series.y, mask=[[True, False], [True, False], [False, False], [False, False]]
        ),
        axes_series=series.field.axes_series,
    )
    averaged = series.calc_average_spectra(2, field=field)
    assert np.ma.isMaskedArray(averaged.y)
    np.testing.assert_array_equal(
        np.ma.getmaskarray(averaged.y), [[True, False], [False, False]]
    )
    assert averaged.y[0, 1] == 2
    np.testing.assert_allclose(averaged.y[1], [5, 6])


@pytest.mark.parametrize("average_number", [0, -1, 1.5, "2", 5])
def test_average_spectra_rejects_invalid_group_sizes(average_number):
    with pytest.raises(ValueError, match="average_number"):
        make_spectrum_series().calc_average_spectra(average_number)


def test_indexing_ms_spectrum_series_returns_a_spectrum():
    """Indexing a series whose technique is registered to a SpectrumSeries class
    returns a single Spectrum."""
    x = np.linspace(1, 50, 50)
    spectrum_list = [
        Spectrum.from_data(x, x * i, tstamp=1e9 + i, x_name="m/z", y_name="signal")
        for i in range(3)
    ]
    series = MSSpectrumSeries.from_spectrum_list(spectrum_list, technique="MS_spectra")

    spectrum = series[1]

    assert not isinstance(spectrum, SpectrumSeries)
    assert np.array_equal(spectrum.y, x * 1)
